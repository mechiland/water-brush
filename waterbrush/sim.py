"""Three-layer watercolour simulation after Curtis et al. 1997.

Layers, per grid cell:
  shallow-water layer  w (water depth, mm), M (wet mask), G (suspended pigment)
  pigment-deposition   d (pigment settled on / in the paper surface)
  capillary layer      s (water held inside the paper, mm)

Pigment amounts G and d are layer thicknesses in Kubelka-Munk units, so the
renderer can use d directly.

Where this differs from Curtis 1997 (see notes/research.md, "M1 model"):
the shallow layer carries an explicit depth w and moves by thin-film
(lubrication) flow driven by gravity, paper relief and tilt, instead of an
incompressible velocity field. Curtis's FlowOutward step, which removes water
near the edge of the wet area, is replaced by its physical cause: evaporation
is stronger near a pinned contact line (Deegan et al. 1997), and the water
lost there is replaced by flow from the interior that carries pigment along.
"""
import numpy as np

from ._kernels import thin_film_flow
from .paper import Paper

# Physical constants (mm, s).
RHO_G_OVER_3MU = 3270.0   # rho*g/(3*mu) for water at 20 C, 1/(mm*s)
CAPILLARY_LENGTH_MM = 2.7  # sqrt(sigma/(rho*g)) for water at 20 C


class Params:
    """Simulation constants. Sources and estimates are in notes/research.md."""

    def __init__(self, **kw):
        # Shallow-water flow
        self.paper_drag = 0.1        # fraction of free-film mobility left on paper (estimate)
        self.edge_scale_mm = 1.0     # length over which the film thins at a pinned edge (estimate)
        self.w_dry = 2e-3            # mm; below this a cell's surface water breaks up (estimate)
        self.lift_ref_depth_mm = 0.1 # depth at which lifting runs at the settling rate (estimate)
        # Evaporation
        self.evap_mm_s = 1.0e-4      # free-surface evaporation, indoor still air
        self.evap_boundary_mm = 5.0  # vapour boundary-layer length L in J = J0*sqrt(L/x) (estimate)
        self.paper_evap_mm_s = 1.0e-4  # evaporation from damp paper with no surface water
        # Capillary layer
        self.cap_diffusion_mm2_s = 0.05  # in-plane wicking inside the paper (estimate)
        self.cap_min_frac = 0.05         # Curtis epsilon: no wicking out of nearly dry fibres
        # Numerics
        self.dt_max = 0.5
        self.cfl = 0.2
        for k, v in kw.items():
            if not hasattr(self, k):
                raise TypeError(f"unknown parameter {k!r}")
            setattr(self, k, v)


class Sim:
    def __init__(self, width_units=1000, height_units=700, cells_per_unit=0.5,
                 mm_per_unit=0.3, paper=None, pigments=(), seed=0, params=None):
        self.paper = paper or Paper()
        self.pigments = list(pigments)
        self.params = params or Params()
        self.cells_per_unit = cells_per_unit
        self.mm_per_unit = mm_per_unit
        self.dx = mm_per_unit / cells_per_unit
        self.H = int(round(height_units * cells_per_unit))
        self.W = int(round(width_units * cells_per_unit))
        shape = (self.H, self.W)
        n = len(self.pigments)

        self.h = self.paper.height_field(shape, self.dx, seed)
        self.relief = self.paper.relief_mm * self.h          # surface elevation, mm
        cap = self.paper.capacity_mm
        # [CGW §4.1]: c = c_min + h (c_max - c_min), c_min 0.3, c_max 0.7 (research.md §2),
        # scaled so the highest cells hold the paper's full capacity.
        self.c = cap * (0.3 + 0.4 * self.h) / 0.7
        self.tilt = (0.0, 0.0)                               # slope dz/dx, dz/dy

        self.w = np.zeros(shape)
        self.M = np.zeros(shape, dtype=bool)
        self.G = np.zeros((n,) + shape)
        self.d = np.zeros((n,) + shape)
        self.s = np.zeros(shape)
        self.t = 0.0
        self._edge = None

    # ------------------------------------------------------------------ input
    def apply_wash(self, mask, water_mm, amounts):
        """Lay water and pigment on the cells in mask.

        water_mm: depth of water added. amounts: pigment thickness carried in
        that water, one value per pigment.
        """
        mask = np.asarray(mask, dtype=bool)
        self.w[mask] += water_mm
        for k, a in enumerate(amounts):
            self.G[k][mask] += a
        self.M |= mask
        self._edge = None

    # ------------------------------------------------------------------ steps
    def _surface(self):
        tx, ty = self.tilt
        if tx == 0.0 and ty == 0.0:
            return self.relief
        yy, xx = np.mgrid[0:self.H, 0:self.W]
        return self.relief + self.dx * (tx * xx + ty * yy)

    def _edge_factor(self):
        """Evaporation enhancement near the contact line.

        For a thin wet film with a pinned edge, vapour diffuses away faster at
        the rim; the flux goes as x^(-1/2) with x the distance to the contact
        line (Deegan et al. 1997, small contact angle). We use
        J/J0 = max(1, sqrt(L/x)).
        """
        if self._edge is None:
            L = self.params.evap_boundary_mm
            n = int(np.ceil(L / self.dx)) + 1
            dist = np.full(self.M.shape, n, dtype=np.int32)
            inside = self.M.copy()
            for k in range(1, n + 1):
                # Cells of `inside` that touch a non-inside cell are k cells in.
                p = np.pad(inside, 1, constant_values=False)
                core = (p[1:-1, 1:-1] & p[:-2, 1:-1] & p[2:, 1:-1] & p[1:-1, :-2] & p[1:-1, 2:]
                        & p[:-2, :-2] & p[:-2, 2:] & p[2:, :-2] & p[2:, 2:])
                dist[inside & ~core] = k
                inside = core
                if not inside.any():
                    break
            x = (dist - 0.5) * self.dx
            self._edge = np.maximum(1.0, np.sqrt(L / x))
        return self._edge

    def _mobility(self):
        p = self.params
        # Capillary pressure adds to gravity at scales below the capillary
        # length; at the edge scale ell it is (l_c/ell)^2 times stronger.
        cap = 1.0 + (CAPILLARY_LENGTH_MM / p.edge_scale_mm) ** 2
        return RHO_G_OVER_3MU * p.paper_drag * cap

    def _flow(self, dt, z):
        """Thin-film flow of the shallow layer with upwind pigment transport.

        Sub-stepped: the explicit scheme needs dt <= cfl * dx^2 / D with
        D = mobility * w^3, which is much shorter than the other processes need.
        """
        mob = self._mobility()
        t = 0.0
        while t < dt:
            wmax = float(self.w.max())
            D = mob * wmax ** 3
            h = dt - t if D <= 0 else min(dt - t, self.params.cfl * self.dx ** 2 / D)
            thin_film_flow(self.w, self.G, self.M, z, mob, self.dx, h)
            t += h

    def _evaporate_and_absorb(self, dt):
        p = self.params
        M = self.M
        # Evaporation from the free surface, stronger near the contact line.
        e = p.evap_mm_s * self._edge_factor() * dt
        self.w = np.where(M, np.maximum(self.w - e, 0.0), self.w)
        # Absorption into the paper (Lucas-Washburn): ds/dt = A^2 / (2 s),
        # slowing as the paper fills and stopping at capacity.
        A = self.paper.sorptivity
        s_eff = np.maximum(self.s, A * np.sqrt(dt))
        a = np.minimum(A * A / (2.0 * s_eff) * dt, np.maximum(self.c - self.s, 0.0))
        a = np.where(M, np.minimum(a, self.w), 0.0)
        frac = a / np.maximum(self.w, 1e-12)
        # Pigment carried down with absorbed water is caught by the fibres.
        caught = self.G * frac
        self.d += caught
        self.G -= caught
        self.w -= a
        self.s += a

    def _transfer_pigment(self, dt):
        """Curtis 1997 sec. 4.5: settling and lifting, modulated by paper height."""
        M = self.M
        h = self.h
        w = np.maximum(self.w, 1e-6)
        for k, pg in enumerate(self.pigments):
            # Settling: flux = speed * concentration, so the deposited fraction
            # per step is v*dt/w. Lifting uses the same rate at a reference
            # depth, divided by staining power (Curtis's form).
            rate_down = np.minimum(pg.settle_mm_s * dt / w, 1.0)
            rate_up = pg.settle_mm_s * dt / self.params.lift_ref_depth_mm
            down = self.G[k] * rate_down * (1.0 - h * pg.granulation)
            up = self.d[k] * rate_up * (1.0 + (h - 1.0) * pg.granulation) / pg.staining
            down = np.where(M, np.clip(down, 0.0, self.G[k]), 0.0)
            up = np.where(M, np.clip(up, 0.0, self.d[k]), 0.0)
            self.G[k] += up - down
            self.d[k] += down - up

    def _dry_cells(self):
        """Cells whose surface water is used up: remaining pigment settles."""
        gone = self.M & (self.w < self.params.w_dry)
        if gone.any():
            self.d[:, gone] += self.G[:, gone]
            self.G[:, gone] = 0.0
            self.s[gone] = np.minimum(self.s[gone] + self.w[gone], self.c[gone])
            self.w[gone] = 0.0
            self.M[gone] = False
            self._edge = None

    def _capillary(self, dt):
        p = self.params
        s, c, dx = self.s, self.c, self.dx
        D = p.cap_diffusion_mm2_s
        active = s > p.cap_min_frac * c
        sat = s / c

        def flux(sa, sb, aa, ab):
            f = -D * (sb - sa) / dx               # driven by saturation difference
            # Water only leaves cells wet enough to wick (Curtis epsilon).
            return np.where(f > 0, f * aa, f * ab)

        fx = flux(sat[:, :-1], sat[:, 1:], active[:, :-1], active[:, 1:]) * np.minimum(c[:, :-1], c[:, 1:])
        fy = flux(sat[:-1, :], sat[1:, :], active[:-1, :], active[1:, :]) * np.minimum(c[:-1, :], c[1:, :])
        ds = np.zeros_like(s)
        ds[:, :-1] -= fx
        ds[:, 1:] += fx
        ds[:-1, :] -= fy
        ds[1:, :] += fy
        s = np.clip(s + dt / dx * ds, 0.0, c)
        # Damp paper not covered by surface water dries.
        s = np.where(self.M, s, np.maximum(s - p.paper_evap_mm_s * dt, 0.0))
        self.s = s

    def _dt(self):
        p = self.params
        return min(p.dt_max, p.cfl * self.dx ** 2 / p.cap_diffusion_mm2_s)

    def step(self, dt=None):
        dt = self._dt() if dt is None else dt
        if self.M.any():
            self._flow(dt, self._surface())
            self._evaporate_and_absorb(dt)
            self._transfer_pigment(dt)
            self._dry_cells()
        self._capillary(dt)
        self.t += dt
        return dt

    def run(self, seconds):
        """Advance simulated time by `seconds`."""
        end = self.t + seconds
        while self.t < end - 1e-9:
            self.step(min(self._dt(), end - self.t))

    def is_dry(self, damp_frac=0.02):
        return (not self.M.any()) and bool(np.all(self.s <= damp_frac * self.c))

    def run_until_dry(self, max_seconds=7200, chunk=5.0):
        start = self.t
        while not self.is_dry():
            if self.t - start > max_seconds:
                raise RuntimeError("painting did not dry within max_seconds")
            self.run(chunk)
        return self.t - start
