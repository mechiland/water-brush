# Research notes

Physical mechanisms and numeric sources for the waterbrush engine.

## Sources

- **[CGW]** C. J. Curtis, S. E. Anderson, J. E. Seims, K. W. Fleischer,
  D. H. Salesin, "Computer-Generated Watercolor", SIGGRAPH 97 Proceedings,
  pp. 421-430. PDF: https://grail.cs.washington.edu/wp-content/uploads/2015/08/curtis-1997-cgw.pdf
  Citations below give section and proceedings page, e.g. [CGW §4.3.1 p.424].
  (PDF page N = proceedings page 420 + N.) Checked against the PDF in round 2.
- **[Deegan96]** R. D. Deegan et al., "Contact line deposits in an
  evaporating drop", U. Chicago preprint 1996. Cited by [CGW §4.3.3] as the
  physical cause of edge darkening (coffee-ring effect).
- **[KM31]** P. Kubelka, F. Munk, "Ein Beitrag zur Optik der Farbanstriche",
  Z. Tech. Physik 12, 1931. Used through the equations in [CGW §5].
- **[Perlin85]** K. Perlin, "An image synthesizer", SIGGRAPH 85.
  **[Worley96]** S. Worley, "A cellular texturing basis function",
  SIGGRAPH 96. The two paper-height generators [CGW §4.1] cites.
- **[Szabo74]** Z. Szabo, *Creative Watercolor Techniques*, 1974. Cited by
  [CGW §4.6] for when backruns occur.
- **估计 (estimate)**: no source; a value chosen by us. Each one says what it
  is meant to achieve so it can be calibrated later.

Tags: `[CGW ...]` means the statement or number is printed in the paper.
`[ours]` means a project decision or an inference not printed in the paper.

## 1. Model overview [CGW §3-4 p.423]

A painting is an ordered set of **glazes** (washes). Each glaze is simulated
independently and then composited optically with KM [CGW §3 p.423].
Each glaze uses three layers [CGW §4 p.423, Fig. 3]:

| Layer | State variables | Meaning |
|---|---|---|
| Shallow-water layer | `u, v` velocity (staggered grid), `p` pressure, `g_k` suspended pigment k, `M` wet-area mask (0/1) | Water and pigment flowing above the paper. Water moves only where `M = 1`. |
| Pigment-deposition layer | `d_k` deposited pigment k | Pigment adsorbed onto the paper; it can be desorbed again. |
| Capillary layer | `s` saturation, `c` fluid capacity | Water absorbed into the paper and diffused by capillary action. [CGW] use this layer **only** for backruns [CGW §4 p.423]. |

Main loop [CGW §4.2 p.424]:

```
for each time step:
    MoveWater(M, u, v, p)          # UpdateVelocities, RelaxDivergence, FlowOutward
    MovePigment(M, u, v, g1..gn)
    TransferPigment(g1..gn, d1..dn)
    SimulateCapillaryFlow(M, s)
```

Scale reference: [CGW Fig. 10 caption p.429] report a painting of 11 glazes,
2750 simulator iterations in total, at 640 × 480 px. The discretization is
resolution-dependent [CGW §4.3.1 p.424], so all constants below are per grid
cell and per step, not physical units.

**Grid for this project [ours].** Canvas 1000 units wide, render 1000-2000 px.
Because [CGW] is resolution-dependent, the simulation grid is fixed in
physical units and the render is resampled. The M1 implementation sets the
units and grid (1 unit = 0.3 mm; see §11 on the M1 branch).

## 2. Paper [CGW §4.1 p.424]

- Height field `h`, `0 < h < 1`, from pseudo-random processes
  [Perlin85, Worley96] [CGW §4.1 p.424].
- Slope `∇h` perturbs velocities (streaks along the flow) [CGW §4.3 p.424].
- Fluid capacity `c = h · (c_max − c_min) + c_min` [CGW §4.1 p.424].
  [CGW] give no values for `c_min, c_max`. 估计: `c_min = 0.3`,
  `c_max = 0.7` (so capacity is a fraction of cell volume, in line with the
  definition "fraction of volume not occupied by paper fibers" [CGW §4 p.423]).
- Sizing slows water absorption and diffusion [CGW §2.1 p.422]. The paper
  has no sizing parameter.

Project paper parameters (M4 requires texture, sizing, weight) [ours]:

| Parameter | Mapping | Source |
|---|---|---|
| Texture: hot press / cold press / rough | Amplitude and feature size of `h` (sum of Worley and Perlin octaves). Hot press: low amplitude, fine; rough: high amplitude, coarse. | Generators from [CGW §4.1]; amplitudes 估计. |
| Sizing | Scales capillary absorption `α` down and raises `δ` (dry paper resists receiving water). Matches the role of sizing in [CGW §2.1]. | Direction from [CGW §2.1]; numbers 估计. |
| Weight: 190 / 300 / 640 gsm | Scales `c_max` and the time the sheet stays damp. | 190/300/640 gsm (90/140/300 lb) are the standard commercial weights; the scaling is 估计. |

## 3. Moving water [CGW §4.3]

[CGW §4.3 p.424] lists six conditions the flow must satisfy: stay inside the
mask; surplus water flows outward; damped (no oscillating waves); perturbed
by paper texture; local changes have global effects; outward flow towards the
edges for edge darkening.

### 3.1 Shallow-water equations [CGW §4.3 eq. 1-2 p.424]

```
∂u/∂t = −(∂u²/∂x + ∂(uv)/∂y) + μ∇²u − ∂p/∂x
∂v/∂t = −(∂v²/∂y + ∂(uv)/∂x) + μ∇²v − ∂p/∂y
```

(The printed equations write the advection denominators as `∂x²`, `∂y²`; that
is a typesetting slip, the terms are first derivatives of `u²` and `uv`.)

### 3.2 UpdateVelocities [CGW §4.3.1 p.424]

Staggered grid: `u` at `(i+½, j)`, `v` at `(i, j+½)`, everything else at cell
centres. A value at a point not stored directly is the average of its two
neighbours, e.g. `u(i,j) = (u(i−½,j) + u(i+½,j))/2`.

```
(u, v) ← (u, v) − ∇h
Δt ← 1 / ceil(max_ij(|u|, |v|))
for t = 0 to 1 by Δt:
    for all cells (i, j):
        A ← u²(i,j) − u²(i+1,j) + (uv)(i+½, j−½) − (uv)(i+½, j+½)
        B ← u(i+1.5,j) + u(i−½,j) + u(i+½,j+1) + u(i+½,j−1) − 4u(i+½,j)
        u'(i+½,j) ← u(i+½,j) + Δt (A − μB + p(i,j) − p(i+1,j) − κ u(i+½,j))
        A ← v²(i,j) − v²(i,j+1) + (uv)(i−½, j+½) − (uv)(i+½, j+½)
        B ← v(i+1,j+½) + v(i−1,j+½) + v(i,j+1.5) + v(i,j−½) − 4v(i,j+½)
        v'(i,j+½) ← v(i,j+½) + Δt (A − μB + p(i,j) − p(i,j+1) − κ v(i,j+½))
    (u, v) ← (u', v')
    EnforceBoundaryConditions(M, u, v)   # zero velocity on the boundary of any cell with M = 0
```

- `μ = 0.1`, `κ = 0.01` "in all of our examples" [CGW §4 p.423].
- **Sign of the viscosity term.** The pseudocode prints `A − μB`, while
  eq. 1 has `+μ∇²u`, and `B` as written is the discrete Laplacian. With
  `−μB` the term is anti-diffusive and amplifies velocity differences.
  [ours]: implement `+μB`, which matches eq. 1 and condition 3 (damping).
  If M1/M2 show the flow too uniform we will test the printed sign and log
  the comparison.
- `Δt` is chosen so no velocity exceeds one cell per step [CGW §4.3.1].

### 3.3 RelaxDivergence [CGW §4.3.2 p.424-425]

Following Foster & Metaxas, relax `∂u/∂x + ∂v/∂y` after each step:

```
t ← 0
repeat
    (u', v') ← (u, v);  δmax ← 0
    for all cells (i, j):
        δ ← ξ (u(i+½,j) − u(i−½,j) + v(i,j+½) − v(i,j−½))
        p(i,j) ← p(i,j) + δ
        u'(i+½,j) += δ;  u'(i−½,j) −= δ
        v'(i,j+½) += δ;  v'(i,j−½) −= δ
        δmax ← max(|δ|, δmax)
    (u, v) ← (u', v');  t ← t + 1
until δmax ≤ τ or t ≥ N
```

`N = 50`, `τ = 0.01`, `ξ = 0.1` [CGW §4.3.2 p.424].

Determinism note [ours]: the loop reads `u, v` and writes `u', v'`, so it is
a Jacobi-style sweep and order-independent; vectorised NumPy gives the same
result as the loop.

### 3.4 FlowOutward: edge darkening [CGW §4.3.3 p.425]

Physical cause: in a wet-on-dry stroke the contact line is pinned by surface
tension; liquid evaporating near the boundary is replenished from the
interior, so there is outward flow that carries pigment to the edge
[Deegan96, via CGW §4.3.3].

Model: lower the pressure near the edge of the mask.

```
M' ← GaussianBlur(M, K × K)
p ← p − η (1 − M') M            (eq. 3)
```

`K = 10`, `0.01 ≤ η ≤ 0.05` [CGW §4.3.3 p.425]. Kernel σ is not given;
[ours] σ = K/3 (so the K × K kernel covers ±1.5σ). 估计.

The paper describes this as "removes at each time step an amount of water
from each cell" [CGW §4.3.3]: in [CGW] pressure `p` doubles as the water
amount in a cell.

### 3.5 Evaporation and drying

[CGW] have **no evaporation routine**. Evaporation is named as a physical
process [CGW §2.1 p.422] and as the cause of edge darkening [CGW §4.3.3],
and FlowOutward removes water near edges, but nothing dries the interior,
and absorption into the capillary layer (§6 below) does not take water out
of `p`. A glaze in [CGW] ends when the user-chosen number of iterations
runs out.

[ours] The easel needs real drying (M1 "wait until dry", M2 "half-dry",
M4 `wait(秒)` and `look` costing 15 s). We add:

- uniform evaporation from the surface water `p` at rate `e_s` per second;
- absorption of surface water into the capillary layer that does subtract
  from `p` (water conserved between layers);
- evaporation from `s` at rate `e_c` per second;
- a cell leaves the mask (`M ← 0`) when its surface water reaches zero;
  a glaze is "dry" when `p = 0` and `s < δ` everywhere.

Values 估计: chosen so a medium-load cold-press wash becomes touch-dry in
roughly 5-15 minutes of simulated time, the usual drying time painters
work with (估计; no measured source found). Calibration happens in M1.

**M1 acceptance link.** "Edge density ≥ 1.3 × interior after drying" comes
from the project goal. [CGW] show edge darkening only in figures
(Fig. 2b, 5) and give no ratio. 估计.

## 4. Moving pigment [CGW §4.4 p.425]

Pigment is advected by the water velocity, moving from each cell to its
neighbours in proportion to the outflow:

```
Δt ← 1 / ceil(max_ij(|u|, |v|))
for each pigment k:
    for t = 0 to 1 by Δt:
        g' ← g ← g^k
        for all cells (i, j):
            g'(i+1,j) += max(0,  u(i+½,j) g(i,j))
            g'(i−1,j) += max(0, −u(i−½,j) g(i,j))
            g'(i,j+1) += max(0,  v(i,j+½) g(i,j))
            g'(i,j−1) += max(0, −v(i,j−½) g(i,j))
            g'(i,j)   −= the sum of those four amounts
        g^k ← g'
```

The printed pseudocode does not multiply by `Δt` inside the sub-step loop;
with `Δt < 1` that would move too much pigment per sub-step. [ours]: scale
each amount by `Δt`. Mass is conserved either way.

[CGW] add no separate pigment diffusion term. Wet-in-wet softness in [CGW]
comes only from the velocity field. [ours] If M2's "wet-in-wet transition
≥ 15 px" fails with advection alone, add a small explicit diffusion of `g`
inside the mask (coefficient 估计) and record that it is not in [CGW].

Pigment particle sizes are 0.05-0.5 µm [CGW §2.1 p.422]; lighter pigments
stay suspended longer and travel further [CGW §2.1].

## 5. Pigment adsorption and desorption [CGW §4.5 p.425]

Per pigment k: density `ρ`, staining power `ω`, granulation `γ`.

```
for each pigment k, for all cells with M(i,j) = 1:
    δdown ← g (1 − h γ) ρ
    δup   ← d (1 + (h − 1) γ) ρ / ω
    if d + δdown > 1: δdown ← max(0, 1 − d)
    if g + δup  > 1: δup   ← max(0, 1 − g)
    d ← d + δdown − δup
    g ← g + δup − δdown
```

- **Granulation**: with large `γ`, low cells (small `h`) adsorb more and
  desorb less, so pigment collects in the paper's hollows [CGW §2.2 p.422,
  §4.5]. [CGW §2.2] note granulation is strongest when the paper is very wet.
- **Separation**: denser pigments settle first [CGW §2.2], from `ρ`.
- **Staining**: `δup` is divided by `ω`, so a high-`ω` pigment barely lifts
  once deposited [CGW §4.5]. This is the mechanism behind M3's "staining
  pigments cannot be lifted, non-staining ones partly can".
- **Glazing**: paints with low resolubility let a new wash go over a dried
  one without disturbing it [CGW §2.2 p.422]. In [CGW] each glaze is a
  separate simulation, so a dry lower glaze never re-wets. [ours] Our
  engine keeps that: once a glaze is dry its `d` becomes a frozen KM layer;
  `lift` with a wet brush re-wets only the current top glaze, and frozen
  layers below are lifted in proportion to `1/ω` (估计).
- **Lift and blot** [ours]: `lift` = add water to the region (joins `M`),
  run TransferPigment so `δup` moves `d` back into `g`, then the brush
  removes a fraction of `p` and `g`. `blot` removes most of `p` and `g`
  without adding water. Pickup fractions 估计.

Values for [CGW]'s 12 example pigments are in §8 below.

## 6. Capillary flow and backruns [CGW §4.6 p.425]

Backruns happen "only when a puddle of water spreads slowly into a region
that is drying but still damp" [CGW §4.6, citing Szabo74]. In a damp region
the only water is in the paper pores, so flow is dominated by capillary
action, not momentum. As the water spreads it pushes pigment, giving
branching shapes with strongly darkened edges [CGW §2.2 p.422].

```
for all cells with M(i,j) > 0:
    s(i,j) ← s(i,j) + max(0, min(α, c(i,j) − s(i,j)))
s' ← s
for all cells (i, j), for each 4-neighbour (k, l):
    if s(i,j) > ε and s(i,j) > s(k,l) and s(k,l) > δ:
        Δs ← max(0, min(s(i,j) − s(k,l), c(k,l) − s(k,l)) / 4)
        s'(i,j) ← s'(i,j) − Δs
        s'(k,l) ← s'(k,l) + Δs
s ← s'
for all cells: if s(i,j) > σ then M(i,j) ← 1
```

| Symbol | Meaning [CGW §4.6] | Value |
|---|---|---|
| `α` | absorption rate from the shallow-water layer | not given; 估计 0.05 per step |
| `ε` | minimum saturation before a cell can give water | not given; 估计 |
| `δ` | saturation below which a cell will not receive water | not given; 估计 |
| `σ` | saturation above which the cell joins the wet-area mask | not given; 估计 |

Irregular branching comes from cell-to-cell variation in `c`, which comes
from the paper height [CGW §4.6].

How the pieces give the three M2 cases [ours, derived from the above]:

- **Wet-on-dry, hard edge**: dry neighbours have `s ≤ δ`, so no capillary
  water reaches them and the mask stays at the brush footprint; FlowOutward
  pushes pigment to that fixed boundary.
- **Wet-in-wet, soft edge**: the paper is already inside `M`, so new pigment
  is carried by the velocity field across the existing wet area.
- **Water into a half-dry wash, backrun**: the damp area has
  `δ < s ≤ σ` and is outside `M`. Capillary water reaches it, cells cross
  `σ` one by one in an irregular order set by `c`, the mask grows slowly,
  and the moving edge with FlowOutward and advection pushes pigment ahead
  of it: pale centre, dark scalloped margin. So we need `δ < σ` and a
  damp region whose `s` lies between them. The paper does not state this
  inequality; it follows from the pseudocode.

[CGW "Future Work: Generalization" p.429] say their model treats backruns and
wet-in-wet as two separate processes, while in real watercolor they are
two ends of one continuum set by paper wetness, and that unifying them
"would be a significant improvement". Our explicit drying (§3.5) is a step
in that direction: the same capillary layer runs all the time, and the
case that results depends on how wet the paper is.

## 7. Drybrush [CGW §4.7 p.426]

A nearly dry brush at a grazing angle wets only the high points of the
paper [CGW §2.2]. Model: exclude from the wet-area mask any cell whose
height is below a user-defined threshold [CGW §4.7]. In the interactive
application the threshold is set by "the depth to which the brush is
allowed to touch the paper" [CGW §6.1 p.427].

[ours] Threshold derived from the brush's remaining water load and
pressure. Mapping 估计.

## 8. Optics: Kubelka-Munk [CGW §5 p.426-427]

Each pigment has absorption `K` and scattering `S`, three values each (RGB)
[CGW §5 p.426].

### 8.1 Pigment from appearance [CGW §5.1 p.426]

Given the user's chosen RGB colour of a unit-thickness layer over white
`R_w` and over black `R_b`, per channel, with `0 < R_b < R_w < 1`:

```
a = ½ (R_w + (R_b − R_w + 1) / R_b)
b = sqrt(a² − 1)
S = (1/b) · arccoth( (b² − (a − R_w)(a − 1)) / (b (1 − R_w)) )
K = S (a − 1)
```

Values of `K` or `S` above 1 can come out for very saturated colours; [CGW]
allow them and report no ill effects [CGW §5.1].

### 8.2 Layer reflectance and transmittance [CGW §5.2 p.426]

```
c = a sinh(b S x) + b cosh(b S x)
R = sinh(b S x) / c
T = b / c
```

Several pigments in one layer: `K` and `S` of each pigment are weighted by
its relative thickness `x_k / x`, with `x = Σ x_k` [CGW §5.2]. A pigment's
thickness is `x_k = g_k + d_k` [CGW §5.2 p.426].

### 8.3 Compositing two layers [CGW §5.2 p.426]

```
R = R1 + T1² R2 / (1 − R1 R2)
T = T1 T2 / (1 − R1 R2)
```

Repeated for each glaze; the final `R` is the pixel colour [CGW §5.2].
[ours] The bottom layer is the paper: white, `R ≈ 0.95` per channel 估计
(watercolor paper is not a perfect reflector).

This is what M3's "glazing after drying matches KM" acceptance tests: the
engine's render of two dried glazes must equal these formulas applied to
the two layers' `K, S, x`.

[CGW §5.4 p.427] list the KM assumptions watercolor violates (refractive
index mismatch at the air and paper boundaries, non-random particle
orientation for metallic paints, non-diffuse light, single wavelength,
clumping) and conclude results are plausible approximations, not accurate.

### 8.4 [CGW] pigment table [CGW Fig. 5 p.427]

| Pigment | K_r | K_g | K_b | S_r | S_g | S_b | ρ | ω | γ |
|---|---|---|---|---|---|---|---|---|---|
| Quinacridone Rose | 0.22 | 1.47 | 0.57 | 0.05 | 0.003 | 0.03 | 0.02 | 5.5 | 0.81 |
| Indian Red | 0.46 | 1.07 | 1.50 | 1.28 | 0.38 | 0.21 | 0.05 | 7.0 | 0.40 |
| Cadmium Yellow | 0.10 | 0.36 | 3.45 | 0.97 | 0.65 | 0.007 | 0.05 | 3.4 | 0.81 |
| Hookers Green | 1.62 | 0.61 | 1.64 | 0.01 | 0.012 | 0.003 | 0.09 | 1.0 | 0.41 |
| Cerulean Blue | 1.52 | 0.32 | 0.25 | 0.06 | 0.26 | 0.40 | 0.01 | 1.0 | 0.31 |
| Burnt Umber | 0.74 | 1.54 | 2.10 | 0.09 | 0.09 | 0.004 | 0.09 | 9.3 | 0.90 |
| Cadmium Red | 0.14 | 1.08 | 1.68 | 0.77 | 0.015 | 0.018 | 0.02 | 1.0 | 0.63 |
| Brilliant Orange | 0.13 | 0.81 | 3.45 | 0.005 | 0.009 | 0.007 | 0.01 | 1.0 | 0.14 |
| Hansa Yellow | 0.06 | 0.21 | 1.78 | 0.50 | 0.88 | 0.009 | 0.06 | 1.0 | 0.08 |
| Phthalo Green | 1.55 | 0.47 | 0.63 | 0.01 | 0.05 | 0.035 | 0.02 | 1.0 | 0.12 |
| French Ultramarine | 0.86 | 0.86 | 0.06 | 0.005 | 0.005 | 0.09 | 0.01 | 3.1 | 0.91 |
| Interference Lilac | 0.08 | 0.11 | 0.07 | 1.25 | 0.42 | 1.43 | 0.06 | 1.0 | 0.08 |

[CGW §5.3] say these colours came from "fairly casual observations" of real
paints over black and white, so the table is a reasonable starting point,
not a measurement.

### 8.5 Our 8 tube pigments (M3) [ours]

| Ours | Starting point |
|---|---|
| Ultramarine (群青) | [CGW] French Ultramarine row, used as is. |
| Burnt sienna (熟赭) | No [CGW] row. Derive K, S from target `R_w, R_b` via §8.1; ρ, ω, γ 估计 (granulating, non-staining). |
| Raw sienna (生赭) | Same route; 估计. |
| Cadmium yellow (镉黄) | [CGW] Cadmium Yellow row, used as is. |
| Alizarin / madder (茜草红) | Same route; staining, so high ω 估计. Nearest [CGW] row: Quinacridone Rose. |
| Prussian blue (普蓝) | Same route; staining, transparent, high ω 估计. |
| Ivory black (象牙黑) | Same route; low γ, non-staining 估计. |
| Raw umber (生褐) | Start from [CGW] Burnt Umber K, S shifted greener; ρ, γ 估计. |

M3 requires ultramarine and burnt sienna to granulate visibly on rough
paper (so high γ, as in [CGW]'s ultramarine 0.91), and staining pigments to
resist lifting. One inconsistency to watch: [CGW] give Burnt Umber
ω = 9.3, the most staining in their table, while real earth pigments are
usually classed non-staining. For our earths we will use low ω and record
that as a departure from [CGW].

## 9. Time, units and determinism [ours]

- [CGW] time is in simulation iterations. Our easel needs seconds. One
  outer step = a fixed `Δt_phys` seconds of painting time; all rates in
  §3.5 are per second and multiplied by `Δt_phys`. Value 估计, set in M1.
- All randomness (paper height, bristle jitter) comes from one seeded
  `numpy.random.Generator`. No wall-clock time, no thread-order dependence.
  Numba kernels must use a fixed loop order so float sums reproduce bit for
  bit; replay tests compare raw arrays, not only PNGs.

## 10. Where we depart from [CGW] (summary)

Written at M0 against the [CGW] model. The M1 implementation (§11 on the
M1 branch) replaces the [CGW] velocity field and FlowOutward with thin-film
flow and edge-enhanced evaporation; where §11 and this table disagree, §11
describes the engine.

| Item | [CGW] | Ours | Why |
|---|---|---|---|
| Viscosity sign | pseudocode `−μB` | `+μB` | Matches eq. 1 and the damping requirement. |
| Pigment advection | no `Δt` factor in sub-step | multiply by `Δt` | Sub-steps otherwise over-transport. |
| Evaporation | none | explicit, per second | Easel needs real drying and wait/look timing. |
| Absorption | does not remove surface water | conserves water | Needed for drying to be meaningful. |
| Capillary layer | backruns only | always on | One model for wet-in-wet through backrun. |
| Earth pigment ω | Burnt Umber 9.3 | low | Real earths are non-staining; M3 needs that. |
| Pigment diffusion | none | only if M2 needs it | Will be logged if added. |
