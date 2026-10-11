# Research notes

Physical mechanisms and numeric sources for the waterbrush engine.

## Source status (read first)

Primary source: C. J. Curtis, S. E. Anderson, J. E. Seims, K. W. Fleischer,
D. H. Salesin, "Computer-Generated Watercolor", SIGGRAPH 97, pp. 421-430.
Below it is cited as **[CGW]** with section numbers.

**The PDF could not be opened in round 1.** The cloud environment's network
policy refused both `grail.cs.washington.edu` and `www.cs.princeton.edu`
(proxy answered 403 to CONNECT). Everything attributed to [CGW] in this file
was written from recall of the paper and is tagged:

- `[CGW §x, recalled]`: content I am confident about, but which has not been
  checked against the PDF text in this repository's history.
- `[CGW §x, recalled, uncertain]`: I remember the mechanism but not the exact
  form or number. Must be checked before code depends on it.
- `估计 (estimate)`: no source; a value chosen by us.

M0 stays open until each `recalled` tag has been checked against the PDF and
either confirmed (tag changed to `[CGW §x, p.N]`) or corrected. See STATE.md.

Other sources referred to:

- **[KM31]** P. Kubelka, F. Munk, "Ein Beitrag zur Optik der Farbanstriche",
  Z. Tech. Physik 12, 1931. Two-flux model. Used here only through the
  closed-form equations restated in [CGW §5].
- **[Haase92]** C. S. Haase, G. W. Meyer, "Modeling pigmented materials for
  realistic image synthesis", ACM TOG 11(4), 1992. KM for paint mixing; cited
  by [CGW] for the KM approach.
- **[Smith03]** R. Smith, *The Artist's Handbook*, DK. Practical properties of
  watercolor pigments (granulation, staining). Used for qualitative pigment
  identity only.
- **[CRC]** CRC Handbook of Chemistry and Physics, water properties table.

## 1. Model overview [CGW §3-4, recalled]

A painting is an ordered list of **glazes**. Each glaze is simulated
independently on a grid and then composited optically with KM ([CGW §5]).
Inside one glaze there are three layers:

| Layer | State variables | Meaning |
|---|---|---|
| Shallow-water layer | `u, v` (velocity, staggered MAC grid), `p` (pressure), `g_k` (suspended pigment k), `M` (wet-area mask) | Water and pigment moving over the paper surface. Water only moves where `M = 1`. |
| Pigment-deposition layer | `d_k` (deposited pigment k) | Pigment that has settled onto the paper and can be picked up again. |
| Capillary layer | `s` (water saturation), `c` (fluid capacity) | Water absorbed into the paper fibres, diffusing through them. |

Paper is a height field `h ∈ [0,1]` per cell. The slope `∇h` pushes surface
water downhill; `h` also sets capillary capacity and drives granulation.

Main loop per time step [CGW §4.2, recalled]:

```
for each time step:
    MoveWater(u, v, p, M)     # 4.3
    MovePigment(u, v, g, M)   # 4.4
    TransferPigment(g, d, h, M)  # 4.5
    CapillaryFlow(s, c, p, M)    # 4.6
```

**Grid for this project (estimate).** Canvas is 1000 units wide, rendered at
1000-2000 px. Simulation cell = 1 unit by default, so 1 cell = 1 px at the
1000 px render. 估计.

## 2. Paper [CGW §4.1, recalled]

- Height field `h` from a pseudo-random process. [CGW] mention Perlin noise
  and Worley cellular textures as generators. [CGW §4.1, recalled, uncertain
  about which generator was used for which figure]
- Fluid capacity: `c = h · (c_max − c_min) + c_min`. [CGW §4.1, recalled]
  Values of `c_min, c_max`: not remembered. 估计 until checked.
- Paper slope also enters velocity update as `u ← u − ∂h/∂x`,
  `v ← v − ∂h/∂y`. [CGW §4.3.1, recalled]

Project paper parameters (M4 requires texture, sizing, weight):

| Parameter | Mapping (proposal) | Source |
|---|---|---|
| Texture: hot press / cold press / rough | Amplitude and feature size of `h`. Hot press: low amplitude, fine; rough: high amplitude, coarse (several mm features). | Qualitative description of paper surfaces, common knowledge of paper making. Numbers 估计. |
| Sizing (施胶) | Scales capillary absorption rate `α` down and raises `δ`/`σ` thresholds; heavily sized paper keeps water on the surface longer. | Qualitative; 估计 for numbers. |
| Weight (克重): 190 / 300 / 640 gsm | Scales `c_max` (thicker sheet holds more water). | 190/300/640 gsm are the standard commercial weights (90/140/300 lb). Proportional scaling 估计. |

## 3. Moving water [CGW §4.3]

### 3.1 Update velocities [CGW §4.3.1, recalled]

Shallow-water equations on a staggered grid; `u` lives at `(i+½, j)`,
`v` at `(i, j+½)`, `p` at cell centres.

```
u ← u − ∇h ;  v ← v − ∇h           # paper slope
Δt ← 1 / ceil(max(|u|, |v|))       # sub-step so nothing moves > 1 cell
for t in 0..1 step Δt:
    for each u-face (i+½, j):
        A ← u²(i,j) − u²(i+1,j) + (uv)(i+½, j−½) − (uv)(i+½, j+½)
        B ← u(i+1.5,j) + u(i−½,j) + u(i+½,j+1) + u(i+½,j−1) − 4u(i+½,j)
        u'(i+½,j) ← u(i+½,j) + Δt (A + μB + p(i,j) − p(i+1,j) − κ u(i+½,j))
    (symmetric update for v)
    u ← u' ; v ← v'
    EnforceBoundaryConditions   # velocity = 0 on faces where M = 0
```

Values at half-integer points are averages of neighbours; values at integer
points of `u` are averages of the two adjacent faces.

- `A`: advection. `B`: discrete Laplacian (viscous diffusion).
  **Sign of `μB`**: physically viscosity adds `+μ∇²u`. I recall the printed
  pseudocode writing `A − μB`; whether that is a typo or a different sign
  convention for `B` must be checked. [CGW §4.3.1, recalled, uncertain]
- Viscosity `μ = 0.1`, viscous drag `κ = 0.01`. [CGW §4.3.1, recalled,
  uncertain about exact values]. These are dimensionless per-step constants,
  not SI viscosity. For reference, water at 20 °C has dynamic viscosity
  1.0 mPa·s [CRC]; there is no direct conversion because [CGW] time and
  length units are grid units.
- Boundary condition: zero velocity on any face touching a dry cell
  (`M = 0`), so water does not leave the wet-area mask. [CGW §4.3.1, recalled]

### 3.2 Relax divergence [CGW §4.3.2, recalled]

Iterative pressure projection so the velocity field is (nearly)
divergence-free:

```
t ← 0
repeat:
    δmax ← 0
    for each cell (i,j):
        δ ← ξ (u(i+½,j) − u(i−½,j) + v(i,j+½) − v(i,j−½))
        p(i,j) += δ
        u(i+½,j) += δ ; u(i−½,j) −= δ
        v(i,j+½) += δ ; v(i,j−½) −= δ
        δmax ← max(|δ|, δmax)
    t += 1
until δmax ≤ τ or t ≥ N
```

Values: `N = 50`, `τ = 0.01`, `ξ = 0.1`. [CGW §4.3.2, recalled]

### 3.3 Flow outward: edge darkening [CGW §4.3.3, recalled]

Mechanism. Near the edge of a wet region the water layer is thin and
evaporates faster, so water flows from the interior towards the edge and
carries pigment with it; pigment accumulates in a dark ring when the area
dries (the "coffee ring" effect). [CGW §2, §4.3.3, recalled]

[CGW] model this by lowering pressure near the mask boundary:

```
M' ← GaussianBlur(M, K × K)
p ← p − η (1 − M') M
```

So pressure drops only inside the wet area, and most where the blurred mask
is far from 1, i.e. within about `K/2` cells of the edge. Values: `η = 0.05`,
`K = 10`. [CGW §4.3.3, recalled]

Evaporation. [CGW] do not run a physically separate evaporation equation
for the surface layer; drying happens by water leaving the shallow layer into
the capillary layer and by the capillary layer losing water at a constant
rate. [CGW §4.6, recalled, uncertain]. Our engine needs a clock in seconds
(M4: `wait(秒)`, `look` costs 15 s), so we add an explicit evaporation rate
from the surface and capillary layers. Value 估计; to be calibrated so a
cold-press flat wash of medium load becomes touch-dry in roughly 5-15
minutes at room conditions, which matches common painting practice (估计,
no measured source found yet).

**Project acceptance link (M1).** M1 requires edge density ≥ 1.3 × interior
after drying. This threshold comes from the project goal, not from [CGW].
[CGW] show edge darkening qualitatively (figures) and give no measured ratio.
Measured ratios for real watercolor washes: no source found yet. 估计.

## 4. Moving pigment [CGW §4.4, recalled]

Suspended pigment `g_k` is advected by the surface velocity, with the same
sub-step `Δt` as in 3.1. Upwind transfer between neighbours:

```
for each pigment k:
    g' ← g
    for each cell (i,j):
        g'(i+1,j) += max(0,  u(i+½,j) · g(i,j))
        g'(i−1,j) += max(0, −u(i−½,j) · g(i,j))
        g'(i,j+1) += max(0,  v(i,j+½) · g(i,j))
        g'(i,j−1) += max(0, −v(i,j−½) · g(i,j))
        g'(i,j)   −= (sum of the four amounts above)
    g ← g'
```

(Each amount is multiplied by `Δt` inside the sub-step loop.) Pigment mass
is conserved by construction. [CGW §4.4, recalled]

Diffusion. Pigment spreading in wet-in-wet comes from the velocity field
plus the viscous Laplacian term; [CGW] do not add a separate pigment
diffusion term in §4.4. [CGW §4.4, recalled, uncertain]. If wet-in-wet
transitions are too sharp in M2 we may add a small explicit diffusion of
`g` (Brownian/turbulent mixing); its coefficient would be 估计.

## 5. Transferring pigment: deposition, lifting, granulation, staining [CGW §4.5, recalled]

Each pigment k has three parameters:

| Symbol | Name | Effect |
|---|---|---|
| `ρ_k` | density | Heavy pigments settle faster. |
| `ω_k` | staining power | High `ω` resists being picked back up. |
| `γ_k` | granularity | How strongly paper height modulates settling and lifting. |

```
for each cell in M, for each pigment k:
    δdown ← g_k (1 − h·γ_k) ρ_k
    δup   ← d_k (1 + (h − 1)·γ_k) ρ_k / ω_k
    if d_k + δdown > 1: δdown ← max(0, 1 − d_k)
    if g_k + δup  > 1: δup   ← max(0, 1 − g_k)
    d_k += δdown − δup
    g_k += δup − δdown
```

- Granulation: with `γ > 0`, valleys (small `h`) receive more deposit and
  lose less, so pigment collects in the paper's low points. [CGW §4.5,
  recalled]
- Staining: dividing `δup` by `ω` means staining pigments, once deposited,
  hardly lift. This is the mechanism behind the M3 acceptance "staining
  pigments cannot be lifted, non-staining ones partly can". [CGW §4.5,
  recalled]
- Lifting with a wet brush (M4 `lift`) is then: add water so the cell is in
  `M`, let `δup` move pigment from `d` back to `g`, remove water and the
  suspended `g` with the brush. Mechanism from [CGW §4.5]; the brush
  pickup model is ours, 估计.

Per-pigment values: [CGW] give `ρ, ω, γ` per pigment in the pigment table of
§5 (alongside K and S). I do not recall those numbers reliably. 估计 until
checked. Qualitative ordering for our 8 pigments from [Smith03] and common
manufacturer pigment charts (granulating: ultramarine, burnt sienna, raw
umber, raw sienna; staining: Prussian blue, alizarin/madder, phthalo-type
colours; non-staining: cadmiums, earths, ivory black). Exact numbers 估计.

## 6. Capillary flow and backruns [CGW §4.6, recalled]

Water from the shallow layer soaks into the paper and diffuses through it:

```
for each cell in M:
    s ← s + max(0, min(α, c − s))          # absorption from surface layer
s' ← s
for each cell (i,j), for each 4-neighbour (k,l):
    if s(i,j) > ε and s(i,j) > s(k,l) and s(k,l) > δ:
        Δs ← max(0, min(s(i,j) − s(k,l), c(k,l) − s(k,l)) / 4)
        s'(i,j) −= Δs
        s'(k,l) += Δs
s ← s'
for each cell: if s > σ then M ← 1          # wet area grows
```

| Symbol | Meaning |
|---|---|
| `α` | absorption rate from surface to paper |
| `ε` | minimum saturation a cell needs to give water |
| `δ` | minimum saturation a cell needs to receive water (dry cells, `s ≤ δ`, receive nothing, so wet areas do not bleed into dry paper: this gives the hard "wet-on-dry" edge) |
| `σ` | saturation above which a cell joins the wet-area mask |

[CGW §4.6, recalled]. Numeric values of `α, ε, δ, σ`: not recalled. 估计.

**Backruns (回水花).** When water is added next to or into a region that is
damp but not wet (`δ < s < σ`), capillary water creeps slowly into the damp
paper. The wet-area mask expands by the `s > σ` rule, so the expanding front
is irregular and slow, and the edge-darkening outflow (3.3) keeps pushing
pigment to that moving front. The result is a pale centre with a dark,
branching ("cauliflower") margin. [CGW] obtain this by requiring `δ < σ` so
that water can move through damp cells without them becoming part of the
flowing region. [CGW §2, §4.6, recalled, uncertain about the exact
inequality] This is the mechanism M2 tests with "add water to a half-dry
area".

Evaporation of capillary water: constant rate removal from `s` per step.
[CGW §4.6, recalled, uncertain]. Our seconds-based rate: 估计 (see 3.3).

## 7. Drybrush [CGW §4.7, recalled]

A brush with little water only wets the high points of the paper: the
wet-area mask for the stroke is restricted to cells with `h` above a
threshold that rises as the brush dries. Gives broken, textured strokes.
[CGW §4.7, recalled]. Threshold mapping from brush water load: 估计.

## 8. Wet-in-wet, wet-on-dry: what produces each edge

| Situation | Mechanism | Expected edge |
|---|---|---|
| Wet-in-wet | New pigment enters cells already in `M`; velocity field and advection spread it across the existing wet region. | Soft; M2 wants transition ≥ 15 px. |
| Wet-on-dry | Mask is exactly the brush footprint; dry neighbours have `s ≤ δ` and do not receive capillary water; edge darkening concentrates pigment at the boundary. | Hard; M2 wants ≤ 3 px. |
| Water into half-dry | See backruns (§6). | Pale centre, dark scalloped margin. |

The 15 px and 3 px thresholds are from the project goal. Source in the
watercolor literature: none found. 估计 (they are plausible for a 1000 px
canvas: 15 px = 1.5 % of width).

## 9. Optics: Kubelka-Munk [CGW §5, recalled]

Each glaze stores, per cell and per pigment, a thickness `x` (taken from
deposited plus suspended pigment after drying). Each pigment has absorption
`K` and scattering `S` per RGB channel.

Single layer reflectance and transmittance [CGW §5.1, recalled; equations
from KM theory as in [Haase92]]:

```
a = 1 + K/S
b = sqrt(a² − 1)
c = a·sinh(b S x) + b·cosh(b S x)
R = sinh(b S x) / c
T = b / c
```

Mixing pigments in one layer: `K` and `S` are the concentration-weighted
sums of each pigment's `K`, `S`. [CGW §5.1, recalled]

Compositing two layers (layer 1 over layer 2) [CGW §5.2, recalled]:

```
R = R1 + T1² R2 / (1 − R1 R2)
T = T1 T2 / (1 − R1 R2)
```

Glazes are composited top-down this way and finally over the paper's
reflectance (white, `R ≈ 1`, possibly slightly tinted; 估计).

This is what M3's "glazing after drying matches KM" acceptance tests
against: two dried glazes composited with these formulas must equal the
rendered result.

Specifying a pigment from its look [CGW §5.1, recalled]: given the user's
chosen colour of a unit-thickness layer over white `R_w` and over black
`R_b`, per channel:

```
a = ½ (R_w + (R_b − R_w + 1) / R_b)
b = sqrt(a² − 1)
S = (1/b) · arccoth( (b² − (a − R_w)(a − 1)) / (b (1 − R_w)) )
K = S (a − 1)
```

This lets us define our 8 tube pigments (ultramarine, burnt sienna, raw
sienna, cadmium yellow, alizarin/madder, Prussian blue, ivory black, raw
umber) from target swatch colours if the [CGW] table does not include them.

[CGW] Figure/table of 12 sample pigments with K and S (quinacridone rose,
indian red, cadmium yellow, hooker's green, cerulean blue, burnt umber,
cadmium red, brilliant orange, hansa yellow, phthalo green, french
ultramarine, interference lilac). I recall the list of names but do not
trust my memory of the 72 numbers; they must be copied from the PDF. Of our
8 pigments only cadmium yellow, burnt umber (as a stand-in neighbour of raw
umber) and ultramarine overlap. The others will be specified via the
`R_w, R_b` route; swatch targets 估计 (to be taken from manufacturer charts
if reachable).

## 10. Time and units (project-specific)

[CGW] time is in simulation steps. The easel needs seconds (`wait(秒)`,
`look` = 15 s). Proposal: fixed `Δt_phys` seconds per outer step, chosen so
the drying time of §3.3 holds. 估计. Determinism: all randomness (paper
height, brush bristle jitter) from a single seeded `numpy.random.Generator`;
no wall-clock or thread-order dependence. Numba kernels must use fixed
iteration order so float sums reproduce bit for bit.

## Open items for M0 verification

1. Check every `recalled` tag against the PDF, record page numbers.
2. Sign of `μB` in UpdateVelocities.
3. Values of `μ, κ` (§4.3.1) and whether [CGW] give `α, ε, δ, σ`.
4. Per-pigment `ρ, ω, γ` and the K/S table.
5. Exact backrun threshold relation.
6. How [CGW] handle evaporation.

## 11. M1 model as implemented (flat wash and drying)

This section records what the engine actually does for M1 and where it
departs from §3-6 above. Code: `waterbrush/sim.py`, `waterbrush/paper.py`,
`waterbrush/pigment.py`, `waterbrush/render.py`.

### 11.1 Departure from [CGW]: thin-film flow instead of a velocity field

[CGW] move surface water with an incompressible velocity field (§3.1-3.2)
and get edge darkening by lowering pressure near the mask edge (§3.3,
FlowOutward). Their water layer has no depth, so it cannot run out, and
drying is not tied to a clock in seconds. The project needs real drying
times (`wait`, `look` costs 15 s) and a wash that dries by losing water, so
the shallow layer here carries a water depth `w` (mm) and moves by
lubrication (thin-film) flow:

```
q = −(ρ g / 3μ) · f_paper · (1 + (l_c / ℓ)²) · w³ · ∇(w + z)
∂w/∂t = −∇·q − E − A
```

- Thin-film flux `q = −(w³ / 3μ) ∇p`: standard lubrication result for a
  film on a no-slip substrate. Source: A. Oron, S. H. Davis, S. G. Bankoff,
  "Long-scale evolution of thin liquid films", Rev. Mod. Phys. 69 (1997)
  931-980, eq. 2.29 form. ρ = 998 kg/m³, μ = 1.0 mPa·s at 20 °C [CRC], so
  ρg/3μ = 3270 mm⁻¹s⁻¹.
- `z`: paper surface elevation = relief × h, plus tilt (M4).
- `f_paper = 0.1`: fraction of free-film mobility left when the film sits in
  the paper's fibre mat. 估计.
- `(1 + (l_c/ℓ)²)`: capillary pressure. Below the capillary length
  `l_c = sqrt(σ/ρg) = 2.7 mm` (σ = 72.8 mN/m [CRC]) surface tension, not
  gravity, drives levelling; at a feature of size ℓ its pressure is
  `(l_c/ℓ)²` times the hydrostatic one. We use ℓ = 1 mm, the width over which
  a pinned edge thins (估计). A full fourth-order surface-tension term would
  need an implicit solver; this is the cheap substitute.
- Pigment moves with the water at the upwind cell's concentration, as in
  [CGW §4.4]; water and pigment are both conserved by construction.
- Wet mask `M`: water only moves between two wet cells, as in [CGW
  §4.3.1]. A cell leaves `M` when its depth falls below `w_dry = 2 µm`
  (估计), and its remaining suspended pigment settles there.

What is kept from [CGW]: the three layers and their roles, the mask
boundary condition, upwind pigment advection, the transfer step with
density/staining/granulation, capillary capacity following `h`, and KM
rendering.

### 11.2 Edge darkening from evaporation at the contact line

Instead of [CGW]'s FlowOutward pressure term, M1 uses its physical cause.
Deegan et al. showed that a liquid film with a pinned contact line loses
more vapour near the edge, the flux diverging as `x^(−λ)` with
`λ = (π − 2θ)/(2π − 2θ)`, which is ½ for a small contact angle θ. Water lost
at the rim is replaced by flow from the interior, which carries suspended
pigment to the rim. Source: R. D. Deegan, O. Bakajin, T. F. Dupont, G. Huber,
S. R. Nagel, T. A. Witten, "Capillary flow as the cause of ring stains from
dried liquid drops", Nature 389 (1997) 827-829.

```
E(x) = E₀ · max(1, sqrt(L / x))      x = distance to the wet-area edge
```

- `E₀ = 1.0e−4 mm/s` (0.36 mm/h): free-surface evaporation in still indoor
  air. 估计; chosen in the range of indoor pan evaporation and so that the
  M1 wash is touch-dry in about 15-20 minutes, which matches the 5-15+
  minute practice figure in §3.3.
- `L = 5 mm`: length over which the edge enhancement acts, of the order of
  the vapour boundary layer above a sheet in still air. 估计.
- Paper with no surface water loses capillary water at `E₀` too. 估计.

### 11.3 Absorption: Lucas-Washburn, calibrated to the Cobb test

```
dA/dt = S² / (2 s)        (s = water already in the paper; stops at capacity c)
```

Lucas-Washburn uptake (absorbed depth ∝ √t). Source: E. W. Washburn,
"The dynamics of capillary flow", Phys. Rev. 17 (1921) 273. Sorptivity `S`
comes from the Cobb60 value (ISO 535: grams of water taken up per m² in
60 s; 1 g/m² = 1 µm of depth): `S = Cobb60 / √60`.

- Cobb60 = 15 + 135 (1 − sizing)² g/m². Hard-sized ≈ 15 g/m², unsized
  ≈ 150 g/m². 估计 (typical ranges quoted for sized writing papers and
  blotting papers; no watercolour-paper measurement found).
- Default sizing 0.8 (Cobb60 ≈ 20 g/m²): artists' papers are internally and
  surface sized. 估计.
- Capacity `c = 0.3 × thickness`, thickness 0.5 mm at 300 gsm, scaled with
  weight; `c` varies ±15 % with `h` as in [CGW §4.1]. 0.5 mm is the
  commonly quoted thickness of 300 gsm cotton paper; porosity 0.3 is 估计.
- Pigment carried down with absorbed water is caught by the fibres and
  added to `d`.
- Capillary wicking inside the paper: diffusion of saturation `s/c` with
  `D = 0.05 mm²/s`, only out of cells above 5 % saturation ([CGW] ε rule).
  估计. No mask growth from capillary water yet; that is M2 (backruns).

### 11.4 Settling and lifting

[CGW §4.5] form, with the per-step density ρ replaced by a Stokes settling
speed `v` so it works in seconds and depends on depth:

```
δdown = G · min(1, v Δt / w) · (1 − h γ)
δup   = d · (v Δt / w_ref) · (1 + (h − 1) γ) / ω        (only where wet)
```

- Stokes: `v = (2/9) Δρ g r² / μ`. Burnt sienna is iron oxide (hematite,
  5.26 g/cm³ [CRC]); for particle radius 0.15-0.25 µm, v ≈ 0.2-0.6 µm/s.
  We use 0.3 µm/s = 3e−4 mm/s. Particle size 估计.
- `w_ref = 0.1 mm`, `ω = 2`, `γ = 0.1` for the M1 pigment: 估计. M3 sets
  the full palette.

### 11.5 Paper relief

| Texture | Tooth size | Relief (peak to valley) |
|---|---|---|
| hot press | 0.6 mm | 0.01 mm |
| cold press | 1.5 mm | 0.05 mm |
| rough | 2.5 mm | 0.10 mm |

All 估计. Height field: two octaves of Gaussian-filtered white noise from
the seed, normalised to [0, 1].

### 11.6 Units and grid

1 canvas unit = 0.3 mm (canvas 1000 units = 300 mm, a quarter imperial
sheet is 280 × 380 mm). Default grid 0.5 cells per unit (0.6 mm cells).
Rendering evaluates KM per cell and resamples the reflectance image to the
output width. 估计.

### 11.7 M1 measurement

Edge band: cells within 1 mm of the wash boundary (distance to the nearest
dry cell, measured to the cell centre). Interior: cells 10 mm or more inside.
Ratio = mean `d` in the edge band / mean `d` in the interior. The 1.3 bar is
the project's; the 1 mm band is our reading of "边缘一圈": the dark line of a
dried wash is under a millimetre wide on cold-press paper (估计, from
looking at real washes, no measured source).
