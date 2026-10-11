# STATE

Current milestone: **M2 Wetness and edges** (not started)

Consecutive failures: 0

## Milestones

- [x] M0 Research: mechanisms and parameters of Curtis 1997 §4-5 in `notes/research.md` (verified against the PDF, round 2)
- [x] M1 Flat wash (round 3; `tests/test_m1_flat_wash.py`, `renders/m1_flat_wash.png`)
  - [x] Dried edge band (outer 1 mm) density >= 1.3 x interior (>= 10 mm in): measured 1.57
  - [x] Hard edge from drying visible in the render
- [ ] M2 Wetness and edges
- [ ] M3 Pigment identity
- [ ] M4 Easel
- [ ] M5 Paintings and blind review

## Next step

M2, smallest item: wet-on-dry transition width <= 3 px. Then wet-in-wet
(needs pigment mixing in the water) and backruns (needs the wet mask to
grow from capillary water, the [CGW §4.6] σ rule, not built yet).

The engine departs from [CGW] in how surface water moves: a thin film with
depth, and edge darkening from evaporation at the contact line instead of
FlowOutward. research.md §10 and §11 record this.

## Known issues

- [CGW] pseudocode prints `A − μB` in UpdateVelocities, which contradicts
  eq. 1 (research.md §3.2). Not used by the current engine, which has no
  velocity field.
- [CGW] give no values for `c_min, c_max, α, ε, δ, σ` and have no
  evaporation model. All M1 constants except water's own properties are
  estimates (research.md §11); the edge ratio depends most on the
  evaporation boundary length L, paper mobility and pigment settling speed.
- Numba is required in practice; the flow kernel falls back to plain
  Python without it, which is far too slow to use.
- M1 tests take about 100 s (one wash simulated for 18 minutes of paint
  time on a 500 x 250 grid).
