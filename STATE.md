# STATE

Current milestone: **M1 Flat wash** (not started)

Consecutive failures: 0

## Milestones

- [x] M0 Research: mechanisms and parameters of Curtis 1997 §4-5 in `notes/research.md` (verified against the PDF, round 2)
- [ ] M1 Flat wash
- [ ] M2 Wetness and edges
- [ ] M3 Pigment identity
- [ ] M4 Easel
- [ ] M5 Paintings and blind review

## Next step

M1, smallest item: paper height field (Worley + Perlin, seeded) with
`c = h (c_max − c_min) + c_min`, plus a render of the three textures.
Then the shallow-water velocity update and divergence relaxation on a
static mask, then FlowOutward, pigment transport and transfer, explicit
drying, and the flat-wash edge/interior density test.

## Known issues

- [CGW] pseudocode prints `A − μB` in UpdateVelocities, which contradicts
  eq. 1. We implement `+μB`; see research.md §3.2.
- [CGW] give no values for `c_min, c_max, α, ε, δ, σ` and have no
  evaporation model. Our values are estimates to be calibrated in M1/M2.
- Numba is not preinstalled; `pip install numba` works.
