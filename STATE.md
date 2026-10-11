# STATE

Current milestone: **M1 Flat wash** (built in a separate thread, draft PR #1, stacked on this M0 branch)

Consecutive failures: 0

## Milestones

- [x] M0 Research: mechanisms and parameters of Curtis 1997 §4-5 in `notes/research.md` (verified against the PDF, round 2)
- [ ] M1 Flat wash
- [ ] M2 Wetness and edges
- [ ] M3 Pigment identity
- [ ] M4 Easel
- [ ] M5 Paintings and blind review

## Next step

M1 is implemented in draft PR #1 (thin-film surface water, edge-enhanced
evaporation; research.md §11 on that branch). When this PR merges, PR #1
retargets onto main and resolves STATE.md / research.md against this
version. Section numbers changed in the verified research.md: FlowOutward
is now §3.4 and drying §3.5 (the round 1 draft had both in §3.3).

## Known issues

- [CGW] pseudocode prints `A − μB` in UpdateVelocities, which contradicts
  eq. 1. We implement `+μB`; see research.md §3.2.
- [CGW] give no values for `c_min, c_max, α, ε, δ, σ` and have no
  evaporation model. Our values are estimates to be calibrated in M1/M2.
- Numba is not preinstalled; `pip install numba` works.
