# STATE

Current milestone: **M0 Research** (draft written, PDF check pending) and **M1 Flat wash** (acceptance passed, stacked on the M0 branch)

Consecutive failures: 0

## Milestones

- [ ] M0 Research: mechanisms and parameters of Curtis 1997 §4-5 in `notes/research.md`
  - [x] Repository skeleton (package, tests, notes, state files)
  - [x] Draft of flow, diffusion, deposition, evaporation, backruns, edge darkening, KM optics
  - [ ] Verify the draft against the PDF and replace every `recalled` tag with a page reference
- [x] M1 Flat wash (passed 2026-10-11; `tests/test_m1_flat_wash.py`, `renders/m1_flat_wash.png`)
  - [x] Dried edge band (outer 1 mm) density >= 1.3 x interior (>= 10 mm in): measured 1.56
  - [x] Hard edge from drying visible in the render
- [ ] M2 Wetness and edges
- [ ] M3 Pigment identity
- [ ] M4 Easel
- [ ] M5 Paintings and blind review

## Next step

Get the Curtis 1997 PDF into the session (either the environment's network
policy allows `grail.cs.washington.edu`, or the PDF is uploaded to the
project files). Then work through "Open items for M0 verification" at the end
of `notes/research.md`. When all are resolved, mark M0 done.

M1 was built before M0 closed (the owner started it in parallel). Its model
departs from [CGW] in one place, recorded in `notes/research.md` §11: the
surface layer is a thin film with depth, and edge darkening comes from
evaporation at the contact line instead of FlowOutward. If the PDF check in
M0 changes anything §11 depends on, re-run M1.

M2 next step (smallest item): wet-on-dry transition width <= 3 px. Then
wet-in-wet (needs pigment mixing in the water) and backruns (needs the
wet mask to grow from capillary water, [CGW §4.6] σ rule, not built yet).

## Known issues

- The paper PDF is unreachable from this environment: the egress proxy
  refuses `grail.cs.washington.edu` and `www.cs.princeton.edu` (403). The
  research draft is written from recall and is tagged as such.
- Numba is not preinstalled; `pip install numba` works. The flow kernel
  falls back to plain Python without it, which is far too slow to use.
- M1 tests take about 90 s (one wash simulated for 18 minutes of paint time on a 500 x 250 grid).
- All M1 physical constants except water's own properties are estimates
  (research.md §11); the edge ratio depends most on the evaporation
  boundary length L, paper mobility and the pigment settling speed.
