# STATE

Current milestone: **M0 Research** (in progress, draft written)

Consecutive failures: 0

## Milestones

- [ ] M0 Research: mechanisms and parameters of Curtis 1997 §4-5 in `notes/research.md`
  - [x] Repository skeleton (package, tests, notes, state files)
  - [x] Draft of flow, diffusion, deposition, evaporation, backruns, edge darkening, KM optics
  - [ ] Verify the draft against the PDF and replace every `recalled` tag with a page reference
- [ ] M1 Flat wash
- [ ] M2 Wetness and edges
- [ ] M3 Pigment identity
- [ ] M4 Easel
- [ ] M5 Paintings and blind review

## Next step

Get the Curtis 1997 PDF into the session (either the environment's network
policy allows `grail.cs.washington.edu`, or the PDF is uploaded to the
project files). Then work through "Open items for M0 verification" at the end
of `notes/research.md`. When all are resolved, mark M0 done and start M1 with
the smallest item: paper height field plus a static wet-area mask and the
shallow-water velocity update.

## Known issues

- The paper PDF is unreachable from this environment: the egress proxy
  refuses `grail.cs.washington.edu` and `www.cs.princeton.edu` (403). The
  research draft is written from recall and is tagged as such.
- Numba is not preinstalled; `pip install numba` works.
