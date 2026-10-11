# LOG

## Round 1 (2026-10-11)

Done:
- Created repository skeleton: `pyproject.toml`, `waterbrush/` package,
  `tests/` with an import smoke test, `renders/`, `notes/research.md`,
  `notes/friction.md`, `STATE.md`, `LOG.md`, `README.md`.
- Tried to download Curtis et al. 1997 from grail.cs.washington.edu and
  the Princeton course mirror. Both refused by the environment's egress
  proxy (403).
- Wrote `notes/research.md` from recall of the paper: three-layer model,
  velocity update, divergence relaxation, flow outward (edge darkening),
  pigment advection, transfer (deposition, lifting, granulation, staining),
  capillary flow and backruns, drybrush, KM optics and compositing, and
  project-specific mappings (paper parameters, time units, determinism).
  Every [CGW] item is tagged `recalled` or `recalled, uncertain`; values
  without a source are marked 估计.

Verification:
- `python3 -m pytest`: 1 passed.
- M0 acceptance (mechanisms and parameters with sources) is not met yet,
  because the sources are recalled, not read. M0 stays open.

Next:
- Obtain the PDF, verify each tagged item, record page numbers, fix
  anything recalled wrongly. Then mark M0 done and begin M1.
