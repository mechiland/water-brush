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

## Round 2 (2026-10-11)

Done:
- The network policy now allows grail.cs.washington.edu. Downloaded the
  PDF (10 pages) and read §3-5 plus the related parts of §2, §6.1 and
  Future Work, checking the pseudocode and the Fig. 5 table on rendered
  page images, since the text extraction drops Greek letters and minus signs.
- Rewrote `notes/research.md` with section and page citations.
  Corrections to the round 1 recall:
  - η is a range, 0.01 ≤ η ≤ 0.05, not a single 0.05.
  - The printed pseudocode really has `A − μB`; eq. 1 has `+μ∇²u`. We
    choose `+μB` and say why.
  - [CGW] have no evaporation routine, and capillary absorption does not
    remove surface water. Our explicit drying is a departure, listed.
  - [CGW] give no values for α, ε, δ, σ, c_min, c_max. Marked 估计.
  - The backrun inequality δ < σ is our inference from the pseudocode,
    not stated in the paper.
  - MovePigment as printed has no Δt factor; we add it.
  - The capillary layer is used only for backruns in [CGW].
- Added the full Fig. 5 pigment table (K, S, ρ, ω, γ for 12 pigments) and
  a starting-point plan for our 8 tube pigments.

Verification:
- `python3 -m pytest`: 1 passed.
- M0 acceptance: every mechanism (flow, diffusion, deposition,
  evaporation, backruns, edge darkening) has a mechanism description and
  parameters with a source or an 估计 tag. Passed. M0 marked done.

Next:
- M1: paper height field and capillary capacity.
