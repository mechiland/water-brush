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

## Round 3 (2026-10-11): M1 flat wash

Done:
- Built the three-layer simulation (`waterbrush/sim.py`): surface water
  with depth, wet mask and suspended pigment; deposited pigment; capillary
  water in the paper. Thin-film flow (Numba kernel in `_kernels.py`),
  evaporation that is stronger near the wet edge (Deegan 1997), Lucas-
  Washburn absorption calibrated to Cobb60, Curtis-style settling and
  lifting with a Stokes settling speed, capillary wicking.
- Paper height field for hot/cold/rough with sizing and weight
  (`paper.py`); KM pigment from reflectance over white and black, one M1
  pigment, burnt sienna (`pigment.py`); KM renderer (`render.py`).
- Acceptance scene and measurements in `waterbrush/acceptance.py`; tests
  in `tests/test_m1_flat_wash.py`; comparison image from
  `scripts/m1_flat_wash.py` to `renders/m1_flat_wash.png`.
- research.md §11 records the model, the departure from [CGW] (thin film
  and evaporation-driven edge flow instead of FlowOutward) and every
  constant's source or 估计 tag.

Tries before passing:
1. Gravity-only thin film, fast absorption (Cobb60 ~ 60 g/m²), settling
   as a fixed rate: ratio about 1.05. Water soaked in within 2-3 minutes,
   so no time for edge flow; paper relief 0.08 mm made blotchy puddles.
2. Washburn absorption with hard sizing (Cobb60 20 g/m²), Stokes
   settling, relief 0.05 mm: outermost cell 1.36, 1 mm band about 1.25.
   The edge dried first because the thin rim could not be refilled.
3. Evaporation profile sqrt(L/x) instead of a blurred-mask boost: 1.31
   in the outermost cells, 1 mm band about 1.2.
4. Added capillary pressure to the film mobility (factor 1 + (l_c/ℓ)²):
   1 mm band 1.56. Passed.
5. Merged the verified M0 notes. Capillary capacity now follows [CGW]'s
   `c = c_min + h (c_max − c_min)` (higher on peaks; it was the other way
   round). §11 moved on top of the verified research.md with section
   references updated (FlowOutward §3.4, drying §3.5), and three rows
   added to the §10 departure table. Re-ran: 1.57.

Verification:
- `python3 -m pytest`: 6 passed (about 100 s). Edge band (outer 1 mm) /
  interior (>= 10 mm) = 1.57; wash touch-dry at t = 1110 s.
- Grid check: the same wash at 0.3 mm cells (1 cell per unit, 1000 x 500
  grid) gives 1.74 for the 1 mm band (1.38 for a 2 mm band), dry at
  t = 1075 s, 10 minutes of wall time. The effect is not a grid artefact; the
  test keeps 0.6 mm cells for speed.
- Looked at `renders/m1_flat_wash.png`: a thin darker line follows the
  whole wash border, sharp on the outside and fading inwards over about
  2 mm; interior shows mild cold-press mottling; nothing outside the brushed
  rectangle.

Next:
- M2 smallest item: wet-on-dry transition
  width <= 3 px.
