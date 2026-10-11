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

## Round 2 (2026-10-11): M1 flat wash

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

Verification:
- `python3 -m pytest`: 6 passed (85 s). Edge band (outer 1 mm) / interior
  (>= 10 mm) = 1.56; wash touch-dry at t = 1080 s.
- Grid check: the same wash at 0.3 mm cells (1 cell per unit, 1000 x 500
  grid) gives 1.74 for the 1 mm band (1.38 for a 2 mm band), dry at
  t = 1075 s, 10 minutes of wall time. The effect is not a grid artefact; the
  test keeps 0.6 mm cells for speed.
- Looked at `renders/m1_flat_wash.png`: a thin darker line follows the
  whole wash border, sharp on the outside and fading inwards over about
  2 mm; interior shows mild cold-press mottling; nothing outside the brushed
  rectangle.

Next:
- M0 still needs the PDF check. M2 smallest item: wet-on-dry transition
  width <= 3 px.
