# water-brush

A physics-based watercolor simulator with an "easel" API that an AI painter
drives by running code blocks. Every visible effect (edge darkening, backruns,
granulation, wet-in-wet soft edges, glazing) must come out of the simulation;
no filters, overlays, textures or image models.

- Physics: Curtis et al. 1997, *Computer-Generated Watercolor* (three-layer model).
- Optics: Kubelka-Munk.
- Python 3 + NumPy, optional Numba. Deterministic: a seed plus `painting.log` replays pixel for pixel.

## Layout

- `waterbrush/`: the engine package.
- `tests/`: acceptance and unit tests (`python3 -m pytest`).
- `renders/`: comparison images produced by acceptance runs.
- `notes/research.md`: physical mechanisms and the source of every parameter.
- `notes/friction.md`: problems the painter hits while using the easel.
- `STATE.md`: current milestone, next step, known issues, consecutive failures.
- `LOG.md`: one entry per loop round.
