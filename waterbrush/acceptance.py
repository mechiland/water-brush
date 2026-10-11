"""Scenes and measurements used by the milestone acceptance tests.

Kept in the package so the tests and the comparison-image scripts run the
exact same scene.
"""
import numpy as np

from .paper import Paper
from .pigment import BURNT_SIENNA
from .sim import Sim


def flat_wash(seed=7, cells_per_unit=0.5, params=None):
    """M1 scene: one flat wash of burnt sienna on cold-press paper.

    Canvas 1000 x 500 units (300 x 150 mm). The wash covers units
    x 200-800, y 100-400 (180 x 90 mm) with 0.12 mm of water, a juicy but
    not puddled load (estimate). Returns the sim and the wash mask.
    """
    sim = Sim(width_units=1000, height_units=500, cells_per_unit=cells_per_unit,
              paper=Paper("cold"), pigments=[BURNT_SIENNA], seed=seed, params=params)
    c = cells_per_unit
    mask = np.zeros((sim.H, sim.W), dtype=bool)
    mask[int(100 * c):int(400 * c), int(200 * c):int(800 * c)] = True
    sim.apply_wash(mask, water_mm=0.12, amounts=[0.8])
    return sim, mask


def distance_from_edge_mm(mask, dx_mm):
    """For cells inside mask: distance to the nearest outside cell, in mm.

    Chebyshev distance in cells, measured to the cell centre (edge cells are
    0.5 cell in). Cells outside the mask get 0.
    """
    dist = np.zeros(mask.shape)
    inside = mask.copy()
    k = 0
    while inside.any():
        k += 1
        p = np.pad(inside, 1, constant_values=False)
        core = (p[1:-1, 1:-1] & p[:-2, 1:-1] & p[2:, 1:-1] & p[1:-1, :-2] & p[1:-1, 2:]
                & p[:-2, :-2] & p[:-2, 2:] & p[2:, :-2] & p[2:, 2:])
        dist[inside & ~core] = (k - 0.5) * dx_mm
        inside = core
    return dist


def edge_profile(d, mask, dx_mm, max_mm=15.0):
    """Mean deposit as a function of distance from the wash edge."""
    dist = distance_from_edge_mm(mask, dx_mm)
    centres = np.unique(dist[mask & (dist <= max_mm)])
    return centres, np.array([d[mask & (dist == c)].mean() for c in centres])


def edge_ratio(d, mask, dx_mm, edge_mm=1.0, interior_mm=10.0):
    """Mean deposit within edge_mm of the edge over the mean beyond interior_mm."""
    dist = distance_from_edge_mm(mask, dx_mm)
    edge = mask & (dist <= edge_mm)
    interior = mask & (dist >= interior_mm)
    return d[edge].mean() / d[interior].mean()
