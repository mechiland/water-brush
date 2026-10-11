"""M1 acceptance: a flat wash on cold-press paper, left to dry.

1. After drying, mean pigment density within 1 mm of the wash edge is at
   least 1.3 times the mean density more than 10 mm inside.
2. The dried edge is hard: density peaks at the outermost cells and falls
   off inwards, and there is no pigment outside the brushed area.

The comparison image is made by scripts/m1_flat_wash.py.
"""
import numpy as np
import pytest

from waterbrush.acceptance import edge_profile, edge_ratio, flat_wash


@pytest.fixture(scope="module")
def dried():
    sim, mask = flat_wash()
    total = sim.G.sum()
    sim.run_until_dry()
    return sim, mask, total


def test_dries(dried):
    sim, mask, _ = dried
    assert sim.is_dry()
    assert not sim.M.any()
    assert sim.G.sum() == 0.0
    # Touch-dry in a time a painter would recognise (5-25 minutes, estimate).
    assert 300 <= sim.t <= 1500


def test_pigment_conserved(dried):
    sim, _, total = dried
    assert sim.d.sum() == pytest.approx(total, rel=1e-9)


def test_edge_darkening_ratio(dried):
    sim, mask, _ = dried
    ratio = edge_ratio(sim.d[0], mask, sim.dx, edge_mm=1.0, interior_mm=10.0)
    assert ratio >= 1.3, f"edge/interior = {ratio:.3f}"


def test_hard_edge(dried):
    sim, mask, _ = dried
    d = sim.d[0]
    assert d[~mask].max() == 0.0
    dist, prof = edge_profile(d, mask, sim.dx, max_mm=5.0)
    # The darkest band is the outermost one, and density falls off inwards.
    assert prof[0] == prof.max()
    assert np.all(np.diff(prof[:4]) < 0)
    assert prof[0] >= 1.3 * prof[-1]


def test_replay_is_identical():
    a, _ = flat_wash(seed=3)
    b, _ = flat_wash(seed=3)
    a.run(120.0)
    b.run(120.0)
    assert np.array_equal(a.w, b.w)
    assert np.array_equal(a.d, b.d)
    assert np.array_equal(a.s, b.s)
