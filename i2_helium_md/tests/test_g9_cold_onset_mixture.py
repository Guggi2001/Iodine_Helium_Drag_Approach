"""G9 cold-onset mixture helpers (THESIS_FINAL_GOALS G9; scorer
``scripts/post_processing/tier2_g9_cold_onset_table.py``).

The mixture is exact pooling of seed-s1 E0 columns: both fragments of a
molecule come from one column, weights sum to 1, the main-blend solve puts
the expected n̄ at the target, and the assignment partitions the molecules.
"""

from __future__ import annotations

import numpy as np
import pytest

from scripts.post_processing.tier2_g9_cold_onset_table import (
    COLD_SHIFT_NODES, assign_molecules, cold_weights, expected_nbar,
    fragment_mask, mixture_weights, solve_main_blend,
)

NODES = list(range(-4, 3))


@pytest.mark.parametrize("shape", ["P", "U"])
@pytest.mark.parametrize("w", [0.0, 0.2, 0.5])
def test_mixture_weights_sum_to_one(shape, w):
    p = mixture_weights(0, 0.3, w, shape)
    assert sum(p.values()) == pytest.approx(1.0, abs=1e-12)
    assert all(v > 0 for v in p.values())


def test_cold_shapes():
    assert cold_weights("P") == {-COLD_SHIFT_NODES: 1.0}
    u = cold_weights("U")
    assert sorted(u) == [-4, -3, -2, -1, 0]
    assert sum(u.values()) == pytest.approx(1.0)
    assert u[0] == u[-4] == pytest.approx(1 / 8)
    with pytest.raises(ValueError):
        cold_weights("X")


def test_point_shape_moves_w_four_nodes_down():
    p = mixture_weights(1, 0.0, 0.25, "P")
    assert p == pytest.approx({1: 0.75, -3: 0.25})


def test_expected_nbar_is_weighted_ratio():
    n_sum = {0: 4000.0, 1: 3000.0}
    n_scored = {0: 1000.0, 1: 2000.0}
    # (0.5*4000 + 0.5*3000) / (0.5*1000 + 0.5*2000) = 3500/1500
    assert expected_nbar({0: 0.5, 1: 0.5}, n_sum, n_scored) == pytest.approx(3500 / 1500)


def test_solve_main_blend_hits_target_on_stub():
    # stub columns: n̄_k = 4 - 0.25 k, equal scored counts
    n_scored = {k: 1800.0 for k in NODES}
    n_sum = {k: (4.0 - 0.25 * k) * 1800.0 for k in NODES}
    for shape in ("P", "U"):
        for w in (0.0, 0.1, 0.3):
            j, f = solve_main_blend(w, shape, NODES, n_sum, n_scored, target=3.89)
            p = mixture_weights(j, f, w, shape)
            assert expected_nbar(p, n_sum, n_scored) == pytest.approx(3.89, abs=1e-9)
            assert min(p) >= NODES[0] and max(p) <= NODES[-1]


def test_solve_main_blend_off_grid_returns_none():
    n_scored = {k: 1000.0 for k in NODES}
    n_sum = {k: 4.0 * 1000.0 for k in NODES}      # flat: no segment brackets 3.0
    assert solve_main_blend(0.2, "P", NODES, n_sum, n_scored, target=3.0) is None


def test_assignment_partitions_molecules():
    p = {-4: 0.1, 0: 0.55, 1: 0.35}
    groups = assign_molecules(p, 1000)
    allm = np.concatenate(list(groups.values()))
    assert np.array_equal(np.sort(allm), np.arange(1000))
    assert {k: g.size for k, g in groups.items()} == {-4: 100, 0: 550, 1: 350}
    # fixed permutation: the same call is reproducible
    again = assign_molecules(p, 1000)
    assert all(np.array_equal(groups[k], again[k]) for k in groups)


def test_fragment_mask_keeps_pairs_together():
    m = fragment_mask(np.array([0, 3]), 5)
    assert m.tolist() == [True, False, False, True, False,
                          True, False, False, True, False]
