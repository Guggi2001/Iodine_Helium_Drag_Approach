"""Pure helpers of scripts/post_processing/tier2_g13_shed_timing.py.

Synthetic arrays only (no run directories): the exit-column rule of
THESIS_DISCUSSION §12.2 (last outward crossing, first stored column after
it) and the one-event-per-unit-change reading of the n_shell history.
"""

import numpy as np
import pytest

from scripts.post_processing.tier2_g13_shed_timing import (
    exit_columns,
    ke_log_split,
    window_events,
)


def test_ke_log_split_terms_sum_to_log_ke_ratio():
    m0, m1 = np.array([2.0, 3.0]), np.array([1.5, 3.0])
    v0, v1 = np.array([4.0, 2.0]), np.array([4.0, 2.5])
    mass_term, speed_term = ke_log_split(m0, v0, m1, v1)
    np.testing.assert_allclose(mass_term + speed_term,
                               np.log((m1 * v1 ** 2) / (m0 * v0 ** 2)), rtol=1e-14)
    np.testing.assert_allclose(speed_term[0], 0.0)    # co-moving shed: mass term only
    np.testing.assert_allclose(mass_term[1], 0.0)     # pure acceleration: speed term only


def test_ke_log_split_rejects_non_positive_input():
    with pytest.raises(ValueError):
        ke_log_split(np.array([1.0]), np.array([0.0]), np.array([1.0]), np.array([1.0]))


def test_exit_column_is_first_column_after_last_outward_crossing():
    depth = np.array([
        [-3.0, -1.0, 0.5, 2.0, 4.0],     # single crossing between 1 and 2
        [-3.0, 1.0, -0.5, 0.0, 2.0],     # leaves, falls back, leaves again (d = 0 is inside)
        [1.0, 2.0, 3.0, 4.0, 5.0],       # born outside, never inside
        [-3.0, 1.0, 2.0, -1.0, -2.0],    # inside at the end: no exit
    ])
    np.testing.assert_array_equal(exit_columns(depth), [2, 4, 0, -1])


def test_exit_columns_rejects_one_dimensional_input():
    with pytest.raises(ValueError):
        exit_columns(np.zeros(5))


def test_window_events_one_event_per_unit_change():
    n = np.array([
        [21, 21, 19, 19, 18],            # two sheds at col 2, one at col 4
        [21, 22, 22, 20, 20],            # pickup at col 1, two sheds at col 3
    ], dtype=float)
    si, sc, pi, pc = window_events(n)
    np.testing.assert_array_equal(si, [0, 0, 0, 1, 1])
    np.testing.assert_array_equal(sc, [2, 2, 4, 3, 3])
    np.testing.assert_array_equal(pi, [1])
    np.testing.assert_array_equal(pc, [1])
    # net change is exact: n_end = n_start + pickups - sheds per ion
    for k in range(2):
        assert n[k, -1] == n[k, 0] + np.sum(pi == k) - np.sum(si == k)
