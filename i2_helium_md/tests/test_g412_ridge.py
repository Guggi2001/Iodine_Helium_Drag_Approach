"""THESIS_FINAL_GOALS G4.12: ridge-point force arithmetic and exit-bound scoring.

Synthetic inputs only (no MD, no reference files). The fixed-c fit itself is
covered on the stub landscape in test_extraction_trajectory_matching.py.
"""

import numpy as np
import pytest

from scripts.extraction.method_b_g412_ridge import (
    AT_BOUND_APS,
    exit_bound,
    lq_force,
    point_label,
    rmse_on_grid,
    transverse_fraction,
    window_mean,
)


def test_lq_force_matches_draft_projection():
    # F(3) = 67 fixed -> a = 22.33 - 3c; F(9.7) = 216.6 + 65c (G4.12 "Why").
    for c in (0.0, 1.0, 3.1):
        a = (67.0 - 9.0 * c) / 3.0
        assert lq_force(a, c, 3.0) == pytest.approx(67.0, rel=1e-12)
        # exact: 9.7 * 67/3 + (9.7**2 - 9.7 * 3) c = 216.633 + 64.99 c
        assert lq_force(a, c, 9.7) == pytest.approx(9.7 * 67.0 / 3.0 + 64.99 * c, rel=1e-12)


def test_point_label():
    assert point_label(0.0) == "c0p000"
    assert point_label(7.224) == "c7p224"


def test_window_mean_closed_and_empty():
    t = np.arange(0.0, 15.0, 0.5)
    assert window_mean(t, t, (13.0, 14.0)) == pytest.approx(13.5)
    with pytest.raises(ValueError, match="no samples"):
        window_mean(t, t, (20.0, 21.0))


def test_transverse_fraction_radial_and_tangential():
    pos = np.array([[0.0, 0.0, 10.0], [0.0, 10.0, 0.0], [3.0, 0.0, 4.0]])
    vel = np.array([[0.0, 0.0, 2.0], [1.0, 0.0, 0.0], [3.0, 0.0, 4.0]])
    np.testing.assert_allclose(transverse_fraction(pos, vel), [0.0, 1.0, 0.0], atol=1e-15)
    # 45 degrees to radial -> sin 45.
    f = transverse_fraction(np.array([[0.0, 0.0, 5.0]]), np.array([[1.0, 0.0, 1.0]]))
    assert f[0] == pytest.approx(np.sqrt(0.5), rel=1e-14)


def test_rmse_on_grid_interpolates_reference():
    t_md = np.linspace(0.0, 10.0, 101)
    t_ref = np.linspace(0.0, 10.0, 11)
    assert rmse_on_grid(t_md, 2.0 * t_md, t_ref, 2.0 * t_ref, (2.0, 8.0)) == pytest.approx(0.0, abs=1e-12)
    assert rmse_on_grid(t_md, t_md + 0.3, t_ref, t_ref, (2.0, 8.0)) == pytest.approx(0.3, rel=1e-12)


def test_exit_bound_sign_and_flag():
    assert exit_bound(4.3, 4.2)["pass"] and not exit_bound(4.3, 4.2)["at_bound"]
    b = exit_bound(4.18, 4.2)
    assert not b["pass"] and b["at_bound"] and b["margin_Aps"] == pytest.approx(-0.02)
    assert exit_bound(4.2, 4.2)["pass"]   # equal speed passes (">=")
    assert AT_BOUND_APS == 0.05
