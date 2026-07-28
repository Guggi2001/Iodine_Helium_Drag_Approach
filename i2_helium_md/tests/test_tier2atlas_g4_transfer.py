"""Atlas G4 Step 1 Block 0 (plan §3.5e) — pure-function unit tests.

Synthetic input only: no run directories, no reference CSVs, no figures. The
report's oracles and joins are exercised by running the script itself (plan
Task 1 Step 5), not here.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.post_processing.tier2atlas_g4_transfer import (  # noqa: E402
    joint_score,
    nbar_bias_model,
    predict_md_nbar,
    rank_license,
)

PROVISIONAL = {
    "w1_ref": 0.571,
    "midhot_ln": float(np.log(1.15)),
    "deepke_ln": float(np.log(1.15)),
}


def test_joint_score_reproduces_the_preregistered_standing_value():
    """finc1v725 pooled (W1 0.571, midHot 1.011, deepKE 0.631) -> S = 4.37.

    The plan §3.5e table is the pre-registration; the score must reproduce it
    under the provisional normalizers or the design's target (S < 4.37) means
    something different from what was registered.
    """
    assert joint_score(0.571, 1.011, 0.631, PROVISIONAL) == pytest.approx(
        4.37, abs=0.01
    )


def test_joint_score_is_zero_at_a_perfect_cell():
    assert joint_score(0.0, 1.0, 1.0, PROVISIONAL) == pytest.approx(0.0)


def test_joint_score_is_symmetric_in_the_ke_ratio():
    """A ratio and its reciprocal are equally wrong — the log form ensures it."""
    assert joint_score(0.0, 1.25, 1.0, PROVISIONAL) == pytest.approx(
        joint_score(0.0, 1 / 1.25, 1.0, PROVISIONAL)
    )


def test_joint_score_drops_an_unlicensed_axis():
    """Block 0 may license only some axes; S must then omit the rest."""
    full = joint_score(0.571, 1.011, 0.631, PROVISIONAL)
    without_deep = joint_score(
        0.571, 1.011, 0.631, PROVISIONAL, axes=("w1", "midhot")
    )
    assert without_deep < full
    assert without_deep == pytest.approx(1.0 + abs(np.log(1.011)) / np.log(1.15))


def test_joint_score_is_nan_when_a_scored_axis_is_nan():
    assert np.isnan(joint_score(0.571, np.nan, 0.631, PROVISIONAL))


def test_nbar_bias_model_recovers_a_planted_line():
    rows = [
        {"twin_nbar": x, "d_nbar": -(0.1 + 0.2 * x)}
        for x in (4.0, 5.0, 6.0, 7.0, 17.0)
    ]
    model = nbar_bias_model(rows)
    assert model["a"] == pytest.approx(-0.1, abs=1e-9)
    assert model["b"] == pytest.approx(-0.2, abs=1e-9)
    assert model["resid_sd"] == pytest.approx(0.0, abs=1e-9)
    assert model["r2"] == pytest.approx(1.0, abs=1e-9)
    assert model["n"] == 5


def test_predict_md_nbar_inverts_the_fit():
    model = {"a": -0.1, "b": -0.2, "resid_sd": 0.0, "r2": 1.0, "n": 5}
    # MD = twin + (a + b*twin) = 5.0 + (-0.1 - 1.0)
    assert predict_md_nbar(5.0, model) == pytest.approx(3.9, abs=1e-9)


def test_nbar_bias_model_needs_at_least_three_points():
    with pytest.raises(ValueError, match="at least 3"):
        nbar_bias_model([{"twin_nbar": 4.0, "d_nbar": -0.5},
                         {"twin_nbar": 5.0, "d_nbar": -0.7}])


def test_rank_license_is_a_hard_threshold_at_0p7():
    twin = [1.0, 2.0, 3.0, 4.0, 5.0]
    up = rank_license(twin, [1.0, 2.0, 3.0, 4.0, 5.0])
    down = rank_license(twin, [5.0, 4.0, 3.0, 2.0, 1.0])
    assert up["licensed"] is True
    assert up["rho"] == pytest.approx(1.0)
    assert up["n"] == 5
    assert down["licensed"] is False


def test_rank_license_ignores_nan_pairs():
    twin = [1.0, 2.0, np.nan, 4.0, 5.0]
    md = [1.0, 2.0, 3.0, np.nan, 5.0]
    assert rank_license(twin, md)["n"] == 3


def test_rank_license_reports_a_ci_that_brackets_rho():
    rng = np.random.default_rng(7)
    twin = np.arange(14, dtype=float)
    md = twin + rng.normal(0.0, 2.0, size=14)
    lic = rank_license(list(twin), list(md), n_boot=200, seed=1)
    assert lic["ci_lo"] <= lic["rho"] <= lic["ci_hi"]
    assert -1.0 <= lic["ci_lo"] <= 1.0
    assert -1.0 <= lic["ci_hi"] <= 1.0
