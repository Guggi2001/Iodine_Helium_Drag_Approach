"""G4.14 (THESIS_FINAL_GOALS, frozen 7726e43): seed block, Step-0 nbar
inversion / filter, surface fitter + residual criterion, feasible-set selection.

All synthetic; no MD, no twin integration.
"""

import math

import numpy as np
import pytest

from scripts.gen_tier2_detfix_battery import MEMBER_SEEDS
from scripts.gen_tier2_g414_tauline import (
    NBAR_TARGETS, OTHER_TAU_PS, PRIMARY_TAU_PS, SEED_BASE, cell_label, plan_cells,
    verify_seed_block,
)
from scripts.post_processing.tier2_g414_table import (
    N1_MIN, c_status, de_at_nbar, fit_observable, fit_surface, predict,
    select_candidates, tradeoff_curve,
)
from scripts.tier2_h2b_forward_model import (
    g414_e0_at_nbar, g414_interp_at, g414_md_prediction, g414_midhot_filter,
)


# --- seed block / plan --------------------------------------------------------
def _line(cs, taus):
    return [{"c": str(c), "tau_ps": str(t), "nbar_target_md": str(n),
             "E0_eV": str(0.3 + 0.01 * t + 0.001 * n)}
            for c in cs for t in taus for n in NBAR_TARGETS]


def test_primary_tau_set_is_frozen():
    assert PRIMARY_TAU_PS == (4.0, 4.4, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6)
    assert OTHER_TAU_PS == (4.8, 5.2, 5.6, 6.0, 6.4)


def test_plan_seed_block_and_order():
    taus = sorted(set(PRIMARY_TAU_PS) | set(OTHER_TAU_PS))
    filt = [{"c": "2.0", "passes_midhot_filter": "1"},
            {"c": "2.5", "passes_midhot_filter": "1"},
            {"c": "1.8", "passes_midhot_filter": "1"},
            {"c": "1.6", "passes_midhot_filter": "0"}]
    cells = plan_cells(_line([1.6, 1.8, 2.0, 2.5], taus), filt)
    assert len(cells) == 24 + 10 + 10
    assert [c.seed for c in cells] == list(range(SEED_BASE, SEED_BASE + 44))
    assert not {c.seed for c in cells} & set(MEMBER_SEEDS.values())
    assert [c.c for c in cells[:24]] == [2.0] * 24          # primary first
    assert [c.c for c in cells[24:]] == [1.8] * 10 + [2.5] * 10
    assert cells[0].label == cell_label(2.0, 4.0, 3.85) == "c200t40n385"
    assert cells[1].nbar_target == 4.05
    assert len({c.label for c in cells}) == len(cells)


def test_plan_rejects_missing_e0():
    rows = _line([2.0], PRIMARY_TAU_PS)
    rows[3]["E0_eV"] = "nan"
    with pytest.raises(ValueError, match="no Step-0 E0"):
        plan_cells(rows, [{"c": "2.0", "passes_midhot_filter": "1"}])


def test_seed_block_guard():
    verify_seed_block([SEED_BASE, SEED_BASE + 1])
    with pytest.raises(AssertionError):
        verify_seed_block([SEED_BASE, SEED_BASE + 2])
    with pytest.raises(AssertionError):
        verify_seed_block([SEED_BASE + 1])


# --- Step 0 helpers -------------------------------------------------------------
def test_e0_at_nbar_linear_stub():
    e0 = np.round(np.arange(0.30, 0.40001, 0.0025), 4)
    nbar = 6.0 - 25.0 * (e0 - 0.30)          # -0.25 per 0.01 eV
    e_star, i = g414_e0_at_nbar(e0, nbar, 5.0)
    assert e_star == pytest.approx(0.34, abs=1e-12)
    assert g414_interp_at(e0, nbar, e_star, i) == pytest.approx(5.0, abs=1e-12)
    assert math.isnan(g414_e0_at_nbar(e0, nbar, 9.0)[0])
    with pytest.raises(ValueError):
        g414_e0_at_nbar(e0[::-1], nbar[::-1], 5.0)


def test_md_prediction_and_filter():
    assert g414_md_prediction(1.0, 0.945, "ratio") == pytest.approx(0.945)
    assert g414_md_prediction(0.6, 0.02, "diff") == pytest.approx(0.62)
    rows = [{"tau_ps": 4.0, "md_pred_midHot": 1.0},       # in band, out of tau range
            {"tau_ps": 5.0, "md_pred_midHot": 1.20}]
    assert not g414_midhot_filter(rows)
    rows.append({"tau_ps": 6.6, "md_pred_midHot": 1.15})  # closed edges
    assert g414_midhot_filter(rows)


# --- fitter / residual criterion ------------------------------------------------
def _grid():
    tau = np.repeat(np.round(np.arange(4.0, 6.61, 0.2), 1), 2)
    de = np.tile([0.0, -0.8], tau.size // 2)
    return tau, de


def test_fitter_recovers_synthetic_surface():
    tau, de = _grid()
    t = tau - 5.6
    truth = 0.20 + 0.01 * t - 0.004 * t ** 2 + 0.001 * t ** 3 + 0.006 * de + 0.002 * t * de
    rng = np.random.default_rng(1)
    y = truth + rng.normal(0.0, 0.002, tau.size)
    res = fit_observable(tau, de, y, 3, False, seed_sd=0.002)
    assert res["status"] == "pass"
    assert res["fit"].beta[4] == pytest.approx(0.006, abs=0.002)
    assert np.max(np.abs(predict(res["fit"], tau, de) - truth)) < 0.003


def test_residual_criterion_adds_tau_order_on_structure():
    tau, de = _grid()
    t = tau - 5.6
    y = 0.9 + 0.05 * t ** 4 + 0.01 * de        # quartic: a cubic misses it
    res = fit_observable(tau, de, y, 3, False, seed_sd=0.001)
    assert res["p_tau"] < 0.05
    assert res["order"] == 4 and res["status"] == "pass_after_order"


def test_residual_criterion_white_noise_asks_second_seed():
    tau, de = _grid()
    rng = np.random.default_rng(7)
    y = 4.0 + 0.3 * de + rng.normal(0.0, 0.2, tau.size)   # noise >> 1.5 x sd
    res = fit_observable(tau, de, y, 3, False, seed_sd=0.07)
    assert res["status"] == "second_seed"
    assert c_status(["pass", "second_seed"]) == "second_seed"
    assert c_status(["pass", "unresolved", "second_seed"]) == "unresolved"
    assert c_status(["pass", "pass_after_order"]) == "resolved"


def test_fit_needs_dof():
    with pytest.raises(ValueError):
        fit_surface([5.0, 5.2], [0.0, 1.0], [1.0, 2.0], 1, False)


# --- selection -------------------------------------------------------------------
def _synthetic_fits():
    tau, de = _grid()
    t = tau - 5.6
    obs = {
        "nbar_det": 3.85 - 0.25 * de,                          # dE = -0.16 at 3.89
        "n1_solv": 0.195 + 0.02 * t + 0.0 * de,                # >= 0.200 for tau >= 5.85
        "KE1_mean": 0.60 + 0.03 * t,                           # rises with tau
        "w1_solv": 0.90 + 0.10 * (t - 0.0) ** 2,               # min at 5.6
        "midHot": 1.05 + 0.0 * t,
    }
    return {k: fit_surface(tau, de, v, 3, k == "n1_solv") for k, v in obs.items()}


def test_de_at_nbar_inverts_linear_nbar():
    fits = _synthetic_fits()
    assert de_at_nbar(fits["nbar_det"], [5.0, 6.0]) == pytest.approx([-0.16, -0.16])


def test_tradeoff_and_candidates():
    fits = _synthetic_fits()
    curve = tradeoff_curve(fits, (4.0, 6.6), (-0.8, 0.0))
    feas = [r for r in curve if r["feasible"]]
    assert min(r["tau_ps"] for r in feas) == pytest.approx(5.85, abs=0.011)
    assert all(r["pred_n1_solv"] >= N1_MIN for r in feas)
    cands = select_candidates(curve)
    roles = {r["role"]: r["tau_ps"] for r in cands}
    assert roles["KE"] == pytest.approx(6.6)                 # KE1 max at the top
    assert roles["W"] == pytest.approx(5.85, abs=0.011)     # W1 min inside feasible
    # out-of-sampled dE -> nothing feasible
    assert select_candidates(tradeoff_curve(fits, (4.0, 6.6), (0.0, 0.5))) == []


def test_candidates_merge_within_0p2():
    curve = [{"tau_ps": 5.0, "feasible": 1, "pred_KE1_mean": 0.6, "pred_w1_solv": 0.9},
             {"tau_ps": 5.15, "feasible": 1, "pred_KE1_mean": 0.5, "pred_w1_solv": 0.8}]
    out = select_candidates(curve)
    assert len(out) == 1 and out[0]["role"] == "KE+W" and out[0]["tau_ps"] == 5.0
