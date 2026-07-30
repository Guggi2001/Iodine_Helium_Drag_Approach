"""Tests for the free-form linear sweep twin extension (arms 1/2).

Covers `scripts/tier2_h2b_forward_model.py`
(`TIER2_FREEFORM_LINEAR_TWIN_SWEEP_PLAN.md` §1/§7/§9):

- ``integrate_pairs`` form-switch validation — the counterfactual forms
  must fail loudly on wrong parameter combinations (rule 4);
- the §7.4 seam: ``linq`` at c = 0 is bit-identical to ``lin`` (toy
  ensemble, shortened window — the identity is per-step arithmetic, so it
  holds at any t_end);
- drag-direction sanity: the lin form dissipates (v_inf below the
  ballistic bracket) and more a means more dissipation;
- the Φ / sub-bare / R2-class helpers (plan §1/§2/§5 arithmetic).

Toy ensembles only — no npz caches, no corrected master, no figures
(testing rules). Tolerances: bit-exact where the seam guarantees identity;
strict inequalities for the direction checks.
"""

from __future__ import annotations

import numpy as np
import pytest

from scripts.tier2_h2b_forward_model import (
    PHI_DENOM,
    integrate_pairs,
    linsweep_f97,
    linsweep_phi_class,
    linsweep_subbare,
)

# A 3-molecule toy: mid-droplet births in a 40 A droplet, mixed cosines.
TOY = dict(
    r0=np.array([5.0, 15.0, 25.0]),
    mu=np.array([0.3, -0.5, 0.9]),
    R_drop=np.array([40.0, 40.0, 40.0]),
    mass_amu=np.array([203.0, 203.0, 203.0]),
)
TOY_T_END = 5.0  # ps — short window, enough steps to accumulate arithmetic


class TestFormSwitchValidation:
    def test_unknown_form_rejected(self):
        with pytest.raises(ValueError, match="unknown drag_form"):
            integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="quadratic")

    def test_cubic_rejects_lin_params(self):
        with pytest.raises(ValueError, match="lin_a/linq_c"):
            integrate_pairs(**TOY, t_end=TOY_T_END, lin_a=40.0)

    def test_lin_requires_positive_a(self):
        with pytest.raises(ValueError, match="requires lin_a > 0"):
            integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="lin")
        with pytest.raises(ValueError, match="requires lin_a > 0"):
            integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="lin",
                            lin_a=-1.0)

    def test_lin_forbids_cap_and_c(self):
        with pytest.raises(ValueError, match="v_c/p_tail forbidden"):
            integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="lin",
                            lin_a=40.0, v_c=5.5, p_tail=-1.0)
        with pytest.raises(ValueError, match="linq_c requires"):
            integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="lin",
                            lin_a=40.0, linq_c=2.0)

    def test_linq_requires_nonnegative_c(self):
        with pytest.raises(ValueError, match="requires linq_c >= 0"):
            integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="linq",
                            lin_a=40.0)
        with pytest.raises(ValueError, match="requires linq_c >= 0"):
            integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="linq",
                            lin_a=40.0, linq_c=-0.5)


class TestSeamAndDirection:
    def test_linq_c0_bit_identical_to_lin(self):
        res_lin = integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="lin",
                                  lin_a=40.0)
        res_c0 = integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="linq",
                                 lin_a=40.0, linq_c=0.0)
        for key in ("K", "v_inf", "v_peak", "depth_end", "trapped"):
            assert np.array_equal(res_lin[key], res_c0[key]), key
        assert np.array_equal(res_lin["t_exit"], res_c0["t_exit"],
                              equal_nan=True)

    def test_lin_dissipates_below_ballistic_and_monotone_in_a(self):
        ball = integrate_pairs(**TOY, t_end=TOY_T_END, drag_on=False)
        soft = integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="lin",
                               lin_a=20.0)
        hard = integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="lin",
                               lin_a=40.0)
        assert (soft["v_inf"] < ball["v_inf"]).all()
        assert (hard["v_inf"] < soft["v_inf"]).all()

    def test_linq_c_adds_drag(self):
        lin = integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="lin",
                              lin_a=20.0)
        linq = integrate_pairs(**TOY, t_end=TOY_T_END, drag_form="linq",
                               lin_a=20.0, linq_c=6.13)
        assert (linq["v_inf"] < lin["v_inf"]).all()


class TestPhiHelpers:
    def test_f97_and_phi_anchors(self):
        # Plan §1: a ~ 43.1 -> Phi = 1; a ~ 30.2 -> Phi = 0.7.
        assert linsweep_f97(43.14) / PHI_DENOM == pytest.approx(1.0, abs=5e-3)
        assert linsweep_f97(30.2) / PHI_DENOM == pytest.approx(0.7, abs=5e-3)
        # linq: F(9.7) = (a + c*9.7)*9.7.
        assert linsweep_f97(20.0, 2.0) == pytest.approx((20 + 19.4) * 9.7)

    def test_subbare_band(self):
        assert linsweep_subbare(170.0) == 2  # below the whole band
        assert linsweep_subbare(200.0) == 1  # inside 175-245
        assert linsweep_subbare(300.0) == 0  # above

    def test_r2_classes(self):
        assert linsweep_phi_class(0.70) == "sub_plateau"
        assert linsweep_phi_class(0.75) == "intermediate"
        assert linsweep_phi_class(1.0) == "plateau_convergent"
        assert linsweep_phi_class(1.2) == "above_window"
