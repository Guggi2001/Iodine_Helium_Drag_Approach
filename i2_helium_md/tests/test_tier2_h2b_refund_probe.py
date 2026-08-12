"""Focused tests for the RQ12 refund-geometry probe (twin `refundscan`).

The probe asks whether the cashed fraction of the E_bind exit toll,
``c = -dKE/dE_bind``, depends on droplet radius. If it does, the Tier-0
co-extraction of E_bind at R = 9/18 A mis-transfers to production at
R ~ 48 A by ``T = c(R_prod)/c(R_cal)``.

These tests pin the frozen grid + prediction bands and, critically, that
the new ``rho_steepness`` override is **byte-inert by default** and
decouples the *density* width from the *solvation-potential* width (the
twin-side mirror of production's `erf_independent` gate).
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

SCRIPT = (
    Path(__file__).resolve().parents[1] / "scripts" / "tier2_h2b_forward_model.py"
)


@pytest.fixture(scope="module")
def twin():
    spec = importlib.util.spec_from_file_location(
        "tier2_h2b_refund_probe_under_test", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def args(twin):
    return (np.array([0.0]), np.array([1.0]), np.array([47.8]),
            np.array([float(twin.complex_mass_amu(twin.N_STAR))]))


class TestRhoSteepnessOverride:
    def test_default_is_byte_inert(self, twin, args):
        """Omitting the override, passing None, and passing the shared width
        must all reproduce the existing path bit-for-bit."""
        kw = dict(r0_sep=twin.R0_SEP_PROD_A, drag_on=True, v_c=5.5,
                  p_tail=-1.0)
        base = twin.integrate_pairs(*args, **kw)
        none_ = twin.integrate_pairs(*args, rho_steepness=None, **kw)
        tied = twin.integrate_pairs(*args, rho_steepness=twin.STEEP_A, **kw)
        for key in ("K", "v_inf", "depth_end"):
            assert np.array_equal(base[key], none_[key])
            assert np.array_equal(base[key], tied[key])

    def test_rejects_nonpositive(self, twin, args):
        with pytest.raises(ValueError, match="rho_steepness must be > 0"):
            twin.integrate_pairs(*args, r0_sep=twin.R0_SEP_PROD_A,
                                 drag_on=True, rho_steepness=0.0)

    def test_touches_density_only_not_the_well(self, twin, args):
        """With drag OFF the width must be inert (the well force stays on
        STEEP_A); with drag ON it must bite."""
        off_wide = twin.integrate_pairs(
            *args, r0_sep=twin.R0_SEP_PROD_A, drag_on=False,
            rho_steepness=14.2)
        off_sharp = twin.integrate_pairs(
            *args, r0_sep=twin.R0_SEP_PROD_A, drag_on=False,
            rho_steepness=3.14)
        assert np.array_equal(off_wide["v_inf"], off_sharp["v_inf"])

        kw = dict(r0_sep=twin.R0_SEP_PROD_A, drag_on=True, v_c=5.5,
                  p_tail=-1.0)
        on_wide = twin.integrate_pairs(*args, rho_steepness=14.2, **kw)
        on_sharp = twin.integrate_pairs(*args, rho_steepness=3.14, **kw)
        assert not np.array_equal(on_wide["v_inf"], on_sharp["v_inf"])


class TestFrozenGrid:
    def test_covers_calibration_and_production_radii(self, twin):
        """The transfer factor is meaningless without both ends."""
        for R in (9.0, 18.0, 26.6, 47.8):
            assert R in twin.REFUND_RADII_A

    def test_steepness_grid_spans_standing_to_literature(self, twin):
        assert twin.REFUND_S_RHO_A[0] == 14.2      # standing shared width
        assert twin.REFUND_S_RHO_A[-1] == 3.14     # Harms DFT 5.7 A / 1.8124
        assert list(twin.REFUND_S_RHO_A) == sorted(
            twin.REFUND_S_RHO_A, reverse=True)

    def test_wells_match_the_measured_step(self, twin):
        """Same well pair §6.5 measured, so c is directly comparable."""
        assert twin.REFUND_WELLS_EV == (0.0482, twin.E_BIND_ION_EV)

    def test_prediction_bands_are_frozen(self, twin):
        assert twin.RF_P1_MIN_SPREAD == 0.05
        assert twin.RF_P4_ENSEMBLE_C == 0.5308     # committed atlas_ebind_md
        assert twin.RF_P4_PROBE_TOL == 0.10
        assert twin.RF_P5_MIN_RATIO == 1.3
        assert abs(twin.RF_FITTED_EBIND_RATIO - 0.154 / 0.071) < 1e-12


class TestCellArithmetic:
    def test_cell_returns_a_sane_cashed_fraction(self, twin):
        m21 = float(twin.complex_mass_amu(twin.N_STAR))
        row = twin._refund_cell(47.8, 14.2, "center", 0.0, m21)
        assert 0.0 < row["c_mean"] < 1.5
        assert abs(row["c_A"] - row["c_B"]) < 1e-9   # symmetric at r0 = 0
        assert abs(row["refund_frac"] - (1.0 - row["c_mean"])) < 1e-9

    def test_deeper_well_always_costs_kinetic_energy(self, twin):
        m21 = float(twin.complex_mass_amu(twin.N_STAR))
        for R in (9.0, 47.8):
            row = twin._refund_cell(R, 14.2, "center", 0.0, m21)
            assert row["KE_lo_A_eV"] > row["KE_hi_A_eV"]
            assert row["c_mean"] > 0.0

    def test_decomp_bands_frozen(self, twin):
        assert twin.DC_P1_MIN_CASCADE == 0.15
        assert twin.DC_P2_MAX_TRAJ == 1.00
        assert twin.DC_P4_MAX_INTERACTION == 0.10
        assert twin.EBDECOMP_WELLS == (("eb0482", 0.0482),
                                       ("eb1168", twin.E_BIND_ION_EV))

    def test_bin1_mean_ke_matches_the_g3_score_rule(self, twin):
        """Same bin rule as the committed scorer: n == 1, not trapped, and
        a minimum count below which the bin is NaN rather than noisy."""
        n_det = np.array([1] * 25 + [2] * 10)
        trapped = np.zeros(35, dtype=bool)
        v = np.concatenate([np.full(25, 8.0), np.full(10, 3.0)])
        ke, cnt = twin._bin1_mean_ke(n_det, trapped, v)
        assert cnt == 25
        expected = float(twin.kinetic_energy_eV(
            twin.complex_mass_amu(np.array([1])), np.array([8.0]))[0])
        assert ke == pytest.approx(expected)

    def test_bin1_mean_ke_excludes_trapped_and_thin_bins(self, twin):
        n_det = np.array([1] * 25)
        trapped = np.zeros(25, dtype=bool)
        trapped[:10] = True            # 15 left -> below the 20 floor
        ke, cnt = twin._bin1_mean_ke(n_det, trapped, np.full(25, 8.0))
        assert cnt == 15 and np.isnan(ke)

    def test_decomposition_closure_is_exact_by_construction(self, twin):
        """c_traj + c_casc must equal c_total for BOTH orderings — the
        identity the DC-P3 check enforces on the real run."""
        ke_lolo, ke_hihi, ke_lohi, ke_hilo = 0.6735, 0.6371, 0.6290, 0.6810
        d = 0.0686
        c_total = (ke_lolo - ke_hihi) / d
        for traj, casc in (((ke_lohi - ke_hihi) / d, (ke_lolo - ke_lohi) / d),
                           ((ke_lolo - ke_hilo) / d, (ke_hilo - ke_hihi) / d)):
            assert traj + casc == pytest.approx(c_total, abs=1e-12)

    def test_refund_by_slices_one_arm_and_width(self, twin):
        rows = [
            {"arm": "center", "s_rho_A": 14.2, "R_A": 9.0, "c_mean": 0.4},
            {"arm": "center", "s_rho_A": 14.2, "R_A": 47.8, "c_mean": 0.5},
            {"arm": "center", "s_rho_A": 3.14, "R_A": 9.0, "c_mean": 0.9},
            {"arm": "frac027", "s_rho_A": 14.2, "R_A": 9.0, "c_mean": 0.1},
        ]
        assert twin._refund_by(rows, "center", 14.2) == {9.0: 0.4, 47.8: 0.5}
