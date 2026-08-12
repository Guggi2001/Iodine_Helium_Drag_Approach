"""Focused tests for the atlas §6.5 Step-3 E_bind MD cell + its scorer.

Config-construction + band-arithmetic tests only — no run, no trajectory
(the gen_tier2atlas_linring test pattern). These lock the one thing the
experiment's validity rests on: the new cell differs from its committed
CRN partner in **exactly** the well (plus the declared §6.5.1 hatch), and
nothing else.
"""

from __future__ import annotations

import types
import warnings

import numpy as np
import pytest

from scripts.gen_tier2atlas_ebindmd import (
    EBIND_EV,
    EXPECTED_DIFF_KEYS,
    LABEL,
    MDP1_SLOPE_BAND,
    MDP2_SLOPE_BAND,
    MDP3_KE2_MAX_GAP,
    MDP4_TRAP_MAX,
    MDP5_GRADING_BAND,
    MDP6_KE1_CEILING,
    N,
    PARTNER_RUN,
    SEED,
    TWIN_SLOPE_CAPPED,
    WELL_TAG,
    build_shallow_well_cell,
    h405_spec,
    md_run_dir_name,
    verify_partner_pairing,
)
from scripts.post_processing.tier2atlas_ebindmd_table import band_mean_ke_eV


@pytest.fixture(scope="module")
def cfg():
    with warnings.catch_warnings():
        # the two documented pairing warnings (§6.5.1 hatch + §6.6 mass)
        warnings.simplefilter("ignore", RuntimeWarning)
        return build_shallow_well_cell()


class TestCellConstruction:
    def test_inherits_the_committed_h405_pins(self, cfg):
        spec = h405_spec()
        assert (spec.label, spec.a, spec.v_c) == ("h405p", None, 5.5)
        assert (spec.tau_ps, spec.e0_eV) == (4.4, 0.405)
        assert cfg.drag_form == "capped_cubic"
        assert cfg.seed == SEED and N == 500

    def test_well_moved_to_the_shallow_provenance_value(self, cfg):
        assert cfg.binding_energy_I_ion_eV == EBIND_EV[WELL_TAG] == 0.0482

    def test_drag_bundle_stamp_is_held_not_moved(self, cfg):
        """Holding the stamp is what makes this a *declared* pairing
        exception rather than a silently self-consistent re-pairing."""
        stamp = cfg.drag_coefficients.effective_binding_energy_I_ion_eV
        assert stamp == EBIND_EV["eb1168"]
        assert stamp != cfg.binding_energy_I_ion_eV
        assert cfg.allow_unvalidated_binding_pairing is True

    def test_expected_diff_is_exactly_the_well_and_its_hatch(self):
        assert EXPECTED_DIFF_KEYS == {
            "binding_energy_I_ion_eV", "allow_unvalidated_binding_pairing"}

    def test_cfg_differs_from_the_crn_partner_in_exactly_those_keys(self, cfg):
        """The CRN validity condition, checked against the committed run."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            diff = verify_partner_pairing(cfg)
        assert set(diff) == EXPECTED_DIFF_KEYS

    def test_run_dir_is_in_the_atlas_namespace_and_names_the_partner(self):
        name = md_run_dir_name()
        assert LABEL in name and "tier2atlas" in name
        assert "_tier2_" not in name and "tier2probe" not in name
        assert PARTNER_RUN.startswith("9A_drag_shared_pure_cubic_N500_")


class TestFrozenPreRegistration:
    def test_prediction_bands_are_frozen(self):
        assert MDP1_SLOPE_BAND == (-0.60, -0.44)
        assert MDP2_SLOPE_BAND == (-0.55, -0.47)
        assert MDP3_KE2_MAX_GAP == 0.10
        assert MDP4_TRAP_MAX == 0.025
        assert MDP5_GRADING_BAND == (0.05, 0.25)
        assert MDP6_KE1_CEILING == 0.75

    def test_md_p2_brackets_the_twin_forecast_scaled_by_the_lin_transfer(self):
        """MD-P2 is 0.94 x the twin slope (the lin arm's measured MD/twin
        ratio), so the band must contain that value and exclude the raw
        twin slope's neighbourhood only if it truly sits outside."""
        predicted = 0.9415 * TWIN_SLOPE_CAPPED
        assert MDP2_SLOPE_BAND[0] <= predicted <= MDP2_SLOPE_BAND[1]


class TestBandArithmetic:
    @staticmethod
    def _read(pairs):
        n, ke = [], []
        for value, count in pairs:
            n.extend([value[0]] * count)
            ke.extend([value[1]] * count)
        return types.SimpleNamespace(
            n_scored=np.array(n), ke_scored_eV=np.array(ke, dtype=float))

    def test_geometric_and_arithmetic_bands(self):
        read = self._read([((2, 0.8), 30), ((4, 0.2), 30),
                           ((10, 0.3), 30), ((12, 0.1), 30)])
        mid, mid_n = band_mean_ke_eV(read, 2, 8, "geometric")
        deep, deep_n = band_mean_ke_eV(read, 10, 17, "arithmetic")
        assert (mid_n, deep_n) == (2, 2)
        assert mid == pytest.approx(np.sqrt(0.8 * 0.2))
        assert deep == pytest.approx(0.2)

    def test_thin_bins_are_dropped_at_the_committed_min_count(self):
        """A bin with < 20 fragments must not enter the band (the
        committed scorer rule) — otherwise a sparse deep bin would swing
        the slope."""
        read = self._read([((2, 0.8), 30), ((3, 5.0), 19)])
        mid, mid_n = band_mean_ke_eV(read, 2, 8, "geometric")
        assert mid_n == 1
        assert mid == pytest.approx(0.8)

    def test_empty_band_is_nan_not_an_error(self):
        read = self._read([((1, 1.0), 40)])
        val, count = band_mean_ke_eV(read, 10, 17, "arithmetic")
        assert count == 0 and np.isnan(val)
