"""Slice 3 -- SimConfig drag surface, config-load guard, coefficient loader.

All checks are config-construction + loader tests: **no run, no trajectory**
(SLICE3_GOALS_config_and_guard.md §6). Tolerances are exact/tight throughout --
these are declarative and validation checks, not numerical physics.

Covers the §6 acceptance table: inert defaults, existing-preset stability, the
two guard branches (§3.3 dissipativity + §6.5 mass<->coefficient consistency),
the loader's exact-match provenance check (DISTINCT from the guard's ~8 amu
physics band), loader failure modes, and enum completeness.
"""

import json
import typing
import warnings

import pytest

from i2_helium_md import (
    SimConfig,
    REFERENCE_DRAG_ROOT,
    load_drag_coefficients,
    single_pulse_N2000,
    single_pulse_N2000_18Angst,
    single_pulse_N2000_drag,
    single_pulse_N2000_18Angst_drag,
)
from i2_helium_md.config import (
    DragForm,
    DragSpatialGate,
    MassScenario,
    NoiseForm,
    NoiseCalibration,
    NoiseGeometry,
    NoiseLowVBehavior,
    PickupRateForm,
    ValidationHistogramMetric,
    check_drag_config,
    _MASS_COEFFICIENT_CONSISTENCY_TOL_AMU,
)
from i2_helium_md.physics.drag import (
    CAPPED_CUBIC,
    DragCoefficients,
    LINEAR_CUBIC,
    LINEAR_QUADRATIC,
    POWER_LAW,
    PURE_LINEAR,
)

# The frozen 9 A extraction mass; both cases share it (verified on disk).
_M_EFF = 202.953908


def _constant_coeffs(*, a=13.86, b=2.58, m_eff=_M_EFF, binding=0.3) -> DragCoefficients:
    """A constant-mass linear_cubic bundle for guard tests.

    ``binding`` defaults to the ``SimConfig`` default ``binding_energy_I_ion_eV``
    (0.3 eV) so the §6.5.1 pairing arm is satisfied and tests of the OTHER guard
    arms stay focused; pass ``binding=None`` to exercise the unstamped-legacy
    refusal explicitly (TestBindingPairingGuard).
    """
    return DragCoefficients(
        form=LINEAR_CUBIC,
        coefficients={"a": a, "b": b},
        extraction_mass_model="constant",
        extraction_mass_amu=m_eff,
        effective_binding_energy_I_ion_eV=binding,
    )


def _form_coeffs(form: str, coefficients: dict, *, binding=0.3) -> DragCoefficients:
    """A constant-mass bundle of an arbitrary realised form (METHOD_B §10)."""
    return DragCoefficients(
        form=form,
        coefficients=coefficients,
        extraction_mass_model="constant",
        extraction_mass_amu=_M_EFF,
        effective_binding_energy_I_ion_eV=binding,
    )


# ---------------------------------------------------------------------------
# Inert defaults (§2, §6) -- a default config must do nothing surprising
# ---------------------------------------------------------------------------
class TestInertDefaults:
    def test_drag_surface_inert_by_default(self):
        cfg = SimConfig()
        assert cfg.drag_form == "linear_cubic"
        assert cfg.drag_coefficients is None
        assert cfg.mass_scenario == "fixed"
        assert cfg.noise_form == "none"
        # form-selector deferred fields default to their inert member
        assert cfg.drag_spatial_gate == "density_proportional"
        assert cfg.pickup_rate_form == "density_only"
        assert cfg.validation_histogram_metric == "wasserstein"

    def test_mass_initial_defaults_to_m_eff(self):
        cfg = SimConfig()
        assert cfg.mass_initial_amu == cfg.m_eff_amu

    def test_gate_steepness_defaults_to_potential_steepness(self):
        cfg = SimConfig()
        assert cfg.drag_gate_steepness == cfg.potential_steepness == 14.2

    def test_default_config_validates_without_drag_checks(self):
        # drag_coefficients is None -> the guard is a no-op; validate() passes.
        SimConfig().validate()  # must not raise


# ---------------------------------------------------------------------------
# Existing presets unchanged (hard-sphere path bit-identical)
# ---------------------------------------------------------------------------
class TestExistingPresetsUnaffected:
    @pytest.mark.parametrize(
        "preset", [single_pulse_N2000, single_pulse_N2000_18Angst]
    )
    def test_existing_presets_have_no_drag_coefficients(self, preset):
        cfg = preset()
        assert cfg.drag_coefficients is None
        assert cfg.mass_scenario == "fixed"
        assert cfg.noise_form == "none"

    def test_existing_preset_core_fields_unchanged(self):
        cfg = single_pulse_N2000()
        # spot-check the hard-sphere/ensemble fields that Slice 3 must not touch
        assert cfg.hard_sphere_collision_mode == 3
        assert cfg.num_molecules == 2000
        assert cfg.R0_GS_angstrom == 9.0
        assert cfg.mass_attach_probability == 0.09
        cfg.validate()  # collision path config still valid


# ---------------------------------------------------------------------------
# Guard: per-form dissipativity (§3.3, relaxed to a >= 0 per METHOD_B §9.5)
# ---------------------------------------------------------------------------
class TestDissipativityGuard:
    def test_linear_cubic_positive_a_passes(self):
        cfg = SimConfig(drag_coefficients=_constant_coeffs(a=13.86, b=2.58))
        check_drag_config(cfg)  # b > 0 => no real turnover => passes

    def test_linear_cubic_negative_a_refused(self):
        cfg = SimConfig(drag_coefficients=_constant_coeffs(a=-1.0, b=2.58))
        with pytest.raises(ValueError, match="a >= 0"):
            check_drag_config(cfg)

    def test_linear_cubic_zero_a_positive_b_passes(self):
        # METHOD_B §9.5 pure-cubic variant: gamma = g*b*v^2 >= 0 vanishes only
        # at v = 0 where no energy can be added -> strictly dissipative.
        cfg = SimConfig(drag_coefficients=_constant_coeffs(a=0.0, b=2.58))
        check_drag_config(cfg)

    def test_linear_cubic_zero_a_nonpositive_b_refused(self):
        # a = 0 is dissipative ONLY while b > 0; a = b = 0 (no drag at all)
        # and a = 0, b < 0 (energy-adding) are both refused.
        for b in (0.0, -1.0):
            cfg = SimConfig(drag_coefficients=_constant_coeffs(a=0.0, b=b))
            with pytest.raises(ValueError, match="requires b > 0"):
                check_drag_config(cfg)

    def test_turnover_assert_and_skip_b_positive(self):
        # b > 0 => v_dagger = sqrt(-a/b) is imaginary => vacuously dissipative.
        cfg = SimConfig(drag_coefficients=_constant_coeffs(a=1.0, b=5.0))
        check_drag_config(cfg)  # passes with no max-speed ceiling needed


# ---------------------------------------------------------------------------
# Guard: §10.3 dissipativity arms for the form-phase families (METHOD_B §10)
# ---------------------------------------------------------------------------
class TestFormPhaseDissipativityGuard:
    def _cfg(self, form, coefficients):
        return SimConfig(
            drag_form=form,
            drag_coefficients=_form_coeffs(form, coefficients),
        )

    # --- linear_quadratic: a >= 0, c >= 0, a + c > 0 ---
    def test_lq_positive_pair_passes(self):
        check_drag_config(self._cfg(LINEAR_QUADRATIC, {"a": 4.0, "c": 11.0}))

    def test_lq_pure_quadratic_passes(self):
        # a = 0, c > 0: the pure-quadratic variant (the Method-A n~+2
        # hypothesis as a closed form) -- strictly dissipative for v > 0.
        check_drag_config(self._cfg(LINEAR_QUADRATIC, {"a": 0.0, "c": 11.0}))

    def test_lq_pure_linear_passes(self):
        # c = 0, a > 0 degenerates to pure Stokes drag -- still dissipative.
        check_drag_config(self._cfg(LINEAR_QUADRATIC, {"a": 4.0, "c": 0.0}))

    def test_lq_negative_a_refused(self):
        with pytest.raises(ValueError, match="a >= 0 and c >= 0"):
            check_drag_config(self._cfg(LINEAR_QUADRATIC, {"a": -1.0, "c": 11.0}))

    def test_lq_negative_c_refused(self):
        with pytest.raises(ValueError, match="a >= 0 and c >= 0"):
            check_drag_config(self._cfg(LINEAR_QUADRATIC, {"a": 4.0, "c": -1.0}))

    def test_lq_zero_drag_refused(self):
        with pytest.raises(ValueError, match=r"a \+ c > 0"):
            check_drag_config(self._cfg(LINEAR_QUADRATIC, {"a": 0.0, "c": 0.0}))

    # --- power_law: C > 0, n >= 1 ---
    def test_pl_interior_passes(self):
        check_drag_config(self._cfg(POWER_LAW, {"C": 10.36, "n": 2.056}))

    def test_pl_n_equal_one_boundary_passes(self):
        # n = 1 is the locked lower fit bound (§10.4.1); gamma -> g*C at rest.
        check_drag_config(self._cfg(POWER_LAW, {"C": 3.0, "n": 1.0}))

    def test_pl_nonpositive_C_refused(self):
        for C in (0.0, -2.0):
            with pytest.raises(ValueError, match="C > 0"):
                check_drag_config(self._cfg(POWER_LAW, {"C": C, "n": 2.0}))

    def test_pl_n_below_one_refused(self):
        # n < 1 diverges at rest (would need the inert §3.8 floor) -> refused.
        with pytest.raises(ValueError, match="n >= 1"):
            check_drag_config(self._cfg(POWER_LAW, {"C": 10.0, "n": 0.5}))

    # --- pure_linear: a > 0 (the whole dissipativity condition; free-form
    #     linear ring instrument, TIER2_FREEFORM_LINEAR_TWIN_SWEEP_PLAN §6.1) ---
    def test_pure_linear_positive_a_passes(self):
        check_drag_config(self._cfg(PURE_LINEAR, {"a": 27.5}))

    def test_pure_linear_nonpositive_a_refused(self):
        for a in (0.0, -1.0):
            with pytest.raises(ValueError, match="a > 0"):
                check_drag_config(self._cfg(PURE_LINEAR, {"a": a}))


# ---------------------------------------------------------------------------
# Guard: capped_cubic dissipativity + adjudicated-tail restriction
# (Tier-2 Addendum I §I.10 Slice T1: b > 0, v_c > 0; tail range widened to
# -4 <= p_tail <= 0 by the §3.5h low-n-KE-axis adjudication, 2026-07-29 —
# the Step-1c set {0, -1} is a subset)
# ---------------------------------------------------------------------------
class TestCappedCubicGuard:
    def _cfg(self, coefficients):
        return SimConfig(
            drag_form=CAPPED_CUBIC,
            drag_coefficients=_form_coeffs(CAPPED_CUBIC, coefficients),
        )

    @pytest.mark.parametrize("p_tail", [0.0, -1.0, -1.5, -2.0, -3.0, -4.0])
    def test_adjudicated_tails_pass(self, p_tail):
        # {0, -1} = the Step-1c set; the continuous extension to [-4, 0] is
        # the §3.5h low-n-KE-axis adjudication (2026-07-29).
        check_drag_config(
            self._cfg({"b": 2.5153509, "v_c": 6.5, "p_tail": p_tail})
        )

    def test_infinite_vc_passes(self):
        # v_c = inf is the byte-identity limit (== the locked pure cubic
        # everywhere); admissible for regression configs.
        check_drag_config(
            self._cfg({"b": 2.5153509, "v_c": float("inf"), "p_tail": 0.0})
        )

    def test_nonpositive_b_refused(self):
        for b in (0.0, -1.0):
            with pytest.raises(ValueError, match="b > 0"):
                check_drag_config(self._cfg({"b": b, "v_c": 6.5, "p_tail": 0.0}))

    def test_nonpositive_vc_refused(self):
        for v_c in (0.0, -3.0):
            with pytest.raises(ValueError, match="v_c > 0"):
                check_drag_config(
                    self._cfg({"b": 2.5, "v_c": v_c, "p_tail": 0.0})
                )

    @pytest.mark.parametrize("p_tail", [1.0, 0.5, -4.5])
    def test_unadjudicated_tail_exponent_refused(self, p_tail):
        # Positive exponents (tail force growing faster than the cap force)
        # and anything below -4 stay refused: the §3.5h adjudication covers
        # dissipative softening only, -4 <= p_tail <= 0.
        with pytest.raises(ValueError, match="p_tail"):
            check_drag_config(
                self._cfg({"b": 2.5, "v_c": 6.5, "p_tail": p_tail})
            )


# ---------------------------------------------------------------------------
# Guard: drag_form typo-recovery + drag_form<->coeffs.form consistency
# (the runtime recovery for the typo-catching given up by Literal vs Enum)
# ---------------------------------------------------------------------------
class TestDragFormGuard:
    def test_unknown_drag_form_rejected_on_default_config(self):
        # Fires before the drag_coefficients-is-None early return, so a typo is
        # caught even with no coefficients loaded.
        cfg = SimConfig(drag_form="linear_cubicc")  # typo
        assert cfg.drag_coefficients is None
        with pytest.raises(ValueError, match="unknown drag_form"):
            check_drag_config(cfg)

    def test_unknown_drag_form_rejected_via_validate(self):
        cfg = SimConfig(drag_form="not_a_form")
        with pytest.raises(ValueError, match="unknown drag_form"):
            cfg.validate()

    def test_unknown_drag_form_rejected_with_coefficients(self):
        cfg = SimConfig(drag_form="bogus", drag_coefficients=_constant_coeffs())
        with pytest.raises(ValueError, match="unknown drag_form"):
            check_drag_config(cfg)

    def test_drag_form_coeffs_mismatch_rejected(self):
        # coeffs are linear_cubic but the field claims threshold -> mismatch.
        cfg = SimConfig(
            drag_form="threshold", drag_coefficients=_constant_coeffs()
        )
        with pytest.raises(ValueError, match="does not match coefficient form"):
            check_drag_config(cfg)

    def test_valid_matching_drag_form_passes(self):
        cfg = SimConfig(
            drag_form="linear_cubic", drag_coefficients=_constant_coeffs()
        )
        check_drag_config(cfg)  # member + matching coeffs form -> ok


# ---------------------------------------------------------------------------
# Guard: mass <-> coefficient consistency (§6.5)
# ---------------------------------------------------------------------------
class TestMassConsistencyGuard:
    def test_fixed_with_matching_constant_passes(self):
        cfg = SimConfig(
            mass_scenario="fixed",
            m_eff_amu=_M_EFF,
            drag_coefficients=_constant_coeffs(m_eff=_M_EFF),
        )
        check_drag_config(cfg)

    def test_fixed_within_band_passes(self):
        # within the ~8 amu band -> still consistent
        cfg = SimConfig(
            mass_scenario="fixed",
            m_eff_amu=_M_EFF,
            drag_coefficients=_constant_coeffs(m_eff=_M_EFF + 5.0),
        )
        check_drag_config(cfg)

    def test_fixed_beyond_band_refused(self):
        cfg = SimConfig(
            mass_scenario="fixed",
            m_eff_amu=_M_EFF,
            drag_coefficients=_constant_coeffs(m_eff=_M_EFF + 10.0),  # > 8 amu
        )
        with pytest.raises(ValueError, match="inconsistent mass"):
            check_drag_config(cfg)

    def test_allow_inconsistent_downgrades_to_warning(self):
        cfg = SimConfig(
            mass_scenario="fixed",
            m_eff_amu=_M_EFF,
            drag_coefficients=_constant_coeffs(m_eff=_M_EFF + 10.0),
            allow_inconsistent_mass_pairing=True,
        )
        with pytest.warns(RuntimeWarning, match="inconsistent mass"):
            check_drag_config(cfg)

    def test_evolving_scenario_with_constant_coeffs_refused(self):
        cfg = SimConfig(
            mass_scenario="anchored_discrete",
            drag_coefficients=_constant_coeffs(),  # constant, not time_resolved
        )
        with pytest.raises(ValueError, match="time_resolved"):
            check_drag_config(cfg)

    def test_evolving_scenario_inconsistency_can_be_overridden(self):
        cfg = SimConfig(
            mass_scenario="biphasic",
            drag_coefficients=_constant_coeffs(),
            allow_inconsistent_mass_pairing=True,
        )
        with pytest.warns(RuntimeWarning):
            check_drag_config(cfg)


# ---------------------------------------------------------------------------
# Guard: drag <-> binding consistency (§6.5.1) -- the coupled pair
# ---------------------------------------------------------------------------
class TestBindingPairingGuard:
    def test_matching_stamped_binding_passes(self):
        cfg = SimConfig(
            binding_energy_I_ion_eV=0.21,
            drag_coefficients=_constant_coeffs(binding=0.21),
        )
        check_drag_config(cfg)

    def test_mismatched_stamped_binding_refused(self):
        cfg = SimConfig(
            binding_energy_I_ion_eV=0.30,
            drag_coefficients=_constant_coeffs(binding=0.21),
        )
        with pytest.raises(ValueError, match="6.5.1"):
            check_drag_config(cfg)

    def test_unstamped_legacy_bundle_refused(self):
        # A Method-A bundle has no jointly-validated binding: hard refuse.
        cfg = SimConfig(drag_coefficients=_constant_coeffs(binding=None))
        with pytest.raises(ValueError, match="not jointly validated"):
            check_drag_config(cfg)

    def test_escape_hatch_downgrades_to_warning(self):
        cfg = SimConfig(
            drag_coefficients=_constant_coeffs(binding=None),
            allow_unvalidated_binding_pairing=True,
        )
        with pytest.warns(RuntimeWarning, match="6.5.1"):
            check_drag_config(cfg)

    def test_escape_hatch_downgrades_mismatch_too(self):
        cfg = SimConfig(
            binding_energy_I_ion_eV=0.30,
            drag_coefficients=_constant_coeffs(binding=0.21),
            allow_unvalidated_binding_pairing=True,
        )
        with pytest.warns(RuntimeWarning, match="mismatch"):
            check_drag_config(cfg)

    def test_no_coefficients_means_no_binding_check(self):
        # Inert default: no drag -> the binding field is the collision-era
        # parameter and is not guarded.
        cfg = SimConfig(binding_energy_I_ion_eV=0.123)
        check_drag_config(cfg)

    def test_rewired_drag_presets_validate_silently(self):
        # Method-B re-wiring (METHOD_B §10, 2026-06-12): the presets carry the
        # stamped shared_pure_cubic bundle and wire its jointly calibrated
        # binding into binding_energy_I_ion_eV, so the §6.5.1 pairing passes
        # with no escape hatch and no warning (warnings escalate to errors
        # here to catch a regression to the transitional state).
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            cfg = single_pulse_N2000_drag()
            cfg.validate()
        assert cfg.allow_unvalidated_binding_pairing is False


# ---------------------------------------------------------------------------
# Coefficient loader (§4)
# ---------------------------------------------------------------------------
class TestLoader:
    @pytest.mark.parametrize("case", ["9A", "18A"])
    def test_loads_real_case_to_valid_bundle(self, case):
        coeffs = load_drag_coefficients(
            REFERENCE_DRAG_ROOT / case / "linear_and_cubic",
            expected_m_eff_amu=_M_EFF,
        )
        assert coeffs.form == LINEAR_CUBIC
        assert set(coeffs.coefficients) == {"a", "b"}
        assert coeffs.extraction_mass_model == "constant"

    def test_stamps_extraction_mass_from_json(self):
        coeffs = load_drag_coefficients(
            REFERENCE_DRAG_ROOT / "9A" / "linear_and_cubic",
            expected_m_eff_amu=_M_EFF,
        )
        # stamped from the JSON, not the caller's expected value
        assert coeffs.extraction_mass_amu == _M_EFF

    def test_provenance_mismatch_refused(self):
        with pytest.raises(ValueError, match="provenance mismatch"):
            load_drag_coefficients(
                REFERENCE_DRAG_ROOT / "9A" / "linear_and_cubic",
                expected_m_eff_amu=_M_EFF + 1.0,
            )

    def test_missing_file_names_path(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="fit_parameters.json"):
            load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)

    def test_malformed_json_refused(self, tmp_path):
        (tmp_path / "fit_parameters.json").write_text("{not valid json", encoding="utf-8")
        with pytest.raises(ValueError, match="malformed"):
            load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)

    def test_missing_keys_refused(self, tmp_path):
        (tmp_path / "fit_parameters.json").write_text(
            json.dumps({"a": 1.0, "b": 2.0}), encoding="utf-8"  # no errs / meff
        )
        with pytest.raises(ValueError, match="missing required keys"):
            load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)


# ---------------------------------------------------------------------------
# Method-B (trajectory_matching) loader mode -- legacy files load unchanged,
# Method-B files must carry the joint-calibration provenance (METHOD_B §6)
# ---------------------------------------------------------------------------
_LEGACY_JSON = {
    "a": 1.0, "b": 2.0, "a_err": 0.1, "b_err": 0.1, "meff_amu": _M_EFF,
}
_METHOD_B_JSON = {
    **_LEGACY_JSON,
    "extraction_method": "trajectory_matching",
    "extraction_mass_model": "constant",
    "effective_binding_energy_I_ion_eV": 0.21,
    "t_start": 2.67,
    "t_end": 14.0,
    "reference_file": "data/reference/drag/9A/velocity_smoothed/cleaned_data_long.csv",
}


class TestLoaderMethodB:
    def _write(self, tmp_path, payload):
        (tmp_path / "fit_parameters.json").write_text(
            json.dumps(payload), encoding="utf-8"
        )

    def test_legacy_file_loads_as_force_balance_with_no_binding(self, tmp_path):
        self._write(tmp_path, _LEGACY_JSON)
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.extraction_method == "force_balance"
        assert coeffs.effective_binding_energy_I_ion_eV is None

    def test_trajectory_matching_file_stamps_binding(self, tmp_path):
        self._write(tmp_path, _METHOD_B_JSON)
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.extraction_method == "trajectory_matching"
        assert coeffs.effective_binding_energy_I_ion_eV == 0.21
        assert coeffs.extraction_mass_model == "constant"
        assert coeffs.extraction_mass_amu == _M_EFF

    @pytest.mark.parametrize(
        "dropped",
        [
            "extraction_mass_model",
            "effective_binding_energy_I_ion_eV",
            "t_start",
            "t_end",
            "reference_file",
        ],
    )
    def test_trajectory_matching_missing_provenance_refused(self, tmp_path, dropped):
        payload = dict(_METHOD_B_JSON)
        del payload[dropped]
        self._write(tmp_path, payload)
        with pytest.raises(ValueError, match="Method-B keys"):
            load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)

    def test_unknown_extraction_method_refused(self, tmp_path):
        self._write(tmp_path, {**_LEGACY_JSON, "extraction_method": "hand_tuned"})
        with pytest.raises(ValueError, match="unknown\\s+extraction_method"):
            load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)

    def test_trajectory_matching_meff_mismatch_still_refused(self, tmp_path):
        self._write(tmp_path, _METHOD_B_JSON)
        with pytest.raises(ValueError, match="provenance mismatch"):
            load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF + 1.0)


# ---------------------------------------------------------------------------
# Form-generic loader (METHOD_B §10.5): per-form keys, legacy contract intact
# ---------------------------------------------------------------------------
_METHOD_B_PROVENANCE = {
    k: _METHOD_B_JSON[k]
    for k in (
        "extraction_method",
        "extraction_mass_model",
        "effective_binding_energy_I_ion_eV",
        "t_start",
        "t_end",
        "reference_file",
    )
}
_LQ_JSON = {
    "form": "linear_quadratic",
    "a": 4.0, "c": 11.0, "a_err": 0.5, "c_err": 0.8, "meff_amu": _M_EFF,
    **_METHOD_B_PROVENANCE,
}
_PL_JSON = {
    "form": "power_law",
    "C": 10.36, "n": 2.056, "C_err": 1.2, "n_err": 0.3, "meff_amu": _M_EFF,
    **_METHOD_B_PROVENANCE,
}


class TestLoaderFormGeneric:
    def _write(self, tmp_path, payload):
        (tmp_path / "fit_parameters.json").write_text(
            json.dumps(payload), encoding="utf-8"
        )

    def test_linear_quadratic_round_trip(self, tmp_path):
        self._write(tmp_path, _LQ_JSON)
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.form == LINEAR_QUADRATIC
        assert coeffs.coefficients == {"a": 4.0, "c": 11.0}
        assert coeffs.extraction_method == "trajectory_matching"
        assert coeffs.effective_binding_energy_I_ion_eV == 0.21

    def test_power_law_round_trip(self, tmp_path):
        self._write(tmp_path, _PL_JSON)
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.form == POWER_LAW
        assert coeffs.coefficients == {"C": 10.36, "n": 2.056}

    def test_capped_cubic_round_trip(self, tmp_path):
        self._write(
            tmp_path,
            {
                "form": "capped_cubic",
                "b": 2.5153509, "v_c": 7.5, "p_tail": -1.0,
                "b_err": 0.1, "v_c_err": 0.5, "p_tail_err": 0.0,
                "meff_amu": _M_EFF,
            },
        )
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.form == CAPPED_CUBIC
        assert coeffs.coefficients == {"b": 2.5153509, "v_c": 7.5, "p_tail": -1.0}

    def test_explicit_linear_cubic_form_key_accepted(self, tmp_path):
        # A new-style file may also stamp form="linear_cubic" explicitly.
        self._write(tmp_path, {**_LEGACY_JSON, "form": "linear_cubic"})
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.form == LINEAR_CUBIC

    def test_reserved_threshold_form_refused(self, tmp_path):
        self._write(
            tmp_path,
            {"form": "threshold", "F_sat": 1.0, "v0": 1.0, "meff_amu": _M_EFF},
        )
        with pytest.raises(ValueError, match="reserved"):
            load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)

    def test_unknown_form_refused(self, tmp_path):
        self._write(tmp_path, {**_LEGACY_JSON, "form": "exotic"})
        with pytest.raises(ValueError, match="declares form"):
            load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)

    def test_per_form_keys_enforced(self, tmp_path):
        # A linear_quadratic file carrying linear_cubic's {a, b} keys must be
        # refused for the missing {c, c_err}, naming the form.
        bad = dict(_LQ_JSON)
        del bad["c"], bad["c_err"]
        bad["b"], bad["b_err"] = 2.0, 0.1
        self._write(tmp_path, bad)
        with pytest.raises(ValueError, match="linear_quadratic.*missing"):
            load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)

    def test_method_b_provenance_required_for_new_forms_too(self, tmp_path):
        payload = dict(_PL_JSON)
        del payload["effective_binding_energy_I_ion_eV"]
        self._write(tmp_path, payload)
        with pytest.raises(ValueError, match="Method-B keys"):
            load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)


# ---------------------------------------------------------------------------
# Loader-tolerance vs guard-tolerance are DISTINCT (§4.4) -- do not collapse
# ---------------------------------------------------------------------------
class TestSeparateTolerances:
    def test_guard_band_is_eight_amu(self):
        assert _MASS_COEFFICIENT_CONSISTENCY_TOL_AMU == 8.0

    def test_small_mass_diff_passes_guard_but_fails_loader(self, tmp_path):
        # A 0.001 amu difference is WELL within the guard's 8 amu physics band
        # but FAR beyond the loader's exact-match provenance identity.
        delta = 0.001
        # guard: passes (within band)
        cfg = SimConfig(
            mass_scenario="fixed",
            m_eff_amu=_M_EFF,
            drag_coefficients=_constant_coeffs(m_eff=_M_EFF + delta),
        )
        check_drag_config(cfg)
        # loader: refuses the same difference (provenance plumbing identity)
        (tmp_path / "fit_parameters.json").write_text(
            json.dumps(
                {"a": 1.0, "b": 2.0, "a_err": 0.1, "b_err": 0.1, "meff_amu": _M_EFF}
            ),
            encoding="utf-8",
        )
        with pytest.raises(ValueError, match="provenance mismatch"):
            load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF + delta)


# ---------------------------------------------------------------------------
# Drag-enabled presets wire valid, consistent coefficients
# ---------------------------------------------------------------------------
class TestDragPresets:
    @pytest.mark.parametrize(
        "preset,r0",
        [(single_pulse_N2000_drag, 9.0), (single_pulse_N2000_18Angst_drag, 18.0)],
    )
    def test_drag_preset_validates_and_carries_coefficients(self, preset, r0):
        cfg = preset()
        assert cfg.drag_coefficients is not None
        assert cfg.drag_coefficients.form == LINEAR_CUBIC
        # Re-wired (METHOD_B §10): both presets carry the one shared
        # pure-cubic Method-B bundle (a = 0 exactly) with its stamped binding
        # wired into the config -- the §6.5.1 identity holds by construction.
        assert cfg.drag_coefficients.extraction_method == "trajectory_matching"
        assert cfg.drag_coefficients.coefficients["a"] == 0.0
        assert cfg.drag_coefficients.coefficients["b"] > 0.0
        assert cfg.binding_energy_I_ion_eV == (
            cfg.drag_coefficients.effective_binding_energy_I_ion_eV
        )
        assert cfg.m_eff_amu == _M_EFF
        assert cfg.R0_GS_angstrom == r0
        cfg.validate()  # fixed + matching constant coeffs -> consistent

    def test_drag_preset_keeps_hard_sphere_path(self):
        # Slice 3: collision path still runs; hard-sphere fields inherited.
        cfg = single_pulse_N2000_drag()
        base = single_pulse_N2000()
        assert cfg.hard_sphere_collision_mode == base.hard_sphere_collision_mode
        assert cfg.geometric_scattering_crosssection_Iplus == (
            base.geometric_scattering_crosssection_Iplus
        )


# ---------------------------------------------------------------------------
# Enum completeness (§6) -- all members declared on each Literal alias
# ---------------------------------------------------------------------------
class TestEnumCompleteness:
    def test_all_members_present(self):
        assert set(typing.get_args(DragForm)) == {
            "linear_cubic", "linear_quadratic", "threshold", "power_law",
            "capped_cubic", "capped_linear_quadratic", "pure_linear",
        }
        assert set(typing.get_args(DragSpatialGate)) == {
            "density_proportional", "erf_tied", "erf_independent", "sharp"
        }
        assert set(typing.get_args(MassScenario)) == {
            "fixed", "biphasic", "anchored_discrete"
        }
        assert set(typing.get_args(NoiseForm)) == {
            "none", "multiplicative_local_fdt", "empirical_residual"
        }
        assert set(typing.get_args(NoiseCalibration)) == {
            "hard_sphere_variance", "tddft_residual", "strict_fdt_bath"
        }
        assert set(typing.get_args(NoiseGeometry)) == {
            "longitudinal", "isotropic", "anisotropic"
        }
        assert set(typing.get_args(NoiseLowVBehavior)) == {"vanish", "blend_to_isotropic"}
        assert set(typing.get_args(PickupRateForm)) == {
            "density_only", "sweeping", "dwell_time"
        }
        assert set(typing.get_args(ValidationHistogramMetric)) == {
            "wasserstein", "chi2", "ks"
        }
