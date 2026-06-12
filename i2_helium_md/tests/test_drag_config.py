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
    MassRateForm,
    ValidationHistogramMetric,
    check_drag_config,
    _MASS_COEFFICIENT_CONSISTENCY_TOL_AMU,
)
from i2_helium_md.physics.drag import DragCoefficients, LINEAR_CUBIC

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
        assert cfg.mass_rate_form == "density_only"
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
            mass_scenario="scenario_B_stripping",
            drag_coefficients=_constant_coeffs(),  # constant, not time_resolved
        )
        with pytest.raises(ValueError, match="time_resolved"):
            check_drag_config(cfg)

    def test_evolving_scenario_inconsistency_can_be_overridden(self):
        cfg = SimConfig(
            mass_scenario="scenario_A_accretion",
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

    def test_transitional_drag_presets_warn_not_refuse(self):
        # Until the Method-B re-wiring, the drag presets pair a legacy bundle
        # with a hand-set binding under the documented escape hatch: loud
        # warning, not refusal. The re-wiring slice flips them to stamped
        # bundles and removes the hatch (and this test's expectation).
        with pytest.warns(RuntimeWarning, match="6.5.1"):
            cfg = single_pulse_N2000_drag()
            cfg.validate()


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
            "linear_cubic", "linear_quadratic", "threshold", "power_law"
        }
        assert set(typing.get_args(DragSpatialGate)) == {
            "density_proportional", "erf_tied", "erf_independent", "sharp"
        }
        assert set(typing.get_args(MassScenario)) == {
            "fixed", "scenario_A_accretion", "scenario_B_stripping", "biphasic"
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
        assert set(typing.get_args(MassRateForm)) == {
            "density_only", "sweeping", "dwell_time"
        }
        assert set(typing.get_args(ValidationHistogramMetric)) == {
            "wasserstein", "chi2", "ks"
        }
