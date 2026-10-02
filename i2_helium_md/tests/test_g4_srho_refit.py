"""THESIS_FINAL_GOALS G4: Tier-0 refit surface at the sourced density width.

Covers the G4.3 code changes, all without production MD:

1. gate pass-through into the extraction setup (``build_case_setup(cfg_overrides=)``);
2. the escape criterion reads the *potential* width, not the drag-gate width;
3. the ``shared_linear`` joint variant and the ``shared_cubic`` /
   ``shared_linear_cubic`` aliases (stub landscape, as in
   ``test_extraction_trajectory_matching``);
4. bundles stamp the gate / bundle id; the loader round-trips them;
5. the §6.5.2 gate<->bundle guard and the bundle-id guard at config load;
6. the ``drag_bundle`` selector: catalog <-> config enum identity,
   ``build_drag_cfg`` brings a stamped bundle's gate along, export folders;
7. ``cfg_diff_vs_reference`` tolerates the new nested stamps at their default;
8. the driver's pre-registered winner rule (pure function).
"""

from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np
import pytest

from i2_helium_md import SimConfig, single_pulse_N2000_drag
from i2_helium_md import config as config_mod
from i2_helium_md import presets as presets_mod
from i2_helium_md.config import drag_gate_steepness
from i2_helium_md.extraction.trajectory_matching import (
    CaseSetup,
    SHARED_3PARAM,
    SHARED_CUBIC,
    SHARED_LINEAR,
    SHARED_LINEAR_CUBIC,
    SHARED_PURE_CUBIC,
    build_case_setup,
    evaluate_objective,
    fit_shared_trajectory_matching,
    joint_sensitivity_halfwidths,
    write_shared_fit_parameters,
)
from i2_helium_md.physics.drag import DragCoefficients
from i2_helium_md.presets import (
    DRAG_BUNDLE_CATALOG,
    load_drag_bundle,
    load_drag_coefficients,
    resolve_drag_bundle_dir,
)

from .test_extraction_trajectory_matching import (  # type: ignore
    _TRUE,
    _joint_setups,
    _make_stub_run_fn,
)

_M_EFF = 202.953908
_GATE_3P5 = {"drag_spatial_gate": "erf_independent", "drag_gate_steepness": 3.5}


@pytest.fixture
def grid():
    t_grid = np.linspace(2.0, 12.0, 50)
    return t_grid, np.full(50, 3.5)


def _setup(cfg_base, t_grid, ref_speed, case="18A"):
    return CaseSetup(
        case=case, cfg_base=cfg_base, neutral=case,
        t_ref_ps=t_grid, ref_speed_Aps=ref_speed,
        window=(float(t_grid[0]), float(t_grid[-1])),
        a0=_TRUE[0], b0=_TRUE[1], reference_path=Path(f"synthetic_{case}.csv"),
    )


def _stamped_coeffs(**kw) -> DragCoefficients:
    base = dict(
        form="linear_cubic", coefficients={"a": 0.0, "b": 2.5},
        extraction_mass_model="constant", extraction_mass_amu=_M_EFF,
        extraction_method="trajectory_matching",
        effective_binding_energy_I_ion_eV=0.1,
    )
    base.update(kw)
    return DragCoefficients(**base)


def _drag_cfg(coeffs, **kw) -> SimConfig:
    return SimConfig(
        drag_form="linear_cubic", drag_coefficients=coeffs,
        binding_energy_I_ion_eV=coeffs.effective_binding_energy_I_ion_eV, **kw,
    )


# ---------------------------------------------------------------------------
# 1. Gate pass-through (real neutral stage, tiny N)
# ---------------------------------------------------------------------------
class TestGatePassThrough:
    def test_default_is_the_preset_and_overrides_reach_cfg_base(self):
        plain = build_case_setup("9A", num_molecules=2, ion_time_ps=20.0)
        gated = build_case_setup(
            "9A", num_molecules=2, ion_time_ps=20.0, cfg_overrides=_GATE_3P5,
        )
        # Default: cfg_base is the drag preset unchanged.
        preset = single_pulse_N2000_drag(
            num_molecules=2, ion_simulation_time=20.0, dt_ion=0.01, seed=20260604,
        )
        assert plain.cfg_base == preset
        # Overrides reach cfg_base; the potential width is untouched.
        assert drag_gate_steepness(gated.cfg_base) == 3.5
        assert gated.cfg_base.potential_steepness == preset.potential_steepness
        # The placeholder bundle is re-stamped so the gate guard accepts it.
        assert (gated.cfg_base.drag_coefficients.extraction_gate_steepness_angstrom
                == 3.5)
        # The gate is an ion-stage quantity: the neutral stage is identical.
        import dataclasses

        compared = 0
        for f in dataclasses.fields(plain.neutral):
            a, b = getattr(plain.neutral, f.name), getattr(gated.neutral, f.name)
            if isinstance(a, np.ndarray):
                np.testing.assert_array_equal(a, b, err_msg=f.name)
                compared += 1
            else:
                assert a == b, f.name
        assert compared >= 2


# ---------------------------------------------------------------------------
# 2. Escape criterion on the potential width
# ---------------------------------------------------------------------------
class TestEscapeOnPotentialWidth:
    def test_ion_outside_sharp_density_gate_but_inside_well_is_not_escaped(
        self, grid,
    ):
        # The stub parks escaped ions at depth 3*drag_gate_steepness: 10.5 A at
        # a 3.5 A gate -- outside the density gate, inside 2*14.2 A of well.
        t_grid, ref = grid
        setup = _setup(SimConfig(**_GATE_3P5), t_grid, ref)
        res = evaluate_objective(_TRUE, setup, run_fn=_make_stub_run_fn(t_grid, ref))
        assert res.escape_fraction == 0.0
        assert res.penalty_Aps > 0.0

    def test_legacy_gate_unchanged(self, grid):
        t_grid, ref = grid
        setup = _setup(SimConfig(), t_grid, ref)
        res = evaluate_objective(_TRUE, setup, run_fn=_make_stub_run_fn(t_grid, ref))
        assert res.escape_fraction == 1.0

    def test_candidate_bundle_stamped_with_run_gate(self, grid):
        t_grid, ref = grid
        seen = {}
        stub = _make_stub_run_fn(t_grid, ref)

        def run_fn(cfg, neutral):
            seen["stamp"] = cfg.drag_coefficients.extraction_gate_steepness_angstrom
            return stub(cfg, neutral)

        evaluate_objective(_TRUE, _setup(SimConfig(**_GATE_3P5), t_grid, ref),
                           run_fn=run_fn)
        assert seen["stamp"] == 3.5


# ---------------------------------------------------------------------------
# 3. shared_linear variant + aliases
# ---------------------------------------------------------------------------
class TestFamilies:
    def test_shared_linear_fixes_b_to_zero_and_recovers_a(self, grid):
        t_grid, ref = grid
        true = (14.0, 0.0, 0.20)
        setups, run_fn = _joint_setups(t_grid, ref, {"18A": true, "9A": true})
        fit = fit_shared_trajectory_matching(
            setups, variant=SHARED_LINEAR, anchors=(14.0, 2.0),
            run_fn=run_fn, maxfev_per_start=400, extra_starts=False,
            prescan_e_bind_eV=[0.15, 0.2, 0.25],
        )
        assert fit.variant == SHARED_LINEAR
        assert fit.best.b == 0.0  # fixed, not fitted
        assert fit.best.a == pytest.approx(true[0], rel=0.01)
        assert set(fit.best.e_bind_eV.values()).pop() == pytest.approx(true[2], abs=0.005)
        # b fixed at 0 -> half-width 0 by convention.
        _, hw_b, _ = joint_sensitivity_halfwidths(setups, fit.best, run_fn=run_fn)
        assert hw_b == 0.0

    @pytest.mark.parametrize(
        "alias, engine",
        [(SHARED_CUBIC, SHARED_PURE_CUBIC), (SHARED_LINEAR_CUBIC, SHARED_3PARAM)],
    )
    def test_alias_runs_engine_unchanged_but_keeps_its_name(self, grid, alias, engine):
        t_grid, ref = grid
        setups, run_fn = _joint_setups(t_grid, ref, {"18A": _TRUE, "9A": _TRUE})
        kw = dict(anchors=(14.0, 2.0), run_fn=run_fn, maxfev_per_start=60,
                  extra_starts=True, prescan_e_bind_eV=[0.15, 0.2, 0.25])
        fa = fit_shared_trajectory_matching(setups, variant=alias, **kw)
        fe = fit_shared_trajectory_matching(setups, variant=engine, **kw)
        assert fa.variant == alias and fe.variant == engine
        assert fa.n_evaluations == fe.n_evaluations
        for ra, re_ in zip(fa.starts, fe.starts):
            assert (ra.a, ra.b, ra.e_bind_eV, ra.objective) == (
                re_.a, re_.b, re_.e_bind_eV, re_.objective)

    def test_gate_provenance_recorded_and_mismatch_refused(self, grid):
        t_grid, ref = grid
        setups, run_fn = _joint_setups(t_grid, ref, {"18A": _TRUE, "9A": _TRUE})
        for s in setups.values():
            s.cfg_base = SimConfig(**_GATE_3P5)
        fit = fit_shared_trajectory_matching(
            setups, variant=SHARED_CUBIC, anchors=(14.0, 2.0), run_fn=run_fn,
            maxfev_per_start=5, extra_starts=False, prescan_e_bind_eV=[0.2],
        )
        assert (fit.drag_spatial_gate, fit.gate_steepness_angstrom,
                fit.potential_steepness_angstrom) == ("erf_independent", 3.5, 14.2)
        setups["9A"].cfg_base = SimConfig()
        with pytest.raises(ValueError, match="spatial gate"):
            fit_shared_trajectory_matching(
                setups, variant=SHARED_CUBIC, anchors=(14.0, 2.0), run_fn=run_fn,
            )


# ---------------------------------------------------------------------------
# 4. Bundle stamps: writer -> loader round trip
# ---------------------------------------------------------------------------
class TestBundleStamps:
    def test_writer_stamps_and_loader_round_trips(self, grid, tmp_path):
        t_grid, ref = grid
        setups, stub = _joint_setups(t_grid, ref, {"18A": _TRUE, "9A": _TRUE})
        for s in setups.values():
            s.cfg_base = SimConfig(**_GATE_3P5)

        def run_fn(cfg, neutral):
            # The stub parks ions at 3 x gate width (10.5 A here), inside the
            # 14.2 A well; move them past 2 x potential width so they escape.
            ion = stub(cfg, neutral)
            ion.positions_x = ion.positions_x + 3.0 * cfg.potential_steepness
            return ion

        fit = fit_shared_trajectory_matching(
            setups, variant=SHARED_CUBIC, anchors=(14.0, 2.0), run_fn=run_fn,
            maxfev_per_start=5, extra_starts=False, prescan_e_bind_eV=[0.2],
        )
        fit.a_err, fit.b_err, fit.e_bind_err_eV, fit.uncertainty_model = 0, 0.1, 0.01, "x"
        out = write_shared_fit_parameters(
            fit, tmp_path, bundle_id="srho3p5A_shared_cubic",
            density_width_source="src", held_out_axis="none (stage-2 only)",
        )
        raw = json.loads(out.read_text(encoding="utf-8"))
        assert raw["variant"] == SHARED_CUBIC
        assert raw["drag_spatial_gate"] == "erf_independent"
        assert raw["drag_gate_steepness_angstrom"] == 3.5
        assert raw["potential_steepness_angstrom"] == 14.2
        assert raw["bundle_id"] == "srho3p5A_shared_cubic"
        assert raw["held_out_axis"] == "none (stage-2 only)"
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.extraction_gate_steepness_angstrom == 3.5
        assert coeffs.bundle_id == "srho3p5A_shared_cubic"

    def test_legacy_bundle_loads_unstamped(self):
        coeffs = load_drag_bundle("shared_pure_cubic", "9A", expected_m_eff_amu=_M_EFF)
        assert coeffs.extraction_gate_steepness_angstrom is None
        assert coeffs.bundle_id is None

    @pytest.mark.parametrize("bad", [0.0, -1.0])
    def test_non_positive_gate_stamp_refused(self, bad):
        with pytest.raises(ValueError, match="extraction_gate_steepness"):
            _stamped_coeffs(extraction_gate_steepness_angstrom=bad)

    def test_empty_bundle_id_refused(self):
        with pytest.raises(ValueError, match="bundle_id"):
            _stamped_coeffs(bundle_id="")


# ---------------------------------------------------------------------------
# 5. Config-load guards
# ---------------------------------------------------------------------------
class TestGateGuard:
    def test_legacy_bundle_at_legacy_gate_passes(self):
        single_pulse_N2000_drag(num_molecules=2).validate()

    def test_legacy_bundle_reads_as_14p2_and_refuses_3p5(self):
        cfg = single_pulse_N2000_drag(num_molecules=2, **_GATE_3P5)
        with pytest.raises(ValueError, match="6.5.2"):
            cfg.validate()

    def test_hatch_downgrades_to_warning(self):
        cfg = single_pulse_N2000_drag(
            num_molecules=2, allow_inconsistent_gate_pairing=True, **_GATE_3P5,
        )
        with pytest.warns(RuntimeWarning, match="6.5.2"):
            cfg.validate()

    def test_matching_stamp_passes_and_mismatch_refused(self):
        coeffs = _stamped_coeffs(extraction_gate_steepness_angstrom=3.5)
        _drag_cfg(coeffs, **_GATE_3P5).validate()
        with pytest.raises(ValueError, match="6.5.2"):
            _drag_cfg(coeffs).validate()  # default gate resolves to 14.2

    def test_resolver_reexported_from_ion(self):
        from i2_helium_md.simulation import ion

        assert ion.drag_gate_steepness is config_mod.drag_gate_steepness


class TestBundleIdGuard:
    def test_matching_id_passes(self):
        coeffs = _stamped_coeffs(extraction_gate_steepness_angstrom=3.5,
                                 bundle_id="srho3p5A_shared_cubic")
        _drag_cfg(coeffs, drag_bundle="srho3p5A_shared_cubic", **_GATE_3P5).validate()

    def test_mismatched_id_refused(self):
        coeffs = _stamped_coeffs(extraction_gate_steepness_angstrom=3.5,
                                 bundle_id="srho3p5A_shared_cubic")
        with pytest.raises(ValueError, match="bundle_id"):
            _drag_cfg(coeffs, **_GATE_3P5).validate()  # drag_bundle = legacy default

    def test_unstamped_coefficients_with_new_selector_refused(self):
        with pytest.raises(ValueError, match="no bundle_id stamp"):
            _drag_cfg(_stamped_coeffs(), drag_bundle="srho3p5A_shared_cubic").validate()

    def test_unknown_selector_refused_even_without_drag(self):
        with pytest.raises(ValueError, match="drag_bundle"):
            SimConfig(drag_bundle="nope").validate()


# ---------------------------------------------------------------------------
# 6. Selector resolution
# ---------------------------------------------------------------------------
class TestSelector:
    def test_catalog_matches_config_enum(self):
        assert set(DRAG_BUNDLE_CATALOG) == set(config_mod._KNOWN_DRAG_BUNDLES)
        assert set(config_mod._UNSTAMPED_LEGACY_DRAG_BUNDLES) < set(DRAG_BUNDLE_CATALOG)

    @pytest.mark.parametrize("bundle_id", sorted(DRAG_BUNDLE_CATALOG))
    @pytest.mark.parametrize("case", ["9A", "18A"])
    def test_every_id_resolves_with_geometry_guard(self, bundle_id, case):
        path = resolve_drag_bundle_dir(bundle_id, case)
        other = "18A" if case == "9A" else "9A"
        if "{case}" in DRAG_BUNDLE_CATALOG[bundle_id]:
            assert case in path.parts and other not in path.parts

    def test_unknown_id_and_case_refused(self):
        with pytest.raises(ValueError, match="drag bundle"):
            resolve_drag_bundle_dir("nope", "9A")
        with pytest.raises(ValueError, match="case"):
            resolve_drag_bundle_dir("shared_pure_cubic", "12A")

    def _write_stamped_bundle(self, root: Path, bundle_id: str, stamp_id: str):
        legacy = json.loads(
            (presets_mod.REFERENCE_DRAG_ROOT / DRAG_BUNDLE_CATALOG["shared_pure_cubic"]
             / "fit_parameters.json").read_text(encoding="utf-8")
        )
        legacy.update(drag_gate_steepness_angstrom=3.5, bundle_id=stamp_id)
        out = root / DRAG_BUNDLE_CATALOG[bundle_id]
        out.mkdir(parents=True)
        (out / "fit_parameters.json").write_text(json.dumps(legacy), encoding="utf-8")

    def test_build_drag_cfg_brings_the_stamped_gate(self, tmp_path, monkeypatch):
        from scripts.tier0_common import build_drag_cfg

        bid = "srho3p5A_shared_cubic"
        self._write_stamped_bundle(tmp_path, bid, bid)
        monkeypatch.setattr(presets_mod, "REFERENCE_DRAG_ROOT", tmp_path)
        cfg = build_drag_cfg("9A", bid, num_molecules=2, ion_time_ps=1.0,
                             dt_ion_ps=0.01, seed=1)
        assert cfg.drag_bundle == bid
        assert cfg.drag_spatial_gate == "erf_independent"
        assert drag_gate_steepness(cfg) == 3.5
        assert cfg.potential_steepness == 14.2
        cfg.validate()

    def test_build_drag_cfg_legacy_unchanged(self):
        from scripts.tier0_common import build_drag_cfg

        cfg = build_drag_cfg("9A", "shared_pure_cubic", num_molecules=2,
                             ion_time_ps=1.0, dt_ion_ps=0.01, seed=1)
        assert cfg.drag_bundle == "shared_pure_cubic"
        assert cfg.drag_spatial_gate == SimConfig().drag_spatial_gate
        cfg.validate()

    def test_misplaced_bundle_refused(self, tmp_path, monkeypatch):
        self._write_stamped_bundle(tmp_path, "srho3p5A_shared_cubic",
                                   "srho3p5A_shared_linear")
        monkeypatch.setattr(presets_mod, "REFERENCE_DRAG_ROOT", tmp_path)
        with pytest.raises(ValueError, match="stamped"):
            load_drag_bundle("srho3p5A_shared_cubic", "9A", expected_m_eff_amu=_M_EFF)

    def test_export_folder_names(self):
        from scripts.tier0_common import reference_export_dir_name

        assert reference_export_dir_name("tier0", "shared_pure_cubic") == "tier0"
        assert reference_export_dir_name("tier1a", "percase_lq") == "tier1a"
        assert (reference_export_dir_name("tier0", "srho3p5A_shared_cubic")
                == "tier0_srho3p5A")
        with pytest.raises(ValueError):
            reference_export_dir_name("tier0", "nope")


# ---------------------------------------------------------------------------
# 7. Reference diff tolerates the new nested stamps at default
# ---------------------------------------------------------------------------
class TestCfgDiffNested:
    def _legacy_reference(self, cfg, tmp_path) -> Path:
        ref = json.loads(json.dumps(asdict(cfg)))
        for key in ("extraction_gate_steepness_angstrom", "bundle_id"):
            ref["drag_coefficients"].pop(key)
        for key in ("drag_bundle", "allow_inconsistent_gate_pairing"):
            ref.pop(key)
        path = tmp_path / "cfg.json"
        path.write_text(json.dumps(ref), encoding="utf-8")
        return path

    def test_pre_g4_reference_gives_no_spurious_diff(self, tmp_path):
        from scripts.tier2_common import cfg_diff_vs_reference

        cfg = single_pulse_N2000_drag(num_molecules=2)
        assert cfg_diff_vs_reference(cfg, self._legacy_reference(cfg, tmp_path)) == []

    def test_legacy_selector_is_a_label_stamped_selector_is_strict(self, tmp_path):
        from scripts.tier2_common import cfg_diff_vs_reference

        cfg = single_pulse_N2000_drag(num_molecules=2)
        ref_path = self._legacy_reference(cfg, tmp_path)
        # Unstamped legacy id: carries no information beyond drag_coefficients.
        assert cfg_diff_vs_reference(
            replace(cfg, drag_bundle="shared_lq"), ref_path) == []
        # Stamped (post-G4) id: not witnessable by a pre-G4 reference.
        with pytest.raises(AssertionError, match="drag_bundle"):
            cfg_diff_vs_reference(
                replace(cfg, drag_bundle="srho3p5A_shared_cubic"), ref_path)

    def test_non_default_nested_stamp_refused(self, tmp_path):
        from scripts.tier2_common import cfg_diff_vs_reference

        cfg = single_pulse_N2000_drag(num_molecules=2)
        ref_path = self._legacy_reference(cfg, tmp_path)
        stamped = replace(cfg, drag_coefficients=replace(
            cfg.drag_coefficients, extraction_gate_steepness_angstrom=3.5))
        with pytest.raises(AssertionError, match="extraction_gate_steepness"):
            cfg_diff_vs_reference(stamped, ref_path)


# ---------------------------------------------------------------------------
# 8. The pre-registered winner rule (driver, pure function)
# ---------------------------------------------------------------------------
class TestWinnerRule:
    @staticmethod
    def _sel(obj_lin, obj_cub, obj_lc, passes=(True, True, True)):
        from scripts.extraction.method_b_srho3p5_refit import select_winner

        objectives = {SHARED_LINEAR: obj_lin, SHARED_CUBIC: obj_cub,
                      SHARED_LINEAR_CUBIC: obj_lc}
        return select_winner(objectives, dict(zip(
            (SHARED_LINEAR, SHARED_CUBIC, SHARED_LINEAR_CUBIC), passes)))

    def test_both_equivalent_cubic_wins(self):
        assert self._sel(0.104, 0.103, 0.100)["winner"] == SHARED_CUBIC

    def test_only_cubic_equivalent(self):
        assert self._sel(0.200, 0.104, 0.100)["winner"] == SHARED_CUBIC

    def test_only_linear_equivalent(self):
        assert self._sel(0.104, 0.200, 0.100)["winner"] == SHARED_LINEAR

    def test_neither_equivalent_full_law(self):
        assert self._sel(0.200, 0.200, 0.100)["winner"] == SHARED_LINEAR_CUBIC

    def test_threshold_is_inclusive_at_0p005(self):
        sel = self._sel(0.300, 0.1050, 0.1000)
        assert sel["T_a0_Aps"] == pytest.approx(0.005)
        # 0.105 - 0.100 is 0.005 up to float rounding; inclusive by rule.
        assert sel["winner"] in (SHARED_CUBIC, SHARED_LINEAR_CUBIC)

    def test_equivalent_but_band_failing_family_is_out(self):
        sel = self._sel(0.200, 0.101, 0.100, passes=(True, False, True))
        assert sel["winner"] == SHARED_LINEAR_CUBIC

    def test_negative_difference_noted_and_counts_equivalent(self):
        sel = self._sel(0.300, 0.099, 0.100)
        assert sel["winner"] == SHARED_CUBIC
        assert sel["T_a0_Aps"] < 0 and sel["notes"]

    def test_no_admissible_family(self):
        sel = self._sel(0.3, 0.3, 0.1, passes=(False, False, False))
        assert sel["winner"] is None
