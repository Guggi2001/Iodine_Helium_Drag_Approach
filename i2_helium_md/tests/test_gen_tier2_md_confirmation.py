"""Coverage for the Slice-T3 MD-confirmation generator (plan §I.10).

Harness tests per ``TIER2_STAIRCASE_PROBE_PLAN.md`` §I.10 Slice T3: the frozen
C1-C4 matrix builds at production kinematics (R0_GS = 2.666 A,
E_coulomb_scale = 1.0, 2.70 eV/fragment) with the common pins (s_eff = 8,
density_scaled cooling, off-center births, relaxation 1000 ps + detection at
the Sourced 8.53 us), the rq4graded / floor1 rung tables are constructed
generator-side from ``d0_of_n`` x the taper multipliers (no taper physics in
the package -- the Slice-T2 boundary), the conf dirs stay out of the F3
campaign scoreboard, and a tiny-N end-to-end run writes all five artifacts
(cfg.json, neutral, ion, relaxation, detection). No production scale, no
figures.
"""

from __future__ import annotations

from dataclasses import replace
from fnmatch import fnmatch

import pytest

from scripts import gen_tier2_md_confirmation as script
from scripts.post_processing import tier2_size_distribution_table as f3
from scripts.tier2_common import (
    EXPERIMENTAL_RELAXATION_TIME_PS,
    PRODUCTION_BUDGET_EV,
    build_biphasic_cfg,
    tier2_confirmation_run_dir_name,
)
from i2_helium_md.physics.constants import N_STAR
from i2_helium_md.physics.dissociation_ladder import d0_of_n, ladder_cumsum


# ---------------------------------------------------------------------------
# Ladder tables (generator-side construction -- the T2 boundary)
# ---------------------------------------------------------------------------


def test_form_u_rungs_match_d0_of_n_and_pin_sigma21():
    """The base table is the mixture kappa=1 Form-U ladder read through
    d0_of_n(1..N_STAR+11) -- the Form-U cache-height convention (Slice T2).
    Numerical anchor: Sigma(21) = 0.18783720 eV (the probe-program wiring
    oracle, findings §4j)."""
    rungs = script.form_u_rungs_eV()
    assert len(rungs) == N_STAR + 11 == 32
    for i, rung in enumerate(rungs, start=1):
        assert rung == pytest.approx(
            d0_of_n(i, picture="statistical_mixture", kappa=1.0)
        )
    assert sum(rungs[:N_STAR]) == pytest.approx(0.18783720, abs=1e-7)
    assert sum(rungs[:N_STAR]) == pytest.approx(
        ladder_cumsum(N_STAR, picture="statistical_mixture", kappa=1.0)
    )


def test_rq4graded_table_tapers_first_three_rungs():
    """rq4graded = the RQ4-graded diagnostic taper 2.2 : 1.5 : 1.3 on rungs
    1-3 of the flat mixture bottom, rest = Form-U verbatim (findings I51)."""
    base = script.form_u_rungs_eV()
    graded = script.rq4graded_rungs_eV()
    assert len(graded) == len(base)
    assert script.RQ4GRADED_TAPER_MULTIPLIERS == (2.2, 1.5, 1.3)
    for i, mult in enumerate(script.RQ4GRADED_TAPER_MULTIPLIERS):
        assert graded[i] == pytest.approx(mult * base[i])
    assert graded[3:] == base[3:]


def test_floor1_table_sets_rung1_absolute():
    """floor1 = the knob-free X2/3Pi transition-at-n=1 variant: rung 1 is the
    [IHe05] 13.3 meV first-rung depth (absolute, not a multiplier), rest =
    Form-U verbatim (plan L3 / findings §4j)."""
    base = script.form_u_rungs_eV()
    floor1 = script.floor1_rungs_eV()
    assert len(floor1) == len(base)
    assert floor1[0] == pytest.approx(script.FLOOR1_RUNG1_EV) == pytest.approx(0.0133)
    assert floor1[1:] == base[1:]


def test_ladder_key_dispatch():
    """'form_u' rides the config default (byte-inert: no selector, no table);
    the tabulated keys return ('tabulated', <their table>); unknown keys are
    refused loudly."""
    assert script.ladder_selector_and_table("form_u") == (None, None)
    sel, table = script.ladder_selector_and_table("rq4graded")
    assert sel == "tabulated" and table == script.rq4graded_rungs_eV()
    sel, table = script.ladder_selector_and_table("floor1")
    assert sel == "tabulated" and table == script.floor1_rungs_eV()
    with pytest.raises(ValueError, match="ladder"):
        script.ladder_selector_and_table("slid2")


# ---------------------------------------------------------------------------
# The frozen C1-C4 matrix
# ---------------------------------------------------------------------------


def test_matrix_is_the_frozen_c1_c4():
    """S2-D3: 2 targets + 2 controls, exactly as frozen in plan §I.10."""
    specs = {spec.label: spec for spec in script.CONFIRMATION_MATRIX}
    assert list(specs) == ["c1", "c2", "c3", "c4"]

    assert specs["c1"].drag_form == "capped_cubic"
    assert specs["c1"].drag_coefficient_overrides == {"v_c": 7.5, "p_tail": -1.0}
    assert specs["c1"].tau_ps == 3.8
    assert specs["c1"].ladder_key == "rq4graded"
    assert specs["c1"].e0_eV == 0.25

    assert specs["c2"].drag_form == "capped_cubic"
    assert specs["c2"].drag_coefficient_overrides == {"v_c": 6.5, "p_tail": 0.0}
    assert specs["c2"].tau_ps == 4.0
    assert specs["c2"].ladder_key == "rq4graded"
    assert specs["c2"].e0_eV == 0.24

    # C3 = baseline control: the current law and the flat Form-U ladder.
    assert specs["c3"].drag_form is None
    assert specs["c3"].drag_coefficient_overrides is None
    assert specs["c3"].tau_ps == 6.55
    assert specs["c3"].ladder_key == "form_u"
    assert specs["c3"].e0_eV == 0.25

    assert specs["c4"].drag_form == "capped_cubic"
    assert specs["c4"].drag_coefficient_overrides == {"v_c": 7.5, "p_tail": -1.0}
    assert specs["c4"].tau_ps == 4.4
    assert specs["c4"].ladder_key == "floor1"
    assert specs["c4"].e0_eV == 0.23


def test_build_confirmation_stamps_common_pins(tmp_path, monkeypatch):
    """Every C-config carries the §I.10 common pins: production kinematics,
    the 2.70 eV budget with f_int = E0/2.70 stamped EXACTLY (never rounded),
    s_eff = 8, gated cooling, off-center births, 1000 ps relaxation cap, and
    the Sourced detector time. Pinned at leg 'a' (the frozen T3 record; the
    A' leg's 3000 ps adequacy cap is asserted in the leg-A' test)."""
    monkeypatch.setattr(script, "LEG", "a")
    scheduled = script.build_confirmation(tmp_path)
    assert len(scheduled) == 4
    for (label, cfg, run_dir), spec in zip(scheduled, script.CONFIRMATION_MATRIX):
        assert spec.label in label
        cfg.validate()
        assert cfg.mass_scenario == "biphasic"
        assert cfg.coulomb_available_eV == PRODUCTION_BUDGET_EV == 2.70
        assert cfg.R0_GS_angstrom == pytest.approx(2.666)
        assert cfg.E_coulomb_scale == pytest.approx(1.0)
        assert cfg.single_initial_position is False
        assert cfg.evap_rrk_dof == pytest.approx(8.0)
        assert cfg.cooling_spatial_gate == "density_scaled"
        assert cfg.pickup_rate_coefficient == pytest.approx(0.9)
        assert cfg.internal_energy_retained_fraction == pytest.approx(0.1)
        # f_int is the exact quotient -- 0.0926 (c1/c3), 0.0889, 0.0852 never
        # collapse onto a two-decimal grid (the C1/C2 aliasing hazard).
        assert cfg.internal_energy_partition_fraction == spec.e0_eV / 2.70
        assert cfg.internal_energy_cooling_tau_ps == pytest.approx(spec.tau_ps)
        assert cfg.relaxation_stage_enabled is True
        assert cfg.relaxation_time_ps == pytest.approx(1000.0)
        assert cfg.detection_stage_enabled is True
        assert cfg.detection_time_ps == pytest.approx(EXPERIMENTAL_RELAXATION_TIME_PS)
        assert cfg.num_molecules == script.N
        assert cfg.seed == script.SEED


def test_build_confirmation_per_config_knobs(tmp_path):
    """Drag form/coefficients and ladder per spec; C3 is byte-inert on both
    axes (current linear_cubic bundle untouched, Form-U ladder, no table)."""
    by_label = {
        label.split()[-1]: cfg for label, cfg, _ in script.build_confirmation(tmp_path)
    }
    control = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        num_molecules=script.N,
        ion_time_ps=script.ION_TIME_PS,
        dt_ion_ps=script.DT_ION_PS,
        seed=script.SEED,
    )

    for lab in ("c1", "c4"):
        assert by_label[lab].drag_form == "capped_cubic"
        assert by_label[lab].drag_coefficients.coefficients["v_c"] == pytest.approx(7.5)
        assert by_label[lab].drag_coefficients.coefficients["p_tail"] == pytest.approx(-1.0)
    assert by_label["c2"].drag_form == "capped_cubic"
    assert by_label["c2"].drag_coefficients.coefficients["v_c"] == pytest.approx(6.5)
    assert by_label["c2"].drag_coefficients.coefficients["p_tail"] == pytest.approx(0.0)
    # b rides the locked shared_pure_cubic bundle on every capped arm.
    for lab in ("c1", "c2", "c4"):
        assert (
            by_label[lab].drag_coefficients.coefficients["b"]
            == control.drag_coefficients.coefficients["b"]
        )

    assert by_label["c3"].drag_form == "linear_cubic"
    assert by_label["c3"].drag_coefficients == control.drag_coefficients
    assert by_label["c3"].dissociation_ladder == "form_u"
    assert by_label["c3"].tabulated_ladder_rungs_eV is None

    assert by_label["c1"].dissociation_ladder == "tabulated"
    assert by_label["c1"].tabulated_ladder_rungs_eV == pytest.approx(
        script.rq4graded_rungs_eV()
    )
    assert by_label["c2"].tabulated_ladder_rungs_eV == pytest.approx(
        script.rq4graded_rungs_eV()
    )
    assert by_label["c4"].dissociation_ladder == "tabulated"
    assert by_label["c4"].tabulated_ladder_rungs_eV == pytest.approx(
        script.floor1_rungs_eV()
    )


def test_build_confirmation_refuses_non_production_budget(tmp_path, monkeypatch):
    """The confirmation build exists only at the production channel (S2-D
    scope); mirroring the probe generator's budget guard in the opposite
    direction."""
    monkeypatch.setattr(script, "BUDGET_EV", 0.80)
    with pytest.raises(ValueError, match="BUDGET_EV"):
        script.build_confirmation(tmp_path)


# ---------------------------------------------------------------------------
# Namespace / run dirs
# ---------------------------------------------------------------------------


def test_run_dirs_unique_and_in_conf_namespace(tmp_path, monkeypatch):
    """Leg A (the delivered T3 configuration) keeps its exact run-dir names —
    the byte-identity lock for the four on-disk pilot dirs."""
    monkeypatch.setattr(script, "LEG", "a")
    dirs = [run_dir for _, _, run_dir in script.build_confirmation(tmp_path)]
    assert len({d.name for d in dirs}) == 4
    for d, spec in zip(dirs, script.CONFIRMATION_MATRIX):
        assert d.name == tier2_confirmation_run_dir_name(
            script.CASE,
            script.VARIANT,
            script.N,
            config_label=spec.label,
            budget_eV=script.BUDGET_EV,
        )
        assert "_tier2probe_conf270_" in d.name
        assert not fnmatch(d.name, "*_tier2_*")


# ---------------------------------------------------------------------------
# T9 leg A' (plan §I.11; the active USER SETTING)
# ---------------------------------------------------------------------------


def test_leg_a_is_byte_inert_on_the_birth_surface(tmp_path, monkeypatch):
    """Under LEG='a' the cfgs carry the T7 byte-inert defaults — the four
    delivered T3 dirs remain reproducible from this generator."""
    monkeypatch.setattr(script, "LEG", "a")
    for _, cfg, _ in script.build_confirmation(tmp_path):
        assert cfg.birth_position_law == "boltzmann"
        assert cfg.initial_position_margin_angstrom == 0.0


def test_leg_aprime_flips_exactly_one_lever(tmp_path, monkeypatch):
    """Leg A' = leg A + the uniform_volume birth law at margin 3 A; every
    other stamped field is identical, and the run dirs carry the ap_ prefix
    (distinct from the delivered T3 dirs, same conf namespace)."""
    monkeypatch.setattr(script, "LEG", "aprime")
    aprime = script.build_confirmation(tmp_path)
    monkeypatch.setattr(script, "LEG", "a")
    a = script.build_confirmation(tmp_path)

    import dataclasses

    for (_, cfg_ap, dir_ap), (_, cfg_a, dir_a), spec in zip(
        aprime, a, script.CONFIRMATION_MATRIX
    ):
        assert cfg_ap.birth_position_law == "uniform_volume"
        assert cfg_ap.initial_position_margin_angstrom == pytest.approx(3.0)
        # The retained policy and the longer relaxation cap are consequences
        # the lever forces (well-trapped ions never decouple; quasi-bound
        # captures need time to resolve), not second physics levers.
        assert cfg_ap.detection_droplet_retained_policy == "exclude"
        assert cfg_ap.relaxation_time_ps == pytest.approx(8000.0)
        # A' stays on the delivered cold default (the A'' lever, not A's).
        assert cfg_ap.evaporation_shed_convention == "cold"
        cfg_ap.validate()
        # exactly one physics lever (+ its two bookkeeping/adequacy
        # consequences): every other field equal
        diff = {
            f.name
            for f in dataclasses.fields(type(cfg_ap))
            if getattr(cfg_ap, f.name) != getattr(cfg_a, f.name)
        }
        assert diff == {
            "birth_position_law",
            "initial_position_margin_angstrom",
            "detection_droplet_retained_policy",
            "relaxation_time_ps",
        }
        # distinct dirs, conf namespace, no collision with the T3 names
        assert dir_ap.name != dir_a.name
        assert f"_tier2probe_conf270_ap{spec.label}" in dir_ap.name
        assert not fnmatch(dir_ap.name, "*_tier2_*")


def test_leg_aprime_cm_flips_exactly_one_lever_vs_aprime(tmp_path, monkeypatch):
    """Leg A'' = leg A' + the co-moving shed convention (the OQ-J working-
    convention adjudication, 2026-07-17) -- exactly one field differs vs A',
    and the run dirs carry the apcm_ prefix (distinct from both the T3 and
    the A' dirs, same conf namespace)."""
    # LEG is a mutable USER SETTING (it rotates per oracle-chain leg), so the
    # test pins it via monkeypatch like every sibling -- asserting the ambient
    # value would turn routine leg rotation into a test failure.
    monkeypatch.setattr(script, "LEG", "aprime_cm")
    acm = script.build_confirmation(tmp_path)
    monkeypatch.setattr(script, "LEG", "aprime")
    aprime = script.build_confirmation(tmp_path)

    import dataclasses

    for (_, cfg_cm, dir_cm), (_, cfg_ap, dir_ap), spec in zip(
        acm, aprime, script.CONFIRMATION_MATRIX
    ):
        assert cfg_cm.evaporation_shed_convention == "co_moving"
        cfg_cm.validate()
        diff = {
            f.name
            for f in dataclasses.fields(type(cfg_cm))
            if getattr(cfg_cm, f.name) != getattr(cfg_ap, f.name)
        }
        assert diff == {"evaporation_shed_convention"}
        assert dir_cm.name != dir_ap.name
        assert f"_tier2probe_conf270_apcm{spec.label}" in dir_cm.name
        assert not fnmatch(dir_cm.name, "*_tier2_*")


def test_unknown_leg_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(script, "LEG", "b")
    with pytest.raises(ValueError, match="LEG"):
        script.build_confirmation(tmp_path)


def test_f3_discovery_ignores_complete_conf_dir(tmp_path):
    """A complete conf run dir must not enter the F3 campaign scoreboard."""
    conf_dir = tmp_path / tier2_confirmation_run_dir_name(
        "9A", "shared_pure_cubic", 50, config_label="c1", budget_eV=2.70
    )
    conf_dir.mkdir()
    for name in ("cfg.json", "ion.npz", "relaxation.npz", "detection.npz"):
        (conf_dir / name).write_text("x", encoding="utf-8")

    assert f3.discover_run_dirs(tmp_path) == []


# ---------------------------------------------------------------------------
# End-to-end (tiny N -- all five artifacts through all four stages)
# ---------------------------------------------------------------------------


def _tiny_capped_tabulated_cfg():
    """A tiny production-kinematics config through the full T3 surface:
    capped_cubic drag + rq4graded table + gated cooling + relaxation +
    detection, small windows so pytest stays cheap."""
    cfg = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        num_molecules=2,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=123,
        f_int=0.25 / 2.70,
        tau_ps=3.8,
        evap_rrk_dof=8.0,
        coulomb_available_eV=PRODUCTION_BUDGET_EV,
        # 100 ps cap: long enough that every Coulomb-pushed fragment is far
        # outside the ~28 A droplet at handover (the detection stage's P1-P3
        # decoupling guard is strict, EPS_DRAIN = 1e-6), short enough for
        # pytest. The production pilots use the 1000 ps §I.10 cap.
        relaxation_time_ps=100.0,
        cooling_spatial_gate="density_scaled",
        drag_form="capped_cubic",
        drag_coefficient_overrides={"v_c": 7.5, "p_tail": -1.0},
        dissociation_ladder="tabulated",
        tabulated_ladder_rungs_eV=script.rq4graded_rungs_eV(),
        R0_GS_angstrom=2.666,
        E_coulomb_scale=1.0,
        single_initial_position=False,
        detection_time_ps=200.0,  # beyond the 0.02 + 100 ps nominal window
    )
    return replace(cfg, t_max_neutral=0.1, dt_neutral=0.01)


def test_run_one_writes_all_five_artifacts(tmp_path):
    run_dir = tmp_path / "run"
    script._run_one("smoke c1", _tiny_capped_tabulated_cfg(), run_dir)
    for artifact in script._REQUIRED_ARTIFACTS:
        assert (run_dir / artifact).exists(), f"missing {artifact}"
    assert "detection.npz" in script._REQUIRED_ARTIFACTS


def test_run_one_refuses_partial_and_skips_complete(tmp_path, monkeypatch):
    monkeypatch.setattr(script, "OVERWRITE_EXISTING_RUN", False)
    monkeypatch.setattr(script, "SKIP_COMPLETED_RUNS", True)
    cfg = _tiny_capped_tabulated_cfg()

    partial = tmp_path / "partial"
    partial.mkdir()
    (partial / "cfg.json").write_text("{}", encoding="utf-8")
    with pytest.raises(FileExistsError):
        script._run_one("partial", cfg, partial)

    complete = tmp_path / "complete"
    complete.mkdir()
    for name in script._REQUIRED_ARTIFACTS:
        (complete / name).write_text("x", encoding="utf-8")
    before = {name: (complete / name).read_text() for name in script._REQUIRED_ARTIFACTS}
    script._run_one("complete", cfg, complete)  # must skip, not clobber
    after = {name: (complete / name).read_text() for name in script._REQUIRED_ARTIFACTS}
    assert after == before


def test_main_schedules_four(tmp_path, monkeypatch):
    scheduled: list = []
    monkeypatch.setattr(script, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(
        script, "_run_one", lambda label, cfg, run_dir: scheduled.append(run_dir)
    )
    assert script.main() == 0
    assert len(scheduled) == 4
    assert scheduled == [run_dir for _, _, run_dir in script.build_confirmation(tmp_path)]
