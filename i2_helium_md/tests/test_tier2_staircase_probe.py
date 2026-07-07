"""Coverage for the Tier-2 staircase capability probe generator (pre-F5).

Harness tests per ``TIER2_STAIRCASE_PROBE_PLAN.md`` §3/§5: the probe tag lives
in its own run-dir namespace (never swept into the F3 campaign scoreboard), the
45-point kappa x picture x tau grid builds with the bridge-pinned timing knobs
(f_int = 0.5 -- the timing-faithful choice, not the Stage-1 floor), the
(kappa=1, mixture, tau=6.55) grid point reproduces the Phase-D bridge config
(the wiring oracle), and a tiny-N end-to-end run writes all four artifacts.
No production scale, no figures.
"""

from __future__ import annotations

from dataclasses import replace
from fnmatch import fnmatch

import pytest

from scripts import gen_tier2_staircase_probe as script
from scripts.post_processing import tier2_size_distribution_table as f3
from scripts.tier2_common import (
    VALIDATION_BUDGET_EV,
    build_biphasic_cfg,
    tier2_probe_run_dir_name,
    tier2_probe_run_tag,
    tier2_run_dir_name,
    tier2_run_tag,
)


# ---------------------------------------------------------------------------
# Probe tag / namespace
# ---------------------------------------------------------------------------


def test_probe_tag_encodes_knobs_with_two_decimals():
    tag = tier2_probe_run_tag(
        picture="statistical_mixture",
        kappa=2.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=6.55,
        budget_eV=0.80,
    )
    assert tag.startswith("tier2probe_b080_mix")
    for token in ("_k2.00", "_l0.90", "_fi0.50", "_fr0.10", "_tau6.55"):
        assert token in tag, f"missing {token} in {tag}"


def test_probe_tag_rejects_unknown_picture():
    with pytest.raises(ValueError, match="picture"):
        tier2_probe_run_tag(
            picture="not_a_picture",
            kappa=1.0,
            lambda0_per_ps=0.9,
            f_int=0.5,
            f_ret=0.1,
            tau_ps=6.55,
            budget_eV=0.80,
        )


def test_probe_dir_namespace_disjoint_from_campaign():
    """Namespace lock (plan §3): the F3 campaign glob must never match a probe
    dir, and the probe glob must never match a campaign dir."""
    probe_name = tier2_probe_run_dir_name(
        "9A",
        "shared_pure_cubic",
        50,
        picture="x2_only",
        kappa=4.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=2.6,
        budget_eV=0.80,
    )
    campaign_name = tier2_run_dir_name(
        "9A",
        "shared_pure_cubic",
        500,
        picture="x2_only",
        kappa=4.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=2.6,
        budget_eV=0.80,
    )
    assert "_tier2probe_" in probe_name
    assert not fnmatch(probe_name, "*_tier2_*")  # the F3 campaign glob
    assert not fnmatch(campaign_name, "*_tier2probe_*")  # the probe glob


def test_f3_discovery_ignores_complete_probe_dir(tmp_path):
    """A *complete* probe run dir must not enter the campaign scoreboard."""
    probe_dir = tmp_path / tier2_probe_run_dir_name(
        "9A",
        "shared_pure_cubic",
        50,
        picture="statistical_mixture",
        kappa=1.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=6.55,
        budget_eV=0.80,
    )
    probe_dir.mkdir()
    for name in ("cfg.json", "ion.npz", "relaxation.npz"):
        (probe_dir / name).write_text("x", encoding="utf-8")

    assert f3.discover_run_dirs(tmp_path) == []


# ---------------------------------------------------------------------------
# Grid + config builds
# ---------------------------------------------------------------------------


def test_probe_grid_default_is_total_strip_ab():
    """The committed USER SETTINGS are the total-strip capability A/B probe: kappa
    and picture bridge-pinned, tau x s_eff swept, cooling_gate A/B'd ->
    1 x 1 x 3 x 4 x 2 = 24 unique points."""
    points = script.probe_grid_points()
    assert len(points) == (
        len(script.PICTURE_LIST)
        * len(script.KAPPA_GRID)
        * len(script.TAU_GRID_PS)
        * len(script.S_EFF_GRID)
        * len(script.COOLING_GATE_GRID)
    )
    assert len(points) == 24
    assert len(set(points)) == 24  # unique (picture, kappa, tau, s_eff, gate) tuples
    # the low-s_eff sweep (s=1 = max-kinetics bound) and the cooling A/B axis
    assert {s_eff for _, _, _, s_eff, _ in points} == {1, 2, 3, 5}
    assert {gate for _, _, _, _, gate in points} == {"none", "density_scaled"}


def test_full_grid_back_compat_none_s_eff(tmp_path, monkeypatch):
    """Addendum A.3 promise: the delivered 45-run probe path is byte-reproducible
    -- the full kappa x picture x tau grid at S_EFF_GRID=[None] and the ungated
    cooling arm yields 45 points, and a probe tag with s_eff/gate unset is
    byte-identical to the pre-addendum tag."""
    monkeypatch.setattr(script, "KAPPA_GRID", [0.5, 1.0, 2.0, 4.0, 8.0])
    monkeypatch.setattr(
        script, "PICTURE_LIST", ["statistical_mixture", "x2_only", "cooling_relaxed"]
    )
    monkeypatch.setattr(script, "TAU_GRID_PS", [2.6, 6.55, 16.5])
    monkeypatch.setattr(script, "S_EFF_GRID", [None])
    monkeypatch.setattr(script, "COOLING_GATE_GRID", ["none"])

    points = script.probe_grid_points()
    assert len(points) == 45
    assert all(s_eff is None for _, _, _, s_eff, _ in points)
    assert all(gate == "none" for _, _, _, _, gate in points)

    # s_eff=None and gate="none" append no suffix -> identical to the delivered tag.
    assert tier2_probe_run_tag(
        picture="statistical_mixture",
        kappa=1.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=6.55,
        budget_eV=0.80,
    ) == "tier2probe_b080_mix_k1.00_l0.90_fi0.50_fr0.10_tau6.55"


def test_build_probe_pins_bridge_timing_knobs(tmp_path):
    """Every probe cfg carries the bridge-pinned f_int/f_ret/lambda0 and the
    0.80 eV validation budget; relaxation (E2) is enabled at the experimental
    cap; the s_eff override is stamped (incl. the None per-n arm); dirs unique."""
    scheduled = script.build_probe(tmp_path)
    assert len(scheduled) == 24
    for _, cfg, _ in scheduled:
        assert cfg.mass_scenario == "biphasic"
        assert cfg.internal_energy_partition_fraction == pytest.approx(0.5)
        assert cfg.internal_energy_retained_fraction == pytest.approx(0.1)
        assert cfg.pickup_rate_coefficient == pytest.approx(0.9)
        assert cfg.coulomb_available_eV == pytest.approx(VALIDATION_BUDGET_EV)
        assert cfg.relaxation_stage_enabled is True
        assert cfg.relaxation_time_ps == pytest.approx(1000.0)  # finite total-strip cap
        assert cfg.ladder_steepness == pytest.approx(1.0)
        assert cfg.ladder_electronic_picture == "statistical_mixture"
    # the low-s_eff sweep and the cooling A/B are stamped onto the configs
    assert {cfg.evap_rrk_dof for _, cfg, _ in scheduled} == {1, 2, 3, 5}
    assert {cfg.cooling_spatial_gate for _, cfg, _ in scheduled} == {
        "none", "density_scaled"
    }
    run_dirs = [run_dir for _, _, run_dir in scheduled]
    assert len(set(run_dirs)) == 24
    assert all("_tier2probe_" in run_dir.name for run_dir in run_dirs)


def test_bridge_wiring_oracle_point_matches_bridge_config(tmp_path, monkeypatch):
    """The (kappa=1, statistical_mixture, tau=6.55, s_eff=None, gate=none) grid point
    IS the Phase-D bridge configuration (plan §2): its generative knobs must equal
    the pure bridge-default build, so its ion stage reproduces the bridge numbers.
    The active total-strip A/B grid drops the s_eff=None / tau=6.55 bridge point, so
    the grid is monkeypatched back to it for this wiring oracle (the pipeline check
    is grid-independent; only the active USER SETTINGS moved)."""
    monkeypatch.setattr(script, "KAPPA_GRID", [1.0])
    monkeypatch.setattr(script, "PICTURE_LIST", ["statistical_mixture"])
    monkeypatch.setattr(script, "TAU_GRID_PS", [6.55])
    monkeypatch.setattr(script, "S_EFF_GRID", [None])
    monkeypatch.setattr(script, "COOLING_GATE_GRID", ["none"])

    bridge_cfg = build_biphasic_cfg(
        script.CASE,
        script.VARIANT,
        num_molecules=script.N,
        ion_time_ps=script.ION_TIME_PS,
        dt_ion_ps=script.DT_ION_PS,
        seed=script.SEED,
    )
    matches = [
        cfg
        for _, cfg, _ in script.build_probe(tmp_path)
        if cfg.ladder_steepness == pytest.approx(1.0)
        and cfg.ladder_electronic_picture == "statistical_mixture"
        and cfg.internal_energy_cooling_tau_ps == pytest.approx(6.55)
        and cfg.evap_rrk_dof is None
        and cfg.cooling_spatial_gate == "none"
    ]
    assert len(matches) == 1
    probe_cfg = matches[0]
    for field in (
        "num_molecules",
        "coulomb_available_eV",
        "pickup_rate_coefficient",
        "internal_energy_partition_fraction",
        "internal_energy_retained_fraction",
        "ladder_electronic_picture",
        "ladder_steepness",
        "internal_energy_cooling_tau_ps",
        "evap_rrk_dof",
        "cooling_spatial_gate",
    ):
        assert getattr(probe_cfg, field) == getattr(bridge_cfg, field), field


def test_build_probe_rejects_non_validation_budget(tmp_path, monkeypatch):
    """The staircase anchor exists only at 9 A / 0.80 eV -- any other budget
    (including the sanctioned production 2.70) is refused by the generator."""
    for bad_budget in (2.70, 8.0):
        monkeypatch.setattr(script, "BUDGET_EV", bad_budget)
        with pytest.raises(ValueError, match="0.80"):
            script.build_probe(tmp_path)


# ---------------------------------------------------------------------------
# _run_one guards + tiny-N end-to-end
# ---------------------------------------------------------------------------


def test_run_one_refuses_partial_and_skips_complete(tmp_path, monkeypatch):
    monkeypatch.setattr(script, "OVERWRITE_EXISTING_RUN", False)
    monkeypatch.setattr(script, "SKIP_COMPLETED_RUNS", True)

    partial = tmp_path / "partial"
    partial.mkdir()
    (partial / "cfg.json").write_text("{}", encoding="utf-8")
    with pytest.raises(FileExistsError, match="OVERWRITE_EXISTING_RUN"):
        script._run_one("partial", object(), partial)

    complete = tmp_path / "complete"
    complete.mkdir()
    for name in script._REQUIRED_ARTIFACTS:
        (complete / name).write_text("x", encoding="utf-8")
    # object() as cfg would explode on any propagation call -> proves early return.
    script._run_one("complete", object(), complete)


def test_main_schedules_expected_grid(tmp_path, monkeypatch):
    scheduled = []
    monkeypatch.setattr(script, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(script, "_run_one", lambda label, cfg, run_dir: scheduled.append(run_dir))

    assert script.main() == 0

    expected_dirs = [
        tmp_path / "data" / "runs" / tier2_probe_run_dir_name(
            script.CASE,
            script.VARIANT,
            script.N,
            picture=picture,
            kappa=kappa,
            lambda0_per_ps=script.LAMBDA0_PER_PS,
            f_int=script.F_INT,
            f_ret=script.F_RET,
            tau_ps=tau_ps,
            budget_eV=script.BUDGET_EV,
            evap_rrk_dof=s_eff,
            cooling_spatial_gate=cooling_gate,
        )
        for picture, kappa, tau_ps, s_eff, cooling_gate in script.probe_grid_points()
    ]
    assert scheduled == expected_dirs


def test_run_one_writes_valid_run_dir(tmp_path):
    """Tiny-N end-to-end: neutral -> ion -> relaxation writes the four artifacts."""
    cfg = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        num_molecules=2,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=123,
        picture="statistical_mixture",
        kappa=1.0,
        tau_ps=6.55,
        f_int=0.5,
        f_ret=0.1,
        coulomb_available_eV=0.80,
        relaxation_time_ps=1.0,  # small cap: bounds the relaxation loop for the test
    )
    cfg = replace(cfg, t_max_neutral=0.1, dt_neutral=0.01)

    run_dir = tmp_path / "run"
    script._run_one("smoke", cfg, run_dir)

    for artifact in ("cfg.json", "neutral.npz", "ion.npz", "relaxation.npz"):
        assert (run_dir / artifact).exists(), f"missing {artifact}"


# ---------------------------------------------------------------------------
# s_eff (RRK-dof) override -- Addendum A
# ---------------------------------------------------------------------------


def test_probe_tag_appends_and_omits_s_eff_suffix():
    """A set s_eff appends _sNN.NN (two decimals, staying in the probe namespace);
    the None (per-n) default appends nothing (byte-identity with the delivered tag)."""
    base_kwargs = dict(
        picture="statistical_mixture",
        kappa=1.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=6.55,
        budget_eV=0.80,
    )
    assert tier2_probe_run_tag(**base_kwargs) == tier2_probe_run_tag(
        **base_kwargs, evap_rrk_dof=None
    )
    tag_s = tier2_probe_run_tag(**base_kwargs, evap_rrk_dof=12.0)
    assert tag_s == tier2_probe_run_tag(**base_kwargs) + "_s12.00"
    assert tag_s.startswith("tier2probe_")
    assert not fnmatch(tag_s, "*_tier2_*")  # still disjoint from the F3 campaign glob


def test_build_biphasic_cfg_stamps_evap_rrk_dof():
    """evap_rrk_dof is a None-sentinel pass-through: set -> stamped and validated;
    None -> rides the config default (None = per-n)."""
    cfg_set = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        num_molecules=2,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=1,
        picture="statistical_mixture",
        kappa=1.0,
        tau_ps=6.55,
        evap_rrk_dof=8.0,
    )
    assert cfg_set.evap_rrk_dof == pytest.approx(8.0)

    cfg_default = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        num_molecules=2,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=1,
    )
    assert cfg_default.evap_rrk_dof is None  # config default, not injected


def test_s_eff_points_yield_distinct_run_dirs():
    """Two grid points differing only in s_eff must not collide (F3 attributes by
    dir name; a collision would abort a run or mis-score it)."""
    common = dict(
        picture="statistical_mixture",
        kappa=1.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=6.55,
        budget_eV=0.80,
    )
    names = {
        tier2_probe_run_dir_name("9A", "shared_pure_cubic", 50, **common, evap_rrk_dof=s)
        for s in (None, 30, 20, 12, 8, 5)
    }
    assert len(names) == 6


# ---------------------------------------------------------------------------
# cooling_spatial_gate (total-strip A/B) override -- probe-scoped
# ---------------------------------------------------------------------------


def test_probe_tag_appends_and_omits_cooling_gate_suffix():
    """density_scaled appends _cgds (staying in the probe namespace); none / None
    append nothing (byte-identity with the pre-arm probe tag). The suffix follows
    any s_eff suffix, so an A/B pair at fixed s_eff differs only by _cgds."""
    base = dict(
        picture="statistical_mixture",
        kappa=1.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=6.55,
        budget_eV=0.80,
    )
    # gate default / "none" -> no suffix
    assert tier2_probe_run_tag(**base) == tier2_probe_run_tag(
        **base, cooling_spatial_gate=None
    )
    assert tier2_probe_run_tag(**base) == tier2_probe_run_tag(
        **base, cooling_spatial_gate="none"
    )
    # density_scaled -> _cgds, appended after the s_eff suffix
    tag_ds = tier2_probe_run_tag(**base, evap_rrk_dof=8.0, cooling_spatial_gate="density_scaled")
    assert tag_ds == tier2_probe_run_tag(**base, evap_rrk_dof=8.0) + "_cgds"
    assert tag_ds.startswith("tier2probe_")
    assert not fnmatch(tag_ds, "*_tier2_*")  # still disjoint from the F3 campaign glob


def test_cooling_gate_points_yield_distinct_run_dirs():
    """The A/B arms (none vs density_scaled) at a fixed knob point must map to
    distinct run dirs so both are generated and scored separately."""
    common = dict(
        picture="statistical_mixture",
        kappa=1.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=6.55,
        budget_eV=0.80,
        evap_rrk_dof=5.0,
    )
    none_dir = tier2_probe_run_dir_name(
        "9A", "shared_pure_cubic", 50, **common, cooling_spatial_gate="none"
    )
    ds_dir = tier2_probe_run_dir_name(
        "9A", "shared_pure_cubic", 50, **common, cooling_spatial_gate="density_scaled"
    )
    assert none_dir != ds_dir
    assert ds_dir.endswith("_cgds")


def test_build_biphasic_cfg_stamps_cooling_spatial_gate():
    """cooling_spatial_gate is a None-sentinel pass-through: set -> stamped and
    validated; None -> rides the config default ('none', ungated cooling)."""
    cfg_set = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        num_molecules=2,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=1,
        cooling_spatial_gate="density_scaled",
    )
    assert cfg_set.cooling_spatial_gate == "density_scaled"

    cfg_default = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        num_molecules=2,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=1,
    )
    assert cfg_default.cooling_spatial_gate == "none"  # config default, not injected


def test_probe_tag_knob_body_matches_campaign_encoder():
    """The probe encoder promises to 'mirror' the campaign encoder's knob body;
    lock that so a future change to one (added knob / changed precision) can't
    silently diverge the two formats. (Review Finding 4: the two remain parallel
    functions -- not refactored to share, to avoid touching the frozen campaign
    encoder -- so this parity test is the divergence guard.)"""
    knobs = dict(
        picture="statistical_mixture",
        kappa=2.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=6.55,
        budget_eV=0.80,
    )
    campaign = tier2_run_tag(**knobs)              # tier2_b080_mix_k2.00_..._tau6.55
    probe = tier2_probe_run_tag(**knobs)           # tier2probe_b080_mix_k2.00_..._tau6.55
    # same knob body; the only intended difference is the namespace prefix.
    assert probe == campaign.replace("tier2_", "tier2probe_", 1)
    # ...including the s_eff suffix (campaign encoder gained it at the 2026-07-06
    # promotion; the parity lock now guards that dimension too).
    knobs_s = dict(knobs, evap_rrk_dof=8.0)
    assert tier2_probe_run_tag(**knobs_s) == tier2_run_tag(**knobs_s).replace(
        "tier2_", "tier2probe_", 1
    )
    # ...and the cooling-gate suffix (the campaign encoder gained _cgds so a gated
    # campaign cannot collide run-dir names; the parity lock now guards it too).
    knobs_cg = dict(knobs_s, cooling_spatial_gate="density_scaled")
    assert tier2_probe_run_tag(**knobs_cg) == tier2_run_tag(**knobs_cg).replace(
        "tier2_", "tier2probe_", 1
    )


def test_campaign_tag_appends_and_omits_cooling_gate_suffix():
    """tier2_run_tag mirrors the probe encoder: density_scaled -> _cgds (after the
    s_eff suffix, before any _totalstrip); none / None append nothing (campaigns not
    sweeping the gate stay byte-identical)."""
    base = dict(
        picture="statistical_mixture",
        kappa=1.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=6.55,
        budget_eV=0.80,
    )
    assert tier2_run_tag(**base) == tier2_run_tag(**base, cooling_spatial_gate=None)
    assert tier2_run_tag(**base) == tier2_run_tag(**base, cooling_spatial_gate="none")
    # density_scaled -> _cgds
    assert (
        tier2_run_tag(**base, cooling_spatial_gate="density_scaled")
        == tier2_run_tag(**base) + "_cgds"
    )
    # after the s_eff suffix
    assert (
        tier2_run_tag(**base, evap_rrk_dof=8.0, cooling_spatial_gate="density_scaled")
        == tier2_run_tag(**base, evap_rrk_dof=8.0) + "_cgds"
    )
    # before the _totalstrip suffix
    assert (
        tier2_run_tag(
            **base, cooling_spatial_gate="density_scaled", total_strip=True
        )
        == tier2_run_tag(**base) + "_cgds_totalstrip"
    )
