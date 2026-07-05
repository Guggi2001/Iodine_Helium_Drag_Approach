"""Smoke coverage for the Tier-2 campaign run-matrix generator (Slice F2).

Two harness tests (plan §2/§5): the grid enumeration schedules the expected 9 A
run dirs (no propagation, via a captured ``_run_one``), and a tiny-N single grid
point runs neutral -> ion -> relaxation end-to-end and writes a valid run dir.
No production scale, no figures.
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from scripts import gen_tier2_runs as script
from scripts.tier2_common import build_biphasic_cfg, tier2_run_dir_name
from i2_helium_md.physics.internal_energy_budget import f_int_floor


def test_campaign_grid_applies_floor_f_int_when_unset():
    """F_INT=None -> every grid point pins f_int at its picture/kappa floor."""
    assert script.F_INT is None  # the delivered Stage-1 default
    points = script.campaign_grid_points()
    assert len(points) == len(script.PICTURE_LIST) * len(script.KAPPA_GRID)
    for picture, kappa, f_int in points:
        assert f_int == pytest.approx(
            f_int_floor(e_avail_eV=script.BUDGET_EV, picture=picture, kappa=kappa)
        )


def test_resolve_f_int_literal_override(monkeypatch):
    """A literal F_INT overrides the floor (conscious timing comparison)."""
    monkeypatch.setattr(script, "F_INT", 0.42)
    assert script.resolve_f_int("statistical_mixture", 3.0) == pytest.approx(0.42)
    # ...and every grid point then carries that literal, not the floor.
    assert all(f_int == pytest.approx(0.42) for _, _, f_int in script.campaign_grid_points())


def test_build_campaign_uses_experimental_relaxation_cap(tmp_path):
    """RELAXATION_TIME_PS=None -> the 8530 ns experimental cap is stamped on
    every grid point, and the E2 stage is enabled."""
    assert script.RELAXATION_TIME_PS is None  # delivered default
    scheduled = script.build_campaign(tmp_path)
    for _, cfg, _ in scheduled:
        assert cfg.relaxation_stage_enabled is True
        assert cfg.relaxation_time_ps == pytest.approx(8.53e6)
        assert cfg.relaxation_forces == script.RELAXATION_FORCES


def test_build_campaign_honours_relaxation_cap_override(tmp_path, monkeypatch):
    monkeypatch.setattr(script, "RELAXATION_TIME_PS", 250.0)
    scheduled = script.build_campaign(tmp_path)
    assert all(cfg.relaxation_time_ps == pytest.approx(250.0) for _, cfg, _ in scheduled)


def test_main_schedules_expected_grid(tmp_path, monkeypatch):
    scheduled = []

    def capture_run(label, cfg, run_dir):
        scheduled.append((label, cfg, run_dir))

    monkeypatch.setattr(script, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(script, "_run_one", capture_run)

    assert script.main() == 0

    expected_dirs = [
        tmp_path / "data" / "runs" / tier2_run_dir_name(
            script.CASE,
            script.VARIANT,
            script.N,
            picture=picture,
            kappa=kappa,
            lambda0_per_ps=script.LAMBDA0_PER_PS,
            f_int=f_int,
            f_ret=script.F_RET,
            tau_ps=script.TAU_PS,
            budget_eV=script.BUDGET_EV,
            total_strip=script.TOTAL_STRIP,
        )
        for picture, kappa, f_int in script.campaign_grid_points()
    ]

    assert [run_dir for _, _, run_dir in scheduled] == expected_dirs
    assert all(cfg.mass_scenario == "biphasic" for _, cfg, _ in scheduled)
    assert all(cfg.relaxation_stage_enabled for _, cfg, _ in scheduled)
    assert all(cfg.coulomb_available_eV == script.BUDGET_EV for _, cfg, _ in scheduled)
    # distinct run dir per grid point
    assert len({run_dir for _, _, run_dir in scheduled}) == len(expected_dirs)


def test_build_campaign_rejects_unsanctioned_budget(tmp_path, monkeypatch):
    """Review finding #3: a mistyped BUDGET_EV is caught before it poisons runs."""
    monkeypatch.setattr(script, "BUDGET_EV", 8.0)  # typo for 0.80
    with pytest.raises(ValueError, match="not a sanctioned budget"):
        script.build_campaign(tmp_path)


def test_run_one_refuses_partial_existing_output(tmp_path, monkeypatch):
    """A partial (incomplete) run dir trips the overwrite guard, never silently kept."""
    monkeypatch.setattr(script, "OVERWRITE_EXISTING_RUN", False)
    monkeypatch.setattr(script, "SKIP_COMPLETED_RUNS", True)
    run_dir = tmp_path / "partial"
    run_dir.mkdir()
    (run_dir / "cfg.json").write_text("{}", encoding="utf-8")  # only 1 of 4 artifacts

    with pytest.raises(FileExistsError, match="OVERWRITE_EXISTING_RUN"):
        script._run_one("partial", object(), run_dir)


def test_run_one_skips_completed_run(tmp_path, monkeypatch):
    """Review finding #5: a fully-written run dir is skipped on resume (no re-run)."""
    monkeypatch.setattr(script, "SKIP_COMPLETED_RUNS", True)
    run_dir = tmp_path / "complete"
    run_dir.mkdir()
    for name in script._REQUIRED_ARTIFACTS:
        (run_dir / name).write_text("x", encoding="utf-8")

    # object() as cfg would explode on any propagation call -> proves we returned early.
    script._run_one("complete", object(), run_dir)


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
        relaxation_time_ps=1.0,   # small cap: bounds the relaxation loop for the test
    )
    # Shrink the neutral stage so the smoke test stays cheap (fixture precedent).
    cfg = replace(cfg, t_max_neutral=0.1, dt_neutral=0.01)

    run_dir = tmp_path / "run"
    script._run_one("smoke", cfg, run_dir)

    for artifact in ("cfg.json", "neutral.npz", "ion.npz", "relaxation.npz"):
        assert (run_dir / artifact).exists(), f"missing {artifact}"
