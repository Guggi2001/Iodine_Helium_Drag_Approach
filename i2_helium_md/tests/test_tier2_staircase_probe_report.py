"""Coverage for the Tier-2 staircase probe report (pre-F5 scorer).

Per ``TIER2_STAIRCASE_PROBE_PLAN.md`` §4/§5: the staircase metrics are
hand-oracled on constructed ``n_shell`` arrays (an ion that walks the anchored
staircase exactly must score a perfect match), one genuine tiny-N probe run is
scored into a row with the expected columns, discovery honours the probe
namespace, and the CSV round-trips. No figures in pytest.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from scripts import gen_tier2_staircase_probe as gen_script
from scripts.post_processing import tier2_staircase_probe_report as report
from scripts.tier2_common import build_biphasic_cfg, tier2_probe_run_dir_name
from i2_helium_md.physics.shell_schedule import (
    ANCHOR_N_END,
    ANCHOR_N_START,
    NUM_SHED_EVENTS,
    build_shell_schedule,
)


# ---------------------------------------------------------------------------
# staircase_metrics: hand oracles on constructed n_shell
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def schedule():
    return build_shell_schedule(5.0)


@pytest.fixture(scope="module")
def time_ps():
    # 0..20 ps at 0.05 ps: every anchored shed time (multiples of 0.05) lies
    # exactly on the grid, so first-shed detection has no discretization slop.
    return np.arange(0.0, 20.0 + 1e-12, 0.05)


def test_metrics_exact_on_the_anchored_staircase(schedule, time_ps):
    """One ion that walks the anchored 21->19->14 staircase exactly scores a
    perfect match: 7 sheds, terminal 14, zero trajectory MAD, first shed at
    the schedule's first event time."""
    n_shell = schedule.n_of_t(time_ps).astype(float)[None, :]  # (1, Nt)

    m = report.staircase_metrics(time_ps, n_shell, schedule)

    assert m["dn_mean_window"] == pytest.approx(float(NUM_SHED_EVENTS))
    assert m["n_ion_end_mean"] == pytest.approx(float(ANCHOR_N_END))
    assert m["n_traj_mad"] == pytest.approx(0.0)
    assert m["t_first_shed_ps"] == pytest.approx(schedule.events[0].time_ps)
    assert m["frac_ions_shed"] == pytest.approx(1.0)


def test_metrics_mixed_ensemble(schedule, time_ps):
    """Two ions -- one frozen at 21, one anchored -- give the ensemble means:
    dn = 3.5, terminal mean 17.5, half the ions shed, and the first-shed time
    averages only over the ions that shed (nanmean semantics)."""
    anchored = schedule.n_of_t(time_ps).astype(float)
    frozen = np.full_like(anchored, float(ANCHOR_N_START))
    n_shell = np.stack([frozen, anchored])

    m = report.staircase_metrics(time_ps, n_shell, schedule)

    assert m["dn_mean_window"] == pytest.approx(NUM_SHED_EVENTS / 2.0)
    assert m["n_ion_end_mean"] == pytest.approx((ANCHOR_N_START + ANCHOR_N_END) / 2.0)
    assert m["n_traj_mad"] > 0.0
    assert m["t_first_shed_ps"] == pytest.approx(schedule.events[0].time_ps)
    assert m["frac_ions_shed"] == pytest.approx(0.5)


def test_metrics_no_shed_ensemble(schedule, time_ps):
    """An ensemble that never sheds: zero sheds, NaN first-shed time (never a
    crash), zero shed fraction."""
    n_shell = np.full((3, time_ps.size), float(ANCHOR_N_START))

    m = report.staircase_metrics(time_ps, n_shell, schedule)

    assert m["dn_mean_window"] == pytest.approx(0.0)
    assert np.isnan(m["t_first_shed_ps"])
    assert m["frac_ions_shed"] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Genuine tiny-N probe run -> one scored row
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def tiny_probe_run(tmp_path_factory):
    """One genuine tiny-N probe run dir (module-scoped: built once)."""
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
        relaxation_time_ps=1.0,
    )
    cfg = replace(cfg, t_max_neutral=0.1, dt_neutral=0.01)
    run_dir = tmp_path_factory.mktemp("probe_runs") / tier2_probe_run_dir_name(
        "9A",
        "shared_pure_cubic",
        2,
        picture="statistical_mixture",
        kappa=1.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=6.55,
        budget_eV=0.80,
    )
    gen_script._run_one("tiny-probe", cfg, run_dir)
    return run_dir


def test_score_probe_run_row(tiny_probe_run, schedule):
    row = report.score_probe_run(tiny_probe_run, schedule)

    assert set(row) == set(report.PROBE_TABLE_COLUMNS)
    # knob columns echo the authoritative cfg.json
    assert row["case"] == "9A"
    assert row["budget_eV"] == pytest.approx(0.80)
    assert row["picture"] == "statistical_mixture"
    assert row["kappa"] == pytest.approx(1.0)
    assert row["tau_ps"] == pytest.approx(6.55)
    assert row["s_eff"] is None  # per-n convention (evap_rrk_dof unset)
    assert row["cooling_gate"] == "none"  # config default (ungated cooling)
    assert row["f_int"] == pytest.approx(0.5)
    assert row["f_ret"] == pytest.approx(0.1)
    assert row["N"] == 2
    # staircase metrics: finite (t_first may legitimately be NaN on a tiny run)
    assert np.isfinite(row["dn_mean_window"])
    assert np.isfinite(row["n_ion_end_mean"])
    assert np.isfinite(row["n_traj_mad"])
    assert 0.0 <= row["frac_ions_shed"] <= 1.0
    # flexibility / total-strip metrics from the relaxed checkpoint
    assert 0.0 <= row["n_relaxed_mean"] <= float(ANCHOR_N_START)
    assert row["n_relaxed_spread"] >= 0.0
    assert 0 <= row["n_relaxed_min"] <= int(ANCHOR_N_START)
    assert row["n_relaxed_min"] <= row["n_relaxed_mean"] + 1e-9  # min <= mean
    assert 0.0 <= row["frac_frozen"] <= 1.0
    # wiring diagnostic
    assert np.isfinite(row["ledger_max_resid_eV"])


# ---------------------------------------------------------------------------
# Discovery, table, CSV, main
# ---------------------------------------------------------------------------


def _fake_complete_probe_dir(root, **kwargs):
    d = root / tier2_probe_run_dir_name("9A", "shared_pure_cubic", 50, **kwargs)
    d.mkdir(parents=True)
    for name in ("cfg.json", "ion.npz", "relaxation.npz"):
        (d / name).write_text("x", encoding="utf-8")
    return d


_TAG_KWARGS = dict(
    picture="statistical_mixture",
    kappa=1.0,
    lambda0_per_ps=0.9,
    f_int=0.5,
    f_ret=0.1,
    tau_ps=6.55,
    budget_eV=0.80,
)


def test_discovery_honours_probe_namespace(tmp_path):
    complete = _fake_complete_probe_dir(tmp_path, **_TAG_KWARGS)

    incomplete = tmp_path / complete.name.replace("tau6.55", "tau2.60")
    incomplete.mkdir()
    (incomplete / "cfg.json").write_text("x", encoding="utf-8")

    campaign = tmp_path / "9A_drag_shared_pure_cubic_N500_tier2_b080_mix_k1.00"
    campaign.mkdir()
    for name in ("cfg.json", "ion.npz", "relaxation.npz"):
        (campaign / name).write_text("x", encoding="utf-8")

    assert report.discover_probe_run_dirs(tmp_path) == [complete]


def test_scored_row_survives_csv_round_trip(tiny_probe_run, schedule, tmp_path):
    row = report.score_probe_run(tiny_probe_run, schedule)
    path = report.write_rows_csv(tmp_path / "probe.csv", [row])

    import csv

    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 1
    assert rows[0]["picture"] == "statistical_mixture"
    assert float(rows[0]["kappa"]) == pytest.approx(1.0)


def test_score_probe_run_reads_s_eff_override(tmp_path_factory, schedule):
    """A run built with a set evap_rrk_dof reports it in the s_eff column
    (read from the authoritative cfg.json, not the tag)."""
    cfg = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        num_molecules=2,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=7,
        picture="statistical_mixture",
        kappa=1.0,
        tau_ps=6.55,
        f_int=0.5,
        f_ret=0.1,
        coulomb_available_eV=0.80,
        relaxation_time_ps=1.0,
        evap_rrk_dof=8.0,
    )
    cfg = replace(cfg, t_max_neutral=0.1, dt_neutral=0.01)
    run_dir = tmp_path_factory.mktemp("probe_seff") / tier2_probe_run_dir_name(
        "9A",
        "shared_pure_cubic",
        2,
        picture="statistical_mixture",
        kappa=1.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=6.55,
        budget_eV=0.80,
        evap_rrk_dof=8.0,
    )
    gen_script._run_one("tiny-seff", cfg, run_dir)

    row = report.score_probe_run(run_dir, schedule)
    assert row["s_eff"] == pytest.approx(8.0)
    assert "s_eff" in report.PROBE_TABLE_COLUMNS


def test_score_probe_run_reads_cooling_gate(tmp_path_factory, schedule):
    """A run built with the density_scaled cooling arm reports it in the
    cooling_gate column (read from the authoritative cfg.json, not the tag)."""
    cfg = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        num_molecules=2,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=7,
        picture="statistical_mixture",
        kappa=1.0,
        tau_ps=6.55,
        f_int=0.5,
        f_ret=0.1,
        coulomb_available_eV=0.80,
        relaxation_time_ps=1.0,
        evap_rrk_dof=5.0,
        cooling_spatial_gate="density_scaled",
    )
    cfg = replace(cfg, t_max_neutral=0.1, dt_neutral=0.01)
    run_dir = tmp_path_factory.mktemp("probe_cg") / tier2_probe_run_dir_name(
        "9A",
        "shared_pure_cubic",
        2,
        picture="statistical_mixture",
        kappa=1.0,
        lambda0_per_ps=0.9,
        f_int=0.5,
        f_ret=0.1,
        tau_ps=6.55,
        budget_eV=0.80,
        evap_rrk_dof=5.0,
        cooling_spatial_gate="density_scaled",
    )
    gen_script._run_one("tiny-cg", cfg, run_dir)

    row = report.score_probe_run(run_dir, schedule)
    assert row["cooling_gate"] == "density_scaled"
    assert "cooling_gate" in report.PROBE_TABLE_COLUMNS
    assert "n_relaxed_min" in row and "frac_frozen" in row


def test_format_table_empty():
    assert "no staircase-probe rows" in report.format_table([])


def test_format_headline_reports_range_and_argmin():
    # The mini-probe sweeps s_eff at a pinned (kappa, picture, tau); the argmin
    # line MUST name s_eff or the operator cannot read the winning RRK-dof.
    rows = [
        {
            "picture": "statistical_mixture",
            "kappa": 1.0,
            "tau_ps": 6.55,
            "s_eff": None,
            "n_relaxed_mean": 20.3,
            "n_traj_mad": 5.1,
        },
        {
            "picture": "statistical_mixture",
            "kappa": 1.0,
            "tau_ps": 6.55,
            "s_eff": 8.0,
            "n_relaxed_mean": 13.9,
            "n_traj_mad": 0.4,
        },
        {
            "picture": "statistical_mixture",
            "kappa": 1.0,
            "tau_ps": 16.5,
            "s_eff": 20.0,
            "n_relaxed_mean": 17.0,
            "n_traj_mad": 2.0,
        },
    ]
    text = report.format_headline(rows)
    # per-picture reachable relaxed terminal-n range
    assert "statistical_mixture" in text
    assert "13.9" in text and "20.3" in text
    # factual argmin of the trajectory MAD (reported, not adjudicated) -- s_eff
    # is the swept lever, so it must appear and disambiguate the winner.
    assert "tau=6.55" in text and "s_eff=8" in text
    assert "not auto-adjudicated" in text


def test_format_headline_argmin_names_per_n_s_eff():
    """When the per-n (s_eff=None) control wins, the argmin renders it 'per-n'."""
    rows = [
        {
            "picture": "statistical_mixture",
            "kappa": 1.0,
            "tau_ps": 6.55,
            "s_eff": None,
            "n_relaxed_mean": 14.0,
            "n_traj_mad": 0.3,
        },
        {
            "picture": "statistical_mixture",
            "kappa": 1.0,
            "tau_ps": 6.55,
            "s_eff": 8.0,
            "n_relaxed_mean": 20.0,
            "n_traj_mad": 4.0,
        },
    ]
    assert "s_eff=per-n" in report.format_headline(rows)


# ---------------------------------------------------------------------------
# best-case selection (drives the single-panel overlay figure; the render
# itself is gated out of pytest, but the pure selection is testable)
# ---------------------------------------------------------------------------


def test_select_best_case_row_picks_lowest_mad():
    """The best case (the one the overlay figure draws) is the lowest-trajectory-
    MAD grid point -- the same 'best' the text headline names, so figure and
    headline agree by construction. Here that point lands n_end ~ 14."""
    rows = [
        {"picture": "statistical_mixture", "kappa": 1.0, "tau_ps": 6.55,
         "s_eff": None, "n_ion_end_mean": 20.3, "n_traj_mad": 5.1},
        {"picture": "statistical_mixture", "kappa": 1.0, "tau_ps": 6.55,
         "s_eff": 8.0, "n_ion_end_mean": 13.6, "n_traj_mad": 1.0},
        {"picture": "statistical_mixture", "kappa": 1.0, "tau_ps": 16.5,
         "s_eff": 20.0, "n_ion_end_mean": 15.0, "n_traj_mad": 1.9},
    ]
    best = report.select_best_case_row(rows)
    assert best is not None
    assert best["s_eff"] == 8.0
    assert best["n_ion_end_mean"] == pytest.approx(13.6)


def test_select_best_case_row_matches_headline_argmin():
    """select_best_case_row is the single source of 'best' -- the headline's
    argmin line must name the same point it returns."""
    rows = [
        {"picture": "statistical_mixture", "kappa": 1.0, "tau_ps": 6.55,
         "s_eff": None, "n_relaxed_mean": 20.3, "n_ion_end_mean": 20.3,
         "n_traj_mad": 5.1},
        {"picture": "statistical_mixture", "kappa": 1.0, "tau_ps": 6.55,
         "s_eff": 8.0, "n_relaxed_mean": 13.9, "n_ion_end_mean": 13.6,
         "n_traj_mad": 0.4},
    ]
    best = report.select_best_case_row(rows)
    assert best["s_eff"] == 8.0
    assert f"s_eff={report._s_eff_label(best['s_eff'])}" in report.format_headline(rows)


def test_select_best_case_row_empty_and_nonfinite():
    """No rows -> None; all-NaN MAD (no ion trajectory scored) -> None, never a
    crash on min() over an empty finite set."""
    assert report.select_best_case_row([]) is None
    rows = [
        {"picture": "x2_only", "kappa": 1.0, "tau_ps": 6.55, "s_eff": 8.0,
         "n_ion_end_mean": 14.0, "n_traj_mad": float("nan")},
    ]
    assert report.select_best_case_row(rows) is None


def _matcher_cfg(gate=None):
    return build_biphasic_cfg(
        "9A", "shared_pure_cubic", num_molecules=2, ion_time_ps=0.02,
        dt_ion_ps=0.01, seed=1, cooling_spatial_gate=gate,
    )


def _matcher_row(cfg, **extra):
    row = {
        "picture": cfg.ladder_electronic_picture,
        "kappa": float(cfg.ladder_steepness),
        "tau_ps": float(cfg.internal_energy_cooling_tau_ps),
        "budget_eV": float(cfg.coulomb_available_eV),
        "s_eff": None if cfg.evap_rrk_dof is None else float(cfg.evap_rrk_dof),
    }
    row.update(extra)
    return row


def test_cfg_matches_row_distinguishes_cooling_gate():
    """The A/B arms share every other knob, so without cooling_gate in the match the
    best-case overlay would re-load the wrong arm's ion checkpoint (Finding 2)."""
    cfg_none = _matcher_cfg("none")
    cfg_ds = _matcher_cfg("density_scaled")
    row_none = _matcher_row(cfg_none, cooling_gate="none")
    assert report._cfg_matches_row(cfg_none, row_none)
    assert not report._cfg_matches_row(cfg_ds, row_none)


def test_cfg_matches_row_back_compat_missing_cooling_gate():
    """A pre-arm row (no cooling_gate key) still matches a default-none cfg."""
    cfg = _matcher_cfg()  # rides the config default 'none'
    assert report._cfg_matches_row(cfg, _matcher_row(cfg))  # no cooling_gate key


def test_main_returns_zero_on_empty_root(tmp_path, monkeypatch):
    monkeypatch.setattr(report, "RUNS_ROOT", tmp_path)
    assert report.main() == 0
