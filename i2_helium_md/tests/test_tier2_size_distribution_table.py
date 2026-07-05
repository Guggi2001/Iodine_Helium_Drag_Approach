"""Harness coverage for the Tier-2 size-distribution scoreboard (Slice F3).

F3 is a *pure scorer*: it composes the delivered E1 (terminal size distribution,
matched-time via ``relaxation.npz`` + sim-end upper bound via ``ion.npz``), E4
(integer-support Wasserstein vs the experimental abundance reference), E5
(``reconstruct_diagnostics`` -> t*, regime label, total-strip reachability), and
the 5-term ``ion_ledger_closure`` residual into one scoreboard row per finished
campaign run dir. It runs nothing (F2 already wrote the artifacts).

These tests build **one** genuine tiny-N biphasic run dir (module-scoped, so the
neutral -> ion -> relaxation pipeline runs once) and reuse it by copying, plus a
few pure-unit checks (freeze-side keying, CSV round-trip, table columns). No
production scale, no figures.
"""

from __future__ import annotations

from dataclasses import replace
import json
import shutil

import numpy as np
import pytest

from i2_helium_md.physics.constants import N_STAR
from i2_helium_md.postprocess import ShellDistribution, load_he_abundance_reference
from i2_helium_md.simulation.ion import run_ion_propagation
from i2_helium_md.simulation.neutral import run_neutral_propagation
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage
from i2_helium_md.simulation.run_directory import RunDirectory
from scripts.tier2_common import (
    PRODUCTION_BUDGET_EV,
    VALIDATION_BUDGET_EV,
    build_biphasic_cfg,
)

import scripts.post_processing.tier2_size_distribution_table as f3mod
from scripts.post_processing.tier2_size_distribution_table import (
    TIER2_TABLE_COLUMNS,
    Tier2RunRecord,
    _case_from_run_dir,
    _distribution_moments,
    _freeze_side_expected,
    _is_total_strip_run,
    collect_tier2_records,
    collect_tier2_rows,
    discover_run_dirs,
    format_table,
    score_tier2_run,
    write_rows_csv,
)


PROJECT_ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]
REFERENCE_PATH = PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"

# The 16 pinned plan-F3 columns + the three confirmed advisory columns.
_PINNED_COLUMNS = [
    "case", "budget_eV", "picture", "kappa", "f_int", "f_ret", "tau_ps", "N",
    "W1_matched", "W1_simend_upper", "t_cross_ps", "regime_label",
    "total_strip_reachable", "ledger_max_resid_eV", "n_terminal_mean",
    "n_terminal_spread",
]
_ADVISORY_COLUMNS = ["t_cross_spread_ps", "t_cross_all_agree", "sanity_flags"]


def _write_tiny_run(run_dir, *, budget_eV=VALIDATION_BUDGET_EV):
    """Build one genuine tiny-N biphasic run dir (neutral -> ion -> relaxation)."""
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
        coulomb_available_eV=budget_eV,
        relaxation_time_ps=1.0,  # small cap: bounds the relaxation loop for the test
    )
    cfg = replace(cfg, t_max_neutral=0.1, dt_neutral=0.01)  # cheap neutral stage
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    run_relaxation_stage(ion, cfg, save_path=run.root / "relaxation.npz")
    return run_dir


@pytest.fixture(scope="module")
def canonical_run(tmp_path_factory):
    """One real tiny-N campaign run dir, built once for the whole module.

    Its basename mirrors the campaign convention (``9A_drag_...`` prefix) so the
    ``case`` column can be recovered from the dir name.
    """
    root = tmp_path_factory.mktemp("tier2_f3")
    run_dir = root / "9A_drag_shared_pure_cubic_N2_tier2_b080_mix_k1.00"
    _write_tiny_run(run_dir)
    return run_dir


@pytest.fixture(scope="module")
def reference():
    return load_he_abundance_reference(REFERENCE_PATH)


def _copy_run(canonical_run, dest):
    shutil.copytree(canonical_run, dest)
    return dest


# ---------------------------------------------------------------------------
# Pure-unit checks (no run pipeline)
# ---------------------------------------------------------------------------
def test_freeze_side_expected_true_at_validation_budget():
    """Decision #2: freeze-side expected at 0.80 eV (Pi>1 is a wiring smell)."""
    cfg = build_biphasic_cfg(
        "9A", "shared_pure_cubic", num_molecules=1, ion_time_ps=0.02,
        dt_ion_ps=0.01, seed=1, coulomb_available_eV=VALIDATION_BUDGET_EV,
    )
    assert _freeze_side_expected(cfg) is True


def test_freeze_side_expected_false_at_production_budget():
    """Decision #2: NOT freeze-side at 2.70 eV (Pi>1 is legitimate physics)."""
    cfg = build_biphasic_cfg(
        "9A", "shared_pure_cubic", num_molecules=1, ion_time_ps=0.02,
        dt_ion_ps=0.01, seed=1, coulomb_available_eV=PRODUCTION_BUDGET_EV,
    )
    assert _freeze_side_expected(cfg) is False


def test_table_columns_are_the_pinned_plus_advisory_set():
    assert TIER2_TABLE_COLUMNS == _PINNED_COLUMNS + _ADVISORY_COLUMNS


def test_write_rows_csv_round_trips(tmp_path):
    import csv

    rows = [{col: (1.5 if col not in ("case", "picture", "regime_label",
                                      "sanity_flags") else "x")
             for col in TIER2_TABLE_COLUMNS}]
    path = write_rows_csv(tmp_path / "board.csv", rows)
    with path.open(newline="", encoding="utf-8") as fh:
        read = list(csv.DictReader(fh))
    assert list(read[0].keys()) == TIER2_TABLE_COLUMNS
    assert len(read) == 1


def test_format_table_handles_empty_rows():
    assert isinstance(format_table([]), str)


# ---------------------------------------------------------------------------
# Scoring one real run
# ---------------------------------------------------------------------------
def test_score_tier2_run_returns_row_with_all_columns(canonical_run, reference):
    row = score_tier2_run(canonical_run, reference)
    assert set(row.keys()) == set(TIER2_TABLE_COLUMNS)


def test_score_tier2_run_values_are_sane(canonical_run, reference):
    row = score_tier2_run(canonical_run, reference)

    # both distances present, finite, non-negative (integer-support W1)
    assert np.isfinite(row["W1_matched"]) and row["W1_matched"] >= 0.0
    assert np.isfinite(row["W1_simend_upper"]) and row["W1_simend_upper"] >= 0.0

    # ledger residual is a finite wiring diagnostic
    assert np.isfinite(row["ledger_max_resid_eV"]) and row["ledger_max_resid_eV"] >= 0.0

    # regime label is one of the two documented outcomes
    assert row["regime_label"] in ("shell_retaining", "total_strip")
    assert isinstance(row["total_strip_reachable"], bool)

    # terminal-n moments come from the matched (relaxed) distribution
    assert 0.0 <= row["n_terminal_mean"] <= float(N_STAR)
    assert row["n_terminal_spread"] >= 0.0

    # knob columns echo the run's cfg
    assert row["case"] == "9A"
    assert row["N"] == 2
    assert row["budget_eV"] == pytest.approx(VALIDATION_BUDGET_EV)
    assert row["picture"] == "statistical_mixture"
    assert row["kappa"] == pytest.approx(1.0)
    assert row["f_ret"] == pytest.approx(0.1)


def test_score_tier2_run_accepts_run_directory_object(canonical_run, reference):
    row = score_tier2_run(RunDirectory(canonical_run), reference)
    assert set(row.keys()) == set(TIER2_TABLE_COLUMNS)


# ---------------------------------------------------------------------------
# Discovery + collection over a runs root
# ---------------------------------------------------------------------------
def test_discover_finds_complete_tier2_runs(canonical_run, tmp_path):
    runs_root = tmp_path / "runs"
    runs_root.mkdir()
    _copy_run(canonical_run, runs_root / canonical_run.name)
    _copy_run(canonical_run, runs_root / (canonical_run.name + "_b"))

    found = discover_run_dirs(runs_root)
    assert len(found) == 2
    assert all((d / "ion.npz").exists() for d in found)


def test_discover_skips_incomplete_run(canonical_run, tmp_path):
    runs_root = tmp_path / "runs"
    runs_root.mkdir()
    partial = _copy_run(canonical_run, runs_root / canonical_run.name)
    (partial / "relaxation.npz").unlink()  # drop a required artifact

    assert discover_run_dirs(runs_root) == []


def test_discover_ignores_non_tier2_dirs(canonical_run, tmp_path):
    runs_root = tmp_path / "runs"
    runs_root.mkdir()
    _copy_run(canonical_run, runs_root / canonical_run.name)
    # a Tier-0/1a-style dir without the tier2 marker must be ignored
    _copy_run(canonical_run, runs_root / "9A_drag_shared_pure_cubic_N50")

    found = discover_run_dirs(runs_root)
    assert len(found) == 1
    assert "tier2" in found[0].name


def test_collect_scores_every_discovered_run(canonical_run, tmp_path):
    runs_root = tmp_path / "runs"
    runs_root.mkdir()
    _copy_run(canonical_run, runs_root / canonical_run.name)
    _copy_run(canonical_run, runs_root / (canonical_run.name + "_b"))

    records, ref = collect_tier2_records(
        runs_root=runs_root, reference_path=REFERENCE_PATH
    )
    assert len(records) == 2
    assert all(isinstance(r, Tier2RunRecord) for r in records)
    assert all(set(r.row.keys()) == set(TIER2_TABLE_COLUMNS) for r in records)

    rows = collect_tier2_rows(runs_root=runs_root, reference_path=REFERENCE_PATH)
    assert len(rows) == 2


def test_collect_budget_filter_selects_matching_runs(canonical_run, tmp_path):
    runs_root = tmp_path / "runs"
    runs_root.mkdir()
    _copy_run(canonical_run, runs_root / canonical_run.name)

    kept = collect_tier2_rows(
        runs_root=runs_root, reference_path=REFERENCE_PATH,
        budget_eV=VALIDATION_BUDGET_EV,
    )
    assert len(kept) == 1

    dropped = collect_tier2_rows(
        runs_root=runs_root, reference_path=REFERENCE_PATH,
        budget_eV=PRODUCTION_BUDGET_EV,
    )
    assert dropped == []


# ---------------------------------------------------------------------------
# Review findings (2026-07-05)
# ---------------------------------------------------------------------------
def test_freeze_side_expected_false_for_total_strip_at_validation():
    """Finding 1: a 0.80 eV total_strip run sheds by construction, so Pi>1 is
    legitimate physics there -- it must NOT be flagged freeze-side."""
    cfg = build_biphasic_cfg(
        "9A", "shared_pure_cubic", num_molecules=1, ion_time_ps=0.02,
        dt_ion_ps=0.01, seed=1, coulomb_available_eV=VALIDATION_BUDGET_EV,
    )
    assert _freeze_side_expected(cfg, total_strip=True) is False
    # ...while the plain validation run stays freeze-side:
    assert _freeze_side_expected(cfg, total_strip=False) is True


def test_discover_excludes_reserved_bridge_run(canonical_run, tmp_path):
    """Finding 2: the reserved Phase-D bridge dir must never enter the campaign
    scoreboard, even if it happens to carry all three artifacts."""
    runs_root = tmp_path / "runs"
    runs_root.mkdir()
    _copy_run(canonical_run, runs_root / canonical_run.name)  # a real campaign run
    # a bridge dir with the reserved tag, complete with all artifacts:
    _copy_run(
        canonical_run,
        runs_root / "9A_drag_shared_pure_cubic_N50_tier2_bridge_biphasic",
    )

    found = discover_run_dirs(runs_root)
    assert len(found) == 1
    assert "bridge" not in found[0].name


# ---------------------------------------------------------------------------
# Coverage locks (regression locks on already-correct behavior)
# ---------------------------------------------------------------------------
def test_is_total_strip_run_reads_suffix():
    from pathlib import Path

    assert _is_total_strip_run(Path("9A_drag_x_N5_tier2_b080_mix_k1.00_totalstrip"))
    assert not _is_total_strip_run(Path("9A_drag_x_N5_tier2_b080_mix_k1.00"))


def test_case_from_run_dir_recovers_case_and_rejects_malformed():
    from pathlib import Path

    assert _case_from_run_dir(Path("18A_drag_shared_pure_cubic_N500_tier2_b270")) == "18A"
    with pytest.raises(ValueError, match="_drag_"):
        _case_from_run_dir(Path("not-a-campaign-dir"))


def test_distribution_moments_hand_oracle():
    # A two-spike distribution at n=0 and n=2 (each 1/2): mean 1, std 1.
    dist = ShellDistribution(
        n_values=np.array([0, 1, 2]),
        counts=np.array([1, 0, 1]),
        fraction=np.array([0.5, 0.0, 0.5]),
        source="relaxed",
    )
    mean, spread = _distribution_moments(dist)
    assert mean == pytest.approx(1.0)
    assert spread == pytest.approx(1.0)


def test_score_wires_relaxed_to_matched_and_simend_to_upper(
    canonical_run, reference, monkeypatch
):
    """The E-review's stated F3 hazard: matched W1 must score the *relaxed*
    distribution and W1_simend the *sim_end* one. Spy on the compare call and
    lock the source tag each W1 sees, in order (matched first, then sim-end)."""
    seen: list[str] = []

    def spy(sim, ref, *, metric):
        seen.append(sim.source)
        return float(len(seen))  # distinct finite values so the row still fills

    monkeypatch.setattr(f3mod, "compare_size_distributions", spy)
    row = score_tier2_run(canonical_run, reference)

    assert seen == ["relaxed", "sim_end"]
    assert row["W1_matched"] == pytest.approx(1.0)      # first compare (relaxed)
    assert row["W1_simend_upper"] == pytest.approx(2.0)  # second compare (sim_end)


def test_real_scored_row_csv_round_trips(canonical_run, reference, tmp_path):
    """A real row (t_cross may be NaN) survives write_rows_csv unbroken."""
    import csv

    row = score_tier2_run(canonical_run, reference)
    path = write_rows_csv(tmp_path / "board.csv", [row])
    with path.open(newline="", encoding="utf-8") as fh:
        read = list(csv.DictReader(fh))
    assert list(read[0].keys()) == TIER2_TABLE_COLUMNS
    assert len(read) == 1


def test_collect_budget_filter_selects_between_two_budgets(canonical_run, tmp_path):
    """The budget filter keys off each run's cfg stamp, not the dir name: two
    runs at different budgets are separated by budget_eV."""
    runs_root = tmp_path / "runs"
    runs_root.mkdir()
    _copy_run(canonical_run, runs_root / canonical_run.name)  # 0.80 eV (as built)
    prod = _copy_run(canonical_run, runs_root / (canonical_run.name + "_prod"))
    # Restamp the copy's budget to 2.70 in its cfg.json (the filter reads cfg).
    cfg_path = prod / "cfg.json"
    payload = json.loads(cfg_path.read_text(encoding="utf-8"))
    payload["coulomb_available_eV"] = PRODUCTION_BUDGET_EV
    cfg_path.write_text(json.dumps(payload), encoding="utf-8")

    val_rows = collect_tier2_rows(
        runs_root=runs_root, reference_path=REFERENCE_PATH,
        budget_eV=VALIDATION_BUDGET_EV,
    )
    prod_rows = collect_tier2_rows(
        runs_root=runs_root, reference_path=REFERENCE_PATH,
        budget_eV=PRODUCTION_BUDGET_EV,
    )
    assert len(val_rows) == 1 and val_rows[0]["budget_eV"] == pytest.approx(VALIDATION_BUDGET_EV)
    assert len(prod_rows) == 1 and prod_rows[0]["budget_eV"] == pytest.approx(PRODUCTION_BUDGET_EV)


def test_main_runs_on_empty_root(tmp_path, monkeypatch):
    """The operator entry returns 0 and reports zero runs on an empty root."""
    monkeypatch.setattr(f3mod, "RUNS_ROOT", tmp_path / "empty")
    monkeypatch.setattr(f3mod, "BUDGET_EV", None)
    assert f3mod.main() == 0
