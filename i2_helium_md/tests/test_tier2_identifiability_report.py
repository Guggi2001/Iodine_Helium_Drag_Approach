"""Harness coverage for the Tier-2 identifiability + regime report (Slice F4).

F4 turns the F3 scoreboard rows into a **per-quantity identifiability statement**
(the campaign's central deliverable), not a point fit. It is a pure assembler --
it reads scoreboard rows and emits the CALIBRATION_MAP "Tier-2 load-bearing"
checklist: which of the 8+ quantities the single size-distribution observable can
actually separate.

These tests drive the assessment logic on **constructed synthetic scoreboard
rows** (the plan's stated F4 test surface) -- no run pipeline, no figures. The
9 Aa-only scope decision (2026-07-03) makes ``f_ret`` and ``lambda_attach``
documented non-deliverables; the tests assert F4 *reports* them as
not-identifiable rather than attempting the (unavailable) 9/18 Aa contrast.
"""

from __future__ import annotations

import csv

import pytest

import scripts.post_processing.tier2_identifiability_report as f4mod
from scripts.post_processing.tier2_identifiability_report import (
    IDENTIFIABILITY_COLUMNS,
    LOAD_BEARING_QUANTITIES,
    STATUS_HELD_FIXED,
    STATUS_IDENTIFIABLE,
    STATUS_INSENSITIVE,
    STATUS_INSUFFICIENT,
    STATUS_NOT_BRACKETED,
    STATUS_NOT_IDENTIFIABLE,
    STATUS_SENSITIVE,
    assess_kappa_landscape,
    assess_picture_separation,
    assess_regime,
    assess_tau_sensitivity,
    build_identifiability_rows,
    build_kappa_landscapes,
    build_tau_sweeps,
    format_report,
    write_report_csv,
)


def _row(**over):
    """One synthetic F3 scoreboard row with campaign-sane defaults."""
    base = dict(
        case="9A",
        budget_eV=0.80,
        picture="statistical_mixture",
        kappa=1.0,
        f_int=0.3,
        f_ret=0.1,
        tau_ps=6.55,
        N=500,
        W1_matched=1.0,
        W1_simend_upper=1.0,
        t_cross_ps=3.0,
        regime_label="shell_retaining",
        total_strip_reachable=False,
        ledger_max_resid_eV=1e-6,
        n_terminal_mean=15.0,
        n_terminal_spread=1.0,
        t_cross_spread_ps=0.1,
        t_cross_all_agree=True,
        sanity_flags="",
    )
    base.update(over)
    return base


# ---------------------------------------------------------------------------
# kappa landscape assessment
# ---------------------------------------------------------------------------
def test_kappa_landscape_identifiable_interior_min():
    # clear interior minimum at kappa=2
    points = [(0.5, 3.0), (1.0, 2.0), (2.0, 1.0), (4.0, 2.5), (8.0, 3.5)]
    status, argmin_kappa, w1_min, _ = assess_kappa_landscape(points)
    assert status == STATUS_IDENTIFIABLE
    assert argmin_kappa == pytest.approx(2.0)
    assert w1_min == pytest.approx(1.0)


def test_kappa_landscape_boundary_min_not_bracketed():
    # monotone decreasing toward the cliff end -> argmin at the grid edge
    points = [(0.5, 5.0), (1.0, 4.0), (2.0, 3.0), (4.0, 2.0), (8.0, 1.0)]
    status, argmin_kappa, _, notes = assess_kappa_landscape(points)
    assert status == STATUS_NOT_BRACKETED
    assert argmin_kappa == pytest.approx(8.0)
    assert "edge" in notes.lower() or "bracket" in notes.lower()


def test_kappa_landscape_flat_not_identifiable():
    # W1 essentially constant across kappa -> degenerate landscape
    points = [(0.5, 1.00), (1.0, 1.01), (2.0, 1.00), (4.0, 1.02), (8.0, 1.01)]
    status, _, _, _ = assess_kappa_landscape(points)
    assert status == STATUS_NOT_IDENTIFIABLE


def test_kappa_landscape_single_point_insufficient():
    status, _, _, _ = assess_kappa_landscape([(1.0, 2.0)])
    assert status == STATUS_INSUFFICIENT


# ---------------------------------------------------------------------------
# picture separation
# ---------------------------------------------------------------------------
def test_picture_separation_identifiable():
    best_w1_by_picture = {"statistical_mixture": 1.0, "x2_only": 3.0}
    status, best_picture, _ = assess_picture_separation(best_w1_by_picture)
    assert status == STATUS_IDENTIFIABLE
    assert best_picture == "statistical_mixture"


def test_picture_separation_flat_not_identifiable():
    best_w1_by_picture = {"statistical_mixture": 1.00, "x2_only": 1.01}
    status, _, _ = assess_picture_separation(best_w1_by_picture)
    assert status == STATUS_NOT_IDENTIFIABLE


def test_picture_separation_single_picture_insufficient():
    status, _, _ = assess_picture_separation({"statistical_mixture": 1.0})
    assert status == STATUS_INSUFFICIENT


# ---------------------------------------------------------------------------
# tau sensitivity (terminal-n insensitivity confirm/flag)
# ---------------------------------------------------------------------------
def test_tau_sensitivity_sensitive():
    # terminal n moves by ~3 He across the tau band -> live dimension
    points = [(2.6, 13.0), (6.55, 15.0), (16.5, 16.0)]
    status, spread, _ = assess_tau_sensitivity(points)
    assert status == STATUS_SENSITIVE
    assert spread == pytest.approx(3.0)


def test_tau_sensitivity_insensitive():
    # terminal n flat across tau -> insensitivity confirmed
    points = [(2.6, 15.0), (6.55, 15.1), (16.5, 15.05)]
    status, spread, _ = assess_tau_sensitivity(points)
    assert status == STATUS_INSENSITIVE
    assert spread < f4mod.TAU_SENSITIVITY_THRESHOLD_HE


def test_tau_sensitivity_single_point_insufficient():
    status, _, _ = assess_tau_sensitivity([(6.55, 15.0)])
    assert status == STATUS_INSUFFICIENT


# ---------------------------------------------------------------------------
# grouping helpers
# ---------------------------------------------------------------------------
def test_build_kappa_landscapes_groups_by_picture():
    rows = [
        _row(picture="statistical_mixture", kappa=0.5, W1_matched=3.0),
        _row(picture="statistical_mixture", kappa=2.0, W1_matched=1.0),
        _row(picture="x2_only", kappa=0.5, W1_matched=4.0),
        _row(picture="x2_only", kappa=2.0, W1_matched=2.0),
    ]
    landscapes = build_kappa_landscapes(rows)
    assert len(landscapes) == 2
    pictures = {ls.picture for ls in landscapes}
    assert pictures == {"statistical_mixture", "x2_only"}


def test_build_kappa_landscapes_skips_single_kappa_cell():
    rows = [_row(kappa=1.0)]  # only one kappa -> not a landscape
    assert build_kappa_landscapes(rows) == []


def test_build_tau_sweeps_needs_two_taus():
    rows = [
        _row(tau_ps=2.6, n_terminal_mean=13.0),
        _row(tau_ps=16.5, n_terminal_mean=16.0),
    ]
    sweeps = build_tau_sweeps(rows)
    assert len(sweeps) == 1
    assert len(sweeps[0].points) == 2


def test_build_tau_sweeps_single_tau_is_no_sweep():
    assert build_tau_sweeps([_row(tau_ps=6.55)]) == []


# ---------------------------------------------------------------------------
# regime determination
# ---------------------------------------------------------------------------
def test_assess_regime_reports_best_w1_regime():
    rows = [
        _row(regime_label="shell_retaining", W1_matched=1.0),
        _row(regime_label="total_strip", W1_matched=3.0),
        _row(regime_label="total_strip", W1_matched=4.0),
    ]
    counts, best_regime, _ = assess_regime(rows)
    assert counts["shell_retaining"] == 1
    assert counts["total_strip"] == 2
    assert best_regime == "shell_retaining"  # lowest W1_matched


# ---------------------------------------------------------------------------
# the checklist
# ---------------------------------------------------------------------------
def _cofit_grid():
    """A kappa x picture co-fit grid with a clean interior-min best cell.

    statistical_mixture has an interior kappa min at 2.0 (W1=1.0) and is the
    globally best picture; x2_only is uniformly worse and well separated.
    """
    rows = []
    for k, w in [(0.5, 3.0), (1.0, 2.0), (2.0, 1.0), (4.0, 2.5), (8.0, 3.5)]:
        rows.append(_row(picture="statistical_mixture", kappa=k, W1_matched=w))
    for k, w in [(0.5, 6.0), (1.0, 5.5), (2.0, 5.0), (4.0, 5.5), (8.0, 6.0)]:
        rows.append(_row(picture="x2_only", kappa=k, W1_matched=w))
    return rows


def test_build_identifiability_rows_covers_all_load_bearing_once():
    rows = build_identifiability_rows(_cofit_grid())
    quantities = [r["quantity"] for r in rows]
    for q in LOAD_BEARING_QUANTITIES:
        assert quantities.count(q) == 1
    assert set(r["budget_eV"] for r in rows) == {0.80}


def test_build_identifiability_rows_kappa_and_picture_identifiable():
    rows = build_identifiability_rows(_cofit_grid())
    by_q = {r["quantity"]: r for r in rows}
    assert by_q["ladder_shape_kappa"]["status"] == STATUS_IDENTIFIABLE
    assert by_q["electronic_picture"]["status"] == STATUS_IDENTIFIABLE
    # the best cell picks kappa*=2.0 under the best picture
    assert "2" in str(by_q["ladder_shape_kappa"]["value_or_range"])
    assert "statistical_mixture" in str(by_q["electronic_picture"]["value_or_range"])


def test_build_identifiability_rows_static_nondeliverables():
    """9 Aa-only scope: f_ret and lambda_attach are reported not-identifiable with
    their distinct rationale; f_int is not separated by construction; nu/s held."""
    rows = build_identifiability_rows(_cofit_grid())
    by_q = {r["quantity"]: r for r in rows}

    assert by_q["f_ret"]["status"] == STATUS_NOT_IDENTIFIABLE
    assert "9" in by_q["f_ret"]["notes"] and "18" in by_q["f_ret"]["notes"]

    assert by_q["lambda_attach"]["status"] == STATUS_NOT_IDENTIFIABLE
    assert "pickup" in by_q["lambda_attach"]["notes"].lower()

    assert by_q["f_int"]["status"] == STATUS_NOT_IDENTIFIABLE
    assert "timing" in by_q["f_int"]["notes"].lower()

    assert by_q["nu"]["status"] == STATUS_HELD_FIXED
    assert by_q["s"]["status"] == STATUS_HELD_FIXED
    # s carries the s<->kappa coupling caveat (prose, not computed)
    assert "kappa" in by_q["s"]["notes"].lower() or "κ" in by_q["s"]["notes"]


def test_build_identifiability_rows_tau_sensitive_when_swept():
    grid = _cofit_grid()
    # add a tau sweep at the best cell (statistical_mixture, kappa=2.0)
    grid += [
        _row(picture="statistical_mixture", kappa=2.0, tau_ps=2.6,
             W1_matched=1.0, n_terminal_mean=13.0),
        _row(picture="statistical_mixture", kappa=2.0, tau_ps=16.5,
             W1_matched=1.0, n_terminal_mean=16.0),
    ]
    rows = build_identifiability_rows(grid)
    by_q = {r["quantity"]: r for r in rows}
    assert by_q["tau"]["status"] == STATUS_SENSITIVE


def test_build_identifiability_rows_tau_insufficient_single_tau():
    rows = build_identifiability_rows(_cofit_grid())  # all tau=6.55
    by_q = {r["quantity"]: r for r in rows}
    assert by_q["tau"]["status"] == STATUS_INSUFFICIENT


# ---------------------------------------------------------------------------
# per-budget sectioning
# ---------------------------------------------------------------------------
def test_build_identifiability_rows_sections_by_budget():
    rows = _cofit_grid() + [
        _row(budget_eV=2.70, picture="statistical_mixture", kappa=k, W1_matched=w)
        for k, w in [(0.5, 3.0), (2.0, 1.0)]
    ]
    checklist = build_identifiability_rows(rows)
    budgets = set(r["budget_eV"] for r in checklist)
    assert budgets == {0.80, 2.70}
    # each budget carries the full load-bearing checklist
    for b in budgets:
        qs = [r["quantity"] for r in checklist if r["budget_eV"] == b]
        assert set(qs) == set(LOAD_BEARING_QUANTITIES)


def test_format_report_has_headline_and_sections():
    text = format_report(_cofit_grid())
    assert "0.8" in text  # budget section header
    assert "identifiab" in text.lower()
    # headline reports an "K of M" separation count
    assert "of" in text.lower()


def test_format_report_empty_rows():
    assert isinstance(format_report([]), str)


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------
def test_write_report_csv_round_trips(tmp_path):
    checklist = build_identifiability_rows(_cofit_grid())
    path = write_report_csv(tmp_path / "identifiability.csv", checklist)
    with path.open(newline="", encoding="utf-8") as fh:
        read = list(csv.DictReader(fh))
    assert list(read[0].keys()) == IDENTIFIABILITY_COLUMNS
    assert len(read) == len(checklist)


def test_main_runs_on_empty_root(tmp_path, monkeypatch):
    monkeypatch.setattr(f4mod, "RUNS_ROOT", tmp_path / "empty")
    monkeypatch.setattr(f4mod, "BUDGET_EV", None)
    assert f4mod.main() == 0


# ===========================================================================
# Review-hardening (workflow high-effort review, 2026-07-05)
# ===========================================================================

# --- Finding 2: assess_kappa_landscape empty-points guard -------------------
def test_kappa_landscape_empty_is_insufficient_not_crash():
    """An empty landscape degrades to insufficient_data, never raises on min()."""
    status, argmin_kappa, w1_min, _ = assess_kappa_landscape([])
    assert status == STATUS_INSUFFICIENT
    assert argmin_kappa != argmin_kappa  # NaN sentinel
    assert w1_min != w1_min


# --- Finding 7: low-edge bracketing + all-zero / zero-min flatness ----------
def test_kappa_landscape_low_edge_not_bracketed():
    # monotone increasing -> argmin at index 0 (the low grid edge)
    points = [(0.5, 1.0), (1.0, 2.0), (2.0, 3.0), (4.0, 4.0), (8.0, 5.0)]
    status, argmin_kappa, _, _ = assess_kappa_landscape(points)
    assert status == STATUS_NOT_BRACKETED
    assert argmin_kappa == pytest.approx(0.5)


def test_kappa_landscape_all_zero_w1_is_flat():
    status, _, _, _ = assess_kappa_landscape([(0.5, 0.0), (1.0, 0.0), (2.0, 0.0)])
    assert status == STATUS_NOT_IDENTIFIABLE


def test_kappa_landscape_zero_min_interior_is_identifiable():
    # a perfect (W1=0) interior match against non-zero neighbours is a real min
    status, argmin_kappa, _, _ = assess_kappa_landscape(
        [(0.5, 1.0), (1.0, 0.0), (2.0, 1.0)]
    )
    assert status == STATUS_IDENTIFIABLE
    assert argmin_kappa == pytest.approx(1.0)


# --- Finding 8: assess_picture_separation zero-best-W1 branches --------------
def test_picture_separation_zero_best_is_identifiable():
    # one picture scores an exact 0 (perfect) against a >0 picture -> separated
    status, best, _ = assess_picture_separation(
        {"statistical_mixture": 0.0, "x2_only": 2.0}
    )
    assert status == STATUS_IDENTIFIABLE
    assert best == "statistical_mixture"


def test_picture_separation_all_zero_is_not_identifiable():
    status, _, _ = assess_picture_separation(
        {"statistical_mixture": 0.0, "x2_only": 0.0}
    )
    assert status == STATUS_NOT_IDENTIFIABLE


# --- Finding 6: NaN guarding on the W1 / terminal-n consumers ----------------
def test_assess_regime_skips_nan_w1_for_best_regime():
    rows = [
        _row(regime_label="shell_retaining", W1_matched=float("nan")),
        _row(regime_label="total_strip", W1_matched=1.0),
    ]
    counts, best_regime, _ = assess_regime(rows)
    assert counts == {"shell_retaining": 1, "total_strip": 1}
    assert best_regime == "total_strip"  # the NaN run is not the argmin


def test_build_kappa_landscapes_skips_nan_w1():
    rows = [
        _row(kappa=1.0, W1_matched=float("nan")),  # skipped
        _row(kappa=1.0, W1_matched=2.0),
        _row(kappa=2.0, W1_matched=1.0),
    ]
    landscapes = build_kappa_landscapes(rows)
    assert len(landscapes) == 1
    assert dict(landscapes[0].points) == {1.0: 2.0, 2.0: 1.0}


def test_build_tau_sweeps_skips_nan_terminal_n():
    rows = [
        _row(kappa=2.0, tau_ps=2.6, n_terminal_mean=float("nan")),  # skipped
        _row(kappa=2.0, tau_ps=2.6, n_terminal_mean=13.0),
        _row(kappa=2.0, tau_ps=16.5, n_terminal_mean=16.0),
    ]
    sweeps = build_tau_sweeps(rows)
    assert len(sweeps) == 1
    assert dict(sweeps[0].points) == {2.6: 13.0, 16.5: 16.0}


# --- Finding 5: duplicate knob-cell collision warns, never silently drops ----
def test_build_kappa_landscapes_warns_on_duplicate_cell():
    rows = [
        _row(kappa=1.0, W1_matched=2.0),
        _row(kappa=1.0, W1_matched=3.0),  # same cell+kappa, different W1
        _row(kappa=2.0, W1_matched=1.0),
    ]
    with pytest.warns(RuntimeWarning, match="duplicate kappa"):
        landscapes = build_kappa_landscapes(rows)
    assert dict(landscapes[0].points)[1.0] == 3.0  # last wins, but loudly


def test_build_kappa_landscapes_silent_on_identical_duplicate():
    rows = [
        _row(kappa=1.0, W1_matched=2.0),
        _row(kappa=1.0, W1_matched=2.0),  # identical -> harmless replay, no warn
        _row(kappa=2.0, W1_matched=1.0),
    ]
    import warnings as _w

    with _w.catch_warnings():
        _w.simplefilter("error")  # any warning would fail the test
        build_kappa_landscapes(rows)


# --- Finding 1 (+ contamination): picture co-fit excludes non-landscape rows -
def test_picture_best_is_not_contaminated_by_non_cofit_run():
    """A total_strip-style run (distinct f_int/tau -> single-kappa cell, not a
    landscape) with a very low W1 must NOT lower its picture's co-fit best or flip
    the reported best picture."""
    rows = _cofit_grid()  # mix best=1.0, x2_only best=5.0 (both kappa-swept)
    rows.append(
        _row(picture="x2_only", kappa=2.0, f_int=0.9, tau_ps=16.5, W1_matched=0.01)
    )
    checklist = build_identifiability_rows(rows)
    by_q = {r["quantity"]: r for r in checklist}
    # best picture stays statistical_mixture; x2_only's co-fit best stays 5, not 0.01
    assert by_q["electronic_picture"]["value_or_range"] == "statistical_mixture"
    assert "x2_only:W1=5" in by_q["electronic_picture"]["evidence"]
    # kappa verdict is attributed to the actual best picture
    assert "picture=statistical_mixture" in by_q["ladder_shape_kappa"]["evidence"]


def test_kappa_verdict_never_falls_back_to_a_non_best_picture():
    """best_picture is landscape-derived, so the kappa entry always assesses the
    best picture's own landscape (Finding 1)."""
    rows = _cofit_grid()
    checklist = build_identifiability_rows(rows)
    by_q = {r["quantity"]: r for r in checklist}
    best_pic = by_q["electronic_picture"]["value_or_range"]
    assert f"picture={best_pic}" in by_q["ladder_shape_kappa"]["evidence"]


# --- Finding 4: tau is read at the co-fit optimum kappa cell -----------------
def test_tau_read_at_best_kappa_not_insertion_order():
    """Two tau sweeps under the best picture: an insensitive one at kappa=1.0
    (inserted first) and a sensitive one at the co-fit optimum kappa=2.0. The tau
    verdict must come from kappa=2.0, not the first-inserted sweep."""
    grid = _cofit_grid()  # best cell = statistical_mixture, kappa*=2.0
    # off-optimum sweep at kappa=1.0 (flat terminal n), inserted first
    grid += [
        _row(picture="statistical_mixture", kappa=1.0, tau_ps=2.6,
             W1_matched=2.0, n_terminal_mean=15.0),
        _row(picture="statistical_mixture", kappa=1.0, tau_ps=16.5,
             W1_matched=2.0, n_terminal_mean=15.05),
    ]
    # optimum sweep at kappa=2.0 (terminal n moves 3 He)
    grid += [
        _row(picture="statistical_mixture", kappa=2.0, tau_ps=2.6,
             W1_matched=1.0, n_terminal_mean=13.0),
        _row(picture="statistical_mixture", kappa=2.0, tau_ps=16.5,
             W1_matched=1.0, n_terminal_mean=16.0),
    ]
    by_q = {r["quantity"]: r for r in build_identifiability_rows(grid)}
    assert by_q["tau"]["status"] == STATUS_SENSITIVE
    assert "kappa=2" in by_q["tau"]["evidence"]


# --- Finding 9: format_report accepts a pre-built checklist ------------------
def test_format_report_prebuilt_checklist_matches_from_rows():
    rows = _cofit_grid()
    checklist = build_identifiability_rows(rows)
    assert format_report(rows, checklist=checklist) == format_report(rows)


# --- collect wrapper wiring (F3 typed records -> checklist) ------------------
def test_collect_identifiability_rows_wraps_f3_records(monkeypatch):
    import types

    rows = _cofit_grid()
    fake = [types.SimpleNamespace(row=r) for r in rows]
    monkeypatch.setattr(
        f4mod, "collect_tier2_records", lambda **kw: (fake, object())
    )
    scoreboard, checklist = f4mod.collect_identifiability_rows()
    assert scoreboard == rows
    assert checklist == build_identifiability_rows(rows)
