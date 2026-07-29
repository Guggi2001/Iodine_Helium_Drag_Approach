"""Tests for the §3.5g low-n KE retro-scan helpers.

Covers `scripts/post_processing/tier2atlas_ke_lown_scan.py`:

- ``fit_slope`` — least-squares slope; ValueError on degenerate input;
- ``max_ingate_move`` — headroom arithmetic of the R5 reachability verdict
  (binding observable wins; zero slopes do not constrain; out-of-band
  current values raise);
- ``lowke_counts`` — per-bin populations;
- ``ke_oracle_mismatches`` — the O2 committed-CSV comparison (pass, numeric
  mismatch, missing label).

Synthetic inputs only — no run dirs, no figures (testing rules). Tolerances:
exact arithmetic identities, so ``pytest.approx`` defaults.
"""

from __future__ import annotations

import csv
from types import SimpleNamespace

import numpy as np
import pytest

from scripts.post_processing.tier2atlas_ke_lown_scan import (
    fit_slope,
    ke_oracle_mismatches,
    lowke_counts,
    max_ingate_move,
)


class TestFitSlope:
    def test_exact_line(self):
        assert fit_slope([(0.0, 1.0), (2.0, 5.0)]) == pytest.approx(2.0)

    def test_least_squares_over_three_points(self):
        # y = 3x with one symmetric perturbation: slope stays 3.
        assert fit_slope([(0, 0.1), (1, 3.0), (2, 5.9)]) == pytest.approx(2.9)

    def test_single_point_raises(self):
        with pytest.raises(ValueError, match=">= 2 points"):
            fit_slope([(1.0, 1.0)])

    def test_identical_x_raises(self):
        with pytest.raises(ValueError, match="all x identical"):
            fit_slope([(1.0, 1.0), (1.0, 2.0)])


class TestMaxIngateMove:
    BANDS = {"a": (0.0, 1.0), "b": (-1.0, 1.0)}

    def test_binding_observable_limits_move(self):
        # a: headroom up 0.2 at slope 0.1 -> 2.0; b: headroom up 0.5 at
        # slope 1.0 -> 0.5 binds.
        move = max_ingate_move({"a": 0.8, "b": 0.5},
                               {"a": 0.1, "b": 1.0}, self.BANDS, +1.0)
        assert move == pytest.approx(0.5)

    def test_negative_direction_uses_lower_headroom(self):
        move = max_ingate_move({"a": 0.8, "b": 0.5},
                               {"a": 0.1, "b": 1.0}, self.BANDS, -1.0)
        # a: down-headroom 0.8 / 0.1 = 8; b: 1.5 / 1 = 1.5 binds.
        assert move == pytest.approx(1.5)

    def test_negative_slope_consumes_lower_band(self):
        # slope < 0 with direction +1 walks the observable DOWN.
        move = max_ingate_move({"a": 0.2, "b": 0.0},
                               {"a": -1.0, "b": 0.0}, self.BANDS, +1.0)
        assert move == pytest.approx(0.2)

    def test_zero_slopes_unbounded(self):
        move = max_ingate_move({"a": 0.5, "b": 0.0},
                               {"a": 0.0, "b": 0.0}, self.BANDS, +1.0)
        assert move == float("inf")

    def test_out_of_band_value_raises(self):
        with pytest.raises(ValueError, match="outside band"):
            max_ingate_move({"a": 1.2, "b": 0.0},
                            {"a": 1.0, "b": 0.0}, self.BANDS, +1.0)

    def test_bad_direction_raises(self):
        with pytest.raises(ValueError, match="direction"):
            max_ingate_move({"a": 0.5, "b": 0.0},
                            {"a": 1.0, "b": 0.0}, self.BANDS, 0.5)


class TestLowkeCounts:
    def test_counts_bins_one_and_two(self):
        read = SimpleNamespace(n_scored=np.array([0, 1, 1, 2, 3, 1, 2]))
        assert lowke_counts(read) == {"KE1_n": 3, "KE2_n": 2}

    def test_empty_bins(self):
        read = SimpleNamespace(n_scored=np.array([0, 3, 4]))
        assert lowke_counts(read) == {"KE1_n": 0, "KE2_n": 0}


def _write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


class TestKeOracleMismatches:
    COLS = ("KE1_mean", "KE2_mean")

    def test_passes_on_exact_match(self, tmp_path):
        p = tmp_path / "committed.csv"
        _write_csv(p, [{"label": "s1", "KE1_mean": "0.6435239722333034",
                        "KE2_mean": "0.5457367035721171"}])
        rows = {"s1": {"KE1_mean": 0.6435239722333034,
                       "KE2_mean": 0.5457367035721171}}
        assert ke_oracle_mismatches(rows, p, cols=self.COLS) == []

    def test_flags_numeric_mismatch(self, tmp_path):
        p = tmp_path / "committed.csv"
        _write_csv(p, [{"label": "s1", "KE1_mean": "0.6435",
                        "KE2_mean": "0.5457"}])
        rows = {"s1": {"KE1_mean": 0.6437, "KE2_mean": 0.5457}}
        problems = ke_oracle_mismatches(rows, p, cols=self.COLS)
        assert len(problems) == 1
        assert "s1.KE1_mean" in problems[0]

    def test_flags_missing_label(self, tmp_path):
        p = tmp_path / "committed.csv"
        _write_csv(p, [{"label": "pooled", "KE1_mean": "0.64",
                        "KE2_mean": "0.55"}])
        problems = ke_oracle_mismatches({}, p, cols=self.COLS)
        assert problems == ["pooled: missing from rescored rows"]

    def test_within_rtol_passes(self, tmp_path):
        p = tmp_path / "committed.csv"
        _write_csv(p, [{"label": "s1", "KE1_mean": "0.6435",
                        "KE2_mean": "0.5457"}])
        rows = {"s1": {"KE1_mean": 0.6435 * (1 + 1e-10),
                       "KE2_mean": 0.5457}}
        assert ke_oracle_mismatches(rows, p, cols=self.COLS) == []
