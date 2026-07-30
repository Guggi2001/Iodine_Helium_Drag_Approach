"""Tests for the free-form linear sweep Step-0 helpers (``ke1auth``).

Covers `scripts/tier2_h2b_forward_model.py`:

- ``spearman_rho`` — rank correlation with average-rank ties; loud failure
  on degenerate input (a licensure verdict must never ride on NaN);
- ``ke1auth_select_rows`` — the frozen replay-set rule of
  `TIER2_FREEFORM_LINEAR_TWIN_SWEEP_PLAN.md` §2 (group filter, pooled
  preference, NaN-KE1 drop, duplicate-pin dedupe to the largest-N group);
- ``ke1auth_verdict`` — the pre-registered licensure bands.

Synthetic inputs only — no run dirs, no chord integrations, no figures
(testing rules). Tolerances: exact arithmetic identities where possible,
``pytest.approx`` defaults for the rank-correlation values.
"""

from __future__ import annotations

import numpy as np
import pytest

from scripts.tier2_h2b_forward_model import (
    ke1auth_select_rows,
    ke1auth_verdict,
    spearman_rho,
)


def _row(group, label, v_c="5.5", well="eb1168", tau="4.8", E0="0.31",
         ke1="0.8", ke2="0.7"):
    return {"group": group, "label": label, "v_c": v_c, "well": well,
            "tau": tau, "E0": E0, "KE1_mean": ke1, "KE2_mean": ke2}


class TestSpearman:
    def test_perfect_monotonic_is_one(self):
        assert spearman_rho([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
        # Non-linear but monotonic: rank correlation stays exactly 1.
        assert spearman_rho([1, 2, 3, 4], [1, 8, 27, 64]) == pytest.approx(1.0)

    def test_reversed_is_minus_one(self):
        assert spearman_rho([1, 2, 3], [5, 4, 3]) == pytest.approx(-1.0)

    def test_ties_use_average_ranks(self):
        # x = [1, 2, 2, 3] -> ranks [1, 2.5, 2.5, 4]; compare against the
        # Pearson correlation of the hand-built rank vectors.
        expected = np.corrcoef([1.0, 2.5, 2.5, 4.0], [1, 2, 3, 4])[0, 1]
        assert spearman_rho([1, 2, 2, 3], [1, 2, 3, 4]) == pytest.approx(expected)

    def test_degenerate_inputs_fail_loudly(self):
        with pytest.raises(ValueError):
            spearman_rho([1.0], [2.0])  # size < 2
        with pytest.raises(ValueError):
            spearman_rho([1, 2], [3])  # shape mismatch
        with pytest.raises(ValueError):
            spearman_rho([1, 1, 1], [1, 2, 3])  # constant input
        with pytest.raises(ValueError):
            spearman_rho([1, np.nan, 3], [1, 2, 3])  # non-finite


class TestVerdict:
    def test_bands(self):
        assert ke1auth_verdict(0.95) == "LICENSED"
        assert ke1auth_verdict(0.8) == "LICENSED"  # boundary inclusive
        assert ke1auth_verdict(0.79) == "DIRECTIONAL_ONLY"
        assert ke1auth_verdict(0.5) == "DIRECTIONAL_ONLY"
        assert ke1auth_verdict(0.49) == "UNLICENSED"
        assert ke1auth_verdict(-0.3) == "UNLICENSED"

    def test_nan_fails_loudly(self):
        with pytest.raises(ValueError):
            ke1auth_verdict(float("nan"))


class TestSelectRows:
    def test_non_corrected_groups_excluded(self):
        kept, excluded = ke1auth_select_rows(
            [_row("incumbent", "s1"), _row("geogrid", "r1l3"),
             _row("g3ring", "a031")])
        assert [r["label"] for r in kept] == ["a031"]
        assert {r["label"] for r, _ in excluded} == {"s1", "r1l3"}

    def test_battery_members_dropped_for_pooled(self):
        kept, excluded = ke1auth_select_rows(
            [_row("h405bat", "s1", E0="0.405"),
             _row("h405bat", "pooled", E0="0.405")])
        assert [(r["group"], r["label"]) for r in kept] == [("h405bat", "pooled")]
        assert excluded[0][1] == "battery member: pooled row preferred"

    def test_nan_ke1_dropped_with_reason(self):
        kept, excluded = ke1auth_select_rows([_row("g3ring", "f725", ke1="nan")])
        assert kept == []
        assert excluded[0][1] == "MD KE1 unmeasured (empty n=1 bin)"

    def test_duplicate_pins_prefer_larger_n_group_either_order(self):
        ring = _row("g3ring", "a037", E0="0.37")
        finals = _row("g4finals", "a037", E0="0.37")
        for order in ([ring, finals], [finals, ring]):
            kept, excluded = ke1auth_select_rows(order)
            assert [(r["group"], r["label"]) for r in kept] == [("g4finals", "a037")]
            assert excluded[0][0]["group"] == "g3ring"
            assert "duplicate pins" in excluded[0][1]

    def test_distinct_pins_both_kept_and_ordered_by_preference(self):
        kept, _ = ke1auth_select_rows(
            [_row("g3ring", "a031", E0="0.31"),
             _row("h405bat", "pooled", tau="4.4", E0="0.405"),
             _row("g4finals", "f1", E0="0.365")])
        assert [r["group"] for r in kept] == ["h405bat", "g4finals", "g3ring"]
