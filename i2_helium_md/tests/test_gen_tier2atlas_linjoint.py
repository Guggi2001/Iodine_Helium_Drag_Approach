"""Focused tests for the (a, τ) joint ring (plan §6.6).

Config-construction, guard and pure-helper tests only — **no run, no
trajectory**. These lock the frozen §6.6 design: the six-cell table, the
shared-seed CRN posture (which here forbids `seed` *and* `num_molecules`
from every diff, unlike §6.5), the free-form provenance inherited from the
ring, the LJ-P1 committed-artifact oracle, and the scorer's Pareto helper.
"""

from __future__ import annotations

import dataclasses
import warnings

import pytest

from scripts.gen_tier2atlas_linjoint import (
    EB_TAG,
    JOINT_MATRIX,
    JOINT_TWIN_COLS,
    JOINT_TWIN_ROWS,
    as_ring_cell,
    build_joint_cell,
    joint_run_dir_name,
    verify_crn_pairing,
    verify_joint_cell,
    verify_joint_preregistration,
)
from scripts.gen_tier2atlas_linring import N, SEED


def _build(label: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return build_joint_cell(label)


class TestFrozenMatrix:
    def test_six_cells_exactly_as_planned(self):
        got = {(c.label, c.a, c.tau_ps, c.e0_eV) for c in JOINT_MATRIX}
        assert got == {
            ("s37", 37.5, 4.8, 0.38),
            ("s40", 40.0, 4.8, 0.39),
            ("s425", 42.5, 4.8, 0.39),
            ("s45", 45.0, 4.8, 0.40),
            ("d2", 35.0, 6.4, 0.30),
            ("j1", 40.0, 6.4, 0.31),
        }

    def test_the_curve_and_the_contrast_are_both_present(self):
        # Four cells on the tau 4.8 a-curve (lr6 supplies a = 35) and two
        # on the tau 6.4 contrast (the clone supplies a = 42.5).
        assert sorted(c.a for c in JOINT_MATRIX if c.tau_ps == 4.8) == [
            37.5, 40.0, 42.5, 45.0]
        assert sorted(c.a for c in JOINT_MATRIX if c.tau_ps == 6.4) == [
            35.0, 40.0]

    def test_twin_rows_cover_every_cell(self):
        assert set(JOINT_TWIN_ROWS) == {c.label for c in JOINT_MATRIX}
        for row in JOINT_TWIN_ROWS.values():
            assert len(row) == len(JOINT_TWIN_COLS)

    def test_j1_twin_gate_is_zero_by_the_twin_side_floor(self):
        # Documented in the plan: twin nbar 4.392 < the TWIN-side floor 4.4,
        # while MD-side it lands inside [3.77, 4.37]. Locked so the flag is
        # never re-read as a predicted failure.
        row = dict(zip(JOINT_TWIN_COLS, JOINT_TWIN_ROWS["j1"]))
        assert row["gate"] == "0"
        assert float(row["nbar_det"]) == pytest.approx(4.392, abs=1e-9)
        assert float(row["nbar_det"]) - 0.55 < 4.37

    def test_namespace(self):
        for cell in JOINT_MATRIX:
            name = joint_run_dir_name(cell.label)
            assert "_tier2_" not in name and "tier2probe" not in name
            assert f"N{N}" in name and name.endswith(f"linj{cell.label}")


class TestCellBuild:
    @pytest.mark.parametrize("cell", JOINT_MATRIX, ids=lambda c: c.label)
    def test_builds_at_the_shared_crn_seed_with_its_own_knobs(self, cell):
        cfg = _build(cell.label)
        assert cfg.seed == SEED               # the whole-family pairing
        assert cfg.num_molecules == N
        assert cfg.internal_energy_cooling_tau_ps == cell.tau_ps
        assert dict(cfg.drag_coefficients.coefficients) == {"a": cell.a}

    def test_free_form_posture_inherited(self):
        cfg = _build("s425")
        coeffs = cfg.drag_coefficients
        assert cfg.drag_form == "pure_linear"
        assert coeffs.extraction_method == "free_form"
        assert coeffs.effective_binding_energy_I_ion_eV is None
        assert cfg.allow_unvalidated_binding_pairing is True

    def test_all_cells_share_one_well(self):
        assert {as_ring_cell(c).eb_tag for c in JOINT_MATRIX} == {EB_TAG}

    def test_unknown_cell_rejected(self):
        with pytest.raises(KeyError):
            build_joint_cell("nope")


class TestCrnGuard:
    """Shared-seed posture: neither key may differ (plan §6.6 oracle 2)."""

    @pytest.mark.parametrize("cell", JOINT_MATRIX, ids=lambda c: c.label)
    def test_seed_and_n_absent_from_every_diff(self, cell):
        diff = set(verify_crn_pairing(_build(cell.label), cell.label))
        assert "seed" not in diff and "num_molecules" not in diff
        assert diff  # the cell does differ from the partner in its own knobs

    def test_a_changed_seed_breaks_the_guard(self):
        cfg = dataclasses.replace(_build("s40"), seed=SEED + 1)
        with pytest.raises(AssertionError, match="CRN guard FAILED"):
            verify_crn_pairing(cfg, "s40")

    def test_a_changed_n_breaks_the_guard(self):
        cfg = dataclasses.replace(_build("s40"), num_molecules=250)
        with pytest.raises(AssertionError, match="CRN guard FAILED"):
            verify_crn_pairing(cfg, "s40")

    def test_full_cell_verification_passes(self):
        for cell in JOINT_MATRIX:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                cfg, diff = verify_joint_cell(cell.label)
            assert cfg.seed == SEED and diff


class TestLjP1Oracle:
    def test_frozen_rows_reproduce_from_the_committed_csv(self, capsys):
        verify_joint_preregistration()
        out = capsys.readouterr().out
        assert "LJ-P1 joint oracle PASSED" in out
        assert "LR-P1 twin oracle PASSED" in out


class TestParetoHelper:
    def test_front_minimises_w1_and_maximises_ke1(self):
        from scripts.post_processing.tier2atlas_linjoint_table import (
            pareto_front,
        )
        pts = [
            ("good_both", 0.70, 0.70),     # dominates everything below
            ("dominated", 0.80, 0.65),
            ("cheap_w1", 0.60, 0.50),      # best W1, worse KE -> on front
            ("hot_ke", 0.95, 0.90),        # best KE, worst W1 -> on front
        ]
        front = {lab for lab, _, _ in pareto_front(pts)}
        assert front == {"good_both", "cheap_w1", "hot_ke"}
        assert "dominated" not in front

    def test_non_finite_points_are_dropped(self):
        from scripts.post_processing.tier2atlas_linjoint_table import (
            pareto_front,
        )
        front = pareto_front([("a", 0.7, 0.7), ("b", float("nan"), 0.9)])
        assert [lab for lab, _, _ in front] == ["a"]
