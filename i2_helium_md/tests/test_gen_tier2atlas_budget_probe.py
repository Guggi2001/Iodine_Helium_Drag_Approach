"""Tests for the §3.5k budget-slope probe generator (pure config construction).

Covers `scripts/gen_tier2atlas_budget_probe.py` without touching any run dir
(testing rules): the built cell differs from the h405 finals config in
exactly the registered budget fields, the kinematic unit oracle reproduces
the registered budget/E0 arithmetic to 1e-9, and the run-dir namespace stays
in the atlas convention. The cfg-diff oracle (which reads the committed
run's ``cfg.json``) is exercised by ``--dry-run``, not here.
"""

from __future__ import annotations

import dataclasses

import pytest

from scripts.gen_tier2atlas_g4finals import _SPEC_BY_LABEL as FINALS_SPECS
from scripts.gen_tier2atlas_g4finals import build_cell
from scripts.gen_tier2atlas_budget_probe import (
    E0_PIN_EV,
    E_FRAG_REF_EV,
    F_INT_STANDING,
    RING,
    budget_run_dir_name,
    build_budget_cell,
    verify_budget_unit_oracle,
)


class TestCellArithmetic:
    def test_reference_budget_is_the_production_point_coulomb(self):
        # 14.39964548 / 2.666 / 2 — the §3.5k registered 2.7006 eV/fragment.
        assert E_FRAG_REF_EV == pytest.approx(2.700609, abs=5e-7)

    @pytest.mark.parametrize("cell", RING, ids=lambda c: c.label)
    def test_scale_roundtrips_to_registered_budget(self, cell):
        assert cell.scale * E_FRAG_REF_EV == pytest.approx(
            cell.budget_eV, abs=1e-12)

    def test_pinned_cells_hold_the_h405_onset(self):
        for cell in RING:
            if cell.pin_E0:
                assert cell.E0_eV == pytest.approx(E0_PIN_EV, abs=1e-12)

    def test_wiring_cell_scales_the_onset_proportionally(self):
        (wiring,) = [c for c in RING if not c.pin_E0]
        assert wiring.f_int == F_INT_STANDING
        assert wiring.E0_eV == pytest.approx(0.6165, abs=1e-12)

    def test_ring_is_the_registered_three(self):
        assert [(c.label, c.budget_eV, c.pin_E0) for c in RING] == [
            ("bud226k", 2.26, True),
            ("bud411k", 4.11, True),
            ("bud411", 4.11, False),
        ]


class TestBuildBudgetCell:
    @pytest.mark.parametrize("cell", RING, ids=lambda c: c.label)
    def test_differs_from_h405_only_in_registered_fields(self, cell):
        base = build_cell(FINALS_SPECS["h405"])
        cfg = build_budget_cell(cell)
        assert cfg.E_coulomb_scale == cell.scale
        assert cfg.coulomb_available_eV == cell.budget_eV
        assert cfg.internal_energy_partition_fraction == cell.f_int
        # The rest of the config is the h405 cell verbatim (CRN pairing:
        # same seed, same geometry, same drag knobs, same bundle).
        assert dataclasses.replace(
            cfg,
            E_coulomb_scale=base.E_coulomb_scale,
            coulomb_available_eV=base.coulomb_available_eV,
            internal_energy_partition_fraction=(
                base.internal_energy_partition_fraction),
        ) == base

    @pytest.mark.parametrize("cell", RING, ids=lambda c: c.label)
    def test_kinematic_unit_oracle_passes(self, cell):
        verify_budget_unit_oracle(build_budget_cell(cell), cell)

    def test_unit_oracle_rejects_geometry_moved_budget(self):
        # The registered guard: the budget must move through E_coulomb_scale,
        # never through the birth separation.
        cell = RING[0]
        cfg = dataclasses.replace(build_budget_cell(cell), R0_GS_angstrom=3.0)
        with pytest.raises(AssertionError, match="birth geometry"):
            verify_budget_unit_oracle(cfg, cell)

    def test_seed_is_the_finals_seed(self):
        assert build_budget_cell(RING[0]).seed == 20260729


class TestRunDirName:
    @pytest.mark.parametrize("cell", RING, ids=lambda c: c.label)
    def test_atlas_namespace(self, cell):
        name = budget_run_dir_name(cell.label)
        assert "tier2atlas" in name and cell.label in name
        assert "_tier2_" not in name and "tier2probe" not in name
