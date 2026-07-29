"""Tests for the §3.5i s(n) probe generator (pure config construction).

Covers `scripts/gen_tier2atlas_sn_probe.py` without touching any run dir
(testing rules): the built cell differs from the h405 finals config in
exactly the three coupling fields, the unit oracle reproduces the
registered s-table, and the run-dir namespace stays in the atlas
convention. The cfg-diff oracle (which reads the committed run's
``cfg.json``) is exercised by ``--dry-run``, not here.
"""

from __future__ import annotations

import dataclasses

import pytest

from scripts.gen_tier2atlas_g4finals import _SPEC_BY_LABEL as FINALS_SPECS
from scripts.gen_tier2atlas_g4finals import build_cell
from scripts.gen_tier2atlas_sn_probe import (
    R_CORE_ANGSTROM,
    RING,
    build_sn_cell,
    sn_run_dir_name,
    verify_sn_unit_oracle,
)


class TestBuildSnCell:
    @pytest.mark.parametrize("cell", RING, ids=lambda c: c.label)
    def test_differs_from_h405_only_in_coupling_fields(self, cell):
        base = build_cell(FINALS_SPECS["h405"])
        cfg = build_sn_cell(cell)
        assert cfg.drag_state_coupling == "shell_area"
        assert cfg.state_coupling_R_core_angstrom == R_CORE_ANGSTROM
        assert cfg.state_coupling_rho_shell_per_A3 == cell.rho_shell_per_A3
        # The rest of the config is the h405 cell verbatim (CRN pairing:
        # same seed, same geometry, same knobs, same bundle).
        assert dataclasses.replace(
            cfg,
            drag_state_coupling=base.drag_state_coupling,
            state_coupling_R_core_angstrom=base.state_coupling_R_core_angstrom,
            state_coupling_rho_shell_per_A3=base.state_coupling_rho_shell_per_A3,
        ) == base

    @pytest.mark.parametrize("cell", RING, ids=lambda c: c.label)
    def test_unit_oracle_reproduces_registered_table(self, cell):
        verify_sn_unit_oracle(build_sn_cell(cell), cell)

    def test_seed_is_the_finals_seed(self):
        assert build_sn_cell(RING[0]).seed == 20260729

    def test_rho_values_are_the_registered_bounded_range(self):
        assert [c.rho_shell_per_A3 for c in RING] == [0.0218, 0.030, 0.0436]


class TestRunDirName:
    @pytest.mark.parametrize("cell", RING, ids=lambda c: c.label)
    def test_atlas_namespace(self, cell):
        name = sn_run_dir_name(cell.label)
        assert "tier2atlas" in name and cell.label in name
        assert "_tier2_" not in name and "tier2probe" not in name

    def test_ring_labels_are_the_registered_three(self):
        from scripts.gen_tier2atlas_sn_probe import _SPEC_BY_LABEL
        assert sorted(_SPEC_BY_LABEL) == ["sa22", "sa30", "sa44"]
