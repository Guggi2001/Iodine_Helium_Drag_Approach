"""Tests for the §3.5l (C) probe generator (pure config construction).

Covers `scripts/gen_tier2atlas_ce_probe.py` without touching any run dir
(testing rules): each built cell differs from the h405 finals config in
exactly its (C) surfaces, the frozen §8 registration inputs are stamped
verbatim, the STOPPED s(n) axis stays off, the unit oracle reproduces the
registered landmark values, and the run-dir namespace stays in the atlas
convention. The cfg-diff oracle (which reads the committed run's
``cfg.json``) is exercised by ``--dry-run``, not here.
"""

from __future__ import annotations

import dataclasses

import pytest

from scripts.gen_tier2atlas_ce_probe import (
    F_INT_Q2,
    F_INT_Q3_FULL,
    F_INT_Q3_PINNED,
    F_INT_SINGLE,
    FROZEN_E_SINGLE,
    FROZEN_F,
    FROZEN_SIGMAS,
    FROZEN_WEIGHTS,
    RING,
    STRIP_BOX_CENTER,
    build_ce_cell,
    ce_run_dir_name,
    verify_ce_unit_oracle,
)
from scripts.gen_tier2atlas_g4finals import _SPEC_BY_LABEL as FINALS_SPECS
from scripts.gen_tier2atlas_g4finals import build_cell


def _reset_ce_fields(cfg, base):
    """Copy of ``cfg`` with every (C) surface reset to ``base``'s values."""
    return dataclasses.replace(
        cfg,
        ce_channel_mode=base.ce_channel_mode,
        ce_channel_weights=base.ce_channel_weights,
        ce_fraction_f=base.ce_fraction_f,
        ce_channel_sigma_eV=base.ce_channel_sigma_eV,
        ce_single_ker_eV=base.ce_single_ker_eV,
        ce_q3_partner_mask=base.ce_q3_partner_mask,
        ce_internal_energy_partition_fractions=(
            base.ce_internal_energy_partition_fractions
        ),
        exit_strip_mode=base.exit_strip_mode,
        exit_strip_v_ref=base.exit_strip_v_ref,
        exit_strip_exponent=base.exit_strip_exponent,
        exit_strip_protect_j0=base.exit_strip_protect_j0,
        exit_strip_width_rungs=base.exit_strip_width_rungs,
        exit_strip_carry_eV=base.exit_strip_carry_eV,
    )


class TestBuildCeCell:
    @pytest.mark.parametrize("cell", RING, ids=lambda c: c.label)
    def test_differs_from_h405_only_in_ce_fields(self, cell):
        base = build_cell(FINALS_SPECS["h405"])
        cfg = build_ce_cell(cell)
        # The rest of the config is the h405 cell verbatim (PC-5: A-only is
        # CRN-paired; the STOPPED s(n) axis stays off).
        assert cfg.drag_state_coupling == "off"
        assert _reset_ce_fields(cfg, base) == base

    def test_cfull_stamps_the_frozen_registration_inputs(self):
        cfg = build_ce_cell([c for c in RING if c.label == "cfull"][0])
        assert cfg.ce_channel_mode == "sampled"
        assert tuple(cfg.ce_channel_weights) == FROZEN_WEIGHTS == (0.30, 0.50, 0.20)
        assert cfg.ce_fraction_f == FROZEN_F == 0.80
        assert tuple(cfg.ce_channel_sigma_eV) == FROZEN_SIGMAS == (0.42, 0.31, 0.55)
        assert cfg.ce_single_ker_eV == FROZEN_E_SINGLE == 0.53
        assert cfg.ce_q3_partner_mask is True
        assert tuple(cfg.ce_internal_energy_partition_fractions) == (
            F_INT_SINGLE, F_INT_Q2, F_INT_Q3_PINNED,
        ) == (0.15, 0.15, 0.0985)
        assert cfg.exit_strip_mode == "depth_graded"
        for name, value in STRIP_BOX_CENTER.items():
            assert getattr(cfg, name) == value

    def test_aonly_keeps_the_scalar_budget(self):
        cfg = build_ce_cell([c for c in RING if c.label == "aonly"][0])
        assert cfg.ce_channel_mode == "off"
        assert cfg.exit_strip_mode == "depth_graded"
        assert cfg.coulomb_available_eV == 2.70
        assert cfg.internal_energy_partition_fraction == pytest.approx(
            0.405 / 2.70,
        )

    def test_bonly_has_no_strip(self):
        cfg = build_ce_cell([c for c in RING if c.label == "bonly"][0])
        assert cfg.ce_channel_mode == "sampled"
        assert cfg.exit_strip_mode == "off"

    def test_coupling_arm_at_the_full_end(self):
        cfg = build_ce_cell([c for c in RING if c.label == "cq3hi"][0])
        assert cfg.ce_internal_energy_partition_fractions[2] \
            == F_INT_Q3_FULL == 0.15

    def test_seed_is_the_finals_seed(self):
        assert build_ce_cell(RING[0]).seed == 20260729

    def test_unit_oracle_passes(self):
        verify_ce_unit_oracle()


class TestRunDirName:
    @pytest.mark.parametrize("cell", RING, ids=lambda c: c.label)
    def test_atlas_namespace(self, cell):
        name = ce_run_dir_name(cell.label)
        assert "tier2atlas" in name and cell.label in name
        assert "_tier2_" not in name and "tier2probe" not in name

    def test_ring_labels_are_the_registered_four(self):
        assert [c.label for c in RING] == ["cfull", "aonly", "bonly", "cq3hi"]
