"""Tests for the §3.5h p_tail ring generator (pure config construction).

Covers `scripts/gen_tier2atlas_ptail_ring.py` without touching any run dir
(testing rules): the built cell differs from the h405 finals config in
exactly the tail exponent, the widened guard accepts the ring values, and
the run-dir namespace stays in the atlas convention.
"""

from __future__ import annotations

import dataclasses

import pytest

from scripts.gen_tier2atlas_g4finals import _SPEC_BY_LABEL as FINALS_SPECS
from scripts.gen_tier2atlas_g4finals import build_cell
from scripts.gen_tier2atlas_ptail_ring import (
    RING,
    build_ptail_cell,
    ptail_run_dir_name,
)


class TestBuildPtailCell:
    @pytest.mark.parametrize("cell", RING, ids=lambda c: c.label)
    def test_differs_from_h405_only_in_p_tail(self, cell):
        base = build_cell(FINALS_SPECS["h405"])
        cfg = build_ptail_cell(cell)
        assert cfg.drag_coefficients.coefficients["p_tail"] == cell.p_tail
        # Everything else in the bundle is bit-identical.
        base_c = dict(base.drag_coefficients.coefficients)
        got_c = dict(cfg.drag_coefficients.coefficients)
        base_c.pop("p_tail"), got_c.pop("p_tail")
        assert got_c == base_c
        assert dataclasses.replace(
            cfg.drag_coefficients, coefficients=base.drag_coefficients.coefficients
        ) == base.drag_coefficients
        # And the rest of the config is the h405 cell verbatim (CRN pairing:
        # same seed, same geometry, same knobs).
        assert dataclasses.replace(
            cfg, drag_coefficients=base.drag_coefficients) == base

    @pytest.mark.parametrize("cell", RING, ids=lambda c: c.label)
    def test_widened_guard_accepts_ring_values(self, cell):
        # cfg.validate() runs inside build_ptail_cell — reaching here means
        # the §3.5h guard admitted the exponent; assert it explicitly too.
        cfg = build_ptail_cell(cell)
        assert -4.0 <= cfg.drag_coefficients.coefficients["p_tail"] <= 0.0

    def test_seed_is_the_finals_seed(self):
        assert build_ptail_cell(RING[0]).seed == 20260729


class TestRunDirName:
    @pytest.mark.parametrize("cell", RING, ids=lambda c: c.label)
    def test_atlas_namespace(self, cell):
        name = ptail_run_dir_name(cell.label)
        assert "tier2atlas" in name and cell.label in name
        assert "_tier2_" not in name and "tier2probe" not in name

    def test_ring_labels_are_the_registered_three(self):
        from scripts.gen_tier2atlas_ptail_ring import _SPEC_BY_LABEL
        assert sorted(_SPEC_BY_LABEL) == ["pt15", "pt20", "pt30"]
