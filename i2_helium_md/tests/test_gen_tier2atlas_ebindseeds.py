"""Unit coverage for scripts/gen_tier2atlas_ebindseeds.py (§6.7 item-2 firm-up).

Config-build + guard wiring only. The qccbigs-diff test needs the item-1
battery run dirs (``data/runs`` uncommitted) and skips when absent.
"""

from __future__ import annotations

import pytest

from scripts.gen_tier2atlas_ebindseeds import (
    ALLOWED_CFG_DIFF_KEYS,
    CUBIC_PAIR_E_BIND_EV,
    EBIND_MATRIX,
    N,
    PROJECT_ROOT,
    SHARED_LQ_E_BIND_EV,
    V_C_APS,
    WELLS,
    atlas_run_dir_name,
    build_cell,
    verify_against_bigs,
)


class TestMatrix:
    def test_six_cells_two_wells_three_seeds(self):
        assert len(EBIND_MATRIX) == 6
        assert sorted(set(c.well_label for c in EBIND_MATRIX)) == ["eb1168", "eb154"]
        assert sorted(set(c.seed for c in EBIND_MATRIX)) == [
            20260722, 20260723, 20260724]

    def test_wells_and_exact_cubic_pair(self):
        assert WELLS["eb1168"] == CUBIC_PAIR_E_BIND_EV == 0.11675778353879479
        assert WELLS["eb154"] == 0.154

    def test_namespace_locks(self):
        for cell in EBIND_MATRIX:
            name = atlas_run_dir_name(cell.label)
            assert "_tier2atlas_" in name
            assert "_tier2_" not in name and "tier2probe" not in name
            assert f"_N{N}_" in name

    def test_pairs_map_seed_to_matched_bigs(self):
        for cell in EBIND_MATRIX:
            assert cell.paired_bigs.endswith(f"bigs{[20260722, 20260723, 20260724].index(cell.seed) + 1}")
            assert f"_N{N}_" in cell.paired_bigs


class TestCellBuild:
    @pytest.mark.parametrize("cell", EBIND_MATRIX, ids=lambda c: c.label)
    def test_builds_validates_well_override(self, cell):
        cfg = build_cell(cell)
        cfg.validate()
        assert cfg.binding_energy_I_ion_eV == cell.binding_eV
        assert cfg.allow_unvalidated_binding_pairing is True
        # drag law + provenance stamp unchanged (the well is NOT re-paired)
        assert cfg.drag_coefficients.effective_binding_energy_I_ion_eV == \
            pytest.approx(SHARED_LQ_E_BIND_EV)
        assert cfg.binding_energy_I_ion_eV != pytest.approx(SHARED_LQ_E_BIND_EV)
        assert cfg.drag_coefficients.coefficients["v_c"] == V_C_APS
        assert cfg.num_molecules == N

    def test_same_well_across_seeds_only_seed_differs(self):
        eb1168 = [c for c in EBIND_MATRIX if c.well_label == "eb1168"]
        cfgs = [build_cell(c) for c in eb1168]
        assert len({c.seed for c in cfgs}) == 3
        for c in cfgs[1:]:
            assert c.binding_energy_I_ion_eV == cfgs[0].binding_energy_I_ion_eV


def _bigs_present() -> bool:
    return all(
        (PROJECT_ROOT / "data" / "runs" / c.paired_bigs / "cfg.json").exists()
        for c in EBIND_MATRIX
    )


@pytest.mark.skipif(not _bigs_present(),
                    reason="item-1 battery run dirs absent (data/runs uncommitted)")
class TestBigsDiff:
    @pytest.mark.parametrize("cell", EBIND_MATRIX, ids=lambda c: c.label)
    def test_diff_vs_paired_bigs_is_well_plus_hatch(self, cell):
        verify_against_bigs(build_cell(cell), cell)

    def test_allowed_diff_set(self):
        assert ALLOWED_CFG_DIFF_KEYS == {
            "binding_energy_I_ion_eV", "allow_unvalidated_binding_pairing"}
