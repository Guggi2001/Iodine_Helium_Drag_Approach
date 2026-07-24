"""Unit coverage for scripts/gen_tier2atlas_ebindscan.py (atlas §6.7 item 2).

Config-build and guard wiring only — no propagation. The qcc-diff test needs
the on-disk qcc cell (``data/runs`` is uncommitted) and skips when absent.
"""

from __future__ import annotations

import pytest

from scripts.gen_tier2atlas_ebindscan import (
    ALLOWED_CFG_DIFF_KEYS,
    CUBIC_PAIR_E_BIND_EV,
    EBIND_MATRIX,
    N,
    PROJECT_ROOT,
    QCC_RUN_DIR_NAME,
    SHARED_LQ_E_BIND_EV,
    SEED,
    V_C_APS,
    atlas_run_dir_name,
    build_cell,
    verify_against_qcc,
)

_QCC_CFG = PROJECT_ROOT / "data" / "runs" / QCC_RUN_DIR_NAME / "cfg.json"


class TestNamespace:
    def test_atlas_tag_locks(self):
        for spec in EBIND_MATRIX:
            name = atlas_run_dir_name(spec.label)
            assert "_tier2atlas_" in name
            assert "_tier2_" not in name
            assert "tier2probe" not in name
            assert "qcc_" in name

    def test_matrix_values(self):
        """Two cells: the exact cubic-pair well and the 0.154 spread top."""
        by = {s.label: s.binding_eV for s in EBIND_MATRIX}
        assert by["eb1168"] == CUBIC_PAIR_E_BIND_EV
        assert by["eb154"] == 0.154
        # eb1168 must be byte-matched to the standing cubic pair for the
        # form-at-matched-well comparison against finc1v725.
        assert CUBIC_PAIR_E_BIND_EV == 0.11675778353879479


class TestCellBuild:
    @pytest.mark.parametrize("spec", EBIND_MATRIX, ids=lambda s: s.label)
    def test_cell_builds_and_validates(self, spec):
        cfg = build_cell(spec)
        cfg.validate()  # passes only because the hatch is set
        # The well is overridden; the drag law + its provenance stamp are NOT.
        assert cfg.binding_energy_I_ion_eV == spec.binding_eV
        assert cfg.allow_unvalidated_binding_pairing is True
        assert cfg.drag_coefficients.effective_binding_energy_I_ion_eV == \
            pytest.approx(SHARED_LQ_E_BIND_EV)
        assert cfg.drag_form == "capped_linear_quadratic"
        assert cfg.drag_coefficients.coefficients["v_c"] == V_C_APS
        assert cfg.seed == SEED
        assert cfg.num_molecules == N

    def test_well_actually_differs_from_stamp(self):
        """The whole point: the climbed well ≠ the co-extracted stamp."""
        for spec in EBIND_MATRIX:
            cfg = build_cell(spec)
            assert cfg.binding_energy_I_ion_eV != pytest.approx(
                cfg.drag_coefficients.effective_binding_energy_I_ion_eV
            )

    def test_landau_arm_and_mass_flags(self):
        cfg = build_cell(EBIND_MATRIX[0])
        assert cfg.relaxation_dissipation == "landau_gated_drag"
        assert cfg.v_limit_m_per_s == 58.0
        assert cfg.mass_scenario == "biphasic"
        assert cfg.allow_inconsistent_mass_pairing is True


@pytest.mark.skipif(
    not _QCC_CFG.exists(),
    reason="on-disk qcc cell not present (data/runs is uncommitted)",
)
class TestQccDiff:
    @pytest.mark.parametrize("spec", EBIND_MATRIX, ids=lambda s: s.label)
    def test_diff_vs_qcc_is_exactly_well_plus_hatch(self, spec):
        verify_against_qcc(build_cell(spec), spec)

    def test_allowed_diff_set(self):
        assert ALLOWED_CFG_DIFF_KEYS == {
            "binding_energy_I_ion_eV",
            "allow_unvalidated_binding_pairing",
        }
