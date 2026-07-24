"""Unit coverage for scripts/gen_tier2atlas_spotcheck.py (atlas §6.6 cells).

Config-build and guard wiring only — no propagation (the expensive stages
are exercised by the generator itself, mirroring the tier0/tier2 generator
test convention). The cfg-diff test needs the standing finc1v725 run dir
(``data/runs`` is not committed) and skips when it is absent.
"""

from __future__ import annotations

import pytest

from scripts.gen_tier2atlas_spotcheck import (
    ALLOWED_CFG_DIFF_KEYS,
    FINC_RUN_DIR_NAME,
    PROJECT_ROOT,
    SPOTCHECK_MATRIX,
    atlas_run_dir_name,
    build_cell,
    verify_against_standing,
)

_FINC_CFG = PROJECT_ROOT / "data" / "runs" / FINC_RUN_DIR_NAME / "cfg.json"


class TestNamespace:
    def test_atlas_tag_locks(self):
        """No atlas dir name may match the campaign or probe globs (§1.4)."""
        for spec in SPOTCHECK_MATRIX:
            name = atlas_run_dir_name(spec.label)
            assert "_tier2atlas_" in name
            assert "_tier2_" not in name
            assert "tier2probe" not in name

    def test_labels_are_tag_charset(self):
        for spec in SPOTCHECK_MATRIX:
            assert spec.label.isascii() and spec.label.isalnum()
            assert spec.label == spec.label.lower()


class TestCellBuild:
    @pytest.mark.parametrize("spec", SPOTCHECK_MATRIX, ids=lambda s: s.label)
    def test_cell_builds_and_validates(self, spec):
        cfg = build_cell(spec)
        cfg.validate()  # raises on any guard failure
        dc = cfg.drag_coefficients
        assert cfg.drag_form == "capped_linear_quadratic"
        assert dc.form == "capped_linear_quadratic"
        # a, c merged from the shared_lq bundle; tail from the spec.
        assert dc.coefficients["a"] == pytest.approx(9.805022771384936e-05)
        assert dc.coefficients["c"] == pytest.approx(12.792201727123915)
        assert dc.coefficients["v_c"] == spec.v_c_aps
        assert dc.coefficients["p_tail"] == -1.0
        assert cfg.internal_energy_cooling_tau_ps == spec.tau_ps

    def test_binding_pairing_is_consistent(self):
        """The lq cells run the jointly-extracted §6.5.1 pair: the cfg binding
        equals the shared_lq bundle's stamped effective binding (0.0482 eV),
        so no unvalidated-binding escape hatch is involved."""
        cfg = build_cell(SPOTCHECK_MATRIX[0])
        assert cfg.binding_energy_I_ion_eV == pytest.approx(
            cfg.drag_coefficients.effective_binding_energy_I_ion_eV
        )
        assert cfg.binding_energy_I_ion_eV == pytest.approx(
            0.048236582655347665
        )

    def test_landau_arm_and_mass_pairing_flags(self):
        cfg = build_cell(SPOTCHECK_MATRIX[0])
        assert cfg.relaxation_dissipation == "landau_gated_drag"
        assert cfg.v_limit_m_per_s == 58.0
        # biphasic + constant-mass coefficients: the §6.6 documented exception.
        assert cfg.mass_scenario == "biphasic"
        assert cfg.allow_inconsistent_mass_pairing is True


@pytest.mark.skipif(
    not _FINC_CFG.exists(),
    reason="standing finc1v725 run dir not present (data/runs is uncommitted)",
)
class TestStandingDiff:
    @pytest.mark.parametrize("spec", SPOTCHECK_MATRIX, ids=lambda s: s.label)
    def test_diff_vs_finc1v725_is_exactly_preregistered(self, spec):
        verify_against_standing(build_cell(spec), spec)

    def test_allowed_diff_set_is_the_preregistered_four(self):
        assert ALLOWED_CFG_DIFF_KEYS == {
            "drag_form",
            "drag_coefficients",
            "binding_energy_I_ion_eV",
            "internal_energy_cooling_tau_ps",
        }
