"""Unit coverage for scripts/gen_tier2atlas_lqbattery.py (atlas §6.7 battery).

Config-build and guard wiring only — no propagation (the expensive stages
are exercised by the generator itself, mirroring the tier0/tier2 generator
test convention). The paired-diff test needs the cubic battery run dirs
(``data/runs`` is not committed) and skips when they are absent.
"""

from __future__ import annotations

import pytest

from scripts.gen_tier2atlas_lqbattery import (
    ALLOWED_CFG_DIFF_KEYS,
    BATTERY_MATRIX,
    N,
    PROJECT_ROOT,
    SHARED_LQ_E_BIND_EV,
    TAU_PS,
    V_C_APS,
    atlas_run_dir_name,
    build_cell,
    verify_against_paired,
)


class TestNamespace:
    def test_atlas_tag_locks(self):
        """No atlas dir name may match the campaign or probe globs (§1.4)."""
        for spec in BATTERY_MATRIX:
            name = atlas_run_dir_name(spec.label)
            assert "_tier2atlas_" in name
            assert "_tier2_" not in name
            assert "tier2probe" not in name
            assert f"_N{N}_" in name

    def test_labels_are_tag_charset(self):
        for spec in BATTERY_MATRIX:
            assert spec.label.isascii() and spec.label.isalnum()
            assert spec.label == spec.label.lower()

    def test_battery_is_five_paired_seeds(self):
        """The five cubic-battery seeds, one member each, paired 1:1."""
        assert len(BATTERY_MATRIX) == 5
        assert [s.seed for s in BATTERY_MATRIX] == [
            20260722, 20260723, 20260724, 20260725, 20260726
        ]
        for k, spec in enumerate(BATTERY_MATRIX, start=1):
            assert spec.paired_cubic_run_dir.endswith(f"bigc1v725s{k}")
            assert f"_N{N}_" in spec.paired_cubic_run_dir


class TestCellBuild:
    @pytest.mark.parametrize("spec", BATTERY_MATRIX, ids=lambda s: s.label)
    def test_cell_builds_and_validates(self, spec):
        cfg = build_cell(spec)
        cfg.validate()  # raises on any guard failure
        dc = cfg.drag_coefficients
        assert cfg.drag_form == "capped_linear_quadratic"
        assert dc.form == "capped_linear_quadratic"
        # a, c merged from the shared_lq bundle (same as the §6.6 spot-check).
        assert dc.coefficients["a"] == pytest.approx(9.805022771384936e-05)
        assert dc.coefficients["c"] == pytest.approx(12.792201727123915)
        assert dc.coefficients["v_c"] == V_C_APS
        assert dc.coefficients["p_tail"] == -1.0
        assert cfg.internal_energy_cooling_tau_ps == TAU_PS
        assert cfg.num_molecules == N
        assert cfg.seed == spec.seed

    def test_binding_pairing_is_consistent(self):
        """The lq battery runs the jointly-extracted §6.5.1 pair: the cfg
        binding equals the shared_lq bundle's stamped effective binding
        (0.0482 eV), so no unvalidated-binding escape hatch is involved."""
        cfg = build_cell(BATTERY_MATRIX[0])
        assert cfg.binding_energy_I_ion_eV == pytest.approx(
            cfg.drag_coefficients.effective_binding_energy_I_ion_eV
        )
        assert cfg.binding_energy_I_ion_eV == pytest.approx(SHARED_LQ_E_BIND_EV)

    def test_landau_arm_and_mass_pairing_flags(self):
        cfg = build_cell(BATTERY_MATRIX[0])
        assert cfg.relaxation_dissipation == "landau_gated_drag"
        assert cfg.v_limit_m_per_s == 58.0
        # biphasic + constant-mass coefficients: the §6.6 documented exception.
        assert cfg.mass_scenario == "biphasic"
        assert cfg.allow_inconsistent_mass_pairing is True

    def test_members_differ_only_by_seed(self):
        """Every member is the qcc cell; only the seed varies across the five."""
        cfgs = [build_cell(s) for s in BATTERY_MATRIX]
        seeds = {c.seed for c in cfgs}
        assert len(seeds) == 5
        for c in cfgs[1:]:
            assert c.drag_coefficients.coefficients == \
                cfgs[0].drag_coefficients.coefficients
            assert c.internal_energy_cooling_tau_ps == \
                cfgs[0].internal_energy_cooling_tau_ps
            assert c.binding_energy_I_ion_eV == cfgs[0].binding_energy_I_ion_eV


def _paired_cfgs_present() -> bool:
    return all(
        (PROJECT_ROOT / "data" / "runs" / s.paired_cubic_run_dir / "cfg.json").exists()
        for s in BATTERY_MATRIX
    )


@pytest.mark.skipif(
    not _paired_cfgs_present(),
    reason="cubic battery run dirs not present (data/runs is uncommitted)",
)
class TestPairedDiff:
    @pytest.mark.parametrize("spec", BATTERY_MATRIX, ids=lambda s: s.label)
    def test_diff_vs_paired_cubic_is_exactly_preregistered(self, spec):
        verify_against_paired(build_cell(spec), spec)

    def test_allowed_diff_set_is_the_preregistered_four(self):
        assert ALLOWED_CFG_DIFF_KEYS == {
            "drag_form",
            "drag_coefficients",
            "binding_energy_I_ion_eV",
            "internal_energy_cooling_tau_ps",
        }
