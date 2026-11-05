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
    atlas_run_dir_name,
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


def test_allowed_diff_set_is_the_preregistered_four():
    assert ALLOWED_CFG_DIFF_KEYS == {
        "drag_form",
        "drag_coefficients",
        "binding_energy_I_ion_eV",
        "internal_energy_cooling_tau_ps",
    }
