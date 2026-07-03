"""Tests for postprocess/abundance_loader.py (Tier-2 Phase E, Slice E3).

The E3 loader reads the experimental I+He_n abundance reference
(``data/reference/integrated_i_he_abundance.csv`` -- the Tier-2 arbiter) under
the same validate-early contract as the HeDFT loader: order-independent header
membership match (missing *and* extra columns raise), loud FileNotFoundError /
ValueError, ``source_path = p.resolve()``. Real CSV loaded once; every
validation arm fired on a synthetic corruption.
"""

from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.postprocess.abundance_loader import (
    HeAbundanceReference,
    load_he_abundance_reference,
)

_REAL_CSV = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "reference"
    / "integrated_i_he_abundance.csv"
)

_HEADER = (
    "n,label,massCenter_u_per_e,massWindowLower_u_per_e,"
    "massWindowUpper_u_per_e,ionCounts,ionPercent"
)

# A tiny, self-consistent 3-row reference (ionPercent sums to 100).
_GOOD_ROWS = [
    ("0", "I^+", "127", "126", "128", "0.5", "50.0"),
    ("1", "I^+He_1", "131.0026", "130.0026", "132.0026", "0.3", "30.0"),
    ("2", "I^+He_2", "135.0052", "134.0052", "136.0052", "0.2", "20.0"),
]


def _write_csv(path: Path, rows, *, header: str = _HEADER) -> Path:
    lines = [header] + [",".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Real reference
# ---------------------------------------------------------------------------
class TestRealReference:
    def test_loads_real_csv(self):
        ref = load_he_abundance_reference(_REAL_CSV)
        assert isinstance(ref, HeAbundanceReference)
        np.testing.assert_array_equal(ref.n, np.arange(21))
        assert ref.label[0] == "I^+"
        assert ref.label[20] == "I^+He_20"
        assert ref.mass_center_u[0] == pytest.approx(127.0)
        assert ref.mass_center_u[1] == pytest.approx(131.0026)
        # Spot-checks from the plan (ionPercent -> fraction is /100 here since
        # the file sums to exactly 100).
        assert ref.ion_fraction[0] == pytest.approx(0.435234834342411)
        assert ref.ion_fraction[20] == pytest.approx(0.00264536882147342)
        assert ref.ion_counts[0] == pytest.approx(0.081412829827488)
        assert ref.source_path == _REAL_CSV.resolve()

    def test_ion_fraction_sums_to_one(self):
        ref = load_he_abundance_reference(_REAL_CSV)
        assert ref.ion_fraction.sum() == pytest.approx(1.0)

    def test_n_is_integer_dtype(self):
        ref = load_he_abundance_reference(_REAL_CSV)
        assert ref.n.dtype.kind == "i"


# ---------------------------------------------------------------------------
# Synthetic happy path + order independence
# ---------------------------------------------------------------------------
class TestSyntheticHappyPath:
    def test_loads_synthetic(self, tmp_path):
        ref = load_he_abundance_reference(_write_csv(tmp_path / "a.csv", _GOOD_ROWS))
        np.testing.assert_array_equal(ref.n, np.arange(3))
        np.testing.assert_allclose(ref.ion_fraction, [0.5, 0.3, 0.2])
        assert ref.ion_fraction.sum() == pytest.approx(1.0)

    def test_column_order_independent(self, tmp_path):
        # Reorder columns: membership match must still succeed and map by name.
        header = (
            "ionPercent,n,ionCounts,label,massWindowLower_u_per_e,"
            "massCenter_u_per_e,massWindowUpper_u_per_e"
        )
        rows = [
            ("50.0", "0", "0.5", "I^+", "126", "127", "128"),
            ("30.0", "1", "0.3", "I^+He_1", "130.0026", "131.0026", "132.0026"),
            ("20.0", "2", "0.2", "I^+He_2", "134.0052", "135.0052", "136.0052"),
        ]
        ref = load_he_abundance_reference(_write_csv(tmp_path / "r.csv", rows, header=header))
        np.testing.assert_array_equal(ref.n, np.arange(3))
        np.testing.assert_allclose(ref.ion_fraction, [0.5, 0.3, 0.2])
        assert ref.label[1] == "I^+He_1"


# ---------------------------------------------------------------------------
# Fail-loud validation arms
# ---------------------------------------------------------------------------
class TestValidation:
    def test_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_he_abundance_reference(tmp_path / "nope.csv")

    def test_missing_column(self, tmp_path):
        header = (
            "n,label,massCenter_u_per_e,massWindowLower_u_per_e,"
            "massWindowUpper_u_per_e,ionCounts"  # ionPercent dropped
        )
        rows = [(r[0], r[1], r[2], r[3], r[4], r[5]) for r in _GOOD_ROWS]
        with pytest.raises(ValueError, match="column"):
            load_he_abundance_reference(_write_csv(tmp_path / "m.csv", rows, header=header))

    def test_extra_column(self, tmp_path):
        header = _HEADER + ",foo"
        rows = [r + ("1.0",) for r in _GOOD_ROWS]
        with pytest.raises(ValueError, match="column"):
            load_he_abundance_reference(_write_csv(tmp_path / "e.csv", rows, header=header))

    def test_non_contiguous_n(self, tmp_path):
        rows = [
            ("0", "I^+", "127", "126", "128", "0.5", "50.0"),
            ("1", "I^+He_1", "131.0026", "130.0026", "132.0026", "0.3", "30.0"),
            ("3", "I^+He_3", "139.0078", "138.0078", "140.0078", "0.2", "20.0"),
        ]
        with pytest.raises(ValueError, match="contiguous|ascending"):
            load_he_abundance_reference(_write_csv(tmp_path / "nc.csv", rows))

    def test_n_not_from_zero(self, tmp_path):
        rows = [
            ("1", "I^+He_1", "131.0026", "130.0026", "132.0026", "0.5", "50.0"),
            ("2", "I^+He_2", "135.0052", "134.0052", "136.0052", "0.3", "30.0"),
            ("3", "I^+He_3", "139.0078", "138.0078", "140.0078", "0.2", "20.0"),
        ]
        with pytest.raises(ValueError, match="contiguous|ascending|zero"):
            load_he_abundance_reference(_write_csv(tmp_path / "nz.csv", rows))

    def test_negative_counts(self, tmp_path):
        rows = [
            ("0", "I^+", "127", "126", "128", "-0.1", "50.0"),
            ("1", "I^+He_1", "131.0026", "130.0026", "132.0026", "0.3", "30.0"),
            ("2", "I^+He_2", "135.0052", "134.0052", "136.0052", "0.2", "20.0"),
        ]
        with pytest.raises(ValueError, match="ionCounts|negative|>= 0|>=0"):
            load_he_abundance_reference(_write_csv(tmp_path / "neg.csv", rows))

    def test_blank_percent_cell_rejected(self, tmp_path):
        # np.genfromtxt fills a blank numeric cell with NaN, and NaN
        # comparisons are False -- a NaN-blind sum guard would load this
        # silently and hand E4 an all-NaN ion_fraction.
        rows = [
            ("0", "I^+", "127", "126", "128", "0.5", "50.0"),
            ("1", "I^+He_1", "131.0026", "130.0026", "132.0026", "0.3", ""),
            ("2", "I^+He_2", "135.0052", "134.0052", "136.0052", "0.2", "50.0"),
        ]
        with pytest.raises(ValueError, match="finite|NaN"):
            load_he_abundance_reference(_write_csv(tmp_path / "bpc.csv", rows))

    def test_blank_counts_cell_rejected(self, tmp_path):
        rows = [
            ("0", "I^+", "127", "126", "128", "", "50.0"),
            ("1", "I^+He_1", "131.0026", "130.0026", "132.0026", "0.3", "30.0"),
            ("2", "I^+He_2", "135.0052", "134.0052", "136.0052", "0.2", "20.0"),
        ]
        with pytest.raises(ValueError, match="finite|NaN"):
            load_he_abundance_reference(_write_csv(tmp_path / "bcc.csv", rows))

    def test_blank_mass_center_cell_rejected(self, tmp_path):
        rows = [
            ("0", "I^+", "", "126", "128", "0.5", "50.0"),
            ("1", "I^+He_1", "131.0026", "130.0026", "132.0026", "0.3", "30.0"),
            ("2", "I^+He_2", "135.0052", "134.0052", "136.0052", "0.2", "20.0"),
        ]
        with pytest.raises(ValueError, match="finite|NaN"):
            load_he_abundance_reference(_write_csv(tmp_path / "bmc.csv", rows))

    def test_negative_percent_rejected(self, tmp_path):
        # Sums to exactly 100, so only a sign guard catches the corruption.
        rows = [
            ("0", "I^+", "127", "126", "128", "0.5", "110.0"),
            ("1", "I^+He_1", "131.0026", "130.0026", "132.0026", "0.3", "-5.0"),
            ("2", "I^+He_2", "135.0052", "134.0052", "136.0052", "0.2", "-5.0"),
        ]
        with pytest.raises(ValueError, match="ionPercent.*negative|negative.*ionPercent"):
            load_he_abundance_reference(_write_csv(tmp_path / "np.csv", rows))

    def test_bad_percent_sum(self, tmp_path):
        rows = [
            ("0", "I^+", "127", "126", "128", "0.5", "50.0"),
            ("1", "I^+He_1", "131.0026", "130.0026", "132.0026", "0.3", "30.0"),
            ("2", "I^+He_2", "135.0052", "134.0052", "136.0052", "0.2", "10.0"),  # sum 90
        ]
        with pytest.raises(ValueError, match="100|sum"):
            load_he_abundance_reference(_write_csv(tmp_path / "bp.csv", rows))


class TestFrozen:
    def test_frozen(self):
        ref = load_he_abundance_reference(_REAL_CSV)
        with pytest.raises(Exception):
            ref.n = np.arange(5)  # frozen dataclass
