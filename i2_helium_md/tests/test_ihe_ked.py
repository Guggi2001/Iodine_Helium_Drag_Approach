"""Tests for i2_helium_md/postprocess/ihe_ked.py.

The loader tests run against the real frozen reference under
``data/reference/ihe_ked/`` (small CSVs, committed) plus synthetic
corrupted files in tmp_path for the failure paths.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.postprocess.ihe_ked import (
    IHeKedReference,
    load_ihe_ked_reference,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
IHE_KED_DIR = PROJECT_ROOT / "data" / "reference" / "ihe_ked"
REFERENCE_CSV = IHE_KED_DIR / "IHe_KED_reference.csv"


class TestLoadIHeKedReference:
    def test_real_file_contract(self):
        ref = load_ihe_ked_reference(REFERENCE_CSV)
        assert isinstance(ref, IHeKedReference)
        np.testing.assert_array_equal(ref.n, np.arange(18))
        # Frozen reference value (README: I+ scale anchor 3.706 eV).
        assert ref.mean_KE_eV[0] == pytest.approx(3.70569, abs=1e-5)
        # Mass model column matches m(n) = 126.90 + 4.0026 n.
        np.testing.assert_allclose(
            ref.mass_center_u, 126.90 + 4.0026 * ref.n, atol=1e-4,
        )
        # noiseLimited is 0 for all 18 fragments (README).
        assert not ref.noise_limited.any()
        assert ref.source_path == REFERENCE_CSV.resolve()

    def test_gold_mask_is_calib_limited_set(self):
        ref = load_ihe_ked_reference(REFERENCE_CSV)
        expected_gold = {0, 3, 4, 5, 6, 7, 8, 9, 10, 12}
        assert set(ref.n[ref.gold_mask].tolist()) == expected_gold

    def test_point_err_is_stat_sys_quadrature(self):
        ref = load_ihe_ked_reference(REFERENCE_CSV)
        np.testing.assert_allclose(
            ref.point_err_eV,
            np.sqrt(
                ref.stat_err_mean_KE_eV ** 2 + ref.sys_err_mean_KE_eV ** 2
            ),
            rtol=1e-12,  # tight: pure arithmetic identity
        )

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_ihe_ked_reference(tmp_path / "nope.csv")

    def test_missing_column_raises(self, tmp_path):
        p = tmp_path / "bad.csv"
        p.write_text("n,label\n0,I+He_0\n", encoding="ascii")
        with pytest.raises(ValueError):
            load_ihe_ked_reference(p)

    def test_non_contiguous_n_raises(self, tmp_path):
        header = (
            "n,label,massCenter_u_per_e,N_counts,N_eff,meanKE_eV,"
            "modeKE_eV,medianKE_eV,sigmaKE_eV,statErr_meanKE_eV,"
            "sysErr_meanKE_eV,calibSyst_frac,conditionSyst_frac,"
            "bgOffShift_eV,dominantError,noiseLimited\n"
        )
        row = "5,I+He_5,146.9,100,90,1.0,0.9,0.95,0.5,0.01,0.02,0.04,0.06,0.0,calib,0\n"
        p = tmp_path / "bad.csv"
        p.write_text(header + row, encoding="ascii")
        with pytest.raises(ValueError, match="contiguous"):
            load_ihe_ked_reference(p)

    def test_unknown_dominant_error_raises(self, tmp_path):
        header = (
            "n,label,massCenter_u_per_e,N_counts,N_eff,meanKE_eV,"
            "modeKE_eV,medianKE_eV,sigmaKE_eV,statErr_meanKE_eV,"
            "sysErr_meanKE_eV,calibSyst_frac,conditionSyst_frac,"
            "bgOffShift_eV,dominantError,noiseLimited\n"
        )
        row = "0,I+He_0,126.9,100,90,1.0,0.9,0.95,0.5,0.01,0.02,0.04,0.06,0.0,vibes,0\n"
        p = tmp_path / "bad.csv"
        p.write_text(header + row, encoding="ascii")
        with pytest.raises(ValueError, match="dominantError"):
            load_ihe_ked_reference(p)
