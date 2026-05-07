"""Tests for i2_helium_md/postprocess/hedft_loader.py."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.postprocess.hedft_loader import (
    HedftTrajectory,
    load_hedft_trajectory,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REF_9A = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"
REF_18A = PROJECT_ROOT / "data" / "reference" / "18A_All_Data.csv"


# ===========================================================================
# Helpers
# ===========================================================================
_HEADER = "Time_ps,V1_mag,V2_mag,V1_z,V2_z,V1_x,V2_x,R_distance"


def _write_synthetic_csv(path: Path, *, t: np.ndarray, r: np.ndarray) -> None:
    """Write a tiny CSV with the 8 expected columns (zeros except t and r)."""
    n = t.size
    z = np.zeros(n)
    rows = [_HEADER]
    for i in range(n):
        rows.append(
            f"{t[i]},{z[i]},{z[i]},{z[i]},{z[i]},{z[i]},{z[i]},{r[i]}"
        )
    path.write_text("\n".join(rows) + "\n", encoding="ascii")


# ===========================================================================
# Synthetic round-trip tests
# ===========================================================================
class TestSyntheticRoundTrip:
    def test_load_round_trip_synthetic(self, tmp_path):
        """Tiny CSV: every column must round-trip to the dataclass."""
        t = np.array([0.0, 0.1, 0.2, 0.3])
        r = np.array([9.0, 9.1, 9.4, 10.0])
        path = tmp_path / "9A_All_Data.csv"
        _write_synthetic_csv(path, t=t, r=r)

        traj = load_hedft_trajectory(path)

        assert isinstance(traj, HedftTrajectory)
        np.testing.assert_array_equal(traj.time_ps, t)
        np.testing.assert_array_equal(traj.distance_A, r)
        np.testing.assert_array_equal(traj.v1_magnitude_Aps, np.zeros(4))
        np.testing.assert_array_equal(traj.v2_magnitude_Aps, np.zeros(4))
        assert traj.droplet_radius_A == 9.0
        assert traj.source_path == path.resolve()

    def test_droplet_radius_inferred_from_18A_prefix(self, tmp_path):
        path = tmp_path / "18A_All_Data.csv"
        _write_synthetic_csv(
            path, t=np.array([0.0, 0.1]), r=np.array([18.0, 18.0])
        )
        traj = load_hedft_trajectory(path)
        assert traj.droplet_radius_A == 18.0

    def test_droplet_radius_override_wins(self, tmp_path):
        """Explicit argument overrides filename inference."""
        path = tmp_path / "9A_All_Data.csv"
        _write_synthetic_csv(
            path, t=np.array([0.0, 0.1]), r=np.array([9.0, 9.0])
        )
        traj = load_hedft_trajectory(path, droplet_radius_A=12.0)
        assert traj.droplet_radius_A == 12.0

    def test_droplet_radius_inference_fails_raises(self, tmp_path):
        path = tmp_path / "no_prefix.csv"
        _write_synthetic_csv(
            path, t=np.array([0.0, 0.1]), r=np.array([9.0, 9.0])
        )
        with pytest.raises(ValueError, match="Cannot infer droplet radius"):
            load_hedft_trajectory(path)

    def test_droplet_radius_inference_fallback_to_arg(self, tmp_path):
        path = tmp_path / "no_prefix.csv"
        _write_synthetic_csv(
            path, t=np.array([0.0, 0.1]), r=np.array([9.0, 9.0])
        )
        traj = load_hedft_trajectory(path, droplet_radius_A=18.0)
        assert traj.droplet_radius_A == 18.0


# ===========================================================================
# Validation tests
# ===========================================================================
class TestValidation:
    def test_missing_file_raises(self, tmp_path):
        bogus = tmp_path / "does_not_exist.csv"
        with pytest.raises(FileNotFoundError, match="HeDFT reference file"):
            load_hedft_trajectory(bogus)

    def test_bad_header_missing_column_raises(self, tmp_path):
        path = tmp_path / "9A_All_Data.csv"
        # Drop the R_distance column
        path.write_text(
            "Time_ps,V1_mag,V2_mag,V1_z,V2_z,V1_x,V2_x\n"
            "0.0,0,0,0,0,0,0\n"
            "0.1,0,0,0,0,0,0\n",
            encoding="ascii",
        )
        with pytest.raises(ValueError, match="missing=.*R_distance"):
            load_hedft_trajectory(path)

    def test_bad_header_extra_column_raises(self, tmp_path):
        path = tmp_path / "9A_All_Data.csv"
        path.write_text(
            "Time_ps,V1_mag,V2_mag,V1_z,V2_z,V1_x,V2_x,R_distance,Extra\n"
            "0.0,0,0,0,0,0,0,9,0\n"
            "0.1,0,0,0,0,0,0,9,0\n",
            encoding="ascii",
        )
        with pytest.raises(ValueError, match="unexpected=.*Extra"):
            load_hedft_trajectory(path)

    def test_non_monotonic_time_raises(self, tmp_path):
        path = tmp_path / "9A_All_Data.csv"
        # Shuffle the time column so it is not strictly increasing.
        _write_synthetic_csv(
            path,
            t=np.array([0.0, 0.2, 0.1, 0.3]),
            r=np.array([9.0, 9.1, 9.2, 9.3]),
        )
        with pytest.raises(ValueError, match="non-monotonic Time_ps"):
            load_hedft_trajectory(path)

    def test_too_few_samples_raises(self, tmp_path):
        path = tmp_path / "9A_All_Data.csv"
        path.write_text(
            _HEADER + "\n0.0,0,0,0,0,0,0,9.0\n",
            encoding="ascii",
        )
        with pytest.raises(ValueError, match="at least 2 time samples"):
            load_hedft_trajectory(path)


# ===========================================================================
# Real reference data
# ===========================================================================
@pytest.mark.skipif(
    not REF_9A.exists(), reason="data/reference/9A_All_Data.csv not present"
)
class TestReal9A:
    def test_load_real_9A_reference_shape(self):
        traj = load_hedft_trajectory(REF_9A)
        assert traj.droplet_radius_A == 9.0
        assert traj.time_ps.shape == traj.distance_A.shape
        assert traj.time_ps.shape[0] == 14082
        assert traj.time_ps[0] == 0.0
        assert traj.distance_A[0] == pytest.approx(9.0)
        # All eight series share the time grid length.
        for arr in (
            traj.v1_magnitude_Aps,
            traj.v2_magnitude_Aps,
            traj.v1_z_Aps,
            traj.v2_z_Aps,
            traj.v1_x_Aps,
            traj.v2_x_Aps,
        ):
            assert arr.shape == traj.time_ps.shape


@pytest.mark.skipif(
    not REF_18A.exists(), reason="data/reference/18A_All_Data.csv not present"
)
class TestReal18A:
    def test_load_real_18A_reference_shape(self):
        traj = load_hedft_trajectory(REF_18A)
        assert traj.droplet_radius_A == 18.0
        assert traj.time_ps.shape == traj.distance_A.shape
        assert traj.time_ps.shape[0] == 14769
        assert traj.time_ps[0] == 0.0
        assert traj.distance_A[0] == pytest.approx(18.0)
