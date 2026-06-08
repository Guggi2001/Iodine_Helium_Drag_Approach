"""Tests for i2_helium_md/postprocess/hedft_loader.py."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.postprocess.hedft_loader import (
    HedftTrajectory,
    SmoothedSpeedReference,
    load_hedft_trajectory,
    load_smoothed_speed_reference,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REF_9A = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"
REF_18A = PROJECT_ROOT / "data" / "reference" / "18A_All_Data.csv"
CLEANED_9A = (
    PROJECT_ROOT / "data" / "reference" / "drag" / "9A"
    / "velocity_smoothed" / "cleaned_data.csv"
)


# ===========================================================================
# Helpers
# ===========================================================================
# Full 16-column reference header: velocities (mags + 3D components),
# separation, and per-atom 3D positions.
_HEADER = (
    "Time_ps,V1_mag,V2_mag,V1_x,V1_y,V1_z,V2_x,V2_y,V2_z,R_distance,"
    "X1,Y1,Z1,X2,Y2,Z2"
)


def _write_synthetic_csv(path: Path, *, t: np.ndarray, r: np.ndarray) -> None:
    """Write a tiny CSV with the 16 expected columns (zeros except t and r).

    Time_ps and R_distance carry the supplied arrays; every velocity component
    and position column is zero. Sufficient for loader round-trip / validation
    tests, which do not exercise the velocity or position values.
    """
    n = t.size
    z = np.zeros(n)

    def _row(i: int) -> str:
        # Time_ps, V1_mag, V2_mag, V1_x, V1_y, V1_z, V2_x, V2_y, V2_z,
        # R_distance, X1, Y1, Z1, X2, Y2, Z2
        vals = [t[i], z[i], z[i], z[i], z[i], z[i], z[i], z[i], z[i],
                r[i], z[i], z[i], z[i], z[i], z[i], z[i]]
        return ",".join(str(v) for v in vals)

    rows = [_HEADER] + [_row(i) for i in range(n)]
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
        # 3D velocity components and per-atom positions round-trip too.
        for arr in (
            traj.v1_x_Aps, traj.v1_y_Aps, traj.v2_x_Aps, traj.v2_y_Aps,
            traj.x1_A, traj.y1_A, traj.z1_A, traj.x2_A, traj.y2_A, traj.z2_A,
        ):
            np.testing.assert_array_equal(arr, np.zeros(4))
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
        # Drop the R_distance column from the otherwise-complete header.
        header_no_r = (
            "Time_ps,V1_mag,V2_mag,V1_x,V1_y,V1_z,V2_x,V2_y,V2_z,"
            "X1,Y1,Z1,X2,Y2,Z2"
        )  # 15 columns: the full 16 minus R_distance
        rest = ",".join(["0"] * 14)
        path.write_text(
            f"{header_no_r}\n0.0,{rest}\n0.1,{rest}\n",
            encoding="ascii",
        )
        with pytest.raises(ValueError, match="missing=.*R_distance"):
            load_hedft_trajectory(path)

    def test_bad_header_extra_column_raises(self, tmp_path):
        path = tmp_path / "9A_All_Data.csv"
        row = ",".join(["0"] * 16) + ",0"  # 16 expected + 1 extra
        path.write_text(
            f"{_HEADER},Extra\n{row}\n{row}\n",
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
        # One data row: R_distance is column index 9, the rest zero.
        row = "0.0," + ",".join(["0"] * 8) + ",9.0," + ",".join(["0"] * 6)
        path.write_text(
            _HEADER + "\n" + row + "\n",
            encoding="ascii",
        )
        with pytest.raises(ValueError, match="at least 2 time samples"):
            load_hedft_trajectory(path)


# ===========================================================================
# Smoothed-speed reference (Tier-0 same-smoothed comparison)
# ===========================================================================
_SMOOTHED_HEADER = "time,cleaned_SG"


def _write_smoothed_csv(path: Path, *, t: np.ndarray, v: np.ndarray) -> None:
    rows = [_SMOOTHED_HEADER]
    for ti, vi in zip(t, v):
        rows.append(f"{ti},{vi}")
    path.write_text("\n".join(rows) + "\n", encoding="ascii")


class TestSmoothedSpeedReference:
    def test_load_round_trip(self, tmp_path):
        t = np.array([2.67, 2.671, 2.672, 2.673])
        v = np.array([4.93, 4.92, 4.91, 4.90])
        path = tmp_path / "cleaned_data.csv"
        _write_smoothed_csv(path, t=t, v=v)

        ref = load_smoothed_speed_reference(path)

        assert isinstance(ref, SmoothedSpeedReference)
        np.testing.assert_array_equal(ref.time_ps, t)
        np.testing.assert_array_equal(ref.speed_Aps, v)
        assert ref.source_path == path.resolve()

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Smoothed-speed reference"):
            load_smoothed_speed_reference(tmp_path / "nope.csv")

    def test_bad_header_raises(self, tmp_path):
        path = tmp_path / "cleaned_data.csv"
        path.write_text("t,v\n2.67,4.9\n2.68,4.8\n", encoding="ascii")
        with pytest.raises(ValueError, match="unexpected columns"):
            load_smoothed_speed_reference(path)

    def test_non_monotonic_time_raises(self, tmp_path):
        path = tmp_path / "cleaned_data.csv"
        _write_smoothed_csv(
            path,
            t=np.array([2.67, 2.69, 2.68, 2.70]),
            v=np.array([4.9, 4.8, 4.7, 4.6]),
        )
        with pytest.raises(ValueError, match="non-monotonic time"):
            load_smoothed_speed_reference(path)

    def test_too_few_samples_raises(self, tmp_path):
        path = tmp_path / "cleaned_data.csv"
        path.write_text(_SMOOTHED_HEADER + "\n2.67,4.9\n", encoding="ascii")
        with pytest.raises(ValueError, match="at least 2 samples"):
            load_smoothed_speed_reference(path)


@pytest.mark.skipif(
    not CLEANED_9A.exists(),
    reason="cleaned 9A velocity reference not present",
)
class TestRealCleaned9A:
    def test_load_real_cleaned_9A(self):
        ref = load_smoothed_speed_reference(CLEANED_9A)
        assert ref.time_ps[0] == pytest.approx(2.67)
        assert ref.time_ps[-1] == pytest.approx(6.0)
        assert ref.time_ps.shape == ref.speed_Aps.shape
        assert np.all(np.diff(ref.time_ps) > 0.0)
        assert np.all(ref.speed_Aps > 0.0)


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
        # All component/position series share the time grid length.
        for arr in (
            traj.v1_magnitude_Aps,
            traj.v2_magnitude_Aps,
            traj.v1_z_Aps,
            traj.v2_z_Aps,
            traj.v1_x_Aps,
            traj.v2_x_Aps,
            traj.v1_y_Aps,
            traj.v2_y_Aps,
            traj.x1_A, traj.y1_A, traj.z1_A,
            traj.x2_A, traj.y2_A, traj.z2_A,
        ):
            assert arr.shape == traj.time_ps.shape
        # Per-atom positions are consistent with the stored separation.
        sep = np.sqrt(
            (traj.x1_A - traj.x2_A) ** 2
            + (traj.y1_A - traj.y2_A) ** 2
            + (traj.z1_A - traj.z2_A) ** 2
        )
        np.testing.assert_allclose(sep, traj.distance_A, rtol=1e-6, atol=1e-4)


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
