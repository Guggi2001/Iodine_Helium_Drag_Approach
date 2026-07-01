"""Tests for i2_helium_md/postprocess/compare_trajectories.py."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.simulation.checkpoint import IonCheckpoint
from i2_helium_md.postprocess.compare_trajectories import (
    TrajectoryComparison,
    compare_distance,
    compare_speed_to_reference,
    compare_velocity_magnitude,
)
from i2_helium_md.postprocess.hedft_loader import HedftTrajectory


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "9A_hedft_comparison"
ION_NPZ = RUN_DIR / "ion.npz"
REF_9A = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"


# ===========================================================================
# Helpers
# ===========================================================================
def _make_ion_with_distance(
    *,
    t_md: np.ndarray,
    distance: np.ndarray,
    num_molecules: int = 2,
    speed_i1: float = 0.3,
    speed_i2: float = 0.7,
) -> IonCheckpoint:
    """Build a minimal IonCheckpoint whose mean I-I distance matches ``distance``.

    All molecules are placed identically along x: I1 at +R/2, I2 at -R/2,
    with y=z=0. Each atom carries a constant velocity in ``+x`` so that the
    velocity-magnitude comparison has a known answer per atom group.
    """
    n = num_molecules
    t_md = np.asarray(t_md, dtype=float)
    distance = np.asarray(distance, dtype=float)
    if t_md.shape != distance.shape:
        raise AssertionError("test helper: t_md and distance must align")

    t_steps = t_md.size
    half = distance / 2.0  # (T,)

    pos_x = np.zeros((2 * n, t_steps))
    pos_x[:n] = +half[None, :]   # atom I1 indices [0, n)
    pos_x[n:] = -half[None, :]   # atom I2 indices [n, 2n)

    pos_y = np.zeros_like(pos_x)
    pos_z = np.zeros_like(pos_x)

    vel_x = np.zeros_like(pos_x)
    vel_x[:n] = speed_i1
    vel_x[n:] = speed_i2
    vel_y = np.zeros_like(pos_x)
    vel_z = np.zeros_like(pos_x)

    diag_zero = np.zeros((2 * n, t_steps))
    return IonCheckpoint(
        num_molecules=n,
        time_ps=t_md,
        positions_x=pos_x,
        positions_y=pos_y,
        positions_z=pos_z,
        velocities_x=vel_x,
        velocities_y=vel_y,
        velocities_z=vel_z,
        positions_final_x=np.zeros(2 * n),
        positions_final_y=np.zeros(2 * n),
        positions_final_z=np.zeros(2 * n),
        velocities_final_x=np.zeros(2 * n),
        velocities_final_y=np.zeros(2 * n),
        velocities_final_z=np.zeros(2 * n),
        mass_kg=np.full(2 * n, 127 * 1.66e-27),
        mass_final_kg=np.full(2 * n, 127 * 1.66e-27),
        mass_history_kg=np.full((2 * n, t_steps), 127 * 1.66e-27),
        droplet_radii_angstrom=np.full(2 * n, 30.0),
        E_kin_eV=diag_zero,
        E_pot_eV=diag_zero,
        E_dissip_eV=diag_zero,
        E_mass_transfer_eV=diag_zero,
        E_int_eV=diag_zero,
        n_shell=diag_zero,
        b_ion_outside=np.zeros(n, dtype=bool),
        relative_loss_per_ps=diag_zero,
        number_of_collisions=np.zeros((2 * n, t_steps), dtype=int),
        temperature_diagnostic=np.full((t_steps, 3), np.nan, dtype=float),
    )


def _make_hedft(
    *,
    t_ps: np.ndarray,
    distance: np.ndarray | None = None,
    v1_mag: np.ndarray | None = None,
    v2_mag: np.ndarray | None = None,
) -> HedftTrajectory:
    """Build a HedftTrajectory directly (no CSV roundtrip)."""
    t = np.asarray(t_ps, dtype=float)
    z = np.zeros_like(t)
    return HedftTrajectory(
        time_ps=t,
        v1_magnitude_Aps=z if v1_mag is None else np.asarray(v1_mag, float),
        v2_magnitude_Aps=z if v2_mag is None else np.asarray(v2_mag, float),
        v1_z_Aps=z,
        v2_z_Aps=z,
        v1_x_Aps=z,
        v2_x_Aps=z,
        v1_y_Aps=z,
        v2_y_Aps=z,
        distance_A=z if distance is None else np.asarray(distance, float),
        x1_A=z, x2_A=z, y1_A=z, y2_A=z, z1_A=z, z2_A=z,
        droplet_radius_A=9.0,
        source_path=Path("synthetic"),
    )


# ===========================================================================
# Distance comparison
# ===========================================================================
class TestCompareDistance:
    def test_zero_when_md_matches_hedft(self):
        """If MD distance == HeDFT distance on the same grid, RMSE is 0."""
        t = np.linspace(0.0, 1.0, 11)
        r = 9.0 + 0.5 * t
        ion = _make_ion_with_distance(t_md=t, distance=r, num_molecules=3)
        hedft = _make_hedft(t_ps=t, distance=r)

        result = compare_distance(ion, hedft)

        assert isinstance(result, TrajectoryComparison)
        assert result.quantity == "distance_A"
        assert result.rmse == pytest.approx(0.0, abs=1e-12)
        assert result.mean_ratio == pytest.approx(1.0, abs=1e-12)
        assert result.num_overlap_points == 11
        assert result.overlap_t_min_ps == 0.0
        assert result.overlap_t_max_ps == 1.0

    def test_known_constant_offset(self):
        """A constant +0.5 angstrom offset must give rmse == 0.5 exactly."""
        t = np.linspace(0.0, 1.0, 21)
        r_hedft = 9.0 + 0.3 * t
        r_md = r_hedft + 0.5
        ion = _make_ion_with_distance(t_md=t, distance=r_md)
        hedft = _make_hedft(t_ps=t, distance=r_hedft)

        result = compare_distance(ion, hedft)

        assert result.rmse == pytest.approx(0.5, abs=1e-12)
        assert result.mean_ratio > 1.0

    def test_resampling_recovers_analytic_value(self):
        """HeDFT grid coarser than MD grid: linear interp should be exact
        for a linear underlying signal."""
        t_md = np.linspace(0.0, 1.0, 101)
        t_hedft = np.linspace(0.0, 1.0, 11)
        r_func = lambda x: 9.0 + 2.0 * x  # noqa: E731

        ion = _make_ion_with_distance(t_md=t_md, distance=r_func(t_md))
        hedft = _make_hedft(t_ps=t_hedft, distance=r_func(t_hedft))

        result = compare_distance(ion, hedft)

        assert result.num_overlap_points == 11
        np.testing.assert_allclose(
            result.md_on_hedft_grid, r_func(t_hedft), atol=1e-12
        )
        assert result.rmse == pytest.approx(0.0, abs=1e-12)

    def test_partial_overlap_truncates_window(self):
        """HeDFT goes to 5 ps but MD only to 3 ps -> overlap caps at 3 ps."""
        t_md = np.linspace(0.0, 3.0, 31)
        t_hedft = np.linspace(0.0, 5.0, 51)
        r_md = 9.0 + 0.0 * t_md
        r_hedft = 9.0 + 0.0 * t_hedft
        ion = _make_ion_with_distance(t_md=t_md, distance=r_md)
        hedft = _make_hedft(t_ps=t_hedft, distance=r_hedft)

        result = compare_distance(ion, hedft)

        assert result.overlap_t_min_ps == 0.0
        assert result.overlap_t_max_ps == 3.0
        # HeDFT samples at t in [0, 3] inclusive: 31 of them.
        assert result.num_overlap_points == 31
        assert result.t_overlap_ps[-1] == pytest.approx(3.0)

    def test_no_overlap_raises(self):
        t_md = np.linspace(0.0, 1.0, 11)
        t_hedft = np.linspace(2.0, 3.0, 11)
        r = np.full_like(t_md, 9.0)
        ion = _make_ion_with_distance(t_md=t_md, distance=r)
        hedft = _make_hedft(t_ps=t_hedft, distance=np.full_like(t_hedft, 9.0))

        with pytest.raises(ValueError, match="do not overlap"):
            compare_distance(ion, hedft)

    def test_single_overlap_sample_raises(self):
        """The overlap must contain at least 2 samples for a meaningful interp."""
        t_md = np.linspace(0.0, 1.0, 11)
        t_hedft = np.array([1.0, 2.0, 3.0])  # only t=1.0 lies in MD range
        ion = _make_ion_with_distance(
            t_md=t_md, distance=np.full_like(t_md, 9.0)
        )
        hedft = _make_hedft(
            t_ps=t_hedft, distance=np.full_like(t_hedft, 9.0)
        )
        with pytest.raises(ValueError, match="do not overlap"):
            compare_distance(ion, hedft)


# ===========================================================================
# Velocity comparison
# ===========================================================================
class TestCompareVelocityMagnitude:
    def test_atom_selection_picks_correct_indices(self):
        """I1 atoms have speed 0.3; I2 atoms have speed 0.7. Verify both pick correctly."""
        t = np.linspace(0.0, 1.0, 11)
        ion = _make_ion_with_distance(
            t_md=t,
            distance=np.full_like(t, 9.0),
            num_molecules=4,
            speed_i1=0.3,
            speed_i2=0.7,
        )
        hedft = _make_hedft(
            t_ps=t,
            v1_mag=np.full_like(t, 0.3),
            v2_mag=np.full_like(t, 0.7),
        )

        r1 = compare_velocity_magnitude(ion, hedft, atom="I1")
        r2 = compare_velocity_magnitude(ion, hedft, atom="I2")

        assert r1.quantity == "v1_magnitude_Aps"
        assert r2.quantity == "v2_magnitude_Aps"
        assert r1.rmse == pytest.approx(0.0, abs=1e-12)
        assert r2.rmse == pytest.approx(0.0, abs=1e-12)
        assert r1.mean_ratio == pytest.approx(1.0, abs=1e-12)
        assert r2.mean_ratio == pytest.approx(1.0, abs=1e-12)

    def test_invalid_atom_raises(self):
        t = np.linspace(0.0, 1.0, 11)
        ion = _make_ion_with_distance(t_md=t, distance=np.full_like(t, 9.0))
        hedft = _make_hedft(t_ps=t)
        with pytest.raises(ValueError, match="atom must be"):
            compare_velocity_magnitude(ion, hedft, atom="I3")  # type: ignore[arg-type]

    def test_zero_reference_denominator_excluded_from_ratio(self):
        """HeDFT velocity is exactly 0 at the first sample (real data
        starts at v=0). Mean ratio must skip that sample, not return inf."""
        t = np.linspace(0.0, 1.0, 11)
        ion = _make_ion_with_distance(
            t_md=t,
            distance=np.full_like(t, 9.0),
            speed_i1=0.5,
            speed_i2=0.5,
        )
        v_ref = np.full_like(t, 0.5)
        v_ref[0] = 0.0
        hedft = _make_hedft(t_ps=t, v1_mag=v_ref, v2_mag=v_ref)

        result = compare_velocity_magnitude(ion, hedft, atom="I1")

        assert np.isfinite(result.mean_ratio)
        # The remaining 10 samples have md == ref == 0.5, so ratio == 1.
        assert result.mean_ratio == pytest.approx(1.0, abs=1e-12)
        assert result.num_overlap_points == 11  # raw overlap unchanged


# ===========================================================================
# Windowing (Tier-0 in-window scoring)
# ===========================================================================
class TestWindow:
    def test_none_is_bit_identical_to_default(self):
        """window=None must reproduce the whole-overlap path exactly."""
        t = np.linspace(0.0, 10.0, 101)
        r_hedft = 9.0 + 0.5 * t
        r_md = r_hedft + 0.2 * np.sin(t)
        ion = _make_ion_with_distance(t_md=t, distance=r_md)
        hedft = _make_hedft(t_ps=t, distance=r_hedft)

        base = compare_distance(ion, hedft)
        explicit_none = compare_distance(ion, hedft, window=None)

        assert explicit_none.rmse == base.rmse
        assert explicit_none.mean_ratio == base.mean_ratio
        assert explicit_none.num_overlap_points == base.num_overlap_points
        assert explicit_none.overlap_t_min_ps == base.overlap_t_min_ps
        assert explicit_none.overlap_t_max_ps == base.overlap_t_max_ps

    def test_window_restricts_scored_samples(self):
        """A window must drop out-of-window samples from the scored RMSE.

        Construct a distance series that matches the reference exactly inside
        [3, 7] but is offset by +5 angstrom outside it. The whole-overlap RMSE
        is dominated by the offset region; the windowed RMSE is ~0.
        """
        t = np.linspace(0.0, 10.0, 101)
        r_hedft = np.full_like(t, 9.0)
        r_md = r_hedft.copy()
        outside = (t < 3.0) | (t > 7.0)
        r_md[outside] += 5.0
        ion = _make_ion_with_distance(t_md=t, distance=r_md)
        hedft = _make_hedft(t_ps=t, distance=r_hedft)

        full = compare_distance(ion, hedft)
        windowed = compare_distance(ion, hedft, window=(3.0, 7.0))

        assert full.rmse > 1.0  # dominated by the out-of-window offset
        assert windowed.rmse == pytest.approx(0.0, abs=1e-12)
        # The scored bounds collapse to the window.
        assert windowed.overlap_t_min_ps == pytest.approx(3.0)
        assert windowed.overlap_t_max_ps == pytest.approx(7.0)
        assert windowed.num_overlap_points < full.num_overlap_points
        assert np.all(windowed.t_overlap_ps >= 3.0)
        assert np.all(windowed.t_overlap_ps <= 7.0)

    def test_window_applies_to_velocity(self):
        """Windowing flows through compare_velocity_magnitude identically."""
        t = np.linspace(0.0, 10.0, 101)
        v_ref = np.full_like(t, 0.5)
        ion = _make_ion_with_distance(
            t_md=t, distance=np.full_like(t, 9.0),
            speed_i1=0.5, speed_i2=0.5,
        )
        hedft = _make_hedft(t_ps=t, v1_mag=v_ref, v2_mag=v_ref)

        windowed = compare_velocity_magnitude(
            ion, hedft, atom="I1", window=(2.67, 8.5)
        )
        assert windowed.rmse == pytest.approx(0.0, abs=1e-12)
        assert windowed.overlap_t_min_ps == pytest.approx(2.67)
        assert windowed.overlap_t_max_ps == pytest.approx(8.5)

    def test_window_intersects_with_available_overlap(self):
        """A window wider than the data clamps to the available overlap."""
        t_md = np.linspace(0.0, 5.0, 51)
        t_hedft = np.linspace(0.0, 5.0, 51)
        r = np.full_like(t_md, 9.0)
        ion = _make_ion_with_distance(t_md=t_md, distance=r)
        hedft = _make_hedft(t_ps=t_hedft, distance=np.full_like(t_hedft, 9.0))

        windowed = compare_distance(ion, hedft, window=(-1.0, 100.0))
        full = compare_distance(ion, hedft)
        assert windowed.overlap_t_min_ps == full.overlap_t_min_ps
        assert windowed.overlap_t_max_ps == full.overlap_t_max_ps
        assert windowed.num_overlap_points == full.num_overlap_points

    def test_too_narrow_window_raises(self):
        """A window admitting <2 samples raises the existing overlap error."""
        t = np.linspace(0.0, 10.0, 11)  # samples at integer ps
        ion = _make_ion_with_distance(t_md=t, distance=np.full_like(t, 9.0))
        hedft = _make_hedft(t_ps=t, distance=np.full_like(t, 9.0))
        # Window [4.1, 4.9] contains no HeDFT sample.
        with pytest.raises(ValueError, match="do not overlap"):
            compare_distance(ion, hedft, window=(4.1, 4.9))


# ===========================================================================
# Speed-vs-arbitrary-reference (Tier-0 same-smoothed consistency comparison)
# ===========================================================================
class TestCompareSpeedToReference:
    """compare_speed_to_reference scores MD |v1|/|v2| against any speed curve.

    Used by the Tier-0 same-smoothed comparison to score MD |v2| against the
    CEEMDAN+SG-denoised |v2| reference (which is a 2-column CSV, not a full
    HedftTrajectory).
    """

    def test_exact_recovery_against_matching_series(self):
        t = np.linspace(0.0, 1.0, 11)
        ion = _make_ion_with_distance(
            t_md=t, distance=np.full_like(t, 9.0),
            speed_i1=0.3, speed_i2=0.7,
        )
        ref = np.full_like(t, 0.7)
        result = compare_speed_to_reference(
            ion, atom="I2", t_ref_ps=t, ref_speed_Aps=ref
        )
        assert isinstance(result, TrajectoryComparison)
        assert result.quantity == "v2_magnitude_Aps"
        assert result.rmse == pytest.approx(0.0, abs=1e-12)
        assert result.mean_ratio == pytest.approx(1.0, abs=1e-12)

    def test_atom_selection(self):
        """I1 -> speed 0.3 series; I2 -> speed 0.7 series; quantities differ."""
        t = np.linspace(0.0, 1.0, 11)
        ion = _make_ion_with_distance(
            t_md=t, distance=np.full_like(t, 9.0),
            num_molecules=4, speed_i1=0.3, speed_i2=0.7,
        )
        r1 = compare_speed_to_reference(
            ion, atom="I1", t_ref_ps=t, ref_speed_Aps=np.full_like(t, 0.3)
        )
        r2 = compare_speed_to_reference(
            ion, atom="I2", t_ref_ps=t, ref_speed_Aps=np.full_like(t, 0.7)
        )
        assert r1.quantity == "v1_magnitude_Aps"
        assert r2.quantity == "v2_magnitude_Aps"
        assert r1.rmse == pytest.approx(0.0, abs=1e-12)
        assert r2.rmse == pytest.approx(0.0, abs=1e-12)

    def test_known_constant_offset(self):
        t = np.linspace(0.0, 1.0, 21)
        ion = _make_ion_with_distance(
            t_md=t, distance=np.full_like(t, 9.0), speed_i2=0.7,
        )
        ref = np.full_like(t, 0.7 - 0.2)  # MD is 0.2 A/ps above reference
        result = compare_speed_to_reference(
            ion, atom="I2", t_ref_ps=t, ref_speed_Aps=ref
        )
        assert result.rmse == pytest.approx(0.2, abs=1e-12)

    def test_window_restricts_scoring(self):
        t = np.linspace(0.0, 10.0, 101)
        ion = _make_ion_with_distance(
            t_md=t, distance=np.full_like(t, 9.0), speed_i2=0.5,
        )
        ref = np.full_like(t, 0.5)
        outside = (t < 3.0) | (t > 7.0)
        ref[outside] += 5.0  # only matches MD inside [3, 7]
        full = compare_speed_to_reference(
            ion, atom="I2", t_ref_ps=t, ref_speed_Aps=ref
        )
        windowed = compare_speed_to_reference(
            ion, atom="I2", t_ref_ps=t, ref_speed_Aps=ref, window=(3.0, 7.0)
        )
        assert full.rmse > 1.0
        assert windowed.rmse == pytest.approx(0.0, abs=1e-12)
        assert windowed.overlap_t_min_ps == pytest.approx(3.0)
        assert windowed.overlap_t_max_ps == pytest.approx(7.0)

    def test_reference_window_subset_of_md(self):
        """Reference covering only [2.67, 6.0] (like the cleaned 9A data) while
        the MD run spans [0, 20] scores on the reference's own extent."""
        t_md = np.linspace(0.0, 20.0, 2001)
        ion = _make_ion_with_distance(
            t_md=t_md, distance=np.full_like(t_md, 9.0), speed_i2=0.7,
        )
        t_ref = np.linspace(2.67, 6.0, 334)
        ref = np.full_like(t_ref, 0.7)
        result = compare_speed_to_reference(
            ion, atom="I2", t_ref_ps=t_ref, ref_speed_Aps=ref
        )
        assert result.overlap_t_min_ps == pytest.approx(2.67)
        assert result.overlap_t_max_ps == pytest.approx(6.0)
        assert result.rmse == pytest.approx(0.0, abs=1e-12)

    def test_invalid_atom_raises(self):
        t = np.linspace(0.0, 1.0, 11)
        ion = _make_ion_with_distance(t_md=t, distance=np.full_like(t, 9.0))
        with pytest.raises(ValueError, match="atom must be"):
            compare_speed_to_reference(
                ion, atom="I3",  # type: ignore[arg-type]
                t_ref_ps=t, ref_speed_Aps=t,
            )

    def test_mismatched_reference_shapes_raise(self):
        t = np.linspace(0.0, 1.0, 11)
        ion = _make_ion_with_distance(t_md=t, distance=np.full_like(t, 9.0))
        with pytest.raises(ValueError, match="equal length"):
            compare_speed_to_reference(
                ion, atom="I2",
                t_ref_ps=t, ref_speed_Aps=np.linspace(0.0, 1.0, 5),
            )


# ===========================================================================
# End-to-end smoke test with real reference + real run (gated)
# ===========================================================================
@pytest.mark.skipif(
    not (ION_NPZ.exists() and REF_9A.exists()),
    reason="real run or 9A reference not present",
)
class TestEndToEndReal:
    def test_compare_real_ion_against_9A_reference(self):
        from i2_helium_md.simulation.run_directory import RunDirectory
        from i2_helium_md.postprocess.hedft_loader import load_hedft_trajectory

        ion = RunDirectory(RUN_DIR).load_ion()
        hedft = load_hedft_trajectory(REF_9A)
        result = compare_distance(ion, hedft)

        assert np.isfinite(result.rmse)
        assert np.isfinite(result.mean_ratio)
        assert result.num_overlap_points >= 2
        # Sanity: MD I-I separation should be roughly the same order as
        # the HeDFT trajectory.
        assert 0.5 < result.mean_ratio < 2.0
        assert result.overlap_t_min_ps >= 0.0
        assert result.overlap_t_max_ps > result.overlap_t_min_ps
