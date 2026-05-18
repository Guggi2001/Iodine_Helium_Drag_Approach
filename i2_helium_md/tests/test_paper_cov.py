"""Tests for the paper-cov post-processing helpers."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.io import savemat

from i2_helium_md.physics.constants import U as U_KG
from i2_helium_md.postprocess._smoothing import moving_mean
from i2_helium_md.postprocess.paper_cov import (
    PAPER_COV_PHI_BINS,
    PAPER_COV_PHI_SMOOTHING_WINDOW,
    PAPER_COV_SMOOTHING_WINDOW,
    PAPER_COV_TRACE_SMOOTHING_WINDOW,
    PAPER_COV_VELOCITY_BINS,
    PAPER_COV_VELOCITY_MAX_APS,
    PAPER_COV_VELOCITY_TRACE_MAX_APS,
    PAPER_COV_VELOCITY_TRACE_MIN_APS,
    PhiAngularDistribution,
    covariance_axis_sum_normalised,
    load_paper_cov_experimental_reference,
    radial_covariance_trace,
    radial_pair_speed_covariance,
    simulated_phi_distribution,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint


def _make_ion(
    *,
    vx: np.ndarray,
    vy: np.ndarray,
    masses_amu: np.ndarray,
    b_outside: np.ndarray | None = None,
) -> IonCheckpoint:
    n2 = vx.size
    n = n2 // 2
    if b_outside is None:
        b_outside = np.ones(n, dtype=bool)
    masses_kg = masses_amu.astype(float) * U_KG
    steps = 3
    zeros = np.zeros((n2, steps))
    return IonCheckpoint(
        num_molecules=n,
        time_ps=np.linspace(0.0, 1.0, steps),
        positions_x=zeros.copy(),
        positions_y=zeros.copy(),
        positions_z=zeros.copy(),
        velocities_x=zeros.copy(),
        velocities_y=zeros.copy(),
        velocities_z=zeros.copy(),
        positions_final_x=np.zeros(n2),
        positions_final_y=np.zeros(n2),
        positions_final_z=np.zeros(n2),
        velocities_final_x=vx.astype(float),
        velocities_final_y=vy.astype(float),
        velocities_final_z=np.zeros(n2),
        mass_kg=masses_kg.copy(),
        mass_final_kg=masses_kg,
        mass_history_kg=np.broadcast_to(masses_kg[:, None], (n2, steps)).copy(),
        droplet_radii_angstrom=np.full(n2, 30.0),
        E_kin_eV=zeros.copy(),
        E_pot_eV=zeros.copy(),
        E_dissip_eV=zeros.copy(),
        E_mass_attach_defect_eV=zeros.copy(),
        b_ion_outside=b_outside.astype(bool),
        relative_loss_per_ps=zeros.copy(),
        number_of_collisions=np.zeros((n2, steps), dtype=int),
        temperature_diagnostic=np.full((steps, 3), np.nan),
    )


def test_radial_pair_speed_covariance_diagonal_zeroed_and_pair_count():
    # Two molecules: atom_a speeds 3 and 6, atom_b speeds 4 and 8.
    # Pair 0: speeds (sqrt(9+16)=5, sqrt(0+0)... no, build via vx/vy)
    # Pair 0: a=(3,0)->5? sqrt(9)=3. b=(0,4)->4. Both pair speeds in bin.
    # Pair 1: a=(6,0)->6, b=(0,8)->8.
    ion = _make_ion(
        vx=np.array([3.0, 6.0, 0.0, 0.0]),
        vy=np.array([0.0, 0.0, 4.0, 8.0]),
        masses_amu=np.full(4, 131.0),
    )

    cov = radial_pair_speed_covariance(
        ion,
        mass_amu=131.0,
        n_velocity_bins=10,
        v_max_Aps=10.0,
        smoothing_window=0,
        remove_diagonal=True,
    )

    assert cov.counts.shape == (10, 10)
    assert cov.num_pairs_used == 2
    assert cov.mass_amu == 131.0
    # Diagonal zeroed.
    np.testing.assert_array_equal(np.diag(cov.counts), np.zeros(10))
    # Each pair lands in a distinct off-diagonal cell.
    assert int(cov.counts.sum()) == 2
    # Expected bins (edges 0, 1, 2, ..., 10).
    # Pair 0: (3, 4) -> bin (3, 4); pair 1: (6, 8) -> bin (6, 8).
    assert cov.counts[3, 4] == 1.0
    assert cov.counts[6, 8] == 1.0


def test_radial_pair_speed_covariance_mass_filter_drops_other_pairs():
    # Pair 0: mass 131 -> kept. Pair 1: mass 127 -> dropped.
    ion = _make_ion(
        vx=np.array([3.0, 6.0, 0.0, 0.0]),
        vy=np.array([0.0, 0.0, 4.0, 8.0]),
        masses_amu=np.array([131.0, 127.0, 131.0, 127.0]),
    )

    cov = radial_pair_speed_covariance(
        ion,
        mass_amu=131.0,
        n_velocity_bins=10,
        v_max_Aps=10.0,
        smoothing_window=0,
        remove_diagonal=False,
    )

    assert cov.num_pairs_used == 1
    assert int(cov.counts.sum()) == 1
    # Pair 0 only: (3, 4) -> bin (3, 4).
    assert cov.counts[3, 4] == 1.0


def test_radial_pair_speed_covariance_smoothing_matches_axiswise_movmean():
    # Three molecules so the histogram has off-diagonal mass to smooth.
    ion = _make_ion(
        vx=np.array([2.0, 5.0, 8.0, 0.0, 0.0, 0.0]),
        vy=np.array([0.0, 0.0, 0.0, 3.0, 6.0, 9.0]),
        masses_amu=np.full(6, 131.0),
    )

    smoothed = radial_pair_speed_covariance(
        ion,
        mass_amu=131.0,
        n_velocity_bins=10,
        v_max_Aps=10.0,
        smoothing_window=2,
        remove_diagonal=True,
    )
    raw = radial_pair_speed_covariance(
        ion,
        mass_amu=131.0,
        n_velocity_bins=10,
        v_max_Aps=10.0,
        smoothing_window=0,
        remove_diagonal=True,
    )

    expected = np.apply_along_axis(lambda c: moving_mean(c, 2), 0, raw.counts)
    expected = np.apply_along_axis(lambda r: moving_mean(r, 2), 1, expected)
    np.testing.assert_allclose(smoothed.counts, expected)


def test_radial_pair_speed_covariance_outside_filter_drops_pair_if_either_atom_inside():
    # Pair 0: both outside -> kept. Pair 1: atom_b inside -> dropped.
    ion = _make_ion(
        vx=np.array([3.0, 6.0, 0.0, 0.0]),
        vy=np.array([0.0, 0.0, 4.0, 8.0]),
        masses_amu=np.full(4, 131.0),
        b_outside=np.array([True, False]),
    )

    cov = radial_pair_speed_covariance(
        ion,
        mass_amu=131.0,
        n_velocity_bins=10,
        v_max_Aps=10.0,
        smoothing_window=0,
    )

    assert cov.num_pairs_used == 1
    assert int(cov.counts.sum()) == 1
    assert cov.counts[3, 4] == 1.0


def test_radial_covariance_trace_sums_band_and_halves():
    # 4x4 cov with velocity centers [1, 2, 3, 4]; band keeps rows 1-2.
    cov = np.array(
        [
            [10.0, 20.0, 30.0, 40.0],
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [100.0, 100.0, 100.0, 100.0],
        ]
    )
    v = np.array([1.0, 2.0, 3.0, 4.0])
    # Band (>1 & <4) selects rows at indices 1 and 2.
    trace = radial_covariance_trace(cov, v, v_min_Aps=1.0, v_max_Aps=4.0)
    expected = (cov[1] + cov[2]) / 2.0
    np.testing.assert_allclose(trace, expected)


def test_radial_covariance_trace_empty_band_returns_zeros():
    cov = np.ones((3, 3))
    v = np.array([1.0, 2.0, 3.0])
    trace = radial_covariance_trace(cov, v, v_min_Aps=10.0, v_max_Aps=20.0)
    np.testing.assert_array_equal(trace, np.zeros(3))


def test_radial_covariance_trace_uses_legacy_window_by_default():
    # Construct cov with v in [0, 25] so the default 4..22 window applies.
    n = 26
    v = np.arange(n, dtype=float)
    cov = np.tile(v, (n, 1))  # row i has value v[j]
    trace = radial_covariance_trace(cov, v)
    band = (v > PAPER_COV_VELOCITY_TRACE_MIN_APS) & (v < PAPER_COV_VELOCITY_TRACE_MAX_APS)
    expected = cov[band, :].sum(axis=0) / 2.0
    np.testing.assert_allclose(trace, expected)


def test_load_paper_cov_experimental_reference_round_trip_mat(tmp_path):
    n_theta = 8
    n_v = 6
    cov_angular = np.arange(n_theta * n_theta, dtype=float).reshape(n_theta, n_theta)
    cov_radial = np.arange(n_v * n_v, dtype=float).reshape(n_v, n_v)
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    v_mps = np.linspace(0.0, 2500.0, n_v)
    path = tmp_path / "iplus_he_covariance.mat"
    savemat(
        path,
        {
            "cov_angular": cov_angular,
            "cov_radial": cov_radial,
            "theta_centers_rad": theta,
            "velocity_centers_mps": v_mps,
        },
    )

    ref = load_paper_cov_experimental_reference(path)

    np.testing.assert_allclose(ref.cov_angular, cov_angular)
    np.testing.assert_allclose(ref.cov_radial, cov_radial)
    np.testing.assert_allclose(ref.theta_centers_rad, theta)
    np.testing.assert_allclose(ref.velocity_centers_mps, v_mps)
    np.testing.assert_allclose(ref.velocity_centers_Aps, v_mps / 100.0)
    assert ref.metadata == {}
    assert ref.source_path == path.resolve()


def test_load_paper_cov_experimental_reference_legacy_aps_key(tmp_path):
    n = 5
    cov = np.eye(n)
    theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    v_Aps = np.linspace(0.0, 25.0, n)
    path = tmp_path / "ref.npz"
    np.savez(
        path,
        cov_angular=cov,
        cov_radial=cov,
        theta_centers_rad=theta,
        velocity_centers_Aps=v_Aps,
    )

    ref = load_paper_cov_experimental_reference(path)
    np.testing.assert_allclose(ref.velocity_centers_Aps, v_Aps)
    np.testing.assert_allclose(ref.velocity_centers_mps, v_Aps * 100.0)


def test_load_paper_cov_experimental_reference_reads_sidecar_json(tmp_path):
    n = 4
    cov = np.zeros((n, n))
    path = tmp_path / "ref.npz"
    np.savez(
        path,
        cov_angular=cov,
        cov_radial=cov,
        theta_centers_rad=np.zeros(n),
        velocity_centers_mps=np.zeros(n),
    )
    sidecar = path.with_suffix(".json")
    sidecar.write_text('{"measurement_ids": [45668, 45662, 45667]}', encoding="utf-8")

    ref = load_paper_cov_experimental_reference(path)
    assert ref.metadata == {"measurement_ids": [45668, 45662, 45667]}


def test_load_paper_cov_experimental_reference_rejects_non_square_matrix(tmp_path):
    path = tmp_path / "ref.npz"
    np.savez(
        path,
        cov_angular=np.zeros((4, 5)),
        cov_radial=np.zeros((4, 4)),
        theta_centers_rad=np.zeros(4),
        velocity_centers_mps=np.zeros(4),
    )
    with pytest.raises(ValueError, match="cov_angular must be a square"):
        load_paper_cov_experimental_reference(path)


def test_radial_pair_speed_covariance_default_smoothing_window_is_two():
    # Sanity check: defaults pull from PAPER_COV_SMOOTHING_WINDOW.
    assert PAPER_COV_SMOOTHING_WINDOW == 2
    assert PAPER_COV_VELOCITY_BINS == 90
    assert PAPER_COV_VELOCITY_MAX_APS == 30.0
    assert PAPER_COV_PHI_BINS == 126
    assert PAPER_COV_PHI_SMOOTHING_WINDOW == 15
    assert PAPER_COV_TRACE_SMOOTHING_WINDOW == 3


def test_simulated_phi_distribution_unsmoothed_bins_match_manual_recipe():
    # Four mass-131 atoms placed at distinct azimuthal directions so each
    # falls into a separate phi bin (atan2(vy, vx) + pi wrapped to [0, 2*pi)).
    # atom 0: (vx=+1, vy= 0) -> phi = pi
    # atom 1: (vx= 0, vy=+1) -> phi = pi/2 + pi   = 3*pi/2
    # atom 2: (vx=-1, vy= 0) -> phi = pi + pi mod 2pi = 0
    # atom 3: (vx= 0, vy=-1) -> phi = -pi/2 + pi = pi/2
    ion = _make_ion(
        vx=np.array([1.0, 0.0, -1.0, 0.0]),
        vy=np.array([0.0, 1.0, 0.0, -1.0]),
        masses_amu=np.full(4, 131.0),
    )

    n_bins = 8  # bin width = 2*pi/8 = pi/4
    dist = simulated_phi_distribution(
        ion, mass_amu=131.0, n_phi_bins=n_bins, smoothing_window=0
    )

    assert isinstance(dist, PhiAngularDistribution)
    assert dist.phi_centers_rad.shape == (n_bins,)
    assert dist.signal_normalised.shape == (n_bins,)
    assert dist.num_samples_used == 4
    assert dist.smoothing_window == 0

    # Each pi/4 bin should contain exactly one atom; after normalisation by
    # max each occupied bin is 1.0 (since no smoothing) and the empty bins
    # are 0.0.
    expected_occupied_bins = {0, 2, 4, 6}  # phi = 0, pi/2, pi, 3*pi/2
    occupied = {int(i) for i in np.where(dist.signal_normalised > 0)[0]}
    assert occupied == expected_occupied_bins
    np.testing.assert_allclose(
        dist.signal_normalised[sorted(expected_occupied_bins)], 1.0
    )


def test_simulated_phi_distribution_applies_movmean_then_max_normalise():
    # Six mass-131 atoms (3 molecule pairs to satisfy the 2*N atom layout
    # the b_outside concat in _paper_cov_atom_selection expects) all at
    # vx > 0, vy = 0 -> phi = pi for every atom. With smoothing_window = 3
    # the central bin smooths into its neighbours giving a 3-bin peak;
    # max-normalisation sends all three bins to 1.0.
    n = 6
    ion = _make_ion(
        vx=np.full(n, 1.0),
        vy=np.zeros(n),
        masses_amu=np.full(n, 131.0),
    )
    n_bins = 8
    dist = simulated_phi_distribution(
        ion, mass_amu=131.0, n_phi_bins=n_bins, smoothing_window=3
    )

    # Manual expectation: histogram has 6 in the pi bin (index 4 of 8) and
    # 0 everywhere else; movmean(3) spreads it to bins {3, 4, 5} all
    # holding 6/3 = 2.0; max-normalise sends those three to 1.0.
    expected = np.zeros(n_bins, dtype=float)
    expected[3:6] = 1.0
    np.testing.assert_allclose(dist.signal_normalised, expected, atol=1e-12)
    assert dist.smoothing_window == 3


def test_simulated_phi_distribution_respects_outside_filter():
    # Same four atoms as the basic case, but mark the first molecule
    # b_outside=False so atoms 0 and 2 are dropped.
    ion = _make_ion(
        vx=np.array([1.0, 0.0, -1.0, 0.0]),
        vy=np.array([0.0, 1.0, 0.0, -1.0]),
        masses_amu=np.full(4, 131.0),
        b_outside=np.array([False, True]),
    )
    dist = simulated_phi_distribution(
        ion, mass_amu=131.0, n_phi_bins=8, smoothing_window=0
    )
    assert dist.num_samples_used == 2


def test_covariance_axis_sum_normalised_matches_manual_recipe():
    cov = np.array(
        [
            [1.0, 2.0, 3.0, 4.0],
            [0.0, 1.0, 0.0, 0.0],
            [4.0, 3.0, 2.0, 1.0],
            [1.0, 1.0, 1.0, 1.0],
        ]
    )

    trace = covariance_axis_sum_normalised(cov, axis=0, smoothing_window=3)

    raw_sum = cov.sum(axis=0)  # [6, 7, 6, 6]
    # MATLAB-style movmean(., 3): centred 3-wide window with shortened ends.
    expected_movmean = np.array(
        [
            np.mean(raw_sum[:2]),                       # [6, 7]
            np.mean(raw_sum[:3]),                       # [6, 7, 6]
            np.mean(raw_sum[1:4]),                      # [7, 6, 6]
            np.mean(raw_sum[2:]),                       # [6, 6]
        ]
    )
    expected = expected_movmean - expected_movmean.min()
    expected = expected / expected.max()
    np.testing.assert_allclose(trace, expected, atol=1e-12)


def test_covariance_axis_sum_normalised_axis_1_transposes_input():
    cov = np.array(
        [
            [1.0, 2.0, 3.0, 4.0],
            [0.0, 1.0, 0.0, 0.0],
            [4.0, 3.0, 2.0, 1.0],
            [1.0, 1.0, 1.0, 1.0],
        ]
    )
    trace_axis0 = covariance_axis_sum_normalised(cov, axis=0, smoothing_window=0)
    trace_axis1 = covariance_axis_sum_normalised(cov, axis=1, smoothing_window=0)
    np.testing.assert_allclose(
        trace_axis1,
        covariance_axis_sum_normalised(cov.T, axis=0, smoothing_window=0),
    )
    # Symmetric matrix would give the same trace on both axes; asymmetric
    # here so they must differ.
    assert not np.allclose(trace_axis0, trace_axis1)


def test_covariance_axis_sum_normalised_zero_input_returns_zeros():
    cov = np.zeros((5, 5))
    trace = covariance_axis_sum_normalised(cov, axis=0, smoothing_window=3)
    np.testing.assert_array_equal(trace, np.zeros(5))
