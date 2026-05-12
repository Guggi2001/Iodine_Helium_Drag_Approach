"""Tests for the paper-v4 post-processing helpers."""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.constants import U as U_KG
from i2_helium_md.postprocess._smoothing import moving_mean
from i2_helium_md.postprocess.paper_v4 import (
    PAPER_V4_COVARIANCE_BINS,
    PAPER_V4_VELOCITY_BIN_WIDTH_APS,
    PAPER_V4_VELOCITY_MAX_APS,
    PAPER_V4_VELOCITY_SMOOTHING_WINDOW,
    load_paper_v4_radial_references,
    paper_v4_angular_pair_covariance,
    paper_v4_velocity_curve,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint


def _make_ion(
    *,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
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
        velocities_final_z=vz.astype(float),
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


def test_radial_reference_loader_reads_paper_v4_directory_and_labels(tmp_path):
    ref_dir = tmp_path / "paper_v4"
    ref_dir.mkdir()
    (ref_dir / "iplus_he_300mw_43563_radial.csv").write_text(
        "v_mps,signal_arb\n0.0,1.0\n10.0,0.5\n",
        encoding="ascii",
    )
    (ref_dir / "iplus_gas_160mw_43555_radial.csv").write_text(
        "v_mps,signal_arb\n0.0,2.0\n10.0,1.0\n",
        encoding="ascii",
    )

    refs = load_paper_v4_radial_references(ref_dir)

    assert [ref.label for ref in refs] == [
        "I+ gas 160 mW (43555)",
        "I+He 300 mW (43563)",
    ]
    np.testing.assert_allclose(refs[1].velocity_mps, [0.0, 10.0])
    np.testing.assert_allclose(refs[1].signal_arb, [1.0, 0.5])


def test_velocity_curve_matches_v4_projected_speed_bins_and_normalisation():
    ion = _make_ion(
        vx=np.array([1.0, 3.0, 99.0, 4.0]),
        vy=np.array([0.0, 4.0, 99.0, 0.0]),
        vz=np.array([100.0, 100.0, 100.0, 100.0]),
        masses_amu=np.array([127.0, 127.0, 131.0, 127.0]),
        b_outside=np.array([True, False]),
    )

    curve = paper_v4_velocity_curve(ion, mass_amu=127.0)

    expected_edges = np.arange(
        0.0,
        PAPER_V4_VELOCITY_MAX_APS + PAPER_V4_VELOCITY_BIN_WIDTH_APS,
        PAPER_V4_VELOCITY_BIN_WIDTH_APS,
    )
    np.testing.assert_allclose(curve.bin_edges_Aps, expected_edges)
    np.testing.assert_allclose(curve.bin_centers_mps, curve.bin_centers_Aps * 100.0)
    assert curve.smoothing_window == PAPER_V4_VELOCITY_SMOOTHING_WINDOW
    assert curve.num_atoms_used == 1

    expected_counts, _ = np.histogram([1.0], bins=expected_edges)
    expected_smoothed = moving_mean(expected_counts, PAPER_V4_VELOCITY_SMOOTHING_WINDOW)
    expected_shifted = expected_smoothed - expected_smoothed.min()
    expected_normalised = expected_shifted / expected_shifted.max()
    np.testing.assert_array_equal(curve.counts, expected_counts)
    np.testing.assert_allclose(curve.normalised, expected_normalised)


def test_velocity_curve_raises_when_no_mass_channel_passes():
    ion = _make_ion(
        vx=np.array([1.0, 1.0]),
        vy=np.array([0.0, 0.0]),
        vz=np.array([0.0, 0.0]),
        masses_amu=np.array([127.0, 127.0]),
    )

    with pytest.raises(ValueError, match="mass=131"):
        paper_v4_velocity_curve(ion, mass_amu=131.0)


def test_angular_pair_covariance_is_literal_v4_no_diagonal_removal():
    # Molecule 0 has both fragments in the same theta bin, proving the
    # v4 helper keeps diagonal counts. Molecule 1 is outside-filtered out.
    ion = _make_ion(
        vx=np.array([1.0, 2.0, 1.0, 2.0]),
        vy=np.array([0.0, 0.0, 0.0, 0.0]),
        vz=np.zeros(4),
        masses_amu=np.full(4, 131.0),
        b_outside=np.array([True, False]),
    )

    cov = paper_v4_angular_pair_covariance(ion, mass_amu=131.0)

    assert cov.counts.shape == (PAPER_V4_COVARIANCE_BINS, PAPER_V4_COVARIANCE_BINS)
    assert cov.num_pairs_used == 1
    assert cov.counts.sum() == 1
    assert np.trace(cov.counts) == 1
    assert cov.theta_pairs_rad.shape == (1, 2)
    expected_theta = np.mod(np.arctan2(1.0, 0.0) + np.pi, 2.0 * np.pi)
    np.testing.assert_allclose(cov.theta_pairs_rad[0], [expected_theta, expected_theta])
