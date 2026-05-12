"""Tests for the literal post_process_single_pulse_paper_v3 helpers."""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.constants import U as U_KG
from i2_helium_md.postprocess.paper_v3 import (
    PAPER_V3_PHI_SMOOTHING_WINDOW,
    PAPER_V3_VELOCITY_SMOOTHING_WINDOW,
    load_paper_v3_phi_reference,
    load_paper_v3_radial_reference,
    matlab_max_normalise,
    paper_v3_phi_curve,
    paper_v3_velocity_curve,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint


def _make_ion(
    *,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    masses_amu: np.ndarray,
    b_outside: np.ndarray | None = None,
    num_steps: int = 4,
) -> IonCheckpoint:
    n2 = vx.size
    if n2 % 2:
        raise AssertionError("Need 2N atoms")
    n = n2 // 2
    if b_outside is None:
        b_outside = np.ones(n, dtype=bool)
    masses_kg = masses_amu.astype(float) * U_KG
    zeros = np.zeros((n2, num_steps))
    return IonCheckpoint(
        num_molecules=n,
        time_ps=np.linspace(0.0, 1.0, num_steps),
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
        mass_history_kg=np.broadcast_to(masses_kg[:, None], (n2, num_steps)).copy(),
        droplet_radii_angstrom=np.full(n2, 30.0),
        E_kin_eV=zeros.copy(),
        E_pot_eV=zeros.copy(),
        E_dissip_eV=zeros.copy(),
        E_mass_attach_defect_eV=zeros.copy(),
        b_ion_outside=np.asarray(b_outside, dtype=bool),
        relative_loss_per_ps=zeros.copy(),
        number_of_collisions=np.zeros((n2, num_steps), dtype=int),
        temperature_diagnostic=np.full((num_steps, 3), np.nan),
    )


def test_velocity_curve_uses_projected_detector_speed_and_exact_v3_bins():
    ion = _make_ion(
        vx=np.array([1.0, 3.0, 5.0, 7.0]),
        vy=np.array([2.0, 4.0, 0.0, 0.0]),
        vz=np.array([100.0, 100.0, 100.0, 100.0]),
        masses_amu=np.array([131.0, 131.0, 135.0, 127.0]),
        b_outside=np.array([True, True]),
    )

    curve = paper_v3_velocity_curve(ion, mass_amu=131.0)

    expected_projected = np.sqrt(np.array([1.0**2 + 2.0**2, 3.0**2 + 4.0**2]))
    expected_edges = np.arange(0.0, 35.0 + 0.05, 0.05)
    expected_counts, _ = np.histogram(expected_projected, bins=expected_edges)
    np.testing.assert_allclose(curve.bin_edges_Aps, expected_edges)
    np.testing.assert_array_equal(curve.counts, expected_counts)
    np.testing.assert_allclose(curve.bin_centers_mps, curve.bin_centers_Aps * 100.0)
    assert curve.num_atoms_used == 2
    assert curve.smoothing_window == PAPER_V3_VELOCITY_SMOOTHING_WINDOW


def test_velocity_curve_requires_molecule_outside_and_rounded_mass_match():
    ion = _make_ion(
        vx=np.array([1.0, 2.0, 3.0, 4.0]),
        vy=np.zeros(4),
        vz=np.zeros(4),
        masses_amu=np.array([130.6, 131.4, 131.0, 135.0]),
        b_outside=np.array([True, False]),
    )

    curve = paper_v3_velocity_curve(ion, mass_amu=131.0)

    assert curve.num_atoms_used == 2
    assert curve.counts.sum() == 2


def test_phi_curve_uses_matlab_atan2_convention_and_edges():
    ion = _make_ion(
        vx=np.array([1.0, 0.0, -1.0, 0.0]),
        vy=np.array([0.0, 1.0, 0.0, -1.0]),
        vz=np.zeros(4),
        masses_amu=np.full(4, 131.0),
    )

    curve = paper_v3_phi_curve(ion, mass_amu=131.0)

    expected_phi = np.mod(np.arctan2(ion.velocities_final_y, ion.velocities_final_x) + np.pi, 2.0 * np.pi)
    expected_edges = np.arange(0.0, 2.0 * np.pi + 0.05, 0.05)
    expected_counts, _ = np.histogram(expected_phi, bins=expected_edges)
    np.testing.assert_allclose(curve.bin_edges_rad, expected_edges)
    np.testing.assert_array_equal(curve.counts, expected_counts)
    assert curve.num_atoms_used == 4
    assert curve.smoothing_window == PAPER_V3_PHI_SMOOTHING_WINDOW


def test_matlab_max_normalise_does_not_subtract_baseline():
    values = np.array([2.0, 4.0, 8.0])
    np.testing.assert_allclose(matlab_max_normalise(values), [0.25, 0.5, 1.0])
    np.testing.assert_allclose(matlab_max_normalise(np.zeros(3)), np.zeros(3))


def test_paper_v3_radial_reference_loader_accepts_multiple_signal_columns(tmp_path):
    path = tmp_path / "radial.csv"
    path.write_text(
        "v_mps,signal_arb,signal_timescan\n"
        "0.0,0.1,0.2\n"
        "100.0,0.3,0.4\n",
        encoding="ascii",
    )

    ref = load_paper_v3_radial_reference(path)

    np.testing.assert_allclose(ref.velocity_mps, [0.0, 100.0])
    assert ref.signal.shape == (2, 2)
    assert ref.signal_labels == ("signal_arb", "signal_timescan")


def test_paper_v3_phi_reference_loader_enforces_columns(tmp_path):
    path = tmp_path / "phi.csv"
    path.write_text("phi_rad,signal_arb\n0.0,1.0\n3.14,0.5\n", encoding="ascii")

    ref = load_paper_v3_phi_reference(path)

    np.testing.assert_allclose(ref.phi_rad, [0.0, 3.14])
    np.testing.assert_allclose(ref.signal_arb, [1.0, 0.5])

    bad = tmp_path / "bad.csv"
    bad.write_text("angle,signal\n0.0,1.0\n", encoding="ascii")
    with pytest.raises(ValueError, match="phi_rad.*signal_arb"):
        load_paper_v3_phi_reference(bad)
