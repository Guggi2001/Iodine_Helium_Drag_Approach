"""Tests for i2_helium_md/postprocess/polar_velocity.py."""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.constants import U as U_KG
from i2_helium_md.postprocess.polar_velocity import (
    AnisotropyFit,
    BetaCurve,
    PolarHistogram,
    anisotropy_fit,
    beta_of_velocity,
    polar_velocity_histogram,
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
    if n2 % 2 != 0:
        raise AssertionError("velocity array length must be 2 * num_molecules")
    n = n2 // 2
    if b_outside is None:
        b_outside = np.ones(n, dtype=bool)
    diag_zero = np.zeros((2 * n, num_steps))
    masses_kg = masses_amu.astype(float) * U_KG
    return IonCheckpoint(
        num_molecules=n,
        time_ps=np.linspace(0.0, 1.0, num_steps),
        positions_x=np.zeros((2 * n, num_steps)),
        positions_y=np.zeros((2 * n, num_steps)),
        positions_z=np.zeros((2 * n, num_steps)),
        velocities_x=np.zeros((2 * n, num_steps)),
        velocities_y=np.zeros((2 * n, num_steps)),
        velocities_z=np.zeros((2 * n, num_steps)),
        positions_final_x=np.zeros(2 * n),
        positions_final_y=np.zeros(2 * n),
        positions_final_z=np.zeros(2 * n),
        velocities_final_x=vx.astype(float),
        velocities_final_y=vy.astype(float),
        velocities_final_z=vz.astype(float),
        mass_kg=masses_kg.copy(),
        mass_final_kg=masses_kg,
        mass_history_kg=np.broadcast_to(
            masses_kg[:, None], (2 * n, num_steps)
        ).copy(),
        droplet_radii_angstrom=np.full(2 * n, 30.0),
        E_kin_eV=diag_zero,
        E_pot_eV=diag_zero,
        E_dissip_eV=diag_zero,
        E_mass_attach_defect_eV=diag_zero,
        b_ion_outside=b_outside,
        relative_loss_per_ps=diag_zero,
        number_of_collisions=np.zeros((2 * n, num_steps), dtype=int),
        temperature_diagnostic=np.full((num_steps, 3), np.nan),
    )


class TestPolarVelocityHistogram:
    def test_counts_conservation(self):
        n_atoms = 200
        rng = np.random.default_rng(0)
        speeds = rng.uniform(1.0, 20.0, size=n_atoms)
        phi = rng.uniform(0.0, 2.0 * np.pi, size=n_atoms)
        ion = _make_ion(
            vx=speeds * np.cos(phi),
            vy=speeds * np.sin(phi),
            vz=np.zeros(n_atoms),
            masses_amu=np.full(n_atoms, 131.0),
        )
        hist = polar_velocity_histogram(
            ion, n_v_bins=40, n_phi_bins=36, mass_amu=131.0,
            v_max_Aps=25.0,
        )
        assert isinstance(hist, PolarHistogram)
        assert hist.counts.shape == (40, 36)
        assert hist.num_atoms_used == n_atoms
        assert hist.counts.sum() == n_atoms

    def test_phi_marginal_matches_1d_phi(self):
        n_atoms = 100
        rng = np.random.default_rng(1)
        speeds = rng.uniform(1.0, 15.0, size=n_atoms)
        phi = rng.uniform(0.0, 2.0 * np.pi, size=n_atoms)
        ion = _make_ion(
            vx=speeds * np.cos(phi),
            vy=speeds * np.sin(phi),
            vz=np.zeros(n_atoms),
            masses_amu=np.full(n_atoms, 131.0),
        )
        hist = polar_velocity_histogram(
            ion, n_v_bins=20, n_phi_bins=18, mass_amu=131.0,
            v_max_Aps=20.0,
        )
        # Marginalising over v should equal a 1-D phi histogram.
        from i2_helium_md.postprocess.energy_balance import phi_histogram
        phi_1d = phi_histogram(
            ion, bin_width_rad=2.0 * np.pi / 18, mass_amu=131.0,
        )
        np.testing.assert_array_equal(
            hist.counts.sum(axis=0).astype(int), phi_1d.counts,
        )

    def test_mass_filter_excludes(self):
        masses = np.array([131, 131, 135, 135], dtype=float)
        ion = _make_ion(
            vx=np.array([5.0, 5.0, 5.0, 5.0]),
            vy=np.zeros(4),
            vz=np.zeros(4),
            masses_amu=masses,
        )
        hist = polar_velocity_histogram(
            ion, n_v_bins=10, n_phi_bins=8, mass_amu=131.0, v_max_Aps=10.0,
        )
        assert hist.num_atoms_used == 2

    def test_bad_args(self):
        ion = _make_ion(
            vx=np.array([1.0, 1.0]),
            vy=np.zeros(2),
            vz=np.zeros(2),
            masses_amu=np.array([131.0, 131.0]),
        )
        with pytest.raises(ValueError, match="n_v_bins"):
            polar_velocity_histogram(ion, n_v_bins=0, n_phi_bins=4)
        with pytest.raises(ValueError, match="n_phi_bins"):
            polar_velocity_histogram(ion, n_v_bins=4, n_phi_bins=0)
        with pytest.raises(ValueError, match="v_max_Aps"):
            polar_velocity_histogram(ion, v_max_Aps=0.0)


class TestAnisotropyFit:
    def test_recovers_known_beta_for_cos2(self):
        # Synthesise a hist where N(phi) = a + b * cos(phi)^2 with known beta.
        a_true = 100.0
        b_true = 200.0  # beta = 2*200 / (2*100 + 200) = 1.0
        n_phi = 72
        phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False) + np.pi / n_phi
        signal = a_true + b_true * np.cos(phi) ** 2
        hist = PolarHistogram(
            counts=signal[None, :],
            v_centers_Aps=np.array([5.0]),
            v_edges_Aps=np.array([4.5, 5.5]),
            phi_centers_rad=phi,
            phi_edges_rad=np.linspace(0.0, 2.0 * np.pi, n_phi + 1),
            mass_amu=131.0,
            num_atoms_used=int(signal.sum()),
        )
        fit = anisotropy_fit(hist)
        assert isinstance(fit, AnisotropyFit)
        assert fit.success
        assert abs(fit.beta - 1.0) < 0.05
        assert abs(fit.a - a_true) < 5.0
        assert abs(fit.b - b_true) < 5.0

    def test_empty_input_returns_nan(self):
        hist = PolarHistogram(
            counts=np.zeros((4, 8)),
            v_centers_Aps=np.linspace(0.5, 3.5, 4),
            v_edges_Aps=np.linspace(0.0, 4.0, 5),
            phi_centers_rad=np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False),
            phi_edges_rad=np.linspace(0.0, 2.0 * np.pi, 9),
            mass_amu=131.0,
            num_atoms_used=0,
        )
        fit = anisotropy_fit(hist)
        assert not fit.success
        assert np.isnan(fit.beta)


class TestBetaOfVelocity:
    def test_low_count_bins_invalid(self):
        n_v = 6
        n_phi = 36
        phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
        counts = np.zeros((n_v, n_phi))
        # Two bins have a clear cos^2 signal; the rest are empty / sparse.
        counts[2] = 100.0 + 200.0 * np.cos(phi) ** 2
        counts[4] = 100.0 + 200.0 * np.cos(phi) ** 2
        counts[5, 0] = 1.0  # one stray count, below threshold
        hist = PolarHistogram(
            counts=counts,
            v_centers_Aps=np.linspace(0.5, 5.5, n_v),
            v_edges_Aps=np.linspace(0.0, 6.0, n_v + 1),
            phi_centers_rad=phi,
            phi_edges_rad=np.linspace(0.0, 2.0 * np.pi, n_phi + 1),
            mass_amu=131.0,
            num_atoms_used=int(counts.sum()),
        )
        curve = beta_of_velocity(hist, min_counts_per_v_bin=50)
        assert isinstance(curve, BetaCurve)
        assert curve.valid[2] and curve.valid[4]
        assert not curve.valid[0]
        assert not curve.valid[5]
        assert np.isnan(curve.beta[5])
