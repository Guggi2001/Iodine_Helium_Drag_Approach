"""Tests for i2_helium_md/physics/leapfrog.py."""

import numpy as np
import pytest

from i2_helium_md import single_pulse_N2000
from i2_helium_md.physics.constants import MASS_I_AMU, U
from i2_helium_md.physics.leapfrog import (
    make_ion_step,
    make_neutral_step,
    velocity_verlet_step,
)


# ---------------------------------------------------------------------------
# Helper: make a molecule of two iodine atoms aligned along x
# ---------------------------------------------------------------------------
def _make_molecule(R: float, N: int = 1, droplet_R: float = 50.0):
    """Create N identical molecules at the origin, aligned along x, separation R."""
    x = np.concatenate([np.full(N, +R / 2), np.full(N, -R / 2)])
    y = np.zeros(2 * N)
    z = np.zeros(2 * N)
    vx = np.zeros(2 * N)
    vy = np.zeros(2 * N)
    vz = np.zeros(2 * N)
    mass = np.full(2 * N, MASS_I_AMU * U)
    droplet_radii = np.full(2 * N, droplet_R)
    return (x, y, z), (vx, vy, vz), mass, droplet_radii


# ===========================================================================
# velocity_verlet_step -- the pure algorithm
# ===========================================================================
class TestVelocityVerlet:
    def test_zero_force_is_straight_line(self):
        """With zero acceleration, positions advance as x + v*dt."""
        x0 = np.array([0.0, 0.0])
        y0 = np.array([0.0, 0.0])
        z0 = np.array([0.0, 0.0])
        vx0 = np.array([1.0, 2.0])
        vy0 = np.array([0.0, 0.0])
        vz0 = np.array([0.0, 0.0])

        def zero_accel(p):
            return (np.zeros_like(p[0]),
                    np.zeros_like(p[0]),
                    np.zeros_like(p[0])), np.zeros(p[0].shape[0] // 2)

        dt = 0.1
        new_pos, new_vel, E_pot = velocity_verlet_step(
            (x0, y0, z0), (vx0, vy0, vz0), zero_accel, dt,
        )

        np.testing.assert_allclose(new_pos[0], [0.1, 0.2])
        np.testing.assert_allclose(new_vel[0], [1.0, 2.0])

    def test_constant_force_kinematic(self):
        """With constant acceleration, match analytical kinematics exactly.

        For constant a, velocity-Verlet is *exact* because Taylor series
        truncates at second order for position and first order for velocity.
        """
        x0 = np.array([0.0])
        y0 = np.array([0.0])
        z0 = np.array([0.0])
        v0 = np.array([1.0])
        zero = np.zeros_like(x0)

        a_const = 3.0

        def const_accel(p):
            ax = np.full_like(p[0], a_const)
            return (ax, zero.copy(), zero.copy()), np.zeros(p[0].shape[0] // 2)

        dt = 0.1
        new_pos, new_vel, _ = velocity_verlet_step(
            (x0, y0, z0), (v0, zero.copy(), zero.copy()), const_accel, dt,
        )

        # analytical: x = v0*dt + 0.5*a*dt^2,  v = v0 + a*dt
        np.testing.assert_allclose(new_pos[0][0], 1.0 * dt + 0.5 * a_const * dt ** 2)
        np.testing.assert_allclose(new_vel[0][0], 1.0 + a_const * dt)


# ===========================================================================
# Neutral integrator
# ===========================================================================
class TestNeutralStep:
    def test_atoms_at_rest_at_R_e_barely_move(self):
        """Molecule at rest at R_e feels ~ no force -> should stay almost still."""
        R_e = 2.666
        pos, vel, mass, droplet_radii = _make_molecule(R_e, N=1)
        cfg = single_pulse_N2000(Xdip_active=False)
        step = make_neutral_step(cfg, mass, droplet_radii)

        dt = 0.001
        new_pos, new_vel, _ = step(pos, vel, dt)

        # Displacement should be tiny (FD residual force only)
        dx = new_pos[0][0] - pos[0][0]
        assert abs(dx) < 1e-5, f"drift {dx} at t=dt seems too large"

    def test_dissociation_atoms_fly_apart(self):
        """With R<<R_e and zero initial velocity, atoms should accelerate apart."""
        R_compressed = 2.0    # compressed -- repulsive Morse force
        pos, vel, mass, droplet_radii = _make_molecule(R_compressed, N=1)
        cfg = single_pulse_N2000(Xdip_active=False)
        step = make_neutral_step(cfg, mass, droplet_radii)

        new_pos, new_vel, _ = step(pos, vel, dt=0.01)

        # atom 1 at +R/2 should have moved in +x and gained +vx
        assert new_pos[0][0] > pos[0][0]
        assert new_vel[0][0] > 0
        # atom 2 at -R/2 should have moved in -x and gained -vx
        assert new_pos[0][1] < pos[0][1]
        assert new_vel[0][1] < 0

    def test_energy_conservation_short_run(self):
        """For a pair oscillating around R_e, energy should be conserved
        over ~one vibration period to ~0.1% with a small dt.

        This is the core guarantee of velocity-Verlet: symplectic integration
        preserves a "shadow" Hamiltonian, bounding energy drift to O(dt^2).
        """
        # Start slightly off-equilibrium to induce a stable oscillation.
        R0 = 2.8   # just past R_e = 2.666
        pos, vel, mass, droplet_radii = _make_molecule(R0, N=1)
        cfg = single_pulse_N2000(
            Xdip_active=False,          # cleanest Morse potential
            partner_interaction=True,
        )
        step = make_neutral_step(cfg, mass, droplet_radii)

        # I2 X-state vibrational period ~ 160 fs -> run 200 fs = 0.2 ps
        dt = 0.001                     # 1 fs -- plenty of resolution
        n_steps = 200
        from i2_helium_md.physics.potentials import morse_X
        from i2_helium_md.physics.constants import EV

        def energy(pos, vel):
            x, y, z = pos
            vx, vy, vz = vel
            # kinetic
            v2 = vx**2 + vy**2 + vz**2  # (A/ps)^2
            # convert A/ps to m/s (factor 100) then KE in J then eV
            ke = 0.5 * mass * (v2 * (100.0)**2) / EV  # eV
            # potential
            dr = np.abs(x[0] - x[1])
            pe = morse_X(np.array([dr]), cfg)[0]  # per pair
            return ke.sum() + pe

        E0 = energy(pos, vel)
        E_hist = [E0]
        for _ in range(n_steps):
            pos, vel, _ = step(pos, vel, dt)
            E_hist.append(energy(pos, vel))

        E_hist = np.array(E_hist)
        drift = (E_hist.max() - E_hist.min()) / abs(E0)
        # velocity-Verlet + FD forces -> bounded drift, expect << 1%
        assert drift < 0.01, f"energy drift {drift:.2%} over {n_steps} steps"


# ===========================================================================
# Ion integrator
# ===========================================================================
class TestIonStep:
    def test_coulomb_explosion(self):
        """Two I+ ions at R=3A with zero velocity should fly apart."""
        R = 3.0
        pos, vel, mass, droplet_radii = _make_molecule(R, N=1, droplet_R=100.0)
        charge = np.ones(2, dtype=int)
        cfg = single_pulse_N2000()
        step = make_ion_step(cfg, mass, droplet_radii, charge)

        new_pos, new_vel, _ = step(pos, vel, dt=0.01)
        # atom 1 at +R/2 should accelerate in +x
        assert new_vel[0][0] > 0
        assert new_pos[0][0] > pos[0][0]

    def test_missing_charge_raises(self):
        """make_ion_step requires a charge array."""
        pos, vel, mass, droplet_radii = _make_molecule(3.0, N=1)
        cfg = single_pulse_N2000()
        # build a step with explicit None charge -- call should error
        from i2_helium_md.physics.leapfrog import _ion_accel_fn, _StepContext
        bad_ctx = _StepContext(mass=mass, droplet_radii=droplet_radii, charge=None)
        with pytest.raises(ValueError):
            _ion_accel_fn(pos, bad_ctx, cfg)


# ===========================================================================
# Droplet force effect
# ===========================================================================
class TestDropletForce:
    def test_atom_at_surface_gets_pulled_inward(self):
        """A stationary atom at the droplet surface should be pushed toward
        the interior (negative radial direction)."""
        # Single atom (2*N=2 for compatibility) sitting just outside droplet
        droplet_R = 30.0
        N = 1
        # Place atom 1 at r = 30 A along +x; atom 2 at origin (inside droplet).
        x = np.array([30.5, 0.0])
        y = np.zeros(2)
        z = np.zeros(2)
        vx = vy = vz = np.zeros(2)
        mass = np.full(2, MASS_I_AMU * U)
        droplet_radii = np.full(2, droplet_R)

        cfg = single_pulse_N2000(partner_interaction=False)
        step = make_neutral_step(cfg, mass, droplet_radii)

        new_pos, new_vel, _ = step((x, y, z), (vx, vy, vz), dt=0.5)
        # atom 1 was at x=30.5 and should be pushed inward (toward smaller x)
        assert new_vel[0][0] < 0, (
            f"atom at surface should be pulled inward; got vx={new_vel[0][0]}"
        )
