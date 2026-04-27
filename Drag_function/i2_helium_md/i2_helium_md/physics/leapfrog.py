"""Velocity-Verlet time integrator.

Ports ``frog_step_neutral.m`` and ``frog_step_ion.m`` from the legacy repo.

Despite the legacy name ``frog_step``, the algorithm is actually the
velocity-Verlet integrator (kick-drift-kick form):

    x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
    v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt

This is symplectic and second-order accurate in dt, the standard choice for
long MD runs where energy conservation matters.

Design
------
The MATLAB version was a monolithic function that handled several force
sources (droplet + partner + charged-droplet) with many optional branches.
Here we separate concerns:

* :func:`velocity_verlet_step`  -- the pure algorithm, takes an ``accel_fn``
  callable and knows nothing about physics.
* :func:`neutral_acceleration` and :func:`ion_acceleration` -- the force
  assemblers, one per simulation stage.
* :func:`make_neutral_step` and :func:`make_ion_step` -- convenience
  factories that bind a config and return a single-argument step function.

Units
-----
* Positions in Angstrom
* Velocities in Angstrom/picosecond
* Accelerations in Angstrom/picosecond^2
* Masses in kg (matches ``cfg`` and ``interactions.py``)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from ..config import SimConfig
from .constants import U
from .interactions import (
    partner_interaction_ion,
    partner_interaction_neutral,
)
from .potentials import droplet_force as droplet_force_fn


# ===========================================================================
# Type aliases for readability
# ===========================================================================
Positions = tuple[np.ndarray, np.ndarray, np.ndarray]          # (x, y, z)  length 2N each
Velocities = tuple[np.ndarray, np.ndarray, np.ndarray]         # (vx, vy, vz)
Accelerations = tuple[np.ndarray, np.ndarray, np.ndarray]      # (ax, ay, az)

# An AccelFn takes positions and returns (accelerations, potential_energy_per_pair).
# The E_pot return matches MATLAB's convention; it's length N (per pair).
AccelFn = Callable[[Positions], tuple[Accelerations, np.ndarray]]


# ===========================================================================
# Conversion factors derived once and memoized as module constants
# ===========================================================================
#   droplet_force returns dU/dr in eV/Angstrom.
#   To get acceleration in Angstrom/ps^2 from mass in kg:
#       a[m/s^2]   = F[N] / m[kg]      with  F[N] = F[eV/A] * 1.602e-9
#       a[A/ps^2] = a[m/s^2] * 1e-14   (since 1 m/s^2 = 1e10 A/m * (1e-12 s/ps)^2)
#   Combined: a[A/ps^2] = F[eV/A] / m[kg] * 1.602e-9 * 1e-14
_DROPLET_FORCE_NEWTON_PER_EV_PER_ANGSTROM: float = 1.602e-9
_MPS2_TO_A_PER_PS2: float = 1e-14


# ===========================================================================
# The pure algorithm
# ===========================================================================
def velocity_verlet_step(
    pos: Positions,
    vel: Velocities,
    acc_fn: AccelFn,
    dt: float,
) -> tuple[Positions, Velocities, np.ndarray]:
    """One velocity-Verlet step (kick-drift-kick).

    Knows nothing about forces — takes an ``acc_fn`` callable that maps
    positions to accelerations. This makes the integrator reusable for
    any force model (neutral, ion, future extensions).

    Parameters
    ----------
    pos : tuple of np.ndarray
        Current positions ``(x, y, z)`` each of shape (2N,) in Angstrom.
    vel : tuple of np.ndarray
        Current velocities ``(vx, vy, vz)`` each of shape (2N,) in Angstrom/ps.
    acc_fn : callable
        Function mapping ``(x, y, z) -> ((ax, ay, az), E_pot)`` where
        accelerations have shape (2N,) in Angstrom/ps^2 and E_pot has shape (N,)
        in eV.
    dt : float
        Timestep in picoseconds.

    Returns
    -------
    new_pos : tuple of np.ndarray
        Updated positions.
    new_vel : tuple of np.ndarray
        Updated velocities.
    E_pot_end : np.ndarray
        Potential energy (length N) evaluated at the new positions. Useful
        for energy-conservation diagnostics.
    """
    x0, y0, z0 = pos
    vx0, vy0, vz0 = vel

    # --- step 1: acceleration at current position ---
    (ax0, ay0, az0), _ = acc_fn((x0, y0, z0))

    # --- step 2: drift to new position ---
    x1 = x0 + dt * vx0 + 0.5 * ax0 * dt ** 2
    y1 = y0 + dt * vy0 + 0.5 * ay0 * dt ** 2
    z1 = z0 + dt * vz0 + 0.5 * az0 * dt ** 2

    # --- step 3: acceleration at new position ---
    (ax1, ay1, az1), E_pot_end = acc_fn((x1, y1, z1))

    # --- step 4: kick velocity using the average acceleration ---
    vx1 = vx0 + 0.5 * (ax0 + ax1) * dt
    vy1 = vy0 + 0.5 * (ay0 + ay1) * dt
    vz1 = vz0 + 0.5 * (az0 + az1) * dt

    return (x1, y1, z1), (vx1, vy1, vz1), E_pot_end


# ===========================================================================
# Force assembly: droplet potential contribution
# ===========================================================================
def _droplet_acceleration(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    mass: np.ndarray,
    droplet_radii: np.ndarray,
    cfg: SimConfig,
    *,
    use_ion_binding: bool,
) -> Accelerations:
    """Acceleration due to the droplet solvation potential.

    The atom experiences a force -dU/dr pushing it toward the interior
    of the droplet (i.e. along -r_hat, where r_hat is its radial unit vector
    relative to the droplet center).

    Parameters
    ----------
    x, y, z : np.ndarray
        Atom positions in Angstrom, shape (2N,).
    mass : np.ndarray
        Atom masses in kg, shape (2N,).
    droplet_radii : np.ndarray
        Per-atom droplet radius (MATLAB replicates this for both atoms of
        a molecule), shape (2N,).
    cfg : SimConfig
        Simulation config (uses ``potential_steepness`` + binding energy).
    use_ion_binding : bool
        True for ions (uses ``binding_energy_I_ion_eV``),
        False for neutrals (uses ``binding_energy_I_atom_eV``).

    Returns
    -------
    ax, ay, az : np.ndarray
        Acceleration components in Angstrom/ps^2, shape (2N,).
    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    # guard against r=0 (shouldn't happen in physics but would produce NaNs)
    r_safe = np.where(r > 0, r, 1.0)

    depth = r - droplet_radii

    binding = (
        cfg.binding_energy_I_ion_eV if use_ion_binding
        else cfg.binding_energy_I_atom_eV
    )
    # dU/dr in eV/Angstrom
    dU_dr = droplet_force_fn(
        depth,
        steepness=cfg.potential_steepness,
        binding_energy=binding,
    )

    # F = -dU/dr acts radially. In MATLAB the sign was absorbed by
    # multiplying acceleration by the *negated* radial unit vector.
    # We do the same explicitly:
    #   F[eV/A] * 1.602e-9 = F[N];  a[m/s^2] = F/m;  a[A/ps^2] = a * 1e-14
    a_mag = (
        -dU_dr
        * _DROPLET_FORCE_NEWTON_PER_EV_PER_ANGSTROM
        / mass
        * _MPS2_TO_A_PER_PS2
    )  # (2N,), Angstrom/ps^2

    # project along the radial unit vector
    r_hat_x = x / r_safe
    r_hat_y = y / r_safe
    r_hat_z = z / r_safe

    return a_mag * r_hat_x, a_mag * r_hat_y, a_mag * r_hat_z


# ===========================================================================
# Force assembly: the neutral and ion acceleration functions
# ===========================================================================
@dataclass
class _StepContext:
    """Stateful inputs bound once per simulation (mass, radii, charge...)."""
    mass: np.ndarray
    droplet_radii: np.ndarray
    # optional:
    charge: np.ndarray | None = None
    state_ids: np.ndarray | None = None


def _neutral_accel_fn(
    pos: Positions,
    ctx: _StepContext,
    cfg: SimConfig,
) -> tuple[Accelerations, np.ndarray]:
    """Total acceleration for neutral atoms.

    Sums two contributions:
    1. Droplet solvation force (always on).
    2. I-I partner interaction (if ``cfg.partner_interaction``).
    """
    x, y, z = pos

    ax, ay, az = _droplet_acceleration(
        x, y, z, ctx.mass, ctx.droplet_radii, cfg, use_ion_binding=False,
    )

    if cfg.partner_interaction:
        ax_p, ay_p, az_p, E_pot = partner_interaction_neutral(
            x, y, z, ctx.mass, cfg,
        )
        ax = ax + ax_p
        ay = ay + ay_p
        az = az + az_p
    else:
        # half-length (per pair) = N = 2N // 2
        E_pot = np.zeros(x.shape[0] // 2)

    return (ax, ay, az), E_pot


def _ion_accel_fn(
    pos: Positions,
    ctx: _StepContext,
    cfg: SimConfig,
) -> tuple[Accelerations, np.ndarray]:
    """Total acceleration for ions.

    Sums two contributions:
    1. Droplet solvation force (with ion binding energy).
    2. Ion-ion partner interaction (pure Coulomb, or Morse+Coulomb if
       ``cfg.single_charge_ionization_allowed``).
    """
    x, y, z = pos
    if ctx.charge is None:
        raise ValueError("ion step requires ctx.charge to be set")

    ax, ay, az = _droplet_acceleration(
        x, y, z, ctx.mass, ctx.droplet_radii, cfg, use_ion_binding=True,
    )

    ax_p, ay_p, az_p, E_pot_per_atom = partner_interaction_ion(
        x, y, z, ctx.mass, ctx.charge, cfg, state_ids=ctx.state_ids,
    )
    ax = ax + ax_p
    ay = ay + ay_p
    az = az + az_p

    # For consistency with the neutral branch, return per-pair energy (length N).
    # partner_interaction_ion gives per-atom (length 2N), each worth half
    # the pair energy, so summing two halves gives the full pair energy.
    N = x.shape[0] // 2
    E_pot_per_pair = E_pot_per_atom[:N] + E_pot_per_atom[N:]

    return (ax, ay, az), E_pot_per_pair


# ===========================================================================
# Public convenience factories
# ===========================================================================
def make_neutral_step(
    cfg: SimConfig,
    mass: np.ndarray,
    droplet_radii: np.ndarray,
) -> Callable[[Positions, Velocities, float], tuple[Positions, Velocities, np.ndarray]]:
    """Build a single-call step function for neutral propagation.

    Parameters
    ----------
    cfg : SimConfig
    mass : np.ndarray, shape (2N,)
        Atom masses in kg.
    droplet_radii : np.ndarray, shape (2N,)
        Per-atom droplet radius in Angstrom.

    Returns
    -------
    step : callable
        ``step(pos, vel, dt) -> (new_pos, new_vel, E_pot_per_pair)``
    """
    ctx = _StepContext(mass=mass, droplet_radii=droplet_radii)

    def acc_fn(p: Positions) -> tuple[Accelerations, np.ndarray]:
        return _neutral_accel_fn(p, ctx, cfg)

    def step(pos, vel, dt):
        return velocity_verlet_step(pos, vel, acc_fn, dt)

    return step


def make_ion_step(
    cfg: SimConfig,
    mass: np.ndarray,
    droplet_radii: np.ndarray,
    charge: np.ndarray,
    state_ids: np.ndarray | None = None,
) -> Callable[[Positions, Velocities, float], tuple[Positions, Velocities, np.ndarray]]:
    """Build a single-call step function for ion propagation.

    Parameters
    ----------
    cfg : SimConfig
    mass : np.ndarray, shape (2N,)
        Atom masses in kg.
    droplet_radii : np.ndarray, shape (2N,)
        Per-atom droplet radius in Angstrom.
    charge : np.ndarray, shape (2N,)
        Per-atom integer charge (0 or 1).
    state_ids : np.ndarray, shape (N,), optional
        Per-molecule I2+ electronic state (0..3). Only used when
        ``cfg.single_charge_ionization_allowed`` is True.

    Returns
    -------
    step : callable
        ``step(pos, vel, dt) -> (new_pos, new_vel, E_pot_per_pair)``
    """
    ctx = _StepContext(
        mass=mass,
        droplet_radii=droplet_radii,
        charge=charge,
        state_ids=state_ids,
    )

    def acc_fn(p: Positions) -> tuple[Accelerations, np.ndarray]:
        return _ion_accel_fn(p, ctx, cfg)

    def step(pos, vel, dt):
        return velocity_verlet_step(pos, vel, acc_fn, dt)

    return step
