"""Pair interactions between the two atoms of each molecule.

Ports three MATLAB files:

* ``atom_interaction_potential.m``  -> :func:`atom_interaction_potential`
* ``ion_interaction_potential.m``   -> :func:`ion_interaction_potential`
* ``add_partner_interaction.m``     -> :func:`partner_interaction_neutral`
* ``add_partner_interaction_ion.m`` -> :func:`partner_interaction_ion`

Array layout convention
-----------------------
To mirror the MATLAB code, coordinate arrays have **length 2N** for N
molecules:

* indices ``0 .. N-1``     - first  atom of each molecule ("atom 1")
* indices ``N .. 2N-1``    - second atom of each molecule ("atom 2", the twin)

Returned acceleration arrays have the same 2N layout.

Units
-----
* Distances in Angstrom
* Energies  in eV
* Accelerations returned in Angstrom / picosecond^2
"""

from __future__ import annotations

import numpy as np

from ..config import SimConfig
from .constants import EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2
from .potentials import morse_I2plus_state_select, morse_X


# ---------------------------------------------------------------------------
# The unit-conversion factor for "force in eV/A on a mass in kg
# -> acceleration in A/ps^2" lives in physics/constants.py as
# EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2. We import and use it here so
# this module and leapfrog.py share a single source of truth.
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Pair potential energies (vectorised over molecules)
# ---------------------------------------------------------------------------
def atom_interaction_potential(
    dr: np.ndarray,
    cfg: SimConfig,
) -> np.ndarray:
    """I-I interaction energy for neutral atoms in the X ground state.

    Thin wrapper around :func:`morse_X` -- the original MATLAB function of
    the same name was itself just a wrapper. Kept as a separate function
    because the rest of the code uses this name (and ion_interaction_potential
    is a genuinely different thing).

    Parameters
    ----------
    dr : np.ndarray, shape (N,)
        Pairwise distances in Angstrom.
    cfg : SimConfig
        Simulation config (used for `Xdip_active`).

    Returns
    -------
    U : np.ndarray, shape (N,)
        Potential energy in eV.
    """
    return morse_X(dr, cfg)


def ion_interaction_potential(
    dr: np.ndarray,
    q1: np.ndarray,
    q2: np.ndarray,
    cfg: SimConfig,
    *,
    state_ids: np.ndarray | None = None,
) -> np.ndarray:
    """Ion-ion interaction energy.

    Direct port of ``ion_interaction_potential.m``. Implements a pure
    Coulomb term plus, in the asymmetric "one-charged, one-neutral" case,
    an optional Morse I2+ state-select term.

    The Coulomb constant ``14.39964548 eV*A`` is exactly e^2 / (4*pi*eps_0)
    converted to these units.

    Parameters
    ----------
    dr : np.ndarray, shape (N,)
        Pairwise distances in Angstrom.
    q1, q2 : np.ndarray, shape (N,)
        Per-molecule integer charges on atom 1 and atom 2.
    cfg : SimConfig
        Uses ``cfg.E_coulomb_scale`` and ``cfg.single_charge_ionization_allowed``.
    state_ids : np.ndarray or None, optional
        Per-molecule I2+ electronic state (0..3) -- only required when
        ``cfg.single_charge_ionization_allowed`` is True. Must have the same
        shape as ``dr``.

    Returns
    -------
    U : np.ndarray, shape (N,)
        Potential energy in eV.

    Notes
    -----
    In the standard I+-I+ ionic case (q1 = q2 = 1) this reduces to
    ``cfg.E_coulomb_scale * 14.4 / dr`` -- see docs/physics_background.md.
    """
    # pure Coulomb term, scaled by the empirical knob E_coulomb_scale
    U_pot = cfg.E_coulomb_scale * q1 * q2 * 14.39964548 / dr

    if cfg.single_charge_ionization_allowed:
        if state_ids is None:
            raise ValueError(
                "state_ids must be provided when "
                "cfg.single_charge_ionization_allowed is True"
            )
        # mask: exactly one atom of the pair is charged
        singly_ionized = (q1 + q2) == 1
        # morse term only applies to those pairs
        morse_term = morse_I2plus_state_select(dr, state_ids)
        U_pot = U_pot + singly_ionized * morse_term

    return U_pot


# ---------------------------------------------------------------------------
# Force calculators (the MATLAB "add_partner_interaction_*" pair)
# ---------------------------------------------------------------------------
def _split_pair_coordinates(
    x: np.ndarray, y: np.ndarray, z: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Split 2N-sized coordinate arrays into atom-1 / atom-2 halves.

    Returns (r1_x, r1_y, r1_z, r2_x, r2_y, r2_z), each of size N.
    Mirrors the MATLAB slicing ``x(1:num_particles/2)`` etc.
    """
    two_N = x.shape[0]
    if two_N % 2 != 0:
        raise ValueError(
            f"coordinate arrays must have even length (got {two_N})"
        )
    N = two_N // 2
    return x[:N], y[:N], z[:N], x[N:], y[N:], z[N:]


def _pair_geometry(
    x: np.ndarray, y: np.ndarray, z: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Compute pair separation magnitudes and unit vectors.

    Parameters
    ----------
    x, y, z : np.ndarray, shape (2N,)
        Cartesian coordinates in the standard 2N layout.

    Returns
    -------
    dr : np.ndarray, shape (N,)
        Scalar pair distance, always positive.
    dr_unit : np.ndarray, shape (N, 3)
        Unit vector pointing from atom 2 to atom 1 of each molecule.
    """
    r1x, r1y, r1z, r2x, r2y, r2z = _split_pair_coordinates(x, y, z)
    dvec = np.stack([r1x - r2x, r1y - r2y, r1z - r2z], axis=-1)  # (N, 3)
    dr = np.linalg.norm(dvec, axis=-1)                            # (N,)
    dr_unit = dvec / dr[:, None]                                  # (N, 3)
    return dr, dr_unit


def _force_from_potential_fd(
    potential_fn,
    dr: np.ndarray,
    h: float = 1e-4,
) -> np.ndarray:
    """Finite-difference force ``F = -dU/dr``, mirroring MATLAB's form.

    MATLAB did ``F = (U(r) - U(r+h))/h`` which gives the force along the
    ``+r`` direction (positive for repulsive potentials at short r).

    We preserve this exactly for byte-compatibility with the legacy code.

    Parameters
    ----------
    potential_fn : callable
        Function mapping an array of distances to an array of energies [eV].
    dr : np.ndarray, shape (N,)
        Distances at which to evaluate the force.
    h : float, optional
        Finite-difference step (default 1e-4 Angstrom, matching MATLAB).

    Returns
    -------
    F : np.ndarray, shape (N,)
        Force in eV/Angstrom.
    """
    return (potential_fn(dr) - potential_fn(dr + h)) / h


def _acceleration_from_force(
    F: np.ndarray,
    dr_unit: np.ndarray,
    mass: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert per-pair scalar force to per-atom (2N, 3) acceleration.

    The force on atom 1 is along ``+dr_unit``; on atom 2 it's along
    ``-dr_unit`` (Newton's third law).

    Parameters
    ----------
    F : np.ndarray, shape (N,)
        Scalar pair force in eV/Angstrom.
    dr_unit : np.ndarray, shape (N, 3)
        Unit vectors from atom 2 toward atom 1.
    mass : np.ndarray, shape (2N,)
        Atomic masses in kg. (MATLAB stored masses in 2N layout.)

    Returns
    -------
    ax, ay, az : np.ndarray, shape (2N,)
        Acceleration components in Angstrom/ps^2.
    """
    N = F.shape[0]

    # Duplicate scalar force for atom-2 (same magnitude) -> (2N,)
    F_full = np.concatenate([F, F])

    # Convert F [eV/A] / mass [kg] -> a [A/ps^2] using the shared constant.
    # See physics/constants.py for the derivation.
    a_magnitude = F_full / mass * EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2

    # Direction: +dr_unit for atom 1, -dr_unit for atom 2
    dr_unit_full = np.concatenate([dr_unit, -dr_unit], axis=0)  # (2N, 3)

    a_vec = a_magnitude[:, None] * dr_unit_full
    return a_vec[:, 0], a_vec[:, 1], a_vec[:, 2]


def partner_interaction_neutral(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    mass: np.ndarray,
    cfg: SimConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Accelerations and potential energies from I-I X-state Morse force.

    Port of ``add_partner_interaction.m``. The finite-difference form of
    the force is retained to match MATLAB bit-for-bit; an analytical
    version is possible and may be added in the future.

    Parameters
    ----------
    x, y, z : np.ndarray, shape (2N,)
        Cartesian coordinates in Angstrom, standard 2N layout.
    mass : np.ndarray, shape (2N,)
        Atomic masses in kg.
    cfg : SimConfig
        Simulation config (passed through to ``morse_X`` for Xdip switch).

    Returns
    -------
    ax, ay, az : np.ndarray, shape (2N,)
        Accelerations in Angstrom/ps^2.
    E_pot : np.ndarray, shape (N,)
        Per-pair potential energy in eV (one value per molecule, *not* per atom).
    """
    dr, dr_unit = _pair_geometry(x, y, z)

    def U_fn(r):
        return atom_interaction_potential(r, cfg)

    E_pot = U_fn(dr)
    F = _force_from_potential_fd(U_fn, dr)
    ax, ay, az = _acceleration_from_force(F, dr_unit, mass)
    return ax, ay, az, E_pot


def partner_interaction_ion(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    mass: np.ndarray,
    charge: np.ndarray,
    cfg: SimConfig,
    *,
    state_ids: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Accelerations and potential energies from ion-ion interaction.

    Port of ``add_partner_interaction_ion.m``.

    Parameters
    ----------
    x, y, z : np.ndarray, shape (2N,)
        Cartesian coordinates in Angstrom, standard 2N layout.
    mass : np.ndarray, shape (2N,)
        Atomic masses in kg.
    charge : np.ndarray, shape (2N,)
        Per-atom charges (0 or 1). The function splits it into q1, q2 halves.
    cfg : SimConfig
        Uses E_coulomb_scale and single_charge_ionization_allowed.
    state_ids : np.ndarray, shape (N,), optional
        Required when single-charge ionization is active.

    Returns
    -------
    ax, ay, az : np.ndarray, shape (2N,)
        Accelerations in Angstrom/ps^2.
    E_pot : np.ndarray, shape (2N,)
        Per-atom potential energy in eV. Matches MATLAB convention of
        duplicating E_pot to both atoms and dividing by 2 so the total
        energy per pair equals the bare pair potential.
    """
    dr, dr_unit = _pair_geometry(x, y, z)

    N = dr.shape[0]
    if charge.shape[0] != 2 * N:
        raise ValueError(
            f"charge must have length 2*N = {2*N}, got {charge.shape[0]}"
        )
    q1 = charge[:N]
    q2 = charge[N:]

    def U_fn(r):
        return ion_interaction_potential(r, q1, q2, cfg, state_ids=state_ids)

    E_pot_pair = U_fn(dr)
    F = _force_from_potential_fd(U_fn, dr)
    ax, ay, az = _acceleration_from_force(F, dr_unit, mass)

    # MATLAB: Epot = [Epot; Epot]/2  -- split per-pair energy across the two atoms
    E_pot_per_atom = np.concatenate([E_pot_pair, E_pot_pair]) / 2.0

    return ax, ay, az, E_pot_per_atom
