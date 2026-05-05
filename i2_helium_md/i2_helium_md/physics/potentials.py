"""Analytical potentials used in the simulation.

Ports four MATLAB functions:

* ``droplet_potential.m``            -> :func:`droplet_potential`
* ``get_morse_potential_X.m``        -> :func:`morse_X`
* ``get_morse_potential_I2plus.m``   -> :func:`morse_I2plus`
* ``morse_potential_I2plus_state_select.m`` -> :func:`morse_I2plus_state_select`

Unit convention
---------------
* Distances in Angstrom
* Energies in eV
* All functions are fully vectorized over NumPy arrays for r.

Reference: J. Chem. Phys. 107, 9046 (1997); doi: 10.1063/1.475194
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import erf

from ..config import SimConfig
from .constants import EV, U


# ===========================================================================
# Private helpers
# ===========================================================================
def _gaussian(mu: float, sig: float, x: np.ndarray) -> np.ndarray:
    """Unit-height Gaussian ``exp(-(x-mu)**2 / (2*sig**2))``.

    Port of MATLAB's ``g(mu, sig, x)``.
    """
    return np.exp(-((x - mu) ** 2) / (2.0 * sig ** 2))


def _morse_a(D_e_eV: float, omega_e_cm: float) -> float:
    """Compute Morse width parameter ``a`` from spectroscopic constants.

    Formula identical to MATLAB sources::

        a = omega_e * (c * 1E2) * 2*pi / sqrt(2 * D_e[J] / mu) * 1E-10

    Parameters
    ----------
    D_e_eV : float
        Well depth in eV.
    omega_e_cm : float
        Vibrational frequency in cm^-1.

    Returns
    -------
    a : float
        Morse width parameter in 1/Angstrom.
    """
    c_m_per_s = 299_792_458.0
    mu_kg = 127.0 / 2.0 * U       # reduced mass of I2 in kg
    D_e_J = D_e_eV * EV
    return (
        omega_e_cm * (c_m_per_s * 1e2) * 2.0 * np.pi
        / np.sqrt(2.0 * D_e_J / mu_kg)
        * 1e-10
    )


# ===========================================================================
# Droplet solvation potential
# ===========================================================================
def droplet_potential(
    r: np.ndarray,
    steepness: float,
    binding_energy: float,
    offset: float = 0.0,
) -> np.ndarray:
    """Erf-based droplet solvation potential.

    Direct port of ``droplet_potential.m``::

        y = ((erf((x - beta3) / beta1) + 1) / 2) * beta2

    Parameters
    ----------
    r : np.ndarray
        Radial coordinate relative to droplet surface, in Angstrom.
        (Pass ``r_atom - droplet_radius`` so that r=0 is the surface.)
    steepness : float
        Width of the erf transition at the droplet surface (beta(1) in MATLAB).
    binding_energy : float
        Asymptotic well depth in eV -- potential goes from 0 inside droplet
        (r << 0) to ``binding_energy`` outside (r >> 0). (beta(2) in MATLAB.)
    offset : float, optional
        Shift along r, default 0 (beta(3) in MATLAB).

    Returns
    -------
    V : np.ndarray
        Solvation energy in eV, same shape as ``r``.
    """
    return ((erf((r - offset) / steepness) + 1.0) / 2.0) * binding_energy


def droplet_force(
    r: np.ndarray,
    steepness: float,
    binding_energy: float,
    offset: float = 0.0,
) -> np.ndarray:
    """Analytical derivative of :func:`droplet_potential` w.r.t. r.

    Replaces MATLAB's finite-difference form
    ``(U(r+h) - U(r)) / h``  (line 192 of vmi_sim_3d_neutral_propa_HeDFT_mimic.m)
    with the exact analytical derivative, avoiding numerical noise
    at small step sizes.

    d/dr [ (erf((r-c)/s) + 1) / 2 * E ]
        = E / (s * sqrt(pi)) * exp(-((r-c)/s)^2)
    """
    s = steepness
    return (
        binding_energy
        / (s * np.sqrt(np.pi))
        * np.exp(-(((r - offset) / s) ** 2))
    )


# ===========================================================================
# I2 X-state Morse potential
# ===========================================================================
@dataclass(frozen=True)
class MorseParams:
    """Spectroscopic Morse parameters for a diatomic electronic state."""
    D_e: float         # well depth [eV]
    omega_e: float     # vibrational frequency [cm^-1]
    omega_e_x_e: float # anharmonicity [cm^-1]
    R_e: float         # equilibrium distance [Angstrom]

    @property
    def a(self) -> float:
        """Morse width parameter [1/Angstrom]."""
        return _morse_a(self.D_e, self.omega_e)


# Canonical ground-state I2 parameters (J. Chem. Phys. 107, 9046)
I2_X_STATE = MorseParams(D_e=1.556, omega_e=214.5, omega_e_x_e=0.65, R_e=2.666)

# Lowest four I2+ potential energy curves (from morse_potential_I2plus_state_select.m)
I2PLUS_STATES: tuple[MorseParams, ...] = (
    MorseParams(D_e=2.70, omega_e=240.0, omega_e_x_e=0.69, R_e=2.61),
    MorseParams(D_e=2.03, omega_e=230.0, omega_e_x_e=0.29, R_e=2.61),
    MorseParams(D_e=1.26, omega_e=141.0, omega_e_x_e=0.32, R_e=2.95),
    MorseParams(D_e=0.56, omega_e=117.0, omega_e_x_e=0.38, R_e=2.95),
)

# Ionization-potential offsets for the four I2+ states [eV]
I2PLUS_IP_REL: tuple[float, ...] = (0.0, 0.63, 1.68, 2.44)
I2PLUS_IP_0: float = 9.36          # baseline IP of I2 [eV]


def morse_X(r: np.ndarray, cfg: SimConfig) -> np.ndarray:
    """Ground-state (X) Morse potential of neutral I2, with optional Xdip.

    Port of ``get_morse_potential_X.m``.

    Parameters
    ----------
    r : np.ndarray
        I-I separation in Angstrom.
    cfg : SimConfig
        Simulation config; only ``cfg.Xdip_active`` is read.

    Returns
    -------
    U : np.ndarray
        Potential energy in eV.
    """
    p = I2_X_STATE
    morse = p.D_e * (1.0 - np.exp(-p.a * (r - p.R_e))) ** 2
    if cfg.Xdip_active:
        # MATLAB: - 0.9 * g(9, 0.3, r) * Xdip_active
        morse = morse - 0.9 * _gaussian(mu=9.0, sig=0.3, x=r)
    return morse


def morse_I2plus(r: np.ndarray, params: MorseParams) -> np.ndarray:
    """Morse potential of a single I2+ electronic state.

    Port of ``get_morse_potential_I2plus.m`` / ``morse_potential_I2plus.m``.
    No Xdip correction here -- only the neutral X-state has that.

    Parameters
    ----------
    r : np.ndarray
        I-I separation in Angstrom.
    params : MorseParams
        Spectroscopic parameters for the target electronic state.

    Returns
    -------
    U : np.ndarray
        Potential energy in eV.
    """
    return params.D_e * (1.0 - np.exp(-params.a * (r - params.R_e))) ** 2


def morse_I2plus_state_select(
    r: np.ndarray,
    state_ids: np.ndarray,
) -> np.ndarray:
    """Per-molecule I2+ Morse potential, one state per molecule.

    Port of ``morse_potential_I2plus_state_select.m``.

    For each molecule ``i``, evaluate the Morse potential of state
    ``state_ids[i]`` at ``r[i]``, shifted to match the absolute energy:

    .. math::

        U_i(r_i) = V_{\\text{morse}}(r_i; \\text{state}_i)
                 - V_{\\text{morse}}(2.666; \\text{state}_i)
                 + \\text{IP}_0 + \\text{IP}_{\\text{rel}}(\\text{state}_i)

    Parameters
    ----------
    r : np.ndarray, shape (N,)
        I-I separation for each molecule in Angstrom.
    state_ids : np.ndarray, shape (N,)
        Integer state index for each molecule, in {0, 1, 2, 3}.

        .. note::
           MATLAB used 1-based indices {1, 2, 3, 4}.
           Python uses 0-based indices {0, 1, 2, 3}.

    Returns
    -------
    U : np.ndarray, shape (N,)
        Absolute ionic potential energy in eV.
    """
    r = np.asarray(r, dtype=float)
    state_ids = np.asarray(state_ids, dtype=int)

    if r.shape != state_ids.shape:
        raise ValueError(
            f"r and state_ids must have the same shape, "
            f"got {r.shape} and {state_ids.shape}"
        )
    if state_ids.min() < 0 or state_ids.max() >= len(I2PLUS_STATES):
        raise ValueError(
            f"state_ids must be in [0, {len(I2PLUS_STATES) - 1}], "
            f"got range [{state_ids.min()}, {state_ids.max()}]"
        )

    # Evaluate vectorised per-state, then pick with fancy indexing.
    # This is O(N * n_states) but n_states=4, so it's essentially O(N).
    per_state = np.stack(
        [morse_I2plus(r, p) - morse_I2plus(np.array(2.666), p) for p in I2PLUS_STATES],
        axis=0,
    )  # shape (n_states, N)

    ip_rel = np.asarray(I2PLUS_IP_REL)      # shape (n_states,)

    idx = np.arange(len(r))
    return per_state[state_ids, idx] + I2PLUS_IP_0 + ip_rel[state_ids]
