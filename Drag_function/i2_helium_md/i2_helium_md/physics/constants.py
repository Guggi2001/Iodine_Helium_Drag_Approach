"""Physical constants used throughout the simulation.

Direct port of ``physical_constants.m``.

Units
-----
Distances:    Angstrom  [A]
Times:        picoseconds [ps]
Velocities:   Angstrom per picosecond [A/ps]   (1 A/ps = 100 m/s)
Energies:    electron volts [eV]  (internal scalar only, SI used for forces)
Masses:      atomic mass units [u] for inputs, kg internally where needed
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Fundamental constants (SI unless noted)
#
# Values are taken from CODATA 2022 / SI 2019. Where a constant is exact
# by definition under the 2019 SI redefinition (E_CHARGE, EV, K_B), we
# carry the full digits. Where measured (U, EPSILON_0), we carry the
# CODATA 2022 best value.
#
# This is a deliberate departure from the legacy MATLAB ``physical_constants.m``,
# which used 4-significant-figure approximations (eV = 1.602e-19, etc.).
# Those caused a ~10 ppm drift in computed energies that was visible in
# our regression tests against the literal MATLAB output. We accept that
# small drift in exchange for physical accuracy: the legacy values are
# wrong by ~100 ppm, our values are correct to all displayed digits.
# ---------------------------------------------------------------------------
E_CHARGE: float = -1.602176634e-19          # C           -- electron charge (exact, 2019 SI)
EPSILON_0: float = 8.8541878188e-12         # F/m         -- vacuum permittivity (CODATA 2022)
U: float = 1.66053906892e-27                # kg          -- atomic mass unit (CODATA 2022)
EV: float = 1.602176634e-19                 # J           -- 1 eV in Joules (exact, 2019 SI)
K_B: float = 1.380649e-23                   # J/K         -- Boltzmann constant (exact, 2019 SI)
HC: float = 1239.841984                     # eV*nm       -- Planck constant * speed of light

# Useful conversion
EV_PER_WAVENUMBER: float = 1.0 / 8065.543937  # eV per cm^-1  (CODATA 2022)


# ---------------------------------------------------------------------------
# Helium droplet properties
#   Source: Phys. Rev. B 58, 3341 (10.1103/PhysRevB.58.3341)
# ---------------------------------------------------------------------------
BULK_DENSITY_HELIUM: float = 0.0219         # atoms / A^3  (liquid He bulk)
DENSITY_DROPLET: float = 0.8 * BULK_DENSITY_HELIUM  # effective density inside droplet


def droplet_radius_bulk_angstrom(N):
    """Droplet radius in A from helium-atom count N, using bulk He density.

    Returns ``(3 N / (4 pi n_he))^(1/3)`` evaluated with the bulk
    density ``n_he = 0.0219 atoms/A^3``, which gives
    ``2.2173 * N^(1/3)`` (numerically).

    Used in the propagation code where the droplet radius defines
    the boundary of the solvation potential and the geometric
    scattering cross-section.

    Difference from legacy MATLAB
    -----------------------------
    The legacy MATLAB code hardcodes ``R = 2.22 * N^(1/3)``, a
    3-significant-figure rounding of the same formula. Using the
    rounded prefactor gives radii ~1200 ppm larger than the exact
    value at N=2000 (28.07 A instead of 27.97 A). We compute the
    exact value here, consistent with our principle of not
    preserving legacy approximations.

    The legacy MATLAB code uses two different droplet-density
    conventions in different files:

    * `vmi_sim_3d_neutral_propa_HeDFT_mimic.m` uses the **bulk**
      density 0.0219 atoms/A^3 (this function).
    * `generate_droplet_sizes.m` (the pickup-cell sampler) uses
      ``0.8 * 0.0219`` instead.

    We mirror this faithfully: the propagation uses bulk density
    (this function), and the pickup sampler uses the 0.8x density
    (the helper inside `sampling/droplet_sizes.py`). Whether the
    legacy choice was physically motivated or a copy-paste accident
    is not known.
    """
    return (3.0 * np.asarray(N) / (4.0 * np.pi * BULK_DENSITY_HELIUM)) ** (1.0 / 3.0)


# ---------------------------------------------------------------------------
# Iodine-specific
# ---------------------------------------------------------------------------
MASS_I_AMU: float = 127.0                   # u           -- iodine atomic mass


# ---------------------------------------------------------------------------
# Coulomb helpers (distance input in Angstrom)
# ---------------------------------------------------------------------------
def coulomb_energy(r_angstrom: np.ndarray | float) -> np.ndarray | float:
    """Coulomb interaction energy in Joules for two unit charges separated by r [A]."""
    return E_CHARGE ** 2 / (4.0 * np.pi * EPSILON_0 * r_angstrom * 1e-10)


def coulomb_velocity(
    r_angstrom: np.ndarray | float,
    mass_kg: np.ndarray | float,
) -> np.ndarray | float:
    """Velocity equivalent of Coulomb energy, sqrt(E/m), for mass in kg."""
    return np.sqrt(coulomb_energy(r_angstrom) / mass_kg)


# ---------------------------------------------------------------------------
# Force-to-acceleration unit conversion
# ---------------------------------------------------------------------------
# Given a force in eV/Angstrom and a mass in kg, the acceleration in
# Angstrom/picosecond^2 is:
#
#     a[A/ps^2] = F[eV/A] / mass[kg] * EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2
#
# Derivation:
#     F [N]       = F [eV/A] * EV [J/eV] / 1e-10 [m/A]
#     a [m/s^2]   = F [N] / m [kg]
#     a [A/ps^2]  = a [m/s^2] * 1e10 [A/m] * (1e-12 [s/ps])^2
#                 = a [m/s^2] * 1e-14
#     => a [A/ps^2] = F [eV/A] / m [kg] * EV * 1e10 * 1e-14
#                   = F [eV/A] / m [kg] * EV * 1e-4
#
# Numerical value: 1.602176634e-19 * 1e-4 = 1.602176634e-23.
# This is the unique source of truth for the conversion -- both the
# droplet_force code in leapfrog.py and the partner_interaction code in
# interactions.py should use this constant.
EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2: float = EV * 1e-4    # = 1.602176634e-23
