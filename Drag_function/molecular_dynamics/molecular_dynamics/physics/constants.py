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
# ---------------------------------------------------------------------------
E_CHARGE: float = -1.602e-19                # C           -- electron charge (signed, as in MATLAB source)
EPSILON_0: float = 8.85418781e-12           # F/m         -- vacuum permittivity
U: float = 1.66053907e-27                   # kg          -- atomic mass unit
EV: float = 1.602e-19                       # J           -- 1 eV in Joules
K_B: float = 1.380649e-23                   # J/K         -- Boltzmann constant
HC: float = 1240.0                          # eV*nm       -- Planck constant * speed of light

# Useful conversion
EV_PER_WAVENUMBER: float = 1.0 / 8065.54429  # eV per cm^-1


# ---------------------------------------------------------------------------
# Helium droplet properties
#   Source: Phys. Rev. B 58, 3341 (10.1103/PhysRevB.58.3341)
# ---------------------------------------------------------------------------
BULK_DENSITY_HELIUM: float = 0.0219         # atoms / A^3  (liquid He bulk)
DENSITY_DROPLET: float = 0.8 * BULK_DENSITY_HELIUM  # effective density inside droplet


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
