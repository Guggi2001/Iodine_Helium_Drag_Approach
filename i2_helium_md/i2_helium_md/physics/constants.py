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
# Helium / iodine-ion masses for the Tier-1a He-shell schedule
#
# These two values are the mass reference for the drag-port Tier-1a anchored
# kinematic schedule (TIER1A_IMPLEMENTATION_PLAN.md §2, §10). The complex mass
# at shell count n is ``MASS_I_ION_AMU + n * MASS_HE_AMU``.
#
# MASS_I_ION_AMU is intentionally 126.90, NOT the rounded ``MASS_I_AMU = 127.0``
# used by the neutral / ion MD. The Tier-1a §10 golden oracle (pre-shed masses
# 210.955 ... 182.936, kick factors, the 1.1532 telescoping invariant) is
# authored at this precision; 127.0 would shift the absolute shell masses by
# ~0.1 amu and break the oracle. Keep both: 127.0 stays the MD iodine mass,
# 126.90 is the shell-schedule reference. Do not "unify" them.
# ---------------------------------------------------------------------------
MASS_HE_AMU: float = 4.0026                 # u           -- helium atomic mass
MASS_I_ION_AMU: float = 126.90              # u           -- I+ mass (Tier-1a shell schedule)


# ---------------------------------------------------------------------------
# Tier-2 Phase-A dissociation ladder (Slice L)
#
# The Form-U single-rung dissociation cost D_0(n) for I+He_n interpolates a
# picture-keyed first rung D_0(1) down to the bulk-He floor across a sigmoid
# cliff centred at n* + 1/2. The first rung and the floor are *sourced* anchors
# (cm^-1, as the IHe05 EPAPS / MASS doc give them); the module converts to eV
# via EV_PER_WAVENUMBER. These are new sourced anchors for the drag-port phase,
# NOT edits to the MD physical-constants table.
#   Provenance: docs/drag_port/Tier2/TIER2_PHASE_A_IMPLEMENTATION_PLAN.md §2.1;
#   MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md R3 / A10.
# ---------------------------------------------------------------------------
N_STAR: int = 21                            # count -- first-shell I+ cation occupancy (MASS OQ8)

D0_1_X2_WAVENUMBER: float = 106.9           # cm^-1 -- IHe05 EPAPS exact J=0 ZPE (MASS R3, +/-3 cm^-1)
D0_1_MIX_WAVENUMBER: float = 74.4           # cm^-1 -- statistical SO mixture (X2+I1+I0)/3 (MASS A10)
# Provisional cooling_relaxed first rung -- a live arm whose *value* is the
# rule-2 carry (the concrete relaxation-weighted blend is pinned at Phase F).
# The source leaves it unpinned ("e.g. a relaxation-weighted blend", strictly
# between mix and X2; MASS rev-2026-06-21 #3); Phase A ships the arithmetic mean
# of the two pinned pictures so the ordering mix < cooling_relaxed < X2 holds.
# Do NOT treat as a 4-figure oracle.
D0_1_COOLING_RELAXED_WAVENUMBER: float = 0.5 * (D0_1_X2_WAVENUMBER + D0_1_MIX_WAVENUMBER)  # 90.65 cm^-1
D_FLOOR_WAVENUMBER: float = 4.97            # cm^-1 -- bulk-He chemical potential |mu_He^bulk| (MASS R3)

# Collective first-shell solvation magnitude |S| for the I+ cation (Tier-2 Phase-A
# Slice K). The cooling asymptote is E_inf(N) = -|S(N)| with the occupancy-resolved
# |S(N)| = |S| * Sigma(N)/Sigma(n*). Sourced eV-primary (the DFT source is native eV);
# = 2484 cm^-1 total, 118 cm^-1/atom (170 K) at n* = 21. A *sourced* anchor (CALIBRATION
# row 12, never tuned); the config field ``solv_struct_asymptote_eV`` defaults to it.
#   Provenance: TIER2_PHASE_A_IMPLEMENTATION_PLAN.md §2.1; MASS K2.
S_ABS_EV: float = 0.308                     # eV -- |S|, DFT first-shell solvation of I+ (MASS K2)


# ---------------------------------------------------------------------------
# Tier-2 Phase-B evaporation channel (Slice Q)
#
# The RRK unimolecular shed prefactor nu for the energy-gated evaporation
# channel: k(E_int, n) = nu * (1 - D_0(n)/E_int)^(s-1) for n >= 2 (saturating,
# bounded to [0, nu)); k = nu for the n=1 direct-dissociation branch. Unlike the
# ladder anchors (sourced in cm^-1 and converted to eV), nu is a *rate*, sourced
# and pinned directly in ps^-1 -- a rate-primary float, NOT wavenumber-converted.
# The config field ``evap_rate_prefactor_per_ps`` defaults to this anchor (the
# Slice-K |S| single-source precedent). Sourced/pinned ([IHe05] curvature, MASS
# A11); never tuned. Provenance: TIER2_PHASE_B_IMPLEMENTATION_PLAN.md §2.1;
# MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md A11.
# ---------------------------------------------------------------------------
NU_EVAP_PER_PS: float = 2.42                 # ps^-1 -- RRK evaporation prefactor (MASS A11, pinned)


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
