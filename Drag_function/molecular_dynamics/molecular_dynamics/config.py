"""Simulation configuration.

This replaces the ~36 MATLAB ``global`` variables scattered across
``run_simulation.m``, ``physical_constants.m`` and the various ``inputfiles_*/*.m``
preset scripts with a single strongly-typed dataclass.

Rule of thumb: every physical parameter lives here. Nothing else in the
code should have a tunable magic number.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional

from .physics.constants import EV, K_B


# ---------------------------------------------------------------------------
# Enumerations for readability
# ---------------------------------------------------------------------------
CollisionMode = Literal[1, 2, 3]
# 1: constant scattering probability per timestep
# 2: scatter after traveling one mean free path
# 3: scatter w/ probability sigma * dR * rho_droplet  (default, "sigma mode")


@dataclass
class SimConfig:
    """All tunable parameters for a single-pulse I2-in-He-droplet simulation.

    Defaults reproduce ``inputfiles_dft_comparison/single_pulse_N2000.m``
    combined with the constants set in ``run_simulation.m``.
    """

    # ------------------------------------------------------------------
    # Reproducibility
    # ------------------------------------------------------------------
    seed: Optional[int] = None          # None -> fresh randomness each run

    # ------------------------------------------------------------------
    # Simulation mode flags
    # ------------------------------------------------------------------
    single_pulse: bool = True
    effusive_dynamics: bool = False
    debug: bool = False                 # was global DEBUG

    # ------------------------------------------------------------------
    # Time grid
    # ------------------------------------------------------------------
    t_max_neutral: float = 200.0        # ps  (hardcoded in neutral script)
    dt_neutral: float = 0.01            # ps

    ion_simulation_time: float = 20.0   # ps  (from ion script)
    dt_ion: float = 0.01                # ps

    # number of timesteps stored from neutral run for export
    num_neutral_export_timesteps: int = 40

    # ------------------------------------------------------------------
    # Laser / pump
    # ------------------------------------------------------------------
    lambda_pump_nm: float = 630.0       # nm
    fwhm_lambda_nm: float = 33.0        # nm
    E_diss_eV: float = 1.556            # dissociation energy of I2 X state

    # ------------------------------------------------------------------
    # Molecule ensemble
    # ------------------------------------------------------------------
    num_molecules: int = 2000
    R0_GS_angstrom: float = 9.0         # ground-state equilibrium distance
    deltaR0_angstrom: float = 0.0       # width of initial R distribution
    T_particles_K: float = 0.4          # translational temperature in droplet
    single_initial_position: bool = True    # all at droplet center
    partner_interaction: bool = True    # include I-I X potential

    # ------------------------------------------------------------------
    # Droplet parameters
    # ------------------------------------------------------------------
    use_single_droplet_size: bool = True
    single_droplet_size: int = 2000     # number of He atoms per droplet
    p_source_mbar: float = 40.0         # nozzle pressure (only for size dist)
    T_source_K: float = 14.0            # nozzle temperature

    # droplet solvation potential
    potential_steepness: float = 14.2                # atoms
    potential_steepness_molecule: float = 14.3324    # from DFT fit
    binding_energy_I_atom_K: float = 318.43          # K -> converted below
    binding_energy_molecule_K: float = 573.3         # K -> converted below

    # Xdip: additional Gaussian dip in X ground-state potential
    Xdip_active: bool = True

    # ------------------------------------------------------------------
    # Landau / minimum velocity
    # ------------------------------------------------------------------
    v_limit_m_per_s: float = 40.0       # Landau velocity, m/s

    # ------------------------------------------------------------------
    # Hard-sphere collisions
    # ------------------------------------------------------------------
    hard_sphere_collision_mode: CollisionMode = 3
    scattering_probability: float = 0.004              # mode 1 only
    geometric_scattering_crosssection_I: float = 30.0  # A^2, neutral I
    geometric_scattering_crosssection_Iplus: float = 2500.0  # A^2, ion I+
    scatter_mass_neutral_amu: float = 4.0              # He mass
    scatter_mass_ion_amu: float = 4.0
    sigma_dependent_on_v: bool = True   # v-dependent cross section for ions
    sigma_ion_exponent: float = -2.0    # sigma ~ v^sigma_ion_exponent

    neutral_scatter_angle_std_deg: float = 0.0
    ion_scatter_angle_std_deg: float = 0.0

    # ------------------------------------------------------------------
    # Ion-specific
    # ------------------------------------------------------------------
    binding_energy_I_ion_eV: float = 0.3
    mass_attach_probability: float = 0.09
    single_charge_ionization_allowed: bool = False
    additional_droplet_charges: int = 0
    highly_charged_iodine: bool = False
    E_coulomb_scale: float = 1.0        # scaling factor for Coulomb potential

    # ------------------------------------------------------------------
    # HeDFT "custom start" mimic
    # ------------------------------------------------------------------
    custom_DFT_start: bool = False

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    output_dir: str = "results"
    data_dir: str = "data/reference"

    # ==================================================================
    # Derived quantities (do not set by hand)
    # ==================================================================
    @property
    def num_timesteps_neutral(self) -> int:
        """Number of leapfrog steps in the neutral propagation stage."""
        import math
        return math.ceil(self.t_max_neutral / self.dt_neutral)

    @property
    def v_limit_angstrom_per_ps(self) -> float:
        """Landau velocity converted to A/ps (1 A/ps = 100 m/s)."""
        return self.v_limit_m_per_s / 100.0

    @property
    def binding_energy_I_atom_eV(self) -> float:
        """binding_energy_I_atom_K converted to eV (MATLAB: 318.43 * k_B / eV)."""
        return self.binding_energy_I_atom_K * K_B / EV

    @property
    def binding_energy_molecule_meV(self) -> float:
        """binding_energy_molecule_K converted to meV (MATLAB: 573.3 * k_B / eV * 1000)."""
        return self.binding_energy_molecule_K * K_B / EV * 1000.0

    @property
    def E_min_eV(self) -> float:
        """Minimum allowed kinetic energy (Landau cutoff), eV.

        MATLAB: E_min = (127*u) * v_limit^2 / 2 / eV, with v_limit in A/ps.
        """
        from .physics.constants import U, EV, MASS_I_AMU
        v = self.v_limit_angstrom_per_ps * 100.0          # back to m/s
        return (MASS_I_AMU * U) * v ** 2 / 2.0 / EV

    # ==================================================================
    # Validation
    # ==================================================================
    def validate(self) -> None:
        """Sanity checks — fail fast rather than produce garbage."""
        if self.E_min_eV > self.binding_energy_I_atom_eV:
            # matches the MATLAB warning "all neutrals will escape!"
            import warnings
            warnings.warn(
                f"E_min ({self.E_min_eV*1000:.2f} meV) > binding energy "
                f"({self.binding_energy_I_atom_eV*1000:.2f} meV): all neutrals "
                "will escape the droplet.",
                RuntimeWarning,
            )
        if self.num_molecules <= 0:
            raise ValueError("num_molecules must be positive")
        if self.dt_neutral <= 0 or self.dt_ion <= 0:
            raise ValueError("timesteps must be positive")
        if self.hard_sphere_collision_mode not in (1, 2, 3):
            raise ValueError("hard_sphere_collision_mode must be 1, 2, or 3")
