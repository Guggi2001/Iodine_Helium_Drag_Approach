"""Preset ``SimConfig`` builders.

Each preset corresponds to one of the old ``inputfiles_*/*.m`` scripts.
Start from a preset and override whichever fields you need:

    >>> cfg = single_pulse_N2000(num_molecules=500, seed=123)
    >>> cfg = single_pulse_N2000_18Angst(seed=123)
    >>> cfg = single_pulse_droplet_distribution(seed=123)
"""

from __future__ import annotations

from dataclasses import replace

from .config import SimConfig


def single_pulse_N2000(**overrides) -> SimConfig:
    """Reproduces ``inputfiles_dft_comparison/single_pulse_N2000.m``.

    This is the canonical preset for He-DFT comparison at R0 = 9 A.

    Parameters
    ----------
    **overrides
        Any ``SimConfig`` field to override from the preset default.
    """
    cfg = SimConfig(
        # --- from single_pulse_N2000.m ---
        R0_GS_angstrom=9.0,
        E_coulomb_scale=1.0,
        single_initial_position=True,
        custom_DFT_start=False,
        deltaR0_angstrom=0.0,
        T_particles_K=0.4,
        sigma_dependent_on_v=True,
        single_pulse=True,
        partner_interaction=True,
        additional_droplet_charges=0,
        highly_charged_iodine=False,
        num_molecules=2000,
        effusive_dynamics=False,
        hard_sphere_collision_mode=3,
        scattering_probability=0.004,
        geometric_scattering_crosssection_I=30.0,
        scatter_mass_neutral_amu=4.0,
        scatter_mass_ion_amu=4.0,
        geometric_scattering_crosssection_Iplus=2500.0,
        binding_energy_I_ion_eV=0.3,
        neutral_scatter_angle_std_deg=0.0,
        ion_scatter_angle_std_deg=0.0,
        mass_attach_probability=0.09,
        single_charge_ionization_allowed=False,
        use_single_droplet_size=True,
        single_droplet_size=2000,
        p_source_mbar=40.0,
        T_source_K=14.0,
        # --- from run_simulation.m ---
        Xdip_active=True,
        debug=False,
        num_neutral_export_timesteps=40,
        v_limit_m_per_s=40.0,
        sigma_ion_exponent=-2.0,
        lambda_pump_nm=630.0,
        E_diss_eV=1.556,
    )
    return replace(cfg, **overrides)


def single_pulse_N2000_18Angst(**overrides) -> SimConfig:
    """Reproduces ``inputfiles_dft_comparison/single_pulse_N2000_18Angst.m``.

    This is the fixed-droplet 18 A He-DFT comparison input file. It shares
    most settings with :func:`single_pulse_N2000`, but the MATLAB file has
    active assignments for the larger initial I-I distance, smaller ensemble,
    weaker I+ hard-sphere cross section, weaker ion binding, and lower helium
    attachment probability.

    Parameters
    ----------
    **overrides
        Any ``SimConfig`` field to override from the preset default.
    """
    cfg = single_pulse_N2000(
        # --- active differences from single_pulse_N2000_18Angst.m ---
        R0_GS_angstrom=18.0,
        num_molecules=2000,
        geometric_scattering_crosssection_Iplus= 1600.0,
        binding_energy_I_ion_eV=0.05,
        mass_attach_probability=0.005,
    )
    return replace(cfg, **overrides)


def single_pulse_droplet_distribution(**overrides) -> SimConfig:
    """Reproduces ``inputfiles_dft_comparison/single_pulse_droplet_distribution.m``.

    This preset keeps the same single-pulse ion/neutral physics as
    :func:`single_pulse_N2000`, but switches from the fixed 2000-atom
    droplet to the source-condition droplet-size sampler. It also uses
    the ground-state I2 equilibrium distance and a larger ensemble,
    matching the active assignments in the MATLAB input file.

    Parameters
    ----------
    **overrides
        Any ``SimConfig`` field to override from the preset default.
    """
    cfg = single_pulse_N2000(
        # --- active differences from single_pulse_droplet_distribution.m ---
        R0_GS_angstrom=2.666,
        E_coulomb_scale=0.8,
        single_initial_position=False,
        num_molecules=8000,
        use_single_droplet_size=False,
    )
    return replace(cfg, **overrides)
