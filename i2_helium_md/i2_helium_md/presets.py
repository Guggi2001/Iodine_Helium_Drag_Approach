"""Preset ``SimConfig`` builders.

Each preset corresponds to one of the old ``inputfiles_*/*.m`` scripts.
Start from a preset and override whichever fields you need:

    >>> cfg = single_pulse_N2000(num_molecules=500, seed=123)
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
