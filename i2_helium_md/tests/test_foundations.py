"""Sanity tests for the foundational modules (Steps 1-3)."""

import math

import pytest

from i2_helium_md import (
    SimConfig,
    single_pulse_N2000,
    single_pulse_N2000_18Angst,
    single_pulse_droplet_distribution,
)
from i2_helium_md.physics import EV, K_B, MASS_I_AMU, U


# ---------------------------------------------------------------------------
# Constants -- cross-check against MATLAB physical_constants.m
# ---------------------------------------------------------------------------
class TestConstants:
    def test_elementary_values(self):
        assert U == pytest.approx(1.66053906892e-27)
        assert EV == pytest.approx(1.602176634e-19)
        assert K_B == pytest.approx(1.380649e-23)
        assert MASS_I_AMU == 127.0

    def test_binding_energy_matlab_match(self):
        """MATLAB: binding_energy_I_atom = 318.43 * k_B / eV."""
        expected = 318.43 * K_B / EV
        cfg = SimConfig()
        assert cfg.binding_energy_I_atom_eV == pytest.approx(expected)

    def test_molecule_binding_energy_match(self):
        """MATLAB: binding_energy_molecule = 573.3 * k_B / eV * 1000  (in meV)."""
        expected = 573.3 * K_B / EV * 1000.0
        cfg = SimConfig()
        assert cfg.binding_energy_molecule_meV == pytest.approx(expected)


class TestDropletRadiusUtility:
    """The droplet_radius_bulk_angstrom helper used by the propagation."""

    def test_scalar_input(self):
        from i2_helium_md.physics.constants import droplet_radius_bulk_angstrom
        # Derived value: (3 / (4 pi 0.0219))^(1/3) = 2.2173
        # So R(N=2000) = 2.2173 * 2000^(1/3) = 27.94 A
        R = droplet_radius_bulk_angstrom(2000)
        assert R == pytest.approx(27.94, abs=0.01)

    def test_vector_input(self):
        import numpy as np
        from i2_helium_md.physics.constants import droplet_radius_bulk_angstrom
        Ns = np.array([1000, 2000, 5000])
        Rs = droplet_radius_bulk_angstrom(Ns)
        assert Rs.shape == (3,)
        # monotonically increasing
        assert Rs[0] < Rs[1] < Rs[2]

    def test_cube_law(self):
        """R^3 should be linear in N (constant density implication)."""
        from i2_helium_md.physics.constants import droplet_radius_bulk_angstrom
        R1 = droplet_radius_bulk_angstrom(1000)
        R8 = droplet_radius_bulk_angstrom(8000)
        # R(8N)/R(N) = 8^(1/3) = 2 exactly
        assert R8 / R1 == pytest.approx(2.0, rel=1e-12)

    def test_legacy_matlab_difference(self):
        """The legacy MATLAB hardcoded 2.22 * N^(1/3); we use the precise
        derived value. Document the ~1200 ppm offset at N=2000.
        """
        from i2_helium_md.physics.constants import droplet_radius_bulk_angstrom
        N = 2000
        R_ours = droplet_radius_bulk_angstrom(N)
        R_legacy = 2.22 * N ** (1.0 / 3.0)
        rel_diff = abs(R_ours - R_legacy) / R_legacy
        assert rel_diff < 0.002  # <0.2%, the legacy was rounded


# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------
class TestDerivedQuantities:
    def test_num_timesteps(self):
        cfg = SimConfig(t_max_neutral=200.0, dt_neutral=0.01)
        assert cfg.num_timesteps_neutral == 20000

    def test_v_limit_units(self):
        cfg = SimConfig(v_limit_m_per_s=40.0)
        # 40 m/s == 0.4 A/ps
        assert cfg.v_limit_angstrom_per_ps == pytest.approx(0.4)

    def test_e_min_matlab_formula(self):
        """MATLAB: E_min = (127*u) * v_limit^2 / 2 / eV, v in m/s."""
        cfg = SimConfig(v_limit_m_per_s=40.0)
        v = 40.0
        expected = MASS_I_AMU * U * v ** 2 / 2.0 / EV
        assert cfg.E_min_eV == pytest.approx(expected)


# ---------------------------------------------------------------------------
# Preset
# ---------------------------------------------------------------------------
class TestSinglePulsePreset:
    def test_preset_reproduces_matlab_inputfile(self):
        cfg = single_pulse_N2000()
        # spot-check every value that's explicitly set in single_pulse_N2000.m
        assert cfg.R0_GS_angstrom == 9.0
        assert cfg.E_coulomb_scale == 1.0
        assert cfg.single_initial_position is True
        assert cfg.custom_DFT_start is False
        assert cfg.deltaR0_angstrom == 0.0
        assert cfg.T_particles_K == 0.4
        assert cfg.sigma_dependent_on_v is True
        assert cfg.single_pulse is True
        assert cfg.partner_interaction is True
        assert cfg.num_molecules == 2000
        assert cfg.effusive_dynamics is False
        assert cfg.hard_sphere_collision_mode == 3
        assert cfg.geometric_scattering_crosssection_I == 30.0
        assert cfg.geometric_scattering_crosssection_Iplus == 2500.0
        assert cfg.binding_energy_I_ion_eV == 0.3
        assert cfg.mass_attach_probability == 0.09
        assert cfg.use_single_droplet_size is True
        assert cfg.single_droplet_size == 2000
        assert cfg.p_source_mbar == 40.0
        assert cfg.T_source_K == 14.0

    def test_preset_overrides_work(self):
        cfg = single_pulse_N2000(num_molecules=500, seed=42)
        assert cfg.num_molecules == 500
        assert cfg.seed == 42
        # others remain at preset defaults
        assert cfg.R0_GS_angstrom == 9.0

    def test_validate_passes_for_preset(self):
        # default preset should validate cleanly (may emit a RuntimeWarning
        # about E_min vs binding energy, but should not raise)
        single_pulse_N2000().validate()

    def test_validate_catches_bad_num_molecules(self):
        cfg = single_pulse_N2000(num_molecules=0)
        with pytest.raises(ValueError):
            cfg.validate()

    def test_validate_catches_bad_collision_mode(self):
        cfg = single_pulse_N2000()
        cfg.hard_sphere_collision_mode = 5  # type: ignore[assignment]
        with pytest.raises(ValueError):
            cfg.validate()


class TestSinglePulseDropletDistributionPreset:
    def test_preset_reproduces_matlab_inputfile(self):
        cfg = single_pulse_droplet_distribution()
        # Spot-check every active assignment in
        # inputfiles_dft_comparison/single_pulse_droplet_distribution.m.
        assert cfg.R0_GS_angstrom == pytest.approx(2.666)
        assert cfg.E_coulomb_scale == pytest.approx(0.8)
        assert cfg.custom_DFT_start is False
        assert cfg.single_initial_position is False
        assert cfg.deltaR0_angstrom == 0.0
        assert cfg.T_particles_K == 0.4
        assert cfg.sigma_dependent_on_v is True
        assert cfg.single_pulse is True
        assert cfg.partner_interaction is True
        assert cfg.additional_droplet_charges == 0
        assert cfg.highly_charged_iodine is False
        assert cfg.num_molecules == 8000
        assert cfg.effusive_dynamics is False
        assert cfg.hard_sphere_collision_mode == 3
        assert cfg.scattering_probability == pytest.approx(0.004)
        assert cfg.geometric_scattering_crosssection_I == pytest.approx(30.0)
        assert cfg.scatter_mass_neutral_amu == pytest.approx(4.0)
        assert cfg.scatter_mass_ion_amu == pytest.approx(4.0)
        assert cfg.geometric_scattering_crosssection_Iplus == pytest.approx(2500.0)
        assert cfg.binding_energy_I_ion_eV == pytest.approx(0.3)
        assert cfg.neutral_scatter_angle_std_deg == pytest.approx(0.0)
        assert cfg.ion_scatter_angle_std_deg == pytest.approx(0.0)
        assert cfg.mass_attach_probability == pytest.approx(0.09)
        assert cfg.single_charge_ionization_allowed is False
        assert cfg.use_single_droplet_size is False
        assert cfg.p_source_mbar == pytest.approx(40.0)
        assert cfg.T_source_K == pytest.approx(14.0)

    def test_preset_overrides_work(self):
        cfg = single_pulse_droplet_distribution(num_molecules=50, seed=7)
        assert cfg.num_molecules == 50
        assert cfg.seed == 7
        assert cfg.use_single_droplet_size is False
        assert cfg.R0_GS_angstrom == pytest.approx(2.666)

    def test_validate_passes_for_preset(self):
        single_pulse_droplet_distribution().validate()


class TestSinglePulseN200018AngstPreset:
    def test_preset_reproduces_matlab_inputfile(self):
        cfg = single_pulse_N2000_18Angst()
        # Spot-check every active assignment in
        # inputfiles_dft_comparison/single_pulse_N2000_18Angst.m.
        assert cfg.R0_GS_angstrom == pytest.approx(18.0)
        assert cfg.E_coulomb_scale == pytest.approx(1.0)
        assert cfg.single_initial_position is True
        assert cfg.custom_DFT_start is False
        assert cfg.deltaR0_angstrom == pytest.approx(0.0)
        assert cfg.T_particles_K == pytest.approx(0.4)
        assert cfg.sigma_dependent_on_v is True
        assert cfg.single_pulse is True
        assert cfg.partner_interaction is True
        assert cfg.additional_droplet_charges == 0
        assert cfg.highly_charged_iodine is False
        assert cfg.num_molecules == 2000
        assert cfg.effusive_dynamics is False
        assert cfg.hard_sphere_collision_mode == 3
        assert cfg.scattering_probability == pytest.approx(0.004)
        assert cfg.geometric_scattering_crosssection_I == pytest.approx(30.0)
        assert cfg.scatter_mass_neutral_amu == pytest.approx(4.0)
        assert cfg.scatter_mass_ion_amu == pytest.approx(4.0)
        assert cfg.geometric_scattering_crosssection_Iplus == pytest.approx(1600.0)
        assert cfg.binding_energy_I_ion_eV == pytest.approx(0.05)
        assert cfg.neutral_scatter_angle_std_deg == pytest.approx(0.0)
        assert cfg.ion_scatter_angle_std_deg == pytest.approx(0.0)
        assert cfg.mass_attach_probability == pytest.approx(0.005)
        assert cfg.single_charge_ionization_allowed is False
        assert cfg.use_single_droplet_size is True
        assert cfg.single_droplet_size == 2000
        assert cfg.p_source_mbar == pytest.approx(40.0)
        assert cfg.T_source_K == pytest.approx(14.0)

    def test_preset_overrides_work(self):
        cfg = single_pulse_N2000_18Angst(num_molecules=500, seed=42)
        assert cfg.num_molecules == 500
        assert cfg.seed == 42
        assert cfg.R0_GS_angstrom == pytest.approx(18.0)
