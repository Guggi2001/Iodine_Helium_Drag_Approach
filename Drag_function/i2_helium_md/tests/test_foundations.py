"""Sanity tests for the foundational modules (Steps 1-3)."""

import math

import pytest

from i2_helium_md import SimConfig, single_pulse_N2000
from i2_helium_md.physics.constants import EV, K_B, MASS_I_AMU, U


# ---------------------------------------------------------------------------
# Constants -- cross-check against MATLAB physical_constants.m
# ---------------------------------------------------------------------------
class TestConstants:
    def test_elementary_values(self):
        assert U == pytest.approx(1.66053907e-27)
        assert EV == pytest.approx(1.602e-19)
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
