"""Slice 4 focused tests: drag-branch per-step accounting + gate assembly.

Covers the pieces unique to Slice 4 (not the BAOAB physics itself, which is
``test_baoab.py``):

* the amu*A^2/ps^2 -> eV dissipation conversion (the §4 baseline idiom),
* the Tier-0 checkpoint fills (defect/collisions/temperature),
* the §5.5 gate collapse (density_proportional == erf_tied today),
* per-step closure rebuild == build-once at fixed mass.
"""

from dataclasses import replace
from functools import partial

import numpy as np
import pytest

from i2_helium_md.physics.baoab import make_ion_baoab_step
from i2_helium_md.physics.constants import EV, U
from i2_helium_md.physics.drag import drag_gamma
from i2_helium_md.physics.leapfrog import make_ion_accel_fn
from i2_helium_md.presets import single_pulse_N2000_drag
from i2_helium_md.simulation.ion import drag_gate_steepness
from i2_helium_md.simulation.ion_propagation_step import (
    IonStepState,
    baoab_propagation_step,
)


def _make_state(n_atoms=4, mass_amu=200.0):
    z = np.zeros(n_atoms)
    return IonStepState(
        x=np.linspace(-1.0, 1.0, n_atoms), y=z.copy(), z=z.copy(),
        vx=np.full(n_atoms, 2.0), vy=z.copy(), vz=z.copy(),
        mass_kg=np.full(n_atoms, mass_amu * U),
        E_kin_eV=z.copy(), E_pot_eV=z.copy(),
        E_dissip_eV=z.copy(),
        E_mass_transfer_eV=z.copy(),
        E_int_eV=z.copy(),
        number_of_collisions=np.zeros(n_atoms, dtype=int),
        time_ps=0.0,
    )


class TestEvConversion:
    def test_dissipation_amu_to_eV_matches_half_m_v2(self):
        """dE_dissip [amu*A^2/ps^2] -> eV via U / 100^2 / EV.

        The conversion must equal 1/2 m v^2 computed in SI eV (m in kg, v in
        m/s) for a dE_dissip that is itself 1/2 m_amu v^2 -- catching a factor
        off by U (amu->kg) or 100^2 (A/ps->m/s).
        """
        n = 4
        m_amu, v = 200.0, 3.0
        dE_dissip_amu = np.full(n, 0.5 * m_amu * v ** 2)  # amu*A^2/ps^2

        def fake_step(pos, vel, dt):
            return pos, vel, np.zeros(n // 2), dE_dissip_amu

        state = _make_state(n_atoms=n, mass_amu=m_amu)
        cfg = single_pulse_N2000_drag(
            num_molecules=n // 2, ion_simulation_time=0.1, dt_ion=0.01,
        )
        new = baoab_propagation_step(
            state, step=fake_step, cfg=cfg, droplet_radii=np.full(n, 30.0),
        )

        expected_eV = dE_dissip_amu * U * 100.0 ** 2 / EV
        np.testing.assert_allclose(new.E_dissip_eV, expected_eV, rtol=1e-12)

        ke_eV = 0.5 * (m_amu * U) * (v * 100.0) ** 2 / EV
        np.testing.assert_allclose(new.E_dissip_eV, ke_eV, rtol=1e-12)

    def test_dissipation_accumulates(self):
        """E_dissip_eV is cumulative: it adds to the carried value."""
        n = 4
        prior = np.full(n, 0.5)

        def fake_step(pos, vel, dt):
            return pos, vel, np.zeros(n // 2), np.ones(n)

        state = replace(_make_state(n), E_dissip_eV=prior.copy())
        cfg = single_pulse_N2000_drag(num_molecules=n // 2, dt_ion=0.01)
        new = baoab_propagation_step(
            state, step=fake_step, cfg=cfg, droplet_radii=np.full(n, 30.0),
        )
        np.testing.assert_allclose(new.E_dissip_eV, prior + U * 100.0 ** 2 / EV)


class TestTierZeroFills:
    def test_defect_collisions_temperature_mass(self):
        n = 4

        def fake_step(pos, vel, dt):
            return pos, vel, np.zeros(n // 2), np.ones(n)

        state = _make_state(n)
        cfg = single_pulse_N2000_drag(num_molecules=n // 2, dt_ion=0.01)
        new = baoab_propagation_step(
            state, step=fake_step, cfg=cfg, droplet_radii=np.full(n, 30.0),
        )
        np.testing.assert_array_equal(new.E_mass_transfer_eV, 0.0)
        np.testing.assert_array_equal(new.number_of_collisions, 0)
        np.testing.assert_array_equal(new.mass_kg, state.mass_kg)  # fixed mass
        assert new.temperature_diagnostic.shape == (3,)
        assert np.all(np.isnan(new.temperature_diagnostic))


class TestGateCollapse:
    def test_density_proportional_equals_erf_tied(self):
        base = single_pulse_N2000_drag(num_molecules=2)
        s_dens = drag_gate_steepness(
            replace(base, drag_spatial_gate="density_proportional")
        )
        s_erf = drag_gate_steepness(replace(base, drag_spatial_gate="erf_tied"))
        assert s_dens == s_erf == base.potential_steepness

    def test_erf_independent_uses_its_own_steepness(self):
        base = single_pulse_N2000_drag(num_molecules=2, drag_gate_steepness=7.5)
        got = drag_gate_steepness(replace(base, drag_spatial_gate="erf_independent"))
        assert got == 7.5

    def test_sharp_gate_rejected(self):
        base = single_pulse_N2000_drag(num_molecules=2)
        with pytest.raises(NotImplementedError, match="sharp"):
            drag_gate_steepness(replace(base, drag_spatial_gate="sharp"))


class TestPerStepRebuild:
    def test_rebuild_equals_build_once_at_fixed_mass(self):
        """At fixed mass the per-step closure rebuild is bit-identical to a
        build-once reference (the driver rebuilds every step; Tier-1-ready)."""
        cfg = single_pulse_N2000_drag(num_molecules=2)
        n = 4
        radii = np.full(n, 30.0)
        charge = np.ones(n)
        m_kg = np.full(n, cfg.m_eff_amu * U)
        gamma_fn = partial(
            drag_gamma, coeffs=cfg.drag_coefficients, steepness=cfg.potential_steepness,
        )
        pos0 = (
            np.array([1.3, 0.0, -1.3, 0.0]),
            np.array([0.0, 1.3, 0.0, -1.3]),
            np.array([0.0, 3.0, 0.0, 3.0]),
        )
        vel0 = (np.zeros(n), np.zeros(n), np.zeros(n))
        dt = 0.01

        # Build once.
        acc = make_ion_accel_fn(cfg, m_kg, radii, charge)
        step_once = make_ion_baoab_step(m_kg / U, radii, acc, gamma_fn, T_eff=0.0)
        p1, v1 = pos0, vel0
        for _ in range(5):
            p1, v1, _, _ = step_once(p1, v1, dt)

        # Rebuild each step.
        p2, v2 = pos0, vel0
        for _ in range(5):
            acc_r = make_ion_accel_fn(cfg, m_kg, radii, charge)
            step_r = make_ion_baoab_step(m_kg / U, radii, acc_r, gamma_fn, T_eff=0.0)
            p2, v2, _, _ = step_r(p2, v2, dt)

        for a, b in zip(p1, p2):
            np.testing.assert_array_equal(a, b)
        for a, b in zip(v1, v2):
            np.testing.assert_array_equal(a, b)
