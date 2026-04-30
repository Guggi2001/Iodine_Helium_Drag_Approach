"""Tests for i2_helium_md/simulation/ion_propagation_step.py."""

from dataclasses import replace

import numpy as np
import pytest

from i2_helium_md.physics.constants import EV, U
from i2_helium_md.presets import single_pulse_N2000
from i2_helium_md.simulation.ion_initial_state import build_initial_ion_state
from i2_helium_md.simulation.ion_propagation_step import (
    IonStepState,
    ion_propagation_step,
)
from i2_helium_md.simulation.neutral import run_neutral_propagation


# ---------------------------------------------------------------------------
# Fixtures: build a real neutral checkpoint + initial ion state once.
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def setup_ion_state():
    """Run a tiny neutral, then a tiny ion init -- shared by all tests."""
    cfg = single_pulse_N2000(num_molecules=8, seed=42)
    cfg = replace(cfg, t_max_neutral=0.1, dt_neutral=0.01, dt_ion=0.01)

    neutral = run_neutral_propagation(cfg, verbose=False)
    ion = build_initial_ion_state(cfg, neutral, num_steps_ion=20)

    N = cfg.num_molecules
    two_N = 2 * N
    state = IonStepState(
        x=ion.positions_x[:, 0].copy(),
        y=ion.positions_y[:, 0].copy(),
        z=ion.positions_z[:, 0].copy(),
        vx=ion.velocities_x[:, 0].copy(),
        vy=ion.velocities_y[:, 0].copy(),
        vz=ion.velocities_z[:, 0].copy(),
        mass_kg=ion.mass_kg.copy(),
        E_kin_eV=ion.E_kin_eV[:, 0].copy(),
        E_pot_eV=ion.E_pot_eV[:, 0].copy(),
        E_dissip_eV=np.zeros(two_N),
        number_of_collisions=np.zeros(two_N, dtype=int),
        time_ps=0.0,
    )
    charge = np.ones(two_N)
    droplet_radii = ion.droplet_radii_angstrom.copy()
    return cfg, state, charge, droplet_radii, ion


# ===========================================================================
# Basic shapes and types
# ===========================================================================
class TestApi:
    def test_returns_ion_step_state(self, setup_ion_state):
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        rng = np.random.default_rng(0)
        new_state = ion_propagation_step(
            state, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
            prev_distance_angstrom=None, rng=rng,
        )
        assert isinstance(new_state, IonStepState)

    def test_does_not_mutate_input(self, setup_ion_state):
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        x_before = state.x.copy()
        v_before = state.vx.copy()
        m_before = state.mass_kg.copy()
        rng = np.random.default_rng(0)
        ion_propagation_step(
            state, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
            prev_distance_angstrom=None, rng=rng,
        )
        np.testing.assert_array_equal(state.x, x_before)
        np.testing.assert_array_equal(state.vx, v_before)
        np.testing.assert_array_equal(state.mass_kg, m_before)

    def test_advances_time_by_dt(self, setup_ion_state):
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        rng = np.random.default_rng(0)
        new_state = ion_propagation_step(
            state, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
            prev_distance_angstrom=None, rng=rng,
        )
        assert new_state.time_ps == pytest.approx(state.time_ps + cfg.dt_ion)

    def test_output_shapes_preserve_2N(self, setup_ion_state):
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        two_N = state.x.shape[0]
        rng = np.random.default_rng(0)
        new_state = ion_propagation_step(
            state, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
            prev_distance_angstrom=None, rng=rng,
        )
        for arr in (new_state.x, new_state.y, new_state.z,
                    new_state.vx, new_state.vy, new_state.vz,
                    new_state.mass_kg,
                    new_state.E_kin_eV, new_state.E_pot_eV, new_state.E_dissip_eV,
                    new_state.number_of_collisions):
            assert arr.shape == (two_N,), f"got shape {arr.shape}"


# ===========================================================================
# First-step behavior (prev_distance is None -> no collisions)
# ===========================================================================
class TestFirstStep:
    def test_no_collisions_first_step(self, setup_ion_state):
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        rng = np.random.default_rng(0)
        new_state = ion_propagation_step(
            state, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
            prev_distance_angstrom=None, rng=rng,
        )
        # No collisions because prev_distance is None
        np.testing.assert_array_equal(new_state.number_of_collisions, 0)
        np.testing.assert_array_equal(new_state.E_dissip_eV, 0.0)

    def test_no_mass_change_first_step(self, setup_ion_state):
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        rng = np.random.default_rng(0)
        new_state = ion_propagation_step(
            state, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
            prev_distance_angstrom=None, rng=rng,
        )
        # No collisions -> no mass attachment
        np.testing.assert_array_equal(new_state.mass_kg, state.mass_kg)


# ===========================================================================
# Reproducibility
# ===========================================================================
class TestReproducibility:
    def test_same_seed_same_result(self, setup_ion_state):
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        prev_distance = np.full(state.x.shape[0], 0.05)  # enough to allow collisions
        s1 = ion_propagation_step(
            state, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
            prev_distance_angstrom=prev_distance,
            rng=np.random.default_rng(123),
        )
        s2 = ion_propagation_step(
            state, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
            prev_distance_angstrom=prev_distance,
            rng=np.random.default_rng(123),
        )
        np.testing.assert_array_equal(s1.x, s2.x)
        np.testing.assert_array_equal(s1.vx, s2.vx)
        np.testing.assert_array_equal(s1.mass_kg, s2.mass_kg)
        np.testing.assert_array_equal(
            s1.number_of_collisions, s2.number_of_collisions
        )

    def test_different_seed_different_result(self, setup_ion_state):
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        prev_distance = np.full(state.x.shape[0], 0.05)
        s1 = ion_propagation_step(
            state, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
            prev_distance_angstrom=prev_distance,
            rng=np.random.default_rng(1),
        )
        s2 = ion_propagation_step(
            state, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
            prev_distance_angstrom=prev_distance,
            rng=np.random.default_rng(2),
        )
        # With collisions sampled differently, at least one of these
        # should differ. (The leapfrog itself is deterministic given
        # the same state; only collision/attachment depend on rng.)
        differ = (
            not np.array_equal(s1.vx, s2.vx)
            or not np.array_equal(
                s1.number_of_collisions, s2.number_of_collisions
            )
        )
        assert differ, "expected RNG-dependent divergence"


# ===========================================================================
# Energy bookkeeping
# ===========================================================================
class TestEnergyBookkeeping:
    def test_E_kin_uses_post_attachment_mass(self, setup_ion_state):
        """E_kin = ½ * m_new * v² (with NEW mass after attachment).

        Match MATLAB line 761 which uses ``mass_i(:, t_id+1)``.
        """
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        # Mass attachment requires collisions, which require prev_distance.
        # Use a generous prev_distance so most atoms collide.
        prev_distance = np.full(state.x.shape[0], 0.5)
        rng = np.random.default_rng(7)
        new_state = ion_propagation_step(
            state, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
            prev_distance_angstrom=prev_distance, rng=rng,
        )
        # Manually recompute E_kin with the NEW mass
        v_sq = new_state.vx ** 2 + new_state.vy ** 2 + new_state.vz ** 2
        E_kin_expected = 0.5 * new_state.mass_kg * v_sq * 100.0 ** 2 / EV
        np.testing.assert_allclose(
            new_state.E_kin_eV, E_kin_expected, rtol=1e-12
        )

    def test_E_dissip_is_cumulative(self, setup_ion_state):
        """E_dissip should monotonically grow when collisions occur."""
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        prev_distance = np.full(state.x.shape[0], 0.1)
        rng = np.random.default_rng(0)
        s = state
        prev_E_dissip = state.E_dissip_eV.copy()
        for _ in range(10):
            s = ion_propagation_step(
                s, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
                prev_distance_angstrom=prev_distance, rng=rng,
            )
            assert np.all(s.E_dissip_eV >= prev_E_dissip - 1e-10), (
                "E_dissip decreased -- it should be monotonically cumulative"
            )
            prev_E_dissip = s.E_dissip_eV.copy()

    def test_n_collisions_is_cumulative(self, setup_ion_state):
        """number_of_collisions should monotonically grow."""
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        prev_distance = np.full(state.x.shape[0], 0.1)
        rng = np.random.default_rng(0)
        s = state
        prev_n = state.number_of_collisions.copy()
        for _ in range(10):
            s = ion_propagation_step(
                s, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
                prev_distance_angstrom=prev_distance, rng=rng,
            )
            assert np.all(s.number_of_collisions >= prev_n)
            prev_n = s.number_of_collisions.copy()


# ===========================================================================
# Mass attachment dynamics
# ===========================================================================
class TestMassAttachment:
    def test_mass_only_increases(self, setup_ion_state):
        """Mass is only ever added (4 amu per attachment), never removed."""
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        prev_distance = np.full(state.x.shape[0], 0.1)
        rng = np.random.default_rng(0)
        s = state
        prev_mass = state.mass_kg.copy()
        for _ in range(20):
            s = ion_propagation_step(
                s, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
                prev_distance_angstrom=prev_distance, rng=rng,
            )
            assert np.all(s.mass_kg >= prev_mass), "mass decreased"
            prev_mass = s.mass_kg.copy()

    def test_mass_change_is_multiple_of_4u(self, setup_ion_state):
        """Each attachment adds exactly 4 amu (one He atom)."""
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        prev_distance = np.full(state.x.shape[0], 0.1)
        rng = np.random.default_rng(0)
        s = state
        for _ in range(20):
            s = ion_propagation_step(
                s, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
                prev_distance_angstrom=prev_distance, rng=rng,
            )
        delta_amu = (s.mass_kg - state.mass_kg) / U
        # Each delta should be a non-negative integer * 4
        n_attach = delta_amu / 4.0
        np.testing.assert_allclose(n_attach, np.round(n_attach), atol=1e-9)
        assert np.all(n_attach >= 0)

    def test_mass_attach_probability_zero_no_attachment(self, setup_ion_state):
        """With p=0, mass should never change even with many collisions."""
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        cfg2 = replace(cfg, mass_attach_probability=0.0)
        prev_distance = np.full(state.x.shape[0], 0.5)
        rng = np.random.default_rng(0)
        s = state
        for _ in range(20):
            s = ion_propagation_step(
                s, cfg=cfg2, droplet_radii=droplet_radii, charge=charge,
                prev_distance_angstrom=prev_distance, rng=rng,
            )
        np.testing.assert_array_equal(s.mass_kg, state.mass_kg)
        # And there must HAVE been collisions, otherwise the test is empty
        assert s.number_of_collisions.sum() > 0


# ===========================================================================
# Energy conservation regression test
# ===========================================================================
class TestEnergyConservation:
    def test_drift_small_without_attachment(self, setup_ion_state):
        """E_kin + E_pot + E_dissip should be conserved within
        leapfrog symplectic-error budget (~ ppm per step) when mass
        attachment is disabled.

        Mass attachment is NOT energy-conserving in our model (the
        attached helium's kinetic energy isn't tracked), so we test
        conservation only with attachment disabled.
        """
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        cfg2 = replace(cfg, mass_attach_probability=0.0)

        E_total_init = (
            state.E_kin_eV.sum()
            + state.E_pot_eV.sum()
            + state.E_dissip_eV.sum()
        )

        rng = np.random.default_rng(1)
        s = state
        prev_s = None
        for _ in range(30):
            if prev_s is None:
                pd = None
            else:
                pd = np.sqrt(
                    (s.x - prev_s.x) ** 2 + (s.y - prev_s.y) ** 2
                    + (s.z - prev_s.z) ** 2
                )
            new_s = ion_propagation_step(
                s, cfg=cfg2, droplet_radii=droplet_radii, charge=charge,
                prev_distance_angstrom=pd, rng=rng,
            )
            prev_s = s
            s = new_s

        E_total_end = (
            s.E_kin_eV.sum() + s.E_pot_eV.sum() + s.E_dissip_eV.sum()
        )
        rel_drift = (E_total_end - E_total_init) / abs(E_total_init)
        assert abs(rel_drift) < 0.005, (
            f"energy drift = {rel_drift*100:.4f}% over 30 steps; "
            f"expected < 0.5%"
        )


# ===========================================================================
# Scope checks
# ===========================================================================
class TestScopeChecks:
    def test_collision_mode_other_than_3_raises(self, setup_ion_state):
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        cfg_bad = replace(cfg, hard_sphere_collision_mode=1)
        with pytest.raises(ValueError, match="mode 3"):
            ion_propagation_step(
                state, cfg=cfg_bad, droplet_radii=droplet_radii, charge=charge,
                prev_distance_angstrom=None, rng=np.random.default_rng(0),
            )

    def test_effusive_dynamics_raises(self, setup_ion_state):
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        cfg_bad = replace(cfg, effusive_dynamics=True)
        with pytest.raises(ValueError, match="effusive_dynamics"):
            ion_propagation_step(
                state, cfg=cfg_bad, droplet_radii=droplet_radii, charge=charge,
                prev_distance_angstrom=None, rng=np.random.default_rng(0),
            )

    def test_single_charge_ionization_raises(self, setup_ion_state):
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        cfg_bad = replace(cfg, single_charge_ionization_allowed=True)
        with pytest.raises(ValueError, match="single_charge_ionization"):
            ion_propagation_step(
                state, cfg=cfg_bad, droplet_radii=droplet_radii, charge=charge,
                prev_distance_angstrom=None, rng=np.random.default_rng(0),
            )

    def test_additional_droplet_charges_raises(self, setup_ion_state):
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        cfg_bad = replace(cfg, additional_droplet_charges=2)
        with pytest.raises(ValueError, match="additional_droplet_charges"):
            ion_propagation_step(
                state, cfg=cfg_bad, droplet_radii=droplet_radii, charge=charge,
                prev_distance_angstrom=None, rng=np.random.default_rng(0),
            )


# ===========================================================================
# Velocity-dependent vs constant cross section
# ===========================================================================
class TestVelocityDependentSigma:
    def test_uses_v_dependent_sigma_when_enabled(self, setup_ion_state):
        """With sigma_dependent_on_v=True, very-slow ions should collide
        much more often than fast ones.

        We construct a state with two atoms: one very slow, one moderate.
        Run a single step with prev_distance=0.05. The slow atom should
        always collide (sigma=inf for v=0); the fast atom rarely.
        """
        cfg, state, charge, droplet_radii, _ = setup_ion_state
        # Modify state: atom 0 stops, atom 1 has speed 50 A/ps
        state2 = IonStepState(
            x=state.x.copy(), y=state.y.copy(), z=state.z.copy(),
            vx=state.vx.copy(), vy=state.vy.copy(), vz=state.vz.copy(),
            mass_kg=state.mass_kg.copy(),
            E_kin_eV=state.E_kin_eV.copy(),
            E_pot_eV=state.E_pot_eV.copy(),
            E_dissip_eV=state.E_dissip_eV.copy(),
            number_of_collisions=state.number_of_collisions.copy(),
            time_ps=state.time_ps,
        )
        # We can't mutate frozen dataclass; reconstruct
        n = state2.x.shape[0]
        vx_new = state2.vx.copy()
        vy_new = state2.vy.copy()
        vz_new = state2.vz.copy()
        vx_new[0] = 0.0; vy_new[0] = 0.0; vz_new[0] = 0.0  # slow
        vx_new[1] = 50.0; vy_new[1] = 0.0; vz_new[1] = 0.0  # fast
        state2 = replace(state2, vx=vx_new, vy=vy_new, vz=vz_new)

        # Place atoms inside droplet so collisions are physically possible
        x_new = state2.x.copy(); y_new = state2.y.copy(); z_new = state2.z.copy()
        x_new[0] = 0.0; y_new[0] = 0.0; z_new[0] = 0.0
        x_new[1] = 0.0; y_new[1] = 0.0; z_new[1] = 0.0
        state2 = replace(state2, x=x_new, y=y_new, z=z_new)

        prev_distance = np.full(n, 0.05)

        # Run many trials. Keep the same state, vary only seed.
        slow_collisions = 0
        fast_collisions = 0
        n_trials = 50
        for seed in range(n_trials):
            rng = np.random.default_rng(seed)
            new_s = ion_propagation_step(
                state2, cfg=cfg, droplet_radii=droplet_radii, charge=charge,
                prev_distance_angstrom=prev_distance, rng=rng,
            )
            slow_collisions += int(new_s.number_of_collisions[0])
            fast_collisions += int(new_s.number_of_collisions[1])

        # With sigma_dependent_on_v + sigma~v^-2 + sigma_0=2500:
        #   slow atom (v=0): sigma=inf -> always collides
        #   fast atom (v=50): sigma=1 -> p_scatter ~ 0.05*1*0.0175 ~ 0.001
        #     -> rarely collides
        assert slow_collisions == n_trials, (
            f"slow atom (v=0) should always collide; got {slow_collisions}/{n_trials}"
        )
        assert fast_collisions < n_trials // 2, (
            f"fast atom should rarely collide; got {fast_collisions}/{n_trials}"
        )
