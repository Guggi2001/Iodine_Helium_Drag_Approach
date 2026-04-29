"""Tests for i2_helium_md/simulation/propagation_step.py (pure-function API)."""

from dataclasses import replace

import numpy as np
import pytest

from i2_helium_md.presets import single_pulse_N2000
from i2_helium_md.simulation.initial_state import build_initial_state
from i2_helium_md.simulation.propagation_step import (
    NeutralStepState,
    neutral_propagation_step,
    state_from_checkpoint_column,
    write_state_to_checkpoint_column,
)


def _bootstrap(num_molecules=5, num_steps=10, *, seed=42, R0=9.0,
               single_pulse=True):
    cfg = single_pulse_N2000(num_molecules=num_molecules, seed=seed,
                              R0_GS_angstrom=R0)
    if not single_pulse:
        cfg = replace(cfg, single_pulse=False)
    ckpt = build_initial_state(cfg, num_steps=num_steps,
                                rng=np.random.default_rng(seed))
    state0 = state_from_checkpoint_column(ckpt, 0)
    return cfg, ckpt, state0


def _step_distance(prev: NeutralStepState, new: NeutralStepState) -> np.ndarray:
    return np.sqrt(
        (new.x - prev.x) ** 2
        + (new.y - prev.y) ** 2
        + (new.z - prev.z) ** 2
    )


# ===========================================================================
class TestValidation:
    def test_collision_mode_not_3_raises(self):
        cfg, ckpt, state0 = _bootstrap()
        cfg = replace(cfg, hard_sphere_collision_mode=1)
        with pytest.raises(ValueError, match="mode 3"):
            neutral_propagation_step(
                state0, cfg=cfg, mass_kg=ckpt.mass_kg,
                droplet_radii=ckpt.droplet_radii,
                prev_distance_angstrom=None,
                rng=np.random.default_rng(0),
            )


class TestPurity:
    def test_input_state_unchanged(self):
        cfg, ckpt, state0 = _bootstrap()
        x_before = state0.x.copy()
        E_dissip_before = state0.E_dissip_eV.copy()
        _ = neutral_propagation_step(
            state0, cfg=cfg, mass_kg=ckpt.mass_kg,
            droplet_radii=ckpt.droplet_radii,
            prev_distance_angstrom=None,
            rng=np.random.default_rng(0),
        )
        np.testing.assert_array_equal(state0.x, x_before)
        np.testing.assert_array_equal(state0.E_dissip_eV, E_dissip_before)


class TestTimeAdvance:
    def test_time_advances_by_dt(self):
        cfg, ckpt, state0 = _bootstrap()
        s1 = neutral_propagation_step(
            state0, cfg=cfg, mass_kg=ckpt.mass_kg,
            droplet_radii=ckpt.droplet_radii,
            prev_distance_angstrom=None,
            rng=np.random.default_rng(0),
        )
        assert s1.time_ps == pytest.approx(state0.time_ps + cfg.dt_neutral)


class TestNoCollisionsInFirstStep:
    def test_E_dissip_zero_when_prev_distance_None(self):
        cfg, ckpt, state0 = _bootstrap(num_molecules=20, single_pulse=False)
        s1 = neutral_propagation_step(
            state0, cfg=cfg, mass_kg=ckpt.mass_kg,
            droplet_radii=ckpt.droplet_radii,
            prev_distance_angstrom=None,
            rng=np.random.default_rng(0),
        )
        np.testing.assert_array_equal(
            s1.E_dissip_eV, np.zeros(2 * ckpt.num_molecules),
        )


class TestEnergyConservationNoCollisions:
    def test_total_energy_drift_small(self):
        cfg, ckpt, state0 = _bootstrap(num_molecules=5, single_pulse=True)
        E0 = (state0.E_kin_eV + state0.E_pot_eV).sum()
        state = state0
        prev_dist = None
        for t in range(9):
            new = neutral_propagation_step(
                state, cfg=cfg, mass_kg=ckpt.mass_kg,
                droplet_radii=ckpt.droplet_radii,
                prev_distance_angstrom=prev_dist,
                rng=np.random.default_rng(t),
            )
            prev_dist = _step_distance(state, new)
            state = new
        if state.E_dissip_eV.sum() == 0:
            E_end = (state.E_kin_eV + state.E_pot_eV).sum()
            assert abs(E_end - E0) / abs(E0) < 0.01


class TestCumulativeBookkeeping:
    def test_E_dissip_monotone(self):
        cfg, ckpt, state0 = _bootstrap(num_molecules=20, num_steps=30,
                                        R0=2.666, single_pulse=False)
        state = state0
        prev_dist = None
        prev_dissip = state0.E_dissip_eV
        for t in range(29):
            new = neutral_propagation_step(
                state, cfg=cfg, mass_kg=ckpt.mass_kg,
                droplet_radii=ckpt.droplet_radii,
                prev_distance_angstrom=prev_dist,
                rng=np.random.default_rng(t),
            )
            assert np.all(new.E_dissip_eV >= prev_dissip - 1e-12)
            prev_dist = _step_distance(state, new)
            prev_dissip = new.E_dissip_eV
            state = new

    def test_L_droplet_monotone(self):
        cfg, ckpt, state0 = _bootstrap(num_molecules=10, num_steps=10)
        state = state0
        prev_dist = None
        prev_L = state0.L_droplet_eV_ps
        for t in range(9):
            new = neutral_propagation_step(
                state, cfg=cfg, mass_kg=ckpt.mass_kg,
                droplet_radii=ckpt.droplet_radii,
                prev_distance_angstrom=prev_dist,
                rng=np.random.default_rng(t),
            )
            assert np.all(new.L_droplet_eV_ps >= prev_L - 1e-12)
            prev_dist = _step_distance(state, new)
            prev_L = new.L_droplet_eV_ps
            state = new


class TestReproducibility:
    def test_same_rng_same_trajectory(self):
        cfg, ckpt, state0 = _bootstrap(num_molecules=5, num_steps=10,
                                        R0=2.666, single_pulse=False)

        def run():
            state = state0
            prev_dist = None
            xs = []
            for t in range(9):
                new = neutral_propagation_step(
                    state, cfg=cfg, mass_kg=ckpt.mass_kg,
                    droplet_radii=ckpt.droplet_radii,
                    prev_distance_angstrom=prev_dist,
                    rng=np.random.default_rng(t + 1000),
                )
                prev_dist = _step_distance(state, new)
                state = new
                xs.append(state.x.copy())
            return xs

        a_xs = run()
        b_xs = run()
        for ax, bx in zip(a_xs, b_xs):
            np.testing.assert_array_equal(ax, bx)


class TestEnergyBookkeeping:
    def test_E_kin_E_pot_E_dissip_balanced(self):
        cfg, ckpt, state0 = _bootstrap(num_molecules=20, num_steps=50,
                                        R0=2.666, single_pulse=False)
        E0 = (state0.E_kin_eV + state0.E_pot_eV + state0.E_dissip_eV).sum()
        state = state0
        prev_dist = None
        for t in range(49):
            new = neutral_propagation_step(
                state, cfg=cfg, mass_kg=ckpt.mass_kg,
                droplet_radii=ckpt.droplet_radii,
                prev_distance_angstrom=prev_dist,
                rng=np.random.default_rng(t),
            )
            prev_dist = _step_distance(state, new)
            state = new
        E_end = (state.E_kin_eV + state.E_pot_eV + state.E_dissip_eV).sum()
        assert abs(E_end - E0) / abs(E0) < 0.05


class TestStateHelpers:
    def test_extract_then_write_round_trip(self):
        cfg, ckpt, _ = _bootstrap(num_molecules=5, num_steps=3)
        x_before = ckpt.positions_x[:, 0].copy()
        s = state_from_checkpoint_column(ckpt, 0)
        write_state_to_checkpoint_column(s, ckpt, 1)
        np.testing.assert_array_equal(ckpt.positions_x[:, 1], x_before)

    def test_state_extraction_copies_arrays(self):
        cfg, ckpt, _ = _bootstrap(num_molecules=3, num_steps=3)
        s = state_from_checkpoint_column(ckpt, 0)
        s.x[0] = 999.0
        assert ckpt.positions_x[0, 0] != 999.0
