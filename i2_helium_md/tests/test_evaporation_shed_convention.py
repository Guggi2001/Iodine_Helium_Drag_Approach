"""OQ-J shed-convention enum (Tier-2 plan SI.11.2 item 1 adjudication, 2026-07-17).

``evaporation_shed_convention`` selects the velocity/momentum convention of the
generative evaporation channel:

* ``"cold"`` (byte-inert default) -- the delivered momentum-conserving reset:
  the shed He is left at rest in the lab frame, the complex keeps its full
  momentum (``v -> m/(m - m_He) * v``), and the ledger books the (negative)
  reduced-mass KE injection. Retained as the diagnostic bound arm.
* ``"co_moving"`` -- the physical zeroth order for thermal evaporation from a
  moving complex (findings §4n / I59): the He leaves co-moving, the velocity
  is unchanged, and the ledger books the (positive) ``+0.5*m_He*|v|^2`` the He
  carries away.

The convention threads through all three shed surfaces from the ONE config
field: ``evaporation_step`` / ``evaporation_step_components`` (ion stage +,
via ``biphasic_step`` under the relaxation cfg view, the E2 stage) and the
detection-stage event loop. The (E_int, n) jump chain reads neither v nor m,
so under a fixed seed the shed *sequence* is convention-invariant -- only the
velocity/ledger bookkeeping moves. Exactness notes: the co-moving ledger term
``+0.5*m_He*|v|^2`` is exact (the He leaves at exactly ``v``), so closure
asserts are machine-precision.
"""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np
import pytest

from i2_helium_md.config import SimConfig, check_evaporation_config
from i2_helium_md.physics.constants import MASS_HE_AMU, U
from i2_helium_md.physics.evaporation import (
    evaporation_step,
    evaporation_step_components,
)
from i2_helium_md.physics.mass_jump import (
    continuous_velocity_shed,
    continuous_velocity_shed_components,
)
from i2_helium_md.simulation.detection_stage import run_detection_stage
from i2_helium_md.simulation.ion_propagation_step import (
    _amu_ang2_ps2_to_eV,
    biphasic_step,
)
from i2_helium_md.simulation.run_directory import RunDirectory

from scripts.tier2_common import build_biphasic_cfg
from tests.test_biphasic_step import _biphasic_cfg, _biphasic_state
from tests.test_detection_stage import _detect_cfg, _far_seed

# Force-fire evaporation: k*dt >> 1 at the in-band (n=21, E_int=0.15 eV) point.
_NU_FORCE = 1.0e5


# ---------------------------------------------------------------------------
# Config surface: default, guard, back-compat, round-trip
# ---------------------------------------------------------------------------
class TestConfigSurface:
    def test_default_is_cold_and_validates(self):
        cfg = SimConfig()
        assert cfg.evaporation_shed_convention == "cold"
        cfg.validate()

    def test_unknown_convention_rejected(self):
        cfg = SimConfig(evaporation_shed_convention="warm")
        with pytest.raises(ValueError, match="evaporation_shed_convention"):
            check_evaporation_config(cfg)
        with pytest.raises(ValueError, match="evaporation_shed_convention"):
            cfg.validate()

    def test_pre_enum_cfg_json_loads_with_default(self, tmp_path):
        # Slice-DS/T7 precedent: a cfg.json written before the field existed
        # must load with the byte-inert default.
        run = RunDirectory(tmp_path / "9A_drag_tier2probe_backcompat_shed")
        run.save_cfg(SimConfig())
        payload = json.loads(run.cfg_path.read_text(encoding="utf-8"))
        assert payload.pop("evaporation_shed_convention") == "cold"
        run.cfg_path.write_text(json.dumps(payload), encoding="utf-8")
        assert run.load_cfg().evaporation_shed_convention == "cold"

    def test_cfg_json_roundtrip(self, tmp_path):
        run = RunDirectory(tmp_path / "9A_drag_tier2probe_shed_roundtrip")
        run.save_cfg(SimConfig(evaporation_shed_convention="co_moving"))
        assert run.load_cfg().evaporation_shed_convention == "co_moving"


# ---------------------------------------------------------------------------
# mass_jump: per-ion generalization of the co-moving operator
# ---------------------------------------------------------------------------
class TestContinuousShedComponents:
    def test_hand_oracle_and_positive_ledger(self):
        vx = np.array([3.0, -1.0])
        vy = np.array([0.0, 2.0])
        vz = np.array([4.0, 0.0])
        m = np.array([210.9546, 130.9071])  # per-ion masses (fires diverged)
        vx_p, vy_p, vz_p, m_p, dE = continuous_velocity_shed_components(
            vx, vy, vz, m, m_he_amu=MASS_HE_AMU,
        )
        np.testing.assert_array_equal(vx_p, vx)  # velocity strictly unchanged
        np.testing.assert_array_equal(vy_p, vy)
        np.testing.assert_array_equal(vz_p, vz)
        np.testing.assert_allclose(m_p, m - MASS_HE_AMU)
        expected = 0.5 * MASS_HE_AMU * (vx**2 + vy**2 + vz**2)
        np.testing.assert_allclose(dE, expected)
        assert np.all(dE > 0.0)  # KE carried away, booked positive

    def test_parity_with_single_atom_oracle(self):
        rng = np.random.default_rng(3)
        vx, vy, vz = (rng.normal(size=5) for _ in range(3))
        m = 130.0 + 40.0 * rng.random(size=5)
        vx_p, vy_p, vz_p, m_p, dE = continuous_velocity_shed_components(
            vx, vy, vz, m, m_he_amu=MASS_HE_AMU,
        )
        for i in range(5):
            res = continuous_velocity_shed(
                np.array([vx[i], vy[i], vz[i]]), float(m[i]),
                m_he_amu=MASS_HE_AMU,
            )
            np.testing.assert_allclose(
                [vx_p[i], vy_p[i], vz_p[i]], res.v_plus
            )
            assert m_p[i] == pytest.approx(res.m_plus_amu)
            assert dE[i] == pytest.approx(res.dE_mass_transfer)

    def test_underweight_ion_fails_loud(self):
        v = np.zeros(2)
        with pytest.raises(ValueError, match="pre-shed mass"):
            continuous_velocity_shed_components(
                v, v, v, np.array([210.0, 3.0]), m_he_amu=MASS_HE_AMU,
            )


# ---------------------------------------------------------------------------
# The evaporation channel: convention dispatch, RNG contract, ledger signs
# ---------------------------------------------------------------------------
def _fire_kwargs(**over):
    kw = dict(
        E_int_eV=0.15, n=21, nu=_NU_FORCE, kappa=1.0, dt_ps=0.01,
    )
    kw.update(over)
    return kw


class TestEvaporationStepConvention:
    V = np.array([3.0, 0.0, 4.0])  # |v|^2 = 25
    M = 210.9546

    def test_co_moving_fire_keeps_velocity_books_positive(self):
        res = evaporation_step(
            rng=np.random.default_rng(0), v=self.V, m_amu=self.M,
            shed_convention="co_moving", **_fire_kwargs(),
        )
        assert res.fired
        np.testing.assert_array_equal(res.v_plus, self.V)
        assert res.m_plus_amu == pytest.approx(self.M - MASS_HE_AMU)
        assert res.dE_mass_transfer == pytest.approx(0.5 * MASS_HE_AMU * 25.0)

    def test_cold_explicit_equals_omitted(self):
        res_default = evaporation_step(
            rng=np.random.default_rng(0), v=self.V, m_amu=self.M,
            **_fire_kwargs(),
        )
        res_cold = evaporation_step(
            rng=np.random.default_rng(0), v=self.V, m_amu=self.M,
            shed_convention="cold", **_fire_kwargs(),
        )
        assert res_default.fired and res_cold.fired
        np.testing.assert_array_equal(res_default.v_plus, res_cold.v_plus)
        assert res_default.dE_mass_transfer == res_cold.dE_mass_transfer
        # and the cold kick is the delivered momentum-conserving reset
        kick = self.M / (self.M - MASS_HE_AMU)
        np.testing.assert_allclose(res_cold.v_plus, kick * self.V)
        assert res_cold.dE_mass_transfer < 0.0

    def test_unknown_convention_raises(self):
        with pytest.raises(ValueError, match="shed_convention"):
            evaporation_step(
                rng=np.random.default_rng(0), v=self.V, m_amu=self.M,
                shed_convention="warm", **_fire_kwargs(),
            )

    def test_components_same_seed_same_fires_opposite_ledgers(self):
        rngs = [np.random.default_rng(11) for _ in range(2)]
        vx = np.array([3.0, 3.0, 3.0])
        vz = np.zeros(3)
        n = np.array([21.0, 21.0, 21.0])
        m = np.full(3, self.M)
        e = np.array([0.15, 0.15, 0.15])
        out = {}
        for conv, rng in zip(("cold", "co_moving"), rngs):
            out[conv] = evaporation_step_components(
                rng=rng, E_int_eV=e, n=n, vx=vx, vy=vz, vz=vz, m_amu=m,
                shed_convention=conv, **{k: v for k, v in _fire_kwargs().items()
                                         if k not in ("E_int_eV", "n")},
            )
        n_c, m_c, vx_c, _, _, dEi_c, dmt_c, fired_c = out["cold"]
        n_m, m_m, vx_m, _, _, dEi_m, dmt_m, fired_m = out["co_moving"]
        # RNG contract: identical draw stream -> identical fire pattern and
        # identical (E_int, n) chain increments; convention moves only v / ledger.
        np.testing.assert_array_equal(fired_c, fired_m)
        np.testing.assert_array_equal(n_c, n_m)
        np.testing.assert_array_equal(dEi_c, dEi_m)
        np.testing.assert_allclose(m_c, m_m)
        assert fired_m.all()
        np.testing.assert_array_equal(vx_m, vx)          # co-moving: unchanged
        np.testing.assert_allclose(vx_c, vx * self.M / (self.M - MASS_HE_AMU))
        # exact sign/magnitude relation: cold books the reduced-mass injection
        # -(m/(m - m_He)) times the co-moving carried KE.
        np.testing.assert_allclose(
            dmt_c, -dmt_m * self.M / (self.M - MASS_HE_AMU)
        )


# ---------------------------------------------------------------------------
# Driver threading: biphasic_step reads cfg.evaporation_shed_convention
# ---------------------------------------------------------------------------
class TestBiphasicStepThreading:
    def _run(self, convention):
        cfg = _biphasic_cfg(
            evap_rate_prefactor_per_ps=_NU_FORCE,
            evaporation_shed_convention=convention,
        )
        state = _biphasic_state(e_int=0.15)  # in the shed band -> fires at nu=1e5
        two_N = state.x.shape[0]
        return state, cfg, biphasic_step(
            state, rng=np.random.default_rng(5), cfg=cfg,
            droplet_radii=np.full(two_N, 30.0), gate_steepness=14.2,
        )

    def test_co_moving_keeps_velocity_and_books_positive(self):
        state, _, new = self._run("co_moving")
        np.testing.assert_array_equal(new.n_shell, state.n_shell - 1.0)
        np.testing.assert_array_equal(new.vx, state.vx)
        np.testing.assert_array_equal(new.vy, state.vy)
        np.testing.assert_array_equal(new.vz, state.vz)
        np.testing.assert_allclose(
            new.mass_kg, state.mass_kg - MASS_HE_AMU * U
        )
        speed_sq = state.vx**2 + state.vy**2 + state.vz**2
        np.testing.assert_allclose(
            new.E_mass_transfer_eV,
            _amu_ang2_ps2_to_eV(0.5 * MASS_HE_AMU * speed_sq),
        )

    def test_conventions_share_the_jump_chain(self):
        # Same seed -> identical (n, E_int) chains; only v / E_mass_transfer move.
        state, _, cold = self._run("cold")
        _, _, com = self._run("co_moving")
        np.testing.assert_array_equal(cold.n_shell, com.n_shell)
        np.testing.assert_allclose(cold.E_int_eV, com.E_int_eV)
        m_amu = state.mass_kg[0] / U
        np.testing.assert_allclose(
            cold.vx, com.vx * m_amu / (m_amu - MASS_HE_AMU)
        )
        assert np.all(cold.E_mass_transfer_eV < 0.0)
        assert np.all(com.E_mass_transfer_eV > 0.0)


# ---------------------------------------------------------------------------
# Detection stage: the event loop rides the same config field
# ---------------------------------------------------------------------------
class TestDetectionStageConvention:
    def _detect(self, convention):
        # n=3 / E_int=0.02 eV: two fires (3->2->1), then frozen below D_0(1);
        # every event is far outside the droplet (free flight).
        cfg = _detect_cfg(evaporation_shed_convention=convention)
        seed = _far_seed(cfg, n_shell=3, E_int_eV=0.02, num_molecules=2)
        return seed, run_detection_stage(seed, cfg)

    def test_co_moving_keeps_velocity_through_events(self):
        seed, res = self._detect("co_moving")
        assert res.event_time_ps.size > 0  # the fixture must actually fire
        np.testing.assert_array_equal(res.vx_detected, seed.velocities_x[:, -1])
        np.testing.assert_array_equal(res.vy_detected, seed.velocities_y[:, -1])
        np.testing.assert_array_equal(res.vz_detected, seed.velocities_z[:, -1])
        assert np.all(res.event_dE_mass_transfer_eV > 0.0)

    def test_cold_kicks_velocity_up(self):
        seed, res = self._detect("cold")
        assert res.event_time_ps.size > 0
        shed_ions = res.n_detected < seed.n_shell[:, -1]
        assert shed_ions.any()
        v0 = np.abs(seed.velocities_x[:, -1][shed_ions])
        assert np.all(np.abs(res.vx_detected[shed_ions]) > v0)
        assert np.all(res.event_dE_mass_transfer_eV < 0.0)

    def test_five_term_closure_holds_co_moving(self):
        # Per ion: E_kin + E_pot + E_dissip + E_mass_transfer + E_int is
        # conserved across the stage under the co-moving convention (the
        # +0.5*m_He*|v|^2 ledger term is exact, so machine precision).
        seed, res = self._detect("co_moving")
        total_seed = (
            seed.E_kin_eV[:, -1] + seed.E_pot_eV[:, -1]
            + seed.E_dissip_eV[:, -1] + seed.E_mass_transfer_eV[:, -1]
            + seed.E_int_eV[:, -1]
        )
        total_det = (
            res.E_kin_detected_eV + res.E_pot_detected_eV
            + res.E_dissip_detected_eV + res.E_mass_transfer_detected_eV
            + res.E_int_detected_eV
        )
        np.testing.assert_allclose(total_det, total_seed, rtol=0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# Script surface: build_biphasic_cfg None-sentinel pass-through
# ---------------------------------------------------------------------------
class TestBuildBiphasicCfgKwarg:
    _ARGS = ("9A", "shared_pure_cubic")
    _COMMON = dict(num_molecules=2, ion_time_ps=0.02, dt_ion_ps=0.01, seed=123)

    def test_none_sentinel_rides_default(self):
        cfg = build_biphasic_cfg(*self._ARGS, **self._COMMON)
        assert cfg.evaporation_shed_convention == "cold"

    def test_explicit_value_stamped(self):
        cfg = build_biphasic_cfg(
            *self._ARGS, evaporation_shed_convention="co_moving", **self._COMMON,
        )
        assert cfg.evaporation_shed_convention == "co_moving"
