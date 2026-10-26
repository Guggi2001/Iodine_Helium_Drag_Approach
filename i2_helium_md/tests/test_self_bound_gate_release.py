"""G13 boil-off of suppressed ions in free flight (THESIS_FINAL_GOALS G13.5-G13.7).

``self_bound_gate_release="decoupled"`` lifts the self-bound gate wherever
the ion is outside the nominal droplet surface (depth > 0): a net
self-unbound complex runs the delivered RRK cascade (epsilon = 0) to bare I+.
Tested on tiny synthetic seeds: the endpoint and kinematics (n = 0, speed
unchanged under co_moving), the 5-term ledger, the time-exhausted lane, the
rate switch, the Stage-I wiring, the detection-stage main loop, and the
config-load surface.
"""

from __future__ import annotations

from dataclasses import fields, replace

import numpy as np
import pytest

from i2_helium_md.config import check_detection_config
from i2_helium_md.physics.constants import MASS_HE_AMU, MASS_I_ION_AMU, U
from i2_helium_md.physics.dissociation_ladder import d0_of_n, ladder_cumsum
from i2_helium_md.physics.evaporation import effective_dof, rrk_rate, shed_probability
from i2_helium_md.simulation.ion_propagation_step import biphasic_step
from i2_helium_md.simulation.detection_stage import DetectionResult, run_detection_stage
from i2_helium_md.simulation.ion_propagation_step import _E_kin_eV

from tests.test_biphasic_step import _biphasic_cfg, _biphasic_state
from tests.test_detection_coulomb_closure import _five_term
from tests.test_detection_stage import _detect_cfg, _far_seed, _ladder


def _cfg(policy="decoupled", **kw):
    base = dict(seed=11, evaporation_shed_convention="co_moving",
                self_bound_gate_release=policy)
    base.update(kw)
    return _detect_cfg(**base)


def _sigma21(cfg) -> float:
    picture, kappa = _ladder(cfg)
    return float(ladder_cumsum(21, picture=picture, kappa=kappa))


def _mixed_seed(cfg):
    """N = 2 ejected ions at n = 21: ions 0, 2 suppressed (G = +2 % of
    Sigma(21)), ions 1, 3 in-band (G = -2 %), so both lanes are live."""
    s21 = _sigma21(cfg)
    seed = _far_seed(cfg, n_shell=21, E_int_eV=1.02 * s21, speed=8.0)
    seed.E_int_eV[[1, 3], 0] = 0.98 * s21
    return seed


def _speed(res):
    return np.sqrt(res.vx_detected ** 2 + res.vy_detected ** 2 + res.vz_detected ** 2)


def _assert_results_equal(a: DetectionResult, b: DetectionResult) -> None:
    for f in fields(DetectionResult):
        np.testing.assert_array_equal(np.asarray(getattr(a, f.name)),
                                      np.asarray(getattr(b, f.name)), err_msg=f.name)


class TestBoilOffEndpoint:
    def test_suppressed_ions_boil_off_to_bare_at_unchanged_speed(self):
        cfg = _cfg(evap_rrk_dof=2.0)
        s21 = _sigma21(cfg)
        seed = _far_seed(cfg, n_shell=21, E_int_eV=1.02 * s21, speed=8.0)
        res = run_detection_stage(seed, cfg)
        assert list(res.state_reason) == ["frozen"] * 4
        np.testing.assert_array_equal(res.n_detected, np.zeros(4))
        np.testing.assert_array_equal(np.diff(res.event_offsets), np.full(4, 21))
        np.testing.assert_allclose(res.mass_detected_kg, MASS_I_ION_AMU * U, rtol=1e-12)
        # co_moving: v never changes, so KE is the handover speed at the bare mass.
        np.testing.assert_allclose(_speed(res), 8.0, rtol=0, atol=1e-12)
        np.testing.assert_allclose(
            res.E_kin_detected_eV,
            _E_kin_eV(np.full(4, MASS_I_ION_AMU * U), np.full(4, 8.0),
                      np.zeros(4), np.zeros(4)), rtol=1e-12)
        # epsilon = 0: every shed drains exactly D_0, so E_int ends at G.
        np.testing.assert_allclose(res.E_int_detected_eV, 0.02 * s21, rtol=1e-9)
        t = res.event_time_ps
        assert np.all((t > res.t_handover_ps) & (t <= res.detection_time_ps))
        for i in range(4):
            sl = slice(res.event_offsets[i], res.event_offsets[i + 1])
            assert np.all(np.diff(t[sl]) >= 0)
            np.testing.assert_array_equal(res.event_pre_shed_n[sl], np.arange(21, 0, -1))

    @pytest.mark.parametrize("shed", ["co_moving", "cold"])
    def test_five_term_invariant_closes(self, shed):
        cfg = _cfg(evap_rrk_dof=2.0, evaporation_shed_convention=shed)
        seed = _mixed_seed(cfg)
        res = run_detection_stage(seed, cfg)
        before, after = _five_term(seed, res)
        np.testing.assert_allclose(after, before, rtol=0, atol=1e-11)

    def test_flight_too_short_leaves_time_exhausted(self):
        # nu ~ 1e-12 /ps: no wait fits the flight, so the model -- not a
        # convention -- leaves the complex intact (time_exhausted at n = 21).
        cfg = _cfg(evap_rate_prefactor_per_ps=1e-12)
        seed = _far_seed(cfg, n_shell=21, E_int_eV=1.02 * _sigma21(cfg))
        res = run_detection_stage(seed, cfg)
        assert list(res.state_reason) == ["time_exhausted"] * 4
        np.testing.assert_array_equal(res.n_detected, np.full(4, 21.0))
        assert res.event_time_ps.size == 0


class TestIsolation:
    def test_intact_default_keeps_suppressed(self):
        cfg = _cfg(policy="never", evap_rrk_dof=2.0)
        res = run_detection_stage(_mixed_seed(cfg), cfg)
        assert list(res.state_reason[[0, 2]]) == ["suppressed"] * 2


class TestRefusals:
    def test_ungated_cooling_still_refused_by_the_handover_guard(self):
        cfg = _cfg(cooling_spatial_gate="none")
        seed = _far_seed(cfg, n_shell=21, E_int_eV=1.5 * _sigma21(cfg))
        with pytest.raises(ValueError, match="P1-P3"):
            run_detection_stage(seed, cfg)


class TestConfig:
    def test_default_is_intact(self):
        assert _detect_cfg().self_bound_gate_release == "never"

    def test_unknown_value_rejected(self):
        with pytest.raises(ValueError, match="self_bound_gate_release"):
            check_detection_config(replace(_detect_cfg(),
                                           self_bound_gate_release="shatter"))

    def test_release_with_gate_onset_override_refused(self):
        with pytest.raises(ValueError, match="gate_onset_override_eV"):
            check_detection_config(replace(
                _detect_cfg(), self_bound_gate_release="decoupled",
                gate_onset_override_eV=0.3))


# ---------------------------------------------------------------------------
# decoupled: the rate switch, Stage-I wiring, detection-stage main loop
# ---------------------------------------------------------------------------
class TestRateSwitch:
    def _kw(self, cfg):
        picture, kappa = _ladder(cfg)
        return dict(nu=cfg.evap_rate_prefactor_per_ps, picture=picture, kappa=kappa)

    def test_lifted_gate_gives_the_rrk_bracket_for_a_self_unbound_ion(self):
        cfg = _detect_cfg()
        kw = self._kw(cfg)
        s21 = _sigma21(cfg)
        E = 1.5 * s21
        assert rrk_rate(E, 21, **kw) == 0.0                         # gated
        d0 = float(d0_of_n(21, picture=kw["picture"], kappa=kw["kappa"]))
        s = float(effective_dof(21))
        np.testing.assert_allclose(rrk_rate(E, 21, gate_lifted=True, **kw),
                                   kw["nu"] * (1.0 - d0 / E) ** (s - 1.0), rtol=1e-12)

    def test_lifted_gate_changes_nothing_in_band_or_at_n_le_1(self):
        cfg = _detect_cfg()
        kw = self._kw(cfg)
        s21 = _sigma21(cfg)
        E = np.array([0.98 * s21, 0.5, 0.5, 0.0])
        n = np.array([21, 1, 0, 5])
        np.testing.assert_array_equal(rrk_rate(E, n, gate_lifted=True, **kw),
                                      rrk_rate(E, n, **kw))

    def test_per_ion_mask(self):
        cfg = _detect_cfg()
        kw = self._kw(cfg)
        E = np.full(4, 1.5 * _sigma21(cfg))
        k = rrk_rate(E, np.full(4, 21), gate_lifted=np.array([0, 1, 0, 1], bool), **kw)
        assert k[0] == k[2] == 0.0 and k[1] == k[3] > 0.0


def _stage_one(release, *, steps=1, seed=3):
    """N = 3 biphasic state at E_int 0.24 eV > Sigma(21) (net self-unbound);
    ions 0-2 deep inside a 30 A droplet, ions 3-5 outside (R = 0.1 A < |x|).
    Cooling ungated and pickup off, so only the evaporation gate differs."""
    cfg = _biphasic_cfg(cooling_spatial_gate="none", self_bound_gate_release=release,
                        internal_energy_cooling_tau_ps=1e9)
    state = _biphasic_state(e_int=0.24)
    droplet = np.array([30.0, 30.0, 30.0, 0.1, 0.1, 0.1])
    rng = np.random.default_rng(seed)
    for _ in range(steps):
        state = biphasic_step(state, rng=rng, cfg=cfg, droplet_radii=droplet,
                              gate_steepness=14.2)
    return cfg, state


class TestStageIDecoupled:
    def test_one_step_outside_ions_on_the_lifted_rate(self):
        cfg_n, never = _stage_one("never")
        cfg_d, dec = _stage_one("decoupled")
        assert np.all(never.n_shell == 21)                           # gated everywhere
        inside, out = np.arange(3), np.arange(3, 6)
        for f in ("n_shell", "E_int_eV", "mass_kg", "vx", "E_dissip_eV"):
            np.testing.assert_array_equal(getattr(dec, f)[inside], getattr(never, f)[inside])
        # never fires nothing, so its E_int is the post-cooling input of the draw
        picture, kappa = _ladder(cfg_d)
        k = rrk_rate(never.E_int_eV[out], np.full(3, 21), nu=cfg_d.evap_rate_prefactor_per_ps,
                     picture=picture, kappa=kappa, gate_lifted=True)
        assert np.all(k > 0.0)
        draws = np.random.default_rng(3).random(6)[out]   # first draw = evaporation
        fired = draws < shed_probability(k, cfg_d.dt_ion)
        np.testing.assert_array_equal(dec.n_shell[out], 21 - fired)

    def test_outside_ions_shed_inside_ions_stay_gated(self):
        _, never = _stage_one("never", steps=500)
        _, dec = _stage_one("decoupled", steps=500)
        assert np.all(never.n_shell == 21)
        assert np.all(dec.n_shell[:3] == 21)
        assert dec.n_shell[3:].sum() < 3 * 21                        # sheds happened
        np.testing.assert_allclose(
            dec.mass_kg / U, MASS_I_ION_AMU + dec.n_shell * MASS_HE_AMU, rtol=0, atol=1e-6)


class TestDetectionDecoupled:
    def test_free_suppressed_ions_boil_off_in_the_main_loop(self):
        cfg = _cfg(policy="decoupled", evap_rrk_dof=2.0)
        seed = _mixed_seed(cfg)
        res = run_detection_stage(seed, cfg)
        assert "suppressed" not in set(np.asarray(res.state_reason).astype(str))
        np.testing.assert_array_equal(res.n_detected[[0, 2]], [0.0, 0.0])
        np.testing.assert_allclose(_speed(res)[[0, 2]], 8.0, rtol=0, atol=1e-12)
        before, after = _five_term(seed, res)
        np.testing.assert_allclose(after, before, rtol=0, atol=1e-11)

    def test_in_band_ions_identical_to_never(self):
        # Lifting the gate does not change an in-band rate, and the main loop
        # draws in ion order -- so ions before the first suppressed one agree.
        cfg_d = _cfg(policy="decoupled", evap_rrk_dof=2.0)
        s21 = _sigma21(cfg_d)
        seed = _far_seed(cfg_d, n_shell=21, E_int_eV=0.98 * s21, speed=8.0)
        a = run_detection_stage(seed, replace(cfg_d, self_bound_gate_release="never"))
        b = run_detection_stage(seed, cfg_d)
        _assert_results_equal(a, b)


class TestReleaseConfig:
    def test_release_requires_biphasic(self):
        from i2_helium_md.config import check_self_bound_gate_release
        cfg = replace(_detect_cfg(), self_bound_gate_release="decoupled",
                      mass_scenario="fixed")
        with pytest.raises(ValueError, match="biphasic"):
            check_self_bound_gate_release(cfg)

    def test_decoupled_does_not_need_the_detection_stage(self):
        check_detection_config(replace(_detect_cfg(), detection_stage_enabled=False,
                                       self_bound_gate_release="decoupled"))
