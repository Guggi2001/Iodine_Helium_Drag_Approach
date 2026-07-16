"""Tests for simulation/detection_stage.py (Tier-2 Slice DS).

Follows the design's seven-step validation plan (TIER2_DETECTION_STAGE_DESIGN.md
§5): the seed-fixed waiting-time formula, the detection.npz round-trip + ragged
offsets, the k = 0 permanent lanes (zero events, zero draws), the multi-event
occupancy against the **hypoexponential closed form** (the analytic oracle --
computed in-test from the delivered ``rrk_rate`` chain, never a production
arm), the 5-term closure across the stage boundary, the statistical
equivalence against a fixed-dt E2 extension over a short horizon, and a
tiny-N driver/report smoke. Plus the Slice-DS additions: the skip-path
seeding lanes (design §1 item 5) and the Wave-7 config back-compat criterion
(probe plan Addendum C.2). Tiny synthetic checkpoints + one tiny real
biphasic run; no figures.
"""

from __future__ import annotations

import json
import math
from dataclasses import replace

import numpy as np
import pytest

from i2_helium_md.config import check_detection_config
from i2_helium_md.physics.constants import MASS_HE_AMU, MASS_I_ION_AMU, U
from i2_helium_md.physics.dissociation_ladder import d0_of_n, ladder_cumsum
from i2_helium_md.physics.evaporation import rrk_rate
from i2_helium_md.postprocess.size_distribution import (
    compute_terminal_shell_distribution,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint, load_ion_checkpoint
from i2_helium_md.simulation.detection_stage import (
    DETECTION_STREAM_KEY,
    EPS_DRAIN,
    DetectionResult,
    load_detection_result,
    run_detection_stage,
    save_detection_result,
)
from i2_helium_md.simulation.ion import run_ion_propagation
from i2_helium_md.simulation.ion_propagation_step import _E_kin_eV
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage
from i2_helium_md.simulation.run_directory import RunDirectory

from tests.test_biphasic_step import _biphasic_driver_cfg, _tiny_neutral


def _detect_cfg(**overrides):
    """A biphasic driver cfg with the detection stage enabled.

    Defaults to the gated production arm (``density_scaled``) because the P3
    handover guard exempts only permanent-state ions under ``none`` (ungated
    cooling acts everywhere): every live-ion fixture needs the gate + an
    ejected seed. ``detection_time_ps`` defaults to the Sourced 8.53e6 ps.
    """
    base = dict(
        detection_stage_enabled=True,
        detection_time_ps=8.53e6,
        cooling_spatial_gate="density_scaled",
    )
    base.update(overrides)
    return _biphasic_driver_cfg(**base)


def _far_seed(cfg, *, n_shell=21, E_int_eV=0.15, num_molecules=2,
              time_ps_start=20.0, speed=0.2, inside=False):
    """A one-column v7 biphasic seed with **self-consistent E_kin**.

    Ions sit at |z| ~ 1000 A, far outside the 30 A droplet (erfc density
    underflows to exactly 0 -> the P3 guard passes for live ions under the
    gated arm), unless ``inside=True`` (all ions at the droplet centre --
    the P3 violator fixture). ``E_kin_eV`` is computed from ``(m, v)`` by the
    delivered idiom so the 5-term sum is checkable across the stage boundary;
    ``E_pot``/``E_dissip``/``E_mass_transfer`` start at 0.
    """
    N = int(num_molecules)
    two_n = 2 * N
    T = 1
    m_amu = MASS_I_ION_AMU + n_shell * MASS_HE_AMU
    mass_kg = np.full(two_n, m_amu * U)

    z_mag = 0.0 if inside else 1000.0
    x = np.tile(np.array([1.3, -1.3]), N).reshape(two_n, 1)
    y = np.zeros((two_n, T))
    z = np.full((two_n, T), z_mag)

    vx = np.full((two_n, T), float(speed))
    vy = np.zeros((two_n, T))
    vz = np.zeros((two_n, T))
    e_kin = _E_kin_eV(mass_kg, vx[:, 0], vy[:, 0], vz[:, 0]).reshape(two_n, 1)

    def col(val):
        return np.full((two_n, T), float(val))

    return IonCheckpoint(
        num_molecules=N,
        time_ps=np.full(T, float(time_ps_start)),
        positions_x=x, positions_y=y, positions_z=z,
        velocities_x=vx, velocities_y=vy, velocities_z=vz,
        positions_final_x=np.zeros(two_n), positions_final_y=np.zeros(two_n),
        positions_final_z=np.zeros(two_n),
        velocities_final_x=np.zeros(two_n), velocities_final_y=np.zeros(two_n),
        velocities_final_z=np.zeros(two_n),
        mass_kg=mass_kg.copy(), mass_final_kg=mass_kg.copy(),
        mass_history_kg=np.tile(mass_kg.reshape(two_n, 1), (1, T)),
        droplet_radii_angstrom=np.full(two_n, 30.0),
        E_kin_eV=e_kin, E_pot_eV=col(0.0), E_dissip_eV=np.zeros((two_n, T)),
        E_mass_transfer_eV=np.zeros((two_n, T)), E_int_eV=col(E_int_eV),
        n_shell=col(n_shell),
        b_ion_outside=np.zeros(N, dtype=bool),
        relative_loss_per_ps=np.zeros((two_n, T)),
        number_of_collisions=np.zeros((two_n, T), dtype=int),
        temperature_diagnostic=np.full((T, 3), np.nan),
        mass_scenario="biphasic",
    )


def _ladder(cfg):
    return cfg.ladder_electronic_picture, float(cfg.ladder_steepness)


# ---------------------------------------------------------------------------
# Step 1 -- waiting-time formula + RNG contract
# ---------------------------------------------------------------------------
class TestFormulaAndRNG:
    def test_first_event_time_matches_hand_formula(self):
        # n = 1 direct dissociation: k = nu exactly, one fire -> n = 0 frozen.
        # Each ion's single event time is t_h - ln(u_i)/nu with u_i the i-th
        # sequential draw of the stage stream (ion-major order) -- the design
        # §5 step-1 hand oracle.
        cfg = _detect_cfg(seed=0)
        picture, kappa = _ladder(cfg)
        d0_1 = float(d0_of_n(1, picture=picture, kappa=kappa))
        seed = _far_seed(cfg, n_shell=1, E_int_eV=5.0 * d0_1)

        result = run_detection_stage(seed, cfg)

        nu = cfg.evap_rate_prefactor_per_ps
        stream = np.random.default_rng(
            np.random.SeedSequence((0, DETECTION_STREAM_KEY))
        )
        expected = 20.0 - np.log(stream.random(size=4)) / nu
        np.testing.assert_allclose(result.event_time_ps, expected,
                                   rtol=0, atol=1e-12)
        np.testing.assert_array_equal(result.n_detected, np.zeros(4))
        assert list(result.state_reason) == ["frozen"] * 4
        np.testing.assert_array_equal(result.event_pre_shed_n, np.ones(4))
        np.testing.assert_allclose(result.event_dE_int_eV, -d0_1,
                                   rtol=0, atol=1e-15)
        np.testing.assert_allclose(result.event_dE_bind_fold_eV, d0_1,
                                   rtol=0, atol=1e-15)

    def test_default_rng_derivation_pinned(self):
        # rng=None MUST derive SeedSequence((cfg.seed, DETECTION_STREAM_KEY)) --
        # seeding from cfg.seed directly would correlate with the ion stream.
        cfg = _detect_cfg(seed=777)
        picture, kappa = _ladder(cfg)
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        cfg = replace(cfg, evap_rrk_dof=2.0)
        seed = _far_seed(cfg, n_shell=21, E_int_eV=0.98 * sigma21)

        r_default = run_detection_stage(seed, cfg)
        r_explicit = run_detection_stage(
            seed, cfg,
            rng=np.random.default_rng(
                np.random.SeedSequence((777, DETECTION_STREAM_KEY))
            ),
        )
        assert result_events_nonzero(r_default)   # draws mattered (non-vacuous)
        np.testing.assert_array_equal(r_default.n_detected, r_explicit.n_detected)
        np.testing.assert_array_equal(r_default.event_time_ps,
                                      r_explicit.event_time_ps)

    def test_zero_draws_on_permanent_seed(self):
        # k = 0 lanes consume NO draws (design §5 step 3): the stage stream
        # state is untouched by a frozen input.
        cfg = _detect_cfg()
        picture, kappa = _ladder(cfg)
        d0_21 = float(d0_of_n(21, picture=picture, kappa=kappa))
        seed = _far_seed(cfg, n_shell=21, E_int_eV=0.5 * d0_21)

        rng = np.random.default_rng(123)
        state_before = rng.bit_generator.state
        run_detection_stage(seed, cfg, rng=rng)
        assert rng.bit_generator.state == state_before

    def test_stream_key_is_stable(self):
        # An accidental key change would silently re-seed every detection run.
        assert DETECTION_STREAM_KEY == 0xD5_2026


def result_events_nonzero(result: DetectionResult) -> bool:
    return result.event_time_ps.size > 0


# ---------------------------------------------------------------------------
# Step 3 -- permanent lanes (k = 0): zero events, reason recorded, n preserved
# ---------------------------------------------------------------------------
class TestPermanentLanes:
    def _run_lane(self, *, n_shell, E_int_eV, gate="none"):
        # Energetically-frozen ions are exempt from the cooling bound (and
        # far-out ions pass the exposure bound), so the frozen lanes run on
        # the ungated arm too -- the ungated per-ion no-op check. Suppressed
        # ions are cooling-reversible and need the gated arm (review-hardened
        # guard, 2026-07-07).
        cfg = _detect_cfg(cooling_spatial_gate=gate)
        seed = _far_seed(cfg, n_shell=n_shell, E_int_eV=E_int_eV)
        return cfg, run_detection_stage(seed, cfg)

    def test_frozen_lane(self):
        cfg = _detect_cfg(cooling_spatial_gate="none")
        picture, kappa = _ladder(cfg)
        d0_21 = float(d0_of_n(21, picture=picture, kappa=kappa))
        _, result = self._run_lane(n_shell=21, E_int_eV=0.5 * d0_21)
        np.testing.assert_array_equal(result.n_detected, np.full(4, 21.0))
        assert list(result.state_reason) == ["frozen"] * 4
        assert result.event_time_ps.size == 0
        np.testing.assert_array_equal(result.event_offsets, np.zeros(5, dtype=int))

    def test_suppressed_lane(self):
        # Gated arm: an ejected suppressed complex has zero cooling (erfc
        # density underflows), so suppression is genuinely permanent and the
        # ion rides to the detector at its handover n.
        cfg = _detect_cfg()
        picture, kappa = _ladder(cfg)
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        _, result = self._run_lane(
            n_shell=21, E_int_eV=1.5 * sigma21, gate="density_scaled",
        )
        np.testing.assert_array_equal(result.n_detected, np.full(4, 21.0))
        assert list(result.state_reason) == ["suppressed"] * 4
        assert result.event_time_ps.size == 0

    def test_ungated_suppressed_rejected(self):
        # Ungated, suppression is cooling-reversible with rho_cool = 1: over
        # the flight K2 WOULD drain E_int below Sigma and re-open the gate,
        # so the guard must refuse rather than record 'suppressed' forever
        # (the review-confirmed Tier-2 observable bias).
        cfg = _detect_cfg(cooling_spatial_gate="none")
        picture, kappa = _ladder(cfg)
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        seed = _far_seed(cfg, n_shell=21, E_int_eV=1.5 * sigma21)
        with pytest.raises(ValueError, match="P1-P3"):
            run_detection_stage(seed, cfg)

    def test_bare_ion_lane(self):
        _, result = self._run_lane(n_shell=0, E_int_eV=0.0)
        np.testing.assert_array_equal(result.n_detected, np.zeros(4))
        assert list(result.state_reason) == ["frozen"] * 4

    def test_rrk_underflow_band_classified_frozen(self):
        # E_int within ~1e-7 relative of D_0(21): strictly in-band (E > D_0,
        # gate margin < 0) but the RRK bracket (1e-7)^59 underflows to k == 0.
        # Review fix: this reachable lane (E2's freeze mask is strict) must
        # classify 'frozen', not crash the whole run.
        cfg = _detect_cfg()
        picture, kappa = _ladder(cfg)
        d0_21 = float(d0_of_n(21, picture=picture, kappa=kappa))
        seed = _far_seed(cfg, n_shell=21, E_int_eV=d0_21 * (1.0 + 1e-7))
        result = run_detection_stage(seed, cfg)
        assert list(result.state_reason) == ["frozen"] * 4
        np.testing.assert_array_equal(result.n_detected, np.full(4, 21.0))
        assert result.event_time_ps.size == 0

    def test_time_exhausted_lane(self):
        # A live in-band ion with an almost-immediate detection time: the
        # first waiting time (~1/k with k ~ 0.36/ps at the per-n s = 60)
        # exceeds the 1e-6 ps horizon with probability 1 - 4e-7 per ion --
        # deterministic under the fixed seed.
        cfg = _detect_cfg(detection_time_ps=20.0 + 1e-6, seed=3)
        picture, kappa = _ladder(cfg)
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        seed = _far_seed(cfg, n_shell=21, E_int_eV=0.98 * sigma21)

        result = run_detection_stage(seed, cfg)
        np.testing.assert_array_equal(result.n_detected, np.full(4, 21.0))
        assert list(result.state_reason) == ["time_exhausted"] * 4
        assert result.event_time_ps.size == 0
        # E_int untouched: constant between sheds under P1-P3.
        np.testing.assert_allclose(result.E_int_detected_eV, 0.98 * sigma21,
                                   rtol=0, atol=1e-15)


# ---------------------------------------------------------------------------
# Input / handover guards (incl. the P3 violator lanes)
# ---------------------------------------------------------------------------
class TestGuards:
    def test_disabled_raises(self):
        cfg = _detect_cfg(detection_stage_enabled=False)
        seed = _far_seed(cfg)
        with pytest.raises(ValueError, match="detection_stage_enabled"):
            run_detection_stage(seed, cfg)

    def test_non_biphasic_seed_rejected(self):
        cfg = _detect_cfg()
        seed = replace(_far_seed(cfg), mass_scenario="fixed")
        with pytest.raises(ValueError, match="mass_scenario"):
            run_detection_stage(seed, cfg)

    def test_stale_seed_column_rejected(self):
        cfg = _detect_cfg()
        seed = _far_seed(cfg)
        seed = replace(seed, mass_final_kg=seed.mass_final_kg + MASS_HE_AMU * U)
        with pytest.raises(ValueError, match="final"):
            run_detection_stage(seed, cfg)

    def test_detection_time_before_handover_rejected(self):
        # Passes the config-load floor (25 > ion window 0.2) but lies before
        # the realized handover time t_h = 30 -> the stage-entry check fires.
        cfg = _detect_cfg(detection_time_ps=25.0)
        seed = _far_seed(cfg, time_ps_start=30.0)
        with pytest.raises(ValueError, match="handover"):
            run_detection_stage(seed, cfg)

    def test_guard_rejects_live_in_bubble_ion(self):
        # Gated arm, live in-band ions at the droplet centre: rho ~ 1, the
        # neglected drain over the 8.53e6 ps flight is macroscopic -> the
        # guard must raise with the violator list and the remedy.
        cfg = _detect_cfg()
        picture, kappa = _ladder(cfg)
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        seed = _far_seed(cfg, n_shell=21, E_int_eV=0.98 * sigma21, inside=True)
        with pytest.raises(ValueError, match="P1-P3"):
            run_detection_stage(seed, cfg)

    def test_guard_rejects_live_ungated_ion_anywhere(self):
        # Ungated cooling acts everywhere (rho_cool = 1): a live ion can
        # never hand over, even fully ejected -- it must freeze in E2.
        cfg = _detect_cfg(cooling_spatial_gate="none")
        picture, kappa = _ladder(cfg)
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        seed = _far_seed(cfg, n_shell=21, E_int_eV=0.98 * sigma21)
        with pytest.raises(ValueError, match="P1-P3"):
            run_detection_stage(seed, cfg)

    def test_guard_rejects_frozen_in_bubble_ion(self):
        # The exposure bound (review fix): a FROZEN ion still inside helium
        # has live pickup (re-heats E_int by f_ret*D_0 per capture -- can
        # unfreeze the cascade) and live drag; the frozen exemption applies
        # to the cooling bound only, so this must fail loud, on both arms.
        cfg = _detect_cfg()   # lambda0 = 0.5 from the driver cfg
        picture, kappa = _ladder(cfg)
        d0_21 = float(d0_of_n(21, picture=picture, kappa=kappa))
        seed = _far_seed(cfg, n_shell=21, E_int_eV=0.5 * d0_21, inside=True)
        with pytest.raises(ValueError, match="P1-P3"):
            run_detection_stage(seed, cfg)

    def test_eps_drain_is_defensive_constant(self):
        # EPS_DRAIN is a documented module constant (order 1e-6), not a knob.
        assert EPS_DRAIN == 1e-6


# ---------------------------------------------------------------------------
# Config-load guard arms + the Wave-7 back-compat criterion
# ---------------------------------------------------------------------------
class TestConfigValidation:
    def test_disabled_is_noop(self):
        check_detection_config(_biphasic_driver_cfg())  # no raise

    def test_non_biphasic_raises(self):
        from i2_helium_md.presets import single_pulse_N2000_drag

        cfg = single_pulse_N2000_drag(
            num_molecules=2, ion_simulation_time=0.2, dt_ion=0.01, seed=0,
        )
        cfg = replace(cfg, mass_scenario="fixed", detection_stage_enabled=True,
                      detection_time_ps=8.53e6)
        with pytest.raises(ValueError, match="biphasic"):
            check_detection_config(cfg)

    def test_missing_time_raises(self):
        cfg = _detect_cfg(detection_time_ps=None)
        with pytest.raises(ValueError, match="detection_time_ps"):
            check_detection_config(cfg)

    @pytest.mark.parametrize("bad_time", [0.0, -1.0])
    def test_zero_or_negative_time_raises(self, bad_time):
        cfg = _detect_cfg(detection_time_ps=bad_time)
        with pytest.raises(ValueError, match="detection_time_ps"):
            check_detection_config(cfg)

    def test_time_inside_ion_window_raises(self):
        # ion_simulation_time = 0.2 on the tiny cfg; 0.1 lies inside it.
        cfg = _detect_cfg(detection_time_ps=0.1)
        with pytest.raises(ValueError, match="seed window"):
            check_detection_config(cfg)

    def test_time_inside_relaxation_window_raises(self):
        # With E2 enabled the floor is ion window + relaxation cap.
        cfg = _detect_cfg(detection_time_ps=0.5, relaxation_stage_enabled=True,
                          relaxation_time_ps=1.0)
        with pytest.raises(ValueError, match="seed window"):
            check_detection_config(cfg)

    def test_relaxation_stage_not_required(self):
        # The skip path (design §1 item 5): detection enabled with the
        # relaxation stage disabled is a VALID config.
        cfg = _detect_cfg()
        assert not cfg.relaxation_stage_enabled
        check_detection_config(cfg)  # no raise

    def test_zero_nu_raises(self):
        # nu = 0 merely WARNS at the biphasic guard (a legal pickup-only
        # diagnostic) but must RAISE here: the detection taxonomy is unsound
        # without a live evaporation channel.
        cfg = _detect_cfg(evap_rate_prefactor_per_ps=0.0)
        with pytest.raises(ValueError, match="evap_rate_prefactor_per_ps"):
            check_detection_config(cfg)

    def test_pre_ds_cfg_json_loads_with_defaults(self, tmp_path):
        # The Wave-7 acceptance criterion (probe plan Addendum C.2): a
        # cfg.json written before the two detection fields existed must load
        # with the defaults (stage disabled) -- no stale-artifact firing.
        run = RunDirectory(tmp_path / "9A_drag_tier2probe_backcompat")
        run.save_cfg(_biphasic_driver_cfg())
        payload = json.loads(run.cfg_path.read_text(encoding="utf-8"))
        assert payload.pop("detection_stage_enabled") is False
        assert payload.pop("detection_time_ps") is None
        run.cfg_path.write_text(json.dumps(payload), encoding="utf-8")

        cfg = run.load_cfg()
        assert cfg.detection_stage_enabled is False
        assert cfg.detection_time_ps is None
        check_detection_config(cfg)  # disabled -> no-op


# ---------------------------------------------------------------------------
# Step 4/5 -- event-loop physics: hypoexponential oracle, cold-shed kick,
# 5-term closure across the stage boundary
# ---------------------------------------------------------------------------
class TestEventLoopPhysics:
    def _rate_chain(self, cfg, n0, E0):
        """The deterministic (rate, rung) chain the seed walks, via the
        delivered rrk_rate -- k values until the permanent state."""
        picture, kappa = _ladder(cfg)
        ks = []
        n, E = int(n0), float(E0)
        while True:
            k = float(rrk_rate(
                E, n, nu=cfg.evap_rate_prefactor_per_ps, picture=picture,
                kappa=kappa, evap_rrk_dof=cfg.evap_rrk_dof,
                gate_onset_eV=cfg.gate_onset_override_eV,
            ))
            if k == 0.0:
                return ks
            ks.append(k)
            E += float(d0_of_n(n, picture=picture, kappa=kappa)) * -1.0
            n -= 1

    @staticmethod
    def _hypoexp_shed_probs(ks, T):
        """P(exactly j sheds by T) for a chain of exponential stages with
        distinct rates ks -- the analytic pytest oracle (design §1 fork 2).
        F_j(T) = 1 - sum_i [prod_{l != i} k_l/(k_l - k_i)] e^{-k_i T}."""
        J = len(ks)
        F = [1.0]
        for j in range(1, J + 1):
            k = np.asarray(ks[:j], dtype=float)
            coef = np.empty(j)
            for i in range(j):
                others = np.delete(k, i)
                coef[i] = np.prod(others / (others - k[i]))
            F.append(float(1.0 - np.sum(coef * np.exp(-k * T))))
        return np.array([F[j] - F[j + 1] for j in range(J)] + [F[J]])

    def test_multi_event_occupancy_matches_hypoexponential(self):
        # 400 iid ions walking a >= 3-stage chain with pairwise-distinct RRK
        # rates; empirical shed-count occupancy at a mid-chain horizon vs the
        # closed form. Tolerance: per-state se = sqrt(p(1-p)/400) <= 0.025;
        # atol 0.08 is > 3 sigma (fixed seed -> deterministic realization).
        cfg = _detect_cfg(seed=42, evap_rrk_dof=3.0)
        picture, kappa = _ladder(cfg)
        n0 = 5
        E0 = 0.6 * float(ladder_cumsum(n0, picture=picture, kappa=kappa))
        ks = self._rate_chain(cfg, n0, E0)
        assert len(ks) >= 3                                # non-vacuous chain
        pairwise = np.abs(np.subtract.outer(ks, ks))[np.triu_indices(len(ks), 1)]
        assert np.all(pairwise > 1e-6)                     # closed form stable

        T = float(np.sum(1.0 / np.asarray(ks[: max(1, len(ks) // 2)])))
        cfg = replace(cfg, detection_time_ps=20.0 + T)
        seed = _far_seed(cfg, n_shell=n0, E_int_eV=E0, num_molecules=200)

        result = run_detection_stage(seed, cfg)
        sheds = (n0 - result.n_detected).astype(int)
        emp = np.bincount(sheds, minlength=len(ks) + 1) / sheds.size

        probs = self._hypoexp_shed_probs(ks, T)
        np.testing.assert_allclose(emp[: len(probs)], probs, rtol=0, atol=0.08)
        # Mean shed count: se ~ std/20; 4 sigma band.
        mean_analytic = float(np.sum(np.arange(len(probs)) * probs))
        std_analytic = float(np.sqrt(
            np.sum((np.arange(len(probs)) - mean_analytic) ** 2 * probs)
        ))
        assert abs(float(np.mean(sheds)) - mean_analytic) < 4 * std_analytic / 20.0

    def test_cold_shed_kick_and_mass_lockstep(self):
        # Sequential cold-shed kicks compose to |v_final| = |v_0| * m_0/m_final
        # (momentum conservation, He leaving cold) -- the as-built correction
        # of the design's "v unchanged" phrasing. Mass stays in m<->n lockstep.
        cfg = _detect_cfg(seed=7, evap_rrk_dof=2.0)
        picture, kappa = _ladder(cfg)
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        seed = _far_seed(cfg, n_shell=21, E_int_eV=0.98 * sigma21, speed=0.3)

        result = run_detection_stage(seed, cfg)
        assert result_events_nonzero(result)

        m0_amu = MASS_I_ION_AMU + 21 * MASS_HE_AMU
        m_final_amu = result.mass_detected_kg / U
        np.testing.assert_allclose(
            m_final_amu, MASS_I_ION_AMU + result.n_detected * MASS_HE_AMU,
            rtol=0, atol=1e-9,
        )
        v_final = np.sqrt(result.vx_detected ** 2 + result.vy_detected ** 2
                          + result.vz_detected ** 2)
        np.testing.assert_allclose(
            v_final, 0.3 * m0_amu / m_final_amu, rtol=1e-12, atol=0,
        )

    def test_five_term_closure_across_stage_boundary(self):
        # The 5-term sum E_kin + E_pot + E_dissip + E_mass_transfer + E_int is
        # invariant across the stage (design §5 step 5): per fire the K1 drain
        # cancels the e_bind fold and the cold-shed dE_kin cancels the booked
        # defect. Pure arithmetic chain (no integrator drift) -> atol 1e-9 eV.
        cfg = _detect_cfg(seed=11, evap_rrk_dof=2.0)
        picture, kappa = _ladder(cfg)
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        seed = _far_seed(cfg, n_shell=21, E_int_eV=0.98 * sigma21, speed=0.5)

        result = run_detection_stage(seed, cfg)
        assert result_events_nonzero(result)               # non-vacuous

        before = (seed.E_kin_eV[:, -1] + seed.E_pot_eV[:, -1]
                  + seed.E_dissip_eV[:, -1] + seed.E_mass_transfer_eV[:, -1]
                  + seed.E_int_eV[:, -1])
        after = (result.E_kin_detected_eV + result.E_pot_detected_eV
                 + result.E_dissip_detected_eV
                 + result.E_mass_transfer_detected_eV
                 + result.E_int_detected_eV)
        np.testing.assert_allclose(after, before, rtol=0, atol=1e-9)

        # Per-event ledger booking sums equal the per-ion field deltas.
        for i in range(4):
            ev = result.events_for_ion(i)
            assert result.E_int_detected_eV[i] - seed.E_int_eV[i, -1] == (
                pytest.approx(float(np.sum(result.event_dE_int_eV[ev])),
                              abs=1e-12)
            )
            assert result.E_pot_detected_eV[i] - seed.E_pot_eV[i, -1] == (
                pytest.approx(float(np.sum(result.event_dE_bind_fold_eV[ev])),
                              abs=1e-12)
            )
            assert (result.E_mass_transfer_detected_eV[i]
                    - seed.E_mass_transfer_eV[i, -1]) == (
                pytest.approx(
                    float(np.sum(result.event_dE_mass_transfer_eV[ev])),
                    abs=1e-12,
                )
            )
        # E_dissip is carried through unchanged (P3: zero cooling drain).
        np.testing.assert_array_equal(result.E_dissip_detected_eV,
                                      seed.E_dissip_eV[:, -1])


# ---------------------------------------------------------------------------
# Step 6 -- equivalence with the fixed-dt E2 chain + the skip path
# ---------------------------------------------------------------------------
class TestEquivalenceAndSkipPath:
    def test_statistical_agreement_with_fixed_dt_e2(self):
        # Same handover state, short horizon T = 1 ps, gated far-out seed (the
        # pure post-ejection jump-chain regime): the event-driven stage and a
        # fixed-dt E2 extension must agree in distribution (design §2.5). 400
        # ions, independent streams. Tolerance: sampling se of the mean-shed
        # difference ~ sqrt(2)*std/20 ~ 0.1 at std ~ 1.4, plus the O(k*dt)
        # Bernoulli bias ~ 0.02 -> 0.45 is a > 4 sigma band (fixed seeds).
        T = 1.0
        cfg = _detect_cfg(seed=5, evap_rrk_dof=2.0,
                          detection_time_ps=20.0 + T)
        picture, kappa = _ladder(cfg)
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        seed = _far_seed(cfg, n_shell=21, E_int_eV=0.98 * sigma21,
                         num_molecules=200)

        detected = run_detection_stage(seed, cfg)

        relax_cfg = replace(cfg, relaxation_stage_enabled=True,
                            relaxation_time_ps=T, relaxation_dt_ps=0.01,
                            relaxation_forces="free_flight")
        relaxed = run_relaxation_stage(seed, relax_cfg)

        d_mean = float(np.mean(detected.n_detected)) - float(
            np.mean(relaxed.terminal_n)
        )
        assert np.any(detected.n_detected < 21)            # non-vacuous
        assert np.any(relaxed.terminal_n < 21)
        assert abs(d_mean) < 0.45

    def test_skip_path_matches_relaxation_seeded_noop(self):
        # Design §5 step 3 skip-path lane: an ion.npz-seeded call equals a
        # relaxation.npz-seeded call whose E2 window did nothing to (n, E_int)
        # dynamics -- realizable exactly on a frozen seed (K2 cools a frozen
        # ion's E_int but never unfreezes it; n is untouched).
        cfg = _detect_cfg(cooling_spatial_gate="none")
        picture, kappa = _ladder(cfg)
        d0_21 = float(d0_of_n(21, picture=picture, kappa=kappa))
        ion_seed = _far_seed(cfg, n_shell=21, E_int_eV=0.5 * d0_21)

        # Skip path: detection straight from the ion checkpoint.
        direct = run_detection_stage(ion_seed, cfg)

        # E2 path: relax first (freezes immediately), then detect from the
        # relaxation checkpoint.
        relax_cfg = replace(cfg, relaxation_stage_enabled=True,
                            relaxation_time_ps=0.05,
                            relaxation_forces="free_flight")
        relaxed = run_relaxation_stage(ion_seed, relax_cfg)
        chained = run_detection_stage(relaxed.checkpoint, cfg)

        np.testing.assert_array_equal(direct.n_detected, chained.n_detected)
        assert list(direct.state_reason) == list(chained.state_reason) == (
            ["frozen"] * 4
        )
        assert direct.event_time_ps.size == chained.event_time_ps.size == 0

    def test_skip_path_rejects_real_in_bubble_ion_checkpoint(self):
        # The real tiny biphasic driver output (relaxation never run): the
        # onset E_int = f_int*E = 0.24 eV > Sigma(21) leaves every ion
        # suppressed *while still inside the droplet* at the 0.2 ps sim-end.
        # Review-hardened guard: suppression is cooling-reversible and pickup
        # is live in-bubble, so this genuine ion.npz-shaped seed is exactly
        # the wrong input the skip path's guard exists to refuse loudly
        # (the pre-review code silently rode these ions to the detector at
        # n = 21 -- the confirmed Tier-2 observable bias).
        cfg = _detect_cfg()
        assert not cfg.relaxation_stage_enabled
        ion = run_ion_propagation(cfg, _tiny_neutral())

        with pytest.raises(ValueError, match="P1-P3"):
            run_detection_stage(ion, cfg)


# ---------------------------------------------------------------------------
# Step 2/7 -- artifact round-trip, fail-loud load, E1 admission, report smoke
# ---------------------------------------------------------------------------
class TestArtifactAndIntegration:
    def _result(self, tmp_path=None):
        cfg = _detect_cfg(seed=9, evap_rrk_dof=2.0)
        picture, kappa = _ladder(cfg)
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        seed = _far_seed(cfg, n_shell=21, E_int_eV=0.98 * sigma21)
        save = None if tmp_path is None else tmp_path / "detection.npz"
        return run_detection_stage(seed, cfg, save_path=save)

    def test_round_trip(self, tmp_path):
        result = self._result(tmp_path)
        path = tmp_path / "detection.npz"
        assert path.exists()

        loaded = load_detection_result(path)
        assert loaded.num_molecules == result.num_molecules
        assert loaded.t_handover_ps == result.t_handover_ps
        assert loaded.detection_time_ps == result.detection_time_ps
        for name in ("n_detected", "E_int_detected_eV", "mass_detected_kg",
                     "vx_detected", "E_kin_detected_eV", "E_pot_detected_eV",
                     "event_offsets", "event_time_ps", "event_pre_shed_n",
                     "event_dE_int_eV", "event_dE_bind_fold_eV",
                     "event_dE_mass_transfer_eV"):
            np.testing.assert_array_equal(getattr(loaded, name),
                                          getattr(result, name), err_msg=name)
        assert list(loaded.state_reason) == list(result.state_reason)

    def test_load_rejects_wrong_schema_version(self, tmp_path):
        result = self._result()
        path = save_detection_result(result, tmp_path / "detection.npz")
        with np.load(path, allow_pickle=False) as npz:
            raw = {k: npz[k] for k in npz.files}
        raw["detection_schema_version"] = np.asarray(99)
        np.savez_compressed(path, **raw)
        with pytest.raises(ValueError, match="detection_schema_version"):
            load_detection_result(path)

    def test_load_rejects_missing_field(self, tmp_path):
        result = self._result()
        path = save_detection_result(result, tmp_path / "detection.npz")
        with np.load(path, allow_pickle=False) as npz:
            raw = {k: npz[k] for k in npz.files}
        raw.pop("event_time_ps")
        np.savez_compressed(path, **raw)
        with pytest.raises(ValueError, match="missing"):
            load_detection_result(path)

    def test_load_rejects_corrupt_offsets(self, tmp_path):
        result = self._result()
        path = save_detection_result(result, tmp_path / "detection.npz")
        with np.load(path, allow_pickle=False) as npz:
            raw = {k: npz[k] for k in npz.files}
        raw["event_offsets"] = raw["event_offsets"][:-1]   # wrong shape (2N,)
        np.savez_compressed(path, **raw)
        with pytest.raises(ValueError, match="event_offsets"):
            load_detection_result(path)

    def test_e1_admission_detected_source(self):
        result = self._result()
        dist = compute_terminal_shell_distribution(result)
        assert dist.source == "detected"
        assert int(dist.counts.sum()) == 4
        for n_val in np.rint(result.n_detected).astype(int):
            assert dist.counts[n_val] >= 1

    def test_probe_report_reads_detection_artifact(self, tmp_path_factory):
        # Tiny-N driver/report smoke (design §5 step 7): a genuine probe run
        # scores "-" in the n_detect_* columns before detection.npz exists and
        # real values after -- the Wave-7 read path end-to-end.
        from scripts import gen_tier2_staircase_probe as gen_script
        from scripts.post_processing import tier2_staircase_probe_report as report
        from scripts.tier2_common import (
            build_biphasic_cfg,
            tier2_probe_run_dir_name,
        )
        from i2_helium_md.physics.shell_schedule import build_shell_schedule

        # 40 ps ungated relaxation: the cascade cools through the band and
        # freezes with the ions well outside the droplet -- a P1-P3-clean
        # handover (a 1 ps cap would leave suppressed ions in-bubble, which
        # the review-hardened guard rightly refuses).
        cfg = build_biphasic_cfg(
            "9A", "shared_pure_cubic", num_molecules=2, ion_time_ps=0.02,
            dt_ion_ps=0.01, seed=123, picture="statistical_mixture", kappa=1.0,
            tau_ps=6.55, f_int=0.5, f_ret=0.1, coulomb_available_eV=0.80,
            relaxation_time_ps=40.0,
        )
        cfg = replace(cfg, t_max_neutral=0.1, dt_neutral=0.01)
        run_dir = tmp_path_factory.mktemp("ds_probe") / tier2_probe_run_dir_name(
            "9A", "shared_pure_cubic", 2, picture="statistical_mixture",
            kappa=1.0, lambda0_per_ps=0.9, f_int=0.5, f_ret=0.1, tau_ps=6.55,
            budget_eV=0.80,
        )
        gen_script._run_one("ds-smoke", cfg, run_dir)
        schedule = build_shell_schedule(5.0)

        row_before = report.score_probe_run(run_dir, schedule)
        assert row_before["n_detect_mean"] is None
        assert row_before["frac_det_suppressed"] is None

        det_cfg = replace(RunDirectory(run_dir).load_cfg(),
                          detection_stage_enabled=True,
                          detection_time_ps=8.53e6)
        relax = load_ion_checkpoint(run_dir / "relaxation.npz")
        run_detection_stage(relax, det_cfg,
                            save_path=run_dir / "detection.npz")

        row = report.score_probe_run(run_dir, schedule)
        assert set(row) == set(report.PROBE_TABLE_COLUMNS)
        # Ungated 40 ps relaxation converges (frac_frozen = 1); the frozen
        # ions no-op through the detection stage, so the detected read equals
        # the relaxed read exactly -- the Wave-7 ungated wiring oracle.
        assert row["frac_frozen"] == pytest.approx(1.0)
        assert row["frac_det_frozen"] == pytest.approx(1.0)
        assert row["frac_det_suppressed"] == pytest.approx(0.0)
        assert row["frac_det_time_exhausted"] == pytest.approx(0.0)
        assert row["n_detect_mean"] == pytest.approx(row["n_relaxed_mean"])
        assert row["n_detect_min"] == row["n_relaxed_min"]

        # Skip-path scoring (review fix): a dir WITHOUT relaxation.npz must
        # still be discovered and scored -- relaxed columns "-", detected
        # columns populated.
        import shutil

        skip_root = tmp_path_factory.mktemp("ds_probe_skip")
        skip_dir = skip_root / run_dir.name
        shutil.copytree(run_dir, skip_dir)
        (skip_dir / "relaxation.npz").unlink()

        assert report.discover_probe_run_dirs(skip_root) == [skip_dir]
        skip_row = report.score_probe_run(skip_dir, schedule)
        assert skip_row["n_relaxed_mean"] is None
        assert skip_row["frac_frozen"] is None
        assert skip_row["n_detect_mean"] == pytest.approx(row["n_detect_mean"])
        assert skip_row["frac_det_frozen"] == pytest.approx(1.0)

        # Stale-artifact coherence guard (review fix): a detection.npz whose
        # ensemble size disagrees with cfg.json fails loud instead of scoring
        # the wrong run's detected read.
        stale_cfg = _detect_cfg()
        picture, kappa = _ladder(stale_cfg)
        d0_21 = float(d0_of_n(21, picture=picture, kappa=kappa))
        stale = run_detection_stage(
            _far_seed(stale_cfg, n_shell=21, E_int_eV=0.5 * d0_21,
                      num_molecules=3),
            stale_cfg,
        )
        save_detection_result(stale, skip_dir / "detection.npz")
        with pytest.raises(ValueError, match="stale detection.npz"):
            report.score_probe_run(skip_dir, schedule)


# ---------------------------------------------------------------------------
# Slice T2 (§I.10): tabulated-ladder wiring through the detection stage.
# ---------------------------------------------------------------------------
class TestTabulatedLadderSliceT2:
    def _rungs(self, cfg, length=32):
        from tests.ladder_feeds import form_u_rungs_for

        return form_u_rungs_for(cfg, length=length)

    def test_form_u_fed_table_detection_byte_identical(self):
        # The §I.10 equivalence oracle at the detector: a tabulated cfg fed the
        # Form-U rungs reproduces the form_u detected read bitwise (same
        # cfg.seed-derived stage stream; the waiting-time draws see identical
        # rates).
        cfg_u = _detect_cfg()
        cfg_t = _detect_cfg(
            dissociation_ladder="tabulated",
            tabulated_ladder_rungs_eV=self._rungs(_detect_cfg()),
        )
        res_u = run_detection_stage(_far_seed(cfg_u, E_int_eV=0.15), cfg_u)
        res_t = run_detection_stage(_far_seed(cfg_t, E_int_eV=0.15), cfg_t)
        assert res_u.event_offsets[-1] > 0    # self-check: the cascade fired
        np.testing.assert_array_equal(res_t.n_detected, res_u.n_detected)
        np.testing.assert_array_equal(
            res_t.E_int_detected_eV, res_u.E_int_detected_eV
        )
        np.testing.assert_array_equal(res_t.state_reason, res_u.state_reason)
        np.testing.assert_array_equal(res_t.event_time_ps, res_u.event_time_ps)
        np.testing.assert_array_equal(
            res_t.event_dE_int_eV, res_u.event_dE_int_eV
        )

    def test_cascade_to_bare_ion_tabulated_matches_form_u(self):
        # Review fix 2026-07-16: the n=1 direct fire leaves a bare ion and the
        # cascade loop re-enters rrk_rate at n=0 -- under a tabulated ladder
        # that lane must ride the Form-U k=0 path (permanent "frozen"), not
        # raise the table-range ValueError. This is the bare-I+ peak path, so
        # every production tabulated run exercises it.
        cfg_u = _detect_cfg()
        picture, kappa = _ladder(cfg_u)
        d0_1 = float(d0_of_n(1, picture=picture, kappa=kappa))
        cfg_t = _detect_cfg(
            dissociation_ladder="tabulated",
            tabulated_ladder_rungs_eV=self._rungs(_detect_cfg()),
        )
        res_u = run_detection_stage(
            _far_seed(cfg_u, n_shell=1, E_int_eV=2.0 * d0_1), cfg_u
        )
        res_t = run_detection_stage(
            _far_seed(cfg_t, n_shell=1, E_int_eV=2.0 * d0_1), cfg_t
        )
        assert int(np.min(res_u.n_detected)) == 0    # self-check: bare I+ reached
        np.testing.assert_array_equal(res_t.n_detected, res_u.n_detected)
        np.testing.assert_array_equal(res_t.state_reason, res_u.state_reason)
        np.testing.assert_array_equal(
            res_t.E_int_detected_eV, res_u.E_int_detected_eV
        )

    def test_tabulated_without_rungs_fails_loud(self):
        # Slice T2 supersedes the DS-review NotImplementedError refusal: the
        # 'tabulated' arm is wired, so the lazy refusal moves to the data path.
        cfg = _detect_cfg(dissociation_ladder="tabulated")
        with pytest.raises(ValueError, match="tabulated_ladder_rungs_eV"):
            run_detection_stage(_far_seed(cfg), cfg)
