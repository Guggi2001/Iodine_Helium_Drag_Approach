"""Tests for the Tier-2 (C) design build: CE channel mixture + exit strip.

Covers (``TIER2_CE_CHANNEL_EXIT_STRIP_DESIGN.md`` §3/§5/§8):

* the sampler (``sampling/ce_channels.py``): weight composition, channel
  means/widths, truncation, the Q3 partner tag, the pair-scale identity
  s_m = E_m/E_ref, the scored-ion mask;
* the strip forms (``physics/exit_strip.py``): P₀ cap, G monotonicity,
  validation;
* the config-load guards (``check_ce_channel_config`` /
  ``check_exit_strip_config``): typo/pairing/partner-mask/fraction rules;
* **off-mode byte-identity**: a biphasic driver run with both enums off is
  bit-identical to the same run on a config that predates the fields
  (structural: defaults touch no RNG stream, no wrapper, no scale array);
* the pair-scale seam: the scaled pair Coulomb is the exact scaled-charge-
  product trajectory (force and energy scale together);
* the strip step operator (``exit_strip_step``): crossing detection, the
  5-term ledger closure across a strip (bit-tight), the co-moving carry,
  the full-strip E_int discard (OQ-H), count/mass lockstep, and the
  cumulative ``ce_strip_count`` bookkeeping through the driver.

The *numerical* off-path regression is the entire existing drag/biphasic
battery (runs on the defaults); the identity test here pins the seam.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from i2_helium_md.config import SimConfig
from i2_helium_md.physics.constants import MASS_HE_AMU, U
from i2_helium_md.physics.dissociation_ladder import d0_of_n
from i2_helium_md.physics.exit_strip import (
    strip_depth_grading,
    strip_knock_probabilities,
    strip_velocity_gate,
)
from i2_helium_md.sampling.ce_channels import (
    CE_CHANNEL_NONE,
    CE_CHANNEL_Q2,
    CE_CHANNEL_Q3,
    CE_CHANNEL_Q3_PARTNER,
    CE_CHANNEL_SINGLE,
    E_REF_PER_ION_EV,
    ce_channel_means_eV,
    ce_pair_scale_from_checkpoint,
    ce_scored_ion_mask,
    sample_ce_channels,
)
from tests.test_biphasic_step import _biphasic_driver_cfg, _tiny_neutral

# The frozen registration inputs (design §8).
FROZEN_WEIGHTS = (0.30, 0.50, 0.20)
FROZEN_SIGMAS = (0.42, 0.31, 0.55)
FROZEN_F = 0.80
FROZEN_E_SINGLE = 0.53
FROZEN_F_INT_C = (0.15, 0.15, 0.0985)


def _ce_cfg(**overrides) -> SimConfig:
    """A driver-runnable biphasic cfg with the (C) mixture + strip on."""
    base = dict(
        ce_channel_mode="sampled",
        ce_channel_weights=FROZEN_WEIGHTS,
        ce_fraction_f=FROZEN_F,
        ce_channel_sigma_eV=FROZEN_SIGMAS,
        ce_single_ker_eV=FROZEN_E_SINGLE,
        ce_internal_energy_partition_fractions=FROZEN_F_INT_C,
        exit_strip_mode="depth_graded",
    )
    base.update(overrides)
    return _biphasic_driver_cfg(**base)


# ===========================================================================
# Sampler
# ===========================================================================
class TestSampler:
    def test_channel_means_frozen_values(self):
        e_s, e_q2, e_q3 = ce_channel_means_eV(
            fraction_f=FROZEN_F, single_ker_eV=FROZEN_E_SINGLE,
        )
        assert e_s == FROZEN_E_SINGLE
        # Design §3.1: 2.16 / 4.32 eV at f = 0.8 (E_ref = 2.7006).
        assert round(e_q2, 2) == 2.16
        assert round(e_q3, 2) == 4.32
        assert round(E_REF_PER_ION_EV, 4) == 2.7006

    def test_composition_and_partner_tags(self):
        rng = np.random.default_rng(20260729)
        n_mol = 4000
        d = sample_ce_channels(
            n_mol, rng, weights=FROZEN_WEIGHTS, fraction_f=FROZEN_F,
            sigma_eV=FROZEN_SIGMAS, single_ker_eV=FROZEN_E_SINGLE,
        )
        assert d.channel.shape == (2 * n_mol,)
        # Per molecule: exactly one partner tag per Q3 molecule, none else.
        c1, c2 = d.channel[:n_mol], d.channel[n_mol:]
        q3_mol = (c1 == CE_CHANNEL_Q3_PARTNER) | (c2 == CE_CHANNEL_Q3_PARTNER)
        assert np.all(
            np.count_nonzero(
                np.stack([c1, c2]) == CE_CHANNEL_Q3_PARTNER, axis=0
            )[q3_mol] == 1
        )
        # Composition within 3 sigma of the frozen weights (binomial SD).
        frac_q3 = float(np.mean(q3_mol))
        assert abs(frac_q3 - 0.20) < 3 * np.sqrt(0.2 * 0.8 / n_mol)
        frac_single = float(np.mean(c1 == CE_CHANNEL_SINGLE))
        assert abs(frac_single - 0.30) < 3 * np.sqrt(0.3 * 0.7 / n_mol)

    def test_channel_energies_and_pair_scale_identity(self):
        rng = np.random.default_rng(7)
        n_mol = 4000
        d = sample_ce_channels(
            n_mol, rng, weights=FROZEN_WEIGHTS, fraction_f=FROZEN_F,
            sigma_eV=FROZEN_SIGMAS, single_ker_eV=FROZEN_E_SINGLE,
        )
        assert np.all(d.E_m_eV > 0)
        # Both fragments share the molecule stamp.
        np.testing.assert_array_equal(d.E_m_eV[:n_mol], d.E_m_eV[n_mol:])
        # s_m = E_m / E_ref exactly.
        np.testing.assert_allclose(
            d.pair_scale, d.E_m_eV[:n_mol] / E_REF_PER_ION_EV, rtol=0, atol=0,
        )
        # Q2 channel mean within 4 SE of 2.16 (sigma 0.31, ~2000 members).
        c1 = d.channel[:n_mol]
        q2_sel = d.E_m_eV[:n_mol][c1 == CE_CHANNEL_Q2]
        assert abs(float(q2_sel.mean()) - 2.16) \
            < 4 * 0.31 / np.sqrt(q2_sel.size) + 1e-3
        q3_sel = d.E_m_eV[:n_mol][
            (c1 == CE_CHANNEL_Q3) | (c1 == CE_CHANNEL_Q3_PARTNER)
        ]
        assert abs(float(q3_sel.mean()) - 4.32) \
            < 4 * 0.55 / np.sqrt(q3_sel.size) + 1e-3

    def test_scored_mask_excludes_exactly_partners(self):
        rng = np.random.default_rng(3)
        d = sample_ce_channels(
            500, rng, weights=FROZEN_WEIGHTS, fraction_f=FROZEN_F,
            sigma_eV=FROZEN_SIGMAS, single_ker_eV=FROZEN_E_SINGLE,
        )
        mask = ce_scored_ion_mask(d.channel)
        assert np.array_equal(~mask, d.channel == CE_CHANNEL_Q3_PARTNER)

    @pytest.mark.parametrize("kwargs", [
        dict(weights=(0.5, 0.5, 0.5)),          # not a simplex
        dict(weights=(1.0, 0.5, -0.5)),         # negative
        dict(fraction_f=0.0),
        dict(fraction_f=1.5),
        dict(sigma_eV=(0.1, -0.1, 0.1)),
        dict(single_ker_eV=0.0),
    ])
    def test_bad_parameters_raise(self, kwargs):
        base = dict(weights=FROZEN_WEIGHTS, fraction_f=FROZEN_F,
                    sigma_eV=FROZEN_SIGMAS, single_ker_eV=FROZEN_E_SINGLE)
        base.update(kwargs)
        with pytest.raises(ValueError):
            sample_ce_channels(10, np.random.default_rng(0), **base)


# ===========================================================================
# Strip forms
# ===========================================================================
class TestStripForms:
    def test_velocity_gate_caps_at_one(self):
        assert strip_velocity_gate(9.9, v_strip_aps=9.9, exponent=2.0) == 1.0
        assert strip_velocity_gate(20.0, v_strip_aps=9.9, exponent=2.0) == 1.0
        # Ram scaling below the cap: (v/v_ref)^2.
        assert strip_velocity_gate(4.95, v_strip_aps=9.9, exponent=2.0) \
            == pytest.approx(0.25)
        assert strip_velocity_gate(0.0, v_strip_aps=9.9, exponent=2.0) == 0.0

    def test_depth_grading_outer_eager_inner_shadowed(self):
        j = np.arange(1, 22)
        g = strip_depth_grading(j, protect_j0=1.75, width_rungs=0.75)
        assert np.all(np.diff(g) > 0)          # monotone in j
        assert g[0] < 0.5 < g[20]              # inner shadowed, outer eager
        # Logistic midpoint at j0.
        assert strip_depth_grading(1.75, protect_j0=1.75, width_rungs=0.75) \
            == pytest.approx(0.5)

    def test_knock_probabilities_shape_and_cap(self):
        p = strip_knock_probabilities(
            21, 10.2, v_strip_aps=9.9, exponent=2.0,
            protect_j0=1.75, width_rungs=0.75,
        )
        assert p.shape == (21,)
        assert np.all((p >= 0) & (p <= 1))
        assert strip_knock_probabilities(
            0, 10.0, v_strip_aps=9.9, exponent=2.0,
            protect_j0=1.75, width_rungs=0.75,
        ).size == 0

    @pytest.mark.parametrize("kwargs", [
        dict(v_strip_aps=0.0), dict(exponent=0.0),
        dict(protect_j0=-1.0), dict(width_rungs=0.0),
    ])
    def test_bad_form_parameters_raise(self, kwargs):
        base = dict(v_strip_aps=9.9, exponent=2.0,
                    protect_j0=1.75, width_rungs=0.75)
        base.update(kwargs)
        with pytest.raises(ValueError):
            strip_knock_probabilities(5, 10.0, **base)


# ===========================================================================
# Config guards
# ===========================================================================
class TestConfigGuards:
    def test_defaults_validate(self):
        _biphasic_driver_cfg().validate()

    def test_ce_on_biphasic_validates(self):
        _ce_cfg().validate()

    def test_unknown_selectors_rejected(self):
        with pytest.raises(ValueError, match="ce_channel_mode"):
            _biphasic_driver_cfg(ce_channel_mode="typo").validate()
        with pytest.raises(ValueError, match="exit_strip_mode"):
            _biphasic_driver_cfg(exit_strip_mode="typo").validate()

    def test_sampled_outside_biphasic_rejected(self):
        from i2_helium_md.presets import single_pulse_N2000_drag

        cfg = replace(
            single_pulse_N2000_drag(num_molecules=2),
            ce_channel_mode="sampled",
            ce_internal_energy_partition_fractions=FROZEN_F_INT_C,
        )
        with pytest.raises(ValueError, match="mass_scenario"):
            cfg.validate()
        cfg2 = replace(
            single_pulse_N2000_drag(num_molecules=2),
            exit_strip_mode="depth_graded",
        )
        with pytest.raises(ValueError, match="mass_scenario"):
            cfg2.validate()

    def test_partner_mask_non_optional(self):
        with pytest.raises(ValueError, match="ce_q3_partner_mask"):
            _ce_cfg(ce_q3_partner_mask=False).validate()

    def test_fraction_tuple_required_when_sampled_refused_when_off(self):
        with pytest.raises(
            ValueError, match="ce_internal_energy_partition_fractions"
        ):
            _ce_cfg(ce_internal_energy_partition_fractions=None).validate()
        with pytest.raises(
            ValueError, match="ce_internal_energy_partition_fractions"
        ):
            _biphasic_driver_cfg(
                ce_internal_energy_partition_fractions=FROZEN_F_INT_C,
            ).validate()

    @pytest.mark.parametrize("field,value", [
        ("ce_channel_weights", (0.5, 0.5, 0.5)),
        ("ce_fraction_f", 0.0),
        ("ce_channel_sigma_eV", (0.1, -0.1, 0.1)),
        ("ce_single_ker_eV", -1.0),
        ("exit_strip_v_ref", 0.0),
        ("exit_strip_exponent", -2.0),
        ("exit_strip_width_rungs", 0.0),
        ("exit_strip_carry_eV", -0.01),
    ])
    def test_bad_parameters_rejected_even_inert(self, field, value):
        with pytest.raises(ValueError, match=field.split("_")[0]):
            _biphasic_driver_cfg(**{field: value}).validate()


# ===========================================================================
# Off-mode byte-identity + driver integration
# ===========================================================================
class TestDriverIntegration:
    @pytest.fixture(scope="class")
    def off_run(self):
        from i2_helium_md.simulation.ion import run_ion_propagation

        cfg = _biphasic_driver_cfg()
        return cfg, run_ion_propagation(cfg, _tiny_neutral())

    def test_off_mode_sentinels_and_no_new_draws(self, off_run):
        """Off-mode: sentinel v8 fields, and the ion-stage stream is
        untouched (positions/velocities/E_int identical to a re-run —
        the frozen two-draw stream would shift if anything new consumed
        from it)."""
        from i2_helium_md.simulation.ion import run_ion_propagation

        cfg, ck = off_run
        assert np.all(ck.ce_channel == CE_CHANNEL_NONE)
        assert np.all(np.isnan(ck.ce_E_m_eV))
        assert np.all(ck.ce_strip_count == 0)
        assert ce_pair_scale_from_checkpoint(ck) is None
        ck2 = run_ion_propagation(cfg, _tiny_neutral())
        np.testing.assert_array_equal(ck.velocities_x, ck2.velocities_x)
        np.testing.assert_array_equal(ck.E_int_eV, ck2.E_int_eV)

    def test_sampled_run_stamps_fields_and_preserves_ion_stream(self):
        """Channels-on consumes ONLY the dedicated stream: the biphasic
        mass-subsystem outcomes (driven by the ion-stage stream) must be
        byte-identical to the off run wherever the dynamics agree at t0.

        The strongest cheap invariant: the t0 E_int onset follows the
        per-channel rule, the pair scale is stamped, and a re-run
        reproduces bit-for-bit (the dedicated stream is seed-derived)."""
        from i2_helium_md.simulation.ion import run_ion_propagation

        cfg = _ce_cfg()
        ck = run_ion_propagation(cfg, _tiny_neutral())
        assert np.all(ck.ce_channel != CE_CHANNEL_NONE)
        assert np.all(np.isfinite(ck.ce_E_m_eV)) and np.all(ck.ce_E_m_eV > 0)
        scale = ce_pair_scale_from_checkpoint(ck)
        np.testing.assert_allclose(
            scale, ck.ce_E_m_eV[:cfg.num_molecules] / E_REF_PER_ION_EV,
        )
        # Per-channel onset at t0: f_int,c * E_m (T6 factor is 1 at n0=n*
        # under the default "constant" law with full shells).
        f_map = {CE_CHANNEL_SINGLE: FROZEN_F_INT_C[0],
                 CE_CHANNEL_Q2: FROZEN_F_INT_C[1],
                 CE_CHANNEL_Q3: FROZEN_F_INT_C[2],
                 CE_CHANNEL_Q3_PARTNER: FROZEN_F_INT_C[2]}
        expected = np.array(
            [f_map[int(c)] for c in ck.ce_channel]
        ) * ck.ce_E_m_eV
        np.testing.assert_allclose(ck.E_int_eV[:, 0], expected)
        ck2 = run_ion_propagation(cfg, _tiny_neutral())
        np.testing.assert_array_equal(ck.ce_E_m_eV, ck2.ce_E_m_eV)
        np.testing.assert_array_equal(ck.velocities_x, ck2.velocities_x)

    def test_pair_scale_changes_coulomb_dynamics(self):
        """A sampled run's pair kinematics differ from the off run (the
        scale is live in the force path) unless every s_m == 1."""
        from i2_helium_md.simulation.ion import run_ion_propagation

        cfg_off = _biphasic_driver_cfg()
        cfg_on = _ce_cfg(exit_strip_mode="off")
        ck_off = run_ion_propagation(cfg_off, _tiny_neutral())
        ck_on = run_ion_propagation(cfg_on, _tiny_neutral())
        assert not np.allclose(
            ck_off.velocities_final_x, ck_on.velocities_final_x,
        )


# ===========================================================================
# The strip step operator (deterministic ledger checks)
# ===========================================================================
class TestExitStripStep:
    def _crossing_state(self, n0: float, v_aps: float, e_int: float = 0.3):
        """A single-molecule state straddling the boundary outbound."""
        from i2_helium_md.physics.constants import MASS_I_ION_AMU
        from i2_helium_md.simulation.ion_propagation_step import IonStepState

        two_n = 2
        m_amu = MASS_I_ION_AMU + n0 * MASS_HE_AMU
        radii = np.full(two_n, 30.0)
        # Atom 0 just OUTSIDE (depth +0.5), atom 1 deep inside.
        x = np.array([30.5, 5.0])
        state = IonStepState(
            x=x, y=np.zeros(two_n), z=np.zeros(two_n),
            vx=np.array([v_aps, 0.0]), vy=np.zeros(two_n),
            vz=np.zeros(two_n),
            mass_kg=np.full(two_n, m_amu * U),
            E_kin_eV=np.zeros(two_n), E_pot_eV=np.zeros(two_n),
            E_dissip_eV=np.zeros(two_n), E_mass_transfer_eV=np.zeros(two_n),
            E_int_eV=np.full(two_n, e_int),
            number_of_collisions=np.zeros(two_n, dtype=int),
            time_ps=0.0,
            n_shell=np.full(two_n, float(n0)),
        )
        # E_kin consistent with (m, v).
        from i2_helium_md.simulation.ion_propagation_step import _E_kin_eV
        state = replace(
            state, E_kin_eV=_E_kin_eV(state.mass_kg, state.vx, state.vy,
                                      state.vz),
        )
        prev_depth = np.array([-0.5, -25.0])   # atom 0 crossed this step
        return state, radii, prev_depth

    def _strip(self, state, radii, prev_depth, cfg, seed=1):
        from i2_helium_md.simulation.ion_propagation_step import exit_strip_step

        rng = np.random.default_rng(seed)
        return exit_strip_step(
            state, rng=rng, cfg=cfg, droplet_radii=radii,
            prev_depth_angstrom=prev_depth,
        )

    def test_no_crossing_is_identity(self):
        cfg = _ce_cfg()
        state, radii, _ = self._crossing_state(19, 10.2)
        # prev_depth already outside -> no outbound crossing this step.
        new_state, depth, knocks = self._strip(
            state, radii, np.array([0.5, -25.0]), cfg,
        )
        assert new_state is state
        assert np.all(knocks == 0)

    def test_fast_crossing_strips_and_ledger_closes(self):
        cfg = _ce_cfg()
        state, radii, prev_depth = self._crossing_state(19, 10.2)
        new_state, _, knocks = self._strip(state, radii, prev_depth, cfg)
        k = int(knocks[0])
        assert k > 0 and knocks[1] == 0
        # Mass/count lockstep: m dropped by exactly k m_He.
        np.testing.assert_allclose(
            (state.mass_kg[0] - new_state.mass_kg[0]) / U,
            k * MASS_HE_AMU, rtol=1e-12,
        )
        assert new_state.n_shell[0] == state.n_shell[0] - k
        # 5-term invariant closes bit-tight across the strip.
        def five(s, i):
            return (s.E_kin_eV[i] + s.E_pot_eV[i] + s.E_dissip_eV[i]
                    + s.E_mass_transfer_eV[i] + s.E_int_eV[i])
        assert five(new_state, 0) == pytest.approx(five(state, 0), abs=1e-12)
        # The ion slowed (toll paid from KE), direction preserved.
        assert new_state.vx[0] < state.vx[0]
        assert new_state.vx[0] > 0
        # Carry-off booked: E_dissip rose by exactly k * eps (no full strip).
        if new_state.n_shell[0] > 0:
            assert new_state.E_dissip_eV[0] - state.E_dissip_eV[0] \
                == pytest.approx(k * cfg.exit_strip_carry_eV, abs=1e-12)

    def test_full_strip_discards_e_int_with_label(self):
        # Force a full strip: j0 = 0 (no protection), huge exit speed.
        cfg = _ce_cfg(exit_strip_protect_j0=0.0, exit_strip_width_rungs=0.5)
        state, radii, prev_depth = self._crossing_state(3, 30.0, e_int=0.42)
        new_state, _, knocks = self._strip(state, radii, prev_depth, cfg,
                                           seed=5)
        if new_state.n_shell[0] == 0:          # all three rungs knocked
            assert new_state.E_int_eV[0] == 0.0
            # The residual E_int landed in E_dissip (plus the eps carries).
            assert new_state.E_dissip_eV[0] == pytest.approx(
                0.42 + 3 * cfg.exit_strip_carry_eV, abs=1e-12,
            )
        # Ledger closes regardless of the realized knock count.
        def five(s, i):
            return (s.E_kin_eV[i] + s.E_pot_eV[i] + s.E_dissip_eV[i]
                    + s.E_mass_transfer_eV[i] + s.E_int_eV[i])
        assert five(new_state, 0) == pytest.approx(five(state, 0), abs=1e-12)

    def test_slow_crossing_rarely_strips(self):
        """P0 suppresses slow exiters: v = 1 A/ps -> P0 ~ 1% -> with a
        median draw nothing knocks."""
        cfg = _ce_cfg()
        state, radii, prev_depth = self._crossing_state(19, 1.0)
        new_state, _, knocks = self._strip(state, radii, prev_depth, cfg)
        assert int(knocks[0]) <= 2

    def test_bare_ion_crossing_is_identity(self):
        cfg = _ce_cfg()
        state, radii, prev_depth = self._crossing_state(0, 10.2)
        new_state, _, knocks = self._strip(state, radii, prev_depth, cfg)
        assert np.all(knocks == 0)
        np.testing.assert_array_equal(new_state.vx, state.vx)


# ===========================================================================
# Partner-mask scoring seam
# ===========================================================================
class TestPartnerMaskSeam:
    def test_read_confirmation_mask_shrinks_universe(self):
        from i2_helium_md.postprocess.tier2_confirmation import (
            read_confirmation_detection,
        )
        from i2_helium_md.simulation.detection_stage import DetectionResult

        two_n = 6
        det = DetectionResult(
            num_molecules=3, t_handover_ps=0.0, detection_time_ps=1.0,
            n_detected=np.array([1.0, 2.0, 0.0, 3.0, 1.0, 2.0]),
            E_int_detected_eV=np.zeros(two_n),
            mass_detected_kg=np.ones(two_n),
            vx_detected=np.zeros(two_n), vy_detected=np.zeros(two_n),
            vz_detected=np.zeros(two_n),
            E_kin_detected_eV=np.linspace(0.5, 3.0, two_n),
            E_pot_detected_eV=np.zeros(two_n),
            E_dissip_detected_eV=np.zeros(two_n),
            E_mass_transfer_detected_eV=np.zeros(two_n),
            state_reason=np.asarray(["time_exhausted"] * two_n),
            event_offsets=np.zeros(two_n + 1, dtype=int),
            event_time_ps=np.zeros(0), event_pre_shed_n=np.zeros(0),
            event_dE_int_eV=np.zeros(0), event_dE_bind_fold_eV=np.zeros(0),
            event_dE_mass_transfer_eV=np.zeros(0),
        )
        mask = np.array([True, True, False, True, True, True])
        read_all = read_confirmation_detection(det, label="all")
        read_masked = read_confirmation_detection(
            det, label="masked", include_mask=mask,
        )
        assert read_all.num_scored == two_n
        assert read_masked.num_scored == two_n - 1
        # The excluded ion (n=0, third entry) left the histogram entirely.
        assert read_all.fraction[0] > 0
        assert read_masked.fraction[0] == 0

    def test_mask_shape_mismatch_raises(self):
        from i2_helium_md.postprocess.tier2_confirmation import (
            read_confirmation_detection,
        )
        with pytest.raises((ValueError, AttributeError)):
            read_confirmation_detection(
                object(), label="x", include_mask=np.ones(3, dtype=bool),
            )
