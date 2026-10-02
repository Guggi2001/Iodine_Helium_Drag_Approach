"""Residual pair-Coulomb closure at the detection handover (§5.4 D1-D4, D8).

TIER2_DetectorStageFix.md §2e / §5.4. Validation order (CLAUDE.md): the pure
closure maths against a direct high-accuracy ODE integration of the same
point-charge problem (the asymptote must be a constant of the motion and the
integrated velocity must approach it), then the stage wiring on tiny synthetic
checkpoints (partner classes, D2 per-ion ledger, D3 well assertion, D4
marginal-partner safeguard), then the config-load surface.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from i2_helium_md.config import check_detection_config
from i2_helium_md.physics.coulomb_closure import (
    fixed_centre_coulomb_asymptote,
    two_body_coulomb_asymptote,
)
from i2_helium_md.physics.constants import U
from i2_helium_md.physics.interactions import ion_interaction_potential
from i2_helium_md.sampling.ce_channels import CE_CHANNEL_Q2, E_REF_PER_ION_EV
from i2_helium_md.simulation import detection_stage as ds
from i2_helium_md.simulation.detection_stage import run_detection_stage
from i2_helium_md.simulation.ion_propagation_step import (
    _E_kin_eV,
    _eV_to_amu_ang2_ps2,
    _amu_ang2_ps2_to_eV,
)

from tests.test_detection_stage import _detect_cfg, _far_seed, _ladder
from i2_helium_md.physics.dissociation_ladder import d0_of_n, ladder_cumsum

# A realistic coupling: 14.4 eV*A in amu*A^3/ps^2 (~1.39e5).
K_MD = float(_eV_to_amu_ang2_ps2(14.39964548))
M1, M2 = 131.0, 139.0          # amu, fragment-like and unequal


def _two_body_rhs(m1, m2, k):
    def rhs(_t, s):
        d = s[:3] - s[3:6]
        f = k * d / np.linalg.norm(d) ** 3
        return np.concatenate([s[6:9], s[9:], f / m1, -f / m2])
    return rhs


def _fixed_rhs(m, k):
    def rhs(_t, s):
        d = s[:3]
        return np.concatenate([s[3:], k * d / np.linalg.norm(d) ** 3 / m])
    return rhs


def _random_pair(rng):
    return (rng.normal(0, 30, 3), rng.normal(0, 5, 3),
            rng.normal(0, 30, 3), rng.normal(0, 5, 3))


# ---------------------------------------------------------------------------
# Closure maths
# ---------------------------------------------------------------------------
class TestTwoBodyAsymptote:
    def test_asymptote_is_a_constant_of_the_motion(self):
        # Integrate the exact two-body problem (DOP853, rtol 1e-12) and close
        # the END state: an exact asymptote must not depend on where along
        # the orbit it is evaluated. Tolerance 1e-9 A/ps (~1e-10 relative) is
        # the integrator's accumulated error, not the closure's.
        rng = np.random.default_rng(7)
        for _ in range(4):
            r1, v1, r2, v2 = _random_pair(rng)
            sol = solve_ivp(_two_body_rhs(M1, M2, K_MD), (0.0, 2.0e3),
                            np.concatenate([r1, r2, v1, v2]), method="DOP853",
                            rtol=1e-12, atol=1e-10)
            s = sol.y[:, -1]
            a1, a2 = two_body_coulomb_asymptote(
                r1[None], v1[None], [M1], r2[None], v2[None], [M2], [K_MD])
            b1, b2 = two_body_coulomb_asymptote(
                s[None, :3], s[None, 6:9], [M1], s[None, 3:6], s[None, 9:],
                [M2], [K_MD])
            np.testing.assert_allclose(a1, b1, rtol=0, atol=1e-9)
            np.testing.assert_allclose(a2, b2, rtol=0, atol=1e-9)

    def test_integrated_velocity_approaches_the_asymptote(self):
        # Magnitude AND direction: after a long flight the integrated velocity
        # differs from the asymptote only by the residual k/r not yet spent
        # (|dv| ~ alpha/(w r) ~ 1e-4 A/ps at r ~ 1e5 A).
        rng = np.random.default_rng(11)
        r1, v1, r2, v2 = _random_pair(rng)
        sol = solve_ivp(_two_body_rhs(M1, M2, K_MD), (0.0, 5.0e4),
                        np.concatenate([r1, r2, v1, v2]), method="DOP853",
                        rtol=1e-12, atol=1e-10)
        s = sol.y[:, -1]
        a1, a2 = two_body_coulomb_asymptote(
            r1[None], v1[None], [M1], r2[None], v2[None], [M2], [K_MD])
        np.testing.assert_allclose(a1[0], s[6:9], rtol=0, atol=1e-3)
        np.testing.assert_allclose(a2[0], s[9:], rtol=0, atol=1e-3)

    def test_momentum_conserved_and_energy_gain_is_pair_energy(self):
        rng = np.random.default_rng(3)
        P = 50
        r1, r2 = rng.normal(0, 200, (P, 3)), rng.normal(0, 200, (P, 3))
        v1, v2 = rng.normal(0, 8, (P, 3)), rng.normal(0, 8, (P, 3))
        m1 = rng.uniform(127, 220, P)
        m2 = rng.uniform(127, 220, P)
        k = np.full(P, K_MD)
        a1, a2 = two_body_coulomb_asymptote(r1, v1, m1, r2, v2, m2, k)
        np.testing.assert_allclose(m1[:, None] * a1 + m2[:, None] * a2,
                                   m1[:, None] * v1 + m2[:, None] * v2,
                                   rtol=1e-12, atol=1e-9)
        dke = (0.5 * m1 * (a1 ** 2).sum(1) + 0.5 * m2 * (a2 ** 2).sum(1)
               - 0.5 * m1 * (v1 ** 2).sum(1) - 0.5 * m2 * (v2 ** 2).sum(1))
        np.testing.assert_allclose(dke, k / np.linalg.norm(r1 - r2, axis=1),
                                   rtol=1e-9)

    def test_radial_orbit_needs_no_angular_momentum(self):
        # h = 0: head-on (inward) and outward motion both end along r_hat.
        r1, r2 = np.array([[10.0, 0, 0]]), np.array([[-10.0, 0, 0]])
        for sign in (+1.0, -1.0):
            v1 = np.array([[sign * 2.0, 0, 0]])
            v2 = -v1
            a1, a2 = two_body_coulomb_asymptote(r1, v1, [M1], r2, v2, [M1], [K_MD])
            assert np.all(np.isfinite(a1))
            assert a1[0, 0] > 0 and a2[0, 0] < 0
            np.testing.assert_allclose(a1[0, 1:], 0.0, atol=1e-12)

    def test_zero_coupling_is_identity(self):
        rng = np.random.default_rng(0)
        r1, v1, r2, v2 = (x[None] for x in _random_pair(rng))
        a1, a2 = two_body_coulomb_asymptote(r1, v1, [M1], r2, v2, [M2], [0.0])
        np.testing.assert_array_equal(a1, v1)
        np.testing.assert_array_equal(a2, v2)

    def test_mixed_coupling_batch_masks_per_pair(self):
        # The `live` mask: k = 0 rows verbatim, k > 0 rows equal to the same
        # pair closed alone (elementwise numpy ops -> round-off level only).
        rng = np.random.default_rng(13)
        P = 6
        r1, r2 = rng.normal(0, 50, (P, 3)), rng.normal(0, 50, (P, 3))
        v1, v2 = rng.normal(0, 5, (P, 3)), rng.normal(0, 5, (P, 3))
        m1, m2 = rng.uniform(127, 220, P), rng.uniform(127, 220, P)
        k = np.where(np.arange(P) % 2 == 0, K_MD, 0.0)
        a1, a2 = two_body_coulomb_asymptote(r1, v1, m1, r2, v2, m2, k)
        dead = k == 0
        np.testing.assert_array_equal(a1[dead], v1[dead])
        np.testing.assert_array_equal(a2[dead], v2[dead])
        for p in np.flatnonzero(~dead):
            b1, b2 = two_body_coulomb_asymptote(
                r1[[p]], v1[[p]], m1[[p]], r2[[p]], v2[[p]], m2[[p]], k[[p]])
            np.testing.assert_allclose(a1[p], b1[0], rtol=1e-14, atol=0)
            np.testing.assert_allclose(a2[p], b2[0], rtol=1e-14, atol=0)

    def test_near_radial_orbit_is_continuous_with_the_radial_limit(self):
        # h -> 0 must not blow up (u_hat is written with e x h, not h_hat).
        # A transverse kick eps gives c = w|h|/alpha ~ 3e-10 here, so the
        # asymptote may move by O(w*c) ~ 1e-8 A/ps; 1e-6 bounds that with
        # margin while catching any 1/|h| instability (which would be O(1)).
        r1, r2 = np.array([[10.0, 0, 0]]), np.array([[-10.0, 0, 0]])
        v1 = np.array([[2.0, 0, 0]])
        a1_rad, a2_rad = two_body_coulomb_asymptote(
            r1, v1, [M1], r2, -v1, [M2], [K_MD])
        for eps in (1e-9, 1e-12):
            v1e = v1 + np.array([[0.0, eps, 0.0]])
            a1, a2 = two_body_coulomb_asymptote(r1, v1e, [M1], r2, -v1e, [M2], [K_MD])
            np.testing.assert_allclose(a1, a1_rad, rtol=0, atol=1e-6)
            np.testing.assert_allclose(a2, a2_rad, rtol=0, atol=1e-6)

    def test_invalid_inputs_raise(self):
        z = np.zeros((1, 3))
        o = np.ones((1, 3))
        with pytest.raises(ValueError, match="coincident"):
            two_body_coulomb_asymptote(o, z, [M1], o, z, [M2], [K_MD])
        with pytest.raises(ValueError, match="k must be"):
            two_body_coulomb_asymptote(o, z, [M1], z, z, [M2], [-1.0])
        with pytest.raises(ValueError, match="masses"):
            two_body_coulomb_asymptote(o, z, [0.0], z, z, [M2], [K_MD])
        with pytest.raises(ValueError, match="shape"):
            two_body_coulomb_asymptote(np.ones(3), z, [M1], z, z, [M2], [K_MD])


class TestFixedCentreAsymptote:
    def test_speed_is_full_pair_energy(self):
        rng = np.random.default_rng(5)
        P = 20
        r, v, c = (rng.normal(0, 100, (P, 3)), rng.normal(0, 6, (P, 3)),
                   rng.normal(0, 100, (P, 3)))
        m = rng.uniform(127, 220, P)
        a = fixed_centre_coulomb_asymptote(r, v, m, c, np.full(P, K_MD))
        np.testing.assert_allclose(
            (a ** 2).sum(1),
            (v ** 2).sum(1) + 2.0 * K_MD / (m * np.linalg.norm(r - c, axis=1)),
            rtol=1e-12,
        )

    def test_asymptote_is_a_constant_of_the_fixed_centre_motion(self):
        rng = np.random.default_rng(9)
        for _ in range(3):
            r, v = rng.normal(0, 30, 3), rng.normal(0, 5, 3)
            s = solve_ivp(_fixed_rhs(M1, K_MD), (0.0, 2.0e3),
                          np.concatenate([r, v]), method="DOP853",
                          rtol=1e-12, atol=1e-10).y[:, -1]
            a = fixed_centre_coulomb_asymptote(r[None], v[None], [M1],
                                               np.zeros((1, 3)), [K_MD])
            b = fixed_centre_coulomb_asymptote(s[None, :3], s[None, 3:], [M1],
                                               np.zeros((1, 3)), [K_MD])
            np.testing.assert_allclose(a, b, rtol=0, atol=1e-9)

    def test_is_the_infinite_partner_mass_limit_of_two_body(self):
        rng = np.random.default_rng(2)
        r1, v1, r2, _ = _random_pair(rng)
        a = fixed_centre_coulomb_asymptote(r1[None], v1[None], [M1], r2[None], [K_MD])
        b1, _ = two_body_coulomb_asymptote(r1[None], v1[None], [M1], r2[None],
                                           np.zeros((1, 3)), [1e12], [K_MD])
        np.testing.assert_allclose(a, b1, rtol=1e-8)   # O(m1/m2) = 1e-10

    def test_zero_coupling_rows_are_identity_in_a_mixed_batch(self):
        rng = np.random.default_rng(17)
        r, v, c = (rng.normal(0, 50, (4, 3)), rng.normal(0, 5, (4, 3)),
                   rng.normal(0, 50, (4, 3)))
        k = np.array([K_MD, 0.0, K_MD, 0.0])
        a = fixed_centre_coulomb_asymptote(r, v, np.full(4, M1), c, k)
        np.testing.assert_array_equal(a[k == 0], v[k == 0])
        assert np.all(np.linalg.norm(a[k > 0], axis=1)
                      > np.linalg.norm(v[k > 0], axis=1))

    def test_invalid_inputs_raise(self):
        z = np.zeros((1, 3))
        o = np.ones((1, 3))
        with pytest.raises(ValueError, match="at the centre"):
            fixed_centre_coulomb_asymptote(o, z, [M1], o, [K_MD])
        with pytest.raises(ValueError, match="k must be"):
            fixed_centre_coulomb_asymptote(o, z, [M1], z, [-1.0])
        with pytest.raises(ValueError, match="masses"):
            fixed_centre_coulomb_asymptote(o, z, [0.0], z, [K_MD])
        with pytest.raises(ValueError, match="shape"):
            fixed_centre_coulomb_asymptote(np.ones(3), z, [M1], z, [K_MD])


# ---------------------------------------------------------------------------
# Stage wiring
# ---------------------------------------------------------------------------
def _cfg(**kw):
    base = dict(seed=5, evaporation_shed_convention="co_moving",
                detection_coulomb_closure="partner_aware")
    base.update(kw)
    return _detect_cfg(**base)


def _free_pair_seed(cfg, **kw):
    """N = 2, all four ions far outside (plateau), pairs separated.

    ``_far_seed`` puts each ion on top of its partner; move them apart so the
    pair energy is finite (~0.02-0.07 eV) and the fragments are not symmetric.
    """
    seed = _far_seed(cfg, **kw)
    seed.positions_x[:, 0] = [300.0, -150.0, -400.0, 500.0]
    seed.positions_y[:, 0] = [20.0, 100.0, -60.0, 0.0]
    seed.positions_z[:, 0] = [1000.0, 1200.0, 900.0, -1100.0]
    seed.velocities_x[:, 0] = [3.0, -1.0, -2.5, 0.5]
    seed.velocities_y[:, 0] = [0.5, 2.0, -1.0, 4.0]
    seed.velocities_z[:, 0] = [1.0, 0.2, -0.3, -2.0]
    seed.E_kin_eV[:, 0] = _E_kin_eV(seed.mass_kg, seed.velocities_x[:, 0],
                                    seed.velocities_y[:, 0],
                                    seed.velocities_z[:, 0])
    return seed


def _pair_energy_eV(seed, cfg):
    x, y, z = (seed.positions_x[:, 0], seed.positions_y[:, 0],
               seed.positions_z[:, 0])
    n = x.size // 2
    r = np.sqrt((x[:n] - x[n:]) ** 2 + (y[:n] - y[n:]) ** 2 + (z[:n] - z[n:]) ** 2)
    return np.asarray(ion_interaction_potential(r, np.ones(n), np.ones(n), cfg))


def _v(res):
    return np.stack([res.vx_detected, res.vy_detected, res.vz_detected], 1)


def _five_term(seed, res):
    before = (seed.E_kin_eV[:, -1] + seed.E_pot_eV[:, -1] + seed.E_dissip_eV[:, -1]
              + seed.E_mass_transfer_eV[:, -1] + seed.E_int_eV[:, -1])
    after = (res.E_kin_detected_eV + res.E_pot_detected_eV
             + res.E_dissip_detected_eV + res.E_mass_transfer_detected_eV
             + res.E_int_detected_eV)
    return before, after


class TestStageClosure:
    def test_default_none_leaves_velocities_and_E_pot_verbatim(self):
        cfg = _cfg(detection_coulomb_closure="none")
        seed = _free_pair_seed(cfg, E_int_eV=0.0)          # frozen: no events
        res = run_detection_stage(seed, cfg)
        np.testing.assert_array_equal(res.vx_detected, seed.velocities_x[:, 0])
        np.testing.assert_array_equal(res.vz_detected, seed.velocities_z[:, 0])
        np.testing.assert_array_equal(res.E_pot_detected_eV, seed.E_pot_eV[:, 0])

    def test_both_free_pairs_take_the_two_body_asymptote(self):
        cfg = _cfg()
        seed = _free_pair_seed(cfg, E_int_eV=0.0)
        res = run_detection_stage(seed, cfg)
        pos = np.stack([seed.positions_x[:, 0], seed.positions_y[:, 0],
                        seed.positions_z[:, 0]], 1)
        vel = np.stack([seed.velocities_x[:, 0], seed.velocities_y[:, 0],
                        seed.velocities_z[:, 0]], 1)
        m = seed.mass_kg / U
        E_c = _pair_energy_eV(seed, cfg)
        r = np.linalg.norm(pos[:2] - pos[2:], axis=1)
        a1, a2 = two_body_coulomb_asymptote(
            pos[:2], vel[:2], m[:2], pos[2:], vel[2:], m[2:],
            _eV_to_amu_ang2_ps2(E_c * r))
        np.testing.assert_array_equal(_v(res), np.concatenate([a1, a2]))
        # D2 ledger: per-ion invariant exact; pair E_pot drops by E_c.
        before, after = _five_term(seed, res)
        np.testing.assert_allclose(after, before, rtol=0, atol=1e-12)
        dE_pot = res.E_pot_detected_eV - seed.E_pot_eV[:, 0]
        np.testing.assert_allclose(dE_pot[:2] + dE_pot[2:], -E_c, rtol=1e-9)
        # The momentum-fixed shares are NOT 1/2 each (unequal state).
        assert not np.allclose(dE_pot[:2], -0.5 * E_c, rtol=1e-3)

    def test_ledger_stays_exact_with_events_after_the_closure(self):
        cfg = _cfg(evap_rrk_dof=2.0)
        picture, kappa = _ladder(cfg)
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        seed = _free_pair_seed(cfg, n_shell=21, E_int_eV=0.98 * sigma21)
        res = run_detection_stage(seed, cfg)
        assert res.event_offsets[-1] > 0                    # non-vacuous
        before, after = _five_term(seed, res)
        np.testing.assert_allclose(after, before, rtol=0, atol=1e-11)

    def test_retained_partner_gives_full_pair_energy_to_the_free_ion(self):
        cfg = _cfg(detection_droplet_retained_policy="exclude")
        seed = _free_pair_seed(cfg, E_int_eV=0.0)
        # ion 2 (partner of 0): at rest at the droplet centre -> bound.
        seed.positions_x[2, 0] = seed.positions_y[2, 0] = seed.positions_z[2, 0] = 0.0
        seed.velocities_x[2, 0] = seed.velocities_y[2, 0] = seed.velocities_z[2, 0] = 0.0
        seed.E_kin_eV[2, 0] = 0.0
        res = run_detection_stage(seed, cfg)
        assert res.state_reason[2] == "droplet_retained"
        E_c = _pair_energy_eV(seed, cfg)
        dke = res.E_kin_detected_eV - seed.E_kin_eV[:, 0]
        np.testing.assert_allclose(dke[0], E_c[0], rtol=1e-9)
        # retained partner verbatim (state and ledger)
        assert res.vx_detected[2] == 0.0 and res.vz_detected[2] == 0.0
        assert res.E_pot_detected_eV[2] == seed.E_pot_eV[2, 0]
        # pair sum again exact: the free ion's E_pot carries the whole -E_c
        np.testing.assert_allclose(res.E_pot_detected_eV[0] - seed.E_pot_eV[0, 0],
                                   -E_c[0], rtol=1e-9)
        before, after = _five_term(seed, res)
        np.testing.assert_allclose(after, before, rtol=0, atol=1e-12)

    def test_retained_first_half_partner_closes_the_second_half_ion(self):
        # The mirror branch (~f1 & f2): the free ion sits in the [N:] half,
        # its retained partner in the [:N] half; molecule 1 stays two-body.
        cfg = _cfg(detection_droplet_retained_policy="exclude")
        seed = _free_pair_seed(cfg, E_int_eV=0.0)
        # ion 0 (partner of 2): at rest at the droplet centre -> bound.
        seed.positions_x[0, 0] = seed.positions_y[0, 0] = seed.positions_z[0, 0] = 0.0
        seed.velocities_x[0, 0] = seed.velocities_y[0, 0] = seed.velocities_z[0, 0] = 0.0
        seed.E_kin_eV[0, 0] = 0.0
        res = run_detection_stage(seed, cfg)
        assert res.state_reason[0] == "droplet_retained"
        pos = np.stack([seed.positions_x[:, 0], seed.positions_y[:, 0],
                        seed.positions_z[:, 0]], 1)
        vel = np.stack([seed.velocities_x[:, 0], seed.velocities_y[:, 0],
                        seed.velocities_z[:, 0]], 1)
        m = seed.mass_kg / U
        E_c = _pair_energy_eV(seed, cfg)
        r = np.linalg.norm(pos[:2] - pos[2:], axis=1)
        k = _eV_to_amu_ang2_ps2(E_c * r)
        want2 = fixed_centre_coulomb_asymptote(pos[[2]], vel[[2]], m[[2]],
                                               pos[[0]], k[[0]])
        np.testing.assert_array_equal(_v(res)[2], want2[0])
        a1, a3 = two_body_coulomb_asymptote(pos[[1]], vel[[1]], m[[1]],
                                            pos[[3]], vel[[3]], m[[3]], k[[1]])
        np.testing.assert_array_equal(_v(res)[[1, 3]], np.concatenate([a1, a3]))
        dke = res.E_kin_detected_eV - seed.E_kin_eV[:, 0]
        np.testing.assert_allclose(dke[2], E_c[0], rtol=1e-9)
        assert dke[0] == 0.0 and res.E_pot_detected_eV[0] == seed.E_pot_eV[0, 0]
        before, after = _five_term(seed, res)
        np.testing.assert_allclose(after, before, rtol=0, atol=1e-12)

    def test_ce_pair_scale_scales_the_closed_pair_energy(self):
        # Tier-2 (C) mixture: the closure's k must carry the per-molecule CE
        # scale s_m from the checkpoint, exactly as the MD pair force did.
        cfg = _cfg()
        seed = _free_pair_seed(cfg, E_int_eV=0.0)
        s_m = np.array([0.8, 1.6])
        seed.ce_channel = np.full(4, CE_CHANNEL_Q2, dtype=int)
        seed.ce_E_m_eV = np.tile(s_m * E_REF_PER_ION_EV, 2)
        res = run_detection_stage(seed, cfg)
        E_c_unscaled = _pair_energy_eV(seed, cfg)
        dke = res.E_kin_detected_eV - seed.E_kin_eV[:, 0]
        np.testing.assert_allclose(dke[:2] + dke[2:], s_m * E_c_unscaled, rtol=1e-9)
        before, after = _five_term(seed, res)
        np.testing.assert_allclose(after, before, rtol=0, atol=1e-12)

    def test_both_retained_pair_is_untouched(self):
        cfg = _cfg(detection_droplet_retained_policy="exclude_all_coupled")
        seed = _free_pair_seed(cfg, E_int_eV=0.0)
        # molecule 1 (ions 1, 3) trapped together inside, fast -> marginal pair
        seed.positions_x[[1, 3], 0] = [5.0, -5.0]
        seed.positions_y[[1, 3], 0] = 0.0
        seed.positions_z[[1, 3], 0] = 0.0
        res = run_detection_stage(seed, cfg)
        # Unbound on the full-credit split (E_c(10 A) ~ 1.4 eV >> E_bind) but
        # helium-coupled -> both marginal, i.e. D4's "trapped pair" shape.
        assert list(res.state_reason[[1, 3]]) == [ds.RETAINED_MARGINAL_REASON] * 2
        for i in (1, 3):
            assert res.vx_detected[i] == seed.velocities_x[i, 0]
            assert res.vy_detected[i] == seed.velocities_y[i, 0]
            assert res.E_pot_detected_eV[i] == seed.E_pot_eV[i, 0]

    def test_free_ion_off_the_well_plateau_refuses(self, monkeypatch):
        # D3: the closure omits the well climb by assertion. Only the closure
        # evaluates the well here (no violators -> escape_energetics unused).
        cfg = _cfg()
        seed = _free_pair_seed(cfg, E_int_eV=0.0)
        monkeypatch.setattr(
            ds, "droplet_potential",
            lambda d, steepness, binding_energy: np.full_like(
                np.asarray(d, dtype=float), binding_energy - 1e-3),
        )
        with pytest.raises(ValueError, match="well plateau"):
            run_detection_stage(seed, cfg)

    def test_unit_conversion_round_trips(self):
        x = np.array([1e-3, 14.39964548, 7.0e2])
        np.testing.assert_allclose(_amu_ang2_ps2_to_eV(_eV_to_amu_ang2_ps2(x)),
                                   x, rtol=1e-15)


class TestMarginalPartnerSafeguard:
    def _escaped_partner_seed(self, cfg):
        e_int = 0.5 * float(d0_of_n(21, picture=cfg.ladder_electronic_picture,
                                    kappa=cfg.ladder_steepness))
        seed = _far_seed(cfg, n_shell=21, E_int_eV=e_int, inside=True, speed=5.0)
        seed.positions_z[1:, 0] = 1.0e5        # ion 0's partner (2) escaped
        return seed

    @pytest.mark.parametrize("closure,shed", [("none", "cold"),
                                              ("partner_aware", "co_moving")])
    def test_marginal_ion_with_escaped_partner_refuses(self, closure, shed):
        cfg = _detect_cfg(seed=3,
                          detection_droplet_retained_policy="exclude_all_coupled",
                          detection_coulomb_closure=closure,
                          evaporation_shed_convention=shed)
        with pytest.raises(ValueError, match="marginal-partner safeguard"):
            run_detection_stage(self._escaped_partner_seed(cfg), cfg)

    def test_trapped_marginal_pair_passes(self):
        cfg = _detect_cfg(seed=3,
                          detection_droplet_retained_policy="exclude_all_coupled")
        seed = self._escaped_partner_seed(cfg)
        seed.positions_z[2, 0] = 0.0            # partner back inside
        res = run_detection_stage(seed, cfg)
        assert res.state_reason[0] == "droplet_retained_marginal"
        assert res.state_reason[2] in ds.RETAINED_REASONS


# ---------------------------------------------------------------------------
# Config surface
# ---------------------------------------------------------------------------
class TestClosureConfig:
    def test_default_is_none(self):
        assert _detect_cfg().detection_coulomb_closure == "none"

    def test_unknown_value_rejected(self):
        cfg = _detect_cfg()
        with pytest.raises(ValueError, match="detection_coulomb_closure"):
            check_detection_config(replace(cfg, detection_coulomb_closure="half"))

    def test_partner_aware_requires_co_moving(self):
        cfg = _detect_cfg()
        with pytest.raises(ValueError, match="co_moving"):
            check_detection_config(replace(
                cfg, detection_coulomb_closure="partner_aware",
                evaporation_shed_convention="cold"))
        check_detection_config(replace(
            cfg, detection_coulomb_closure="partner_aware",
            evaporation_shed_convention="co_moving"))

    def test_partner_aware_with_stage_disabled_refused(self):
        cfg = _detect_cfg()
        with pytest.raises(ValueError, match="silently inert"):
            check_detection_config(replace(
                cfg, detection_stage_enabled=False,
                detection_coulomb_closure="partner_aware",
                evaporation_shed_convention="co_moving"))
