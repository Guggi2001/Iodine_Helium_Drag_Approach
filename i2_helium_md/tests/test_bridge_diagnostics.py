"""Tests for the Phase-D bridge reconstruction helpers (Slice Z).

Covers ``postprocess/bridge_diagnostics.py`` — the two thin post-hoc
reconstructions (`t×`, ``Π(t)``) plus the mean-``n(t)`` reduction — against
the plan §2.2 / §3 oracles (`TIER2_PHASE_D_IMPLEMENTATION_PLAN.md`):

* ``crossing_time_ps`` returns the exact per-ion first-crossing time on a
  synthetic ramp, NaN where the gate never opens in-window, and matches the
  §2.2 analytic Newton-cooling closed form ``τ·ln(E_int(0)/Σ(21))`` to stored-
  step (dt) discreteness on a cooling-only synthetic.
* ``regime_parameter`` reproduces ``λ(n)·f_ret·τ`` on known inputs, with the
  ``ρ_ratio`` re-derivation matching a direct ``helium_density`` call on
  hand-built positions; ``Π → 0`` outside the droplet and ``Π = 0`` at
  ``n = n*`` (Langmuir cap); fail-loud on unresolvable config.
* ``mean_shell_count`` matches a hand average on a tiny synthetic ``n_shell``.
* A **few-step** generative driver smoke (the only live-driver check in
  pytest, plan §6): the 5-term invariant closes (reuse the Phase-C
  ``ion_ledger_closure`` gate) and all three helpers run on the real v7
  checkpoint — per-ion ``t×`` all-agree (pre-crossing dynamics are
  deterministic, §2.2).

No figures, no production-sized runs, no reference-data reads (plan §6).
Tolerances: exact identities at machine precision; the analytic-oracle bound
is one stored dt (sampling discreteness of a continuous crossing); the smoke
closure bound reuses the Phase-C 1e-2 Verlet-drift gate.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from i2_helium_md.config import SimConfig
from i2_helium_md.physics.dissociation_ladder import ladder_cumsum
from i2_helium_md.physics.helium_density import rho_he_ratio
from i2_helium_md.physics.pickup import lambda_attach
from i2_helium_md.postprocess.bridge_diagnostics import (
    crossing_time_ps,
    mean_shell_count,
    regime_parameter,
)

from tests.test_biphasic_step import _biphasic_driver_cfg, _tiny_neutral


PICTURE = "statistical_mixture"
KAPPA = 1.0


# =============================================================================
# mean_shell_count
# =============================================================================

class TestMeanShellCount:
    def test_matches_hand_average(self):
        n_shell = np.array([
            [21.0, 20.0, 19.0],
            [21.0, 21.0, 20.0],
            [21.0, 19.0, 18.0],
            [21.0, 20.0, 20.0],
        ])
        expected = np.array([21.0, 20.0, 19.25])  # column means by hand
        result = mean_shell_count(n_shell)
        np.testing.assert_array_equal(result, expected)

    def test_returns_1d_float_of_length_T(self):
        n_shell = np.full((6, 5), 21)
        result = mean_shell_count(n_shell)
        assert isinstance(result, np.ndarray)
        assert result.shape == (5,)
        assert result.dtype == np.float64

    def test_rejects_non_2d(self):
        with pytest.raises(ValueError, match="2-D"):
            mean_shell_count(np.array([21.0, 20.0]))


# =============================================================================
# crossing_time_ps
# =============================================================================

def _sigma21() -> float:
    return float(ladder_cumsum(21, picture=PICTURE, kappa=KAPPA))


class TestCrossingTimeRamp:
    """Exact first-crossing index per ion on synthetic E_int(t) ramps."""

    def test_exact_first_crossing_per_ion(self):
        time_ps = np.linspace(0.0, 1.0, 11)  # dt = 0.1
        sigma = _sigma21()
        n_shell = np.full((3, 11), 21.0)
        E = np.empty((3, 11))
        # Ion 0 drops below Sigma(21) first at index 3; ion 1 at index 7;
        # ion 2 stays strictly above throughout (never opens).
        E[0] = np.where(np.arange(11) >= 3, sigma - 0.01, sigma + 0.01)
        E[1] = np.where(np.arange(11) >= 7, sigma - 0.01, sigma + 0.01)
        E[2] = sigma + 0.01

        t_x = crossing_time_ps(
            E, n_shell, time_ps, picture=PICTURE, kappa=KAPPA
        )
        assert t_x.shape == (3,)
        assert t_x[0] == time_ps[3]
        assert t_x[1] == time_ps[7]
        assert np.isnan(t_x[2])  # gate never opens in-window -> NaN

    def test_crossing_is_strict_inequality(self):
        # E_int == Sigma exactly is still self-bound (t_x = min{t: E_int < Sigma}).
        time_ps = np.array([0.0, 0.1, 0.2])
        sigma = _sigma21()
        n_shell = np.full((1, 3), 21.0)
        E = np.array([[sigma + 0.01, sigma, sigma - 1e-9]])
        t_x = crossing_time_ps(
            E, n_shell, time_ps, picture=PICTURE, kappa=KAPPA
        )
        assert t_x[0] == time_ps[2]

    def test_already_open_at_t0(self):
        time_ps = np.array([0.0, 0.1])
        n_shell = np.full((1, 2), 21.0)
        E = np.full((1, 2), 0.5 * _sigma21())
        t_x = crossing_time_ps(
            E, n_shell, time_ps, picture=PICTURE, kappa=KAPPA
        )
        assert t_x[0] == time_ps[0]

    def test_threshold_follows_n_of_t(self):
        # The gate threshold is Sigma(n(t)), not a constant: an E_int that sits
        # between Sigma(14) and Sigma(21) is self-unbound at n = 14 but drops
        # below the threshold the moment occupancy grows to 21 (pickup raises
        # Sigma past it) -- the crossing is n-driven, at fixed E.
        time_ps = np.array([0.0, 0.1, 0.2, 0.3])
        sigma_21 = _sigma21()
        sigma_14 = float(ladder_cumsum(14, picture=PICTURE, kappa=KAPPA))
        assert sigma_14 < sigma_21  # self-check on the ladder
        n_shell = np.array([[14.0, 14.0, 21.0, 21.0]])
        E = np.full((1, 4), 0.5 * (sigma_14 + sigma_21))  # between the two
        t_x = crossing_time_ps(
            E, n_shell, time_ps, picture=PICTURE, kappa=KAPPA
        )
        assert t_x[0] == time_ps[2]

    def test_gate_onset_override_threads_to_threshold(self):
        # A diagnostic gate_onset_override_eV run replaces Sigma(n) with the
        # fixed override (evaporation._gate_threshold_eV); the reconstruction
        # must resolve the same gate or it silently misrepresents such a run.
        # E sits below Sigma(21) from t0 but crosses the (lower) override
        # only at index 2 -- the two gates give different crossings.
        time_ps = np.array([0.0, 0.1, 0.2])
        sigma = _sigma21()
        override_eV = 0.5 * sigma
        n_shell = np.full((1, 3), 21.0)
        E = np.array([[0.8 * sigma, 0.8 * sigma, 0.4 * sigma]])

        t_x_default = crossing_time_ps(
            E, n_shell, time_ps, picture=PICTURE, kappa=KAPPA
        )
        t_x_override = crossing_time_ps(
            E, n_shell, time_ps, picture=PICTURE, kappa=KAPPA,
            gate_onset_eV=override_eV,
        )
        assert t_x_default[0] == time_ps[0]   # below Sigma(21) immediately
        assert t_x_override[0] == time_ps[2]  # below the override only at t2

    def test_shape_mismatch_rejected(self):
        time_ps = np.linspace(0.0, 1.0, 5)
        with pytest.raises(ValueError, match="shape"):
            crossing_time_ps(
                np.zeros((2, 4)), np.full((2, 5), 21.0), time_ps,
                picture=PICTURE, kappa=KAPPA,
            )
        with pytest.raises(ValueError, match="shape"):
            crossing_time_ps(
                np.zeros((2, 5)), np.full((2, 5), 21.0), np.zeros(4),
                picture=PICTURE, kappa=KAPPA,
            )


class TestCrossingTimeAnalyticOracle:
    """Plan §2.2 closed form: pure Newton decay at fixed n = 21 crosses at
    t_x = tau*ln(E_int(0)/Sigma(21)), reconstructed to dt discreteness."""

    def test_newton_cooling_closed_form(self):
        tau_ps = 6.55
        e0_eV = 0.40           # f_int * E_avail at the pinned point (0.5 * 0.80)
        dt_ps = 0.01
        time_ps = np.arange(0.0, 8.0, dt_ps)
        sigma = _sigma21()
        assert e0_eV > sigma   # gate shut at t0 (self-check)

        E_row = e0_eV * np.exp(-time_ps / tau_ps)
        E = np.tile(E_row, (4, 1))
        n_shell = np.full((4, time_ps.size), 21.0)

        t_x = crossing_time_ps(
            E, n_shell, time_ps, picture=PICTURE, kappa=KAPPA
        )
        t_exact = tau_ps * np.log(e0_eV / sigma)
        # The continuous crossing lies in (t_k - dt, t_k]; the reconstruction
        # returns the first *sample* below Sigma, so it overshoots by < dt.
        assert np.all(t_x >= t_exact - 1e-12)
        assert np.all(t_x < t_exact + dt_ps + 1e-12)
        # Identical rows -> identical reconstruction (the all-agree wiring check).
        assert np.all(t_x == t_x[0])


# =============================================================================
# regime_parameter
# =============================================================================

def _stub_ckpt(positions_r, n_shell, droplet_radius=30.0):
    """Duck-typed checkpoint stub with only the fields regime_parameter reads.

    ``positions_r`` (2N, T) radial distances are laid along x (y = z = 0) so the
    intended depth r - R_droplet is exact by construction.
    """
    r = np.asarray(positions_r, dtype=float)
    two_N = r.shape[0]
    return SimpleNamespace(
        positions_x=r.copy(),
        positions_y=np.zeros_like(r),
        positions_z=np.zeros_like(r),
        droplet_radii_angstrom=np.full(two_N, float(droplet_radius)),
        n_shell=np.asarray(n_shell, dtype=float),
    )


def _bridge_cfg(**overrides):
    """Minimal cfg carrying the knobs regime_parameter reads (no validate())."""
    base = dict(
        mass_scenario="biphasic",
        internal_energy_partition_fraction=0.5,
        internal_energy_retained_fraction=0.1,
        internal_energy_cooling_tau_ps=6.55,
        pickup_rate_coefficient=0.9,
    )
    base.update(overrides)
    return SimConfig(**base)


class TestRegimeParameter:
    def test_reproduces_lambda_fret_tau_on_known_inputs(self):
        # Radii spanning deep-inside -> surface -> outside a 30 A droplet.
        positions_r = np.array([
            [5.0, 20.0, 45.0],
            [10.0, 30.0, 60.0],
        ])
        n_shell = np.array([
            [21.0, 20.0, 14.0],
            [21.0, 21.0, 18.0],
        ])
        ckpt = _stub_ckpt(positions_r, n_shell)
        cfg = _bridge_cfg()

        result = regime_parameter(ckpt, cfg)

        # The rho_ratio re-derivation must match a direct helium_density call on
        # the hand-built positions (depth = r - R_droplet), using the same
        # steepness surface the driver resolves (default drag gate collapses to
        # cfg.potential_steepness -- ion._drag_gate_steepness section-5.5).
        depth = positions_r - 30.0
        rho = rho_he_ratio(depth, steepness=cfg.potential_steepness)
        lam = lambda_attach(
            rho, n_shell,
            lambda0=cfg.pickup_rate_coefficient,
            p=cfg.pickup_occupancy_exponent,
            cap=cfg.pickup_occupancy_cap,
            pickup_rate_form=cfg.pickup_rate_form,
        )
        expected = (
            lam
            * cfg.internal_energy_retained_fraction
            * cfg.internal_energy_cooling_tau_ps
        )
        np.testing.assert_array_equal(result, expected)
        assert result.shape == positions_r.shape

    def test_zero_at_full_shell(self):
        # n = n* = 21 zeroes the Langmuir cap: Pi = 0 exactly at gate-open.
        ckpt = _stub_ckpt(np.full((3, 2), 5.0), np.full((3, 2), 21.0))
        result = regime_parameter(ckpt, _bridge_cfg())
        np.testing.assert_array_equal(result, np.zeros((3, 2)))

    def test_vanishes_outside_droplet(self):
        # Far outside (depth >> steepness) rho -> 0, so Pi -> 0: the exit-driven
        # termination route (plan §2.2).
        ckpt = _stub_ckpt(np.full((2, 2), 500.0), np.full((2, 2), 14.0))
        result = regime_parameter(ckpt, _bridge_cfg())
        np.testing.assert_allclose(result, 0.0, atol=1e-12)

    def test_freeze_side_at_pinned_point(self):
        # Plan §2.2 corrected oracle: at the pinned point even the deep-bulk,
        # n = 14 corner gives Pi ~ 0.9 * (1/3) * 0.1 * 6.55 ~ 0.2 < 1.
        ckpt = _stub_ckpt(np.full((1, 1), 5.0), np.full((1, 1), 14.0))
        result = regime_parameter(ckpt, _bridge_cfg())
        assert 0.0 < result[0, 0] < 1.0

    def test_resolves_independent_gate_steepness(self):
        # The section-2.2 contract is "the same resolved steepness the run
        # used": with drag_spatial_gate='erf_independent' the resolver returns
        # cfg.drag_gate_steepness, NOT cfg.potential_steepness (ion.
        # _drag_gate_steepness section-5.5). A divergent pair discriminates a
        # wrong hard-coding that the default (coincident) gate cannot.
        positions_r = np.array([[20.0, 28.0, 36.0]])  # near-surface: erf sensitive
        n_shell = np.full((1, 3), 14.0)
        ckpt = _stub_ckpt(positions_r, n_shell)
        cfg = _bridge_cfg(
            drag_spatial_gate="erf_independent", drag_gate_steepness=5.0
        )
        assert cfg.drag_gate_steepness != cfg.potential_steepness  # divergent pair

        result = regime_parameter(ckpt, cfg)

        depth = positions_r - 30.0
        expected_lam = lambda_attach(
            rho_he_ratio(depth, steepness=cfg.drag_gate_steepness), n_shell,
            lambda0=cfg.pickup_rate_coefficient,
            p=cfg.pickup_occupancy_exponent,
            cap=cfg.pickup_occupancy_cap,
            pickup_rate_form=cfg.pickup_rate_form,
        )
        expected = (
            expected_lam
            * cfg.internal_energy_retained_fraction
            * cfg.internal_energy_cooling_tau_ps
        )
        np.testing.assert_array_equal(result, expected)
        # ... and differs from the potential_steepness reconstruction, so the
        # assert above actually discriminates the two resolutions.
        wrong_rho = rho_he_ratio(depth, steepness=cfg.potential_steepness)
        assert not np.allclose(
            rho_he_ratio(depth, steepness=cfg.drag_gate_steepness), wrong_rho
        )

    def test_rejects_unset_f_ret(self):
        ckpt = _stub_ckpt(np.full((1, 1), 5.0), np.full((1, 1), 21.0))
        with pytest.raises(ValueError, match="internal_energy_retained_fraction"):
            regime_parameter(ckpt, _bridge_cfg(internal_energy_retained_fraction=None))

    def test_rejects_non_positive_tau(self):
        # tau <= 0 is unphysical (mirrors the check_solvation_cooling_config
        # load guard); a hand-built cfg bypassing validate() must still fail
        # loud here, not return a sign-flipped Pi.
        ckpt = _stub_ckpt(np.full((1, 1), 5.0), np.full((1, 1), 21.0))
        for bad_tau in (None, 0.0, -6.55):
            with pytest.raises(ValueError, match="internal_energy_cooling_tau_ps"):
                regime_parameter(
                    ckpt, _bridge_cfg(internal_energy_cooling_tau_ps=bad_tau)
                )

    def test_rejects_tabulated_density_profile(self):
        # The re-derivation assumes the analytic erf-complement gate; a
        # 'tabulated' run could not be reconstructed this way -- mirror the
        # driver's point-of-use refusal instead of silently substituting erf.
        ckpt = _stub_ckpt(np.full((1, 1), 5.0), np.full((1, 1), 21.0))
        with pytest.raises(NotImplementedError, match="helium_density_profile"):
            regime_parameter(ckpt, _bridge_cfg(helium_density_profile="tabulated"))


# =============================================================================
# Few-step generative driver smoke (the only live-driver check; plan §6)
# =============================================================================

class TestBridgeDriverSmoke:
    """Few-step biphasic run: 5-term invariant closes; helpers consume the
    real v7 checkpoint; per-ion t_x all-agree (deterministic pre-crossing)."""

    @pytest.fixture(scope="class")
    def smoke_run(self):
        from i2_helium_md.simulation.ion import run_ion_propagation

        # Fast-cooling override pulls t_x inside the 0.2 ps smoke window
        # (t_x = tau*ln(0.24/Sigma(21)) ~ 0.014 ps at tau = 0.05); every other
        # knob is the shared Phase-C smoke cfg.
        cfg = _biphasic_driver_cfg(internal_energy_cooling_tau_ps=0.05)
        ck = run_ion_propagation(cfg, _tiny_neutral())
        return cfg, ck

    def test_five_term_invariant_closes(self, smoke_run):
        from i2_helium_md.postprocess.energy_balance import ion_ledger_closure

        _, ck = smoke_run
        closure = ion_ledger_closure(ck)
        e0 = abs(closure.E_system_eV[0])
        assert closure.max_abs_residual_eV / e0 < 1e-2

    def test_mean_shell_count_on_real_checkpoint(self, smoke_run):
        _, ck = smoke_run
        mean_n = mean_shell_count(ck.n_shell)
        assert mean_n.shape == (ck.time_ps.size,)
        assert np.all(np.isfinite(mean_n))
        assert mean_n[0] == 21.0  # biphasic runs seed the full shell

    def test_crossing_time_per_ion_all_agree(self, smoke_run):
        cfg, ck = smoke_run
        t_x = crossing_time_ps(
            ck.E_int_eV, ck.n_shell, ck.time_ps,
            picture=cfg.ladder_electronic_picture,
            kappa=cfg.ladder_steepness,
        )
        assert t_x.shape == (2 * ck.num_molecules,)
        # Pre-crossing dynamics are deterministic and identical across ions
        # (same S2 onset, same cooling, pickup dead at n = n*): disagreement is
        # a wiring flag, not physics (plan §2.2).
        assert np.all(np.isfinite(t_x))
        assert np.all(t_x == t_x[0])

    def test_regime_parameter_on_real_checkpoint(self, smoke_run):
        cfg, ck = smoke_run
        pi = regime_parameter(ck, cfg)
        assert pi.shape == ck.n_shell.shape
        assert np.all(np.isfinite(pi))
        assert np.all(pi >= 0.0)
        # Column 0 is gate-open at the full shell: Pi = 0 exactly.
        np.testing.assert_array_equal(pi[:, 0], 0.0)
