"""Tests for the s(n) drag-shell state coupling (atlas §3.5i; design S1/S2/S4).

Covers the geometric closure module (`physics/state_coupling.py`), the
config-load guard (`check_drag_state_coupling_config`), and the driver seam
(BC-1..BC-3): the `off` default never builds a wrapper (structural
bit-identity — the *numerical* off-path regression is the entire existing
drag/biphasic test battery, which runs on the default and predates the seam),
a live `shell_area` coupling changes the dynamics in both the ion stage and
the E2 landau_gated_drag arm, and one BAOAB step under a live s matches the
analytic O-step damping exactly.

Registered numbers: the s-table asserted here is the design §3/§8 table
after the 2026-07-29 pre-build erratum (exact closure arithmetic).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from i2_helium_md.physics.constants import U
from i2_helium_md.physics.drag import drag_gamma
from i2_helium_md.physics.state_coupling import (
    N_REF_INTEGER_TOL,
    apply_state_factor,
    derive_n_ref_amu,
    effective_radius_angstrom,
    shell_area_state_factor,
)
from tests.test_biphasic_step import _biphasic_driver_cfg, _tiny_neutral

# The standing shared bundle's extraction-mass stamp (I+ + 19 He).
STANDING_STAMP_AMU = 202.953908

# Design §3 prior table (R_core 3.2, rho 0.030, n_ref 19; post-erratum).
PRIOR_S_TABLE = {21: 1.057, 19: 1.000, 14: 0.850, 8: 0.650,
                 5: 0.538, 2: 0.412, 1: 0.366, 0: 0.317}

_PRIOR_KWARGS = dict(R_core_angstrom=3.2, rho_shell_per_A3=0.030)


class TestGeometricClosure:
    def test_registered_prior_table_to_3_decimals(self):
        for n, want in PRIOR_S_TABLE.items():
            got = float(shell_area_state_factor(
                float(n), n_ref=19, **_PRIOR_KWARGS))
            assert round(got, 3) == want, f"s({n}) = {got}"

    def test_s_at_n_ref_is_exactly_one(self):
        # Identical arithmetic on both sides of the ratio -> exact, not approx.
        s = shell_area_state_factor(19.0, n_ref=19, **_PRIOR_KWARGS)
        assert float(s) == 1.0

    def test_monotone_increasing_in_n_and_positive(self):
        n = np.arange(0, 30, dtype=float)
        s = shell_area_state_factor(n, n_ref=19, **_PRIOR_KWARGS)
        assert np.all(s > 0)
        assert np.all(np.diff(s) > 0)

    def test_bare_radius_is_the_core_radius(self):
        r0 = effective_radius_angstrom(0.0, **_PRIOR_KWARGS)
        # (R_core**3)**(1/3) round-trips through floating cube/cube-root.
        np.testing.assert_allclose(float(r0), 3.2, rtol=1e-12)

    def test_units_shell_volume_term(self):
        # One He at shell density rho adds exactly 1/rho of volume [A^3].
        r1 = effective_radius_angstrom(1.0, **_PRIOR_KWARGS)
        vol_added = (4.0 / 3.0) * np.pi * (float(r1) ** 3 - 3.2 ** 3)
        np.testing.assert_allclose(vol_added, 1.0 / 0.030, rtol=1e-12)

    @pytest.mark.parametrize("bad", [
        dict(R_core_angstrom=0.0, rho_shell_per_A3=0.030),
        dict(R_core_angstrom=-3.2, rho_shell_per_A3=0.030),
        dict(R_core_angstrom=3.2, rho_shell_per_A3=0.0),
        dict(R_core_angstrom=3.2, rho_shell_per_A3=np.inf),
    ])
    def test_bad_parameters_raise(self, bad):
        with pytest.raises(ValueError):
            effective_radius_angstrom(1.0, **bad)

    def test_negative_n_raises(self):
        with pytest.raises(ValueError, match="n must be >= 0"):
            shell_area_state_factor(
                np.array([1.0, -0.5]), n_ref=19, **_PRIOR_KWARGS)

    def test_bad_n_ref_raises(self):
        with pytest.raises(ValueError, match="n_ref"):
            shell_area_state_factor(1.0, n_ref=0, **_PRIOR_KWARGS)


class TestDeriveNRef:
    def test_standing_stamp_derives_19(self):
        assert derive_n_ref_amu(STANDING_STAMP_AMU) == 19

    def test_non_integer_stamp_raises(self):
        # Half a He off any integer shell state (way past the tolerance).
        from i2_helium_md.physics.constants import MASS_HE_AMU, MASS_I_ION_AMU
        bad = MASS_I_ION_AMU + 18.5 * MASS_HE_AMU
        with pytest.raises(ValueError, match="n_ref"):
            derive_n_ref_amu(bad)

    def test_non_positive_stamp_raises(self):
        with pytest.raises(ValueError, match="positive"):
            derive_n_ref_amu(0.0)

    def test_tolerance_is_tight(self):
        # The constants.py mass-table slack (~0.0011) must pass, but the guard
        # must stay an order tighter than half a He.
        assert N_REF_INTEGER_TOL < 0.5 / 4.0
        assert N_REF_INTEGER_TOL >= 0.002


class TestApplyStateFactor:
    def test_elementwise_scaling(self):
        s = np.array([0.5, 1.0, 2.0])

        def base(speed, depth):
            return np.asarray(speed) * 10.0

        wrapped = apply_state_factor(base, s)
        speed = np.array([1.0, 2.0, 3.0])
        depth = np.zeros(3)
        np.testing.assert_allclose(
            wrapped(speed, depth), s * base(speed, depth), rtol=0.0, atol=0.0)


class TestConfigGuard:
    def test_default_off_validates(self):
        cfg = _biphasic_driver_cfg()
        assert cfg.drag_state_coupling == "off"
        cfg.validate()

    def test_shell_area_on_biphasic_validates(self):
        replace(_biphasic_driver_cfg(), drag_state_coupling="shell_area").validate()

    def test_unknown_selector_rejected(self):
        with pytest.raises(ValueError, match="unknown drag_state_coupling"):
            replace(_biphasic_driver_cfg(),
                    drag_state_coupling="shell_are").validate()

    def test_shell_area_outside_biphasic_rejected(self):
        from i2_helium_md.presets import single_pulse_N2000_drag
        cfg = single_pulse_N2000_drag(num_molecules=2)
        with pytest.raises(ValueError, match="mass_scenario='biphasic'"):
            replace(cfg, drag_state_coupling="shell_area").validate()

    @pytest.mark.parametrize("field,value", [
        ("state_coupling_R_core_angstrom", 0.0),
        ("state_coupling_R_core_angstrom", -1.0),
        ("state_coupling_rho_shell_per_A3", 0.0),
        ("state_coupling_rho_shell_per_A3", float("nan")),
    ])
    def test_bad_coefficients_rejected_even_inert(self, field, value):
        # Rule 2: fail-early on unphysical geometry even under "off".
        with pytest.raises(ValueError, match="finite and > 0"):
            replace(_biphasic_driver_cfg(), **{field: value}).validate()

    def test_shell_area_without_bundle_rejected(self):
        cfg = replace(
            _biphasic_driver_cfg(),
            drag_state_coupling="shell_area", drag_coefficients=None,
        )
        with pytest.raises(ValueError, match="drag_coefficients"):
            cfg.validate()

    def test_n_ref_stamp_consistency_checked_at_load(self):
        cfg = _biphasic_driver_cfg()
        bad_bundle = replace(cfg.drag_coefficients,
                             extraction_mass_amu=200.0)  # not I+ + k*He
        with pytest.raises(ValueError, match="n_ref"):
            replace(cfg, drag_state_coupling="shell_area",
                    drag_coefficients=bad_bundle).validate()


class TestOneStepDeterministic:
    """One BAOAB step with a live s: the O-step damping is exactly
    exp(-s*gamma*dt/m) (zero conservative forces, so B/A are inert on v)."""

    def test_o_step_decay_carries_s(self):
        from i2_helium_md.physics.baoab import make_ion_baoab_step

        two_n = 4
        m_amu = np.full(two_n, 202.953908)
        radii = np.full(two_n, 30.0)
        s = np.array([0.32, 0.65, 1.0, 1.06])
        cfg = _biphasic_driver_cfg()
        coeffs = cfg.drag_coefficients
        steep = cfg.potential_steepness

        def zero_acc(pos):
            x, y, z = pos
            return (np.zeros_like(x), np.zeros_like(y), np.zeros_like(z)), \
                np.zeros(two_n // 2)

        def base_gamma(speed, depth):
            return drag_gamma(speed, depth, coeffs=coeffs, steepness=steep)

        pos = (np.zeros(two_n), np.zeros(two_n), np.zeros(two_n))  # center: deep
        v0 = 3.0
        vel = (np.full(two_n, v0), np.zeros(two_n), np.zeros(two_n))
        dt = 0.01

        step = make_ion_baoab_step(
            m_amu, radii, zero_acc, apply_state_factor(base_gamma, s),
            T_eff=0.0,
        )
        _, (vx, vy, vz), _, _ = step(pos, vel, dt)

        # BAOAB = B A O A B: with zero acceleration the O-step sees the
        # unchanged speed v0 but the position after the first half-drift,
        # so the gate depth is |v0*dt/2| - R, not -R.
        depth = np.abs(v0 * dt / 2.0) - radii
        expected_decay = np.exp(-s * base_gamma(np.full(two_n, v0), depth)
                                * dt / m_amu)
        np.testing.assert_allclose(vx, v0 * expected_decay, rtol=1e-12)


class TestDriverSeam:
    def test_off_never_builds_a_wrapper(self, monkeypatch):
        """Structural bit-identity of the off path: apply_state_factor is
        unreachable, so the base gamma_fn partial goes to BAOAB verbatim."""
        import i2_helium_md.simulation.ion as ion_mod

        def _boom(*a, **k):  # pragma: no cover - failure path
            raise AssertionError("apply_state_factor called under 'off'")

        monkeypatch.setattr(ion_mod, "apply_state_factor", _boom)
        from i2_helium_md.simulation.ion import run_ion_propagation
        run_ion_propagation(_biphasic_driver_cfg(), _tiny_neutral())

    def test_shell_area_changes_ion_stage_dynamics(self):
        """Strong coupling (bulk rho) must move the trajectories off the
        s == 1 baseline while preserving finiteness."""
        from i2_helium_md.simulation.ion import run_ion_propagation

        base = run_ion_propagation(_biphasic_driver_cfg(), _tiny_neutral())
        coupled_cfg = replace(
            _biphasic_driver_cfg(),
            drag_state_coupling="shell_area",
            state_coupling_rho_shell_per_A3=0.0218,
        )
        coupled_cfg.validate()
        coupled = run_ion_propagation(coupled_cfg, _tiny_neutral())
        assert np.all(np.isfinite(coupled.velocities_final_x))
        assert not np.allclose(
            base.velocities_final_x, coupled.velocities_final_x)

    def test_shell_area_scales_e2_landau_arm(self):
        """BC-3: the E2 landau_gated_drag damping carries s. Same seed and
        identical (translation-independent) shed sequence, so any difference
        is the drag scaling."""
        from tests.test_relaxation_stage import _relax_cfg, _seed_checkpoint
        from i2_helium_md.simulation.relaxation_stage import (
            run_relaxation_stage,
        )

        def _run(**coupling):
            cfg = _relax_cfg(
                relaxation_dissipation="landau_gated_drag",
                relaxation_time_ps=0.5,
                **coupling,
            )
            cfg.validate()
            seed = _seed_checkpoint(cfg, n_shell=5, E_int_eV=0.0)
            # Super-Landau speeds (v_limit ~ 0.58 A/ps) so the gated drag arm
            # is live; E_int = 0 freezes the mass subsystem (n stays 5, well
            # off n_ref = 19, so s != 1 on the coupled leg).
            for arr in (seed.velocities_x,):
                arr[:] = 2.0
            return run_relaxation_stage(seed, cfg)

        base = _run()
        coupled = _run(drag_state_coupling="shell_area",
                       state_coupling_rho_shell_per_A3=0.0218)
        vb = base.checkpoint.velocities_final_x
        vc = coupled.checkpoint.velocities_final_x
        assert np.all(np.isfinite(vc))
        assert not np.allclose(vb, vc)
        # n = 5 -> s < 1 -> weaker drag -> the coupled ions stay faster.
        assert np.all(np.abs(vc) >= np.abs(vb))
