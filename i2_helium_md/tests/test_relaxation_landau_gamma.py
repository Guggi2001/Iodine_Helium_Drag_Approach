"""E2 Landau-gated drag arm: gamma gate + integration (plan §I.11.2 item 2, arm (c)).

The arm replaces the relaxation coulomb-translate's ``_zero_gamma`` with a speed-
gated friction coefficient::

    gamma(v, depth) = 0                            for v <= v_L   (superfluid)
    gamma(v, depth) = g(depth) * b * v**2          for v >  v_L   (locked pure cubic)

with ``v_L = cfg.v_limit_angstrom_per_ps`` (the existing Landau cutoff) and ``b``
drawn from the config's production drag bundle (single source, CLAUDE.md rule 1).
The gate is on **speed** (mass-independent Landau criterion; preserves the drag
module's mass-agnostic contract). Below the cutoff the arm is byte-identical to
the delivered zero-gamma closure; above it, drag threads dissipation into
``E_dissip`` through the same BAOAB machinery that closes the ion-stage ledger, so
the 5-term invariant closes.
"""

from dataclasses import replace

import numpy as np

from i2_helium_md.physics.drag import drag_gamma
from i2_helium_md.postprocess.energy_balance import ion_ledger_closure
from i2_helium_md.simulation.ion import drag_gate_steepness, run_ion_propagation
from i2_helium_md.simulation.relaxation_stage import (
    _make_relaxation_gamma_fn,
    _zero_gamma,
    run_relaxation_stage,
)

from tests.test_biphasic_step import _biphasic_driver_cfg, _tiny_neutral
from tests.test_relaxation_stage import _relax_cfg, _seed_checkpoint

_DEEP = -20.0   # depth (r_atom - r_droplet) deep inside the droplet -> g ~ 1


def _landau_cfg(**overrides):
    return _relax_cfg(relaxation_dissipation="landau_gated_drag",
                      relaxation_forces="coulomb", **overrides)


def _seed_at_speed(cfg, speed, **kw):
    """A `_seed_checkpoint` with the seed velocity set to `speed` along x."""
    ion = _seed_checkpoint(cfg, **kw)
    ion.velocities_x[:] = float(speed)
    ion.velocities_y[:] = 0.0
    ion.velocities_z[:] = 0.0
    return ion


# ---------------------------------------------------------------------------
# gamma gate: formula / shape / spatial gate  (Tier 0-1 of the 7-step order)
# ---------------------------------------------------------------------------
def test_zero_gamma_arm_is_the_delivered_closure():
    # Byte-inert: the default arm reuses the *same* function object, so the
    # translate path is literally unchanged.
    cfg = _relax_cfg()
    fn = _make_relaxation_gamma_fn(cfg, drag_gate_steepness(cfg))
    assert fn is _zero_gamma


def test_landau_below_cutoff_is_exactly_zero():
    cfg = _landau_cfg()
    fn = _make_relaxation_gamma_fn(cfg, drag_gate_steepness(cfg))
    speed = np.array([0.5 * cfg.v_limit_angstrom_per_ps])
    assert fn(speed, np.array([_DEEP]))[0] == 0.0


def test_landau_above_cutoff_matches_production_drag_gamma():
    cfg = _landau_cfg()
    s = drag_gate_steepness(cfg)
    fn = _make_relaxation_gamma_fn(cfg, s)
    speed = np.array([3.0])   # > v_L (0.4 A/ps)
    expected = drag_gamma(speed, np.array([_DEEP]),
                          coeffs=cfg.drag_coefficients, steepness=s)
    np.testing.assert_array_equal(fn(speed, np.array([_DEEP])), expected)


def test_landau_boundary_belongs_to_the_superfluid_side():
    # Heaviside is strict >: exactly at v_L there is no friction.
    cfg = _landau_cfg()
    fn = _make_relaxation_gamma_fn(cfg, drag_gate_steepness(cfg))
    v_L = cfg.v_limit_angstrom_per_ps
    assert fn(np.array([v_L]), np.array([_DEEP]))[0] == 0.0


def test_landau_spatial_gate_kills_drag_outside_droplet():
    # Above v_L but well outside the droplet: the erf gate g -> 0 removes drag.
    cfg = _landau_cfg()
    fn = _make_relaxation_gamma_fn(cfg, drag_gate_steepness(cfg))
    assert fn(np.array([5.0]), np.array([+200.0]))[0] < 1e-9


def test_landau_is_elementwise_over_the_array():
    cfg = _landau_cfg()
    s = drag_gate_steepness(cfg)
    fn = _make_relaxation_gamma_fn(cfg, s)
    v_L = cfg.v_limit_angstrom_per_ps
    depth = np.full(3, _DEEP)
    out = fn(np.array([0.1, v_L, 3.0]), depth)
    assert out[0] == 0.0
    assert out[1] == 0.0
    expected_hi = drag_gamma(np.array([3.0]), np.array([_DEEP]),
                             coeffs=cfg.drag_coefficients, steepness=s)[0]
    assert out[2] == expected_hi


# ---------------------------------------------------------------------------
# integration: deterministic multi-step + energy bookkeeping
# ---------------------------------------------------------------------------
def test_below_cutoff_run_matches_zero_gamma():
    # One step from a sub-cutoff seed: with all speeds <= v_L the Landau arm must
    # reproduce the zero-gamma closure bit-for-bit (the gate never opens).
    cfg0 = _relax_cfg(relaxation_time_ps=0.01)          # one dt_ion step
    ck0 = run_relaxation_stage(_seed_at_speed(cfg0, 0.2), cfg0).checkpoint
    ckL = run_relaxation_stage(
        _seed_at_speed(cfg0, 0.2),
        replace(cfg0, relaxation_dissipation="landau_gated_drag"),
    ).checkpoint

    # Precondition: the window really stayed sub-cutoff (else the test is vacuous).
    v = np.sqrt(ckL.velocities_x**2 + ckL.velocities_y**2 + ckL.velocities_z**2)
    assert np.all(v <= cfg0.v_limit_angstrom_per_ps)

    np.testing.assert_array_equal(ck0.E_kin_eV, ckL.E_kin_eV)
    np.testing.assert_array_equal(ck0.E_dissip_eV, ckL.E_dissip_eV)
    np.testing.assert_array_equal(ck0.velocities_x, ckL.velocities_x)


def test_above_cutoff_dissipates_and_isolates_the_drag_channel():
    # A supercritical seed: the mass subsystem (n, E_int) is drag-independent, so
    # it is identical between the arms; only the drag channel differs -> the
    # Landau run books strictly more E_dissip and ends colder.
    cfg0 = _relax_cfg(relaxation_time_ps=1.0)
    ck0 = run_relaxation_stage(_seed_at_speed(cfg0, 3.0), cfg0).checkpoint
    ckL = run_relaxation_stage(
        _seed_at_speed(cfg0, 3.0),
        replace(cfg0, relaxation_dissipation="landau_gated_drag"),
    ).checkpoint

    np.testing.assert_array_equal(ck0.n_shell, ckL.n_shell)
    np.testing.assert_array_equal(ck0.E_int_eV, ckL.E_int_eV)
    assert np.all(ckL.E_dissip_eV[:, -1] > ck0.E_dissip_eV[:, -1])   # drag fired
    assert np.all(ckL.E_kin_eV[:, -1] < ck0.E_kin_eV[:, -1])         # braked


def test_five_term_invariant_closes_under_landau_arm():
    # The correctness gate: with drag engaged the ledger still closes (the BAOAB
    # step threads drag dE_dissip into E_dissip, as the ion stage does).
    cfg = _biphasic_driver_cfg()
    ion = run_ion_propagation(cfg, _tiny_neutral())
    cfg = replace(cfg, relaxation_stage_enabled=True, relaxation_time_ps=3.0,
                  relaxation_forces="coulomb",
                  relaxation_dissipation="landau_gated_drag",
                  evap_rrk_dof=2.0, internal_energy_cooling_tau_ps=1.0)
    ck = run_relaxation_stage(ion, cfg).checkpoint
    closure = ion_ledger_closure(ck)
    e0 = abs(closure.E_system_eV[0])
    assert closure.max_abs_residual_eV / e0 < 1e-4


# ---------------------------------------------------------------------------
# Finding-1: the freeze-all early-exit must not truncate the drag arm for the
# droplet-retained class it exists to relax (evaporation freeze != translational
# equilibrium once drag is on).
# ---------------------------------------------------------------------------
def test_relaxation_converged_requires_sub_landau_under_the_drag_arm():
    from i2_helium_md.physics.dissociation_ladder import resolve_ladder
    from i2_helium_md.simulation.ion_propagation_step import (
        ion_state_from_checkpoint_column,
    )
    from i2_helium_md.simulation.relaxation_stage import (
        _freeze_mask,
        _relaxation_converged,
    )

    cfg = _landau_cfg()
    ladder = resolve_ladder(cfg.dissociation_ladder, cfg.tabulated_ladder_rungs_eV)
    # E_int = 0 -> evaporation-frozen; speed 5 A/ps -> super-Landau (v_L = 0.4).
    st = ion_state_from_checkpoint_column(
        _seed_at_speed(cfg, 5.0, E_int_eV=0.0), -1,
    )
    frozen = _freeze_mask(st, picture=cfg.ladder_electronic_picture,
                          kappa=cfg.ladder_steepness, ladder=ladder)
    assert bool(frozen.all())   # the seed really is evaporation-frozen

    v_L = cfg.v_limit_angstrom_per_ps
    # zero_gamma: the evaporation freeze alone ends the loop (delivered behaviour).
    assert _relaxation_converged(st, frozen, dissipation="zero_gamma", v_limit=v_L)
    # landau_gated_drag: still super-Landau -> NOT converged (drag must keep going).
    assert not _relaxation_converged(
        st, frozen, dissipation="landau_gated_drag", v_limit=v_L,
    )
    # ...and once damped sub-Landau, the drag arm does converge.
    st_slow = replace(st, vx=np.zeros_like(st.vx), vy=np.zeros_like(st.vy),
                      vz=np.zeros_like(st.vz))
    assert _relaxation_converged(
        st_slow, frozen, dissipation="landau_gated_drag", v_limit=v_L,
    )


def test_landau_arm_damps_past_the_evaporation_freeze():
    # End-to-end: a seed evaporation-frozen from t0 (E_int=0) but super-Landau.
    # The freeze-all early-exit alone would stop the loop at the first step; the
    # drag arm must keep braking past it, so it runs many more steps than the
    # zero_gamma arm (which correctly stops at the freeze).
    cfg = _landau_cfg(relaxation_time_ps=5.0)
    dt = cfg.dt_ion if cfg.relaxation_dt_ps is None else cfg.relaxation_dt_ps

    resL = run_relaxation_stage(_seed_at_speed(cfg, 5.0, E_int_eV=0.0), cfg)
    res0 = run_relaxation_stage(
        _seed_at_speed(cfg, 5.0, E_int_eV=0.0),
        replace(cfg, relaxation_dissipation="zero_gamma"),
    )

    # zero_gamma stops at the freeze (~one step); the drag arm runs well past it.
    assert res0.time_relaxed_ps <= 2.0 * dt
    assert resL.time_relaxed_ps > 5.0 * dt
    assert resL.time_relaxed_ps > res0.time_relaxed_ps
    # Drag fired and braked the retained ion below the seed speed.
    ckL = resL.checkpoint
    assert np.all(ckL.E_dissip_eV[:, -1] > 0.0)
    vT = np.sqrt(ckL.velocities_x[:, -1] ** 2 + ckL.velocities_y[:, -1] ** 2
                 + ckL.velocities_z[:, -1] ** 2)
    assert np.all(vT < 5.0)
