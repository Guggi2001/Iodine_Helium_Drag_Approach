"""Post-ejection relaxation stage (Tier-2 Phase E, Slice E2).

Propagates the biphasic **mass subsystem** past the 20 ps ion stage under
fixed-dt integration -- the only correct integrator **while K2 cooling is
live** (cooling makes ``E_int`` decay continuously between sheds, which the
event-driven detection stage structurally cannot represent).

Contract (re-framed at Slice DS -- design §1 item 5; supersedes the original
R5 "reach the terminal by matched-time integration" mission, which Wave 5
showed is unachievable on the gated arm: the terminal is flight-time
dependent there and the physical ~8.5e6 ps flight is beyond any tractable
fixed-dt cap):

* **Ungated (``cooling_spatial_gate == "none"``): the terminal solver.**
  Cooling acts everywhere and reliably drives ``E_int`` below ``D_0(n)``, so
  the stage runs to freeze-out (the early exit) and its relaxed read IS the
  converged terminal; the detection stage then no-ops per ion.
* **Gated (``"density_scaled"``): a handover bridge.** The stage carries the
  ions from sim-end to a P3-clean state (ejected, cooling erfc-suppressed);
  the *detection stage* (``simulation/detection_stage.py``, Slice DS) owns
  the microsecond tail exactly. The cap can therefore be short (~10 ps
  class); the DS P1-P3 handover guard fails loud if it is shortened past
  ejection/decoupling. Gated runs may also skip this stage entirely
  (``relaxation_stage_enabled=False`` -> DS seeds from ``ion.npz``).

The **relaxed read is a convergence diagnostic** (with ``frac_frozen``), not
the Tier-2 arbitration observable -- that is the Slice-DS *detected* read at
the Sourced ``t_detect`` (TIER2_DETECTION_STAGE_DESIGN.md §3.5).

Composition (decision R1 -- reuse, not re-composition)
------------------------------------------------------
The stage is **not** new physics. Per relaxation step it calls the delivered
:func:`~i2_helium_md.simulation.ion_propagation_step.biphasic_step` **verbatim**
under a *relaxation view* of the run config
(``replace(cfg, pickup_rate_coefficient=0.0, dt_ion=dt_relax)``, constructed
internally and **not** re-``validate()``-d -- the lambda_0 = 0 advisory in
``config.check_biphasic_config`` is the sanctioned evaporation-only limit). That
runs the locked K2 -> gate -> RRK -> K1 sequence; the Poisson pickup channel runs
too but is **structurally inert** (lambda_0 = 0 => P_attach == 0) while its RNG
draw is still consumed, so the Slice-X frozen two-draw stream is preserved
byte-identically.

Translation is drag-off (a **zero-gamma** BAOAB closure: decay = 1,
dE_dissip == 0, the delivered bookkeeping intact) with two arms:

* ``"coulomb"`` (default) -- the delivered ``make_ion_accel_fn`` conservative
  field (residual inter-ion Coulomb + droplet solvation, the latter vanishing
  outside the droplet), for ledger/asymptotic-state fidelity. Mirrors the ion
  driver's biphasic branch exactly, minus drag.
* ``"free_flight"`` -- ballistic positions (zero conservative acceleration), the
  MD potential *held* between sheds. Exact for the observable, not an
  approximation of it (the decoupling fact below).

Both arms fold the same ``e_bind_pair(n)`` term into ``E_pot`` (rule-1 reuse), so
a shed shifts ``E_pot`` by exactly ``+D_0(n)`` and the 5-term invariant closes in
both.

Decoupling fact (audit)
-----------------------
With lambda_0 = 0 and gamma = 0 the mass subsystem (K2 -> gate -> RRK -> K1) is
**independent of translation** -- positions enter only through the (inert at
lambda_0 = 0) pickup density gate. So ``coulomb`` and ``free_flight`` produce the
**identical shed sequence** ``n(t)`` / ``E_int(t)`` under the same seed; only the
translational split (x, v, E_kin, E_pot) differs. Translation is retained solely
for ledger and asymptotic-state fidelity.

Termination
-----------
Runs to ``relaxation_time_ps`` **or** stops early when every fragment is frozen:
``n_i == 0`` or ``E_int_i < D_0(n_i)`` (then the evaporation rate is
``k == 0`` forever, since ``E_int`` is monotone non-increasing under K2 + K1 with
pickup off -- verified against ``physics.evaporation.rrk_rate``). Each step
consumes both draws before the freeze check, so a cold input still consumes the
frozen two-draw stream on its terminating step.

The early-exit is **guaranteed** only while K2 actively drains ``E_int`` toward 0
(the ``cooling_spatial_gate == "none"`` default): it is K2 -- not evaporation --
that reliably pushes ``E_int`` below ``D_0(n)``. Under
``cooling_spatial_gate == "density_scaled"`` an ejected fragment (``rho_He -> 0``)
stops cooling, so ``E_int`` is merely held (monotone non-increasing, but not
driven to 0) and can sit above ``D_0(n)`` indefinitely -- evaporation's RRK rate
vanishes near threshold. Freeze is then not guaranteed and only
``relaxation_time_ps`` bounds the run; the ``config.check_relaxation_config``
step-budget guard rejects a runaway cap for that arm, and ``frac_frozen`` reports
per-run freeze completeness.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from functools import partial
from pathlib import Path
from typing import Callable, Optional

import numpy as np

from ..config import SimConfig, check_relaxation_config
from ..physics.constants import U
from ..physics.baoab import make_ion_baoab_step
from ..physics.dissociation_ladder import d0_of_n, resolve_ladder
from ..physics.drag import drag_gamma
from ..physics.leapfrog import make_ion_accel_fn
from ..physics.solvation_cooling import e_bind_pair_eV
from .checkpoint import (
    IonCheckpoint,
    check_biphasic_seed_checkpoint,
    save_ion_checkpoint,
    stage_stream_rng,
)
from .ion import DEFAULT_MAX_CHECKPOINT_BYTES_ION, _decide_stride_ion, drag_gate_steepness
from .ion_propagation_step import (
    IonStepState,
    baoab_propagation_step,
    biphasic_step,
    ion_state_from_checkpoint_column,
    write_ion_state_to_checkpoint_column,
)

# Fixed module-level key that spawns the relaxation stage's own PCG64 stream via
# SeedSequence((cfg.seed, RELAXATION_STREAM_KEY)) -- the ion-stage stream is never
# touched (the Slice-X freeze is *extended by a new stage*, never re-ordered).
RELAXATION_STREAM_KEY: int = 0xE2_2026


@dataclass(frozen=True)
class RelaxationResult:
    """Result of one relaxation-stage run.

    Attributes
    ----------
    checkpoint : IonCheckpoint
        A bona fide v7 ``IonCheckpoint`` covering the relaxation window (column 0
        = the seed = the ion stage's final column). Round-trips through
        ``load_ion_checkpoint``; ``ion_ledger_closure`` applies unchanged.
    freeze_flags : np.ndarray, shape (2N,), bool
        Per-ion frozen status at the final relaxed state (``n == 0`` or
        ``E_int < D_0(n)``).
    time_relaxed_ps : float
        **Window-relative** duration actually relaxed: the realized full window
        (``ceil(relaxation_time_ps / dt_relax) * dt_relax``, which exceeds
        ``relaxation_time_ps`` by less than one ``dt_relax``) or the earlier
        all-frozen time. Compare against ``cfg.relaxation_time_ps`` to detect an
        early freeze. The checkpoint's ``time_ps`` axis stays **absolute** (it
        continues the ion-stage axis from the seed column).
    terminal_n : np.ndarray, shape (2N,)
        The matched-time terminal shell count -- the E1 relaxed-input hook
        (``compute_terminal_shell_distribution`` duck-types on this attribute and
        reports ``source="relaxed"``).
    """

    checkpoint: IonCheckpoint
    freeze_flags: np.ndarray
    time_relaxed_ps: float
    terminal_n: np.ndarray


def _zero_gamma(speed: np.ndarray, depth: np.ndarray) -> np.ndarray:
    """Zero friction coefficient -> BAOAB O-step decay = 1, dE_dissip = 0 (drag off)."""
    return np.zeros_like(np.asarray(speed, dtype=float))


def _make_relaxation_gamma_fn(
    cfg: SimConfig, gate_steepness: float,
) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    """Resolve the coulomb-translate friction coefficient for the dissipation arm.

    ``zero_gamma`` (default) returns the delivered :func:`_zero_gamma` closure
    **verbatim** -- same function object, so the translate path is byte-identical.

    ``landau_gated_drag`` (plan §I.11.2 item 2, arm (c)) returns a **speed-gated**
    coefficient: ``gamma = 0`` for ``speed <= cfg.v_limit_angstrom_per_ps`` (the
    Landau cutoff -- sub-Landau superfluid motion is dissipationless) and the
    locked pure-cubic drag ``g(depth)*b*speed**2`` above it. ``b`` and the erf
    spatial gate come from the config's production drag bundle via
    :func:`~i2_helium_md.physics.drag.drag_gamma` (single source, CLAUDE.md rule
    1), so the arm never re-derives the coefficient. The gate is on **speed**, not
    kinetic energy: the Landau critical velocity is mass-independent, and drag.py
    is mass-agnostic by contract. Heaviside is strict ``>`` -- exactly at the
    cutoff there is no friction.

    Parameters
    ----------
    cfg : SimConfig
        Carries ``relaxation_dissipation``, the drag bundle, and the cutoff. The
        ``landau_gated_drag`` arm requires ``drag_coefficients is not None`` (the
        biphasic scenario guarantees it; validated at config-load).
    gate_steepness : float
        The drag erf-gate width [Angstrom], resolved once by
        :func:`~i2_helium_md.simulation.ion.drag_gate_steepness`.

    Returns
    -------
    Callable[[speed, depth], np.ndarray]
        A ``gamma_fn(speed, depth) -> gamma [amu/ps]`` matching the
        :func:`~i2_helium_md.physics.baoab.make_ion_baoab_step` contract.
    """
    if cfg.relaxation_dissipation == "zero_gamma":
        return _zero_gamma
    if cfg.relaxation_dissipation == "landau_gated_drag":
        base_gamma = partial(
            drag_gamma, coeffs=cfg.drag_coefficients, steepness=gate_steepness,
        )
        v_limit = cfg.v_limit_angstrom_per_ps

        def _landau_gated_gamma(speed: np.ndarray, depth: np.ndarray) -> np.ndarray:
            speed = np.asarray(speed, dtype=float)
            # base_gamma is regular everywhere (evaluated for all samples); the
            # gate zeroes the sub-Landau side so decay = 1 there (drag off).
            return np.where(speed > v_limit, base_gamma(speed, depth), 0.0)

        return _landau_gated_gamma
    # check_relaxation_config rejects unknown arms at config-load; defensive.
    raise ValueError(
        f"unknown relaxation_dissipation {cfg.relaxation_dissipation!r}"
    )


def _zero_accel_fn(num_molecules: int):
    """Conservative-force-free acceleration closure for the free-flight arm.

    Matches the ``AccelFn`` contract ``pos -> ((ax, ay, az), E_pot_per_pair(N,))``
    with zero acceleration and zero per-pair potential; the BAOAB step then
    advances positions ballistically (``x -> x + v*dt``) and leaves velocity
    unchanged. E_pot is overwritten by the *held* value in
    :func:`_free_flight_translate`.
    """

    def acc(pos):
        x, y, z = pos
        zero = np.zeros_like(x)
        return (zero, np.zeros_like(y), np.zeros_like(z)), np.zeros(num_molecules)

    return acc


def _freeze_mask(state: IonStepState, *, picture: str, kappa: float,
                 ladder=None) -> np.ndarray:
    """Per-ion frozen mask: ``n == 0`` or ``E_int < D_0(n)`` (the k==0 floor).

    ``D_0(0)`` is undefined on the ladder, so ``n`` is clamped to >= 1 for the
    rung lookup; the ``n == 0`` clause dominates there regardless.
    """
    n = np.asarray(state.n_shell, dtype=float)
    E_int = np.asarray(state.E_int_eV, dtype=float)
    d0 = np.asarray(
        d0_of_n(np.maximum(n, 1.0), picture=picture, kappa=kappa, ladder=ladder)
    )
    return (n <= 0.0) | (E_int < d0)


def _coulomb_translate(
    state: IonStepState,
    *,
    relax_cfg: SimConfig,
    droplet_radii: np.ndarray,
    charge: np.ndarray,
    picture: str,
    kappa: float,
    gamma_fn,
    ladder=None,
) -> IonStepState:
    """One conservative BAOAB translation step + the e_bind_pair fold.

    Byte-for-byte the ion driver's biphasic translation. ``gamma_fn`` selects the
    dissipation arm (:func:`_make_relaxation_gamma_fn`): the default
    :func:`_zero_gamma` (decay = 1, dE_dissip = 0) or the Landau-gated pure-cubic
    drag. Either way the BAOAB step threads its ``dE_dissip`` into ``E_dissip`` on
    top of the K2 cooling drain, so the 5-term invariant closes.
    """
    acc_fn = make_ion_accel_fn(relax_cfg, state.mass_kg, droplet_radii, charge)
    step = make_ion_baoab_step(
        state.mass_kg / U, droplet_radii, acc_fn, gamma_fn, T_eff=0.0,
    )
    new_state = baoab_propagation_step(
        state, step=step, cfg=relax_cfg, droplet_radii=droplet_radii,
    )
    return replace(
        new_state,
        E_pot_eV=new_state.E_pot_eV
        + e_bind_pair_eV(new_state.n_shell, picture=picture, kappa=kappa,
                         ladder=ladder),
    )


def _free_flight_translate(
    state: IonStepState,
    *,
    relax_cfg: SimConfig,
    droplet_radii: np.ndarray,
    held_md_pot: np.ndarray,
    num_molecules: int,
    picture: str,
    kappa: float,
    ladder=None,
) -> IonStepState:
    """One ballistic (force-free) translation step; MD potential held, e_bind tracks n.

    Reuses the BAOAB machinery with a zero acceleration + zero gamma (so positions
    advance ballistically and E_kin is recomputed from the post-jump (m, v) by the
    tested path), then **replaces** E_pot with ``held_md_pot + e_bind_pair(n)`` --
    no conservative work is done under free flight, so the MD potential is
    constant and only the discrete-binding fold moves (by ``+D_0(n)`` per shed).
    """
    step = make_ion_baoab_step(
        state.mass_kg / U, droplet_radii, _zero_accel_fn(num_molecules), _zero_gamma,
        T_eff=0.0,
    )
    new_state = baoab_propagation_step(
        state, step=step, cfg=relax_cfg, droplet_radii=droplet_radii,
    )
    return replace(
        new_state,
        E_pot_eV=held_md_pot
        + e_bind_pair_eV(new_state.n_shell, picture=picture, kappa=kappa,
                         ladder=ladder),
    )


def _build_relaxation_checkpoint(
    ion: IonCheckpoint,
    num_molecules: int,
    stored: list[IonStepState],
) -> IonCheckpoint:
    """Assemble a v7 ``IonCheckpoint`` from the stored relaxation states.

    Trajectory arrays are written column-by-column with the delivered
    :func:`write_ion_state_to_checkpoint_column` (no new I/O). Static /
    pass-through fields: ``droplet_radii_angstrom`` and ``b_ion_outside`` from the
    input ion checkpoint; ``number_of_collisions`` zero; ``mass_scenario`` biphasic.
    """
    two_n = 2 * num_molecules
    T = len(stored)

    def _z2nt() -> np.ndarray:
        return np.zeros((two_n, T), dtype=float)

    ckpt = IonCheckpoint(
        num_molecules=num_molecules,
        time_ps=np.zeros(T, dtype=float),
        positions_x=_z2nt(), positions_y=_z2nt(), positions_z=_z2nt(),
        velocities_x=_z2nt(), velocities_y=_z2nt(), velocities_z=_z2nt(),
        positions_final_x=np.zeros(two_n), positions_final_y=np.zeros(two_n),
        positions_final_z=np.zeros(two_n),
        velocities_final_x=np.zeros(two_n), velocities_final_y=np.zeros(two_n),
        velocities_final_z=np.zeros(two_n),
        mass_kg=stored[0].mass_kg.copy(),           # initial (seed) mass
        mass_final_kg=np.zeros(two_n),
        mass_history_kg=_z2nt(),
        droplet_radii_angstrom=ion.droplet_radii_angstrom.copy(),   # pass-through
        E_kin_eV=_z2nt(), E_pot_eV=_z2nt(), E_dissip_eV=_z2nt(),
        E_mass_transfer_eV=_z2nt(), E_int_eV=_z2nt(), n_shell=_z2nt(),
        b_ion_outside=np.asarray(ion.b_ion_outside).copy(),         # pass-through
        relative_loss_per_ps=_z2nt(),               # not booked on the relaxation path
        number_of_collisions=np.zeros((two_n, T), dtype=int),        # no collisions
        temperature_diagnostic=np.full((T, 3), np.nan, dtype=float),
        mass_scenario="biphasic",
    )

    for t_id, st in enumerate(stored):
        write_ion_state_to_checkpoint_column(st, ckpt, t_id, mass_scenario="biphasic")

    last = stored[-1]
    ckpt.positions_final_x[:] = last.x
    ckpt.positions_final_y[:] = last.y
    ckpt.positions_final_z[:] = last.z
    ckpt.velocities_final_x[:] = last.vx
    ckpt.velocities_final_y[:] = last.vy
    ckpt.velocities_final_z[:] = last.vz
    ckpt.mass_final_kg[:] = last.mass_kg
    return ckpt


def run_relaxation_stage(
    ion: IonCheckpoint,
    cfg: SimConfig,
    *,
    rng: Optional[np.random.Generator] = None,
    save_path: Optional[str | Path] = None,
) -> RelaxationResult:
    """Propagate the biphasic mass subsystem while K2 cooling is live (Slice E2).

    Ungated: runs to freeze-out -- the converged terminal solver. Gated: a
    short handover bridge to a P3-clean state for the detection stage (see
    the module docstring's re-framed contract, Slice DS).

    Parameters
    ----------
    ion : IonCheckpoint
        A finished ``biphasic`` ion-stage checkpoint (v7). The relaxation seeds
        from its **final** stored column.
    cfg : SimConfig
        The run config. Must have ``relaxation_stage_enabled=True`` and pass
        :func:`~i2_helium_md.config.check_relaxation_config` (biphasic scenario,
        ``relaxation_time_ps`` set and > 0, the ``nu*dt <= 0.1`` guard, the forces
        enum). Internally viewed as ``replace(cfg, pickup_rate_coefficient=0.0,
        dt_ion=dt_relax)`` -- not re-validated.
    rng : np.random.Generator, optional
        The stage's generator. Default ``None`` -> a fresh stream seeded
        ``SeedSequence((cfg.seed, RELAXATION_STREAM_KEY))`` (or non-reproducible
        when ``cfg.seed is None``). The ion-stage stream is never touched.
    save_path : str or Path, optional
        If given, the relaxation-window checkpoint is written there via
        ``save_ion_checkpoint`` (Phase F passes ``<run_dir>/relaxation.npz``).

    Returns
    -------
    RelaxationResult
        The relaxation-window checkpoint, per-ion freeze flags, the time reached,
        and the matched-time terminal ``n`` (the E1 relaxed-input hook).

    Raises
    ------
    ValueError
        If ``cfg.relaxation_stage_enabled`` is False, any
        ``check_relaxation_config`` arm fails, the seed checkpoint is not a
        ``biphasic``-scenario checkpoint, or the seed's final stored column is
        not the ion stage's true final state (stride truncation; see below).
    """
    if not cfg.relaxation_stage_enabled:
        raise ValueError(
            "run_relaxation_stage requires cfg.relaxation_stage_enabled=True "
            "(the stage is opt-in; enable it and set relaxation_time_ps)."
        )
    check_relaxation_config(cfg)   # fires the enabled-path checks (fail-loud)

    # Seed-coherence guards (review 2026-07-03; hoisted to the shared
    # checkpoint helper at the Slice-DS review, 2026-07-07 -- one source for
    # both continuation stages): biphasic scenario + the mass stride guard.
    check_biphasic_seed_checkpoint(ion, stage="run_relaxation_stage")

    dt_relax = cfg.dt_ion if cfg.relaxation_dt_ps is None else cfg.relaxation_dt_ps
    # Relaxation view: pickup inert (lambda_0=0), relaxation dt. NOT re-validated
    # (the lambda_0=0 advisory is the sanctioned evaporation-only limit).
    relax_cfg = replace(cfg, pickup_rate_coefficient=0.0, dt_ion=dt_relax)

    if rng is None:
        rng = stage_stream_rng(cfg.seed, RELAXATION_STREAM_KEY)

    num_molecules = int(ion.num_molecules)
    two_n = 2 * num_molecules
    droplet_radii = np.asarray(ion.droplet_radii_angstrom, dtype=float)
    charge = np.ones(two_n, dtype=float)
    picture = cfg.ladder_electronic_picture
    kappa = cfg.ladder_steepness
    # Slice T2 (§I.10): resolve the ladder once for the stage-local consumers
    # (freeze mask, e_bind folds); the pre-T2 stage had NO ladder guard and
    # would have silently run Form-U under a tabulated cfg. biphasic_step
    # resolves its own injection from the same (carried) cfg fields.
    ladder = resolve_ladder(cfg.dissociation_ladder, cfg.tabulated_ladder_rungs_eV)
    forces = cfg.relaxation_forces
    gate_steepness = drag_gate_steepness(cfg)
    # Dissipation arm (§I.11.2 item 2, arm (c)): zero_gamma (default) or the
    # Landau-gated pure-cubic drag. Built once from relax_cfg (carries the drag
    # bundle + cutoff); only the coulomb translation consumes it (the config-load
    # guard forbids landau_gated_drag under free_flight).
    gamma_fn = _make_relaxation_gamma_fn(relax_cfg, gate_steepness)

    seed = ion_state_from_checkpoint_column(ion, -1)

    # Free-flight holds the MD conservative potential constant (no work under free
    # flight); recover it once by removing the seed's e_bind fold.
    held_md_pot = None
    if forces == "free_flight":
        held_md_pot = seed.E_pot_eV - e_bind_pair_eV(
            seed.n_shell, picture=picture, kappa=kappa, ladder=ladder,
        )

    num_internal_steps = math.ceil(cfg.relaxation_time_ps / dt_relax)
    stride, _ = _decide_stride_ion(
        num_molecules, num_internal_steps + 1, DEFAULT_MAX_CHECKPOINT_BYTES_ION,
    )

    stored: list[IonStepState] = [seed]
    state = seed
    freeze_flags = _freeze_mask(seed, picture=picture, kappa=kappa, ladder=ladder)

    for internal_id in range(1, num_internal_steps + 1):
        # Mass subsystem (verbatim biphasic_step under the lambda_0=0 view): draws
        # the frozen two-draw stream (evaporation then the inert pickup) every step.
        state = biphasic_step(
            state, rng=rng, cfg=relax_cfg, droplet_radii=droplet_radii,
            gate_steepness=gate_steepness,
        )
        # Translation (drag off).
        if forces == "coulomb":
            state = _coulomb_translate(
                state, relax_cfg=relax_cfg, droplet_radii=droplet_radii,
                charge=charge, picture=picture, kappa=kappa, gamma_fn=gamma_fn,
                ladder=ladder,
            )
        else:
            state = _free_flight_translate(
                state, relax_cfg=relax_cfg, droplet_radii=droplet_radii,
                held_md_pot=held_md_pot, num_molecules=num_molecules,
                picture=picture, kappa=kappa, ladder=ladder,
            )

        if internal_id % stride == 0:
            stored.append(state)

        freeze_flags = _freeze_mask(state, picture=picture, kappa=kappa,
                                    ladder=ladder)
        if bool(freeze_flags.all()):
            break

    if stored[-1] is not state:
        stored.append(state)

    ckpt = _build_relaxation_checkpoint(ion, num_molecules, stored)

    if save_path is not None:
        save_ion_checkpoint(ckpt, save_path)

    return RelaxationResult(
        checkpoint=ckpt,
        freeze_flags=freeze_flags,
        # Window-relative (the seed enters at the ion stage's absolute end time,
        # ~20 ps in production; the checkpoint's time_ps axis stays absolute).
        time_relaxed_ps=float(state.time_ps) - float(seed.time_ps),
        terminal_n=np.asarray(state.n_shell, dtype=float).copy(),
    )
