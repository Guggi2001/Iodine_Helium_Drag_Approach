"""Post-ejection relaxation stage (Tier-2 Phase E, Slice E2 -- the R5 mitigation).

Propagates the biphasic **mass subsystem** past the 20 ps ion stage to the
experimental timescale, so the terminal I+He_n size distribution (E1) is read at
*matched time* rather than at the truncated sim-end upper bound (MASS doc §R5:
the energy-gated cascade continues for hundreds of ps as loosely-bound outer He
keep evaporating).

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
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional

import numpy as np

from ..config import SimConfig, check_relaxation_config
from ..physics.constants import U
from ..physics.baoab import make_ion_baoab_step
from ..physics.dissociation_ladder import d0_of_n
from ..physics.leapfrog import make_ion_accel_fn
from ..physics.solvation_cooling import e_bind_pair_eV
from .checkpoint import IonCheckpoint, save_ion_checkpoint
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
        The time actually reached (``relaxation_time_ps`` or the earlier
        all-frozen time).
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


def _freeze_mask(state: IonStepState, *, picture: str, kappa: float) -> np.ndarray:
    """Per-ion frozen mask: ``n == 0`` or ``E_int < D_0(n)`` (the k==0 floor).

    ``D_0(0)`` is undefined on the ladder, so ``n`` is clamped to >= 1 for the
    rung lookup; the ``n == 0`` clause dominates there regardless.
    """
    n = np.asarray(state.n_shell, dtype=float)
    E_int = np.asarray(state.E_int_eV, dtype=float)
    d0 = np.asarray(d0_of_n(np.maximum(n, 1.0), picture=picture, kappa=kappa))
    return (n <= 0.0) | (E_int < d0)


def _coulomb_translate(
    state: IonStepState,
    *,
    relax_cfg: SimConfig,
    droplet_radii: np.ndarray,
    charge: np.ndarray,
    picture: str,
    kappa: float,
) -> IonStepState:
    """One conservative (zero-gamma BAOAB) translation step + the e_bind_pair fold.

    Byte-for-byte the ion driver's biphasic translation, with the drag ``gamma_fn``
    replaced by :func:`_zero_gamma` (decay = 1, dE_dissip = 0).
    """
    acc_fn = make_ion_accel_fn(relax_cfg, state.mass_kg, droplet_radii, charge)
    step = make_ion_baoab_step(
        state.mass_kg / U, droplet_radii, acc_fn, _zero_gamma, T_eff=0.0,
    )
    new_state = baoab_propagation_step(
        state, step=step, cfg=relax_cfg, droplet_radii=droplet_radii,
    )
    return replace(
        new_state,
        E_pot_eV=new_state.E_pot_eV
        + e_bind_pair_eV(new_state.n_shell, picture=picture, kappa=kappa),
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
        + e_bind_pair_eV(new_state.n_shell, picture=picture, kappa=kappa),
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
    """Propagate the biphasic mass subsystem to the experimental timescale (Slice E2).

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
        If ``cfg.relaxation_stage_enabled`` is False, or any
        ``check_relaxation_config`` arm fails.
    """
    if not cfg.relaxation_stage_enabled:
        raise ValueError(
            "run_relaxation_stage requires cfg.relaxation_stage_enabled=True "
            "(the stage is opt-in; enable it and set relaxation_time_ps)."
        )
    check_relaxation_config(cfg)   # fires the enabled-path checks (fail-loud)

    dt_relax = cfg.dt_ion if cfg.relaxation_dt_ps is None else cfg.relaxation_dt_ps
    # Relaxation view: pickup inert (lambda_0=0), relaxation dt. NOT re-validated
    # (the lambda_0=0 advisory is the sanctioned evaporation-only limit).
    relax_cfg = replace(cfg, pickup_rate_coefficient=0.0, dt_ion=dt_relax)

    if rng is None:
        if cfg.seed is None:
            rng = np.random.default_rng()
        else:
            rng = np.random.default_rng(
                np.random.SeedSequence((int(cfg.seed), RELAXATION_STREAM_KEY))
            )

    num_molecules = int(ion.num_molecules)
    two_n = 2 * num_molecules
    droplet_radii = np.asarray(ion.droplet_radii_angstrom, dtype=float)
    charge = np.ones(two_n, dtype=float)
    picture = cfg.ladder_electronic_picture
    kappa = cfg.ladder_steepness
    forces = cfg.relaxation_forces
    gate_steepness = drag_gate_steepness(cfg)

    seed = ion_state_from_checkpoint_column(ion, -1)

    # Free-flight holds the MD conservative potential constant (no work under free
    # flight); recover it once by removing the seed's e_bind fold.
    held_md_pot = None
    if forces == "free_flight":
        held_md_pot = seed.E_pot_eV - e_bind_pair_eV(
            seed.n_shell, picture=picture, kappa=kappa,
        )

    num_internal_steps = math.ceil(cfg.relaxation_time_ps / dt_relax)
    stride, _ = _decide_stride_ion(
        num_molecules, num_internal_steps + 1, DEFAULT_MAX_CHECKPOINT_BYTES_ION,
    )

    stored: list[IonStepState] = [seed]
    state = seed
    freeze_flags = _freeze_mask(seed, picture=picture, kappa=kappa)

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
                charge=charge, picture=picture, kappa=kappa,
            )
        else:
            state = _free_flight_translate(
                state, relax_cfg=relax_cfg, droplet_radii=droplet_radii,
                held_md_pot=held_md_pot, num_molecules=num_molecules,
                picture=picture, kappa=kappa,
            )

        if internal_id % stride == 0:
            stored.append(state)

        freeze_flags = _freeze_mask(state, picture=picture, kappa=kappa)
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
        time_relaxed_ps=float(state.time_ps),
        terminal_n=np.asarray(state.n_shell, dtype=float).copy(),
    )
