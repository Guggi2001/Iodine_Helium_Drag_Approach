"""Top-level driver for ion propagation.

Public API:

* :func:`run_ion_propagation` -- run a full ion-stage simulation from a
  finished neutral checkpoint and return an :class:`IonCheckpoint`.
  Optionally save it via a :class:`RunDirectory`.

The driver mirrors :func:`run_neutral_propagation` (``simulation/neutral.py``)
in structure:

1. Refuses out-of-scope configs up front (``cfg.single_pulse=False`` uses
   MATLAB dt-switching which is pump-probe scope; not implemented).
2. Decides how many internal timesteps to run from
   ``cfg.ion_simulation_time / cfg.dt_ion``.
3. Decides the storage stride: if a full-resolution checkpoint would
   exceed ``max_bytes`` (default 1 GB), only every K-th internal step
   is stored. Internal steps still happen at ``cfg.dt_ion`` because the
   Mode-3 collision sampler requires the previous-internal-step
   displacement.
4. Builds the t=0 state with :func:`build_initial_ion_state` sized for
   the *stored* number of steps.
5. Runs the inner loop using :func:`ion_propagation_step`, storing
   every K-th step into the checkpoint (always including the final
   state in the last reachable column so the trajectory ends at the
   actual end-time rather than a pre-stride snapshot).
6. Populates the post-loop "final state" fields
   (``positions_final_*``, ``velocities_final_*``, ``mass_final_kg``,
   ``b_ion_outside``).
7. Returns the checkpoint and, if a :class:`RunDirectory` was provided,
   saves it (auto-saving cfg if not already present).

Reference MATLAB implementation: ``vmi_sim_3d_ion_propa.m``. Two known
MATLAB t=0 bookkeeping bugs (E_kin missing v_z and squared, E_pot
using 2D radius and missing partner Coulomb) are fixed by
``build_initial_ion_state`` and must NOT be reintroduced here.
"""

from __future__ import annotations

import math
from dataclasses import replace
from functools import partial
from typing import Optional

import numpy as np

from ..config import SimConfig
from ..physics.baoab import make_ion_baoab_step
from ..physics.constants import U
from ..physics.drag import drag_gamma
from ..physics.leapfrog import make_ion_accel_fn
from ..physics.shell_schedule import build_onset_strip_schedule, build_shell_schedule
from ..physics.state_coupling import (
    apply_state_factor,
    derive_n_ref_amu,
    shell_area_state_factor,
)
from ..physics.dissociation_ladder import resolve_ladder
from ..physics.solvation_cooling import e_bind_pair_eV
from ..sampling.ce_channels import ce_pair_scale_from_checkpoint
from .checkpoint import IonCheckpoint, NeutralCheckpoint, stage_stream_rng
from .ion_initial_state import build_initial_ion_state
from .ion_propagation_step import (
    IonStepState,
    baoab_propagation_step,
    biphasic_step,
    exit_strip_step,
    ion_propagation_step,
    ion_state_from_checkpoint_column,
    shed_step,
    write_ion_state_to_checkpoint_column,
    _check_drag_scope,
    _depth,
)
from .run_directory import RunDirectory

# Fixed module-level key for the ion driver's exit-strip Bernoulli stream via
# SeedSequence((cfg.seed, EXIT_STRIP_ION_STREAM_KEY)) -- dedicated, appended
# after all existing streams (the CE design §5 draw-order contract: the
# ion-stage stream is never touched; strip-off consumes nothing). The
# relaxation stage derives its own strip stream from its own key.
EXIT_STRIP_ION_STREAM_KEY: int = 0xCE2_2026


# ===========================================================================
# Memory-budget settings
# ===========================================================================
#: Maximum approximate ion checkpoint size before auto-stride kicks in.
#: Same default as the neutral driver (1 GB).
DEFAULT_MAX_CHECKPOINT_BYTES_ION: int = 1000000000

#: Number of (2N, num_steps) trajectory/diagnostic arrays in an IonCheckpoint:
#: 6 positions/velocities + E_kin + E_pot + E_dissip + E_mass_transfer
#: + E_int_eV (schema v7) + n_shell + relative_loss_per_ps
#: + number_of_collisions + mass_history_kg = 15.
#: All are 8 bytes per cell (number_of_collisions is int64). Must track the
#: checkpoint schema (guarded by test_ion.py's schema-count test): an
#: undercount silently grants stride 1 past the byte budget.
_NUM_2N_T_ARRAYS_ION: int = 15


# ===========================================================================
# Public API
# ===========================================================================
def run_ion_propagation(
    cfg: SimConfig,
    neutral_ckpt: NeutralCheckpoint,
    *,
    rng: np.random.Generator | None = None,
    run_dir: Optional[RunDirectory] = None,
    max_bytes: int = DEFAULT_MAX_CHECKPOINT_BYTES_ION,
    verbose: bool = False,
) -> IonCheckpoint:
    """Run ion propagation and return the resulting checkpoint.

    Parameters
    ----------
    cfg : SimConfig
        Simulation config. Reads ``num_molecules``, ``dt_ion``,
        ``ion_simulation_time``, ``single_pulse``, plus everything used
        by ``build_initial_ion_state`` and ``ion_propagation_step``.
    neutral_ckpt : NeutralCheckpoint
        End-state of the neutral stage. The ion stage starts from its
        last column (matches the production single-pulse use case).
    rng : np.random.Generator, optional
        Reproducible RNG. If None, built from ``cfg.seed`` (same
        convention as the neutral driver).
    run_dir : RunDirectory, optional
        If provided, the resulting checkpoint is saved via
        ``run_dir.save_ion(...)``. The cfg is also saved if not
        already present.
    max_bytes : int, optional
        Approximate cap on the size of the stored checkpoint. If a
        full-resolution checkpoint would exceed this, the stride is
        increased so only every K-th internal step is stored. Default 1 GB.
    verbose : bool, optional
        If True, print progress messages (number of steps, stride,
        memory estimate).

    Returns
    -------
    IonCheckpoint
        With shape ``(2N, num_stored_steps)`` for trajectory arrays and
        all final-state fields populated from the last internal step.

    Raises
    ------
    NotImplementedError
        If ``cfg.single_pulse`` is False (the MATLAB dt-switching
        non-single-pulse branch is pump-probe scope and out of scope
        for this stage of the project).
    """
    if rng is None:
        rng = np.random.default_rng(cfg.seed)

    # 1. Driver-level scope check.
    _check_scope_ion_driver(cfg)

    # 2. Decide internal step count.
    num_internal_steps = _internal_step_count_ion(cfg)

    # 3. Decide stride.
    stride, num_stored_steps = _decide_stride_ion(
        cfg.num_molecules, num_internal_steps, max_bytes,
    )
    if verbose:
        est_mb = _estimate_checkpoint_bytes_ion(
            cfg.num_molecules, num_stored_steps,
        ) / 1e6
        print(
            f"ion propagation: "
            f"{num_internal_steps} internal steps, "
            f"stride={stride}, {num_stored_steps} stored steps, "
            f"~{est_mb:.1f} MB checkpoint"
        )

    # 4. Build initial state (allocates the storage-sized arrays).
    ckpt = build_initial_ion_state(
        cfg, neutral_ckpt,
        num_steps_ion=num_stored_steps,
        start_id=-1,
        rng=rng,
    )

    # 5. Inner loop.
    if num_internal_steps <= 1:
        # No propagation at all; only column 0 is meaningful.
        # Still populate the final-state fields from column 0 so callers
        # see a self-consistent checkpoint.
        state = ion_state_from_checkpoint_column(ckpt, 0)
        _write_final_state(state, ckpt, cfg)
        if run_dir is not None:
            _save_with_cfg_ion(ckpt, cfg, run_dir)
        return ckpt

    # Allocate constants used inside the loop once.
    charge = np.ones(2 * cfg.num_molecules, dtype=float)
    droplet_radii = ckpt.droplet_radii_angstrom
    dt = cfg.dt_ion   # ps; the internal step (the schedule shed-window is (t, t+dt])

    # Dispatch once per run (D3): a validated drag-coefficient bundle routes the
    # ion stage onto the BAOAB drag path; otherwise the hard-sphere collision
    # path runs unchanged. The predicate is both *necessary* (BAOAB cannot run
    # without coefficients) and *already validated* (a non-None bundle passed the
    # Slice 3 check_drag_config: consistency + dissipativity + form agreement).
    use_drag = cfg.drag_coefficients is not None
    gamma_fn = None
    schedule = None
    next_shed_idx = 0
    ladder = None
    # Tier-2 (C) surfaces (both None/off = byte-identical): the per-molecule
    # CE Coulomb scale from the t0 channel draw (stored in the v8 fields by
    # build_initial_ion_state), and the depth-graded exit strip.
    pair_scale = ce_pair_scale_from_checkpoint(ckpt)
    strip_live = cfg.exit_strip_mode == "depth_graded"
    strip_rng = None
    strip_counts = None
    prev_depth = None
    if strip_live:
        strip_rng = stage_stream_rng(cfg.seed, EXIT_STRIP_ION_STREAM_KEY)
        strip_counts = np.zeros(2 * cfg.num_molecules, dtype=int)
    if use_drag:
        if cfg.mass_scenario == "biphasic":
            # Slice T2 (§I.10): resolve the ladder once for the per-step
            # e_bind_pair E_pot fold below (biphasic_step resolves its own
            # injection from the same cfg fields; fail-loud on a tabulated
            # selector without a table).
            ladder = resolve_ladder(
                cfg.dissociation_ladder, cfg.tabulated_ladder_rungs_eV
            )
        # Pass the *realized* initial ion mass (ckpt.mass_kg, downstream of the
        # build_initial_ion_state m_eff override) so the scope guard's mass
        # trip-wire checks what the stepper will actually integrate, not a
        # config field.
        _check_drag_scope(cfg, ckpt.mass_kg)
        gate_steepness = drag_gate_steepness(cfg)
        # gamma_fn(speed, depth) -> gamma [amu/ps]; the erf gate (§5.5) lives
        # inside drag_gamma via the steepness arg (hard FDT coupling, §5.2).
        gamma_fn = partial(
            drag_gamma, coeffs=cfg.drag_coefficients, steepness=gate_steepness,
        )
        # Tier-2 atlas §3.5i s(n) state coupling (design S2, BC-1/BC-2): under
        # "shell_area" (biphasic-only, config-load-guarded) each step's closure
        # rebuild wraps gamma_fn with the per-ion factor s(n_shell) evaluated
        # on the post-event shell state -- jump-then-O, the same state the
        # rebuild reads for m(t). Under "off" the base gamma_fn is passed
        # verbatim below (no wrapper object; structurally bit-identical).
        state_coupling_live = cfg.drag_state_coupling == "shell_area"
        n_ref_coupling = None
        if state_coupling_live:
            n_ref_coupling = derive_n_ref_amu(
                cfg.drag_coefficients.extraction_mass_amu
            )
        # Tier-1a: the anchored He-shell schedule drives the variable mass m(t).
        # Built once; queried per step (shed_step). Absent under `fixed`, so the
        # fixed-mass path below is byte-for-byte the Tier-0 path (regression guard).
        if cfg.mass_scenario == "anchored_discrete":
            if cfg.anchor_mode == "time":
                schedule = build_shell_schedule(cfg.t_star_ps)
            elif cfg.anchor_mode == "onset_strip":
                schedule = build_onset_strip_schedule(
                    t_strip_ps=cfg.t_star_ps,
                    n_final=cfg.anchor_n_final,
                )
            else:
                raise NotImplementedError(f"unsupported anchor_mode={cfg.anchor_mode!r}")

    state = ion_state_from_checkpoint_column(ckpt, 0)
    prev_dist: np.ndarray | None = None
    next_storage_idx = 1
    if strip_live:
        prev_depth = _depth(state.x, state.y, state.z, droplet_radii)

    for internal_id in range(1, num_internal_steps):
        if use_drag:
            if cfg.mass_scenario == "biphasic":
                # Tier-2 generative pre-step (Slice G): K2 cooling + <=1 mass event
                # (evaporation-then-pickup) + the E_int/E_mass_transfer/E_dissip
                # bookings, at the step seam BEFORE the closure rebuild -- so the
                # rebuild reads the post-jump mass m+ for both the conservative kicks
                # and the O-step (jump-then-O). Threads rng (the channel draws).
                state = biphasic_step(
                    state, rng=rng, cfg=cfg, droplet_radii=droplet_radii,
                    gate_steepness=gate_steepness,
                )
            elif schedule is not None:
                # SQ2/SQ3 (anchored_discrete only): apply at most one scheduled cold
                # shed BEFORE rebuilding the closure, so the rebuild reads the
                # post-shed mass m+ for both the conservative kicks and the O-step.
                state, next_shed_idx = shed_step(
                    state, schedule, next_shed_idx, dt,
                )
            # Rebuild the BAOAB closure every step, matching the make_ion_step
            # rebuild pattern. Mass enters here in amu (kg -> amu via U); under
            # `fixed` it is constant, under `anchored_discrete` it follows the
            # shed schedule, under `biphasic` the generative channels. Noise dormant
            # (T_eff=0).
            acc_fn = make_ion_accel_fn(
                cfg, state.mass_kg, droplet_radii, charge,
                pair_scale=pair_scale,
            )
            step_gamma_fn = gamma_fn
            if state_coupling_live:
                step_gamma_fn = apply_state_factor(
                    gamma_fn,
                    shell_area_state_factor(
                        state.n_shell,
                        R_core_angstrom=cfg.state_coupling_R_core_angstrom,
                        rho_shell_per_A3=cfg.state_coupling_rho_shell_per_A3,
                        n_ref=n_ref_coupling,
                    ),
                )
            step = make_ion_baoab_step(
                state.mass_kg / U, droplet_radii, acc_fn, step_gamma_fn,
                T_eff=0.0,
            )
            new_state = baoab_propagation_step(
                state, step=step, cfg=cfg, droplet_radii=droplet_radii,
            )
            if cfg.mass_scenario == "biphasic":
                # Fold the pair-binding potential E_bind^pair(n) = -Sigma(n) into the
                # freshly-recomputed MD E_pot (a pure function of the post-event n), so
                # the +/- D_0 pickup/shed bookings close the 5-term invariant (MASS §6;
                # the t0 seed lives in build_initial_ion_state).
                new_state = replace(
                    new_state,
                    E_pot_eV=new_state.E_pot_eV + e_bind_pair_eV(
                        new_state.n_shell,
                        picture=cfg.ladder_electronic_picture,
                        kappa=cfg.ladder_steepness,
                        ladder=ladder,
                    ),
                )
                # Tier-2 (C) exit strip at this step's outbound crossings
                # (design §3.3; AFTER the fold -- the operator books its own
                # fold delta for the knocked rungs). Dedicated stream; off =
                # structurally absent.
                if strip_live:
                    new_state, prev_depth, knocks = exit_strip_step(
                        new_state, rng=strip_rng, cfg=cfg,
                        droplet_radii=droplet_radii,
                        prev_depth_angstrom=prev_depth, ladder=ladder,
                    )
                    strip_counts += knocks
        else:
            new_state = ion_propagation_step(
                state,
                cfg=cfg,
                droplet_radii=droplet_radii,
                charge=charge,
                prev_distance_angstrom=prev_dist,
                rng=rng,
            )
            prev_dist = _state_step_distance_ion(state, new_state)
        state = new_state

        # Store every stride-th internal step.
        if internal_id % stride == 0 and next_storage_idx < num_stored_steps:
            write_ion_state_to_checkpoint_column(
                state, ckpt, next_storage_idx, mass_scenario=cfg.mass_scenario,
            )
            if state.temperature_diagnostic is not None:
                ckpt.temperature_diagnostic[next_storage_idx, :] = (
                    state.temperature_diagnostic
                )
            next_storage_idx += 1

    # If we ended up with fewer stored steps than allocated (rare; happens
    # when num_internal_steps - 1 is not divisible by stride), make sure the
    # last reachable column holds the final state. This keeps the trajectory
    # ending at the actual end-time rather than at a pre-stride snapshot.
    if next_storage_idx < num_stored_steps:
        write_ion_state_to_checkpoint_column(
            state, ckpt, next_storage_idx, mass_scenario=cfg.mass_scenario,
        )
        if state.temperature_diagnostic is not None:
            ckpt.temperature_diagnostic[next_storage_idx, :] = (
                state.temperature_diagnostic
            )

    # 6. Final-state fields, taken from the actual last internal step.
    _write_final_state(state, ckpt, cfg)
    if strip_counts is not None:
        # Cumulative v8 strip counter (build_initial_ion_state zero-seeds it).
        ckpt.ce_strip_count[:] = ckpt.ce_strip_count + strip_counts

    # 7. Save if run_dir given.
    if run_dir is not None:
        _save_with_cfg_ion(ckpt, cfg, run_dir)

    return ckpt


# ===========================================================================
# Internal helpers
# ===========================================================================
def _check_scope_ion_driver(cfg: SimConfig) -> None:
    """Refuse to run with cfg flags that need driver-level support
    we haven't implemented.

    The other ion-stage scope flags (``effusive_dynamics``,
    ``single_charge_ionization_allowed``, ``additional_droplet_charges``,
    ``highly_charged_iodine``, ``hard_sphere_collision_mode != 3``) are
    not re-checked here -- ``build_initial_ion_state`` and
    ``ion_propagation_step`` already raise on those, and the failure
    surfaces before any expensive stepping.
    """
    if cfg.mass_scenario == "biphasic" and cfg.drag_coefficients is None:
        # Slice-G review fix: biphasic is a drag-path scenario. Without a bundle
        # the loop below would dispatch onto the hard-sphere collision path, where
        # E_int is never evolved, the E_pot binding fold is lost after column 0,
        # and the writer would store a frozen n_shell while attachment grows the
        # mass. config.check_biphasic_config refuses this at load; this driver
        # check covers configs built without validate().
        raise NotImplementedError(
            "mass_scenario='biphasic' requires a drag_coefficients bundle: the "
            "generative biphasic mechanism runs only on the BAOAB drag path, and "
            "drag_coefficients=None would dispatch onto the hard-sphere collision "
            "path."
        )
    if not cfg.single_pulse:
        raise NotImplementedError(
            "run_ion_propagation requires cfg.single_pulse=True. The "
            "non-single-pulse MATLAB branch uses dt-switching "
            "(dt_fine/dt_coarse + switchtime), which is pump-probe "
            "scope and not implemented in this stage."
        )


def drag_gate_steepness(cfg: SimConfig) -> float:
    """Resolve the spatial-gate steepness for the drag branch (§5.5 collapse).

    The drag gate stays G2: ``density_proportional`` (the default) collapses to
    the erf complement (G2) -- identical to ``erf_tied`` -- and both use
    ``cfg.potential_steepness``. ``erf_independent`` (G3) uses its own
    ``cfg.drag_gate_steepness`` (which itself defaults to ``potential_steepness``).
    (The Tier-2 ``cfg.helium_density_profile`` is the *pickup* occupancy gate, a
    separate G2 surface-density quantity -- **not** a drag-gate G4 promotion; the
    Slice-rho density gate reuses this resolver's steepness so the two share one
    surface.) The discarded sharp boolean gate (G1) has no
    continuous implementation in ``physics/drag.py`` and is rejected here so the
    default preset cannot silently fall through with no gate.
    """
    gate = cfg.drag_spatial_gate
    if gate in ("density_proportional", "erf_tied"):
        return cfg.potential_steepness
    if gate == "erf_independent":
        return cfg.drag_gate_steepness
    raise NotImplementedError(
        f"drag_spatial_gate={gate!r} is not supported on the Tier-0 drag path; "
        "use 'density_proportional' (default), 'erf_tied', or 'erf_independent'. "
        "'sharp' (G1) is discarded -- a discontinuous force breaks the BAOAB "
        "O-step (DRAG_PORT_DESIGN_DECISIONS.md §5.4)."
    )


def _internal_step_count_ion(cfg: SimConfig) -> int:
    """Number of leapfrog steps to run for the ion stage.

    Equals ``ceil(cfg.ion_simulation_time / cfg.dt_ion)``. Matches the
    MATLAB ``ion_timesteps = ceil(ion_simulation_time/dt)`` line in
    ``vmi_sim_3d_ion_propa.m``. For the production single-pulse case
    (20 ps / 0.01 ps) this is 2000.
    """
    return math.ceil(cfg.ion_simulation_time / cfg.dt_ion)


def _estimate_checkpoint_bytes_ion(num_molecules: int, num_steps: int) -> int:
    """Estimate the in-memory size of an IonCheckpoint.

    Counts:

    * ``_NUM_2N_T_ARRAYS_ION`` (15) trajectory/diagnostic arrays of shape
      (2N, num_steps) at 8 bytes/cell
    * Static (2N,) arrays for mass_kg, mass_final_kg, droplet_radii_angstrom,
      the six positions/velocities final placeholders, and the three v8
      (C)-design per-ion fields (12 arrays).
    * Static (N,) array for b_ion_outside (1 byte, but we count 8 for slack).
    * (num_steps,) for time_ps.
    """
    n_atoms = 2 * num_molecules
    bytes_2N_T = _NUM_2N_T_ARRAYS_ION * n_atoms * num_steps * 8
    bytes_static_2N = 12 * n_atoms * 8
    bytes_static_N = num_molecules * 8
    bytes_time = num_steps * 8
    return bytes_2N_T + bytes_static_2N + bytes_static_N + bytes_time


def _decide_stride_ion(
    num_molecules: int,
    num_internal_steps: int,
    max_bytes: int,
) -> tuple[int, int]:
    """Return ``(stride, num_stored_steps)`` for the given budget.

    Storage stride is the smallest integer K such that storing every
    K-th step keeps the checkpoint under ``max_bytes``. Internal steps
    still happen at ``cfg.dt_ion``; only storage is downsampled.

    The stored trajectory always includes step 0 (set by
    ``build_initial_ion_state``), so
    ``num_stored_steps = ceil(num_internal_steps / stride)``.
    """
    full_size = _estimate_checkpoint_bytes_ion(num_molecules, num_internal_steps)
    if full_size <= max_bytes:
        return 1, num_internal_steps

    bytes_per_step = _NUM_2N_T_ARRAYS_ION * 2 * num_molecules * 8
    if bytes_per_step <= 0:
        return 1, num_internal_steps
    max_stored = max(2, int(max_bytes // bytes_per_step))
    stride = max(1, math.ceil(num_internal_steps / max_stored))
    num_stored = math.ceil(num_internal_steps / stride)
    return stride, num_stored


def _state_step_distance_ion(
    prev: IonStepState,
    new: IonStepState,
) -> np.ndarray:
    """Per-atom distance traveled between two consecutive ion states."""
    return np.sqrt(
        (new.x - prev.x) ** 2
        + (new.y - prev.y) ** 2
        + (new.z - prev.z) ** 2
    )


def _write_final_state(
    state: IonStepState,
    ckpt: IonCheckpoint,
    cfg: SimConfig,
) -> None:
    """Populate the post-loop final-state fields of the checkpoint.

    ``b_ion_outside`` (shape (N,)) uses the per-molecule OR rule: the
    molecule counts as having an ion outside the droplet if either of
    its two atoms has ``depth > 0`` at the final state.
    """
    ckpt.positions_final_x[:] = state.x
    ckpt.positions_final_y[:] = state.y
    ckpt.positions_final_z[:] = state.z
    ckpt.velocities_final_x[:] = state.vx
    ckpt.velocities_final_y[:] = state.vy
    ckpt.velocities_final_z[:] = state.vz
    ckpt.mass_final_kg[:] = state.mass_kg

    N = cfg.num_molecules
    depth_final = (
        np.sqrt(state.x ** 2 + state.y ** 2 + state.z ** 2)
        - ckpt.droplet_radii_angstrom
    )
    atom1_outside = depth_final[:N] > 0
    atom2_outside = depth_final[N:] > 0
    ckpt.b_ion_outside[:] = atom1_outside | atom2_outside


def _save_with_cfg_ion(
    ckpt: IonCheckpoint,
    cfg: SimConfig,
    run_dir: RunDirectory,
) -> None:
    """Save the ion checkpoint and the cfg in the run directory."""
    if not run_dir.has_cfg():
        run_dir.save_cfg(cfg)
    run_dir.save_ion(ckpt)
