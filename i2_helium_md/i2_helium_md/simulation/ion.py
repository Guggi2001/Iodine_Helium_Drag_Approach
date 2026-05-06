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
from typing import Optional

import numpy as np

from ..config import SimConfig
from .checkpoint import IonCheckpoint, NeutralCheckpoint
from .ion_initial_state import build_initial_ion_state
from .ion_propagation_step import (
    IonStepState,
    ion_propagation_step,
    ion_state_from_checkpoint_column,
    write_ion_state_to_checkpoint_column,
)
from .run_directory import RunDirectory


# ===========================================================================
# Memory-budget settings
# ===========================================================================
#: Maximum approximate ion checkpoint size before auto-stride kicks in.
#: Same default as the neutral driver (1 GB).
DEFAULT_MAX_CHECKPOINT_BYTES_ION: int = 1000000000

#: Number of (2N, num_steps) trajectory/diagnostic arrays in an IonCheckpoint:
#: 6 positions/velocities + E_kin + E_pot + E_dissip + E_mass_attach_defect
#: + relative_loss_per_ps + number_of_collisions + mass_history_kg = 13.
#: All are 8 bytes per cell (number_of_collisions is int64).
_NUM_2N_T_ARRAYS_ION: int = 13


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

    state = ion_state_from_checkpoint_column(ckpt, 0)
    prev_dist: np.ndarray | None = None
    next_storage_idx = 1

    for internal_id in range(1, num_internal_steps):
        new_state = ion_propagation_step(
            state,
            cfg=cfg,
            droplet_radii=ckpt.droplet_radii_angstrom,
            charge=charge,
            prev_distance_angstrom=prev_dist,
            rng=rng,
        )
        prev_dist = _state_step_distance_ion(state, new_state)
        state = new_state

        # Store every stride-th internal step.
        if internal_id % stride == 0 and next_storage_idx < num_stored_steps:
            write_ion_state_to_checkpoint_column(state, ckpt, next_storage_idx)
            next_storage_idx += 1

    # If we ended up with fewer stored steps than allocated (rare; happens
    # when num_internal_steps - 1 is not divisible by stride), make sure the
    # last reachable column holds the final state. This keeps the trajectory
    # ending at the actual end-time rather than at a pre-stride snapshot.
    if next_storage_idx < num_stored_steps:
        write_ion_state_to_checkpoint_column(state, ckpt, next_storage_idx)

    # 6. Final-state fields, taken from the actual last internal step.
    _write_final_state(state, ckpt, cfg)

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
    if not cfg.single_pulse:
        raise NotImplementedError(
            "run_ion_propagation requires cfg.single_pulse=True. The "
            "non-single-pulse MATLAB branch uses dt-switching "
            "(dt_fine/dt_coarse + switchtime), which is pump-probe "
            "scope and not implemented in this stage."
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

    * 12 trajectory/diagnostic arrays of shape (2N, num_steps) at 8 bytes/cell
    * Static (2N,) arrays for mass_kg, mass_final_kg, droplet_radii_angstrom,
      and the six positions/velocities final placeholders (9 arrays).
    * Static (N,) array for b_ion_outside (1 byte, but we count 8 for slack).
    * (num_steps,) for time_ps.
    """
    n_atoms = 2 * num_molecules
    bytes_2N_T = _NUM_2N_T_ARRAYS_ION * n_atoms * num_steps * 8
    bytes_static_2N = 9 * n_atoms * 8
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
