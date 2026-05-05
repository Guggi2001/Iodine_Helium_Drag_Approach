"""Top-level driver for neutral propagation.

Public API:

* :func:`run_neutral_propagation` -- run a full neutral simulation and
  return a :class:`NeutralCheckpoint`. Optionally save it via a
  :class:`RunDirectory`.

The driver:

1. Decides how many internal timesteps to run (single-pulse mode runs
   only ``2*dt`` total, matching MATLAB).
2. Decides the storage stride: if a full-resolution checkpoint would
   exceed ``MAX_CHECKPOINT_BYTES`` (default 300 MB), only every
   K-th internal step is stored. Internal steps still happen at
   ``cfg.dt_neutral``.
3. Builds the initial state with :func:`build_initial_state` sized
   for the *stored* number of steps.
4. Optionally calls a (currently stubbed) DFT pre-fill if
   ``cfg.custom_DFT_start`` is True.
5. Runs the inner loop using :func:`neutral_propagation_step`,
   storing every K-th step into the checkpoint.
6. Returns the checkpoint and, if a :class:`RunDirectory` was provided,
   saves it.

The auto-stride lets us scale from short single-pulse runs (~1 MB)
all the way to long runs (200 ps × 20,000 atoms internal × stride 22)
without changing the user-facing API.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from ..config import SimConfig
from .checkpoint import NeutralCheckpoint
from .initial_state import build_initial_state
from .propagation_step import (
    NeutralStepState,
    neutral_propagation_step,
    state_from_checkpoint_column,
    write_state_to_checkpoint_column,
)
from .run_directory import RunDirectory


# ===========================================================================
# Memory-budget settings
# ===========================================================================
#: Maximum approximate checkpoint size before auto-stride kicks in.
#: Set to 300 MB by default; can be passed in or modified for special runs.
DEFAULT_MAX_CHECKPOINT_BYTES: int = 300_000_000

#: Number of (2N, num_steps) trajectory arrays in a NeutralCheckpoint.
#: 6 positions+velocities + 4 energy/L_droplet diagnostics.
_NUM_2N_T_ARRAYS: int = 10


# ===========================================================================
# Public API
# ===========================================================================
def run_neutral_propagation(
    cfg: SimConfig,
    *,
    rng: np.random.Generator | None = None,
    run_dir: Optional[RunDirectory] = None,
    max_bytes: int = DEFAULT_MAX_CHECKPOINT_BYTES,
    verbose: bool = False,
) -> NeutralCheckpoint:
    """Run neutral propagation and return the resulting checkpoint.

    Parameters
    ----------
    cfg : SimConfig
        Simulation config. Reads ``num_molecules``, ``dt_neutral``,
        ``t_max_neutral``, ``single_pulse``, ``custom_DFT_start``,
        plus everything used by ``build_initial_state`` and
        ``neutral_propagation_step``.
    rng : np.random.Generator, optional
        Reproducible RNG. If None, built from ``cfg.seed``.
    run_dir : RunDirectory, optional
        If provided, the resulting checkpoint is saved via
        ``run_dir.save_neutral(...)``. The cfg is also saved if not
        already present.
    max_bytes : int, optional
        Approximate cap on the size of the stored checkpoint. If a
        full-resolution checkpoint would exceed this, the stride
        is increased so only every K-th step is stored. Default 300 MB.
    verbose : bool, optional
        If True, print progress messages (number of steps, stride,
        memory estimate).

    Returns
    -------
    NeutralCheckpoint
        With shape ``(2N, num_stored_steps)`` for trajectory arrays.

    Raises
    ------
    NotImplementedError
        If ``cfg.custom_DFT_start`` is True (the DFT pre-fill is not
        yet implemented).
    """
    if rng is None:
        rng = np.random.default_rng(cfg.seed)

    # 1. Decide internal step count.
    num_internal_steps = _internal_step_count(cfg)

    # 2. Decide stride.
    stride, num_stored_steps = _decide_stride(
        cfg.num_molecules, num_internal_steps, max_bytes,
    )
    if verbose:
        est_mb = _estimate_checkpoint_bytes(
            cfg.num_molecules, num_stored_steps,
        ) / 1e6
        print(
            f"neutral propagation: "
            f"{num_internal_steps} internal steps, "
            f"stride={stride}, {num_stored_steps} stored steps, "
            f"~{est_mb:.1f} MB checkpoint"
        )

    # 3. Build initial state (allocates the storage-sized arrays).
    ckpt = build_initial_state(cfg, num_steps=num_stored_steps, rng=rng)

    # 4. Optional DFT pre-fill.
    if cfg.custom_DFT_start:
        raise NotImplementedError(
            "custom_DFT_start (TD-HeDFT initial conditions) is not "
            "yet implemented. Set cfg.custom_DFT_start=False to run "
            "without DFT pre-fill, or implement apply_dft_prefill()."
        )

    # 5. Inner loop.
    if num_internal_steps <= 1:
        # No propagation at all; only column 0 is meaningful.
        if run_dir is not None:
            _save_with_cfg(ckpt, cfg, run_dir)
        return ckpt

    state = state_from_checkpoint_column(ckpt, 0)
    prev_dist: np.ndarray | None = None
    next_storage_idx = 1

    for internal_id in range(1, num_internal_steps):
        new_state = neutral_propagation_step(
            state,
            cfg=cfg,
            mass_kg=ckpt.mass_kg,
            droplet_radii=ckpt.droplet_radii,
            prev_distance_angstrom=prev_dist,
            rng=rng,
        )
        prev_dist = _state_step_distance(state, new_state)
        state = new_state

        # Store every stride-th internal step.
        if internal_id % stride == 0 and next_storage_idx < num_stored_steps:
            write_state_to_checkpoint_column(state, ckpt, next_storage_idx)
            next_storage_idx += 1

    # If we ended up with fewer stored steps than allocated (rare; happens
    # when num_internal_steps - 1 is not divisible by stride), make sure the
    # last reachable column holds the final state. This keeps the trajectory
    # ending at the actual end-time rather than at a pre-stride snapshot.
    if next_storage_idx < num_stored_steps:
        write_state_to_checkpoint_column(state, ckpt, next_storage_idx)

    # 6. Save if run_dir given.
    if run_dir is not None:
        _save_with_cfg(ckpt, cfg, run_dir)

    return ckpt


# ===========================================================================
# Internal helpers
# ===========================================================================
def _internal_step_count(cfg: SimConfig) -> int:
    """Number of leapfrog steps to run.

    For ``single_pulse=True`` the MATLAB sets ``t_max = dt * 2`` -- the
    propagation freezes time near the moment of photoexcitation, so
    barely any time evolution happens. We mirror that here: 2 internal
    steps means we go from t=0 to t=dt (one column 0, one column 1).

    For ``single_pulse=False`` we use ``cfg.num_timesteps_neutral``
    (= ceil(t_max / dt)).
    """
    if cfg.single_pulse:
        return 2
    return cfg.num_timesteps_neutral


def _estimate_checkpoint_bytes(num_molecules: int, num_steps: int) -> int:
    """Estimate the in-memory size of a NeutralCheckpoint.

    Counts:

    * 10 trajectory arrays of shape (2N, num_steps) at 8 bytes/cell
    * Static (2N,) arrays for mass and droplet radii
    * Static (N,) arrays for r0 and E_initial
    * (num_steps,) for time_ps

    Pickle/npz overhead is small enough that this rough estimate is
    fine for stride decisions.
    """
    n_atoms = 2 * num_molecules
    bytes_2N_T = _NUM_2N_T_ARRAYS * n_atoms * num_steps * 8
    bytes_static_2N = 2 * n_atoms * 8                  # mass_kg, droplet_radii
    bytes_static_N = 2 * num_molecules * 8             # r0, E_initial_eV
    bytes_time = num_steps * 8
    return bytes_2N_T + bytes_static_2N + bytes_static_N + bytes_time


def _decide_stride(
    num_molecules: int,
    num_internal_steps: int,
    max_bytes: int,
) -> tuple[int, int]:
    """Return ``(stride, num_stored_steps)`` for the given budget.

    Storage stride is the smallest integer K such that storing every
    K-th step keeps the checkpoint under ``max_bytes``. Internal steps
    still happen at ``cfg.dt_neutral``; only storage is downsampled.

    The stored trajectory always includes step 0 (set by
    ``build_initial_state``), so ``num_stored_steps = ceil(num_internal_steps / stride)``.
    """
    full_size = _estimate_checkpoint_bytes(num_molecules, num_internal_steps)
    if full_size <= max_bytes:
        return 1, num_internal_steps

    # Roughly, num_stored * bytes_per_step + overhead <= max_bytes
    # bytes_per_step = 10 * 2N * 8
    bytes_per_step = _NUM_2N_T_ARRAYS * 2 * num_molecules * 8
    if bytes_per_step <= 0:
        return 1, num_internal_steps
    max_stored = max(2, int(max_bytes // bytes_per_step))
    stride = max(1, math.ceil(num_internal_steps / max_stored))
    num_stored = math.ceil(num_internal_steps / stride)
    return stride, num_stored


def _state_step_distance(
    prev: NeutralStepState,
    new: NeutralStepState,
) -> np.ndarray:
    """Per-atom distance traveled between two consecutive states."""
    return np.sqrt(
        (new.x - prev.x) ** 2
        + (new.y - prev.y) ** 2
        + (new.z - prev.z) ** 2
    )


def _save_with_cfg(
    ckpt: NeutralCheckpoint,
    cfg: SimConfig,
    run_dir: RunDirectory,
) -> None:
    """Save the checkpoint and the cfg in the run directory."""
    if not run_dir.has_cfg():
        run_dir.save_cfg(cfg)
    run_dir.save_neutral(ckpt)
