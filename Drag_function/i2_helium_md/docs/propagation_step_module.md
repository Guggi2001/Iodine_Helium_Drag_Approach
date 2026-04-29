# The `simulation/propagation_step.py` module

## What problem does this file solve?

This module implements the **per-timestep advance** for the neutral
propagation, as a **pure function**: input state → output state, no
side effects.

A `NeutralStepState` is the minimum-sufficient state to advance:
positions, velocities, cumulative energy/path bookkeeping, and time.
The function reads one and returns a new one.

The driver (`simulation/neutral.py`) orchestrates: it calls the step
function in a loop, threading the previous-step distance for the
collision sampler, and copies the resulting state into the
`NeutralCheckpoint` storage at the appropriate stride.

This split (pure step + driver) lets us **decouple internal stepping
from storage stride**, which is essential for memory-bounded runs
(see `simulation/neutral.py` for the auto-stride logic).

It replaces the body of the `while t_id < num_timesteps` loop in
`vmi_sim_3d_neutral_propa_HeDFT_mimic.m` (lines ~536-913) excluding
the `attach_he` mass-attachment branch (out of scope for the neutral
stage).

## Public API

```python
from i2_helium_md.simulation.propagation_step import (
    NeutralStepState,
    neutral_propagation_step,
    state_from_checkpoint_column,
    write_state_to_checkpoint_column,
)

# Bootstrap from a checkpoint:
state = state_from_checkpoint_column(ckpt, 0)

# Advance one step:
new_state = neutral_propagation_step(
    state,
    cfg=cfg,
    mass_kg=ckpt.mass_kg,
    droplet_radii=ckpt.droplet_radii,
    prev_distance_angstrom=None,   # None for the first step
    rng=np.random.default_rng(42),
)

# Optionally persist:
write_state_to_checkpoint_column(new_state, ckpt, t_id=1)
```

## What's inside (step by step)

### Step 1: Leapfrog integration
Uses `physics.leapfrog.make_neutral_step(cfg, mass, droplet_radii)` to
get a closure that handles per-pair Morse + per-atom droplet
acceleration. Returns candidate `(x1, y1, z1, vx1, vy1, vz1)` and the
per-pair Morse potential `E_pot_partner` (shape `(N,)`, eV).

### Step 2: Depth into droplet
`r1 = sqrt(x1² + y1² + z1²)`, `depth = r1 - droplet_radii`.
`depth < 0` means the atom is inside the droplet.

### Step 3: Pre-collision energy E0
Computed from post-leapfrog velocities. Used by the collision sampler
to enforce the Landau cutoff `E_min_eV`.

### Step 4: Collision event sampling (Mode 3)
**Important:** the Mode 3 probability uses the *previous* step's
displacement, passed in as `prev_distance_angstrom`.

If `prev_distance_angstrom is None` (first step of a run), no
collisions are sampled. This matches MATLAB's `if t_id > 1`.

The collision sampler also rejects atoms outside the droplet
(`depth >= 0`) and below the Landau cutoff (`E0 < E_min`).

### Step 5: Apply collisions
For atoms with `b_collision=True`, replace `(vx1, vy1, vz1)` with
elastically-scattered velocities and accumulate the per-atom energy
loss `dE = E0 - E1`. Non-colliders pass through unchanged.

### Step 6: Energy diagnostics
- **E_kin**: from post-collision velocities, per atom in eV.
- **E_pot**: per-atom droplet solvation + half of per-pair Morse.
  The 50/50 split between atoms of a pair matches MATLAB.

### Step 7: Cumulative bookkeeping
- **E_dissip**: prev + per-atom dE.
- **L_droplet**: prev + per-atom step length (only if ended inside
  the droplet, matching MATLAB).

### Step 8: Assemble new state
A new frozen `NeutralStepState` is returned. Inputs are not mutated.

## NeutralStepState dataclass

| Field | Shape | Units | Notes |
|---|---|---|---|
| `x, y, z` | `(2N,)` | Å | Cartesian positions |
| `vx, vy, vz` | `(2N,)` | Å/ps | Velocities |
| `E_kin_eV` | `(2N,)` | eV | Per-atom KE |
| `E_pot_eV` | `(2N,)` | eV | Per-atom droplet + half partner |
| `E_dissip_eV` | `(2N,)` | eV | Cumulative per-atom dissipation |
| `L_droplet_eV_ps` | `(2N,)` | Å | Cumulative path inside droplet (legacy MATLAB field name) |
| `time_ps` | scalar | ps | Time at this state |

The dataclass is frozen, so the function signature documents the
purity contract. Underlying NumPy arrays are still mutable in
principle; the function nonetheless does not mutate them.

## Helper functions

`state_from_checkpoint_column(ckpt, t_id)`: extract a
`NeutralStepState` by **copying** column `t_id` arrays from a
`NeutralCheckpoint`. Mutations on the returned state cannot affect
the checkpoint.

`write_state_to_checkpoint_column(state, ckpt, t_id)`: copy a
`NeutralStepState` into column `t_id` of a `NeutralCheckpoint`.

## Mode 3 only

The `cfg.hard_sphere_collision_mode` field exists for backward
compatibility, but only Mode 3 is implemented. Modes 1 and 2 raise
`ValueError`. They were never used in production runs and we
deliberately did not port them.

## Tests

- **11 pytest tests** in `test_propagation_step.py`: validation,
  purity, time advance, no-collisions-on-first-step, atoms-outside,
  cumulative monotonicity, reproducibility, energy bookkeeping,
  state-helper round trip.
- **14 sandbox checks** in `smoke_test_propagation_step.py`.

## Departures from MATLAB

1. **Pure function instead of in-place on a 2D array.** Decouples the
   physics step from the storage stride.
2. **Mode 1, 2 not ported.** Raise `ValueError`.
3. **No `attach_he`, `effusive_dynamics`, `relative_energy_loss`, `DEBUG`.**
   All deliberately not ported (out of scope or unused in production).
4. **`prev_distance` is an explicit argument**, not derived from a
   trajectory array. The driver tracks it as it loops.
