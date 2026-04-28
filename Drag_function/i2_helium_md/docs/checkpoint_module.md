# The `checkpoint.py` module — a walkthrough

## What problem does this file solve?

The MATLAB pipeline runs in two stages: **neutral propagation** (atoms before
the laser pulse) and **ion propagation** (after ionization). The two stages
live in separate scripts and pass data via disk:

```matlab
% end of neutral script:
save('neutral_propagation_checkpoint', 'positions', 'velocities', ...)

% start of ion script:
load('neutral_propagation_checkpoint');
```

We need a Python equivalent that:
1. Saves the neutral-stage state to disk.
2. Loads it cleanly into the ion-stage driver.
3. Saves the final ion-stage state for postprocessing.

This module provides exactly that.

## Position in the dependency chain

```
config.py
   ↓
physics/* (all)
sampling/* (all)
   ↓
simulation/checkpoint.py    ← THIS MODULE
   ↓
simulation/neutral.py       (saves a NeutralCheckpoint)
simulation/ion.py           (loads a NeutralCheckpoint, saves IonCheckpoint)
   ↓
postprocess/* (loads IonCheckpoint)
```

`checkpoint.py` is the I/O boundary between simulation stages.

## Public API

```python
from i2_helium_md.simulation.checkpoint import (
    NeutralCheckpoint,
    IonCheckpoint,
    save_neutral_checkpoint,
    load_neutral_checkpoint,
    save_ion_checkpoint,
    load_ion_checkpoint,
)
```

### Typical usage

```python
# At end of neutral propagation
ckpt = NeutralCheckpoint(
    num_molecules=cfg.num_molecules,
    time_ps=time_axis,
    positions_x=x_traj,
    ...
)
save_neutral_checkpoint(ckpt, "data/runs/run01/neutral.npz")

# At start of ion propagation
ckpt = load_neutral_checkpoint("data/runs/run01/neutral.npz", cfg=cfg)
# ckpt.positions_x[:, -1] gives the final neutral positions, etc.
```

## What's saved

### `NeutralCheckpoint`

The minimum sufficient set to:
- Continue with ion propagation (positions, velocities, masses, droplet radii)
- Reproduce energy diagnostics (E_kin, E_pot, E_initial, E_dissip, L_droplet)
- Run postprocess on the neutral trajectory alone (time axis, r0)

### `IonCheckpoint`

The minimum sufficient set to:
- Generate VMI images and momentum spectra (`positions_final_*`,
  `velocities_final_*`, `mass_final_kg`)
- Reproduce ion-stage energy diagnostics
- Track collisional history (`number_of_collisions`,
  `relative_loss_per_ps`, `b_ion_outside`)

## What is **not** saved

- **Physical constants** (`eV`, `u`, `K_B`, ...) — they live in `constants.py`
  and are imported by readers.
- **Config flags** (`single_pulse`, `effusive_dynamics`, `Xdip_active`, ...) —
  these come from `cfg`. The loader can validate against a passed `cfg`.
- **Mode-specific scalars** (`he_direction_scattering`, `binding_energy_*`,
  ...) — derived from `cfg`.
- **Plot-only diagnostics** that can be recomputed from positions/velocities.

This is a deliberate design choice: a checkpoint is **state**, not
**configuration**. Reproducing constants from `constants.py` and recovering
config from `cfg` keeps the file small and prevents drift between stages
(e.g. an old checkpoint with a stale `eV` constant).

## File format: `.npz`

NumPy's native `.savez_compressed`. Pros:
- One file per checkpoint, with all named arrays inside.
- No extra dependencies (no HDF5, no pickle).
- Forward-compatible: loaders use explicit field names and can ignore extras.
- ZIP-deflate compression, which is effective for smooth physical
  trajectories (typically 30–60% size reduction).

The loader uses `allow_pickle=False` for security — refuses to load files
that contain Python pickles.

## Schema versioning

Every checkpoint carries a `schema_version: int` field. The loader checks it
on read and refuses to load incompatible versions.

**Rule for bumping the version:**
- Adding a field: backward-compatible — bump *not required*.
- Removing or renaming a field: bump version.
- Changing the meaning or units of an existing field: bump version.

When bumping, update `_NEUTRAL_SCHEMA_VERSION` or `_ION_SCHEMA_VERSION` in
`checkpoint.py` and document the change in `migration_log.md`.

## Validation against `SimConfig`

`load_*_checkpoint(path, cfg=...)` cross-checks the checkpoint against a
config:
- `num_molecules` must match.
- 2N-shaped arrays (`mass_kg`, `droplet_radii`) must have shape `(2*cfg.num_molecules,)`.

If `cfg=None`, no validation is performed against config (useful for
inspection scripts that don't want to bring up a full SimConfig).

## Departures from MATLAB

1. **Two clean dataclasses instead of bare variable dumps.** MATLAB's `save`
   serializes whatever is in scope. Our dataclasses make the contract
   explicit.

2. **No constants saved.** MATLAB saved `eV`, `u`, `mass`, `m`, ... — all
   recoverable from current code.

3. **No flags saved.** MATLAB saved `effusive_dynamics`, `hard_sphere_collision_mode`,
   ... — these belong to `cfg`, not state.

4. **Path is explicit.** MATLAB hardcoded `'neutral_propagation_checkpoint'`
   and used `cd()` to change directories. Python takes a `Path` argument and
   creates parent directories as needed.

5. **Schema version + validation.** MATLAB had no version field; an old
   checkpoint silently produced wrong results if the variable list changed.
   We fail loudly.

6. **`.npz` instead of `.mat`.** Avoids a `scipy.io` dependency and integrates
   better with NumPy.

## Common pitfalls

1. **Don't pass `mass` in amu.** It's `mass_kg`, in kilograms. The dataclass
   field name is the documentation.

2. **Don't save trajectories without sub-sampling.** A 100-fs simulation at
   dt=1 fs has 100,000 timesteps, and `(2N, num_steps)` arrays scale fast.
   Decide whether to save full trajectories or only every Nth timestep
   *before* allocating the dataclass.

3. **Don't write to the project's source tree.** The convention is
   `data/runs/<run_name>/neutral.npz` (or `ion.npz`). Add `data/runs/`
   to `.gitignore`.

4. **The loader's `cfg` is optional but recommended.** Skipping it means a
   shape mismatch only manifests downstream as a confusing broadcasting
   error.

## Future extensions

1. **Resumable mid-stage checkpointing.** Currently a checkpoint is "end of
   stage". For long ion runs we may want to save every K timesteps and resume
   from there. Adding a `step_index` field is the natural extension.

2. **Streaming format.** For very large runs the dataclass approach loads
   everything into memory. A future version could expose iterator access:
   `for batch in stream_neutral_checkpoint(path, batch_size=100): ...`.

3. **Metadata block.** A `metadata: dict[str, Any]` field could carry the
   `cfg` itself (serialized) so the checkpoint is fully self-describing
   without an external `cfg`. Trade-off: tightly couples the format to the
   `SimConfig` dataclass.
