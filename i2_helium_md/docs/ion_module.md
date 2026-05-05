# The `simulation/ion.py` module

## What problem does this file solve?

This is the top-level driver for ion propagation. Given a `SimConfig`
and a finished `NeutralCheckpoint`, it runs a complete ion-stage
simulation and produces an `IonCheckpoint`. It orchestrates the three
pieces built in Steps 11a/b/c:

1. **`build_initial_ion_state`** — sets up the ion-stage t=0 column
   from the last column of the neutral checkpoint and pre-allocates
   the ion trajectory arrays.
2. **`ion_propagation_step`** (pure) — advances one ion timestep
   (leapfrog + collisions + mass attachment + energy bookkeeping).
3. **`RunDirectory`** — optional checkpoint serialization.

It also handles three run-time concerns the lower-level pieces don't:

- **Auto-stride** for memory-bounded checkpoints: if the full-resolution
  trajectory would exceed ~1 GB, only every K-th internal step is
  stored. Internal stepping still happens at `cfg.dt_ion` so the
  physics (collision rate, leapfrog stability) is unchanged.
- **Final-state population**: the post-loop `positions_final_*`,
  `velocities_final_*`, `mass_final_kg`, and `b_ion_outside` fields are
  filled from the *actual last internal step*, not from the last stored
  column (under stride > 1 these can differ).
- **Driver-level scope check**: refuses `cfg.single_pulse=False` up
  front, since that branch needs the MATLAB dt-switching scheme
  (pump-probe, out of scope for this stage).

## Public API

```python
from i2_helium_md.simulation.ion import run_ion_propagation
from i2_helium_md.simulation.run_directory import RunDirectory

# Simple call: returns an IonCheckpoint
ion_ckpt = run_ion_propagation(cfg, neutral_ckpt, rng=np.random.default_rng(42))

# With saving:
rd = RunDirectory("runs/my-run-1")
ion_ckpt = run_ion_propagation(cfg, neutral_ckpt, run_dir=rd, verbose=True)

# Custom memory budget (default 1 GB):
ion_ckpt = run_ion_propagation(cfg, neutral_ckpt, max_bytes=100_000_000)
```

The `neutral_ckpt` argument is required (positional). The driver always
starts from its last column (matches the production single-pulse use
case where ionization happens at the end of the neutral stage).

### Two-process pipeline

The expected end-to-end pipeline is two scripts sharing a `RunDirectory`:

```python
# Script A: run_neutral.py
from i2_helium_md.presets import single_pulse_N2000
from i2_helium_md.simulation.neutral import run_neutral_propagation
from i2_helium_md.simulation.run_directory import RunDirectory

cfg = single_pulse_N2000(num_molecules=500, seed=42)
rd = RunDirectory("data/runs/test01")
run_neutral_propagation(cfg, run_dir=rd)
```

```python
# Script B: run_ion.py (later)
from i2_helium_md.simulation.ion import run_ion_propagation
from i2_helium_md.simulation.run_directory import RunDirectory

rd = RunDirectory("data/runs/test01")
cfg = rd.load_cfg()
neutral = rd.load_neutral(cfg=cfg)
run_ion_propagation(cfg, neutral, run_dir=rd)
# data/runs/test01/ion.npz now exists alongside cfg.json and neutral.npz
```

## Single-pulse vs non-single-pulse

The `cfg.single_pulse` flag controls how the ion stage runs — but the
semantics are the **opposite** of the neutral stage:

- **`single_pulse=True`** (default and the only supported mode here):
  the ion stage runs the full
  `num_internal = ceil(cfg.ion_simulation_time / cfg.dt_ion)`
  internal steps at fixed `dt_ion`. For the production preset
  (20 ps / 0.01 ps) this is 2000 steps.
- **`single_pulse=False`**: the driver raises `NotImplementedError`.
  The MATLAB non-single-pulse branch uses a dt-switching scheme
  (`dt_fine`/`dt_coarse` + `switchtime`) to handle the long-time
  evolution after a pump-probe sequence. That scheme is out of scope
  for this stage of the project.

This matches the MATLAB:

```matlab
if single_pulse
    ion_simulation_time = 20;
    dt = 0.01;
    ion_timesteps = ceil(ion_simulation_time/dt);
else
    % dt_fine / dt_coarse / switchtime branch -- not ported
end
```

## Auto-stride for memory budget

If the full-resolution checkpoint would exceed `max_bytes`, the driver
chooses an integer storage stride `K` such that storing every K-th
internal step keeps the checkpoint within the budget:

```
stride = ceil(num_internal_steps / max_storable_steps)
num_stored_steps = ceil(num_internal_steps / stride)
```

The internal loop **always** runs at full resolution (`dt_ion`).
Two reasons this matters more for ion than neutral:

1. The Mode-3 collision sampler relies on the previous internal step's
   per-atom displacement; skipping internal steps would corrupt the
   collision rate.
2. The mass-attachment RNG draw happens **every** internal step. With
   stride > 1 the stored `mass_history_kg` columns skip intermediate
   attachments, but the underlying step state still accumulates them
   correctly — so the mass at the last stored column is the true mass
   at that time, not "the mass after K MATLAB steps assuming none of
   the skipped attachments fired".

The state itself is carried as an `IonStepState` dataclass between
internal calls; only every K-th state is copied into the checkpoint
arrays.

### Memory estimates

Per-step bytes are higher than neutral because `IonCheckpoint` carries
12 `(2N, T)` arrays vs. 10 for neutral (`mass_history_kg`,
`relative_loss_per_ps`, `number_of_collisions` are ion-only).

| Run | N | t_ion | Internal steps | Full size | Stride | Stored | Saved size |
|---|---|---|---|---|---|---|---|
| Single-pulse demo | 4 | 0.1 ps | 10 | <0.1 MB | 1 | 10 | <0.1 MB |
| HeDFT 9 Å | 500 | 20 ps | 2000 | ~190 MB | 1 | 2000 | ~190 MB |
| HeDFT 9 Å | 2000 | 20 ps | 2000 | ~770 MB | 3 | 667 | ~256 MB |
| Long single-pulse | 2000 | 50 ps | 5000 | ~1920 MB | 7 | 715 | ~275 MB |

## Final-state fields

After the inner loop, the driver populates the post-loop fields of
`IonCheckpoint` from the *actual last internal step*:

```python
ckpt.positions_final_x[:] = state.x
ckpt.positions_final_y[:] = state.y
ckpt.positions_final_z[:] = state.z
ckpt.velocities_final_x[:] = state.vx
ckpt.velocities_final_y[:] = state.vy
ckpt.velocities_final_z[:] = state.vz
ckpt.mass_final_kg[:]      = state.mass_kg

# Per-molecule "did either ion of this pair end up outside the droplet?"
depth_final = sqrt(x^2 + y^2 + z^2) - droplet_radii_angstrom
ckpt.b_ion_outside[i] = (depth_final[atom1_i] > 0) | (depth_final[atom2_i] > 0)
```

When `stride == 1` the last stored column equals these final fields.
Under `stride > 1` they can differ: `positions_final_*` is from the
true end of the run, while `positions_x[:, -1]` is the most recent
strided snapshot.

`relative_loss_per_ps` is **not** populated by the driver; it stays as
zeros from `build_initial_ion_state`. Its physical meaning under our
hard-sphere mode-3 scope is undefined (the legacy MATLAB only uses it
in a different energy-loss model that we don't implement). It is
reserved for postprocessing.

## RunDirectory integration

If `run_dir` is provided:

1. The cfg is saved to `cfg.json` if not already present.
2. The checkpoint is saved to `ion.npz`.

This is the same auto-saving cfg behavior the neutral driver has, so
the two-process pipeline above works without manual cfg duplication.
The driver does NOT re-validate the cfg before saving — that's the
job of `load_ion_checkpoint(path, cfg=...)` on read.

## Scope checks

The driver fails fast if cfg requests an unsupported feature:

| Flag | Where it raises | Exception |
|---|---|---|
| `single_pulse=False` | `_check_scope_ion_driver` (driver) | `NotImplementedError` |
| `effusive_dynamics=True` | `build_initial_ion_state` (build) | `NotImplementedError` |
| `single_charge_ionization_allowed=True` | `build_initial_ion_state` | `NotImplementedError` |
| `additional_droplet_charges > 0` | `build_initial_ion_state` | `NotImplementedError` |
| `highly_charged_iodine=True` | `build_initial_ion_state` | `NotImplementedError` |
| `hard_sphere_collision_mode != 3` | `ion_propagation_step` (first step) | `ValueError` |

All driver-level and build-time checks fire before any expensive
stepping. The collision-mode check fires on the first inner-loop call
to `ion_propagation_step`.

## What's inside

```
run_ion_propagation(cfg, neutral_ckpt, *, rng, run_dir, max_bytes, verbose)
├─ default_rng(cfg.seed) if rng is None
├─ _check_scope_ion_driver(cfg)           -- single_pulse=False -> raise
├─ _internal_step_count_ion(cfg)          -- ceil(ion_sim_time / dt_ion)
├─ _decide_stride_ion(N, internal_steps, cap)
├─ build_initial_ion_state(cfg, neutral_ckpt, num_steps_ion=stored, start_id=-1)
├─ inner loop:
│    state  = ion_state_from_checkpoint_column(ckpt, 0)
│    charge = ones(2N)                    -- allocated once
│    for internal_id in 1..num_internal_steps-1:
│        new = ion_propagation_step(state, ..., prev_distance, rng)
│        prev_distance = |new - state|
│        if internal_id % stride == 0:
│            write_ion_state_to_checkpoint_column(new, ckpt, next_storage_idx++)
│        state = new
│    if next_storage_idx < num_stored_steps:    -- defensive; unreachable
│        write_ion_state_to_checkpoint_column(state, ckpt, next_storage_idx)
├─ _write_final_state(state, ckpt, cfg)
│    ├─ positions_final_*, velocities_final_*, mass_final_kg <- state
│    └─ b_ion_outside <- (atom1_outside | atom2_outside)
├─ _save_with_cfg_ion(ckpt, cfg, run_dir)  if run_dir is not None
└─ return ckpt
```

The two `ion_state_from_checkpoint_column` /
`write_ion_state_to_checkpoint_column` helpers live in
`simulation/ion_propagation_step.py` next to `IonStepState`, mirroring
the neutral helpers' placement next to `NeutralStepState` in
`simulation/propagation_step.py`.

## Tests

**25 pytest tests** in `test_ion.py` covering:

- API contract (returns `IonCheckpoint`, output shapes match schema,
  time axis consistency, final state stored when `stride=1`)
- Mass history (column 0 equals neutral mass, monotonic non-decreasing,
  `mass_final_kg == mass_history_kg[:, -1]` for `stride=1`)
- Reproducibility (same explicit RNG seed → bit-identical checkpoint)
- Stride math (small/large cases, end-to-end strided run with tight
  `max_bytes` produces a checkpoint smaller than the cap and still
  ending at a non-initial state)
- RunDirectory round trip and "do not overwrite existing cfg"
- Scope checks (all six unsupported flags + collision mode)
- Energy bookkeeping (dissipation non-decreasing; with
  `mass_attach_probability=0.0`, total energy `E_kin + E_pot + E_dissip`
  is conserved to <1% over a short run)
- Internal helpers (step count, byte estimate scales linearly)

The most informative regression guard is **energy bookkeeping with
attachment disabled**: with `mass_attach_probability=0.0` the only
energy sink is collisional dissipation captured in `E_dissip_eV`, so
the sum should be conserved up to leapfrog symplectic error
(~ppm/step). The 1% tolerance leaves headroom for the driver test's
short run length.

## Departures from MATLAB

1. **Auto-stride is new.** MATLAB stores every internal ion step and
   only downsamples after the loop with `reduction_timesteps`. The
   Python driver decides the stride up front, allocates the
   appropriately sized checkpoint, and never holds the full-resolution
   trajectory in memory.
2. **`single_pulse=False` raises rather than running pump-probe.** The
   MATLAB dt-switching branch is intentionally unimplemented; raising
   loudly is preferable to silently using a single `dt_ion` for a
   different physics scenario.
3. **t=0 energy bookkeeping fixes are preserved.** `build_initial_ion_state`
   already corrects two MATLAB t=0 bugs (E_kin missing v_z and squared,
   E_pot using 2D radius and missing partner Coulomb). The driver does
   not touch these — it just runs `ion_propagation_step` from t=0
   onward, where the MATLAB code already used the correct formulas.
4. **`b_ion_outside` shape change.** MATLAB stores per-atom `(2N,)`;
   the Python `IonCheckpoint` schema declares per-molecule `(N,)`. The
   driver reduces with a per-molecule **OR** rule: the molecule counts
   as having an ion outside if either of its two atoms has
   `depth > 0` at the final state.
5. **`relative_loss_per_ps` left as zeros.** MATLAB only ever saves a
   scalar `relative_ion_energy_loss_per_ps` config value, but the
   Python schema declares `(2N, T)`. Under our hard-sphere mode-3
   scope this field has no defined per-step semantics, so the driver
   leaves it zero and defers any meaningful population to postprocess.
6. **Final state always stored.** Same guarantee as the neutral driver:
   `positions_final_*` etc. are populated from the actual last internal
   step, even if it's not aligned with the storage stride.