# The `run_directory.py` module — a walkthrough

## What problem does this file solve?

`checkpoint.py` gave us low-level save/load functions, but each one takes a
raw `Path`. That's fine for one-off scripts but creates two real problems
for a real pipeline:

1. **Filename invention.** Every script picks its own conventions —
   `neutral.npz`, `Neutral.npz`, `neut.npz`, ... — and now two scripts
   that should agree don't.
2. **Path remembering.** The neutral-stage script writes to one path; the
   ion-stage script (possibly run minutes or weeks later, possibly by a
   different person) needs to find that same file.

`RunDirectory` solves both by treating one folder as the unit of a
simulation run. The user names the folder; we pick the filenames inside.

This addresses the exact pain point the legacy MATLAB had: `save('neutral_propagation_checkpoint',...)` followed by `cd()` calls to switch directories,
with no real concept of a "run" as a first-class artifact.

## Layout

```
<root>/
    cfg.json         the SimConfig that produced this run (JSON)
    neutral.npz      NeutralCheckpoint
    ion.npz          IonCheckpoint
    figures/         (optional) postprocess plots
    logs/            (optional) text logs
```

A run directory is **self-describing**: anyone can open it six months
later and see the config that produced the data.

## Public API

```python
from i2_helium_md.simulation.run_directory import RunDirectory

run = RunDirectory("data/runs/my_run")     # does NOT create the dir yet

# Save things (creates dir lazily)
run.save_cfg(cfg)
run.save_neutral(neutral_ckpt)            # convention-named: neutral.npz
run.save_ion(ion_ckpt)                    # convention-named: ion.npz

# Inspect
run.exists()                              # is the dir there?
run.has_cfg(); run.has_neutral(); run.has_ion()

# Load
cfg = run.load_cfg()
neutral = run.load_neutral()              # auto-validates against cfg.json
ion = run.load_ion()
```

## Auto-saved cfg

When you call `run.save_neutral(ckpt, cfg=cfg)` with a cfg argument, the
cfg is also saved if `cfg.json` doesn't already exist. This is convenient
for the typical first-write case:

```python
run = RunDirectory("data/runs/run01")
run.save_neutral(neutral_ckpt, cfg=cfg)   # saves both neutral.npz and cfg.json
```

The cfg is **not overwritten** if it already exists — we never silently
mutate the config of a running pipeline. This matters: if the ion-stage
script is run with a slightly different cfg by accident, we want loud
failure, not silent corruption.

## Validation chain

When you call `run.load_neutral()` with no args, the loader:

1. Looks for `cfg.json` in the run directory.
2. If found, parses it into a `SimConfig`.
3. Passes that cfg to `load_neutral_checkpoint()` for shape validation.

If you pass an explicit `cfg=...` to `load_neutral()`, that cfg takes
precedence — useful when you want to validate against a cfg that's
different from the one saved (e.g. you renamed a field and want to check
old runs against the new schema).

## Two-process pipeline (the original motivating use case)

```python
# ----- script A: run_neutral.py -----
cfg = single_pulse_N2000(num_molecules=500, seed=42)
run = RunDirectory("data/runs/test01")
run.save_cfg(cfg)
neutral_ckpt = run_neutral_propagation(cfg)
run.save_neutral(neutral_ckpt)

# ----- script B: run_ion.py (later, possibly different process) -----
run = RunDirectory("data/runs/test01")
cfg = run.load_cfg()                        # picks up the saved config
neutral = run.load_neutral()                # validated against cfg.json
ion_ckpt = run_ion_propagation(cfg, neutral)
run.save_ion(ion_ckpt)

# ----- script C: postprocess.py -----
run = RunDirectory("data/runs/test01")
ion = run.load_ion()
make_vmi_image(ion)
```

All three scripts only need to agree on the **directory path** —
`"data/runs/test01"`. Filenames are an implementation detail.

## Departures from MATLAB

1. **No global current-working-directory dependency.** MATLAB `save('foo')`
   wrote to wherever `cd()` had last been called. We always use absolute or
   explicit relative paths, never CWD.

2. **No file-name branching by config flag.** MATLAB used different
   filenames for `effusive_dynamics` (`ion_propagation_checkpoint_gas`)
   vs `single_droplet_size && single_initial_position`
   (`ion_propagation_checkpoint_hedft`). We use **the same filenames**
   inside different run directories. The variant lives in `cfg.json`,
   not in the filename.

3. **Self-describing runs.** MATLAB had no equivalent of `cfg.json`; if you
   came back to a checkpoint file later, you couldn't reconstruct what
   produced it. Our run directories contain that information.

## Implementation choices worth highlighting

### Lazy directory creation

```python
run = RunDirectory("data/runs/foo")    # filesystem untouched
```

The directory is created when the first `save_*` is called. This is so
that:
- Building a `RunDirectory` for inspection of an existing run is
  side-effect free.
- We don't litter the filesystem with empty directories from typoed paths.

### Cfg unknown-field check

When loading `cfg.json`, any field in the JSON that is not a `SimConfig`
field raises `ValueError` with the offending field names. This catches
the case where a run was produced by a different version of the code
and the schema has drifted. The alternative (silently dropping the
field) would hide real bugs.

### `save_neutral` accepts `cfg=None`

The signature is `save_neutral(ckpt, cfg=None)` rather than requiring cfg,
because:
- For unit tests and inspection scripts you may not have a cfg handy.
- For the common case where you've already called `save_cfg()` separately,
  the second cfg argument would be redundant.

The convenience pattern (`save_neutral(ckpt, cfg=cfg)`) is for scripts
that want to do everything in one call.

## Future extensions

1. **Variant runs (pump-probe).** Add `save_pump(ckpt)` and `save_probe(ckpt)`
   methods. The layout becomes
   `<root>/{cfg.json, pump.npz, probe.npz, ion.npz}`.

2. **Postprocess outputs.** Add `figures_dir`, `logs_dir` properties that
   return `<root>/figures/`, `<root>/logs/` and create them on demand.

3. **Run hashing / naming.** A helper `RunDirectory.from_cfg(cfg, base="data/runs")`
   that hashes the cfg into a deterministic name. Two runs with identical
   cfgs would land in the same directory; this could either be a feature
   (idempotency) or a footgun (silent overwrites). Worth thinking through
   before adding.

4. **Run discovery.** `RunDirectory.list_runs("data/runs")` to enumerate
   all existing runs in a parent directory. Useful for batch postprocess
   scripts.
