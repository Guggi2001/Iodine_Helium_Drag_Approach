# Tier-0 scripts

Tier 0 runs one law at fixed mass, without noise, and compares it with the TDDFT
traces. The laws are the three thesis bundles: `LQ` (B2), `LQ_eff` (B3) and `CC`
(h405). See `docs/drag_port/Tier0/TIER0_FINDINGS.md` for the history of the
verdict and `CLAUDE.md` for the thesis result.

```
tier0_common.py                  -> cfg builders, run naming, scoring windows
gen_tier0_runs.py                -> one run for a CASE / LAW (data/runs, gitignored)
tier0_drag_comparison.py         -> score one run, show its figures
tier0_calibration_figures.py     -> regenerate and save the five thesis figures
tier0_thesis_figure.py           -> series + scores CSVs of B2 / B3 / h405
tier0_law_figure.py              -> the law and ridge figures from those CSVs
```

## `scripts/tier0_common.py`

- `build_tier0_cfg(case, law, ...)`: the Tier-0 config of a law. `CC` comes capped
  at the 3.5 Å gate, as in the thesis, under `allow_inconsistent_gate_pairing`.
  `ridge_c` takes a and E_bind from a `calibration_ridge` point instead (`LQ` only).
- `build_drag_cfg`: the layer below, with free coefficient, binding and mass
  overrides. The Tier-2 generators use it too.
- `run_dir_name`: the run directory both sides agree on. A renamed bundle keeps
  its old id in the name, so existing runs are still found.
- `TIER0_WINDOWS_PS`: scoring window per droplet case (9 Å 2.67–14.081 ps,
  18 Å 4.54–14.768 ps).

## `scripts/gen_tier0_runs.py`

Edit USER SETTINGS (`CASE`, `LAW`, `RIDGE_C`, `N`, `ION_TIME_PS`, `E_BIND_OVERRIDE`,
`M_EFF_AMU`, `RUN_TAG`), then run it. The default mass is I⁺ + 21 He (G11.2). Give
a run with changed settings its own `RUN_TAG`, otherwise it replaces the earlier
run of that name.

## `scripts/post_processing/tier0_drag_comparison.py`

Set the same `CASE`, `VARIANT`, `N` and `RUN_TAG` as the generator. It prints the
in-window scores against the HeDFT reference and the smoothed trace, and shows the
figures in a window. It does not save them.

## `scripts/post_processing/tier0_calibration_figures.py`

No settings to edit: its `RUNS` table holds the three runs and five figures. A run
is reused if its `cfg.json` matches and regenerated otherwise. PNGs are half the text width, in thesis style,
and go to `data/figures/thesis_microscopic_calibration/`; pass another folder as
the first argument to try something out without touching them.
