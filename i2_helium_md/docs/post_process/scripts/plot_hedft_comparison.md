# `scripts/plot_hedft_comparison.py` — a walkthrough

## What does this script do?

Reproduces the matplotlib equivalent of the two legacy MATLAB
post-processing scripts:

```
legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/
    simulation_image.m                    -- full output (3 plots)
    simulation_image_only_trajectories.m  -- subset (2 plots)
```

It opens two figures in interactive matplotlib windows for one MD
ion-stage run plus a HeDFT reference. Nothing is written to disk; save
manually from the matplotlib UI if needed.

## Usage

```
python scripts/plot_hedft_comparison.py
```

Edit the `USER SETTINGS` block at the top of the file — same
convention as `scripts/run_single_pulse.py`. There is no argparse
front-end on purpose: the script is meant to be tweaked in the
editor and rerun, not driven from the command line.

| Setting | Default | Effect |
|---|---|---|
| `RUN_DIR` | `data/runs/single_pulse_N_2000` | which MD run to read |
| `HEDFT_PATH` | `data/reference/9A_All_Data.csv` | HeDFT reference (use `18A_All_Data.csv` for the 18 Å droplet) |
| `VMI_HE_PATH`, `VMI_GAS_PATH` | the two existing references | only loaded when `SHOW_VMI_TILE=True` |
| `SHOW_VMI_TILE` | `True` | `False` mirrors `simulation_image_only_trajectories.m` (no bottom tile) |
| `MD_DISTANCE_STRIDE` | 50 | every Nth molecule on Figure 1 |
| `MD_VELOCITY_STRIDE` | 15 | every Nth atom on Figure 2 top tile (matches MATLAB) |
| `MASS_I_HE_AMU`, `MASS_I_HE2_AMU` | 131, 135 | mass filter for the simulation curves |
| `HIST_NUM_BINS`, `HIST_V_MAX_APS` | 120, 28 | bottom-tile histogram resolution |

## Figure 1 — distance trajectories

Single axis. Reproduces `simulation_image.m:54-91`.

- **HeDFT** R(t) — solid black, linewidth 1.5
- **MD trajectories** — every `MD_DISTANCE_STRIDE`th molecule's I-I
  separation, plotted with `color=(1.0, 0.2, 0.6, 0.1)` magenta-with-
  alpha. Only the first MD line is in the legend (matches MATLAB
  `HandleVisibility` toggle).
- xlim=[0, 6] ps, ylim=[8, 40] Å.

## Figure 2 — velocity panels

`GridSpec(2, 1)` when `SHOW_VMI_TILE=True`; otherwise a single axis.
Reproduces `simulation_image.m:93-268`.

### Top tile — velocity vs time (subplot label "a")
- **HeDFT** |v| from `hedft.v1_magnitude_Aps` — dotted, linewidth 2
- **MD velocity** — every `MD_VELOCITY_STRIDE`th atom's
  `sqrt(vx²+vy²+vz²)` over time, alpha-blue
- **mean MD velocity** — population mean across all atoms, dashed
- xlim=[0, 12] ps

### Bottom tile — velocity distribution at the detector (label "b")
Reproduces `simulation_image.m:159-262`. Four series, each
self-normalised:

| Series | Source | Style |
|---|---|---|
| I₂:I⁺ (gas) | `data/reference/vmi_iplus_gas.csv` | solid |
| I₂He_N:I⁺He (droplet) | `data/reference/vmi_iplus_he.csv` | dotted |
| simulation I⁺He (mass=131 amu) | `compute_final_velocity_histogram(ion, mass_amu=131.0)` | dashed |
| simulation I⁺He₂ (mass=135 amu) | `compute_final_velocity_histogram(ion, mass_amu=135.0)` | dash-dot |

Colors come from `plt.colormaps["viridis"]` sampled at five points —
a pragmatic stand-in for MATLAB's `colorcet('L08', 'N', 5)`.

xlim=[0, 28] Å/ps, ylim=[0, 1.1].

## Limitations / known differences from MATLAB

- **Display only.** No PDF or PNG export by default; the MATLAB
  `exportgraphics(..., 'simulation_image_inv.pdf')` is intentionally
  omitted.
- **Simulation traces are unsmoothed.** The MATLAB pipeline applies
  `movmean(..., 15)` to the simulated histograms; we plot the raw
  density. This is visible as bin-edge noise on the dashed/dash-dot
  curves at low atom counts.
- **No Bayesian uncertainty.** The MATLAB `bayes_hist` adds
  uncertainty bars; we skip them because the figure does not display
  them.
- **Color palette is viridis-based, not `colorcet('L08')`.** The
  ordering and meaning of the four colors are preserved; the exact
  hue mapping is not.
- **No mass/outside-flag warnings.** If a target mass yields zero
  matching atoms, `compute_final_velocity_histogram` raises
  `ValueError` and the script aborts; tweak `MASS_I_HE_AMU` /
  `MASS_I_HE2_AMU` (or the `mass_tolerance_amu` argument inside
  `velocity_distribution.py`) to recover.

## Smoke test

`tests/test_plot_hedft_comparison_smoke.py` drives the script via
`importlib`, monkeypatches `plt.show`, and asserts that two figures
are created with both `SHOW_VMI_TILE=True` and `False`. The test is
gated on the real run + reference files being present.
