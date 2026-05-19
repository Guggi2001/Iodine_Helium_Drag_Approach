# `plot_paper_v3.py` paper-v3 VMI comparison

This script ports the active non-effusive droplet branch of
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_v3.m`.
It is a focused paper-figure reproduction, not a general run-summary
diagnostic.

The goal is to compare the simulation to the same VMI observables used by the
legacy MATLAB paper script:

1. detector-plane radial velocity distribution,
2. detector-plane azimuthal angle distribution,
3. final ion mass population as a separate diagnostic.

## Outputs

`plot_paper_v3.py` writes two groups of files under
`<run>/figures/`:

| Output | Purpose |
|---|---|
| `compare_simulation_and_measurement.pdf` / `.png` | Main two-panel VMI comparison. |
| `ion_mass_histogram.pdf` / `.png` | Separate final-mass diagnostic, matching the extra figure opened at the end of the MATLAB script. |

The main figure has two stacked panels. The top panel compares radial velocity
curves. The bottom panel compares angular `phi` curves.

## Top panel: radial velocity distribution

The top panel asks:

> Do the simulated ions produce the same detector-plane speed distribution as
> the experimental I+He VMI measurement?

The x-axis is velocity in `m/s`. The y-axis is normalized signal in arbitrary
units, so this is a shape comparison rather than an absolute-yield comparison.

Curves:

| Curve | Meaning |
|---|---|
| `I2:I+He TS (296:297)` | Experimental timescan radial VMI reference from MATLAB `mean_timescan_2d_VMI([296:297], false, [524.5297 380.8430], false)`. |
| `I2:I+He (43563)` | Experimental high-SNR processed I+He radial VMI reference from the MATLAB `res_sum` path; 43563 is the I+He measurement assignment kept in the v3 provenance. |
| `simulated v.distr. m=127` | Simulated bare iodine ion channel, approximately `I+`. |
| `simulated v.distr. m=131` | Simulated `I+He` channel. |
| `simulated v.distr. m=135` | Simulated `I+He2` channel. |

The simulation recipe follows the unified Strategy A pipeline:

```text
round(mass / u) == mass_select
b_ion_outside == true
v_projected = sqrt(vx^2 + vy^2)
edges_velocity = 0:0.05:35  A/ps
x_plot = centers_velocity * 100  m/s
smooth = movmean(histogram_counts, 15)
y_plot = smooth / max(smooth)
```

The important point is the projection: `sqrt(vx² + vy²)` intentionally
ignores `vz`, because the paper-v3 figure compares to a VMI detector
projection. See `post_processing_strategy.md` §3 for the Strategy A
framing and §5 for the cross-script recipe table.

## Bottom panel: phi distribution

The bottom panel asks:

> Do the simulated ions leave in the same azimuthal directions as the
> processed experimental VMI image?

Curves:

| Curve | Meaning |
|---|---|
| `I2:I+He (43563)` | Experimental angular distribution from the processed polar VMI image. MATLAB averages `res.image_polar` over the selected radial range and normalizes the result. |
| `simulation m=127` | Simulated phi histogram for the `I+` mass channel. |
| `simulation m=131` | Simulated phi histogram for the `I+He` mass channel. |
| `simulation m=135` | Simulated phi histogram for the `I+He2` mass channel. |

The simulation-side angular recipe is:

```text
round(mass / u) == mass_select
b_ion_outside == true
phi = atan2(vy, vx) + pi
edges_phi = 0:0.05:2*pi
smooth = movmean(histogram_counts, 15)
y_plot = smooth / max(smooth)
```

This panel is about directionality and anisotropy. It is not a speed
distribution.

## Why this differs from `plot_experimental_comparison.py`

`plot_paper_v3.py` follows Strategy A (2-D projected speed vs raw VMI);
`plot_experimental_comparison.py` follows Strategy B (3-D speed vs
Abel-inverted reference). Peak positions and curve shapes therefore
differ by design, not by porting regression. See
`post_processing_strategy.md` §3-§4 for the full strategy comparison and
§5 Table for the per-script recipe row.

## Provenance and limits

The experimental curves are not re-derived directly from raw lab inputs in
Python. They are exported by the legacy MATLAB processing path into small CSV
reference files:

| CSV | Role |
|---|---|
| `data/reference/paper_v3/iplus_he_300mw_43563_radial.csv` | High-SNR I+He radial reference. |
| `data/reference/paper_v3/timescan_296_297_radial.csv` | Timescan radial reference. |
| `data/reference/paper_v3/iplus_he_300mw_43563_phi.csv` | I+He angular reference. |

The MATLAB exporter precedent is
`data/reference/scripts/export_paper_v3_reference_data.m`. It documents the
measurement IDs, centers, velocity factors, MATLAB functions, MAT-file inputs,
and external VMI toolbox requirements used to create those CSVs.

The Python port uses the final-state checkpoint vector from the selected
`RunDirectory`. Project-wide scope rules (Abel inversion, raw VMI image
interpretation, effusive branch, MATLAB multi-start matrix behavior all
out-of-scope) live in `CLAUDE.md` §"Current Scope".
