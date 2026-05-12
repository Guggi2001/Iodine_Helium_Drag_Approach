# `plot_paper_figure.py` paper-v3 VMI comparison

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

`plot_paper_figure.py` writes two groups of files under
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
| `I2:I+He TS` | Experimental timescan radial VMI reference from MATLAB `mean_timescan_2d_VMI([296:297], false, [524.5297 380.8430], false)`. |
| `I2:I+He` | Experimental high-SNR processed I+He radial VMI reference from the MATLAB `res_sum` path. |
| `simulated v.distr. m=127` | Simulated bare iodine ion channel, approximately `I+`. |
| `simulated v.distr. m=131` | Simulated `I+He` channel. |
| `simulated v.distr. m=135` | Simulated `I+He2` channel. |

The simulation recipe follows MATLAB v3 literally:

```text
round(mass / u) == mass_select
b_ion_outside == true
v_projected = sqrt(vx^2 + vy^2)
edges_velocity = 0:0.05:35  A/ps
x_plot = centers_velocity * 100  m/s
smooth = movmean(histogram_counts, 20)
y_plot = smooth / max(smooth)
```

The important point is the projection:

```text
sqrt(vx^2 + vy^2)
```

This is the detector-plane projected speed. It intentionally ignores `vz`,
because the paper-v3 figure is comparing to a VMI detector projection.

## Bottom panel: phi distribution

The bottom panel asks:

> Do the simulated ions leave in the same azimuthal directions as the
> processed experimental VMI image?

Curves:

| Curve | Meaning |
|---|---|
| `I2:I+He` | Experimental angular distribution from the processed polar VMI image. MATLAB averages `res.image_polar` over the selected radial range and normalizes the result. |
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

`plot_paper_figure.py` and `plot_experimental_comparison.py` both compare
simulation to experimental VMI-derived references, but they do not compare the
same observable.

| Detail | `plot_experimental_comparison.py` | `plot_paper_figure.py` |
|---|---|---|
| MATLAB source | `simulation_image.m`-style comparison | `post_process_single_pulse_paper_v3.m` active droplet branch |
| Experimental CSVs | `data/reference/vmi_iplus_*.csv` | `data/reference/paper_v3_*.csv` |
| Velocity units | `A/ps` | `m/s` |
| Simulation speed | Full 3D `sqrt(vx^2 + vy^2 + vz^2)` | Detector-plane `sqrt(vx^2 + vy^2)` |
| Mass channels | `131`, `135` | `127`, `131`, `135` |
| Angular panel | No | Yes |

Different-looking curves are therefore expected. The 3D speed histogram
includes motion along `vz`, while the paper-v3 plot compares the detector-plane
projection. The reference data also comes from different MATLAB processing
recipes.

Use `plot_paper_figure.py` when reproducing the paper-v3 MATLAB figure. Use
`plot_experimental_comparison.py` when reproducing the older 1-D experimental
velocity overlay from the `simulation_image.m` workflow.

## Provenance and limits

The experimental curves are not re-derived directly from raw lab inputs in
Python. They are exported by the legacy MATLAB processing path into small CSV
reference files:

| CSV | Role |
|---|---|
| `data/reference/paper_v3_iplus_he_radial.csv` | High-SNR I+He radial reference. |
| `data/reference/paper_v3_timescan_radial.csv` | Timescan radial reference. |
| `data/reference/paper_v3_iplus_he_phi.csv` | I+He angular reference. |

The MATLAB exporter precedent is
`data/reference/scripts/export_paper_v3_reference_data.m`. It documents the
measurement IDs, centers, velocity factors, MATLAB functions, MAT-file inputs,
and external VMI toolbox requirements used to create those CSVs.

Out of scope for this Python port:

- Abel inversion,
- raw experimental VMI image interpretation,
- full experimental image-processing expansion,
- the `effusive_dynamics` branch,
- MATLAB's full multi-start `vx_total(:, start_id)` matrix behavior.

The Python port uses the final-state checkpoint vector from the selected
`RunDirectory` and reproduces the plotted MATLAB recipe first.
