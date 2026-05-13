# `plot_paper_v4.py` paper-v4 VMI comparison

This script ports the active non-effusive droplet branch of
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_v4.m`.
It is a focused paper-v4 reproduction, separate from the v3 paper figure and
from the consolidated run-summary diagnostics.

## Outputs

`plot_paper_v4.py` writes three figures under `<run>/figures/`:

| Output | Purpose |
|---|---|
| `compare_simulation_and_measurement_simpler.png` | One-panel radial VMI comparison. |
| `paper_v4_angular_pair_covariance.png` | Simulated I+He angular pair-covariance heatmap with scatter overlay. |
| `paper_v4_ion_mass_histogram.png` | Final ion mass diagnostic. |

## Experimental radial references

The experimental I+He curves are loaded from `data/reference/paper_v4/`. They
are created by `data/reference/scripts/export_paper_v4_reference_data.m`, which
uses the legacy MATLAB VMI processing path as the source of truth.

The exporter writes one radial CSV for each measurement declared in v4 lines
26-34:

| CSV | Curve label |
|---|---|
| `iplus_gas_160mw_43555_radial.csv` | `I+ gas 160 mW (43555)` |
| `iplus_drop_160mw_43554_radial.csv` | `I+ drop 160 mW (43554)` |
| `iplus_gas_600mw_43568_radial.csv` | `I+ gas 600 mW (43568)` |
| `iplus_drop_600mw_43567_radial.csv` | `I+ drop 600 mW (43567)` |
| `iplus_he_160mw_43556_radial.csv` | `I+He 160 mW (43556)` |
| `iplus_he_300mw_43563_radial.csv` | `I+He 300 mW (43563)` |

Each CSV has columns:

```text
v_mps,signal_arb
```

Python normalizes each plotted experimental signal by its maximum for plotting.

## Simulated radial curves

The simulation overlays match the active v4 MATLAB loop over masses `127` and
`131` amu:

```text
round(mass / u) == mass_select
b_ion_outside == true
v_projected = sqrt(vx^2 + vy^2)
edges_velocity = 0:0.05:35  A/ps
smooth = movmean(histogram_counts, 40)
shifted = smooth - min(smooth)
y_plot = shifted / max(shifted)
x_plot = centers_velocity * 100  m/s
```

As in the MATLAB paper script, this is a detector-plane projected-speed
comparison, not a full 3D speed distribution.

## Angular pair covariance

The separate covariance figure reproduces the simulated v4 block after the
radial comparison. It uses the `131 u` channel and keeps only molecule pairs
where both fragments pass the mass and outside-droplet selection.

The angular convention is literal v4:

```text
theta = atan2(vx, vy) + pi
numbins = 90
```

Unlike the older consolidated run-summary covariance helper, this focused v4
plot does not remove diagonal counts. It also overlays the selected
`(theta1, theta2)` points as scatter markers, matching the MATLAB figure.

## Limits

Only the active non-effusive branch is in scope. The effusive branch, raw VMI
interpretation, Abel inversion, and experimental covariance-matrix generation
remain deferred unless explicitly requested.
