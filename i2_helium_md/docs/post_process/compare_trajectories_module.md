# The `compare_trajectories.py` module — a walkthrough

## What problem does this file solve?

Once the MD ion stage has produced an `ion.npz`, we need a numerical
answer to *how closely does my MD pipeline reproduce the published
HeDFT/TDDFT trajectory?* The legacy MATLAB block that does this lives in
`legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/
simulation_image_only_trajectories.m:97-118`:

```matlab
dx = data_ion.x_ci(1:Nmol,:) - data_ion.x_ci(1+Nmol:end,:);
dy = data_ion.y_ci(1:Nmol,:) - data_ion.y_ci(1+Nmol:end,:);
dz = data_ion.z_ci(1:Nmol,:) - data_ion.z_ci(1+Nmol:end,:);
dR_mean = mean(sqrt(dx.^2+dy.^2+dz.^2), 1);

tmax = min(max(tR), max(data_ion.time_i));
% mask both axes to [tmin, tmax]
dR_mean_on_tR = interp1(t_md, d_md, tR_use, 'linear');

good = isfinite(dR_mean_on_tR) & isfinite(R_use);
rmse  = sqrt(mean((dR_mean_on_tR(good) - R_use(good)).^2));
ratio = mean(dR_mean_on_tR(good) ./ R_use(good));
```

This module ports that contract to Python and applies the same
overlap-and-resample recipe to the I1/I2 velocity-magnitude curves
(which the MATLAB script plots but does not reduce to a number).

The full `simulation_image.m` script adds Abel inversion, mass filtering,
Bayesian histograms, and PDF export — all out of scope here per CLAUDE.md.

## Public API

```python
from i2_helium_md.simulation.run_directory import RunDirectory
from i2_helium_md.postprocess import (
    load_hedft_trajectory,
    compare_distance,
    compare_velocity_magnitude,
)

ion = RunDirectory("data/runs/single_pulse_N_2000").load_ion()
hedft = load_hedft_trajectory("data/reference/9A_All_Data.csv")

r = compare_distance(ion, hedft)
print(r.rmse, r.mean_ratio, r.overlap_t_min_ps, r.overlap_t_max_ps)

r1 = compare_velocity_magnitude(ion, hedft, atom="I1")
r2 = compare_velocity_magnitude(ion, hedft, atom="I2")
```

All three functions return the same `TrajectoryComparison` dataclass:

```python
@dataclass(frozen=True)
class TrajectoryComparison:
    quantity:           str        # "distance_A", "v1_magnitude_Aps", "v2_magnitude_Aps"
    overlap_t_min_ps:   float
    overlap_t_max_ps:   float
    num_overlap_points: int        # length of t_overlap_ps
    rmse:               float      # angstrom (distance) / angstrom-per-ps (velocity)
    mean_ratio:         float      # MD / HeDFT, dimensionless; nan if no nonzero refs
    t_overlap_ps:       np.ndarray
    md_on_hedft_grid:   np.ndarray
    hedft_on_overlap:   np.ndarray
```

## What gets compared

| Function | MD series | HeDFT series |
|---|---|---|
| `compare_distance` | mean over molecules of `sqrt(dx^2+dy^2+dz^2)` from `positions_x/y/z` | `distance_A` (the `R_distance` column) |
| `compare_velocity_magnitude(atom="I1")` | mean over molecules of `sqrt(vx^2+vy^2+vz^2)` for atoms `[0, N)` from `velocities_x/y/z` | `v1_magnitude_Aps` |
| `compare_velocity_magnitude(atom="I2")` | same for atoms `[N, 2N)` | `v2_magnitude_Aps` |

Particle indexing follows the rest of the codebase: `[0, num_molecules)`
are the I1 atoms, `[num_molecules, 2 * num_molecules)` are the I2 atoms.

## The overlap recipe

Both functions use the same private `_compare_series` helper:

1. Compute `t_min = max(ion.time_ps[0], hedft.time_ps[0])` and
   `t_max = min(ion.time_ps[-1], hedft.time_ps[-1])`. This is symmetric
   to the MATLAB `tmax = min(...)` and avoids extrapolation on either
   end.
2. Keep HeDFT samples inside `[t_min, t_max]` -> `t_overlap_ps`.
3. `np.interp(t_overlap_ps, ion.time_ps, md_series)` -> the MD series on
   the HeDFT grid (linear interpolation, matching MATLAB
   `interp1(..., 'linear')`).
4. RMSE is computed over samples where both interpolated MD and HeDFT
   are finite.
5. `mean_ratio` is computed over samples where both are finite *and*
   the HeDFT denominator is non-zero. The HeDFT velocities start at
   exactly 0, so without this guard the velocity ratio would be `inf`.
6. Empty overlap (`num_overlap_points < 2`) raises `ValueError`.

## Units

| Quantity | Units |
|---|---|
| `t_overlap_ps`, `overlap_t_min_ps`, `overlap_t_max_ps` | picoseconds |
| `rmse` for `quantity == "distance_A"` | angstrom |
| `rmse` for `quantity == "v1_magnitude_Aps"` / `v2_magnitude_Aps` | angstrom / picosecond |
| `mean_ratio` | dimensionless |

No conversion happens — both sides come in with these units already.

## Smoke check on the existing run

Against `data/runs/single_pulse_N_2000` (2000 molecules, 20 ps) and the
9 angstrom reference (14082 samples, 14.08 ps):

| Comparison | RMSE | mean_ratio | overlap |
|---|---|---|---|
| `compare_distance` | 26.19 A | 1.508 | [0.0, 14.08] ps |
| `compare_velocity_magnitude(atom="I1")` | 1.68 A/ps | 0.959 | [0.0, 14.08] ps |
| `compare_velocity_magnitude(atom="I2")` | 1.39 A/ps | 1.167 | [0.0, 14.08] ps |

The MD trajectory diverges from HeDFT at long times (MD ions reach
~112 A while HeDFT predicts ~84 A), so a mean distance ratio above 1
and a sizeable RMSE are expected; the velocity-magnitude ratios stay
within roughly 5-15% of unity over the same window.

## Out of scope

This module is numerical comparison only: no plotting, no CLI, no mass
filtering, no Bayesian histogramming. Plotting is done by
`scripts/post_processing/plot_hedft_comparison.py`. For project-wide
scope rules (Abel inversion, pump-probe, experimental VMI image
interpretation) see `CLAUDE.md` §"Current Scope".
