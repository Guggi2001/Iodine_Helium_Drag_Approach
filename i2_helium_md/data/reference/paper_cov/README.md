# paper_cov reference data

Experimental pair-covariance reference for the Python port of
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_IplusHe_comparison_cov.m`.

## Files

| File | Purpose |
|---|---|
| `iplus_he_covariance.mat` | Experimental angular + radial pair covariance, exported by the MATLAB script in `../scripts/export_paper_cov_reference_data.m`. Required by `scripts/post_processing/plot_paper_cov.py` for Fig 3 and Fig 4. |
| `iplus_he_covariance.json` | Sidecar metadata: source MATLAB function, measurement IDs, center, velocity factor, mass correction, axis equations, processing steps, export timestamp. |
| `iplus_he_cov_angular_preview.png` | Quick-look preview of the angular pair-covariance matrix. Not loaded by Python. |
| `iplus_he_cov_radial_preview.png` | Quick-look preview of the radial pair-covariance matrix. Not loaded by Python. |

## On-disk schema (`.mat`)

| Field | Shape | Units | Notes |
|---|---|---|---|
| `cov_angular` | `(N_theta, N_theta)` | counts | Diagonal already zeroed (no smoothing applied). |
| `cov_radial` | `(N_v, N_v)` | counts | Diagonal already zeroed AND 2 x 2 `movmean` applied along each axis. |
| `theta_centers_rad` | `(N_theta,)` | radian | Range `[-pi, pi]` (from `result.theta`). |
| `velocity_centers_mps` | `(N_v,)` | m/s | Mass-corrected by `mass_correction_factor` (`sqrt(127/131)`). The Python loader divides by 100 to recover A/ps. |

## Generating / regenerating

From MATLAB with the legacy VMI toolbox + single-pulse repository on the path:

```matlab
cd('data/reference/scripts');
export_paper_cov_reference_data;
```

This overwrites `iplus_he_covariance.mat` and `.json` in place. The two preview PNGs are also regenerated.

## Measurement provenance

| Item | Value |
|---|---|
| Measurement IDs | `[45668, 45662, 45667]` |
| Date | 2024-12-03 |
| Probe power | 300 mW |
| Sample | I+He droplet, low doping |
| Center calibration | `autocenter_from_extended_data([45668, 45662, 45667])` |

## Covariance call arguments (literal from `_cov.m:362`)

```
generate_VMI_covariance_matrices(
    [45668, 45662, 45667],   % measurement IDs
    [0, 600],                % velocity range in pixels
    center,                  % auto-detected
    [90, 90],                % theta x velocity bins
    false,                   % apply_angular_filter
    true,                    % event_filter
    pi,                      % theta_target (back-scattering)
    40/180*pi                % theta_range (40 degrees)
);
```

## See also

- `i2_helium_md/postprocess/paper_cov.py` — Python loader and simulated counterparts.
- `i2_helium_md/scripts/post_processing/plot_paper_cov.py` — driver script.
- `legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_IplusHe_comparison_cov.m` — original MATLAB source.
