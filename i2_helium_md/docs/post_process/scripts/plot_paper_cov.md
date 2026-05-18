# `plot_paper_cov.py` — pair-covariance comparison

`scripts/post_processing/plot_paper_cov.py` is the Python port of
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_IplusHe_comparison_cov.m`.
It ports the **active non-effusive** branch only.

The MATLAB script extends the I+He comparison figure with two extra
diagnostics: a 2-D pair-angular covariance `(θ_a, θ_b)` and a 2-D
pair-radial-speed covariance `(v_r_a, v_r_b)`. The Python port splits
the result into six standalone figures.

## Outputs

Written to `<run>/figures/`:

| File | Content |
|---|---|
| `paper_cov_vmi_comparison.png` | exp + sim 2-D VMI map (side-by-side) |
| `paper_cov_radial_distribution.png` | 1-D velocity-distribution overlay (panel c) |
| `paper_cov_phi_distribution.png` | exp + sim 1-D phi(angle) distribution overlay (panel f) |
| `paper_cov_angular_pair_cov.png` | exp + sim angular pair covariance |
| `paper_cov_radial_pair_cov.png` | exp + sim radial pair-speed covariance |
| `paper_cov_pair_cov_traces.png` | 1-D axis-sum traces of the angular and radial pair covariance (sim vs exp), panels (g) and (h) |

The three cov-derived figures (`paper_cov_angular_pair_cov.png`,
`paper_cov_radial_pair_cov.png`, `paper_cov_pair_cov_traces.png`) are
produced only when the experimental covariance reference
(`data/reference/paper_cov/iplus_he_covariance.mat`) is present.

The phi distribution figure always renders. The experimental overlay
curve is loaded from the precomputed
`data/reference/paper_cov/iplus_he_phi.csv` (literal MATLAB output of
`mean(res_Iplus_He.image_polar(:, b_r), 2) / max(...)`, computed by the
exporter from the same three measurement IDs as the covariance
reference). If that CSV is missing, only the simulated curve is drawn.

> **Why a precomputed CSV instead of polar-binning the Cartesian image?**
> An earlier port tried to derive the experimental phi curve in Python
> by polar-binning the existing 2-D Cartesian VMI image
> (`data/reference/paper_v2/images/iplus_he_high_snr_vmi_image.mat`).
> That approach was wrong for two reasons: (1) the high-SNR `res_sum`
> source feeding the paper_v2 reference is a different averaged
> measurement than the three IDs `_cov.m` itself uses, and (2) Python's
> on-the-fly polar binning of a Cartesian grid cannot reproduce the
> closed-source VMI-toolbox polar transform (centring refinement,
> radial interpolation) that builds `res.image_polar`. Loading the
> precomputed CSV avoids both issues — the curve on disk is exactly
> what `_cov.m` plots.

## Panel (c) overlays

Panel (c) is the 1-D velocity-distribution comparison. The Python port
deliberately deviates from the MATLAB recipe here: the MATLAB script
overlays the I⁺ gas (43632) and averaged 300 mW I+He (45668/45662/45667)
curves; the Python script swaps both for **one cleaner experimental
reference** (the high-SNR `res_sum` I+He, also used by paper_v2) and
adds a **simulated** radial curve as a fine line. Gas is dropped as
uninformative for the I+He discussion.

The Python panel (c) plots four traces:

| Trace | Source | Style |
|---|---|---|
| I+He high-SNR | `data/reference/paper_v2/iplus_he_high_snr_radial.csv` via `load_paper_v2_radial_reference` | solid, lw 1.6, blue |
| exp v-cov trace | `radial_covariance_trace(cov_ref.cov_radial, cov_ref.velocity_centers_Aps)` | dashed, lw 1.2, black |
| sim radial | `paper_v2_velocity_curve(ion, mass_amu=131.0).normalised` vs `bin_centers_mps` | solid, lw 0.9, orange |
| sim v-cov trace | `radial_covariance_trace(sim_radial.counts, sim_radial.velocity_centers_Aps)` | red scatter dots, s 18 |

All four traces are max-normalised. The high-SNR reference is loaded
with the mass-correction already applied at export time
(`sqrt(127/131)`); the simulated curves use the package mass selection
on the Python ion checkpoint.

If the experimental covariance reference (`iplus_he_covariance.mat`) is
missing, the exp v-cov dashed trace is skipped but the rest of the
figure still renders. If the high-SNR reference CSV is missing, that
line is skipped with a warning.

The legacy `_cov.m` overlay (gas + averaged I+He) is documented here for
provenance. To reproduce it manually, plot the high-SNR file alongside
the gas/averaged-I+He CSVs produced by `plot_processed_VMI` — those
CSVs are not produced by the paper-cov exporter; they exist under
`data/reference/paper_v2/` only for the gas curve (43562) and the
power-scan I+He references.

The mass correction `sqrt(127/131)` maps the measured I⁺He velocity
back to a single-iodine equivalent so the gas and droplet curves can
share an axis (relevant for any direct comparison). `vf_single = 8.6178`
is the VMI toolbox detector calibration in metres-per-second per
pixel-per-time-of-flight (so `vf_single / 100` is Å/ps per pixel).

## Reference files written by the MATLAB exporter

`data/reference/scripts/export_paper_cov_reference_data.m` writes:

| File | Purpose |
|---|---|
| `iplus_he_covariance.mat` | `cov_angular`, `cov_radial`, `theta_centers_rad`, `velocity_centers_mps` (see below for processing) |
| `iplus_he_covariance.json` | provenance sidecar (lists the phi CSV under `companion_files`) |
| `iplus_he_phi.csv` | `phi_rad,signal_arb` — 1-D phi distribution from `mean(res_Iplus_He.image_polar(:, b_r), 2) / max(...)` |
| `iplus_he_cov_angular_preview.png` | quick-look preview (not loaded by Python) |
| `iplus_he_cov_radial_preview.png` | quick-look preview (not loaded by Python) |

The high-SNR I+He radial CSV (`iplus_he_high_snr_radial.csv`) lives
under `data/reference/paper_v2/` and is produced by the paper_v2
exporter (`export_paper_v2_reference_data.m`).

## Experimental pair covariance: how it is calculated

The experimental covariance matrices come from the external VMI MATLAB
toolbox function `generate_VMI_covariance_matrices`. The legacy script
calls it once with the three I+He droplet measurement IDs:

```matlab
center = autocenter_from_extended_data([45668, 45662, 45667]);
result = generate_VMI_covariance_matrices( ...
    [45668, 45662, 45667],     % measurement IDs
    [0, 600],                  % velocity range in detector pixels
    center,                    % auto-detected image centre
    [90, 90],                  % bins:  N_theta x N_velocity
    false,                     % apply_angular_filter
    true,                      % event_filter
    pi,                        % theta_target (back-scattering)
    40/180*pi);                % theta_range (40 degrees)
```

Although the implementation lives in a separate toolbox, the observable
outputs and conventional usage tell us what it is doing:

### Inputs

The function takes raw VMI detector data — every laser shot in the
listed measurement IDs contains zero, one, two, or more ion hits. Each
hit has a pixel position `(x, y)` and is converted to centred
velocity-space coordinates via the calibrated `(x − c_x, y − c_y) *
vf_single`. Only hits within the velocity-radius window
`[0, 600] px → [0, 5170 m/s ≈ 51.7 A/ps]` are kept.

### Per-shot pair binning

For each shot with two or more surviving hits, the function takes every
ordered pair `(hit_a, hit_b)` (with `a ≠ b`) and bins it into two
matrices:

- **Angular covariance** (`cov_angular`):
  `(θ_a, θ_b)` with `θ = atan2(y, x)`, binned into an `N_θ × N_θ`
  histogram on the range `[−π, π]`.
- **Radial covariance** (`cov_radial`):
  `(r_a, r_b)` with `r = sqrt(x² + y²)`, binned into an `N_r × N_r`
  histogram on the range `[0, 600] px`.

The `event_filter = true` flag restricts the sum to "valid" shots — in
the toolbox convention this means shots that survive the standard
detector quality filters (TOF gate, ROI, etc.). The
`apply_angular_filter = false` flag means no extra `theta_target` /
`theta_range` restriction is applied even though those parameters are
passed (they are kept for symmetry with other callsites).

The result is an **integer count matrix** whose `(i, j)` entry is the
total number of (a, b) hit-pair observations across all qualifying
shots that fell in bin `i` for `a` and bin `j` for `b`. This is **not**
a normalised correlation coefficient. It is essentially a 2-D joint
histogram of two single-particle observables, treated as a covariance
because every diagonal-removed off-diagonal entry encodes the joint
event rate that single-particle distributions cannot.

Why this is called "covariance": at fixed total number of pairs, the
2-D joint histogram is mathematically the **second-moment matrix** of
the per-shot single-particle histogram. After mean-subtraction (or, as
in this script, after diagonal removal — which absorbs the
single-particle marginals where the strongest self-correlation lives)
what remains is the inter-particle correlation structure.

### Diagonal removal

Both matrices have their diagonals zeroed:

```matlab
cov_angular = cov_angular - diag(diag(cov_angular));
cov_radial  = cov_radial  - diag(diag(cov_radial));
```

The diagonals are dominated by self-binning (the same single-particle
distribution overlaid against itself). Removing them isolates the
genuine pair-correlation structure.

### Smoothing (radial only)

The radial matrix is smoothed by a 2-sample moving mean along each
axis:

```matlab
cov_radial = movmean(cov_radial, 2, 1);
cov_radial = movmean(cov_radial, 2, 2);
```

The angular matrix is left unsmoothed.

### Axes on disk

- `theta_centers_rad` is `result.theta` in radians (range `[−π, π)`).
- `velocity_centers_mps` is `result.r * vf_single *
  sqrt(127/131)`. The `/100` from the legacy `result.r * vf_single /
  100 * mass_correction_factor` (which produces A/ps) cancels with the
  canonical on-disk choice of m/s, so the exporter just drops it. The
  Python loader divides by 100 to recover A/ps.

## Simulated pair covariance: how it is calculated

The simulated counterparts come from two Python helpers — both operate
on the **final-state** `IonCheckpoint` for the run directory and apply
the same atom selection used by paper_v4:

- mass `round(mass / U) == 131` (I⁺He),
- `b_ion_outside == True` for both atoms of a molecule (one atom
  failing the filter drops the whole pair),
- pair indexing `[0, N) ↔ [N, 2N)` (the molecule-pair convention used
  everywhere in this package).

### Angular pair covariance (`paper_v4_angular_pair_covariance`)

```python
theta_a = (arctan2(vx_a, vy_a) + pi) mod 2*pi
theta_b = (arctan2(vx_b, vy_b) + pi) mod 2*pi
edges   = linspace(0, 2*pi, n_bins + 1)
counts  = histogram2d(theta_a, theta_b, bins=(edges, edges))
```

Note the simulation wraps θ into `[0, 2π)` whereas the experimental
matrix uses `[−π, π]`. The driver `plot_paper_cov.py` rolls the sim
matrix by `N/2` along each axis (and shifts the labelled axis by `−π`)
so the side-by-side comparison shares an origin.

This helper does **not** zero the diagonal — for paper-v4 we wanted to
see the bin-matching events. The driver script for paper_cov is fine
with that because the diagonal in a 90-bin sim matrix from thousands of
pairs is not particularly hot.

### Radial pair-speed covariance (`radial_pair_speed_covariance`)

The new helper in `i2_helium_md/postprocess/paper_cov.py`:

```python
v_r_a   = sqrt(vx_a**2 + vy_a**2)            # A/ps
v_r_b   = sqrt(vx_b**2 + vy_b**2)            # A/ps
edges   = linspace(0, v_max_Aps, n_bins + 1) # default: 0..30, 90 bins
counts  = histogram2d(v_r_a, v_r_b, bins=(edges, edges))
counts[diag] = 0                              # remove_diagonal=True
counts  = movmean(counts, 2, axis=0)
counts  = movmean(counts, 2, axis=1)          # MATLAB-style 2x2 smoothing
```

Defaults match the legacy figure: 90 bins, diagonal zeroed, 2×2 moving
mean. Both can be turned off via keyword arguments for testing.

The package-level `moving_mean` in `i2_helium_md/postprocess/_smoothing.py`
implements MATLAB's `movmean(x, k)`: a centred sliding window with
shortened endpoint averages. Applying it along axis 0 then axis 1
matches `movmean(X, k, 1); movmean(X, k, 2)`.

## The 1-D `v-cov` trace on panel (c)

`_cov.m` line 419–429 builds a 1-D trace from the experimental
`cov_radial` and overlays it on the radial distribution:

```matlab
v_radial_corr = result.r * velocity_factor / 100 * mass_correction_factor;
b_v = v_radial_corr > vmin & v_radial_corr < vmax;       % vmin=4, vmax=22 A/ps
vd_radial_corr = sum(cov_radial(b_v, :), 1) / 2;          % divide by 2: symmetric
scatter(v_radial_corr, vd_radial_corr / max(vd_radial_corr));
```

Interpretation: pick the velocity rows where the *first* fragment lands
in the I+He fragmentation channel (roughly the bright lobes at 4–22
A/ps), then sum those rows along axis 0. The result is a 1-D
distribution of *partner* velocities given that the first fragment is
in the chosen band — equivalently, the marginal of the joint
distribution conditioned on the active radial band. Dividing by 2
accounts for the matrix being symmetric (each unordered pair is counted
twice).

The Python equivalent is
`i2_helium_md.postprocess.paper_cov.radial_covariance_trace`. It is
called twice in the driver:

- once on `cov_ref.cov_radial` (experimental, loaded from the MATLAB
  exporter) → black dashed line "exp v-cov trace",
- once on the simulated `radial_pair_speed_covariance.counts` (computed
  from the ion checkpoint) → red dots "sim v-cov trace".

Both are max-normalised before plotting so the shapes can be compared
on a common axis.

## Why the v-cov trace is small at low velocity

The `[4, 22] A/ps` band excludes the slow ion population near
`v ≲ 4 A/ps`. The trace is the conditional marginal: it is mechanically
zero below 4 A/ps because no row contributes to the sum at those bin
indices. This is the intended MATLAB behaviour and not a bug.

## Color clipping convention

Both 2-D covariance heatmaps use the legacy convention
`clim([0, 0.7 * max])`. The driver implements this via a
`matplotlib.colors.Normalize(vmin=0, vmax=0.7 * vmax)`. The VMI panel
(a + b) keeps the paper_v2 convention `[0, 0.8 * max]` so it matches
the existing I+He comparison figure.

## Conventions preserved verbatim from `_cov.m`

| Item | Value |
|---|---|
| Velocity factor `vf_single` | `8.6178` (m/s per pixel-per-ps) |
| Mass correction | `sqrt(127 / 131)` |
| `[4, 22] A/ps` band for the v-cov trace | hard-coded |
| Simulated VMI bins | `−35 : 0.5 : 35` A/ps (both axes) |
| Angular cov bins | 90 (matches `[90, 90]` MATLAB call) |
| Radial cov smoothing | `movmean(cov, 2, 1); movmean(cov, 2, 2)` |
| Diagonal removal | `cov - diag(diag(cov))` for both matrices |
| Color clipping | `[0, 0.7 * max]` for cov panels |
| Pair indexing | `[0, N) ↔ [N, 2N)` |

## Out of scope

- Effusive branch of `_cov.m` (out of project scope).
- Re-implementing `generate_VMI_covariance_matrices` in Python. The
  experimental matrices are treated as frozen MATLAB reference data,
  loaded from `.mat`. If the toolbox math is ever needed in Python, it
  would be a separate task and would need to reproduce the toolbox
  event filtering, calibration, and centring exactly.
- Abel inversion, pump-probe, full experimental VMI interpretation.

## Related files

- `i2_helium_md/postprocess/paper_cov.py` — Python helpers
  (`load_paper_cov_experimental_reference`,
  `radial_pair_speed_covariance`, `radial_covariance_trace`,
  `simulated_phi_distribution`,
  `covariance_axis_sum_normalised`).
- `i2_helium_md/postprocess/paper_v2.py` — reused for
  `load_paper_v2_phi_reference` (generic `phi_rad,signal_arb` CSV
  loader; consumed by `plot_paper_cov.py` to read the new
  `data/reference/paper_cov/iplus_he_phi.csv`).
- `i2_helium_md/postprocess/paper_v4.py` — reused for
  `paper_v4_angular_pair_covariance`.
- `i2_helium_md/postprocess/paper_v2.py` — reused for
  `paper_v2_velocity_map` (simulated 2-D VMI), `paper_v2_velocity_curve`
  (sim radial line), and `load_paper_v2_radial_reference` (high-SNR CSV).
- `data/reference/scripts/export_paper_cov_reference_data.m` — MATLAB
  exporter.
- `data/reference/paper_cov/README.md` — provenance for the experimental
  reference files.
- `tests/test_paper_cov.py`, `tests/test_plot_paper_cov_smoke.py` —
  unit and smoke tests.
