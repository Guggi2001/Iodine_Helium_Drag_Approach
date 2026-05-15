# VMI summary references (Abel-inverted)

This folder holds the Abel-inverted, image-smoothed, mass-corrected
experimental VMI references consumed by
`scripts/post_processing/plot_experimental_comparison.py` and the
consolidated `scripts/post_processing/plot_run_summary.py` driver.

These files are **not** directly comparable to the paper-era radial
references under `paper_v2/`, `paper_v3/`, or `paper_v4/`. The reduction
pipeline is documented below.

## Files

```text
vmi_iplus_he.csv   columns: v_mps,signal_arb
vmi_iplus_gas.csv  columns: v_mps,signal_arb
```

Velocity is stored in m/s on disk. Python loaders convert to Å/ps
internally so downstream plotting code keeps the documented Å/ps binning
conventions (CLAUDE.md "Known Plotting Conventions").

## Pipeline (MATLAB)

Source script: `data/reference/scripts/export_vmi_reference_data.m`.

I+He channel:

1. Average four raw measurements: 45668, 45662, 45667, 45686.
2. Process each with `plot_processed_VMI(fn, true, [524.5297 380.8430], true)`.
3. Floor negative pixels to zero.
4. Apply `movmean(image, 3, 1)` then `movmean(image, 3, 2)` to smooth the
   2-D image.
5. Run `abel_invert_processed_VMI()` on the smoothed image.
6. Velocity axis: `r * vf_single * sqrt(127/131)` with `vf_single = 8.6178`.
7. Radial signal: `movmean(radial_distribution, 1)` (window size 1, effectively
   no-op; preserved for parity with the legacy script).

Gas channel:

1. Single raw measurement: 43632.
2. Process with `plot_processed_VMI(fn, true, [482.9299 392.4866], true)`.
3. Run `abel_invert_processed_VMI()`.
4. Velocity axis: `r * vf_single` with `vf_single = 8.6178` (no mass correction).
5. Radial signal: `radial_distribution` directly.

## Why this is different from paper_v2/v3/v4

- **Different raw measurements.** The paper-era exports use 43554, 43555,
  43556, 43562, 43563, 43567, 43568, 43569 (depending on script). This
  folder uses 45668/45662/45667/45686 (I+He average) and 43632 (gas).
- **Abel inversion.** Paper exports plot the raw VMI radial distribution
  directly; this folder inverts to recover the 3-D velocity distribution.
- **Image smoothing.** Paper exports do not smooth the 2-D image before
  radial extraction.
- **Mass correction.** Paper exports do not apply the `sqrt(127/131)`
  isotope correction to the I+He radial axis.

If you want to compare these curves against any paper-era export, you
must keep the pipeline differences in mind. They are not the same
observable computed at a different time; they are a different physical
quantity (Abel-inverted 3-D distribution vs raw 2-D projection).
