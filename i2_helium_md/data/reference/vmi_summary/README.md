# VMI summary references (Abel-inverted 3-D speed distributions)

This folder holds the Abel-inverted, image-smoothed, mass-corrected
experimental VMI references consumed by
`scripts/post_processing/plot_experimental_comparison.py` and the
consolidated `scripts/post_processing/plot_run_summary.py` driver.

> **Convention fix (2026-07):** `signal_arb` is now pyabel's speed
> distribution `I(v)` from `distr.rIbeta()` — the properly v²- and
> sin θ-weighted **3-D speed distribution**, the like-for-like counterpart
> of the simulation's histogram of |v| (`compute_final_velocity_histogram`).
> The previous export shipped `radial_distribution` (the 2-D slice radial
> sum), which is one factor of v short of a speed distribution and gave the
> MD overlay a systematic ~v tilt; its low-v Abel center spike also
> dominated the gas curve. Regenerate old plots before comparing against
> archived figures.

These files are **not** directly comparable to the paper-era radial
references under `paper_v2/`, `paper_v3/`, or `paper_v4/`. The reduction
pipeline is documented below.

## Files

```text
vmi_iplus_he.csv           columns: v_mps,signal_arb
vmi_iplus_gas.csv          columns: v_mps,signal_arb
vmi_iplus_he_high_snr.csv  columns: v_mps,signal_arb
```

Velocity is stored in m/s on disk. Python loaders convert to Å/ps
internally so downstream plotting code keeps the documented Å/ps binning
conventions (CLAUDE.md "Known Plotting Conventions").

> **⚠️ The two I⁺He files are byte-identical (noted 2026-07-14).**
> `vmi_iplus_he.csv` and `vmi_iplus_he_high_snr.csv` are the same file: the
> "high-SNR" `res_sum` MAT is a pre-average of the same four measurements
> (45668/45662/45667/45686, ~23k counts total) that the I+He pipeline
> averages itself, so the two pipelines below produce identical output.
> Treat them as ONE reference, not two independent measurements. That single
> reference is a different-campaign (45xxx) measurement: its core agrees
> with the 17.10.24 `ihe_ked/` I⁺He reference (~63k counts), while its ⟨E⟩
> sits ~16% lower (condition-dependent tail weight).

## Pipeline (MATLAB)

Source script: `data/reference/scripts/export_vmi_reference_data.m`.

I+He channel:

1. Average four raw measurements: 45668, 45662, 45667, 45686.
2. Process each with `plot_processed_VMI(fn, true, [524.5297 380.8430], true)`.
3. Floor negative pixels to zero.
4. Apply `movmean(image, 3, 1)` then `movmean(image, 3, 2)` to smooth the
   2-D image.
5. Run `abel_invert_processed_VMI()` on the smoothed image.
6. Velocity axis: `speed_r * vf_single * sqrt(127/131)` with `vf_single = 8.6178`
   (`speed_r` is pyabel's radial grid converted to the same pixel units as `r`).
7. Signal: `speed_distribution` — pyabel `distr.rIbeta()` `I(v)`.

Gas channel:

1. Single raw measurement: 43632.
2. Process with `plot_processed_VMI(fn, true, [482.9299 392.4866], true)`.
3. Floor negative pixels to zero, then `movmean(image, 3, 1)` and
   `movmean(image, 3, 2)` — same smoothing as the I+He channels (added
   2026-07 together with the speed-distribution fix; the legacy script left
   the gas image raw).
4. Run `abel_invert_processed_VMI()`.
5. Velocity axis: `speed_r * vf_single` with `vf_single = 8.6178` (no mass correction).
6. Signal: `speed_distribution` — pyabel `distr.rIbeta()` `I(v)`.

High-SNR I+He channel:

1. Load the pre-averaged `res_sum` struct from the same MAT file used by
   `data/reference/scripts/export_paper_v2_reference_data.m`
   (`...\high_snr\ressumI2HeNI^+He[.mat]` on the legacy MATLAB share).
2. Floor negative pixels to zero.
3. Apply `movmean(image, 3, 1)` then `movmean(image, 3, 2)` to smooth the
   2-D image.
4. Run `abel_invert_processed_VMI()` on the smoothed image.
5. Velocity axis: `speed_r * vf_single * sqrt(127/131)` with `vf_single = 8.6178`.
6. Signal: `speed_distribution` — pyabel `distr.rIbeta()` `I(v)`.

The high-SNR file shares its source MAT with the paper-v2 high-SNR radial
export but applies a different pipeline: the paper-v2 export reads
`res.radial_distribution` raw (2-D projection), while this file is the
Abel-inverted, image-smoothed 3-D speed distribution `I(v)`.

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

## See also

The high-SNR `res_sum` struct also carries the full 2-D polar VMI image
`res.image_polar`. The polar image is exported separately by the
paper-v2 pipeline as
`data/reference/paper_v2/images/iplus_he_high_snr_vmi_polar_image.mat`
(fields `phi_rad`, `v_radius_mps`, `intensity_polar`) for the
`(phi, v)` side-by-side comparison panel in `plot_paper_v2.py`. It is
intentionally kept under `paper_v2/` rather than here because, like the
other paper-v2 outputs, it is the raw 2-D projection rather than an
Abel-inverted distribution.
