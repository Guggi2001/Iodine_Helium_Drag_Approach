# Reference data layout

This directory holds normalized, Python-consumable copies of the
experimental VMI and HeDFT trajectory references used by the
post-processing scripts.

## Subfolders

```text
paper_v2/       <- post_process_single_pulse_paper_IplusHe_comparison.m
                   (export_paper_v2_reference_data.m)
paper_v3/       <- post_process_single_pulse_paper_v3.m
                   (export_paper_v3_reference_data.m)
paper_v4/       <- post_process_single_pulse_paper_v4.m
                   (export_paper_v4_reference_data.m)
vmi_summary/    <- export_vmi_reference_data.m (legacy Abel-inverted summary
                   consumed by plot_experimental_comparison.py and the
                   consolidated plot_run_summary.py)
scripts/        <- the MATLAB exporters that generate the above CSV/MAT files,
                   plus optional Python verification helpers.
```

Top-level files:

```text
9A_All_Data.csv   <- HeDFT trajectory reference (9 angstrom droplet)
18A_All_Data.csv  <- HeDFT trajectory reference (normalized 18 angstrom droplet)
```

## On-disk units convention

All reference CSVs use **m/s** as the canonical velocity unit. CSV
columns are named with explicit unit suffixes (`v_mps`, `vx_mps`,
`vy_mps`).

Python loaders convert m/s → Å/ps (divide by 100) at load time so
downstream plotting and binning code can stay in the documented Å/ps
conventions (CLAUDE.md "Known Plotting Conventions"). Until the MATLAB
exporters have been re-run, some on-disk files may still carry the
legacy `v_Aps` column; the loaders detect this and convert without
warning.

## Known issues

- **43556 image center disagreement.** `paper_v2/iplus_he_160mw_43556_radial.csv`
  uses VMI center `[524.5297, 380.8430]`; `paper_v4/iplus_he_160mw_43556_radial.csv`
  uses `[509.3664, 387.6409]`. Same raw measurement, two different radial
  profiles. We do not currently have a reliable VMI-centering procedure
  to choose between the two. Treat overlays of `paper_v2/43556` and
  `paper_v4/43556` as approximate until centering is verified.
- **`paper_v3` I+He curves are the high-SNR co-add, not 43563.** Files
  there are named `iplus_he_high_snr_{radial,phi}.csv` to reflect this.
  If you need the standalone 43563 radial, use
  `paper_v4/iplus_he_300mw_43563_radial.csv` (which is also not the same
  measurement as the high-SNR co-add).
- **`vmi_summary/` is a different reduction pipeline.** See
  `vmi_summary/README.md`. It runs Abel inversion, image smoothing, and a
  mass correction not applied in the paper-era exports, and it averages
  raw measurements that no paper-era export uses.

## Forbidden

Per CLAUDE.md, do not delete reference data or commit large checkpoints
into this directory. Keep reference CSVs small and inspectable.
