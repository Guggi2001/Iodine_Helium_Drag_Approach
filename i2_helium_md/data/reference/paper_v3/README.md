# Paper V3 Reference CSVs

This directory is the target for
`data/reference/scripts/export_paper_v3_reference_data.m`.

Generated files are MATLAB-exported references for
`scripts/post_processing/plot_paper_v3.py`:

```text
iplus_he_high_snr_radial.csv
timescan_296_297_radial.csv
iplus_he_high_snr_phi.csv
```

Each CSV uses columns:

```text
v_mps,signal_arb          (or v_mps,signal_296,signal_297 for the timescan)
phi_rad,signal_arb        (for the phi file)
```

## Provenance note

The I+He files are named `iplus_he_high_snr_*`, **not**
`iplus_he_300mw_43563_*`, because the MATLAB exporter loads the
high-SNR co-added `res_sum` MAT file (path documented in the script
header) and overwrites the freshly computed standalone-43563 result
before export. The exported curves therefore correspond to the
co-added dataset, not the standalone 43563 measurement.

If you need the standalone 43563 radial, use
`data/reference/paper_v4/iplus_he_300mw_43563_radial.csv`. The two
files are not directly comparable.
