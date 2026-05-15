# Paper V4 Reference CSVs

This directory is the target for
`data/reference/scripts/export_paper_v4_reference_data.m`.

Generated files are small MATLAB-exported radial VMI references used by
`scripts/post_processing/plot_paper_v4.py`. They are named by channel,
probe power, measurement ID, and observable, for example:

```text
iplus_he_300mw_43563_radial.csv
```

Each CSV has columns:

```text
v_mps,signal_arb
```

## Known issues

- **43556 image center disagrees with `paper_v2`.** This folder's
  `iplus_he_160mw_43556_radial.csv` is exported with VMI image center
  `[509.3664, 387.6409]`. The `paper_v2/iplus_he_160mw_43556_radial.csv`
  for the same raw measurement uses center `[524.5297, 380.8430]`. The
  centers shift the radial-pixel ring assignment, so the two radial profiles
  differ in shape even though the underlying raw frames are identical. We
  do not currently have a reliable VMI-centering procedure to determine
  which value is correct. Treat overlays of `paper_v2/43556` and
  `paper_v4/43556` as approximate until centering is verified.
- **43563 center is the same `[509.3664, 387.6409]` used here.** The
  `paper_v3` folder's `iplus_he_high_snr_radial.csv` (renamed from the
  misleadingly-named `iplus_he_300mw_43563_radial.csv`) is the high-SNR
  co-add, not the standalone 43563 measurement; comparing it to this
  folder's 43563 export is not apples-to-apples.
