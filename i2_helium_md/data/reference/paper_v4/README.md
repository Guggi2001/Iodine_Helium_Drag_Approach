# Paper V4 Reference CSVs

This directory is the target for
`data/reference/scripts/export_paper_v4_reference_data.m`.

Generated files are small MATLAB-exported radial VMI references used by
`scripts/post_processing/plot_paper_v4_figure.py`. They are named by channel,
probe power, measurement ID, and observable, for example:

```text
iplus_he_300mw_43563_radial.csv
```

Each CSV has columns:

```text
v_mps,signal_arb
```
