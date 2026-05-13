# Paper v2 reference exports

This folder is for MATLAB-exported references used by
`scripts/post_processing/plot_paper_v2.py`. The Python name `paper_v2` refers
to the legacy MATLAB script
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_IplusHe_comparison.m`.

Expected radial CSV files use columns `v_Aps,signal_arb`. The first two are
the curves actually plotted by the MATLAB paper-v2 radial panel; the final two
are optional I+He power-scan comparison curves:

- `iplus_gas_300mw_43562_radial.csv`
- `iplus_he_high_snr_radial.csv`
- `iplus_he_160mw_43556_radial.csv`
- `iplus_he_600mw_43569_radial.csv`

The optional phi comparison CSV uses columns `phi_rad,signal_arb`:

- `iplus_he_high_snr_phi.csv`

Processed 2-D VMI image references live under `images/` as MATLAB `.mat`
matrix data with a JSON sidecar. Python also accepts `.npz` files with the same
fields for manually converted references. The expected fields are `vx_Aps`,
`vy_Aps`, and `intensity`. For MATLAB exports, `vx_Aps` and `vy_Aps` are full
2-D coordinate grids matching `intensity`, normalized for Matplotlib
`pcolormesh(X, Y, C)`: `vx_Aps` is the plot x-grid and `vy_Aps` is the plot
y-grid. These image references are not raw detector data; they are the
processed MATLAB/VMI-toolbox result that the legacy script plotted.
