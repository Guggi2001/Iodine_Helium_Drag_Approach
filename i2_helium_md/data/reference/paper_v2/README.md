# Paper v2 reference exports

This folder is for MATLAB-exported references used by
`scripts/post_processing/plot_paper_v2.py`. The Python name `paper_v2` refers
to the legacy MATLAB script
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_IplusHe_comparison.m`.

Expected radial CSV files use columns `v_mps,signal_arb` (after the on-disk
unit normalisation; the Python loader converts to Å/ps internally). The first
two are the curves actually plotted by the MATLAB paper-v2 radial panel; the
final two are optional I+He power-scan comparison curves:

- `iplus_gas_300mw_43562_radial.csv`
- `iplus_he_high_snr_radial.csv`
- `iplus_he_160mw_43556_radial.csv`
- `iplus_he_600mw_43569_radial.csv`

The optional phi comparison CSV uses columns `phi_rad,signal_arb`:

- `iplus_he_high_snr_phi.csv`

Processed 2-D VMI image references live under `images/` as MATLAB `.mat`
matrix data with a JSON sidecar. Python also accepts `.npz` files with the
same fields for manually converted references. The expected fields are
`vx_mps`, `vy_mps`, and `intensity`. For MATLAB exports, `vx_mps` and `vy_mps`
are full 2-D coordinate grids matching `intensity`, normalized for Matplotlib
`pcolormesh(X, Y, C)`: `vx_mps` is the plot x-grid and `vy_mps` is the plot
y-grid. These image references are not raw detector data; they are the
processed MATLAB/VMI-toolbox result that the legacy script plotted.

## Known issues

- **43556 image center disagrees with `paper_v4`.** This folder's
  `iplus_he_160mw_43556_radial.csv` is exported with VMI image center
  `[524.5297, 380.8430]`. The `paper_v4/iplus_he_160mw_43556_radial.csv`
  for the same raw measurement uses center `[509.3664, 387.6409]`. The
  centers shift the radial-pixel ring assignment, so the two radial profiles
  differ in shape even though the underlying raw frames are identical. We
  do not currently have a reliable VMI-centering procedure to determine
  which value is correct. Treat overlays of `paper_v2/43556` and
  `paper_v4/43556` as approximate until centering is verified.
