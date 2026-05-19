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

## Provenance: I+He high-SNR (`iplus_he_high_snr_*`)

The high-SNR I+He radial CSV, phi CSV, and 2-D VMI image MAT all derive
from a single averaged `res_sum` produced by
`data/reference/scripts/generate_high_snr_iplus_he_mat.m` and saved to
`i2_helium_md/old_data/ressumI2HeNI^+He.mat`. The averaging recipe is:

- measurement IDs: `[45668, 45662, 45667, 45686]` (I+He droplet, low doping,
  300 mW, 03.12.24) — same set as `legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/simulation_image.m:153`,
  and a superset of the `paper_cov` triplet by the addition of 45686
- center: `[524.5297, 380.8430]` (hardcoded; matches `simulation_image.m:159`
  and the `paper_cov` phi pipeline)
- velocity factor: `vf_single = 8.6178`
- global background subtraction: OFF by default (the previous 17.10.24
  bg frame 43655 is from a different session; a 03.12.24 bg frame has
  not been identified)

To regenerate, open `generate_high_snr_iplus_he_mat.m` in MATLAB, run once
with `INSPECT_CENTERS_ONLY = true` to verify the per-file centers, then
flip the flag to `false` and rerun.

The I+He2 high-SNR exports (`iplus_he2_high_snr_*`) still derive from the
external VMI_matlab high-SNR directory until a corresponding generator is
written.

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
