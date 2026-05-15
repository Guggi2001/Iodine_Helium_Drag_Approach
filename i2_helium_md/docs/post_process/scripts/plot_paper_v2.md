# `plot_paper_v2.py` I+He comparison and 2-D VMI images

`scripts/post_processing/plot_paper_v2.py` is the Python consistency name for
the legacy MATLAB script
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_IplusHe_comparison.m`.
It ports the active non-effusive I+He comparison branch, not a historical
script literally named `paper_v2`.

## Outputs

`plot_paper_v2.py` writes two figures under `<run>/figures/`:

- `paper_v2_compare_simulation_and_measurement.png`
- `paper_v2_phi_comparison.png`

The main figure follows MATLAB's 2x2 layout:

- top-left: processed experimental 2-D VMI image exported from MATLAB
- top-right: simulated `vx/vy` detector-plane velocity map
- bottom row: radial VMI comparison

The phi comparison remains a separate figure because the MATLAB script opens a
new figure for the angular distribution.

## Experimental references

Radial experimental references live in `data/reference/paper_v2/` as CSV files
with columns `v_mps,signal_arb` (velocity in m/s; the Python loader converts
to A/ps internally). Legacy files with `v_Aps,signal_arb` are also accepted
during the transition. The first two are the curves actually plotted in
MATLAB's radial panel:

- `iplus_gas_300mw_43562_radial.csv`: `res_Iplus_gas`, measurement 43562,
  center `[482.9299 392.4866]`.
- `iplus_he_high_snr_radial.csv`: `res_Iplus_He` after the MATLAB script
  overwrites the measurement load with high-SNR `res_sum`.

The exporter also writes optional I+He power-scan comparison curves:

- `iplus_he_160mw_43556_radial.csv`
- `iplus_he_600mw_43569_radial.csv`

The optional high-SNR phi reference is `iplus_he_high_snr_phi.csv` with columns
`phi_rad,signal_arb`. When present, it is overlaid in the separate phi
comparison figure; when absent, the script still writes the simulated phi
diagnostic and prints a skip message.

The processed 2-D VMI image reference lives under
`data/reference/paper_v2/images/` as MATLAB `.mat` plus a JSON sidecar. Python
also accepts `.npz` files with the same fields for manually converted
references. The image fields are:

- `vx_mps`
- `vy_mps`
- `intensity`

For MATLAB exports, `vx_mps` and `vy_mps` are full 2-D coordinate grids
matching `intensity`, normalized for Matplotlib's `pcolormesh(X, Y, C)`. The
Python loader divides them by 100 at load time so the resulting axes
(`PaperV2VMIImageReference.vx_Aps`, `.vy_Aps`) are in A/ps. Legacy files using
`vx_Aps,vy_Aps,intensity` axes are also accepted unchanged.

This is intentionally not raw detector data. The MATLAB exporter writes the
processed image that the legacy script plots after the VMI toolbox and
high-SNR MAT-file processing have already happened.

Exporter:
`data/reference/scripts/export_paper_v2_reference_data.m`.

## Simulation recipe

The simulated image panel is a literal port of the MATLAB block:

- mass selection: `round(mass/u) == 131`
- outside-droplet filter: `b_ion_outside`
- detector-plane components: `vx_total(:,1)` and `vy_total(:,1)`
- image bins: `-35:0.2:35` A/ps
- counting rule: nearest bin center via MATLAB-style `get_closest_index`
- plotting orientation: the helper stores `counts[vx_index, vy_index]`, then
  the Python figure transposes for Matplotlib so `v_x` is horizontal and `v_y`
  is vertical.

The radial simulation curve uses projected detector speed
`sqrt(vx^2 + vy^2)`, velocity edges `0:0.05:35` A/ps, moving mean window `10`,
and max normalization. The phi curve uses `atan2(vy, vx) + pi`, edges
`0:0.05:2*pi`, moving mean window `15`, and max normalization.

## Limits

Python does not perform raw VMI extraction, Abel inversion, or full image
interpretation. The experimental image is a processed MATLAB reference artifact
used for comparison. Experimental covariance images and effusive/gas-phase
branches remain out of scope unless explicitly requested.
