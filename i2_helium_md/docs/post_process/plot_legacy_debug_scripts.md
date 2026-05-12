# `scripts/post_processing/plot_*` — legacy-debug walkthrough

Five scripts under `scripts/post_processing/` reproduce the legacy
MATLAB live-debug figures and the simulation-side panels of
`post_process_single_pulse_paper_v3.m` / `v4.m`. All run post-hoc
from a finished `RunDirectory`; none of them touch the simulation modules.

| Script | Reproduces (legacy MATLAB) | Output |
|---|---|---|
| `plot_neutral_energy_balance.py` | `vmi_sim_3d_neutral_propa_HeDFT_mimic.m:965` | `<run>/figures/neutral_energy_balance.png` |
| `plot_ion_energy_balance.py` | `vmi_sim_3d_ion_propa.m:898` | `<run>/figures/ion_energy_balance.png` |
| `plot_ion_temperature_diagnostic.py` | `vmi_sim_3d_ion_propa.m:683` / `:883` | `<run>/figures/ion_temperature_diagnostic.png` |
| `plot_paper_figure.py` | `post_process_single_pulse_paper_v3.m` active droplet branch | `<run>/figures/compare_simulation_and_measurement.{pdf,png}`, `<run>/figures/ion_mass_histogram.{pdf,png}` |
| `plot_paper_v4_figure.py` | `post_process_single_pulse_paper_v4.m` active droplet branch | `<run>/figures/compare_simulation_and_measurement_simpler.{pdf,png}`, `<run>/figures/paper_v4_angular_pair_covariance.{pdf,png}`, `<run>/figures/paper_v4_ion_mass_histogram.{pdf,png}` |

## Common conventions

- Each script has a USER SETTINGS block at the top: edit `RUN_DIR` to
  point at a different run directory and rerun. `plot_paper_figure.py`
  and `plot_paper_v4_figure.py` also accept CLI path overrides and
  `--no-show` for headless checks.
- Output is interactive (`plt.show()`); a copy is also written under
  `<run>/figures/`. The directory is created on demand.
- All scripts load via `RunDirectory(...).load_neutral()` /
  `.load_ion()`. The ion scripts therefore inherit the schema-v5
  loader's "regenerate older runs" rule.

## `plot_neutral_energy_balance.py`

Single-axes figure with four traces: sums of `E_kin_eV`, `E_pot_eV`,
`E_dissip_eV` over atoms, plus `E_system = E_kin + E_pot + E_dissip`.
Title `Energy balance neutral atoms` matches the MATLAB.

Sanity check: `E_system` should be ~flat (Verlet drift) and
`E_dissip` monotone non-decreasing. `E_kin + E_pot` falls with
`E_dissip` rising in lockstep.

## `plot_ion_energy_balance.py`

Single-axes figure with five traces: per-molecule sums of `E_kin`,
`E_pot`, `E_dissip`, `E_mass_attach_defect`, total `E_system =
E_kin + E_pot + E_dissip + E_mass_attach_defect`. Title `energy
balance ions`.

The `E_mass_attach_defect` curve is typically negative (or small
positive) because helium-attaches at non-zero velocity overstate
the true post-attachment kinetic energy.

## `plot_ion_temperature_diagnostic.py`

Twin-y-axis figure. Left axis: actual `<T'/T>` and mass-ratio
`<T'/T>` from `temperature_diagnostic[:, 0:2]`. Right axis:
`<theta_lab>` in degrees from `temperature_diagnostic[:, 2] * 180/pi`.
NaN rows (no collision in that stored step) are masked before
plotting.

The third column is the **lab-frame** scattering angle, not the
COM-frame one. For heavy projectile on light target the lab cone is
very narrow.

Sanity check for I+ on He (`rho = 127/4 ≈ 31.75`):
`<T'/T>_from_mass_ratio` should cluster near
`(1 + rho^2) / (1 + rho)^2 ≈ 0.943`. `<theta_lab>` should sit between
~0 and `asin(1/rho) ≈ 1.81 deg`, matching the legacy MATLAB plot's
y-axis range of roughly 0.4 to 1.8.

If the loaded checkpoint has no valid (non-NaN) rows the script
prints a warning and exits with code 1; this happens on pre-v5
checkpoints or runs with no collisions.

## `plot_paper_figure.py`

Two-panel figure saved as `compare_simulation_and_measurement.pdf`,
matching the active non-effusive branch of
`post_process_single_pulse_paper_v3.m`. The top panel overlays v3
experimental radial VMI references with simulated mass-selected
detector-plane projected-velocity curves for 127, 131, and 135 amu.
The bottom panel overlays the v3 experimental phi reference with
simulation phi histograms for the same mass channels.

The final ion mass histogram from MATLAB line 397 is preserved as a
separate `ion_mass_histogram.{pdf,png}` output because the legacy script
opens it after exporting the main paper figure.

Detailed curve-by-curve explanation:
[`docs/post_process/scripts/plot_paper_figure_paper_v3.md`](scripts/plot_paper_figure_paper_v3.md).

## `plot_paper_v4_figure.py`

One-panel radial velocity comparison saved as
`compare_simulation_and_measurement_simpler.pdf`, plus separate
`paper_v4_angular_pair_covariance.{pdf,png}` and
`paper_v4_ion_mass_histogram.{pdf,png}` outputs. The radial panel loads all
available v4 experimental radial CSVs from `data/reference/paper_v4/` and
overlays simulated projected-velocity curves for masses 127 and 131 amu.

Detailed curve-by-curve explanation:
[`docs/post_process/scripts/plot_paper_v4_figure.md`](scripts/plot_paper_v4_figure.md).

## Tests

`tests/test_plot_legacy_debug_smoke.py` imports each script as a
module, monkeypatches `plt.show` and `plt.savefig` (so figures are
not written during test runs), invokes `main()`, and asserts the
expected figure count. The whole module is gated on the
realistic-experimental run + VMI references being present, so a CI
without those artifacts skips cleanly.

The temperature-diagnostic test additionally requires a v5 ion
checkpoint with at least one non-NaN row; otherwise it skips with a
clear message.
