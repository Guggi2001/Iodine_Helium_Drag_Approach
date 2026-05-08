# `scripts/post_processing/plot_*` — legacy-debug walkthrough

Four scripts under `scripts/post_processing/` reproduce the legacy
MATLAB live-debug figures and the simulation-side panels of
`post_process_single_pulse_paper_v3.m`. All run post-hoc from a
finished `RunDirectory`; none of them touch the simulation modules.

| Script | Reproduces (legacy MATLAB) | Output |
|---|---|---|
| `plot_neutral_energy_balance.py` | `vmi_sim_3d_neutral_propa_HeDFT_mimic.m:965` | `<run>/figures/neutral_energy_balance.png` |
| `plot_ion_energy_balance.py` | `vmi_sim_3d_ion_propa.m:898` | `<run>/figures/ion_energy_balance.png` |
| `plot_ion_temperature_diagnostic.py` | `vmi_sim_3d_ion_propa.m:683` / `:883` | `<run>/figures/ion_temperature_diagnostic.png` |
| `plot_paper_figure.py` | `post_process_single_pulse_paper_v3.m` (minimal) | `<run>/figures/compare_simulation_and_measurement.{pdf,png}` |

## Common conventions

- Each script has a USER SETTINGS block at the top: edit `RUN_DIR` to
  point at a different run directory and rerun. There is no argparse
  front-end on purpose.
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

1×3 figure saved as `compare_simulation_and_measurement.pdf` (the
filename matches the legacy MATLAB `print()` target on line 391):

- **Tile A** -- experimental I+ gas / I+He droplet 1-D references
  overlaid with simulated mass-selected final-velocity histograms
  (I+He = 131 amu, I+He2 = 135 amu). Reuses
  `compute_final_velocity_histogram` and the smoothing helpers from
  `postprocess._smoothing` (15-bin moving mean, unit-max
  normalisation, `xlim=(0, 28) Å/ps`, `0:0.04:26` bin edges -- all
  matching the MATLAB).
- **Tile B** -- simulated azimuthal phi histogram from
  `phi_histogram(ion, mass_amu=131.0)`, smoothed and normalised to
  unit max, `xlim=(0, 2*pi)`.
- **Tile C** -- final ion mass spectrum from
  `mass_spectrum(ion, bin_width_amu=1.0)` plotted as a bar histogram.

Out of scope: the polar-VMI panels of the legacy figure (cos^2
angular anisotropy fit on `Tile 2`, `beta(v)` function, 3-D `surf`
of polar VMI image) require a 2-D polar VMI image not present in
`data/reference/`. They remain deferred per CLAUDE.md scope.

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
