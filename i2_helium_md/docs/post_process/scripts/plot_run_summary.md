# `plot_run_summary.py` — consolidated post-processing summary

## What this script is

Per `CLAUDE.md` "Current Working Mode", this is the **preferred
single-entry-point diagnostic** for a finished run directory. It
consolidates every in-scope MATLAB post-processing figure into one
multi-page PDF plus per-panel PNGs.

The script has no CLI: configure it by editing the `USER SETTINGS`
block at the top, then run it from PyCharm's run button or a shell:

```bash
python scripts/post_processing/plot_run_summary.py
```

It is not a replacement for the focused scripts. Those scripts are
still the right tool when iterating on a single figure (faster cycle
time, more knobs exposed). `plot_run_summary.py` locks the panels to
the CLAUDE.md conventions and emits the bundle that documents a
finished run.

## Legacy MATLAB scripts consolidated here

```text
vmi_sim_3d_neutral_propa_HeDFT_mimic.m            → neutral energy balance
vmi_sim_3d_ion_propa.m                            → ion energy balance + temperature
simulation_image_only_trajectories.m              → HeDFT R(t), v(t)
post_process_single_pulse_paper_v3.m              → 1-D / 2-D polar VMI panels
post_process_single_pulse_paper_IplusHe_comparison_cov.m
                                                  → covariance-paper overlays
post_process_single_pulse_paper.m                 → bimodal Gaussian fit
post_process_single_pulse.m                       → 2-D (vx, vy) histogram
post_process_compare_radial_distributions.m       → time-resolved radial,
                                                    interatomic distance,
                                                    Boltzmann reference
compare_neutral_dynamics_to_HeDFT.m               → neutral cumtrapz r(t)
```

Out of scope (deferred per `CLAUDE.md`): Abel inversion, pump-probe,
effusive / gas-phase comparison, live-debug 3D animations.

## USER SETTINGS

The block at the top of the file controls one run and one set of
reference inputs. Three sentinel modes cover the in-scope workflows:

### 1. Experimental-condition droplet summary

```python
RUN_DIR                = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet_long"
HEDFT_REF_PATH         = None
VMI_REF_HE_PATH        = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_he.csv"
VMI_REF_GAS_PATH       = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_gas.csv"
VMI_REF_HE_HIGH_SNR_PATH = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_he_high_snr.csv"
PAPER_V2_REFERENCE_DIR   = PROJECT_ROOT / "data" / "reference" / "paper_v2"
PAPER_COV_REFERENCE_DIR  = PROJECT_ROOT / "data" / "reference" / "paper_cov"
```

This is the default checked into the file. HeDFT-trajectory panels are
skipped because experimental runs don't have a matching HeDFT
reference; every VMI panel is drawn.

### 2. 9 Å HeDFT comparison summary

```python
RUN_DIR        = PROJECT_ROOT / "data" / "runs" / "9A_hedft_comparison"
HEDFT_REF_PATH = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"
VMI_REF_HE_PATH         = None
VMI_REF_GAS_PATH        = None
VMI_REF_HE_HIGH_SNR_PATH = None
PAPER_V2_REFERENCE_DIR   = None
PAPER_COV_REFERENCE_DIR  = None
```

HeDFT R(t) and v(t) panels are drawn; VMI and paper-figure panels are
skipped. This matches `CLAUDE.md` "Keep these workflows separate."

### 3. Bare run, no references

```python
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "<run_name>"
HEDFT_REF_PATH         = None
VMI_REF_HE_PATH        = None
VMI_REF_GAS_PATH       = None
VMI_REF_HE_HIGH_SNR_PATH = None
PAPER_V2_REFERENCE_DIR   = None
PAPER_COV_REFERENCE_DIR  = None
```

Every reference-only panel is skipped silently (not an error). Useful
for new runs before reference data has been prepared.

### Other settings

| Name | Purpose |
|---|---|
| `OUT_DIR` | Output directory. `None` → `<RUN_DIR>/figures/`. |
| `SHOW_FIGURES` | Open figures interactively after writing. `False` is the right setting for PyCharm / headless / smoke tests. |
| `EXPERIMENTAL_NOISE_FLOOR` | Fraction of experimental panel max intensity below which pixels clip to background. Affects only experimental panels; simulated panels are unaffected. |
| `PAPER_V2_MASS_AMU` | Mass channel used for the paper-v2 simulated curves (default 131.0 amu = I⁺He). |

### Locked plot-tuning constants

These are *not* user settings — they are pinned to the CLAUDE.md
"Known Plotting Conventions":

```python
MASS_I               = 127.0     # bare I+
MASS_I_HE_AMU        = 131.0     # I+He
MASS_I_HE2_AMU       = 135.0     # I+He2
HIST_BIN_WIDTH_APS   = 0.04      # = 4 m/s on display
HIST_EDGE_MAX_APS    = 26.0      # = 2600 m/s on display
HIST_SMOOTHING_WINDOW = 15       # MATLAB movmean(..., 15)
VELOCITY_PLOT_V_MAX_MPS = 2800.0 # display range
```

To override any of these for a focused figure, use the matching
focused script (e.g. `plot_paper_v3.py`) rather than editing this
driver.

## Outputs

A single multi-page PDF plus one PNG per panel. Default layout under
`<RUN_DIR>/figures/`:

```text
run_summary.pdf
metadata.png
neutral_energy_balance.png
ion_energy_balance.png
ion_temperature_diagnostic.png
mass_spectrum.png
radial_velocity_with_vmi.png
paper_v2_vmi_comparison.png
paper_cov_radial_distribution.png
paper_cov_phi_distribution.png
paper_v2_polar_image_comparison.png
mass_resolved_velocities.png
radial_evolution_heatmap.png
interparticle_distance_histogram.png
paper_cov_angular_pair_cov.png
paper_cov_radial_pair_cov.png
paper_cov_pair_cov_traces.png
hedft_neutral_comparison.png            (only if HEDFT_REF_PATH set)
hedft_ion_comparison.png                (only if HEDFT_REF_PATH set)
boltzmann_overlay_initial.png           (only if cfg.json + ion checkpoint present)
```

## Panel guide

| Section label | Helper module(s) | What it shows |
|---|---|---|
| `metadata` | (cfg dump) | Run parameters from `cfg.json`, plus checkpoint sizes and the resolved reference paths. |
| `neutral_energy_balance` | `energy_balance.neutral_energy_totals` | Summed-over-atoms `E_kin / E_pot / E_dissip / E_system` traces from the neutral checkpoint. |
| `ion_energy_balance` | `energy_balance.ion_energy_totals` | Per-molecule energy traces for the ion stage, including `E_mass_attach_defect_eV`. |
| `ion_temperature_diagnostic` | `energy_balance` temperature path | The (T, 3) `temperature_diagnostic` field captured during ion propagation (schema v5). |
| `mass_spectrum` | `energy_balance.mass_spectrum` | Final-mass histogram on 1-amu bins, peaks centred on 127 / 131 / 135 amu. |
| `radial_velocity_with_vmi` | `velocity_distribution`, `_smoothing` | 4-curve overlay matching the legacy `simulation_image.m` figure (experimental I⁺He, experimental I⁺ gas, simulated 131 amu, simulated 135 amu). Display in m/s; 15-bin moving mean. |
| `paper_v2_vmi_comparison` | `paper_v2` loaders, `paper_v2_plotting` | Paper-v2 detector-plane VMI overlay. |
| `paper_cov_radial_distribution` | `paper_cov` loaders, `paper_cov_plotting` | Covariance-paper radial-distribution panel. |
| `paper_cov_phi_distribution` | `paper_cov` loaders | Covariance-paper azimuthal distribution. |
| `paper_v2_polar_image_comparison` | `paper_v2_plotting` | Paper-v2 polar (\|v\|, φ) image overlay. |
| `mass_resolved_velocities` | `velocity_distribution` | Separate final-velocity histograms per mass channel (127 / 131 / 135 amu), normalised. |
| `radial_evolution_heatmap` | `time_resolved.radial_distribution_evolution` | Per-atom `|r|` heat-map over uniformly-sub-sampled time (60 slices × 100 r-bins by default). |
| `interparticle_distance_histogram` | `pair_correlation.interparticle_distance_histogram` | Final I–I separation distribution. |
| `paper_cov_angular_pair_cov` | `pair_correlation.angular_pair_covariance`, `paper_cov_plotting` | 2-D `(θ_a, θ_b)` covariance heat-map (diagonal removed). |
| `paper_cov_radial_pair_cov` | `paper_cov.radial_pair_speed_covariance` | Radial pair-speed covariance against the covariance-paper reference. |
| `paper_cov_pair_cov_traces` | `paper_cov_plotting` | 1-D traces extracted from the covariance matrix. |
| `hedft_neutral_comparison` | `compare_trajectories.compare_neutral_to_hedft` | Neutral-stage `r(t)` overlap against the HeDFT trajectory (R only — no v before ionisation). |
| `hedft_ion_comparison` | `compare_trajectories.compare_distance / compare_velocity_magnitude` | Ion-stage distance + I1 / I2 velocity-magnitude overlap. |
| `boltzmann_overlay_initial` | `boltzmann_overlay.boltzmann_population` | Analytic `exp(-V/kT)` curve over `droplet_potential` overlaid on the initial-state radial histogram. |

Mass channels 127 / 131 / 135 amu, the 0.04 A/ps histogram bin width,
the 15-bin moving mean, and the m/s display axis up to 2800 m/s are
all pinned to `CLAUDE.md` "Known Plotting Conventions".

## Reference gating

Each section that needs optional reference data is wrapped in a
build-time guard. When the reference path is `None` (or the file is
missing), the section is skipped with a `[run_summary] skip <label>:
<reason>` print rather than aborting the whole figure. Concretely:

- `HEDFT_REF_PATH = None` → `hedft_neutral_comparison` and
  `hedft_ion_comparison` are skipped.
- `VMI_REF_HE_PATH = None` / `VMI_REF_GAS_PATH = None` →
  `radial_velocity_with_vmi` falls back to simulated-only curves or
  is skipped, depending on which VMI references are missing.
- `PAPER_V2_REFERENCE_DIR = None` → every `paper_v2_*` section is
  skipped.
- `PAPER_COV_REFERENCE_DIR = None` → every `paper_cov_*` section is
  skipped.
- Missing `cfg.json` → `boltzmann_overlay_initial` and the metadata
  page's cfg dump are skipped.

This keeps the HeDFT-comparison and experimental-VMI workflows
separable per `CLAUDE.md` "Post-Processing Workflow".

## How it differs from the focused scripts

| Focused script | Consolidated panel | When to still reach for the focused script |
|---|---|---|
| `plot_neutral_energy_balance.py` | `neutral_energy_balance` | Debugging an energy-bookkeeping discrepancy where you need separate `E_kin / E_pot / E_dissip` curves rather than the summed view. |
| `plot_ion_energy_balance.py` | `ion_energy_balance` | Same. |
| `plot_ion_temperature_diagnostic.py` | `ion_temperature_diagnostic` | When tuning the legacy MATLAB temperature recipe against `temperature_diagnostic` columns. |
| `plot_hedft_comparison.py` | `hedft_*_comparison` | When you want larger axes or to override the `~30` overlay-trace cap on the v(t) panel. |
| `plot_experimental_comparison.py` | `radial_velocity_with_vmi` | When iterating on the Strategy B (Abel-inverted 3-D) recipe with different mass tolerances. |
| `plot_paper_v2.py` | `paper_v2_*` | When iterating on paper-v2 specifically. |
| `plot_paper_v3.py` | `paper_cov_*` panels share the recipe | When iterating on paper-v3 specifically. |
| `plot_paper_v4.py` | `paper_cov_angular_pair_cov` reuses the pair recipe | When iterating on paper-v4 specifically. |
| `plot_paper_cov.py` | `paper_cov_*` | When iterating on paper-cov reference loading or the covariance traces. |

The focused scripts let you change settings the consolidated driver
locks in (mass tolerance, smoothing window, axis ranges, mass channel
of the headline curve). The consolidated driver is for "produce the
canonical bundle from a finished run."

## Provenance & scope

Out-of-scope per `CLAUDE.md`: Abel inversion, pump-probe, effusive /
gas-phase comparison, full experimental VMI image interpretation,
broad refactors. For the Strategy A vs Strategy B framing of the
detector-plane vs Abel-inverted comparison, see
`post_processing_strategy.md` §3–§4.

For the per-helper contract and units behind each panel, see:

- `compare_trajectories_module.md`
- `energy_balance_module.md`
- `hedft_loader_module.md`
- `velocity_distribution_module.md`
- `smaller_helper_modules.md` (`_smoothing`, `pair_correlation`,
  `polar_velocity`, `velocity_2d`, `time_resolved`,
  `boltzmann_overlay`)

For the paper-specific reference loaders and their CSV/MAT inputs,
see the matching script docs in this directory.
