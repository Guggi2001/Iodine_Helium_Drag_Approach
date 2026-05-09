# Consolidated MATLAB Post-Processing Port Plan

**Status:** implemented. All in-scope helper modules, the consolidated
`plot_run_summary.py` driver, and focused tests are merged. Full pytest
suite (422 tests) passes. End-to-end PDFs are produced from
`data/runs/9A_hedft_comparison` (HeDFT-comparison run, 17 sections) and
`data/runs/single_pulse_droplet` (experimental-conditions run, 15
sections). Visual comparison vs. the legacy MATLAB figures is the only
remaining verification step and is left to human review of
`data/runs/*/figures/run_summary.pdf`.

## Context

A subset of the legacy MATLAB single-pulse post-processing has already been
ported into Python (energy balance, temperature, HeDFT trajectory comparison,
1D radial / phi velocity histograms, mass spectrum, mass-resolved velocity
histograms with experimental VMI overlay). The legacy MATLAB tree contains
~30 post-processing scripts.

The goal of this plan is to survey all of those legacy scripts and design a
single consolidated Python script that contains every *useful* post-processing
operation, so that one entry point produces every legacy diagnostic from a
finished run directory.

Scope is bounded by `CLAUDE.md`:

- **In:** post-processing of finished `RunDirectory` outputs against reference
  data already in `data/reference/`.
- **Out:** Abel inversion, full experimental VMI interpretation, pump-probe,
  effusive dynamics, new physics, broad refactors.
- Different simulation runs (HeDFT-comparison vs. experimental-conditions vs.
  droplet-distribution) must stay separable.

## MATLAB Inventory and Verdict

Status legend: **DONE** = already ported · **PORT** = port useful pieces ·
**SKIP** = redundant, debug-only, or out-of-scope per `CLAUDE.md`.

### Already ported (do not re-do)

| MATLAB | Python | Status |
| ------ | ------ | ------ |
| `vmi_sim_3d_neutral_propa_HeDFT_mimic.m` (energy panel) | `plot_neutral_energy_balance.py` | DONE |
| `vmi_sim_3d_ion_propa.m` (energy + temperature panels) | `plot_ion_energy_balance.py`, `plot_ion_temperature_diagnostic.py` | DONE |
| `single_pulse_simulation/HeDFT_comparison/simulation_image_only_trajectories.m` | `plot_hedft_comparison.py` + `compare_trajectories.py` | DONE |
| `single_pulse_simulation/HeDFT_comparison/simulation_image.m` (velocity-overlay panel) | `plot_experimental_comparison.py` + `velocity_distribution.py` | DONE |
| `single_pulse_simulation/post_process_single_pulse_paper_v3.m` (radial v + phi + mass spectrum) | `plot_paper_figure.py` | DONE (1D panels only) |

### Useful — port unique operations into the consolidated script

| MATLAB | Unique post-processing operation | Notes |
| ------ | -------------------------------- | ----- |
| `post_process_single_pulse_paper_v3.m` (deferred panels) | 2D polar VMI histogram (v, phi); cos² anisotropy fit; β(v); 3D polar surface | Doable from existing 3D velocities — **no Abel inversion needed**; we have full velocity vectors, not a projected experimental image. Was deferred only because the polar 2D *experimental* image is missing; the simulation side does not require Abel. |
| `post_process_single_pulse_paper_v4.m` | Angular pair covariance matrix θ₁ × θ₂ | Pair correlation of I-I ion-pair angles. |
| `post_process_single_pulse_paper.m` (older) | Bimodal Gaussian fit to radial velocity distribution | Two-population decomposition (high-v tail vs. main peak). |
| `post_process_single_pulse.m` | 2D velocity density (vₓ, vᵧ) | Pure rebin of existing 3D velocities. |
| `post_process_compare_radial_distributions.m` | Time-resolved radial distribution heatmap; final inter-particle distance histogram; Boltzmann reference overlay on initial r₀ distribution | Boltzmann curve uses the *existing* `droplet_potential` already in the package — analytic overlay only, no new physics. |
| `HeDFT_MD_comparison_neutral/compare_neutral_dynamics_to_HeDFT.m` | Cumulative-integral trajectory reconstruction (`cumtrapz(v) → r`) overlaying HeDFT v(t) reference; total-energy conservation check across MD vs HeDFT grids | Neutral-side counterpart to existing ion HeDFT comparison. |
| `post_process_single_pulse_paper_IplusHe_comparison.m` | Multi-reference comparative overlay (I⁺ gas vs I⁺He droplet vs simulation) | The existing `plot_paper_figure.py` already does the equivalent overlay; only the explicit "comparison filtering" (`v_proj > VMIN`) is novel — fold as an option. |

### Skip — out of scope or duplicate

| MATLAB | Reason |
| ------ | ------ |
| `single_pulse_simulation/extra_figures.m` | Abel inversion + experimental VMI image processing — out per `CLAUDE.md`. |
| `single_pulse_simulation/HeDFT_comparison/simulation_image_pumpprobe_comp.m` | Pump-probe — out per `CLAUDE.md`. |
| `post_process_single_pulse_paper_gase_phase_comp.m` | Gas-phase / effusive comparison — out per `CLAUDE.md`. |
| `vmi_sim_visualize_distributions.m`, `vmi_sim_visualize_ensemble.m`, `vmi_sim_visualize_particle_paths.m`, `vmi_sim_visualize_trajectory.m` | Live-debug visualizers (3D point clouds, animations); no post-hoc metric. |
| `generate_post_process_figures.m`, `post_process_function.m`, `vmi_sim_post_process.m` | Dispatcher / scaffold; replaced by the new Python driver. |
| `single_pulse_simulation/HeDFT_comparison/analyze_trajectories.m`, `DFT_results_visualization.m`, `plot_HeDFT_neutral_data.m` | Pre-sim mechanics inspection / HeDFT-only viz; not run-output post-processing. |
| `HeDFT_results/*.m`, `I2_potential_energy_curve_check/*.m`, `I2+_potential_energy_curve_check/*.m`, `energy_estimations.m` | Pre-simulation reference / Morse fitting / literature PEC plots — not post-processing of run output. |

## Consolidated Script Design

### Entry point

`scripts/post_processing/plot_run_summary.py`

- CLI: `python plot_run_summary.py <run_dir> [--hedft-ref PATH] [--vmi-ref-he PATH] [--vmi-ref-gas PATH] [--out-dir PATH] [--no-show]`.
- Loads `RunDirectory(run_dir)`; uses `has_neutral()` / `has_ion()` to gate sections.
- References are optional; sections that need a missing reference are skipped with a log line, preserving the "keep different runs separate" rule from `CLAUDE.md` (the user passes only the references appropriate to that run).
- Output: one multi-page PDF `<run_dir>/figures/run_summary.pdf` plus per-figure PNGs alongside it. Existing single-figure scripts (`plot_*_energy_balance.py`, `plot_paper_figure.py`, etc.) stay as-is for focused workflows; the new script just consolidates them and adds the new panels.
- `main(argv=...)` is callable from tests so smoke coverage can drive the
  whole pipeline without touching the CLI parser.

### Panel sections (in PDF order)

1. **Run metadata**: text page with `cfg.json` summary (preset name, N, droplet R, pulse params, run timestamp).
2. **Neutral energy balance** (if `neutral.npz`): reuse `energy_balance.neutral_energy_totals`.
3. **Ion energy balance + temperature diagnostic** (if `ion.npz`): reuse `energy_balance.ion_energy_totals` + the existing temperature plot helper.
4. **Mass spectrum**: reuse `energy_balance.mass_spectrum`.
5. **Final radial velocity histogram** (1D), optional VMI overlay: reuse `velocity_distribution.compute_final_velocity_histogram` + `load_vmi_reference`. Add **bimodal Gaussian fit** overlay when fit converges (new helper, see below).
6. **Phi (azimuthal) histogram** (1D): reuse `energy_balance.phi_histogram`.
7. **2D polar velocity histogram (v, φ)** *(new)*: surface / pcolormesh.
8. **cos² anisotropy fit** *(new)*: fit `a + b·cos(φ)^β` per v-bin to row-wise polar histogram. Also produce **β(v)** line plot.
9. **2D velocity density (vₓ, vᵧ)** *(new)*: 2D histogram of final velocity components in the lab frame.
10. **Mass-resolved final-velocity histograms** for `I⁺`, `I⁺He`, `I⁺He₂`: reuse mass-filtered `compute_final_velocity_histogram`.
11. **Time-resolved radial distribution heatmap** *(new)*: |r| of each ion vs. time, binned into `(n_time_slices, n_r_bins)`.
12. **Final inter-particle distance histogram** *(new)*: per-molecule |rₐ − r_b| at last time step.
13. **Angular pair covariance** *(new)*: 2D histogram of `(θₐ, θ_b)` for each molecular pair, with self-diagonal removed (matches v4 output structure).
14. **HeDFT comparison** (if `--hedft-ref` provided):
    - Ion side: reuse the existing `compare_trajectories` panel.
    - Neutral side *(new)*: add cumulative-integral trajectory reconstruction and total-energy conservation check from `compare_neutral_dynamics_to_HeDFT.m`.
15. **Initial population vs. Boltzmann reference overlay** *(new)*: histogram of `r0` (already on checkpoints) overlaid with `exp(−V(r) / k_B T)` using the package's existing `droplet_potential`. Analytic overlay only.

### New helper modules (under `i2_helium_md/postprocess/`)

| File | Public API | Purpose |
| ---- | ---------- | ------- |
| `polar_velocity.py` | `polar_velocity_histogram(ckpt, n_v_bins, n_phi_bins, mass_amu=None) -> PolarHistogram`, `anisotropy_fit(polar) -> AnisotropyFit`, `beta_of_velocity(polar) -> BetaCurve` | Polar (v, φ) histogram and cos² anisotropy fit using `scipy.optimize.curve_fit`. |
| `velocity_2d.py` | `velocity_density_2d(ckpt, axes=('x','y'), n_bins=200, v_max=...) -> Velocity2DHistogram` | 2D vₓ/vᵧ histogram. |
| `pair_correlation.py` | `interparticle_distance_histogram(ckpt, bins) -> DistanceHistogram`, `angular_pair_covariance(ckpt, n_theta_bins) -> CovarianceMatrix` | Per-molecule pair-distance and pair-angle correlations. |
| `time_resolved.py` | `radial_distribution_evolution(ckpt, n_time_slices, n_r_bins) -> RadialEvolution` | Time-binned \|r\| distribution, used for the heatmap panel. |
| `velocity_distribution.py` (extend) | `bimodal_gaussian_fit(hist: FinalVelocityHistogram) -> BimodalGaussianFit` | Two-Gaussian decomposition; returns `np.nan` parameters on non-convergence. |
| `compare_trajectories.py` (extend) | `compare_neutral_to_hedft(neutral_ckpt: NeutralCheckpoint, hedft: HedftTrajectory) -> NeutralComparison` | Cumulative-integral reconstruction + energy conservation (mirrors the legacy neutral comparison). |
| `boltzmann_overlay.py` | `boltzmann_population(droplet_radius, temperature_K, r_grid) -> BoltzmannCurve` | Analytic Boltzmann curve using `physics.droplet_potential`; no new physics. |

All new dataclasses follow the pattern already used in `velocity_distribution.py` and `compare_trajectories.py` (frozen dataclasses with numpy arrays).

### Re-used existing pieces (no changes)

- `RunDirectory` (`load_cfg`, `load_neutral`, `load_ion`).
- `energy_balance.{ion_energy_totals, neutral_energy_totals, phi_histogram, mass_spectrum}`.
- `velocity_distribution.{load_vmi_reference, compute_final_velocity_histogram}`.
- `compare_trajectories.{compare_distance, compare_velocity_magnitude}` and `hedft_loader.load_hedft_trajectory`.
- `_smoothing.{moving_mean, normalise_trace}`.
- The existing single-figure scripts in `scripts/post_processing/` (kept for focused workflows; the consolidated script imports their underlying helpers, not the scripts themselves).

### Files to be added

- `i2_helium_md/postprocess/polar_velocity.py`
- `i2_helium_md/postprocess/velocity_2d.py`
- `i2_helium_md/postprocess/pair_correlation.py`
- `i2_helium_md/postprocess/time_resolved.py`
- `i2_helium_md/postprocess/boltzmann_overlay.py`
- `scripts/post_processing/plot_run_summary.py`

### Files to be modified

- `i2_helium_md/postprocess/velocity_distribution.py` (add `bimodal_gaussian_fit`).
- `i2_helium_md/postprocess/compare_trajectories.py` (add `compare_neutral_to_hedft`).
- `i2_helium_md/postprocess/__init__.py` (export new symbols).

### Tests to be added (focused, synthetic checkpoints only)

- `tests/test_polar_velocity.py` — polar histogram normalization; cos² fit recovers known β on synthetic anisotropic distribution; β(v) returns NaN on empty bins.
- `tests/test_velocity_2d.py` — counts conservation; symmetry of synthetic isotropic input.
- `tests/test_pair_correlation.py` — pair-distance histogram on hand-built 2-molecule checkpoint; covariance diagonal removed.
- `tests/test_time_resolved.py` — heatmap shape and that summing across time slices recovers per-particle radial histogram.
- `tests/test_bimodal_gaussian_fit.py` — recovers two known peaks within tolerance; fails gracefully on flat input.
- `tests/test_compare_neutral_to_hedft.py` — cumulative-integral reconstruction on a synthetic constant-velocity case.
- `tests/test_boltzmann_overlay.py` — analytic Boltzmann curve normalises to 1 on a uniform potential.
- `tests/test_plot_run_summary_smoke.py` — non-interactive matplotlib smoke test using a tiny synthetic run dir.

No production-sized checkpoints; no figure files committed.

## Verification

End-to-end manual checks before reporting complete:

1. `pytest tests/test_*` — all new tests + existing post-processing tests pass. **Done:** 422 tests green.
2. `python scripts/post_processing/plot_run_summary.py data/runs/9A_hedft_comparison --hedft-ref data/reference/9A_All_Data.csv` — produces `figures/run_summary.pdf` containing HeDFT panels and skipping VMI overlay. **Done:** 17-section PDF (~860 KB).
3. `python scripts/post_processing/plot_run_summary.py data/runs/single_pulse_droplet --vmi-ref-he data/reference/vmi_iplus_he.csv --vmi-ref-gas data/reference/vmi_iplus_gas.csv` — produces VMI overlay, anisotropy, polar histogram, mass-resolved velocity panels, and skips HeDFT panels. **Done:** 15-section PDF (~860 KB).
4. Visually compare the new PDF panels against the legacy MATLAB figures in
   `legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_v3.m`,
   `..._v4.m`,
   `post_process_compare_radial_distributions.m`,
   `HeDFT_MD_comparison_neutral/compare_neutral_dynamics_to_HeDFT.m`
   for qualitative agreement (peak positions, overall shape). **Pending human review.**

## Deferred / Explicitly Out of Scope

- Abel inversion and reconstruction of an experimental polar VMI image (`extra_figures.m`).
- Pump-probe variants (`simulation_image_pumpprobe_comp.m`).
- Effusive / gas-phase MD comparison (`post_process_single_pulse_paper_gase_phase_comp.m`).
- Reference / pre-simulation PEC plots and Morse fits (`HeDFT_results/`, `I2_potential_energy_curve_check/`, `I2+_potential_energy_curve_check/`, `energy_estimations.m`).
- Live-debug 3D scatter / animation viewers (`vmi_sim_visualize_*`).

These remain blocked behind explicit user request per `CLAUDE.md`.

## Progress Tracker

Update this table as items are implemented. Status values: **TODO**, **WIP**, **DONE**, **SKIPPED**.

### New helper modules

| Item | File | Status | Notes |
| ---- | ---- | ------ | ----- |
| Polar (v, φ) histogram | `i2_helium_md/postprocess/polar_velocity.py` | DONE | |
| cos² anisotropy fit | `i2_helium_md/postprocess/polar_velocity.py` | DONE | `anisotropy_fit` |
| β(v) curve | `i2_helium_md/postprocess/polar_velocity.py` | DONE | `beta_of_velocity` |
| 2D velocity density (vₓ, vᵧ) | `i2_helium_md/postprocess/velocity_2d.py` | DONE | |
| Inter-particle distance histogram | `i2_helium_md/postprocess/pair_correlation.py` | DONE | |
| Angular pair covariance | `i2_helium_md/postprocess/pair_correlation.py` | DONE | |
| Time-resolved radial heatmap | `i2_helium_md/postprocess/time_resolved.py` | DONE | |
| Boltzmann overlay | `i2_helium_md/postprocess/boltzmann_overlay.py` | DONE | |
| Bimodal Gaussian fit | `i2_helium_md/postprocess/velocity_distribution.py` (extend) | DONE | |
| Neutral-side HeDFT comparison | `i2_helium_md/postprocess/compare_trajectories.py` (extend) | DONE | `compare_neutral_to_hedft` |
| Package exports | `i2_helium_md/postprocess/__init__.py` | DONE | re-export new symbols |

### Driver script

| Item | File | Status | Notes |
| ---- | ---- | ------ | ----- |
| Consolidated run-summary driver | `scripts/post_processing/plot_run_summary.py` | DONE | CLI + multi-page PDF |
| Section: run metadata | (driver) | DONE | |
| Section: neutral energy balance | (driver) | DONE | reuse existing helpers |
| Section: ion energy balance + temperature | (driver) | DONE | reuse existing helpers |
| Section: mass spectrum | (driver) | DONE | reuse `energy_balance.mass_spectrum` |
| Section: 1D radial velocity + VMI overlay + bimodal fit | (driver) | DONE | |
| Section: 1D phi histogram | (driver) | DONE | |
| Section: 2D polar (v, φ) histogram | (driver) | DONE | |
| Section: cos² anisotropy fit + β(v) | (driver) | DONE | |
| Section: 2D velocity density (vₓ, vᵧ) | (driver) | DONE | |
| Section: mass-resolved final velocities | (driver) | DONE | |
| Section: time-resolved radial heatmap | (driver) | DONE | |
| Section: final inter-particle distance histogram | (driver) | DONE | |
| Section: angular pair covariance | (driver) | DONE | |
| Section: HeDFT comparison (ion + neutral) | (driver) | DONE | reference-gated |
| Section: Boltzmann initial-population overlay | (driver) | DONE | |

### Tests

| Item | File | Status | Notes |
| ---- | ---- | ------ | ----- |
| Polar velocity / anisotropy | `tests/test_polar_velocity.py` | DONE | 7 tests |
| 2D velocity density | `tests/test_velocity_2d.py` | DONE | 4 tests |
| Pair correlation | `tests/test_pair_correlation.py` | DONE | 5 tests |
| Time-resolved radial | `tests/test_time_resolved.py` | DONE | 3 tests |
| Bimodal Gaussian fit | `tests/test_bimodal_gaussian_fit.py` | DONE | 2 tests |
| Neutral HeDFT comparison | `tests/test_compare_neutral_to_hedft.py` | DONE | 4 tests |
| Boltzmann overlay | `tests/test_boltzmann_overlay.py` | DONE | 5 tests |
| Run-summary smoke test | `tests/test_plot_run_summary_smoke.py` | DONE | 2 tests, matplotlib `Agg` |

### End-to-end verification

| Step | Status | Notes |
| ---- | ------ | ----- |
| `pytest tests/` all green | DONE | 422 passed |
| Run on `data/runs/9A_hedft_comparison` with HeDFT ref | DONE | 17 sections incl. HeDFT, VMI overlay skipped |
| Run on `data/runs/single_pulse_droplet` with VMI refs | DONE | 15 sections incl. VMI, HeDFT comparison skipped |
| Visual comparison vs. legacy MATLAB figures | TODO | Awaiting human review of `data/runs/*/figures/run_summary.pdf` |
