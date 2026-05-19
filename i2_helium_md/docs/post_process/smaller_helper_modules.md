# Smaller post-processing helper modules

## What this file covers

The four large post-processing helpers each have their own walkthrough:

- `compare_trajectories.py` — see `compare_trajectories_module.md`
- `energy_balance.py` — see `energy_balance_module.md`
- `hedft_loader.py` — see `hedft_loader_module.md`
- `velocity_distribution.py` — see `velocity_distribution_module.md`

This document covers the **six remaining helper modules** under
`i2_helium_md/postprocess/`. They are smaller (one or two public
entry points each) but every in-scope plotting script reaches for at
least one of them, so they need a single place where their contract,
units, and MATLAB legacy origin are pinned.

| Module | One-line purpose |
|---|---|
| `_smoothing.py` | MATLAB-style `movmean` + trace baseline/unit-max normaliser |
| `pair_correlation.py` | I-I distance histogram + angular pair covariance |
| `polar_velocity.py` | (\|v\|, φ) polar histogram + cos² anisotropy fit |
| `velocity_2d.py` | 2-D (vx, vy) heat-map of final ion velocities |
| `time_resolved.py` | Radial-distribution heat-map over simulation time |
| `boltzmann_overlay.py` | Analytic `exp(-V/kT)` reference curve over `droplet_potential` |

Strategy-level framing (Strategy A vs Strategy B, projection conventions,
why the recipes differ between scripts) lives in
`post_processing_strategy.md`. Project-wide scope rules live in
`CLAUDE.md`.

## Position in the dependency chain

```
IonCheckpoint / NeutralCheckpoint
        │
        ▼
postprocess/    _smoothing.py
                pair_correlation.py
                polar_velocity.py        ← THIS FILE
                velocity_2d.py
                time_resolved.py
                boltzmann_overlay.py
        │
        ▼
scripts/post_processing/
        plot_experimental_comparison.py
        plot_paper_v2.py / v3 / v4
        plot_paper_cov.py
        plot_run_summary.py
```

All six modules are pure functions over the saved checkpoints (and, for
the Boltzmann curve, over `cfg.json`). No plotting, no I/O.

---

## `_smoothing.py` — MATLAB `movmean` + trace normalisation

### What problem does this file solve?

The legacy paper figures apply a 15-bin centred moving mean to the
simulation curves before normalising them to unit maximum. Both
`plot_experimental_comparison.py` and `plot_paper_v3.py` (and the
consolidated `plot_run_summary.py`) need exactly that recipe. Keeping
it in a single module gives unit tests one place to pin the
convention and avoids the script-to-script imports the legacy MATLAB
code resorts to.

Mirrors MATLAB `movmean(..., 15)` (centred window, shortened endpoint
windows) plus the implicit `(x - min(x)) / max(x - min(x))` baseline
normalisation used on the paper plots.

### Public API

```python
from i2_helium_md.postprocess._smoothing import (
    moving_mean,
    normalise_trace,
)

smoothed     = moving_mean(counts, window=15)   # 1-D ndarray, same shape as input
normalised   = normalise_trace(smoothed)        # 1-D ndarray, in [0, 1]
```

Both functions accept array-like input, return `np.ndarray` of dtype
`float`, and operate on the input as a 1-D sequence.

### Contract

| Concern | Behaviour |
|---|---|
| `window <= 1` | Returns a copy of the input (no smoothing) |
| Empty input | Returns an empty `np.ndarray` |
| Centred window with shortened endpoints | Matches MATLAB `movmean` (slice `[max(0, i - half_left) : min(N, i + half_right + 1)]`) |
| Constant input → `normalise_trace` | Returns zeros (no divide-by-zero) |
| Baseline | `normalise_trace` subtracts `min(values)` before dividing by `max(shifted)` |

The companion `energy_balance_module.md` references this module from
its "Companion module: `_smoothing.py`" section; the contract is the
same.

### Used by

- `scripts/post_processing/plot_experimental_comparison.py`
- `scripts/post_processing/plot_paper_v3.py`
- `scripts/post_processing/plot_run_summary.py`

---

## `pair_correlation.py` — distance histogram + angular pair covariance

### What problem does this file solve?

Two final-state pair-correlation diagnostics for an ion run:

1. **Inter-particle distance histogram** — mirrors the
   `histogram_data_interatomic_distance` block of
   `legacy_matlab_repository/post_process_compare_radial_distributions.m`.
2. **Angular pair covariance** — mirrors the
   `simulated_angular_covariance` block of
   `legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_v4.m`.
   The MATLAB script bins `atan2(vx, vy)` for each ion pair into a
   2-D histogram and removes the diagonal so the covariance structure
   shows inter-atom angular correlation rather than trivial
   self-binning.

The MATLAB atom-layout convention is reused exactly: indices `[0, N)`
are the first iodine atom of each molecule, indices `[N, 2N)` the
second, so the per-molecule pair is `(atom i, atom i + N)`.

### Public API

```python
from i2_helium_md.simulation.run_directory import RunDirectory
from i2_helium_md.postprocess import (
    interparticle_distance_histogram,
    angular_pair_covariance,
)

ion = RunDirectory("data/runs/single_pulse_droplet").load_ion()

dist = interparticle_distance_histogram(ion, num_bins=100)
cov  = angular_pair_covariance(
    ion,
    n_theta_bins=50,
    mass_amu=131.0,         # optional mass channel filter
    require_outside=True,   # default
    remove_diagonal=True,   # default — matches MATLAB cov_angular - diag(...)
)
```

### Returned dataclasses

```python
@dataclass(frozen=True)
class DistanceHistogram:
    bin_centers_A:  np.ndarray   # shape (B,)  angstrom
    bin_edges_A:    np.ndarray   # shape (B+1,)
    counts:         np.ndarray   # shape (B,)  int
    num_pairs_used: int

@dataclass(frozen=True)
class CovarianceMatrix:
    counts:            np.ndarray  # shape (n, n), float
    theta_centers_rad: np.ndarray  # shape (n,)   radians, in [0, 2 pi)
    theta_edges_rad:   np.ndarray  # shape (n+1,)
    num_pairs_used:    int
```

### Contract and conventions

| Concern | Behaviour |
|---|---|
| Distance source | `positions_final_{x,y,z}` (asymptotic position, not the last trajectory sample). |
| `max_distance_A=None` | Auto-sets to `1.05 × nanmax(|r_a - r_b|)`. |
| Mass filter (covariance) | Atom-level mask `|round(mass / U) - mass_amu| ≤ tolerance` applied to both atoms of each molecule; a single failing atom drops the pair. |
| `require_outside=True` | Both atoms must satisfy `b_ion_outside` for the molecule to enter the sample. |
| φ convention | `theta = arctan2(vx, vy) + π`, wrapped into `[0, 2π)` — same bin layout as `energy_balance.phi_histogram`. |
| `remove_diagonal=True` | Sets `counts.diagonal() = 0`. Set to `False` to recover the raw 2-D histogram for debugging. |

### Used by

- `scripts/post_processing/plot_paper_v3.py`
- `scripts/post_processing/plot_paper_v4.py`
- `scripts/post_processing/plot_paper_cov.py`
- `scripts/post_processing/plot_run_summary.py`

---

## `polar_velocity.py` — (\|v\|, φ) histogram + cos² anisotropy

### What problem does this file solve?

Reproduces the simulation side of the polar VMI panels of
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_v3.m`
(lines 96–111): an angular distribution plus a `nlinfit` cos² model
plus a per-velocity `beta(v)` curve.

The legacy MATLAB code computed an anisotropy fit on top of an
Abel-inverted 2-D experimental VMI image. We have the full 3-D
simulated velocities, so no Abel inversion is needed — we bin the
lab-frame final velocities directly and fit

```text
f(phi) = a + b * cos(phi - phi0)^2
```

(equivalent to a Legendre-P2 expansion). The conventional
photodissociation anisotropy parameter is recovered as

```text
beta = 2 * b / (2 * a + b)         # range [-1, 2]; +2 = pure cos^2, -1 = pure sin^2
```

### Public API

```python
from i2_helium_md.postprocess import (
    polar_velocity_histogram,
    anisotropy_fit,
    beta_of_velocity,
)

polar = polar_velocity_histogram(
    ion,
    n_v_bins=80,
    n_phi_bins=72,
    v_max_Aps=28.0,
    mass_amu=131.0,            # optional
    mass_tolerance_amu=0.5,
    require_outside=True,
)

fit   = anisotropy_fit(polar)                                 # one fit over all |v|
betas = beta_of_velocity(polar, min_counts_per_v_bin=50)      # per-|v| beta curve
```

### Returned dataclasses

```python
@dataclass(frozen=True)
class PolarHistogram:
    counts:          np.ndarray  # shape (n_v, n_phi)
    v_centers_Aps:   np.ndarray
    v_edges_Aps:     np.ndarray
    phi_centers_rad: np.ndarray  # in [0, 2 pi)
    phi_edges_rad:   np.ndarray
    mass_amu:        float       # NaN when mass_amu was None
    num_atoms_used:  int

@dataclass(frozen=True)
class AnisotropyFit:
    a:        float
    b:        float
    phi0_rad: float
    beta:     float    # NaN if (2 a + b) == 0 or fit failed
    residual: float
    success:  bool

@dataclass(frozen=True)
class BetaCurve:
    v_centers_Aps:    np.ndarray
    beta:             np.ndarray  # NaN where the per-bin fit was skipped or failed
    beta_uncertainty: np.ndarray  # residual norm, NaN where invalid
    valid:            np.ndarray  # bool mask, True where the fit converged
```

### Conventions

- φ = `arctan2(vy_final, vx_final) + π`, wrapped into `[0, 2π)` so the
  bin layout matches `energy_balance.phi_histogram` and the angular
  covariance from `pair_correlation.py`.
- 1-D radial recovery: `polar.counts.sum(axis=1)`.
- 1-D φ recovery: `polar.counts.sum(axis=0)`.
- `min_counts_per_v_bin=50` is the CLAUDE.md default for `beta(v)`;
  bins with fewer counts get `NaN` and `valid=False`. Lower it
  cautiously on small samples.
- The `cos²` fit uses `scipy.optimize.curve_fit` with
  `p0 = (min(counts), max(counts) - min(counts), phi[argmax(counts)])`
  and a `maxfev=5000` limit. On convergence failure it returns an
  `AnisotropyFit` with `success=False` and all-NaN parameters rather
  than raising.

### Used by

- `scripts/post_processing/plot_paper_v3.py`
- `scripts/post_processing/plot_paper_v4.py`
- `scripts/post_processing/plot_paper_cov.py`
- `scripts/post_processing/plot_run_summary.py`

---

## `velocity_2d.py` — 2-D (vx, vy) heat-map

### What problem does this file solve?

Reproduces the (vx, vy) heat-map panel of
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse.m`
(`velocity_bins = -22:0.5:22`, `hist2(vx_total, vy_total)`). A plain
rebin of the existing 3-D simulated final velocities — no new physics,
no Abel inversion.

### Public API

```python
from i2_helium_md.postprocess import velocity_density_2d

heat = velocity_density_2d(
    ion,
    axes=("x", "y"),           # any two distinct axes of {"x", "y", "z"}
    n_bins=200,
    v_max_Aps=22.0,            # symmetric ±v_max grid
    mass_amu=131.0,            # optional
    mass_tolerance_amu=0.5,
    require_outside=True,
)
```

### Returned dataclass

```python
@dataclass(frozen=True)
class Velocity2DHistogram:
    counts:           np.ndarray  # shape (n_bins, n_bins)
    bin_centers_a_Aps: np.ndarray
    bin_edges_a_Aps:   np.ndarray
    bin_centers_b_Aps: np.ndarray
    bin_edges_b_Aps:   np.ndarray
    axis_a:           str          # "x" / "y" / "z"
    axis_b:           str
    mass_amu:         float        # NaN when mass_amu was None
    num_atoms_used:   int
```

### Contract

| Concern | Behaviour |
|---|---|
| `axes[0] == axes[1]` | `ValueError`. |
| `v_max_Aps <= 0` | `ValueError`. |
| `n_bins < 1` | `ValueError`. |
| Default grid | Matches legacy `-22:0.5:22` (`n_bins=88` if you want the exact bin width; `n_bins=200` is the in-package default and is the choice the plotting scripts use). |
| Mass + outside selection | Identical to `velocity_distribution.compute_final_velocity_histogram` and `pair_correlation.angular_pair_covariance`. |

### Used by

- `scripts/post_processing/plot_paper_v2.py`
- `scripts/post_processing/plot_paper_v3.py`
- `scripts/post_processing/plot_run_summary.py`

---

## `time_resolved.py` — radial distribution heat-map over time

### What problem does this file solve?

Reproduces the time-evolution radial-distribution panel of
`legacy_matlab_repository/post_process_compare_radial_distributions.m`
(`histogram_data_radius(i, :)` looped over time indices).

For each of `n_time_slices` uniformly-spaced time samples, bins per-atom
lab-frame `|r|` into `n_r_bins`. Time slices are sub-sampled across the
stored time axis, not per stored simulation step, so the heat-map stays
affordable on long runs.

### Public API

```python
from i2_helium_md.postprocess import radial_distribution_evolution

# Works on either checkpoint type — only reads positions_{x,y,z} + time_ps.
heat = radial_distribution_evolution(
    ion,                       # or neutral checkpoint
    n_time_slices=50,
    n_r_bins=100,
    r_max_A=None,              # autoscale to 1.05 * nanmax(|r|)
)
```

### Returned dataclass

```python
@dataclass(frozen=True)
class RadialEvolution:
    counts:          np.ndarray  # shape (n_time_slices, n_r_bins)
    time_centers_ps: np.ndarray  # shape (n_time_slices,)
    time_indices:    np.ndarray  # indices into ckpt.time_ps that were sampled
    r_centers_A:     np.ndarray  # shape (n_r_bins,)
    r_edges_A:       np.ndarray  # shape (n_r_bins + 1,)
```

### Contract

| Concern | Behaviour |
|---|---|
| `n_time_slices > num_steps` | Clamped silently to `num_steps`. |
| `n_time_slices < 1` or `n_r_bins < 1` | `ValueError`. |
| Empty time axis | `ValueError("checkpoint has no time samples")`. |
| `r_max_A=None` | Auto-sets to `1.05 × nanmax(|r|)` over the selected slices; falls back to `1.0` if no finite radii exist. |
| Accepts both `IonCheckpoint` and `NeutralCheckpoint` | Only the position arrays and the time axis are read; nothing ion-specific. |

### Used by

- `scripts/post_processing/plot_run_summary.py`

(The MATLAB origin script — the radial-distribution-evolution figure —
is consolidated into `plot_run_summary.py` rather than getting its own
focused Python script.)

---

## `boltzmann_overlay.py` — analytic Boltzmann reference curve

### What problem does this file solve?

Reproduces the `p_boltzmann = exp(-V/(k_B*T))` overlay block of
`legacy_matlab_repository/post_process_compare_radial_distributions.m`.

The potential `V(r)` is the simulation's own
`i2_helium_md.physics.potentials.droplet_potential` — the *same*
`V(r)` the sampler in `i2_helium_md/sampling/` and the neutral
propagation already use. This module just normalises
`exp(-V / k_B T)` onto a user-chosen radial grid so it can be plotted
on top of the initial-population histogram as an analytic reference
curve. **No new physics.**

### Public API

```python
from i2_helium_md.postprocess import boltzmann_population

curve = boltzmann_population(
    droplet_radius_A=18.0,
    temperature_K=0.37,
    steepness_A=0.5,
    binding_energy_eV=2.6e-4,
    r_grid_A=None,           # default: linspace(0, 2*R, n_points)
    n_points=400,
)
```

### Returned dataclass

```python
@dataclass(frozen=True)
class BoltzmannCurve:
    r_grid_A:        np.ndarray
    density:         np.ndarray  # trapezoidal integral over r_grid_A == 1
    unnormalised:    np.ndarray  # bare exp(-V / k_B T)
    droplet_radius_A: float
    temperature_K:    float
```

### Contract

| Concern | Behaviour |
|---|---|
| `V(r)` argument | `droplet_potential(r - R, steepness=..., binding_energy=...)`; `r = R` is the droplet surface. |
| Energy unit conversion | `k_B * T / eV` happens inside the function (uses `physics.constants.K_B`, `EV`). |
| Normalisation | Trapezoidal integral on `r_grid_A` set to 1. Falls back to zeros if the integral is zero (e.g. all-zero exponent). |
| `temperature_K <= 0` or `droplet_radius_A <= 0` | `ValueError`. |
| User-supplied `r_grid_A` with `< 2` samples | `ValueError`. |

### Convention note (per CLAUDE.md)

For the initial-state overlay use `cfg.potential_steepness_molecule`
for `steepness_A` and convert `cfg.binding_energy_molecule_K` to eV
(`E_eV = k_B * T_K / eV`). Normalise the histogram by trapezoidal
integration on the same `r_grid_A` so the two curves share a y-axis
scale.

### Used by

- `scripts/post_processing/plot_run_summary.py`
- `scripts/post_processing/plot_paper_cov.py`

(via the package `__init__` re-exports.)

---

## Deferred / out of scope for this file

- **Paper reference loaders** — `paper_v2.py`, `paper_v3.py`,
  `paper_v4.py`, `paper_cov.py`. These are CSV/MAT loaders +
  curve-extraction helpers tied to specific paper figures. Their
  contract, expected file layout, and MATLAB provenance are documented
  alongside the figures they feed:
  - `scripts/plot_paper_v2.md`
  - `scripts/plot_paper_v3.md`
  - `scripts/plot_paper_v4.md`
  - `scripts/plot_paper_cov.md`
- **Paper plotting helpers** — `paper_v2_plotting.py`,
  `paper_cov_plotting.py`. These are matplotlib glue with a single
  consumer each (`plot_paper_v2.py` / `plot_run_summary.py`,
  `plot_paper_cov.py` / `plot_run_summary.py`). Their layout
  decisions are discussed in the same script docs above.
- **Project-wide scope rules** — Abel inversion, pump-probe, effusive
  comparison, and full experimental VMI image interpretation remain
  out of scope per `CLAUDE.md`.
