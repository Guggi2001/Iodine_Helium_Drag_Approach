# The `droplet_sizes.py` module — a walkthrough

## What problem does this file solve?

Every simulated I₂ molecule in this project is embedded in a helium droplet.
Real helium nanodroplet experiments produce droplets in a **distribution of
sizes**, set by the source nozzle conditions (pressure, temperature, nozzle
diameter). Realistic simulations should sample droplet sizes from this
distribution.

This module provides:

1. **An empirical formula** for the mean droplet size given source conditions.
2. **A log-normal sampler** that draws individual droplet sizes from the
   distribution.
3. **A pickup-cell simulator** that models how I₂ pickup biases and reshapes
   the size distribution before the experiment proper starts.

## Position in the dependency chain

```
config.py
   ↓
sampling/droplet_sizes.py    ← THIS MODULE
   ↓
simulation/neutral.py        (uses N to set per-molecule droplet radii)
```

It's the first member of the **sampling** subpackage. Sampling is its own
layer because it sits between configuration ("what kind of droplets do
we want?") and simulation ("what is each droplet's size?").

## When this is actually used

For the **`single_pulse_N2000` preset** (our default scope), the simulation
sets `use_single_droplet_size = True` and assigns every molecule the same
droplet size (2000 atoms). The full sampler is bypassed.

We port the sampler anyway because:

1. The config exposes `use_single_droplet_size = False` for realistic runs.
2. Future thesis work and parameter scans need realistic distributions.
3. Tests of the sampler give us a clean entry point into the sampling layer.

## Public API

```python
from i2_helium_md.sampling.droplet_sizes import (
    mean_droplet_size,        # empirical correlation: returns (N_mean, R_mean)
    droplet_radius_from_N,    # convert atom count to radius in Angstrom
    sample_droplet_sizes,     # the main sampler
)
```

Typical use:

```python
from i2_helium_md import single_pulse_N2000
from i2_helium_md.sampling.droplet_sizes import sample_droplet_sizes

cfg = single_pulse_N2000(num_molecules=500, seed=42, use_single_droplet_size=False)

N = sample_droplet_sizes(cfg, mode="post_pickup")    # shape (500,)
# N now contains 500 droplet sizes drawn from a realistic post-pickup distribution.
# Pass this to the simulation init code (next step) to set droplet_radii.
```

## The two modes

### `mode="raw"`

Pure log-normal sampling. The mean and width come from:

- **Mean N**: empirical formula `N = k1 * p^k2 * T^k3 * d^k4` with constants
  `k1=4e5, k2=0.97, k3=-3.88, k4=2`. From the Lackner thesis and literature.
- **Log-normal width δ = 0.625**: a fixed value from Kornilov (Mol. Phys. 2009)
  that fits experimental size distributions across a range of conditions.
- **Log-normal location μ**: derived from the relation `mu = log(N_mean) - δ²/2`,
  which makes the **arithmetic mean** of the log-normal equal to N_mean.

Use this mode when you want the **theoretical** distribution as it leaves
the nozzle, before any I₂ pickup.

### `mode="post_pickup"` (default)

The realistic mode: also simulates I₂ pickup and resulting evaporation.

**Why this matters physically:**

When a fast-moving I₂ molecule enters a droplet, it deposits its kinetic
energy (~26 meV at 300 K) plus its solvation energy (~14 meV). This energy
is carried away by evaporating helium atoms — about 50 He atoms per pickup
event. So the droplets that actually contain I₂ are smaller than the raw
sampled droplets.

**However**, the pickup probability scales with the geometric cross section
of the droplet (`σ ∝ R² ∝ N^(2/3)`), so larger droplets are preferentially
"selected." The post-pickup distribution is the result of:

- Larger droplets having higher pickup probability (selection bias toward N).
- Evaporation removing ~50 atoms per pickup (slight bias toward lower N).

Net effect: the **mean** post-pickup N can be larger than the raw mean.
The **shape** of the distribution is also altered.

In our smoke run with default conditions (p=40 mbar, T=14 K):

| Mode | Mean N |
|---|---|
| raw          | ~12,700 |
| post_pickup  | ~16,300 |

The pickup simulation is iterative: each round, droplets either pick up a
molecule (and then evaporate atoms) or pass through. Droplets with exactly
one pickup are kept; multi-pickup events and destroyed droplets (N≤0 after
evaporation) are discarded. The function generates **10× the requested
number** of raw samples to ensure enough one-pickup droplets remain after
selection.

## Internal walkthrough

### `mean_droplet_size(p_mbar, T_K, d_um=5)`

Direct port of `get_dropletsize.m`. One-line empirical formula. Returns
both the mean atom count and the mean droplet radius.

### `droplet_radius_from_N(N)`

A standalone helper for converting a single atom count to a droplet radius,
assuming the same density (0.8 of bulk He, ~0.0175 atoms/Å³). Used in
several places in the simulation init.

### `sample_droplet_sizes(cfg, mode, ...)`

The user-facing function. Generates either raw log-normal samples or
post-pickup samples, with a reproducible RNG.

**Reproducibility** comes from `cfg.seed` if set, otherwise a fresh
`np.random.default_rng()` is used. You can also pass an `rng` argument
to share a single RNG with other samplers (the radial position sampler,
for example).

### `_simulate_pickup(N_raw, num_target, ...)` (private)

The pickup-cell simulator. Mutates working arrays in a while-loop until
the maximum pickup count (1) is reached or all droplets have been
processed. Returns a list of "completed" droplets that ended with exactly
one pickup.

**Termination guards:**

- Loop exits when no more droplets can pick up another molecule.
- A `max_iters = 100` safety guard prevents pathological cases from
  hanging.
- If too few one-pickup droplets are found (`one_pickup.size < num_target`),
  raises `ValueError` with an actionable message instead of returning
  silently truncated data — a deliberate departure from the legacy MATLAB
  function which would have returned a shorter array via index error.

### `_droplet_radius_m`, `_impact_parameter_threshold_m`,
###  `_pickup_cross_section_m2`, `_evaporation_per_pickup` (private)

These reproduce the four anonymous-function helpers from the MATLAB code
verbatim, with names that explain their physical meaning. Each is
independently testable.

## Departures from MATLAB

These are intentional and documented:

1. **No plotting code.** The legacy function generated histograms inline.
   Visualisation is the job of post-processing modules, not samplers.

2. **Single function instead of two.** MATLAB had `generate_droplet_sizes.m`
   and `generate_droplet_sizes_simpler.m` (which was just the early-exit
   path of the full version). We unify them into one function with a
   `mode` argument.

3. **Hard error on insufficient samples, with auto-retry.** MATLAB would
   `error()` only at the very end. We do the same, but with auto-retry up
   to `max_retries` rounds (each at 4× the previous oversample factor).
   Only then do we raise, with a clearer message that points to the fix.

4. **`num_evap` returns negative `dN`** instead of MATLAB's positive
   "evaporated count." This way the same expression `samples + dN`
   updates droplet sizes correctly in both modes -- matches the MATLAB
   behaviour but with clearer intent.

5. **Reproducibility via `np.random.Generator`.** MATLAB had no clean way
   to seed `lognrnd` without using a global RNG. Python's `Generator` is
   per-call and thread-safe.

## Thesis figure 3.2 reproduction

The supervisor's thesis includes a figure (3.2) showing post-pickup
droplet size distributions at T = 12, 15, 18 K with both normal and
reduced cross sections.

**The figure was produced by a separate analytical script, not by our
Monte Carlo sampler.** We have two different code paths for the two
different physical models:

* :func:`sample_droplet_sizes` -- our Monte Carlo with explicit
  evaporation. Used to set initial conditions for the actual MD
  simulation. Default ``E_solv = 14 meV`` (matches the legacy production
  MATLAB code).
* :func:`conditional_size_distributions_analytical` -- literal port of
  Treiber's didactic script
  (``conditional_droplet_size_distribution_simplified.m``). Returns
  **both** the normal-σ and reduced-σ conditional distributions as a
  tuple. Used only for thesis-figure reproduction and analytical
  comparison.

Two implementation details that matter for the analytical port (both
caught after the user pushed back on a wrong reproduction; see
`migration_log.md` for the full debugging story):

1. **Density convention.** Treiber's script uses **bulk** helium density
   ``n_he = 2.18e28`` for ``R(N) = (3N / 4π n_he)^(1/3)``, NOT the
   ``0.8 × n_he`` droplet density used elsewhere. Our analytical function
   exposes ``n_he_per_m3`` as a parameter with Treiber's default.

2. **Normalisation convention.** Both branches share the **same**
   normalisation (``p_k_normalization`` from the normal-σ Poisson
   convolution). The reduced-σ distribution therefore integrates to
   **less than 1** -- its integral is the fraction of one-pickup events
   that come from droplets above the kinetic-energy threshold. This is
   why the dashed peaks in the thesis figure can be smaller than the
   solid peaks at high T (where the threshold rejects most droplets) but
   taller at low T (where the threshold passes most droplets).

The diagnostic helper :func:`droplet_sizes_diagnostics.plot_thesis_figure_3_2`
reproduces the thesis figure visually, peak-by-peak. The reproduction
matches pixel-for-pixel.

## Testing

The test file `tests/test_droplet_sizes.py` covers:

- **Formula port**: spot-check that `mean_droplet_size(40, 23)` matches a
  hand calculation, and that the dependence on p, T, d goes the right way.
- **Sampler shape and reproducibility**: same seed → same output.
- **Statistical properties**: with 50,000 samples, the arithmetic mean is
  within 5% of the target N_mean, and `log(N)` is approximately normal
  (|skewness| < 0.1).
- **Post-pickup distribution**: returns positive samples in a sensible
  range, errors loudly on impossible requests.
- **Internal helper sanity**: radius monotonically increasing, reduced
  cross section ≤ geometric, evaporation always negative, etc.

## Future improvements

1. **Match-to-data calibration.** The empirical constants `k1...k4` come
   from a single source. A future improvement would be to expose them as
   parameters and provide a fit-to-data utility for re-calibration.

2. **Support multi-pickup post-distributions.** `max_pickup = 1` is
   currently hardcoded. Adding multi-pickup support is straightforward
   (the loop already tracks it) and would enable simulating cluster
   experiments.

3. **Vectorised iteration with NumPy ufuncs only.** The current
   implementation uses boolean masks and array reassignment in a
   `while` loop. For `num_molecules > 100,000` this could become slow.
   Switching to a single vectorised pickup probability per droplet
   (without iteration) is possible if we accept some loss of physical
   fidelity.

## References

- Lackner, F., dissertation, TU Graz — nozzle correlation derivation.
- Kornilov, O. & Toennies, J.P., *Mol. Phys.* 107, 2071 (2009) — log-normal width.
- Yang et al., 10.1007/978-3-030-94896-2_1 — He droplet chemical potential.
