# The `droplet_sizes_diagnostics.py` module

## What problem does this file solve?

The legacy MATLAB `generate_droplet_sizes.m` had a `debug_pickup_plot=true`
flag that, when enabled, would interleave histograms and `fprintf` calls
into the production sampler. We split that responsibility into a dedicated
diagnostics module so the production `droplet_sizes.py` stays free of any
plotting code (per project principle: "samplers don't import matplotlib").

## Public API

```python
from i2_helium_md.sampling.droplet_sizes_diagnostics import diagnose_pickup

cfg = single_pulse_N2000(num_molecules=2000, p_source_mbar=25, T_source_K=15)
diagnostics, fig = diagnose_pickup(cfg, reduced_crosssection=False)
fig.savefig("pickup_diagnostics.png")
```

Returns:
- A **dict** of per-round physics + raw / one-pickup / no-pickup arrays
- A **4-panel matplotlib Figure** (or None if `plot=False`)

Also prints a per-round summary table to stdout (suppressible with
`print_report=False`).

## What the figure shows

The figure has four panels:

| Panel | Content | MATLAB equivalent |
|---|---|---|
| Top-left | Raw log-normal distribution + formula `<N>` | (new — provides reference) |
| **Top-right** | **Initial vs one-pickup distribution (raw counts)** | **MATLAB lines 226-229** |
| Bottom-left | Per-round yield (alive, picked, destroyed) | MATLAB per-round histograms |
| **Bottom-right** | **Bayesian-smoothed density of `new_samples` per round** | **MATLAB lines 191-196 (`bayes_hist` + `plot`)** |

The top-right panel is the most informative: it overlays the initial
log-normal distribution with the post-pickup one-molecule distribution
**as raw counts** (not density), exactly like MATLAB's
``histogram(samples_initial); histogram(samples_1_pickup)`` default. The
one-pickup histogram appears smaller in absolute count because most
droplets don't end up with exactly one pickup, and biased toward
larger N because the geometric cross section grows with droplet size.

If you instead want a normalised shape comparison (each distribution
integrating to 1), bin the arrays yourself with `density=True`:

```python
import matplotlib.pyplot as plt
diag, _ = diagnose_pickup(cfg, plot=False)
plt.hist(diag["N_raw"],        bins=80, density=True, alpha=0.5)
plt.hist(diag["N_one_pickup"], bins=80, density=True, alpha=0.5)
```

But for verification against the MATLAB reference plot, raw counts is
the right comparison.

### Panel 4 — Bayesian-smoothed per-round densities

Panel 4 plots the **density of `new_samples` after each round** of the
pickup loop, using Laplace-smoothed binning (Bayesian histogram):

    p_i = (n_i + 1) / (N + nbins + 1)
    h_i = p_i / binwidth

This is a literal port of MATLAB's `bayes_hist()` helper. Smoothing
ensures no bin has zero density even when only a few samples land in
it -- useful when the pickup loop kills many droplets and per-round
counts are sparse.

For the production case (`max_pickup=1`) there's only one round, so
Panel 4 shows a single curve. The y-axis units are 1/N (probability
density), and the curve integrates to ~1 modulo Laplace smoothing.

The shape of the curve answers: "after one round of Poisson pickup +
evaporation, what does the distribution of droplets that picked up
exactly one molecule look like?" Compared to Panel 2's count plot,
Panel 4 normalises away the absolute number (which can be small if
yield is low) and lets you see the **shape** of the conditional
distribution -- e.g. whether evaporation has shifted the mode left
or whether the cross-section weighting has shifted it right.

## What the printed report shows

Mirrors the MATLAB `fprintf` style:

```
Pickup-cell simulation diagnostics
  source: p=25 mbar, T=15 K, d=5.0 um
  E_solv = 14.0 meV, reduced sigma = False
  formula <N> = 6205

  raw samples :    20000, mean N = 6237
  one-pickup  :     3050, mean N = 8063
  no-pickup   :    16950

  round | n_alive  | <p_pickup> | <dN evap> | picked | destroyed
  --------------------------------------------------------------
      1 |   20000  |   0.1535   |    -44    |  3050  |     0
```

The `mean p_pickup` and `mean evaporated atoms` rows directly correspond to
the MATLAB `fprintf` calls at lines 153-154 of the legacy code.

## Implementation notes

### Tolerance to low yields

Unlike the production sampler, the diagnostic function does **not** raise
on insufficient single-pickup yield. The whole point of running diagnostics
is often *to investigate why* yield is low, so we'd be defeating the
purpose by raising. Internally this is achieved by:

```python
_simulate_pickup(..., return_diagnostics=True)
```

When `return_diagnostics=True`, the count check is skipped and whatever
samples were produced are returned alongside the diagnostics dict.

### Dependency on the production sampler

This module imports `_simulate_pickup` from `droplet_sizes.py` directly,
so the diagnostic and production paths share **exactly the same physics**
— there's no second copy of the pickup loop to drift out of sync.

### Optional matplotlib import

`matplotlib.pyplot` is imported **inside** `_make_pickup_figure`, not at
module load time. This way the diagnostics module can be imported in
headless contexts (e.g. CI without matplotlib) and only the plotting
function will fail if matplotlib isn't installed.

## Thesis figure 3.2 reproduction

The diagnostics module also exposes :func:`plot_thesis_figure_3_2`, which
reproduces the post-pickup distribution figure from chapter 3 of the
supervisor's thesis using the literal Poisson-convolution formula
:func:`droplet_sizes.conditional_size_distributions_analytical`.

```python
from i2_helium_md.sampling.droplet_sizes_diagnostics import plot_thesis_figure_3_2

fig = plot_thesis_figure_3_2()
fig.savefig("thesis_3_2_reproduction.png")
```

The defaults match Treiber's didactic script
(``conditional_droplet_size_distribution_simplified.m``) exactly:
``p_source = 40 mbar``, ``d = 5 μm``, ``E_solv = 30 meV``,
``E_kin_thermal = 38.78 meV``, ``p_pickup_gas = 8.7e-5 mbar``,
``T_pickup_gas = 293 K``, and crucially ``n_he = 2.18e28`` (bulk, NOT
the 0.8x droplet density used in the production sampler).

These parameters and conventions are **different** from
:func:`sample_droplet_sizes`'s defaults. Both are valid in their
respective contexts:

* The production sampler is a Monte-Carlo simulation with evaporation,
  feeding initial conditions for actual MD runs.
* The analytical reproduction is a closed-form didactic formula with no
  evaporation, used only for thesis-figure regression testing.

See ``migration_log.md`` for the full debugging story (this took
four wrong attempts before settling on a literal MATLAB transliteration
that matches the figure pixel-for-pixel).

The reproduction is a regression target -- the regression tests in
``test_droplet_sizes_diagnostics.py`` lock in the specific peak-height
ordering at each temperature:

* ``T=12 K``: dashed (reduced) peak HIGHER than solid (normal) peak
* ``T=18 K``: dashed peak much LOWER than solid peak (< 40%)
* ``T=18 K``: solid peak position in [1500, 4000] atoms

These signatures are visible in the thesis figure and would have caught
my earlier wrong attempts immediately if I had encoded them as tests
from the start.

## Use cases

- **Reproducing thesis figure 3.2**: call :func:`plot_thesis_figure_3_2`
  with no arguments. Compare side-by-side to the thesis figure.
- **Investigating why a parameter sweep fails**: the per-round table from
  :func:`diagnose_pickup` shows whether the issue is low pickup
  probability (small droplets), evaporation destroying droplets, or a
  mismatch between the formula ``<N>`` and what the pickup process
  actually selects.
- **Verifying physics changes**: when modifying any of the constants
  (``E_solv_meV``, ``_MEAN_FREE_PATH_M``,
  ``_REL_ENERGY_LOSS_PER_COLLISION``, ``_PICKUP_REGION_LENGTH_M``),
  running diagnostics at known conditions gives a quick visual sanity
  check.
