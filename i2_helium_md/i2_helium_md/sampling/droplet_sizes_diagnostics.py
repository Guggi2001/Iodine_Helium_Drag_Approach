"""Diagnostic plots and reports for droplet-size sampling.

These functions are debugging aids -- they re-run the sampler with
instrumentation enabled and produce histograms / console reports for
understanding what's happening inside the pickup simulation.

They are kept in a separate module so the production sampler in
``droplet_sizes.py`` stays free of plotting code, matching the project
principle that samplers should not import matplotlib.

Replaces the ``debug_pickup_plot`` flag from the legacy MATLAB code
(``generate_droplet_sizes.m``).

Two kinds of diagnostics are exposed:

* :func:`diagnose_pickup`              -- runs our Monte-Carlo
  pickup simulator (``_simulate_pickup``) and reports per-round physics
  + plots the raw vs one-pickup distribution.
* :func:`plot_thesis_figure_3_2`       -- reproduces the thesis figure 3.2
  using the closed-form analytical formula
  (``conditional_size_distributions_analytical``). This is a regression
  target: the figure should look identical to the thesis if our analytical
  code is correct.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..config import SimConfig
from .droplet_sizes import (
    _DEFAULT_E_SOLV_meV,
    _DEFAULT_NOZZLE_DIAMETER_UM,
    _LOGNORMAL_DELTA,
    _simulate_pickup,
    conditional_size_distributions_analytical,
    mean_droplet_size,
)

if TYPE_CHECKING:
    import matplotlib.figure


def diagnose_pickup(
    cfg: SimConfig,
    *,
    reduced_crosssection: bool = False,
    nozzle_diameter_um: float = _DEFAULT_NOZZLE_DIAMETER_UM,
    E_solv_meV: float = _DEFAULT_E_SOLV_meV,
    oversample_factor: int = 10,
    rng: np.random.Generator | None = None,
    plot: bool = True,
    print_report: bool = True,
) -> tuple[dict, "matplotlib.figure.Figure | None"]:
    """Run the pickup simulation with instrumentation and return diagnostics.

    Equivalent to running ``sample_droplet_sizes`` with the legacy MATLAB
    ``debug_pickup_plot=true`` flag enabled. Produces:

    * A printed per-round report (mean pickup probability, mean evaporated
      atoms, number of droplets picked up vs destroyed each round).
    * A 4-panel matplotlib figure overlaying the raw, no-pickup, and
      one-pickup distributions plus per-round statistics (optional).

    Parameters
    ----------
    cfg : SimConfig
        Used for ``num_molecules``, ``p_source_mbar``, ``T_source_K``, ``seed``.
    reduced_crosssection : bool, optional
        Whether to use the kinetic-energy-dissipation-reduced cross section.
    nozzle_diameter_um, E_solv_meV : float, optional
        Override the corresponding sampler defaults.
    oversample_factor : int, optional
        How many raw samples to generate per requested molecule. Default 10.
    rng : np.random.Generator, optional
        Reproducible RNG. If None, build one from ``cfg.seed``.
    plot : bool, optional
        If True (default), build a matplotlib figure. If False, return None
        for the figure.
    print_report : bool, optional
        If True (default), print a per-round summary table.

    Returns
    -------
    diagnostics : dict
        Returned by :func:`_simulate_pickup` with ``return_diagnostics=True``.
        Keys:
        - ``"N_raw"``           -- shape (oversample,)
        - ``"N_one_pickup"``    -- all droplets with exactly one pickup
        - ``"N_no_pickup"``     -- droplets that never picked anything up
        - ``"per_round"``       -- list of dicts with per-iteration stats
    fig : matplotlib.figure.Figure or None
        4-panel figure if ``plot=True``, else None.
    """
    if rng is None:
        rng = np.random.default_rng(cfg.seed)

    N_mean, _ = mean_droplet_size(
        cfg.p_source_mbar, cfg.T_source_K, nozzle_diameter_um,
    )
    mu = np.log(N_mean) - _LOGNORMAL_DELTA ** 2 / 2.0

    oversample = max(cfg.num_molecules * oversample_factor, 10_000)
    N_raw = rng.lognormal(mean=mu, sigma=_LOGNORMAL_DELTA, size=oversample)

    _, diagnostics = _simulate_pickup(
        N_raw,
        num_target=min(cfg.num_molecules, len(N_raw)),
        reduced_crosssection=reduced_crosssection,
        rng=rng,
        E_solv_meV=E_solv_meV,
        return_diagnostics=True,
    )

    # also stash the source conditions for reproducibility
    diagnostics["source_conditions"] = {
        "p_source_mbar": cfg.p_source_mbar,
        "T_source_K": cfg.T_source_K,
        "nozzle_diameter_um": nozzle_diameter_um,
        "E_solv_meV": E_solv_meV,
        "reduced_crosssection": reduced_crosssection,
        "N_mean_formula": N_mean,
    }

    if print_report:
        _print_pickup_report(diagnostics)

    fig = _make_pickup_figure(diagnostics) if plot else None

    return diagnostics, fig


# ===========================================================================
# Internal helpers
# ===========================================================================
def _print_pickup_report(diag: dict) -> None:
    """Print a per-round summary table to stdout."""
    src = diag["source_conditions"]
    print(
        f"\nPickup-cell simulation diagnostics\n"
        f"  source: p={src['p_source_mbar']} mbar, "
        f"T={src['T_source_K']} K, d={src['nozzle_diameter_um']} um\n"
        f"  E_solv = {src['E_solv_meV']} meV, "
        f"reduced sigma = {src['reduced_crosssection']}\n"
        f"  formula <N> = {src['N_mean_formula']:.0f}\n"
    )
    print(f"  raw samples : {diag['N_raw'].size:>8d}, mean N = {diag['N_raw'].mean():.0f}")
    if diag['N_one_pickup'].size > 0:
        print(f"  one-pickup  : {diag['N_one_pickup'].size:>8d}, mean N = {diag['N_one_pickup'].mean():.0f}")
    else:
        print(f"  one-pickup  : 0  (no droplets ended with exactly one pickup)")
    print(f"  no-pickup   : {diag['N_no_pickup'].size:>8d}")
    print()
    print(f"  {'round':>5} | {'n_alive':>9} | {'<p_pickup>':>10} | {'<dN evap>':>10} | "
          f"{'picked':>7} | {'destroyed':>9}")
    print("  " + "-" * 65)
    for r in diag["per_round"]:
        print(
            f"  {r['iter']:>5} | {r['n_alive_before']:>9} | {r['mean_p_pickup']:>10.4f} | "
            f"{r['mean_evap']:>10.0f} | {r['n_picked_up']:>7} | "
            f"{r.get('n_destroyed', 0):>9}"
        )
    print()


def _make_pickup_figure(diag: dict):
    """Build the 4-panel debugging figure."""
    import matplotlib.pyplot as plt   # local import: avoid hard dep at module load

    src = diag["source_conditions"]
    N_raw = diag["N_raw"]
    N_one = diag["N_one_pickup"]
    N_no = diag["N_no_pickup"]

    bin_max = float(np.percentile(N_raw, 99.5))
    bins = np.linspace(0, bin_max, 80)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Panel 1: raw histogram
    ax = axes[0, 0]
    ax.hist(N_raw, bins=bins, density=True, alpha=0.7, color="C0",
            label=f"raw  (mean={N_raw.mean():.0f})")
    ax.axvline(src["N_mean_formula"], color="k", ls="--",
               label=f"<N> formula = {src['N_mean_formula']:.0f}")
    ax.set_xlabel("droplet size N")
    ax.set_ylabel("p(N)")
    ax.set_title("Raw log-normal distribution")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # Panel 2: raw vs one-pickup -- raw COUNTS like MATLAB's histogram() default,
    # so the visual height of "one pickup" reflects how many droplets actually
    # ended up with one pickup (smaller than the raw set, biased toward large N).
    # This is the diagnostic that matches generate_droplet_sizes.m line 226-229.
    ax = axes[0, 1]
    ax.hist(N_raw, bins=bins, density=False, alpha=0.5, color="C0",
            label=f"initial droplet sizes (n={N_raw.size})")
    if N_one.size > 0:
        ax.hist(N_one, bins=bins, density=False, alpha=0.6, color="C1",
                label=f"after 1 pickup (n={N_one.size}, mean={N_one.mean():.0f})")
    ax.set_xlabel("droplet size N")
    ax.set_ylabel("count")
    title_extra = " (reduced σ)" if src["reduced_crosssection"] else " (normal σ)"
    ax.set_title("Initial vs one-pickup distribution (counts)" + title_extra)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # Panel 3: yield breakdown across rounds
    ax = axes[1, 0]
    rounds = [r["iter"] for r in diag["per_round"]]
    n_alive = [r["n_alive_before"] for r in diag["per_round"]]
    n_picked = [r["n_picked_up"] for r in diag["per_round"]]
    n_destroyed = [r.get("n_destroyed", 0) for r in diag["per_round"]]
    width = 0.4
    x = np.arange(len(rounds))
    ax.bar(x - width / 2, n_alive, width, label="alive", color="C0")
    ax.bar(x + width / 2, n_picked, width, label="picked up", color="C3")
    if any(n_destroyed):
        ax.bar(x + width / 2, n_destroyed, width, bottom=n_picked,
                label="destroyed", color="C1", alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(rounds)
    ax.set_xlabel("round")
    ax.set_ylabel("droplet count")
    ax.set_title("Per-round yield")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y")

    # Panel 4: per-round new_samples density (Bayesian-smoothed histogram).
    # This is a direct port of MATLAB's debug_pickup_plot inner-loop figure
    # (generate_droplet_sizes.m lines 191-196), which uses bayes_hist() to
    # plot the distribution of droplets that picked up at least one molecule
    # in the current round, after evaporation. For max_pickup=1 there's only
    # one round, so this shows a single curve. For multi-pickup runs the
    # panel would show one line per round, color-coded.
    ax = axes[1, 1]
    N_max_for_bins = float(N_raw.max())
    panel4_bins = np.linspace(1, N_max_for_bins, 60)
    cmap = plt.get_cmap("viridis")
    n_rounds = len(diag["per_round"])
    for round_info in diag["per_round"]:
        ns = round_info.get("new_samples")
        if ns is None or ns.size == 0:
            continue
        centers, h, sigma_h = _bayes_hist(ns, panel4_bins)
        # color by round index, matching MATLAB's bar_colors cycling
        c = cmap((round_info["iter"] - 1) / max(n_rounds, 1))
        ax.plot(centers, h, "-", lw=1.2, color=c,
                label=f"{round_info['iter']} pickup round")
        # error band like MATLAB's errorbar(centers, h, sigma_h)
        ax.fill_between(centers, h - sigma_h, h + sigma_h,
                        color=c, alpha=0.2, linewidth=0)
    ax.set_xlabel("N")
    ax.set_ylabel("p(N)")
    ax.set_title("Distributions after pickup events")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    title_str = (
        f"Pickup diagnostics: p={src['p_source_mbar']} mbar, "
        f"T={src['T_source_K']} K, d={src['nozzle_diameter_um']} um, "
        f"E_solv={src['E_solv_meV']} meV"
        + (", reduced σ" if src["reduced_crosssection"] else ", normal σ")
    )
    fig.suptitle(title_str)
    fig.tight_layout()
    return fig


def _bayes_hist(
    samples: np.ndarray,
    bin_edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bayesian-smoothed histogram density (port of MATLAB ``bayes_hist.m``).

    Uses Laplace smoothing to avoid zero-probability bins:

        p_i      = (n_i + 1) / (N + nbins + 1)
        sigma_p2 = (n_i + 2) / (N + nbins + 2) * p_i  -  p_i^2

    where ``n_i`` is the raw count in bin ``i``, ``N`` is the total sample
    count and ``nbins`` is the number of bins. Returns ``(centers, h,
    sigma_h)`` with ``h = p / binwidth`` and ``sigma_h = sigma_p / binwidth``.

    This is what MATLAB's diagnostic plot uses for the per-round new_samples
    distribution; we keep it byte-equivalent for verification.
    """
    n_counts, edges = np.histogram(samples, bins=bin_edges)
    nbins = n_counts.size
    binwidth = edges[1] - edges[0]
    centers = edges[:-1] + binwidth / 2.0
    N = samples.size

    p = (n_counts + 1.0) / (N + nbins + 1.0)
    sigma_p2 = (n_counts + 2.0) / (N + nbins + 2.0) * p - p ** 2
    sigma_p2 = np.maximum(sigma_p2, 0.0)
    sigma_p = np.sqrt(sigma_p2)

    h = p / binwidth
    sigma_h = sigma_p / binwidth
    return centers, h, sigma_h


# ===========================================================================
# Thesis figure 3.2 reproduction (analytical)
# ===========================================================================
def plot_thesis_figure_3_2(
    *,
    p_source_mbar: float = 40.0,
    nozzle_diameter_um: float = 5.0,
    temperatures_K: tuple[float, ...] = (12.0, 15.0, 18.0),
    N_max: int = 100_000,
):
    """Reproduce thesis figure 3.2 using the analytical formula.

    Returns a matplotlib Figure with two panels matching the thesis:

    * (a) Mean conditional droplet size vs source temperature, for both
      normal and reduced cross sections.
    * (b) Conditional droplet-size distributions at the requested
      temperatures, normal and reduced.

    Uses the parameters from
    ``conditional_droplet_size_distribution_simplified.m`` (the didactic
    script that produced the thesis figure). All these parameters live
    in the underlying :func:`conditional_size_distributions_analytical`
    with their thesis-script defaults; we just pass through here.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    N_grid = np.arange(1, N_max + 1, dtype=float)
    common = dict(
        p_source_mbar=p_source_mbar,
        nozzle_diameter_um=nozzle_diameter_um,
    )

    # purple -> blue -> pink palette like the thesis
    colors = {
        T: c for T, c in zip(
            sorted(temperatures_K),
            ["#1f3a8a", "#7a2682", "#d72660"],
        )
    }

    fig, (ax_top, ax_main) = plt.subplots(
        2, 1, figsize=(9, 6), sharex=True,
        gridspec_kw={"height_ratios": [1, 4]},
    )

    # -- panel a: mean N vs T --
    T_grid = np.linspace(min(temperatures_K), max(temperatures_K), 25)
    mean_normal = np.zeros_like(T_grid)
    mean_reduced = np.zeros_like(T_grid)
    for i, T in enumerate(T_grid):
        p_n, p_r = conditional_size_distributions_analytical(
            N_grid, T_source_K=T, **common,
        )
        # Normal: <N> = trapz(N * p_normal)  (p_normal integrates to 1)
        mean_normal[i] = np.trapezoid(N_grid * p_n, N_grid)
        # Reduced: same shape as MATLAB line ~95 -- mean weighted by the
        # *renormalised* reduced distribution. We renormalise locally.
        Z_r = np.trapezoid(p_r, N_grid)
        if Z_r > 0:
            mean_reduced[i] = np.trapezoid(N_grid * p_r, N_grid) / Z_r
        else:
            mean_reduced[i] = np.nan
    ax_top.plot(mean_normal, T_grid, "k-", lw=1.6, label="normal σ")
    ax_top.plot(mean_reduced, T_grid, "k--", lw=1.6, label="reduced σ")
    ax_top.set_ylabel("T / K")
    ax_top.set_xlim(0, 30000)
    ax_top.set_ylim(min(temperatures_K), max(temperatures_K))
    ax_top.legend(loc="upper right", fontsize=8)
    ax_top.grid(alpha=0.3)

    # -- panel b: distributions --
    for T in sorted(temperatures_K):
        p_n, p_r = conditional_size_distributions_analytical(
            N_grid, T_source_K=T, **common,
        )
        ax_main.plot(N_grid, p_n, "-", color=colors[T], lw=1.6,
                      label=f"T = {T:.0f} K")
        ax_main.plot(N_grid, p_r, "--", color=colors[T], lw=1.6,
                      label=f"T = {T:.0f} K, reduced σ")
        # tick marks at the means: normal uses simple trapz, reduced uses
        # the locally-renormalised value (matches MATLAB line 95).
        Nm_n = np.trapezoid(N_grid * p_n, N_grid)
        Z_r = np.trapezoid(p_r, N_grid)
        Nm_r = np.trapezoid(N_grid * p_r, N_grid) / Z_r if Z_r > 0 else np.nan
        ax_main.plot([Nm_n, Nm_n], [2.45e-4, 2.55e-4],
                      color=colors[T], lw=1.6)
        if not np.isnan(Nm_r):
            ax_main.plot([Nm_r, Nm_r], [2.45e-4, 2.55e-4],
                          color=colors[T], lw=1.6, ls="--")

    ax_main.set_xlabel("droplet size N")
    ax_main.set_ylabel("probability density")
    ax_main.set_xlim(0, 30000)
    ax_main.set_ylim(0, 2.7e-4)
    ax_main.legend(loc="upper right", fontsize=8)
    ax_main.grid(alpha=0.3)

    fig.suptitle(
        f"Conditional size distribution P(N | k=1)  "
        f"(p={p_source_mbar} mbar, d={nozzle_diameter_um} μm, "
        f"thesis-script defaults)"
    )
    fig.tight_layout()
    return fig
