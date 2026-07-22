"""Shared section builders for the run-summary family of scripts.

The detector-facing figure recipes consolidated out of
``plot_run_summary.py`` (2026-07-22, `TIER2_DETECTION_SUMMARY_SPEC.md`
§2.1) so that ``plot_run_summary.py`` (MD-window document) and
``plot_detection_summary.py`` (detection-stage document) render every
panel through the **same** implementation — rule 1, no duplicated
recipes. All numeric recipes are verbatim moves; behavior on legacy run
dirs is unchanged.

The mass-gated builders accept any object exposing the checkpoint gate
surface (``mass_final_kg``, ``velocities_final_{x,y,z}``,
``b_ion_outside``, ``num_molecules``) — an ``IonCheckpoint`` or a
``DetectedEnsembleView``.

Window-only sections (metadata, energy balances, temperature, radial
evolution, pair distance, Boltzmann, HeDFT comparisons) stay in
``plot_run_summary.py``: they are about the MD trajectory, not the
detector.
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    compute_final_velocity_histogram,
    fragment_gate_counts,
    fragment_mean_kinetic_energy,
    load_ihe_ked_curve,
    mass_spectrum,
    paper_v2_velocity_curve,
    paper_v2_velocity_map,
    paper_v4_angular_pair_covariance,
    radial_pair_speed_covariance,
    speed_mps_of_energy_eV,
)
from i2_helium_md.physics.shell_schedule import complex_mass_amu  # noqa: E402
from i2_helium_md.postprocess._smoothing import (  # noqa: E402
    moving_mean,
    normalise_trace,
)
from i2_helium_md.postprocess.paper_v2_plotting import (  # noqa: E402
    build_polar_image_figure,
    build_vmi_figure,
    load_optional_image,
    load_optional_polar_image,
    polar_histogram_matched_to_reference,
)
from i2_helium_md.postprocess.paper_cov_plotting import (  # noqa: E402
    build_angular_cov_figure as build_paper_cov_angular_cov_figure,
    build_pair_cov_traces_figure as build_paper_cov_pair_cov_traces_figure,
    build_phi_distribution_figure as build_paper_cov_phi_distribution_figure,
    build_radial_cov_figure as build_paper_cov_radial_cov_figure,
    build_radial_distribution_figure as build_paper_cov_radial_distribution_figure,
    load_optional_cov_reference,
    load_optional_high_snr_radial,
    load_optional_phi_reference as load_optional_paper_cov_phi_reference,
)
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    load_confirmation_run,
    score_histogram_vs_reference,
    score_ke_curve_vs_reference,
    solvated_renormalized,
)


# =============================================================================
# Plot-tuning constants (kept aligned with the existing focused scripts)
# =============================================================================
MASS_I = 127.0
MASS_I_HE_AMU = 131.0
MASS_I_HE2_AMU = 135.0
MASS_SPECTRUM_MAX_AMU = 127.0 + 5.0 * 4.0

HIST_BIN_WIDTH_APS = 0.04
HIST_EDGE_MAX_APS = 26.0
HIST_NUM_BINS = int(round(HIST_EDGE_MAX_APS / HIST_BIN_WIDTH_APS))
HIST_SMOOTHING_WINDOW = 15
VELOCITY_PLOT_V_MAX_APS = 28.0

# ihe_ked overlay panels: cover the full detector range (crop edge
# 354 px = 3052 m/s at m = 127; reference README).
IHE_KED_HIST_V_MAX_APS = 31.0
IHE_KED_HIST_NUM_BINS = int(round(IHE_KED_HIST_V_MAX_APS / HIST_BIN_WIDTH_APS))
IHE_KED_PLOT_V_MAX_MPS = 3100.0

# Display smoothing for the experimental curve overlays: the exporter's own
# mask/normalization window ("movmean": 10, omitnan, in
# IHe_KED_reference.provenance.json). Using the same window means the
# smoothed envelope peaks at ~1 by construction of the reference
# normalization.
IHE_KED_EXP_SMOOTHING_WINDOW = 10

# Mean-energy figure: left panel covers the trusted-curve fragments
# (n < split), right panel the evaporative tail (n >= split), each on
# its own linear scale -- log-y compressed the low-n comparison where
# the reference is strongest.
IHE_KED_MEAN_ENERGY_SPLIT_N = 5

# Detection-stage sections (Tier-2 confirmation runs only; adjudication D4
# option (a), log entry "PRODUCTION POINT ADJUDICATED" 2026-07-21): scorer
# kwargs mirror scripts/post_processing/tier2_confirmation_score.py — the
# solvated branch n >= 1, thin bins need >= 2 ions, sim-side SE widens the
# per-point sigma, and the n = 1 anchor is the I77 operative "median"
# (the mean-legacy chi^2 is annotated alongside, never substituted).
DETECTED_KE_N_MIN = 1
DETECTED_KE_MIN_BIN_COUNT = 2
DETECTED_KE_INCLUDE_SIM_SE = True
DETECTED_KE_N1_ANCHOR = "median"


class _SectionSkipped(RuntimeError):
    pass


def _section_mass_spectrum(ion, abundance) -> plt.Figure:
    """Final ion mass spectrum; side-by-side with the experimental
    I+He_n abundance when the reference is configured.

    Sim fractions use the same mass gates as the velocity diagnostics
    (m(n) +/- 0.5 amu, outside ions only), so a detected-fraction bar is
    directly the mass-gated ensemble the other ihe_ked sections draw from.
    """
    if abundance is None:
        spec = mass_spectrum(ion, bin_width_amu=1.0)
        fig, ax = plt.subplots(figsize=(8.0, 4.0), constrained_layout=True)
        ax.bar(spec.bin_centers_amu, spec.counts, width=0.9,
               edgecolor="black", linewidth=0.5)
        ax.set(title="Final ion mass spectrum",
               xlabel="m / u", ylabel="count")
        ax.set_xlim(left=MASS_I - 1, right=MASS_SPECTRUM_MAX_AMU + 1)
        ax.set_xticks(np.arange(MASS_I, MASS_SPECTRUM_MAX_AMU + 1, 4))
        return fig

    counts = fragment_gate_counts(ion, abundance.n)
    total = counts.sum()
    if total == 0:
        raise _SectionSkipped("no outside ions in any I+He_n mass gate")
    sim_percent = 100.0 * counts / total
    exp_percent = 100.0 * abundance.ion_fraction

    fig, ax = plt.subplots(figsize=(9.5, 4.5), constrained_layout=True)
    width = 0.4
    ax.bar(abundance.n - width / 2, exp_percent, width=width,
           color="tab:blue", label="experiment (ionPercent)")
    ax.bar(abundance.n + width / 2, sim_percent, width=width,
           color="tab:red", label=f"simulation (N={int(total)} atoms)")
    ax.set(title="Final I$^+$He$_n$ size distribution vs experimental abundance",
           xlabel="n (attached He atoms)", ylabel="fraction / %")
    ax.set_xticks(abundance.n)
    ax.legend(frameon=False)
    return fig


def _section_ihe_ked_mean_energy(ion, ked_ref) -> plt.Figure:
    """Mean kinetic energy per fragment: sim vs experiment, mean-to-mean.

    Error model per the reference README: per-point error is
    sqrt(stat^2 + sys^2); the calibration and condition bands are
    correlated (they shift the whole experimental curve coherently) and
    are drawn as shaded envelopes, never folded into point errors. Gold
    points are calib-limited: a disagreement there is real physics beyond
    the two bands.

    Rendered as two side-by-side linear-scale panels split at
    ``IHE_KED_MEAN_ENERGY_SPLIT_N``: a single log-y panel compresses the
    visual distance exactly where absolute sim-vs-exp discrepancies matter
    most (low n, the trusted-curve / gold-point regime) and flatters
    near-zero sim values at high n. The left panel covers the
    trusted-curve fragments (n < split); the right panel the evaporative
    tail (n >= split); each on its own linear scale.
    """
    if ked_ref is None:
        raise _SectionSkipped("IHE_KED_REFERENCE_DIR is None")

    sim_points = []
    for n in ked_ref.n:
        try:
            sim_points.append(fragment_mean_kinetic_energy(ion, int(n)))
        except ValueError:
            continue

    split = IHE_KED_MEAN_ENERGY_SPLIT_N
    mean = ked_ref.mean_KE_eV
    gold = ked_ref.gold_mask
    n_max = int(ked_ref.n.max())

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), constrained_layout=True)
    panel_specs = (
        (axes[0], ked_ref.n < split, lambda n: n < split,
         f"n = 0-{split - 1} (trusted-curve fragments)", True),
        (axes[1], ked_ref.n >= split, lambda n: n >= split,
         f"n = {split}-{n_max}", False),
    )

    for ax, subset, in_panel, panel_title, show_legend in panel_specs:
        if not np.any(subset):
            ax.set_axis_off()
            continue

        for frac, label, alpha in (
            (ked_ref.calib_syst_frac, "calibration band (correlated)", 0.20),
            (ked_ref.condition_syst_frac, "condition band (correlated)", 0.12),
        ):
            ax.fill_between(
                ked_ref.n[subset],
                mean[subset] * (1.0 - frac[subset]),
                mean[subset] * (1.0 + frac[subset]),
                color="tab:blue", alpha=alpha, linewidth=0, label=label,
            )
        ax.errorbar(
            ked_ref.n[subset], mean[subset],
            yerr=ked_ref.point_err_eV[subset], fmt="o",
            color="tab:blue", markersize=4, capsize=2,
            label=r"experiment $\langle E\rangle$ (stat $\oplus$ sys)")
        gold_subset = gold & subset
        ax.plot(ked_ref.n[gold_subset], mean[gold_subset], "o", markersize=10,
                markerfacecolor="none", markeredgecolor="goldenrod",
                markeredgewidth=1.5, label="gold points (calib-limited)")

        panel_sim = [p for p in sim_points if in_panel(p.n)]
        if panel_sim:
            ax.errorbar(
                [p.n for p in panel_sim],
                [p.mean_KE_eV for p in panel_sim],
                yerr=[p.stat_err_mean_KE_eV for p in panel_sim],
                fmt="s", linestyle="none", color="tab:red", markersize=5,
                capsize=2, label=r"simulation $\langle E\rangle$",
            )
            for p in panel_sim:
                ax.annotate(f"N={p.num_atoms_used}",
                            (p.n, p.mean_KE_eV),
                            textcoords="offset points", xytext=(0, 7),
                            fontsize=7, color="tab:red", ha="center")

        ax.set_ylim(bottom=0.0)
        ax.set(title=panel_title,
               xlabel="n (attached He atoms)",
               ylabel=r"$\langle E\rangle$ / eV")
        ax.set_xticks(ked_ref.n[subset])
        if show_legend:
            ax.legend(frameon=False, fontsize=8)

    fig.suptitle(r"I$^+$He$_n$ mean kinetic energy (mean-to-mean)")
    return fig


def _nan_aware_moving_mean(values: np.ndarray, window: int) -> np.ndarray:
    """Centered moving mean with MATLAB ``movmean(..., 'omitnan')`` semantics
    that additionally preserves NaN positions.

    Each finite entry becomes the mean of the finite entries inside its
    window; entries that are NaN in the input stay NaN in the output, so the
    reference's deliberate NaN cuts (low-E Abel exclusions in the 3-D
    columns) render as gaps and smoothing never bleeds across them.

    Parameters
    ----------
    values : (N,) float array, may contain NaN
    window : positive int, boxcar length in samples

    Returns
    -------
    (N,) float array with NaN exactly where ``values`` has NaN.
    """
    if window < 1:
        raise ValueError(f"window must be >= 1, got {window}")
    v = np.asarray(values, dtype=float)
    finite = np.isfinite(v)
    kernel = np.ones(int(window))
    num = np.convolve(np.where(finite, v, 0.0), kernel, mode="same")
    den = np.convolve(finite.astype(float), kernel, mode="same")
    out = np.divide(num, den, out=np.full_like(v, np.nan), where=den > 0.0)
    out[~finite] = np.nan
    return out


def _section_ihe_ked_curves(ion, ked_dir, ked_ref, representation) -> plt.Figure:
    """Per-fragment speed-distribution overlays for n = 0..3.

    n = 4 is dropped from this display: it is the weakest trusted curve
    (N_eff ~= 11.7k vs 193.7k at n = 0) and forced a 2x3 grid with a dead
    sixth axis. The n = 4 export stays on disk and loadable
    (``load_ihe_ked_curve`` still supports n up to ``CURVE_N_MAX`` = 4) --
    only this figure's panel count changed.

    representation: "3d" overlays the sim 3-D |v| histogram on the
    reconstructed ``signal_3d_Pv``; "2d" overlays the sim in-plane
    projected speed sqrt(vx^2+vy^2) on the detector projection
    ``signal_2d_Pv``. Reference curves are peak-normalized (smoothed
    envelope = 1): compare shapes, never amplitudes. Vertical markers sit
    at v(<E>) -- not <v> -- on both sides; the reference NaN cuts render
    as gaps. The experimental curve itself is drawn as faint raw points
    plus a movmean-smoothed line (``IHE_KED_EXP_SMOOTHING_WINDOW``) so the
    underlying scatter stays visible under the smoothed envelope.
    """
    if ked_dir is None or ked_ref is None:
        raise _SectionSkipped("IHE_KED_REFERENCE_DIR is None")
    if representation not in ("2d", "3d"):
        raise ValueError(f"representation must be '2d' or '3d', "
                         f"got {representation!r}")

    titles = {
        "3d": ("3-D speed distributions P(v) vs I$^+$He$_n$ reference "
               "(peak-normalized, shapes only)"),
        "2d": ("2-D detector projections (in-plane speed) vs I$^+$He$_n$ "
               "reference (peak-normalized, shapes only)"),
    }
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.0),
                             constrained_layout=True)
    axes = axes.ravel()
    fig.suptitle(titles[representation])

    for n in range(4):
        ax = axes[n]
        curve = load_ihe_ked_curve(ked_dir, n)
        ref_signal = getattr(curve, f"signal_{representation}_Pv")
        # Raw export as faint points joined by a fine line; matplotlib
        # breaks the line at the NaN cuts, so gaps stay gaps.
        ax.plot(curve.v_mps, ref_signal, ".-", color="tab:blue",
                markersize=2.5, linewidth=0.5, alpha=0.30,
                label="experiment (raw)")
        ax.plot(curve.v_mps,
                _nan_aware_moving_mean(ref_signal, IHE_KED_EXP_SMOOTHING_WINDOW),
                color="tab:blue", linewidth=1.4,
                label="experiment (movmean 10)")

        mass_amu = float(complex_mass_amu(n))
        sim_note = None
        try:
            hist = compute_final_velocity_histogram(
                ion, mass_amu=mass_amu,
                num_bins=IHE_KED_HIST_NUM_BINS,
                v_max_Aps=IHE_KED_HIST_V_MAX_APS,
                projected=(representation == "2d"),
            )
        except ValueError:
            sim_note = "sim: no atoms in gate"
        else:
            sim_density = normalise_trace(
                moving_mean(hist.density, HIST_SMOOTHING_WINDOW)
            )
            ax.plot(hist.bin_centers_mps, sim_density, "--",
                    color="tab:red", linewidth=1.4,
                    label=f"simulation (N={hist.num_atoms_used})")

        # v(<E>) markers (labelled v(<E>), not <v>): experiment from the
        # reference table, simulation from the mass-gated ensemble.
        v_exp = speed_mps_of_energy_eV(
            float(ked_ref.mean_KE_eV[n]), mass_amu,
        )
        ax.axvline(v_exp, color="tab:blue", linestyle=":", linewidth=1.0,
                   label=r"exp $v(\langle E\rangle)$")
        try:
            sim_mean = fragment_mean_kinetic_energy(ion, n)
        except ValueError:
            pass
        else:
            ax.axvline(sim_mean.v_of_mean_E_mps, color="tab:red",
                       linestyle=":", linewidth=1.0,
                       label=r"sim $v(\langle E\rangle)$")

        if sim_note:
            ax.annotate(sim_note, (0.97, 0.9), xycoords="axes fraction",
                        ha="right", fontsize=8, color="tab:red")
        ax.set(title=f"n = {n}", xlim=(0.0, IHE_KED_PLOT_V_MAX_MPS),
               ylim=(0.0, 1.25))
        ax.set_xlabel("v / m/s")
        ax.set_ylabel("signal / arb. units")
        if n == 0:
            ax.legend(frameon=False, fontsize=7)
    return fig


def _section_paper_v2_vmi(ion, reference_dir, noise_floor, *, mass_amu) -> plt.Figure:
    if reference_dir is None:
        raise _SectionSkipped("PAPER_V2_REFERENCE_DIR is None")
    image_ref = load_optional_image(reference_dir, log_prefix="[run_summary]")
    velocity_map = paper_v2_velocity_map(ion, mass_amu=mass_amu)
    return build_vmi_figure(
        image_ref=image_ref,
        velocity_map=velocity_map,
        experimental_noise_floor=noise_floor,
    )


def _section_paper_cov_radial_distribution(
    ion, reference_dir, paper_v2_reference_dir, *, mass_amu,
) -> plt.Figure:
    if reference_dir is None:
        raise _SectionSkipped("PAPER_COV_REFERENCE_DIR is None")
    high_snr_ref = None
    if paper_v2_reference_dir is not None:
        high_snr_ref = load_optional_high_snr_radial(
            paper_v2_reference_dir, log_prefix="[run_summary]",
        )
    cov_ref = load_optional_cov_reference(reference_dir, log_prefix="[run_summary]")
    try:
        sim_radial_curve = paper_v2_velocity_curve(
            ion, mass_amu=mass_amu,
        )
        sim_radial = radial_pair_speed_covariance(
            ion, mass_amu=mass_amu,
        )
    except ValueError as exc:
        raise _SectionSkipped(str(exc))
    return build_paper_cov_radial_distribution_figure(
        high_snr_ref=high_snr_ref,
        sim_radial_curve=sim_radial_curve,
        sim_radial=sim_radial,
        cov_ref=cov_ref,
        title = 'run_summary'
    )


def _section_paper_cov_phi_distribution(ion, reference_dir, *, mass_amu) -> plt.Figure:
    if reference_dir is None:
        raise _SectionSkipped("PAPER_COV_REFERENCE_DIR is None")
    phi_ref = load_optional_paper_cov_phi_reference(
        reference_dir, log_prefix="[run_summary]",
    )
    try:
        return build_paper_cov_phi_distribution_figure(
            phi_ref=phi_ref,
            ion=ion,
            mass_amu=mass_amu,
        )
    except ValueError as exc:
        raise _SectionSkipped(str(exc))


def _section_paper_cov_angular(ion, reference_dir, *, mass_amu) -> plt.Figure:
    cov_ref = _required_paper_cov_reference(reference_dir)
    try:
        sim_angular = paper_v4_angular_pair_covariance(
            ion, mass_amu=mass_amu,
        )
    except ValueError as exc:
        raise _SectionSkipped(str(exc))
    return build_paper_cov_angular_cov_figure(
        cov_ref=cov_ref, sim_angular=sim_angular,
    )


def _section_paper_cov_radial(ion, reference_dir, *, mass_amu) -> plt.Figure:
    cov_ref = _required_paper_cov_reference(reference_dir)
    try:
        sim_radial = radial_pair_speed_covariance(
            ion, mass_amu=mass_amu,
        )
    except ValueError as exc:
        raise _SectionSkipped(str(exc))
    return build_paper_cov_radial_cov_figure(
        cov_ref=cov_ref, sim_radial=sim_radial,
    )


def _section_paper_cov_traces(ion, reference_dir, *, mass_amu) -> plt.Figure:
    cov_ref = _required_paper_cov_reference(reference_dir)
    try:
        sim_angular = paper_v4_angular_pair_covariance(
            ion, mass_amu=mass_amu,
        )
        sim_radial = radial_pair_speed_covariance(
            ion, mass_amu=mass_amu,
        )
    except ValueError as exc:
        raise _SectionSkipped(str(exc))
    return build_paper_cov_pair_cov_traces_figure(
        cov_ref=cov_ref,
        sim_angular=sim_angular,
        sim_radial=sim_radial,
    )


def _required_paper_cov_reference(reference_dir):
    if reference_dir is None:
        raise _SectionSkipped("PAPER_COV_REFERENCE_DIR is None")
    cov_ref = load_optional_cov_reference(reference_dir, log_prefix="[run_summary]")
    if cov_ref is None:
        raise _SectionSkipped("no paper-cov experimental covariance reference found")
    return cov_ref


def _section_paper_v2_polar(ion, reference_dir, noise_floor, *, mass_amu) -> plt.Figure:
    if reference_dir is None:
        raise _SectionSkipped("PAPER_V2_REFERENCE_DIR is None")
    polar_ref = load_optional_polar_image(
        reference_dir, log_prefix="[run_summary]"
    )
    if polar_ref is None:
        raise _SectionSkipped("no polar VMI image reference found")
    polar_hist = polar_histogram_matched_to_reference(
        ion, polar_ref, mass_amu=mass_amu,
    )
    return build_polar_image_figure(
        polar_ref, polar_hist,
        experimental_noise_floor=noise_floor,
    )


def _section_mass_resolved(ion) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(9.0, 4.5), constrained_layout=True)
    drew_anything = False
    for mass, label, style in (
        (127.0, r"$I^+$", "-"),
        (MASS_I_HE_AMU, r"$I^+ He$", "--"),
        (MASS_I_HE2_AMU, r"$I^+ He_2$", "-."),
    ):
        try:
            h = compute_final_velocity_histogram(
                ion, mass_amu=mass,
                num_bins=HIST_NUM_BINS, v_max_Aps=HIST_EDGE_MAX_APS,
            )
        except ValueError:
            continue
        smoothed = normalise_trace(
            moving_mean(h.density, HIST_SMOOTHING_WINDOW)
        )
        ax.plot(h.bin_centers_Aps, smoothed, style,
                label=f"{label} (n={h.num_atoms_used})")
        drew_anything = True
    if not drew_anything:
        raise _SectionSkipped("no mass channel produced a histogram")
    ax.set(title="Mass-resolved final-velocity histograms",
           xlim=(0.0, VELOCITY_PLOT_V_MAX_APS), ylim=(0.0, 1.1),
           xlabel=r"v / $\mathrm{\AA}/\mathrm{ps}$",
           ylabel="signal / arb. units")
    ax.legend(frameon=False)
    return fig


# =============================================================================
# Detection-stage sections (Tier-2 confirmation runs; adjudication D4 (a))
# =============================================================================
def _load_detection_read(run_root: Path):
    """The detected read for a Tier-2 confirmation run dir, or ``None``.

    Gating contract: a legacy run dir (no ``detection.npz``) returns
    ``None`` and the summary renders exactly as before. All scoring
    conventions (retained-excluded, suppressed at n = 0) live in
    ``i2_helium_md.postprocess.tier2_confirmation``.
    """
    run_root = Path(run_root)
    if not (run_root / "detection.npz").exists():
        return None
    return load_confirmation_run(run_root)


def _section_detected_size_distribution(read, abundance) -> plt.Figure:
    """Detected I$^+$He$_n$ size distribution at t_detect (the Tier-2
    scored surface: ``droplet_retained`` excluded, suppressed ions in the
    n = 0 bare bin). Left: the full scored histogram. Right: the solvated
    (n >= 1 renormalized, RQ8) branch vs the experimental abundance
    reference with the W1/n1/ratio scores annotated; sim-only single
    panel when no reference is configured (mass-spectrum precedent)."""
    n_vals = read.n_values
    frac = read.fraction

    ncols = 2 if abundance is not None else 1
    fig, axes = plt.subplots(
        1, ncols, figsize=(6.0 * ncols + 0.5, 4.5), constrained_layout=True,
        squeeze=False,
    )
    ax = axes[0][0]
    ax.bar(n_vals, 100.0 * frac, width=0.8, color="tab:red",
           edgecolor="black", linewidth=0.4)
    ax.set(title="Detected size distribution (scored read)",
           xlabel="n (attached He atoms; suppressed at n = 0)",
           ylabel="scored fraction / %")
    fig.suptitle(f"Detected I$^+$He$_n$ at t_detect — {read.label}",
                 fontsize=11)
    ax.text(
        0.97, 0.95,
        f"scored {read.num_scored}/{read.num_ions} ions\n"
        f"retained excluded: {read.num_ions - read.num_scored}\n"
        f"suppressed (bare): {100.0 * read.suppressed_frac:.1f} %\n"
        rf"$\bar{{n}}_\mathrm{{det}}$ = {read.n_mean:.3f}",
        transform=ax.transAxes, ha="right", va="top", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    if abundance is None:
        return fig

    score = score_histogram_vs_reference(read, abundance)
    sim_n, sim_f = solvated_renormalized(n_vals, frac, n_min=score.n_min)
    ref_n, ref_f = solvated_renormalized(
        abundance.n, abundance.ion_fraction, n_min=score.n_min
    )
    ax = axes[0][1]
    width = 0.4
    ax.bar(ref_n - width / 2, 100.0 * ref_f, width=width, color="tab:blue",
           label="experiment (solvated)")
    ax.bar(sim_n + width / 2, 100.0 * sim_f, width=width, color="tab:red",
           label="simulation (solvated)")
    ratio = score.ratio_n1_over_n2
    if score.ref_n2 > 0.0:
        ref_ratio = score.ref_n1 / score.ref_n2
    else:
        # scorer convention (score_histogram_vs_reference): inf when only
        # the numerator carries mass, nan when both are empty
        ref_ratio = float("inf") if score.ref_n1 > 0.0 else float("nan")
    ax.text(
        0.97, 0.95,
        rf"$W_1$(solv) = {score.w1_solvated:.3f} bins" "\n"
        rf"$n_1$ = {score.sim_n1:.3f} (ref {score.ref_n1:.3f})" "\n"
        rf"$n_1/n_2$ = {ratio:.2f} (ref "
        rf"{ref_ratio:.2f})",
        transform=ax.transAxes, ha="right", va="top", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )
    ax.set(title="Solvated branch (n ≥ 1 renormalized) vs abundance",
           xlabel="n (attached He atoms)", ylabel="fraction / %")
    ax.legend(frameon=False)
    return fig


def _detected_ke_scores(read, ked_ref):
    """The committed KE-curve score under the operative median anchor,
    plus the mean-legacy chi^2 (both always reported side by side)."""
    kwargs = dict(
        n_min=DETECTED_KE_N_MIN,
        min_count=DETECTED_KE_MIN_BIN_COUNT,
        include_sim_se=DETECTED_KE_INCLUDE_SIM_SE,
    )
    med = score_ke_curve_vs_reference(
        read, ked_ref, n1_anchor=DETECTED_KE_N1_ANCHOR, **kwargs
    )
    legacy = (
        med if DETECTED_KE_N1_ANCHOR == "mean"
        else score_ke_curve_vs_reference(
            read, ked_ref, n1_anchor="mean", **kwargs
        )
    )
    return med, legacy


def _section_detected_ked_mean_energy(read, ked_ref) -> plt.Figure:
    """Detected per-n mean KE vs the ihe_ked reference under the committed
    error model (per-point stat (+) sys; correlated bands as envelopes),
    with the n = 1 median anchor (I77) marked and both chi^2 columns
    annotated. Same split-panel layout as the ion-stage section."""
    if ked_ref is None:
        raise _SectionSkipped("IHE_KED_REFERENCE_DIR is None")

    med, legacy = _detected_ke_scores(read, ked_ref)
    ke = read.ke_by_n(min_count=DETECTED_KE_MIN_BIN_COUNT)

    split = IHE_KED_MEAN_ENERGY_SPLIT_N
    mean = ked_ref.mean_KE_eV
    n_max = int(ked_ref.n.max())

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0),
                             constrained_layout=True)
    panel_specs = (
        (axes[0], ked_ref.n < split,
         f"n = 1–{split - 1} (detected read)", True),
        (axes[1], ked_ref.n >= split,
         f"n = {split}–{n_max} (detected read)", False),
    )
    for ax, subset, panel_title, show_legend in panel_specs:
        subset = subset & (ked_ref.n >= DETECTED_KE_N_MIN)
        if not np.any(subset):
            ax.set_axis_off()
            continue
        for frac, label, alpha in (
            (ked_ref.calib_syst_frac, "calibration band (correlated)", 0.20),
            (ked_ref.condition_syst_frac, "condition band (correlated)", 0.12),
        ):
            ax.fill_between(
                ked_ref.n[subset],
                mean[subset] * (1.0 - frac[subset]),
                mean[subset] * (1.0 + frac[subset]),
                color="tab:blue", alpha=alpha, linewidth=0, label=label,
            )
        ax.errorbar(
            ked_ref.n[subset], mean[subset],
            yerr=ked_ref.point_err_eV[subset], fmt="o",
            color="tab:blue", markersize=4, capsize=2,
            label=r"experiment $\langle E\rangle$ (stat $\oplus$ sys)")
        anchor_hit = subset & (ked_ref.n == 1)
        if np.any(anchor_hit):
            ax.plot(ked_ref.n[anchor_hit], ked_ref.median_KE_eV[anchor_hit],
                    "D", markersize=7, markerfacecolor="none",
                    markeredgecolor="tab:purple", markeredgewidth=1.5,
                    label="n = 1 median anchor (I77, operative)")
        in_panel = (ke.n < split) if ax is axes[0] else (ke.n >= split)
        in_panel = in_panel & (ke.n >= DETECTED_KE_N_MIN)  # unscored bare bin stays off
        if np.any(in_panel):
            ax.errorbar(
                ke.n[in_panel], ke.mean_eV[in_panel],
                yerr=ke.se_eV[in_panel], fmt="s", linestyle="none",
                color="tab:red", markersize=5, capsize=2,
                label=r"detected simulation $\langle E\rangle$",
            )
        ax.set(title=panel_title, xlabel="n (attached He atoms)",
               ylabel="mean kinetic energy / eV")
        if show_legend:
            ax.legend(frameon=False, fontsize=8)
    axes[0].text(
        0.03, 0.05,
        rf"$\chi^2$_med = {med.chi2_profiled:.1f}"
        rf" (mean-legacy {legacy.chi2_profiled:.1f});"
        f" {med.n_points} bins",
        transform=axes[0].transAxes, ha="left", va="bottom", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )
    fig.suptitle(
        f"Detected mean KE vs ihe_ked reference ({read.label})", fontsize=11
    )
    return fig


def _section_detected_ke_anatomy(read, ked_ref) -> plt.Figure:
    """Per-bin profiled residuals z of the detected KE curve (the §4z
    anatomy read): positive = simulation hot vs reference, negative =
    cold; the two correlated bands are already profiled out."""
    if ked_ref is None:
        raise _SectionSkipped("IHE_KED_REFERENCE_DIR is None")

    med, legacy = _detected_ke_scores(read, ked_ref)
    z = med.z_profiled
    colors = ["tab:red" if v > 0 else "tab:blue" for v in z]

    fig, ax = plt.subplots(figsize=(9.0, 4.5), constrained_layout=True)
    ax.bar(med.n, z, width=0.8, color=colors, edgecolor="black",
           linewidth=0.4)
    ax.axhline(0.0, color="black", linewidth=0.8)
    for lvl in (-2.0, 2.0):
        ax.axhline(lvl, color="gray", linewidth=0.8, linestyle="--")
    ax.set(
        title="Detected KE per-bin anatomy (profiled residual z)",
        xlabel="n (attached He atoms)",
        ylabel="z (sim − ref, per-point $\\sigma$)",
    )
    ax.text(
        0.03, 0.05,
        rf"$\chi^2$_med = {med.chi2_profiled:.1f}"
        rf" (mean-legacy {legacy.chi2_profiled:.1f}), {med.n_points} bins"
        "\n"
        f"a_calib = {med.a_calib:.2f}, b_cond = {med.b_condition:.2f}",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )
    return fig
