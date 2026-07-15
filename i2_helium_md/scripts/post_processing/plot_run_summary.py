"""Consolidated post-processing summary for one finished run directory.

Produces every numerical diagnostic the legacy MATLAB post-processing
scripts produced *that is in scope per CLAUDE.md*, in a single
multi-page PDF plus per-figure PNGs. Sections that need optional
reference data (HeDFT trajectory, experimental VMI) are gated on the
matching ``*_REF_PATH`` user setting so each run type stays separable
(per CLAUDE.md "Keep these workflows separate.").

The following legacy scripts are consolidated here:

* ``vmi_sim_3d_neutral_propa_HeDFT_mimic.m``    -> neutral energy balance
* ``vmi_sim_3d_ion_propa.m``                   -> ion energy balance + temperature
* ``simulation_image_only_trajectories.m``     -> HeDFT R(t), v(t)
* ``post_process_single_pulse_paper_v3.m``     -> 1D and 2D polar VMI panels
* ``post_process_single_pulse_paper_IplusHe_comparison_cov.m``
                                                 -> covariance-paper overlays
* ``post_process_single_pulse_paper.m``        -> bimodal Gaussian fit
* ``post_process_single_pulse.m``              -> 2D (vx, vy) histogram
* ``post_process_compare_radial_distributions.m`` -> time-resolved radial,
                                                     interatomic distance,
                                                     Boltzmann reference
* ``compare_neutral_dynamics_to_HeDFT.m``      -> neutral cumtrapz r(t)

The legacy ``simulation_image.m`` velocity overlay against the
``vmi_summary`` CSVs was retired 2026-07-14 in favor of the frozen
``data/reference/ihe_ked/`` per-fragment reference (mean-KE table +
trusted 2-D/3-D curves for n = 0..4); see
``docs/superpowers/specs/2026-07-14-ihe-ked-run-summary-design.md``.

Out of scope (deferred per CLAUDE.md): Abel inversion, pump-probe,
effusive / gas-phase comparison, live-debug 3D animations.

Edit the ``USER SETTINGS`` block below and run the script (e.g. from
PyCharm)::

    python scripts/post_processing/plot_run_summary.py
"""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import sys
from types import SimpleNamespace


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    boltzmann_population,
    compare_distance,
    compare_neutral_to_hedft,
    compare_velocity_magnitude,
    compute_final_velocity_histogram,
    fragment_gate_counts,
    fragment_mean_kinetic_energy,
    interparticle_distance_histogram,
    ion_energy_totals,
    load_he_abundance_reference,
    load_hedft_trajectory,
    load_ihe_ked_curve,
    load_ihe_ked_reference,
    mass_spectrum,
    neutral_energy_totals,
    paper_v2_velocity_curve,
    paper_v2_velocity_map,
    paper_v4_angular_pair_covariance,
    radial_distribution_evolution,
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
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


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

PAIR_DIST_NUM_BINS = 100
TIME_HEATMAP_N_SLICES = 60
TIME_HEATMAP_N_R_BINS = 100


# =============================================================================
# USER SETTINGS -- edit these and run the script (e.g. from PyCharm)
# =============================================================================
# Path to the run directory holding cfg.json + neutral.npz + ion.npz.
RUN_DIR: Path = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet_long"

# Optional reference CSVs. Set to ``None`` to skip the matching section.
# Typical experimental-droplet configuration:
HEDFT_REF_PATH: Path | None = None

# Directory holding the frozen I+He_n kinetic-energy reference
# (IHe_KED_reference.csv + IHe_KED_curves_n{0..4}.csv). ``None`` skips
# the ihe_ked mean-energy and curve-overlay sections.
IHE_KED_REFERENCE_DIR: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "ihe_ked"
)

# Experimental I+He_n abundance CSV for the mass-spectrum side-by-side.
# ``None`` keeps the plain simulated mass spectrum.
ABUNDANCE_REF_PATH: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"
)

# Directory holding paper-v2 reference CSVs and images/. Set to ``None`` to
# skip the paper-v2 VMI and polar sections.
PAPER_V2_REFERENCE_DIR: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "paper_v2"
)

# Directory holding paper-cov covariance references. Set to ``None`` to skip
# the covariance-paper radial, phi, and pair-covariance replacement sections.
PAPER_COV_REFERENCE_DIR: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "paper_cov"
)

# Fraction of experimental panel max intensity below which pixels clip to
# background. Raise to suppress more low-level noise; 0 disables. Only
# experimental panels use this; simulated panels are unaffected.
EXPERIMENTAL_NOISE_FLOOR: float = 0.20

# Mass channel used for the paper-v2 simulated curves.
PAPER_V2_MASS_AMU: float = 131.0

# Typical 9 A HeDFT comparison configuration (uncomment and comment out the
# experimental block above to switch):
# RUN_DIR = PROJECT_ROOT / "data" / "runs" / "9A_hedft_comparison"
# HEDFT_REF_PATH = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"
# IHE_KED_REFERENCE_DIR = None
# ABUNDANCE_REF_PATH = None

# Output directory for the PDF and per-panel PNGs. ``None`` -> <RUN_DIR>/figures.
OUT_DIR: Path | None = None

# Show figures interactively after writing them. False is the right setting
# for PyCharm / headless / smoke-test runs.
SHOW_FIGURES: bool = False


def main() -> int:
    run_dir = Path(RUN_DIR)
    hedft_ref = Path(HEDFT_REF_PATH) if HEDFT_REF_PATH else None
    ihe_ked_dir = Path(IHE_KED_REFERENCE_DIR) if IHE_KED_REFERENCE_DIR else None
    abundance_ref_path = Path(ABUNDANCE_REF_PATH) if ABUNDANCE_REF_PATH else None
    paper_v2_ref_dir = (
        Path(PAPER_V2_REFERENCE_DIR) if PAPER_V2_REFERENCE_DIR else None
    )
    paper_cov_ref_dir = (
        Path(PAPER_COV_REFERENCE_DIR) if PAPER_COV_REFERENCE_DIR else None
    )
    args = SimpleNamespace(
        run_dir=run_dir,
        hedft_ref=hedft_ref,
        ihe_ked_dir=ihe_ked_dir,
        abundance_ref_path=abundance_ref_path,
        paper_v2_ref_dir=paper_v2_ref_dir,
        paper_cov_ref_dir=paper_cov_ref_dir,
    )

    run = RunDirectory(run_dir)

    out_dir = Path(OUT_DIR) if OUT_DIR else (run.root / "figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "run_summary.pdf"

    cfg = run.load_cfg() if run.has_cfg() else None
    neutral = run.load_neutral() if run.has_neutral() else None
    ion = run.load_ion() if run.has_ion() else None

    hedft = load_hedft_trajectory(hedft_ref) if hedft_ref else None
    ked_ref = (
        load_ihe_ked_reference(ihe_ked_dir / "IHe_KED_reference.csv")
        if ihe_ked_dir else None
    )
    abundance = (
        load_he_abundance_reference(abundance_ref_path)
        if abundance_ref_path else None
    )

    print(f"[run_summary] writing {pdf_path}")
    with PdfPages(pdf_path) as pdf:
        sections = [
            ("metadata", lambda: _section_metadata(cfg, ion, neutral, args)),
        ]
        if neutral is not None:
            sections += [
                ("neutral_energy_balance",
                 lambda: _section_neutral_energy(neutral)),
            ]
        if ion is not None:
            sections += [
                ("ion_energy_balance",
                 lambda: _section_ion_energy(ion)),
                ("ion_temperature_diagnostic",
                 lambda: _section_temperature(ion)),
                ("mass_spectrum",
                 lambda: _section_mass_spectrum(ion, abundance)),
                ("ihe_ked_mean_energy",
                 lambda: _section_ihe_ked_mean_energy(ion, ked_ref)),
                ("ihe_ked_curves_3d",
                 lambda: _section_ihe_ked_curves(
                     ion, ihe_ked_dir, ked_ref, "3d")),
                ("ihe_ked_curves_2d",
                 lambda: _section_ihe_ked_curves(
                     ion, ihe_ked_dir, ked_ref, "2d")),
                ("paper_v2_vmi_comparison",
                 lambda: _section_paper_v2_vmi(
                     ion, paper_v2_ref_dir, EXPERIMENTAL_NOISE_FLOOR)),
                ("paper_cov_radial_distribution",
                 lambda: _section_paper_cov_radial_distribution(
                     ion, paper_cov_ref_dir, paper_v2_ref_dir)),
                ("paper_cov_phi_distribution",
                 lambda: _section_paper_cov_phi_distribution(
                     ion, paper_cov_ref_dir)),
                ("paper_v2_polar_image_comparison",
                 lambda: _section_paper_v2_polar(
                     ion, paper_v2_ref_dir, EXPERIMENTAL_NOISE_FLOOR)),
                ("mass_resolved_velocities",
                 lambda: _section_mass_resolved(ion)),
                ("radial_evolution_heatmap",
                 lambda: _section_radial_evolution(ion)),
                ("interparticle_distance_histogram",
                 lambda: _section_pair_distance(ion)),
                ("paper_cov_angular_pair_cov",
                 lambda: _section_paper_cov_angular(ion, paper_cov_ref_dir)),
                ("paper_cov_radial_pair_cov",
                 lambda: _section_paper_cov_radial(ion, paper_cov_ref_dir)),
                ("paper_cov_pair_cov_traces",
                 lambda: _section_paper_cov_traces(ion, paper_cov_ref_dir)),
            ]
        if neutral is not None and hedft is not None:
            sections += [
                ("hedft_neutral_comparison",
                 lambda: _section_hedft_neutral(neutral, hedft)),
            ]
        if ion is not None and hedft is not None:
            sections += [
                ("hedft_ion_comparison",
                 lambda: _section_hedft_ion(ion, hedft)),
            ]
        if cfg is not None and ion is not None:
            sections += [
                ("boltzmann_overlay_initial",
                 lambda: _section_boltzmann(cfg, ion)),
            ]

        for label, builder in sections:
            try:
                fig = builder()
            except _SectionSkipped as skip:
                print(f"[run_summary] skip {label}: {skip}")
                continue
            if fig is None:
                continue
            pdf.savefig(fig)
            fig.savefig(out_dir / f"{label}.png", dpi=150)
            plt.close(fig)
            print(f"[run_summary] wrote {label}")

    if SHOW_FIGURES:
        plt.show()
    return 0


# =============================================================================
# Section builders
# =============================================================================
class _SectionSkipped(RuntimeError):
    pass


def _section_metadata(cfg, ion, neutral, args) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8.5, 6.0), constrained_layout=True)
    ax.axis("off")
    ax.set_title(f"Run summary -- {args.run_dir.name}")

    lines: list[str] = [f"run_dir = {args.run_dir}"]
    if cfg is not None:
        cfg_dict = asdict(cfg)
        for key in (
            "num_molecules", "T_particles_K", "T_source_K", "dt_ps",
            "potential_steepness", "potential_steepness_molecule",
            "binding_energy_I_atom_K", "binding_energy_molecule_K",
            "binding_energy_I_ion_eV",
        ):
            if key in cfg_dict:
                lines.append(f"  {key} = {cfg_dict[key]}")
    else:
        lines.append("  no cfg.json found")
    if neutral is not None:
        lines.append(
            f"neutral: N={neutral.num_molecules}, "
            f"t=[{neutral.time_ps[0]:.3f},{neutral.time_ps[-1]:.3f}] ps "
            f"({neutral.time_ps.size} steps)"
        )
    if ion is not None:
        lines.append(
            f"ion:     N={ion.num_molecules}, "
            f"t=[{ion.time_ps[0]:.3f},{ion.time_ps[-1]:.3f}] ps "
            f"({ion.time_ps.size} steps)"
        )
    refs = []
    if args.hedft_ref:
        refs.append(f"HeDFT: {args.hedft_ref}")
    if args.ihe_ked_dir:
        refs.append(f"IHe KED (mean-KE + curves): {args.ihe_ked_dir}")
    if args.abundance_ref_path:
        refs.append(f"I+He_n abundance: {args.abundance_ref_path}")
    if args.paper_v2_ref_dir:
        refs.append(f"paper-v2 ref dir: {args.paper_v2_ref_dir}")
    if args.paper_cov_ref_dir:
        refs.append(f"paper-cov ref dir: {args.paper_cov_ref_dir}")
    if refs:
        lines.append("references:")
        lines.extend(f"  {r}" for r in refs)

    ax.text(
        0.0, 1.0, "\n".join(lines),
        family="monospace", fontsize=9,
        verticalalignment="top",
    )
    return fig


def _section_neutral_energy(neutral) -> plt.Figure:
    totals = neutral_energy_totals(neutral)
    fig, ax = plt.subplots(figsize=(8.0, 4.5), constrained_layout=True)
    ax.plot(totals.time_ps, totals.E_kin_eV, label=r"$E_{kin}$")
    ax.plot(totals.time_ps, totals.E_pot_eV, label=r"$E_{pot}$")
    ax.plot(totals.time_ps, totals.E_dissip_eV, label=r"$E_{dissip}$")
    ax.plot(totals.time_ps, totals.E_system_eV, "k", label=r"$E_{system}$")
    ax.set(title="Neutral energy balance", xlabel="t / ps", ylabel="E / eV")
    ax.legend(frameon=False)
    return fig


def _section_ion_energy(ion) -> plt.Figure:
    totals = ion_energy_totals(ion)
    fig, ax = plt.subplots(figsize=(8.0, 4.5), constrained_layout=True)
    ax.plot(totals.time_ps, totals.E_kin_eV, label=r"$E_{kin}$")
    ax.plot(totals.time_ps, totals.E_pot_eV, label=r"$E_{pot}$")
    ax.plot(totals.time_ps, totals.E_dissip_eV, label=r"$E_{dissip}$")
    if totals.E_mass_transfer_eV is not None:
        ax.plot(totals.time_ps, totals.E_mass_transfer_eV,
                label=r"$E_{mass\,transfer}$")
    if totals.E_int_eV is not None:
        ax.plot(totals.time_ps, totals.E_int_eV, label=r"$E_{int}$")
    ax.plot(totals.time_ps, totals.E_system_eV, "k", label=r"$E_{system}$")
    ax.set(title="Ion energy balance (per molecule)",
           xlabel="t / ps", ylabel="E / eV")
    ax.legend(frameon=False, loc="best")
    return fig


def _section_temperature(ion) -> plt.Figure:
    td = ion.temperature_diagnostic
    valid = np.isfinite(td[:, 0])
    if valid.sum() == 0:
        raise _SectionSkipped("no collision rows in temperature_diagnostic")
    t = ion.time_ps[valid]
    td = td[valid]
    fig, ax_left = plt.subplots(figsize=(8.5, 4.5), constrained_layout=True)
    ax_right = ax_left.twinx()
    ax_left.plot(t, td[:, 0], label=r"actual $\langle T'/T\rangle$")
    ax_left.plot(t, td[:, 1], "--",
                 label=r"$\langle T'/T\rangle$ from mass ratio")
    ax_right.plot(t, td[:, 2] * 180.0 / np.pi, color="tab:red",
                  label=r"$\langle\theta\rangle$")
    ax_left.set(title="Ion temperature diagnostic",
                xlabel="t / ps", ylabel=r"$\langle T'/T\rangle$")
    ax_right.set_ylabel(r"$\theta_{lab}$ / deg")
    h1, l1 = ax_left.get_legend_handles_labels()
    h2, l2 = ax_right.get_legend_handles_labels()
    ax_left.legend(h1 + h2, l1 + l2, loc="center left",
                   bbox_to_anchor=(0.0, 0.5), frameon=False)
    return fig


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


def _section_paper_v2_vmi(ion, reference_dir, noise_floor) -> plt.Figure:
    if reference_dir is None:
        raise _SectionSkipped("PAPER_V2_REFERENCE_DIR is None")
    image_ref = load_optional_image(reference_dir, log_prefix="[run_summary]")
    velocity_map = paper_v2_velocity_map(ion, mass_amu=PAPER_V2_MASS_AMU)
    return build_vmi_figure(
        image_ref=image_ref,
        velocity_map=velocity_map,
        experimental_noise_floor=noise_floor,
    )


def _section_paper_cov_radial_distribution(
    ion, reference_dir, paper_v2_reference_dir,
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
            ion, mass_amu=PAPER_V2_MASS_AMU,
        )
        sim_radial = radial_pair_speed_covariance(
            ion, mass_amu=PAPER_V2_MASS_AMU,
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


def _section_paper_cov_phi_distribution(ion, reference_dir) -> plt.Figure:
    if reference_dir is None:
        raise _SectionSkipped("PAPER_COV_REFERENCE_DIR is None")
    phi_ref = load_optional_paper_cov_phi_reference(
        reference_dir, log_prefix="[run_summary]",
    )
    try:
        return build_paper_cov_phi_distribution_figure(
            phi_ref=phi_ref,
            ion=ion,
            mass_amu=PAPER_V2_MASS_AMU,
        )
    except ValueError as exc:
        raise _SectionSkipped(str(exc))


def _section_paper_cov_angular(ion, reference_dir) -> plt.Figure:
    cov_ref = _required_paper_cov_reference(reference_dir)
    try:
        sim_angular = paper_v4_angular_pair_covariance(
            ion, mass_amu=PAPER_V2_MASS_AMU,
        )
    except ValueError as exc:
        raise _SectionSkipped(str(exc))
    return build_paper_cov_angular_cov_figure(
        cov_ref=cov_ref, sim_angular=sim_angular,
    )


def _section_paper_cov_radial(ion, reference_dir) -> plt.Figure:
    cov_ref = _required_paper_cov_reference(reference_dir)
    try:
        sim_radial = radial_pair_speed_covariance(
            ion, mass_amu=PAPER_V2_MASS_AMU,
        )
    except ValueError as exc:
        raise _SectionSkipped(str(exc))
    return build_paper_cov_radial_cov_figure(
        cov_ref=cov_ref, sim_radial=sim_radial,
    )


def _section_paper_cov_traces(ion, reference_dir) -> plt.Figure:
    cov_ref = _required_paper_cov_reference(reference_dir)
    try:
        sim_angular = paper_v4_angular_pair_covariance(
            ion, mass_amu=PAPER_V2_MASS_AMU,
        )
        sim_radial = radial_pair_speed_covariance(
            ion, mass_amu=PAPER_V2_MASS_AMU,
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


def _section_paper_v2_polar(ion, reference_dir, noise_floor) -> plt.Figure:
    if reference_dir is None:
        raise _SectionSkipped("PAPER_V2_REFERENCE_DIR is None")
    polar_ref = load_optional_polar_image(
        reference_dir, log_prefix="[run_summary]"
    )
    if polar_ref is None:
        raise _SectionSkipped("no polar VMI image reference found")
    polar_hist = polar_histogram_matched_to_reference(
        ion, polar_ref, mass_amu=PAPER_V2_MASS_AMU,
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


def _section_radial_evolution(ion) -> plt.Figure:
    ev = radial_distribution_evolution(
        ion,
        n_time_slices=TIME_HEATMAP_N_SLICES,
        n_r_bins=TIME_HEATMAP_N_R_BINS,
    )
    t_c = ev.time_centers_ps
    if t_c.size >= 2:
        t_edges = np.concatenate((
            [t_c[0] - 0.5 * (t_c[1] - t_c[0])],
            0.5 * (t_c[:-1] + t_c[1:]),
            [t_c[-1] + 0.5 * (t_c[-1] - t_c[-2])],
        ))
    else:
        t_edges = np.array([t_c[0] - 0.5, t_c[0] + 0.5])
    fig, ax = plt.subplots(figsize=(9.0, 5.0), constrained_layout=True)
    pcm = ax.pcolormesh(
        t_edges, ev.r_edges_A, ev.counts.T,
        shading="flat", cmap="viridis",
    )
    ax.set(title="Radial distribution evolution |r|(t)",
           xlabel="t / ps", ylabel=r"|r| / $\mathrm{\AA}$")
    fig.colorbar(pcm, ax=ax, label="count")
    return fig


def _section_pair_distance(ion) -> plt.Figure:
    h = interparticle_distance_histogram(ion, num_bins=PAIR_DIST_NUM_BINS)
    fig, ax = plt.subplots(figsize=(8.0, 4.0), constrained_layout=True)
    ax.plot(h.bin_centers_A, h.counts)
    ax.set(title="Final per-molecule I-I separation",
           xlabel=r"|r$_a$ - r$_b$| / $\mathrm{\AA}$",
           ylabel="count")
    return fig


def _section_hedft_neutral(neutral, hedft) -> plt.Figure:
    cmp1 = compare_neutral_to_hedft(neutral, hedft, atom="I1")
    cmp2 = compare_neutral_to_hedft(neutral, hedft, atom="I2")
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5),
                             constrained_layout=True)
    ax_v, ax_z = axes
    ax_v.plot(cmp1.t_overlap_ps, cmp1.v_md_z_Aps, label="MD I1")
    ax_v.plot(cmp1.t_overlap_ps, cmp1.v_hedft_z_Aps, "--", label="HeDFT I1")
    ax_v.plot(cmp2.t_overlap_ps, cmp2.v_md_z_Aps, label="MD I2")
    ax_v.plot(cmp2.t_overlap_ps, cmp2.v_hedft_z_Aps, "--", label="HeDFT I2")
    ax_v.set(title=("Neutral v_z(t)  rmse(I1)="
                    f"{cmp1.rmse_velocity_Aps:.3f} A/ps"),
             xlabel="t / ps", ylabel=r"$v_z$ / $\mathrm{\AA}/\mathrm{ps}$")
    ax_v.legend(frameon=False)

    ax_z.plot(cmp1.t_overlap_ps, cmp1.z_md_A, label="MD I1")
    ax_z.plot(cmp1.t_overlap_ps, cmp1.z_hedft_A, "--", label="HeDFT I1")
    ax_z.plot(cmp2.t_overlap_ps, cmp2.z_md_A, label="MD I2")
    ax_z.plot(cmp2.t_overlap_ps, cmp2.z_hedft_A, "--", label="HeDFT I2")
    ax_z.set(title=(f"z(t) from cumtrapz(v_z)  rmse(I1)="
                    f"{cmp1.rmse_position_A:.3f} A"),
             xlabel="t / ps", ylabel=r"z / $\mathrm{\AA}$")
    ax_z.legend(frameon=False)
    return fig


def _section_hedft_ion(ion, hedft) -> plt.Figure:
    dist = compare_distance(ion, hedft)
    v1 = compare_velocity_magnitude(ion, hedft, atom="I1")
    v2 = compare_velocity_magnitude(ion, hedft, atom="I2")
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5),
                             constrained_layout=True)
    ax_d, ax_v = axes
    ax_d.plot(dist.t_overlap_ps, dist.md_on_hedft_grid, label="MD")
    ax_d.plot(dist.t_overlap_ps, dist.hedft_on_overlap, "--", label="HeDFT")
    ax_d.set(title=f"I-I distance, RMSE={dist.rmse:.3f} A",
             xlabel="t / ps", ylabel=r"R / $\mathrm{\AA}$")
    ax_d.legend(frameon=False)

    ax_v.plot(v1.t_overlap_ps, v1.md_on_hedft_grid, label="MD |v_I1|")
    ax_v.plot(v1.t_overlap_ps, v1.hedft_on_overlap, "--", label="HeDFT |v_I1|")
    ax_v.plot(v2.t_overlap_ps, v2.md_on_hedft_grid, label="MD |v_I2|")
    ax_v.plot(v2.t_overlap_ps, v2.hedft_on_overlap, "--", label="HeDFT |v_I2|")
    ax_v.set(title=(f"|v|(t), RMSE I1={v1.rmse:.3f}, "
                    f"I2={v2.rmse:.3f} A/ps"),
             xlabel="t / ps", ylabel=r"|v| / $\mathrm{\AA}/\mathrm{ps}$")
    ax_v.legend(frameon=False)
    return fig


def _section_boltzmann(cfg, ion) -> plt.Figure:
    if not hasattr(ion, "droplet_radii_angstrom"):
        raise _SectionSkipped("ion checkpoint has no droplet_radii_angstrom")
    radii = np.asarray(ion.droplet_radii_angstrom)
    if radii.size == 0:
        raise _SectionSkipped("droplet_radii_angstrom is empty")
    droplet_radius = float(np.median(radii))

    from i2_helium_md.physics.constants import EV, K_B
    binding_energy_eV = cfg.binding_energy_molecule_K * K_B / EV
    curve = boltzmann_population(
        droplet_radius_A=droplet_radius,
        temperature_K=cfg.T_particles_K,
        steepness_A=cfg.potential_steepness_molecule,
        binding_energy_eV=binding_energy_eV,
    )

    n = ion.num_molecules
    r_init = np.sqrt(
        ion.positions_x[:n, 0] ** 2
        + ion.positions_y[:n, 0] ** 2
        + ion.positions_z[:n, 0] ** 2
    )
    fig, ax = plt.subplots(figsize=(8.0, 4.5), constrained_layout=True)
    if r_init.size > 0:
        ax.hist(r_init, bins=40, density=True, alpha=0.4,
                label="initial |r| histogram")
    ax.plot(curve.r_grid_A, curve.density,
            label=f"Boltzmann T={cfg.T_particles_K:.2f} K, "
                  f"R={droplet_radius:.1f} A")
    ax.set(title="Initial population vs Boltzmann reference",
           xlabel=r"|r| / $\mathrm{\AA}$", ylabel=r"density / $1/\mathrm{\AA}$")
    ax.legend(frameon=False)
    return fig


if __name__ == "__main__":
    raise SystemExit(main())
