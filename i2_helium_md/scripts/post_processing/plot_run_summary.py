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
    interparticle_distance_histogram,
    ion_energy_totals,
    load_he_abundance_reference,
    load_hedft_trajectory,
    load_ihe_ked_reference,
    neutral_energy_totals,
    radial_distribution_evolution,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

# Shared section builders (rule 1: one implementation per recipe, used by
# this script and plot_detection_summary.py alike). The module lives next
# to this script, which is not automatically importable when this file is
# loaded via importlib (the test idiom) -- hence the explicit path insert.
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from summary_sections import (  # noqa: E402, F401  -- constants and
    # _nan_aware_moving_mean are re-exported: the thin interactive
    # comparison scripts (plot_speed_distribution_comparison.py,
    # plot_mean_kinetic_comparison.py, plot_histogram_comparison.py) and
    # the pin tests consume them through this module's namespace.
    DETECTED_KE_INCLUDE_SIM_SE,
    DETECTED_KE_MIN_BIN_COUNT,
    DETECTED_KE_N1_ANCHOR,
    DETECTED_KE_N_MIN,
    HIST_BIN_WIDTH_APS,
    HIST_EDGE_MAX_APS,
    HIST_NUM_BINS,
    HIST_SMOOTHING_WINDOW,
    IHE_KED_EXP_SMOOTHING_WINDOW,
    IHE_KED_HIST_NUM_BINS,
    IHE_KED_HIST_V_MAX_APS,
    IHE_KED_MEAN_ENERGY_SPLIT_N,
    IHE_KED_PLOT_V_MAX_MPS,
    MASS_I,
    MASS_I_HE2_AMU,
    MASS_I_HE_AMU,
    MASS_SPECTRUM_MAX_AMU,
    VELOCITY_PLOT_V_MAX_APS,
    _SectionSkipped,
    _load_detection_read,
    _nan_aware_moving_mean,
    _section_detected_ke_anatomy,
    _section_detected_ked_mean_energy,
    _section_detected_size_distribution,
    _section_ihe_ked_curves,
    _section_ihe_ked_mean_energy,
    _section_mass_resolved,
    _section_mass_spectrum,
    _section_paper_cov_angular,
    _section_paper_cov_phi_distribution,
    _section_paper_cov_radial,
    _section_paper_cov_radial_distribution,
    _section_paper_cov_traces,
    _section_paper_v2_polar,
    _section_paper_v2_vmi,
)


# =============================================================================
# Plot-tuning constants for the window-only sections (the detector-facing
# recipe constants moved to summary_sections.py with their builders)
# =============================================================================
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
    # Tier-2 confirmation runs only: None on a legacy run dir (no
    # detection.npz), which keeps the section list — and every legacy
    # figure — exactly as before.
    detection_read = _load_detection_read(run.root)

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
                     ion, paper_v2_ref_dir, EXPERIMENTAL_NOISE_FLOOR,
                     mass_amu=PAPER_V2_MASS_AMU)),
                ("paper_cov_radial_distribution",
                 lambda: _section_paper_cov_radial_distribution(
                     ion, paper_cov_ref_dir, paper_v2_ref_dir,
                     mass_amu=PAPER_V2_MASS_AMU)),
                ("paper_cov_phi_distribution",
                 lambda: _section_paper_cov_phi_distribution(
                     ion, paper_cov_ref_dir, mass_amu=PAPER_V2_MASS_AMU)),
                ("paper_v2_polar_image_comparison",
                 lambda: _section_paper_v2_polar(
                     ion, paper_v2_ref_dir, EXPERIMENTAL_NOISE_FLOOR,
                     mass_amu=PAPER_V2_MASS_AMU)),
                ("mass_resolved_velocities",
                 lambda: _section_mass_resolved(ion)),
                ("radial_evolution_heatmap",
                 lambda: _section_radial_evolution(ion)),
                ("interparticle_distance_histogram",
                 lambda: _section_pair_distance(ion)),
                ("paper_cov_angular_pair_cov",
                 lambda: _section_paper_cov_angular(
                     ion, paper_cov_ref_dir, mass_amu=PAPER_V2_MASS_AMU)),
                ("paper_cov_radial_pair_cov",
                 lambda: _section_paper_cov_radial(
                     ion, paper_cov_ref_dir, mass_amu=PAPER_V2_MASS_AMU)),
                ("paper_cov_pair_cov_traces",
                 lambda: _section_paper_cov_traces(
                     ion, paper_cov_ref_dir, mass_amu=PAPER_V2_MASS_AMU)),
            ]
        if detection_read is not None:
            sections += [
                ("detected_size_distribution",
                 lambda: _section_detected_size_distribution(
                     detection_read, abundance)),
                ("detected_ihe_ked_mean_energy",
                 lambda: _section_detected_ked_mean_energy(
                     detection_read, ked_ref)),
                ("detected_ke_anatomy",
                 lambda: _section_detected_ke_anatomy(
                     detection_read, ked_ref)),
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
# Section builders (window-only; detector-facing builders live in
# summary_sections.py)
# =============================================================================
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
