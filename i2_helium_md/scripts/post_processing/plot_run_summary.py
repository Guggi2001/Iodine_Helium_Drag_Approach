"""Consolidated post-processing summary for one finished run directory.

Produces every numerical diagnostic the legacy MATLAB post-processing
scripts produced *that is in scope per CLAUDE.md*, in a single
multi-page PDF plus per-figure PNGs. Sections that need optional
reference data (HeDFT trajectory, experimental VMI) are gated on the
corresponding ``--*-ref`` arguments so each run type stays separable
(per CLAUDE.md "Keep these workflows separate.").

The following legacy scripts are consolidated here:

* ``vmi_sim_3d_neutral_propa_HeDFT_mimic.m``    -> neutral energy balance
* ``vmi_sim_3d_ion_propa.m``                   -> ion energy balance + temperature
* ``simulation_image_only_trajectories.m``     -> HeDFT R(t), v(t)
* ``post_process_single_pulse_paper_v3.m``     -> 1D and 2D polar VMI panels
* ``post_process_single_pulse_paper_v4.m``     -> angular pair covariance
* ``post_process_single_pulse_paper.m``        -> bimodal Gaussian fit
* ``post_process_single_pulse.m``              -> 2D (vx, vy) histogram
* ``post_process_compare_radial_distributions.m`` -> time-resolved radial,
                                                     interatomic distance,
                                                     Boltzmann reference
* ``compare_neutral_dynamics_to_HeDFT.m``      -> neutral cumtrapz r(t)

Out of scope (deferred per CLAUDE.md): Abel inversion, pump-probe,
effusive / gas-phase comparison, live-debug 3D animations.

Usage::

    python scripts/post_processing/plot_run_summary.py <run_dir> \
        [--hedft-ref PATH] [--vmi-ref-he PATH] [--vmi-ref-gas PATH] \
        [--out-dir PATH] [--no-show]
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path
import sys


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
    angular_pair_covariance,
    anisotropy_fit,
    beta_of_velocity,
    bimodal_gaussian_fit,
    boltzmann_population,
    compare_distance,
    compare_neutral_to_hedft,
    compare_velocity_magnitude,
    compute_final_velocity_histogram,
    interparticle_distance_histogram,
    ion_energy_totals,
    load_hedft_trajectory,
    load_vmi_reference,
    mass_spectrum,
    neutral_energy_totals,
    phi_histogram,
    polar_velocity_histogram,
    radial_distribution_evolution,
    velocity_density_2d,
)
from i2_helium_md.postprocess._smoothing import (  # noqa: E402
    moving_mean,
    normalise_trace,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


# =============================================================================
# Plot-tuning constants (kept aligned with the existing focused scripts)
# =============================================================================
MASS_I_HE_AMU = 131.0
MASS_I_HE2_AMU = 135.0

HIST_BIN_WIDTH_APS = 0.04
HIST_EDGE_MAX_APS = 26.0
HIST_NUM_BINS = int(round(HIST_EDGE_MAX_APS / HIST_BIN_WIDTH_APS))
HIST_SMOOTHING_WINDOW = 15
VELOCITY_PLOT_V_MAX_APS = 28.0

PHI_BIN_WIDTH_RAD = 0.05
PHI_SMOOTHING_WINDOW = 15

POLAR_N_V_BINS = 80
POLAR_N_PHI_BINS = 72

VEL2D_N_BINS = 200
VEL2D_V_MAX_APS = 22.0

PAIR_DIST_NUM_BINS = 100
PAIR_COV_N_THETA_BINS = 50

TIME_HEATMAP_N_SLICES = 60
TIME_HEATMAP_N_R_BINS = 100


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build a consolidated post-processing PDF for one run."
    )
    p.add_argument("run_dir", type=Path, help="Path to a run directory.")
    p.add_argument("--hedft-ref", type=Path, default=None,
                   help="Optional HeDFT reference CSV (e.g. 9A_All_Data.csv).")
    p.add_argument("--vmi-ref-he", type=Path, default=None,
                   help="Optional experimental I+He VMI reference CSV.")
    p.add_argument("--vmi-ref-gas", type=Path, default=None,
                   help="Optional experimental I+ gas-phase VMI reference CSV.")
    p.add_argument("--out-dir", type=Path, default=None,
                   help="Output directory; default is <run_dir>/figures.")
    p.add_argument("--no-show", action="store_true",
                   help="Skip plt.show() (useful in headless smoke tests).")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    run = RunDirectory(args.run_dir)

    out_dir = args.out_dir or (run.root / "figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "run_summary.pdf"

    cfg = run.load_cfg() if run.has_cfg() else None
    neutral = run.load_neutral() if run.has_neutral() else None
    ion = run.load_ion() if run.has_ion() else None

    hedft = load_hedft_trajectory(args.hedft_ref) if args.hedft_ref else None
    vmi_he = load_vmi_reference(args.vmi_ref_he) if args.vmi_ref_he else None
    vmi_gas = load_vmi_reference(args.vmi_ref_gas) if args.vmi_ref_gas else None

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
                 lambda: _section_mass_spectrum(ion)),
                ("radial_velocity_with_vmi",
                 lambda: _section_radial_velocity(ion, vmi_he, vmi_gas)),
                ("phi_histogram",
                 lambda: _section_phi(ion)),
                ("polar_v_phi_histogram",
                 lambda: _section_polar(ion)),
                ("anisotropy_fit",
                 lambda: _section_anisotropy(ion)),
                ("velocity_density_2d",
                 lambda: _section_velocity_2d(ion)),
                ("mass_resolved_velocities",
                 lambda: _section_mass_resolved(ion)),
                ("radial_evolution_heatmap",
                 lambda: _section_radial_evolution(ion)),
                ("interparticle_distance_histogram",
                 lambda: _section_pair_distance(ion)),
                ("angular_pair_covariance",
                 lambda: _section_pair_covariance(ion)),
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

    if not args.no_show:
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
    if args.vmi_ref_he:
        refs.append(f"VMI(I+He): {args.vmi_ref_he}")
    if args.vmi_ref_gas:
        refs.append(f"VMI(gas): {args.vmi_ref_gas}")
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
    if totals.E_mass_attach_defect_eV is not None:
        ax.plot(totals.time_ps, totals.E_mass_attach_defect_eV,
                label=r"$E_{mass\,defect}$")
    ax.plot(totals.time_ps, totals.E_system_eV, "k", label=r"$E_{system}$")
    ax.set(title="Ion energy balance (per molecule)",
           xlabel="t / ps", ylabel="E / eV")
    ax.legend(frameon=False)
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
    ax_left.legend(h1 + h2, l1 + l2, frameon=False)
    return fig


def _section_mass_spectrum(ion) -> plt.Figure:
    spec = mass_spectrum(ion, bin_width_amu=1.0)
    fig, ax = plt.subplots(figsize=(8.0, 4.0), constrained_layout=True)
    ax.bar(spec.bin_centers_amu, spec.counts, width=0.9,
           edgecolor="black", linewidth=0.5)
    ax.set(title="Final ion mass spectrum",
           xlabel="m / u", ylabel="count")
    return fig


def _section_radial_velocity(ion, vmi_he, vmi_gas) -> plt.Figure:
    try:
        sim_he = compute_final_velocity_histogram(
            ion, mass_amu=MASS_I_HE_AMU,
            num_bins=HIST_NUM_BINS, v_max_Aps=HIST_EDGE_MAX_APS,
        )
    except ValueError as exc:
        raise _SectionSkipped(str(exc))
    try:
        sim_he2 = compute_final_velocity_histogram(
            ion, mass_amu=MASS_I_HE2_AMU,
            num_bins=HIST_NUM_BINS, v_max_Aps=HIST_EDGE_MAX_APS,
        )
    except ValueError:
        sim_he2 = None

    fig, ax = plt.subplots(figsize=(9.0, 4.5), constrained_layout=True)
    if vmi_gas is not None:
        mask = vmi_gas.velocity_Aps > 4.0
        max_gas = float(vmi_gas.signal_arb[mask].max())
        ax.plot(vmi_gas.velocity_Aps, vmi_gas.signal_arb / max_gas,
                label=r"exp. $I_2$:$I^+$")
    if vmi_he is not None:
        max_he = float(vmi_he.signal_arb.max())
        ax.plot(vmi_he.velocity_Aps, vmi_he.signal_arb / max_he,
                ":", label=r"exp. $I_2 He_N$:$I^+ He$")
    sim_he_smooth = normalise_trace(
        moving_mean(sim_he.density, HIST_SMOOTHING_WINDOW)
    )
    ax.plot(sim_he.bin_centers_Aps, sim_he_smooth, "--",
            label=r"sim. $I^+ He$")
    if sim_he2 is not None:
        sim_he2_smooth = normalise_trace(
            moving_mean(sim_he2.density, HIST_SMOOTHING_WINDOW)
        )
        ax.plot(sim_he2.bin_centers_Aps, sim_he2_smooth, "-.",
                label=r"sim. $I^+ He_2$")

    bim = bimodal_gaussian_fit(sim_he)
    if bim.success:
        norm = sim_he_smooth.max() / max(sim_he.density.max(), 1e-12)
        ax.plot(sim_he.bin_centers_Aps, bim.fitted_curve * norm,
                color="tab:gray", linewidth=1.0,
                label=(f"bimodal fit "
                       f"(mu1={bim.mean_1_Aps:.1f}, mu2={bim.mean_2_Aps:.1f})"))

    ax.set(title="(a) radial velocity with experimental VMI overlay",
           xlim=(0.0, VELOCITY_PLOT_V_MAX_APS), ylim=(0.0, 1.1),
           xlabel=r"v / $\mathrm{\AA}/\mathrm{ps}$",
           ylabel="signal / arb. units")
    ax.legend(frameon=False, fontsize=9)
    return fig


def _section_phi(ion) -> plt.Figure:
    try:
        ph = phi_histogram(
            ion, bin_width_rad=PHI_BIN_WIDTH_RAD, mass_amu=MASS_I_HE_AMU,
        )
    except ValueError as exc:
        raise _SectionSkipped(str(exc))
    fig, ax = plt.subplots(figsize=(8.0, 4.0), constrained_layout=True)
    if ph.density.sum() > 0:
        smoothed = normalise_trace(
            moving_mean(ph.density, PHI_SMOOTHING_WINDOW)
        )
    else:
        smoothed = np.zeros_like(ph.density, dtype=float)
    ax.plot(ph.bin_centers_rad, smoothed)
    ax.set(title=f"Azimuthal phi distribution, m={MASS_I_HE_AMU:.0f} u "
                 f"(n={ph.num_atoms_used})",
           xlabel=r"$\varphi$ / rad", ylabel="signal / arb. units",
           xlim=(0.0, 2.0 * np.pi), ylim=(0.0, 1.1))
    return fig


def _section_polar(ion) -> plt.Figure:
    pol = polar_velocity_histogram(
        ion, n_v_bins=POLAR_N_V_BINS, n_phi_bins=POLAR_N_PHI_BINS,
        v_max_Aps=VELOCITY_PLOT_V_MAX_APS, mass_amu=MASS_I_HE_AMU,
    )
    if pol.num_atoms_used == 0:
        raise _SectionSkipped("no atoms passed mass+outside filter")
    fig, ax = plt.subplots(figsize=(8.5, 4.5), constrained_layout=True)
    pcm = ax.pcolormesh(
        pol.phi_edges_rad, pol.v_edges_Aps, pol.counts,
        shading="auto", cmap="magma",
    )
    ax.set(title=f"Polar velocity histogram, m={MASS_I_HE_AMU:.0f} u "
                 f"(n={pol.num_atoms_used})",
           xlabel=r"$\varphi$ / rad",
           ylabel=r"|v| / $\mathrm{\AA}/\mathrm{ps}$")
    fig.colorbar(pcm, ax=ax, label="count")
    return fig


def _section_anisotropy(ion) -> plt.Figure:
    pol = polar_velocity_histogram(
        ion, n_v_bins=POLAR_N_V_BINS, n_phi_bins=POLAR_N_PHI_BINS,
        v_max_Aps=VELOCITY_PLOT_V_MAX_APS, mass_amu=MASS_I_HE_AMU,
    )
    if pol.num_atoms_used == 0:
        raise _SectionSkipped("no atoms passed mass+outside filter")
    fit = anisotropy_fit(pol)
    curve = beta_of_velocity(pol, min_counts_per_v_bin=50)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5),
                             constrained_layout=True)
    ax_phi, ax_beta = axes
    phi_marg = pol.counts.sum(axis=0)
    ax_phi.plot(pol.phi_centers_rad, phi_marg, label="counts")
    if fit.success:
        from i2_helium_md.postprocess.polar_velocity import _cos2_model
        ax_phi.plot(
            pol.phi_centers_rad,
            _cos2_model(pol.phi_centers_rad, fit.a, fit.b, fit.phi0_rad),
            "--", label=fr"cos$^2$ fit, $\beta=${fit.beta:.2f}",
        )
    ax_phi.set(
        title="Marginal phi + cos^2 fit",
        xlabel=r"$\varphi$ / rad", ylabel="count",
        xlim=(0.0, 2.0 * np.pi),
    )
    ax_phi.legend(frameon=False)

    ax_beta.plot(curve.v_centers_Aps[curve.valid],
                 curve.beta[curve.valid], "o-")
    ax_beta.axhline(0.0, color="grey", linewidth=0.5)
    ax_beta.set(
        title=r"$\beta(v)$",
        xlabel=r"v / $\mathrm{\AA}/\mathrm{ps}$",
        ylabel=r"$\beta$",
    )
    return fig


def _section_velocity_2d(ion) -> plt.Figure:
    h = velocity_density_2d(
        ion, axes=("x", "y"), n_bins=VEL2D_N_BINS, v_max_Aps=VEL2D_V_MAX_APS,
        mass_amu=MASS_I_HE_AMU,
    )
    fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
    pcm = ax.pcolormesh(
        h.bin_edges_a_Aps, h.bin_edges_b_Aps, h.counts.T,
        shading="auto", cmap="magma",
    )
    ax.set_aspect("equal")
    ax.set(title=f"2D velocity density, axes=({h.axis_a},{h.axis_b}), "
                 f"m={MASS_I_HE_AMU:.0f} u (n={h.num_atoms_used})",
           xlabel=fr"$v_{h.axis_a}$ / $\mathrm{{\AA}}/\mathrm{{ps}}$",
           ylabel=fr"$v_{h.axis_b}$ / $\mathrm{{\AA}}/\mathrm{{ps}}$")
    fig.colorbar(pcm, ax=ax, label="count")
    return fig


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


def _section_pair_covariance(ion) -> plt.Figure:
    cov = angular_pair_covariance(
        ion, n_theta_bins=PAIR_COV_N_THETA_BINS, mass_amu=MASS_I_HE_AMU,
    )
    if cov.num_pairs_used == 0:
        raise _SectionSkipped("no pairs passed mass+outside filter")
    fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
    pcm = ax.pcolormesh(
        cov.theta_edges_rad, cov.theta_edges_rad, cov.counts.T,
        shading="auto", cmap="magma",
    )
    ax.set_aspect("equal")
    ax.set(
        title=f"Angular pair covariance (n_pairs={cov.num_pairs_used})",
        xlabel=r"$\theta_a$ / rad",
        ylabel=r"$\theta_b$ / rad",
    )
    fig.colorbar(pcm, ax=ax, label="count")
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
