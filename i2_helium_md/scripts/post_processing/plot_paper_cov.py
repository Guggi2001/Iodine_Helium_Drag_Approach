"""Paper-cov port of post_process_single_pulse_paper_IplusHe_comparison_cov.m.

The legacy MATLAB script extends the I+He comparison figure with two
pair-covariance diagnostics computed from the experimental I+He droplet
VMI measurements (IDs 45668, 45662, 45667). This Python port splits the
output into six standalone figures:

- ``paper_cov_vmi_comparison.png``            experimental + simulated 2-D VMI
- ``paper_cov_radial_distribution.png``       1-D velocity distribution overlay
                                              (I+ gas, I+He, exp v-cov trace,
                                              sim v-cov trace)
- ``paper_cov_angular_pair_cov.png``          experimental + simulated angular
                                              pair covariance
- ``paper_cov_radial_pair_cov.png``           experimental + simulated radial
                                              pair-speed covariance
- ``paper_cov_phi_distribution.png``          experimental + simulated 1-D
                                              phi(angle) distribution overlay
- ``paper_cov_pair_cov_traces.png``           1-D axis-sum traces of the
                                              angular and radial pair
                                              covariance (sim vs exp overlay)

The experimental pair-covariance matrices are loaded from the frozen
MATLAB reference under ``data/reference/paper_cov/``; the simulated
counterparts are computed from the Python ion checkpoint. If the
covariance reference is not present the two heatmap covariance figures
and the 1-D trace figure are skipped (and the 1-D radial figure drops
the experimental v-cov trace). The phi-distribution figure always
renders; the experimental overlay curve is loaded from the precomputed
``data/reference/paper_cov/iplus_he_phi.csv`` (literal MATLAB output of
``mean(res_Iplus_He.image_polar(:, b_r), 2) / max(...)`` for the same
three measurement IDs as the covariance reference). If that CSV is
absent, only the simulated curve is drawn.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# USER SETTINGS / DEFAULTS
# =============================================================================
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet_long"
REFERENCE_DIR = PROJECT_ROOT / "data" / "reference" / "paper_cov"
PAPER_V2_REFERENCE_DIR = PROJECT_ROOT / "data" / "reference" / "paper_v2"
MASS_SELECTION_AMU = 131.0
EXPERIMENTAL_NOISE_FLOOR = 0.20

# Single experimental reference overlaid on panel (c): the high-SNR I+He
# radial. The MATLAB _cov.m script overlays (gas 43632 + averaged 300 mW
# I+He) instead, but the Python port substitutes the higher-SNR
# `res_sum` curve so the experimental shape is cleaner and drops the gas
# trace as uninformative for the I+He covariance discussion.
HIGH_SNR_RADIAL_FILENAME = "iplus_he_high_snr_radial.csv"

# 1-D phi(angle) experimental reference for panel (f). Produced by the
# extended `export_paper_cov_reference_data.m` from the same averaged
# res_Iplus_He as the covariance reference; literal MATLAB recipe
# (mean over polar-image radii, max-normalised).
EXPERIMENTAL_PHI_FILENAME = "iplus_he_phi.csv"


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    covariance_axis_sum_normalised,
    load_paper_cov_experimental_reference,
    load_paper_v2_phi_reference,
    load_paper_v2_radial_reference,
    load_paper_v2_vmi_image_reference,
    paper_v2_velocity_curve,
    paper_v2_velocity_map,
    paper_v4_angular_pair_covariance,
    radial_covariance_trace,
    radial_pair_speed_covariance,
    simulated_phi_distribution,
)
from i2_helium_md.postprocess.paper_cov import (  # noqa: E402
    PAPER_COV_COLOR_CLIP_FRACTION,
)
from i2_helium_md.postprocess.paper_v2 import max_normalise  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Reproduce the active droplet branch of "
        "post_process_single_pulse_paper_IplusHe_comparison_cov.m."
    )
    p.add_argument("--run-dir", type=Path, default=RUN_DIR,
                   help="RunDirectory containing cfg.json, neutral.npz, and ion.npz.")
    p.add_argument("--reference-dir", type=Path, default=REFERENCE_DIR,
                   help="Directory containing iplus_he_covariance.mat (paper_cov references).")
    p.add_argument("--paper-v2-reference-dir", type=Path, default=PAPER_V2_REFERENCE_DIR,
                   help="Directory containing the paper-v2 radial CSV references.")
    p.add_argument("--no-show", action="store_true",
                   help="Skip plt.show(), useful for tests and headless runs.")
    p.add_argument("--noise-floor", type=float, default=EXPERIMENTAL_NOISE_FLOOR,
                   help="Fraction of max intensity below which experimental "
                        "VMI pixels clip to background (default: %(default)s; "
                        "set 0 to disable). Simulated panels are unaffected.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    run = RunDirectory(args.run_dir)
    ion = run.load_ion()

    # ----- references -----
    high_snr_ref = _load_optional_high_snr_radial(args.paper_v2_reference_dir)
    image_ref = _load_optional_vmi_image(args.paper_v2_reference_dir)
    cov_ref = _load_optional_cov_reference(args.reference_dir)
    phi_ref = _load_optional_paper_cov_phi_reference(args.reference_dir)

    # ----- simulated diagnostics -----
    velocity_map = paper_v2_velocity_map(ion, mass_amu=MASS_SELECTION_AMU)
    sim_radial_curve = paper_v2_velocity_curve(ion, mass_amu=MASS_SELECTION_AMU)
    sim_angular = paper_v4_angular_pair_covariance(ion, mass_amu=MASS_SELECTION_AMU)
    sim_radial = radial_pair_speed_covariance(ion, mass_amu=MASS_SELECTION_AMU)

    # ----- figures -----
    figures: dict[str, plt.Figure] = {}

    figures["paper_cov_vmi_comparison.png"] = _build_vmi_figure(
        image_ref=image_ref,
        velocity_map=velocity_map,
        experimental_noise_floor=args.noise_floor,
    )
    figures["paper_cov_radial_distribution.png"] = _build_radial_distribution_figure(
        high_snr_ref=high_snr_ref,
        sim_radial_curve=sim_radial_curve,
        sim_radial=sim_radial,
        cov_ref=cov_ref,
    )
    figures["paper_cov_phi_distribution.png"] = _build_phi_distribution_figure(
        phi_ref=phi_ref, ion=ion
    )
    if cov_ref is not None:
        figures["paper_cov_angular_pair_cov.png"] = _build_angular_cov_figure(
            cov_ref=cov_ref, sim_angular=sim_angular
        )
        figures["paper_cov_radial_pair_cov.png"] = _build_radial_cov_figure(
            cov_ref=cov_ref, sim_radial=sim_radial
        )
        figures["paper_cov_pair_cov_traces.png"] = _build_pair_cov_traces_figure(
            cov_ref=cov_ref,
            sim_angular=sim_angular,
            sim_radial=sim_radial,
        )

    out_dir = run.root / "figures"
    out_dir.mkdir(exist_ok=True)
    for filename, fig in figures.items():
        fig.savefig(out_dir / filename, dpi=150)
        print(f"Saved {filename} to {out_dir / filename}")

    if not args.no_show:
        plt.show()
    return 0


# -----------------------------------------------------------------------------
# Optional reference loaders
# -----------------------------------------------------------------------------
def _load_optional_high_snr_radial(reference_dir: Path):
    path = reference_dir / HIGH_SNR_RADIAL_FILENAME
    if not path.exists():
        print(
            f"[paper_cov] optional high-SNR I+He radial reference not found: {path}"
        )
        return None
    return load_paper_v2_radial_reference(path)


def _load_optional_vmi_image(reference_dir: Path):
    image_dir = reference_dir / "images"
    for ext in (".npz", ".mat"):
        path = image_dir / f"iplus_he_high_snr_vmi_image{ext}"
        if path.exists():
            return load_paper_v2_vmi_image_reference(path)
    print(f"[paper_cov] optional 2-D VMI image reference not found in {image_dir}")
    return None


def _load_optional_cov_reference(reference_dir: Path):
    for filename in ("iplus_he_covariance.mat", "iplus_he_covariance.npz"):
        path = reference_dir / filename
        if path.exists():
            return load_paper_cov_experimental_reference(path)
    print(
        f"[paper_cov] experimental covariance reference not found in {reference_dir}; "
        "run data/reference/scripts/export_paper_cov_reference_data.m to generate it. "
        "Skipping the two covariance figures."
    )
    return None


def _load_optional_paper_cov_phi_reference(reference_dir: Path):
    path = reference_dir / EXPERIMENTAL_PHI_FILENAME
    if not path.exists():
        print(
            f"[paper_cov] experimental phi reference not found: {path}. "
            "Re-run data/reference/scripts/export_paper_cov_reference_data.m "
            "to generate it. The phi figure will draw the simulated curve only."
        )
        return None
    return load_paper_v2_phi_reference(path)


# -----------------------------------------------------------------------------
# Figure builders
# -----------------------------------------------------------------------------
def _build_vmi_figure(
    *, image_ref, velocity_map, experimental_noise_floor: float
) -> plt.Figure:
    fig, (ax_exp, ax_sim) = plt.subplots(
        1, 2, figsize=(10.0, 4.8), constrained_layout=True
    )
    _draw_experimental_vmi(
        ax_exp, image_ref, noise_floor_fraction=experimental_noise_floor
    )
    _draw_simulated_vmi(ax_sim, velocity_map)
    return fig


def _build_radial_distribution_figure(
    *, high_snr_ref, sim_radial_curve, sim_radial, cov_ref
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
    handles, labels = [], []

    if high_snr_ref is not None:
        (line,) = ax.plot(
            high_snr_ref.velocity_mps,
            max_normalise(high_snr_ref.signal_arb),
            lw=1.6,
            color="tab:blue",
        )
        handles.append(line)
        labels.append("I+He high-SNR")

    if cov_ref is not None:
        exp_trace = radial_covariance_trace(
            cov_ref.cov_radial, cov_ref.velocity_centers_Aps
        )
        (line,) = ax.plot(
            cov_ref.velocity_centers_mps,
            max_normalise(exp_trace),
            lw=1.2,
            ls="--",
            color="black",
        )
        handles.append(line)
        labels.append("exp v-cov trace")

    (line,) = ax.plot(
        sim_radial_curve.bin_centers_mps,
        sim_radial_curve.normalised,
        lw=0.9,
        color="tab:orange",
    )
    handles.append(line)
    labels.append("sim radial")

    sim_trace = radial_covariance_trace(
        sim_radial.counts, sim_radial.velocity_centers_Aps
    )
    scatter = ax.scatter(
        sim_radial.velocity_centers_mps,
        max_normalise(sim_trace),
        s=18,
        color="tab:red",
        marker="o",
        zorder=3,
    )
    handles.append(scatter)
    labels.append("sim v-cov trace")

    ax.set_xlabel("v / m/s")
    ax.set_ylabel("signal / arb. units")
    ax.set_xlim(0.0, 2800.0)
    ax.set_ylim(0.0, 1.05)
    ax.set_title("(c) velocity distribution comparison")
    ax.legend(handles, labels, fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)
    return fig


def _build_angular_cov_figure(*, cov_ref, sim_angular) -> plt.Figure:
    fig, (ax_exp, ax_sim) = plt.subplots(
        1, 2, figsize=(10.0, 4.8), constrained_layout=True
    )
    _draw_cov_heatmap(
        ax_exp,
        cov_ref.cov_angular,
        x=cov_ref.theta_centers_rad,
        y=cov_ref.theta_centers_rad,
        xlabel="angle / rad",
        ylabel="angle / rad",
        title="(d) experimental angular pair covariance",
    )
    # Roll the simulated [0, 2pi) matrix into [-pi, pi] for visual parity.
    centers_2pi = sim_angular.theta_centers_rad
    counts = np.asarray(sim_angular.counts, dtype=float)
    counts_centered, theta_centered = _roll_to_centered_theta(counts, centers_2pi)
    _draw_cov_heatmap(
        ax_sim,
        counts_centered,
        x=theta_centered,
        y=theta_centered,
        xlabel="angle / rad",
        ylabel="angle / rad",
        title=f"(d') simulated angular pair covariance, "
              f"{sim_angular.num_pairs_used} pairs",
    )
    return fig


def _build_radial_cov_figure(*, cov_ref, sim_radial) -> plt.Figure:
    fig, (ax_exp, ax_sim) = plt.subplots(
        1, 2, figsize=(10.0, 4.8), constrained_layout=True
    )
    _draw_cov_heatmap(
        ax_exp,
        cov_ref.cov_radial,
        x=cov_ref.velocity_centers_Aps,
        y=cov_ref.velocity_centers_Aps,
        xlabel=r"v / $\AA$/ps",
        ylabel=r"v / $\AA$/ps",
        title="(e) experimental radial pair covariance",
    )
    _draw_cov_heatmap(
        ax_sim,
        sim_radial.counts,
        x=sim_radial.velocity_centers_Aps,
        y=sim_radial.velocity_centers_Aps,
        xlabel=r"v / $\AA$/ps",
        ylabel=r"v / $\AA$/ps",
        title=f"(e') simulated radial pair covariance, "
              f"{sim_radial.num_pairs_used} pairs",
    )
    return fig


def _build_phi_distribution_figure(*, phi_ref, ion) -> plt.Figure:
    """Overlay of experimental + simulated 1-D phi distribution.

    Mirrors the MATLAB overlay at lines 187-205 (experimental) and 240-252
    (simulated) of ``_cov.m``. The experimental curve is loaded from the
    precomputed CSV ``data/reference/paper_cov/iplus_he_phi.csv`` (the
    literal MATLAB output ``mean(res_Iplus_He.image_polar(:, b_r), 2) /
    max(...)`` computed in the exporter from the same three measurement
    IDs as the covariance reference). If that CSV is missing, the figure
    renders the simulated curve alone.
    """

    fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
    handles, labels = [], []

    if phi_ref is not None:
        (line,) = ax.plot(
            phi_ref.phi_rad,
            max_normalise(phi_ref.signal_arb),
            lw=1.6,
            color="tab:blue",
        )
        handles.append(line)
        labels.append("I+He (exp)")

    sim_phi = simulated_phi_distribution(ion, mass_amu=MASS_SELECTION_AMU)
    (line,) = ax.plot(
        sim_phi.phi_centers_rad,
        sim_phi.signal_normalised,
        lw=1.1,
        color="tab:orange",
    )
    handles.append(line)
    labels.append(f"I+He (sim), {sim_phi.num_samples_used} atoms")

    ax.set_xlabel(r"$\phi$ / radian")
    ax.set_ylabel("signal / arb. units")
    ax.set_xlim(0.0, 2.0 * np.pi)
    ax.set_ylim(0.0, 1.05)
    ax.set_title("(f) phi angular distribution")
    ax.legend(handles, labels, fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)
    return fig


def _build_pair_cov_traces_figure(*, cov_ref, sim_angular, sim_radial) -> plt.Figure:
    """Two side-by-side 1-D pair-covariance axis-sum traces.

    Mirrors lines 486-522 of ``_cov.m`` (nexttiles 2 and 3 of the
    second MATLAB figure): for both the angular and radial pair
    covariance matrices, sum along axis 0, smooth with ``movmean(., 3)``,
    subtract minimum, and normalise by maximum. Sim is solid, exp is
    dashed. The simulated angular matrix is rolled into the
    ``[-pi, pi)`` axis range first so it lines up with the experimental
    ``theta_centers_rad`` axis (the same rolling step already used by
    the angular covariance heatmap figure).
    """

    fig, (ax_ang, ax_rad) = plt.subplots(
        1, 2, figsize=(10.0, 4.5), constrained_layout=True
    )

    sim_ang_counts = np.asarray(sim_angular.counts, dtype=float)
    sim_ang_centered, theta_centered = _roll_to_centered_theta(
        sim_ang_counts, sim_angular.theta_centers_rad
    )
    sim_ang_trace = covariance_axis_sum_normalised(sim_ang_centered, axis=0)
    exp_ang_trace = covariance_axis_sum_normalised(cov_ref.cov_angular, axis=0)
    ax_ang.plot(
        theta_centered, sim_ang_trace,
        lw=1.4, color="tab:orange", label="sim",
    )
    ax_ang.plot(
        cov_ref.theta_centers_rad, exp_ang_trace,
        lw=1.2, ls="--", color="black", label="exp",
    )
    ax_ang.set_xlabel(r"$\theta$ / radian")
    ax_ang.set_ylabel("normalised trace / arb. units")
    ax_ang.set_xlim(-np.pi, np.pi)
    ax_ang.set_ylim(0.0, 1.05)
    ax_ang.set_title("(g) angular pair-cov trace")
    ax_ang.legend(fontsize=9, loc="upper right")
    ax_ang.grid(True, alpha=0.3)

    sim_rad_trace = covariance_axis_sum_normalised(sim_radial.counts, axis=0)
    exp_rad_trace = covariance_axis_sum_normalised(cov_ref.cov_radial, axis=0)
    ax_rad.plot(
        sim_radial.velocity_centers_Aps, sim_rad_trace,
        lw=1.4, color="tab:orange", label="sim",
    )
    ax_rad.plot(
        cov_ref.velocity_centers_Aps, exp_rad_trace,
        lw=1.2, ls="--", color="black", label="exp",
    )
    ax_rad.set_xlabel(r"$v_r$ / $\AA$/ps")
    ax_rad.set_ylabel("normalised trace / arb. units")
    ax_rad.set_ylim(0.0, 1.05)
    ax_rad.set_title("(h) radial pair-cov trace")
    ax_rad.legend(fontsize=9, loc="upper right")
    ax_rad.grid(True, alpha=0.3)

    return fig


# -----------------------------------------------------------------------------
# Drawing helpers
# -----------------------------------------------------------------------------
def _draw_experimental_vmi(ax, image_ref, *, noise_floor_fraction: float) -> None:
    cmap_floor = plt.get_cmap("magma")(0.0)
    if image_ref is None:
        ax.text(
            0.5, 0.5, "experimental VMI image not exported",
            ha="center", va="center",
        )
        ax.set_xlim(-3500.0, 3500.0)
        ax.set_ylim(-3500.0, 3500.0)
    else:
        mesh = ax.pcolormesh(
            image_ref.vx_mps,
            image_ref.vy_mps,
            image_ref.intensity,
            shading="auto",
            cmap="magma",
            norm=_color_norm(
                image_ref.intensity, noise_floor_fraction=noise_floor_fraction
            ),
        )
        plt.colorbar(mesh, ax=ax, label="signal / arb. units")
    ax.set_title("(a) experimental I+He VMI")
    ax.set_xlabel("v_x / m/s")
    ax.set_ylabel("v_y / m/s")
    ax.set_aspect("equal")
    ax.set_xlim(-3500.0, 3500.0)
    ax.set_ylim(-3500.0, 3500.0)
    ax.set_facecolor(cmap_floor)


def _draw_simulated_vmi(ax, velocity_map) -> None:
    bins_mps = velocity_map.velocity_bins_Aps * 100.0
    mesh = ax.pcolormesh(
        bins_mps,
        bins_mps,
        velocity_map.counts.T,
        shading="auto",
        cmap="magma",
        norm=_color_norm(velocity_map.counts),
    )
    plt.colorbar(mesh, ax=ax, label="counts")
    ax.set_title(
        f"(b) simulated I+He VMI, {velocity_map.mass_amu:.0f} u, "
        f"{velocity_map.num_atoms_used} atoms"
    )
    ax.set_xlabel("v_x / m/s")
    ax.set_ylabel("v_y / m/s")
    ax.set_aspect("equal")
    ax.set_xlim(-3500.0, 3500.0)
    ax.set_ylim(-3500.0, 3500.0)
    ax.set_facecolor(plt.get_cmap("magma")(0.0))


def _draw_cov_heatmap(
    ax, counts, *, x, y, xlabel: str, ylabel: str, title: str
) -> None:
    counts_arr = np.asarray(counts, dtype=float)
    finite = counts_arr[np.isfinite(counts_arr)]
    vmax = float(np.nanmax(finite)) if finite.size else 0.0
    if vmax <= 0.0:
        norm = None
    else:
        norm = mcolors.Normalize(vmin=0.0, vmax=PAPER_COV_COLOR_CLIP_FRACTION * vmax)
    mesh = ax.pcolormesh(
        x, y, counts_arr, shading="auto", cmap="magma", norm=norm
    )
    plt.colorbar(mesh, ax=ax, label="covariance / arb. units")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_aspect("equal")


def _color_norm(values, *, noise_floor_fraction: float = 0.0):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return None
    vmax = float(np.nanmax(finite))
    if vmax <= 0.0:
        return None
    if noise_floor_fraction > 0.0:
        return mcolors.Normalize(
            vmin=float(noise_floor_fraction) * vmax,
            vmax=vmax,
        )
    return mcolors.Normalize(vmin=0.0, vmax=0.8 * vmax)


def _roll_to_centered_theta(
    counts: np.ndarray, theta_centers_2pi: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Roll a covariance matrix from [0, 2pi) bins into [-pi, pi) for plotting.

    The simulated angular pair covariance uses
    ``theta = arctan2(vx, vy) + pi`` wrapped into ``[0, 2 pi)``; the MATLAB
    experimental matrix uses ``[-pi, pi]``. Shifting both axes by half the
    bin count realigns the centre of the figure so the side-by-side
    comparison shares an origin.
    """
    n = theta_centers_2pi.size
    shift = n // 2
    rolled = np.roll(np.roll(counts, shift, axis=0), shift, axis=1)
    theta_centered = theta_centers_2pi - np.pi
    return rolled, theta_centered


if __name__ == "__main__":
    raise SystemExit(main())
