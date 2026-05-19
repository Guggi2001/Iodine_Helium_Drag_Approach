"""Shared plotting helpers for the paper-cov comparison panels."""

from __future__ import annotations

from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

from .paper_cov import (
    PAPER_COV_COLOR_CLIP_FRACTION,
    covariance_axis_sum_normalised,
    load_paper_cov_experimental_reference,
    radial_covariance_trace,
    simulated_phi_distribution,
)
from .paper_v2 import (
    load_paper_v2_phi_reference,
    load_paper_v2_radial_reference,
    load_paper_v2_vmi_image_reference,
    max_normalise,
)


HIGH_SNR_RADIAL_FILENAME = "iplus_he_high_snr_radial.csv"
EXPERIMENTAL_PHI_FILENAME = "iplus_he_phi.csv"


def load_optional_high_snr_radial(
    reference_dir: Path, *, log_prefix: str = "[paper_cov]",
):
    path = reference_dir / HIGH_SNR_RADIAL_FILENAME
    if not path.exists():
        print(f"{log_prefix} optional high-SNR I+He radial reference not found: {path}")
        return None
    return load_paper_v2_radial_reference(path)


def load_optional_vmi_image(reference_dir: Path, *, log_prefix: str = "[paper_cov]"):
    image_dir = reference_dir / "images"
    for ext in (".npz", ".mat"):
        path = image_dir / f"iplus_he_high_snr_vmi_image{ext}"
        if path.exists():
            return load_paper_v2_vmi_image_reference(path)
    print(f"{log_prefix} optional 2-D VMI image reference not found in {image_dir}")
    return None


def load_optional_cov_reference(
    reference_dir: Path, *, log_prefix: str = "[paper_cov]",
):
    for filename in ("iplus_he_covariance.mat", "iplus_he_covariance.npz"):
        path = reference_dir / filename
        if path.exists():
            return load_paper_cov_experimental_reference(path)
    print(
        f"{log_prefix} experimental covariance reference not found in {reference_dir}; "
        "run data/reference/scripts/export_paper_cov_reference_data.m to generate it."
    )
    return None


def load_optional_phi_reference(
    reference_dir: Path, *, log_prefix: str = "[paper_cov]",
):
    path = reference_dir / EXPERIMENTAL_PHI_FILENAME
    if not path.exists():
        print(f"{log_prefix} experimental phi reference not found: {path}")
        return None
    return load_paper_v2_phi_reference(path)


def build_vmi_figure(
    *, image_ref, velocity_map, experimental_noise_floor: float
) -> plt.Figure:
    fig, (ax_exp, ax_sim) = plt.subplots(
        1, 2, figsize=(10.0, 4.8), constrained_layout=True
    )
    draw_experimental_vmi(
        ax_exp, image_ref, noise_floor_fraction=experimental_noise_floor
    )
    draw_simulated_vmi(ax_sim, velocity_map)
    return fig


def build_radial_distribution_figure(
    *, high_snr_ref, sim_radial_curve, sim_radial, cov_ref, title = 'standard',
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
    if title == 'standard':
        ax.set_title("Velocity distribution comparison")
    elif title == 'run_summary':
        ax.set_title('2-D detector plane speed using experimental VMI data \n and simulation with corresponding v-cov trace')
    ax.legend(handles, labels, fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)
    return fig


def build_angular_cov_figure(*, cov_ref, sim_angular) -> plt.Figure:
    fig, (ax_exp, ax_sim) = plt.subplots(
        1, 2, figsize=(10.0, 4.8), constrained_layout=True
    )
    draw_cov_heatmap(
        ax_exp,
        cov_ref.cov_angular,
        x=cov_ref.theta_centers_rad,
        y=cov_ref.theta_centers_rad,
        xlabel="angle / rad",
        ylabel="angle / rad",
        title="(a) Experimental angular pair covariance",
    )
    counts_centered, theta_centered = roll_to_centered_theta(
        np.asarray(sim_angular.counts, dtype=float),
        sim_angular.theta_centers_rad,
    )
    draw_cov_heatmap(
        ax_sim,
        counts_centered,
        x=theta_centered,
        y=theta_centered,
        xlabel="angle / rad",
        ylabel="angle / rad",
        title=(
            f"(b) Simulated angular pair covariance, "
            f"{sim_angular.num_pairs_used} pairs"
        ),
    )
    return fig


def build_radial_cov_figure(*, cov_ref, sim_radial) -> plt.Figure:
    fig, (ax_exp, ax_sim) = plt.subplots(
        1, 2, figsize=(10.0, 4.8), constrained_layout=True
    )
    draw_cov_heatmap(
        ax_exp,
        cov_ref.cov_radial,
        x=cov_ref.velocity_centers_Aps,
        y=cov_ref.velocity_centers_Aps,
        xlabel=r"v / $\AA$/ps",
        ylabel=r"v / $\AA$/ps",
        title="(c) Experimental radial pair covariance",
    )
    draw_cov_heatmap(
        ax_sim,
        sim_radial.counts,
        x=sim_radial.velocity_centers_Aps,
        y=sim_radial.velocity_centers_Aps,
        xlabel=r"v / $\AA$/ps",
        ylabel=r"v / $\AA$/ps",
        title=(
            f"(d) Simulated radial pair covariance, "
            f"{sim_radial.num_pairs_used} pairs"
        ),
    )
    return fig


def build_phi_distribution_figure(*, phi_ref, ion, mass_amu: float) -> plt.Figure:
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
        labels.append("I+He high SNR")

    sim_phi = simulated_phi_distribution(ion, mass_amu=mass_amu)
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
    ax.set_title("Phi angular distribution")
    ax.legend(handles, labels, fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)
    return fig


def build_pair_cov_traces_figure(*, cov_ref, sim_angular, sim_radial) -> plt.Figure:
    fig, (ax_ang, ax_rad) = plt.subplots(
        1, 2, figsize=(10.0, 4.5), constrained_layout=True
    )

    sim_ang_centered, theta_centered = roll_to_centered_theta(
        np.asarray(sim_angular.counts, dtype=float),
        sim_angular.theta_centers_rad,
    )
    sim_ang_trace = covariance_axis_sum_normalised(sim_ang_centered, axis=0)
    exp_ang_trace = covariance_axis_sum_normalised(cov_ref.cov_angular, axis=0)
    ax_ang.plot(theta_centered, sim_ang_trace, lw=1.4, color="tab:orange", label="sim")
    ax_ang.plot(
        cov_ref.theta_centers_rad,
        exp_ang_trace,
        lw=1.2,
        ls="--",
        color="black",
        label="exp",
    )
    ax_ang.set_xlabel(r"$\theta$ / radian")
    ax_ang.set_ylabel("normalised trace / arb. units")
    ax_ang.set_xlim(-np.pi, np.pi)
    ax_ang.set_ylim(0.0, 1.05)
    ax_ang.set_title("(e) Angular pair-cov trace")
    ax_ang.legend(fontsize=9, loc="upper right")
    ax_ang.grid(True, alpha=0.3)

    sim_rad_trace = covariance_axis_sum_normalised(sim_radial.counts, axis=0)
    exp_rad_trace = covariance_axis_sum_normalised(cov_ref.cov_radial, axis=0)
    ax_rad.plot(
        sim_radial.velocity_centers_Aps,
        sim_rad_trace,
        lw=1.4,
        color="tab:orange",
        label="sim",
    )
    ax_rad.plot(
        cov_ref.velocity_centers_Aps,
        exp_rad_trace,
        lw=1.2,
        ls="--",
        color="black",
        label="exp",
    )
    ax_rad.set_xlabel(r"$v_r$ / $\AA$/ps")
    ax_rad.set_ylabel("normalised trace / arb. units")
    ax_rad.set_ylim(0.0, 1.05)
    ax_rad.set_title("(f) Radial pair-cov trace")
    ax_rad.legend(fontsize=9, loc="upper right")
    ax_rad.grid(True, alpha=0.3)
    return fig


def draw_experimental_vmi(ax, image_ref, *, noise_floor_fraction: float) -> None:
    cmap_floor = plt.get_cmap("magma")(0.0)
    if image_ref is None:
        ax.text(
            0.5,
            0.5,
            "experimental VMI image not exported",
            ha="center",
            va="center",
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
            norm=color_norm(
                image_ref.intensity, noise_floor_fraction=noise_floor_fraction
            ),
        )
        plt.colorbar(mesh, ax=ax, label="signal / arb. units")
    ax.set_title("(a) Experimental I+He VMI")
    ax.set_xlabel("v_x / m/s")
    ax.set_ylabel("v_y / m/s")
    ax.set_aspect("equal")
    ax.set_xlim(-3500.0, 3500.0)
    ax.set_ylim(-3500.0, 3500.0)
    ax.set_facecolor(cmap_floor)


def draw_simulated_vmi(ax, velocity_map) -> None:
    bins_mps = velocity_map.velocity_bins_Aps * 100.0
    mesh = ax.pcolormesh(
        bins_mps,
        bins_mps,
        velocity_map.counts.T,
        shading="auto",
        cmap="magma",
        norm=color_norm(velocity_map.counts),
    )
    plt.colorbar(mesh, ax=ax, label="counts")
    ax.set_title(
        f"(b) Simulated I+He VMI, {velocity_map.mass_amu:.0f} u, "
        f"{velocity_map.num_atoms_used} atoms"
    )
    ax.set_xlabel("v_x / m/s")
    ax.set_ylabel("v_y / m/s")
    ax.set_aspect("equal")
    ax.set_xlim(-3500.0, 3500.0)
    ax.set_ylim(-3500.0, 3500.0)
    ax.set_facecolor(plt.get_cmap("magma")(0.0))


def draw_cov_heatmap(
    ax, counts, *, x, y, xlabel: str, ylabel: str, title: str
) -> None:
    counts_arr = np.asarray(counts, dtype=float)
    finite = counts_arr[np.isfinite(counts_arr)]
    vmax = float(np.nanmax(finite)) if finite.size else 0.0
    norm = None
    if vmax > 0.0:
        norm = mcolors.Normalize(vmin=0.0, vmax=PAPER_COV_COLOR_CLIP_FRACTION * vmax)
    mesh = ax.pcolormesh(x, y, counts_arr, shading="auto", cmap="magma", norm=norm)
    plt.colorbar(mesh, ax=ax, label="covariance / arb. units")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_aspect("equal")


def color_norm(values, *, noise_floor_fraction: float = 0.0):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return None
    vmax = float(np.nanmax(finite))
    if vmax <= 0.0:
        return None
    if noise_floor_fraction > 0.0:
        return mcolors.Normalize(vmin=float(noise_floor_fraction) * vmax, vmax=vmax)
    return mcolors.Normalize(vmin=0.0, vmax=0.8 * vmax)


def roll_to_centered_theta(
    counts: np.ndarray, theta_centers_2pi: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    n = theta_centers_2pi.size
    shift = n // 2
    rolled = np.roll(np.roll(counts, shift, axis=0), shift, axis=1)
    theta_centered = theta_centers_2pi - np.pi
    return rolled, theta_centered
