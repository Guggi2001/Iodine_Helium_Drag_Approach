"""Shared plotting helpers for the paper-v2 panels.

These helpers were lifted from
``scripts/post_processing/plot_paper_v2.py`` so the same builders can be
reused by ``plot_run_summary.py`` without duplicating ~200 lines of
figure-construction code. The functions are kept thin and matplotlib-only;
all numerical work lives in :mod:`i2_helium_md.postprocess.paper_v2`.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

from i2_helium_md.postprocess import (
    load_paper_v2_phi_reference,
    load_paper_v2_vmi_image_reference,
    load_paper_v2_vmi_polar_image_reference,
    polar_velocity_histogram,
)
from i2_helium_md.postprocess.paper_v2 import max_normalise


def color_norm(values, *, noise_floor_fraction: float = 0.0):
    """Build a Normalize covering the panel's intensity range.

    ``noise_floor_fraction == 0`` (default): legacy contrast,
    ``Normalize(0, 0.8 * vmax)``. Used by simulated panels.

    ``noise_floor_fraction > 0``: ``Normalize(floor * vmax, vmax)``.
    Pixels below the floor clip to the bottom colormap entry; the full
    bright range stays visible at the top. Used by the experimental
    panels to suppress the high-SNR noise band.

    Returns ``None`` when ``vmax <= 0`` so callers fall back to default
    matplotlib autoscaling.
    """
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


def load_optional_image(reference_dir: Path, *, log_prefix: str = "[paper_v2]"):
    image_dir = reference_dir / "images"
    image_paths = []
    if image_dir.exists():
        image_paths = sorted(image_dir.glob("*_vmi_image.npz"))
        image_paths.extend(sorted(image_dir.glob("*_vmi_image.mat")))
    if not image_paths:
        print(f"{log_prefix} no optional 2-D VMI image references found in {image_dir}")
        return None
    if len(image_paths) > 1:
        print(f"{log_prefix} using first 2-D VMI image reference: {image_paths[0].name}")
    return load_paper_v2_vmi_image_reference(image_paths[0])


def load_optional_polar_image(
    reference_dir: Path, *, log_prefix: str = "[paper_v2]",
):
    image_dir = reference_dir / "images"
    image_paths = []
    if image_dir.exists():
        image_paths = sorted(image_dir.glob("*_vmi_polar_image.npz"))
        image_paths.extend(sorted(image_dir.glob("*_vmi_polar_image.mat")))
    if not image_paths:
        print(f"{log_prefix} no optional polar VMI image references found in {image_dir}")
        return None
    if len(image_paths) > 1:
        print(f"{log_prefix} using first polar VMI image reference: {image_paths[0].name}")
    return load_paper_v2_vmi_polar_image_reference(image_paths[0])


def polar_histogram_matched_to_reference(ion, polar_ref, *, mass_amu: float = 131.0):
    """Bin simulated final velocities so the (phi, v) axes align with the
    experimental polar reference. ``polar_velocity_histogram`` builds uniform
    edges via ``np.linspace(0, v_max, n+1)`` and centers at the midpoints, so
    a half-bin shift relative to the experimental sample points is expected
    when the reference axis is not center-aligned. This is left unsmoothed
    by design (no 2-D interpolation per the plan).
    """
    n_phi = int(polar_ref.phi_rad.size)
    n_v = int(polar_ref.v_radius_Aps.size)
    if n_v < 2:
        raise ValueError("polar reference must have at least two v_radius bins")
    dv = float(polar_ref.v_radius_Aps[1] - polar_ref.v_radius_Aps[0])
    v_max = float(polar_ref.v_radius_Aps[-1] + dv)
    return polar_velocity_histogram(
        ion,
        n_v_bins=n_v,
        n_phi_bins=n_phi,
        v_max_Aps=v_max,
        mass_amu=mass_amu,
    )


def load_optional_phi(reference_dir: Path, *, log_prefix: str = "[paper_v2]"):
    phi_path = reference_dir / "iplus_he_high_snr_phi.csv"
    if not phi_path.exists():
        print(f"{log_prefix} optional phi reference not found: {phi_path}")
        return None
    return load_paper_v2_phi_reference(phi_path)


def build_vmi_figure(
    *, image_ref, velocity_map, experimental_noise_floor: float = 0.0,
) -> plt.Figure:
    fig, (ax_exp, ax_sim) = plt.subplots(
        1, 2, figsize=(10.0, 4.8), constrained_layout=True
    )
    draw_experimental_image(
        ax_exp, image_ref, noise_floor_fraction=experimental_noise_floor
    )
    draw_simulated_map(ax_sim, velocity_map)
    return fig


def build_radial_figure(radial_refs, velocity_curve) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
    draw_radial_panel(ax, radial_refs, velocity_curve)
    return fig


def draw_experimental_image(
    ax, image_ref, *, noise_floor_fraction: float = 0.0,
) -> None:
    cmap_floor = plt.get_cmap("magma")(0.0)
    if image_ref is None:
        ax.text(0.5, 0.5, "experimental VMI image not exported", ha="center", va="center")
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
    ax.set_title("(a) experimental high-SNR VMI image")
    ax.set_xlabel("v_x / m/s")
    ax.set_ylabel("v_y / m/s")
    ax.set_aspect("equal")
    ax.set_xlim(-3500.0, 3500.0)
    ax.set_ylim(-3500.0, 3500.0)
    ax.set_facecolor(cmap_floor)


def draw_simulated_map(ax, velocity_map) -> None:
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
    ax.set_title(f"(b) simulated VMI map, {velocity_map.mass_amu:.0f} u")
    ax.set_xlabel("v_x / m/s")
    ax.set_ylabel("v_y / m/s")
    ax.set_aspect("equal")
    ax.set_xlim(-3500.0, 3500.0)
    ax.set_ylim(-3500.0, 3500.0)
    ax.set_facecolor(plt.get_cmap("magma")(0.0))


def draw_radial_panel(ax, radial_refs, velocity_curve) -> None:
    for ref in radial_refs:
        ax.plot(
            ref.velocity_mps,
            max_normalise(ref.signal_arb),
            linewidth=1.4,
            label=ref.label,
        )
    ax.plot(
        velocity_curve.bin_centers_mps,
        velocity_curve.normalised,
        "--",
        color="black",
        linewidth=1.5,
        label=f"simulation I+He {velocity_curve.mass_amu:.0f} u",
    )
    ax.set_title("2-D detector-plane speed vs raw VMI radial profile")
    ax.set_xlabel("v / m/s")
    ax.set_ylabel("signal / arb. units")
    ax.set_xlim(0.0, 3500.0)
    ax.set_ylim(0.0, 1.1)
    ax.legend(frameon=False, fontsize=8)


def build_phi_figure(phi_curve, *, phi_ref) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
    if phi_ref is not None:
        ax.plot(
            phi_ref.phi_rad,
            max_normalise(phi_ref.signal_arb),
            linewidth=1.5,
            label="I+He high-SNR",
        )
    ax.plot(
        phi_curve.bin_centers_rad,
        phi_curve.normalised,
        linewidth=1.5,
        label=f"simulation I+He {phi_curve.mass_amu:.0f} u",
    )
    ax.set_title("paper v2 azimuthal distribution")
    ax.set_xlabel("phi / radian")
    ax.set_ylabel("signal / arb. units")
    ax.set_xlim(0.0, 2.0 * np.pi)
    ax.set_ylim(0.0, 1.1)
    ax.legend(frameon=False, fontsize=9)
    return fig


def build_polar_image_figure(
    polar_ref, polar_hist, *, experimental_noise_floor: float = 0.0,
) -> plt.Figure:
    fig, (ax_exp, ax_sim) = plt.subplots(
        1, 2, figsize=(10.0, 4.5), constrained_layout=True
    )

    mesh_e = ax_exp.pcolormesh(
        polar_ref.phi_rad,
        polar_ref.v_radius_mps,
        polar_ref.intensity.T,
        shading="auto",
        cmap="magma",
        norm=color_norm(
            polar_ref.intensity, noise_floor_fraction=experimental_noise_floor
        ),
    )
    plt.colorbar(mesh_e, ax=ax_exp, label="signal / arb. units")

    mesh_s = ax_sim.pcolormesh(
        polar_hist.phi_centers_rad,
        polar_hist.v_centers_Aps * 100.0,
        polar_hist.counts,
        shading="auto",
        cmap="magma",
        norm=color_norm(polar_hist.counts),
    )
    plt.colorbar(mesh_s, ax=ax_sim, label="counts")

    v_max_mps = float(polar_ref.v_radius_mps[-1])
    panels = (
        (ax_exp, "(a) experimental polar VMI"),
        (
            ax_sim,
            f"(b) simulated polar histogram, "
            f"{polar_hist.mass_amu:.0f} u (3-D |v|)",
        ),
    )
    for ax, title in panels:
        ax.set_xlabel("phi / radian")
        ax.set_ylabel("v / m/s")
        ax.set_xlim(0.0, 2.0 * np.pi)
        ax.set_ylim(0.0, v_max_mps)
        ax.set_title(title)
        ax.set_facecolor(plt.get_cmap("magma")(0.0))
    return fig
