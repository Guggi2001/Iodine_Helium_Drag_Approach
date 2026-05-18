"""Paper-v2 port of post_process_single_pulse_paper_IplusHe_comparison.m.

The Python name ``paper_v2`` is a consistency alias. The legacy source is the
I+He comparison MATLAB script. The output is split into one figure per
diagnostic so each can be inspected independently:

- ``paper_v2_vmi_comparison.png``        experimental + simulated 2-D VMI
- ``paper_v2_radial_comparison.png``     radial velocity distribution
- ``paper_v2_he2_vmi_comparison.png``    experimental + simulated I+He2 VMI
- ``paper_v2_he2_radial_comparison.png`` I+He2 radial velocity distribution
- ``paper_v2_phi_comparison.png``        azimuthal phi distribution
- ``paper_v2_polar_image_comparison.png``  experimental + simulated polar VMI
                                          (only when a polar reference exists)

``EXPERIMENTAL_NOISE_FLOOR`` (or ``--noise-floor``) clips the colormap
floor on the experimental VMI panels to suppress the low-level noise band
so the bright signal lobes stand out against a dark background. Simulated
panels are unaffected.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# USER SETTINGS / DEFAULTS
# =============================================================================
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet"
REFERENCE_DIR = PROJECT_ROOT / "data" / "reference" / "paper_v2"
MASS_SELECTION_AMU = 131.0
MASS_SELECTION_AMU_2 = 135.0
# Fraction of the experimental panel's max intensity below which pixels
# collapse to the bottom (background) colormap entry. Raise to suppress
# more of the noise band; lower (or set to 0.0) to show more of the dim
# background. Only the experimental Cartesian and polar panels use this.
# (--noise-floor on the command line overrides this constant.)
EXPERIMENTAL_NOISE_FLOOR = 0.20


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    load_paper_v2_phi_reference,
    load_paper_v2_he2_radial_references,
    load_paper_v2_radial_references,
    load_paper_v2_vmi_image_reference,
    load_paper_v2_vmi_polar_image_reference,
    paper_v2_phi_curve,
    paper_v2_velocity_curve,
    paper_v2_velocity_map,
    polar_velocity_histogram,
)
from i2_helium_md.postprocess.paper_v2 import max_normalise  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Reproduce the active droplet branch of post_process_single_pulse_paper_IplusHe_comparison.m."
    )
    p.add_argument("--run-dir", type=Path, default=RUN_DIR,
                   help="RunDirectory containing cfg.json, neutral.npz, and ion.npz.")
    p.add_argument("--reference-dir", type=Path, default=REFERENCE_DIR,
                   help="Directory containing paper-v2 reference CSVs and images/.")
    p.add_argument("--no-show", action="store_true",
                   help="Skip plt.show(), useful for tests and headless runs.")
    p.add_argument("--noise-floor", type=float, default=EXPERIMENTAL_NOISE_FLOOR,
                   help="Fraction of max intensity below which experimental "
                        "VMI pixels clip to background (default: %(default)s; "
                        "set 0 to disable). Simulated panels are unaffected.")
    return p.parse_args(argv)


def _color_norm(values, *, noise_floor_fraction: float = 0.0):
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


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    run = RunDirectory(args.run_dir)
    ion = run.load_ion()

    radial_refs = load_paper_v2_radial_references(args.reference_dir)
    if not radial_refs:
        print(f"[paper_v2] no optional radial references found in {args.reference_dir}")
    he2_radial_refs = load_paper_v2_he2_radial_references(args.reference_dir)
    if not he2_radial_refs:
        print(f"[paper_v2] no optional I+He2 radial references found in {args.reference_dir}")

    image_ref = _load_optional_image(args.reference_dir)
    he2_image_ref = _load_optional_he2_image(args.reference_dir)
    polar_ref = _load_optional_polar_image(args.reference_dir)
    phi_ref = _load_optional_phi(args.reference_dir)
    velocity_map = paper_v2_velocity_map(ion, mass_amu=MASS_SELECTION_AMU)
    velocity_map_2 = paper_v2_velocity_map(ion, mass_amu=MASS_SELECTION_AMU_2)
    velocity_curve = paper_v2_velocity_curve(ion, mass_amu=MASS_SELECTION_AMU)
    velocity_curve_2 = paper_v2_velocity_curve(ion, mass_amu=MASS_SELECTION_AMU_2)
    phi_curve = paper_v2_phi_curve(ion, mass_amu=MASS_SELECTION_AMU)

    fig_vmi = _build_vmi_figure(
        image_ref=image_ref,
        velocity_map=velocity_map,
        experimental_noise_floor=args.noise_floor,
    )
    fig_vmi_2 = _build_vmi_figure(
        image_ref=he2_image_ref,
        velocity_map=velocity_map_2,
        experimental_noise_floor=args.noise_floor,
    )
    fig_radial = _build_radial_figure(radial_refs, velocity_curve)
    fig_radial_2 = _build_radial_figure(he2_radial_refs, velocity_curve_2)
    fig_phi = _build_phi_figure(phi_curve, phi_ref=phi_ref)
    fig_polar = None
    if polar_ref is not None:
        polar_hist = _polar_histogram_matched_to_reference(ion, polar_ref)
        fig_polar = _build_polar_image_figure(
            polar_ref, polar_hist, experimental_noise_floor=args.noise_floor
        )

    out_dir = run.root / "figures"
    out_dir.mkdir(exist_ok=True)
    fig_vmi.savefig(out_dir / "paper_v2_vmi_comparison.png", dpi=150)
    fig_vmi_2.savefig(out_dir / "paper_v2_he2_vmi_comparison.png", dpi=150)
    fig_radial.savefig(out_dir / "paper_v2_radial_comparison.png", dpi=150)
    fig_radial_2.savefig(out_dir / "paper_v2_he2_radial_comparison.png", dpi=150)
    fig_phi.savefig(out_dir / "paper_v2_phi_comparison.png", dpi=150)
    print(f"Saved VMI comparison to {out_dir / 'paper_v2_vmi_comparison.png'}")
    print(f"Saved I+He2 VMI comparison to {out_dir / 'paper_v2_he2_vmi_comparison.png'}")
    print(f"Saved radial comparison to {out_dir / 'paper_v2_radial_comparison.png'}")
    print(f"Saved I+He2 radial comparison to {out_dir / 'paper_v2_he2_radial_comparison.png'}")
    print(f"Saved phi comparison to {out_dir / 'paper_v2_phi_comparison.png'}")
    if fig_polar is not None:
        fig_polar.savefig(out_dir / "paper_v2_polar_image_comparison.png", dpi=150)
        print(f"Saved polar image comparison to {out_dir / 'paper_v2_polar_image_comparison.png'}")

    if not args.no_show:
        plt.show()
    return 0


def _load_optional_image(reference_dir: Path):
    return _load_optional_image_stem(
        reference_dir,
        "iplus_he_high_snr_vmi_image",
        "I+He 2-D VMI image",
        load_paper_v2_vmi_image_reference,
    )


def _load_optional_he2_image(reference_dir: Path):
    return _load_optional_image_stem(
        reference_dir,
        "iplus_he2_high_snr_vmi_image",
        "I+He2 2-D VMI image",
        load_paper_v2_vmi_image_reference,
    )


def _load_optional_image_stem(reference_dir: Path, stem: str, description: str, loader):
    image_dir = reference_dir / "images"
    image_paths = [image_dir / f"{stem}.npz", image_dir / f"{stem}.mat"]
    for path in image_paths:
        if path.exists():
            return loader(path)
    print(f"[paper_v2] optional {description} reference not found in {image_dir}")
    return None


def _load_optional_polar_image(reference_dir: Path):
    return _load_optional_image_stem(
        reference_dir,
        "iplus_he_high_snr_vmi_polar_image",
        "I+He polar VMI image",
        load_paper_v2_vmi_polar_image_reference,
    )


def _polar_histogram_matched_to_reference(ion, polar_ref):
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
        mass_amu=MASS_SELECTION_AMU,
    )


def _load_optional_phi(reference_dir: Path):
    phi_path = reference_dir / "iplus_he_high_snr_phi.csv"
    if not phi_path.exists():
        print(f"[paper_v2] optional phi reference not found: {phi_path}")
        return None
    return load_paper_v2_phi_reference(phi_path)


def _build_vmi_figure(
    *, image_ref, velocity_map, experimental_noise_floor: float = 0.0
) -> plt.Figure:
    fig, (ax_exp, ax_sim) = plt.subplots(
        1, 2, figsize=(10.0, 4.8), constrained_layout=True
    )
    _draw_experimental_image(
        ax_exp, image_ref, noise_floor_fraction=experimental_noise_floor
    )
    _draw_simulated_map(ax_sim, velocity_map)
    return fig


def _build_radial_figure(radial_refs, velocity_curve) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
    _draw_radial_panel(ax, radial_refs, velocity_curve)
    return fig


def _draw_experimental_image(
    ax, image_ref, *, noise_floor_fraction: float = 0.0
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
            norm=_color_norm(
                image_ref.intensity, noise_floor_fraction=noise_floor_fraction
            ),
        )
        plt.colorbar(mesh, ax=ax, label="signal / arb. units")
    ax.set_title(_experimental_image_title(image_ref))
    ax.set_xlabel("v_x / m/s")
    ax.set_ylabel("v_y / m/s")
    ax.set_aspect("equal")
    ax.set_xlim(-3500.0, 3500.0)
    ax.set_ylim(-3500.0, 3500.0)
    ax.set_facecolor(cmap_floor)


def _draw_simulated_map(ax, velocity_map) -> None:
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
    label = _simulation_channel_label(velocity_map.mass_amu)
    if label.endswith(" u"):
        title = f"(b) simulated VMI map, {label}"
    else:
        title = f"(b) simulated {label} VMI map, {velocity_map.mass_amu:.0f} u"
    ax.set_title(title)
    ax.set_xlabel("v_x / m/s")
    ax.set_ylabel("v_y / m/s")
    ax.set_aspect("equal")
    ax.set_xlim(-3500.0, 3500.0)
    ax.set_ylim(-3500.0, 3500.0)
    ax.set_facecolor(plt.get_cmap("magma")(0.0))


def _draw_radial_panel(ax, radial_refs, velocity_curve) -> None:
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
        label=_simulation_curve_label(velocity_curve.mass_amu),
    )
    ax.set_title("2-D detector-plane speed vs raw VMI radial profile")
    ax.set_xlabel("v / m/s")
    ax.set_ylabel("signal / arb. units")
    ax.set_xlim(0.0, 3500.0)
    ax.set_ylim(0.0, 1.1)
    ax.legend(frameon=False, fontsize=8)


def _experimental_image_title(image_ref) -> str:
    if image_ref is not None:
        name = image_ref.source_path.name.lower()
        if "he2" in name:
            return "(a) experimental I+He2 high-SNR VMI image"
        if "he_high" in name:
            return "(a) experimental I+He high-SNR VMI image"
    return "(a) experimental high-SNR VMI image"


def _simulation_channel_label(mass_amu: float) -> str:
    if np.isclose(mass_amu, 135.0):
        return "I+He2"
    if np.isclose(mass_amu, 131.0):
        return "I+He"
    return f"{mass_amu:.0f} u"


def _simulation_curve_label(mass_amu: float) -> str:
    label = _simulation_channel_label(mass_amu)
    if label.endswith(" u"):
        return f"simulation {label}"
    return f"simulation {label} {mass_amu:.0f} u"


def _build_phi_figure(phi_curve, *, phi_ref) -> plt.Figure:
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


def _build_polar_image_figure(
    polar_ref, polar_hist, *, experimental_noise_floor: float = 0.0
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
        norm=_color_norm(
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
        norm=_color_norm(polar_hist.counts),
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


if __name__ == "__main__":
    raise SystemExit(main())
