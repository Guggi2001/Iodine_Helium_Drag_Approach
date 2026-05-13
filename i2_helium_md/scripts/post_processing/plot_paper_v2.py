"""Paper-v2 port of post_process_single_pulse_paper_IplusHe_comparison.m.

The Python name ``paper_v2`` is a consistency alias. The legacy source is the
I+He comparison MATLAB script, whose active droplet branch combines a processed
experimental 2-D VMI image, a simulated vx/vy velocity map, a radial VMI
comparison, and a separate phi figure.
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


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    load_paper_v2_phi_reference,
    load_paper_v2_radial_references,
    load_paper_v2_vmi_image_reference,
    paper_v2_phi_curve,
    paper_v2_velocity_curve,
    paper_v2_velocity_map,
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
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    run = RunDirectory(args.run_dir)
    ion = run.load_ion()

    radial_refs = load_paper_v2_radial_references(args.reference_dir)
    if not radial_refs:
        print(f"[paper_v2] no optional radial references found in {args.reference_dir}")

    image_ref = _load_optional_image(args.reference_dir)
    phi_ref = _load_optional_phi(args.reference_dir)
    velocity_map = paper_v2_velocity_map(ion, mass_amu=MASS_SELECTION_AMU)
    velocity_curve = paper_v2_velocity_curve(ion, mass_amu=MASS_SELECTION_AMU)
    phi_curve = paper_v2_phi_curve(ion, mass_amu=MASS_SELECTION_AMU)

    fig_main = _build_main_figure(
        image_ref=image_ref,
        velocity_map=velocity_map,
        radial_refs=radial_refs,
        velocity_curve=velocity_curve,
    )
    fig_phi = _build_phi_figure(phi_curve, phi_ref=phi_ref)

    out_dir = run.root / "figures"
    out_dir.mkdir(exist_ok=True)
    fig_main.savefig(out_dir / "paper_v2_compare_simulation_and_measurement.png", dpi=150)
    fig_phi.savefig(out_dir / "paper_v2_phi_comparison.png", dpi=150)
    print(f"Saved figure to {out_dir / 'paper_v2_compare_simulation_and_measurement.png'}")
    print(f"Saved phi comparison to {out_dir / 'paper_v2_phi_comparison.png'}")

    if not args.no_show:
        plt.show()
    return 0


def _load_optional_image(reference_dir: Path):
    image_dir = reference_dir / "images"
    image_paths = []
    if image_dir.exists():
        image_paths = sorted(image_dir.glob("*_vmi_image.npz"))
        image_paths.extend(sorted(image_dir.glob("*_vmi_image.mat")))
    if not image_paths:
        print(f"[paper_v2] no optional 2-D VMI image references found in {image_dir}")
        return None
    if len(image_paths) > 1:
        print(f"[paper_v2] using first 2-D VMI image reference: {image_paths[0].name}")
    return load_paper_v2_vmi_image_reference(image_paths[0])


def _load_optional_phi(reference_dir: Path):
    phi_path = reference_dir / "iplus_he_high_snr_phi.csv"
    if not phi_path.exists():
        print(f"[paper_v2] optional phi reference not found: {phi_path}")
        return None
    return load_paper_v2_phi_reference(phi_path)


def _build_main_figure(*, image_ref, velocity_map, radial_refs, velocity_curve) -> plt.Figure:
    fig = plt.figure(figsize=(8.0, 8.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.8])
    ax_exp = fig.add_subplot(gs[0, 0])
    ax_sim = fig.add_subplot(gs[0, 1])
    ax_radial = fig.add_subplot(gs[1, :])

    _draw_experimental_image(ax_exp, image_ref)
    _draw_simulated_map(ax_sim, velocity_map)
    _draw_radial_panel(ax_radial, radial_refs, velocity_curve)
    return fig


def _draw_experimental_image(ax, image_ref) -> None:
    if image_ref is None:
        ax.text(0.5, 0.5, "experimental VMI image not exported", ha="center", va="center")
        ax.set_xlim(-35.0, 35.0)
        ax.set_ylim(-35.0, 35.0)
    else:
        vmax = float(np.nanmax(image_ref.intensity)) * 0.8
        mesh = ax.pcolormesh(
            image_ref.vx_Aps,
            image_ref.vy_Aps,
            image_ref.intensity,
            shading="auto",
            cmap="magma",
            vmin=0.0,
            vmax=vmax if vmax > 0.0 else None,
        )
        plt.colorbar(mesh, ax=ax, label="signal / arb. units")
    ax.set_title("(a) experimental VMI image")
    ax.set_xlabel("v_x / A/ps")
    ax.set_ylabel("v_y / A/ps")
    ax.set_aspect("equal")
    ax.set_xlim(-35.0, 35.0)
    ax.set_ylim(-35.0, 35.0)


def _draw_simulated_map(ax, velocity_map) -> None:
    vmax = float(np.nanmax(velocity_map.counts)) * 0.8
    mesh = ax.pcolormesh(
        velocity_map.velocity_bins_Aps,
        velocity_map.velocity_bins_Aps,
        velocity_map.counts.T,
        shading="auto",
        cmap="magma",
        vmin=0.0,
        vmax=vmax if vmax > 0.0 else None,
    )
    plt.colorbar(mesh, ax=ax, label="counts")
    ax.set_title(f"(b) simulated VMI map, {velocity_map.mass_amu:.0f} u")
    ax.set_xlabel("v_x / A/ps")
    ax.set_ylabel("v_y / A/ps")
    ax.set_aspect("equal")
    ax.set_xlim(-35.0, 35.0)
    ax.set_ylim(-35.0, 35.0)


def _draw_radial_panel(ax, radial_refs, velocity_curve) -> None:
    for ref in radial_refs:
        ax.plot(
            ref.velocity_Aps,
            max_normalise(ref.signal_arb),
            linewidth=1.4,
            label=ref.label,
        )
    ax.plot(
        velocity_curve.bin_centers_Aps,
        velocity_curve.normalised,
        "--",
        color="black",
        linewidth=1.5,
        label=f"simulation I+He {velocity_curve.mass_amu:.0f} u",
    )
    ax.set_title("(c) radial velocity distribution")
    ax.set_xlabel("v / A/ps")
    ax.set_ylabel("signal / arb. units")
    ax.set_xlim(0.0, 35.0)
    ax.set_ylim(0.0, 1.1)
    ax.legend(frameon=False, fontsize=8)


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


if __name__ == "__main__":
    raise SystemExit(main())
