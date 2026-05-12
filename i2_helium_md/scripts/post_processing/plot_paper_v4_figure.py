"""Literal droplet-branch port of post_process_single_pulse_paper_v4.m."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# USER SETTINGS / DEFAULTS
# =============================================================================
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet"
REFERENCE_DIR = PROJECT_ROOT / "data" / "reference" / "paper_v4"
MASS_SELECTIONS_AMU = (127.0, 131.0)


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    load_paper_v4_radial_references,
    mass_spectrum,
    paper_v4_angular_pair_covariance,
    paper_v4_velocity_curve,
)
from i2_helium_md.postprocess.paper_v4 import max_normalise  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Reproduce the active droplet branch of post_process_single_pulse_paper_v4.m."
    )
    p.add_argument("--run-dir", type=Path, default=RUN_DIR,
                   help="RunDirectory containing cfg.json, neutral.npz, and ion.npz.")
    p.add_argument("--reference-dir", type=Path, default=REFERENCE_DIR,
                   help="Directory containing v4 *_radial.csv reference exports.")
    p.add_argument("--no-show", action="store_true",
                   help="Skip plt.show(), useful for tests and headless runs.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    run = RunDirectory(args.run_dir)
    ion = run.load_ion()

    references = load_paper_v4_radial_references(args.reference_dir, only_include_iplus_he = True)
    if not references:
        print(f"[paper_v4] no optional radial references found in {args.reference_dir}")

    velocity_curves = []
    for mass in MASS_SELECTIONS_AMU:
        try:
            velocity_curves.append(paper_v4_velocity_curve(ion, mass_amu=mass))
        except ValueError as exc:
            print(f"[paper_v4] skip velocity mass {mass:.0f}: {exc}")

    covariance = paper_v4_angular_pair_covariance(ion, mass_amu=131.0)

    fig_main = _build_radial_figure(references, velocity_curves)
    fig_cov = _build_covariance_figure(covariance)
    fig_mass = _build_mass_figure(ion)

    out_dir = run.root / "figures"
    out_dir.mkdir(exist_ok=True)
    fig_main.savefig(out_dir / "compare_simulation_and_measurement_simpler.png", dpi=150)
    fig_cov.savefig(out_dir / "paper_v4_angular_pair_covariance.png", dpi=150)
    fig_mass.savefig(out_dir / "paper_v4_ion_mass_histogram.png", dpi=150)
    print(f"Saved figure to {out_dir / 'compare_simulation_and_measurement_simpler.png'}")
    print(f"Saved covariance to {out_dir / 'paper_v4_angular_pair_covariance.png'}")
    print(f"Saved mass histogram to {out_dir / 'paper_v4_ion_mass_histogram.png'}")

    if not args.no_show:
        plt.show()
    return 0


def _build_radial_figure(references, velocity_curves) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7.0, 4.6), constrained_layout=True)
    for ref in references:
        ax.plot(
            ref.velocity_mps,
            max_normalise(ref.signal_arb),
            linewidth=1.4,
            label=ref.label,
        )

    styles = {
        127.0: ("--", "simulated, mass = 127 u"),
        131.0: (":", "simulated, mass = 131 u"),
    }
    for curve in velocity_curves:
        linestyle, label = styles.get(curve.mass_amu, ("--", f"simulated, mass = {curve.mass_amu:.0f} u"))
        ax.plot(
            curve.bin_centers_mps,
            curve.normalised,
            linestyle=linestyle,
            linewidth=1.5,
            color=(0.1, 0.1, 0.1),
            label=label,
        )

    ax.set_xlabel("v / m/s")
    ax.set_ylabel("signal / arb. units")
    ax.set_xlim(0.0, 3500.0)
    ax.set_ylim(0.0, 1.1)
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("paper v4 radial velocity comparison")
    return fig


def _build_covariance_figure(covariance) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(6.4, 5.6), constrained_layout=True)
    theta_edges = covariance.theta_edges_rad - 3.141592653589793
    pcm = ax.pcolormesh(
        theta_edges,
        theta_edges,
        covariance.counts.T,
        shading="auto",
        cmap="magma",
    )
    if covariance.theta_pairs_rad.size:
        ax.scatter(
            covariance.theta_pairs_rad[:, 0] - 3.141592653589793,
            covariance.theta_pairs_rad[:, 1] - 3.141592653589793,
            marker="x",
            s=10,
            linewidths=0.5,
            color="black",
        )
    ax.set_aspect("equal")
    ax.set_xlabel("theta / radian")
    ax.set_ylabel("theta / radian")
    ax.set_title(f"paper v4 angular pair covariance (n_pairs={covariance.num_pairs_used})")
    fig.colorbar(pcm, ax=ax, label="counts")
    return fig


def _build_mass_figure(ion) -> plt.Figure:
    masses = mass_spectrum(ion, bin_width_amu=1.0)
    fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
    ax.bar(masses.bin_centers_amu, masses.counts, width=0.9, edgecolor="black", linewidth=0.5)
    ax.set_title("ion mass histogram")
    ax.set_xlabel("m / u")
    ax.set_ylabel("count")
    return fig


if __name__ == "__main__":
    raise SystemExit(main())
