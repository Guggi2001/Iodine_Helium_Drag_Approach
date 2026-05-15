"""Literal droplet-branch port of post_process_single_pulse_paper_v3.m.

This script follows the active non-effusive I+He comparison branch of the
legacy MATLAB paper figure. Experimental curves are optional CSV exports from
the MATLAB recipe; simulation curves are computed from a finished
``RunDirectory`` using the same mass selection, projected velocity, binning,
smoothing, and max-only normalization as the MATLAB script.
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
V3_REFERENCE_DIR = PROJECT_ROOT / "data" / "reference" / "paper_v3"
V3_IHE_RADIAL_PATH = V3_REFERENCE_DIR / "iplus_he_high_snr_radial.csv"
V3_TIMESCAN_RADIAL_PATH = V3_REFERENCE_DIR / "timescan_296_297_radial.csv"
V3_IHE_PHI_PATH = V3_REFERENCE_DIR / "iplus_he_high_snr_phi.csv"

MASS_SELECTIONS_AMU = (127.0, 131.0, 135.0)
VMIN_ANGULAR_DISTR_MPS = 0.0


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    load_paper_v3_phi_reference,
    load_paper_v3_radial_reference,
    mass_spectrum,
    paper_v3_phi_curve,
    paper_v3_velocity_curve,
)
from i2_helium_md.postprocess.paper_v3 import matlab_max_normalise  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Reproduce the active droplet branch of post_process_single_pulse_paper_v3.m."
    )
    p.add_argument("--run-dir", type=Path, default=RUN_DIR,
                   help="RunDirectory containing cfg.json, neutral.npz, and ion.npz.")
    p.add_argument("--ihe-radial-ref", type=Path, default=V3_IHE_RADIAL_PATH,
                   help="Optional v3 I+He radial CSV export (v_mps + signal columns).")
    p.add_argument("--timescan-radial-ref", type=Path, default=V3_TIMESCAN_RADIAL_PATH,
                   help="Optional v3 timescan radial CSV export (v_mps + signal columns).")
    p.add_argument("--ihe-phi-ref", type=Path, default=V3_IHE_PHI_PATH,
                   help="Optional v3 I+He phi CSV export (phi_rad,signal_arb).")
    p.add_argument("--no-show", action="store_true",
                   help="Skip plt.show(), useful for tests and headless runs.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    run = RunDirectory(args.run_dir)
    ion = run.load_ion()

    ihe_radial = _load_optional_radial(args.ihe_radial_ref, "I+He radial")
    timescan_radial = _load_optional_radial(args.timescan_radial_ref, "timescan radial")
    ihe_phi = _load_optional_phi(args.ihe_phi_ref, "I+He phi")

    velocity_curves = []
    phi_curves = []
    for mass in MASS_SELECTIONS_AMU:
        try:
            velocity_curves.append(paper_v3_velocity_curve(ion, mass_amu=mass))
        except ValueError as exc:
            print(f"[paper_v3] skip velocity mass {mass:.0f}: {exc}")
        try:
            phi_curves.append(paper_v3_phi_curve(ion, mass_amu=mass))
        except ValueError as exc:
            print(f"[paper_v3] skip phi mass {mass:.0f}: {exc}")

    fig_main = _build_main_figure(
        ihe_radial=ihe_radial,
        timescan_radial=timescan_radial,
        ihe_phi=ihe_phi,
        velocity_curves=velocity_curves,
        phi_curves=phi_curves,
    )
    fig_mass = _build_mass_figure(ion)

    out_dir = run.root / "figures"
    out_dir.mkdir(exist_ok=True)
    fig_main.savefig(out_dir / "compare_simulation_and_measurement.png", dpi=150)
    fig_mass.savefig(out_dir / "ion_mass_histogram.png", dpi=150)
    print(f"Saved figure to {out_dir / 'compare_simulation_and_measurement.png'}")
    print(f"Saved mass histogram to {out_dir / 'ion_mass_histogram.png'}")

    if not args.no_show:
        plt.show()
    return 0


def _load_optional_radial(path: Path, label: str):
    if path is None or not path.exists():
        print(f"[paper_v3] optional {label} reference not found: {path}")
        return None
    return load_paper_v3_radial_reference(path)


def _load_optional_phi(path: Path, label: str):
    if path is None or not path.exists():
        print(f"[paper_v3] optional {label} reference not found: {path}")
        return None
    return load_paper_v3_phi_reference(path)


def _build_main_figure(
    *,
    ihe_radial,
    timescan_radial,
    ihe_phi,
    velocity_curves,
    phi_curves,
) -> plt.Figure:
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 8.0), constrained_layout=True)
    _draw_velocity_panel(axes[0], ihe_radial, timescan_radial, velocity_curves)
    _draw_phi_panel(axes[1], ihe_phi, phi_curves)
    return fig


def _draw_velocity_panel(ax, ihe_radial, timescan_radial, velocity_curves) -> None:
    if timescan_radial is not None:
        _plot_radial_reference(ax, timescan_radial, "I2:I+He TS (296:297)", linestyle="-", alpha=0.65)
    if ihe_radial is not None:
        _plot_radial_reference(ax, ihe_radial, "I2:I+He (high-SNR)", linestyle=":", alpha=1.0)

    styles = {
        127.0: ("--", "simulated v.distr. m=127"),
        131.0: (":", "simulated v.distr. m=131"),
        135.0: ("-.", "simulated v.distr. m=135"),
    }
    for curve in velocity_curves:
        style, label = styles.get(curve.mass_amu, ("--", f"simulated v.distr. m={curve.mass_amu:.0f}"))
        ax.plot(curve.bin_centers_mps, curve.normalised, linestyle=style, linewidth=1.5, label=label)

    ax.set_xlabel("v / m/s")
    ax.set_ylabel("signal / arb. units")
    ax.set_xlim(0.0, 3500.0)
    ax.set_ylim(0.0, 1.1)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title("(a) radial velocity distribution")


def _plot_radial_reference(ax, ref, label: str, *, linestyle: str, alpha: float) -> None:
    for idx, signal_label in enumerate(ref.signal_labels):
        trace_label = label if idx == 0 else f"{label} {signal_label}"
        ax.plot(
            ref.velocity_mps,
            matlab_max_normalise(ref.signal[:, idx]),
            linestyle=linestyle,
            linewidth=1.5,
            alpha=alpha,
            label=trace_label,
        )


def _draw_phi_panel(ax, ihe_phi, phi_curves) -> None:
    if ihe_phi is not None:
        ax.plot(
            ihe_phi.phi_rad,
            matlab_max_normalise(ihe_phi.signal_arb),
            linewidth=1.5,
            label="I2:I+He (high-SNR)",
        )

    styles = {
        127.0: ("--", "simulation m=127"),
        131.0: (":", "simulation m=131"),
        135.0: ("-.", "simulation m=135"),
    }
    for curve in phi_curves:
        style, label = styles.get(curve.mass_amu, ("--", f"simulation m={curve.mass_amu:.0f}"))
        ax.plot(curve.bin_centers_rad, curve.normalised, linestyle=style, linewidth=1.5, label=label)

    ax.set_xlabel("phi / radian")
    ax.set_ylabel("signal / arb. units")
    ax.set_xlim(0.0, 2.0 * 3.141592653589793)
    ax.set_ylim(0.0, 1.1)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title("(b) azimuthal distribution")


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
