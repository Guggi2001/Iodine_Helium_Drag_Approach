"""Paper-v2 port of post_process_single_pulse_paper_IplusHe_comparison.m.

The Python name ``paper_v2`` is a consistency alias. The legacy source is the
I+He comparison MATLAB script. The output is split into one figure per
diagnostic so each can be inspected independently:

- ``paper_v2_vmi_comparison.png``        experimental + simulated 2-D VMI
- ``paper_v2_radial_comparison.png``     radial velocity distribution
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

import matplotlib.pyplot as plt  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    load_paper_v2_radial_references,
    paper_v2_phi_curve,
    paper_v2_velocity_curve,
    paper_v2_velocity_map,
)
from i2_helium_md.postprocess.paper_v2_plotting import (  # noqa: E402
    build_phi_figure,
    build_polar_image_figure,
    build_radial_figure,
    build_vmi_figure,
    load_optional_image,
    load_optional_phi,
    load_optional_polar_image,
    polar_histogram_matched_to_reference,
)
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


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    run = RunDirectory(args.run_dir)
    ion = run.load_ion()

    radial_refs = load_paper_v2_radial_references(args.reference_dir)
    if not radial_refs:
        print(f"[paper_v2] no optional radial references found in {args.reference_dir}")

    image_ref = load_optional_image(args.reference_dir)
    polar_ref = load_optional_polar_image(args.reference_dir)
    phi_ref = load_optional_phi(args.reference_dir)
    velocity_map = paper_v2_velocity_map(ion, mass_amu=MASS_SELECTION_AMU)
    velocity_curve = paper_v2_velocity_curve(ion, mass_amu=MASS_SELECTION_AMU)
    phi_curve = paper_v2_phi_curve(ion, mass_amu=MASS_SELECTION_AMU)

    fig_vmi = build_vmi_figure(
        image_ref=image_ref,
        velocity_map=velocity_map,
        experimental_noise_floor=args.noise_floor,
    )
    fig_radial = build_radial_figure(radial_refs, velocity_curve)
    fig_phi = build_phi_figure(phi_curve, phi_ref=phi_ref)
    fig_polar = None
    if polar_ref is not None:
        polar_hist = polar_histogram_matched_to_reference(
            ion, polar_ref, mass_amu=MASS_SELECTION_AMU,
        )
        fig_polar = build_polar_image_figure(
            polar_ref, polar_hist, experimental_noise_floor=args.noise_floor
        )

    out_dir = run.root / "figures"
    out_dir.mkdir(exist_ok=True)
    fig_vmi.savefig(out_dir / "paper_v2_vmi_comparison.png", dpi=150)
    fig_radial.savefig(out_dir / "paper_v2_radial_comparison.png", dpi=150)
    fig_phi.savefig(out_dir / "paper_v2_phi_comparison.png", dpi=150)
    print(f"Saved VMI comparison to {out_dir / 'paper_v2_vmi_comparison.png'}")
    print(f"Saved radial comparison to {out_dir / 'paper_v2_radial_comparison.png'}")
    print(f"Saved phi comparison to {out_dir / 'paper_v2_phi_comparison.png'}")
    if fig_polar is not None:
        fig_polar.savefig(out_dir / "paper_v2_polar_image_comparison.png", dpi=150)
        print(f"Saved polar image comparison to {out_dir / 'paper_v2_polar_image_comparison.png'}")

    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
