"""Interactive comparison: I+He_n per-fragment speed distributions (n = 0-3).

Shows BOTH the 2-D and 3-D per-fragment overlay figures -- the interactive
twin of ``plot_run_summary.py``'s ``_section_ihe_ked_curves`` section,
called once per representation (imported from there, not reimplemented).
Intended for interactive IDE use (e.g. PyCharm): edit the
``USER SETTINGS`` block below, then either run this file directly::

    python scripts/post_processing/plot_speed_distribution_comparison.py

or run it from the IDE. Both figures are shown together with a single
``plt.show()`` call.
"""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import matplotlib.pyplot as plt  # noqa: E402

from i2_helium_md.postprocess import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

import plot_run_summary as run_summary  # noqa: E402


# =============================================================================
# USER SETTINGS -- edit these and run the script (e.g. from PyCharm)
# =============================================================================
# Path to the run directory holding cfg.json + neutral.npz + ion.npz.
RUN_DIR: Path = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet"

# Directory holding the frozen I+He_n kinetic-energy reference
# (IHe_KED_reference.csv + IHe_KED_curves_n{0..4}.csv). Required -- this
# script exists specifically to compare against the reference curves, so
# unlike plot_run_summary.py's optional section gating, None is refused
# with a clear error rather than silently skipping the figures.
IHE_KED_REFERENCE_DIR: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "ihe_ked"
)


def main() -> int:
    """Load the configured run's ion checkpoint and the ihe_ked reference,
    build both the 3-D and 2-D per-fragment speed-distribution overlay
    figures via the shared section builder, and show them together.

    Returns
    -------
    int
        0 on success (mirrors ``plot_run_summary.main``'s convention).

    Raises
    ------
    ValueError
        If ``IHE_KED_REFERENCE_DIR`` is ``None`` (this script exists to
        compare against the reference) or if the run directory has no ion
        checkpoint.
    """
    if IHE_KED_REFERENCE_DIR is None:
        raise ValueError(
            "IHE_KED_REFERENCE_DIR is None -- this script exists to "
            "compare against the I+He_n kinetic-energy reference curves; "
            "set it to a valid directory holding IHe_KED_reference.csv "
            "and IHe_KED_curves_n{0..4}.csv."
        )
    ked_dir = Path(IHE_KED_REFERENCE_DIR)

    run = RunDirectory(RUN_DIR)
    if not run.has_ion():
        raise ValueError(
            f"{RUN_DIR} has no ion checkpoint -- this is an ion-stage "
            "comparison (I+He_n speed distributions)."
        )
    ion = run.load_ion()

    ked_ref = load_ihe_ked_reference(ked_dir / "IHe_KED_reference.csv")

    run_summary._section_ihe_ked_curves(ion, ked_dir, ked_ref, "3d")
    run_summary._section_ihe_ked_curves(ion, ked_dir, ked_ref, "2d")

    plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
