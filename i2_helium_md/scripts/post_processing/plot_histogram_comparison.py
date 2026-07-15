"""Interactive comparison: final I+He_n mass spectrum vs experimental abundance.

Shows ONLY the mass-spectrum comparison figure -- the interactive
single-figure twin of ``plot_run_summary.py``'s ``_section_mass_spectrum``
section (imported from there, not reimplemented). Intended for interactive
IDE use (e.g. PyCharm): edit the ``USER SETTINGS`` block below, then either
run this file directly::

    python scripts/post_processing/plot_histogram_comparison.py

or run it from the IDE.
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

from i2_helium_md.postprocess import load_he_abundance_reference  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

import plot_run_summary as run_summary  # noqa: E402


# =============================================================================
# USER SETTINGS -- edit these and run the script (e.g. from PyCharm)
# =============================================================================
# Path to the run directory holding cfg.json + neutral.npz + ion.npz.
RUN_DIR: Path = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet_long"

# Experimental I+He_n abundance CSV for the mass-spectrum side-by-side.
# None keeps the plain simulated mass spectrum -- this is the section's own
# documented fallback (see plot_run_summary._section_mass_spectrum), so
# unlike the other two scripts in this trio, None is allowed here.
ABUNDANCE_REF_PATH: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"
)


def main() -> int:
    """Load the configured run's ion checkpoint, build the mass-spectrum
    comparison figure via the shared section builder, and show it.

    Returns
    -------
    int
        0 on success (mirrors ``plot_run_summary.main``'s convention).

    Raises
    ------
    ValueError
        If the run directory has no ion checkpoint, or if the mass
        spectrum has zero outside ions in every I+He_n gate (the section
        builder's own ``_SectionSkipped`` re-raised here as a clear
        error -- this script exists to show a figure, not to skip one).
    """
    run = RunDirectory(RUN_DIR)
    if not run.has_ion():
        raise ValueError(
            f"{RUN_DIR} has no ion checkpoint -- this is an ion-stage "
            "comparison (final I+He_n mass spectrum)."
        )
    ion = run.load_ion()

    abundance = (
        load_he_abundance_reference(ABUNDANCE_REF_PATH)
        if ABUNDANCE_REF_PATH else None
    )

    try:
        run_summary._section_mass_spectrum(ion, abundance)
    except run_summary._SectionSkipped as skip:
        raise ValueError(f"mass spectrum comparison unavailable: {skip}") from skip

    plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
