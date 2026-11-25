"""The five Tier-0 microscopic-calibration figures of the thesis.

Three LQ (B2) runs at m = I+ + 21 He, each drawn with the builders of
``tier0_drag_comparison.py``:

* ``18A_Stuck_Velocity_R.png``, ``18A_Stuck_Force_Energy.png`` -- 18 A at the
  TDDFT binding 0.308 eV: the ion is recaptured on the longer 50 ps timescale.
* ``18A_B2.png``, ``18A_B2_Force_Energy.png`` -- 18 A at the B2 binding.
* ``9A_B2.png`` -- 9 A at the B2 binding; panels (c), (d), since it sits next
  to ``18A_B2.png`` (a), (b) in the thesis.

The runs (``scripts/tier0_runs.py``) are made by
``scripts/gen_tier0_runs.py``; a missing or different run is refused. The PNGs are half the text width, in thesis style, and go
side by side in pairs at 0.5 \\linewidth, so their sizes must stay identical. The
exception is ``18A_B2_Force_Energy.png``, which stands alone at the full text width.

Usage::

    python scripts/post_processing/tier0_calibration_figures.py
"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import NamedTuple, Optional

# =============================================================================
# USER SETTINGS
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "data" / "figures" / "thesis_microscopic_calibration"
DPI = 200                 # as the other thesis figures
FORCE_ATOM_INDEX = 0      # force / energy figure: fragment 1


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tier0_runs import (  # noqa: E402
    B2_9A,
    B2_18A,
    STUCK_18A,
    CalibrationRun,
    finished_tier0_run,
)


class Run(NamedTuple):
    """The figures drawn from one Tier-0 run."""
    sim: CalibrationRun
    t_max_ps: float             # x-axis limit of the figures [ps]
    velocity_png: str           # distance + speed figure
    force_png: Optional[str]    # force + energy figure, or None
    force_alone: bool = False   # force figure not in a pair: full width, side legend
    velocity_letters: tuple[str, str] = ("a", "b")   # panel letters of the velocity figure
    mark_surface: bool = False  # force figure: dashed line at the droplet surface crossing


RUNS = (
    Run(STUCK_18A, 50.0, "18A_Stuck_Velocity_R.png", "18A_Stuck_Force_Energy.png"),
    Run(B2_18A, 20.0, "18A_B2.png", "18A_B2_Force_Energy.png",
        force_alone=True, mark_surface=True),
    # 9A_B2.png sits next to 18A_B2.png (a, b) in the thesis, so it continues with (c, d).
    Run(B2_9A, 20.0, "9A_B2.png", None, velocity_letters=("c", "d")),
)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from i2_helium_md.postprocess import load_hedft_trajectory  # noqa: E402
from i2_helium_md.presets import tddft_reference_path  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.post_processing.thesis_style import FIG_RC  # noqa: E402
from scripts.post_processing.tier0_drag_comparison import (  # noqa: E402
    build_figure,
    build_force_figure,
)
from scripts.tier0_common import TIER0_WINDOWS_PS  # noqa: E402


def draw(run: Run, directory: RunDirectory, out_dir: Path) -> None:
    """Write the figures of ``run`` into ``out_dir``."""
    ion, cfg = directory.load_ion(), directory.load_cfg()
    hedft = load_hedft_trajectory(tddft_reference_path(run.sim.case))
    window = TIER0_WINDOWS_PS[run.sim.case]
    with plt.rc_context(FIG_RC):              # thesis_style font (Palatino)
        figures = {run.velocity_png: build_figure(
            ion, hedft, window, False, float(ion.droplet_radii_angstrom[0]),
            case=run.sim.case, t_max=run.t_max_ps, letters=run.velocity_letters)}
        if run.force_png is not None:
            figures[run.force_png] = build_force_figure(
                ion, cfg, window, hedft, atom_index=FORCE_ATOM_INDEX, t_max=run.t_max_ps,
                alone=run.force_alone, mark_surface=run.mark_surface, letters=("a", "b") if run.force_alone else ("c", "d"))
        for name, fig in figures.items():
            fig.savefig(out_dir / name, dpi=DPI)
            plt.close(fig)
            print(f"[tier0_figures] wrote {out_dir / name}", flush=True)


def main(out_dir: Path = OUT_DIR) -> int:
    """Write all five figures to ``out_dir``."""
    out_dir.mkdir(parents=True, exist_ok=True)
    directories = [finished_tier0_run(run.sim.name) for run in RUNS]
    for run, directory in zip(RUNS, directories):
        draw(run, directory, out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1]) if len(sys.argv) > 1 else OUT_DIR))
