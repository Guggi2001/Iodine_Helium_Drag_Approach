"""The five Tier-0 microscopic-calibration figures of the thesis.

Three LQ (B2) runs at m = I+ + 21 He, each drawn with the builders of
``tier0_drag_comparison.py``:

* ``18A_Stuck_Velocity_R.png``, ``18A_Stuck_Force_Energy.png`` -- 18 A at the
  TDDFT binding 0.308 eV: the ion is recaptured on the longer 50 ps timescale.
* ``18A_B2.png``, ``18A_B2_Force_Energy.png`` -- 18 A at the B2 binding.
* ``9A_B2.png`` -- 9 A at the B2 binding.

A run is reused when its ``cfg.json`` matches and regenerated otherwise (N = 50,
about a minute each). The PNGs are half the text width, in thesis style, and go
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
N = 50
SEED = 20260604
DT_ION_PS = 0.01
M_EFF_AMU = 126.90 + 21 * 4.0026   # I+ + 21 He (G11.2) [amu]
FORCE_ATOM_INDEX = 0      # force / energy figure: fragment 1


class Run(NamedTuple):
    """One Tier-0 run and the figures drawn from it."""
    case: str                   # droplet geometry
    run_tag: str                # run-dir marker
    ion_time_ps: float          # ion-stage duration [ps]
    t_max_ps: float             # x-axis limit of the figures [ps]
    velocity_png: str           # distance + speed figure
    force_png: Optional[str]    # force + energy figure, or None
    e_bind_eV: Optional[float] = None   # binding [eV]; None = the LQ bundle's
    force_alone: bool = False   # force figure not in a pair: full width, side legend


RUNS = (
    Run("18A", "B2m21stuck", 50.0, 50.0, "18A_Stuck_Velocity_R.png",
        "18A_Stuck_Force_Energy.png", e_bind_eV=0.308),
    Run("18A", "B2m21", 30.0, 20.0, "18A_B2.png", "18A_B2_Force_Energy.png",
        force_alone=True),
    Run("9A", "B2m21", 30.0, 20.0, "9A_B2.png", None),
)

# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.postprocess import load_hedft_trajectory  # noqa: E402
from i2_helium_md.presets import tddft_reference_path  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.gen_tier0_runs import run_tier0  # noqa: E402
from scripts.post_processing.thesis_style import FIG_RC  # noqa: E402
from scripts.post_processing.tier0_drag_comparison import (  # noqa: E402
    build_figure,
    build_force_figure,
)
from scripts.tier0_common import TIER0_WINDOWS_PS, build_tier0_cfg, run_dir_name  # noqa: E402
from scripts.tier2_common import with_current_bundle_ids  # noqa: E402


def build_cfg(run: Run) -> SimConfig:
    """The Tier-0 config of ``run``."""
    return build_tier0_cfg(
        run.case, "LQ", num_molecules=N, ion_time_ps=run.ion_time_ps,
        dt_ion_ps=DT_ION_PS, seed=SEED, m_eff_amu=M_EFF_AMU,
        e_bind_override=run.e_bind_eV,
    )


def ensure_run(run: Run) -> RunDirectory:
    """The finished run directory of ``run``, generated if missing or different."""
    cfg = build_cfg(run)
    path = PROJECT_ROOT / "data" / "runs" / run_dir_name(run.case, "LQ", N, run.run_tag)
    directory = RunDirectory(path)
    if not (directory.has_cfg() and (path / "ion.npz").is_file()
            and with_current_bundle_ids(directory.load_cfg()) == cfg):
        print(f"[tier0_figures] running {path.name} ...", flush=True)
        run_tier0(cfg, path)
    return directory


def draw(run: Run, directory: RunDirectory, out_dir: Path) -> None:
    """Write the figures of ``run`` into ``out_dir``."""
    ion, cfg = directory.load_ion(), directory.load_cfg()
    hedft = load_hedft_trajectory(tddft_reference_path(run.case))
    window = TIER0_WINDOWS_PS[run.case]
    with plt.rc_context(FIG_RC):              # thesis_style font (Palatino)
        figures = {run.velocity_png: build_figure(
            ion, hedft, window, False, float(ion.droplet_radii_angstrom[0]),
            case=run.case, t_max=run.t_max_ps)}
        if run.force_png is not None:
            figures[run.force_png] = build_force_figure(
                ion, cfg, window, hedft, atom_index=FORCE_ATOM_INDEX, t_max=run.t_max_ps,
                alone=run.force_alone, letters=("a", "b") if run.force_alone else ("c", "d"))
        for name, fig in figures.items():
            fig.savefig(out_dir / name, dpi=DPI)
            plt.close(fig)
            print(f"[tier0_figures] wrote {out_dir / name}", flush=True)


def main(out_dir: Path = OUT_DIR) -> int:
    """Generate any missing run, then write all five figures to ``out_dir``."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for run in RUNS:
        draw(run, ensure_run(run), out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1]) if len(sys.argv) > 1 else OUT_DIR))
