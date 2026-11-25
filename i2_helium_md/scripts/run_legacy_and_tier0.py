"""Run the legacy hard-sphere runs and the Tier-0 runs behind the thesis figures.

Set the options below, then: ``python scripts/run_legacy_and_tier0.py``.

A finished run with the same model is skipped; a partial one or one with another
model is refused. The legacy runs are seeded, under new names: the runs behind the
submitted legacy figures were unseeded, so these agree with them statistically only.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.presets import LEGACY_THESIS_RUNS  # noqa: E402
from i2_helium_md.simulation.battery import claimed, is_finished  # noqa: E402
from scripts.extraction.method_b_g412_ridge import (  # noqa: E402
    BOUND_VARIANT,
    G11_M_EFF_AMU,
    point_cfg,
)
from scripts.gen_tier0_runs import run_tier0  # noqa: E402
from scripts.post_processing import tier0_calibration_figures as calibration  # noqa: E402
from scripts.post_processing.tier0_thesis_figure import LQ_LAWS, h405_9a_cfg  # noqa: E402
from scripts.tier0_common import run_dir_name  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================
RUN_LEGACY = True                            # 9 A, 18 A and droplet hard-sphere runs
RUN_TIER0 = True                             # Tier-0 drag runs of the calibration figures
RUN_SIZE = "production"                      # "smoke": N = 4, 0.5 ps, in OUTPUT_ROOT/smoke
OUTPUT_ROOT = PROJECT_ROOT / "data" / "runs"
# =============================================================================

ARTIFACTS = ("cfg.json", "neutral.npz", "ion.npz")
SMOKE = {"num_molecules": 4, "ion_simulation_time": 0.5}     # ion time in ps


def tier0_runs() -> dict[str, SimConfig]:
    """Tier-0 runs the thesis figures read: directory name -> config."""
    runs = {
        run_dir_name(run.case, "LQ", calibration.N, run.run_tag): calibration.build_cfg(run)
        for run in calibration.RUNS
    }
    cfg = h405_9a_cfg()
    runs[run_dir_name("9A", "CC", cfg.num_molecules, "g11fig_h405pp35m21")] = cfg
    for label, (a, c, e_bind, _) in LQ_LAWS.items():
        cfg = point_cfg(c, a, e_bind, mass_amu=G11_M_EFF_AMU)
        runs[run_dir_name("9A", BOUND_VARIANT, cfg.num_molecules, f"g11m21{label}")] = cfg
    return runs


def run_one(cfg: SimConfig, run_dir: Path) -> None:
    """Neutral and ion stage of ``cfg`` into ``run_dir``, unless finished there."""
    if is_finished(run_dir, cfg, ARTIFACTS):
        print(f"[{run_dir.name}] skip (complete)", flush=True)
        return
    cfg.validate()
    with claimed(run_dir, cfg, ARTIFACTS) as start:
        if start:
            print(f"[{run_dir.name}] running, seed {cfg.seed}", flush=True)
            run_tier0(cfg, run_dir)
            print(f"[{run_dir.name}] done", flush=True)


def main() -> int:
    """Run the selected sets."""
    if RUN_SIZE not in ("production", "smoke"):
        raise SystemExit(f"RUN_SIZE must be 'production' or 'smoke', got {RUN_SIZE!r}")
    runs: dict[str, SimConfig] = {}
    if RUN_LEGACY:
        runs |= {name: build() for name, build in LEGACY_THESIS_RUNS.items()}
    if RUN_TIER0:
        runs |= tier0_runs()
    root = OUTPUT_ROOT
    if RUN_SIZE == "smoke":
        root = OUTPUT_ROOT / "smoke"
        runs = {name: replace(cfg, **SMOKE) for name, cfg in runs.items()}
    for name, cfg in runs.items():
        run_one(cfg, root / name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
