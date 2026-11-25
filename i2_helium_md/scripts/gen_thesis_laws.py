"""Run the thesis laws B2 ("LQ" in the thesis), B3 and h405: members, pool, score.

Set the options below, then: ``python scripts/gen_thesis_laws.py``.

A finished run or pool with the same model is skipped; a partial one or one with
another model is refused. To vary a parameter, set ``OVERRIDES`` and a ``TAG``;
the tag keeps the variant's directories apart from the thesis runs.
"""

from __future__ import annotations

import os
import sys
from dataclasses import replace
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_var, "1")     # one thread per member process

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.postprocess.battery_score import (  # noqa: E402
    format_table,
    score_pool,
    write_rows_csv,
)
from i2_helium_md.presets import THESIS_MEMBER_SEEDS, THESIS_PRESETS  # noqa: E402
from i2_helium_md.simulation.battery import (  # noqa: E402
    POOL_ARTIFACTS,
    is_finished,
    pool_battery,
    run_battery,
)
from i2_helium_md.simulation import thesis_runs  # noqa: E402
from i2_helium_md.simulation.thesis_runs import RUNS_ROOT  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================
LAWS = ("B2", "B3", "h405")
MEMBERS = "all"                              # or a list such as ["s1"]
RUN_SIZE = "production"                      # "smoke": N = 20, 1 ps, no pool, no score
CONCURRENCY = 2                              # members at a time
OUTPUT_ROOT = RUNS_ROOT
TAG = ""                                     # e.g. "_rrk10"; required with OVERRIDES
OVERRIDES: dict[str, object] = {}            # e.g. {"evap_rrk_dof": 10.0}
SCORES_CSV: Path | None = None               # None: OUTPUT_ROOT/thesis_scores{TAG}.csv
# =============================================================================

SMOKE = {"num_molecules": 20, "ion_simulation_time": 1.0}    # ion time in ps
SUMMARY_COLUMNS = ("label", "num_scored", "nbar_det", "n1_solv", "w1_solv", "midHot",
                   "deepKE", "KE1_mean", "KE2_mean", "chi2_med", "S", "trap")


def member_dir(law: str, member: str) -> Path:
    """Run directory of one member."""
    if RUN_SIZE == "smoke":
        return OUTPUT_ROOT / "smoke" / f"{law}{TAG}_N{SMOKE['num_molecules']}_{member}"
    return thesis_runs.member_dir(law, member, root=OUTPUT_ROOT, tag=TAG)


def pool_dir(law: str) -> Path:
    """Directory of the six-member pool."""
    return thesis_runs.pool_dir(law, root=OUTPUT_ROOT, tag=TAG)


def scores_csv() -> Path:
    """Path of the score table."""
    return SCORES_CSV or OUTPUT_ROOT / f"thesis_scores{TAG}.csv"


def _build(law: str, member: str) -> SimConfig:
    extra = {**OVERRIDES, **SMOKE} if RUN_SIZE == "smoke" else OVERRIDES
    return THESIS_PRESETS[law](member, **extra)


def _pool_cfg(law: str) -> SimConfig:
    cfg = _build(law, next(iter(THESIS_MEMBER_SEEDS)))
    return replace(cfg, num_molecules=cfg.num_molecules * len(THESIS_MEMBER_SEEDS))


def _check_settings() -> list[str]:
    """Refuse impossible settings; return the selected members."""
    if RUN_SIZE not in ("production", "smoke"):
        raise SystemExit(f"RUN_SIZE must be 'production' or 'smoke', got {RUN_SIZE!r}")
    if OVERRIDES and not TAG:
        raise SystemExit("OVERRIDES needs a TAG, so the variant gets its own directories")
    fixed = sorted({"seed", "num_molecules"} & set(OVERRIDES))
    if fixed:
        raise SystemExit(f"OVERRIDES cannot set {fixed}: the directory name states them")
    for law in LAWS:
        if law not in THESIS_PRESETS:
            raise KeyError(f"unknown law {law!r}; expected one of {list(THESIS_PRESETS)}")
    members = list(THESIS_MEMBER_SEEDS) if MEMBERS == "all" else list(MEMBERS)
    for member in members:
        if member not in THESIS_MEMBER_SEEDS:
            raise KeyError(
                f"unknown member {member!r}; expected one of {list(THESIS_MEMBER_SEEDS)}")
    return members


def main() -> int:
    """Run, pool and score the selected laws."""
    members = _check_settings()
    smoke = RUN_SIZE == "smoke"
    run_battery({member_dir(law, m): _build(law, m) for law in LAWS for m in members},
                concurrency=CONCURRENCY, smoke=smoke)
    if smoke:
        return 0
    if set(members) != set(THESIS_MEMBER_SEEDS):
        print("not pooled: a pool needs all members", flush=True)
        return 0
    for law in LAWS:
        pool_battery([member_dir(law, m) for m in THESIS_MEMBER_SEEDS], pool_dir(law))
    # Every law pooled so far, so that a call for one law keeps the others' rows.
    rows = [score_pool(pool_dir(law), law + TAG) for law in THESIS_PRESETS
            if is_finished(pool_dir(law), _pool_cfg(law), POOL_ARTIFACTS)]
    print(format_table([{c: row[c] for c in SUMMARY_COLUMNS} for row in rows]))
    print(f"scores -> {write_rows_csv(scores_csv(), rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
