"""Add the Stage-I event log to the 18 thesis run dirs (THESIS_DISCUSSION §16).

The thesis members ``data/runs/{B2,B3,h405}_N1000_s1..s6`` predate
``stage1_events.npz``. Recording it is RNG-free, so a re-run must reproduce
every existing artifact exactly. Per member, sequentially:

1. Oracles: the generator's cfg guards, and the production cfg.json must equal
   the rebuilt thesis cfg (refuse otherwise).
2. Back up the production artifacts to ``data/runs/_stage1_events_backup/<dir>``.
3. Re-run neutral -> ion (with the log) -> detection into
   ``data/runs/_stage1_events_regen/<dir>`` (production files are not written).
4. Bit-identity: cfg.json text and every array of neutral / ion / detection
   (key set, dtype, shape, raw bytes) against the backup.
5. On pass: copy ``stage1_events.npz`` into the production dir, delete the
   regen dir and the backup. On fail: keep both, stop with the differing keys.

Usage::

    python scripts/regen_tier2_stage1_events.py                 # all 18
    python scripts/regen_tier2_stage1_events.py --laws B2 --members s1
"""

from __future__ import annotations

import argparse
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import shutil
import sys

for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2_detfix_battery import MEMBER_SEEDS  # noqa: E402
from scripts.gen_tier2_detfix_postpickup import (  # noqa: E402
    THESIS_LQ_JOBS, build_thesis_lq, run_thesis_member, thesis_run_dir, verify_lq_member,
)
from scripts import gen_tier2_h405_pp35 as h405  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from i2_helium_md.simulation.stage1_events import (  # noqa: E402
    STAGE1_EVENTS_FILENAME, load_stage1_events, replay_n_shell,
)

# =============================================================================
# USER SETTINGS
# =============================================================================
LAWS = ("B2", "B3", "h405")
BACKUP_ROOT = PROJECT_ROOT / "data" / "runs" / "_stage1_events_backup"
REGEN_ROOT = PROJECT_ROOT / "data" / "runs" / "_stage1_events_regen"
ARTIFACTS = ("cfg.json", "neutral.npz", "ion.npz", "detection.npz")
MAX_CONCURRENCY = 2
# =============================================================================


def thesis_cfg(law: str, member: str) -> SimConfig:
    """The thesis cfg of ``law`` / ``member`` after the generator's oracles."""
    if law == "h405":
        v_c, tau, e0 = h405.SCAN_POOL_CELL
        h405.verify(member, e0, v_c, tau)
        return h405.build_thesis(member)
    verify_lq_member(member, *THESIS_LQ_JOBS[law])
    return build_thesis_lq(law, member)


def npz_differences(a: Path, b: Path) -> list[str]:
    """Keys whose dtype, shape or raw bytes differ between two npz files
    (plus keys present in only one). Empty list = bit-identical content."""
    with np.load(a, allow_pickle=False) as za, np.load(b, allow_pickle=False) as zb:
        diff = sorted(set(za.files) ^ set(zb.files))
        for k in sorted(set(za.files) & set(zb.files)):
            x, y = za[k], zb[k]
            if x.dtype != y.dtype or x.shape != y.shape or x.tobytes() != y.tobytes():
                diff.append(k)
    return diff


def regen_member(job: tuple[str, str]) -> str:
    """Steps 1-5 of the module docstring for one ``(law, member)``."""
    law, member = job
    label = f"{law} {member}"
    prod = thesis_run_dir(law, member)
    if (prod / STAGE1_EVENTS_FILENAME).exists():
        return f"[{label}] skip: {STAGE1_EVENTS_FILENAME} already present"
    missing = [f for f in ARTIFACTS if not (prod / f).exists()]
    if missing:
        raise FileNotFoundError(f"[{label}] {prod.name} lacks {missing}")
    cfg = thesis_cfg(law, member)
    if RunDirectory(prod).load_cfg() != cfg:
        raise RuntimeError(f"[{label}] {prod.name}/cfg.json differs from the thesis cfg")

    backup, regen = BACKUP_ROOT / prod.name, REGEN_ROOT / prod.name
    for d in (backup, regen):
        if d.exists():
            raise FileExistsError(f"[{label}] {d} exists from an earlier attempt; inspect it")
    backup.mkdir(parents=True)
    for f in ARTIFACTS:
        shutil.copy2(prod / f, backup / f)

    run_thesis_member(cfg, regen, label)

    diffs = {}
    if (regen / "cfg.json").read_bytes() != (backup / "cfg.json").read_bytes():
        diffs["cfg.json"] = ["<text>"]
    for f in ARTIFACTS[1:]:
        d = npz_differences(backup / f, regen / f)
        if d:
            diffs[f] = d
    if diffs:
        raise AssertionError(f"[{label}] NOT bit-identical, kept {backup} and {regen}: {diffs}")

    log = load_stage1_events(regen / STAGE1_EVENTS_FILENAME)
    ion = RunDirectory(regen).load_ion()
    steps = np.rint(np.asarray(ion.time_ps) / cfg.dt_ion).astype(int)
    if not np.array_equal(replay_n_shell(log, ion.n_shell[:, 0], steps),
                          ion.n_shell.astype(int)):
        raise AssertionError(f"[{label}] event log does not replay n_shell; kept {regen}")

    shutil.copy2(regen / STAGE1_EVENTS_FILENAME, prod / STAGE1_EVENTS_FILENAME)
    shutil.rmtree(regen)
    shutil.rmtree(backup)
    return (f"[{label}] bit-identical; {log.event_ion.size} events, "
            f"{log.crossing_ion.size} crossings -> {prod.name}/{STAGE1_EVENTS_FILENAME}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--laws", nargs="*", default=list(LAWS), choices=LAWS)
    parser.add_argument("--members", nargs="*", default=None,
                        help="thesis labels s1..s6 (default all)")
    parser.add_argument("--concurrency", type=int, default=MAX_CONCURRENCY)
    args = parser.parse_args(argv)
    by_label = {("s6" if m == "r6" else m): m for m in MEMBER_SEEDS}
    labels = args.members or list(by_label)
    unknown = sorted(set(labels) - set(by_label))
    if unknown:
        raise SystemExit(f"unknown member(s) {unknown}; expected {list(by_label)}")
    jobs = [(law, by_label[s]) for law in args.laws for s in labels]
    if args.concurrency <= 1:
        for job in jobs:
            print(regen_member(job), flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.concurrency) as pool:
            for msg in pool.map(regen_member, jobs):
                print(msg, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
