"""h405 contrast law at both production settings: ``post_pickup`` + 3.5 A gate.

Question (user, 2026-10-14): the thesis compares B2 (``post_pickup``
sampler, 3.5 A density gate) with the h405 landing, which was run with the
``raw`` sampler and the 14.2 A gate. h405 was only ever run with one of the
two changed (``h405g35``: raw + 3.5 A; h405 ``post_pickup``: 14.2 A, not
re-landed, nbar 4.59). This generator runs h405 with **both** production
settings and E0 re-placed to restore the landing, so B2 vs h405 differs in
the drag law only.

Cell: the h405 production member (``gen_tier2_detfix_battery.build_member``:
capped_cubic v_c 5.5, tau 4.4 ps, E0 0.405 eV) with

* the 3.5 A density gate exactly as ``h405g35`` (G4.13 definition:
  ``erf_independent`` s_rho 3.5 A + ``allow_inconsistent_gate_pairing`` --
  the h405 bundle was extracted at 14.2 A, so this is a production setting,
  not a Tier-0 refit),
* ``droplet_size_sampler_mode="post_pickup"``,
* E0 re-placed (the only re-tuned knob, as for B2 under G5).

E0 placement. Projected nbar shift of the two flips at h405 E0: sampler
+0.717 (``detfix_postpickup_table.csv`` pooled_pp5000 - pooled_raw5000) and
gate -0.033 (``g414_table.csv`` h405g35_pooled6000 - h405_pooled6000);
restoring B2's pooled nbar 3.903 needs ~ -0.67, at dnbar/dE0 ~ -25 /eV
(h405 A/B) ~ +0.027 eV. Step ``--column``: E0 in :data:`COLUMN_E0_EV` at s1.
Step ``--pool``: :data:`POOL_E0_EV` (set from the column, interpolated in
nbar to :data:`NBAR_TARGET`) at s1-s5 + r6.

Oracle 1 (in memory; the h405 member run dirs are no longer on disk): the
cfg differs from ``build_member(member)`` in exactly the h405g35 gate keys
+ the sampler + the E0 key.

Usage::

    python scripts/gen_tier2_h405_pp35.py --column [--dry-run]
    python scripts/gen_tier2_h405_pp35.py --pool [--dry-run]

Scorer: ``scripts/post_processing/tier2_h405_pp35_table.py``; pool target
``build_pooled_detection_container.py --target h405_pp35``. Max two
concurrent MD runs; ion.npz dropped after detection verifies.
"""

from __future__ import annotations

import argparse
import dataclasses
import math
import os
from pathlib import Path
import sys
from typing import Optional

for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2_detfix_battery import (  # noqa: E402
    MAX_BYTES_ION, MEMBER_SEEDS, T_HANDOVER_PS, build_member, verify_handover_alignment,
)
from scripts.gen_tier2_detfix_postpickup import LQ_DROP_REASON  # noqa: E402
from scripts.gen_tier2_g413_ridge_battery import (  # noqa: E402
    H405_GATE_LABEL, _E0_KEY, _H405G35_KEYS, build as g413_build,
)
from scripts.gen_tier2_g414_tauline import drop_ion_npz  # noqa: E402
from scripts.gen_tier2atlas_g4finals import BUDGET_EV, CASE, N  # noqa: E402
from scripts.tier0_common import run_dir_name  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.gen_tier2atlas_linclone import cfg_field_diff  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================

SAMPLER_MODE = "post_pickup"
#: E0 column [eV] at s1 (h405 own E0 0.405; projected landing ~0.432).
COLUMN_E0_EV = (0.425, 0.435, 0.445)
COLUMN_MEMBER = "s1"
#: nbar the pool is placed at: B2's pooled nbar (G5 Outcome).
NBAR_TARGET = 3.903
#: Pool E0 [eV], set from the column (s1, 2026-10-14): nbar 4.034 at 0.425,
#: 3.819 at 0.435 -> linear in nbar to NBAR_TARGET 3.903 -> 0.4311 -> 0.431.
POOL_E0_EV: Optional[float] = 0.431
MAX_CONCURRENCY = 2

# =============================================================================


def run_dir_for(member: str, e0_eV: float) -> Path:
    """e.g. ``…shared_pure_cubic_N1000_detfix_conf270_h405pp35e4350s1_th500``."""
    tag = f"h405pp35e{round(e0_eV * 1e4):04d}{member}"
    return PROJECT_ROOT / "data" / "runs" / run_dir_name(
        CASE, "shared_pure_cubic", N, run_tag=f"detfix_conf270_{tag}_th500")


def build(member: str, e0_eV: float) -> SimConfig:
    """h405 member with the h405g35 gate, ``post_pickup`` and E0 [eV]."""
    cfg = dataclasses.replace(g413_build(H405_GATE_LABEL, member),
                              droplet_size_sampler_mode=SAMPLER_MODE,
                              **{_E0_KEY: float(e0_eV) / BUDGET_EV})
    cfg.validate()
    return cfg


def verify(member: str, e0_eV: float) -> list[str]:
    """Oracle 1: in-memory diff vs the h405 member is exactly the h405g35
    gate keys + the sampler + the E0 key. Returns the keys."""
    diff = sorted(cfg_field_diff(build_member(member), build(member, e0_eV)))
    want = set(_H405G35_KEYS) | {"droplet_size_sampler_mode", _E0_KEY}
    if set(diff) != want:
        raise AssertionError(f"[{member} E0 {e0_eV}] diff vs h405 {diff}; "
                             f"pre-registered {sorted(want)}")
    return diff


def run_one(job: tuple[str, float]) -> str:
    """Guard, neutral -> ion (to t_h) -> detection, drop ion.npz."""
    member, e0 = job
    verify(member, e0)
    cfg = build(member, e0)
    steps, stride = verify_handover_alignment(cfg)
    run_dir = run_dir_for(member, e0)
    tag = f"{member} E0 {e0}"
    if (run_dir / "detection.npz").exists() and (run_dir / "cfg.json").exists():
        drop_ion_npz(run_dir, reason=LQ_DROP_REASON)
        print(f"[{tag}] skip (complete) -> {run_dir.name}", flush=True)
        return tag
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{tag}] seed {cfg.seed}, {steps} steps (stride {stride}) -> {run_dir.name}",
          flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False,
                              max_bytes=MAX_BYTES_ION)
    t_last = float(np.asarray(ion.time_ps)[-1])
    if not math.isclose(t_last, T_HANDOVER_PS, abs_tol=1e-9):
        raise AssertionError(f"[{tag}] ion.npz last column at {t_last} ps, not t_h")
    det = run_detection_stage(ion, cfg, save_path=run.root / "detection.npz")
    del ion, neutral
    drop_ion_npz(run_dir, reason=LQ_DROP_REASON)
    print(f"[{tag}] done: scored {int(det.detected_mask.sum())}/{det.detected_mask.size}",
          flush=True)
    return tag


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--column", action="store_true", help="E0 column at s1")
    group.add_argument("--pool", action="store_true", help="POOL_E0_EV at s1-s5 + r6")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    if args.pool:
        if POOL_E0_EV is None:
            raise SystemExit("POOL_E0_EV not set: run --column and set it first")
        jobs = [(m, POOL_E0_EV) for m in MEMBER_SEEDS]
    else:
        jobs = [(COLUMN_MEMBER, e0) for e0 in COLUMN_E0_EV]
    for m, e0 in jobs:
        diff = verify(m, e0)
        steps, stride = verify_handover_alignment(build(m, e0))
        print(f"[{m} E0 {e0}] guards PASSED: diff {diff}, {steps} steps, stride {stride} "
              f"-> {run_dir_for(m, e0).name}", flush=True)
    if args.dry_run:
        return 0
    from concurrent.futures import ProcessPoolExecutor, as_completed
    failures = 0
    with ProcessPoolExecutor(max_workers=min(MAX_CONCURRENCY, len(jobs))) as pool:
        futures = {pool.submit(run_one, j): j for j in jobs}
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[{futures[fut]}] FAILED: {exc!r}", flush=True)
    print(f"{len(jobs) - failures}/{len(jobs)} runs complete.", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
