"""Atlas G4 Step 2 Block V — the pooled verification battery at h405.

Runs the `TIER2_SENSITIVITY_ATLAS_PLAN` §3.5f Block V: **5 members ×
N = 1000, five fresh seeds**, all at the G4 Step 1 successor candidate
**h405** (v_c 5.5 / τ 4.4 / E₀ 0.405, corrected geometry). The §4cc battery
pattern: one cell, many seeds — verification (GV-P1..P3) plus 5× the
statistics for the Block-D W₁ anatomy (single-seed W₁ SD is 0.095).

**Zero pin duplication:** each member cfg is
`gen_tier2atlas_g4finals.build_cell(h405)` with only the seed swapped, and
an oracle asserts the member cfg diffs against the committed Block-3
`g4fh405` run's `cfg.json` in **exactly** `{"seed"}` — same cell, fresh
seed, no drifted pin.

Pre-registered (plan §3.5f):

* **GV-P1 (gate):** pooled n1_solv ∈ [0.19, 0.30] ∧ pooled nbar_det ∈
  [3.77, 4.37].
* **GV-P2 (floor consistency):** pooled w1_solv ∈ [0.64, 0.78]; per-seed
  SD comparable to the standing battery's 0.095 (informational).
* **GV-P3 (score):** pooled S (licensed axes) < a037's re-measured 1.683,
  near the Block-3 single-seed 1.559.

GV-P1/P2 failure ⇒ h405 not seed-robust ⇒ no adjudication fires. All
three pass ⇒ the G4 adjudications become FIREABLE (user calls) — nothing
adopts here.

**Do not read the launch print as n̄** — `n_detect_mean` lacks the §4r
suppressed→bin-0 convention (measured up to 2.9 He hot on the G3 ring).
Only the Block-V scorer's `nbar_det` is n̄.

Atlas stance: instrument runs — nothing here moves `finc1v725`.

Usage::

    python scripts/gen_tier2atlas_g4step2_battery.py --dry-run   # oracles only
    python scripts/gen_tier2atlas_g4step2_battery.py             # all 5 members
    python scripts/gen_tier2atlas_g4step2_battery.py s1 s4       # selected
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
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


# =============================================================================
# USER SETTINGS
# =============================================================================

# Five fresh seeds (plan §3.5f; verified unused by any committed generator —
# the committed range ends at 20260729).
MEMBER_SEEDS: dict[str, int] = {
    "s1": 20260730,
    "s2": 20260731,
    "s3": 20260732,
    "s4": 20260733,
    "s5": 20260734,
}

DEFAULT_CONCURRENCY = 3
SKIP_COMPLETED_RUNS = True
# Complete members are skipped untouched, so a relaunch only ever rebuilds a
# member that died mid-stage (the G1/G3/G4 safe-relaunch precedent).
OVERWRITE_EXISTING_RUN = True


# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2atlas_g4finals import (  # noqa: E402
    CASE,
    N,
    TWIN_ROWS,
    TWIN_SCAN_CSV,
    VARIANT,
    _SPEC_BY_LABEL,
    _TWIN_COLUMNS,
    build_cell,
    finals_run_dir_name,
    verify_corrected_geometry,
)
from scripts.tier0_common import run_dir_name  # noqa: E402
from scripts.tier2_common import cfg_diff_vs_reference  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

H405 = _SPEC_BY_LABEL["h405"]
COMMITTED_H405_RUN = finals_run_dir_name("h405")   # the Block-3 seed-20260729 run

_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "cfg.json", "neutral.npz", "ion.npz", "relaxation.npz", "detection.npz",
)


def verify_h405_twin_row() -> None:
    """Launch oracle: h405's frozen twin row reproduces string-exact from the
    committed fine-scan record (the finals G4F-P1 oracle, h405 slice)."""
    csv_path = PROJECT_ROOT / TWIN_SCAN_CSV
    key = (str(H405.v_c), H405.eb_tag, str(H405.tau_ps), str(H405.e0_eV))
    with open(csv_path, newline="") as fh:
        idx = {
            (r["v_c"], r["E_bind_tag"], r["tau_ps"], r["E0_eV"]): r
            for r in csv.DictReader(fh)
        }
    row = idx.get(key)
    if row is None:
        raise AssertionError(f"h405 twin oracle FAILED: no scan row at {key}.")
    got = tuple(row[c] for c in _TWIN_COLUMNS)
    if got != TWIN_ROWS["h405"]:
        raise AssertionError(
            f"h405 twin oracle FAILED: committed scan row {got} != frozen "
            f"pre-registration {TWIN_ROWS['h405']}."
        )
    print("h405 twin oracle PASSED: the frozen row reproduces string-exact "
          "from the committed fine-scan record.", flush=True)


def build_member(member: str) -> SimConfig:
    """h405's frozen cell, seed swapped — nothing else may move."""
    cfg = dataclasses.replace(build_cell(H405), seed=MEMBER_SEEDS[member])
    cfg.validate()
    verify_corrected_geometry(cfg, H405)
    return cfg


def verify_member_vs_committed_h405(cfg: SimConfig, member: str) -> list[str]:
    """Oracle: the member cfg diffs vs the committed g4fh405 cfg.json in
    exactly ``{"seed"}`` — same cell, fresh seed, no drifted pin."""
    ref = PROJECT_ROOT / "data" / "runs" / COMMITTED_H405_RUN / "cfg.json"
    diff = cfg_diff_vs_reference(cfg, ref, context=f"g4s2 {member}")
    if set(diff) != {"seed"}:
        raise AssertionError(
            f"[{member}] cfg vs committed h405 differs in {sorted(diff)}; "
            f"pre-registered exactly ['seed']."
        )
    return diff


def member_run_dir_name(member: str) -> str:
    """Run-dir basename in the atlas namespace (no ``_tier2_`` substring)."""
    if member not in MEMBER_SEEDS:
        raise KeyError(f"unknown member {member!r}")
    name = run_dir_name(
        CASE, VARIANT, N, run_tag=f"tier2atlas_conf270_g4s2h405{member}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _run_one(member: str) -> str:
    cfg = build_member(member)
    verify_member_vs_committed_h405(cfg, member)
    run_dir = PROJECT_ROOT / "data" / "runs" / member_run_dir_name(member)
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{member}] skip (complete) -> {run_dir}", flush=True)
            return member
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(f"{run_dir} exists and is not complete.")
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{member}] h405 (v_c={H405.v_c} {H405.eb_tag} tau={H405.tau_ps} "
          f"E0={H405.e0_eV}) seed={cfg.seed}", flush=True)
    print(f"[{member}] neutral ...", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{member}] ion ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    print(f"[{member}] relaxation (E2) ...", flush=True)
    relax = run_relaxation_stage(ion, cfg, save_path=run.root / "relaxation.npz")
    print(f"[{member}] detection ...", flush=True)
    detect = run_detection_stage(
        relax.checkpoint, cfg, save_path=run.root / "detection.npz")
    dm = detect.detected_mask
    nd = float(detect.n_detected[dm].mean()) if dm.any() else float("nan")
    print(f"[{member}] done -> {run_dir} (n_detect_mean={nd:.2f} [NOT nbar — "
          f"§4r], droplet_retained={int(np.count_nonzero(~dm))}/{dm.size})",
          flush=True)
    return member


def _select(labels: list[str]) -> list[str]:
    if not labels:
        return list(MEMBER_SEEDS)
    unknown = sorted(set(labels) - set(MEMBER_SEEDS))
    if unknown:
        raise SystemExit(f"unknown member(s) {unknown}; expected among "
                         f"{list(MEMBER_SEEDS)}")
    return labels


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    members = _select(args.labels)
    print(f"Atlas G4 Step 2 Block V battery: {len(members)} member(s) at "
          f"h405, N={N}, seeds {sorted(MEMBER_SEEDS[m] for m in members)}, "
          f"concurrency={args.concurrency}"
          f"{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument runs -- standing point does not move).", flush=True)

    verify_h405_twin_row()

    if args.dry_run:
        for member in members:
            cfg = build_member(member)
            diff = verify_member_vs_committed_h405(cfg, member)
            print(f"[{member}] cfg OK seed={cfg.seed} -> "
                  f"{member_run_dir_name(member)}", flush=True)
            print(f"[{member}]   diff vs committed h405: {sorted(diff)}",
                  flush=True)
        print("dry run complete: twin oracle + all cfg guards pass.",
              flush=True)
        return 0

    concurrency = max(1, args.concurrency)
    if concurrency == 1 or len(members) == 1:
        for member in members:
            _run_one(member)
        return 0

    from concurrent.futures import ProcessPoolExecutor, as_completed
    workers = min(concurrency, len(members))
    print(f"launching {len(members)} members through a {workers}-slot pool ...",
          flush=True)
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as poolx:
        futures = {poolx.submit(_run_one, m): m for m in members}
        for fut in as_completed(futures):
            member = futures[fut]
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[{member}] FAILED: {exc!r}", flush=True)
    if failures:
        print(f"{failures}/{len(members)} member(s) failed.", flush=True)
        return 1
    print(f"all {len(members)} members complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
