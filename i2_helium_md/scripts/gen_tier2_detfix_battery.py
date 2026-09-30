"""Detector-stage fix — the new thesis reference battery at h405 (production pipeline).

TIER2_DetectorStageFix.md §5.2 / §5.4 D6. The production pipeline:
**Stage I** (full physics, pickup live, no Landau gate) to the fixed handover
``t_h = 500 ps`` → **E2 skipped** (``relaxation_stage_enabled=False``; the
detection stage seeds from ``ion.npz``) → **detection stage** with the exact
partner-aware residual-Coulomb closure at handover
(``detection_coulomb_closure="partner_aware"``) and the marginal-partner
safeguard.

Members:

* ``s1``–``s5`` — seeds 20260730–20260734, each CRN-paired with its E2-era
  twin in ``gen_tier2atlas_g4step2_battery.py`` (``…g4s2h405sK``). Pooled
  N = 5000 = the new thesis reference.
* ``r6`` — seed 20260729, paired with the committed finals ``g4fh405`` and
  the t_h = 502.84 ps validation run (``gen_detstagefix_validation.py``);
  answers "does 500 vs 502.84 ps matter". Not part of the reference bundle.

**Zero pin duplication:** each cfg is ``gen_tier2atlas_g4finals.build_cell(h405)``
with the seed swapped plus the three pipeline keys. Oracles before any MD:

1. ``verify_member_vs_twin`` — with the closure set back to its default, the
   cfg diffs against the twin's committed ``cfg.json`` in exactly
   :data:`TWIN_DIFF_KEYS`, and the closure is the only further in-memory
   difference. (``cfg_diff_vs_reference`` refuses a non-default value in a
   field the older ``cfg.json`` cannot witness, hence the two-step oracle.)
2. ``verify_handover_alignment`` — 50001 internal steps, last step id
   50000 × dt = 500.00 ps; ``max_bytes`` = :data:`MAX_BYTES_ION` gives
   storage stride 16 (smallest divisor of 50000 at or above the default 13),
   so the last stored column — the detection seed — **is** the final state at
   t_h. Re-asserted on the finished ``ion.npz`` before the detection stage.

The pipeline identity is the cfg itself (``relaxation_stage_enabled``,
``detection_coulomb_closure``; §5.4 D5 — no stamp field).

Usage::

    python scripts/gen_tier2_detfix_battery.py --dry-run     # oracles only
    python scripts/gen_tier2_detfix_battery.py               # all 6, 2 at a time
    python scripts/gen_tier2_detfix_battery.py s1 r6         # selected

Max two concurrent MD runs (``--concurrency`` is capped at 2).
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

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2atlas_g4finals import (  # noqa: E402
    CASE, N, VARIANT, build_cell, finals_run_dir_name, verify_corrected_geometry,
)
from scripts.gen_tier2atlas_g4step2_battery import (  # noqa: E402
    H405, MEMBER_SEEDS as E2_BATTERY_SEEDS, member_run_dir_name as e2_member_run_dir_name,
)
from scripts.gen_tier2atlas_linclone import cfg_field_diff  # noqa: E402
from scripts.tier0_common import run_dir_name  # noqa: E402
from scripts.tier2_common import cfg_diff_vs_reference  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import (  # noqa: E402
    _NUM_2N_T_ARRAYS_ION, _decide_stride_ion, _internal_step_count_ion,
    run_ion_propagation,
)
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================

#: s1–s5 = the E2-era battery seeds (CRN twins); r6 = the finals / validation seed.
MEMBER_SEEDS: dict[str, int] = {**E2_BATTERY_SEEDS, "r6": 20260729}

T_HANDOVER_PS = 500.0
# ceil(500.005 / 0.01) = 50001 internal steps -> last step id 50000 = t_h.
ION_TIME_PS = 500.005
STORAGE_STRIDE = 16
# Budget for exactly 3126 stored columns -> _decide_stride_ion gives
# ceil(50001 / 3126) = 16 (3126..3333 columns all map to 16); ~750 MB at N=1000.
MAX_BYTES_ION = _NUM_2N_T_ARRAYS_ION * 2 * N * 8 * 3126

#: The production pipeline keys.
PIPELINE = dict(
    ion_simulation_time=ION_TIME_PS,
    relaxation_stage_enabled=False,
    detection_coulomb_closure="partner_aware",
)
#: Keys the twin's cfg.json can witness (closure reset to its default).
TWIN_DIFF_KEYS = frozenset({"ion_simulation_time", "relaxation_stage_enabled"})

MAX_CONCURRENCY = 2
SKIP_COMPLETED_RUNS = True
_REQUIRED_ARTIFACTS = ("cfg.json", "neutral.npz", "ion.npz", "detection.npz")


def _select(labels: list[str]) -> list[str]:
    if not labels:
        return list(MEMBER_SEEDS)
    unknown = sorted(set(labels) - set(MEMBER_SEEDS))
    if unknown:
        raise SystemExit(f"unknown member(s) {unknown}; expected among "
                         f"{list(MEMBER_SEEDS)}")
    return labels


def member_run_dir_name(member: str) -> str:
    """Run-dir basename; distinct from the E2-era ``tier2atlas_`` namespace."""
    if member not in MEMBER_SEEDS:
        raise KeyError(f"unknown member {member!r}")
    return run_dir_name(CASE, VARIANT, N,
                        run_tag=f"detfix_conf270_h405{member}_th500")


def twin_run_dir(member: str) -> Path:
    """The E2-era run this member is CRN-paired with."""
    name = (finals_run_dir_name("h405") if member == "r6"
            else e2_member_run_dir_name(member))
    return PROJECT_ROOT / "data" / "runs" / name


def build_member(member: str) -> SimConfig:
    """h405's frozen cell, seed swapped, production pipeline keys set."""
    cfg = dataclasses.replace(build_cell(H405), seed=MEMBER_SEEDS[member], **PIPELINE)
    cfg.validate()
    verify_corrected_geometry(cfg, H405)
    return cfg


def verify_member_vs_twin(member: str) -> list[str]:
    """Oracle 1 (see module docstring). Returns the witnessed diff keys."""
    cfg = build_member(member)
    witnessed = dataclasses.replace(cfg, detection_coulomb_closure="none")
    diff = cfg_diff_vs_reference(witnessed, twin_run_dir(member) / "cfg.json",
                                 context=f"detfix {member}")
    if set(diff) != TWIN_DIFF_KEYS:
        raise AssertionError(
            f"[{member}] cfg vs E2-era twin differs in {sorted(diff)}; "
            f"pre-registered exactly {sorted(TWIN_DIFF_KEYS)}."
        )
    extra = cfg_field_diff(witnessed, cfg)
    if set(extra) != {"detection_coulomb_closure"}:
        raise AssertionError(f"[{member}] unexpected in-memory diff {extra}.")
    return diff


def verify_handover_alignment(cfg: SimConfig) -> tuple[int, int]:
    """Oracle 2: the last stored column is the final step at t_h. Returns
    (internal steps, stride)."""
    steps = _internal_step_count_ion(cfg)
    last_id = steps - 1
    if not math.isclose(last_id * cfg.dt_ion, T_HANDOVER_PS, abs_tol=1e-9):
        raise AssertionError(
            f"last step {last_id} * dt = {last_id * cfg.dt_ion} ps != t_h "
            f"{T_HANDOVER_PS} ps.")
    stride, _ = _decide_stride_ion(cfg.num_molecules, steps, MAX_BYTES_ION)
    if stride != STORAGE_STRIDE or last_id % stride:
        raise AssertionError(
            f"stride {stride} (expected {STORAGE_STRIDE}) must divide the last "
            f"step {last_id}: the detection seed column would not be the final "
            "state.")
    return steps, stride


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def run_one(member: str) -> str:
    """Guard, then neutral -> ion (to t_h) -> detection (skip path + closure)."""
    verify_member_vs_twin(member)
    cfg = build_member(member)
    steps, stride = verify_handover_alignment(cfg)
    run_dir = PROJECT_ROOT / "data" / "runs" / member_run_dir_name(member)
    if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
        print(f"[{member}] skip (complete) -> {run_dir}", flush=True)
        return member
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{member}] seed={cfg.seed} ion {steps} steps (stride {stride}) to "
          f"t_h={T_HANDOVER_PS} ps, E2 off, closure partner_aware -> {run_dir}",
          flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{member}] ion ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False,
                              max_bytes=MAX_BYTES_ION)
    t_last = float(np.asarray(ion.time_ps)[-1])
    if not math.isclose(t_last, T_HANDOVER_PS, abs_tol=1e-9):
        raise AssertionError(
            f"[{member}] ion.npz last stored column at {t_last} ps, not t_h "
            f"{T_HANDOVER_PS} ps.")
    print(f"[{member}] detection (seeded from ion.npz at {t_last} ps) ...",
          flush=True)
    detect = run_detection_stage(ion, cfg, save_path=run.root / "detection.npz")
    dm = detect.detected_mask
    print(f"[{member}] done: scored {int(dm.sum())}/{dm.size}, retained "
          f"{int(np.count_nonzero(~dm))}", flush=True)
    return member


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=MAX_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    members = _select(args.labels)

    for member in members:
        diff = verify_member_vs_twin(member)
        steps, stride = verify_handover_alignment(build_member(member))
        print(f"[{member}] guards PASSED: seed {MEMBER_SEEDS[member]}, twin diff "
              f"{diff} + closure, {steps} steps, stride {stride} -> "
              f"{member_run_dir_name(member)}", flush=True)
    if args.dry_run:
        return 0

    workers = max(1, min(args.concurrency, MAX_CONCURRENCY, len(members)))
    if workers == 1:
        for member in members:
            run_one(member)
        return 0
    from concurrent.futures import ProcessPoolExecutor, as_completed
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as poolx:
        futures = {poolx.submit(run_one, m): m for m in members}
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[{futures[fut]}] FAILED: {exc!r}", flush=True)
    if failures:
        print(f"{failures}/{len(members)} member(s) failed.", flush=True)
        return 1
    print(f"all {len(members)} members complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
