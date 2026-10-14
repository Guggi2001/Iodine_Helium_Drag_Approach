"""RRK s_eff scan at the thesis production law B2 (seed s1, N = 1000 per cell).

Question (user, 2026-10-14): does the RRK effective dof ``s_eff``
(``evap_rrk_dof``, production 8) move the detected observables at B2, or only
the shedding clock? Every committed s_eff evidence (staircase probe I10, Wave 7
I20, Wave 9 I32) predates the lq ridge, the ``post_pickup`` sampler and the
absolute-E0 arm; ``TIER2_PARAMETER_INFLUENCE.md`` flags s_eff as transferred.

Cells: the committed B2 s1 member (``gen_tier2_detfix_postpickup.py --cell lq
--b2``, label ``g5lqppc160e3297s1``) with **one field moved**:
``evap_rrk_dof`` ∈ :data:`SEFF_GRID`. Same seed as B2 s1; s_eff = 8 is the
committed B2 s1 run itself (not re-run). The cells share the seed, but the
evaporation draws follow the rates, so downstream streams differ: read
differences against the per-member seed sd, not as paired deltas.

Oracles before any MD:

1. ``verify_seff_cell`` — the cfg differs from the committed B2 s1
   ``cfg.json`` in exactly ``{"evap_rrk_dof"}``; at s_eff = 8 it differs in
   nothing (wiring check of the builder against the committed run).
2. ``verify_handover_alignment`` (re-used) — last stored ion column at t_h.

Usage::

    python scripts/gen_tier2_b2_seff_scan.py --dry-run   # oracles only
    python scripts/gen_tier2_b2_seff_scan.py             # all cells, 2 at a time
    python scripts/gen_tier2_b2_seff_scan.py 12 30       # selected s_eff

**--keep-ion (user, 2026-10-14):** re-runs every cell incl. s_eff = 8 into
``…seff##s1keepion…`` dirs and keeps ``ion.npz`` (the in-window n_shell
history, for the shed-timing read against Albrechtsen et al., Nature 623,
319 (2023): "the vast majority of dissociation events are finished after
400 ps"). Oracle: each detection.npz equals its scan twin's (n_detected,
state_reason, event times) — the runs are deterministic at fixed seed.

    python scripts/gen_tier2_b2_seff_scan.py --keep-ion

**--pool (user, 2026-10-14):** s_eff = :data:`POOL_SEFF` at all six
production members (s1-s5 + r6, the B2 seeds; s1 = the scan cell), ion.npz
dropped; oracle 1 per member against that member's committed B2 cfg.json.
Pool: ``build_pooled_detection_container.py --target g5_lqppb2seff05``.

    python scripts/gen_tier2_b2_seff_scan.py --pool

**--pool-e0 (user, 2026-10-14):** the same six members with E0 re-placed to
:data:`POOL_E0_EV` to restore the B2 pooled nbar (3.903 vs 3.721 at s_eff 5,
E0 0.3297): dnbar/dE0 ~ -32 /eV from the committed B1 E0 column (s1,
post_pickup, c 2.0: 0.3155 -> 0.3355 eV, nbar 4.557 -> 3.910), so
dE0 = -0.182 / 32 = -0.0057 eV. Oracle 1 diff = {evap_rrk_dof, E0 key}.
Pool target ``g5_lqppb2seff05e3240``.

    python scripts/gen_tier2_b2_seff_scan.py --pool-e0

Scorer: ``scripts/post_processing/tier2_b2_seff_table.py``; shed timing:
``scripts/post_processing/tier2_b2_seff_timing.py``.
Max two concurrent MD runs.
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
    MAX_BYTES_ION, MEMBER_SEEDS, T_HANDOVER_PS, verify_handover_alignment,
)
from scripts.gen_tier2_detfix_postpickup import (  # noqa: E402
    B2_C, B2_E0_EV, B2_TAU_PS, LQ_DROP_REASON, _E0_KEY, b1c_e0, build_lq_member, lq_run_dir,
)
from scripts.gen_tier2_g414_tauline import ION_DROPPED_MARKER, drop_ion_npz  # noqa: E402
from scripts.gen_tier2atlas_g4finals import CASE, N  # noqa: E402
from scripts.gen_tier2_g413_ridge_battery import LQ_BUNDLE  # noqa: E402
from scripts.tier0_common import run_dir_name  # noqa: E402
from scripts.tier2_common import cfg_diff_vs_reference  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================

#: s_eff values run [dimensionless]; 8 = production (the committed B2 s1 run).
SEFF_GRID = (5.0, 12.0, 20.0, 30.0)
SEFF_PRODUCTION = 8.0
MEMBER = "s1"
SEFF_DIFF_KEYS = frozenset({"evap_rrk_dof"})
#: s_eff of the six-member pool (--pool).
POOL_SEFF = 5.0
#: E0 [eV] of the nbar-restoring pool (--pool-e0); see module docstring.
POOL_E0_EV = 0.3240
MAX_CONCURRENCY = 2

# =============================================================================


def b2_run_dir(member: str = MEMBER) -> Path:
    """The committed B2 run of ``member`` (the s_eff = 8 cell)."""
    return lq_run_dir(member, B2_E0_EV, B2_TAU_PS, B2_C)


def seff_run_dir(s_eff: float, keep_ion: bool = False, member: str = MEMBER,
                 e0_eV: Optional[float] = None) -> Path:
    """Run dir of an s_eff cell (optionally at E0 ``e0_eV`` [eV]); s_eff = 8
    at the B2 E0 maps to the committed B2 run unless ``keep_ion``."""
    if s_eff == SEFF_PRODUCTION and not keep_ion and e0_eV is None:
        return b2_run_dir(member)
    e0_tag = "" if e0_eV is None else f"e{round(e0_eV * 1e4):04d}"
    tag = (f"g5lqppc160e3297seff{round(s_eff):02d}{e0_tag}{member}"
           + ("keepion" if keep_ion else ""))
    return PROJECT_ROOT / "data" / "runs" / run_dir_name(
        CASE, LQ_BUNDLE, N, run_tag=f"detfix_conf270_{tag}_th500")


def build_seff_cell(s_eff: float, member: str = MEMBER,
                    e0_eV: Optional[float] = None) -> SimConfig:
    """B2 member cfg with ``evap_rrk_dof`` = s_eff [dimensionless, >= 1] and,
    optionally, E0 re-placed to ``e0_eV`` [eV]."""
    if B2_E0_EV != b1c_e0(B2_C):
        raise AssertionError(f"B2_E0_EV {B2_E0_EV} != b1c_e0({B2_C})")
    e0 = B2_E0_EV if e0_eV is None else e0_eV
    cfg = dataclasses.replace(build_lq_member(member, e0, B2_TAU_PS, B2_C),
                              evap_rrk_dof=float(s_eff))
    cfg.validate()
    return cfg


def verify_seff_cell(s_eff: float, member: str = MEMBER,
                     e0_eV: Optional[float] = None) -> list[str]:
    """Oracle 1: diff vs the committed B2 member cfg.json is exactly
    ``{"evap_rrk_dof"}`` (+ the E0 key when re-placed; empty at the
    production value). Returns the keys."""
    diff = cfg_diff_vs_reference(build_seff_cell(s_eff, member, e0_eV),
                                 b2_run_dir(member) / "cfg.json",
                                 context=f"B2 {member} s_eff {s_eff} E0 {e0_eV}")
    want = ((set() if s_eff == SEFF_PRODUCTION else set(SEFF_DIFF_KEYS))
            | (set() if e0_eV is None else {_E0_KEY}))
    if set(diff) != want:
        raise AssertionError(f"[s_eff {s_eff}] cfg vs B2 s1 differs in {sorted(diff)}; "
                             f"pre-registered exactly {sorted(want)}.")
    return diff


def verify_same_detection(s_eff: float) -> None:
    """Oracle (keep-ion): the re-run's detection.npz equals the scan twin's."""
    from i2_helium_md.simulation.detection_stage import load_detection_result
    a = load_detection_result(seff_run_dir(s_eff) / "detection.npz")
    b = load_detection_result(seff_run_dir(s_eff, keep_ion=True) / "detection.npz")
    for f in ("n_detected", "state_reason", "event_offsets", "event_time_ps"):
        if not np.array_equal(np.asarray(getattr(a, f)), np.asarray(getattr(b, f))):
            raise AssertionError(f"[s_eff {s_eff}] keep-ion re-run differs in {f}")


def run_cell(s_eff: float, keep_ion: bool = False, member: str = MEMBER,
             e0_eV: Optional[float] = None) -> float:
    """Guard, neutral -> ion (to t_h) -> detection; drop ion.npz unless keep_ion."""
    verify_seff_cell(s_eff, member, e0_eV)
    cfg = build_seff_cell(s_eff, member, e0_eV)
    steps, stride = verify_handover_alignment(cfg)
    run_dir = seff_run_dir(s_eff, keep_ion, member, e0_eV)
    done = all((run_dir / n).exists() for n in ("cfg.json", "neutral.npz", "detection.npz"))
    if keep_ion and done and (run_dir / "ion.npz").exists():
        verify_same_detection(s_eff)
        print(f"[s_eff {s_eff:g}] skip (complete, ion kept) -> {run_dir.name}", flush=True)
        return s_eff
    done = done and ((run_dir / "ion.npz").exists() or (run_dir / ION_DROPPED_MARKER).exists())
    if done and not keep_ion:
        drop_ion_npz(run_dir, reason=LQ_DROP_REASON)
        print(f"[s_eff {s_eff:g}] skip (complete) -> {run_dir.name}", flush=True)
        return s_eff
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[s_eff {s_eff:g}] seed={cfg.seed} ion {steps} steps (stride {stride}) "
          f"-> {run_dir.name}", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False,
                              max_bytes=MAX_BYTES_ION)
    t_last = float(np.asarray(ion.time_ps)[-1])
    if not math.isclose(t_last, T_HANDOVER_PS, abs_tol=1e-9):
        raise AssertionError(f"[s_eff {s_eff}] ion.npz last column at {t_last} ps, "
                             f"not t_h {T_HANDOVER_PS} ps.")
    det = run_detection_stage(ion, cfg, save_path=run.root / "detection.npz")
    del ion, neutral
    if keep_ion:
        verify_same_detection(s_eff)
    else:
        drop_ion_npz(run_dir, reason=LQ_DROP_REASON)
    print(f"[s_eff {s_eff:g}] done: scored {int(det.detected_mask.sum())}/"
          f"{det.detected_mask.size}", flush=True)
    return s_eff


def run_pool(dry_run: bool, e0_eV: Optional[float] = None) -> int:
    """s_eff = POOL_SEFF at every B2 member (E0 optionally re-placed);
    ion.npz dropped."""
    for m in MEMBER_SEEDS:
        diff = verify_seff_cell(POOL_SEFF, m, e0_eV)
        steps, stride = verify_handover_alignment(build_seff_cell(POOL_SEFF, m, e0_eV))
        print(f"[{m} s_eff {POOL_SEFF:g} E0 {e0_eV}] guards PASSED: diff {diff}, {steps} "
              f"steps, stride {stride} -> "
              f"{seff_run_dir(POOL_SEFF, member=m, e0_eV=e0_eV).name}", flush=True)
    if dry_run:
        return 0
    from concurrent.futures import ProcessPoolExecutor, as_completed
    failures = 0
    with ProcessPoolExecutor(max_workers=MAX_CONCURRENCY) as pool:
        futures = {pool.submit(run_cell, POOL_SEFF, False, m, e0_eV): m
                   for m in MEMBER_SEEDS}
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[{futures[fut]}] FAILED: {exc!r}", flush=True)
    print(f"{len(MEMBER_SEEDS) - failures}/{len(MEMBER_SEEDS)} members complete.", flush=True)
    return 1 if failures else 0


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seff", nargs="*", type=float)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--keep-ion", action="store_true",
                        help="re-run every cell incl. s_eff 8, keep ion.npz")
    parser.add_argument("--pool", action="store_true",
                        help="s_eff POOL_SEFF at all six B2 members (s1-s5 + r6)")
    parser.add_argument("--pool-e0", action="store_true",
                        help="as --pool, with E0 re-placed to POOL_E0_EV")
    args = parser.parse_args(argv)
    if args.pool or args.pool_e0:
        return run_pool(args.dry_run, POOL_E0_EV if args.pool_e0 else None)
    grid = tuple(args.seff) or (tuple(sorted((*SEFF_GRID, SEFF_PRODUCTION)))
                                if args.keep_ion else SEFF_GRID)
    unknown = sorted(set(grid) - set(SEFF_GRID) - ({SEFF_PRODUCTION} if args.keep_ion else set()))
    if unknown:
        raise SystemExit(f"s_eff {unknown} not in SEFF_GRID {SEFF_GRID}")

    verify_seff_cell(SEFF_PRODUCTION)
    print(f"[s_eff 8] builder == committed B2 s1 cfg ({b2_run_dir().name})", flush=True)
    for s in grid:
        diff = verify_seff_cell(s)
        steps, stride = verify_handover_alignment(build_seff_cell(s))
        print(f"[s_eff {s:g}] guards PASSED: diff {diff}, {steps} steps, stride {stride} "
              f"-> {seff_run_dir(s, args.keep_ion).name}", flush=True)
    if args.dry_run:
        return 0

    from concurrent.futures import ProcessPoolExecutor, as_completed
    failures = 0
    with ProcessPoolExecutor(max_workers=min(MAX_CONCURRENCY, len(grid))) as pool:
        futures = {pool.submit(run_cell, s, args.keep_ion): s for s in grid}
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[s_eff {futures[fut]:g}] FAILED: {exc!r}", flush=True)
    print(f"{len(grid) - failures}/{len(grid)} cells complete.", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
