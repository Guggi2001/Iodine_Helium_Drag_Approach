"""G9 cold-onset minority screen — E0 columns at the lq A5 cell (seed s1, CRN).

THESIS_FINAL_GOALS.md G9 (frozen 2026-10-08). The A5 law (c 2.0, a / E_bind
from the committed ridge JSON, 3.5 A gate, tau 6.4, uncapped), production
pipeline, N = 1000, seed s1 (the A5 s1 births). E0 grid step 0.0125 eV
around the A5 value: offsets k = -4 .. +2 (0.2655 .. 0.3405 eV); k = 0 is
the committed A5 s1 member (not re-run) -> 6 new runs. Mixtures are built by
the scorer (``tier2_g9_cold_onset_table.py``), zero MD.

Oracles before any MD (``--dry-run`` runs them only):

1. the rebuilt A5 s1 cfg equals its stored ``cfg.json``;
2. each column cfg diffs from it in exactly
   ``internal_energy_partition_fraction`` (= E0 / 2.70);
3. handover alignment (the production generator's check, reused).

Usage::

    python scripts/gen_tier2_g9_cold_onset.py --dry-run
    python scripts/gen_tier2_g9_cold_onset.py          # 6 runs, 2 at a time
"""

from __future__ import annotations

import argparse
import dataclasses
import math
import os
from pathlib import Path
import sys
from typing import NamedTuple, Optional

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

from scripts.gen_tier2_detfix_battery import MAX_BYTES_ION, verify_handover_alignment  # noqa: E402
from scripts.gen_tier2_g413_ridge_battery import LQ_BUNDLE, _E0_KEY  # noqa: E402
from scripts.gen_tier2_g414_tauline import (  # noqa: E402
    ION_DROPPED_MARKER, build as g414_build, confirm_cells, drop_ion_npz,
    run_dir_for as g414_run_dir_for,
)
from scripts.gen_tier2atlas_g4finals import CASE, N  # noqa: E402
from scripts.gen_tier2atlas_linclone import cfg_field_diff  # noqa: E402
from scripts.tier0_common import run_dir_name  # noqa: E402
from scripts.tier2_common import cfg_diff_vs_reference  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402

# =============================================================================
# USER SETTINGS (pre-registered G9; FIXED)
# =============================================================================
E0_STEP_EV = 0.0125
OFFSETS = (-4, -3, -2, -1, 0, 1, 2)       # grid nodes; 0 = the committed A5 s1 member
PARENT_MEMBER = "s1"
MAX_CONCURRENCY = 2
_REQUIRED_ARTIFACTS = ("cfg.json", "neutral.npz", "detection.npz")


class Column(NamedTuple):
    k: int
    e0_eV: float
    label: str


def parent_cell():
    """The A5 confirmation cell of the s1 member (G4.14)."""
    return next(c for c in confirm_cells() if c.label.endswith(PARENT_MEMBER))


def columns() -> list[Column]:
    """The frozen grid: E0 = E0_A5 + k * 0.0125 eV, k in OFFSETS."""
    e0_a5 = parent_cell().e0_eV
    out = []
    for k in OFFSETS:
        e0 = e0_a5 + k * E0_STEP_EV
        out.append(Column(k, e0, f"g9e{round(e0 * 1e4):04d}{PARENT_MEMBER}"))
    return out


def build(col: Column) -> SimConfig:
    return g414_build(parent_cell()._replace(e0_eV=col.e0_eV))


def run_dir_for(col: Column) -> Path:
    if col.k == 0:
        return g414_run_dir_for(parent_cell())
    name = run_dir_name(CASE, LQ_BUNDLE, N, run_tag=f"detfix_conf270_{col.label}_th500")
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"namespace violated by {name!r}")
    return PROJECT_ROOT / "data" / "runs" / name


def verify(col: Column) -> list[str]:
    """Oracles 1 and 2."""
    parent = g414_build(parent_cell())
    stored = cfg_diff_vs_reference(parent, g414_run_dir_for(parent_cell()) / "cfg.json",
                                   context="G9 parent A5 s1")
    if stored:
        raise AssertionError(f"rebuilt A5 s1 differs from its cfg.json in {stored}")
    diff = cfg_field_diff(parent, build(col))
    want = [] if col.k == 0 else [_E0_KEY]
    if diff != want:
        raise AssertionError(f"[{col.label}] cfg diff {diff}; pre-registered {want}")
    return diff


def _run_is_complete(run_dir: Path) -> bool:
    return (all((run_dir / n).exists() for n in _REQUIRED_ARTIFACTS)
            and ((run_dir / "ion.npz").exists() or (run_dir / ION_DROPPED_MARKER).exists()))


def run_column(col: Column) -> str:
    from i2_helium_md.simulation.detection_stage import run_detection_stage
    from i2_helium_md.simulation.ion import run_ion_propagation
    from i2_helium_md.simulation.neutral import run_neutral_propagation
    from i2_helium_md.simulation.run_directory import RunDirectory

    verify(col)
    cfg = build(col)
    steps, stride = verify_handover_alignment(cfg)
    run_dir = run_dir_for(col)
    if _run_is_complete(run_dir):
        drop_ion_npz(run_dir)
        print(f"[{col.label}] skip (complete) -> {run_dir.name}", flush=True)
        return col.label
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{col.label}] E0 {col.e0_eV:.5f} seed {cfg.seed}, {steps} steps "
          f"(stride {stride}) -> {run_dir.name}", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False,
                              max_bytes=MAX_BYTES_ION)
    t_last = float(np.asarray(ion.time_ps)[-1])
    if not math.isclose(t_last, cfg.ion_simulation_time - 0.005, abs_tol=1e-9):
        raise AssertionError(f"[{col.label}] last stored column at {t_last} ps")
    det = run_detection_stage(ion, cfg, save_path=run.root / "detection.npz")
    del ion, neutral
    drop_ion_npz(run_dir)
    print(f"[{col.label}] done: scored {int(det.detected_mask.sum())}/"
          f"{det.detected_mask.size}", flush=True)
    return col.label


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    jobs = []
    for col in columns():
        diff = verify(col)
        steps, stride = verify_handover_alignment(build(col))
        status = "committed A5 s1" if col.k == 0 else "new"
        print(f"[{col.label}] k {col.k:+d} E0 {col.e0_eV:.5f} ({status}): guards PASSED "
              f"(diff {diff}; {steps} steps, stride {stride})", flush=True)
        if col.k != 0:
            jobs.append(col)
    if args.dry_run:
        return 0
    from concurrent.futures import ProcessPoolExecutor, as_completed
    failures = 0
    with ProcessPoolExecutor(max_workers=MAX_CONCURRENCY) as poolx:
        futures = {poolx.submit(run_column, j): j.label for j in jobs}
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
