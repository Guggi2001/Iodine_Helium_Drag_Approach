"""Detector-stage fix — validation MD (TIER2_DetectorStageFix.md §3b).

Re-runs two committed G4 finalists (**h405**, **b031**) with the proposed
pipeline: Stage I (full physics, pickup live, no Landau gate) extended to the
handover time ``t_h``, **E2 skipped** (``relaxation_stage_enabled=False``, the
detection stage seeds from ``ion.npz``), detection as-is. Everything else is
the finals generator's cell verbatim (``build_cell``), so the committed run is
the paired reference for the §3b sharp checks.

``t_h`` = 502.84 ps is chosen so that

* it equals a stored column of the committed ``relaxation.npz`` (E2 column
  245 = internal step 50284; stride 193 steps from the 29.99 ps seed), giving
  an equal-time per-ion comparison, and
* the Stage-I storage stride at the default byte budget (13) divides the last
  internal step (50284), so the checkpoint's last stored column — the one the
  detection stage seeds from — **is** the final state at exactly ``t_h``.

Both are asserted before any MD (``verify_handover_alignment``), as is the cfg
identity with the committed run (``verify_validation_cfg``: vs the committed
``cfg.json`` the base cell diffs in nothing, the validation cfg in exactly
``VALIDATION_DIFF_KEYS``).

No Coulomb closure is applied here; §2e shows it is exact post hoc under the
``co_moving`` shed convention, so it is evaluated on the finished runs.

Usage::

    python scripts/gen_detstagefix_validation.py --dry-run     # guards only
    python scripts/gen_detstagefix_validation.py h405           # one cell
    python scripts/gen_detstagefix_validation.py                # both, serial

Run the two cells as two separate processes for concurrency (max two
concurrent MD runs).
"""

from __future__ import annotations

import argparse
import dataclasses
import math
from pathlib import Path
import sys
from typing import Optional

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2atlas_g4finals import (  # noqa: E402
    CASE, FINAL_MATRIX, N, VARIANT, build_cell, finals_run_dir_name,
    verify_corrected_geometry,
)
from scripts.gen_tier2atlas_linclone import cfg_field_diff  # noqa: E402
from scripts.tier0_common import run_dir_name  # noqa: E402
from scripts.tier2_common import cfg_diff_vs_reference  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import (  # noqa: E402
    DEFAULT_MAX_CHECKPOINT_BYTES_ION, _decide_stride_ion,
    _internal_step_count_ion, run_ion_propagation,
)
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================

VALIDATION_LABELS = ("h405", "b031")

# Handover: E2 column 245 of the committed runs = internal step 50284.
T_HANDOVER_PS = 502.84
# ceil(502.845 / 0.01) = 50285 internal steps -> last step id 50284 = t_h.
ION_TIME_PS = 502.845
E2_COLUMN = 245

# The only keys allowed to differ from the committed finals cfg.
VALIDATION_DIFF_KEYS = frozenset({"ion_simulation_time", "relaxation_stage_enabled"})

_REQUIRED_ARTIFACTS = ("cfg.json", "neutral.npz", "ion.npz", "detection.npz")
SKIP_COMPLETED_RUNS = True

_SPEC_BY_LABEL = {c.label: c for c in FINAL_MATRIX}


def validation_run_dir_name(label: str) -> str:
    """Run-dir basename; distinct from the committed finals namespace."""
    return run_dir_name(CASE, VARIANT, N, run_tag=f"detfix_conf270_{label}_th503")


def committed_run_dir(label: str) -> Path:
    """The committed finals run this validation cell is paired with."""
    return PROJECT_ROOT / "data" / "runs" / finals_run_dir_name(label)


def build_validation_cfg(label: str) -> SimConfig:
    """The finals cell with Stage I extended to t_h and E2 switched off."""
    cfg = dataclasses.replace(
        build_cell(_SPEC_BY_LABEL[label]),
        ion_simulation_time=ION_TIME_PS,
        relaxation_stage_enabled=False,
    )
    cfg.validate()
    return cfg


def verify_validation_cfg(label: str) -> list[str]:
    """Base cell == committed cfg.json; validation cfg differs only in
    :data:`VALIDATION_DIFF_KEYS`. Returns the validation diff keys."""
    cell = _SPEC_BY_LABEL[label]
    base = build_cell(cell)
    verify_corrected_geometry(base, cell)
    ref = committed_run_dir(label) / "cfg.json"
    base_diff = cfg_diff_vs_reference(base, ref, context=label)
    if base_diff:
        raise AssertionError(
            f"[{label}] finals cell no longer reproduces the committed cfg.json: "
            f"diff {sorted(base_diff)}."
        )
    diff = cfg_field_diff(base, build_validation_cfg(label))
    if set(diff) != VALIDATION_DIFF_KEYS:
        raise AssertionError(
            f"[{label}] validation cfg diff {diff}; expected exactly "
            f"{sorted(VALIDATION_DIFF_KEYS)}."
        )
    return diff


def verify_handover_alignment(cfg: SimConfig, label: str) -> tuple[int, int]:
    """Assert Stage I's last stored column is the final step at t_h and that
    t_h is the committed E2 column :data:`E2_COLUMN`. Returns (steps, stride)."""
    steps = _internal_step_count_ion(cfg)
    last_id = steps - 1
    if not math.isclose(last_id * cfg.dt_ion, T_HANDOVER_PS, abs_tol=1e-9):
        raise AssertionError(
            f"[{label}] last step {last_id} * dt = {last_id * cfg.dt_ion} ps "
            f"!= t_h {T_HANDOVER_PS} ps."
        )
    stride, _ = _decide_stride_ion(
        cfg.num_molecules, steps, DEFAULT_MAX_CHECKPOINT_BYTES_ION,
    )
    if last_id % stride:
        raise AssertionError(
            f"[{label}] stride {stride} does not divide the last step {last_id}: "
            "the detection seed column would not be the final state."
        )
    with np.load(committed_run_dir(label) / "relaxation.npz",
                 allow_pickle=False) as z:
        t_e2 = float(z["time_ps"][E2_COLUMN])
    if not math.isclose(t_e2, T_HANDOVER_PS, abs_tol=1e-6):
        raise AssertionError(
            f"[{label}] committed E2 column {E2_COLUMN} is at {t_e2} ps, "
            f"not t_h {T_HANDOVER_PS} ps."
        )
    return steps, stride


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def run_one(label: str) -> None:
    """Build, guard, and run one validation cell (neutral -> ion -> detection)."""
    verify_validation_cfg(label)
    cfg = build_validation_cfg(label)
    steps, stride = verify_handover_alignment(cfg, label)
    run_dir = PROJECT_ROOT / "data" / "runs" / validation_run_dir_name(label)
    if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
        print(f"[{label}] skip (complete) -> {run_dir}", flush=True)
        return
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{label}] seed={cfg.seed} ion {steps} steps (stride {stride}) to "
          f"t_h={T_HANDOVER_PS} ps, E2 off -> {run_dir}", flush=True)
    print(f"[{label}] neutral ...", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{label}] ion ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    print(f"[{label}] detection (skip path, seeded from ion.npz) ...", flush=True)
    detect = run_detection_stage(ion, cfg, save_path=run.root / "detection.npz")
    dm = detect.detected_mask
    print(f"[{label}] done: scored {int(dm.sum())}/{dm.size}, "
          f"retained {int(np.count_nonzero(~dm))}", flush=True)


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    labels = args.labels or list(VALIDATION_LABELS)
    unknown = sorted(set(labels) - set(VALIDATION_LABELS))
    if unknown:
        raise SystemExit(f"unknown label(s) {unknown}; expected among "
                         f"{list(VALIDATION_LABELS)}")
    for label in labels:
        diff = verify_validation_cfg(label)
        steps, stride = verify_handover_alignment(build_validation_cfg(label), label)
        print(f"[{label}] guards PASSED: cfg diff {diff}, {steps} steps, "
              f"stride {stride}, t_h {T_HANDOVER_PS} ps = E2 column {E2_COLUMN}",
              flush=True)
    if args.dry_run:
        return 0
    for label in labels:
        run_one(label)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
