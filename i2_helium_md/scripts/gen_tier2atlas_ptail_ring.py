"""Atlas §3.5h — the p_tail ring at the h405 pins (3 × N = 1000).

Runs the `TIER2_SENSITIVITY_ATLAS_PLAN` §3.5h cells: the drag-tail exponent
opened as the low-n KE axis. Each cell is the committed `g4fh405` config
**verbatim** — corrected geometry, (v_c 5.5, τ 4.4, E₀ 0.405), seed
20260729 (the finals seed, so every comparison is CRN-paired) — with only
`p_tail` moved off the −1 convention:

* **pt15** (−1.5) — below-target bracket (placement KE₁ 0.889),
* **pt20** (−2.0) — at/above-target bracket (placement 1.062; the
  KE₁ = 1.00 anchor sits at p ≈ −1.80),
* **pt30** (−3.0) — far probe: overshoot + midHot-damage measurement.

**Oracle (§1.4, enforced before any MD and by --dry-run):** each cell's cfg
diffs against the committed `g4fh405` `cfg.json` in **exactly**
`{"drag_coefficients"}`, and inside the bundle only `p_tail` differs
(b, v_c and all provenance stamps bit-identical). The corrected-geometry
stamps are re-checked field by field (the finals guard, reused).

The conditional E₀-recenter cell (plan §3.5h PT-P4) is NOT here — it fires
only per its pre-registered clause, after the scorer has read these three.

Atlas stance: instrument runs — nothing here moves `finc1v725`.

Usage::

    python scripts/gen_tier2atlas_ptail_ring.py --dry-run   # oracles only
    python scripts/gen_tier2atlas_ptail_ring.py             # all 3 cells
    python scripts/gen_tier2atlas_ptail_ring.py pt20        # selected cells
"""

from __future__ import annotations

import argparse
import dataclasses
import json
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
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2atlas_g4finals import (  # noqa: E402
    _SPEC_BY_LABEL as _FINALS_SPECS,
    N,
    SEED,
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

# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"
VARIANT = "shared_pure_cubic"

# The committed h405 run this ring is CRN-paired against (its cfg.json is
# the cfg-diff oracle reference; its seed is inherited via build_cell).
H405_RUN = finals_run_dir_name("h405")

DEFAULT_CONCURRENCY = 3
SKIP_COMPLETED_RUNS = True
OVERWRITE_EXISTING_RUN = True

_REQUIRED_ARTIFACTS = ("cfg.json", "neutral.npz", "ion.npz",
                       "relaxation.npz", "detection.npz")


class PtailCell(NamedTuple):
    """One §3.5h ring cell."""

    label: str
    p_tail: float
    role: str


RING: tuple[PtailCell, ...] = (
    PtailCell("pt15", -1.5, "below-target bracket (placement KE1 0.889)"),
    PtailCell("pt20", -2.0, "at/above-target bracket (placement 1.062)"),
    PtailCell("pt30", -3.0, "far probe: overshoot + midHot damage"),
)
_SPEC_BY_LABEL = {c.label: c for c in RING}


def ptail_run_dir_name(label: str) -> str:
    """Run-dir basename in the atlas namespace (no ``_tier2_`` substring)."""
    name = run_dir_name(
        CASE, VARIANT, N, run_tag=f"tier2atlas_conf270_{label}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_ptail_cell(cell: PtailCell) -> SimConfig:
    """The committed h405 config with only the tail exponent moved."""
    cfg = build_cell(_FINALS_SPECS["h405"])
    bundle = cfg.drag_coefficients
    coefficients = {**bundle.coefficients, "p_tail": cell.p_tail}
    cfg = dataclasses.replace(
        cfg,
        drag_coefficients=dataclasses.replace(bundle,
                                              coefficients=coefficients),
    )
    cfg.validate()
    return cfg


def verify_ptail_oracle(cfg: SimConfig, cell: PtailCell) -> None:
    """§3.5h cfg-diff oracle: exactly {"drag_coefficients"}, p_tail-only."""
    ref_path = PROJECT_ROOT / "data" / "runs" / H405_RUN / "cfg.json"
    diff = cfg_diff_vs_reference(cfg, ref_path, context=f"ptail {cell.label}")
    if set(diff) != {"drag_coefficients"}:
        raise AssertionError(
            f"[{cell.label}] cfg diff vs committed h405 is {sorted(diff)}; "
            f"pre-registered exactly ['drag_coefficients']."
        )
    ref_bundle = json.loads(ref_path.read_text(encoding="utf-8"))[
        "drag_coefficients"]
    got = dataclasses.asdict(cfg.drag_coefficients)
    for key, ref_val in ref_bundle.items():
        if key == "coefficients":
            continue
        if got[key] != ref_val:
            raise AssertionError(
                f"[{cell.label}] bundle field {key!r} drifted: "
                f"{got[key]!r} vs committed {ref_val!r}."
            )
    for name, ref_val in ref_bundle["coefficients"].items():
        if name == "p_tail":
            continue
        if float(got["coefficients"][name]) != float(ref_val):
            raise AssertionError(
                f"[{cell.label}] coefficient {name!r} drifted: "
                f"{got['coefficients'][name]!r} vs committed {ref_val!r}."
            )
    if float(got["coefficients"]["p_tail"]) != cell.p_tail:
        raise AssertionError(
            f"[{cell.label}] p_tail {got['coefficients']['p_tail']!r} != "
            f"cell value {cell.p_tail!r}."
        )


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _run_one(cell: PtailCell) -> str:
    cfg = build_ptail_cell(cell)
    verify_corrected_geometry(cfg, _FINALS_SPECS["h405"])
    verify_ptail_oracle(cfg, cell)
    run_dir = PROJECT_ROOT / "data" / "runs" / ptail_run_dir_name(cell.label)
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{cell.label}] skip (complete) -> {run_dir}", flush=True)
            return cell.label
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(f"{run_dir} exists and is not complete.")
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{cell.label}] p_tail={cell.p_tail} ({cell.role}) "
          f"seed={cfg.seed}", flush=True)
    print(f"[{cell.label}] neutral ...", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{cell.label}] ion ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    print(f"[{cell.label}] relaxation (E2) ...", flush=True)
    relax = run_relaxation_stage(ion, cfg, save_path=run.root / "relaxation.npz")
    print(f"[{cell.label}] detection ...", flush=True)
    detect = run_detection_stage(
        relax.checkpoint, cfg, save_path=run.root / "detection.npz")
    dm = detect.detected_mask
    nd = float(detect.n_detected[dm].mean()) if dm.any() else float("nan")
    print(f"[{cell.label}] done -> {run_dir} (n_detect_mean={nd:.2f}, "
          f"droplet_retained={int(np.count_nonzero(~dm))}/{dm.size})",
          flush=True)
    return cell.label


def _run_one_by_label(label: str) -> str:
    return _run_one(_SPEC_BY_LABEL[label])


def _select(labels: list[str]) -> list[PtailCell]:
    if not labels:
        return list(RING)
    unknown = sorted(set(labels) - set(_SPEC_BY_LABEL))
    if unknown:
        raise SystemExit(f"unknown label(s) {unknown}; expected among "
                         f"{[c.label for c in RING]}")
    return [_SPEC_BY_LABEL[label] for label in labels]


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cells = _select(args.labels)
    print(f"Atlas 3.5h p_tail ring: {len(cells)} cell(s), N={N}, seed={SEED} "
          f"(= the finals seed, CRN-paired), h405 pins, "
          f"concurrency={args.concurrency}"
          f"{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument runs -- standing point does not move).", flush=True)

    if args.dry_run:
        for cell in cells:
            cfg = build_ptail_cell(cell)
            verify_corrected_geometry(cfg, _FINALS_SPECS["h405"])
            verify_ptail_oracle(cfg, cell)
            print(f"[{cell.label}] cfg OK p_tail={cell.p_tail} -> "
                  f"{ptail_run_dir_name(cell.label)}", flush=True)
        print("dry run complete: all cfg oracles pass.", flush=True)
        return 0

    concurrency = max(1, args.concurrency)
    if concurrency == 1 or len(cells) == 1:
        for cell in cells:
            _run_one(cell)
        return 0

    from concurrent.futures import ProcessPoolExecutor, as_completed
    with ProcessPoolExecutor(max_workers=concurrency) as pool:
        futures = {pool.submit(_run_one_by_label, cell.label): cell.label
                   for cell in cells}
        for future in as_completed(futures):
            future.result()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
