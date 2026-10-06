"""G4.14 Step 1 — MD tau lines for the ridge law (production pipeline, independent seeds).

THESIS_FINAL_GOALS.md G4.14 (frozen 7726e43; do NOT re-tune). The law
F = g(rho_hat) * (a*v + c*v^2), each c's (a, E_bind) read from its committed
ridge JSON (``g412_ridge/c*.json``), density gate ``erf_independent`` 3.5 A,
no cap. Everything else is the production pipeline exactly as G4.13
(``build_member`` = h405's frozen cell, Stage I to 500 ps, E2 off,
partner-aware closure, N = 1000).

Members:

* **h405g35 pool** — h405's own cell with the 3.5 A density gate (the G4.13
  ``h405g35`` definition, run through the G4.13 generator) at s2-s5 + r6;
  pooled with the committed s1 member by the scorer. Runs first.
* **Step 1** — c = 2.0: tau in {4.0, 4.4} U {4.8, 5.0, ..., 6.6}; every other
  c that passed the Step-0 midHot filter: tau in {4.8, 5.2, 5.6, 6.0, 6.4};
  each tau at the two Step-0 E0 values (predicted MD nbar 3.85 / 4.05, read
  from the committed ``atlas_g414_twin_line.csv``, never retyped).
  **Seeds:** one fresh seed per run, consecutive from 20261101 in plan order
  (c = 2.0 first, tau ascending, nbar target ascending; then the other c
  ascending), disjoint from the production seeds. No CRN.

Oracles before any MD (``--dry-run`` runs them only):

1. each cfg diffs against the h405 detfix s1 member in exactly the law keys,
   tau, E0 and seed (any other key aborts);
2. handover alignment (the production generator's check, reused);
3. the seed block is consecutive, unique and disjoint from every production
   seed.

Usage::

    python scripts/gen_tier2_g414_tauline.py --dry-run
    python scripts/gen_tier2_g414_tauline.py --h405g35     # reference pool
    python scripts/gen_tier2_g414_tauline.py               # Step 1, 2 at a time
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
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

from scripts.gen_tier2_detfix_battery import (  # noqa: E402
    MAX_BYTES_ION, MEMBER_SEEDS, build_member, verify_handover_alignment,
)
from scripts.gen_tier2_g413_ridge_battery import (  # noqa: E402
    DRAG_FIELDS, H405_GATE_LABEL, LQ_BUNDLE, S_RHO_A, STAGE_B_MEMBERS,
    _E0_KEY, _LAW_KEYS, _TAU_KEY, run_one as g413_run_one,
    verify_cell as g413_verify_cell,
)
from scripts.gen_tier2atlas_g4finals import BUDGET_EV, CASE, N  # noqa: E402
from scripts.gen_tier2atlas_linclone import cfg_field_diff  # noqa: E402
from scripts.tier0_common import build_drag_cfg, run_dir_name  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

# =============================================================================
# USER SETTINGS (pre-registered G4.14; FIXED)
# =============================================================================
RIDGE_DIR = PROJECT_ROOT / "data/reference/drag/18A/trajectory_matching_srho3p5A/g412_ridge"
TWIN_DIR = PROJECT_ROOT / "data/runs/h2b_forward_model"
LINE_CSV = TWIN_DIR / "atlas_g414_twin_line.csv"
FILTER_CSV = TWIN_DIR / "atlas_g414_filter.csv"

PRIMARY_C = 2.0
PRIMARY_TAU_PS = (4.0, 4.4) + tuple(float(t) for t in np.round(np.arange(4.8, 6.6001, 0.2), 1))
OTHER_TAU_PS = (4.8, 5.2, 5.6, 6.0, 6.4)
NBAR_TARGETS = (3.85, 4.05)
SEED_BASE = 20261101
# Amendment A5 (2026-10-06): one KE-led pooled confirmation at the production
# seeds; E0 = the fitted nbar = 3.89 point of the committed trade-off curve.
CONFIRM_C = 2.0
CONFIRM_TAU_PS = 6.4
CONFIRM_NBAR = 3.89
CONFIRM_E0_EV_4DP = 0.3155            # frozen in G4.14 A5; read value must match
TRADEOFF_CSV = TWIN_DIR / "g414_tradeoff.csv"
REFERENCE_MEMBER = "s1"          # the diff oracle's h405 detfix member

MAX_CONCURRENCY = 2
SKIP_COMPLETED_RUNS = True
# Disk (user decision 2026-10-06, T: full): the scorer reads only cfg.json +
# detection.npz, so each Step-1 run's 361 MB ion.npz is deleted once
# detection.npz reloads with scored ions; ION_DROPPED_MARKER records it.
DROP_ION_NPZ = True
ION_DROPPED_MARKER = "ion_npz_dropped.json"
_REQUIRED_ARTIFACTS = ("cfg.json", "neutral.npz", "detection.npz")


class Cell(NamedTuple):
    label: str
    c: float
    tau_ps: float
    nbar_target: float
    e0_eV: float
    seed: int


def cell_label(c: float, tau_ps: float, nbar_target: float) -> str:
    """e.g. (2.0, 4.8, 3.85) -> ``c200t48n385``."""
    return f"c{round(c * 100):03d}t{round(tau_ps * 10):d}n{round(nbar_target * 100):03d}"


def plan_cells(line_rows: list[dict], filter_rows: list[dict]) -> list[Cell]:
    """The frozen Step-1 run list with its seed block (plan order).

    ``line_rows`` / ``filter_rows``: the committed Step-0 CSV rows (strings).
    c = 2.0 always runs; another c runs iff ``passes_midhot_filter == 1``.
    Raises if a planned (c, tau, target) has no finite Step-0 E0.
    """
    passing = sorted(float(f["c"]) for f in filter_rows
                     if int(f["passes_midhot_filter"]) == 1 and float(f["c"]) != PRIMARY_C)
    plan = [(PRIMARY_C, PRIMARY_TAU_PS)] + [(c, OTHER_TAU_PS) for c in passing]
    e0_at = {(float(r["c"]), float(r["tau_ps"]), float(r["nbar_target_md"])): r["E0_eV"]
             for r in line_rows}
    cells, k = [], 0
    for c, taus in plan:
        for tau in taus:
            for target in NBAR_TARGETS:
                raw = e0_at.get((c, tau, target), "nan")
                e0 = float(raw)
                if not math.isfinite(e0):
                    raise ValueError(f"no Step-0 E0 for c={c}, tau={tau}, nbar {target}")
                cells.append(Cell(cell_label(c, tau, target), c, tau, target, e0,
                                  SEED_BASE + k))
                k += 1
    verify_seed_block([cl.seed for cl in cells])
    return cells


def verify_seed_block(seeds: list[int]) -> None:
    """Oracle 3: consecutive from SEED_BASE, unique, no production seed."""
    if seeds != list(range(SEED_BASE, SEED_BASE + len(seeds))):
        raise AssertionError("seed block is not consecutive from SEED_BASE")
    clash = set(seeds) & set(MEMBER_SEEDS.values())
    if clash:
        raise AssertionError(f"seed block reuses production seeds {sorted(clash)}")


def load_plan() -> list[Cell]:
    for p in (LINE_CSV, FILTER_CSV):
        if not p.exists():
            raise FileNotFoundError(f"{p} missing: run G4.14 Step 0 (g414twin) first")
    with open(LINE_CSV, newline="", encoding="utf-8") as fh:
        line_rows = list(csv.DictReader(fh))
    with open(FILTER_CSV, newline="", encoding="utf-8") as fh:
        filter_rows = list(csv.DictReader(fh))
    return plan_cells(line_rows, filter_rows)


def confirm_cells() -> list[Cell]:
    """The A5 confirmation: one cell per production member (s1-s5 + r6)."""
    with open(TRADEOFF_CSV, newline="", encoding="utf-8") as fh:
        row = next(r for r in csv.DictReader(fh)
                   if float(r["c"]) == CONFIRM_C and float(r["tau_ps"]) == CONFIRM_TAU_PS)
    e0 = float(row["E0_eV"])
    if round(e0, 4) != CONFIRM_E0_EV_4DP or abs(float(row["pred_nbar_det"]) - CONFIRM_NBAR) > 1e-9:
        raise AssertionError(f"confirmation E0 {e0!r} / nbar {row['pred_nbar_det']} "
                             f"!= frozen {CONFIRM_E0_EV_4DP} / {CONFIRM_NBAR}")
    tag = cell_label(CONFIRM_C, CONFIRM_TAU_PS, CONFIRM_NBAR).split("n")[0]
    return [Cell(f"{tag}conf{m}", CONFIRM_C, CONFIRM_TAU_PS, CONFIRM_NBAR, e0, seed)
            for m, seed in MEMBER_SEEDS.items()]


def ridge_point(c: float) -> dict:
    """A committed ridge point (read, never retyped); must pass its bands."""
    path = RIDGE_DIR / f"c{c:.3f}".replace(".", "p")
    rec = json.loads(path.with_suffix(".json").read_text(encoding="utf-8"))
    if rec["c_amu_per_A"] != c or not rec["bands_pass"]:
        raise AssertionError(f"{path}.json is not a passing c = {c} point")
    return {"a": float(rec["a_amu_per_ps"]), "c": c, "e_bind_eV": float(rec["e_bind_eV"])}


def run_dir_for(cell: Cell) -> Path:
    name = run_dir_name(CASE, LQ_BUNDLE, N, run_tag=f"detfix_conf270_g414{cell.label}_th500")
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"namespace violated by {name!r}")
    return PROJECT_ROOT / "data" / "runs" / name


def build(cell: Cell) -> SimConfig:
    """h405 production member, seed swapped, with the ridge law / tau / E0."""
    base = dataclasses.replace(build_member(REFERENCE_MEMBER), seed=cell.seed)
    rp = ridge_point(cell.c)
    drag = build_drag_cfg(
        CASE, LQ_BUNDLE, num_molecules=base.num_molecules,
        ion_time_ps=base.ion_simulation_time, dt_ion_ps=base.dt_ion,
        seed=base.seed, coeff_overrides={"a": rp["a"], "c": rp["c"]},
        e_bind_override=rp["e_bind_eV"])
    repl = {k: getattr(drag, k) for k in DRAG_FIELDS}
    if repl["drag_gate_steepness"] != S_RHO_A:
        raise AssertionError(f"bundle gate {repl['drag_gate_steepness']} != {S_RHO_A}")
    repl[_TAU_KEY] = float(cell.tau_ps)
    repl[_E0_KEY] = float(cell.e0_eV) / BUDGET_EV
    cfg = dataclasses.replace(base, **repl)
    cfg.validate()
    return cfg


def expected_diff_keys(cell: Cell) -> frozenset[str]:
    """Frozen diff vs the h405 detfix s1 member: law keys + seed, plus tau / E0
    whenever the cell's value differs from h405's (tau 4.4 is h405's own)."""
    base = build_member(REFERENCE_MEMBER)
    keys = set(_LAW_KEYS)
    if int(cell.seed) != int(base.seed):
        keys.add("seed")
    if float(cell.tau_ps) != float(getattr(base, _TAU_KEY)):
        keys.add(_TAU_KEY)
    if float(cell.e0_eV) / BUDGET_EV != float(getattr(base, _E0_KEY)):
        keys.add(_E0_KEY)
    return frozenset(keys)


def verify_cell(cell: Cell) -> list[str]:
    """Oracle 1: diff vs the h405 detfix s1 member is exactly the frozen set."""
    diff = sorted(cfg_field_diff(build_member(REFERENCE_MEMBER), build(cell)))
    want = expected_diff_keys(cell)
    if set(diff) != want:
        raise AssertionError(f"[{cell.label}] cfg diff vs h405 {REFERENCE_MEMBER} is "
                             f"{diff}; pre-registered {sorted(want)}")
    return diff


def _run_is_complete(run_dir: Path) -> bool:
    return (all((run_dir / n).exists() for n in _REQUIRED_ARTIFACTS)
            and ((run_dir / "ion.npz").exists() or (run_dir / ION_DROPPED_MARKER).exists()))


def drop_ion_npz(run_dir: Path) -> None:
    """Delete ``ion.npz`` after checking ``detection.npz`` reloads with scored
    ions; writes :data:`ION_DROPPED_MARKER` first. No-op if already dropped."""
    from i2_helium_md.simulation.detection_stage import load_detection_result

    ion_path = run_dir / "ion.npz"
    if not ion_path.exists():
        return
    det = load_detection_result(run_dir / "detection.npz")
    n_scored = int(np.asarray(det.detected_mask).sum())
    if n_scored <= 0:
        raise AssertionError(f"{run_dir.name}: detection.npz has no scored ions; ion.npz kept")
    (run_dir / ION_DROPPED_MARKER).write_text(json.dumps(
        {"ion_npz_bytes": ion_path.stat().st_size, "detection_scored": n_scored,
         "reason": "G4.14 Step 1: disk (T: full); scorer reads cfg.json + detection.npz"},
        indent=2), encoding="utf-8")
    ion_path.unlink()


def run_cell(cell: Cell) -> str:
    verify_cell(cell)
    cfg = build(cell)
    steps, stride = verify_handover_alignment(cfg)
    run_dir = run_dir_for(cell)
    if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
        if DROP_ION_NPZ:
            drop_ion_npz(run_dir)
        print(f"[{cell.label}] skip (complete) -> {run_dir.name}", flush=True)
        return cell.label
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{cell.label}] c {cell.c} tau {cell.tau_ps} E0 {cell.e0_eV} seed {cfg.seed}, "
          f"{steps} steps (stride {stride}) -> {run_dir.name}", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False,
                              max_bytes=MAX_BYTES_ION)
    t_last = float(np.asarray(ion.time_ps)[-1])
    if not math.isclose(t_last, cfg.ion_simulation_time - 0.005, abs_tol=1e-9):
        raise AssertionError(f"[{cell.label}] last stored column at {t_last} ps")
    det = run_detection_stage(ion, cfg, save_path=run.root / "detection.npz")
    del ion, neutral
    if DROP_ION_NPZ:
        drop_ion_npz(run_dir)
    print(f"[{cell.label}] done: scored {int(det.detected_mask.sum())}/"
          f"{det.detected_mask.size}", flush=True)
    return cell.label


def _h405g35_job(member: str) -> str:
    return g413_run_one((H405_GATE_LABEL, member))


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cells", nargs="*", help="subset of Step-1 labels")
    parser.add_argument("--h405g35", action="store_true", help="the reference pool")
    parser.add_argument("--confirm", action="store_true", help="A5 pooled confirmation")
    parser.add_argument("--concurrency", type=int, default=MAX_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    if args.h405g35:
        fn, jobs = _h405g35_job, list(STAGE_B_MEMBERS)
        for m in jobs:
            print(f"[{H405_GATE_LABEL}/{m}] guards PASSED: diff "
                  f"{g413_verify_cell(H405_GATE_LABEL, m)}", flush=True)
    else:
        plan = confirm_cells() if args.confirm else load_plan()
        unknown = sorted(set(args.cells) - {c.label for c in plan})
        if unknown:
            raise SystemExit(f"unknown cell(s) {unknown}")
        jobs = [c for c in plan if not args.cells or c.label in args.cells]
        fn = run_cell
        for c in sorted({c.c for c in jobs}):
            rp = ridge_point(c)
            print(f"law c {c}: a {rp['a']!r} amu/ps, E_bind {rp['e_bind_eV']!r} eV, "
                  f"gate {S_RHO_A} A", flush=True)
        for cell in jobs:
            diff = verify_cell(cell)
            steps, stride = verify_handover_alignment(build(cell))
            print(f"[{cell.label}] E0 {cell.e0_eV} seed {cell.seed}: guards PASSED "
                  f"(diff {len(diff)} keys; {steps} steps, stride {stride})", flush=True)
        print(f"{len(jobs)} {'confirmation' if args.confirm else 'Step-1'} runs planned "
              f"(seeds {[j.seed for j in jobs]})", flush=True)
    if args.dry_run:
        return 0
    workers = max(1, min(args.concurrency, MAX_CONCURRENCY, len(jobs)))
    from concurrent.futures import ProcessPoolExecutor, as_completed
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as poolx:
        futures = {poolx.submit(fn, j): j for j in jobs}
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
