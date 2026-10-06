"""G4.13 — MD battery for the G4.12 ridge law c = 2.0, uncapped (production pipeline).

THESIS_FINAL_GOALS.md G4.13 (pre-registered 4e0908c; do NOT re-tune). The law
F = g(rho_hat) * (a*v + c*v^2) with the committed G4.12 ridge point c = 2.0
(a, E_bind read from ``g412_ridge/c2p000.json`` at build time), density gate
``erf_independent`` 3.5 A, no cap. Everything else is the production pipeline
(``gen_tier2_detfix_battery.py``): h405's frozen cell, Stage I to 500 ps, E2
off, partner-aware closure, N = 1000.

Members:

* **Stage A** (seed 20260730 = h405 detfix s1, CRN): the 16-cell (tau, E0)
  grid ring, plus ``h405g35`` = h405's own cell with the 3.5 A density gate.
* **Stage B**: one cell (chosen by the frozen rule in the scorer,
  ``tier2_g413_table.py``) at s2-s5 + r6 -> pooled N = 6000 with its Stage-A
  s1 member. Selected with ``--battery <cell>``.

Oracles before any MD (``--dry-run`` runs them only):

1. each cfg diffs against the h405 detfix member of the same seed in exactly
   :func:`expected_diff_keys` (any other key aborts);
2. handover alignment (the production generator's check, reused).

Usage::

    python scripts/gen_tier2_g413_ridge_battery.py --dry-run
    python scripts/gen_tier2_g413_ridge_battery.py              # Stage A, 2 at a time
    python scripts/gen_tier2_g413_ridge_battery.py --battery t52e370   # Stage B
"""

from __future__ import annotations

import argparse
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
    MAX_BYTES_ION, MEMBER_SEEDS, build_member, member_run_dir_name,
    verify_handover_alignment,
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
# USER SETTINGS (pre-registered G4.13; FIXED)
# =============================================================================
RIDGE_POINT_JSON = (PROJECT_ROOT / "data/reference/drag/18A/"
                    "trajectory_matching_srho3p5A/g412_ridge/c2p000.json")
LQ_BUNDLE = "srho3p5A18A_linear_quadratic"   # stamped gate 3.5 A; {a,c,E} overridden
S_RHO_A = 3.5
STAGE_A_MEMBER = "s1"                        # seed 20260730, CRN with h405 detfix s1
STAGE_B_MEMBERS = ("s2", "s3", "s4", "s5", "r6")

#: (tau [ps], E0 [eV]) grid -- G4.13 Stage-A table.
GRID: tuple[tuple[float, float], ...] = (
    (4.8, 0.39), (4.8, 0.40), (4.8, 0.41), (4.8, 0.42),
    (5.2, 0.36), (5.2, 0.37), (5.2, 0.38), (5.2, 0.39),
    (5.6, 0.335), (5.6, 0.345), (5.6, 0.355), (5.6, 0.365),
    (6.4, 0.30), (6.4, 0.31), (6.4, 0.32), (6.4, 0.33),
)
H405_GATE_LABEL = "h405g35"

MAX_CONCURRENCY = 2
SKIP_COMPLETED_RUNS = True
_REQUIRED_ARTIFACTS = ("cfg.json", "neutral.npz", "ion.npz", "detection.npz")

#: Drag-law fields the c = 2.0 swap replaces (from the Tier-0 drag builder).
DRAG_FIELDS = ("drag_form", "drag_coefficients", "drag_bundle",
               "binding_energy_I_ion_eV", "drag_spatial_gate",
               "drag_gate_steepness")
#: Pre-registered diff vs the h405 detfix member of the same seed.
_LAW_KEYS = frozenset(DRAG_FIELDS)
_GATE_KEYS = frozenset({"drag_spatial_gate", "drag_gate_steepness"})
_H405G35_KEYS = _GATE_KEYS | {"allow_inconsistent_gate_pairing"}
_TAU_KEY = "internal_energy_cooling_tau_ps"
_E0_KEY = "internal_energy_partition_fraction"


class Cell(NamedTuple):
    label: str
    tau_ps: Optional[float]    # None -> h405's own
    e0_eV: Optional[float]     # None -> h405's own
    ridge_law: bool            # False -> h405 law with the 3.5 A gate only


def cell_label(tau_ps: float, e0_eV: float) -> str:
    """e.g. (5.2, 0.37) -> ``t52e370``; (5.6, 0.335) -> ``t56e335``."""
    return f"t{round(tau_ps * 10):d}e{round(e0_eV * 1000):03d}"


CELLS: dict[str, Cell] = {
    **{cell_label(t, e): Cell(cell_label(t, e), t, e, True) for t, e in GRID},
    H405_GATE_LABEL: Cell(H405_GATE_LABEL, None, None, False),
}


def ridge_point() -> dict:
    """The committed G4.12 c = 2.0 ridge point (read, never retyped)."""
    rec = json.loads(RIDGE_POINT_JSON.read_text(encoding="utf-8"))
    if rec["c_amu_per_A"] != 2.0 or not rec["bands_pass"]:
        raise AssertionError(f"{RIDGE_POINT_JSON} is not the passing c = 2.0 point")
    return {"a": float(rec["a_amu_per_ps"]), "c": 2.0,
            "e_bind_eV": float(rec["e_bind_eV"])}


def run_dir_for(cell: str, member: str) -> Path:
    variant = LQ_BUNDLE if CELLS[cell].ridge_law else "shared_pure_cubic"
    name = run_dir_name(CASE, variant, N,
                        run_tag=f"detfix_conf270_g413{cell}{member}_th500")
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"namespace violated by {name!r}")
    return PROJECT_ROOT / "data" / "runs" / name


def build(cell: str, member: str) -> SimConfig:
    """h405 production member of this seed, with the G4.13 law / gate / knobs."""
    spec = CELLS[cell]
    base = build_member(member)
    if spec.ridge_law:
        rp = ridge_point()
        drag = build_drag_cfg(
            CASE, LQ_BUNDLE, num_molecules=base.num_molecules,
            ion_time_ps=base.ion_simulation_time, dt_ion_ps=base.dt_ion,
            seed=base.seed, coeff_overrides={"a": rp["a"], "c": rp["c"]},
            e_bind_override=rp["e_bind_eV"])
        repl = {k: getattr(drag, k) for k in DRAG_FIELDS}
        if repl["drag_gate_steepness"] != S_RHO_A:
            raise AssertionError(f"bundle gate {repl['drag_gate_steepness']} != {S_RHO_A}")
        repl[_TAU_KEY] = float(spec.tau_ps)
        repl[_E0_KEY] = float(spec.e0_eV) / BUDGET_EV
    else:
        # The h405 bundle was extracted at 14.2 A: moving only its density gate
        # trips the §6.5.2 pairing guard by design -> the documented hatch.
        repl = {"drag_spatial_gate": "erf_independent", "drag_gate_steepness": S_RHO_A,
                "allow_inconsistent_gate_pairing": True}
    cfg = dataclasses.replace(base, **repl)
    cfg.validate()
    return cfg


def expected_diff_keys(cell: str) -> frozenset[str]:
    """Pre-registered diff vs the h405 detfix member of the same seed."""
    spec = CELLS[cell]
    if not spec.ridge_law:
        return _H405G35_KEYS
    return _LAW_KEYS | {_TAU_KEY, _E0_KEY}


def verify_cell(cell: str, member: str) -> list[str]:
    """Oracle 1: the in-memory diff vs the h405 member is exactly the
    pre-registered key set."""
    diff = sorted(cfg_field_diff(build_member(member), build(cell, member)))
    if set(diff) != expected_diff_keys(cell):
        raise AssertionError(
            f"[{cell}/{member}] cfg diff vs h405 {member} is {diff}; "
            f"pre-registered {sorted(expected_diff_keys(cell))}")
    return diff


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / n).exists() for n in _REQUIRED_ARTIFACTS)


def run_one(job: tuple[str, str]) -> str:
    cell, member = job
    verify_cell(cell, member)
    cfg = build(cell, member)
    steps, stride = verify_handover_alignment(cfg)
    run_dir = run_dir_for(cell, member)
    tag = f"{cell}/{member}"
    if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
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
    if not math.isclose(t_last, cfg.ion_simulation_time - 0.005, abs_tol=1e-9):
        raise AssertionError(f"[{tag}] last stored column at {t_last} ps")
    det = run_detection_stage(ion, cfg, save_path=run.root / "detection.npz")
    print(f"[{tag}] done: scored {int(det.detected_mask.sum())}/"
          f"{det.detected_mask.size}", flush=True)
    return tag


def jobs_for(args) -> list[tuple[str, str]]:
    if args.battery:
        if args.battery not in CELLS or not CELLS[args.battery].ridge_law:
            raise SystemExit(f"--battery needs a ridge grid cell, got {args.battery!r}")
        return [(args.battery, m) for m in STAGE_B_MEMBERS]
    cells = args.cells or list(CELLS)
    unknown = sorted(set(cells) - set(CELLS))
    if unknown:
        raise SystemExit(f"unknown cell(s) {unknown}")
    return [(c, STAGE_A_MEMBER) for c in cells]


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cells", nargs="*")
    parser.add_argument("--battery", default=None)
    parser.add_argument("--concurrency", type=int, default=MAX_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    jobs = jobs_for(args)
    rp = ridge_point()
    print(f"G4.13 law: a {rp['a']!r} amu/ps, c {rp['c']} amu/A, "
          f"E_bind {rp['e_bind_eV']!r} eV, gate {S_RHO_A} A", flush=True)
    for cell, member in jobs:
        diff = verify_cell(cell, member)
        steps, stride = verify_handover_alignment(build(cell, member))
        print(f"[{cell}/{member}] guards PASSED: diff {diff}; {steps} steps, "
              f"stride {stride}", flush=True)
    if args.dry_run:
        return 0
    workers = max(1, min(args.concurrency, MAX_CONCURRENCY, len(jobs)))
    from concurrent.futures import ProcessPoolExecutor, as_completed
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as poolx:
        futures = {poolx.submit(run_one, j): j for j in jobs}
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
