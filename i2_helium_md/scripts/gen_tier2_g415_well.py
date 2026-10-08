"""Fixed-well MD check of the lq A5 cell (thesis discussion 2026-10-06).

D0 §1.G4.14: the high-n tail is set by the escape well E_bind. The lq ridge
point c = 2.0 refitted on 18 A at the fixed shallow well E_bind = 0.080 eV
(``method_b_g415_fixed_well.py``; band pass, RMSE 0.159 A/ps) replaces the
committed (a, E_bind) of the A5 confirmation cell (c 2.0, tau 6.4). Everything
else is the A5 cell exactly (``gen_tier2_g414_tauline.build`` path: h405
production member, Stage I 500 ps, E2 off, partner-aware closure, N = 1000).

Runs (production seed s1, paired with the committed A5 member
``g414c200t64confs1``):

* ``w080e3155`` -- A5 E0 0.3155 eV;
* ``w080e3255`` -- E0 0.3255 eV (the shallower well raises nbar by about
  +0.27, CALIBRATION_MAP dnbar/dE_bind ~ -7.3 /eV; E0 +0.01 eV ~ -0.25).

Oracles before any MD (``--dry-run`` runs them only):

1. the fixed-well JSON passes its bands and is the c = 2.0 point;
2. each cfg diffs against the committed A5 s1 cfg in exactly the drag
   coefficient / binding keys (plus E0 for the second run) -- any other key
   aborts;
3. handover alignment (the production generator's check, reused).

ion.npz is kept (per-ion exit reads); two runs ~0.7 GB.

Usage::

    python scripts/gen_tier2_g415_well.py --dry-run
    python scripts/gen_tier2_g415_well.py
    python scripts/gen_tier2_g415_well.py --well-only --members s2 s3 s4 s5 r6

The last form completes the well-only arm to the six production members
(pooled N = 6000 per well: ``build_pooled_detection_container.py --target
g415_wo080`` / ``g415_wo048``).
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
from pathlib import Path
import sys
from typing import Optional

for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2_detfix_battery import (  # noqa: E402
    MAX_BYTES_ION, MEMBER_SEEDS, verify_handover_alignment,
)
from scripts.gen_tier2_g413_ridge_battery import DRAG_FIELDS, LQ_BUNDLE, S_RHO_A, _E0_KEY  # noqa: E402
from scripts.gen_tier2_g414_tauline import (  # noqa: E402
    ION_DROPPED_MARKER, Cell, build as a5_build, confirm_cells, ridge_point,
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
# USER SETTINGS
# =============================================================================
WELL_JSON = (PROJECT_ROOT / "data/reference/drag/18A/trajectory_matching_srho3p5A"
             / "g415_fixed_well/c2p000_e0p0800.json")
MEMBER = "s1"   # default member (the committed s1 runs; the scorer reads it)
# Well-only pooled arm (user, 2026-10-08): the remaining production members
# for an N = 6000 pool per well (``--members``); each member pairs with the
# committed A5 member of the same seed.
POOL_MEMBERS = tuple(MEMBER_SEEDS)   # s1-s5 + r6
E0_EV = {"w080e3155": None, "w080e3255": 0.3255}   # None = the A5 cell's own E0
# Informative well-only arm (user, 2026-10-06): the A5 law (committed c2p000
# a, c) with ONLY E_bind lowered -- NOT Tier-0-consistent (the 18 A refit
# moves a); A5 E0. ion.npz dropped after detection verifies (disk).
WELL_ONLY_EV = {"wo080": 0.080, "wo048": 0.0482}
MAX_CONCURRENCY = 2


def well_point() -> dict:
    """The fixed-well 18 A point (read, never retyped); must pass its bands."""
    rec = json.loads(WELL_JSON.read_text(encoding="utf-8"))
    if rec["c_amu_per_A"] != 2.0 or not rec["bands_pass"]:
        raise AssertionError(f"{WELL_JSON.name} is not a passing c = 2.0 point")
    return {"a": float(rec["a_amu_per_ps"]), "c": 2.0, "e_bind_eV": float(rec["e_bind_eV"])}


def a5_cell(member: str = MEMBER) -> Cell:
    return next(c for c in confirm_cells() if c.label.endswith(f"conf{member}"))


def build(label: str, member: str = MEMBER) -> SimConfig:
    """The A5 member cfg with the fixed-well (a, E_bind) and, optionally, E0;
    for a well-only label, the A5 law with only E_bind replaced."""
    base = a5_build(a5_cell(member))
    if label in WELL_ONLY_EV:
        wp = {**ridge_point(2.0), "e_bind_eV": WELL_ONLY_EV[label]}
    else:
        wp = well_point()
    drag = build_drag_cfg(
        CASE, LQ_BUNDLE, num_molecules=base.num_molecules,
        ion_time_ps=base.ion_simulation_time, dt_ion_ps=base.dt_ion,
        seed=base.seed, coeff_overrides={"a": wp["a"], "c": wp["c"]},
        e_bind_override=wp["e_bind_eV"])
    repl = {k: getattr(drag, k) for k in DRAG_FIELDS}
    if repl["drag_gate_steepness"] != S_RHO_A:
        raise AssertionError(f"bundle gate {repl['drag_gate_steepness']} != {S_RHO_A}")
    if E0_EV.get(label) is not None:
        repl[_E0_KEY] = float(E0_EV[label]) / BUDGET_EV
    cfg = dataclasses.replace(base, **repl)
    cfg.validate()
    return cfg


def verify(label: str, member: str = MEMBER) -> list[str]:
    """Oracle 2: the diff vs the committed A5 cfg of the same member is
    coefficient / binding (+E0) only; the seed is that member's production seed."""
    base = a5_build(a5_cell(member))
    cfg = build(label, member)
    diff = sorted(cfg_field_diff(base, cfg))
    want = {"drag_coefficients", "binding_energy_I_ion_eV"} | (
        {_E0_KEY} if E0_EV.get(label) is not None else set())
    if set(diff) != want:
        raise AssertionError(f"[{label}{member}] cfg diff vs A5 {member}: {diff}; "
                             f"expected {sorted(want)}")
    if cfg.seed != MEMBER_SEEDS[member]:
        raise AssertionError(f"[{label}{member}] seed {cfg.seed} != production {member}")
    return diff


def run_dir_for(label: str, member: str = MEMBER) -> Path:
    name = run_dir_name(CASE, LQ_BUNDLE, N, run_tag=f"detfix_conf270_g415{label}{member}_th500")
    return PROJECT_ROOT / "data" / "runs" / name


def run_one(label: str, member: str = MEMBER) -> str:
    verify(label, member)
    cfg = build(label, member)
    steps, stride = verify_handover_alignment(cfg)
    run_dir = run_dir_for(label, member)
    tag = f"{label}{member}"
    if (run_dir / "detection.npz").exists():
        print(f"[{tag}] skip (complete) -> {run_dir.name}", flush=True)
        return tag
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{tag}] E0 {getattr(cfg, _E0_KEY) * BUDGET_EV:.4f} seed {cfg.seed}, "
          f"{steps} steps (stride {stride}) -> {run_dir.name}", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False, max_bytes=MAX_BYTES_ION)
    t_last = float(np.asarray(ion.time_ps)[-1])
    if not math.isclose(t_last, cfg.ion_simulation_time - 0.005, abs_tol=1e-9):
        raise AssertionError(f"[{tag}] last stored column at {t_last} ps")
    det = run_detection_stage(ion, cfg, save_path=run.root / "detection.npz")
    n_scored = int(det.detected_mask.sum())
    print(f"[{tag}] done: scored {n_scored}/{det.detected_mask.size}", flush=True)
    del ion, neutral
    if label in WELL_ONLY_EV and n_scored > 0:
        ion_path = run.root / "ion.npz"
        (run.root / ION_DROPPED_MARKER).write_text(json.dumps(
            {"ion_npz_bytes": ion_path.stat().st_size, "detection_scored": n_scored,
             "reason": "g415 well-only arm: disk (T: full); scorer reads cfg.json + detection.npz"},
            indent=2), encoding="utf-8")
        ion_path.unlink()
    return tag


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--well-only", action="store_true", help="the informative well-only arm")
    parser.add_argument("--members", nargs="+", default=[MEMBER], choices=POOL_MEMBERS,
                        help="production members to run (default s1; the pool: s1 s2 s3 s4 s5 r6)")
    args = parser.parse_args(argv)
    labels = list(WELL_ONLY_EV) if args.well_only else list(E0_EV)
    jobs = [(lab, m) for lab in labels for m in args.members]
    wp = ridge_point(2.0) if args.well_only else well_point()
    print(f"law: a {wp['a']!r} amu/ps, c {wp['c']}; E_bind "
          f"{list(WELL_ONLY_EV.values()) if args.well_only else wp['e_bind_eV']} eV", flush=True)
    for label, member in jobs:
        diff = verify(label, member)
        steps, stride = verify_handover_alignment(build(label, member))
        print(f"[{label}{member}] guards PASSED: diff {diff}; seed {MEMBER_SEEDS[member]}; "
              f"{steps} steps, stride {stride}", flush=True)
    if args.dry_run:
        return 0
    from concurrent.futures import ProcessPoolExecutor, as_completed
    failures = 0
    with ProcessPoolExecutor(max_workers=min(MAX_CONCURRENCY, len(jobs))) as poolx:
        futures = {poolx.submit(run_one, lab, m): f"{lab}{m}" for lab, m in jobs}
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
