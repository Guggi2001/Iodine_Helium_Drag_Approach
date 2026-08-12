"""Atlas free-form linear — the (a, τ) joint ring (plan §6.6).

Runs the `TIER2_FREEFORM_LINEAR_TWIN_SWEEP_PLAN` §6.6 ring: **6 cells ×
N = 500 at seed 20260731** — deliberately the §6.3 ring seed, so every
cell is CRN-paired ion-for-ion with lr1–lr7, the `h405p` partner **and**
the §6.5 clone member s1. All cells `eb0482`, corrected geometry,
standing pins otherwise.

The ring exists because the §6.5 anatomy isolated a mechanism: the
clone's W₁ damage is entirely a **tail** deficit, the tail was
**evaporated rather than trapped** (per-ion CRN flow: heavily-dressed
ions fall from mean n 18.3 → 9.6), E₀ is measured *not* to be the width
knob, and `a` improves W₁ monotonically — leaving **τ** as the owner by
elimination. τ *downward* is unavailable (τ 3.2 ⇒ twin n̄ 6.7–16.0).

Cells (plan §6.6 table; per-cell knobs = (a, τ, E₀), everything else the
standing pins at the corrected geometry):

* `s37`   a 37.5 / τ 4.8 / E₀ 0.38 — a-curve
* `s40`   a 40.0 / τ 4.8 / E₀ 0.39 — a-curve
* `s425`  a 42.5 / τ 4.8 / E₀ 0.39 — **decomposition** (the clone's a at
  the ring's clock)
* `s45`   a 45.0 / τ 4.8 / E₀ 0.40 — a-trend saturation
* `d2`    a 35.0 / τ 6.4 / E₀ 0.30 — **τ mirror** (lr6's a at the clone's
  clock)
* `j1`    a 40.0 / τ 6.4 / E₀ 0.31 — **the joint candidate**

`j1` carries twin `gate = 0` because its twin n̄ 4.392 sits just under the
*twin*-side floor 4.4; under the validated −0.55 transfer it lands inside
the MD-side band [3.77, 4.37]. Recorded so the twin flag is not misread.

**Oracles (enforced before any MD and by --dry-run):**

1. **LJ-P1** — all six frozen twin rows reproduce **string-exact** from
   the committed `atlas_linsweep.csv`, plus the ring's own LR-P1.
2. **CRN guard** — the seed is *shared* with the committed `h405p`
   partner, so `seed` and `num_molecules` must both be **absent** from
   every cell's cfg diff against it; the remaining diff must be exactly
   the pre-registered key set.
3. The standing-battery cfg-diff guard, reused verbatim from the ring.

Cell construction and the corrected-geometry / free-form provenance
guards are **imported** from `gen_tier2atlas_linring.py` (rule 1).

Atlas stance: instrument runs — nothing here moves `finc1v725`/h405.

Usage::

    python scripts/gen_tier2atlas_linjoint.py --dry-run   # oracles + guards
    python scripts/gen_tier2atlas_linjoint.py             # all 6 cells
    python scripts/gen_tier2atlas_linjoint.py s425 d2     # selected cells
"""

from __future__ import annotations

import argparse
import csv
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

from scripts.gen_tier2atlas_linclone import (  # noqa: E402
    CRN_FORBIDDEN_DIFF_KEYS,
    PARTNER_DIFF_KEYS,
    PARTNER_RUN,
)
from scripts.gen_tier2atlas_linring import (  # noqa: E402
    LINSWEEP_CSV,
    LinRingCell,
    N,
    SEED,
    build_cell,
    verify_against_reference,
    verify_corrected_geometry,
    verify_twin_preregistration,
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
# USER SETTINGS — the frozen §6.6 cell table
# =============================================================================

class JointCell(NamedTuple):
    """One §6.6 ring cell (all `pure_linear`, all eb0482)."""

    label: str
    a: float            # amu/ps
    tau_ps: float
    e0_eV: float
    role: str


JOINT_MATRIX: tuple[JointCell, ...] = (
    JointCell("s37", 37.5, 4.8, 0.38, "a-curve"),
    JointCell("s40", 40.0, 4.8, 0.39, "a-curve"),
    JointCell("s425", 42.5, 4.8, 0.39, "decomposition: clone a at ring clock"),
    JointCell("s45", 45.0, 4.8, 0.40, "a-trend saturation"),
    JointCell("d2", 35.0, 6.4, 0.30, "tau mirror: lr6 a at clone clock"),
    JointCell("j1", 40.0, 6.4, 0.31, "the joint candidate"),
)
_SPEC_BY_LABEL = {c.label: c for c in JOINT_MATRIX}

EB_TAG = "eb0482"

# LJ-P1: the frozen twin rows (committed atlas_linsweep.csv strings,
# compared EXACTLY). Columns below.
JOINT_TWIN_COLS = ("trapped_frac", "suppressed_frac", "nbar_det", "n1_solv",
                   "w1_solv", "n1_ke_eV", "ke2_eV", "deepke", "midhot_arith",
                   "phi", "gate")
JOINT_TWIN_ROWS: dict[str, tuple[str, ...]] = {
    "s37": ("0.0215", "0.1527", "4.564", "0.2027", "0.6027", "0.6365",
            "0.567", "1.5087", "1.3384", "0.869", "1"),
    "s40": ("0.0332", "0.1789", "4.506", "0.2029", "0.6657", "0.5748",
            "0.5071", "1.2285", "1.1691", "0.927", "1"),
    "s425": ("0.0476", "0.1644", "4.72", "0.193", "0.8079", "0.5412",
             "0.4737", "1.059", "1.0694", "0.985", "1"),
    "s45": ("0.0661", "0.1873", "4.641", "0.1938", "0.8795", "0.487",
            "0.4227", "0.8788", "0.9346", "1.043", "1"),
    "d2": ("0.0128", "0.026", "4.423", "0.1923", "0.8031", "0.7688",
           "0.669", "1.1753", "1.4372", "0.811", "1"),
    # j1's twin gate is 0 ONLY because twin n̄ 4.392 sits under the
    # twin-side floor 4.4; MD-side it lands in [3.77, 4.37] (plan §6.6).
    "j1": ("0.0332", "0.0518", "4.392", "0.2107", "0.6814", "0.6497",
           "0.5491", "0.7515", "1.1091", "0.927", "0"),
}

DEFAULT_CONCURRENCY = 3
SKIP_COMPLETED_RUNS = True
OVERWRITE_EXISTING_RUN = True

_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "cfg.json", "neutral.npz", "ion.npz", "relaxation.npz", "detection.npz",
)


# =============================================================================
# Oracles and guards
# =============================================================================


def as_ring_cell(cell: JointCell) -> LinRingCell:
    """This ring's cell in the §6.3 constructor's vocabulary (rule 1 — the
    cfg itself is built by `gen_tier2atlas_linring.build_cell`)."""
    return LinRingCell(cell.label, cell.a, None, EB_TAG, cell.tau_ps,
                       cell.e0_eV, cell.role)


def verify_joint_preregistration() -> None:
    """LJ-P1: every frozen twin row reproduces string-exact from the
    committed ``atlas_linsweep.csv``, and the ring's LR-P1 still passes."""
    verify_twin_preregistration()
    with open(PROJECT_ROOT / LINSWEEP_CSV, newline="", encoding="utf-8") as fh:
        idx = {
            (r["arm"], r["a"], r["E_bind_tag"], r["tau_ps"], r["E0_eV"]): r
            for r in csv.DictReader(fh)
        }
    for cell in JOINT_MATRIX:
        key = ("lin", f"{cell.a:.1f}", EB_TAG, str(cell.tau_ps),
               str(cell.e0_eV))
        row = idx.get(key)
        if row is None:
            raise AssertionError(
                f"[{cell.label}] LJ-P1 FAILED: no committed linsweep row at "
                f"{key}."
            )
        got = tuple(row[col] for col in JOINT_TWIN_COLS)
        if got != JOINT_TWIN_ROWS[cell.label]:
            raise AssertionError(
                f"[{cell.label}] LJ-P1 FAILED: committed linsweep row {got} "
                f"!= frozen pre-registration {JOINT_TWIN_ROWS[cell.label]}."
            )
    print(f"LJ-P1 joint oracle PASSED: all {len(JOINT_TWIN_ROWS)} frozen "
          "(a, tau) rows reproduce string-exact from the committed "
          "atlas_linsweep.csv.", flush=True)


def joint_run_dir_name(label: str) -> str:
    """Run-dir basename in the atlas namespace (no ``_tier2_`` substring)."""
    name = run_dir_name(
        "9A", "shared_pure_cubic", N,
        run_tag=f"tier2atlas_conf270_linj{label}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_joint_cell(label: str) -> SimConfig:
    """One ring cell's validated config (seed = the shared §6.3 CRN seed)."""
    if label not in _SPEC_BY_LABEL:
        raise KeyError(f"unknown cell {label!r}; expected among "
                       f"{sorted(_SPEC_BY_LABEL)}")
    cfg = build_cell(as_ring_cell(_SPEC_BY_LABEL[label]))
    if cfg.seed != SEED:
        raise AssertionError(
            f"[{label}] seed {cfg.seed} != the shared CRN seed {SEED} — the "
            "whole-family pairing depends on this."
        )
    cfg.validate()
    return cfg


def verify_crn_pairing(cfg: SimConfig, label: str) -> list[str]:
    """The §6.6 CRN guard.

    The seed is *shared* with the committed `h405p` partner, so unlike
    §6.5 neither `seed` nor `num_molecules` may appear in the diff for any
    cell — their identity is the whole-family pairing.
    """
    ref_path = PROJECT_ROOT / "data" / "runs" / PARTNER_RUN / "cfg.json"
    diff = cfg_diff_vs_reference(cfg, ref_path, context=label)
    hit = sorted(CRN_FORBIDDEN_DIFF_KEYS & set(diff))
    if hit:
        raise AssertionError(
            f"[{label}] CRN guard FAILED: {hit} differ(s) vs {PARTNER_RUN} — "
            "the whole-family CRN pairing is exactly that identity."
        )
    if set(diff) != set(PARTNER_DIFF_KEYS):
        raise AssertionError(
            f"[{label}] cfg diff vs the partner is {sorted(diff)}; "
            f"pre-registered exactly {sorted(PARTNER_DIFF_KEYS)}."
        )
    return diff


def verify_joint_cell(label: str) -> tuple[SimConfig, list[str]]:
    """All guards for one cell; returns the cfg and its partner diff."""
    cell = as_ring_cell(_SPEC_BY_LABEL[label])
    cfg = build_joint_cell(label)
    verify_corrected_geometry(cfg, cell)
    verify_against_reference(cfg, cell)   # vs the standing battery
    return cfg, verify_crn_pairing(cfg, label)


# =============================================================================
# Runner
# =============================================================================


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _run_one(label: str) -> str:
    cell = _SPEC_BY_LABEL[label]
    cfg, _ = verify_joint_cell(label)
    run_dir = PROJECT_ROOT / "data" / "runs" / joint_run_dir_name(label)
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{label}] skip (complete) -> {run_dir}", flush=True)
            return label
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(f"{run_dir} exists and is not complete.")
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{label}] pure_linear a={cell.a} {EB_TAG} tau={cell.tau_ps} "
          f"E0={cell.e0_eV} seed={cfg.seed} ({cell.role})", flush=True)
    print(f"[{label}] neutral ...", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{label}] ion ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    print(f"[{label}] relaxation (E2) ...", flush=True)
    relax = run_relaxation_stage(ion, cfg, save_path=run.root / "relaxation.npz")
    print(f"[{label}] detection ...", flush=True)
    detect = run_detection_stage(
        relax.checkpoint, cfg, save_path=run.root / "detection.npz")
    dm = detect.detected_mask
    nd = float(detect.n_detected[dm].mean()) if dm.any() else float("nan")
    print(f"[{label}] done -> {run_dir} (n_detect_mean={nd:.2f}, "
          f"droplet_retained={int(np.count_nonzero(~dm))}/{dm.size})",
          flush=True)
    return label


def _select(labels: list[str]) -> list[str]:
    if not labels:
        return [c.label for c in JOINT_MATRIX]
    unknown = sorted(set(labels) - set(_SPEC_BY_LABEL))
    if unknown:
        raise SystemExit(f"unknown cell(s) {unknown}; expected among "
                         f"{[c.label for c in JOINT_MATRIX]}")
    return labels


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    labels = _select(args.labels)
    print(f"Atlas (a, tau) joint ring (plan §6.6): {len(labels)} cell(s), "
          f"N={N}, seed={SEED} (the §6.3 ring seed — CRN with lr1-lr7, h405p "
          f"and the §6.5 clone s1), all {EB_TAG}, "
          f"concurrency={args.concurrency}"
          f"{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument runs -- standing point does not move).", flush=True)

    verify_joint_preregistration()

    if args.dry_run:
        for label in labels:
            cell = _SPEC_BY_LABEL[label]
            cfg, diff = verify_joint_cell(label)
            print(f"[{label}] cfg OK a={cell.a} tau={cell.tau_ps} "
                  f"E0={cell.e0_eV} seed={cfg.seed} -> "
                  f"{joint_run_dir_name(label)}", flush=True)
            print(f"[{label}]   diff vs h405p partner: {sorted(diff)}",
                  flush=True)
        print("dry run complete: LJ-P1 + LR-P1 + all cfg/CRN guards pass.",
              flush=True)
        return 0

    concurrency = max(1, args.concurrency)
    if concurrency == 1 or len(labels) == 1:
        for label in labels:
            _run_one(label)
        return 0

    from concurrent.futures import ProcessPoolExecutor, as_completed
    workers = min(concurrency, len(labels))
    print(f"launching {len(labels)} cells through a {workers}-slot pool ...",
          flush=True)
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as poolx:
        futures = {poolx.submit(_run_one, lb): lb for lb in labels}
        for fut in as_completed(futures):
            label = futures[fut]
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[{label}] FAILED: {exc!r}", flush=True)
    if failures:
        print(f"{failures}/{len(labels)} cell(s) failed.", flush=True)
        return 1
    print(f"all {len(labels)} cells complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
