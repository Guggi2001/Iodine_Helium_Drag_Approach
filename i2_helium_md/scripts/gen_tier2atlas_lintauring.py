"""Atlas free-form linear — the short-τ decisive pair (plan §6.8).

Runs the two cells that adjudicate whether the §6.7 τ refinement found a
real joint winner, at **N = 500, seed 20260731** — the §6.3 ring seed, so
both stay CRN-paired ion-for-ion with the whole 14-cell linear family
(lr1–lr7, h405p, the §6.5 clone, the §6.6 joint ring).

* `t1`  a 35.0 / eb0482 / τ 4.0 / E₀ 0.45 — the best predicted **joint**
  winner (predictor: W₁ 0.532, KE₁ 0.659; cost midHot 1.380)
* `t2`  a 42.5 / eb0482 / τ 4.4 / E₀ 0.44 — the **in-band** W₁ optimum
  (predictor: W₁ 0.490 at midHot 0.923; KE₁ 0.506, expected below h405)

Both sit at τ values the committed sweep never sampled, which is why the
§6.7 refinement had to run first. The two available instruments
**disagree by ~0.5** on these cells — twin W₁ says 1.01/1.06, the (n₁,
tail) predictor says 0.53/0.49 — and this ring adjudicates between them
(plan §6.8 TR-P3).

**Oracles (enforced before any MD and by --dry-run):**

1. **TR-P0** — both frozen twin rows reproduce **string-exact** from the
   committed `atlas_lintau.csv` (the §6.7 artifact), plus the §6.3 LR-P1.
2. **CRN guard** — the seed is shared with the committed `h405p`, so
   `seed` and `num_molecules` must both be absent from every cell's diff
   against it (the §6.6 posture).
3. The standing-battery cfg-diff guard, reused from the ring.

Cell construction and the corrected-geometry / free-form provenance
guards are imported from `gen_tier2atlas_linring.py` (rule 1); the
detection-only resume path is imported from `gen_tier2atlas_linjoint.py`.

Atlas stance: instrument runs — nothing here moves `finc1v725`/h405.

Usage::

    python scripts/gen_tier2atlas_lintauring.py --dry-run
    python scripts/gen_tier2atlas_lintauring.py            # both cells
    python scripts/gen_tier2atlas_lintauring.py t1
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
from scripts.gen_tier2atlas_linjoint import (  # noqa: E402
    _PRE_DETECTION_ARTIFACTS,
    _run_detection_only,
)
from scripts.gen_tier2atlas_linring import (  # noqa: E402
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

import zipfile  # noqa: E402


# =============================================================================
# USER SETTINGS — the frozen §6.8 pair
# =============================================================================

class TauCell(NamedTuple):
    label: str
    a: float
    tau_ps: float
    e0_eV: float
    role: str


TAU_MATRIX: tuple[TauCell, ...] = (
    TauCell("t1", 35.0, 4.0, 0.45, "best predicted joint winner"),
    TauCell("t2", 42.5, 4.4, 0.44, "in-band W1 optimum"),
)
_SPEC_BY_LABEL = {c.label: c for c in TAU_MATRIX}

EB_TAG = "eb0482"

# The committed h405p partner's own knobs — needed because `t2` coincides
# with its cooling clock (see `expected_partner_diff`).
PARTNER_TAU_PS = 4.4
PARTNER_E0_EV = 0.405

# The §6.7 twin artifact this ring's cells are frozen against.
LINTAU_CSV = Path("data/runs/h2b_forward_model/atlas_lintau.csv")

TAU_TWIN_COLS = ("trapped_frac", "suppressed_frac", "nbar_det", "n1_solv",
                 "w1_solv", "n1_ke_eV", "ke2_eV", "deepke", "midhot_arith",
                 "twin_tail_ge10", "twin_tail_ge13", "phi", "gate")
TAU_TWIN_ROWS: dict[str, tuple[str, ...]] = {
    "t1": ("0.0128", "0.2349", "4.509", "0.1932", "1.0054", "0.6448",
           "0.5923", "2.2654", "1.5236", "0.21945", "0.11727", "0.811", "1"),
    "t2": ("0.0476", "0.255", "4.418", "0.1943", "1.0615", "0.4914",
           "0.4358", "1.1279", "1.019", "0.225", "0.12032", "0.985", "1"),
}

# Plan §6.8 predicted MD values (the (n1, tail) predictor + the measured
# transfers) — frozen here so the comparison is pre-registered, not fitted
# after the fact. LOO RMSE on the W1 prediction is 0.081.
TAU_PREDICTED: dict[str, dict[str, float]] = {
    "t1": {"w1": 0.532, "KE1": 0.659, "tail": 0.1507, "midHot": 1.380,
           "n1": 0.1940, "nbar": 3.919},
    "t2": {"w1": 0.490, "KE1": 0.506, "tail": 0.1563, "midHot": 0.923,
           "n1": 0.1951, "nbar": 3.828},
}
PREDICTOR_LOO_RMSE = 0.081

DEFAULT_CONCURRENCY = 2          # measured ceiling on this machine
SKIP_COMPLETED_RUNS = True
OVERWRITE_EXISTING_RUN = True

_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "cfg.json", "neutral.npz", "ion.npz", "relaxation.npz", "detection.npz",
)


# =============================================================================
# Oracles and guards
# =============================================================================


def as_ring_cell(cell: TauCell) -> LinRingCell:
    """This ring's cell in the §6.3 constructor's vocabulary (rule 1)."""
    return LinRingCell(cell.label, cell.a, None, EB_TAG, cell.tau_ps,
                       cell.e0_eV, cell.role)


def verify_tau_preregistration() -> None:
    """TR-P0: both frozen twin rows reproduce string-exact from the committed
    §6.7 artifact, and the §6.3 LR-P1 still passes."""
    verify_twin_preregistration()
    with open(PROJECT_ROOT / LINTAU_CSV, newline="", encoding="utf-8") as fh:
        idx = {
            (r["a"], r["E_bind_tag"], r["tau_ps"], r["E0_eV"]): r
            for r in csv.DictReader(fh) if r["arm"] == "lin"
        }
    for cell in TAU_MATRIX:
        key = (f"{cell.a:.1f}", EB_TAG, str(cell.tau_ps), str(cell.e0_eV))
        row = idx.get(key)
        if row is None:
            raise AssertionError(
                f"[{cell.label}] TR-P0 FAILED: no committed lintau row at "
                f"{key}.")
        got = tuple(row[col] for col in TAU_TWIN_COLS)
        if got != TAU_TWIN_ROWS[cell.label]:
            raise AssertionError(
                f"[{cell.label}] TR-P0 FAILED: committed lintau row {got} != "
                f"frozen pre-registration {TAU_TWIN_ROWS[cell.label]}.")
    print(f"TR-P0 oracle PASSED: both frozen §6.8 twin rows reproduce "
          "string-exact from the committed atlas_lintau.csv.", flush=True)


def tau_run_dir_name(label: str) -> str:
    """Run-dir basename in the atlas namespace."""
    name = run_dir_name("9A", "shared_pure_cubic", N,
                        run_tag=f"tier2atlas_conf270_lint{label}")
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_tau_cell(label: str) -> SimConfig:
    """One cell's validated config at the shared whole-family CRN seed."""
    if label not in _SPEC_BY_LABEL:
        raise KeyError(f"unknown cell {label!r}; expected among "
                       f"{sorted(_SPEC_BY_LABEL)}")
    cfg = build_cell(as_ring_cell(_SPEC_BY_LABEL[label]))
    if cfg.seed != SEED:
        raise AssertionError(
            f"[{label}] seed {cfg.seed} != the shared CRN seed {SEED}.")
    cfg.validate()
    return cfg


def expected_partner_diff(cell: TauCell) -> frozenset[str]:
    """The pre-registered diff key set for this cell vs `h405p`.

    `t2` sits at τ 4.4, which is **exactly the partner's own cooling
    clock**, so `internal_energy_cooling_tau_ps` legitimately does not
    differ for it. The guard computes the expected set per cell rather
    than being loosened — a coincidence in a knob value must not become a
    licence for that knob to drift unnoticed on the other cell.
    """
    keys = set(PARTNER_DIFF_KEYS)
    if cell.tau_ps == PARTNER_TAU_PS:
        keys.discard("internal_energy_cooling_tau_ps")
    if cell.e0_eV == PARTNER_E0_EV:
        keys.discard("internal_energy_partition_fraction")
    return frozenset(keys)


def verify_crn_pairing(cfg: SimConfig, label: str) -> list[str]:
    """§6.6 posture: the seed is shared with the committed partner, so
    neither `seed` nor `num_molecules` may appear in the diff."""
    ref_path = PROJECT_ROOT / "data" / "runs" / PARTNER_RUN / "cfg.json"
    diff = cfg_diff_vs_reference(cfg, ref_path, context=label)
    hit = sorted(CRN_FORBIDDEN_DIFF_KEYS & set(diff))
    if hit:
        raise AssertionError(
            f"[{label}] CRN guard FAILED: {hit} differ(s) vs {PARTNER_RUN} — "
            "the whole-family pairing is exactly that identity.")
    expected = expected_partner_diff(_SPEC_BY_LABEL[label])
    if set(diff) != set(expected):
        raise AssertionError(
            f"[{label}] cfg diff vs the partner is {sorted(diff)}; "
            f"pre-registered exactly {sorted(expected)}.")
    return diff


def verify_tau_cell(label: str) -> tuple[SimConfig, list[str]]:
    cell = as_ring_cell(_SPEC_BY_LABEL[label])
    cfg = build_tau_cell(label)
    verify_corrected_geometry(cfg, cell)
    verify_against_reference(cfg, cell)
    return cfg, verify_crn_pairing(cfg, label)


# =============================================================================
# Runner
# =============================================================================


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _run_one(label: str) -> str:
    cell = _SPEC_BY_LABEL[label]
    cfg, _ = verify_tau_cell(label)
    run_dir = PROJECT_ROOT / "data" / "runs" / tau_run_dir_name(label)
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{label}] skip (complete) -> {run_dir}", flush=True)
            return label
        if all((run_dir / n_).exists() for n_ in _PRE_DETECTION_ARTIFACTS):
            try:
                return _run_detection_only(label, cfg, run_dir)
            except (zipfile.BadZipFile, EOFError, OSError, ValueError) as exc:
                print(f"[{label}] stored relaxation.npz unusable ({exc!r}) — "
                      "recomputing from scratch.", flush=True)
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
        return [c.label for c in TAU_MATRIX]
    unknown = sorted(set(labels) - set(_SPEC_BY_LABEL))
    if unknown:
        raise SystemExit(f"unknown cell(s) {unknown}")
    return labels


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    labels = _select(args.labels)
    print(f"Atlas short-τ decisive pair (plan §6.8): {len(labels)} cell(s), "
          f"N={N}, seed={SEED} (whole-family CRN), {EB_TAG}, "
          f"concurrency={args.concurrency}"
          f"{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument runs -- standing point does not move).", flush=True)

    verify_tau_preregistration()

    if args.dry_run:
        for label in labels:
            cell = _SPEC_BY_LABEL[label]
            cfg, diff = verify_tau_cell(label)
            p = TAU_PREDICTED[label]
            print(f"[{label}] cfg OK a={cell.a} tau={cell.tau_ps} "
                  f"E0={cell.e0_eV} -> {tau_run_dir_name(label)}", flush=True)
            print(f"[{label}]   diff vs h405p: {sorted(diff)}", flush=True)
            print(f"[{label}]   PRE-REGISTERED prediction: W1 {p['w1']:.3f} "
                  f"(+-{PREDICTOR_LOO_RMSE} LOO), KE1 {p['KE1']:.3f}, "
                  f"midHot {p['midHot']:.3f}", flush=True)
        print("dry run complete: TR-P0 + LR-P1 + all cfg/CRN guards pass.",
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
