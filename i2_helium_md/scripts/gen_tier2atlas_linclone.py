"""Atlas free-form linear — the h405-clone MD battery (plan §6.5).

Runs the `TIER2_FREEFORM_LINEAR_TWIN_SWEEP_PLAN` §6.5 battery: **one cell,
three seeds, N = 500** at the G2-adopted corrected geometry. The cell is
the §6.4 item-2 "h405 clone" — the committed `atlas_linsweep.csv` row that
reproduces the h405 landing vector from a **kink-free, uncapped** pure
linear law:

* `pure_linear a = 42.5 amu/ps / eb0482 / τ 6.4 ps / E₀ 0.31 eV` (Φ 0.985)

Seeds (plan §6.5): **20260731** (the §6.3 ring seed — the clone thereby
slots into the committed ring table as a same-N, same-seed 8th cell *and*
CRN-pairs for free against the committed `h405p` partner) plus fresh
**20260812 / 20260813**. No h405 partner is re-run: the committed g4
Step-2 battery measured h405's per-seed SD at 0.0024 (KE₁) / 0.0074
(midHot) / 0.0102 (n₁), so the Δ target is a constant (plan §6.5).

The §6.3 ring covered a ∈ [27.5, 35] at τ 4.8 — this cell sits **outside
the measured twin↔MD bias box in both swept axes**, which is why MD is
required and why the CL-P6 box extension is a product regardless of the
verdict.

**Oracles (enforced before any MD and by --dry-run):**

1. **LC-P1** — the frozen clone twin row reproduces string-exact from the
   committed `atlas_linsweep.csv`, and the ring's own LR-P1 (7 lin rows +
   the h405 authority anchor) re-runs unchanged.
2. **CRN guard — the pairing *is* the guard.** The seed-20260731 cell's
   cfg diff vs the committed `h405p` cfg must contain neither `seed` nor
   `num_molecules` (their identity is what makes the pair common-random-
   number), and must otherwise be exactly the pre-registered key set. The
   fresh-seed cells must diff against the 20260731 cell in exactly
   `{"seed"}` (the g4 Step-2 precedent).
3. The standing-battery cfg-diff guard, reused verbatim from the ring.

Cell construction, geometry/provenance guards and the free-form posture
(`extraction_method="free_form"`, no §6.5.1 binding stamp, documented
`allow_unvalidated_binding_pairing` hatch) are **imported** from
`gen_tier2atlas_linring.py` with only the seed replaced — there is no
second copy of the cell construction (rule 1).

Atlas stance: instrument runs — nothing here moves `finc1v725`/h405.

Usage::

    python scripts/gen_tier2atlas_linclone.py --dry-run   # oracles + guards
    python scripts/gen_tier2atlas_linclone.py             # all 3 seeds
    python scripts/gen_tier2atlas_linclone.py s1          # selected seeds
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
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

from scripts.gen_tier2atlas_linring import (  # noqa: E402
    LINSWEEP_CSV,
    LinRingCell,
    N,
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
# USER SETTINGS
# =============================================================================

# The §6.4 item-2 clone cell. N, geometry, ladder and every other pin come
# from the ring module — only (a, E_bind, tau, E0) are this cell's.
CLONE_CELL = LinRingCell(
    "clone", 42.5, None, "eb0482", 6.4, 0.31,
    "h405 clone (plan §6.4 item 2), Phi 0.985",
)

# Three seeds (plan §6.5). s1 is the §6.3 ring seed: same N, same seed =>
# CRN pair against the committed h405p partner at zero extra MD.
MEMBER_SEEDS: dict[str, int] = {
    "s1": 20260731,
    "s2": 20260812,
    "s3": 20260813,
}
CRN_MEMBER = "s1"

# The committed §6.3 ring partner the CRN member pairs against.
PARTNER_RUN = "9A_drag_shared_pure_cubic_N500_tier2atlas_conf270_linrh405p"

# The pre-registered cfg-diff key set vs that partner (form + chord +
# well + the §6.5.1 hatch + this cell's τ/E₀). `seed` and `num_molecules`
# are FORBIDDEN here — their identity is the CRN pairing.
PARTNER_DIFF_KEYS = frozenset({
    "drag_form",
    "drag_coefficients",
    "binding_energy_I_ion_eV",
    "allow_unvalidated_binding_pairing",
    "internal_energy_cooling_tau_ps",
    "internal_energy_partition_fraction",
})
CRN_FORBIDDEN_DIFF_KEYS = frozenset({"seed", "num_molecules"})

# LC-P1: the frozen clone twin row (committed atlas_linsweep.csv strings,
# compared EXACTLY). Columns below.
CLONE_TWIN_COLS = ("trapped_frac", "suppressed_frac", "nbar_det", "n1_solv",
                   "w1_solv", "n1_ke_eV", "ke2_eV", "deepke", "midhot_arith",
                   "phi", "gate")
CLONE_TWIN_ROW = ("0.0476", "0.0456", "4.576", "0.1966", "0.6857", "0.6153",
                  "0.5155", "0.6444", "1.0144", "0.985", "1")

DEFAULT_CONCURRENCY = 3
SKIP_COMPLETED_RUNS = True
OVERWRITE_EXISTING_RUN = True

_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "cfg.json", "neutral.npz", "ion.npz", "relaxation.npz", "detection.npz",
)


# =============================================================================
# Oracles and guards
# =============================================================================


def verify_clone_preregistration() -> None:
    """LC-P1: the frozen clone twin row reproduces string-exact from the
    committed ``atlas_linsweep.csv``, and the ring's LR-P1 still passes.

    Fails loud — a mismatch means the frozen row or a committed CSV drifted,
    and the battery must not launch (gate-on-committed-artifacts rule).
    """
    verify_twin_preregistration()
    key = ("lin", f"{CLONE_CELL.a:.1f}", CLONE_CELL.eb_tag,
           str(CLONE_CELL.tau_ps), str(CLONE_CELL.e0_eV))
    with open(PROJECT_ROOT / LINSWEEP_CSV, newline="", encoding="utf-8") as fh:
        row = next(
            (r for r in csv.DictReader(fh)
             if (r["arm"], r["a"], r["E_bind_tag"], r["tau_ps"], r["E0_eV"])
             == key),
            None,
        )
    if row is None:
        raise AssertionError(
            f"LC-P1 FAILED: no committed linsweep row at {key}."
        )
    got = tuple(row[col] for col in CLONE_TWIN_COLS)
    if got != CLONE_TWIN_ROW:
        raise AssertionError(
            f"LC-P1 FAILED: committed clone row {got} != frozen "
            f"pre-registration {CLONE_TWIN_ROW}."
        )
    print("LC-P1 clone oracle PASSED: the frozen h405-clone twin row "
          "reproduces string-exact from the committed atlas_linsweep.csv.",
          flush=True)


def clone_run_dir_name(member: str) -> str:
    """Run-dir basename in the atlas namespace (no ``_tier2_`` substring)."""
    name = run_dir_name(
        "9A", "shared_pure_cubic", N,
        run_tag=f"tier2atlas_conf270_linclone{member}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_member(member: str) -> SimConfig:
    """The clone cell with only the seed replaced (rule 1: the cell itself is
    built by :func:`gen_tier2atlas_linring.build_cell`)."""
    if member not in MEMBER_SEEDS:
        raise KeyError(f"unknown member {member!r}; expected among "
                       f"{sorted(MEMBER_SEEDS)}")
    cfg = dataclasses.replace(build_cell(CLONE_CELL), seed=MEMBER_SEEDS[member])
    cfg.validate()
    return cfg


def cfg_field_diff(cfg_a: SimConfig, cfg_b: SimConfig) -> list[str]:
    """Sorted field names where two built configs differ (in-memory twin of
    :func:`cfg_diff_vs_reference`, which needs a stored ``cfg.json``)."""
    return sorted(
        f.name for f in dataclasses.fields(SimConfig)
        if getattr(cfg_a, f.name) != getattr(cfg_b, f.name)
    )


def verify_crn_pairing(cfg: SimConfig, member: str) -> list[str]:
    """The §6.5 CRN guard.

    For the CRN member the cfg diff vs the committed `h405p` partner must be
    exactly :data:`PARTNER_DIFF_KEYS`; for the fresh-seed members the same set
    plus ``seed``. ``num_molecules`` may never appear — N identity is half of
    what makes the pair common-random-number, and ``seed`` identity is the
    other half (so it is forbidden on the CRN member specifically).
    """
    ref_path = PROJECT_ROOT / "data" / "runs" / PARTNER_RUN / "cfg.json"
    diff = cfg_diff_vs_reference(cfg, ref_path, context=member)
    is_crn = MEMBER_SEEDS[member] == MEMBER_SEEDS[CRN_MEMBER]
    forbidden = CRN_FORBIDDEN_DIFF_KEYS if is_crn else frozenset({"num_molecules"})
    hit = sorted(forbidden & set(diff))
    if hit:
        raise AssertionError(
            f"[{member}] CRN guard FAILED: {hit} differ(s) vs {PARTNER_RUN} — "
            "the CRN pairing is exactly the identity of those keys."
        )
    expected = PARTNER_DIFF_KEYS if is_crn else PARTNER_DIFF_KEYS | {"seed"}
    if set(diff) != expected:
        raise AssertionError(
            f"[{member}] cfg diff vs the partner is {sorted(diff)}; "
            f"pre-registered exactly {sorted(expected)}."
        )
    return diff


def verify_seed_only_vs_crn_member(cfg: SimConfig, member: str) -> None:
    """Fresh-seed members must differ from the CRN member in exactly the seed
    (the g4 Step-2 precedent: same cell, fresh seed, no drifted pin)."""
    if member == CRN_MEMBER:
        return
    diff = cfg_field_diff(cfg, build_member(CRN_MEMBER))
    if diff != ["seed"]:
        raise AssertionError(
            f"[{member}] differs from {CRN_MEMBER} in {diff}; pre-registered "
            "exactly ['seed']."
        )


def verify_member(member: str) -> tuple[SimConfig, list[str]]:
    """All four guards for one member; returns the cfg and its partner diff."""
    cfg = build_member(member)
    verify_corrected_geometry(cfg, CLONE_CELL)
    verify_against_reference(cfg, CLONE_CELL)   # vs the standing battery
    partner_diff = verify_crn_pairing(cfg, member)
    verify_seed_only_vs_crn_member(cfg, member)
    return cfg, partner_diff


# =============================================================================
# Runner
# =============================================================================


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _run_one(member: str) -> str:
    cfg, _ = verify_member(member)
    run_dir = PROJECT_ROOT / "data" / "runs" / clone_run_dir_name(member)
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{member}] skip (complete) -> {run_dir}", flush=True)
            return member
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(f"{run_dir} exists and is not complete.")
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{member}] pure_linear a={CLONE_CELL.a} {CLONE_CELL.eb_tag} "
          f"tau={CLONE_CELL.tau_ps} E0={CLONE_CELL.e0_eV} seed={cfg.seed} "
          f"({CLONE_CELL.role})", flush=True)
    print(f"[{member}] neutral ...", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{member}] ion ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    print(f"[{member}] relaxation (E2) ...", flush=True)
    relax = run_relaxation_stage(ion, cfg, save_path=run.root / "relaxation.npz")
    print(f"[{member}] detection ...", flush=True)
    detect = run_detection_stage(
        relax.checkpoint, cfg, save_path=run.root / "detection.npz")
    dm = detect.detected_mask
    nd = float(detect.n_detected[dm].mean()) if dm.any() else float("nan")
    print(f"[{member}] done -> {run_dir} (n_detect_mean={nd:.2f}, "
          f"droplet_retained={int(np.count_nonzero(~dm))}/{dm.size})",
          flush=True)
    return member


def _select(labels: list[str]) -> list[str]:
    if not labels:
        return list(MEMBER_SEEDS)
    unknown = sorted(set(labels) - set(MEMBER_SEEDS))
    if unknown:
        raise SystemExit(f"unknown member(s) {unknown}; expected among "
                         f"{sorted(MEMBER_SEEDS)}")
    return labels


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    members = _select(args.labels)
    print(f"Atlas h405-clone MD battery (plan §6.5): {len(members)} seed(s) of "
          f"pure_linear a={CLONE_CELL.a} {CLONE_CELL.eb_tag} "
          f"tau={CLONE_CELL.tau_ps} E0={CLONE_CELL.e0_eV}, N={N}, "
          f"seeds={[MEMBER_SEEDS[m] for m in members]} "
          f"(s1 = the §6.3 ring seed => CRN vs the committed h405p), "
          f"concurrency={args.concurrency}"
          f"{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument runs -- standing point does not move).", flush=True)

    verify_clone_preregistration()

    if args.dry_run:
        for member in members:
            cfg, partner_diff = verify_member(member)
            print(f"[{member}] cfg OK seed={cfg.seed} -> "
                  f"{clone_run_dir_name(member)}", flush=True)
            print(f"[{member}]   diff vs h405p partner: {sorted(partner_diff)}",
                  flush=True)
        print("dry run complete: LC-P1 + LR-P1 + all cfg/CRN guards pass.",
              flush=True)
        return 0

    concurrency = max(1, args.concurrency)
    if concurrency == 1 or len(members) == 1:
        for member in members:
            _run_one(member)
        return 0

    from concurrent.futures import ProcessPoolExecutor, as_completed
    workers = min(concurrency, len(members))
    print(f"launching {len(members)} members through a {workers}-slot pool ...",
          flush=True)
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as poolx:
        futures = {poolx.submit(_run_one, m): m for m in members}
        for fut in as_completed(futures):
            member = futures[fut]
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[{member}] FAILED: {exc!r}", flush=True)
    if failures:
        print(f"{failures}/{len(members)} member(s) failed.", flush=True)
        return 1
    print(f"all {len(members)} members complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
