"""Detector-stage fix — scorer for the new h405 reference battery (§5.4 D9).

Scores the ``gen_tier2_detfix_battery.py`` members (h405, Stage I to
t_h = 500 ps, E2 skipped, partner-aware Coulomb closure) against their
E2-era CRN twins, and ``r6`` against the t_h = 502.84 ps validation run.

Pure scorer — runs no MD and writes only the CSVs below. The observable
vector is the committed battery scorer's (``observable_columns`` +
``lowke_columns`` + ``md_gate`` + ``S``; rule 1). The E2-era twins are closed
**post hoc** by re-running the detection stage on their committed
``relaxation.npz`` with ``detection_coulomb_closure="partner_aware"``: the
stage-private RNG stream and, under ``co_moving``, a velocity-independent
jump chain make the event record identical to the committed one (asserted),
so this is exactly "E2-era + closure" with the package's single closure
implementation.

Order of operations (oracles first; any failure stops before a new number):

1. **Committed-row oracle** — each E2-era member rescored from its committed
   ``detection.npz`` equals its ``atlas_g4step2_battery.csv`` row on
   :data:`ORACLE_COLS` to 4 decimals (r6's twin: the finals-table h405 row).
2. **Pipeline-identity oracle** (§5.4 D5) — the new members' cfg is the skip
   path with the closure, the twins' the E2 path without it; never mixed.
3. **Post-hoc determinism oracle** — the re-run detection on each twin
   reproduces the committed ``n_detected`` / ``state_reason`` / event times
   bit-exactly.
4. **Closure anchor** — the r6 twin (``g4fh405``) closed post hoc reproduces
   the §3d "committed + closure" row, and the th503 validation run closed in
   the stage reproduces the §3d "new + closure" row (4 decimals).

Sharp checks per member (reported PASS/FAIL, §3d precedent):

* ``neutral.npz`` identical to the twin's;
* ``ion.npz`` bit-identical at every common stored step (twin: the 0–30 ps
  Stage I; r6 vs th503: all common steps ≤ 500 ps — same ``default_rng``
  stream throughout Stage I);
* excluded (retained) ion ids identical.

Outputs (committed; every column, incl. KE — memory "scorer CSVs carry full
vector"):

* ``data/runs/h2b_forward_model/detfix_battery_table.csv`` — one row per
  (member, pipeline): new, E2-era raw, E2-era closed; s1–s5 + pooled
  (N = 5000, the paired comparison); r6 rows incl. the th503 validation run
  closed; and ``pooled_new6000`` (s1–s5 + r6) = **the production reference**,
  cross-checked against the pooled figures container.
* ``data/runs/h2b_forward_model/detfix_battery_paired.csv`` — new minus
  E2-era-closed per member + pooled, plus the sharp-check outcomes.

Invocation::

    python scripts/post_processing/tier2_detfix_battery_table.py
"""

from __future__ import annotations

import csv
import dataclasses
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    pool_confirmation_reads,
    read_confirmation_detection,
)
from i2_helium_md.simulation.checkpoint import load_ion_checkpoint  # noqa: E402
from i2_helium_md.simulation.detection_stage import (  # noqa: E402
    load_detection_result,
    run_detection_stage,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.gen_detstagefix_validation import validation_run_dir_name  # noqa: E402
from scripts.gen_tier2_detfix_battery import (  # noqa: E402
    MEMBER_SEEDS,
    member_run_dir_name,
    twin_run_dir,
)
from scripts.post_processing.tier2atlas_g4step2_battery_table import (  # noqa: E402
    FINALS_TABLE_CSV,
    GATE_N1,
    GATE_NBAR,
    SAVE_CSV_PATH as E2_BATTERY_CSV,
    _s,
    lowke_columns,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    RUNS_ROOT,
    observable_columns,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

OUT_DIR = RUNS_ROOT / "h2b_forward_model"
TABLE_CSV = OUT_DIR / "detfix_battery_table.csv"
PAIRED_CSV = OUT_DIR / "detfix_battery_paired.csv"
#: The production reference figures container (build_pooled_detection_container.py
#: --target detfix_h405); cross-checked against the pooled_new6000 row if present.
PRODUCTION_CONTAINER = "9A_drag_shared_pure_cubic_N6000_detfix_conf270_h405pooled_th500"

BATTERY_MEMBERS = ("s1", "s2", "s3", "s4", "s5")
ORACLE_COLS = ("nbar_det", "n1_solv", "w1_solv", "midHot", "deepKE")
PAIRED_COLS = ("nbar_det", "n1_solv", "w1_solv", "midHot", "deepKE",
               "chi2_med", "KE1_mean", "KE2_mean", "S")

#: TIER2_DetectorStageFix.md §3d (committed numbers; 4-decimal anchors).
ANCHOR_TWIN_R6_CLOSED = {"midHot": 0.9569, "deepKE": 0.6039, "nbar_det": 3.9551}
ANCHOR_TH503_CLOSED = {"midHot": 0.9572, "deepKE": 0.6057, "nbar_det": 3.9584}

_ION_ARRAYS = ("positions_x", "positions_y", "positions_z", "velocities_x",
               "velocities_y", "velocities_z", "mass_history_kg", "E_kin_eV",
               "E_pot_eV", "E_dissip_eV", "E_mass_transfer_eV", "E_int_eV",
               "n_shell")

# ---------------------------------------------------------------------------


class OracleFailure(AssertionError):
    """An oracle failed — no new number may be read."""


def _pipeline(cfg) -> str:
    """The pipeline identity derived from cfg (§5.4 D5; no stamp field)."""
    stage = "E2" if cfg.relaxation_stage_enabled else "skip"
    return f"{stage}+closure:{cfg.detection_coulomb_closure}"


def _row(label, member, pipeline, seed, read, ab, ked) -> dict[str, Any]:
    row: dict[str, Any] = {"label": label, "member": member,
                           "pipeline": pipeline, "seed": seed}
    row.update(observable_columns(read, ab, ked))
    row.update(lowke_columns(read))
    row["md_gate"] = int(GATE_N1[0] <= float(row["n1_solv"]) <= GATE_N1[1]
                         and GATE_NBAR[0] <= float(row["nbar_det"]) <= GATE_NBAR[1])
    row["S"] = _s(row)
    return row


def _csv_row(path: Path, label: str) -> dict[str, str]:
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["label"] == label:
                return r
    raise OracleFailure(f"no row {label!r} in {path}")


def _check_4dp(tag: str, got: dict, want: dict, cols) -> None:
    bad = [f"{c}: {float(got[c]):.4f} vs {float(want[c]):.4f}" for c in cols
           if f"{float(got[c]):.4f}" != f"{float(want[c]):.4f}"]
    if bad:
        raise OracleFailure(f"{tag}: " + "; ".join(bad))


def closed_post_hoc(run_dir: Path, committed):
    """Re-run the detection stage with the closure on a finished run's seed
    checkpoint (``relaxation.npz`` if present, else ``ion.npz``); assert the
    event record is identical to the committed one (oracle 3)."""
    cfg = RunDirectory(run_dir).load_cfg()
    cfg = dataclasses.replace(cfg, detection_coulomb_closure="partner_aware")
    seed_file = ("relaxation.npz" if (run_dir / "relaxation.npz").exists()
                 else "ion.npz")
    det = run_detection_stage(load_ion_checkpoint(run_dir / seed_file), cfg)
    for field in ("n_detected", "state_reason", "event_time_ps", "event_offsets"):
        if not np.array_equal(getattr(det, field), getattr(committed, field)):
            raise OracleFailure(
                f"post-hoc determinism: {run_dir.name} {field} differs from the "
                "committed detection.npz")
    return det


def _same_npz(a: Path, b: Path) -> bool:
    with np.load(a, allow_pickle=False) as za, np.load(b, allow_pickle=False) as zb:
        if set(za.files) != set(zb.files):
            return False
        return all(np.array_equal(za[k], zb[k]) for k in za.files)


def common_step_identity(a: Path, b: Path, dt_ps: float) -> tuple[int, bool]:
    """Bit-identity of two ion.npz at every common stored step.

    Steps are matched by ``round(time_ps / dt)``. Returns (number of common
    columns, all identical)."""
    with np.load(a, allow_pickle=False) as za, np.load(b, allow_pickle=False) as zb:
        sa = np.rint(za["time_ps"] / dt_ps).astype(np.int64)
        sb = np.rint(zb["time_ps"] / dt_ps).astype(np.int64)
        common, ia, ib = np.intersect1d(sa, sb, return_indices=True)
        ok = all(np.array_equal(za[k][:, ia], zb[k][:, ib]) for k in _ION_ARRAYS)
    return int(common.size), bool(ok)


def excluded_ids(det) -> np.ndarray:
    return np.flatnonzero(~det.detected_mask)


def sharp_checks(new_dir: Path, ref_dir: Path, new_det, ref_det, dt_ps: float) -> dict:
    n_cols, ion_ok = common_step_identity(new_dir / "ion.npz", ref_dir / "ion.npz", dt_ps)
    return {
        "neutral_identical": _same_npz(new_dir / "neutral.npz", ref_dir / "neutral.npz"),
        "ion_common_cols": n_cols,
        "ion_common_identical": ion_ok,
        "excluded_ids_identical": bool(np.array_equal(excluded_ids(new_det),
                                                      excluded_ids(ref_det))),
        "n_excluded_new": int(excluded_ids(new_det).size),
        "n_excluded_ref": int(excluded_ids(ref_det).size),
    }


def ke_agreement(a, b) -> tuple[int, float, float]:
    """Per-ion |dKE|/KE for ions scored in both with the same detected n.
    Returns (count, median, p99)."""
    both = a.detected_mask & b.detected_mask & (a.n_detected == b.n_detected)
    rel = (np.abs(a.E_kin_detected_eV[both] - b.E_kin_detected_eV[both])
           / b.E_kin_detected_eV[both])
    return int(both.sum()), float(np.median(rel)), float(np.quantile(rel, 0.99))


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    missing = [m for m in MEMBER_SEEDS
               if not (RUNS_ROOT / member_run_dir_name(m) / "detection.npz").exists()]
    if missing:
        print(f"not yet run: {missing} — stopping before any number.")
        return 1

    rows: list[dict[str, Any]] = []
    paired: list[dict[str, Any]] = []
    reads = {"new": [], "e2raw": [], "e2closed": []}

    for member in MEMBER_SEEDS:
        new_dir = RUNS_ROOT / member_run_dir_name(member)
        twin_dir = twin_run_dir(member)
        cfg_new = RunDirectory(new_dir).load_cfg()
        cfg_twin = RunDirectory(twin_dir).load_cfg()
        # Oracle 2: pipeline identity, never mixed.
        if (_pipeline(cfg_new) != "skip+closure:partner_aware"
                or _pipeline(cfg_twin) != "E2+closure:none"):
            raise OracleFailure(f"[{member}] pipelines {_pipeline(cfg_new)} / "
                                f"{_pipeline(cfg_twin)}")
        det_new = load_detection_result(new_dir / "detection.npz")
        det_twin = load_detection_result(twin_dir / "detection.npz")
        read_twin = read_confirmation_detection(det_twin, label=f"{member}_E2")
        raw_row = _row(f"{member}_E2raw", member, _pipeline(cfg_twin),
                       cfg_twin.seed, read_twin, ab, ked)
        # Oracle 1: committed rows.
        if member == "r6":
            _check_4dp("[r6] committed finals h405 row", raw_row,
                       _csv_row(FINALS_TABLE_CSV, "h405"), ORACLE_COLS)
        else:
            _check_4dp(f"[{member}] committed battery row", raw_row,
                       _csv_row(E2_BATTERY_CSV, member), ORACLE_COLS)
        # Oracle 3 (inside closed_post_hoc).
        det_closed = closed_post_hoc(twin_dir, det_twin)
        read_closed = read_confirmation_detection(det_closed, label=f"{member}_E2c")
        closed_row = _row(f"{member}_E2closed", member,
                          "E2+closure:partner_aware(post-hoc)", cfg_twin.seed,
                          read_closed, ab, ked)
        read_new = read_confirmation_detection(det_new, label=member)
        new_row = _row(f"{member}_new", member, _pipeline(cfg_new), cfg_new.seed,
                       read_new, ab, ked)
        rows += [new_row, raw_row, closed_row]

        sharp = sharp_checks(new_dir, twin_dir, det_new, det_twin, cfg_new.dt_ion)
        n_ke, ke_med, ke_p99 = ke_agreement(det_new, det_closed)
        p = {"member": member, "vs": "E2closed"}
        p.update({f"d_{c}": float(new_row[c]) - float(closed_row[c]) for c in PAIRED_COLS})
        p.update(sharp)
        p.update({"ke_same_n_count": n_ke, "ke_rel_med": ke_med, "ke_rel_p99": ke_p99})
        paired.append(p)

        if member == "r6":
            # Oracle 4a: the §3d committed + closure row.
            _check_4dp("[r6] §3d committed+closure anchor", closed_row,
                       ANCHOR_TWIN_R6_CLOSED, ANCHOR_TWIN_R6_CLOSED)
            th_dir = RUNS_ROOT / validation_run_dir_name("h405")
            det_th = closed_post_hoc(th_dir, load_detection_result(th_dir / "detection.npz"))
            th_row = _row("r6_th503closed", member,
                          "skip(502.84ps)+closure:partner_aware(post-hoc)",
                          cfg_new.seed, read_confirmation_detection(det_th, label="th503"),
                          ab, ked)
            # Oracle 4b: the §3d new + closure row.
            _check_4dp("[r6] §3d th503 new+closure anchor", th_row,
                       ANCHOR_TH503_CLOSED, ANCHOR_TH503_CLOSED)
            rows.append(th_row)
            sharp_th = sharp_checks(new_dir, th_dir, det_new, det_th, cfg_new.dt_ion)
            n_ke, ke_med, ke_p99 = ke_agreement(det_new, det_th)
            p = {"member": member, "vs": "th503closed"}
            p.update({f"d_{c}": float(new_row[c]) - float(th_row[c]) for c in PAIRED_COLS})
            p.update(sharp_th)
            p.update({"ke_same_n_count": n_ke, "ke_rel_med": ke_med, "ke_rel_p99": ke_p99})
            paired.append(p)
            read_r6_new = read_new
        else:
            reads["new"].append(read_new)
            reads["e2raw"].append(read_twin)
            reads["e2closed"].append(read_closed)
        print(f"[{member}] scored; sharp: {sharp}", flush=True)

    pooled = {}
    for key, pipe in (("new", "skip+closure:partner_aware"),
                      ("e2raw", "E2+closure:none"),
                      ("e2closed", "E2+closure:partner_aware(post-hoc)")):
        pooled[key] = _row(f"pooled_{key}", "pooled", pipe, "-",
                           pool_confirmation_reads(reads[key], label=f"pooled_{key}"),
                           ab, ked)
        rows.append(pooled[key])
    # Pooled E2raw must equal the committed pooled battery row.
    _check_4dp("committed pooled battery row", pooled["e2raw"],
               _csv_row(E2_BATTERY_CSV, "pooled"), ORACLE_COLS)
    # The production reference: all six members (s1-s5 + r6), N = 6000.
    # Oracle: it equals the read of the pooled figures container (in-memory
    # pooling == the pair-preserving concatenated detection.npz).
    prod = _row("pooled_new6000", "pooled6000", "skip+closure:partner_aware",
                "-", pool_confirmation_reads(reads["new"] + [read_r6_new],
                                             label="pooled_new6000"), ab, ked)
    container = RUNS_ROOT / PRODUCTION_CONTAINER / "detection.npz"
    if container.exists():
        _check_4dp("N=6000 container", prod, _row(
            "container", "-", "-", "-", read_confirmation_detection(
                load_detection_result(container), label="container"), ab, ked),
            ORACLE_COLS + ("chi2_med", "S"))
    rows.append(prod)
    p = {"member": "pooled", "vs": "E2closed"}
    p.update({f"d_{c}": float(pooled["new"][c]) - float(pooled["e2closed"][c])
              for c in PAIRED_COLS})
    paired.append(p)
    keys = list(paired[0].keys())
    paired = [{k: r.get(k, "") for k in keys} for r in paired]

    show = ("label", "pipeline", "num_scored", "trap", "trap_bound", "trap_marg",
            "nbar_det", "n1_solv", "w1_solv", "midHot", "deepKE", "chi2_med",
            "KE1_mean", "KE2_mean", "S")
    print("\n=== detfix battery — observable vector ===")
    print(format_table([{c: r[c] for c in show} for r in rows]))
    print("\n=== paired: new minus reference ===")
    print(format_table(paired))
    print(f"\nCSV -> {write_rows_csv(TABLE_CSV, rows)}")
    print(f"CSV -> {write_rows_csv(PAIRED_CSV, paired)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
