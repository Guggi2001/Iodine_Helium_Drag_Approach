"""G4.13 scorer — c = 2.0 ridge-law MD battery vs the production reference.

THESIS_FINAL_GOALS.md G4.13 (pre-registered 4e0908c). Pure scorer: runs no MD.
Same observable vector as the production battery scorer (``_row`` from
``tier2_detfix_battery_table.py``: observables + low-n KE + MD gate + S).

Order (oracles first; any failure stops before a new number):

1. **Anchor** — the h405 detfix ``s1`` run dir rescored equals the committed
   ``detfix_battery_table.csv`` row ``s1_new`` (4 dp, :data:`ANCHOR_COLS`);
   ``pooled_new6000`` is read from the same committed CSV.
2. **Stage A** — every grid cell + ``h405g35`` at seed s1; the frozen Stage-B
   selection rule (:func:`select_stage_b`) applied to the gated cells.
3. **Stage B** (when its five members exist) — the selected cell pooled over
   s1-s5 + r6 (N = 6000) and the frozen decision rule (:func:`decide`).

Outputs (committed; every column): ``data/runs/h2b_forward_model/
g413_table.csv`` (one row per run + pooled) and ``g413_verdict.csv``.

    python scripts/post_processing/tier2_g413_table.py
"""

from __future__ import annotations

import csv
import sys
import warnings
from pathlib import Path
from typing import Any, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import load_he_abundance_reference  # noqa: E402
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    pool_confirmation_reads,
    read_confirmation_detection,
)
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.gen_tier2_detfix_battery import member_run_dir_name  # noqa: E402
from scripts.gen_tier2_g413_ridge_battery import (  # noqa: E402
    CELLS, STAGE_A_MEMBER, STAGE_B_MEMBERS, run_dir_for,
)
from scripts.post_processing.tier2_detfix_battery_table import (  # noqa: E402
    TABLE_CSV as DETFIX_TABLE_CSV,
    _check_4dp,
    _csv_row,
    _pipeline,
    _row,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    RUNS_ROOT,
)

OUT_DIR = RUNS_ROOT / "h2b_forward_model"
TABLE_CSV = OUT_DIR / "g413_table.csv"
VERDICT_CSV = OUT_DIR / "g413_verdict.csv"
ANCHOR_COLS = ("nbar_det", "n1_solv", "w1_solv", "midHot", "deepKE",
               "KE1_mean", "chi2_med")

# Pre-registered (G4.13; FIXED).
MIDHOT_BAND = (0.85, 1.15)
KE1_H405_REF = 0.640


def select_stage_b(rows: list[dict[str, Any]]) -> Optional[str]:
    """Frozen Stage-B rule over Stage-A ridge rows: among gated cells, midHot in
    band with max KE1_mean; none in band -> midHot nearest the band; ties ->
    lower W1. ``None`` if no cell gates."""
    gated = [r for r in rows if int(r["md_gate"]) == 1]
    if not gated:
        return None
    lo, hi = MIDHOT_BAND

    def dist(r):
        m = float(r["midHot"])
        return max(lo - m, m - hi, 0.0)

    in_band = [r for r in gated if dist(r) == 0.0]
    if in_band:
        best = min(in_band, key=lambda r: (-float(r["KE1_mean"]), float(r["w1_solv"])))
    else:
        best = min(gated, key=lambda r: (dist(r), float(r["w1_solv"])))
    return str(best["cell"])


def decide(pooled: dict[str, Any]) -> str:
    """Frozen decision rule on the Stage-B pooled row."""
    ok = (int(pooled["md_gate"]) == 1
          and float(pooled["KE1_mean"]) >= KE1_H405_REF
          and MIDHOT_BAND[0] <= float(pooled["midHot"]) <= MIDHOT_BAND[1])
    return "matches_or_beats_h405" if ok else "dominated"


def _write(rows: list[dict[str, Any]], path: Path) -> None:
    """CSV over the union of keys (rows differ: ref / Stage A / pooled)."""
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {path}")


def _read(run_dir: Path, label: str):
    det = load_detection_result(run_dir / "detection.npz")
    return read_confirmation_detection(det, label=label)


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warnings
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # 1. anchor
    h405_s1 = RUNS_ROOT / member_run_dir_name("s1")
    cfg = RunDirectory(h405_s1).load_cfg()
    anchor = _row("h405_s1", "s1", _pipeline(cfg), cfg.seed,
                  _read(h405_s1, "h405_s1"), ab, ked)
    _check_4dp("[anchor] h405 detfix s1", anchor,
               _csv_row(DETFIX_TABLE_CSV, "s1_new"), ANCHOR_COLS)
    ref6000 = _csv_row(DETFIX_TABLE_CSV, "pooled_new6000")
    print("ANCHOR PASSED: h405 detfix s1 rescored = committed s1_new (4 dp)")

    rows: list[dict[str, Any]] = []
    anchor.update({"cell": "h405", "stage": "ref"})
    rows.append(anchor)

    # 2. Stage A
    stage_a = []
    for cell in list(CELLS):
        rd = run_dir_for(cell, STAGE_A_MEMBER)
        if not (rd / "detection.npz").exists():
            print(f"  [{cell}] not yet run")
            continue
        c = RunDirectory(rd).load_cfg()
        if _pipeline(c) != "skip+closure:partner_aware":
            raise AssertionError(f"[{cell}] pipeline {_pipeline(c)}")
        r = _row(f"{cell}_s1", STAGE_A_MEMBER, _pipeline(c), c.seed,
                 _read(rd, cell), ab, ked)
        spec = CELLS[cell]
        r.update({"cell": cell, "stage": "A", "tau_ps": spec.tau_ps,
                  "E0_eV": spec.e0_eV, "drag_gate_A": c.drag_gate_steepness})
        rows.append(r)
        if spec.ridge_law:
            stage_a.append(r)
    expected = len(CELLS)
    done = sum(1 for r in rows if r["stage"] == "A")
    verdict: dict[str, Any] = {"stage_a_done": done, "stage_a_expected": expected,
                               "n_gated": sum(int(r["md_gate"]) for r in stage_a)}
    sel = select_stage_b(stage_a) if done == expected else None
    verdict["stage_b_cell"] = sel or ""
    print(f"Stage A: {done}/{expected} runs; gated ridge cells "
          f"{verdict['n_gated']}; Stage-B selection: {sel}")

    # 3. Stage B
    if sel:
        members = [run_dir_for(sel, m) for m in STAGE_B_MEMBERS]
        if all((m / "detection.npz").exists() for m in members):
            reads = [_read(run_dir_for(sel, STAGE_A_MEMBER), f"{sel}_s1")]
            for m, rd in zip(STAGE_B_MEMBERS, members):
                c = RunDirectory(rd).load_cfg()
                rd_read = _read(rd, f"{sel}_{m}")
                reads.append(rd_read)
                r = _row(f"{sel}_{m}", m, _pipeline(c), c.seed, rd_read, ab, ked)
                r.update({"cell": sel, "stage": "B"})
                rows.append(r)
            pooled = _row(f"{sel}_pooled6000", "pooled6000",
                          "skip+closure:partner_aware", "-",
                          pool_confirmation_reads(reads, label=f"{sel}_pooled6000"),
                          ab, ked)
            pooled.update({"cell": sel, "stage": "B_pooled"})
            rows.append(pooled)
            verdict["decision"] = decide(pooled)
            for k in ("nbar_det", "n1_solv", "w1_solv", "midHot", "deepKE",
                      "chi2_med", "KE1_mean", "KE2_mean", "trap", "S", "md_gate"):
                verdict[f"pooled_{k}"] = pooled[k]
                verdict[f"ref6000_{k}"] = ref6000[k]
            print(f"Stage B pooled N=6000 ({sel}): decision {verdict['decision']}")
        else:
            print(f"Stage B ({sel}): members not yet run")

    _write(rows, TABLE_CSV)
    _write([verdict], VERDICT_CSV)
    cols = ("label", "md_gate", "n1_solv", "nbar_det", "KE1_mean", "midHot",
            "w1_solv", "deepKE", "chi2_med", "trap", "S")
    for r in rows:
        print("  " + "  ".join(
            f"{k}={float(r[k]):.4g}" if k != "label" else f"{r[k]}"
            for k in cols))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
