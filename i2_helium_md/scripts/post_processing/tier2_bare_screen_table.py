"""G8 scorer — bare-start channel screen (what leaves, at which n, with what KE).

THESIS_FINAL_GOALS.md G8 (frozen 2026-10-06). Pure scorer: runs no MD.
Reported, not ranked: no weight scan, no pass/fail rule.

Order (anchors first; any failure stops before a new number):

1. **Anchors** — each parent run rescored through the same row builder
   reproduces its committed row (4 dp): h405 detfix s1 = ``s1_new``
   (``tier2_detfix_battery_table``); lq A5 s1 = ``c200t64confs1``
   (``g414_table.csv``).
2. **Rows** (``g8_bare_table.csv``) — per law, the bare run and its intact
   parent (same seed, same births): the production observable vector + band
   columns, plus fates over **all** ions (scored / each state reason) and the
   bare-class shares (n = 0 among scored, n >= 4 among solvated).
3. **By n** (``g8_bare_by_n.csv``) — per law and arm (bare / intact parent),
   per detected n: count, fraction of scored, solvated fraction, KE mean /
   median / IQR; the experiment's solvated fraction and KE mean / median
   beside it. Bare rows add the mechanics from ``bare_summary.npz``: mean
   pickup exposure, exit time, n at handover, max E_int.

    python scripts/post_processing/tier2_bare_screen_table.py
"""

from __future__ import annotations

import csv
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

N_SHOW_MAX = 21            # by-n rows run 0..21 (rows with no ions in any arm skipped)
TABLE_NAME = "g8_bare_table.csv"
BY_N_NAME = "g8_bare_by_n.csv"


def _write(rows: list[dict[str, Any]], path: Path) -> None:
    keys: list[str] = []
    for r in rows:
        keys.extend(k for k in r if k not in keys)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def fate_columns(state_reason: np.ndarray) -> dict[str, float]:
    """Fraction of **all** ions per terminal state reason, plus ``scored``."""
    from i2_helium_md.simulation.detection_stage import RETAINED_REASONS, STATE_REASONS

    state = np.asarray(state_reason)
    out = {f"fate_{s}": float(np.mean(state == s)) for s in STATE_REASONS}
    out["fate_scored"] = float(np.mean(~np.isin(state, sorted(RETAINED_REASONS))))
    return out


def scored_mask_and_n(det) -> tuple[np.ndarray, np.ndarray]:
    """Per-ion scored mask (retained excluded) and scored n (suppressed -> 0),
    the §4r conventions of ``read_confirmation_detection``."""
    from i2_helium_md.simulation.detection_stage import RETAINED_REASONS

    state = np.asarray(det.state_reason)
    keep = ~np.isin(state, sorted(RETAINED_REASONS))
    n = np.asarray(det.n_detected, dtype=float).copy()
    n[state == "suppressed"] = 0.0
    return keep, n[keep]


def ke_stats(ke: np.ndarray) -> dict[str, float]:
    if ke.size == 0:
        return {"KE_mean": np.nan, "KE_median": np.nan, "KE_iqr": np.nan}
    q1, med, q3 = np.percentile(ke, [25, 50, 75])
    return {"KE_mean": float(ke.mean()), "KE_median": float(med), "KE_iqr": float(q3 - q1)}


def main() -> int:
    from i2_helium_md.postprocess.abundance_loader import load_he_abundance_reference
    from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference
    from i2_helium_md.postprocess.tier2_confirmation import (
        read_confirmation_detection, solvated_renormalized,
    )
    from i2_helium_md.simulation.detection_stage import load_detection_result
    from i2_helium_md.simulation.run_directory import RunDirectory
    from scripts.gen_tier2_bare_screen import (
        LAWS, SUMMARY_NAME, parent_run_dir, run_dir_for,
    )
    from scripts.post_processing.tier2_detfix_battery_table import (
        TABLE_CSV as DETFIX_TABLE_CSV, _check_4dp, _csv_row, _pipeline, _row,
    )
    from scripts.post_processing.tier2_g413_table import ANCHOR_COLS
    from scripts.post_processing.tier2_g414_table import band_columns
    from scripts.post_processing.tier2atlas_geometry_table import (
        ABUNDANCE_REFERENCE_CSV, IHE_KED_REFERENCE_CSV, RUNS_ROOT,
    )

    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warnings
    out_dir = RUNS_ROOT / "h2b_forward_model"
    g414_csv = out_dir / "g414_table.csv"
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    exp_n, exp_f = solvated_renormalized(ab.n, ab.ion_fraction, n_min=1)
    exp_solv = dict(zip(np.asarray(exp_n).tolist(), np.asarray(exp_f).tolist()))
    exp_ke = {int(n): (float(m), float(md)) for n, m, md in
              zip(ked.n, ked.mean_KE_eV, ked.median_KE_eV)}
    anchor_rows = {"h405": (DETFIX_TABLE_CSV, "s1_new"), "lq": (g414_csv, "c200t64confs1")}

    rows: list[dict[str, Any]] = []
    by_n: list[dict[str, Any]] = []
    for law in LAWS:
        arms = {"intact": parent_run_dir(law), "bare": run_dir_for(law)}
        if not (arms["bare"] / "detection.npz").exists():
            print(f"[{law}] bare run not yet complete")
            continue
        per_arm = {}
        for arm, rd in arms.items():
            cfg = RunDirectory(rd).load_cfg()
            det = load_detection_result(rd / "detection.npz")
            read = read_confirmation_detection(det, label=f"{law}_{arm}")
            sn, sf = solvated_renormalized(read.n_values, read.fraction, n_min=1)
            row = {**_row(f"g8_{law}_{arm}", "s1", _pipeline(cfg), cfg.seed, read, ab, ked),
                   **band_columns(sn, sf, exp_n, exp_f)}
            if arm == "intact":
                path, label = anchor_rows[law]
                _check_4dp(f"[anchor] {law} parent", row, _csv_row(path, label), ANCHOR_COLS)
                print(f"ANCHOR PASSED: {law} parent rescored = committed {label} (4 dp)")
            keep, n_sc = scored_mask_and_n(det)
            if not np.array_equal(np.sort(n_sc), np.sort(np.asarray(read.n_scored))):
                raise AssertionError(f"[{law}/{arm}] scored n differs from the scorer read")
            solv = n_sc >= 1
            row.update({"law": law, "arm": arm, "initial_shell_model": cfg.initial_shell_model,
                        **fate_columns(det.state_reason),
                        "n0_share_scored": float(np.mean(n_sc == 0)),
                        "n4plus_share_solv": float(np.mean(n_sc[solv] >= 4)),
                        "nbar_solv": float(n_sc[solv].mean())})
            rows.append(row)
            ke_sc = np.asarray(det.E_kin_detected_eV, dtype=float)[keep]
            summ = None
            if arm == "bare":
                with np.load(rd / SUMMARY_NAME) as z:
                    summ = {k: np.asarray(z[k])[keep] for k in z.files}
            per_arm[arm] = (n_sc, ke_sc, summ)
            print(f"[{law}/{arm}] scored {row['fate_scored']:.3f}, retained "
                  f"{row['trap']:.3f}, n=0 {row['n0_share_scored']:.3f}, n1_solv "
                  f"{row['n1_solv']:.3f}, n>=4 {row['n4plus_share_solv']:.3f}, KE1 "
                  f"{row['KE1_mean']:.3f}, KE2 {row['KE2_mean']:.3f}, midHot {row['midHot']:.3f}")
        for n in range(N_SHOW_MAX + 1):
            for arm, (n_sc, ke_sc, summ) in per_arm.items():
                m = n_sc == n
                if not m.any():
                    continue
                e_mean, e_med = exp_ke.get(n, (np.nan, np.nan))
                r = {"law": law, "arm": arm, "n": n, "count": int(m.sum()),
                     "frac_scored": float(m.mean()),
                     "frac_solv": float(m.sum() / max(1, int((n_sc >= 1).sum())))
                     if n >= 1 else np.nan,
                     **ke_stats(ke_sc[m]),
                     "exp_frac_solv": exp_solv.get(n, np.nan) if n >= 1 else np.nan,
                     "exp_KE_mean": e_mean, "exp_KE_median": e_med}
                if summ is not None:
                    r.update({"exposure_ps_mean": float(summ["exposure_ps"][m].mean()),
                              "t_exit_ps_mean": float(np.nanmean(summ["t_exit_ps"][m])),
                              "n_handover_mean": float(summ["n_handover"][m].mean()),
                              "E_int_max_eV_mean": float(summ["E_int_max_eV"][m].mean())})
                by_n.append(r)
    if rows:
        _write(rows, out_dir / TABLE_NAME)
        _write(by_n, out_dir / BY_N_NAME)
        print(f"wrote {TABLE_NAME} ({len(rows)} rows), {BY_N_NAME} ({len(by_n)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
