"""(a, τ) joint ring scorer (plan §6.6 frozen reads JR-P1..JR-P5).

Scores the finished `gen_tier2atlas_linjoint.py` cells (6 × N = 500, all
at the §6.3 ring seed 20260731) together with the committed context rows
they are CRN-paired against — `lr6` (the a-curve's a = 35 anchor), the
§6.5 clone pool, and the `h405p` partner (the bar).

Reads (plan §6.6, pre-registered before launch):

* **JR-P1** — the decomposition. τ owns the clone's W₁ damage ⟺
  `s425` W₁ ≤ 0.85 **and** `d2` W₁ ≥ 0.95. Killed if `s425` ≥ 1.0 (the
  a-trend broke) or `d2` ≤ 0.85 (a × τ interaction). Anything else is
  reported as mixed and claims nothing.
* **JR-P2** — the a-curve at τ 4.8 over a ∈ {35, 37.5, 40, 42.5, 45}.
* **JR-P3** — the joint target: any cell beating `h405p` on **both**
  W₁ and KE₁ with midHot ≤ 1.15.
* **JR-P4** — if JR-P3 is empty, the (W₁, KE₁) Pareto front over the
  whole measured family is the result, naming the boundary.
* **JR-P5** — the tail diagnostic: scored n ≥ 10 fraction, max n, and
  the per-ion CRN shell flow vs `lr6`.

**The hard gate is REPORTED, never selecting** (user call 2026-08-12):
every cell in this family sits at n₁ 0.19–0.21 against a 0.19 floor
because of the shelved (C) source-side deficit that h405 shares, so
gating on it discriminates nothing *within* the family. Selection is on
W₁ + KE₁ with midHot as the cost line. The §6.5 three-way band verdict
is still printed per cell under both retained-policy arms.

Absolute values carry the per-seed SDs the §6.5 battery measured in this
exact family at N = 500 (W₁ 0.031, n₁ 0.0181, n̄ 0.099, KE₁ 0.0022,
midHot 0.0024); **cell-to-cell differences are CRN-paired and far
tighter** — every cell here shares seed 20260731.

Pure scorer — runs nothing, mutates nothing. Full vector to the
committed `atlas_linjoint_table.csv`.

Invocation::

    python scripts/post_processing/tier2atlas_linjoint_table.py
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Optional

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import csv  # noqa: E402

from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    load_confirmation_run,
    pool_confirmation_reads,
)
from i2_helium_md.simulation.detection_stage import (  # noqa: E402
    RETAINED_BOUND_REASON,
    RETAINED_MARGINAL_REASON,
    load_detection_result,
)
from scripts.gen_tier2atlas_linclone import (  # noqa: E402
    MEMBER_SEEDS as CLONE_SEEDS,
    clone_run_dir_name,
)
from scripts.gen_tier2atlas_linjoint import (  # noqa: E402
    EB_TAG,
    JOINT_MATRIX,
    JOINT_TWIN_COLS,
    JOINT_TWIN_ROWS,
    joint_run_dir_name,
    verify_joint_preregistration,
)
from scripts.gen_tier2atlas_linring import ring_run_dir_name  # noqa: E402
from scripts.post_processing.tier2atlas_g4step2_battery_table import (  # noqa: E402
    lowke_columns,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    ORACLE_RUN,
    RUNS_ROOT,
    check_oracle,
    observable_columns,
    observable_row,
)
from scripts.post_processing.tier2atlas_linclone_table import (  # noqa: E402
    arm_b_columns,
    arm_b_read,
    band_verdict,
    committed_row,
)
from scripts.post_processing.tier2atlas_linring_table import (  # noqa: E402
    GATE_N1,
    GATE_NBAR,
    hard_gate,
    ke_path_oracle,
    ke_width_columns,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS — frozen §6.6 reads
# ---------------------------------------------------------------------------

# JR-P1 decomposition thresholds.
JRP1_S425_TAU_OWNS_MAX = 0.85     # s425 at or below => consistent with "tau"
JRP1_S425_ATREND_BROKE = 1.00     # s425 at or above => the a-trend broke
JRP1_D2_TAU_OWNS_MIN = 0.95       # d2 at or above  => consistent with "tau"
JRP1_D2_INTERACTION_MAX = 0.85    # d2 at or below  => a x tau interaction

# JR-P3 cost line.
JRP3_MAX_MIDHOT = 1.15

# Per-seed SDs measured by the §6.5 battery in this family at N = 500.
SEED_SD = {"w1_solv": 0.031, "n1_solv": 0.0181, "nbar_det": 0.099,
           "KE1_mean": 0.0022, "midHot": 0.0024}

LINRING_CSV = RUNS_ROOT / "h2b_forward_model" / "atlas_linring_table.csv"
PARTNER_LABEL = "h405p"
CURVE_ANCHOR = "lr6"              # the a = 35 / tau 4.8 point of the a-curve

# Committed-row anchors re-derived before any context row is used.
ANCHOR_COLS = ("KE1_mean", "n1_solv", "nbar_det", "w1_solv", "midHot",
               "deepKE", "trap", "chi2_med")

SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_linjoint_table.csv"

TWIN_FORECAST = {
    "twin_trap": "trapped_frac", "twin_supp": "suppressed_frac",
    "twin_nbar": "nbar_det", "twin_n1": "n1_solv", "twin_w1": "w1_solv",
    "twin_KE1": "n1_ke_eV", "twin_KE2": "ke2_eV", "twin_deepKE": "deepke",
    "twin_midHot": "midhot_arith", "twin_phi": "phi", "twin_gate": "gate",
}


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------


def twin_columns(label: str) -> dict[str, float]:
    """The frozen §6.6 twin forecast for one cell, as floats."""
    raw = dict(zip(JOINT_TWIN_COLS, JOINT_TWIN_ROWS[label]))
    return {out: float(raw[src]) for out, src in TWIN_FORECAST.items()}


def tail_columns(read) -> dict[str, Any]:
    """JR-P5 tail diagnostic from one scored read."""
    n = np.asarray(read.n_scored)
    if n.size == 0:
        return {"tail_ge10": float("nan"), "tail_ge13": float("nan"),
                "n_max": float("nan")}
    return {
        "tail_ge10": float((n >= 10).mean()),
        "tail_ge13": float((n >= 13).mean()),
        "n_max": int(n.max()),
    }


def shell_flow_vs_anchor(cell_dir: Path, anchor_dir: Path) -> dict[str, float]:
    """JR-P5 per-ion CRN shell flow vs the a-curve anchor (`lr6`).

    Same seed + same N + same master draw => fragment index i is the same
    ion in both runs. Reports what happens to the anchor's heavily-dressed
    population (n >= 10): trapped, or shed down?
    """
    det_c = load_detection_result(cell_dir / "detection.npz")
    det_a = load_detection_result(anchor_dir / "detection.npz")
    retained = (RETAINED_BOUND_REASON, RETAINED_MARGINAL_REASON)
    n_c = np.asarray(det_c.n_detected, dtype=float)
    n_a = np.asarray(det_a.n_detected, dtype=float)
    trap_c = np.isin(np.asarray(det_c.state_reason), retained)
    trap_a = np.isin(np.asarray(det_a.state_reason), retained)

    hi = (n_a >= 10) & ~trap_a
    out: dict[str, float] = {"flow_hi_n": int(hi.sum())}
    if not hi.any():
        return out | {"flow_hi_trapped": 0, "flow_hi_dn": float("nan"),
                      "flow_hi_still": 0}
    surv = hi & ~trap_c
    out["flow_hi_trapped"] = int((hi & trap_c).sum())
    out["flow_hi_dn"] = (float((n_c[surv] - n_a[surv]).mean())
                         if surv.any() else float("nan"))
    out["flow_hi_still"] = int((n_c[surv] >= 10).sum()) if surv.any() else 0
    return out


def anchor_oracle(row: dict[str, Any], label: str) -> list[str]:
    """A committed ring row must re-derive before it is used as context."""
    committed = committed_row(LINRING_CSV, label)
    problems = []
    for col in ANCHOR_COLS:
        want, got = float(committed[col]), float(row[col])
        if not np.isclose(got, want, rtol=1e-9, atol=0.0):
            problems.append(f"{label}.{col}: rescored {got!r} vs committed "
                            f"{want!r}")
    return problems


def pareto_front(points: list[tuple[str, float, float]]
                 ) -> list[tuple[str, float, float]]:
    """Non-dominated (label, W1, KE1) points: minimise W1, maximise KE1."""
    front = []
    for lab, w1, ke in points:
        if not np.isfinite(w1) or not np.isfinite(ke):
            continue
        dominated = any(
            (w2 <= w1 and k2 >= ke) and (w2 < w1 or k2 > ke)
            for lab2, w2, k2 in points
            if lab2 != lab and np.isfinite(w2) and np.isfinite(k2)
        )
        if not dominated:
            front.append((lab, w1, ke))
    return sorted(front, key=lambda t: t[1])


def score_run(label: str, run_dir: Path, abundance_ref, ked_ref,
              read=None) -> dict[str, Any]:
    """One full row: observable vector + KE + width + Arm B + tail."""
    if read is None:
        read = load_confirmation_run(run_dir, label=label)
    row: dict[str, Any] = {"label": label}
    row.update(observable_columns(read, abundance_ref, ked_ref))
    row.update(lowke_columns(read))
    row.update(ke_width_columns(read))
    row.update(tail_columns(read))
    if run_dir is not None:
        b_read, n_marg = arm_b_read(run_dir, label)
        row.update(arm_b_columns(b_read, n_marg, abundance_ref))
    row["md_gate_A"] = hard_gate(float(row["n1_solv"]), float(row["nbar_det"]))
    if "gate_B" in row:
        row["policy_blocked"] = int(row["md_gate_A"] != row["gate_B"])
    return row


# ---------------------------------------------------------------------------


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # ---- Oracles before any new number -----------------------------------
    verify_joint_preregistration()              # LJ-P1 + LR-P1
    oracle_row = observable_row("pooled(oracle)", RUNS_ROOT / ORACLE_RUN,
                                abundance_ref, ked_ref, with_geometry=False)
    drift = check_oracle(oracle_row)
    if drift:
        print("*** SCORER-DRIFT ORACLE FAILED — do not read the ring ***")
        for line in drift:
            print(f"    {line}")
        return
    ke_drift = ke_path_oracle(abundance_ref, ked_ref)
    if ke_drift:
        print("*** KE-PATH ORACLE FAILED — do not read the ring ***")
        for line in ke_drift:
            print(f"    {line}")
        return

    context: dict[str, dict[str, Any]] = {}
    for label in (PARTNER_LABEL, CURVE_ANCHOR):
        d = RUNS_ROOT / ring_run_dir_name(label)
        if not (d / "detection.npz").exists():
            print(f"*** committed context run missing at {d} ***")
            return
        context[label] = score_run(label, d, abundance_ref, ked_ref)
        problems = anchor_oracle(context[label], label)
        if problems:
            print(f"*** ANCHOR ORACLE FAILED for {label} — do not read ***")
            for line in problems:
                print(f"    {line}")
            return
    print("oracles OK (LJ-P1 + LR-P1 + drift + KE-path + committed "
          f"{PARTNER_LABEL}/{CURVE_ANCHOR} anchors).")
    print()

    partner = context[PARTNER_LABEL]
    anchor_dir = RUNS_ROOT / ring_run_dir_name(CURVE_ANCHOR)

    # ---- Ring rows --------------------------------------------------------
    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    for cell in JOINT_MATRIX:
        run_dir = RUNS_ROOT / joint_run_dir_name(cell.label)
        if not (run_dir / "detection.npz").exists():
            missing.append(cell.label)
            continue
        row = score_run(cell.label, run_dir, abundance_ref, ked_ref)
        row.update({"a": cell.a, "well": EB_TAG, "tau": cell.tau_ps,
                    "E0": cell.e0_eV, "role": cell.role})
        row.update(twin_columns(cell.label))
        row.update(shell_flow_vs_anchor(run_dir, anchor_dir))
        row["d_w1_h405"] = float(row["w1_solv"]) - float(partner["w1_solv"])
        row["d_KE1_h405"] = float(row["KE1_mean"]) - float(partner["KE1_mean"])
        rows.append(row)

    if missing:
        print(f"not yet run: {missing}")
    if not rows:
        print("(no finished cells)")
        return

    # Context rows: the clone pool + the committed anchors.
    clone_dirs = [RUNS_ROOT / clone_run_dir_name(m) for m in CLONE_SEEDS]
    if all((d / "detection.npz").exists() for d in clone_dirs):
        pooled = pool_confirmation_reads(
            [load_confirmation_run(d, label=d.name) for d in clone_dirs],
            label="clonepool")
        clone_row = score_run("clone(pool)", None, abundance_ref, ked_ref,
                              read=pooled)
        clone_row.update({"a": 42.5, "well": EB_TAG, "tau": 6.4, "E0": 0.31,
                          "role": "§6.5 clone (context)"})
        context["clone(pool)"] = clone_row

    for lab, row in context.items():
        row.setdefault("role", "context")
        row["d_w1_h405"] = float(row["w1_solv"]) - float(partner["w1_solv"])
        row["d_KE1_h405"] = float(row["KE1_mean"]) - float(partner["KE1_mean"])
    context[CURVE_ANCHOR].update({"a": 35.0, "tau": 4.8, "E0": 0.37})
    context[PARTNER_LABEL].update({"a": "-", "tau": 4.4, "E0": 0.405})

    table_rows = rows + [context[k] for k in context]
    key_set = set().union(*(r.keys() for r in table_rows))
    for r in table_rows:
        for k in key_set:
            r.setdefault(k, "-")

    print("=== (a, tau) joint ring (plan §6.6), N = 500, seed 20260731 "
          "(CRN with the whole family) ===")
    print(format_table(table_rows))
    print()

    by = {r["label"]: r for r in table_rows}

    # ---- JR-P1: the decomposition ----------------------------------------
    print("=== Verdicts (plan §6.6, pre-registered before launch) ===")
    s425 = by.get("s425")
    d2 = by.get("d2")
    if s425 is None or d2 is None:
        print("JR-P1: decomposition cells not both finished — no verdict.")
    else:
        w425, wd2 = float(s425["w1_solv"]), float(d2["w1_solv"])
        tau_owns = (w425 <= JRP1_S425_TAU_OWNS_MAX
                    and wd2 >= JRP1_D2_TAU_OWNS_MIN)
        atrend_broke = w425 >= JRP1_S425_ATREND_BROKE
        interaction = wd2 <= JRP1_D2_INTERACTION_MAX
        print(f"JR-P1 (decomposition): s425 W1 {w425:.4f} "
              f"(tau-owns if <= {JRP1_S425_TAU_OWNS_MAX}, a-trend broke if "
              f">= {JRP1_S425_ATREND_BROKE}); d2 W1 {wd2:.4f} "
              f"(tau-owns if >= {JRP1_D2_TAU_OWNS_MIN}, interaction if "
              f"<= {JRP1_D2_INTERACTION_MAX})")
        if tau_owns:
            print("  => TAU OWNS the clone's W1 damage — confirmed on both "
                  "cells. The a-dial is clean and the clock is the cost.")
        elif interaction and not atrend_broke:
            print("  => a x TAU INTERACTION: tau 6.4 is harmless at a 35, so "
                  "the clone's damage needed the high chord too — j1 is the "
                  "live joint candidate.")
        elif atrend_broke:
            print("  => THE a-TREND BROKE at 42.5 independently of tau; the "
                  "W1 gain from a saturates below the clone's chord.")
        else:
            print("  => MIXED — neither leg of the pre-registered condition "
                  "is met; claims nothing (plan §6.6).")

    # ---- JR-P2: the a-curve at tau 4.8 ------------------------------------
    curve = [(35.0, CURVE_ANCHOR)] + [(c.a, c.label) for c in JOINT_MATRIX
                                      if c.tau_ps == 4.8]
    print("JR-P2 (a-curve at tau 4.8):")
    for a_val, lab in sorted(curve):
        r = by.get(lab)
        if r is None:
            continue
        print(f"  a {a_val:5.1f} ({lab:5s}): W1 {float(r['w1_solv']):.4f}  "
              f"KE1 {float(r['KE1_mean']):.4f}  midHot "
              f"{float(r['midHot']):.3f}  chi2 {float(r['chi2_med']):7.1f}  "
              f"tail>=10 {float(r['tail_ge10']):.4f}  n_max {r['n_max']}  "
              f"n1 {float(r['n1_solv']):.4f}  nbar {float(r['nbar_det']):.3f}")

    # ---- JR-P3 / JR-P4 ----------------------------------------------------
    bar_w1 = float(partner["w1_solv"])
    bar_ke = float(partner["KE1_mean"])
    winners = [r for r in rows
               if float(r["w1_solv"]) < bar_w1
               and float(r["KE1_mean"]) > bar_ke
               and float(r["midHot"]) <= JRP3_MAX_MIDHOT]
    print(f"JR-P3 (joint target: W1 < {bar_w1:.4f} AND KE1 > {bar_ke:.4f} "
          f"AND midHot <= {JRP3_MAX_MIDHOT}):")
    if winners:
        for r in winners:
            print(f"  *** SUCCESS {r['label']}: W1 {float(r['w1_solv']):.4f} "
                  f"({float(r['d_w1_h405']):+.4f}), KE1 "
                  f"{float(r['KE1_mean']):.4f} "
                  f"({float(r['d_KE1_h405']):+.4f}), midHot "
                  f"{float(r['midHot']):.3f}")
    else:
        print("  none — JR-P4 applies: the (W1, KE1) Pareto front IS the "
              "result (the family's measured boundary).")

    pts = [(r["label"], float(r["w1_solv"]), float(r["KE1_mean"]))
           for r in table_rows]
    with open(LINRING_CSV, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["label"] not in {p[0] for p in pts}:
                pts.append((r["label"], float(r["w1_solv"]),
                            float(r["KE1_mean"])))
    print("JR-P4 (Pareto front over the whole measured family, "
          "minimise W1 / maximise KE1):")
    for lab, w1, ke in pareto_front(pts):
        print(f"  {lab:12s} W1 {w1:.4f}  KE1 {ke:.4f}")

    # ---- JR-P5 tail diagnostic + the reported gate ------------------------
    print("JR-P5 (tail diagnostic + CRN shell flow vs lr6):")
    for r in rows:
        print(f"  {r['label']:5s}: tail>=10 {float(r['tail_ge10']):.4f} "
              f"(>=13 {float(r['tail_ge13']):.4f}), n_max {r['n_max']}, "
              f"supp {float(r['supp']):.3f}; of lr6's {int(r['flow_hi_n'])} "
              f"n>=10 ions: {int(r['flow_hi_trapped'])} trapped, mean dn "
              f"{float(r['flow_hi_dn']):+.2f}, {int(r['flow_hi_still'])} "
              "still n>=10")

    print(f"Gate (REPORTED, never selecting — user call 2026-08-12; "
          f"n1 {list(GATE_N1)} x nbar {list(GATE_NBAR)}, three-way band "
          "verdict at the measured per-seed SD):")
    for r in rows:
        v_n1 = band_verdict(float(r["n1_solv"]), GATE_N1, SEED_SD["n1_solv"])
        v_nb = band_verdict(float(r["nbar_det"]), GATE_NBAR,
                            SEED_SD["nbar_det"])
        blocked = r.get("policy_blocked", 0)
        print(f"  {r['label']:5s}: n1 {float(r['n1_solv']):.4f} -> {v_n1}; "
              f"nbar {float(r['nbar_det']):.3f} -> {v_nb}"
              + ("  [policy-blocked]" if blocked == 1 else ""))

    print(f"\nDisclosure: absolute values carry the §6.5 per-seed SDs "
          f"{SEED_SD}; cell-to-cell differences are CRN-paired and tighter.")
    if SAVE_CSV_PATH is not None:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, table_rows)}")


if __name__ == "__main__":
    main()
