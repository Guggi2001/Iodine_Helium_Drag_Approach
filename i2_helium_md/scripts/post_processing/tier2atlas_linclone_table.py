"""h405-clone MD battery scorer (plan §6.5 frozen reads CL-P1..CL-P6).

Scores the finished `gen_tier2atlas_linclone.py` runs — one cell
(`pure_linear a 42.5 / eb0482 / τ 6.4 / E₀ 0.31`, the §6.4 item-2 clone)
at three seeds, N = 500 — plus the committed §6.3 ring partner `h405p`,
which is re-scored here (never read from its CSV) and CRN-pairs with the
seed-20260731 member.

This is an **equivalence** test, not a gain test: the ring's K-KE kill is
inverted here — the clone *wants* ΔKE₁ ≈ 0 against h405.

Frozen reads (plan §6.5, pre-registered before launch):

* **CL-P1** — hard gate on the POOLED ensemble under both retained-policy
  arms: n1_solv ∈ [0.19, 0.30] ∧ nbar_det ∈ [3.77, 4.37]. **Three-way**:
  a pooled value within 1 SE of a band edge returns `gate-marginal`,
  neither pass nor fail, where SE is this battery's own measured
  per-seed SD/√3. No single seed decides.
* **CL-P2** — |ΔKE₁| ≤ 0.05 eV vs h405 (primary reference = the committed
  g4 Step-2 pooled battery, whose per-seed KE₁ SD is 0.0024; the CRN pair
  against `h405p` is the cross-check).
* **CL-P3** — midHot ≤ 1.15: does the §6.3 ring's mid-band overheating
  (1.39–1.98) vanish at the clone, or is it intrinsic to constant γ?
* **CL-P4** — trap reported decomposed, never gating (§4.1); flagged on a
  non-zero marginal class or Δtrap > +0.05 vs h405.
* **CL-P5** — W₁ **reported-only** (user 2026-08-12); a pooled loss
  > +0.10 is an honest cost line and cannot overturn CL-P1..P3.
* **CL-P6** — the twin−MD deltas at a 42.5 / τ 6.4 vs the §6.3
  ring-measured in-family bias ranges. Lands regardless of the verdict.

χ²_med is **N-extensive** and is therefore compared to the N = 500
partner only, never to the pooled N = 1000 battery number (plan §6.5).

Any CL-P1/P2/P3 miss routes to §6.2 item 1 recalibration (bias-corrected
twin re-scan, zero MD), never to a new ring.

Pure scorer — runs nothing, mutates nothing. Per-seed + pooled rows with
the full observable vector, every KE column and the twin-forecast columns
go to the committed `atlas_linclone_table.csv`.

Invocation::

    python scripts/post_processing/tier2atlas_linclone_table.py
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
    read_confirmation_detection,
    score_histogram_vs_reference,
)
from scripts.gen_tier2atlas_linclone import (  # noqa: E402
    CLONE_CELL,
    CLONE_TWIN_COLS,
    CLONE_TWIN_ROW,
    CRN_MEMBER,
    MEMBER_SEEDS,
    clone_run_dir_name,
    verify_clone_preregistration,
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
from scripts.post_processing.tier2atlas_linring_table import (  # noqa: E402
    GATE_N1,
    GATE_NBAR,
    fate_flow,
    hard_gate,
    ke_path_oracle,
    ke_width_columns,
)
from scripts.post_processing.tier2atlas_retained_bracket import (  # noqa: E402
    arm_b_detection,
    marginal_dossier,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS — frozen bands (plan §6.5; GATE_* reused from the ring scorer)
# ---------------------------------------------------------------------------

CLONE_MAX_ABS_DKE1_EV = 0.05     # CL-P2 equivalence band
CLONE_MAX_MIDHOT = 1.15          # CL-P3 soft band top
TRAP_FLAG_DELTA = 0.05           # CL-P4 flag threshold vs h405
W1_COST_DELTA = 0.10             # CL-P5 cost-line threshold
W1_SINGLE_SEED_DISCLOSURE = 0.13

# CL-P6: the §6.3 ring-measured in-family bias ranges, TWIN MINUS MD.
RING_BIAS_TWIN_MINUS_MD = {
    "KE1": (-0.015, -0.010),     # twin uniformly ~1-2 % cold in-family
    "n1": (-0.006, 0.006),       # |transfer| <= 0.006
    "nbar": (0.46, 0.59),        # twin-hot, small end of the bracket
}

# Committed h405 references (read at runtime — gate-on-committed-artifacts).
BATTERY_CSV = RUNS_ROOT / "h2b_forward_model" / "atlas_g4step2_battery.csv"
LINRING_CSV = RUNS_ROOT / "h2b_forward_model" / "atlas_linring_table.csv"
PARTNER_LABEL = "h405p"

# Partner-anchor oracle columns (the committed ring row this scorer must
# re-derive before the partner is used as a Δ target).
PARTNER_ANCHOR_COLS = ("KE1_mean", "KE2_mean", "n1_solv", "nbar_det",
                       "w1_solv", "midHot", "deepKE", "trap", "chi2_med")

SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_linclone_table.csv"

TWIN_FORECAST = {
    "twin_trap": "trapped_frac", "twin_supp": "suppressed_frac",
    "twin_nbar": "nbar_det", "twin_n1": "n1_solv", "twin_w1": "w1_solv",
    "twin_KE1": "n1_ke_eV", "twin_KE2": "ke2_eV", "twin_deepKE": "deepke",
    "twin_midHot": "midhot_arith", "twin_phi": "phi", "twin_gate": "gate",
}


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------


def twin_row() -> dict[str, float]:
    """The frozen §6.4 clone twin forecast as floats (LC-P1 verifies the
    strings against the committed CSV before this is used)."""
    raw = dict(zip(CLONE_TWIN_COLS, CLONE_TWIN_ROW))
    return {out: float(raw[src]) for out, src in TWIN_FORECAST.items()}


def committed_row(csv_path: Path, label: str) -> dict[str, str]:
    """One committed scorer row by label (string fields)."""
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["label"] == label:
                return row
    raise AssertionError(f"no {label!r} row in {csv_path}")


def arm_b_read(run_dir: Path, label: str):
    """The Arm-B (marginal-injection) read for one run, or the Arm-A read when
    the cell has no marginal class (then Arm B is definitionally Arm A).

    Same §3.5b construction the §6.3 ring scorer uses; factored to a *read* so
    the three seeds can be pooled on the Arm-B end as well as the Arm-A end.
    """
    dossier = marginal_dossier(run_dir)
    if dossier is None:
        return load_confirmation_run(run_dir, label=label), 0
    return (read_confirmation_detection(arm_b_detection(dossier),
                                        label=f"{label}B"),
            int(dossier["marginal"].sum()))


def arm_b_columns(read, n_marg: int, abundance_ref) -> dict[str, float]:
    """Arm-B gate surface from an Arm-B read (per-seed or pooled)."""
    hist = score_histogram_vs_reference(read, abundance_ref)
    return {
        "n_marg": n_marg,
        "n1_B": float(hist.sim_n1),
        "nbar_B": float(read.n_mean),
        "w1_B": float(hist.w1_solvated),
        "gate_B": hard_gate(float(hist.sim_n1), float(read.n_mean)),
    }


def band_verdict(value: float, band: tuple[float, float], se: float) -> str:
    """CL-P1's three-way outcome for one observable (plan §6.5).

    ``pass`` inside the band by more than 1 SE; ``fail`` outside by more than
    1 SE; ``gate-marginal`` within 1 SE of either edge — neither pass nor
    fail, forced by the N = 500 cut.
    """
    lo, hi = band
    if not np.isfinite(value):
        return "nan"
    if not np.isfinite(se) or se <= 0.0:
        se = 0.0
    inside = lo <= value <= hi
    edge_gap = min(abs(value - lo), abs(value - hi))
    if edge_gap <= se:
        return "gate-marginal"
    return "pass" if inside else "fail"


def partner_anchor_oracle(row: dict[str, Any]) -> list[str]:
    """The re-scored h405p row must reproduce the committed
    `atlas_linring_table.csv` values before it is used as a Δ target."""
    committed = committed_row(LINRING_CSV, PARTNER_LABEL)
    problems = []
    for col in PARTNER_ANCHOR_COLS:
        want = float(committed[col])
        got = float(row[col])
        if not np.isclose(got, want, rtol=1e-9, atol=0.0):
            problems.append(f"{col}: rescored {got!r} vs committed {want!r}")
    return problems


# ---------------------------------------------------------------------------


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # ---- Oracles, all before any new number ------------------------------
    verify_clone_preregistration()          # LC-P1 + the ring's LR-P1
    oracle_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
    drift = check_oracle(oracle_row)
    if drift:
        print("*** SCORER-DRIFT ORACLE FAILED — do not read the battery ***")
        for line in drift:
            print(f"    {line}")
        return
    ke_drift = ke_path_oracle(abundance_ref, ked_ref)
    if ke_drift:
        print("*** KE-PATH ORACLE FAILED — do not read the battery ***")
        for line in ke_drift:
            print(f"    {line}")
        return
    print("scorer-drift + KE-path oracles OK.")

    partner_dir = RUNS_ROOT / ring_run_dir_name(PARTNER_LABEL)
    if not (partner_dir / "detection.npz").exists():
        print(f"*** committed ring partner missing at {partner_dir} ***")
        return
    partner = observable_row(PARTNER_LABEL, partner_dir, abundance_ref,
                             ked_ref, with_geometry=False)
    partner_read = load_confirmation_run(partner_dir, label=PARTNER_LABEL)
    partner.update(lowke_columns(partner_read))
    partner.update(ke_width_columns(partner_read))
    anchor_drift = partner_anchor_oracle(partner)
    if anchor_drift:
        print("*** PARTNER-ANCHOR ORACLE FAILED — the committed h405p row "
              "does not re-derive; do not read the battery ***")
        for line in anchor_drift:
            print(f"    {line}")
        return
    print("partner-anchor oracle OK (committed h405p row re-derived).")
    print()

    twin = twin_row()
    h405_pooled = committed_row(BATTERY_CSV, "pooled")

    # ---- Rows -------------------------------------------------------------
    rows: list[dict[str, Any]] = []
    reads = []
    arm_b_reads = []
    n_marg_total = 0
    missing: list[str] = []
    for member, seed in MEMBER_SEEDS.items():
        run_dir = RUNS_ROOT / clone_run_dir_name(member)
        if not (run_dir / "detection.npz").exists():
            missing.append(member)
            continue
        read = load_confirmation_run(run_dir, label=member)
        reads.append(read)
        row: dict[str, Any] = {"label": member, "seed": seed}
        row.update(observable_columns(read, abundance_ref, ked_ref))
        row.update(lowke_columns(read))
        row.update(ke_width_columns(read))
        row.update({"KE1_n": int((read.n_scored == 1).sum()),
                    "KE2_n": int((read.n_scored == 2).sum())})
        b_read, n_marg = arm_b_read(run_dir, member)
        arm_b_reads.append(b_read)
        n_marg_total += n_marg
        row.update(arm_b_columns(b_read, n_marg, abundance_ref))
        row["md_gate_A"] = hard_gate(float(row["n1_solv"]),
                                     float(row["nbar_det"]))
        row["policy_blocked"] = int(row["md_gate_A"] != row["gate_B"])
        rows.append(row)

    if missing:
        print(f"not yet run: {missing}")
    if not rows:
        print("(no finished members)")
        return

    pooled_read = pool_confirmation_reads(reads, label="clonepooled")
    pooled: dict[str, Any] = {"label": "pooled", "seed": "-"}
    pooled.update(observable_columns(pooled_read, abundance_ref, ked_ref))
    pooled.update(lowke_columns(pooled_read))
    pooled.update(ke_width_columns(pooled_read))
    pooled.update({"KE1_n": int((pooled_read.n_scored == 1).sum()),
                   "KE2_n": int((pooled_read.n_scored == 2).sum())})
    pooled_b = pool_confirmation_reads(arm_b_reads, label="clonepooledB")
    pooled.update(arm_b_columns(pooled_b, n_marg_total, abundance_ref))
    pooled["md_gate_A"] = hard_gate(float(pooled["n1_solv"]),
                                    float(pooled["nbar_det"]))
    pooled["policy_blocked"] = int(pooled["md_gate_A"] != pooled["gate_B"])

    # Twin join (§6.2 item 6 convention) + the CL-P6 twin-minus-MD deltas.
    for row in rows + [pooled]:
        row.update(twin)
        row.update({
            "tmd_KE1": float(twin["twin_KE1"]) - float(row["KE1_mean"]),
            "tmd_n1": float(twin["twin_n1"]) - float(row["n1_solv"]),
            "tmd_nbar": float(twin["twin_nbar"]) - float(row["nbar_det"]),
        })
        row["d_KE1_h405bat"] = (float(row["KE1_mean"])
                                - float(h405_pooled["KE1_mean"]))
        row["d_KE1_h405p"] = (float(row["KE1_mean"])
                              - float(partner["KE1_mean"]))
        row["d_w1_h405bat"] = (float(row["w1_solv"])
                               - float(h405_pooled["w1_solv"]))
        row["d_trap_h405bat"] = (float(row["trap"])
                                 - float(h405_pooled["trap"]))

    partner_out = dict(partner)
    partner_out.update({
        "label": PARTNER_LABEL, "seed": MEMBER_SEEDS[CRN_MEMBER],
        "KE1_n": int((partner_read.n_scored == 1).sum()),
        "KE2_n": int((partner_read.n_scored == 2).sum()),
        "md_gate_A": hard_gate(float(partner["n1_solv"]),
                               float(partner["nbar_det"])),
    })
    # The partner carries no Arm-B / twin-forecast / Δ columns (it IS the Δ
    # target and its own §6.3 Arm-B read lives in the committed ring CSV):
    # pad so the shared table keeps one column set.
    for key in pooled:
        partner_out.setdefault(key, "-")
    table_rows = rows + [pooled, partner_out]

    print(f"=== h405-clone MD battery (plan §6.5): pure_linear "
          f"a={CLONE_CELL.a} {CLONE_CELL.eb_tag} tau={CLONE_CELL.tau_ps} "
          f"E0={CLONE_CELL.e0_eV}, N=500 x {len(rows)} seeds ===")
    print(format_table(table_rows))
    print()

    # ---- Verdicts (plan §6.5 CL-P1..CL-P6) --------------------------------
    print("=== Verdicts (plan §6.5, pre-registered before launch) ===")

    def per_seed_se(col: str) -> float:
        vals = np.asarray([float(r[col]) for r in rows], dtype=float)
        if vals.size < 2:
            return float("nan")
        return float(vals.std(ddof=1) / np.sqrt(vals.size))

    se_n1 = per_seed_se("n1_solv")
    se_nbar = per_seed_se("nbar_det")
    v_n1 = band_verdict(float(pooled["n1_solv"]), GATE_N1, se_n1)
    v_nbar = band_verdict(float(pooled["nbar_det"]), GATE_NBAR, se_nbar)
    v_n1_B = band_verdict(float(pooled["n1_B"]), GATE_N1, se_n1)
    v_nbar_B = band_verdict(float(pooled["nbar_B"]), GATE_NBAR, se_nbar)
    order = {"fail": 0, "gate-marginal": 1, "pass": 2, "nan": 0}
    cl_p1_A = min((v_n1, v_nbar), key=lambda v: order[v])
    cl_p1_B = min((v_n1_B, v_nbar_B), key=lambda v: order[v])
    print(f"CL-P1 (pooled hard gate, measured SE from this battery):")
    print(f"  Arm A: n1 {float(pooled['n1_solv']):.4f} +-{se_n1:.4f} in "
          f"{list(GATE_N1)} -> {v_n1}; nbar {float(pooled['nbar_det']):.3f} "
          f"+-{se_nbar:.3f} in {list(GATE_NBAR)} -> {v_nbar}  => {cl_p1_A}")
    print(f"  Arm B: n1 {float(pooled['n1_B']):.4f} -> {v_n1_B}; "
          f"nbar {float(pooled['nbar_B']):.3f} -> {v_nbar_B}  => {cl_p1_B}"
          f"   (marginal class n={int(pooled['n_marg'])})")
    if cl_p1_A != cl_p1_B:
        print("  POLICY-BLOCKED: the two retained-policy ends disagree "
              "(§4.1) — escalates to the user.")
    print("  per-seed gate (Arm A): "
          + "; ".join(f"{r['label']} n1={float(r['n1_solv']):.4f} "
                      f"nbar={float(r['nbar_det']):.3f} "
                      f"gate={int(r['md_gate_A'])}" for r in rows))

    d_bat = float(pooled["d_KE1_h405bat"])
    d_crn = float(rows[0]["d_KE1_h405p"]) if rows else float("nan")
    crn_row = next((r for r in rows if r["label"] == CRN_MEMBER), None)
    if crn_row is not None:
        d_crn = float(crn_row["d_KE1_h405p"])
    cl_p2 = abs(d_bat) <= CLONE_MAX_ABS_DKE1_EV
    print(f"CL-P2 (equivalence, NOT a gain test): pooled KE1 "
          f"{float(pooled['KE1_mean']):.4f} vs h405 battery "
          f"{float(h405_pooled['KE1_mean']):.4f} -> dKE1 {d_bat:+.4f} eV, "
          f"|d| <= {CLONE_MAX_ABS_DKE1_EV} => {'PASS' if cl_p2 else 'MISS'}")
    print(f"  CRN cross-check ({CRN_MEMBER} vs the committed h405p, same seed "
          f"+ same N): dKE1 {d_crn:+.4f} eV")

    cl_p3 = float(pooled["midHot"]) <= CLONE_MAX_MIDHOT
    print(f"CL-P3 (the cost read): pooled midHot {float(pooled['midHot']):.4f} "
          f"<= {CLONE_MAX_MIDHOT} => {'PASS' if cl_p3 else 'MISS'}  "
          f"[§6.3 ring lin cells 1.39-1.98; h405 "
          f"{float(h405_pooled['midHot']):.4f}]")
    if not cl_p3:
        print("  => the mid-band overheating is FORM-INTRINSIC to constant "
              "gamma, not an a-dial artifact (plan §6.5).")

    d_trap = float(pooled["d_trap_h405bat"])
    trap_flag = (float(pooled["trap_marg"]) > 0.0) or (d_trap > TRAP_FLAG_DELTA)
    print(f"CL-P4 (trap, reported / never gating): pooled trap "
          f"{float(pooled['trap']):.4f} (bound {float(pooled['trap_bound']):.4f} "
          f"/ marg {float(pooled['trap_marg']):.4f}); dtrap vs h405 "
          f"{d_trap:+.4f}; twin floor {twin['twin_trap']:.4f} => "
          f"{'FLAGGED' if trap_flag else 'no flag'}")

    d_w1 = float(pooled["d_w1_h405bat"])
    print(f"CL-P5 (W1, REPORTED-ONLY): pooled W1 {float(pooled['w1_solv']):.4f} "
          f"vs h405 {float(h405_pooled['w1_solv']):.4f} -> {d_w1:+.4f}"
          + (f"  [COST LINE: loss > +{W1_COST_DELTA}]"
             if d_w1 > W1_COST_DELTA else "")
          + f"  (single-seed disclosure +-{W1_SINGLE_SEED_DISCLOSURE}; "
            "cannot overturn CL-P1..P3)")
    # chi2_med is N-EXTENSIVE: the pooled value sums over 3x the ions and is
    # NOT comparable to anything here (it is committed to the CSV, not read).
    # Only the per-seed values may be compared, and only to the N = 500
    # partner (plan §6.5).
    print(f"  chi2_med (N-EXTENSIVE — per-seed only, vs the N=500 partner "
          f"h405p {float(partner['chi2_med']):.1f}): "
          + "; ".join(f"{r['label']} {float(r['chi2_med']):.1f}" for r in rows)
          + f"   [pooled {float(pooled['chi2_med']):.1f} spans "
            f"{len(rows)}x500 ions and is NOT comparable]")

    print("CL-P6 (authority-box extension at a 42.5 / tau 6.4 — lands "
          "regardless of the verdict), twin minus MD:")
    for key, col in (("KE1", "tmd_KE1"), ("n1", "tmd_n1"),
                     ("nbar", "tmd_nbar")):
        lo, hi = RING_BIAS_TWIN_MINUS_MD[key]
        val = float(pooled[col])
        inside = lo <= val <= hi
        print(f"  {key}: twin-MD {val:+.4f} vs ring range [{lo:+.3f}, "
              f"{hi:+.3f}] => {'inside' if inside else 'OUTSIDE'}")
    if not all(RING_BIAS_TWIN_MINUS_MD[k][0] <= float(pooled[c])
               <= RING_BIAS_TWIN_MINUS_MD[k][1]
               for k, c in (("KE1", "tmd_KE1"), ("n1", "tmd_n1"),
                            ("nbar", "tmd_nbar"))):
        print("  => the in-family bias vector is a-/tau-DEPENDENT: the §6.3 "
              "box does not extend to this corner.")

    # CRN fate flow vs the committed partner (reported, never gating).
    if crn_row is not None:
        flow = fate_flow(RUNS_ROOT / clone_run_dir_name(CRN_MEMBER),
                         partner_dir)
        crn_row.update(flow)
        print(f"CRN fate flow ({CRN_MEMBER} vs h405p, same seed + same N): "
              f"new-trapped {int(flow['flow_new_trap'])}, freed "
              f"{int(flow['flow_freed'])}")

    clone_holds = (cl_p1_A == "pass" and cl_p1_B == "pass" and cl_p2 and cl_p3)
    hard_miss = (cl_p1_A == "fail" or cl_p1_B == "fail"
                 or not cl_p2 or not cl_p3)
    print()
    if not clone_holds and not hard_miss:
        # CL-P1's third outcome. Report what it would COST to resolve, so the
        # "just run more seeds" reflex is priced rather than assumed.
        n1_vals = np.asarray([float(r["n1_solv"]) for r in rows], dtype=float)
        sd = float(n1_vals.std(ddof=1)) if n1_vals.size > 1 else float("nan")
        gap = min(abs(float(pooled["n1_solv"]) - GATE_N1[0]),
                  abs(float(pooled["n1_solv"]) - GATE_N1[1]))
        need = (sd / gap) ** 2 if gap > 0 and np.isfinite(sd) else float("inf")
        print("VERDICT: NOT confirmed and NOT missed — CL-P1 returned "
              "gate-marginal while CL-P2 and CL-P3 both PASS. The clone sits "
              f"ON the n1 gate floor: pooled {float(pooled['n1_solv']):.4f} vs "
              f"{GATE_N1[0]}, gap {gap:.4f}, per-seed SD {sd:.4f} => "
              f"~{need:.0f} seeds of N = 500 would be needed to separate it "
              "from the edge. That is the measurement, not a resolution "
              "shortfall: more MD at this cell buys a decision only at "
              "absurd cost. Routing (recalibration vs accept-marginal vs "
              "re-site the cell) is a USER GATE — the §6.2 item-1 branch is "
              "pre-registered for a MISS, and this is not one.")
    elif clone_holds:
        print("VERDICT: CL-P1 ∧ CL-P2 ∧ CL-P3 all PASS — an uncapped, "
              "kink-free pure_linear cell reproduces the h405 landing in MD. "
              "This OPENS the form-choice on the §6.4 item-3 "
              "deviation-accounting grounds; it does not settle it (Tier-0's "
              "in-band rejection n_hat = 2.927 is untouched). Nothing is "
              "adopted by this run (§0).")
    else:
        print("VERDICT: the clone does NOT hold as specified -> plan §6.2 "
              "item 1 RECALIBRATION, not a kill and not a new ring: apply "
              "the CL-P6 bias vector to the committed atlas_linsweep.csv and "
              "re-scan the twin bias-corrected (ZERO MD) to name the "
              "corrected clone cells. A second MD iteration is a NEW USER "
              "GATE.")

    if SAVE_CSV_PATH is not None:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, table_rows)}")


if __name__ == "__main__":
    main()
