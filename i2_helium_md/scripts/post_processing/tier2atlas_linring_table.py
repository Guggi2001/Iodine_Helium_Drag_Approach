"""Free-form linear MD ring scorer (plan §6.1 bands, §4.1 kill surface,
§6.2 outcome semantics).

Scores the finished `gen_tier2atlas_linring.py` cells (7 `pure_linear`
N = 500 runs + the CRN h405 capped partner) on the committed observable
vector + the low-n KE columns, under BOTH retained-policy bracket ends
(plan §6.2 item 3):

* **Arm A** — `exclude_all_coupled` headline (marginals excluded);
* **Arm B** — marginals injected at handover n with conservative
  asymptotic KE (the §3.5b construction, reused from
  `tier2atlas_retained_bracket.py`).

A cell whose hard-gate verdict differs between the arms is
**policy-blocked** (neither passed nor killed — escalates to the user).

Hard gate (MD-side realization of the §3.5c gate — the g3ring/ke_lown
precedent GATE bands; §6.1's quoted "n̄ ∈ [4.4, 7.1]" is the *twin*-side
band = target 4.07 + the twin-hot bias bracket, and applying it MD-side
would fail the h405 partner itself): n1_solv in [0.19, 0.30] AND
nbar_det in [3.77, 4.37].

Verdict semantics (plan §4.1 + §6.2, pre-registered):

* K-KE kill: CRN-paired ΔKE₁ < +0.05 eV vs h405p at EVERY gated cell.
* Success stamp: gated cell with KE₁ ≥ 0.75 (beats the (A)-ceiling 0.708).
* Survivor (escalation-eligible): gated + not policy-blocked +
  ΔKE₁ ≥ +0.05 (§6.2 item 2 — success is NOT the escalation threshold).
* All lin cells missing the gate = RECALIBRATION, not a kill (§6.2
  item 1): the twin↔MD bias vector below is the product.
* deepKE: [0.6, 1.6] neutral-to-improved; form-questioned ONLY if
  ≥ 2.0 at every lin cell; mixed/1.6–2.0 reported-no-flag (§6.2 item 4).
* Trap is reported decomposed, never gating (§4.1); the CRN fate-flow
  read (who traps under lin that escaped under capped, by stratum) is
  the landing-by-amputation diagnostic — reported, never gating.
* h405p is Δ-only (§6.2 item 5): lin gate verdicts are absolute against
  the frozen bands; a partner band-miss is disclosed as seed noise.

W₁ carries the ±0.13 single-seed disclosure and is nowhere gating.

Pure scorer — runs nothing, mutates nothing. The full observable vector
including every KE column AND the twin forecast columns (§6.2 item 6:
the lin-family authority box / bias vector is a pure CSV read) goes to
the committed `atlas_linring_table.csv`.

Invocation:

    python scripts/post_processing/tier2atlas_linring_table.py
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

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
    read_confirmation_detection,
    score_histogram_vs_reference,
)
from i2_helium_md.simulation.checkpoint import load_ion_checkpoint  # noqa: E402
from i2_helium_md.simulation.detection_stage import (  # noqa: E402
    RETAINED_BOUND_REASON,
    RETAINED_MARGINAL_REASON,
    load_detection_result,
)
from i2_helium_md.simulation.ion_propagation_step import (  # noqa: E402
    ion_state_from_checkpoint_column,
)
from scripts.gen_tier2atlas_linring import (  # noqa: E402
    H405_ANCHOR,
    RING_MATRIX,
    SEED,
    TWIN_ROWS,
    ring_run_dir_name,
    verify_twin_preregistration,
)
from scripts.post_processing.tier2atlas_g3ring_table import (  # noqa: E402
    sampled_geometry_columns,
)
from scripts.post_processing.tier2atlas_g4step2_battery_table import (  # noqa: E402
    lowke_columns,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    ORACLE_RUN,
    RUNS_ROOT,
    check_oracle,
    observable_row,
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
# USER SETTINGS — frozen bands and semantics (plan §6.1/§6.2/§4.1)
# ---------------------------------------------------------------------------

# MD-side realization of the §3.5c hard gate (g3ring/ke_lown precedent; see
# module docstring for why the twin-side [4.4, 7.1] n̄ band does not apply).
GATE_N1 = (0.19, 0.30)
GATE_NBAR = (3.77, 4.37)

KKE_MIN_DELTA_EV = 0.05          # K-KE kill threshold (CRN-paired vs h405p)
SUCCESS_KE1_EV = 0.75            # interpretive stamp: beats (A)-ceiling 0.708
DEEPKE_NEUTRAL = (0.6, 1.6)      # neutral-to-improved band (NB-RQ11-12)
DEEPKE_HOT = 2.0                 # form-questioned only if EVERY lin cell >= 2
W1_SINGLE_SEED_DISCLOSURE = 0.13

# KE-path oracle: the committed incumbent pooled KE row (atlas_ke_lown_scan
# .csv, group=incumbent label=pooled), reproduced at rtol 1e-9.
KE_LOWN_CSV = RUNS_ROOT / "h2b_forward_model" / "atlas_ke_lown_scan.csv"
KE_ORACLE_COLS = ("KE1_mean", "KE1_med", "KE1_mode",
                  "KE2_mean", "KE2_med", "KE2_mode")

SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_linring_table.csv"

_TWIN_JOIN = ("twin_trap", "twin_supp", "twin_nbar", "twin_n1", "twin_w1",
              "twin_KE1", "twin_KE2", "twin_deepKE", "twin_gate")


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------


def ke_width_columns(read) -> dict[str, float]:
    """KE₁ width read (plan §6.1 'the twin predicts nothing here')."""
    sel = read.ke_scored_eV[read.n_scored == 1]
    if sel.size < 5:
        return {"KE1_SD": float("nan"), "above115_n1": float("nan")}
    return {
        "KE1_SD": float(sel.std(ddof=1)),
        "above115_n1": float((sel > 1.15).mean()),
    }


def hard_gate(n1: float, nbar: float) -> int:
    return int(GATE_N1[0] <= n1 <= GATE_N1[1]
               and GATE_NBAR[0] <= nbar <= GATE_NBAR[1])


def arm_b_gate_columns(run_dir: Path, label: str, abundance_ref
                       ) -> dict[str, float]:
    """Arm-B (marginal-injection) gate surface for one cell.

    No marginals => Arm B is definitionally Arm A; the columns then repeat
    the Arm-A values so the policy-flip comparison stays uniform.
    """
    dossier = marginal_dossier(run_dir)
    if dossier is None:
        read = load_confirmation_run(run_dir, label=label)
        n_marg = 0
    else:
        read = read_confirmation_detection(
            arm_b_detection(dossier), label=f"{label}B")
        n_marg = int(dossier["marginal"].sum())
    hist = score_histogram_vs_reference(read, abundance_ref)
    return {
        "n_marg": n_marg,
        "n1_B": float(hist.sim_n1),
        "nbar_B": float(read.n_mean),
        "w1_B": float(hist.w1_solvated),
        "gate_B": hard_gate(float(hist.sim_n1), float(read.n_mean)),
    }


def deepke_class(value: float) -> str:
    if not np.isfinite(value):
        return "nan"
    if value >= DEEPKE_HOT:
        return "hot(>=2.0)"
    if DEEPKE_NEUTRAL[0] <= value <= DEEPKE_NEUTRAL[1]:
        return "neutral[0.6,1.6]"
    if DEEPKE_NEUTRAL[1] < value < DEEPKE_HOT:
        return "strip(1.6,2.0)"
    return "cold(<0.6)"


def fate_flow(cell_dir: Path, partner_dir: Path) -> dict[str, float]:
    """CRN per-ion fate flow vs the h405 partner (§4.1 item 3 — reported,
    never gating). Same master draw => fragment index i is the same ion in
    both runs; 'trapped' = the retained classes (bound + marginal)."""
    det_l = load_detection_result(cell_dir / "detection.npz")
    det_p = load_detection_result(partner_dir / "detection.npz")
    reasons_l = np.asarray(det_l.state_reason)
    reasons_p = np.asarray(det_p.state_reason)
    retained = (RETAINED_BOUND_REASON, RETAINED_MARGINAL_REASON)
    trap_l = np.isin(reasons_l, retained)
    trap_p = np.isin(reasons_p, retained)

    new_trap = trap_l & ~trap_p          # trapped under lin, free under capped
    freed = ~trap_l & trap_p             # the reverse flow

    out: dict[str, float] = {
        "flow_new_trap": int(new_trap.sum()),
        "flow_freed": int(freed.sum()),
    }
    if np.any(new_trap):
        # Stratum of the amputated ions: droplet size, birth depth, v0 —
        # plus what they were under h405 (their detected n), i.e. which part
        # of the histogram the lin landing removed.
        ion_ckpt = load_ion_checkpoint(cell_dir / "ion.npz")
        state0 = ion_state_from_checkpoint_column(ion_ckpt, 0)
        radii = np.asarray(ion_ckpt.droplet_radii_angstrom, dtype=float)
        r0 = np.sqrt(np.asarray(state0.x) ** 2 + np.asarray(state0.y) ** 2
                     + np.asarray(state0.z) ** 2)
        v0 = np.sqrt(np.asarray(state0.vx) ** 2 + np.asarray(state0.vy) ** 2
                     + np.asarray(state0.vz) ** 2)
        n_p = np.asarray(det_p.n_detected, dtype=float)
        out.update({
            "flow_R_med": float(np.median(radii[new_trap])),
            "flow_depth_med": float(np.median((r0 - radii)[new_trap])),
            "flow_v0_med": float(np.median(v0[new_trap])),
            "flow_h405_n_med": float(np.median(n_p[new_trap])),
        })
    else:
        out.update({"flow_R_med": float("nan"), "flow_depth_med": float("nan"),
                    "flow_v0_med": float("nan"),
                    "flow_h405_n_med": float("nan")})
    return out


def ke_path_oracle(abundance_ref, ked_ref) -> list[str]:
    """O3: the pooled-battery KE columns reproduce the committed
    atlas_ke_lown_scan.csv incumbent pooled row at rtol 1e-9."""
    committed = None
    with open(KE_LOWN_CSV, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["group"] == "incumbent" and r["label"] == "pooled":
                committed = r
                break
    if committed is None:
        return ["no incumbent pooled row in the committed KE CSV"]
    read = load_confirmation_run(RUNS_ROOT / ORACLE_RUN, label="oracle")
    got = lowke_columns(read)
    problems = []
    for col in KE_ORACLE_COLS:
        want = float(committed[col])
        if not np.isclose(float(got[col]), want, rtol=1e-9, atol=0.0):
            problems.append(f"{col}: rescored {got[col]!r} vs committed "
                            f"{want!r}")
    return problems


# ---------------------------------------------------------------------------


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # Oracles — all three before any new number is read.
    verify_twin_preregistration()
    oracle_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
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
    print("scorer-drift + KE-path oracles OK (pooled battery reproduces).")
    print()

    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    for cell in RING_MATRIX:
        run_dir = RUNS_ROOT / ring_run_dir_name(cell.label)
        if not (run_dir / "detection.npz").exists():
            missing.append(cell.label)
            continue
        read = load_confirmation_run(run_dir, label=cell.label)
        row = observable_row(cell.label, run_dir, abundance_ref, ked_ref,
                             with_geometry=False)
        row.update(sampled_geometry_columns(run_dir))
        row.update(lowke_columns(read))
        row.update(ke_width_columns(read))
        row.update({"KE1_n": int((read.n_scored == 1).sum()),
                    "KE2_n": int((read.n_scored == 2).sum())})
        row.update(arm_b_gate_columns(run_dir, cell.label, abundance_ref))
        row.update({
            "form": "lin" if cell.a is not None else "capped",
            "a": cell.a if cell.a is not None else "-",
            "v_c": cell.v_c if cell.v_c is not None else "-",
            "well": cell.eb_tag, "tau": cell.tau_ps, "E0": cell.e0_eV,
        })
        # Twin forecast join (§6.2 item 6): frozen strings -> floats.
        if cell.a is not None:
            (t_trap, t_supp, t_nbar, t_n1, t_w1,
             t_ke1, t_ke2, t_deep, t_gate) = TWIN_ROWS[cell.label]
        else:
            _md1, _md2, t_ke1, t_ke2, t_n1, t_nbar, t_w1 = H405_ANCHOR
            t_trap = t_supp = t_deep = "nan"
            t_gate = "1"
        row.update({
            "twin_trap": float(t_trap), "twin_supp": float(t_supp),
            "twin_nbar": float(t_nbar), "twin_n1": float(t_n1),
            "twin_w1": float(t_w1), "twin_KE1": float(t_ke1),
            "twin_KE2": float(t_ke2), "twin_deepKE": float(t_deep),
            "twin_gate": int(t_gate),
            "d_n1_twin": float(row["n1_solv"]) - float(t_n1),
            "d_nbar_twin": float(row["nbar_det"]) - float(t_nbar),
            "d_KE1_twin": float(row["KE1_mean"]) - float(t_ke1),
        })
        row.update({
            "md_gate_A": hard_gate(float(row["n1_solv"]),
                                   float(row["nbar_det"])),
        })
        row["policy_blocked"] = int(row["md_gate_A"] != row["gate_B"])
        rows.append(row)

    print(f"=== Free-form linear MD ring (corrected geometry, N = 500, "
          f"seed {SEED}, CRN) ===")
    print(format_table(rows) if rows else "(no finished cells)")
    if missing:
        print(f"\nnot yet run: {missing}")
    if not rows:
        return
    by = {r["label"]: r for r in rows}

    if "h405p" not in by:
        print("\nh405p partner not finished — no paired reads possible.")
        if SAVE_CSV_PATH is not None:
            print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")
        return

    partner = by["h405p"]
    partner_dir = RUNS_ROOT / ring_run_dir_name("h405p")
    lin_rows = [r for r in rows if r["form"] == "lin"]

    # CRN-paired ΔKE₁ + fate flow, per lin cell.
    for r in lin_rows:
        r["d_KE1_h405"] = float(r["KE1_mean"]) - float(partner["KE1_mean"])
        r["d_KE2_h405"] = float(r["KE2_mean"]) - float(partner["KE2_mean"])
        r.update(fate_flow(RUNS_ROOT / ring_run_dir_name(r["label"]),
                           partner_dir))
        gated = bool(r["md_gate_A"]) and not bool(r["policy_blocked"])
        r["survivor"] = int(gated and r["d_KE1_h405"] >= KKE_MIN_DELTA_EV)
        r["success"] = int(bool(r["md_gate_A"])
                           and float(r["KE1_mean"]) >= SUCCESS_KE1_EV)
    partner["d_KE1_h405"] = 0.0
    partner["d_KE2_h405"] = 0.0
    partner["survivor"] = "-"
    partner["success"] = "-"

    # ---- Verdicts (plan §4.1 + §6.2, pre-registered) ---------------------
    print("\n=== Verdicts (plan §4.1 kill surface + §6.2 semantics) ===")

    # Partner disclosure (§6.2 item 5): Δ-only, gate read is informational.
    p_gate = hard_gate(float(partner["n1_solv"]), float(partner["nbar_det"]))
    print(f"h405p partner (Δ-only): n1={partner['n1_solv']:.3f} "
          f"nbar={partner['nbar_det']:.3f} KE1={partner['KE1_mean']:.4f} "
          f"in-band={bool(p_gate)}"
          + ("" if p_gate else "  [seed-noise disclosure: partner outside a "
             "band; lin verdicts are absolute and unaffected]"))

    gated_cells = [r for r in lin_rows if r["md_gate_A"]]
    blocked = [r["label"] for r in lin_rows if r["policy_blocked"]]
    print(f"hard gate (Arm A, MD-side n1 {list(GATE_N1)} x nbar "
          f"{list(GATE_NBAR)}): gated {[r['label'] for r in gated_cells] or 'NONE'}"
          f"; policy-blocked {blocked or 'none'}")

    if not gated_cells:
        print("ALL lin cells miss the hard gate -> §6.2 item 1: "
              "RECALIBRATION, not a kill. The twin<->MD bias vector "
              "(d_n1_twin / d_nbar_twin / d_KE1_twin columns) is the "
              "product; one bias-corrected twin re-scan is licensed; "
              "family dies only if the corrected re-scan finds no in-grid "
              "basin or the bias vector is incoherent.")
    else:
        clean_gated = [r for r in gated_cells if not r["policy_blocked"]]
        kke_all_below = all(r["d_KE1_h405"] < KKE_MIN_DELTA_EV
                            for r in clean_gated) if clean_gated else False
        for r in gated_cells:
            print(f"  {r['label']}: dKE1_h405={r['d_KE1_h405']:+.4f} eV "
                  f"KE1={r['KE1_mean']:.4f} survivor={r['survivor']} "
                  f"success={r['success']}"
                  + ("  [policy-blocked]" if r["policy_blocked"] else ""))
        if kke_all_below:
            print(f"K-KE KILL FIRES: ΔKE₁ < +{KKE_MIN_DELTA_EV} eV at every "
                  "gated (non-blocked) cell — co-selection wins in real MD; "
                  "the drag-side KE₁ question closes with the full "
                  "instrument chain.")
        else:
            surv = [r["label"] for r in lin_rows if r.get("survivor") == 1]
            succ = [r["label"] for r in lin_rows if r.get("success") == 1]
            print(f"K-KE kill does NOT fire. survivors (escalation-eligible)="
                  f"{surv or 'none'}; success (KE1 >= {SUCCESS_KE1_EV})="
                  f"{succ or 'none'}")
            if succ:
                print("SUCCESS band met: escalation = one N = 1000 x 5-seed "
                      "battery at the single best cell (user gate).")

    # deepKE interpretation band (§6.2 item 4).
    classes = {r["label"]: deepke_class(float(r["deepKE"])) for r in lin_rows}
    print("deepKE band read: " + "; ".join(f"{k}={v}"
                                           for k, v in classes.items()))
    if all(v == "hot(>=2.0)" for v in classes.values()):
        print("deepKE >= 2.0 at EVERY lin cell -> the twin's hot flag is "
              "confirmed: FORM-QUESTIONED, user adjudicates.")
    else:
        print("deepKE: mixed / in-band -> reported-no-flag (§6.2 item 4).")

    # Fate-flow / trap decomposition (reported, never gating).
    print("trap decomposed + CRN fate flow vs h405p (reported, never "
          "gating):")
    for r in lin_rows:
        print(f"  {r['label']}: trap={r['trap']:.3f} "
              f"(bound {r['trap_bound']:.3f} / marg {r['trap_marg']:.3f}); "
              f"new-trapped vs h405p {int(r['flow_new_trap'])} "
              f"(freed {int(r['flow_freed'])}), stratum R_med="
              f"{r['flow_R_med']:.1f} A depth_med={r['flow_depth_med']:.1f} A "
              f"v0_med={r['flow_v0_med']:.2f} A/ps, was n_med="
              f"{r['flow_h405_n_med']:.1f} under h405")

    print(f"W1 disclosure: single-seed +-{W1_SINGLE_SEED_DISCLOSURE}; W1 is "
          "nowhere gating.")

    if SAVE_CSV_PATH is not None and rows:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")


if __name__ == "__main__":
    main()
