"""Atlas G4 Step 1 Block 3 — MD finalists scorer report (plan §3.5e).

Scores the finished `gen_tier2atlas_g4finals.py` cells (6 corrected-geometry
N = 1000 runs) on the committed observable vector, prints each next to its
frozen twin row, and evaluates the pre-registered predictions.

Pure scorer — runs nothing, mutates nothing. Physics conventions live in
`i2_helium_md.postprocess.tier2_confirmation`; the per-run observable row is
imported from `tier2atlas_geometry_table.py` (rule 1), which also carries the
pooled-battery scorer-drift oracle this report re-runs first (§1.4).

**Do not read the generator's launch print as n̄.** `n_detect_mean` in the
launch log is `n_detected[detected_mask].mean()`, without the §4r convention
that forces the suppressed (never-opened) class into bin 0. Measured on the
G3 ring the two differ by up to 2.9 He (a037: print 6.97 vs scored 4.10). Only
the numbers in this report are `nbar_det`.

MD-level acceptance (unchanged from §3.5d, bias-free): n1_solv in
[0.19, 0.30] AND nbar_det in [3.77, 4.37].

Pre-registered predictions (roles frozen in the generator before launch;
the deep-KE ordering clause G4F-P4 was added after three cells finished but
before any scored number existed — the launch prints it saw are the
uninformative raw means described above, disclosed here rather than
presented as a pre-launch registration):

* **G4F-P1** — the frozen twin rows reproduce string-exact (generator oracle),
  and the pooled battery reproduces within the scorer-drift tolerance.
* **G4F-P2** — f1 and f2 (the clean ridge optima) land the MD acceptance.
* **G4F-P3** — the N-scaling test: a037 at N = 1000 vs its N = 500 ring value
  (W1 0.826, midHot 0.939, deepKE 0.429, nbar 4.097, n1 0.211). A W1 drop of
  order the per-seed SD (0.095) would mean the ring's W1 spread was largely
  sampling noise.
* **G4F-P4** — deep-KE accessibility: f1/f2/f3 span twin deepKE
  0.691 / 0.776 / 0.887 at nearly equal histogram quality. If MD deepKE
  tracks that ordering the deep-KE axis is parameter-accessible after all;
  if it stays flat near 0.4 the axis is not reachable on the drag surface and
  Block 2 (p_tail) is authorized.
* **G4F-P5** — x345 misses the acceptance on n1 (twin 0.1555 < 0.19).

Invocation:

    python scripts/post_processing/tier2atlas_g4finals_table.py
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

from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.simulation.checkpoint import load_neutral_checkpoint  # noqa: E402
from scripts.gen_tier2atlas_g4finals import (  # noqa: E402
    FINAL_MATRIX,
    GATE_N1,
    GATE_NBAR,
    TWIN_ROWS,
    finals_run_dir_name,
    verify_twin_preregistration,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    ORACLE_RUN,
    RUNS_ROOT,
    check_oracle,
    observable_row,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)
from scripts.post_processing.tier2atlas_g4_transfer import (  # noqa: E402
    PROVISIONAL_NORM,
    joint_score,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_g4finals_table.csv"

# a037 as measured in the G3 ring (N = 500, seed 20260728) — the G4F-P3
# comparison point, read from the committed ring table.
RING_A037 = {"w1_solv": 0.826, "midHot": 0.939, "deepKE": 0.429,
             "nbar_det": 4.097, "n1_solv": 0.211}
RING_B031 = {"w1_solv": 1.061, "midHot": 0.499, "deepKE": 0.715,
             "nbar_det": 4.314, "n1_solv": 0.203}

# The incumbent, on the Block-0-licensed axes (W1 + midHot), measured.
INCUMBENT = {"w1_solv": 0.5791, "midHot": 1.0108, "deepKE": 0.6327}
LICENSED_AXES = ("w1", "midhot")
W1_SEED_SD = 0.0954     # per-seed SD at N = 1000 (Block 0)

# Twin column order frozen by the generator.
TWIN_COLS = ("trapped_frac", "suppressed_frac", "nbar_det", "pred_md_nbar",
             "n1_solv", "w1_solv", "midhot_geo", "deepke", "gate")

# ---------------------------------------------------------------------------


def sampled_geometry_columns(run_dir: Path) -> dict[str, Any]:
    """Measured sampled-geometry columns (per-molecule sizes at this geometry)."""
    ckpt = load_neutral_checkpoint(run_dir / "neutral.npz")
    radii = np.asarray(ckpt.droplet_radii, dtype=float)
    r0 = np.asarray(ckpt.r0, dtype=float)
    if radii.size == 2 * r0.size:        # per-fragment radii, per-molecule r0
        r0 = np.concatenate([r0, r0])    # exact for the mean/quantile columns
    return {
        "R_q50": float(np.quantile(radii, 0.50)),
        "depth_mean_A": float((radii - r0).mean()),
    }


def _s(row: dict[str, Any]) -> float:
    """Joint score on the licensed axes, provisional norms (plan §3.5e)."""
    return joint_score(float(row["w1_solv"]), float(row["midHot"]),
                       float(row["deepKE"]), PROVISIONAL_NORM,
                       axes=LICENSED_AXES)


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # §1.4 oracles: frozen twin rows, then the scorer-drift oracle.
    verify_twin_preregistration()
    oracle_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
    drift = check_oracle(oracle_row)
    if drift:
        print("*** SCORER-DRIFT ORACLE FAILED — do not read the finalists ***")
        for line in drift:
            print(f"    {line}")
        return
    print("scorer-drift oracle OK (pooled battery reproduces).\n")

    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    for cell in FINAL_MATRIX:
        run_dir = RUNS_ROOT / finals_run_dir_name(cell.label)
        if not (run_dir / "detection.npz").exists():
            missing.append(cell.label)
            continue
        row = observable_row(cell.label, run_dir, abundance_ref, ked_ref,
                             with_geometry=False)
        row.update(sampled_geometry_columns(run_dir))
        twin = dict(zip(TWIN_COLS, TWIN_ROWS[cell.label]))
        row.update({
            "v_c": cell.v_c, "well": cell.eb_tag, "tau": cell.tau_ps,
            "E0": cell.e0_eV, "role": cell.role,
            "twin_nbar": float(twin["nbar_det"]),
            "pred_nbar": float(twin["pred_md_nbar"]),
            "twin_n1": float(twin["n1_solv"]),
            "twin_w1": float(twin["w1_solv"]),
            "twin_midHot": float(twin["midhot_geo"]),
            "twin_deepKE": float(twin["deepke"]),
        })
        row["d_pred_nbar"] = float(row["nbar_det"]) - row["pred_nbar"]
        row["md_gate"] = int(
            GATE_N1[0] <= float(row["n1_solv"]) <= GATE_N1[1]
            and GATE_NBAR[0] <= float(row["nbar_det"]) <= GATE_NBAR[1]
        )
        row["S"] = _s(row)
        rows.append(row)

    print("=== G4 Block 3 MD finalists (corrected geometry, N = 1000, "
          "seed 20260729) ===")
    print(format_table(rows) if rows else "(no finished cells)")
    if missing:
        print(f"\nnot yet run: {missing}")
    if not rows:
        return
    by = {r["label"]: r for r in rows}

    print("\n=== Pre-registered predictions (plan §3.5e Block 3) ===")

    # G4F-P2 — the clean ridge optima land.
    landed = [lbl for lbl in ("f1", "f2") if by.get(lbl, {}).get("md_gate")]
    if {"f1", "f2"} <= set(by):
        print(f"G4F-P2 (clean ridge optima land): {landed or 'NONE'} -> "
              f"{'CONFIRMED' if len(landed) == 2 else 'PARTIAL' if landed else 'REFUTED'}")

    # G4F-P3 — the N-scaling test on the bridge controls.
    for lbl, ring in (("a037", RING_A037), ("b031", RING_B031)):
        if lbl not in by:
            continue
        r = by[lbl]
        dw1 = float(r["w1_solv"]) - ring["w1_solv"]
        print(f"G4F-P3 ({lbl} N=500 -> N=1000): "
              f"W1 {ring['w1_solv']:.3f} -> {float(r['w1_solv']):.3f} "
              f"(Δ {dw1:+.3f}, per-seed SD {W1_SEED_SD:.3f}); "
              f"nbar {ring['nbar_det']:.2f} -> {float(r['nbar_det']):.2f}; "
              f"n1 {ring['n1_solv']:.3f} -> {float(r['n1_solv']):.3f}; "
              f"midHot {ring['midHot']:.3f} -> {float(r['midHot']):.3f}; "
              f"deepKE {ring['deepKE']:.3f} -> {float(r['deepKE']):.3f}")
    if "a037" in by:
        dw1 = float(by["a037"]["w1_solv"]) - RING_A037["w1_solv"]
        verdict = ("largely SAMPLING NOISE" if dw1 <= -W1_SEED_SD
                   else "NOT explained by sample size" if abs(dw1) < W1_SEED_SD
                   else "moved the wrong way")
        print(f"   => the §3.5e W1 confound reads: {verdict} "
              f"(incumbent pooled W1 {INCUMBENT['w1_solv']:.3f})")

    # G4F-P4 — is the deep-KE axis parameter-accessible?
    trio = [by[l] for l in ("f2", "f1", "f3") if l in by]
    if len(trio) == 3:
        tw = [float(r["twin_deepKE"]) for r in trio]
        md = [float(r["deepKE"]) for r in trio]
        monotone = all(b >= a for a, b in zip(md, md[1:]))
        span = max(md) - min(md)
        print(f"G4F-P4 (deep-KE accessibility): twin "
              f"{tw[0]:.3f}/{tw[1]:.3f}/{tw[2]:.3f} (f2/f1/f3) -> MD "
              f"{md[0]:.3f}/{md[1]:.3f}/{md[2]:.3f}; ordering "
              f"{'PRESERVED' if monotone else 'NOT preserved'}, MD span "
              f"{span:.3f}")
        if span < 0.10:
            print("   => the deep-KE axis is NOT reachable on the (v_c, tau, "
                  "E0) surface at fixed form: BLOCK 2 (p_tail) IS AUTHORIZED.")
        else:
            print("   => the deep-KE axis moves with the drag knobs; Block 2 "
                  "(p_tail) is not required by this evidence.")

    # G4F-P5 — the fail control misses.
    if "x345" in by:
        r = by["x345"]
        on_n1 = not (GATE_N1[0] <= float(r["n1_solv"]) <= GATE_N1[1])
        print(f"G4F-P5 (fail control x345): n1 {float(r['n1_solv']):.3f} "
              f"(predicted < {GATE_N1[0]}), nbar {float(r['nbar_det']):.2f}, "
              f"gate {r['md_gate']} -> "
              f"{'CONFIRMED (missed on n1)' if on_n1 and not r['md_gate'] else 'CHECK'}")

    # Successor-point input: the joint score against the incumbent.
    inc_s = joint_score(INCUMBENT["w1_solv"], INCUMBENT["midHot"],
                        INCUMBENT["deepKE"], PROVISIONAL_NORM,
                        axes=LICENSED_AXES)
    print(f"\n=== Joint score S on the licensed axes (lower better; "
          f"incumbent finc1v725 = {inc_s:.3f}) ===")
    for r in sorted(rows, key=lambda r: r["S"]):
        flag = "BEATS INCUMBENT" if r["S"] < inc_s else ""
        print(f"  {r['label']:6s} S = {r['S']:.3f}  "
              f"(W1 {float(r['w1_solv']):.3f}, midHot {float(r['midHot']):.3f}, "
              f"deepKE {float(r['deepKE']):.3f}, gate {r['md_gate']})  {flag}")
    beaters = [r["label"] for r in rows if r["S"] < inc_s and r["md_gate"]]
    print(f"\n  gated cells beating the incumbent: {beaters or 'NONE'}")
    print("  (G4 successor-point adjudication is the user's; a pooled "
          "5 × N = 1000 battery at the winner is the verification step.)")

    if SAVE_CSV_PATH is not None:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")


if __name__ == "__main__":
    main()
