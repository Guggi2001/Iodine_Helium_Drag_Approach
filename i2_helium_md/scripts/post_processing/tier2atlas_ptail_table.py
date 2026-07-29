"""Atlas §3.5h — the p_tail ring scorer (plan §3.5h, pre-registered).

Scores the finished `gen_tier2atlas_ptail_ring.py` cells (pt15/pt20/pt30,
N = 1000, seed 20260729) against the committed h405 finals row (the
CRN-paired p_tail = −1 baseline) and evaluates the pre-registered
predictions PT-P1..PT-P5, on the §3.5g conventions (KE₁ vs the 1.00 eV
peak anchor, KE₂ vs the 0.706 eV reference mean).

Pure scorer — runs nothing, mutates nothing. Physics conventions are
imported from the committed modules (rule 1): reads via
`load_confirmation_run`, observables via `observable_columns`, KE columns
via the Block-V `lowke_columns`, row layout via the §3.5g scan.

Order of operations (§1.4 — oracles first, hard-fail before any new
number):

1. **O1 scorer-drift** — the standing pooled battery reproduces its
   recorded observable row (`check_oracle`).
2. **O2 committed-row** — the committed Block-3 `g4fh405` run rescored
   through this script's own code path equals the committed
   `atlas_g4finals_table.csv` h405 row to 4 decimals.

Invocation:

    python scripts/post_processing/tier2atlas_ptail_table.py
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
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    load_confirmation_run,
)
from scripts.gen_tier2atlas_g4finals import finals_run_dir_name  # noqa: E402
from scripts.gen_tier2atlas_ptail_ring import (  # noqa: E402
    RING,
    ptail_run_dir_name,
)
from scripts.post_processing.tier2atlas_g4step2_battery_table import (  # noqa: E402
    ORACLE_ROW_COLS,
    committed_h405_row,
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
from scripts.post_processing.tier2atlas_ke_lown_scan import (  # noqa: E402
    GATE_N1,
    GATE_NBAR,
    KE1_ANCHOR_EV,
    KE2_REF_EV,
    lowke_counts,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_ptail_ring.csv"

MIDHOT_CEILING = 1.15          # PT-P3 kill criterion band edge
KE1_SUCCESS = (0.95, 1.16)     # plan §3.5h success shape on KE1
DEEPKE_BASE = 0.50             # h405 battery pooled deepKE (PT-P5 center)
DEEPKE_TOL = 0.05              # PT-P5 band half-width

# PT-P1 placement bands (plan §3.5h, placement ± 0.10 eV).
PLACEMENT_BAND = {"pt15": (0.79, 0.99), "pt20": (0.96, 1.16),
                  "pt30": (1.14, 1.34)}

# The CRN-paired baseline row (committed finals h405, p_tail = -1).
H405_BASELINE = {"label": "h405", "p_tail": -1.0}

# ---------------------------------------------------------------------------


def scored_row(label: str, p_tail: float, run_dir: Path,
               abundance_ref, ked_ref) -> dict[str, Any]:
    """One ring row on the §3.5g column conventions."""
    read = load_confirmation_run(run_dir, label=label)
    obs = observable_columns(read, abundance_ref, ked_ref)
    ke = lowke_columns(read)
    row: dict[str, Any] = {"label": label, "p_tail": p_tail}
    row.update({
        "num_scored": obs["num_scored"],
        "trap": obs["trap"], "supp": obs["supp"],
        "nbar_det": obs["nbar_det"], "n1_solv": obs["n1_solv"],
        "w1_solv": obs["w1_solv"], "midHot": obs["midHot"],
        "deepKE": obs["deepKE"], "chi2_med": obs["chi2_med"],
    })
    row.update(ke)
    row.update(lowke_counts(read))
    row["KE1_r100"] = float(ke["KE1_mean"]) / KE1_ANCHOR_EV
    row["KE2_r"] = float(ke["KE2_mean"]) / KE2_REF_EV
    row["md_gate"] = int(
        GATE_N1[0] <= float(obs["n1_solv"]) <= GATE_N1[1]
        and GATE_NBAR[0] <= float(obs["nbar_det"]) <= GATE_NBAR[1]
    )
    return row


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # ---- O1: scorer drift.
    oracle_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
    drift = check_oracle(oracle_row)
    if drift:
        print("*** O1 SCORER-DRIFT ORACLE FAILED — do not read the ring ***")
        for line in drift:
            print(f"    {line}")
        return
    print("O1 scorer-drift oracle OK (pooled standing battery reproduces).")

    # ---- O2: committed finals h405 row through this code path.
    committed = committed_h405_row()
    baseline = scored_row("h405", -1.0,
                          RUNS_ROOT / finals_run_dir_name("h405"),
                          abundance_ref, ked_ref)
    bad = [
        f"{c}: rescored {float(baseline[c]):.4f} vs committed "
        f"{float(committed[c]):.4f}"
        for c in ORACLE_ROW_COLS
        if f"{float(baseline[c]):.4f}" != f"{float(committed[c]):.4f}"
    ]
    if bad:
        print("*** O2 COMMITTED-ROW ORACLE FAILED — do not read the ring ***")
        for line in bad:
            print(f"    {line}")
        return
    print("O2 committed-row oracle OK (finals h405 rescored to 4 decimals).\n")

    # ---- The ring (new numbers — read only after O1/O2).
    rows: list[dict[str, Any]] = [baseline]
    missing: list[str] = []
    for cell in RING:
        run_dir = RUNS_ROOT / ptail_run_dir_name(cell.label)
        if not (run_dir / "detection.npz").exists():
            missing.append(cell.label)
            continue
        rows.append(scored_row(cell.label, cell.p_tail, run_dir,
                               abundance_ref, ked_ref))
    if missing:
        print(f"not yet run: {missing} — ring incomplete, stopping before "
              "any verdict.")
        return

    print("=== §3.5h p_tail ring (h405 pins, seed 20260729, CRN-paired "
          "baseline first) ===")
    print(format_table(rows))

    ring_rows = rows[1:]

    print("\n=== Pre-registered PT verdicts (plan §3.5h) ===")
    # PT-P1: monotone KE1 + placement bands.
    ke1 = [float(r["KE1_mean"]) for r in ring_rows]
    monotone = all(a < b for a, b in zip(
        [float(baseline["KE1_mean"])] + ke1[:-1], ke1))
    in_band = {
        r["label"]: PLACEMENT_BAND[r["label"]][0]
        <= float(r["KE1_mean"]) <= PLACEMENT_BAND[r["label"]][1]
        for r in ring_rows
    }
    print(f"PT-P1: KE1 monotone in |p_tail| -> {monotone}; placement bands "
          + ", ".join(f"{lb} {'IN' if ok else 'OUT'} "
                      f"({PLACEMENT_BAND[lb][0]:.2f}-{PLACEMENT_BAND[lb][1]:.2f})"
                      for lb, ok in in_band.items())
          + f" => {'CONFIRMED' if monotone and all(in_band.values()) else 'REFUTED'}")

    # PT-P2: KE2 rises with the same ordering; KE3 via chi2 is not scored —
    # KE2 only (the pre-registered clause).
    ke2 = [float(r["KE2_mean"]) for r in ring_rows]
    p2 = all(a < b for a, b in zip(
        [float(baseline["KE2_mean"])] + ke2[:-1], ke2))
    print(f"PT-P2: KE2 rises with |p_tail| "
          f"({float(baseline['KE2_mean']):.3f} -> "
          + " -> ".join(f"{v:.3f}" for v in ke2)
          + f", ref {KE2_REF_EV}) => {'CONFIRMED' if p2 else 'REFUTED'}")

    # PT-P3: kill criterion.
    reaching = [r for r in ring_rows if float(r["KE1_mean"]) >= KE1_SUCCESS[0]]
    killed = bool(reaching) and all(
        float(r["midHot"]) > MIDHOT_CEILING for r in reaching)
    if not reaching:
        print("PT-P3: no cell reaches KE1 >= 0.95 — kill criterion not "
              "evaluable on this ring.")
    else:
        print(f"PT-P3 (kill): cells with KE1 >= {KE1_SUCCESS[0]}: "
              + ", ".join(f"{r['label']} midHot {float(r['midHot']):.3f}"
                          for r in reaching)
              + f" => {'KILLED (every reaching cell exceeds midHot 1.15)' if killed else 'NOT killed'}")

    # PT-P4: back-reaction signs + the conditional recenter clause.
    n1 = [float(r["n1_solv"]) for r in ring_rows]
    nbar = [float(r["nbar_det"]) for r in ring_rows]
    p4_sign = (all(a >= b for a, b in zip(
        [float(baseline["n1_solv"])] + n1[:-1], n1))
        and all(a <= b for a, b in zip(
            [float(baseline["nbar_det"])] + nbar[:-1], nbar)))
    print(f"PT-P4: n1 falls / nbar rises with |p_tail| => "
          f"{'CONFIRMED' if p4_sign else 'REFUTED'} "
          f"(n1 {float(baseline['n1_solv']):.3f} -> "
          + " -> ".join(f"{v:.3f}" for v in n1)
          + f"; nbar {float(baseline['nbar_det']):.3f} -> "
          + " -> ".join(f"{v:.3f}" for v in nbar) + ")")
    recenter = [r["label"] for r in ring_rows
                if KE1_SUCCESS[0] <= float(r["KE1_mean"]) <= KE1_SUCCESS[1]
                and float(r["n1_solv"]) < GATE_N1[0]]
    print(f"PT-P4 recenter clause: "
          + (f"FIRES for {recenter} (E0-recenter cell registered)"
             if recenter else "does not fire"))

    # PT-P5: deepKE orthogonality.
    p5 = all(abs(float(r["deepKE"]) - DEEPKE_BASE) <= DEEPKE_TOL
             for r in ring_rows)
    print(f"PT-P5: deepKE within {DEEPKE_BASE} ± {DEEPKE_TOL}: "
          + ", ".join(f"{r['label']} {float(r['deepKE']):.3f}"
                      for r in ring_rows)
          + f" => {'CONFIRMED' if p5 else 'REFUTED'}")

    # Success shape (plan §3.5h): KE1 in band, midHot <= 1.15, gate intact.
    winners = [r["label"] for r in ring_rows
               if KE1_SUCCESS[0] <= float(r["KE1_mean"]) <= KE1_SUCCESS[1]
               and float(r["midHot"]) <= MIDHOT_CEILING
               and int(r["md_gate"]) == 1]
    print(f"\nSuccess-shape cells (KE1 in {KE1_SUCCESS}, midHot <= "
          f"{MIDHOT_CEILING}, gated): {winners if winners else 'NONE'}"
          + (" -> pooled-battery candidate (GV precedent); nothing adopts here."
             if winners else ""))

    if SAVE_CSV_PATH is not None:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")


if __name__ == "__main__":
    main()
