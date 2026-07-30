"""Atlas §3.5l — the (C) probe scorer (CP-1..8, registration frozen §8).

Scores the finished `gen_tier2atlas_ce_probe.py` cells
(cfull/aonly/bonly/cq3hi, N = 1000, seed 20260729) against the committed
h405 finals row and evaluates the frozen CP-1..8 bands
(`TIER2_CE_CHANNEL_EXIT_STRIP_DESIGN.md` §8, user-approved 2026-07-29).

Pure scorer — runs nothing, mutates nothing. Physics conventions imported
from the committed modules (rule 1). **Partner mask (OQ-I):** every
channel-sampled cell is scored through `ce_partner_include_mask` — the
emulated q3-partner fragments leave the scored universe entirely. The
committed-row oracle (O2) and the baseline row are unmasked (no partners
exist there), so the delivered code path is exercised unchanged. PC-5
caveat: CRN pairing vs the committed h405 is clean only for aonly (strip
Bernoullis aside); cfull/bonly/cq3hi break per-molecule pairing at the
source (every E_m differs).

Order of operations (§1.4 — oracles first, hard-fail before any new number):

1. **O1 scorer-drift** — the standing pooled battery reproduces its
   recorded observable row (`check_oracle`).
2. **O2 committed-row** — the committed `g4fh405` run rescored through this
   script's code path equals the committed finals h405 row to 4 decimals.

Invocation:

    python scripts/post_processing/tier2atlas_ce_probe_table.py
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
    ce_partner_include_mask,
    load_confirmation_run,
)
from scripts.gen_tier2atlas_ce_probe import (  # noqa: E402
    RING,
    ce_run_dir_name,
)
from scripts.gen_tier2atlas_g4finals import finals_run_dir_name  # noqa: E402
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
    KE1_ANCHOR_EV,
    KE2_REF_EV,
    lowke_counts,
)
from scripts.post_processing.tier2atlas_sn_table import sn_ke_columns  # noqa: E402
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS — the frozen CP-1..8 magnitudes (design §8; do not tune)
# ---------------------------------------------------------------------------

SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_ce_probe.csv"

CP1_SUPP_CEIL = 0.05
CP2_KE1_SD_BAND = (0.30, 0.70)
CP3_KE1_MEAN_BAND = (0.85, 1.20)
CP3_ABOVE115_BAND = (0.35, 0.55)
CP4_N1_BAND = (0.18, 0.26)
CP5_SLOW_BARE_CEIL = 0.01          # ensemble share: n = 0 AND KE < 1.0 eV
CP5_BARE_LT1_SHARE_CEIL = 0.05     # within the bare row: share below 1.0 eV
CP6_MIDHOT_BAND = (0.85, 1.15)
CP7_W1_CEIL = 0.80

ABOVE_EV = 1.15                    # the n = 1 KED upper-share threshold
SLOW_BARE_EV = 1.0
REF_ABOVE115 = 0.487               # experimental n = 1 KED share above 1.15
REF_N1_KED_MODE = 0.891            # single-mode reference (Tension 1 context)


def _in(band: tuple[float, float], x: float) -> bool:
    return band[0] <= x <= band[1]


def ce_extra_columns(read) -> dict[str, Any]:
    """The (C)-specific observables beyond the committed column stack.

    * ``bare0``          — scored fraction in bin 0 (suppressed + physical
      bare; under (C)-on supp → ~0 so this is ~the physical bare row),
    * ``slow_bare``      — ensemble share with n = 0 AND KE < 1.0 eV (CP-5),
    * ``bare_lt1_share`` — within-bare share below 1.0 eV (CP-5 second
      clause; NaN on an empty bare row),
    * ``KE1_above115``   — share of the n = 1 KED above 1.15 eV (CP-3),
    * ``bareKE_mean/med``— bare-row KE stats (CP-8 qualitative context).
    """
    n = np.asarray(read.n_scored)
    ke = np.asarray(read.ke_scored_eV)
    bare = n == 0
    n1 = n == 1
    out: dict[str, Any] = {
        "bare0": float(np.mean(bare)),
        "slow_bare": float(np.mean(bare & (ke < SLOW_BARE_EV))),
        "bare_lt1_share": (
            float(np.mean(ke[bare] < SLOW_BARE_EV)) if bare.any()
            else float("nan")
        ),
        "KE1_above115": (
            float(np.mean(ke[n1] > ABOVE_EV)) if n1.any() else float("nan")
        ),
        "bareKE_mean": float(ke[bare].mean()) if bare.any() else float("nan"),
        "bareKE_med": (
            float(np.median(ke[bare])) if bare.any() else float("nan")
        ),
    }
    return out


def strip_columns(run_dir: Path) -> dict[str, Any]:
    """Strip bookkeeping from the run's v8 checkpoint fields (lazy read)."""
    for name in ("relaxation.npz", "ion.npz"):
        p = run_dir / name
        if not p.exists():
            continue
        with np.load(p, allow_pickle=False) as npz:
            if "ce_strip_count" not in npz.files:
                return {"stripped_frac": float("nan"),
                        "strip_mean": float("nan")}
            counts = np.asarray(npz["ce_strip_count"], dtype=int)
        return {
            "stripped_frac": float(np.mean(counts > 0)),
            "strip_mean": float(counts[counts > 0].mean())
            if np.any(counts > 0) else 0.0,
        }
    raise FileNotFoundError(f"no checkpoint in {run_dir}")


def scored_row(label: str, run_dir: Path, abundance_ref, ked_ref,
               *, masked: bool) -> dict[str, Any]:
    """One probe row on the committed column conventions + the CE columns."""
    include = ce_partner_include_mask(run_dir) if masked else None
    read = load_confirmation_run(run_dir, label=label, include_mask=include)
    obs = observable_columns(read, abundance_ref, ked_ref)
    ke = lowke_columns(read)
    row: dict[str, Any] = {"label": label}
    row.update({
        "num_scored": obs["num_scored"],
        "trap": obs["trap"], "supp": obs["supp"],
        "nbar_det": obs["nbar_det"], "n1_solv": obs["n1_solv"],
        "w1_solv": obs["w1_solv"], "midHot": obs["midHot"],
        "deepKE": obs["deepKE"], "chi2_med": obs["chi2_med"],
    })
    row.update(ke)
    row.update(lowke_counts(read))
    row.update(sn_ke_columns(read))
    row.update(ce_extra_columns(read))
    row.update(strip_columns(run_dir))
    row["KE1_r100"] = float(ke["KE1_mean"]) / KE1_ANCHOR_EV
    row["KE2_r"] = (
        float(ke["KE2_mean"]) / KE2_REF_EV
        if np.isfinite(ke["KE2_mean"]) else float("nan")
    )
    return row


def cp_verdicts(row: dict[str, Any]) -> list[str]:
    """The frozen CP-1..8 verdict lines for one C-cell row."""
    lines = []
    supp = float(row["supp"])
    lines.append(f"CP-1  supp {supp:.3f} <= {CP1_SUPP_CEIL} -> "
                 f"{'PASS' if supp <= CP1_SUPP_CEIL else 'FAIL'}")
    sd = float(row["KE1_sd"])
    lines.append(f"CP-2  KE1_sd {sd:.3f} in {CP2_KE1_SD_BAND} -> "
                 f"{'PASS' if _in(CP2_KE1_SD_BAND, sd) else 'FAIL'}")
    m = float(row["KE1_mean"])
    a115 = float(row["KE1_above115"])
    ok3 = _in(CP3_KE1_MEAN_BAND, m) and _in(CP3_ABOVE115_BAND, a115)
    lines.append(
        f"CP-3  KE1_mean {m:.3f} in {CP3_KE1_MEAN_BAND} AND above-1.15 "
        f"{a115:.3f} in {CP3_ABOVE115_BAND} (ref {REF_ABOVE115}) -> "
        f"{'PASS' if ok3 else 'FAIL'}"
    )
    n1 = float(row["n1_solv"])
    lines.append(f"CP-4  n1_solv {n1:.3f} in {CP4_N1_BAND} -> "
                 f"{'PASS' if _in(CP4_N1_BAND, n1) else 'FAIL'} "
                 "(forecast band, NOT the 0.31 reference — Tension 2)")
    sb = float(row["slow_bare"])
    blt = float(row["bare_lt1_share"])
    ok5 = sb <= CP5_SLOW_BARE_CEIL and (
        not np.isfinite(blt) or blt <= CP5_BARE_LT1_SHARE_CEIL
    )
    lines.append(
        f"CP-5* slow_bare {sb:.4f} <= {CP5_SLOW_BARE_CEIL} AND bare<1eV "
        f"share {blt:.3f} <= {CP5_BARE_LT1_SHARE_CEIL} -> "
        f"{'PASS' if ok5 else 'KILL FIRED'}"
    )
    mh = float(row["midHot"])
    lines.append(f"CP-6* midHot {mh:.3f} in {CP6_MIDHOT_BAND} -> "
                 f"{'PASS' if _in(CP6_MIDHOT_BAND, mh) else 'KILL FIRED'}")
    w1 = float(row["w1_solv"])
    lines.append(f"CP-7* W1_solv {w1:.3f} <= {CP7_W1_CEIL} -> "
                 f"{'PASS' if w1 <= CP7_W1_CEIL else 'KILL FIRED'}")
    lines.append(
        f"CP-8  (soft, qualitative) bare row: n={row['bare0']:.3f} of "
        f"scored, KE mean {row['bareKE_mean']:.2f} / med "
        f"{row['bareKE_med']:.2f} eV — two-lump read is visual "
        "(report the bare KED; Q4 rider acknowledged)"
    )
    return lines


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
        print("*** O1 SCORER-DRIFT ORACLE FAILED — do not read the probe ***")
        for line in drift:
            print(f"    {line}")
        return
    print("O1 scorer-drift oracle OK (pooled standing battery reproduces).")

    # ---- O2: committed finals h405 row through this code path (unmasked —
    # no partners exist pre-(C), so the delivered path is exercised as-is).
    committed = committed_h405_row()
    baseline = scored_row("h405", RUNS_ROOT / finals_run_dir_name("h405"),
                          abundance_ref, ked_ref, masked=False)
    bad = [
        f"{c}: rescored {float(baseline[c]):.4f} vs committed "
        f"{float(committed[c]):.4f}"
        for c in ORACLE_ROW_COLS
        if f"{float(baseline[c]):.4f}" != f"{float(committed[c]):.4f}"
    ]
    if bad:
        print("*** O2 COMMITTED-ROW ORACLE FAILED — do not read the probe ***")
        for line in bad:
            print(f"    {line}")
        return
    print("O2 committed-row oracle OK (finals h405 rescored to 4 decimals).\n")

    # ---- The probe (new numbers — read only after O1/O2).
    rows: list[dict[str, Any]] = [baseline]
    missing: list[str] = []
    by_label: dict[str, dict[str, Any]] = {"h405": baseline}
    for cell in RING:
        run_dir = RUNS_ROOT / ce_run_dir_name(cell.label)
        if not (run_dir / "detection.npz").exists():
            missing.append(cell.label)
            continue
        row = scored_row(cell.label, run_dir, abundance_ref, ked_ref,
                         masked=cell.channels)
        rows.append(row)
        by_label[cell.label] = row
    if missing:
        print(f"not yet run: {missing} — probe incomplete, stopping before "
              "any verdict.")
        return

    print("=== (C) probe (h405 pins, seed 20260729; PC-5: CRN vs h405 clean "
          "for aonly only; partner-masked where channels on) ===")
    print(format_table(rows))

    print("\n=== Frozen CP-1..8 verdicts on C-full (design §8) ===")
    for line in cp_verdicts(by_label["cfull"]):
        print(line)

    print("\n=== Coupling arm (CP-6 at both bracket ends) ===")
    for lb in ("cfull", "cq3hi"):
        r = by_label[lb]
        print(f"{lb}: f_int,Q3 {'0.0985' if lb == 'cfull' else '0.15'} -> "
              f"midHot {float(r['midHot']):.3f}, supp {float(r['supp']):.3f}, "
              f"KE1 {float(r['KE1_mean']):.3f}, n1 {float(r['n1_solv']):.3f}")

    print("\n=== Attribution reads (measured, not gated) ===")
    a, b = by_label["aonly"], by_label["bonly"]
    print(f"A-only (strip alone, scalar budget): supp {float(a['supp']):.3f}, "
          f"n1 {float(a['n1_solv']):.3f}, KE1 {float(a['KE1_mean']):.3f} "
          f"(the §3.5j bracket ceiling was ~0.71–0.77), "
          f"KE1_sd {float(a['KE1_sd']):.3f}")
    print(f"B-only (mixture alone): supp {float(b['supp']):.3f} "
          f"(§3.5k complementarity expects the fast branch parked), "
          f"n1 {float(b['n1_solv']):.3f}, KE1 {float(b['KE1_mean']):.3f}")

    print("\nTension-1 context: n=1 KED mode reference (single-mode) "
          f"{REF_N1_KED_MODE} eV; the P3 blend was two-lump — inspect the "
          "cfull n=1 KED for lump merging (named discriminator, not a gate).")

    if SAVE_CSV_PATH is not None:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")


if __name__ == "__main__":
    main()
