"""Atlas G4 Step 2 Block D — the W₁ residual anatomy (plan §3.5f).

Signed per-bin CDF-gap decomposition of the solvated-branch W₁ — exact,
since on unit rung spacing ``W₁ = Σ_n |F_sim(n) − F_ref(n)|`` — for:

1. **pooled h405** (the five Block-V battery members, §4cc in-memory pool),
2. **pooled `finc1v725`** (the standing N = 5000 container — what the 0.579
   landing looked like bin-wise at the wrong geometry),
3. **the τ 4.4 E₀ arm** (committed Block-3 runs h405/h410/h415,
   CRN-paired): the finite-difference **bin signature of the E₀ lever** —
   the measured template a knob signature looks like, and the direct read
   of which bins are E₀-inaccessible (those carry the floor).

Pure scorer — runs nothing, mutates nothing. Sign convention:
``gap(n) = F_sim(n) − F_ref(n)`` on the RQ8 solvated branch (n ≥ 1
renormalized) — ``> 0`` means the sim CDF leads (sim mass sits at smaller
n than the reference up to this bin).

Order of operations (§1.4 — oracles first, hard-fail before any new
number): the scorer-drift oracle (standing pooled battery); the
**identity oracle** ``solvated_gap_profile(read).w1 == w1_solvated`` to
1e-12 for every distribution scored here; the **committed-row oracle**
(rescored W₁ of g4fh405/h410/h415 equals `atlas_g4finals_table.csv` to 4
decimals).

**This scorer reports data, not verdicts.** Block F (the mechanism
fingerprints) is frozen in `TIER2_SENSITIVITY_ATLAS_FINDINGS.md`; the
§3.5f match rule is applied there, in the findings, by the session that
reads this output.

Invocation:

    python scripts/post_processing/tier2atlas_g4step2_w1_anatomy.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from types import SimpleNamespace
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
from i2_helium_md.postprocess.distribution_compare import (  # noqa: E402
    CdfGapProfile,
    cdf_gap_profile,
)
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    ConfirmationDetectionRead,
    load_confirmation_run,
    pool_confirmation_reads,
    score_histogram_vs_reference,
    solvated_renormalized,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

RUNS_ROOT = PROJECT_ROOT / "data" / "runs"
SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_g4step2_w1_anatomy.csv"
FINALS_TABLE_CSV = RUNS_ROOT / "h2b_forward_model" / "atlas_g4finals_table.csv"
ABUNDANCE_REFERENCE_CSV = (
    PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"
)

TOP_SHARE = 0.70          # the §3.5f match-rule share
E0_STEP_EV = 0.005        # the arm's rung spacing (E₀-lever normalization)

# The τ 4.4 E₀ arm, ascending E₀ (committed Block-3 runs, CRN-paired).
ARM_LABELS = ("h405", "h410", "h415")
ARM_E0 = {"h405": 0.405, "h410": 0.41, "h415": 0.415}

IDENTITY_TOL = 1e-12

# ---------------------------------------------------------------------------
# helpers (unit-tested in tests/test_tier2atlas_g4step2_anatomy.py)
# ---------------------------------------------------------------------------


def solvated_gap_profile(read: ConfirmationDetectionRead,
                         abundance_ref) -> CdfGapProfile:
    """Gap profile on the RQ8 solvated branch (n >= 1 renormalized) —
    EXACTLY the branch ``w1_solvated`` is scored on (rule 1: the profile
    delegates to the same ``cdf_gap_profile`` the committed W1 uses)."""
    sim_n, sim_f = solvated_renormalized(read.n_values, read.fraction)
    ref_n, ref_f = solvated_renormalized(abundance_ref.n,
                                         abundance_ref.ion_fraction)
    return cdf_gap_profile(
        SimpleNamespace(n_values=sim_n, fraction=sim_f),
        SimpleNamespace(n=ref_n, ion_fraction=ref_f),
    )


def top_share_bins(profile: CdfGapProfile, share: float = TOP_SHARE
                   ) -> tuple[np.ndarray, float]:
    """Minimal bin set carrying >= ``share`` of W1.

    Bins are ranked by ``(-|gap|, n)`` (largest contribution first, ties to
    the lower n) and taken until the cumulative |gap| fraction first
    reaches ``share``. Returns ``(bins ascending, exact share reached)``.

    Raises
    ------
    ValueError
        Unless ``0 < share <= 1``.
    """
    if not 0.0 < share <= 1.0:
        raise ValueError(f"share must be in (0, 1]; got {share}.")
    absg = np.abs(profile.gaps)
    total = float(absg.sum())
    order = sorted(range(absg.size), key=lambda i: (-absg[i], profile.support[i]))
    picked: list[int] = []
    acc = 0.0
    for i in order:
        if absg[i] == 0.0:
            break
        picked.append(i)
        acc += float(absg[i])
        if acc / total >= share:
            break
    bins = np.sort(profile.support[picked])
    return bins, acc / total


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------


def _identity_oracle(label: str, read: ConfirmationDetectionRead,
                     abundance_ref) -> CdfGapProfile:
    """The per-distribution identity oracle: profile.w1 == w1_solvated."""
    prof = solvated_gap_profile(read, abundance_ref)
    committed = score_histogram_vs_reference(read, abundance_ref).w1_solvated
    if abs(prof.w1 - committed) > IDENTITY_TOL:
        raise AssertionError(
            f"[{label}] identity oracle FAILED: profile W1 {prof.w1!r} != "
            f"committed w1_solvated {committed!r}."
        )
    return prof


def _profile_rows(label: str, prof: CdfGapProfile) -> list[dict[str, Any]]:
    top, reached = top_share_bins(prof, TOP_SHARE)
    absg = np.abs(prof.gaps)
    total = float(absg.sum())
    order = np.argsort(-absg, kind="stable")
    cum = np.zeros_like(absg)
    acc = 0.0
    rank = {}
    for i in order:
        acc += float(absg[i])
        cum[i] = acc / total
        rank[int(prof.support[i])] = cum[i]
    rows = []
    for i, n in enumerate(prof.support):
        rows.append({
            "dist": label,
            "n": int(n),
            "f_sim": f"{prof.f_sim[i]:.6f}",
            "f_ref": f"{prof.f_ref[i]:.6f}",
            "gap": f"{prof.gaps[i]:+.6f}",
            "abs_gap": f"{absg[i]:.6f}",
            "share": f"{absg[i] / total:.4f}" if total else "nan",
            "in_top70": int(n in set(int(b) for b in top)),
        })
    return rows


def _print_profile(label: str, prof: CdfGapProfile) -> None:
    top, reached = top_share_bins(prof, TOP_SHARE)
    print(f"\n--- {label}: W₁ = {prof.w1:.4f}; top-{TOP_SHARE:.0%} bins "
          f"{list(int(b) for b in top)} (share reached {reached:.1%}) ---")
    absg = np.abs(prof.gaps)
    total = float(absg.sum())
    print("  n   f_sim     f_ref     gap        share")
    for i, n in enumerate(prof.support):
        if absg[i] / total < 0.005 and prof.f_sim[i] == 0 and prof.f_ref[i] == 0:
            continue
        mark = " *" if int(n) in set(int(b) for b in top) else ""
        print(f"  {int(n):2d}  {prof.f_sim[i]:.4f}    {prof.f_ref[i]:.4f}    "
              f"{prof.gaps[i]:+.4f}    {absg[i] / total:5.1%}{mark}")


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)

    # §1.4 oracles, part 1 — the scorer-drift oracle via the geometry table's
    # committed machinery (imported lazily to keep this module import-light
    # for the unit tests).
    from scripts.post_processing.tier2atlas_geometry_table import (
        IHE_KED_REFERENCE_CSV, ORACLE_RUN, check_oracle, observable_row,
    )
    from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    oracle_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
    drift = check_oracle(oracle_row)
    if drift:
        print("*** SCORER-DRIFT ORACLE FAILED — do not read the anatomy ***")
        for line in drift:
            print(f"    {line}")
        return
    print("scorer-drift oracle OK (pooled standing battery reproduces).")

    # §1.4 oracles, part 2 — committed-row oracle on the arm cells.
    from scripts.gen_tier2atlas_g4finals import finals_run_dir_name
    committed_w1: dict[str, str] = {}
    with open(FINALS_TABLE_CSV, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["label"] in ARM_LABELS:
                committed_w1[row["label"]] = row["w1_solv"]
    arm_reads: dict[str, ConfirmationDetectionRead] = {}
    for label in ARM_LABELS:
        read = load_confirmation_run(
            RUNS_ROOT / finals_run_dir_name(label), label=label)
        arm_reads[label] = read
        rescored = score_histogram_vs_reference(read, abundance_ref).w1_solvated
        if f"{rescored:.4f}" != f"{float(committed_w1[label]):.4f}":
            raise AssertionError(
                f"[{label}] committed-row oracle FAILED: rescored W₁ "
                f"{rescored:.4f} != committed {float(committed_w1[label]):.4f}."
            )
    print("committed-row oracle OK (arm W₁ h405/h410/h415 rescored to 4 "
          "decimals).")

    # Distributions to profile. Pooled h405 needs the Block-V battery.
    from scripts.gen_tier2atlas_g4step2_battery import (
        MEMBER_SEEDS, member_run_dir_name,
    )
    member_reads = []
    missing = []
    for member in MEMBER_SEEDS:
        run_dir = RUNS_ROOT / member_run_dir_name(member)
        if not (run_dir / "detection.npz").exists():
            missing.append(member)
            continue
        member_reads.append(load_confirmation_run(run_dir, label=member))
    if missing:
        print(f"battery members not yet run: {missing} — stopping before any "
              "anatomy number.")
        return

    profiles: dict[str, CdfGapProfile] = {}
    pooled_read = pool_confirmation_reads(member_reads, label="h405pooled")
    profiles["h405pooled"] = _identity_oracle("h405pooled", pooled_read,
                                              abundance_ref)
    standing_read = load_confirmation_run(RUNS_ROOT / ORACLE_RUN,
                                          label="finc1v725pooled")
    profiles["finc1v725pooled"] = _identity_oracle(
        "finc1v725pooled", standing_read, abundance_ref)
    for label in ARM_LABELS:
        profiles[label] = _identity_oracle(label, arm_reads[label],
                                           abundance_ref)
    print(f"identity oracle OK ({len(profiles)} distributions: "
          f"Σ|gap| == w1_solvated to {IDENTITY_TOL}).")

    print("\n=== G4 Step 2 Block D — W₁ residual anatomy "
          "(gap = F_sim − F_ref, solvated branch) ===")
    _print_profile("h405pooled (corrected geometry, 5 × N = 1000)",
                   profiles["h405pooled"])
    _print_profile("finc1v725pooled (standing, WRONG geometry, N = 5000)",
                   profiles["finc1v725pooled"])
    for label in ARM_LABELS:
        _print_profile(f"{label} (arm, E₀ {ARM_E0[label]})", profiles[label])

    # The E₀-lever bin signature: per-bin Δgap per +0.005 eV, both arm steps.
    print("\n=== E₀-lever bin signature (Δgap per +0.005 eV, CRN-paired) ===")
    print("  n   h405→h410   h410→h415   |  pooled-h405 gap (the residual)")
    base = profiles["h405pooled"]
    for pair in (("h405", "h410"), ("h410", "h415")):
        a, b = profiles[pair[0]], profiles[pair[1]]
        if not np.array_equal(a.support, b.support):
            raise AssertionError(f"arm supports differ: {pair}")
    a, b, c = profiles["h405"], profiles["h410"], profiles["h415"]
    for i, n in enumerate(a.support):
        d1 = (b.gaps[i] - a.gaps[i]) / (E0_STEP_EV / 0.005)
        d2 = (c.gaps[i] - b.gaps[i]) / (E0_STEP_EV / 0.005)
        j = np.nonzero(base.support == n)[0]
        resid = f"{base.gaps[j[0]]:+.4f}" if j.size else "   -  "
        if abs(a.gaps[i]) < 5e-4 and abs(d1) < 5e-4 and abs(d2) < 5e-4:
            continue
        print(f"  {int(n):2d}  {d1:+.4f}     {d2:+.4f}     |  {resid}")

    print("\nBlock F is frozen in the findings doc; apply the §3.5f match "
          "rule there. This scorer reports data, not verdicts.")

    if SAVE_CSV_PATH is not None:
        rows: list[dict[str, Any]] = []
        for label, prof in profiles.items():
            rows.extend(_profile_rows(label, prof))
        SAVE_CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
        with open(SAVE_CSV_PATH, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nCSV -> {SAVE_CSV_PATH.resolve()}")


if __name__ == "__main__":
    main()
