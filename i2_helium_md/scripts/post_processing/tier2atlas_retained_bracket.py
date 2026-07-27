"""Atlas Axis A — the ONE-OFF marginal-class bracket diagnostic (plan §3.5b item 8).

**This script is deliberately temporary.** It exists to answer one question:
does excluding the `droplet_retained_marginal` class bias the scored observable
vector of the anchored-radius geometry cells? Per the pre-registered drop rule,
if the bracket comes out **tight** the result is recorded in the findings doc
and this file is deleted; if it comes out **wide** the offending observable is
flagged and the follow-up is the kinetic forward model. Do not build anything
on top of it.

Why a bracket is needed at all
------------------------------
Scoring a cell under `exclude_all_coupled` cannot report its own bias: the
excluded ions are absent from every observable *by construction*. So the cell
is scored twice —

* **Arm A** — marginals excluded. The headline vector (what the grid table
  prints).
* **Arm B** — marginals injected into the scored ensemble at their handover
  ``n`` and their exact conservative **asymptotic** kinetic energy.

Arm B costs nothing new. The E2 relaxation is zero-gamma
(`relaxation_stage.py:43`), so an ion that clears its effective-potential
barrier arrives at infinity with exactly ``E_tot - U(inf)`` — a quantity
`detection_stage.escape_energetics` already computes while classifying the ion.
Nothing is integrated and no new physics is introduced.

What Arm B is NOT
-----------------
It is the **conservative** corner: it holds ``n`` at its handover value and
lets the ion coast out unchanged. Over the ~80–690 ns these ions actually need,
pickup is live (λ₀ρ) and would raise ``n`` and mass-load them — possibly into
the bound class. Arm B therefore brackets the *cold-arrival* direction, and a
tight Arm A/Arm B gap does not by itself retire the pickup question; the
kinetic forward model does. Stated in the printed output too.

Known benign warning
--------------------
Arm B injects ions whose asymptotic KE is exactly 0, so
``tier2_confirmation``'s within-25 %-agreement counter evaluates ``1/ratio``
at ratio = 0 and numpy emits a divide-by-zero RuntimeWarning. It affects only
that diagnostic count, never the observables the bracket tests, and it cannot
occur on a real detected ensemble (a detector arrival has KE > 0). Left
unsuppressed on purpose: a silenced warning in a one-off script is worse than
a documented one.

Pre-registered reading rule (frozen before this ran, plan §3.5b item 8)
-----------------------------------------------------------------------
* Expected direction: Arm B adds low-KE, mid-to-high-n fragments ⇒ deepKE
  **down**, n̄ slightly up, midHot ≈ flat, supp unchanged.
* The **fate split is excluded from the test**: Arm B reclassifies
  marginal → detected *definitionally*, so that Δ is not bias.
* Tight ⇔ |Δ| < 1 seed-SD on ``nbar_det``, ``n1_solv``, ``midHot``, ``deepKE``,
  ``chi2_med`` **and** |Δw1_solv| < 0.15 (the §3.1 single-seed W₁ floor).
* The seed-SDs are **measured here**, not remembered: the five N = 1000
  standing-battery members are scored and their per-observable sample SD is
  scaled to each cell's realized detected count by √(N_member / N_cell).

Invocation:

    python scripts/post_processing/tier2atlas_retained_bracket.py
"""

from __future__ import annotations

import dataclasses
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
    deep_bin_ke_ratio,
    midhot_ratio,
    read_confirmation_detection,
    score_histogram_vs_reference,
    score_ke_curve_vs_reference,
)
from i2_helium_md.simulation.checkpoint import load_ion_checkpoint  # noqa: E402
from i2_helium_md.simulation.detection_stage import (  # noqa: E402
    RETAINED_MARGINAL_REASON,
    escape_energetics,
    load_detection_result,
)
from i2_helium_md.simulation.ion_propagation_step import (  # noqa: E402
    ion_state_from_checkpoint_column,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    CELL_LABELS,
    IHE_KED_REFERENCE_CSV,
    INCLUDE_SIM_SE,
    KE_N_MIN,
    MIN_KE_BIN_COUNT,
    N1_ANCHOR,
    RUN_SELECTION,
    RUNS_ROOT,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

# The five standing-battery members the seed-SD yardstick is measured from.
BATTERY_MEMBERS: tuple[str, ...] = tuple(
    f"9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s{i}"
    for i in range(1, 6)
)

# Observables the tightness test runs on (the fate split is excluded on
# purpose — Arm B moves it by construction).
TEST_COLUMNS: tuple[str, ...] = (
    "nbar_det", "n1_solv", "midHot", "deepKE", "chi2_med",
)
# W1 gets its own absolute floor instead of a seed-SD (plan §3.1).
W1_FLOOR = 0.15

# Reason written onto an injected marginal ion. Any non-suppressed,
# non-retained value scores identically (the scorer only special-cases those
# two); "time_exhausted" is the honest one — the cascade is still live.
INJECTED_REASON = "time_exhausted"

# ---------------------------------------------------------------------------


def observable_vector(read, abundance_ref, ked_ref) -> dict[str, float]:
    """The committed observable vector for one scored read (rule-1 reuse)."""
    hist = score_histogram_vs_reference(read, abundance_ref)
    ke = score_ke_curve_vs_reference(
        read, ked_ref, n_min=KE_N_MIN, min_count=MIN_KE_BIN_COUNT,
        include_sim_se=INCLUDE_SIM_SE, n1_anchor=N1_ANCHOR,
    )
    return {
        "num_scored": float(read.num_scored),
        "supp": read.suppressed_frac,
        "nbar_det": read.n_mean,
        "n1_solv": hist.sim_n1,
        "w1_solv": hist.w1_solvated,
        "midHot": midhot_ratio(read, ked_ref).value,
        "deepKE": deep_bin_ke_ratio(read, ked_ref).value,
        "chi2_med": ke.chi2_profiled,
    }


def measure_seed_sd(abundance_ref, ked_ref) -> tuple[dict[str, float], float]:
    """Per-observable seed SD from the five N = 1000 battery members.

    Returns ``(sd_by_column, mean_scored_count)``. Measured rather than
    remembered so the yardstick cannot drift away from the data it claims to
    come from.
    """
    vectors: list[dict[str, float]] = []
    for name in BATTERY_MEMBERS:
        detection = load_detection_result(RUNS_ROOT / name / "detection.npz")
        read = read_confirmation_detection(detection, label=name[-2:])
        vectors.append(observable_vector(read, abundance_ref, ked_ref))
    sd = {
        col: float(np.std([v[col] for v in vectors], ddof=1))
        for col in TEST_COLUMNS
    }
    sd["w1_solv"] = float(np.std([v["w1_solv"] for v in vectors], ddof=1))
    mean_scored = float(np.mean([v["num_scored"] for v in vectors]))
    return sd, mean_scored


def marginal_dossier(run_dir: Path) -> dict[str, Any] | None:
    """Asymptotic-escape dossier for one cell's marginal class.

    Recomputes the escape energetics from the stored ``relaxation.npz`` — the
    same call the detection stage made when it classified these ions — and
    returns the marginal mask together with the exact conservative asymptotic
    kinetic energies. ``None`` when the cell has no marginal ions.
    """
    detection = load_detection_result(run_dir / "detection.npz")
    reasons = np.asarray(detection.state_reason)
    marginal = reasons == RETAINED_MARGINAL_REASON
    if not np.any(marginal):
        return None

    cfg = RunDirectory(run_dir).load_cfg()
    seed_ckpt = load_ion_checkpoint(run_dir / "relaxation.npz")
    seed = ion_state_from_checkpoint_column(seed_ckpt, -1)
    energetics = escape_energetics(seed, seed_ckpt, cfg)

    radii = np.asarray(seed_ckpt.droplet_radii_angstrom, dtype=float)
    r = np.sqrt(np.asarray(seed.x) ** 2 + np.asarray(seed.y) ** 2
                + np.asarray(seed.z) ** 2)
    ke_asym = np.asarray(energetics.asymptotic_ke_eV)
    n_hand = np.asarray(detection.n_detected, dtype=float)
    return {
        "detection": detection,
        "marginal": marginal,
        "ke_asym_eV": ke_asym,
        "ke_ceiling_eV": np.asarray(energetics.asymptotic_ke_ceiling_eV),
        "n_handover": n_hand,
        "depth_A": (r - radii),          # negative inside the droplet
        "cfg": cfg,
    }


def arm_b_detection(dossier: dict[str, Any]):
    """Arm B: the same detection artifact with marginals injected as arrivals.

    Only two fields move — the state reason (so the scorer stops excluding
    them) and the kinetic energy (handover -> asymptotic). ``n_detected``
    already holds the verbatim handover ``n`` for retained rows, which is
    exactly the conservative-corner value Arm B claims.
    """
    detection = dossier["detection"]
    marginal = dossier["marginal"]
    reasons = np.asarray(detection.state_reason).copy()
    reasons[marginal] = INJECTED_REASON
    ke = np.asarray(detection.E_kin_detected_eV, dtype=float).copy()
    ke[marginal] = dossier["ke_asym_eV"][marginal]
    return dataclasses.replace(
        detection, state_reason=reasons, E_kin_detected_eV=ke
    )


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    print("Atlas §3.5b item 8 — one-off marginal-class bracket diagnostic\n")

    sd, member_scored = measure_seed_sd(abundance_ref, ked_ref)
    print(f"Seed-SD yardstick measured from {len(BATTERY_MEMBERS)} × N = 1000 "
          f"battery members (mean scored {member_scored:.0f} ions):")
    print("  " + "  ".join(f"{k}={v:.4g}" for k, v in sd.items()))
    print()

    dossier_rows: list[dict[str, Any]] = []
    bracket_rows: list[dict[str, Any]] = []
    verdicts: list[tuple[str, list[str]]] = []

    for label in CELL_LABELS:
        run_dir = RUNS_ROOT / RUN_SELECTION[label]
        if not (run_dir / "detection.npz").exists():
            continue
        dossier = marginal_dossier(run_dir)
        if dossier is None:
            continue

        marginal = dossier["marginal"]
        ke_m = dossier["ke_asym_eV"][marginal]
        ke_ceiling_m = dossier["ke_ceiling_eV"][marginal]
        n_m = dossier["n_handover"][marginal]
        depth_m = dossier["depth_A"][marginal]
        ref_mean_at_n = np.interp(
            n_m, ked_ref.n.astype(float), ked_ref.mean_KE_eV
        )
        dossier_rows.append({
            "label": label,
            "n_marg": int(marginal.sum()),
            "frac": float(marginal.mean()),
            "KEasym_med_eV": float(np.median(ke_m)),
            "KEasym_p90_eV": float(np.percentile(ke_m, 90)),
            "KEceil_med_eV": float(np.median(ke_ceiling_m)),
            "n_hand_med": float(np.median(n_m)),
            "depth_med_A": float(np.median(depth_m)),
            "refKE_at_n_med_eV": float(np.median(ref_mean_at_n)),
            "KEasym/refKE": float(np.median(ke_m / ref_mean_at_n)),
        })

        arm_a = observable_vector(
            read_confirmation_detection(dossier["detection"], label=label),
            abundance_ref, ked_ref,
        )
        arm_b = observable_vector(
            read_confirmation_detection(arm_b_detection(dossier), label=label),
            abundance_ref, ked_ref,
        )
        scale = np.sqrt(member_scored / max(arm_a["num_scored"], 1.0))
        row: dict[str, Any] = {"label": label}
        failures: list[str] = []
        for col in (*TEST_COLUMNS, "w1_solv"):
            delta = arm_b[col] - arm_a[col]
            row[f"{col}_A"] = arm_a[col]
            row[f"{col}_B"] = arm_b[col]
            row[f"d_{col}"] = delta
            if not np.isfinite(delta):
                continue
            limit = W1_FLOOR if col == "w1_solv" else sd[col] * scale
            if abs(delta) > limit:
                failures.append(f"{col} |Δ|={abs(delta):.4g} > {limit:.4g}")
        bracket_rows.append(row)
        verdicts.append((label, failures))

    if not dossier_rows:
        print("No cell has a marginal class — nothing to bracket.")
        return

    print("=== Marginal-class dossier (exact conservative asymptotics) ===")
    print(format_table(dossier_rows))
    print(
        "\nKEasym is the half-credit pair split (the physical asymptotic "
        "partition); KEceil is the full-credit upper bound. depth is r − R, "
        "negative inside the droplet. refKE is the experimental mean KE "
        "interpolated at the same n — it is the reference's first moment, NOT "
        "a detector acceptance, so 'KEasym/refKE' says how cold these "
        "fragments would arrive relative to the observed population, not "
        "whether the instrument would register them."
    )
    print()

    print("=== Arm A (excluded) vs Arm B (injected at asymptotic KE) ===")
    print(format_table(bracket_rows))
    print()

    print("=== Pre-registered drop rule ===")
    all_tight = True
    for label, failures in verdicts:
        if failures:
            all_tight = False
            print(f"  {label}: WIDE — " + "; ".join(failures))
        else:
            print(f"  {label}: tight")
    print()
    if all_tight:
        print(
            "VERDICT: TIGHT on every cell and every tested observable. Per the "
            "drop rule, record in the findings doc, strike the bracket columns "
            "from the grid report, and retire the §3.5b item-6 travelling "
            "caveat — but only for the cold-arrival direction: Arm B holds n "
            "at handover, so the pickup/mass-loading question stays open and "
            "belongs to the kinetic forward model."
        )
    else:
        print(
            "VERDICT: NOT TIGHT. The flagged observables above are biased by "
            "the marginal exclusion; keep the bracket columns, carry the "
            "§3.5b item-6 caveat, and escalate to the kinetic forward model "
            "(conservative orbit + live Poisson pickup + RRK evaporation), "
            "which additionally tests whether pickup mass-loads marginals "
            "into the bound class."
        )


if __name__ == "__main__":
    main()
