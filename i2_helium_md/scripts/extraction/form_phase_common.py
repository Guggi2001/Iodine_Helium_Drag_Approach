"""Shared helpers + locked constants for the §10 form-discrimination drivers.

Both sibling drivers (``method_b_form_refit_linear_quadratic.py`` and
``method_b_form_refit_power_law.py``) import from here so there is no
copy-paste between them (project quality rule 1). This module carries:

* the LOCKED §10.4.1 pre-registered constants (first-runs rule: fixed
  2026-06-12 before any form-phase run, not re-tuned after) as named
  constants with provenance,
* the reused §9.4 Stage-2 band check, made form-generic (it scores a
  :class:`JointFormObjectiveResult`),
* a form-generic fit summary,
* the two-threshold ``T_form`` classifier (§10.4.1),
* JSON load/dump + a date/branch stamp,
* a read-modify-write of the combined ``form_comparison_verdict.json``
  (keyed by family) that sits next to the §9 ``verdict.json``.

The §9.4 Stage-2 bands are reused UNCHANGED -- they are form-agnostic
statements about trajectory reproduction, so reusing them avoids any
post-hoc tuning (METHOD_B §10.4).
"""

from __future__ import annotations

import datetime
import json
from pathlib import Path
import sys

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.extraction.trajectory_matching import (  # noqa: E402
    JointFormObjectiveResult,
    SharedFormTrajectoryMatchingFit,
    _git_branch,
)
from i2_helium_md.presets import REFERENCE_DRAG_ROOT  # noqa: E402

# =============================================================================
# LOCKED §10.4.1 pre-registered constants (first-runs rule; do NOT re-tune)
# Provenance per the METHOD_B §10.4.1 table and form_phase_tasks_left.md.
# =============================================================================
CASES = ("18A", "9A")
OUT_ROOT = REFERENCE_DRAG_ROOT / "shared" / "trajectory_matching"

# --- §9 run-comparability (same N / seed / ion run as the shared refit) ------
N = 50                   # ensemble size per forward integration (Tier-0)
SEED = 20260604          # neutral-stage seed (ion drag path is RNG-free)
ION_TIME_PS = 20.0       # covers the [t*, ~14 ps] windows with escape margin
DT_ION_PS = 0.01

# --- reused §9.4 Stage-2 bands (gating; form-agnostic, unchanged) ------------
S2_RMSE_18A_MAX_APS = 0.19   # reused §9.4 band
S2_RMSE_9A_MAX_APS = 0.45    # reused §9.4 band
# Escape fraction must be exactly 1.0 for both cases (reused §9.4).

# --- Stage-1 analog (recorded, NON-gating) -----------------------------------
S1_ANALOG_RMSE_MAX_APS = 0.45  # comparability with pure-cubic's 0.2685

# --- two-threshold T_form classification (§10.4.1) ---------------------------
T_FORM_EQUIV_APS = 0.005     # equivalence-zone half-width (optimizer/T_a0 scale)
T_FORM_BETTER_APS = 0.013    # escalation bar (10%-sensitivity of the incumbent)
# Exact objective_mean_Aps of shared_pure_cubic/fit_parameters.json (the doc's
# 0.1293 is the rounded display); read this session.
INCUMBENT_OBJECTIVE_APS = 0.1292649398514104

# --- lq pure-quadratic equivalence (the §9.4 T_a0 reused per §10.4) ----------
T_A0_ANALOG_APS = 0.005

# --- power_law locked numerics -----------------------------------------------
N_BOUNDS = (1.0, 4.0)        # nesting points 2 and 3 interior; n>=1 -> gamma
                             # finite at rest (drag_low_v_floor stays inert)
V_REF_APS = 3.0              # pivot speed inside both in-window speed ranges

# --- pre-registered ANCHOR constants (read off Method-A bundles ONCE; values
# LOCKED here, NEVER a runtime bundle read -- the re-wired shared bundle's
# a0 = 0 must not condition any fit). 18 A-only, a convention choice (§10.4.1).
A0 = 14.555626399148123      # 18A/linear_and_cubic/fit_parameters.json (lq a-scale)
C0_LQ = 11.016050300970692   # 18A/quadratic/fit_parameters.json key "a" (lq c-scale, d1ca029)
PL_C0 = 10.36139380949775    # 18A/power/fit_parameters.json key "gamma" (legacy key; frozen evidence)
PL_N0 = 2.0557526931077588   # same bundle

ANCHOR_PROVENANCE = (
    "pre-registered §10.4.1 constants (18 A-only, a convention choice): "
    "a0=14.555626399148123 from 18A/linear_and_cubic; c0=11.016050300970692 "
    "from 18A/quadratic (key 'a', Method-A quadratic export d1ca029); "
    "C0=10.36139380949775, n0=2.0557526931077588 from 18A/power (key 'gamma', "
    "legacy key, frozen evidence read once); V_REF_APS=3.0. Read off the "
    "Method-A bundles once and LOCKED in the driver USER SETTINGS -- never a "
    "runtime bundle read (the re-wired shared bundle's a0 = 0 must not "
    "condition a fit)."
)

# Sensitivity-only band (seed sweep omitted per the §8 seed-insensitivity
# finding -- same convention as the §9 shared fit).
UNCERTAINTY_MODEL = "rmse_sensitivity_only_seed_sweep_omitted_per_s8"

VERDICT_PATH = OUT_ROOT / "form_comparison_verdict.json"


# =============================================================================
# JSON I/O + stamp
# =============================================================================
def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def dump_json(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    print(f"  wrote {path}")


def stamp() -> dict:
    return {
        "date": datetime.date.today().isoformat(),
        "branch": _git_branch(PROJECT_ROOT),
    }


# =============================================================================
# §9.4 Stage-2 band check (form-generic) + fit summary
# =============================================================================
def band_checks(joint_best: JointFormObjectiveResult) -> dict:
    """Apply the reused §9.4 Stage-2 bands to one joint form result.

    Form-agnostic: scores ``joint_best.per_case`` (FormObjectiveResult per
    case). Returns the per-band records plus an ``all_pass`` rollup. Only
    these Stage-2 bands carry verdict power (§9.4 / §10.4).
    """
    r18 = joint_best.per_case["18A"]
    r9 = joint_best.per_case["9A"]
    checks = {
        "S2_RMSE_18A": {
            "value_Aps": r18.rmse_Aps, "band_max_Aps": S2_RMSE_18A_MAX_APS,
            "pass": bool(r18.rmse_Aps <= S2_RMSE_18A_MAX_APS),
        },
        "S2_RMSE_9A": {
            "value_Aps": r9.rmse_Aps, "band_max_Aps": S2_RMSE_9A_MAX_APS,
            "pass": bool(r9.rmse_Aps <= S2_RMSE_9A_MAX_APS),
        },
        "S2_ESCAPE": {
            "value": {"18A": r18.escape_fraction, "9A": r9.escape_fraction},
            "band": 1.0,
            "pass": bool(
                r18.escape_fraction == 1.0 and r9.escape_fraction == 1.0
            ),
        },
    }
    checks["all_pass"] = all(
        c["pass"] for k, c in checks.items() if k != "all_pass"
    )
    return checks


def fit_summary(fit: SharedFormTrajectoryMatchingFit) -> dict:
    """Form-generic summary (coefficients as the raw closed-form dict)."""
    best = fit.best
    return {
        "variant": fit.variant,
        "form": fit.form,
        "coefficients": dict(best.coefficients),
        "e_bind_eV": best.e_bind_eV,
        "objective_mean_Aps": best.objective,
        "per_case": {
            c: {
                "rmse_Aps": r.rmse_Aps,
                "escape_fraction": r.escape_fraction,
                "objective_Aps": r.objective,
                "num_overlap_points": r.num_overlap_points,
            }
            for c, r in best.per_case.items()
        },
        "n_evaluations": fit.n_evaluations,
        "converged": fit.converged,
        "start_minima": [
            {
                "coefficients": dict(r.coefficients),
                "e_bind_eV": r.e_bind_eV,
                "objective_Aps": r.objective,
            }
            for r in fit.starts
        ],
    }


def fit_is_trapped(fit: SharedFormTrajectoryMatchingFit) -> bool:
    """True if any case's best result carries a nonzero escape penalty."""
    return any(r.penalty_Aps != 0.0 for r in fit.best.per_case.values())


# =============================================================================
# Two-threshold T_form classification (§10.4.1)
# =============================================================================
def classify_t_form(best_objective_Aps: float, bands_pass: bool) -> dict:
    """Classify a family against the pure-cubic incumbent (§10.4.1).

    ``delta = best_variant_objective - INCUMBENT_OBJECTIVE_APS``. Zones
    (asymmetric, two-threshold):

    * ``delta <= -T_FORM_BETTER_APS``  -> ``genuinely_better`` (escalate, but
      ONLY if it also passes the §9.4 bands -- a family must pass the bands to
      be classified better at all);
    * ``-T_FORM_BETTER_APS < delta < -T_FORM_EQUIV_APS`` -> ``marginally_better``
      (recorded, no escalation; inside the objective's own flatness scale);
    * ``|delta| <= T_FORM_EQUIV_APS`` -> ``equivalent`` (exponent degeneracy
      recorded; pure-cubic stays the candidate; VMI post-Tier-1 arbitrates);
    * ``delta > +T_FORM_EQUIV_APS`` -> ``worse`` (rejected-by-trajectory-
      objective).
    """
    delta = best_objective_Aps - INCUMBENT_OBJECTIVE_APS
    if delta <= -T_FORM_BETTER_APS:
        zone = "genuinely_better"
    elif delta < -T_FORM_EQUIV_APS:
        zone = "marginally_better"
    elif delta <= T_FORM_EQUIV_APS:
        zone = "equivalent"
    else:
        zone = "worse"

    # A family must ALSO pass the §9.4 bands to count as a real improvement.
    escalate = bool(zone == "genuinely_better" and bands_pass)
    if zone == "genuinely_better" and not bands_pass:
        interpretation = (
            "objective beats the incumbent beyond T_FORM_BETTER but the "
            "family FAILS the §9.4 bands -> NOT a competing candidate "
            "(bands gate); recorded rejected-by-bands"
        )
    elif zone == "genuinely_better":
        interpretation = (
            "beats the incumbent beyond T_FORM_BETTER and passes the §9.4 "
            "bands -> competing production candidate; ESCALATE to the user"
        )
    elif zone == "marginally_better":
        interpretation = (
            "better than the incumbent but inside the objective's flatness "
            "scale (-0.013 < delta < -0.005) -> recorded marginally-better, "
            "not escalation-worthy; pure-cubic stays the candidate"
        )
    elif zone == "equivalent":
        interpretation = (
            "indistinguishable from the incumbent within +/-T_FORM_EQUIV -> "
            "exponent degeneracy recorded as the finding; pure-cubic stays "
            "the candidate; VMI post-Tier-1 arbitrates"
        )
    else:
        interpretation = (
            "worse than the incumbent beyond +T_FORM_EQUIV -> family bundle "
            "recorded rejected-by-trajectory-objective; incumbent confirmed"
        )
    return {
        "delta_objective_Aps": delta,
        "incumbent_objective_Aps": INCUMBENT_OBJECTIVE_APS,
        "best_variant_objective_Aps": best_objective_Aps,
        "T_FORM_EQUIV_APS": T_FORM_EQUIV_APS,
        "T_FORM_BETTER_APS": T_FORM_BETTER_APS,
        "bands_pass": bands_pass,
        "zone": zone,
        "escalate_to_user": escalate,
        "interpretation": interpretation,
    }


# =============================================================================
# Combined form_comparison_verdict.json (read-modify-write, keyed by family)
# =============================================================================
def update_form_comparison_verdict(family_key: str, family_record: dict) -> None:
    """Merge one family's record into the combined verdict file.

    Sits next to the §9 ``verdict.json`` in ``OUT_ROOT``; keyed by family so
    the two sibling drivers can each write their own slice without clobbering
    the other.
    """
    record = load_json(VERDICT_PATH) if VERDICT_PATH.is_file() else {}
    record.setdefault("_meta", {})
    record["_meta"].update({
        **stamp(),
        "phase": "METHOD_B §10 alternative drag-form discrimination",
        "incumbent": "shared_pure_cubic (gamma = g*b*v^2, n=3)",
        "incumbent_objective_Aps": INCUMBENT_OBJECTIVE_APS,
        "note": (
            "model selection on SEEN data under a pre-registered protocol "
            "(the §9 cross-case axis was spent by the Stage-2 joint fit); "
            "NOT fresh held-out validation -- the winner's external test is "
            "VMI after Tier 1. Presets are NOT re-wired; Tier 1 not started; "
            "VMI stays pending/post-Tier-1."
        ),
    })
    record[family_key] = family_record
    dump_json(record, VERDICT_PATH)
