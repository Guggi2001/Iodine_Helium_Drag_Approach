"""Method-B shared-form joint refit driver (METHOD_B doc §9).

Two stages, run in order so the held-out content is recorded BEFORE the
joint fit consumes the cross-case axis:

* **Stage 1** -- score the **untouched 9A prediction** from the delivered §8
  18A per-case Method-B bundle (NO refit, per the 2026-06-11 pre-run
  clarification: the §8 bundle is reused as-is). One forward integration +
  scoring; recorded to ``stage1_prediction.json``. **Recorded-only,
  non-gating** -- only the Stage-2 bands carry verdict power.
* **Stage 2** -- the joint fit over both cases: primary ``shared_3param``
  (fully shared {a, b, E_bind}, ``a`` free with lower bound 0) plus the
  ``shared_pure_cubic`` (``a = 0``) variant; the pure-cubic-equivalence
  threshold ``T_a0`` classifies the empirical conclusion on ``a``'s
  identifiability. On a primary band failure the ``diagnostic_4param``
  variant (shared {a, b}, per-case E_bind) localizes drag-vs-binding.

All thresholds below are the PRE-REGISTERED §9.4 bands (first-runs rule:
fixed before any run, not re-tuned after; anchors recorded in §9.4). They
SUPERSEDE the §8 provisional cross-case bands -- in particular the 18A A<->B
a-ratio check is NOT re-applied (2026-06-11 decision; it would auto-fail any
a -> 0 fit and is dead under the weak-a finding).

Verdict mapping (§9.4): Stage-2 primary passes its bands -> the shared
bundle is the production candidate and **Tier 1 ungates** (VMI becomes the
post-Tier-1 final arbiter). Primary fails -> diagnostic localizes; the
per-case bundles stay flagged not-yet-usable.

Both shared-variant bundles are recorded (variant subdirectories under
``data/reference/drag/shared/trajectory_matching/``); the production-
candidate choice between them is deferred to the preset-rewiring decision.
**Presets are NOT re-wired by this script.**

DRY_RUN evaluates the joint objective twice at the anchor point, asserts
bitwise determinism (the drag ion path must be RNG-free), prints timing,
and writes nothing.

Usage: edit USER SETTINGS, then::

    python scripts/extraction/method_b_shared_refit.py
"""

from __future__ import annotations

import datetime
import json
from pathlib import Path
import sys
import time

# =============================================================================
# USER SETTINGS
# =============================================================================
N = 50                  # ensemble size per forward integration (Tier-0 evidence)
SEED = 20260604         # neutral-stage seed (ion drag path is RNG-free)
ION_TIME_PS = 20.0      # covers the [t*, ~14 ps] windows with escape margin
DT_ION_PS = 0.01
DRY_RUN = False         # True: 1 joint eval x2 + determinism assert, no writes
MAXFEV_PER_START = 400
RUN_DIAGNOSTIC_ON_FAIL = True

# --- Pre-registered §9.4 bands (FIXED before any run; do not re-tune) -------
S1_PRED_RMSE_MAX_APS = 0.45   # Stage-1 9A prediction RMSE (recorded, NON-gating)
S2_RMSE_18A_MAX_APS = 0.19    # Stage-2 joint per-case 18A RMSE (gating)
S2_RMSE_9A_MAX_APS = 0.45     # Stage-2 joint per-case 9A RMSE (gating)
T_A0_APS = 0.005              # pure-cubic equivalence: obj(a==0) - obj(a free)
# Escape fraction must be exactly 1.0 everywhere (S1_PRED_ESCAPE / S2_ESCAPE).

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.extraction.trajectory_matching import (  # noqa: E402
    DIAGNOSTIC_4PARAM,
    SHARED_3PARAM,
    SHARED_PURE_CUBIC,
    _git_branch,
    build_case_setup,
    evaluate_joint_objective,
    evaluate_objective,
    fit_shared_trajectory_matching,
    joint_sensitivity_halfwidths,
    write_shared_fit_parameters,
)
from i2_helium_md.presets import REFERENCE_DRAG_ROOT  # noqa: E402

CASES = ("18A", "9A")
OUT_ROOT = REFERENCE_DRAG_ROOT / "shared" / "trajectory_matching"
STAGE1_BUNDLE_PATH = (
    REFERENCE_DRAG_ROOT / "18A" / "trajectory_matching" / "fit_parameters.json"
)

# The §8 seed sweep showed the HeDFT-comparison presets are insensitive to
# the neutral seed (parameter std ~ 0), so the shared refit carries a
# sensitivity-only band and omits the sweep (recorded in the bundle label).
UNCERTAINTY_MODEL = "rmse_sensitivity_only_seed_sweep_omitted_per_s8"


def _load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def _dump_json(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    print(f"  wrote {path}")


def _stamp() -> dict:
    return {
        "date": datetime.date.today().isoformat(),
        "branch": _git_branch(PROJECT_ROOT),
    }


def _band_checks(joint_best) -> dict:
    """Apply the pre-registered Stage-2 §9.4 bands to one joint result."""
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


def _fit_summary(fit) -> dict:
    return {
        "variant": fit.variant,
        "a": fit.best.a,
        "b": fit.best.b,
        "e_bind_eV": fit.best.e_bind_eV,
        "objective_mean_Aps": fit.best.objective,
        "per_case": {
            c: {
                "rmse_Aps": r.rmse_Aps,
                "escape_fraction": r.escape_fraction,
                "objective_Aps": r.objective,
            }
            for c, r in fit.best.per_case.items()
        },
        "n_evaluations": fit.n_evaluations,
        "converged": fit.converged,
        "start_minima": [
            {"a": r.a, "b": r.b, "e_bind_eV": r.e_bind_eV,
             "objective_Aps": r.objective}
            for r in fit.starts
        ],
    }


def dry_run(setups, anchors) -> None:
    a0, b0 = anchors
    e_by_case = {c: 0.2 for c in CASES}
    t0 = time.perf_counter()
    r1 = evaluate_joint_objective(a0, b0, e_by_case, setups)
    t1 = time.perf_counter()
    r2 = evaluate_joint_objective(a0, b0, e_by_case, setups)
    print(f"  joint eval at (a0={a0:.4f}, b0={b0:.4f}, E=0.2 eV), "
          f"{t1 - t0:.2f} s/joint-eval (= 2 ion runs):")
    print(f"    mean objective = {r1.objective:.6f} A/ps; per case: "
          + ", ".join(
              f"{c} rmse {r.rmse_Aps:.4f} escape {r.escape_fraction:.3f}"
              for c, r in r1.per_case.items()
          ))
    if r1.objective != r2.objective:
        raise AssertionError(
            "joint objective is NOT deterministic: repeated evaluation "
            f"differs ({r1.objective!r} vs {r2.objective!r}) -- "
            "stochasticity has leaked into the ion drag path"
        )
    print("    repeated evaluation bitwise identical -> deterministic OK")
    print("  DRY_RUN: nothing written.")


def main() -> int:
    print("===== Method-B shared-form joint refit (METHOD_B §9) =====")
    print(f"  N={N}, seed={SEED}, ion {ION_TIME_PS} ps @ dt={DT_ION_PS} ps")

    # --- per-case contexts (one neutral run each) ---
    setups = {}
    for case in CASES:
        print(f"  building {case} case setup (one neutral run) ...")
        setups[case] = build_case_setup(
            case, num_molecules=N, seed=SEED,
            ion_time_ps=ION_TIME_PS, dt_ion_ps=DT_ION_PS,
        )
        w = setups[case].window
        print(f"    reference {setups[case].reference_path.name}, "
              f"window [{w[0]:.2f}, {w[1]:.2f}] ps")

    # §9.3: anchors from the 18A Method-A bundle ONLY (the 18A drag preset
    # wires exactly that bundle; the 9A Method-A artifact is provenance-
    # broken and must not condition the fit).
    anchors = (setups["18A"].a0, setups["18A"].b0)
    print(f"  anchors (18A Method-A only): a0={anchors[0]:.4f}, "
          f"b0={anchors[1]:.4f}")

    if DRY_RUN:
        dry_run(setups, anchors)
        return 0

    # =====================================================================
    # Stage 1 -- 9A prediction from the delivered §8 18A bundle (no refit)
    # =====================================================================
    print("  --- Stage 1: 9A prediction from the §8 18A bundle (no refit) ---")
    bundle18 = _load_json(STAGE1_BUNDLE_PATH)
    s1_params = (
        float(bundle18["a"]),
        float(bundle18["b"]),
        float(bundle18["effective_binding_energy_I_ion_eV"]),
    )
    print(f"    18A-fit params (reused, 2026-06-11 decision): "
          f"a={s1_params[0]:.4f}, b={s1_params[1]:.4f}, "
          f"E_bind={s1_params[2]:.4f} eV")
    s1 = evaluate_objective(s1_params, setups["9A"])
    s1_rmse_pass = bool(s1.rmse_Aps <= S1_PRED_RMSE_MAX_APS)
    s1_escape_pass = bool(s1.escape_fraction == 1.0)
    print(f"    9A prediction: rmse={s1.rmse_Aps:.4f} A/ps "
          f"(band <= {S1_PRED_RMSE_MAX_APS}) -> "
          f"{'PASS' if s1_rmse_pass else 'FAIL'}; "
          f"escape={s1.escape_fraction:.3f} (band 1.0) -> "
          f"{'PASS' if s1_escape_pass else 'FAIL'}   [recorded, NON-gating]")

    stage1_record = {
        **_stamp(),
        "stage": "stage1_prediction",
        "gating": (
            "recorded_only_non_gating (2026-06-11 decision: only the "
            "Stage-2 bands carry verdict power)"
        ),
        "params_source": {
            "bundle": str(
                STAGE1_BUNDLE_PATH.relative_to(PROJECT_ROOT).as_posix()
            ),
            "note": (
                "delivered §8 18A per-case Method-B fit reused as-is "
                "(old 0.1*a0 lower bound; no refit under the new a >= 0 "
                "bound, per the 2026-06-11 pre-run clarification)"
            ),
            "a": s1_params[0], "b": s1_params[1], "e_bind_eV": s1_params[2],
        },
        "prediction_case": "9A",
        "rmse_Aps": s1.rmse_Aps,
        "escape_fraction": s1.escape_fraction,
        "num_overlap_points": s1.num_overlap_points,
        "checks": {
            "S1_PRED_RMSE_MAX": {
                "value_Aps": s1.rmse_Aps,
                "band_max_Aps": S1_PRED_RMSE_MAX_APS,
                "pass": s1_rmse_pass,
            },
            "S1_PRED_ESCAPE": {
                "value": s1.escape_fraction, "band": 1.0,
                "pass": s1_escape_pass,
            },
        },
        "anchor_note": (
            "0.45 = Tier-0 9A dimensionality residual (~0.39 A/ps) + margin "
            "(METHOD_B §9.4)"
        ),
    }
    _dump_json(stage1_record, OUT_ROOT / "stage1_prediction.json")

    # =====================================================================
    # Stage 2 -- joint fits (primary + pure-cubic variant)
    # =====================================================================
    fits = {}
    for variant in (SHARED_3PARAM, SHARED_PURE_CUBIC):
        print(f"  --- Stage 2: {variant} joint fit ---")
        t0 = time.perf_counter()
        fits[variant] = fit_shared_trajectory_matching(
            setups, variant=variant, anchors=anchors,
            maxfev_per_start=MAXFEV_PER_START,
        )
        t1 = time.perf_counter()
        fit = fits[variant]
        best = fit.best
        print(f"    done in {t1 - t0:.0f} s, {fit.n_evaluations} joint evals "
              f"(x2 ion runs each), converged={fit.converged}")
        for i, r in enumerate(fit.starts):
            print(f"      start {i}: a={r.a:.4f}, b={r.b:.4f}, "
                  f"E={r.e_bind_eV[fit.cases[0]]:.4f}, "
                  f"objective={r.objective:.4f}")
        print(f"    BEST: a={best.a:.4f} amu/ps, b={best.b:.4f} amu*ps/A^2, "
              f"E_bind={best.e_bind_eV[fit.cases[0]]:.4f} eV, "
              f"mean objective={best.objective:.4f} A/ps")
        for c, r in best.per_case.items():
            print(f"      {c}: rmse={r.rmse_Aps:.4f} A/ps, "
                  f"escape={r.escape_fraction:.3f}")

    fit3 = fits[SHARED_3PARAM]
    fitc = fits[SHARED_PURE_CUBIC]

    # --- pure-cubic equivalence (T_a0; sign per the corrected §9.4) ---
    delta_a0 = fitc.best.objective - fit3.best.objective
    pure_cubic_equivalent = bool(delta_a0 <= T_A0_APS)
    print(f"  T_a0: obj(a==0) - obj(a free) = {delta_a0:.6f} A/ps "
          f"(<= {T_A0_APS}) -> "
          + ("EQUIVALENT (pure-cubic recorded as the empirical conclusion "
             "on a's identifiability)" if pure_cubic_equivalent
             else "DISTINGUISHABLE (a carries real information)"))

    # --- pre-registered Stage-2 band checks ---
    checks3 = _band_checks(fit3.best)
    checksc = _band_checks(fitc.best)
    primary_pass = checks3["all_pass"]
    for name, checks in ((SHARED_3PARAM, checks3), (SHARED_PURE_CUBIC, checksc)):
        print(f"  bands [{name}]: "
              + ", ".join(
                  f"{k} {'PASS' if c['pass'] else 'FAIL'}"
                  for k, c in checks.items() if k != "all_pass"
              )
              + f" -> {'ALL PASS' if checks['all_pass'] else 'FAIL'}")

    # --- Stage-1 <-> Stage-2 parameter shift (qualitative, recorded-only) ---
    e3 = fit3.best.e_bind_eV[fit3.cases[0]]
    shift_record = {
        "note": (
            "qualitative, recorded-only (2026-06-11 decision: no "
            "pre-registered band, no gate power); a large shift when 9A "
            "enters the objective flags 9A dragging the law toward its "
            "transverse contamination (METHOD_B §9.2)"
        ),
        "stage1_18A_fit": {
            "a": s1_params[0], "b": s1_params[1], "e_bind_eV": s1_params[2],
        },
        "stage2_shared_3param": {
            "a": fit3.best.a, "b": fit3.best.b, "e_bind_eV": e3,
        },
        "shift": {
            "delta_a": fit3.best.a - s1_params[0],
            "delta_b": fit3.best.b - s1_params[1],
            "delta_e_bind_eV": e3 - s1_params[2],
        },
    }
    print(f"  Stage-1<->Stage-2 shift (recorded): "
          f"da={shift_record['shift']['delta_a']:+.4f}, "
          f"db={shift_record['shift']['delta_b']:+.4f}, "
          f"dE={shift_record['shift']['delta_e_bind_eV']:+.4f} eV")

    # --- diagnostic variant (only on primary failure) ---
    diagnostic_summary = None
    if not primary_pass and RUN_DIAGNOSTIC_ON_FAIL:
        print("  --- primary FAILED its bands: running diagnostic_4param ---")
        fitd = fit_shared_trajectory_matching(
            setups, variant=DIAGNOSTIC_4PARAM, anchors=anchors,
            maxfev_per_start=MAXFEV_PER_START,
        )
        diagnostic_summary = _fit_summary(fitd)
        e_by_case = fitd.best.e_bind_eV
        print(f"    diagnostic: a={fitd.best.a:.4f}, b={fitd.best.b:.4f}, "
              + ", ".join(f"E_{c}={e:.4f}" for c, e in e_by_case.items())
              + f", mean objective={fitd.best.objective:.4f}")
        _dump_json(
            {**_stamp(), **diagnostic_summary,
             "purpose": (
                 "failure localization only (drag law vs binding); never a "
                 "production bundle (METHOD_B §9.2)"
             )},
            OUT_ROOT / "diagnostic_4param.json",
        )

    # --- sensitivity bands + bundle writing (both shared variants) ---
    bundle_paths = {}
    for variant, fit in ((SHARED_3PARAM, fit3), (SHARED_PURE_CUBIC, fitc)):
        trapped = any(
            r.penalty_Aps != 0.0 for r in fit.best.per_case.values()
        )
        if trapped:
            print(f"  NOT writing {variant} bundle: fit traps ions "
                  "(recorded in the verdict).")
            continue
        print(f"  sensitivity scan ({variant}) ...")
        hw_a, hw_b, hw_e = joint_sensitivity_halfwidths(setups, fit.best)
        fit.a_err, fit.b_err, fit.e_bind_err_eV = hw_a, hw_b, hw_e
        fit.uncertainty_model = UNCERTAINTY_MODEL
        out = write_shared_fit_parameters(fit, OUT_ROOT / variant)
        bundle_paths[variant] = str(out.relative_to(PROJECT_ROOT).as_posix())
        print(f"  wrote {out}")

    # --- verdict (§9.4 mapping; 2026-06-11 gating clarifications) ---
    if primary_pass:
        verdict = (
            "PASS: the shared_3param fit passes the pre-registered Stage-2 "
            "bands -> the shared bundle is the production candidate; Tier 1 "
            "UNGATES (gate policy §9.1); VMI = post-Tier-1 final arbiter. "
            "Production-candidate choice between the recorded variants is "
            "deferred to the preset-rewiring decision."
        )
    else:
        verdict = (
            "FAIL: the shared_3param fit violates the pre-registered "
            "Stage-2 bands -> the shared form does not generalize as "
            "fitted; Tier 1 STAYS GATED; per-case bundles stay flagged "
            "not-yet-usable; see diagnostic_4param.json for the "
            "localization; EXTRACTION_FRAME_FIX_milestone.md remains the "
            "recorded contingency."
        )
    print(f"  VERDICT: {verdict}")

    thresholds = {
        "S1_PRED_RMSE_MAX_Aps": S1_PRED_RMSE_MAX_APS,
        "S2_RMSE_18A_MAX_Aps": S2_RMSE_18A_MAX_APS,
        "S2_RMSE_9A_MAX_Aps": S2_RMSE_9A_MAX_APS,
        "escape_band": 1.0,
        "T_a0_Aps": T_A0_APS,
        "provenance": (
            "pre-registered §9.4 bands (fixed 2026-06-11 before any run; "
            "anchors recorded in METHOD_B §9.4); SUPERSEDE the §8 "
            "provisional cross-case bands -- the 18A A<->B a-ratio check "
            "is not re-applied"
        ),
        "gating": (
            "Stage-2 bands only (Stage 1 recorded non-gating; T_a0 is a "
            "classification, not a gate) -- 2026-06-11 decisions"
        ),
    }
    verdict_record = {
        **_stamp(),
        "thresholds": thresholds,
        "stage1": stage1_record,
        "stage2": {
            SHARED_3PARAM: {**_fit_summary(fit3), "bands": checks3},
            SHARED_PURE_CUBIC: {**_fit_summary(fitc), "bands": checksc},
        },
        "pure_cubic_equivalence": {
            "delta_objective_Aps": delta_a0,
            "threshold_Aps": T_A0_APS,
            "equivalent": pure_cubic_equivalent,
            "conclusion": (
                "pure-cubic recorded as the empirical conclusion on a's "
                "identifiability (a not identified by this objective)"
                if pure_cubic_equivalent else
                "a-free and pure-cubic distinguishable: a carries real "
                "information over these windows"
            ),
        },
        "stage1_stage2_parameter_shift": shift_record,
        "diagnostic_4param": diagnostic_summary,
        "primary_pass": primary_pass,
        "bundles_written": bundle_paths,
        "verdict": verdict,
        "vmi": "pending (unreachable until Tier 1; post-Tier-1 final arbiter)",
    }
    _dump_json(verdict_record, OUT_ROOT / "verdict.json")

    # --- update the per-case held_out_validation.json files (§9.6) ---
    for case in CASES:
        path = (
            REFERENCE_DRAG_ROOT / case / "trajectory_matching"
            / "held_out_validation.json"
        )
        record = _load_json(path) if path.is_file() else {}
        record["shared_form_refit"] = {
            **_stamp(),
            "note": (
                "scored against the pre-registered §9.4 bands, which "
                "supersede the provisional first-runs bands above "
                "(2026-06-11 decision; the A<->B a-ratio check is not "
                "re-applied)"
            ),
            "stage1_9A_prediction": {
                "rmse_Aps": s1.rmse_Aps,
                "escape_fraction": s1.escape_fraction,
                "gating": "recorded_only_non_gating",
            },
            "stage2_shared_3param": {
                "per_case_rmse_Aps": {
                    c: r.rmse_Aps for c, r in fit3.best.per_case.items()
                },
                "bands": checks3,
            },
            "primary_pass": primary_pass,
            "verdict": verdict,
            "shared_artifacts": str(
                OUT_ROOT.relative_to(PROJECT_ROOT).as_posix()
            ),
        }
        _dump_json(record, path)

    print("  NOTE: presets are NOT re-wired by this script (gated on the "
          "verdict plus a separate user decision).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
