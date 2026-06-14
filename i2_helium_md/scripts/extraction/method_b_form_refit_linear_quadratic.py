"""Method-B §10 form-discrimination driver -- ``linear_quadratic`` family.

Sibling of the frozen §9 ``method_b_shared_refit.py`` (which stays the
incumbent's record). Runs the ``linear_quadratic`` family
(``F = g*(a*v + c*|v|*v)``, ``gamma = g*(a + c*|v|)``) through the same
shared-form joint-refit machinery against the pure-cubic incumbent, including
the **pure-quadratic** ``a == 0`` variant (the Method-A ``n ~ +2`` hypothesis
as a closed form, METHOD_B §10.2/§10.3).

Stages (same order/discipline as §9):

* **Stage-1 analog** (recorded, NON-gating, honestly weakened -- the 9 A data
  has been seen repeatedly, §10.4): fit the family's FULL variant
  (``lq_shared_3param``) on **18 A only**, then score the untouched 9 A
  prediction vs the same 0.45 A/ps number for comparability with pure-cubic's
  0.2685. Record-only ``stage1_analog_linear_quadratic.json`` -- no bundle
  (the writer refuses single-case fits).
* **Stage-2 shared joint fits**: ``lq_shared_3param`` (``a`` free, lower bound
  0) + ``lq_shared_pure_quadratic`` (``a == 0``), then the ``T_a0``-analog
  equivalence classification between them (threshold 0.005 A/ps,
  Delta = obj(a==0) - obj(a-free) >= 0 nested).

Verdict scoring (all pre-registered, LOCKED §10.4.1 -- see
``form_phase_common.py``): the reused §9.4 Stage-2 bands (gating) plus the
two-threshold ``T_form`` classification against the incumbent. Escalation is
reserved for a family that beats the incumbent beyond ``T_FORM_BETTER`` AND
passes the bands; otherwise the result is recorded and the run stops. This is
**model selection on seen data** -- not fresh held-out validation; the
winner's external test is VMI after Tier 1. **Presets are NOT re-wired.**

DRY_RUN evaluates the joint objective twice at the anchor point, asserts
bitwise determinism (the drag ion path must be RNG-free), prints timing, and
writes nothing.

Usage: edit USER SETTINGS, then::

    python scripts/extraction/method_b_form_refit_linear_quadratic.py
"""

from __future__ import annotations

import time

from form_phase_common import (
    A0,
    ANCHOR_PROVENANCE,
    C0_LQ,
    CASES,
    DT_ION_PS,
    ION_TIME_PS,
    N,
    OUT_ROOT,
    PROJECT_ROOT,
    S1_ANALOG_RMSE_MAX_APS,
    SEED,
    T_A0_ANALOG_APS,
    UNCERTAINTY_MODEL,
    band_checks,
    classify_t_form,
    dump_json,
    fit_is_trapped,
    fit_summary,
    stamp,
    update_form_comparison_verdict,
)
from i2_helium_md.extraction.trajectory_matching import (
    LQ_SHARED_3PARAM,
    LQ_SHARED_PURE_QUADRATIC,
    build_case_setup,
    evaluate_form_objective,
    evaluate_joint_form_objective,
    fit_form_trajectory_matching,
    fit_shared_form_trajectory_matching,
    form_sensitivity_halfwidths,
    write_shared_form_fit_parameters,
)
from i2_helium_md.physics.drag import LINEAR_QUADRATIC

# =============================================================================
# USER SETTINGS
# =============================================================================
DRY_RUN = False          # True: 1 joint eval x2 + determinism assert, no writes
MAXFEV_PER_START = 400   # §9 comparability

FAMILY = "linear_quadratic"
FORM = LINEAR_QUADRATIC
FULL_VARIANT = LQ_SHARED_3PARAM
VARIANTS = (LQ_SHARED_3PARAM, LQ_SHARED_PURE_QUADRATIC)
# Pre-registered anchors (§10.4.1 LOCKED constants; NEVER setup.a0/b0 -- the
# re-wired shared bundle reads a0 = 0).
ANCHORS = {"a0": A0, "c0": C0_LQ}


def dry_run(setups) -> None:
    coeffs = {"a": A0, "c": C0_LQ}
    e_by_case = {c: 0.2 for c in CASES}
    t0 = time.perf_counter()
    r1 = evaluate_joint_form_objective(FORM, coeffs, e_by_case, setups)
    t1 = time.perf_counter()
    r2 = evaluate_joint_form_objective(FORM, coeffs, e_by_case, setups)
    print(f"  joint eval at (a0={A0:.4f}, c0={C0_LQ:.4f}, E=0.2 eV), "
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
    print("===== Method-B §10 form refit: linear_quadratic =====")
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

    print(f"  anchors (§10.4.1 LOCKED, 18A-only): a0={A0:.4f}, c0={C0_LQ:.4f}")

    if DRY_RUN:
        dry_run(setups)
        return 0

    # =====================================================================
    # Stage 1 analog -- fit lq_shared_3param on 18A only, predict 9A
    # =====================================================================
    print("  --- Stage 1 analog: fit lq_shared_3param on 18A, predict 9A ---")
    s1fit = fit_form_trajectory_matching(
        setups["18A"], variant=FULL_VARIANT, anchors=ANCHORS,
        maxfev_per_start=MAXFEV_PER_START,
    )
    s1_coeffs = dict(s1fit.best.coefficients)
    s1_e_bind = s1fit.best.e_bind_eV["18A"]
    print(f"    18A-only fit: a={s1_coeffs['a']:.4f}, c={s1_coeffs['c']:.4f}, "
          f"E_bind={s1_e_bind:.4f} eV, "
          f"objective={s1fit.best.objective:.4f} A/ps")
    s1pred = evaluate_form_objective(FORM, s1_coeffs, s1_e_bind, setups["9A"])
    s1_rmse_pass = bool(s1pred.rmse_Aps <= S1_ANALOG_RMSE_MAX_APS)
    s1_escape_pass = bool(s1pred.escape_fraction == 1.0)
    print(f"    9A prediction: rmse={s1pred.rmse_Aps:.4f} A/ps "
          f"(band <= {S1_ANALOG_RMSE_MAX_APS}) -> "
          f"{'PASS' if s1_rmse_pass else 'FAIL'}; "
          f"escape={s1pred.escape_fraction:.3f} -> "
          f"{'PASS' if s1_escape_pass else 'FAIL'}   [recorded, NON-gating]")

    stage1_record = {
        **stamp(),
        "family": FAMILY,
        "stage": "stage1_analog",
        "gating": (
            "recorded_only_non_gating (§10.4: honestly weakened -- the 9A "
            "data has been seen repeatedly; audit-trail value, no verdict "
            "power; only the Stage-2 bands gate)"
        ),
        "fit_18A_only": {
            "variant": FULL_VARIANT,
            "coefficients": s1_coeffs,
            "e_bind_eV": s1_e_bind,
            "objective_Aps": s1fit.best.objective,
        },
        "prediction_case": "9A",
        "rmse_Aps": s1pred.rmse_Aps,
        "escape_fraction": s1pred.escape_fraction,
        "num_overlap_points": s1pred.num_overlap_points,
        "checks": {
            "S1_ANALOG_RMSE_MAX": {
                "value_Aps": s1pred.rmse_Aps,
                "band_max_Aps": S1_ANALOG_RMSE_MAX_APS,
                "pass": s1_rmse_pass,
            },
            "S1_ANALOG_ESCAPE": {
                "value": s1pred.escape_fraction, "band": 1.0,
                "pass": s1_escape_pass,
            },
        },
        "anchor_note": (
            "0.45 = same number as pure-cubic's Stage-1 (Tier-0 9A "
            "dimensionality residual ~0.39 + margin); compare to pure-cubic's "
            "0.2685 (METHOD_B §9.7 / §10.4.1)"
        ),
    }
    dump_json(stage1_record, OUT_ROOT / "stage1_analog_linear_quadratic.json")

    # =====================================================================
    # Stage 2 -- shared joint fits (a-free + pure-quadratic)
    # =====================================================================
    fits = {}
    for variant in VARIANTS:
        print(f"  --- Stage 2: {variant} joint fit ---")
        t0 = time.perf_counter()
        fits[variant] = fit_shared_form_trajectory_matching(
            setups, variant=variant, anchors=ANCHORS,
            maxfev_per_start=MAXFEV_PER_START,
        )
        t1 = time.perf_counter()
        fit = fits[variant]
        best = fit.best
        print(f"    done in {t1 - t0:.0f} s, {fit.n_evaluations} joint evals "
              f"(x2 ion runs each), converged={fit.converged}")
        for i, r in enumerate(fit.starts):
            print(f"      start {i}: {r.coefficients}, "
                  f"E={r.e_bind_eV[fit.cases[0]]:.4f}, "
                  f"objective={r.objective:.4f}")
        print(f"    BEST: a={best.coefficients['a']:.4f} amu/ps, "
              f"c={best.coefficients['c']:.4f} amu/A, "
              f"E_bind={best.e_bind_eV[fit.cases[0]]:.4f} eV, "
              f"mean objective={best.objective:.4f} A/ps")
        for c, r in best.per_case.items():
            print(f"      {c}: rmse={r.rmse_Aps:.4f} A/ps, "
                  f"escape={r.escape_fraction:.3f}")

    fit3 = fits[LQ_SHARED_3PARAM]
    fitq = fits[LQ_SHARED_PURE_QUADRATIC]

    # --- pure-quadratic equivalence (T_a0 analog; nested so Delta >= 0) ---
    delta_a0 = fitq.best.objective - fit3.best.objective
    pure_quad_equivalent = bool(delta_a0 <= T_A0_ANALOG_APS)
    print(f"  T_a0-analog: obj(a==0) - obj(a free) = {delta_a0:.6f} A/ps "
          f"(<= {T_A0_ANALOG_APS}) -> "
          + ("EQUIVALENT (the linear term carries no information; "
             "pure-quadratic is the empirical reduced form)"
             if pure_quad_equivalent
             else "DISTINGUISHABLE (a carries real information)"))

    # --- §9.4 Stage-2 band checks (gating) ---
    checks = {v: band_checks(fits[v].best) for v in VARIANTS}
    for v in VARIANTS:
        print(f"  bands [{v}]: "
              + ", ".join(
                  f"{k} {'PASS' if c['pass'] else 'FAIL'}"
                  for k, c in checks[v].items() if k != "all_pass"
              )
              + f" -> {'ALL PASS' if checks[v]['all_pass'] else 'FAIL'}")

    # --- family best variant (lowest objective) + T_form classification ---
    best_variant = min(VARIANTS, key=lambda v: fits[v].best.objective)
    best_obj = fits[best_variant].best.objective
    bands_pass = checks[best_variant]["all_pass"]
    t_form = classify_t_form(best_obj, bands_pass)
    print(f"  T_form: best variant = {best_variant} "
          f"(objective {best_obj:.6f} A/ps); "
          f"Delta = {t_form['delta_objective_Aps']:+.6f} A/ps -> "
          f"zone={t_form['zone']}, escalate={t_form['escalate_to_user']}")
    print(f"    {t_form['interpretation']}")

    # --- sensitivity bands + bundle writing (both variants) ---
    bundle_paths = {}
    for variant in VARIANTS:
        fit = fits[variant]
        if fit_is_trapped(fit):
            print(f"  NOT writing {variant} bundle: fit traps ions "
                  "(recorded in the verdict).")
            continue
        print(f"  sensitivity scan ({variant}) ...")
        hw = form_sensitivity_halfwidths(setups, fit.best)
        fit.coeff_errs = {"a": hw["a"], "c": hw["c"]}
        fit.e_bind_err_eV = hw["e_bind_eV"]
        fit.uncertainty_model = UNCERTAINTY_MODEL
        out = write_shared_form_fit_parameters(
            fit, OUT_ROOT / variant, anchor_provenance=ANCHOR_PROVENANCE,
        )
        bundle_paths[variant] = str(out.relative_to(PROJECT_ROOT).as_posix())
        print(f"  wrote {out}")

    # --- verdict string ---
    if t_form["escalate_to_user"]:
        verdict = (
            f"ESCALATE: {FAMILY} (best variant {best_variant}) beats the "
            "pure-cubic incumbent beyond T_FORM_BETTER AND passes the §9.4 "
            "bands -> a competing production candidate. User decision needed "
            "before any preset re-wiring."
        )
    elif t_form["zone"] == "marginally_better":
        verdict = (
            f"RECORDED (marginally better): {FAMILY} edges the incumbent but "
            "within the objective's flatness scale -> not escalation-worthy; "
            "pure-cubic stays the candidate."
        )
    elif t_form["zone"] == "equivalent":
        verdict = (
            f"RECORDED (equivalent): {FAMILY} is indistinguishable from the "
            "pure-cubic incumbent within T_FORM_EQUIV -> exponent degeneracy "
            "recorded as the finding; pure-cubic stays the candidate; VMI "
            "post-Tier-1 arbitrates."
        )
    else:  # worse, or "better"-but-bands-fail
        verdict = (
            f"RECORDED (rejected): {FAMILY} does not improve on the "
            "pure-cubic incumbent under the trajectory objective / §9.4 "
            "bands -> incumbent confirmed."
        )
    print(f"  VERDICT: {verdict}")

    # --- per-family verdict + combined verdict update ---
    family_record = {
        **stamp(),
        "family": FAMILY,
        "form": FORM,
        "variants": list(VARIANTS),
        "anchors": {**ANCHORS, "provenance": ANCHOR_PROVENANCE},
        "stage1_analog": stage1_record,
        "stage2": {v: {**fit_summary(fits[v]), "bands": checks[v]}
                   for v in VARIANTS},
        "pure_quadratic_equivalence": {
            "delta_objective_Aps": delta_a0,
            "threshold_Aps": T_A0_ANALOG_APS,
            "equivalent": pure_quad_equivalent,
            "conclusion": (
                "pure-quadratic recorded as the empirical reduced form (a "
                "not identified by this objective)" if pure_quad_equivalent
                else "a-free and pure-quadratic distinguishable: a carries "
                "real information over these windows"
            ),
        },
        "best_variant": best_variant,
        "t_form": t_form,
        "bundles_written": bundle_paths,
        "verdict": verdict,
        "presets_rewired": False,
        "tier1_started": False,
        "vmi": "pending (post-Tier-1 final arbiter, regardless of form outcome)",
        "methodological_status": (
            "model selection on seen data under a pre-registered protocol "
            "(the §9 cross-case axis was spent); NOT fresh held-out validation"
        ),
    }
    dump_json(family_record, OUT_ROOT / "verdict_linear_quadratic.json")
    update_form_comparison_verdict(FAMILY, family_record)

    print("  NOTE: presets are NOT re-wired; Tier 1 not started; VMI stays "
          "pending/post-Tier-1; model selection on seen data (§10.4).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
