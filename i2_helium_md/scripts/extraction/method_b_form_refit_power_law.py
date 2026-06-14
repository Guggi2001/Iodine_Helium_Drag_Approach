"""Method-B §10 form-discrimination driver -- ``power_law`` family.

Sibling of the frozen §9 ``method_b_shared_refit.py``. Runs the ``power_law``
family (``F = g*C*|v|^n``, ``gamma = g*C*|v|^(n-1)``) through the same
shared-form joint-refit machinery against the pure-cubic incumbent. The free
exponent ``n`` **nests both incumbents** -- ``power_law(n=2, C=c)`` is the
pure-quadratic law and ``power_law(n=3, C=b)`` is the pure-cubic law -- so the
fitted ``n_hat`` (with its sensitivity half-width) is the direct measurement
of exponent identifiability, the sharpest single number this phase produces
(METHOD_B §10.2/§10.3).

The optimizer works in the locked **pivot parameterization** ``(gamma_ref, n)``
with ``gamma_ref = C*v_ref**(n-1)`` at ``V_REF_APS = 3.0`` (§10.4.1) so the
``log C ~ const - n*log v_bar`` matching ridge is axis-aligned and the
``n_hat`` half-width is meaningful; the stamped bundle records raw ``{C, n}``
(the §10.3 closed form -- loader contract unchanged).

Stages (same order/discipline as §9):

* **Stage-1 analog** (recorded, NON-gating, honestly weakened): fit
  ``pl_shared_3param`` on **18 A only**, score the untouched 9 A prediction vs
  the same 0.45 A/ps number (comparability with pure-cubic's 0.2685).
  Record-only ``stage1_analog_power_law.json`` -- no bundle.
* **Stage-2 shared joint fit**: ``pl_shared_3param`` only.

Verdict scoring (all pre-registered, LOCKED §10.4.1 in ``form_phase_common``):
the reused §9.4 Stage-2 bands (gating) plus the two-threshold ``T_form``
classification against the incumbent. Escalation only if the family beats the
incumbent beyond ``T_FORM_BETTER`` AND passes the bands; otherwise record and
stop. **Model selection on seen data; presets NOT re-wired.**

DRY_RUN: joint objective twice at the anchor point + bitwise-determinism
assert + timing, no writes.

Usage: edit USER SETTINGS, then::

    python scripts/extraction/method_b_form_refit_power_law.py
"""

from __future__ import annotations

import time

from form_phase_common import (
    ANCHOR_PROVENANCE,
    CASES,
    DT_ION_PS,
    ION_TIME_PS,
    N,
    N_BOUNDS,
    OUT_ROOT,
    PL_C0,
    PL_N0,
    PROJECT_ROOT,
    S1_ANALOG_RMSE_MAX_APS,
    SEED,
    UNCERTAINTY_MODEL,
    V_REF_APS,
    band_checks,
    classify_t_form,
    dump_json,
    fit_is_trapped,
    fit_summary,
    stamp,
    update_form_comparison_verdict,
)
from i2_helium_md.extraction.trajectory_matching import (
    PL_SHARED_3PARAM,
    build_case_setup,
    evaluate_form_objective,
    evaluate_joint_form_objective,
    fit_form_trajectory_matching,
    fit_shared_form_trajectory_matching,
    form_sensitivity_halfwidths,
    write_shared_form_fit_parameters,
)
from i2_helium_md.physics.drag import POWER_LAW

# =============================================================================
# USER SETTINGS
# =============================================================================
DRY_RUN = False          # True: 1 joint eval x2 + determinism assert, no writes
MAXFEV_PER_START = 400   # §9 comparability

FAMILY = "power_law"
FORM = POWER_LAW
VARIANT = PL_SHARED_3PARAM
# Pre-registered anchors (§10.4.1 LOCKED constants; NEVER setup.a0/b0).
ANCHORS = {"C0": PL_C0, "n0": PL_N0, "v_ref_Aps": V_REF_APS}


def dry_run(setups) -> None:
    coeffs = {"C": PL_C0, "n": PL_N0}
    e_by_case = {c: 0.2 for c in CASES}
    t0 = time.perf_counter()
    r1 = evaluate_joint_form_objective(FORM, coeffs, e_by_case, setups)
    t1 = time.perf_counter()
    r2 = evaluate_joint_form_objective(FORM, coeffs, e_by_case, setups)
    print(f"  joint eval at (C0={PL_C0:.4f}, n0={PL_N0:.4f}, E=0.2 eV), "
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
    print("===== Method-B §10 form refit: power_law =====")
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

    print(f"  anchors (§10.4.1 LOCKED, 18A-only): C0={PL_C0:.4f}, "
          f"n0={PL_N0:.4f}, v_ref={V_REF_APS:.2f} A/ps; n_bounds={N_BOUNDS}")

    if DRY_RUN:
        dry_run(setups)
        return 0

    # =====================================================================
    # Stage 1 analog -- fit pl_shared_3param on 18A only, predict 9A
    # =====================================================================
    print("  --- Stage 1 analog: fit pl_shared_3param on 18A, predict 9A ---")
    s1fit = fit_form_trajectory_matching(
        setups["18A"], variant=VARIANT, anchors=ANCHORS,
        n_bounds=N_BOUNDS, maxfev_per_start=MAXFEV_PER_START,
    )
    s1_coeffs = dict(s1fit.best.coefficients)
    s1_e_bind = s1fit.best.e_bind_eV["18A"]
    print(f"    18A-only fit: C={s1_coeffs['C']:.4f}, n={s1_coeffs['n']:.4f}, "
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
            "variant": VARIANT,
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
    dump_json(stage1_record, OUT_ROOT / "stage1_analog_power_law.json")

    # =====================================================================
    # Stage 2 -- shared joint fit (pl_shared_3param)
    # =====================================================================
    print(f"  --- Stage 2: {VARIANT} joint fit (pivot parameterization) ---")
    t0 = time.perf_counter()
    fit = fit_shared_form_trajectory_matching(
        setups, variant=VARIANT, anchors=ANCHORS,
        n_bounds=N_BOUNDS, maxfev_per_start=MAXFEV_PER_START,
    )
    t1 = time.perf_counter()
    best = fit.best
    print(f"    done in {t1 - t0:.0f} s, {fit.n_evaluations} joint evals "
          f"(x2 ion runs each), converged={fit.converged}")
    for i, r in enumerate(fit.starts):
        print(f"      start {i}: {r.coefficients}, "
              f"E={r.e_bind_eV[fit.cases[0]]:.4f}, "
              f"objective={r.objective:.4f}")
    print(f"    BEST: C={best.coefficients['C']:.4f}, "
          f"n={best.coefficients['n']:.4f}, "
          f"E_bind={best.e_bind_eV[fit.cases[0]]:.4f} eV, "
          f"mean objective={best.objective:.4f} A/ps")
    for c, r in best.per_case.items():
        print(f"      {c}: rmse={r.rmse_Aps:.4f} A/ps, "
              f"escape={r.escape_fraction:.3f}")
    n_hat = best.coefficients["n"]
    print(f"    n_hat = {n_hat:.4f} (nesting points: 2 = pure-quadratic, "
          "3 = pure-cubic incumbent)")

    # --- §9.4 Stage-2 band checks (gating) ---
    checks = band_checks(best)
    print("  bands [pl_shared_3param]: "
          + ", ".join(
              f"{k} {'PASS' if c['pass'] else 'FAIL'}"
              for k, c in checks.items() if k != "all_pass"
          )
          + f" -> {'ALL PASS' if checks['all_pass'] else 'FAIL'}")

    # --- T_form classification (single variant => family best) ---
    bands_pass = checks["all_pass"]
    t_form = classify_t_form(best.objective, bands_pass)
    print(f"  T_form: objective {best.objective:.6f} A/ps; "
          f"Delta = {t_form['delta_objective_Aps']:+.6f} A/ps -> "
          f"zone={t_form['zone']}, escalate={t_form['escalate_to_user']}")
    print(f"    {t_form['interpretation']}")

    # --- sensitivity (pivot space) + bundle writing ---
    bundle_path = None
    n_err = None
    c_err = None
    if fit_is_trapped(fit):
        print(f"  NOT writing {VARIANT} bundle: fit traps ions "
              "(recorded in the verdict).")
    else:
        print(f"  sensitivity scan ({VARIANT}, pivot space) ...")
        hw = form_sensitivity_halfwidths(setups, best, v_ref_Aps=V_REF_APS)
        # Convert the gamma_ref half-width to a C band at fixed n_hat
        # (C_err = gamma_ref_err / v_ref**(n_hat-1); §10.4.1 / task step 6).
        c_err = hw["gamma_ref"] / V_REF_APS ** (n_hat - 1.0)
        n_err = hw["n"]
        print(f"    n_err half-width = {n_err:.4f} (the exponent-"
              "identifiability number; a large value = degeneracy, §10.4.1)")
        fit.coeff_errs = {"C": c_err, "n": n_err}
        fit.e_bind_err_eV = hw["e_bind_eV"]
        fit.uncertainty_model = UNCERTAINTY_MODEL
        out = write_shared_form_fit_parameters(
            fit, OUT_ROOT / VARIANT, anchor_provenance=ANCHOR_PROVENANCE,
        )
        bundle_path = str(out.relative_to(PROJECT_ROOT).as_posix())
        print(f"  wrote {out}")

    # --- verdict string ---
    if t_form["escalate_to_user"]:
        verdict = (
            f"ESCALATE: {FAMILY} beats the pure-cubic incumbent beyond "
            "T_FORM_BETTER AND passes the §9.4 bands -> a competing "
            "production candidate. User decision needed before any preset "
            "re-wiring."
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
            f"recorded (n_hat={n_hat:.3f}); pure-cubic stays the candidate; "
            "VMI post-Tier-1 arbitrates."
        )
    else:
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
        "variants": [VARIANT],
        "anchors": {**ANCHORS, "provenance": ANCHOR_PROVENANCE},
        "stage1_analog": stage1_record,
        "stage2": {VARIANT: {**fit_summary(fit), "bands": checks}},
        "exponent": {
            "n_hat": n_hat,
            "n_err_halfwidth": n_err,
            "C_err": c_err,
            "nesting_note": (
                "n=2 -> pure-quadratic (Method-A n~+2), n=3 -> pure-cubic "
                "incumbent; the n_hat half-width measures exponent "
                "identifiability (§10.2/§10.4.1: a large half-width is a "
                "live expectation -- the degeneracy is itself the finding)"
            ),
        },
        "best_variant": VARIANT,
        "t_form": t_form,
        "bundles_written": {VARIANT: bundle_path} if bundle_path else {},
        "verdict": verdict,
        "presets_rewired": False,
        "tier1_started": False,
        "vmi": "pending (post-Tier-1 final arbiter, regardless of form outcome)",
        "methodological_status": (
            "model selection on seen data under a pre-registered protocol "
            "(the §9 cross-case axis was spent); NOT fresh held-out validation"
        ),
    }
    dump_json(family_record, OUT_ROOT / "verdict_power_law.json")
    update_form_comparison_verdict(FAMILY, family_record)

    print("  NOTE: presets are NOT re-wired; Tier 1 not started; VMI stays "
          "pending/post-Tier-1; model selection on seen data (§10.4).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
