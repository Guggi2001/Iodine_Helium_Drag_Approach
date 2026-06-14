"""Method-B §10.8 per-case (single-case) form fits -- DIAGNOSTIC ONLY.

Fills the missing diagonal: individual per-case (9 A-only, 18 A-only) optimal
fits for the ``linear_quadratic`` and ``power_law`` families, parallel to what
§8 gave ``linear_cubic``. So far these two families have only the §10.7 shared
joint refits; this driver produces the per-case bundles that answer questions
the shared fit cannot (METHOD_B §10.8):

* per-case ``n_hat`` for ``power_law`` -- does 9 A alone vs 18 A alone prefer
  ``n ~ 2`` or ``n ~ 3``, and how poorly is ``n`` identified WITHIN a single
  curve (the §10.4.1-predicted large ``n_err``)?
* does each case independently collapse to pure-quadratic (``a -> 0``) the way
  the shared ``lq`` fit did?

**Methodological status (stamped in every bundle).** These are
calibrated-not-validated, single-case fits: for these forms the §9 cross-case
held-out axis is already spent, and a per-case fit trivially matches its own
curve (the §5/§8 warning); 9 A additionally carries the transverse-contamination
flag. They are DIAGNOSTIC artifacts only, NEVER wired to presets -- the incumbent
stays ``shared_pure_cubic``. No outcome here reopens the form question (the
2026-06-14 strictly-documentary decision).

Sibling of ``method_b_form_refit_*.py``; it reuses ``form_phase_common`` (LOCKED
§10.4.1 anchors, conditioning-only) and the form-generic single-case engine
``fit_form_trajectory_matching``. It does NOT generalize the §8
``method_b_extraction.py`` (which is wired to ``linear_cubic`` three ways:
fit-fn/param-vector, preset-derived anchors, and the linear_cubic writer schema).

DRY_RUN: one single-case form objective eval x2 + bitwise-determinism assert
per case + timing, no writes (mirrors the existing drivers).

Usage: edit USER SETTINGS, then::

    python scripts/extraction/method_b_per_case_form.py
"""

from __future__ import annotations

import time

from form_phase_common import (
    A0,
    ANCHOR_PROVENANCE,
    C0_LQ,
    DT_ION_PS,
    ION_TIME_PS,
    N,
    N_BOUNDS,
    OUT_ROOT,
    PL_C0,
    PL_N0,
    PROJECT_ROOT,
    SEED,
    UNCERTAINTY_MODEL,
    V_REF_APS,
    dump_json,
    fit_is_trapped,
    stamp,
)
from i2_helium_md.extraction.trajectory_matching import (
    LQ_SHARED_3PARAM,
    LQ_SHARED_PURE_QUADRATIC,
    PL_SHARED_3PARAM,
    build_case_setup,
    evaluate_form_objective,
    fit_form_trajectory_matching,
    form_sensitivity_halfwidths,
    write_per_case_form_fit_parameters,
)
from i2_helium_md.physics.drag import LINEAR_QUADRATIC, POWER_LAW
from i2_helium_md.presets import REFERENCE_DRAG_ROOT

# =============================================================================
# USER SETTINGS
# =============================================================================
DRY_RUN = False          # True: 1 single-case eval x2 + determinism, no writes
MAXFEV_PER_START = 400   # §9 comparability

CASES = ("18A", "9A")
VARIANTS = (LQ_SHARED_3PARAM, LQ_SHARED_PURE_QUADRATIC, PL_SHARED_3PARAM)

# Pre-registered anchors per variant (§10.4.1 LOCKED constants; NEVER setup
# a0/b0 -- the re-wired presets carry a0 = 0 and must not condition a fit).
VARIANT_ANCHORS = {
    LQ_SHARED_3PARAM: {"a0": A0, "c0": C0_LQ},
    LQ_SHARED_PURE_QUADRATIC: {"a0": A0, "c0": C0_LQ},
    PL_SHARED_3PARAM: {"C0": PL_C0, "n0": PL_N0, "v_ref_Aps": V_REF_APS},
}
_IS_POWER_LAW = {
    LQ_SHARED_3PARAM: False,
    LQ_SHARED_PURE_QUADRATIC: False,
    PL_SHARED_3PARAM: True,
}


def dry_run(setups) -> None:
    coeffs = {"a": A0, "c": C0_LQ}
    for case, setup in setups.items():
        t0 = time.perf_counter()
        r1 = evaluate_form_objective(LINEAR_QUADRATIC, coeffs, 0.2, setup)
        t1 = time.perf_counter()
        r2 = evaluate_form_objective(LINEAR_QUADRATIC, coeffs, 0.2, setup)
        print(f"  {case}: single-case eval at (a0={A0:.4f}, c0={C0_LQ:.4f}, "
              f"E=0.2 eV), {t1 - t0:.2f} s/eval (= 1 ion run): "
              f"rmse {r1.rmse_Aps:.4f} escape {r1.escape_fraction:.3f}")
        if r1.objective != r2.objective:
            raise AssertionError(
                f"{case} single-case objective is NOT deterministic: "
                f"{r1.objective!r} vs {r2.objective!r} -- stochasticity has "
                "leaked into the ion drag path"
            )
        print("    repeated evaluation bitwise identical -> deterministic OK")
    print("  DRY_RUN: nothing written.")


def _fit_one(case, setup, variant):
    """Run one (case, variant) fit; return (record_dict, bundle_path | None).

    A trapped fit is recorded as a skip (not a failure, §10.8): no bundle is
    written, but the outcome is kept in the summary.
    """
    anchors = VARIANT_ANCHORS[variant]
    is_pl = _IS_POWER_LAW[variant]
    transverse = case == "9A"

    fit = fit_form_trajectory_matching(
        setup, variant=variant, anchors=anchors,
        n_bounds=N_BOUNDS, maxfev_per_start=MAXFEV_PER_START,
    )
    best = fit.best
    r = best.per_case[case]
    coeffs = dict(best.coefficients)
    e_bind = best.e_bind_eV[case]

    base = {
        "case": case,
        "variant": variant,
        "form": fit.form,
        "coefficients": coeffs,
        "e_bind_eV": e_bind,
        "objective_rmse_Aps": r.rmse_Aps,
        "escape_fraction": r.escape_fraction,
        "num_overlap_points": r.num_overlap_points,
        "converged": fit.converged,
        "transverse_contaminated": transverse,
    }

    if fit_is_trapped(fit):
        print(f"    {case}/{variant}: TRAPPED (escape "
              f"{r.escape_fraction:.3f}) -> SKIP, no bundle [recorded]")
        return {**base, "status": "skipped_trapped", "bundle": None}, None

    # --- sensitivity-only band (seed sweep omitted per §8) ---
    if is_pl:
        hw = form_sensitivity_halfwidths(
            {case: setup}, best, v_ref_Aps=V_REF_APS
        )
        n_hat = coeffs["n"]
        # Fixed-n partial band: C_err = gamma_ref_err / v_ref**(n_hat-1)
        # (§10.4.1; NOT a marginal uncertainty, see the bundle pivot note).
        c_err = hw["gamma_ref"] / V_REF_APS ** (n_hat - 1.0)
        fit.coeff_errs = {"C": c_err, "n": hw["n"]}
        print(f"    {case}/{variant}: C={coeffs['C']:.4f}, n_hat={n_hat:.4f} "
              f"+/- {hw['n']:.4f} (exponent identifiability), "
              f"E={e_bind:.4f} eV, rmse={r.rmse_Aps:.4f} A/ps")
    else:
        hw = form_sensitivity_halfwidths({case: setup}, best)
        fit.coeff_errs = {"a": hw["a"], "c": hw["c"]}
        print(f"    {case}/{variant}: a={coeffs['a']:.4f} +/- {hw['a']:.4f} "
              f"(collapse check), c={coeffs['c']:.4f} +/- {hw['c']:.4f}, "
              f"E={e_bind:.4f} eV, rmse={r.rmse_Aps:.4f} A/ps")
    fit.e_bind_err_eV = hw["e_bind_eV"]
    fit.uncertainty_model = UNCERTAINTY_MODEL

    out_dir = REFERENCE_DRAG_ROOT / case / "trajectory_matching" / variant
    out = write_per_case_form_fit_parameters(
        fit, out_dir, transverse_contaminated=transverse,
        anchor_provenance=ANCHOR_PROVENANCE,
        v_ref_Aps=V_REF_APS if is_pl else None,
    )
    bundle_path = str(out.relative_to(PROJECT_ROOT).as_posix())
    print(f"      wrote {out}")
    return (
        {**base, "status": "written", "coeff_errs": dict(fit.coeff_errs),
         "e_bind_err_eV": fit.e_bind_err_eV, "bundle": bundle_path},
        bundle_path,
    )


def main() -> int:
    print("===== Method-B §10.8 per-case form fits (DIAGNOSTIC) =====")
    print(f"  N={N}, seed={SEED}, ion {ION_TIME_PS} ps @ dt={DT_ION_PS} ps")
    print("  NOT preset-wired; incumbent stays shared_pure_cubic; strictly "
          "documentary (2026-06-14).")

    # --- per-case contexts (one neutral run per case, reused across variants)
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

    if DRY_RUN:
        dry_run(setups)
        return 0

    results = []
    for case in CASES:
        print(f"  --- case {case} ---")
        for variant in VARIANTS:
            t0 = time.perf_counter()
            record, _ = _fit_one(case, setups[case], variant)
            record["wall_s"] = time.perf_counter() - t0
            results.append(record)

    summary = {
        **stamp(),
        "phase": "METHOD_B §10.8 per-case alternative-form fits (diagnostic)",
        "status": (
            "calibrated-not-validated single-case fits; cross-case axis spent "
            "for these forms; NEVER preset-wired; incumbent shared_pure_cubic "
            "unchanged; strictly documentary (no outcome reopens the form "
            "question, 2026-06-14)"
        ),
        "settings": {
            "N": N, "seed": SEED, "ion_time_ps": ION_TIME_PS,
            "dt_ion_ps": DT_ION_PS, "maxfev_per_start": MAXFEV_PER_START,
            "n_bounds": list(N_BOUNDS), "v_ref_Aps": V_REF_APS,
            "uncertainty_model": UNCERTAINTY_MODEL,
        },
        "anchors": {**VARIANT_ANCHORS, "provenance": ANCHOR_PROVENANCE},
        "results": results,
    }
    dump_json(summary, OUT_ROOT / "per_case_form_summary.json")

    n_written = sum(1 for r in results if r["status"] == "written")
    print(f"  done: {n_written}/{len(results)} bundles written "
          f"({len(results) - n_written} skipped-trapped).")
    print("  NOTE: presets NOT re-wired; Tier 1 not started; diagnostic only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
