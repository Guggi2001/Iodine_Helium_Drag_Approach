"""G4.11 step 1: 18 A-only trajectory-matching fits of the quadratic families.

THESIS_FINAL_GOALS.md G4.11 (pre-registered 10144c3; do NOT re-tune). The
cross-case drag-law test needs every family's 18 A trajectory-matched law at
the sourced density width s_rho = 3.5 A. ``linear`` and ``cubic`` are reused
from G4.9 (``method_b_18A_only_refit.py``); this driver adds

* ``quadratic``        -- ``lq_shared_pure_quadratic`` (``F = c v^2``),
* ``linear_quadratic`` -- ``lq_shared_3param``         (``F = a v + c v^2``),

on 18 A only, with the G4.9 engine settings (N = 50, seed 20260604, 20 ps at
0.01 ps, the long smoothed reference), the LOCKED §10.4.1 anchors
(a0 14.5556, c0 11.0161; ``form_phase_common``) and sensitivity-only
half-widths.

Bands (as G4.9; a failing family is reported and still scored in step 2):
18 A RMSE <= 0.19 A/ps, escape fraction == 1.

Modes (USER SETTINGS ``MODE``; optional argv override ``python ... fit quadratic``):

* ``"oracle"`` -- (i) the committed G4.9 cubic point reproduces its 18 A RMSE
  (0.1092) at the 3.5 A gate; (ii) the June §10.8 per-case quadratic point
  (c 11.27, E_bind 0.060) reproduces its 18 A RMSE (0.0993) at the 14.2 A
  gate. Bit-exact, writes nothing.
* ``"dry_run"`` -- one 18 A form evaluation twice at 3.5 A: determinism + timing.
* ``"fit"`` -- fit each family in ``FAMILIES``; writes
  ``<family>/fit_summary.json`` and, unless trapped, the per-case form bundle
  ``<family>/fit_parameters.json`` (gate-stamped; diagnostic, not preset-wired).

Output: ``data/reference/drag/18A/trajectory_matching_srho3p5A/{quadratic,
linear_quadratic}/`` -- never the 14.2 A ``18A/trajectory_matching/`` folders.

Usage: edit USER SETTINGS, then::

    python scripts/extraction/method_b_g411_lq_18A_refit.py
"""

from __future__ import annotations

from pathlib import Path
import sys
import time

# =============================================================================
# USER SETTINGS
# =============================================================================
MODE = "oracle"   # "oracle" | "dry_run" | "fit"
FAMILIES = ("quadratic", "linear_quadratic")

# --- Pre-registered G4.9 / G4.11 bands (FIXED; do not re-tune) ---------------
RMSE_18A_MAX_APS = 0.19
# Escape fraction must be exactly 1.0.

MAXFEV_PER_START = 400   # G4.9 / §9 comparability
REFERENCE_FILENAME = "cleaned_data_long.csv"   # the G4.9 18 A reference

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.extraction.trajectory_matching import (  # noqa: E402
    LQ_SHARED_3PARAM,
    LQ_SHARED_PURE_QUADRATIC,
    build_case_setup,
    evaluate_form_objective,
    evaluate_objective,
    fit_form_trajectory_matching,
    form_sensitivity_halfwidths,
    write_per_case_form_fit_parameters,
)
from i2_helium_md.physics.drag import LINEAR_QUADRATIC  # noqa: E402
from i2_helium_md.presets import DRAG_BUNDLE_CATALOG, REFERENCE_DRAG_ROOT  # noqa: E402
from scripts.extraction.form_phase_common import (  # noqa: E402
    A0,
    ANCHOR_PROVENANCE,
    C0_LQ,
    DT_ION_PS,
    ION_TIME_PS,
    N,
    SEED,
    UNCERTAINTY_MODEL,
    dump_json,
    fit_is_trapped,
    load_json,
    stamp,
)
from scripts.extraction.method_b_srho3p5_refit import GATE_OVERRIDES  # noqa: E402

CASE = "18A"
OUT_ROOT = REFERENCE_DRAG_ROOT / CASE / "trajectory_matching_srho3p5A"
LEGACY_GATE_ROOT = REFERENCE_DRAG_ROOT / CASE / "trajectory_matching"
FAMILY_VARIANT = {
    "quadratic": LQ_SHARED_PURE_QUADRATIC,
    "linear_quadratic": LQ_SHARED_3PARAM,
}
ANCHORS = {"a0": A0, "c0": C0_LQ}   # LOCKED §10.4.1 (both lq variants)
HELD_OUT_AXIS = "9A (G4.11 C-X: force-balance drag, no trajectory fit; G4.9)"

# Oracle points (committed artifacts, read at run time -- never retyped).
G49_CUBIC_SUMMARY = OUT_ROOT / "cubic" / "fit_summary.json"
JUNE_QUADRATIC_BUNDLE = LEGACY_GATE_ROOT / LQ_SHARED_PURE_QUADRATIC / "fit_parameters.json"


def bundle_id(family: str) -> str:
    """Catalog id of a family's bundle (must map to ``OUT_ROOT / family``)."""
    bid = f"srho3p5A18A_{family}"
    if bid not in DRAG_BUNDLE_CATALOG:
        raise ValueError(f"bundle id {bid!r} is not in DRAG_BUNDLE_CATALOG")
    if (REFERENCE_DRAG_ROOT / DRAG_BUNDLE_CATALOG[bid]).resolve() != (
        OUT_ROOT / family
    ).resolve():
        raise ValueError(f"catalog maps {bid!r} elsewhere than {OUT_ROOT / family}")
    return bid


def band_checks_18a(rmse_Aps: float, escape_fraction: float) -> dict:
    """Apply the pre-registered G4.9 bands (reused by G4.11) to one 18 A score."""
    checks = {
        "RMSE_18A": {"value_Aps": rmse_Aps, "band_max_Aps": RMSE_18A_MAX_APS,
                     "pass": bool(rmse_Aps <= RMSE_18A_MAX_APS)},
        "ESCAPE_18A": {"value": escape_fraction, "band": 1.0,
                       "pass": bool(escape_fraction == 1.0)},
    }
    checks["all_pass"] = all(c["pass"] for c in checks.values())
    return checks


def _setup(gate_overrides: dict | None):
    label = gate_overrides or "the legacy 14.2 A gate"
    print(f"  building {CASE} case setup (one neutral run) at {label} ...")
    return build_case_setup(
        CASE, num_molecules=N, seed=SEED, ion_time_ps=ION_TIME_PS,
        dt_ion_ps=DT_ION_PS, reference_filename=REFERENCE_FILENAME,
        cfg_overrides=gate_overrides,
    )


def _check(label: str, got: float, committed: float) -> None:
    ok = got == committed
    print(f"  {label}: rmse {got!r} vs committed {committed!r} -> "
          f"{'BIT-EXACT' if ok else 'MISMATCH'}")
    if not ok:
        raise AssertionError(f"oracle failed: {label}")


def run_oracle() -> None:
    g49 = load_json(G49_CUBIC_SUMMARY)
    setup35 = _setup(GATE_OVERRIDES)
    r = evaluate_objective((g49["a"], g49["b"], g49["e_bind_eV"]), setup35)
    _check("(i) G4.9 cubic @ 3.5 A", r.rmse_Aps, g49["per_case"][CASE]["rmse_Aps"])

    june = load_json(JUNE_QUADRATIC_BUNDLE)
    if june.get("calibration_case") != CASE or june.get("variant") != LQ_SHARED_PURE_QUADRATIC:
        raise ValueError(f"{JUNE_QUADRATIC_BUNDLE} is not the June 18 A quadratic bundle")
    setup14 = _setup(None)
    r = evaluate_form_objective(
        LINEAR_QUADRATIC, {"a": june["a"], "c": june["c"]},
        june["effective_binding_energy_I_ion_eV"], setup14,
    )
    _check("(ii) June quadratic @ 14.2 A", r.rmse_Aps, june["objective_rmse_Aps"])
    print("  ORACLE PASS: the G4.11 path reproduces both committed RMSEs.")


def dry_run() -> None:
    setup = _setup(GATE_OVERRIDES)
    coeffs = {"a": A0, "c": C0_LQ}
    t0 = time.perf_counter()
    r1 = evaluate_form_objective(LINEAR_QUADRATIC, coeffs, 0.2, setup)
    t1 = time.perf_counter()
    r2 = evaluate_form_objective(LINEAR_QUADRATIC, coeffs, 0.2, setup)
    print(f"  18A form eval @ 3.5 A: {t1 - t0:.2f} s; rmse {r1.rmse_Aps:.4f}, "
          f"escape {r1.escape_fraction:.3f}")
    if r1.objective != r2.objective:
        raise AssertionError("single-case form objective is NOT deterministic")
    print("    repeated evaluation bitwise identical -> deterministic OK")


def _summary(fit, family: str) -> dict:
    r = fit.best.per_case[CASE]
    return {
        **stamp(),
        "goal": "THESIS_FINAL_GOALS.md G4.11 step 1 (18 A-only, s_rho = 3.5 A)",
        "preregistration": "THESIS_FINAL_GOALS.md G4.11, committed 10144c3 before any run",
        "family": family, "variant": fit.variant, "form": fit.form,
        "bundle_id": bundle_id(family),
        "calibration_cases": fit.cases,
        "coefficients": dict(fit.best.coefficients),
        "e_bind_eV": fit.best.e_bind_eV[CASE],
        "objective_Aps": fit.best.objective,
        "per_case": {CASE: {"rmse_Aps": r.rmse_Aps, "escape_fraction": r.escape_fraction,
                            "objective_Aps": r.objective,
                            "num_overlap_points": r.num_overlap_points}},
        "bands": band_checks_18a(r.rmse_Aps, r.escape_fraction),
        "n_evaluations": fit.n_evaluations, "converged": fit.converged,
        "start_minima": [{"coefficients": dict(s.coefficients),
                          "e_bind_eV": s.e_bind_eV[CASE], "objective_Aps": s.objective}
                         for s in fit.starts],
        "anchors": dict(fit.anchors),
        "gate": {"drag_spatial_gate": fit.drag_spatial_gate,
                 "drag_gate_steepness_angstrom": fit.gate_steepness_angstrom,
                 "potential_steepness_angstrom": fit.potential_steepness_angstrom},
        "reference_file": str(fit.reference_paths[CASE].relative_to(PROJECT_ROOT).as_posix()),
        "n_molecules": fit.n_molecules, "seed": fit.seed, "meff_amu": fit.meff_amu,
        "held_out_axis": HELD_OUT_AXIS,
    }


def run_fits(families) -> None:
    setup = _setup(GATE_OVERRIDES)
    for family in families:
        if family not in FAMILY_VARIANT:
            raise ValueError(f"unknown family {family!r}; expected {sorted(FAMILY_VARIANT)}")
        out_dir = OUT_ROOT / family
        if out_dir.resolve().is_relative_to(LEGACY_GATE_ROOT.resolve()):
            raise AssertionError(f"{out_dir} is inside the 14.2 A folders")
        t0 = time.perf_counter()
        fit = fit_form_trajectory_matching(
            setup, variant=FAMILY_VARIANT[family], anchors=ANCHORS,
            maxfev_per_start=MAXFEV_PER_START,
        )
        if fit.gate_steepness_angstrom != GATE_OVERRIDES["drag_gate_steepness"]:
            raise AssertionError(f"fit ran at gate {fit.gate_steepness_angstrom} A, not 3.5 A")
        b = fit.best
        print(f"    {family}: {time.perf_counter() - t0:.0f} s, {fit.n_evaluations} evals; "
              f"{b.coefficients}, E_bind={b.e_bind_eV[CASE]:.4f} eV, "
              f"rmse={b.per_case[CASE].rmse_Aps:.4f}")
        summary = _summary(fit, family)
        trapped = fit_is_trapped(fit)
        if not trapped:
            hw = form_sensitivity_halfwidths({CASE: setup}, b)
            fit.coeff_errs = {"a": hw["a"], "c": hw["c"]}
            fit.e_bind_err_eV = hw["e_bind_eV"]
            fit.uncertainty_model = UNCERTAINTY_MODEL
            summary["sensitivity_halfwidths"] = hw
            print(f"      half-widths {hw}")
            write_per_case_form_fit_parameters(
                fit, out_dir, transverse_contaminated=False,
                anchor_provenance=ANCHOR_PROVENANCE, bundle_id=bundle_id(family),
            )
        summary["bundle_written"] = not trapped
        dump_json(summary, out_dir / "fit_summary.json")


def main(mode: str = MODE, families=FAMILIES) -> int:
    print(f"===== G4.11 step 1: 18A-only lq fits at s_rho = 3.5 A, MODE={mode} =====")
    if mode == "oracle":
        run_oracle()
    elif mode == "dry_run":
        dry_run()
    elif mode == "fit":
        print(f"  anchors (LOCKED §10.4.1): a0={A0:.4f}, c0={C0_LQ:.4f}")
        run_fits(families)
    else:
        raise ValueError(f"unknown MODE {mode!r}")
    return 0


if __name__ == "__main__":
    # Optional command-line override: ``python ... fit [family ...]``.
    args = sys.argv[1:]
    raise SystemExit(main(args[0] if args else MODE, tuple(args[1:]) or FAMILIES))
