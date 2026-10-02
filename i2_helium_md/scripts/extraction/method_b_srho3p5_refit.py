"""Method-B Stage-2 joint refit at the sourced He density width s_rho = 3.5 A.

THESIS_FINAL_GOALS.md G4 (= atlas plan §6.8 T5). The density/drag-gate width
moves from the solvation-potential width (14.2 A, reused as a density width) to
the literature value s_rho = 3.5 A (Harms, Toennies, Dalfovo, PRB 58, 3341
(1998); 10-90 % width 1.81*s = 6.3 A). ``potential_steepness`` stays 14.2 A.
The extraction mass stays m_eff = 19 He (G3).

**Stage 2 only** (joint fit over the 18 A + 9 A TDDFT traces; no Stage-1
held-out prediction, no per-case refits). Three families, each fitting its
coefficients plus a shared E_bind on the existing ``linear_cubic`` force:

==========================  =================  ==========================
family                      force              free
==========================  =================  ==========================
``shared_linear``           F = a v            a [amu/ps], E_bind
``shared_cubic``            F = b v^3          b [amu ps/A^2], E_bind
``shared_linear_cubic``     F = a v + b v^3    a, b, E_bind
==========================  =================  ==========================

Pre-registered (G4.2, committed b729c28 before any run; do NOT re-tune):

* bands per family (METHOD_B §9.4): RMSE 18 A <= 0.19 A/ps, RMSE 9 A <= 0.45
  A/ps, escape fraction == 1 in both cases; a family failing its bands is out;
* nested-model tests: ``T_a0 = obj(cubic) - obj(linear_cubic) <= 0.005`` ->
  cubic equivalent; ``T_b0 = obj(linear) - obj(linear_cubic) <= 0.005`` ->
  linear equivalent;
* winner: exactly one boundary equivalent -> that family; neither ->
  ``shared_linear_cubic``; both -> ``shared_cubic``.

Modes (USER SETTINGS ``MODE``):

* ``"dry_run"`` -- build the setups at the configured gate, evaluate the joint
  objective twice at the anchor point, assert bitwise determinism, print the
  time per joint evaluation; writes nothing.
* ``"fit"`` -- fit each family in ``FAMILIES``; always writes
  ``<family>/fit_summary.json`` (band checks included, also for a trapped fit)
  and, when the fit permits full escape, the loadable
  ``<family>/fit_parameters.json`` (stamped ``bundle_id``, gate, source,
  ``held_out_axis``). Split ``FAMILIES`` across two processes to stay within
  the two-concurrent-runs rule.
* ``"verdict"`` -- read the three summaries, apply the rule above, run the
  ``diagnostic_4param`` localization if any family fails its bands, write
  ``verdict.json``.

Usage: edit USER SETTINGS, then::

    python scripts/extraction/method_b_srho3p5_refit.py
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
MODE = "dry_run"        # "dry_run" | "fit" | "verdict"
FAMILIES = ("shared_cubic", "shared_linear_cubic", "shared_linear")
N = 50                  # ensemble size per forward integration (Tier-0 evidence)
SEED = 20260604         # neutral-stage seed (ion drag path is RNG-free)
ION_TIME_PS = 20.0      # covers the [t*, ~14 ps] windows with escape margin
DT_ION_PS = 0.01
MAXFEV_PER_START = 400
RUN_DIAGNOSTIC_ON_FAIL = True

# G2: the sourced density width (drag gate only; potential_steepness untouched).
GATE_OVERRIDES = {"drag_spatial_gate": "erf_independent", "drag_gate_steepness": 3.5}
BUNDLE_TAG = "srho3p5A"  # catalog-id prefix + output folder suffix

# Smoothed TDDFT speed references, file under data/reference/drag/<case>/velocity_smoothed/
REFERENCE_FILENAMES = {"18A": "cleaned_data_long.csv", "9A": "cleaned_data_long.csv"}

# --- Pre-registered G4.2 bands (FIXED before any run; do not re-tune) --------
S2_RMSE_18A_MAX_APS = 0.19
S2_RMSE_9A_MAX_APS = 0.45
T_A0_APS = 0.005        # cubic equivalent:  obj(cubic)  - obj(linear_cubic)
T_B0_APS = 0.005        # linear equivalent: obj(linear) - obj(linear_cubic)
# Escape fraction must be exactly 1.0 in both cases.

DENSITY_WIDTH_SOURCE = (
    "Harms, Toennies, Dalfovo, PRB 58, 3341 (1998): 4He droplet 10-90 % "
    "surface width 5.7 A (DFT) / 6-8 A (exp); erf-complement 10-90 % = "
    "1.81 s -> s = 3.5 A gives 6.3 A. Borrowed (atlas plan §6.8.2)."
)
HELD_OUT_AXIS = "none (stage-2 only)"

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.extraction.trajectory_matching import (  # noqa: E402
    DIAGNOSTIC_4PARAM,
    SHARED_CUBIC,
    SHARED_LINEAR,
    SHARED_LINEAR_CUBIC,
    _git_branch,
    build_case_setup,
    evaluate_joint_objective,
    fit_shared_trajectory_matching,
    joint_sensitivity_halfwidths,
    write_shared_fit_parameters,
)
from i2_helium_md.presets import DRAG_BUNDLE_CATALOG, REFERENCE_DRAG_ROOT  # noqa: E402

CASES = ("18A", "9A")
OUT_ROOT = REFERENCE_DRAG_ROOT / "shared" / f"trajectory_matching_{BUNDLE_TAG}"
ANCHOR_BUNDLE = REFERENCE_DRAG_ROOT / "18A" / "linear_and_cubic" / "fit_parameters.json"
LEGACY_BUNDLE = (
    REFERENCE_DRAG_ROOT / "shared" / "trajectory_matching" / "shared_pure_cubic"
    / "fit_parameters.json"
)
ALL_FAMILIES = (SHARED_LINEAR, SHARED_CUBIC, SHARED_LINEAR_CUBIC)
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


def bundle_id(family: str) -> str:
    """Catalog id of a family's bundle (must exist in DRAG_BUNDLE_CATALOG)."""
    bid = f"{BUNDLE_TAG}_{family}"
    if bid not in DRAG_BUNDLE_CATALOG:
        raise ValueError(f"bundle id {bid!r} is not in DRAG_BUNDLE_CATALOG")
    expected = OUT_ROOT / family
    if (REFERENCE_DRAG_ROOT / DRAG_BUNDLE_CATALOG[bid]).resolve() != expected.resolve():
        raise ValueError(
            f"catalog maps {bid!r} to {DRAG_BUNDLE_CATALOG[bid]!r}, not {expected}"
        )
    return bid


def band_checks(per_case: dict[str, dict]) -> dict:
    """Apply the pre-registered G4.2 bands to one family's per-case scores."""
    r18, r9 = per_case["18A"], per_case["9A"]
    checks = {
        "S2_RMSE_18A": {
            "value_Aps": r18["rmse_Aps"], "band_max_Aps": S2_RMSE_18A_MAX_APS,
            "pass": bool(r18["rmse_Aps"] <= S2_RMSE_18A_MAX_APS),
        },
        "S2_RMSE_9A": {
            "value_Aps": r9["rmse_Aps"], "band_max_Aps": S2_RMSE_9A_MAX_APS,
            "pass": bool(r9["rmse_Aps"] <= S2_RMSE_9A_MAX_APS),
        },
        "S2_ESCAPE": {
            "value": {"18A": r18["escape_fraction"], "9A": r9["escape_fraction"]},
            "band": 1.0,
            "pass": bool(
                r18["escape_fraction"] == 1.0 and r9["escape_fraction"] == 1.0
            ),
        },
    }
    checks["all_pass"] = all(c["pass"] for c in checks.values())
    return checks


def select_winner(objectives: dict[str, float], passes: dict[str, bool]) -> dict:
    """The pre-registered G4.2 selection rule (pure function, unit-tested).

    ``objectives`` / ``passes`` are keyed by the three family names. Returns
    the signed nested-model differences, the equivalence flags, and the
    winner (``None`` if no admissible family exists). A family failing its
    bands is out before the rule is applied; a single-term family can only
    win as "equivalent" if it passes its own bands, and ``linear_cubic`` can
    only win (neither boundary equivalent) if it passes.
    """
    lc = objectives[SHARED_LINEAR_CUBIC]
    t_a0 = objectives[SHARED_CUBIC] - lc
    t_b0 = objectives[SHARED_LINEAR] - lc
    cubic_eq = bool(t_a0 <= T_A0_APS)
    linear_eq = bool(t_b0 <= T_B0_APS)
    cubic_ok = cubic_eq and passes[SHARED_CUBIC]
    linear_ok = linear_eq and passes[SHARED_LINEAR]
    if cubic_ok:  # alone, or both equivalent -> cubic (Tier-0 continuity)
        winner, reason = SHARED_CUBIC, (
            "both boundaries equivalent -> cubic (continuity; production cap "
            "is defined on cubic)" if linear_ok else "cubic boundary equivalent"
        )
    elif linear_ok:
        winner, reason = SHARED_LINEAR, "linear boundary equivalent"
    elif passes[SHARED_LINEAR_CUBIC]:
        winner, reason = SHARED_LINEAR_CUBIC, (
            "no admissible equivalent boundary -> full linear_cubic"
        )
    else:
        winner, reason = None, "no family admissible (all band failures)"
    notes = [
        f"{name} = {val:+.6f} A/ps < 0: the linear_cubic optimizer did not "
        "reach the nested optimum (obj(lc) must be <= its boundary in exact "
        "arithmetic); the boundary counts as equivalent"
        for name, val in (("T_a0", t_a0), ("T_b0", t_b0)) if val < 0.0
    ]
    return {
        "T_a0_Aps": t_a0, "T_a0_threshold_Aps": T_A0_APS, "cubic_equivalent": cubic_eq,
        "T_b0_Aps": t_b0, "T_b0_threshold_Aps": T_B0_APS, "linear_equivalent": linear_eq,
        "passes": dict(passes),
        "winner": winner, "reason": reason, "notes": notes,
    }


def build_setups() -> dict:
    setups = {}
    for case in CASES:
        print(f"  building {case} case setup (one neutral run) at {GATE_OVERRIDES} ...")
        setups[case] = build_case_setup(
            case, num_molecules=N, seed=SEED,
            ion_time_ps=ION_TIME_PS, dt_ion_ps=DT_ION_PS,
            reference_filename=REFERENCE_FILENAMES[case],
            cfg_overrides=GATE_OVERRIDES,
        )
        w = setups[case].window
        print(f"    reference {setups[case].reference_path}, "
              f"window [{w[0]:.2f}, {w[1]:.2f}] ps")
    return setups


def anchors() -> tuple[float, float]:
    """§9.3 anchors from the 18 A Method-A bundle ONLY (atlas-oracle precedent;
    the presets now wire the a = 0 shared bundle, which must not condition a fit)."""
    ma = _load_json(ANCHOR_BUNDLE)
    return float(ma["a"]), float(ma["b"])


def dry_run(setups, anch) -> None:
    a0, b0 = anch
    e_by_case = {c: 0.2 for c in CASES}
    t0 = time.perf_counter()
    r1 = evaluate_joint_objective(a0, b0, e_by_case, setups)
    t1 = time.perf_counter()
    r2 = evaluate_joint_objective(a0, b0, e_by_case, setups)
    print(f"  joint eval at (a0={a0:.4f}, b0={b0:.4f}, E=0.2 eV): "
          f"{t1 - t0:.2f} s/joint-eval (= 2 ion runs)")
    print(f"    mean objective = {r1.objective:.6f} A/ps; per case: "
          + ", ".join(f"{c} rmse {r.rmse_Aps:.4f} escape {r.escape_fraction:.3f}"
                      for c, r in r1.per_case.items()))
    if r1.objective != r2.objective:
        raise AssertionError(
            f"joint objective is NOT deterministic ({r1.objective!r} vs "
            f"{r2.objective!r}) -- stochasticity leaked into the drag ion path"
        )
    print("    repeated evaluation bitwise identical -> deterministic OK")
    print("  DRY_RUN: nothing written.")


def _summary(fit, setups) -> dict:
    per_case = {
        c: {"rmse_Aps": r.rmse_Aps, "escape_fraction": r.escape_fraction,
            "objective_Aps": r.objective, "num_overlap_points": r.num_overlap_points}
        for c, r in fit.best.per_case.items()
    }
    return {
        **_stamp(),
        "family": fit.variant,
        "bundle_id": bundle_id(fit.variant),
        "a": fit.best.a, "b": fit.best.b,
        "e_bind_eV": fit.best.e_bind_eV[fit.cases[0]],
        "objective_mean_Aps": fit.best.objective,
        "per_case": per_case,
        "bands": band_checks(per_case),
        "n_evaluations": fit.n_evaluations,
        "converged": fit.converged,
        "start_minima": [
            {"a": r.a, "b": r.b, "e_bind_eV": r.e_bind_eV, "objective_Aps": r.objective}
            for r in fit.starts
        ],
        "anchors": {"a0": fit.anchors[0], "b0": fit.anchors[1],
                    "source": str(ANCHOR_BUNDLE.relative_to(PROJECT_ROOT).as_posix())},
        "gate": {"drag_spatial_gate": fit.drag_spatial_gate,
                 "drag_gate_steepness_angstrom": fit.gate_steepness_angstrom,
                 "potential_steepness_angstrom": fit.potential_steepness_angstrom},
        "reference_files": {
            c: str(setups[c].reference_path.relative_to(PROJECT_ROOT).as_posix())
            for c in CASES
        },
        "n_molecules": fit.n_molecules, "seed": fit.seed, "meff_amu": fit.meff_amu,
    }


def run_fits(setups, anch) -> None:
    for family in FAMILIES:
        if family not in ALL_FAMILIES:
            raise ValueError(f"unknown family {family!r}; expected {ALL_FAMILIES}")
        print(f"  --- Stage 2: {family} joint fit ---")
        t0 = time.perf_counter()
        fit = fit_shared_trajectory_matching(
            setups, variant=family, anchors=anch, maxfev_per_start=MAXFEV_PER_START,
        )
        print(f"    done in {time.perf_counter() - t0:.0f} s, {fit.n_evaluations} "
              f"joint evals, converged={fit.converged}")
        best = fit.best
        print(f"    BEST: a={best.a:.4f} amu/ps, b={best.b:.4f} amu*ps/A^2, "
              f"E_bind={best.e_bind_eV[fit.cases[0]]:.4f} eV, "
              f"objective={best.objective:.6f} A/ps")
        summary = _summary(fit, setups)
        trapped = any(r.penalty_Aps != 0.0 for r in best.per_case.values())
        if trapped:
            print(f"    NOT writing a loadable {family} bundle: fit traps ions.")
        else:
            hw_a, hw_b, hw_e = joint_sensitivity_halfwidths(setups, best)
            fit.a_err, fit.b_err, fit.e_bind_err_eV = hw_a, hw_b, hw_e
            fit.uncertainty_model = UNCERTAINTY_MODEL
            summary["sensitivity_halfwidths"] = {"a": hw_a, "b": hw_b, "e_bind_eV": hw_e}
            write_shared_fit_parameters(
                fit, OUT_ROOT / family,
                bundle_id=bundle_id(family),
                density_width_source=DENSITY_WIDTH_SOURCE,
                held_out_axis=HELD_OUT_AXIS,
            )
        summary["bundle_written"] = not trapped
        _dump_json(summary, OUT_ROOT / family / "fit_summary.json")


def run_verdict() -> None:
    summaries = {f: _load_json(OUT_ROOT / f / "fit_summary.json") for f in ALL_FAMILIES}
    objectives = {f: float(s["objective_mean_Aps"]) for f, s in summaries.items()}
    passes = {f: bool(s["bands"]["all_pass"]) for f, s in summaries.items()}
    selection = select_winner(objectives, passes)
    print(f"  T_a0 = {selection['T_a0_Aps']:+.6f} A/ps -> cubic "
          f"{'EQUIVALENT' if selection['cubic_equivalent'] else 'distinguishable'}")
    print(f"  T_b0 = {selection['T_b0_Aps']:+.6f} A/ps -> linear "
          f"{'EQUIVALENT' if selection['linear_equivalent'] else 'distinguishable'}")
    print(f"  bands: {passes}")
    print(f"  WINNER: {selection['winner']} ({selection['reason']})")

    diagnostic = None
    if not all(passes.values()) and RUN_DIAGNOSTIC_ON_FAIL:
        print("  --- a family FAILED its bands: running diagnostic_4param ---")
        setups = build_setups()
        fitd = fit_shared_trajectory_matching(
            setups, variant=DIAGNOSTIC_4PARAM, anchors=anchors(),
            maxfev_per_start=MAXFEV_PER_START,
        )
        diagnostic = {
            **_stamp(),
            "a": fitd.best.a, "b": fitd.best.b, "e_bind_eV": fitd.best.e_bind_eV,
            "objective_mean_Aps": fitd.best.objective,
            "per_case": {c: {"rmse_Aps": r.rmse_Aps, "escape_fraction": r.escape_fraction}
                         for c, r in fitd.best.per_case.items()},
            "purpose": "failure localization only (drag law vs binding); never a "
                       "production bundle (METHOD_B §9.2)",
        }
        _dump_json(diagnostic, OUT_ROOT / "diagnostic_4param.json")

    legacy = _load_json(LEGACY_BUNDLE)
    verdict = {
        **_stamp(),
        "goal": "THESIS_FINAL_GOALS.md G4 (atlas §6.8 T5): Stage-2 refit at s_rho = 3.5 A",
        "preregistration": "THESIS_FINAL_GOALS.md G4.2, committed b729c28 before any run",
        "thresholds": {
            "S2_RMSE_18A_MAX_Aps": S2_RMSE_18A_MAX_APS,
            "S2_RMSE_9A_MAX_Aps": S2_RMSE_9A_MAX_APS,
            "escape_band": 1.0, "T_a0_Aps": T_A0_APS, "T_b0_Aps": T_B0_APS,
        },
        "families": {f: {k: s[k] for k in ("bundle_id", "a", "b", "e_bind_eV",
                                             "objective_mean_Aps", "per_case", "bands",
                                             "converged", "bundle_written")}
                     for f, s in summaries.items()},
        "selection": selection,
        "diagnostic_4param": diagnostic,
        "legacy_bundle_for_comparison": {
            "path": str(LEGACY_BUNDLE.relative_to(PROJECT_ROOT).as_posix()),
            "b": legacy["b"], "e_bind_eV": legacy["effective_binding_energy_I_ion_eV"],
            "objective_mean_Aps": legacy["objective_mean_Aps"],
            "note": "fit at the 14.2 A gate; recorded for the predictions on record "
                    "(E_bind expected to fall, b shift small), not a gate",
        },
        "gate": summaries[SHARED_CUBIC]["gate"],
        "reference_files": summaries[SHARED_CUBIC]["reference_files"],
    }
    _dump_json(verdict, OUT_ROOT / "verdict.json")


def main() -> int:
    print(f"===== Method-B Stage-2 refit at {GATE_OVERRIDES} (G4), MODE={MODE} =====")
    if MODE == "verdict":
        run_verdict()
        return 0
    anch = anchors()
    print(f"  anchors (18A Method-A only): a0={anch[0]:.4f}, b0={anch[1]:.4f}")
    setups = build_setups()
    if MODE == "dry_run":
        dry_run(setups, anch)
    elif MODE == "fit":
        run_fits(setups, anch)
    else:
        raise ValueError(f"MODE must be 'dry_run', 'fit' or 'verdict', got {MODE!r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
