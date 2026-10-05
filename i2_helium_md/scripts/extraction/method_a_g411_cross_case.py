"""G4.11 step 2: cross-case drag-law test (force balance vs trajectory fits).

THESIS_FINAL_GOALS.md G4.11 (pre-registered 10144c3; criteria, band, decision
rule, E_bind pairing and pipeline FROZEN -- do not re-tune). No MD.

For each family f in {linear, quadratic, linear_quadratic, cubic (control)}:

* ``L18(v)`` -- f's 18 A trajectory-matched law at s_rho = 3.5 A (G4.9 /
  G4.11 step 1 bundles), evaluated ungated through ``physics.drag``;
* force-balance drag ``F(v)`` on the G4.10-PH pipeline (June data handling:
  short references, [4.54, 8.0] / [2.67, 6.0] ps, MAD spline, 0.5 ps trim;
  projected Coulomb, well -dU/dr (r_hat . v_hat), F / rho_hat at 3.5 A,
  m_eff), with the well at f's OWN 18 A E_bind (E_bind pairing);
* C-X (main): ``s_X = sum(F9 L18) / sum(L18^2)`` on the 9 A points;
* C-M (secondary): ``s_M``, same on the 18 A points;
* pass: ``|s - 1| <= 0.19``; f is *consistent* if C-X and C-M pass;
* C-A (recorded, single-coefficient families): force-balance coefficient
  9 A / 18 A; shape (recorded): ``RMS((F9 - s_X L18) / mean(F9))``.

Decision rule: only the cubic consistent -> cubic supported by the TDDFT
data; several consistent -> not discriminated, cubic stays by decision;
cubic not consistent -> open as a decision. The verdict is on the
``primary`` row; the G4.10-PH sensitivity rows are recorded.

Oracle first: with the cubic bundle the pipeline must reproduce every
committed G4.10-PH cubic b (``posthoc_june_pipeline/fits.csv``) bit-exactly.

Writes ``data/reference/drag/force_balance_srho3p5A/g411_cross_case/``:
``scores.csv`` and ``verdict.json``.

Usage::

    python scripts/extraction/method_a_g411_cross_case.py
"""

from __future__ import annotations

import csv
import datetime
import json
from pathlib import Path
import sys

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from i2_helium_md.extraction.force_balance import (  # noqa: E402
    bulk_law_force,
    cross_case_scale,
    fit_drag_form,
    shape_residual_rms,
)
from i2_helium_md.extraction.trajectory_matching import _git_branch  # noqa: E402
from i2_helium_md.presets import (  # noqa: E402
    load_drag_coefficients,
    single_pulse_N2000_18Angst_drag,
)
from scripts.extraction.method_a_force_balance_g410_posthoc import (  # noqa: E402
    OUT_DIR as G410PH_DIR,
    _trace,
    _variants,
    load_case_data,
)

# =============================================================================
# Pre-registered G4.11 constants (FROZEN 10144c3; do not re-tune)
# =============================================================================
BAND = 0.19   # |s - 1| <= BAND; = 2 * 0.2127 / 2.2286 (the G4.10 band), rounded
PRIMARY = "primary"
# G4.10-PH rows that G4.11 records (the decomposition rows are not part of it).
SENSITIVITY_PREFIX_EXCLUDED = "decomp_"

# family -> (18 A trajectory-matching bundle folder, single-coefficient FB form)
TM_ROOT = REPO_ROOT / "data" / "reference" / "drag" / "18A" / "trajectory_matching_srho3p5A"
FAMILIES = {
    "linear": ("linear", "a"),
    "quadratic": ("quadratic", "c"),
    "linear_quadratic": ("linear_quadratic", None),
    "cubic": ("cubic", "b"),
}
CONTROL = "cubic"
OUT_DIR = G410PH_DIR.parent / "g411_cross_case"


def _law_label(coeffs) -> str:
    return ", ".join(f"{k} {v:.4f}" for k, v in coeffs.coefficients.items())


def score_row(family: str, coeffs, fb_by_case: dict) -> dict:
    """C-X / C-M / C-A / shape of one family under one pipeline row."""
    form, coef_name = FAMILIES[family]
    out = {"family": family}
    pts = {c: (fb.speed_Aps[fb.interior], fb.F_drag_bulk[fb.interior])
           for c, fb in fb_by_case.items()}
    for case, key in (("9A", "X"), ("18A", "M")):
        if case not in pts:
            out.update({f"s_{key}": None, f"pass_{key}": None, f"shape_{key}": None})
            continue
        v, F = pts[case]
        L = bulk_law_force(v, coeffs)
        s = cross_case_scale(F, L)
        out.update({f"s_{key}": s, f"pass_{key}": bool(abs(s - 1.0) <= BAND),
                    f"shape_{key}": shape_residual_rms(F, L, s)})
        out[f"n_{case}"] = int(v.size)
        out[f"v_range_{case}"] = [float(v.min()), float(v.max())]
    for case, (v, F) in pts.items():
        fit = fit_drag_form(v, F, form)
        out[f"fb_{case}"] = {"a": fit.a, "b": fit.b, "c": fit.c,
                             "rmse_force": fit.rmse_force}
    if coef_name is not None and {"9A", "18A"} <= set(pts):
        out["C_A_ratio_9A_over_18A"] = out["fb_9A"][coef_name] / out["fb_18A"][coef_name]
    else:
        out["C_A_ratio_9A_over_18A"] = None
    out["consistent"] = (bool(out["pass_X"] and out["pass_M"])
                         if out["pass_X"] is not None and out["pass_M"] is not None
                         else None)
    return out


def decide(consistent: dict[str, bool]) -> dict:
    """The pre-registered G4.11 decision rule on the primary row."""
    passing = sorted(f for f, ok in consistent.items() if ok)
    if not consistent[CONTROL]:
        outcome, reading = "open_decision", (
            "cubic not consistent -> open as a decision")
    elif passing == [CONTROL]:
        outcome, reading = "cubic_supported", (
            "only the cubic consistent -> the cubic form is supported by the "
            "TDDFT data (not only chosen)")
    else:
        outcome, reading = "not_discriminated", (
            f"several consistent ({', '.join(passing)}) -> the TDDFT data do not "
            "discriminate those forms; the cubic stays by decision (continuity, "
            "production cap)")
    return {"consistent_families": passing, "outcome": outcome, "reading": reading}


def run_oracle(data: dict, bundle: dict) -> None:
    """The cubic bundle must reproduce the committed G4.10-PH cubic b (bit-exact)."""
    with (G410PH_DIR / "fits.csv").open(newline="") as fh:
        committed = {(r["variant"], r["case_set"]): float(r["b"])
                     for r in csv.DictReader(fh) if r["form"] == "cubic"}
    n = 0
    for name, s in _variants(bundle).items():
        cases = ("9A",) if name == "mass_anchored_9A" else ("18A", "9A")
        for c in cases:
            fb, _ = _trace(s, data[c][1], data[c][0])
            b = fit_drag_form(fb.speed_Aps[fb.interior], fb.F_drag_bulk[fb.interior], "cubic").b
            if b != committed[(name, c)]:
                raise AssertionError(f"oracle failed at {name}/{c}: {b!r} vs "
                                     f"{committed[(name, c)]!r}")
            n += 1
    print(f"  ORACLE PASS: {n} committed G4.10-PH cubic b reproduced bit-exactly "
          f"(primary 18A {committed[(PRIMARY, '18A')]:.4f}, 9A {committed[(PRIMARY, '9A')]:.4f}).")


def main() -> None:
    m_eff = float(single_pulse_N2000_18Angst_drag().m_eff_amu)
    data = load_case_data()
    bundles, laws, step1 = {}, {}, {}
    for family, (folder, _) in FAMILIES.items():
        d = TM_ROOT / folder
        bundles[family] = json.loads((d / "fit_parameters.json").read_text())
        laws[family] = load_drag_coefficients(d, expected_m_eff_amu=m_eff)
        summary = json.loads((d / "fit_summary.json").read_text())
        step1[family] = {
            "law": _law_label(laws[family]), "form": laws[family].form,
            "e_bind_eV": bundles[family]["effective_binding_energy_I_ion_eV"],
            "e_bind_err_eV": bundles[family]["effective_binding_energy_err_eV"],
            "rmse_18A_Aps": summary["per_case"]["18A"]["rmse_Aps"],
            "bands_18A_pass": summary["bands"]["all_pass"],
            "gate_A": bundles[family]["drag_gate_steepness_angstrom"],
        }
        if step1[family]["gate_A"] != 3.5:
            raise AssertionError(f"{family} bundle is not a 3.5 A-gate fit")
    run_oracle(data, bundles[CONTROL])

    rows = []
    for family in FAMILIES:
        for name, s in _variants(bundles[family]).items():
            if name.startswith(SENSITIVITY_PREFIX_EXCLUDED):
                continue
            cases = ("9A",) if name == "mass_anchored_9A" else ("18A", "9A")
            fbs = {c: _trace(s, data[c][1], data[c][0])[0] for c in cases}
            rows.append({"row": name, "e_bind_eV": s["binding_energy_eV"],
                         **score_row(family, laws[family], fbs)})

    primary = {r["family"]: r for r in rows if r["row"] == PRIMARY}
    decision = decide({f: primary[f]["consistent"] for f in FAMILIES})
    robustness = {
        f: {key: [r["row"] for r in rows if r["family"] == f and r[f"pass_{key}"] is not None
                  and r[f"pass_{key}"] != primary[f][f"pass_{key}"]]
            for key in ("X", "M")}
        for f in FAMILIES
    }
    cubic_rel = 2.0 * bundles[CONTROL]["b_err"] / bundles[CONTROL]["b"]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    flat = []
    for r in rows:
        f = {k: v for k, v in r.items() if not isinstance(v, (dict, list))}
        for case in ("18A", "9A"):
            for k, v in (r.get(f"fb_{case}") or {}).items():
                f[f"fb_{case}_{k}"] = v
            vr = r.get(f"v_range_{case}")
            f[f"v_min_{case}"], f[f"v_max_{case}"] = (vr if vr else (None, None))
        flat.append(f)
    fields = list(dict.fromkeys(k for f in flat for k in f))
    with (OUT_DIR / "scores.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(flat)
    verdict = {
        "date": datetime.date.today().isoformat(),
        "branch": _git_branch(REPO_ROOT),
        "goal": "THESIS_FINAL_GOALS.md G4.11 cross-case drag-law test",
        "preregistration": "THESIS_FINAL_GOALS.md G4.11, committed 10144c3 before any run",
        "band": BAND,
        "band_definition": {"cubic_2hw_over_b": cubic_rel,
                            "note": "frozen value 0.19 applied to every family"},
        "pipeline": "G4.10-PH (June Method-A data handling + corrections)",
        "e_bind_pairing": "force balance for family f at f's own 18 A trajectory-fit E_bind",
        "step1_laws": step1,
        "primary": {f: {k: primary[f][k] for k in (
            "s_X", "pass_X", "s_M", "pass_M", "consistent", "C_A_ratio_9A_over_18A",
            "shape_X", "shape_M", "fb_18A", "fb_9A")} for f in FAMILIES},
        "decision": decision,
        "sensitivity_rows_flipping_primary_pass": robustness,
        "sensitivity_rows": sorted({r["row"] for r in rows}),
    }
    (OUT_DIR / "verdict.json").write_text(json.dumps(verdict, indent=2))
    print(f"  wrote {OUT_DIR / 'scores.csv'} and verdict.json")
    for f in FAMILIES:
        p = primary[f]
        ca = p["C_A_ratio_9A_over_18A"]
        print(f"  {f:17s} [{step1[f]['law']}, E {step1[f]['e_bind_eV']:.4f}] "
              f"s_X {p['s_X']:.3f} ({'pass' if p['pass_X'] else 'FAIL'}), "
              f"s_M {p['s_M']:.3f} ({'pass' if p['pass_M'] else 'FAIL'}), "
              f"C-A {'n/a' if ca is None else f'{ca:.3f}'}, shape_X {p['shape_X']:.3f}"
              f" -> {'CONSISTENT' if p['consistent'] else 'not consistent'}")
        print(f"      sensitivity rows flipping the primary pass: {robustness[f]}")
    print(f"  DECISION: {decision['outcome']} -- {decision['reading']}")


if __name__ == "__main__":
    main()
