"""Atlas D4 Step-1 report: Method-B form table + gamma(v) diagnostics + tails.

Sensitivity-atlas stage 2a (`TIER2_SENSITIVITY_ATLAS_PLAN.md` section 6.2).
Three read-only sections, all compiled from the committed Tier-0 artifacts
under ``data/reference/drag/`` -- **no fits are run and nothing is written**
(optional CSV export excepted, off by default):

1. **Form table** -- the full Tier-0 three-mode protocol per form (shared
   joint fit / held-out 18A->9A prediction / per-case single-curve),
   coefficients + objectives + delta vs the pure-cubic incumbent, from the
   stored ``verdict*.json`` / ``stage1*.json`` / ``per_case_form_summary``
   artifacts, reused verbatim (plan section 6.2: reuse, not rerun).
2. **gamma(v) diagnostic table** -- every form's implied friction
   coefficient gamma(v) [amu/ps] at the plan's checkpoints
   v = 2 / 7.25 / 10.5 A/ps plus the calibrated band edges 2.54 / 4.95,
   computed through :func:`i2_helium_md.physics.drag.drag_gamma` (gate = 1,
   i.e. rho_hat = 1) for every realised form. The two **parked** forms
   (saturating/Pade cubic, subtractive gated cubic -- plan section 6.1) are
   evaluated with locked b and *implied* parameters (v_s from the arbitrated
   cap, v_f from the section-4ee twin window); their formulas live only here
   because no production enum realises them yet (their zero-cost entry gate).
3. **Trace-tail inspection** -- the low-v end of the smoothed 9/18 A TDDFT
   references: window-end speeds, deceleration, local effective power-law
   exponent n_loc(v) = dln(-dv/dt)/dln(v) over the tail, and a
   force-balance-style empirical gamma_emp = m_eff * (-dv/dt) / v at low v
   (qualitative -- the tail mixes drag with the exit-well climb; the
   findings doc carries the interpretation caveats).

Usage: edit USER SETTINGS, then::

    python scripts/extraction/atlas_d4_step1_report.py
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

# =============================================================================
# USER SETTINGS
# =============================================================================
V_GRID_APS = (2.0, 2.54, 4.95, 7.25, 10.5)  # plan checkpoints + band edges
V_C_APS = 7.25          # arbitrated production cap (standing point)
P_TAIL = -1.0           # arbitrated tail exponent (constant-force)
V_S_IMPLIED_APS = 7.25  # Pade v_s implied from the arbitrated cap (plan 6.1)
V_F_IMPLIED_APS = (1.5, 2.0)  # gated-form v_f from the 4ee twin window
TAIL_LAST_PS = 3.0      # tail segment: last TAIL_LAST_PS of each window
TAIL_V_MAX_APS = 3.0    # and separately, all in-window samples below this v
CSV_OUT = None          # e.g. Path("gamma_table.csv"); None = print only

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.physics.drag import (  # noqa: E402
    CAPPED_CUBIC,
    LINEAR_CUBIC,
    LINEAR_QUADRATIC,
    POWER_LAW,
    DragCoefficients,
    drag_gamma,
)
from i2_helium_md.postprocess.hedft_loader import (  # noqa: E402
    load_smoothed_speed_reference,
)
from i2_helium_md.presets import REFERENCE_DRAG_ROOT  # noqa: E402

SHARED_ROOT = REFERENCE_DRAG_ROOT / "shared" / "trajectory_matching"
CASES = ("18A", "9A")

# gate = 1 stand-in: deep inside the droplet the erf gate is exactly 1.0
_DEEP_DEPTH_A = -50.0
_GATE_STEEPNESS_A = 1.0


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


# -----------------------------------------------------------------------------
# gamma(v) evaluation
# -----------------------------------------------------------------------------
def gamma_realised(form: str, coefficients: dict, v_aps) -> np.ndarray:
    """gamma(v) [amu/ps] for a realised form at gate = 1 (rho_hat = 1)."""
    coeffs = DragCoefficients(
        form=form,
        coefficients=coefficients,
        extraction_mass_model="constant",
        extraction_mass_amu=202.953908,
        extraction_method="trajectory_matching",
    )
    v = np.asarray(v_aps, dtype=float)
    return drag_gamma(v, np.full_like(v, _DEEP_DEPTH_A), coeffs,
                      _GATE_STEEPNESS_A)


def gamma_pade(b: float, v_s: float, v_aps) -> np.ndarray:
    """Parked saturating (Pade) cubic: gamma = b*v^2 / (1 + (v/v_s)^3).

    [amu/ps]; reduces to b*v^2 for v << v_s; force asymptote b*v_s^3.
    Not a realised production form -- zero-cost entry-gate arithmetic only.
    """
    v = np.asarray(v_aps, dtype=float)
    return b * v**2 / (1.0 + (v / v_s) ** 3)


def gamma_subtractive(b: float, v_f: float, v_aps) -> np.ndarray:
    """Parked subtractive gated cubic: gamma = max(b*(v^2 - v_f^2), 0).

    [amu/ps]; the ``linear_cubic`` family with a = -b*v_f^2 < 0 plus the
    zero-clamp the production enum does not implement yet -- zero-cost
    entry-gate arithmetic only.
    """
    v = np.asarray(v_aps, dtype=float)
    return np.maximum(b * (v**2 - v_f**2), 0.0)


def build_gamma_rows() -> list[tuple[str, np.ndarray]]:
    """Assemble (label, gamma-at-V_GRID) rows from the committed artifacts."""
    v = np.asarray(V_GRID_APS, dtype=float)
    pc = _load(SHARED_ROOT / "shared_pure_cubic" / "fit_parameters.json")
    lq = _load(SHARED_ROOT / "lq_shared_3param" / "fit_parameters.json")
    pl = _load(SHARED_ROOT / "pl_shared_3param" / "fit_parameters.json")
    cub9 = _load(REFERENCE_DRAG_ROOT / "9A" / "trajectory_matching"
                 / "fit_parameters.json")
    cub18 = _load(REFERENCE_DRAG_ROOT / "18A" / "trajectory_matching"
                  / "fit_parameters.json")
    b = float(pc["b"])

    rows = [
        ("capped_cubic v_c=7.25 (production)",
         gamma_realised(CAPPED_CUBIC,
                        {"b": b, "v_c": V_C_APS, "p_tail": P_TAIL}, v)),
        ("pure cubic shared (uncapped)",
         gamma_realised(LINEAR_CUBIC, {"a": 0.0, "b": b}, v)),
        ("lin+quad shared (lq_shared_3param)",
         gamma_realised(LINEAR_QUADRATIC,
                        {"a": float(lq["a"]), "c": float(lq["c"])}, v)),
        ("power law shared (n=%.3f)" % float(pl["n"]),
         gamma_realised(POWER_LAW,
                        {"C": float(pl["C"]), "n": float(pl["n"])}, v)),
        ("lin+cubic 9A per-case",
         gamma_realised(LINEAR_CUBIC,
                        {"a": float(cub9["a"]), "b": float(cub9["b"])}, v)),
        ("lin+cubic 18A per-case",
         gamma_realised(LINEAR_CUBIC,
                        {"a": float(cub18["a"]), "b": float(cub18["b"])}, v)),
        ("PARKED Pade sat. cubic v_s=%.2f" % V_S_IMPLIED_APS,
         gamma_pade(b, V_S_IMPLIED_APS, v)),
    ]
    for v_f in V_F_IMPLIED_APS:
        rows.append(("PARKED subtractive cubic v_f=%.1f" % v_f,
                     gamma_subtractive(b, v_f, v)))
    return rows


def print_gamma_table() -> list[tuple[str, np.ndarray]]:
    rows = build_gamma_rows()
    v = np.asarray(V_GRID_APS, dtype=float)
    ref = rows[0][1]  # production capped_cubic is the comparison row
    head = "gamma(v) [amu/ps] at gate=1 (rho_hat=1)"
    print("=" * 78)
    print("SECTION 2 -- " + head)
    print("  band 2.54-4.95 = TDDFT-calibrated; 7.25 = cap; 10.5 = prod peak")
    print("-" * 78)
    print(f"{'form':38s}" + "".join(f"  v={x:<5.4g}" for x in v))
    for label, g in rows:
        print(f"{label:38s}" + "".join(f"  {x:7.2f}" for x in g))
    print("-" * 78)
    print("ratio to production capped_cubic:")
    for label, g in rows[1:]:
        with np.errstate(divide="ignore", invalid="ignore"):
            r = g / ref
        print(f"{label:38s}" + "".join(f"  {x:7.3f}" for x in r))
    print("-" * 78)
    print("cap pressure F(10.5)/F(7.25) (arbitration wanted ~1.0, p_tail=-1):")
    for label, g in rows:
        f_lo = g[list(V_GRID_APS).index(7.25)] * 7.25
        f_hi = g[list(V_GRID_APS).index(10.5)] * 10.5
        ratio = f_hi / f_lo if f_lo > 0 else float("nan")
        print(f"  {label:38s} {ratio:6.2f}")
    if CSV_OUT is not None:
        with open(CSV_OUT, "w", encoding="utf-8") as fh:
            fh.write("form," + ",".join(f"v={x}" for x in v) + "\n")
            for label, g in rows:
                fh.write(label.replace(",", ";") + ","
                         + ",".join(f"{x:.6f}" for x in g) + "\n")
        print(f"  wrote {CSV_OUT}")
    return rows


# -----------------------------------------------------------------------------
# Section 1 -- form table from the committed artifacts
# -----------------------------------------------------------------------------
def print_form_table() -> None:
    verdict = _load(SHARED_ROOT / "verdict.json")
    fc = _load(SHARED_ROOT / "form_comparison_verdict.json")
    pcs = _load(SHARED_ROOT / "per_case_form_summary.json")["results"]
    incumbent_obj = float(fc["_meta"]["incumbent_objective_Aps"])

    print("=" * 78)
    print("SECTION 1 -- Method-B three-mode form table (artifact reuse)")
    print(f"  incumbent shared_pure_cubic objective = {incumbent_obj:.5f} A/ps"
          f" ; equivalence band +-0.005 ; held-out band <= 0.45 A/ps")
    print("-" * 78)

    print("mode: SHARED joint fit (both traces, one coefficient set + E_bind)")
    s2 = verdict["stage2"]
    for key in ("shared_pure_cubic", "shared_3param"):
        r = s2[key]
        print(f"  {key:24s} a={float(r['a']):.4g} b={float(r['b']):.5g} "
              f"E_bind={float(r['e_bind_eV']['18A']):.5g} eV "
              f"obj={float(r['objective_mean_Aps']):.5f} "
              f"dObj={float(r['objective_mean_Aps']) - incumbent_obj:+.5f}")
    for fam, key in (("linear_quadratic", "lq_shared_3param"),
                     ("power_law", "pl_shared_3param")):
        r = fc[fam]["stage2"][key]
        cs = ", ".join(f"{k}={float(x):.5g}"
                       for k, x in r["coefficients"].items())
        print(f"  {key:24s} {cs} "
              f"E_bind={float(r['e_bind_eV']['18A']):.5g} eV "
              f"obj={float(r['objective_mean_Aps']):.5f} "
              f"dObj={float(r['objective_mean_Aps']) - incumbent_obj:+.5f}")
        print(f"    verdict: {fc[fam]['verdict']}")

    print("mode: HELD-OUT (fit 18A alone -> predict untouched 9A; "
          "band <= 0.45)")
    s1 = verdict["stage1"]
    print(f"  linear_cubic (sec-8 bundle)  rmse_9A={float(s1['rmse_Aps']):.4f}"
          f"  pass={s1['checks']['S1_PRED_RMSE_MAX']['pass']}")
    for fam in ("linear_quadratic", "power_law"):
        r = fc[fam]["stage1_analog"]
        cs = ", ".join(f"{k}={float(x):.5g}"
                       for k, x in r["fit_18A_only"]["coefficients"].items())
        print(f"  {fam:28s} rmse_9A={float(r['rmse_Aps']):.4f}"
              f"  pass={r['checks']['S1_ANALOG_RMSE_MAX']['pass']}"
              f"  (18A-only fit: {cs})")

    print("mode: PER-CASE single-curve (diagnostic only, never preset-wired)")
    for case in CASES:
        tm = _load(REFERENCE_DRAG_ROOT / case / "trajectory_matching"
                   / "fit_parameters.json")
        print(f"  {case} linear_cubic  a={float(tm['a']):.5g} "
              f"b={float(tm['b']):.5g} "
              f"E_bind={float(tm['effective_binding_energy_I_ion_eV']):.5g}")
    for r in pcs:
        cs = ", ".join(f"{k}={float(x):.5g}"
                       for k, x in r["coefficients"].items())
        print(f"  {r['case']} {r['variant']:26s} {cs} "
              f"E_bind={float(r['e_bind_eV']):.5g} "
              f"obj={float(r['objective_rmse_Aps']):.5f}")


# -----------------------------------------------------------------------------
# Section 3 -- low-v trace-tail inspection
# -----------------------------------------------------------------------------
def tail_inspection(case: str, m_eff_amu: float = 202.953908) -> None:
    ref_path = (REFERENCE_DRAG_ROOT / case / "velocity_smoothed"
                / "cleaned_data_long.csv")
    ref = load_smoothed_speed_reference(ref_path)
    t = np.asarray(ref.time_ps, dtype=float)
    vv = np.asarray(ref.speed_Aps, dtype=float)
    dvdt = np.gradient(vv, t)

    print(f"case {case}: window [{t[0]:.2f}, {t[-1]:.2f}] ps, "
          f"{t.size} samples")
    print(f"  v(start)={vv[0]:.3f}  v(end)={vv[-1]:.3f}  "
          f"min v={vv.min():.3f} A/ps at t={t[np.argmin(vv)]:.2f} ps")

    for tag, mask in (
        (f"last {TAIL_LAST_PS:.0f} ps", t >= t[-1] - TAIL_LAST_PS),
        (f"v < {TAIL_V_MAX_APS:.0f} A/ps", vv < TAIL_V_MAX_APS),
    ):
        sel = mask & (dvdt < 0)
        n_dec = int(np.sum(sel))
        print(f"  tail [{tag}]: {int(np.sum(mask))} samples, "
              f"{n_dec} decelerating")
        if n_dec >= 5:
            ln_v = np.log(vv[sel])
            ln_f = np.log(-dvdt[sel])
            slope, intercept = np.polyfit(ln_v, ln_f, 1)
            print(f"    local exponent n_loc = dln(-dv/dt)/dln v "
                  f"= {slope:.2f}")
            v_med = float(np.median(vv[sel]))
            g_emp = m_eff_amu * float(np.median(-dvdt[sel])) / v_med
            print(f"    gamma_emp(median v={v_med:.2f}) ~ {g_emp:.1f} amu/ps "
                  f"(force-balance read; mixes drag + exit-well climb)")
        else:
            print("    too few decelerating samples for an exponent read")
    # coast check: does the trace flatten at finite v?
    last = t >= t[-1] - 1.0
    print(f"  coast check (last 1 ps): mean dv/dt = "
          f"{float(np.mean(dvdt[last])):+.4f} A/ps^2 at "
          f"mean v = {float(np.mean(vv[last])):.3f} A/ps")


def main() -> int:
    print_form_table()
    print_gamma_table()
    print("=" * 78)
    print("SECTION 3 -- low-v trace-tail inspection (smoothed references)")
    print("-" * 78)
    for case in CASES:
        tail_inspection(case)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
