"""Tier-0 thesis runs (THESIS_FINAL_GOALS G11.4 step 3): B2 / B3 / h405 vs TDDFT.

Runs the Tier-0 integrations and writes their series + scores; the thesis
figure is drawn from those CSVs by ``tier0_law_figure`` (plot-only, called at
the end) as ``data/figures/thesis_results_discussion/tier2_B2*.png``. The three panels:

(a) 18 A: ensemble-mean |v| of I2 (the Method-B scored atom) vs the smoothed
    TDDFT reference (raw TDDFT |v1|, |v2| in the background), scoring window
    shaded, in-window RMSE in the legend.
(b) 9 A (validity check, radial MD): ensemble-mean |v| vs HeDFT atoms 1 / 2,
    exit window shaded, exit margin vs atom 1 in the legend (G4.12 step 2).
(c) The ungated bulk laws F(v), with the speed range the 18 A trace probes.

Laws (production coefficients, not refit):

* B2 -- thesis production law: linear + quadratic, c 1.6 (committed
  ``calibration_ridge/ridge/check_B2.json``), 3.5 A gate, run at the G11.2 extraction mass
  m = I+ + 21 He.
* B3 -- effective high-tail variant: B2 with a 19.8, E_bind 0.0482 eV
  (``LQ_eff`` bundle), same setup; not
  Tier-0-consistent by construction (G10.5), shown as the departure.
* h405 -- the contrast law re-landed at the production settings
  (``gen_tier2_h405_pp35.py``): capped cubic (``CC`` bundle b with its stamped
  cap, p_tail -1, v_c 5.40), 3.5 A gate (``allow_inconsistent_gate_pairing``:
  the bundle was extracted at 14.2 A), run at m = 21 He like B2 / B3. Its
  Tier-2-only knobs (tau 4.4, E0 0.424) do not enter fixed-mass Tier 0.

Anchors (fail loudly): the B2 18 A RMSE equals the committed
``calibration_ridge/ridge/check_B2.json`` bit-exactly and its 9 A margin the committed
``calibration_ridge/exit_bound_9A/B2.json``. B3 and the re-landed h405 have no committed
Tier-0 value (new numbers, reported as measured).

Writes the series ``data/figures/thesis_microscopic_calibration/tier0_B2_B3_h405_series.csv`` +
``tier0_B2_B3_h405_scores.csv`` (the inputs of ``tier0_law_figure``, which
also draws ``tier0_ridge_m21_laws.png`` from the committed ``ridge_m21.csv``).

Usage::

    python scripts/post_processing/tier0_thesis_figure.py
"""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# =============================================================================
# USER SETTINGS
# =============================================================================
OUT_DIR = PROJECT_ROOT / "data" / "figures" / "thesis_microscopic_calibration"
GATE_A = 3.5                # density-gate width of every law here [A]

from i2_helium_md import load_capped_cc, load_drag_bundle, tddft_reference_path  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.extraction.trajectory_matching import (  # noqa: E402
    build_case_setup,
    evaluate_form_objective,
)
from i2_helium_md.physics.drag import CAPPED_CUBIC, LINEAR_QUADRATIC  # noqa: E402
from i2_helium_md.postprocess import load_hedft_trajectory  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.extraction.method_b_g412_ridge import (  # noqa: E402
    BOUND_DT_ION_PS, BOUND_ION_TIME_PS, BOUND_N, BOUND_VARIANT, DT_ION_PS, ION_TIME_PS, N, SEED,
    G11_BOUND_DIR, G11_M_EFF_AMU, G11_RIDGE_DIR, GATE_OVERRIDES,
    HEDFT_9A_PATH, REFERENCE_FILENAME, RUNS_ROOT, _run_point, check_point, load_json,
    _window_9a, score_run,
)
from scripts.post_processing.tier0_drag_comparison import ensemble_mean_series  # noqa: E402
from scripts.tier0_common import M_EFF_AMU, build_tier0_cfg, run_dir_name  # noqa: E402


HEDFT_18A_PATH = tddft_reference_path("18A")
_B2 = check_point("B2")
_B3 = load_drag_bundle("LQ_eff", expected_m_eff_amu=M_EFF_AMU)
#: lq laws at m = 21 / 3.5 A gate: label -> (a [amu/ps], c [amu/A], E_bind [eV], anchored).
LQ_LAWS = {
    "B2": (_B2["a_amu_per_ps"], _B2["c_amu_per_A"], _B2["e_bind_eV"], True),
    "B3": (_B3.coefficients["a"], _B3.coefficients["c"],
           _B3.effective_binding_energy_I_ion_eV, False),
}
STEM = "tier0_B2_B3_h405"


def _scored_run(form, coefficients, e_bind, setup):
    """Engine objective + the ion checkpoint it scored (one integration)."""
    captured = []

    def run_fn(cfg, neutral):
        captured.append((run_ion_propagation(cfg, neutral), cfg.drag_coefficients))
        return captured[-1][0]

    r = evaluate_form_objective(form, coefficients, e_bind, setup, run_fn=run_fn)
    return r, *captured[0]


def run_18a() -> dict:
    """18 A ensemble means + RMSE per law (all at m = 21, 3.5 A gate)."""
    setup21 = build_case_setup(
        "18A", num_molecules=N, seed=SEED, ion_time_ps=ION_TIME_PS, dt_ion_ps=DT_ION_PS,
        reference_filename=REFERENCE_FILENAME,
        cfg_overrides={**GATE_OVERRIDES, "m_eff_amu": G11_M_EFF_AMU,
                       "mass_initial_amu": G11_M_EFF_AMU},
    )
    out = {}
    for label, (a, c, e_bind, anchored) in LQ_LAWS.items():
        r, ion, coeffs = _scored_run(LINEAR_QUADRATIC, {"a": a, "c": c}, e_bind, setup21)
        if anchored:
            anchor = load_json(G11_RIDGE_DIR / f"check_{label}.json")["rmse_18A_Aps"]
            if r.rmse_Aps != anchor:
                raise AssertionError(f"{label}: 18 A rmse {r.rmse_Aps!r} != committed {anchor!r}")
        out[label] = dict(ion=ion, coeffs=coeffs, rmse=r.rmse_Aps, escape=r.escape_fraction,
                          m_amu=G11_M_EFF_AMU, gate_A=GATE_A, e_bind=e_bind,
                          law=f"a {a:.2f}, c {c}, E {e_bind:.3f}")
        print(f"  18A {label}: rmse {r.rmse_Aps:.4f}"
              f"{' (anchor bit-exact)' if anchored else ' (unanchored)'}, escape {r.escape_fraction}")

    # h405 re-landed: the CC bundle's b / E_bind with its stamped cap, scored in
    # the same m = 21 / 3.5 A setup (no committed anchor at 3.5 A).
    h = load_capped_cc(expected_m_eff_amu=M_EFF_AMU)
    e_bind = h.effective_binding_energy_I_ion_eV
    r, ion, coeffs = _scored_run(CAPPED_CUBIC, h.coefficients, e_bind, setup21)
    out["h405"] = dict(ion=ion, coeffs=coeffs, rmse=r.rmse_Aps, escape=r.escape_fraction,
                       m_amu=G11_M_EFF_AMU, gate_A=GATE_A, e_bind=e_bind,
                       law=f"b {h.coefficients['b']:.3f}, v_c {h.coefficients['v_c']:.2f}")
    print(f"  18A h405: rmse {r.rmse_Aps:.4f}, escape {r.escape_fraction}, "
          f"b {h.coefficients['b']:.4f}, E_bind {e_bind:.4f}")
    return out, setup21


def h405_9a_cfg() -> SimConfig:
    """Radial 9 A Tier-0 config of the re-landed h405 (G4.12 bound settings, m = 21, 3.5 A)."""
    return build_tier0_cfg("9A", "CC", num_molecules=BOUND_N,
                           ion_time_ps=BOUND_ION_TIME_PS, dt_ion_ps=BOUND_DT_ION_PS,
                           seed=SEED, m_eff_amu=G11_M_EFF_AMU)


def _h405_9a_run() -> Path:
    """The finished :func:`h405_9a_cfg` run, generated if missing."""
    cfg = h405_9a_cfg()
    run_dir = RUNS_ROOT / run_dir_name("9A", "CC", BOUND_N, "g11fig_h405pp35m21")
    if not (run_dir / "ion.npz").exists():
        run = RunDirectory(run_dir)
        run.save_cfg(cfg)
        neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
        run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    return run_dir


def run_9a() -> tuple[dict, object]:
    hedft, window = load_hedft_trajectory(HEDFT_9A_PATH), _window_9a()
    dirs = {}
    for label, (a, c, e_bind, _) in LQ_LAWS.items():
        d = RUNS_ROOT / run_dir_name("9A", BOUND_VARIANT, BOUND_N, f"g11m21{label}")
        if not (d / "ion.npz").exists():
            d = _run_point(c, a, e_bind, mass_amu=G11_M_EFF_AMU, tag=f"g11m21{label}")
        dirs[label] = d
    dirs["h405"] = _h405_9a_run()
    out = {}
    for label, d in dirs.items():
        row = score_run(d, hedft, window)
        if label in LQ_LAWS and LQ_LAWS[label][3]:
            anchor = load_json(G11_BOUND_DIR / f"{label}.json")["margin_Aps"]
            if row["margin_Aps"] != anchor:
                raise AssertionError(f"{label}: 9 A margin {row['margin_Aps']!r} != {anchor!r}")
        t, _, v1, v2 = ensemble_mean_series(RunDirectory(d).load_ion())
        out[label] = dict(t=t, v=0.5 * (v1 + v2), row=row)
        print(f"  9A {label}: margin {row['margin_Aps']:+.3f} ({'pass' if row['pass'] else 'FAIL'}), "
              f"m {row['m_eff_amu']:.3f}")
    return out, hedft


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("===== Tier-0 thesis runs: B2 / B3 / h405 =====")
    runs18, _ = run_18a()
    runs9, _ = run_9a()

    series_rows = []
    for label, d in runs18.items():
        t, _, _, v2 = ensemble_mean_series(d["ion"])
        series_rows += [("18A", label, ti, vi) for ti, vi in zip(t, v2)]
    for label, d in runs9.items():
        series_rows += [("9A", label, ti, vi) for ti, vi in zip(d["t"], d["v"])]

    with open(OUT_DIR / f"{STEM}_series.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(("case", "law", "time_ps", "mean_speed_Aps"))
        w.writerows(series_rows)
    with open(OUT_DIR / f"{STEM}_scores.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(("law", "drag", "e_bind_eV", "m_amu", "gate_A", "rmse_18A_Aps",
                    "escape_18A", "v_exit_9A_Aps", "v_exit_atom1_Aps", "margin_9A_Aps", "pass_9A"))
        for label, d in runs18.items():
            row = runs9[label]["row"]
            w.writerow((label, d["law"], d["e_bind"], d["m_amu"], d["gate_A"], d["rmse"],
                        d["escape"], row["v_exit_md_Aps"], row["v_exit_atom1_Aps"],
                        row["margin_Aps"], row["pass"]))
    # the thesis figures (law + ridge), from the CSVs just written (lazy:
    # tier0_law_figure imports this module)
    from scripts.post_processing import tier0_law_figure
    return tier0_law_figure.main()


if __name__ == "__main__":
    raise SystemExit(main())
