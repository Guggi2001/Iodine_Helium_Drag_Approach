"""Tier-0 thesis figure (THESIS_FINAL_GOALS G11.4 step 3): B2 / B3 / h405 vs TDDFT.

Three panels:

(a) 18 A: ensemble-mean |v| of I2 (the Method-B scored atom) vs the smoothed
    TDDFT reference (raw TDDFT |v1|, |v2| in the background), scoring window
    shaded, in-window RMSE in the legend.
(b) 9 A (validity check, radial MD): ensemble-mean |v| vs HeDFT atoms 1 / 2,
    exit window shaded, exit margin vs atom 1 in the legend (G4.12 step 2).
(c) The ungated bulk laws F(v), with the speed range the 18 A trace probes.

Laws (production coefficients, not refit):

* B2 -- thesis production law: linear + quadratic, c 1.6 (committed
  ``g412_ridge/c1p600.json``), 3.5 A gate, run at the G11.2 extraction mass
  m = I+ + 21 He.
* B3 -- effective high-tail variant: B2 with a 19.8, E_bind 0.0482 eV
  (``gen_tier2_detfix_postpickup.py`` B3 cell), same setup; not
  Tier-0-consistent by construction (G10.5), shown as the departure.
* h405 -- the contrast law re-landed at the production settings
  (``gen_tier2_h405_pp35.py``): capped cubic (``shared_pure_cubic`` bundle b,
  p_tail -1, cap v_c 5.40), 3.5 A gate (``allow_inconsistent_gate_pairing``:
  the bundle was extracted at 14.2 A), run at m = 21 He like B2 / B3. Its
  Tier-2-only knobs (tau 4.4, E0 0.424) do not enter fixed-mass Tier 0.

Anchors (fail loudly): the B2 18 A RMSE equals the committed
``g11_ridge_m21/check_B2.json`` bit-exactly and its 9 A margin the committed
``g11_exit_bound_m21/B2.json``. B3 and the re-landed h405 have no committed
Tier-0 value (new numbers, reported as measured).

Writes ``data/figures/thesis_tier0/tier0_B2_B3_h405.png`` and the plotted
series ``tier0_B2_B3_h405_series.csv`` + ``tier0_B2_B3_h405_scores.csv``, plus
``tier0_ridge_m21_laws.png``: panel (c) on its own for every ``g11_ridge_m21``
ridge point, the c = 1.6 point replaced by B2.

Usage::

    python scripts/post_processing/tier0_thesis_figure.py
"""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# =============================================================================
# USER SETTINGS
# =============================================================================
OUT_DIR = PROJECT_ROOT / "data" / "figures" / "thesis_tier0"
SHOW_FIGURE = False
H405_V_C_APS = 5.40         # re-landed h405 cap [A/ps] (gen_tier2_h405_pp35 SCAN_POOL_CELL)
H405_P_TAIL = -1.0          # gen_tier2atlas_g4finals.P_TAIL

import matplotlib  # noqa: E402

if not SHOW_FIGURE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from i2_helium_md import load_drag_bundle  # noqa: E402
from i2_helium_md.extraction.trajectory_matching import (  # noqa: E402
    build_case_setup,
    evaluate_form_objective,
)
from i2_helium_md.physics.drag import CAPPED_CUBIC, LINEAR_QUADRATIC, drag_force  # noqa: E402
from i2_helium_md.postprocess import load_hedft_trajectory  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.extraction.form_phase_common import DT_ION_PS, ION_TIME_PS, N, SEED, load_json  # noqa: E402
from scripts.extraction.method_b_g412_ridge import (  # noqa: E402
    BOUND_DT_ION_PS, BOUND_ION_TIME_PS, BOUND_N, BOUND_VARIANT, EXIT_WINDOW_PS,
    G11_BOUND_DIR, G11_M_EFF_AMU, G11_RIDGE_DIR, GATE_OVERRIDES,
    HEDFT_9A_PATH, REFERENCE_FILENAME, RUNS_ROOT, _committed_point, _run_point,
    _window_9a, score_run,
)
from scripts.post_processing.tier0_drag_comparison import ensemble_mean_series  # noqa: E402
from scripts.tier0_common import M_EFF_AMU, build_drag_cfg, run_dir_name  # noqa: E402


scaling_factor = 0.85
HEDFT_18A_PATH = PROJECT_ROOT / "data" / "reference" / "18A_All_Data.csv"
H405_VARIANT = "shared_pure_cubic"
_B2 = _committed_point(1.6)
#: lq laws at m = 21 / 3.5 A gate: label -> (a [amu/ps], c [amu/A], E_bind [eV], anchored).
LQ_LAWS = {
    "B2": (_B2["a_amu_per_ps"], _B2["c_amu_per_A"], _B2["e_bind_eV"], True),
    "B3": (19.8, 1.6, 0.0482, False),   # gen_tier2_detfix_postpickup B3_A / B1W_*
}
STEM = "tier0_B2_B3_h405"
RIDGE_STEM = "tier0_ridge_m21_laws"
STYLE = {"B2": dict(color="#1f77b4", ls="-"),
         "B3": dict(color="#9467bd", ls="--"),
         "h405": dict(color="#d62728", ls="-.")}


def _capped(coeffs):
    """The h405 law: the bundle's cubic b with the production cap."""
    return replace(coeffs, form=CAPPED_CUBIC,
                   coefficients={"b": float(coeffs.coefficients["b"]),
                                 "v_c": H405_V_C_APS, "p_tail": H405_P_TAIL})


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
                          m_amu=G11_M_EFF_AMU, gate_A=3.5, e_bind=e_bind,
                          law=f"a {a:.2f}, c {c}, E {e_bind:.3f}")
        print(f"  18A {label}: rmse {r.rmse_Aps:.4f}"
              f"{' (anchor bit-exact)' if anchored else ' (unanchored)'}, escape {r.escape_fraction}")

    # h405 re-landed: the preset bundle's (shared_pure_cubic) b / E_bind, capped,
    # scored in the same m = 21 / 3.5 A setup (no committed anchor at 3.5 A).
    bundle = load_drag_bundle(H405_VARIANT, "18A", expected_m_eff_amu=M_EFF_AMU)
    h = _capped(bundle)
    e_bind = bundle.effective_binding_energy_I_ion_eV
    r, ion, coeffs = _scored_run(CAPPED_CUBIC, h.coefficients, e_bind, setup21)
    out["h405"] = dict(ion=ion, coeffs=coeffs, rmse=r.rmse_Aps, escape=r.escape_fraction,
                       m_amu=G11_M_EFF_AMU, gate_A=3.5, e_bind=e_bind,
                       law=f"b {h.coefficients['b']:.3f}, v_c {H405_V_C_APS:.2f}")
    print(f"  18A h405: rmse {r.rmse_Aps:.4f}, escape {r.escape_fraction}, "
          f"b {h.coefficients['b']:.4f}, E_bind {e_bind:.4f}")
    return out, setup21


def _h405_9a_run() -> Path:
    """Radial 9 A Tier-0 run of the re-landed h405 (G4.12 bound settings, m = 21, 3.5 A)."""
    cfg = build_drag_cfg("9A", H405_VARIANT, num_molecules=BOUND_N,
                         ion_time_ps=BOUND_ION_TIME_PS, dt_ion_ps=BOUND_DT_ION_PS, seed=SEED,
                         m_eff_amu=G11_M_EFF_AMU)
    cfg = replace(cfg, drag_form=CAPPED_CUBIC, drag_coefficients=_capped(cfg.drag_coefficients),
                  **GATE_OVERRIDES, allow_inconsistent_gate_pairing=True)
    cfg.validate()
    run_dir = RUNS_ROOT / run_dir_name("9A", H405_VARIANT, BOUND_N, "g11fig_h405pp35m21")
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


def plot_ridge_laws(b2_coeffs, v_ref) -> Path:
    """Standalone panel-(c) figure: the bulk F(v) of every m = 21 ridge point.

    The ridge points are the rows of ``g11_ridge_m21/ridge_m21.csv``; the
    c = 1.6 row is replaced by B2 (its a / E_bind, the production point on the
    ridge). ``b2_coeffs`` is B2's scored DragCoefficients (form + gate stamp
    reused for every point); ``v_ref`` [A/ps] the smoothed 18 A trace speeds.
    """
    with open(G11_RIDGE_DIR / "ridge_m21.csv", newline="", encoding="utf-8") as fh:
        # Ridge-grid rows only (labels c<x>p<yyy>); the B2 / A5 check rows are skipped.
        ridge = [(float(r["c_amu_per_A"]), float(r["a_amu_per_ps"]), r["bound_pass"] == "True")
                 for r in csv.DictReader(fh) if r["label"].startswith("c")]
    if not any(c == LQ_LAWS["B2"][1] for c, _, _ in ridge):
        raise AssertionError(f"ridge_m21.csv has no c = {LQ_LAWS['B2'][1]} row to replace by B2")
    c_max = max(c for c, _, _ in ridge)
    cmap = plt.get_cmap("viridis")
    v_lo, v_hi = float(v_ref.min()), float(v_ref.max())
    a_b2, c_b2 = b2_coeffs.coefficients["a"], b2_coeffs.coefficients["c"]

    fig, (ax, axz) = plt.subplots(1, 2, figsize=(10.4*scaling_factor, 4.6*scaling_factor), constrained_layout=True)
    # (b) zoom: the trace-speed band with a 0.25 A/ps margin each side.
    for axis, v in ((ax, np.linspace(0.0, 12.0, 241)),
                    (axz, np.linspace(v_lo - 0.25, v_hi + 0.25, 201))):
        bulk = np.full_like(v, -1e3)
        axis.axvspan(v_lo, v_hi, color="0.92", zorder=0, label="18 Å trace speeds")
        for c, a, bound_pass in ridge:
            if c == LQ_LAWS["B2"][1]:
                continue
            coeffs = replace(b2_coeffs, coefficients={"a": a, "c": c})
            axis.plot(v, drag_force(v, bulk, coeffs, 3.5), lw=1.2, alpha=0.55,
                      ls="-" if bound_pass else "--", color=cmap(0.9 * c / c_max),
                      label=f"c {c:g} (a {a:.2f})" + ("" if bound_pass else ", fails 9 Å"))
        axis.plot(v, drag_force(v, bulk, b2_coeffs, 3.5), lw=2.4, zorder=5,
                  **STYLE["B2"], label=f"B2: c {c_b2:g} (a {a_b2:.2f})")
        axis.set(xlim=(v[0], v[-1]), xlabel="|v| / Å/ps", ylabel="$F_{drag}$ / amu·Å/ps²")
    ax.set(ylim=(0, None), title="(a) 18 Å ridge (m = 21 He): drag laws in bulk He")
    axz.set(title="(b) Zoom: 18 Å trace speeds")
    ax.legend(fontsize=7, loc="upper left")
    png = OUT_DIR / f"{RIDGE_STEM}.png"
    fig.savefig(png, dpi=200)
    return png


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("===== Tier-0 thesis figure: B2 / B3 / h405 =====")
    runs18, setup21 = run_18a()
    runs9, hedft9 = run_9a()
    window18 = setup21.window

    fig, (ax18, ax9, axF) = plt.subplots(1, 3, figsize=(15.0, 4.6), constrained_layout=True)

    # (a) 18 A
    ax18.axvspan(*window18, color="0.92", zorder=0)
    hedft18 = load_hedft_trajectory(HEDFT_18A_PATH)
    ax18.plot(hedft18.time_ps, hedft18.v1_magnitude_Aps, color="0.6", lw=0.9, alpha=0.7,
              label="TDDFT 18 Å |v1|, |v2| (raw)")
    ax18.plot(hedft18.time_ps, hedft18.v2_magnitude_Aps, color="0.6", lw=0.9, alpha=0.7)
    ax18.plot(setup21.t_ref_ps, setup21.ref_speed_Aps, color="k", lw=2.2,
              label="TDDFT 18 Å (smoothed)")
    series_rows = []
    for label, d in runs18.items():
        t, _, _, v2 = ensemble_mean_series(d["ion"])
        ax18.plot(t, v2, lw=1.6, **STYLE[label],
                  label=f"{label} ({d['law']}; RMSE {d['rmse']:.3f} Å/ps)")
        series_rows += [("18A", label, ti, vi) for ti, vi in zip(t, v2)]
    ax18.set(xlim=(0, 16), xlabel="t [ps]", ylabel="mean |v| [Å/ps]",
             title="(a) 18 Å: fitted trace (band RMSE ≤ 0.19)")
    ax18.legend(fontsize=7.5, loc="lower right")

    # (b) 9 A
    ax9.axvspan(*EXIT_WINDOW_PS, color="0.92", zorder=0)
    ax9.plot(hedft9.time_ps, hedft9.v1_magnitude_Aps, color="k", lw=2.0, label="TDDFT 9 Å atom 1")
    ax9.plot(hedft9.time_ps, hedft9.v2_magnitude_Aps, color="0.55", lw=2.0, label="TDDFT 9 Å atom 2")
    for label, d in runs9.items():
        ax9.plot(d["t"], d["v"], lw=1.6, **STYLE[label],
                 label=f"{label} (exit margin {d['row']['margin_Aps']:+.2f} Å/ps)")
        series_rows += [("9A", label, ti, vi) for ti, vi in zip(d["t"], d["v"])]
    ax9.set(xlim=(0, 16), xlabel="t [ps]", ylabel="mean |v| [Å/ps]",
            title="(b) 9 Å: exit bound (radial MD, held out)")
    ax9.legend(fontsize=7.5, loc="lower right")

    # (c) bulk laws
    v = np.linspace(0.0, 12.0, 241)
    v_ref = np.asarray(setup21.ref_speed_Aps)
    axF.axvspan(float(v_ref.min()), float(v_ref.max()), color="0.92", zorder=0,
                label="18 Å trace speeds")
    for label, d in runs18.items():
        F = drag_force(v, np.full_like(v, -1e3), d["coeffs"], 3.5)
        axF.plot(v, F, lw=1.8, **STYLE[label], label=label)
    axF.set(xlim=(0, 12), ylim=(0, None), xlabel="|v| [Å/ps]",
            ylabel="$F_{drag}$ / amu·Å/ps²", title="(c) drag laws in bulk He")
    axF.legend(fontsize=8, loc="upper left")

    png = OUT_DIR / f"{STEM}.png"
    fig.savefig(png, dpi=200)
    ridge_png = plot_ridge_laws(runs18["B2"]["coeffs"], v_ref)

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
    for p in (png, ridge_png):
        print(f"  wrote {p.relative_to(PROJECT_ROOT).as_posix()}")
    if SHOW_FIGURE:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
