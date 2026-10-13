"""G6 chord-map figure (THESIS_FINAL_GOALS G6): chord -> (n, KE) at B2.

Three panels, per-ion at the thesis production law B2 (``post_pickup``, lq
ridge c 1.6, E0 0.3297 eV; pooled N = 6000, s1-s5 + r6), with the old-bundle
map (h405, capped cubic v_c 5.5, 14.2 A gate, ``post_pickup``; binned only --
its run directories are no longer on disk) overlaid from the committed
``detfix_postpickup_chord.csv``:

(a) detected n vs emission-axis chord [A] (per-ion scatter + binned mean n);
(b) detected KE [eV] vs chord, n = 1 ions highlighted, the B2 n = 1 window
    (chord q05-q95 of the n = 1 ions) shaded, experimental KE1 mean / median
    as horizontal lines;
(c) chord distribution per sampler (raw: A5; ``post_pickup``: B2, whose
    births and chords are bit-identical to A5 ``post_pickup`` -- checked):
    geometry reweights the chord axis only; the law does not enter.

Chord = straight path from the birth point along the explosion axis to the
droplet surface (``emission_chords``; ignores deflection and the soft erf
edge). "Scored" = the detection stage's ``detected_mask`` (same selection as
``chord_rows``).

Pure scorer + plot -- runs no MD. Anchor (fails loudly): ``chord_rows`` over
the B2 member dirs reproduces the committed ``lqppb2`` rows of
``g5_postpickup_chord.csv`` to 4 dp.

Also writes the G7 numbers (n = 1 window, KE1 spread, KE slope along the
chord) to ``data/runs/h2b_forward_model/g6_chord_stats.csv`` and the KE
width at fixed n = 1-4 (model vs experiment, IQR / median) to
``g7_ke_width.csv``. Experimental quantiles come from the unsmoothed 3-D
P(E) of ``IHe_KED_curves_n{n}.csv`` (dE-weighted; the n = 1 curve equals
the spectra file). They do not reproduce the reference median / mean
columns exactly (n = 1: median 1.168 vs 1.128, mean 1.372 vs 1.302 eV --
the reference applies cuts not re-done here), so only the relative width
is taken from them; central values stay the reference columns.

Writes ``data/figures/thesis_tier2/g6_chord_map.{png,pdf}`` (not committed).

Usage::

    python scripts/post_processing/tier2_g6_chord_map.py
"""

from __future__ import annotations

import csv
from pathlib import Path
import sys
import warnings

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# =============================================================================
# USER SETTINGS
# =============================================================================
OUT_DIR = PROJECT_ROOT / "data" / "figures" / "thesis_tier2"
SHOW_FIGURE = False
CHORD_RANGE_A = (20.0, 90.0)        # plotted chord axis [A]
HIST_BIN_A = 1.0                    # panel (c) chord bin width [A]
SCATTER_MAX = 6000                  # per-ion points drawn in (a)/(b) (random subset)
SCATTER_SEED = 0

import matplotlib  # noqa: E402

if not SHOW_FIGURE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from scripts.gen_tier2_detfix_battery import MEMBER_SEEDS  # noqa: E402
from scripts.gen_tier2_detfix_postpickup import (  # noqa: E402
    B2_C, B2_E0_EV, B2_TAU_PS, a5_member_cell, lq_run_dir,
)
from scripts.gen_tier2_g414_tauline import run_dir_for as a5_run_dir_for  # noqa: E402
from scripts.post_processing.tier2_detfix_postpickup_table import (  # noqa: E402
    CHORD_CSV as H405_CHORD_CSV,
    chord_rows,
    emission_chords,
)
from scripts.post_processing.tier2_g5_postpickup_table import CHORD_CSV as G5_CHORD_CSV  # noqa: E402
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    IHE_KED_REFERENCE_CSV,
    RUNS_ROOT,
)

STATS_CSV = RUNS_ROOT / "h2b_forward_model" / "g6_chord_stats.csv"
WIDTH_CSV = RUNS_ROOT / "h2b_forward_model" / "g7_ke_width.csv"
KED_CURVE_TEMPLATE = IHE_KED_REFERENCE_CSV.parent / "IHe_KED_curves_n{n}.csv"
WIDTH_N = (1, 2, 3, 4)
ANCHOR_COLS = ("ions", "retained_frac", "n1_frac", "mean_n", "mean_KE_eV")
STYLE = {"B2": dict(color="#1f77b4", ls="-"),
         "h405": dict(color="#d62728", ls="-."),
         "A5 raw": dict(color="0.45", ls="--")}


def _b2_dirs() -> list[Path]:
    return [lq_run_dir(m, B2_E0_EV, B2_TAU_PS, B2_C) for m in MEMBER_SEEDS]


def _per_ion(run_dirs: list[Path]) -> dict[str, np.ndarray]:
    """Pooled per-ion chord [A], detected n, scored mask, detected KE [eV]."""
    c, n, ok, ke = [], [], [], []
    for rd in run_dirs:
        if not (rd / "detection.npz").exists():
            raise FileNotFoundError(f"missing detection.npz in {rd}")
        det = load_detection_result(rd / "detection.npz")
        c.append(emission_chords(rd))
        n.append(det.n_detected)
        ok.append(det.detected_mask)
        ke.append(det.E_kin_detected_eV)
    return {"chord": np.concatenate(c), "n": np.concatenate(n),
            "scored": np.concatenate(ok).astype(bool), "ke": np.concatenate(ke)}


def _read_chord_csv(path: Path, arm: str) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["arm"] == arm]
    if not rows:
        raise ValueError(f"arm {arm!r} not in {path}")
    return rows


def _check_anchor(b2_dirs: list[Path]) -> None:
    """B2 chord_rows vs the committed lqppb2 rows, 4 dp (fail loudly)."""
    fresh = chord_rows("lqppb2", b2_dirs)
    committed = _read_chord_csv(G5_CHORD_CSV, "lqppb2")
    if len(fresh) != len(committed):
        raise RuntimeError("anchor: bin count differs")
    for f, c in zip(fresh, committed):
        for k in ANCHOR_COLS:
            a, b = f[k], c[k]
            if (a == "") != (b == ""):
                raise RuntimeError(f"anchor: bin {c['chord_lo_A']} {k} {a!r} vs {b!r}")
            if a != "" and round(float(a), 4) != round(float(b), 4):
                raise RuntimeError(f"anchor: bin {c['chord_lo_A']} {k} {a} vs {b}")


def _binned(rows: list[dict[str, str]], key: str) -> tuple[np.ndarray, np.ndarray]:
    """Bin midpoints [A] and values of one column (empty cells dropped; the
    open-ended last bin plotted at its lower edge + 10 A)."""
    x, y = [], []
    for r in rows:
        if r[key] == "":
            continue
        lo, hi = float(r["chord_lo_A"]), float(r["chord_hi_A"])
        x.append(0.5 * (lo + hi) if hi - lo < 50 else lo + 10.0)
        y.append(float(r[key]))
    return np.asarray(x), np.asarray(y)


def _g7_stats(d: dict[str, np.ndarray]) -> dict[str, float]:
    """n = 1 window on the chord axis and the KE1 spread it implies."""
    s = d["scored"] & np.isfinite(d["ke"])
    one = s & (d["n"] == 1)
    c1, k1 = d["chord"][one], d["ke"][one]
    slope, _ = np.polyfit(c1, k1, 1)
    q05, q50, q95 = np.quantile(c1, [0.05, 0.5, 0.95])
    resid = k1 - np.polyval(np.polyfit(c1, k1, 1), c1)
    return {"N_ions": int(d["chord"].size), "N_scored": int(s.sum()), "N_n1": int(one.sum()),
            "n1_chord_q05_A": q05, "n1_chord_median_A": q50, "n1_chord_q95_A": q95,
            "n1_window_width_A": q95 - q05,
            "KE1_mean_eV": float(k1.mean()), "KE1_sd_eV": float(k1.std(ddof=1)),
            "KE1_slope_eV_per_A": float(slope),
            "KE1_sd_from_slope_eV": float(abs(slope) * c1.std(ddof=1)),
            "KE1_resid_sd_eV": float(resid.std(ddof=1))}


def _exp_quartiles(n: int) -> tuple[float, float, float]:
    """Experimental KE quartiles [eV] at fixed n from the unsmoothed 3-D P(E)
    (non-uniform E grid -> trapezoid-weighted CDF; NaN regions dropped)."""
    with KED_CURVE_TEMPLATE.with_name(f"IHe_KED_curves_n{n}.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    E = np.array([float(r["E_eV"]) for r in rows])
    P = np.array([float(r["signal_3d_PE"]) for r in rows])
    ok = np.isfinite(E) & np.isfinite(P)
    E, P = E[ok], np.clip(P[ok], 0.0, None)
    cdf = np.r_[0.0, np.cumsum(0.5 * (P[1:] + P[:-1]) * np.diff(E))]
    if cdf[-1] <= 0:
        raise ValueError(f"empty P(E) for n = {n}")
    cdf /= cdf[-1]
    q1, med, q3 = (float(np.interp(q, cdf, E)) for q in (0.25, 0.5, 0.75))
    return q1, med, q3


def _width_rows(d: dict[str, np.ndarray], ked) -> list[dict[str, float]]:
    """KE width at fixed n: model (B2 per ion) vs experiment, IQR / median."""
    s = d["scored"] & np.isfinite(d["ke"])
    rows = []
    for n in WIDTH_N:
        k = d["ke"][s & (d["n"] == n)]
        mq1, mmed, mq3 = np.quantile(k, [0.25, 0.5, 0.75])
        eq1, emed, eq3 = _exp_quartiles(n)
        i = int(np.flatnonzero(ked.n == n)[0])
        rows.append({"n": n, "model_N": int(k.size), "model_median_eV": float(mmed),
                     "model_sd_eV": float(k.std(ddof=1)),
                     "model_iqr_over_med": float((mq3 - mq1) / mmed),
                     "exp_median_ref_eV": float(ked.median_KE_eV[i]),
                     "exp_sigma_ref_eV": float(ked.sigma_KE_eV[i]),
                     "exp_iqr_over_med": (eq3 - eq1) / emed,
                     "ratio_exp_over_model": ((eq3 - eq1) / emed) / float((mq3 - mq1) / mmed)})
    return rows


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    b2_dirs = _b2_dirs()
    _check_anchor(b2_dirs)
    b2 = _per_ion(b2_dirs)
    a5_raw = _per_ion([a5_run_dir_for(a5_member_cell(m)) for m in MEMBER_SEEDS])
    a5_pp = _per_ion([lq_run_dir(m) for m in MEMBER_SEEDS])
    b2_rows = _read_chord_csv(G5_CHORD_CSV, "lqppb2")
    h405_rows = _read_chord_csv(H405_CHORD_CSV, "pp")
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    stats = _g7_stats(b2)
    stats["exp_KE1_sigma_eV"] = _exp_ke1(ked)[2]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with STATS_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["arm", *stats])
        w.writeheader()
        w.writerow({"arm": "B2", **stats})
    print(f"CSV -> {STATS_CSV}")
    width = _width_rows(b2, ked)
    with WIDTH_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(width[0]))
        w.writeheader()
        w.writerows(width)
    print(f"CSV -> {WIDTH_CSV}")
    for r in width:
        print(f"  n={r['n']}: model IQR/med {r['model_iqr_over_med']:.3f}  "
              f"exp {r['exp_iqr_over_med']:.3f}  ratio {r['ratio_exp_over_model']:.1f}")
    for k, v in stats.items():
        print(f"  {k:24s} {v:.4f}" if isinstance(v, float) else f"  {k:24s} {v}")

    s = b2["scored"] & np.isfinite(b2["ke"])
    idx = np.flatnonzero(s)
    rng = np.random.default_rng(SCATTER_SEED)
    if idx.size > SCATTER_MAX:
        idx = np.sort(rng.choice(idx, SCATTER_MAX, replace=False))
    is1 = b2["n"][idx] == 1
    lo, hi = CHORD_RANGE_A
    win = (stats["n1_chord_q05_A"], stats["n1_chord_q95_A"])

    fig, (axn, axk, axc) = plt.subplots(1, 3, figsize=(15.0, 4.6), constrained_layout=True)

    for ax in (axn, axk):
        ax.axvspan(*win, color="0.92", zorder=0, label="B2 n = 1 window (q05–q95)")
    axn.scatter(b2["chord"][idx], b2["n"][idx] + rng.uniform(-0.25, 0.25, idx.size),
                s=3, alpha=0.25, color=STYLE["B2"]["color"], lw=0, rasterized=True)
    for label, rows in (("B2", b2_rows), ("h405", h405_rows)):
        x, y = _binned(rows, "mean_n")
        axn.plot(x, y, marker="o", ms=4, lw=1.8, **STYLE[label],
                 label=f"{label} binned mean n")
    axn.set(xlim=CHORD_RANGE_A, ylim=(-0.5, 22), xlabel="emission-axis chord [Å]",
            ylabel="detected n", title="(a) chord → n")
    axn.legend(fontsize=8, loc="upper right")

    axk.scatter(b2["chord"][idx][~is1], b2["ke"][idx][~is1], s=3, alpha=0.2,
                color=STYLE["B2"]["color"], lw=0, rasterized=True, label="B2 ions, n ≥ 2")
    axk.scatter(b2["chord"][idx][is1], b2["ke"][idx][is1], s=4, alpha=0.6,
                color="#ff7f0e", lw=0, rasterized=True, label="B2 ions, n = 1")
    for label, rows in (("B2", b2_rows), ("h405", h405_rows)):
        x, y = _binned(rows, "mean_KE_eV")
        axk.plot(x, y, marker="o", ms=4, lw=1.8, **STYLE[label],
                 label=f"{label} binned mean KE")
    ke1_mean, ke1_med, _ = _exp_ke1(ked)
    axk.axhline(ke1_mean, color="k", lw=1.2, label=f"exp KE₁ mean {ke1_mean:.2f} eV")
    axk.axhline(ke1_med, color="k", lw=1.2, ls="--", label=f"exp KE₁ median {ke1_med:.2f} eV")
    axk.set(xlim=CHORD_RANGE_A, ylim=(0, 1.6), xlabel="emission-axis chord [Å]",
            ylabel="detected KE [eV]", title="(b) chord → KE")
    axk.legend(fontsize=7.5, loc="upper right")

    if not np.array_equal(a5_pp["chord"], b2["chord"]):
        raise RuntimeError("A5 post_pickup and B2 chords differ (births should be identical)")
    edges = np.arange(lo, hi + HIST_BIN_A, HIST_BIN_A)
    for label, d, text in (("A5 raw", a5_raw, "raw (A5)"),
                           ("B2", b2, "post_pickup (B2 = A5, identical births)")):
        axc.hist(d["chord"], bins=edges, density=True, histtype="step", lw=1.6,
                 color=STYLE[label]["color"], ls=STYLE[label]["ls"],
                 label=f"{text}, N = {d['chord'].size}")
    axc.axvspan(*win, color="0.92", zorder=0)
    axc.set_ylim(0, 1.25 * axc.get_ylim()[1])
    axc.set(xlim=CHORD_RANGE_A, xlabel="emission-axis chord [Å]",
            ylabel="probability density [1/Å]", title="(c) chord distribution per sampler")
    axc.legend(fontsize=8)

    png, pdf = OUT_DIR / "g6_chord_map.png", OUT_DIR / "g6_chord_map.pdf"
    fig.savefig(png, dpi=200)
    fig.savefig(pdf)
    print(f"figure -> {png}")
    if SHOW_FIGURE:
        plt.show()
    return 0


def _exp_ke1(ked) -> tuple[float, float, float]:
    """Experimental n = 1 KE mean, median and sigma [eV] (I+He KED reference)."""
    i = np.flatnonzero(ked.n == 1)
    if i.size != 1:
        raise ValueError("I+He KED reference has no unique n = 1 row")
    k = int(i[0])
    return (float(ked.mean_KE_eV[k]), float(ked.median_KE_eV[k]),
            float(ked.sigma_KE_eV[k]))


if __name__ == "__main__":
    raise SystemExit(main())
