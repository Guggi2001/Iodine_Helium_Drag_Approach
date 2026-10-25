"""G6 chord-map figure (THESIS_FINAL_GOALS G6): chord -> (n, KE) at B2.

Three panels, per-ion at the thesis production law B2 (``post_pickup``, lq
ridge c 1.6, E0 0.3297 eV; pooled N = 6000, s1-s5 + r6), written twice:
``g6_chord_map_b2_only`` (B2 alone) and ``g6_chord_map_with_h405`` (the
re-landed contrast law h405 overlaid in every panel: ``post_pickup``, 3.5 A
gate, v_c 5.40, tau 4.4, E0 0.424 -- ``gen_tier2_h405_pp35.SCAN_POOL_CELL``,
THESIS_DISCUSSION §11; per ion from its s1-s5 + r6 member directories, binned
on the B2 chord bins):

(a) scored n vs emission-axis chord [A] (per-ion scatter + binned mean n),
    in the scorer's convention: ``suppressed`` ions (ejected net-self-unbound
    complexes, raw ``n_detected`` = 21) count as n = 0, the bare-candidate
    class (``tier2_confirmation``: bin 0; excluded from the solvated read,
    RQ8). The committed B2 chord CSV carries the raw n = 21; its mean n is
    converted where the bin holds only n = 1 and n = 21 ions (checked
    exactly; chord < 35 A) and must equal the per-ion binned curve;
(b) detected KE [eV] vs chord, n = 1 ions highlighted, suppressed ions
    marked (KE of the intact complex -- their fragmentation is open, RQ3), the B2 n = 1 window
    (chord q05-q95 of the n = 1 ions) shaded, experimental KE1 mean / median
    as horizontal lines;
(c) chord distribution: geometry reweights the chord axis only; the law
    does not enter (B2 and the re-landed h405 share seeds and sampler, so
    their chords coincide ion by ion). (The raw-sampler A5 overlay and the
    B2 = A5 ``post_pickup`` chord identity check, done 2026-10-13, were
    dropped 2026-10-24: the A5 run directories are no longer on disk.)

Chord = straight path from the birth point along the explosion axis to the
droplet surface (``emission_chords``; ignores deflection and the soft erf
edge). "Scored" = the detection stage's ``detected_mask`` (same selection as
``chord_rows``).

Pure scorer + plot -- runs no MD. Anchors (fail loudly): ``chord_rows`` over
the B2 member dirs reproduces the committed ``lqppb2`` rows of
``g5_postpickup_chord.csv`` to 4 dp; the pooled re-landed h405 read (scored
count, n̄, suppressed share) reproduces its committed ``h405_pp35_table.csv``
row to 4 dp.

Also writes the G7 numbers (n = 1 window, KE1 spread, KE slope along the
chord) to ``data/runs/h2b_forward_model/g6_chord_stats.csv`` and the KE
width at fixed n = 1-6 (model vs experiment, IQR / median; the B2-3 KED
panels) to ``g7_ke_width.csv``. Experimental quantiles come from the
unsmoothed 3-D P(E) of the corrected (v2, G12) ``IHe_KED_curves_n{n}.csv``
(dE-weighted). They do not reproduce the reference median column exactly
(curve medians 4-11 % high, n = 1: 1.094 vs 1.053 eV -- the reference
applies cuts not re-done here), so only the relative width is taken from
them; central values stay the reference columns.

Writes ``data/figures/thesis_tier2/g6_chord_map_b2_only.{png,pdf}`` and
``g6_chord_map_with_h405.{png,pdf}`` (not committed).

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
from scripts.gen_tier2_detfix_postpickup import B2_C, B2_E0_EV, B2_TAU_PS, lq_run_dir  # noqa: E402
from scripts.gen_tier2_h405_pp35 import SCAN_POOL_CELL, run_dir_for  # noqa: E402
from scripts.post_processing.tier2_detfix_postpickup_table import (  # noqa: E402
    chord_rows,
    emission_chords,
)
from scripts.post_processing.tier2_g5_postpickup_table import CHORD_CSV as G5_CHORD_CSV  # noqa: E402
from scripts.post_processing.tier2_h405_pp35_table import OUT_CSV as H405_TABLE_CSV  # noqa: E402
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    IHE_KED_REFERENCE_CSV,
    RUNS_ROOT,
)

STATS_CSV = RUNS_ROOT / "h2b_forward_model" / "g6_chord_stats.csv"
WIDTH_CSV = RUNS_ROOT / "h2b_forward_model" / "g7_ke_width.csv"
KED_CURVE_TEMPLATE = IHE_KED_REFERENCE_CSV.parent / "IHe_KED_curves_n{n}.csv"
WIDTH_N = (1, 2, 3, 4, 5, 6)
SUPPRESSED = "suppressed"
N_START = 21                        # intact shell size of a suppressed ion
SUPP_MAX_CHORD_A = 35.0             # bins entirely below this may hold suppressed ions
ANCHOR_COLS = ("ions", "retained_frac", "n1_frac", "mean_n", "mean_KE_eV")
H405_ANCHOR_COLS = ("num_scored", "nbar_det", "supp")
STYLE = {"B2": dict(color="#1f77b4", ls="-"),
         "h405": dict(color="#d62728", ls="-.")}
NAMES = {"B2": "B2", "h405": "h405 re-landed"}
VARIANTS = (("g6_chord_map_b2_only", False), ("g6_chord_map_with_h405", True))


def _b2_dirs() -> list[Path]:
    return [lq_run_dir(m, B2_E0_EV, B2_TAU_PS, B2_C) for m in MEMBER_SEEDS]


def _h405_dirs() -> list[Path]:
    """Re-landed h405 pool members (``SCAN_POOL_CELL`` = (v_c, tau, E0))."""
    v_c, tau, e0 = SCAN_POOL_CELL
    return [run_dir_for(m, e0, v_c, tau) for m in MEMBER_SEEDS]


def _check_h405_anchor(d: dict[str, np.ndarray]) -> None:
    """Pooled re-landed h405 per-ion read vs its committed
    ``h405_pp35_table.csv`` row, 4 dp (fail loudly)."""
    v_c, tau, _ = SCAN_POOL_CELL
    label = f"h405pp35_v{v_c}_t{tau}_pooled{1000 * len(MEMBER_SEEDS)}"
    with H405_TABLE_CSV.open(newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["label"] == label]
    if len(rows) != 1:
        raise RuntimeError(f"anchor: expected one row {label!r} in {H405_TABLE_CSV.name}")
    s = d["scored"]
    got = {"num_scored": float(s.sum()), "nbar_det": float(d["n_scorer"][s].mean()),
           "supp": float((d["supp"] & s).sum() / s.sum())}
    for k in H405_ANCHOR_COLS:
        if round(got[k], 4) != round(float(rows[0][k]), 4):
            raise RuntimeError(f"anchor: h405 {k} {got[k]:.4f} vs committed {rows[0][k]}")


def _binned_from_ions(d: dict[str, np.ndarray], edges: list[float]
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bin midpoints [A], mean scorer-n and mean detected KE [eV] over scored
    ions per chord bin (empty bins dropped; open last bin plotted at its lower
    edge + 10 A, as in ``_binned``)."""
    s = d["scored"] & np.isfinite(d["ke"])
    x, mn, mk = [], [], []
    for lo_e, hi_e in zip(edges[:-1], edges[1:]):
        sel = s & (d["chord"] >= lo_e) & (d["chord"] < hi_e)
        if not sel.any():
            continue
        x.append(0.5 * (lo_e + hi_e) if hi_e - lo_e < 50 else lo_e + 10.0)
        mn.append(float(d["n_scorer"][sel].mean()))
        mk.append(float(d["ke"][sel].mean()))
    return np.asarray(x), np.asarray(mn), np.asarray(mk)


def _per_ion(run_dirs: list[Path]) -> dict[str, np.ndarray]:
    """Pooled per-ion chord [A], detected n, scored mask, detected KE [eV],
    suppressed mask, and the scorer's n (suppressed -> 0)."""
    c, n, ok, ke, sup = [], [], [], [], []
    for rd in run_dirs:
        if not (rd / "detection.npz").exists():
            raise FileNotFoundError(f"missing detection.npz in {rd}")
        det = load_detection_result(rd / "detection.npz")
        c.append(emission_chords(rd))
        n.append(det.n_detected)
        ok.append(det.detected_mask)
        ke.append(det.E_kin_detected_eV)
        sup.append(np.asarray(det.state_reason).astype(str) == SUPPRESSED)
    d = {"chord": np.concatenate(c), "n": np.concatenate(n),
         "scored": np.concatenate(ok).astype(bool), "ke": np.concatenate(ke),
         "supp": np.concatenate(sup)}
    d["n_scorer"] = np.where(d["supp"], 0.0, d["n"])
    return d


def _binned_n_scorer(rows: list[dict[str, str]]) -> tuple[np.ndarray, np.ndarray]:
    """Binned mean n in the scorer's convention from a committed chord CSV.

    Bins below ``SUPP_MAX_CHORD_A`` must hold only n = 1 and n = N_START
    ions (mean_n = f1 + N_START (1 - f1), checked to 1e-9); their n = 21
    ions are the suppressed class -> mean = f1. Other bins pass unchanged.
    """
    x, y = _binned(rows, "mean_n")
    _, f1 = _binned(rows, "n1_frac")
    hi = np.array([float(r["chord_hi_A"]) for r in rows if r["mean_n"] != ""])
    out = y.copy()
    low = hi <= SUPP_MAX_CHORD_A
    if not np.allclose(y[low], f1[low] + N_START * (1 - f1[low]), atol=1e-9, rtol=0):
        raise RuntimeError("short-chord bin holds n other than 1 / 21: cannot convert")
    out[low] = f1[low]
    return x, out


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
    s_all = b2["scored"]
    if not (np.array_equal(b2["supp"] & s_all, s_all & (b2["n"] == N_START))
            and b2["chord"][b2["supp"]].max() < SUPP_MAX_CHORD_A):
        raise RuntimeError("B2: scored n = 21 is not exactly the suppressed class "
                           "below the conversion edge")
    b2_rows = _read_chord_csv(G5_CHORD_CSV, "lqppb2")
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    edges_csv = [float(r["chord_lo_A"]) for r in b2_rows] + [float(b2_rows[-1]["chord_hi_A"])]
    binned = {"B2": _binned_from_ions(b2, edges_csv)}
    xc, nc = _binned_n_scorer(b2_rows)          # per-ion B2 curves vs committed CSV
    _, kc = _binned(b2_rows, "mean_KE_eV")
    xb, nb, kb = binned["B2"]
    if not (np.array_equal(xb, xc) and np.array_equal(nb.round(4), nc.round(4))
            and np.array_equal(kb.round(4), kc.round(4))):
        raise RuntimeError("per-ion B2 binned curves differ from the committed chord CSV")
    h405 = _per_ion(_h405_dirs())
    _check_h405_anchor(h405)
    binned["h405"] = _binned_from_ions(h405, edges_csv)

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

    ke1 = _exp_ke1(ked)[:2]
    win = (stats["n1_chord_q05_A"], stats["n1_chord_q95_A"])
    for stem, with_h405 in VARIANTS:
        fig = _draw(b2, h405 if with_h405 else None, binned, win, ke1)
        png, pdf = OUT_DIR / f"{stem}.png", OUT_DIR / f"{stem}.pdf"
        fig.savefig(png, dpi=200)
        fig.savefig(pdf)
        print(f"figure -> {png}")
    if SHOW_FIGURE:
        plt.show()
    return 0


def _draw(b2: dict[str, np.ndarray], h405: dict[str, np.ndarray] | None,
          binned: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
          win: tuple[float, float], ke1: tuple[float, float]) -> plt.Figure:
    """The three-panel chord map; ``h405 = None`` leaves the overlay out."""
    arms = ("B2",) if h405 is None else ("B2", "h405")
    s = b2["scored"] & np.isfinite(b2["ke"])
    idx = np.flatnonzero(s)
    rng = np.random.default_rng(SCATTER_SEED)
    if idx.size > SCATTER_MAX:
        idx = np.sort(rng.choice(idx, SCATTER_MAX, replace=False))
    is1 = b2["n"][idx] == 1
    isS = b2["supp"][idx]
    lo, hi = CHORD_RANGE_A

    fig, (axn, axk, axc) = plt.subplots(1, 3, figsize=(15.0, 4.6), constrained_layout=True)

    for ax in (axn, axk):
        ax.axvspan(*win, color="0.92", zorder=0, label="B2 n = 1 window (q05–q95)")
    jit = rng.uniform(-0.25, 0.25, idx.size)
    axn.scatter(b2["chord"][idx][~isS], b2["n_scorer"][idx][~isS] + jit[~isS],
                s=3, alpha=0.25, color=STYLE["B2"]["color"], lw=0, rasterized=True)
    axn.scatter(b2["chord"][idx][isS], b2["n_scorer"][idx][isS] + jit[isS], s=3,
                alpha=0.35, color="0.35", lw=0, rasterized=True,
                label="B2 suppressed → n = 0 (bare candidate)")
    for arm in arms:
        x, y, _ = binned[arm]
        axn.plot(x, y, marker="o", ms=4, lw=1.8, **STYLE[arm],
                 label=f"{NAMES[arm]} binned mean n")
    axn.set(xlim=CHORD_RANGE_A, ylim=(-0.7, 18), xlabel="emission-axis chord [Å]",
            ylabel="scored n", title="(a) chord → n")
    axn.legend(fontsize=8, loc="upper left")

    axk.scatter(b2["chord"][idx][isS], b2["ke"][idx][isS], s=3, alpha=0.35,
                color="0.35", lw=0, rasterized=True,
                label="B2 suppressed (intact-complex KE)")
    rest = ~is1 & ~isS
    axk.scatter(b2["chord"][idx][rest], b2["ke"][idx][rest], s=3, alpha=0.2,
                color=STYLE["B2"]["color"], lw=0, rasterized=True, label="B2 ions, n ≥ 2")
    axk.scatter(b2["chord"][idx][is1], b2["ke"][idx][is1], s=4, alpha=0.6,
                color="#ff7f0e", lw=0, rasterized=True, label="B2 ions, n = 1")
    for arm in arms:
        x, _, y = binned[arm]
        axk.plot(x, y, marker="o", ms=4, lw=1.8, **STYLE[arm],
                 label=f"{NAMES[arm]} binned mean KE")
    ke1_mean, ke1_med = ke1
    axk.axhline(ke1_mean, color="k", lw=1.2, label=f"exp KE₁ mean {ke1_mean:.2f} eV")
    axk.axhline(ke1_med, color="k", lw=1.2, ls="--", label=f"exp KE₁ median {ke1_med:.2f} eV")
    axk.set(xlim=CHORD_RANGE_A, ylim=(0, 1.6), xlabel="emission-axis chord [Å]",
            ylabel="detected KE [eV]", title="(b) chord → KE")
    axk.legend(fontsize=7.5, loc="upper right")

    edges = np.arange(lo, hi + HIST_BIN_A, HIST_BIN_A)
    axc.hist(b2["chord"], bins=edges, density=True, histtype="step", lw=1.6,
             **STYLE["B2"], label=f"B2 (post_pickup), N = {b2['chord'].size}")
    if h405 is not None:
        axc.hist(h405["chord"], bins=edges, density=True, histtype="step", lw=1.6,
                 **STYLE["h405"], label=f"{NAMES['h405']} (post_pickup), N = {h405['chord'].size}")
    axc.axvspan(*win, color="0.92", zorder=0)
    axc.set_ylim(0, 1.25 * axc.get_ylim()[1])
    axc.set(xlim=CHORD_RANGE_A, xlabel="emission-axis chord [Å]",
            ylabel="probability density [1/Å]", title="(c) chord distribution")
    axc.legend(fontsize=8)
    return fig


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
