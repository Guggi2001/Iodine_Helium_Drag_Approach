"""Tier-0 law figure in the thesis style: B2 (and B3 / CC5.4) vs TDDFT.

Plot-only twin of ``tier0_thesis_figure`` (THESIS_FINAL_GOALS G11.4 step 3):
it runs **no MD**. It reads the Tier-0 series and scores that script writes
(committed ``data/figures/thesis_tier0/tier0_B2_B3_h405_{series,scores}.csv``),
the TDDFT references and the drag-law coefficients, and draws three panels at
the thesis text width with the shared ``thesis_style``:

(a) 18 Å: ensemble-mean |v| of the scored atom vs the smoothed TDDFT trace
    (raw |v1|, |v2| behind it); scoring window shaded, in-window RMSE
    annotated.
(b) 9 Å (held out, radial MD): ensemble-mean |v| vs TDDFT atoms 1 / 2 (drawn
    like the raw traces in (a)); exit window shaded, the exit-speed bound
    (atom 1's exit-window mean) dashed, exit margin annotated (G4.12 step 2).
(c) The ungated bulk laws F(v); the 18 Å trace speeds shaded.

Laws as in ``tier0_thesis_figure`` (m = I+ + 21 He, 3.5 Å gate): B2, B3 and
the re-landed h405 under its thesis name CC5.4 (legend names and colours as
``tier2_law_comparison_figure``). Written once per ``VARIANTS`` entry to
``data/figures/thesis_tier2/``: ``tier2_B2_h405`` (B2 vs CC5.4) and
``tier2_B2_B3_h405`` (all three).

Anchors (fail loudly): the B2 18 Å RMSE in the scores CSV equals the committed
``g11_ridge_m21/check_B2.json`` and its 9 Å margin the committed
``g11_exit_bound_m21/B2.json``, bit-exactly; every law of a variant must have
a series in both cases.

Usage::

    python scripts/post_processing/tier0_law_figure.py
"""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT_DIR = Path(__file__).resolve().parent
for _p in (PROJECT_ROOT, _SCRIPT_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from scripts.post_processing.thesis_style import FIG_RC, TEXTWIDTH_IN, VIRIDIS  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================
TIER0_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_tier0"
SERIES_CSV: Path = TIER0_DIR / "tier0_B2_B3_h405_series.csv"
SCORES_CSV: Path = TIER0_DIR / "tier0_B2_B3_h405_scores.csv"
OUT_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_tier2"
# file stem -> laws drawn (draw order)
VARIANTS: dict[str, tuple[str, ...]] = {
    "tier2_B2_h405": ("B2", "h405"),
    "tier2_B2_B3_h405": ("B2", "B3", "h405"),
}
FORMATS: tuple[str, ...] = ("png",)
DPI: int = 200
SHOW_FIGURE: bool = False

FIG_T0: dict = {
    "figsize": (TEXTWIDTH_IN, 3.4),      # default; see height_by_nlaws
    # height [in] by number of laws: just tall enough for the score boxes
    # (one line per law) to sit below the traces
    "height_by_nlaws": {2: 2.9, 3: 3.2},
    "dash": {"B2": "-", "B3": "--", "h405": "-."},
    "linewidth": 1.4,
    "t_xlim": (0.0, 16.0),
    "v_xlim": (0.0, 12.0),
    "shade": "0.90",                     # scoring / exit window, trace-speed band
    "tddft_raw": dict(color="0.65", lw=0.7, alpha=0.8),   # raw |v1|, |v2| in (a) and (b)
    "tddft": dict(color="black", lw=1.5),                  # smoothed 18 Å |v1|
    "exit_bound": dict(color=VIRIDIS[0.00], ls="--", lw=1.2),
    "exit_bound_label": "Exit-speed bound",
    "titles": ("18 Å (fitted)", "9 Å (held out)", "Bulk drag law"),
    "annot_loc": (0.97, 0.04),           # lower right of (a) and (b)
    "legend_ncols": 4,
    "xlabel_t": r"Time $t$ / ps",
    "ylabel_v": r"Mean speed $|\bar{v}|$ / Å ps$^{-1}$",
    "xlabel_v": r"Speed $|v|$ / Å ps$^{-1}$",
    "ylabel_F": r"Drag force $F$ / amu Å ps$^{-2}$",
}

# The ridge figure: bulk F(v) of every m = 21 ridge point (committed
# ``g11_ridge_m21/ridge_m21.csv``), the c = 1.6 point replaced by B2.
RIDGE_PNG: Path = TIER0_DIR / "tier0_ridge_m21_laws.png"
FIG_RIDGE: dict = {
    "figsize": (TEXTWIDTH_IN, 2.9),
    "width_ratios": (1.0, 1.0),
    "v_xlim": (0.0, 12.0),
    "zoom_margin_Aps": 0.25,             # (b): the trace-speed band +- this
    "cmap_span": (0.0, 0.9),             # viridis(0 .. 0.9) over c = 0 .. c_max
    "ridge_line": dict(lw=1.0, alpha=0.55),
    "b2_line": dict(lw=2.0, zorder=5),
    "shade": "0.90",
    "titles": ("Ridge laws in bulk He", "Zoom: 18 Å trace speeds"),
    "legend_loc": "upper left",
    "cbar_label": r"Ridge $c$ / amu Å$^{-1}$",
    "xlabel": r"Speed $|v|$ / Å ps$^{-1}$",
    "ylabel": r"Drag force $F$ / amu Å ps$^{-2}$",
}

# =============================================================================
# Imports (after the settings so SHOW_FIGURE picks the backend)
# =============================================================================
import matplotlib  # noqa: E402

if not SHOW_FIGURE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import ListedColormap, Normalize  # noqa: E402

from i2_helium_md import load_drag_bundle  # noqa: E402
from i2_helium_md.extraction.trajectory_matching import REFERENCE_DRAG_ROOT  # noqa: E402
from i2_helium_md.physics.drag import drag_force  # noqa: E402
from i2_helium_md.postprocess import load_hedft_trajectory  # noqa: E402
from i2_helium_md.postprocess.hedft_loader import load_smoothed_speed_reference  # noqa: E402
from scripts.extraction.form_phase_common import load_json  # noqa: E402
from scripts.extraction.method_b_g412_ridge import (  # noqa: E402
    BOUND_VARIANT, EXIT_WINDOW_PS, G11_BOUND_DIR, G11_RIDGE_DIR, HEDFT_9A_PATH,
    REFERENCE_FILENAME,
)
from scripts.post_processing.tier0_thesis_figure import (  # noqa: E402
    GATE_A, H405_VARIANT, HEDFT_18A_PATH, LQ_LAWS, _capped,
)
from scripts.post_processing.tier2_b2_thesis_figure import (  # noqa: E402
    _annot, _finish_axis, _panel_letters, merged_style,
)
from scripts.post_processing.tier2_law_comparison_figure import LAWS  # noqa: E402
from scripts.tier0_common import M_EFF_AMU  # noqa: E402

BULK_DEPTH_A = -1e3          # deep inside the droplet: gate g = 1 (ungated bulk law)


# =============================================================================
# Data layer
# =============================================================================
def read_series(path: Path = SERIES_CSV) -> dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]:
    """``(case, law) -> (time_ps, mean_speed_Aps)`` from the Tier-0 series CSV."""
    rows: dict[tuple[str, str], list[tuple[float, float]]] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows.setdefault((r["case"], r["law"]), []).append(
                (float(r["time_ps"]), float(r["mean_speed_Aps"])))
    return {k: (np.array([t for t, _ in v]), np.array([s for _, s in v]))
            for k, v in rows.items()}


def read_scores(path: Path = SCORES_CSV) -> dict[str, dict[str, str]]:
    """``law -> row`` of the Tier-0 scores CSV (values as written)."""
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["law"]: r for r in csv.DictReader(fh)}


def law_coefficients() -> dict[str, object]:
    """The bulk drag coefficients per law, built as ``tier0_thesis_figure`` scores
    them: B2 / B3 = the lq bundle with their (a, c); h405 = its bundle capped."""
    lq = load_drag_bundle(BOUND_VARIANT, "18A", expected_m_eff_amu=M_EFF_AMU)
    out = {label: replace(lq, coefficients={"a": a, "c": c})
           for label, (a, c, _, _) in LQ_LAWS.items()}
    out["h405"] = _capped(load_drag_bundle(H405_VARIANT, "18A", expected_m_eff_amu=M_EFF_AMU))
    return out


def check_anchor(scores: dict[str, dict[str, str]],
                 series: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]) -> None:
    """B2 RMSE / margin bit-exact vs the committed G11 JSONs; every law has both series."""
    b2 = scores["B2"]
    for name, value, committed in (
            ("18 Å RMSE", float(b2["rmse_18A_Aps"]),
             load_json(G11_RIDGE_DIR / "check_B2.json")["rmse_18A_Aps"]),
            ("9 Å margin", float(b2["margin_9A_Aps"]),
             load_json(G11_BOUND_DIR / "B2.json")["margin_Aps"])):
        if value != committed:
            raise ValueError(f"B2 {name}: scores CSV {value!r} != committed {committed!r}")
    missing = [(case, law) for keys in VARIANTS.values() for law in keys
               for case in ("18A", "9A") if (case, law) not in series]
    if missing:
        raise ValueError(f"series CSV lacks {missing}")


def read_ridge(b2_coeffs, path: Path = G11_RIDGE_DIR / "ridge_m21.csv"
               ) -> list[tuple[float, float, bool]]:
    """Ridge-grid rows ``(c [amu/A], a [amu/ps], passes 9 Å)`` of ``ridge_m21.csv``.

    Only the grid rows (labels ``c<x>p<yyy>``); the B2 / A5 check rows are
    skipped. Raises unless the grid has a c = 1.6 row (replaced by B2 in the
    figure) and the CSV's B2 check row carries ``b2_coeffs`` bit-exactly.
    """
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    ridge = [(float(r["c_amu_per_A"]), float(r["a_amu_per_ps"]), r["bound_pass"] == "True")
             for r in rows if r["label"].startswith("c")]
    c_b2, a_b2 = b2_coeffs.coefficients["c"], b2_coeffs.coefficients["a"]
    if not any(c == c_b2 for c, _, _ in ridge):
        raise ValueError(f"{path.name} has no c = {c_b2} row to replace by B2")
    b2_rows = [r for r in rows if r["label"] == "B2"]
    if len(b2_rows) != 1 or (float(b2_rows[0]["c_amu_per_A"]),
                             float(b2_rows[0]["a_amu_per_ps"])) != (c_b2, a_b2):
        raise ValueError(f"{path.name}: B2 row does not carry the B2 law (c {c_b2}, a {a_b2})")
    return ridge


def exit_bound_Aps(scores: dict[str, dict[str, str]], keys: tuple[str, ...]) -> float:
    """The 9 Å exit-speed bound [A/ps]: TDDFT atom 1's mean speed over the exit
    window (``v_exit_atom1_Aps``; a law passes if its own exit mean is >= it).
    The same for every law of the scores CSV, checked."""
    values = {float(scores[k]["v_exit_atom1_Aps"]) for k in keys}
    if len(values) != 1:
        raise ValueError(f"laws disagree on the 9 Å exit bound: {sorted(values)}")
    return values.pop()


# =============================================================================
# Figure builder
# =============================================================================
def build_tier0(keys: tuple[str, ...], series, scores, coeffs, ref18, hedft18, hedft9,
                style: dict | None = None) -> plt.Figure:
    """(a) 18 Å trace, (b) 9 Å exit bound, (c) bulk F(v) for the laws ``keys``."""
    st = merged_style(FIG_T0 if style is None else style)
    figsize = (st["figsize"][0], st["height_by_nlaws"].get(len(keys), st["figsize"][1]))
    fig, (ax18, ax9, axF) = plt.subplots(1, 3, figsize=figsize, constrained_layout=True)
    names = {k: LAWS[k][1] for k in keys}
    line = {k: dict(color=LAWS[k][2], ls=st["dash"][k], lw=st["linewidth"]) for k in keys}

    # (a) 18 Å
    window = (float(ref18.time_ps[0]), float(ref18.time_ps[-1]))
    ax18.axvspan(*window, color=st["shade"], zorder=0, lw=0)
    # raw TDDFT |v1|, |v2| drawn alike in (a) and (b): one legend entry for both
    ax18.plot(hedft18.time_ps, hedft18.v1_magnitude_Aps, **st["tddft_raw"],
              label=r"TDDFT $|v_1|$, $|v_2|$")
    ax18.plot(hedft18.time_ps, hedft18.v2_magnitude_Aps, **st["tddft_raw"])
    ax18.plot(ref18.time_ps, ref18.speed_Aps, **st["tddft"], label=r"Smoothed $|v_1|$")
    for k in keys:
        ax18.plot(*series[("18A", k)], **line[k], label=names[k])
    _annot(ax18, "\n".join([r"RMSE / Å ps$^{-1}$"] + [
        f"{names[k]}: {float(scores[k]['rmse_18A_Aps']):.3f}" for k in keys]),
        st, st["annot_loc"], ha="right", va="bottom")

    # (b) 9 Å
    ax9.axvspan(*EXIT_WINDOW_PS, color=st["shade"], zorder=0, lw=0)
    ax9.plot(hedft9.time_ps, hedft9.v1_magnitude_Aps, **st["tddft_raw"])
    ax9.plot(hedft9.time_ps, hedft9.v2_magnitude_Aps, **st["tddft_raw"])
    # the bound every law's exit-window mean speed must reach (atom 1, G4.12 step 2)
    ax9.axhline(exit_bound_Aps(scores, keys), **st["exit_bound"], label=st["exit_bound_label"])
    for k in keys:
        ax9.plot(*series[("9A", k)], **line[k])
    _annot(ax9, "\n".join([r"Exit margin / Å ps$^{-1}$"] + [
        f"{names[k]}: {float(scores[k]['margin_9A_Aps']):+.2f}".replace("-", "−")
        + ("" if scores[k]["pass_9A"] == "True" else " (fails)") for k in keys]),
        st, st["annot_loc"], ha="right", va="bottom")

    # (c) bulk laws
    v = np.linspace(*st["v_xlim"], 241)
    axF.axvspan(float(np.min(ref18.speed_Aps)), float(np.max(ref18.speed_Aps)),
                color=st["shade"], zorder=0, lw=0)
    for k in keys:
        axF.plot(v, drag_force(v, np.full_like(v, BULK_DEPTH_A), coeffs[k], GATE_A), **line[k])

    for ax, title in zip((ax18, ax9, axF), st["titles"]):
        ax.set_title(title, fontsize=st["title_fontsize"])
    ax18.set(xlim=st["t_xlim"], ylim=(0.0, None), xlabel=st["xlabel_t"], ylabel=st["ylabel_v"])
    ax9.set(xlim=st["t_xlim"], ylim=(0.0, None), xlabel=st["xlabel_t"], ylabel=st["ylabel_v"])
    axF.set(xlim=st["v_xlim"], ylim=(0.0, None), xlabel=st["xlabel_v"], ylabel=st["ylabel_F"])
    for ax in (ax18, ax9, axF):
        _finish_axis(ax, st)
    entries: dict[str, object] = {}
    for ax in (ax18, ax9):
        for h, lab in zip(*ax.get_legend_handles_labels()):
            entries.setdefault(lab, h)
    # references first, then the laws in draw order
    labels = [lab for lab in entries if lab not in names.values()] + [names[k] for k in keys]
    handles = [entries[lab] for lab in labels]
    # one row while it fits, else three columns (column-major: references, then laws)
    ncols = len(labels) if len(labels) <= st["legend_ncols"] else 3
    fig.legend(handles, labels, loc="outside upper center", ncols=ncols,
               frameon=st["legend_frameon"], fontsize=st["legend_fontsize"])
    _panel_letters((ax18, ax9, axF), st)
    return fig


def build_ridge(ridge: list[tuple[float, float, bool]], b2_coeffs, ref18,
                style: dict | None = None) -> plt.Figure:
    """(a) bulk F(v) of every ridge point over 0..12 Å/ps, (b) zoom on the 18 Å
    trace speeds; ridge points coloured by c (colour bar), the 9 Å failures
    dashed, B2 on top in its thesis colour."""
    st = merged_style(FIG_RIDGE if style is None else style)
    fig, axes = plt.subplots(1, 2, figsize=st["figsize"], constrained_layout=True,
                             gridspec_kw={"width_ratios": st["width_ratios"]})
    c_max = max(c for c, _, _ in ridge)
    cmap = ListedColormap(plt.get_cmap("viridis")(np.linspace(*st["cmap_span"], 256)))
    norm = Normalize(0.0, c_max)
    v_lo, v_hi = float(np.min(ref18.speed_Aps)), float(np.max(ref18.speed_Aps))
    c_b2 = b2_coeffs.coefficients["c"]
    pad = st["zoom_margin_Aps"]
    for ax, v in ((axes[0], np.linspace(*st["v_xlim"], 241)),
                  (axes[1], np.linspace(v_lo - pad, v_hi + pad, 201))):
        bulk = np.full_like(v, BULK_DEPTH_A)
        ax.axvspan(v_lo, v_hi, color=st["shade"], zorder=0, lw=0, label="18 Å trace speeds")
        failed = False
        for c, a, bound_pass in ridge:
            if c == c_b2:
                continue
            coeffs = replace(b2_coeffs, coefficients={"a": a, "c": c})
            ax.plot(v, drag_force(v, bulk, coeffs, GATE_A), **st["ridge_line"],
                    ls="-" if bound_pass else "--", color=cmap(norm(c)),
                    label=None if bound_pass or failed else "Fails 9 Å bound")
            failed = failed or not bound_pass
        ax.plot(v, drag_force(v, bulk, b2_coeffs, GATE_A), **st["b2_line"],
                color=LAWS["B2"][2], label=rf"{LAWS['B2'][1]} ($c$ = {c_b2:g})")
        ax.set(xlim=(v[0], v[-1]), xlabel=st["xlabel"], ylabel=st["ylabel"])
        _finish_axis(ax, st)
    axes[0].set_ylim(0.0, None)
    axes[1].set_ylabel("")
    for ax, title in zip(axes, st["titles"]):
        ax.set_title(title, fontsize=st["title_fontsize"])
    axes[0].legend(loc=st["legend_loc"], frameon=st["legend_frameon"],
                   fontsize=st["legend_fontsize"])
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=list(axes), pad=0.02)
    cbar.set_label(st["cbar_label"], fontsize=st["label_fontsize"])
    cbar.ax.tick_params(labelsize=st["tick_fontsize"])
    _panel_letters(axes, st)
    return fig


# =============================================================================
# Main
# =============================================================================
def main() -> int:
    series, scores = read_series(), read_scores()
    check_anchor(scores, series)
    print("[tier0_law] anchor OK: B2 RMSE / 9 Å margin vs committed G11 JSONs")
    coeffs = law_coefficients()
    ref18 = load_smoothed_speed_reference(
        REFERENCE_DRAG_ROOT / "18A" / "velocity_smoothed" / REFERENCE_FILENAME)
    hedft18, hedft9 = load_hedft_trajectory(HEDFT_18A_PATH), load_hedft_trajectory(HEDFT_9A_PATH)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(FIG_RC):                  # thesis_style font (Palatino)
        for stem, keys in VARIANTS.items():
            fig = build_tier0(keys, series, scores, coeffs, ref18, hedft18, hedft9)
            for ext in FORMATS:
                path = OUT_DIR / f"{stem}.{ext}"
                fig.savefig(path, dpi=DPI)
                print(f"[tier0_law] wrote {path}")
        fig = build_ridge(read_ridge(coeffs["B2"]), coeffs["B2"], ref18)
        fig.savefig(RIDGE_PNG, dpi=DPI)
        print(f"[tier0_law] wrote {RIDGE_PNG}")
    if SHOW_FIGURE:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
