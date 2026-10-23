"""Tier-2 thesis figures for the production law B2 (detected ensemble).

Three figures, rendered from the pooled N = 6000 B2 container
(``…g5lqppc160e3297pooled_th500``, s1-s5 + r6) and written to
``data/figures/thesis_tier2/``:

* **B2-1 size distribution** -- (a) scored histogram (droplet-retained
  excluded; suppressed ions in the bare n = 0 bin, hatched), (b) solvated
  branch (n >= 1 renormalized, RQ8) vs the experimental abundance, with
  n̄ against the MD landing band, n₁ and W₁ annotated.
* **B2-2 mean KE vs n** -- (a) detected per-n mean KE vs the ihe_ked
  reference (correlated calibration / condition bands as envelopes, n = 1
  median anchor) with the midHot (n 2-8, geometric) and deepKE (n 10-17,
  arithmetic) bands; bins below the scorer's chi^2 min-count are drawn as
  open markers with their ion count -- they enter deepKE (committed
  ``min_count = 1``) but not chi^2; (b) profiled residuals z (the KE
  anatomy) on a shared n axis.
* **B2-3 KED shapes** -- per-n 3-D P(v) (or 2-D detector projection),
  experiment raw + movmean vs simulation, peak-normalized, in Å/ps, with
  v(<E>) markers.

Structure: a **data layer** (``*_data`` functions; numbers only, every
recipe imported -- ``tier2_confirmation`` scorers, the ``postprocess``
histogram / smoothing helpers, the ``summary_sections`` recipe constants;
nothing re-derived), a **draw layer** (``draw_*`` functions; draw into a
given axis, never create figures) and three figure builders. All styling
lives in ``COMMON_STYLE`` plus one override dict per figure
(``FIG_B2_1`` / ``FIG_B2_2`` / ``FIG_B2_3``); an override wins key by key
(nested dicts merged one level deep).

Anchor guard: before drawing, every quoted scalar (n̄, suppressed, n₁,
W₁, midHot, deepKE, chi^2_med, KE₁, number of chi^2 bins, scored count) is
checked against the committed B2 row of
``data/runs/h2b_forward_model/g5_winners_table.csv`` to 4 dp; a mismatch
raises.

Edit the ``USER SETTINGS`` block and run::

    python scripts/post_processing/tier2_b2_thesis_figure.py
"""

from __future__ import annotations

import copy
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT_DIR = Path(__file__).resolve().parent
for _p in (PROJECT_ROOT, _SCRIPT_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

# =============================================================================
# USER SETTINGS
# =============================================================================
RUN_DIR: Path = (
    PROJECT_ROOT / "data" / "runs"
    / "9A_drag_srho3p5A18A_linear_quadratic_N6000_detfix_conf270_g5lqppc160e3297pooled_th500"
)
WINNERS_CSV: Path = PROJECT_ROOT / "data" / "runs" / "h2b_forward_model" / "g5_winners_table.csv"
WINNERS_LABEL: str = "B2_pp_c160"
ANCHOR_CHECK: bool = True
IHE_KED_REFERENCE_DIR: Path = PROJECT_ROOT / "data" / "reference" / "ihe_ked"
ABUNDANCE_REF_PATH: Path = PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"

OUT_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_tier2"
FIGURES: tuple[str, ...] = ("b2_1", "b2_2", "b2_3")
FORMATS: tuple[str, ...] = ("png", "pdf")
DPI: int = 200
SHOW_FIGURE: bool = False
MPLSTYLE: Path | None = None     # optional .mplstyle applied to every figure

# Shared defaults; each FIG_* dict below overrides key by key.
COMMON_STYLE: dict = {
    "exp": dict(color="black"),
    "exp_light": dict(color="0.55"),
    "sim": dict(color="#1f77b4"),        # B2 colour of the Tier-0 / G6 figures
    "band": dict(color="0.80"),          # midHot / deepKE / landing-band shading
    "band_shade_alpha": 0.25,
    "label_fontsize": 10,
    "tick_fontsize": 9,
    "legend_fontsize": 8,
    "annot_fontsize": 8,
    "panel_letters": True,
    "panel_letter_fmt": "({})",
    "legend_frameon": False,
}

FIG_B2_1: dict = {
    "stem": "b2_1_size_distribution",
    "layout": "1x2",                     # "1x2" | "2x1" | "solvated_only"
    "figsize": (10.0, 3.8),
    "n_xlim": (-0.7, 21.0),
    "bar_width_scored": 0.8,
    "bar_width_pair": 0.4,
    "suppressed_hatch": "///",
    "show_landing_band": True,
    "annot_scored": True,
    "annot_solvated": True,
    "annot_loc": (0.97, 0.95),           # axes fraction, top-right anchor
    "legend_loc": "center right",
    "ylabel_scored": "scored fraction / %",
    "ylabel_solvated": "fraction (n ≥ 1) / %",
    "xlabel": "n (attached He atoms)",
}

FIG_B2_2: dict = {
    "stem": "b2_2_mean_ke",
    "figsize": (7.0, 5.6),
    "height_ratios": (2.6, 1.0),
    "show_residual": True,
    "yscale": "log",                     # "log" | "linear"
    "ylim": None,
    "n_xlim": (0.4, 17.6),
    "show_bands": True,                  # correlated calib / condition envelopes
    "band_alphas": (0.35, 0.20),
    "show_median_anchor": True,
    "show_midhot_band": True,
    "show_deepke_band": True,
    "show_thin_bin_counts": True,
    "annot_chi2": True,
    "annot_loc": (0.03, 0.05),
    "legend_loc": "upper right",
    "z_guides": (-2.0, 2.0),
    "z_colors": ("#d62728", "#1f77b4"),  # (sim hot, sim cold)
    "ylabel": r"$\langle E\rangle$ / eV",
    "ylabel_z": r"$z$",
    "xlabel": "n (attached He atoms)",
}

FIG_B2_3: dict = {
    "stem": "b2_3_ked_shapes",
    "representation": "3d",              # "3d" (P(v)) | "2d" (detector projection)
    "n_values": (1, 2, 3, 4),
    "ncols": 2,
    "figsize": (8.0, 6.0),
    "v_xlim_Aps": (0.0, 26.0),
    "ylim": (0.0, 1.25),
    "show_exp_raw": True,
    "show_v_of_mean_E": True,
    "show_counts": True,
    "annot_ke_n": (1,),                  # panels that print <E> sim / exp
    "legend_panel": "figure",            # "figure" (one legend on top) | panel index | None
    "legend_loc": "upper right",         # used for a panel legend only
    "xlabel": r"$v$ / Å ps$^{-1}$",
    "ylabel": "P(v) / peak-normalized",
}

# =============================================================================
# Imports (after the settings so SHOW_FIGURE picks the backend)
# =============================================================================
import matplotlib  # noqa: E402

if not SHOW_FIGURE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

from i2_helium_md.physics.shell_schedule import complex_mass_amu  # noqa: E402
from i2_helium_md.postprocess import (  # noqa: E402
    compute_final_velocity_histogram,
    detected_ensemble_view,
    fragment_mean_kinetic_energy,
    load_he_abundance_reference,
    load_ihe_ked_curve,
    load_ihe_ked_reference,
    speed_mps_of_energy_eV,
)
from i2_helium_md.postprocess._smoothing import moving_mean, normalise_trace  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    DEEP_KE_BAND,
    MIDHOT_BAND,
    deep_bin_ke_ratio,
    midhot_ratio,
    read_confirmation_detection,
    score_histogram_vs_reference,
    score_ke_curve_vs_reference,
    solvated_renormalized,
)
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from scripts.post_processing.tier2_g414_table import NBAR_LANDING  # noqa: E402
from summary_sections import (  # noqa: E402
    DETECTED_KE_INCLUDE_SIM_SE,
    DETECTED_KE_MIN_BIN_COUNT,
    DETECTED_KE_N1_ANCHOR,
    DETECTED_KE_N_MIN,
    HIST_SMOOTHING_WINDOW,
    IHE_KED_EXP_SMOOTHING_WINDOW,
    IHE_KED_HIST_NUM_BINS,
    IHE_KED_HIST_V_MAX_APS,
    _nan_aware_moving_mean,
)

MPS_PER_APS = 100.0   # 1 Å/ps = 100 m/s (exact)


# =============================================================================
# Style helpers
# =============================================================================
def merged_style(override: dict, common: dict = COMMON_STYLE) -> dict:
    """``common`` updated by ``override``; nested dicts merged one level deep."""
    out = copy.deepcopy(common)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = {**out[key], **value}
        else:
            out[key] = copy.deepcopy(value)
    return out


def _finish_axis(ax, st: dict, integer_x: bool = False) -> None:
    if integer_x:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.tick_params(labelsize=st["tick_fontsize"])
    ax.xaxis.label.set_size(st["label_fontsize"])
    ax.yaxis.label.set_size(st["label_fontsize"])


def _panel_letters(axes, st: dict) -> None:
    if not st["panel_letters"]:
        return
    for i, ax in enumerate(axes):
        ax.text(0.0, 1.02, st["panel_letter_fmt"].format(chr(ord("a") + i)),
                transform=ax.transAxes, ha="left", va="bottom",
                fontsize=st["label_fontsize"], fontweight="bold")


def _legend(ax, st: dict, loc: str) -> None:
    ax.legend(loc=loc, frameon=st["legend_frameon"], fontsize=st["legend_fontsize"])


# =============================================================================
# Data layer
# =============================================================================
@dataclass(frozen=True)
class SizeData:
    """Scored + solvated size histograms (fractions, not %) and their scores."""
    n_scored: np.ndarray
    frac_scored: np.ndarray
    n_sim_solv: np.ndarray
    frac_sim_solv: np.ndarray
    n_ref_solv: np.ndarray
    frac_ref_solv: np.ndarray
    n_mean: float
    suppressed_frac: float
    num_scored: int
    num_ions: int
    w1_solvated: float
    sim_n1: float
    ref_n1: float


def size_data(read, abundance) -> SizeData:
    """Scorer read path of ``_section_detected_size_distribution``."""
    score = score_histogram_vs_reference(read, abundance)
    sim_n, sim_f = solvated_renormalized(read.n_values, read.fraction, n_min=score.n_min)
    ref_n, ref_f = solvated_renormalized(abundance.n, abundance.ion_fraction, n_min=score.n_min)
    return SizeData(
        n_scored=np.asarray(read.n_values), frac_scored=np.asarray(read.fraction),
        n_sim_solv=sim_n, frac_sim_solv=sim_f, n_ref_solv=ref_n, frac_ref_solv=ref_f,
        n_mean=float(read.n_mean), suppressed_frac=float(read.suppressed_frac),
        num_scored=int(read.num_scored), num_ions=int(read.num_ions),
        w1_solvated=float(score.w1_solvated), sim_n1=float(score.sim_n1),
        ref_n1=float(score.ref_n1),
    )


@dataclass(frozen=True)
class KEData:
    """Per-n mean KE (all occupied bins, ``min_count = 1``), reference, scores.

    ``scored_mask`` marks the bins that pass the chi^2 min-count
    (``DETECTED_KE_MIN_BIN_COUNT``); the others are display-only thin bins.
    """
    n: np.ndarray
    mean_eV: np.ndarray
    se_eV: np.ndarray
    count: np.ndarray
    scored_mask: np.ndarray
    ref: object
    chi2_score: object
    midhot: object
    deepke: object


def ke_data(read, ked_ref) -> KEData:
    """Same scorer kwargs as ``summary_sections._detected_ke_scores``;
    midHot / deepKE at their committed conventions (geometric / arithmetic,
    ``min_count = 1``) -- the values in ``g5_winners_table.csv``."""
    ke = read.ke_by_n(min_count=1)
    keep = ke.n >= DETECTED_KE_N_MIN
    chi2 = score_ke_curve_vs_reference(
        read, ked_ref, n_min=DETECTED_KE_N_MIN, min_count=DETECTED_KE_MIN_BIN_COUNT,
        include_sim_se=DETECTED_KE_INCLUDE_SIM_SE, n1_anchor=DETECTED_KE_N1_ANCHOR,
    )
    count = np.asarray(ke.count)[keep]
    return KEData(
        n=np.asarray(ke.n)[keep], mean_eV=np.asarray(ke.mean_eV)[keep],
        se_eV=np.asarray(ke.se_eV)[keep], count=count,
        scored_mask=count >= DETECTED_KE_MIN_BIN_COUNT, ref=ked_ref,
        chi2_score=chi2, midhot=midhot_ratio(read, ked_ref),
        deepke=deep_bin_ke_ratio(read, ked_ref),
    )


@dataclass(frozen=True)
class PvData:
    """One n-panel of B2-3; speeds in Å/ps, traces peak-normalized."""
    n: int
    v_exp_Aps: np.ndarray
    exp_raw: np.ndarray
    exp_smooth: np.ndarray
    v_sim_Aps: np.ndarray | None
    sim_trace: np.ndarray | None
    num_sim: int
    v_mean_E_exp_Aps: float
    v_mean_E_sim_Aps: float | None
    mean_E_exp_eV: float
    mean_E_sim_eV: float | None


def pv_data(view, ked_dir: Path, ked_ref, n: int, representation: str) -> PvData:
    """Recipe of ``summary_sections._section_ihe_ked_curves`` for one n,
    converted to Å/ps (display unit only; histogram bins are the Å/ps
    ``IHE_KED_HIST_*`` grid already)."""
    if representation not in ("2d", "3d"):
        raise ValueError(f"representation must be '2d' or '3d', got {representation!r}")
    curve = load_ihe_ked_curve(ked_dir, n)
    raw = np.asarray(getattr(curve, f"signal_{representation}_Pv"), dtype=float)
    mass_amu = float(complex_mass_amu(n))
    ref_idx = np.flatnonzero(np.asarray(ked_ref.n) == n)
    if ref_idx.size != 1:
        raise ValueError(f"ihe_ked reference has no unique row for n = {n}")
    mean_exp = float(ked_ref.mean_KE_eV[ref_idx[0]])
    try:
        hist = compute_final_velocity_histogram(
            view, mass_amu=mass_amu, num_bins=IHE_KED_HIST_NUM_BINS,
            v_max_Aps=IHE_KED_HIST_V_MAX_APS, projected=(representation == "2d"),
        )
    except ValueError:
        hist = None
    try:
        sim_mean = fragment_mean_kinetic_energy(view, n)
    except ValueError:
        sim_mean = None
    return PvData(
        n=n,
        v_exp_Aps=np.asarray(curve.v_mps, dtype=float) / MPS_PER_APS,
        exp_raw=raw,
        exp_smooth=_nan_aware_moving_mean(raw, IHE_KED_EXP_SMOOTHING_WINDOW),
        v_sim_Aps=None if hist is None else hist.bin_centers_Aps,
        sim_trace=None if hist is None else normalise_trace(
            moving_mean(hist.density, HIST_SMOOTHING_WINDOW)),
        num_sim=0 if hist is None else int(hist.num_atoms_used),
        v_mean_E_exp_Aps=float(speed_mps_of_energy_eV(mean_exp, mass_amu)) / MPS_PER_APS,
        v_mean_E_sim_Aps=None if sim_mean is None else sim_mean.v_of_mean_E_mps / MPS_PER_APS,
        mean_E_exp_eV=mean_exp,
        mean_E_sim_eV=None if sim_mean is None else float(sim_mean.mean_KE_eV),
    )


# =============================================================================
# Anchor guard
# =============================================================================
def anchor_values(size: SizeData, ke: KEData) -> dict[str, float]:
    """The scalars the figures quote, keyed by ``g5_winners_table.csv`` column."""
    ke1 = ke.mean_eV[ke.n == 1]
    return {
        "num_scored": size.num_scored,
        "supp": size.suppressed_frac,
        "nbar_det": size.n_mean,
        "n1_solv": size.sim_n1,
        "w1_solv": size.w1_solvated,
        "midHot": ke.midhot.value,
        "midHot_bins": ke.midhot.n_bins,
        "deepKE": ke.deepke.value,
        "deepKE_bins": ke.deepke.n_bins,
        "chi2_med": ke.chi2_score.chi2_profiled,
        "ke_npts": ke.chi2_score.n_points,
        "KE1_mean": float(ke1[0]) if ke1.size else float("nan"),
    }


def check_anchor(values: dict[str, float], csv_path: Path, label: str,
                 container: str | None = None) -> None:
    """Raise ``ValueError`` unless every value matches the committed row to 4 dp."""
    with open(csv_path, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh) if r["label"] == label]
    if len(rows) != 1:
        raise ValueError(f"{csv_path.name}: expected one row {label!r}, found {len(rows)}")
    row = rows[0]
    if container is not None and row.get("container") != container:
        raise ValueError(f"{label!r} row belongs to {row.get('container')!r}, not {container!r}")
    bad = [f"{k}: figure {v:.4f} vs committed {float(row[k]):.4f}"
           for k, v in values.items() if round(float(v), 4) != round(float(row[k]), 4)]
    if bad:
        raise ValueError("anchor mismatch vs committed B2 row:\n  " + "\n  ".join(bad))


# =============================================================================
# Draw layer
# =============================================================================
def _annot(ax, text: str, st: dict, loc, ha="right", va="top") -> None:
    ax.text(*loc, text, transform=ax.transAxes, ha=ha, va=va, fontsize=st["annot_fontsize"],
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="0.7", alpha=0.9))


def draw_size_scored(ax, d: SizeData, st: dict) -> None:
    """Scored histogram; the n = 0 (suppressed / bare-candidate) bar hatched."""
    bare = d.n_scored == 0
    ax.bar(d.n_scored[~bare], 100.0 * d.frac_scored[~bare], width=st["bar_width_scored"],
           edgecolor="black", linewidth=0.4, **st["sim"])
    ax.bar(d.n_scored[bare], 100.0 * d.frac_scored[bare], width=st["bar_width_scored"],
           facecolor="white", edgecolor=st["sim"]["color"], hatch=st["suppressed_hatch"],
           linewidth=0.8, label="suppressed (n = 0)")
    if st["annot_scored"]:
        lo, hi = NBAR_LANDING
        flag = "∈" if lo <= d.n_mean <= hi else "∉"
        _annot(ax, rf"$\bar n$ = {d.n_mean:.2f} {flag} [{lo:.2f}, {hi:.2f}]" "\n"
                   f"suppressed {100.0 * d.suppressed_frac:.1f} %\n"
                   f"scored {d.num_scored} / {d.num_ions} ions", st, st["annot_loc"])
    if st["show_landing_band"]:
        ax.axvspan(*NBAR_LANDING, zorder=0, alpha=2 * st["band_shade_alpha"], linewidth=0,
                   label=r"$\bar n$ landing band", **st["band"])
        ax.axvline(d.n_mean, color=st["sim"]["color"], ls="--", lw=1.0, label=r"$\bar n$")
    ax.set(xlim=st["n_xlim"], xlabel=st["xlabel"], ylabel=st["ylabel_scored"])
    _legend(ax, st, st["legend_loc"])
    _finish_axis(ax, st, integer_x=True)


def draw_size_solvated(ax, d: SizeData, st: dict) -> None:
    """Solvated branch, experiment vs B2, side-by-side bars."""
    w = st["bar_width_pair"]
    ax.bar(d.n_ref_solv - w / 2, 100.0 * d.frac_ref_solv, width=w, label="experiment",
           **st["exp"])
    ax.bar(d.n_sim_solv + w / 2, 100.0 * d.frac_sim_solv, width=w, label="B2", **st["sim"])
    if st["annot_solvated"]:
        _annot(ax, rf"$n_1$ = {d.sim_n1:.3f} (exp {d.ref_n1:.3f})" "\n"
                   rf"$W_1$ = {d.w1_solvated:.3f} bins", st, st["annot_loc"])
    ax.set(xlim=st["n_xlim"], xlabel=st["xlabel"], ylabel=st["ylabel_solvated"])
    _legend(ax, st, st["legend_loc"])
    _finish_axis(ax, st, integer_x=True)


def draw_ke_mean(ax, d: KEData, st: dict) -> None:
    """Mean KE vs n: reference (bands, error bars, median anchor) + B2."""
    ref = d.ref
    sel = np.asarray(ref.n) >= DETECTED_KE_N_MIN
    n_ref = np.asarray(ref.n)[sel]
    mean = np.asarray(ref.mean_KE_eV)[sel]
    if st["show_midhot_band"]:
        ax.axvspan(MIDHOT_BAND[0] - 0.5, MIDHOT_BAND[1] + 0.5, zorder=0,
                   alpha=st["band_shade_alpha"], linewidth=0, **st["band"])
        ax.text(0.5 * sum(MIDHOT_BAND), 1.0, f"midHot {d.midhot.value:.3f}",
                transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=st["annot_fontsize"])
    if st["show_deepke_band"]:
        ax.axvspan(DEEP_KE_BAND[0] - 0.5, DEEP_KE_BAND[1] + 0.5, zorder=0,
                   alpha=st["band_shade_alpha"], linewidth=0, **st["band"])
        ax.text(0.5 * sum(DEEP_KE_BAND), 1.0, f"deepKE {d.deepke.value:.3f}",
                transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=st["annot_fontsize"])
    if st["show_bands"]:
        for frac, alpha in zip((ref.calib_syst_frac, ref.condition_syst_frac), st["band_alphas"]):
            frac = np.asarray(frac)[sel]
            ax.fill_between(n_ref, mean * (1.0 - frac), mean * (1.0 + frac), linewidth=0,
                            alpha=alpha, **st["exp_light"])
    ax.errorbar(n_ref, mean, yerr=np.asarray(ref.point_err_eV)[sel], fmt="o", ms=4, capsize=2,
                label="experiment", **st["exp"])
    if st["show_median_anchor"]:
        hit = n_ref == 1
        ax.plot(n_ref[hit], np.asarray(ref.median_KE_eV)[sel][hit], "D", ms=6,
                mfc="none", mew=1.2, color=st["exp"]["color"], label="experiment median (n = 1)")
    ok = d.scored_mask
    ax.errorbar(d.n[ok], d.mean_eV[ok], yerr=d.se_eV[ok], fmt="s", ms=4.5, capsize=2,
                label="B2", **st["sim"])
    if np.any(~ok):
        ax.plot(d.n[~ok], d.mean_eV[~ok], "s", ms=4.5, mfc="none", mew=1.0,
                color=st["sim"]["color"],
                label=f"B2, < {DETECTED_KE_MIN_BIN_COUNT} ions (not in $\\chi^2$)")
        if st["show_thin_bin_counts"]:
            for n, e, c in zip(d.n[~ok], d.mean_eV[~ok], d.count[~ok]):
                ax.annotate(f"{c}", (n, e), xytext=(0, 5), textcoords="offset points",
                            ha="center", va="bottom", fontsize=st["annot_fontsize"] - 1,
                            color=st["sim"]["color"])
    if st["annot_chi2"]:
        _annot(ax, rf"$\chi^2_\mathrm{{med}}$ = {d.chi2_score.chi2_profiled:.1f}"
                   f" ({d.chi2_score.n_points} bins)", st, st["annot_loc"], ha="left", va="bottom")
    ax.set(yscale=st["yscale"], xlim=st["n_xlim"], ylabel=st["ylabel"])
    if st["ylim"] is not None:
        ax.set_ylim(st["ylim"])
    _legend(ax, st, st["legend_loc"])
    _finish_axis(ax, st, integer_x=True)


def draw_ke_residual(ax, d: KEData, st: dict) -> None:
    """Profiled residuals z (scored chi^2 bins only)."""
    s = d.chi2_score
    hot, cold = st["z_colors"]
    ax.bar(s.n, s.z_profiled, width=0.8, edgecolor="black", linewidth=0.4,
           color=[hot if z > 0 else cold for z in s.z_profiled])
    ax.axhline(0.0, color="black", lw=0.8)
    for lvl in st["z_guides"]:
        ax.axhline(lvl, color="0.5", lw=0.8, ls="--")
    ax.set(xlim=st["n_xlim"], xlabel=st["xlabel"], ylabel=st["ylabel_z"])
    _finish_axis(ax, st, integer_x=True)


def draw_pv(ax, d: PvData, st: dict) -> None:
    """One KED-shape panel (Å/ps)."""
    if st["show_exp_raw"]:
        ax.plot(d.v_exp_Aps, d.exp_raw, ".-", ms=2.0, lw=0.4, alpha=0.35,
                label="experiment (raw)", **st["exp_light"])
    ax.plot(d.v_exp_Aps, d.exp_smooth, lw=1.4,
            label=f"experiment (movmean {IHE_KED_EXP_SMOOTHING_WINDOW})", **st["exp"])
    if d.sim_trace is not None:
        ax.plot(d.v_sim_Aps, d.sim_trace, lw=1.4, label="B2", **st["sim"])
    if st["show_v_of_mean_E"]:
        ax.axvline(d.v_mean_E_exp_Aps, ls=":", lw=1.0, **st["exp"])
        if d.v_mean_E_sim_Aps is not None:
            ax.axvline(d.v_mean_E_sim_Aps, ls=":", lw=1.0, **st["sim"])
    lines = [rf"I$^+$He$_{{{d.n}}}$" if d.n > 1 else r"I$^+$He"]
    if st["show_counts"]:
        lines.append(f"N = {d.num_sim}")
    if d.n in st["annot_ke_n"] and d.mean_E_sim_eV is not None:
        lines.append(rf"$\langle E\rangle$ {d.mean_E_sim_eV:.3f} / exp {d.mean_E_exp_eV:.3f} eV")
    ax.text(0.03, 0.95, "\n".join(lines), transform=ax.transAxes, ha="left", va="top",
            fontsize=st["annot_fontsize"])
    ax.set(xlim=st["v_xlim_Aps"], ylim=st["ylim"])
    _finish_axis(ax, st)


# =============================================================================
# Figure builders
# =============================================================================
def build_b2_1(d: SizeData, style: dict | None = None) -> plt.Figure:
    st = merged_style(FIG_B2_1 if style is None else style)
    layout = st["layout"]
    if layout == "solvated_only":
        fig, ax = plt.subplots(figsize=st["figsize"], constrained_layout=True)
        draw_size_solvated(ax, d, st)
        return fig
    if layout not in ("1x2", "2x1"):
        raise ValueError(f"unknown B2-1 layout {layout!r}")
    shape = (1, 2) if layout == "1x2" else (2, 1)
    fig, axes = plt.subplots(*shape, figsize=st["figsize"], constrained_layout=True)
    draw_size_scored(axes[0], d, st)
    draw_size_solvated(axes[1], d, st)
    _panel_letters(axes, st)
    return fig


def build_b2_2(d: KEData, style: dict | None = None) -> plt.Figure:
    st = merged_style(FIG_B2_2 if style is None else style)
    if not st["show_residual"]:
        fig, ax = plt.subplots(figsize=st["figsize"], constrained_layout=True)
        draw_ke_mean(ax, d, st)
        ax.set_xlabel(st["xlabel"])
        return fig
    fig, (ax, axz) = plt.subplots(2, 1, sharex=True, figsize=st["figsize"],
                                  constrained_layout=True,
                                  gridspec_kw={"height_ratios": st["height_ratios"]})
    draw_ke_mean(ax, d, st)
    draw_ke_residual(axz, d, st)
    _panel_letters((ax, axz), st)
    return fig


def build_b2_3(panels: list[PvData], style: dict | None = None) -> plt.Figure:
    st = merged_style(FIG_B2_3 if style is None else style)
    ncols = st["ncols"]
    nrows = -(-len(panels) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=st["figsize"], sharex=True, sharey=True,
                             constrained_layout=True, squeeze=False)
    flat = axes.ravel()
    for ax, d in zip(flat, panels):
        draw_pv(ax, d, st)
    for ax in flat[len(panels):]:
        ax.set_axis_off()
    for ax in axes[-1]:
        ax.set_xlabel(st["xlabel"])
    for ax in axes[:, 0]:
        ax.set_ylabel(st["ylabel"])
    if st["legend_panel"] == "figure":
        handles, labels = flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside upper center", ncols=len(labels),
                   frameon=st["legend_frameon"], fontsize=st["legend_fontsize"])
    elif st["legend_panel"] is not None:
        _legend(flat[st["legend_panel"]], st, st["legend_loc"])
    _panel_letters(flat[:len(panels)], st)
    return fig


# =============================================================================
# Main
# =============================================================================
def _save(fig: plt.Figure, stem: str) -> None:
    for ext in FORMATS:
        path = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(path, dpi=DPI)
        print(f"[b2_thesis] wrote {path}")


def main() -> int:
    det = load_detection_result(RUN_DIR / "detection.npz")
    read = read_confirmation_detection(det, label=RUN_DIR.name)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_DIR / "IHe_KED_reference.csv")
    abundance = load_he_abundance_reference(ABUNDANCE_REF_PATH)

    size = size_data(read, abundance)
    ke = ke_data(read, ked_ref)
    if ANCHOR_CHECK:
        check_anchor(anchor_values(size, ke), WINNERS_CSV, WINNERS_LABEL, RUN_DIR.name)
        print(f"[b2_thesis] anchor OK vs {WINNERS_CSV.name}:{WINNERS_LABEL}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    style_ctx = plt.style.context(str(MPLSTYLE)) if MPLSTYLE else plt.style.context({})
    with style_ctx:
        if "b2_1" in FIGURES:
            _save(build_b2_1(size), FIG_B2_1["stem"])
        if "b2_2" in FIGURES:
            _save(build_b2_2(ke), FIG_B2_2["stem"])
        if "b2_3" in FIGURES:
            view = detected_ensemble_view(det)
            panels = [pv_data(view, IHE_KED_REFERENCE_DIR, ked_ref, n,
                              FIG_B2_3["representation"]) for n in FIG_B2_3["n_values"]]
            _save(build_b2_3(panels), FIG_B2_3["stem"])
    if SHOW_FIGURE:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
