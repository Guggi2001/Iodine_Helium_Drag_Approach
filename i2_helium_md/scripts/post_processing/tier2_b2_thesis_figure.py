"""Tier-2 thesis figures for the production law B2 (detected ensemble).

Four figures, rendered from the pooled N = 6000 B2 thesis container
(``B2_N6000_pooled``, s1-s6, ``self_bound_gate_release = decoupled``,
THESIS_FINAL_GOALS G13) at the thesis text width
(``TEXTWIDTH_CM``) and written to ``data/figures/thesis_tier2/``. Pipeline
vocabulary (scored / detected read, median anchor provenance, chi^2,
smoothing window) lives in the captions, not on the figures:

* **B2-1 size distribution** -- (a) scored histogram (droplet-retained
  excluded; bare I+ in the n = 0 bin, bare fraction annotated), (b) solvated
  branch (n >= 1 renormalized, RQ8) vs the experimental abundance, with
  n̄ against the MD landing band, n₁, W₁ and the n >= 15 mass annotated.
* **B2-2 mean KE vs n** -- (a) n = 1-4 and (b) n = 5-17, linear, each on
  its own scale -- correlated calibration / condition bands as envelopes,
  per-point stat (+) sys bars, n = 1 median, sim bins with SE (open
  markers: < chi^2 min-count ions; they enter the committed deepKE,
  ``min_count = 1``); (c) per-bin ratio sim / experiment mean
  (``ke_band_ratio``, the midHot / deepKE ingredient) with the midHot
  (n 2-8, geometric) and deepKE (n 10-17, arithmetic) band means as
  segments. ``lower_panel = "z"`` restores the profiled-residual panel.
* **B2-3 KED shapes** -- the detection-summary curve panels for n = 1-6
  (3 x 2), 3-D P(v), experiment raw + movmean vs simulation,
  peak-normalized, v(<E>) markers; speed axis in Å/ps, per-n range.
* **B2-4** -- B2-3 on the 2-D detector projection (discussion figure,
  ``THESIS_DISCUSSION.md`` §13; 2-D is a diagnostic, not a score).

Structure: a **data layer** (``*_data`` functions; numbers only, every
recipe imported -- ``tier2_confirmation`` scorers, the ``postprocess``
histogram / smoothing helpers, the ``summary_sections`` recipe constants;
nothing re-derived), a **draw layer** (``draw_*`` functions; draw into a
given axis, never create figures) and three figure builders. All styling
lives in ``COMMON_STYLE`` plus one override dict per figure
(``FIG_B2_1`` / ``FIG_B2_2`` / ``FIG_B2_3`` / ``FIG_B2_4``, the 2-D twin of
B2-3); an override wins key by key
(nested dicts merged one level deep).

Anchor guard: before drawing, every quoted scalar (n̄, suppressed, n₁,
W₁, n >= 15 mass, midHot, deepKE, chi^2_med, KE₁, number of chi^2 bins,
scored count) is
checked against the decoupled B2 row of
``data/runs/h2b_forward_model/g13_decoupled_table.csv`` (v2 scorer, the
corrected ihe_ked reference, G12) to 4 dp; a mismatch raises.

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
RUN_DIR: Path = PROJECT_ROOT / "data" / "runs" / "B2_N6000_pooled"
# G13: decoupled rows, v2 scorer (corrected ihe_ked reference, G12)
WINNERS_CSV: Path = PROJECT_ROOT / "data" / "runs" / "h2b_forward_model" / "g13_decoupled_table.csv"
WINNERS_LABEL: str = "B2"
ANCHOR_CHECK: bool = True
IHE_KED_REFERENCE_DIR: Path = PROJECT_ROOT / "data" / "reference" / "ihe_ked"
ABUNDANCE_REF_PATH: Path = PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"

OUT_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_tier2"
FIGURES: tuple[str, ...] = ("b2_1", "b2_2", "b2_2_split", "b2_3", "b2_4")
FORMATS: tuple[str, ...] = ("png",)
DPI: int = 200
SHOW_FIGURE: bool = False
MPLSTYLE: Path | None = None     # optional .mplstyle applied to every figure

from scripts.post_processing.thesis_style import (  # noqa: E402  (shared thesis style)
    EXP_COLOR, FIG_RC, FONT, SIM_COLOR, TEXTWIDTH_CM, TEXTWIDTH_IN, VIRIDIS,
)

# Shared defaults; each FIG_* dict below overrides key by key.
COMMON_STYLE: dict = {
    "exp": dict(color=EXP_COLOR),        # thesis_style: experiment warm accent,
    "sim": dict(color=SIM_COLOR),        # simulation viridis blue
    "band": dict(color="0.80"),          # landing-band shading
    "band_shade_alpha": 0.25,
    **FONT,
    "panel_letters": True,
    "panel_letter_fmt": "({})",
    "legend_frameon": False,
}

FIG_B2_1: dict = {
    "stem": "b2_1_size_distribution",
    "layout": "1x2",                     # "1x2" | "2x1" | "solvated_only"
    "figsize": (TEXTWIDTH_IN, 2.8),
    "n_xlim": (-0.7, 21.0),
    "bar_width_scored": 0.8,
    "bar_width_pair": 0.4,
    "show_landing_band": True,
    "annot_scored": True,
    "annot_solvated": True,
    "tail_n": 15,                        # n >= tail_n solvated mass (committed column)
    "annot_loc": (0.97, 0.95),           # axes fraction, top-right anchor
    "legend_loc": "center right",
    "ylabel_scored": "Fraction / %",
    "ylabel_solvated": r"Fraction ($n \geq 1$) / %",
    "xlabel": r"Shell size $n$",
}

FIG_B2_2: dict = {
    # "stacked" (default): (a) <E> vs n = 1..n_max on a log axis, full
    # width, over (b) the lower panel, shared n axis. "split": (a) n =
    # 1..split-1, (b) n = split..n_max (linear, own scales) -- the
    # plot_detection_summary mean-KE section -- over (c). Lower panel:
    # "ratio" (default) or the profiled-residual anatomy "z".
    "stem": "b2_2_mean_ke",
    "layout": "stacked",                 # "stacked" | "split"
    "figsize": (TEXTWIDTH_IN, 5.0),
    "height_ratios": (1.2, 1.0),
    "figsize_split": (TEXTWIDTH_IN, 4.6),
    "height_ratios_split": (1.15, 1.0),
    "anchor": dict(color=EXP_COLOR),     # experimental n = 1 median: experiment colour, own marker
    "split_n": 5,                        # = summary_sections.IHE_KED_MEAN_ENERGY_SPLIT_N
    "yscale_stacked": "log",
    "yscale": "linear",                  # split layout
    "ylim_low": None,
    "ylim_high": None,
    "show_bands": True,                  # correlated calib / condition envelopes
    "band_alphas": (0.20, 0.12),
    "show_median_anchor": True,
    "show_thin_bins": True,              # bins < chi^2 min-count (they enter deepKE only)
    "show_midhot_band": False,
    "show_deepke_band": False,
    "annot_chi2": False,
    "annot_loc": (0.03, 0.05),
    "legend_loc": "upper right",
    "titles": ("$n$ = 1–{hi}", "$n$ = {lo}–{n_max}", None),
    "lower_panel": "ratio",              # "ratio" | "z"
    "ratio_ylim": None,
    "annot_band_means": True,            # midHot / deepKE segments + values (legend)
    "legend_loc_ratio": "upper right",
    "ylabel_ratio": r"$\bar{E}_\mathrm{sim} / \bar{E}_\mathrm{exp}$",
    "z_guides": (-2.0, 2.0),
    "z_colors": (SIM_COLOR, VIRIDIS[0.85]),  # (sim hot, sim cold)
    "annot_z": True,
    "ylabel": "Mean kinetic energy / eV",
    "ylabel_z": r"$z$ (sim − ref, per-point $\sigma$)",
    "xlabel": r"Shell size $n$",
}

# B2-2 split variant: linear n = 1-4 | 5-17 over the same lower panel, at
# the stacked figure's size for a side-by-side comparison.
FIG_B2_2_SPLIT: dict = {
    **copy.deepcopy(FIG_B2_2),
    "stem": "b2_2_mean_ke_split",
    "layout": "split",
    "figsize_split": FIG_B2_2["figsize"],
    "height_ratios_split": FIG_B2_2["height_ratios"],
    "legend_loc_ratio": "lower center",  # the shorter (c) has no room above the 1.0 band
}

FIG_B2_3: dict = {
    # The plot_detection_summary ihe_ked curve section for n = 1..6 (3 x 2)
    # with the speed axis in Å/ps; n = 0 (bare I+) is a different
    # population and is left out.
    "stem": "b2_3_ked_shapes",
    "representation": "3d",              # "3d" (P(v)) | "2d" (detector projection)
    "n_values": (1, 2, 3, 4, 5, 6),      # n 5-7: extended exports, not verified
    "ncols": 2,
    "figsize": (TEXTWIDTH_IN, 7.0),
    "sharex": False,
    "sharey": True,
    "v_xlim_Aps": (0.0, 31.0),           # IHE_KED_PLOT_V_MAX_MPS / 100 (fallback)
    "v_xlim_by_n": {1: (0.0, 26.0), 2: (0.0, 20.0), 3: (0.0, 16.0),
                    4: (0.0, 14.0), 5: (0.0, 14.0), 6: (0.0, 12.0)},
    "ylim": (0.0, 1.35),
    "show_exp_raw": True,
    "exp_raw_alpha": 0.15,               # raw trace stays in the background of the smoothed one
    "show_v_of_mean_E": True,
    "show_counts": True,                 # per-panel N text
    "annot_ke_n": (),                    # panels that print <E> sim / exp
    "titles": "$n = {n}$",
    "suptitle": None,
    "legend_panel": "figure",            # "figure" | panel index | None
    "legend_ncols": 3,                   # figure legend only
    "legend_loc": "best",
    "outer_labels_only": True,
    "xlabel": r"Speed $v$ / Å ps$^{-1}$",
    "ylabel": "Signal / arb. units",
}

# B2-4: B2-3 on the 2-D detector projection (sim: sqrt(vx^2 + vy^2));
# starts as a copy of FIG_B2_3, override keys here.
FIG_B2_4: dict = {
    **copy.deepcopy(FIG_B2_3),
    "stem": "b2_4_ked_shapes_2d",
    "representation": "2d",
}

# =============================================================================
# Imports (after the settings so SHOW_FIGURE picks the backend)
# =============================================================================
import matplotlib  # noqa: E402

if not SHOW_FIGURE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

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
    ke_band_ratio,
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
    ax.title.set_size(st["title_fontsize"])


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
    bare_frac: float
    num_scored: int
    num_ions: int
    w1_solvated: float
    sim_n1: float
    ref_n1: float
    sim_tail_mass: float
    ref_tail_mass: float


def size_data(read, abundance, tail_n: int = FIG_B2_1["tail_n"]) -> SizeData:
    """Scorer read path of ``_section_detected_size_distribution``; the tail
    mass is the solvated share at ``n >= tail_n`` (the
    ``mass_n15plus_solv`` column of ``tier2_g5_winners_table``); the bare
    fraction is the n = 0 share of the scored histogram."""
    score = score_histogram_vs_reference(read, abundance)
    sim_n, sim_f = solvated_renormalized(read.n_values, read.fraction, n_min=score.n_min)
    ref_n, ref_f = solvated_renormalized(abundance.n, abundance.ion_fraction, n_min=score.n_min)
    sim_n, sim_f = np.asarray(sim_n), np.asarray(sim_f, dtype=float)
    ref_n, ref_f = np.asarray(ref_n), np.asarray(ref_f, dtype=float)
    n_scored, frac_scored = np.asarray(read.n_values), np.asarray(read.fraction, dtype=float)
    return SizeData(
        n_scored=n_scored, frac_scored=frac_scored,
        n_sim_solv=sim_n, frac_sim_solv=sim_f, n_ref_solv=ref_n, frac_ref_solv=ref_f,
        n_mean=float(read.n_mean), suppressed_frac=float(read.suppressed_frac),
        bare_frac=float(frac_scored[n_scored == 0].sum()),
        num_scored=int(read.num_scored), num_ions=int(read.num_ions),
        w1_solvated=float(score.w1_solvated), sim_n1=float(score.sim_n1),
        ref_n1=float(score.ref_n1),
        sim_tail_mass=float(sim_f[sim_n >= tail_n].sum()),
        ref_tail_mass=float(ref_f[ref_n >= tail_n].sum()),
    )


@dataclass(frozen=True)
class KEData:
    """Per-n mean KE (all occupied bins, ``min_count = 1``), reference, scores.

    ``scored_mask`` marks the bins that pass the chi^2 min-count
    (``DETECTED_KE_MIN_BIN_COUNT``); the others are display-only thin bins.
    ``ratio`` is the per-bin ``ke_band_ratio`` over every bin (sim mean /
    raw reference mean -- the ingredient of midHot / deepKE), aligned with
    ``n``; ``ratio_err`` is its display error bar, sim SE and reference
    point error in quadrature (NaN where the sim SE is undefined).
    """
    n: np.ndarray
    mean_eV: np.ndarray
    se_eV: np.ndarray
    count: np.ndarray
    scored_mask: np.ndarray
    ratio: np.ndarray
    ratio_err: np.ndarray
    ref: object
    chi2_score: object
    midhot: object
    deepke: object


def ke_data(read, ked_ref) -> KEData:
    """Same scorer kwargs as ``summary_sections._detected_ke_scores``;
    midHot / deepKE at their committed conventions (geometric / arithmetic,
    ``min_count = 1``) -- the values in ``g13_decoupled_table.csv``."""
    ke = read.ke_by_n(min_count=1)
    keep = ke.n >= DETECTED_KE_N_MIN
    chi2 = score_ke_curve_vs_reference(
        read, ked_ref, n_min=DETECTED_KE_N_MIN, min_count=DETECTED_KE_MIN_BIN_COUNT,
        include_sim_se=DETECTED_KE_INCLUDE_SIM_SE, n1_anchor=DETECTED_KE_N1_ANCHOR,
    )
    n = np.asarray(ke.n)[keep]
    mean = np.asarray(ke.mean_eV, dtype=float)[keep]
    se = np.asarray(ke.se_eV, dtype=float)[keep]
    count = np.asarray(ke.count)[keep]
    per_bin = ke_band_ratio(read, ked_ref, n_lo=DETECTED_KE_N_MIN, n_hi=int(n.max()),
                            min_count=1)
    ratio_of = dict(zip(per_bin.n.tolist(), per_bin.ratio.tolist()))
    ratio = np.array([ratio_of.get(int(k), np.nan) for k in n])
    ref_idx = {int(v): i for i, v in enumerate(np.asarray(ked_ref.n))}
    rel_ref = np.array([
        ked_ref.point_err_eV[ref_idx[int(k)]] / ked_ref.mean_KE_eV[ref_idx[int(k)]]
        if int(k) in ref_idx else np.nan for k in n])
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio_err = ratio * np.hypot(se / mean, rel_ref)
    return KEData(
        n=n, mean_eV=mean, se_eV=se, count=count,
        scored_mask=count >= DETECTED_KE_MIN_BIN_COUNT, ratio=ratio, ratio_err=ratio_err,
        ref=ked_ref, chi2_score=chi2, midhot=midhot_ratio(read, ked_ref),
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
    """The scalars the figures quote, keyed by ``g13_decoupled_table.csv`` column."""
    ke1 = ke.mean_eV[ke.n == 1]
    return {
        "num_scored": size.num_scored,
        "supp": size.suppressed_frac,
        "nbar_det": size.n_mean,
        "n1_solv": size.sim_n1,
        "w1_solv": size.w1_solvated,
        "mass_n15plus_solv": size.sim_tail_mass,
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
            bbox=dict(boxstyle="square,pad=0.2", facecolor="white", edgecolor="none", alpha=0.85))


def draw_size_scored(ax, d: SizeData, st: dict) -> None:
    """Scored histogram, n = 0 (bare I+) included; bare fraction annotated."""
    ax.bar(d.n_scored, 100.0 * d.frac_scored, width=st["bar_width_scored"],
           edgecolor="black", linewidth=0.4, **st["sim"])
    if st["annot_scored"]:
        lo, hi = NBAR_LANDING
        flag = r"\in" if lo <= d.n_mean <= hi else r"\notin"
        _annot(ax, rf"$\bar n = {d.n_mean:.2f} {flag} [{lo:.2f}, {hi:.2f}]$" "\n"
                   f"Bare fraction {100.0 * d.bare_frac:.1f} %", st, st["annot_loc"])
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
    ax.bar(d.n_ref_solv - w / 2, 100.0 * d.frac_ref_solv, width=w, label="Experiment",
           **st["exp"])
    ax.bar(d.n_sim_solv + w / 2, 100.0 * d.frac_sim_solv, width=w, label="B2", **st["sim"])
    if st["annot_solvated"]:
        _annot(ax, rf"$n_1 = {d.sim_n1:.3f}$ (exp {d.ref_n1:.3f})" "\n"
                   rf"$W_1 = {d.w1_solvated:.3f}$ bins" "\n"
                   rf"$n \geq {st['tail_n']}$: {d.sim_tail_mass:.4f} (exp {d.ref_tail_mass:.4f})",
               st, st["annot_loc"])
    ax.set(xlim=st["n_xlim"], xlabel=st["xlabel"], ylabel=st["ylabel_solvated"])
    _legend(ax, st, st["legend_loc"])
    _finish_axis(ax, st, integer_x=True)


def draw_ke_mean(ax, d: KEData, st: dict, n_lo: int, n_hi: int, *,
                 show_legend: bool, show_chi2: bool) -> None:
    """Mean KE vs n over ``n_lo..n_hi`` -- the ``_section_detected_ked_mean_energy``
    panel: correlated bands as envelopes, per-point stat (+) sys bars, n = 1
    median anchor, scored sim bins (chi^2 min-count) with SE bars."""
    ref = d.ref
    n_all = np.asarray(ref.n)
    sel = (n_all >= max(n_lo, DETECTED_KE_N_MIN)) & (n_all <= n_hi)
    n_ref = n_all[sel]
    mean = np.asarray(ref.mean_KE_eV)[sel]
    for show, (b_lo, b_hi) in ((st["show_midhot_band"], MIDHOT_BAND),
                               (st["show_deepke_band"], DEEP_KE_BAND)):
        if show and b_lo <= n_hi and b_hi >= n_lo:
            ax.axvspan(max(b_lo, n_lo) - 0.5, min(b_hi, n_hi) + 0.5, zorder=0,
                       alpha=st["band_shade_alpha"], linewidth=0, **st["band"])
    if st["show_bands"]:
        for frac, label, alpha in zip(
                (ref.calib_syst_frac, ref.condition_syst_frac),
                ("Correlated bands (calibration, condition)", None), st["band_alphas"]):
            frac = np.asarray(frac)[sel]
            ax.fill_between(n_ref, mean * (1.0 - frac), mean * (1.0 + frac), linewidth=0,
                            alpha=alpha, label=label, **st["exp"])
    ax.errorbar(n_ref, mean, yerr=np.asarray(ref.point_err_eV)[sel], fmt="o", markersize=4,
                capsize=2, label=r"Experiment $\bar{E}$", **st["exp"])
    hit = n_ref == 1
    if st["show_median_anchor"] and np.any(hit):
        ax.plot(n_ref[hit], np.asarray(ref.median_KE_eV)[sel][hit], "D", markersize=6,
                markerfacecolor="none", markeredgecolor=st["anchor"]["color"],
                markeredgewidth=1.3, label="Experiment median ($n = 1$)")
    in_panel = (d.n >= n_lo) & (d.n <= n_hi)
    ok = in_panel & d.scored_mask
    if np.any(ok):
        ax.errorbar(d.n[ok], d.mean_eV[ok], yerr=d.se_eV[ok], fmt="s", linestyle="none",
                    markersize=4, capsize=2, label=r"Simulation (B2) $\bar{E}$",
                    **st["sim"])
    thin = in_panel & ~d.scored_mask
    if st["show_thin_bins"] and np.any(thin):
        ax.plot(d.n[thin], d.mean_eV[thin], "s", markersize=4, markerfacecolor="none",
                color=st["sim"]["color"],
                label=f"Simulation, < {DETECTED_KE_MIN_BIN_COUNT} ions")
    if show_chi2 and st["annot_chi2"]:
        _annot(ax, rf"$\chi^2_\mathrm{{med}} = {d.chi2_score.chi2_profiled:.1f}$; "
                   f"{d.chi2_score.n_points} bins", st, st["annot_loc"], ha="left", va="bottom")
    ax.set(yscale=st["yscale"], xlabel=st["xlabel"], ylabel=st["ylabel"])
    if show_legend:
        _legend(ax, st, st["legend_loc"])
    _finish_axis(ax, st, integer_x=True)


def draw_ke_ratio(ax, d: KEData, st: dict) -> None:
    """Per-bin sim / experiment mean-KE ratio (``ke_band_ratio``, raw
    reference mean) with the correlated bands as 1 +- frac envelopes and the
    committed midHot / deepKE band means as segments."""
    ref = d.ref
    n_ref = np.asarray(ref.n)
    sel = (n_ref >= DETECTED_KE_N_MIN) & (n_ref <= int(d.n.max()))
    if st["show_bands"]:
        for frac, alpha in zip((ref.calib_syst_frac, ref.condition_syst_frac),
                               st["band_alphas"]):
            frac = np.asarray(frac)[sel]
            ax.fill_between(n_ref[sel], 1.0 - frac, 1.0 + frac, linewidth=0, alpha=alpha,
                            **st["exp"])
    ax.axhline(1.0, linewidth=0.8, **st["exp"])
    ok = d.scored_mask & np.isfinite(d.ratio)
    ax.errorbar(d.n[ok], d.ratio[ok], yerr=d.ratio_err[ok], fmt="s", linestyle="none",
                markersize=4, capsize=2, **st["sim"])
    thin = ~d.scored_mask & np.isfinite(d.ratio)
    if st["show_thin_bins"] and np.any(thin):
        ax.plot(d.n[thin], d.ratio[thin], "s", markersize=4, markerfacecolor="none",
                color=st["sim"]["color"])
    if st["annot_band_means"]:
        for band, name, ls in ((d.midhot, "midHot", "--"), (d.deepke, "deepKE", ":")):
            if band.n_bins == 0:
                continue
            ax.hlines(band.value, band.n.min() - 0.4, band.n.max() + 0.4, linestyles=ls,
                      linewidth=1.4, label=f"{name} ($n$ = {band.n.min()}–{band.n.max()}): "
                                           f"{band.value:.3f}", **st["sim"])
        _legend(ax, st, st["legend_loc_ratio"])
    if st["ratio_ylim"] is not None:
        ax.set_ylim(st["ratio_ylim"])
    ax.set(xlim=(DETECTED_KE_N_MIN - 0.6, int(np.asarray(ref.n).max()) + 0.6),
           xlabel=st["xlabel"], ylabel=st["ylabel_ratio"])
    _finish_axis(ax, st, integer_x=True)


def draw_ke_residual(ax, d: KEData, st: dict) -> None:
    """Profiled residuals z (scored chi^2 bins) -- the ``_section_detected_ke_anatomy`` panel."""
    s = d.chi2_score
    hot, cold = st["z_colors"]
    ax.bar(s.n, s.z_profiled, width=0.8, edgecolor="black", linewidth=0.4,
           color=[hot if z > 0 else cold for z in s.z_profiled])
    ax.axhline(0.0, color="black", linewidth=0.8)
    for lvl in st["z_guides"]:
        ax.axhline(lvl, color="gray", linewidth=0.8, linestyle="--")
    if st["annot_z"]:
        _annot(ax, rf"$\chi^2_\mathrm{{med}} = {s.chi2_profiled:.1f}$, {s.n_points} bins" "\n"
                   rf"$a_\mathrm{{calib}} = {s.a_calib:.2f}$, $b_\mathrm{{cond}} = {s.b_condition:.2f}$",
               st, st["annot_loc"], ha="left", va="bottom")
    ax.set(xlim=(DETECTED_KE_N_MIN - 0.6, int(np.asarray(d.ref.n).max()) + 0.6),
           xlabel=st["xlabel"], ylabel=st["ylabel_z"])
    _finish_axis(ax, st, integer_x=True)


def draw_pv(ax, d: PvData, st: dict, *, show_legend: bool = False) -> None:
    """One KED-shape panel -- the ``_section_ihe_ked_curves`` panel in Å/ps."""
    if st["show_exp_raw"]:
        ax.plot(d.v_exp_Aps, d.exp_raw, ".-", markersize=2.0, linewidth=0.5,
                alpha=st["exp_raw_alpha"],
                label="Experiment (raw)", **st["exp"])
    ax.plot(d.v_exp_Aps, d.exp_smooth, linewidth=1.4, label="Experiment (smoothed)",
            **st["exp"])
    if d.sim_trace is not None:
        ax.plot(d.v_sim_Aps, d.sim_trace, "--", linewidth=1.4, label="Simulation (B2)",
                **st["sim"])
    else:
        ax.annotate("Simulation: no atoms in gate", (0.97, 0.9), xycoords="axes fraction",
                    ha="right", fontsize=st["annot_fontsize"], color=st["sim"]["color"])
    if st["show_v_of_mean_E"]:
        ax.axvline(d.v_mean_E_exp_Aps, linestyle=":", linewidth=1.0,
                   label=r"Experiment $v(\bar{E})$", **st["exp"])
        if d.v_mean_E_sim_Aps is not None:
            ax.axvline(d.v_mean_E_sim_Aps, linestyle=":", linewidth=1.0,
                       label=r"Simulation $v(\bar{E})$", **st["sim"])
    lines = []
    if st["show_counts"]:
        lines.append(f"$N = {d.num_sim}$")
    if d.n in st["annot_ke_n"] and d.mean_E_sim_eV is not None:
        lines.append(rf"$\bar{{E}}$ {d.mean_E_sim_eV:.3f} / exp {d.mean_E_exp_eV:.3f} eV")
    if lines:
        ax.text(0.97, 0.95, "\n".join(lines), transform=ax.transAxes, ha="right", va="top",
                fontsize=st["annot_fontsize"])
    xlim = st.get("v_xlim_by_n", {}).get(d.n, st["v_xlim_Aps"])
    ax.set(title=st["titles"].format(n=d.n), xlim=xlim, ylim=st["ylim"],
           xlabel=st["xlabel"], ylabel=st["ylabel"])
    if show_legend:
        _legend(ax, st, st["legend_loc"])
    _finish_axis(ax, st, integer_x=True)


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
    """``stacked``: (a) all n on a log axis over (b) ratio or z, shared n.
    ``split``: (a) n = 1..split-1 | (b) n = split..n_max on top, (c) below."""
    st = merged_style(FIG_B2_2 if style is None else style)
    split = st["split_n"]
    n_max = int(np.asarray(d.ref.n).max())
    lower = {"ratio": draw_ke_ratio, "z": draw_ke_residual}.get(st["lower_panel"])
    if lower is None:
        raise ValueError(f"lower_panel must be 'ratio' or 'z', got {st['lower_panel']!r}")
    if st["layout"] == "stacked":
        fig, (ax_e, ax_l) = plt.subplots(2, 1, figsize=st["figsize"], sharex=True,
                                         constrained_layout=True,
                                         gridspec_kw={"height_ratios": st["height_ratios"]})
        draw_ke_mean(ax_e, d, st | {"yscale": st["yscale_stacked"]}, DETECTED_KE_N_MIN, n_max,
                     show_legend=True, show_chi2=True)
        lower(ax_l, d, st)
        ax_e.set_xlabel("")
        if st["yscale_stacked"] == "log":
            ax_e.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
        _panel_letters((ax_e, ax_l), st)
        return fig
    if st["layout"] != "split":
        raise ValueError(f"unknown B2-2 layout {st['layout']!r}")
    fig, axd = plt.subplot_mosaic([["a", "b"], ["c", "c"]], figsize=st["figsize_split"],
                                  constrained_layout=True,
                                  gridspec_kw={"height_ratios": st["height_ratios_split"]})
    t_lo, t_hi, t_c = st["titles"]
    draw_ke_mean(axd["a"], d, st, DETECTED_KE_N_MIN, split - 1, show_legend=False, show_chi2=True)
    draw_ke_mean(axd["b"], d, st, split, n_max, show_legend=False, show_chi2=False)
    # one legend for (a) + (b) above the figure: wider than a half-width panel
    handles = {}
    for key in ("a", "b"):
        for h, lab in zip(*axd[key].get_legend_handles_labels()):
            handles.setdefault(lab, h)
    fig.legend(handles.values(), handles.keys(), loc="outside upper center", ncols=3,
               frameon=st["legend_frameon"], fontsize=st["legend_fontsize"])
    lower(axd["c"], d, st)
    axd["a"].set_title(t_lo.format(hi=split - 1), fontsize=st["title_fontsize"])
    axd["b"].set_title(t_hi.format(lo=split, n_max=n_max), fontsize=st["title_fontsize"])
    if t_c:
        axd["c"].set_title(t_c, fontsize=st["title_fontsize"])
    for key, lim in (("a", st["ylim_low"]), ("b", st["ylim_high"])):
        if lim is not None:
            axd[key].set_ylim(lim)
    _panel_letters((axd["a"], axd["b"], axd["c"]), st)
    return fig


def build_b2_3(panels: list[PvData], style: dict | None = None) -> plt.Figure:
    st = merged_style(FIG_B2_3 if style is None else style)
    ncols = st["ncols"]
    nrows = -(-len(panels) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=st["figsize"], sharex=st["sharex"],
                             sharey=st["sharey"], constrained_layout=True, squeeze=False)
    flat = axes.ravel()
    for k, (ax, d) in enumerate(zip(flat, panels)):
        draw_pv(ax, d, st, show_legend=(st["legend_panel"] == k))
        if st.get("outer_labels_only"):
            if k % ncols:
                ax.set_ylabel("")
            if k + ncols < len(panels):
                ax.set_xlabel("")
    for ax in flat[len(panels):]:
        ax.set_axis_off()
    if st["legend_panel"] == "figure":
        handles, labels = flat[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside upper center", ncols=st["legend_ncols"],
                   frameon=st["legend_frameon"], fontsize=st["legend_fontsize"])
    if st["suptitle"]:
        fig.suptitle(st["suptitle"])
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
    with style_ctx, plt.rc_context(FIG_RC):      # thesis_style font (Palatino)
        if "b2_1" in FIGURES:
            _save(build_b2_1(size), FIG_B2_1["stem"])
        if "b2_2" in FIGURES:
            _save(build_b2_2(ke), FIG_B2_2["stem"])
        if "b2_2_split" in FIGURES:
            _save(build_b2_2(ke, FIG_B2_2_SPLIT), FIG_B2_2_SPLIT["stem"])
        view = detected_ensemble_view(det)
        for key, fig_style in (("b2_3", FIG_B2_3), ("b2_4", FIG_B2_4)):
            if key in FIGURES:
                panels = [pv_data(view, IHE_KED_REFERENCE_DIR, ked_ref, n,
                                  fig_style["representation"]) for n in fig_style["n_values"]]
                _save(build_b2_3(panels, fig_style), fig_style["stem"])
    if SHOW_FIGURE:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
