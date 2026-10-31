"""Tier-2 discussion figures: B2, B3 and h405 vs experiment (detected ensemble).

The discussion-chapter twins of the B2 thesis figures
(``tier2_b2_thesis_figure``), keyed to ``THESIS_DISCUSSION.md`` §0
(identifiability), §11 (B2 vs the re-landed h405 under identical settings)
and §15 (final observable table). Every law is read from its pooled N = 6000
thesis container (``{B2,B3,h405}_N6000_pooled``, s1-s6,
``self_bound_gate_release = decoupled``), written to
``data/figures/thesis_tier2/`` -- L-1 and L-2 at half the thesis text width
(side-by-side subfigures, equal height), L-3 to L-5 at the full width --
once per ``FIGURE_VARIANTS`` entry of that figure: ``…_lq`` (B2 alone; the
n = 0 figures L-4 / L-5, written only so), ``…_h405`` (B2 vs the re-landed
h405, the discussion pair) and ``…_h405_b3`` (plus the effective high-tail
variant B3; both L-1 to L-3):

* **L-1 size distribution** -- solvated branch (n >= 1 renormalized, RQ8),
  experimental abundance as bars, each law as a marker line; (a) linear
  over (b) the n >= ``tail_lo`` tail on a log axis (§7, §15.3). n₁ and W₁
  per law annotated (the n >= 15 mass goes in the caption).
* **L-2 mean KE vs n** -- (a) all n on a log axis: experiment with the
  correlated bands, per-point bars and the n = 1 median, each law with SE
  bars (open markers: < chi^2 min-count ions); (b) per-bin ratio sim /
  experiment mean (``ke_band_ratio``) per law; midHot / deepKE per law
  annotated in (a). Legend labels shortened for the half width
  (``legend_short``). Laws are offset in n by ``slot x n_offset`` (fixed per
  law, so h405 sits at the same place in both variants; display only, the
  values are at integer n).
* **L-3 KED shapes** -- the B2-3 panels for n = 1-6 (3-D P(E) on the
  kinetic-energy axis, B2-3's ``x_axis``; peak-normalized), every law's
  trace over the experiment.
* **L-4 abundance with bare I+** -- (a) all detected ions, n = 0 included
  and not renormalised (experiment bars, n = 0 hatched; each law's line runs
  through n = 0, an open marker there), bare shares annotated; (b) the
  solvated branch, L-1 (a). n = 0 is shown, not scored: the experimental
  bare peak is mostly channel-external (RQ8), the simulated n = 0 is the
  boiled-off class under decoupled release (G13).
* **L-5 KED shapes n = 0, 1, 2** -- the L-3 traces side by side, mean
  (filled triangle) and median (open diamond) of the experiment and each law
  as marker rows above the traces (value box for <= 2 series). The n = 0
  panel carries the symmetric point-charge benchmark E_C / 2 (Q = 2)
  derived from the pool's ``cfg.json`` (``coulomb_benchmarks_eV``), guarded
  against the documented 2.70 eV.
  Simulated median = ``np.median`` of the scored-ion KE at n (the committed
  ``KE{n}_med`` recipe, ``scored_ke_median``).

Every recipe is the B2 script's data layer (``size_data``, ``ke_data``,
``pv_data`` -- the ``tier2_confirmation`` scorers); the draw layer is reused
for the experiment and the first law, the further laws are overlaid in their
own colour / marker / dash. Nothing is re-derived.

Anchor guard: per law, the quoted scalars (``anchor_values``) are checked
against that law's decoupled row of
``data/runs/h2b_forward_model/g13_decoupled_table.csv`` to 4 dp (row label
= law, container = pool name), plus KE1 / KE2 mean and median there; the
bare share, the n = 0 mean / median KE and the KE1 / KE2 medians against
``thesis_observables_extra.csv``; a mismatch raises.

Edit the ``USER SETTINGS`` block and run::

    python scripts/post_processing/tier2_law_comparison_figure.py
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT_DIR = Path(__file__).resolve().parent
for _p in (PROJECT_ROOT, _SCRIPT_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from scripts.post_processing.thesis_style import (  # noqa: E402  (shared thesis style)
    B3_COLOR, CONTRAST_COLOR, EXP_COLOR, FIG_RC, LIT_COLOR, SIM_COLOR, TEXTWIDTH_IN,
    THESIS_NAMES,
)

# =============================================================================
# USER SETTINGS
# =============================================================================
RUNS_DIR: Path = PROJECT_ROOT / "data" / "runs"
# law -> (pool, legend name, colour, sim marker, P(v) dash, n-offset slot);
# order = draw order. The slot (x n_offset) is fixed per law so a law sits at
# the same place in every variant; B2 (drawn first, by the B2 draw layer) is 0.
LAWS: dict[str, tuple[str, str, str, str, str, int]] = {
    # legend names = thesis_style.THESIS_NAMES (B2 -> LQ, h405 -> CC5.4)
    "B2": ("B2_N6000_pooled", THESIS_NAMES["B2"], SIM_COLOR, "s", "--", 0),
    "B3": ("B3_N6000_pooled", THESIS_NAMES["B3"], B3_COLOR, "^", "-.", -1),
    "h405": ("h405_N6000_pooled", THESIS_NAMES["h405"], CONTRAST_COLOR, "v", ":", 1),
}
# file-stem suffix -> laws drawn
VARIANTS: dict[str, tuple[str, ...]] = {
    "lq": ("B2",),                        # production law alone
    "h405": ("B2", "h405"),               # discussion: B2 vs the contrast law
    "h405_b3": ("B2", "B3", "h405"),      # with the effective high-tail variant
}
# figure -> the variants it is written for
FIGURE_VARIANTS: dict[str, tuple[str, ...]] = {
    "l_1": ("h405", "h405_b3"),
    "l_2": ("h405", "h405_b3"),
    "l_3": ("h405", "h405_b3"),
    "l_4": ("lq",),                       # n = 0 discussion: production law only
    "l_5": ("lq",),
}
WINNERS_CSV: Path = RUNS_DIR / "h2b_forward_model" / "g13_decoupled_table.csv"
# bare share and the n = 0 / median KE columns (not in the g13 table)
EXTRA_CSV: Path = RUNS_DIR / "h2b_forward_model" / "thesis_observables_extra.csv"
ANCHOR_CHECK: bool = True
IHE_KED_REFERENCE_DIR: Path = PROJECT_ROOT / "data" / "reference" / "ihe_ked"
ABUNDANCE_REF_PATH: Path = PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"

OUT_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_tier2"
FIGURES: tuple[str, ...] = ("l_1", "l_2", "l_3", "l_4", "l_5")
FORMATS: tuple[str, ...] = ("png",)
DPI: int = 200
SHOW_FIGURE: bool = False

# L-1 and L-2 are drawn at half the text width, to sit side by side in the
# thesis (two subfigures of 0.5 \textwidth); same height, so they align.
HALF_FIGSIZE: tuple[float, float] = (TEXTWIDTH_IN / 2, 3.7)

FIG_L_1: dict = {
    "stem": "laws_1_size_distribution",
    "figsize": HALF_FIGSIZE,
    "height_ratios": (1.35, 1.0),        # (a) linear over (b) log tail
    "n_xlim": (0.3, 20.7),
    "tail_lo": 10,                       # (b) starts here
    "tail_ylim": (5e-3, 5.0),            # % ; zero bins are not drawn on the log axis
    "bar_width": 0.8,
    "exp_bar_alpha": 0.45,
    "markersize": 3.0,
    "linewidth": 1.0,
    "annot_loc": (0.97, 0.95),           # n₁ / W₁ box in (a); the n >= 15 mass goes in the caption
    "legend_max_cols": 3,
}

FIG_L_2: dict = {
    "stem": "laws_2_mean_ke",
    "figsize": HALF_FIGSIZE,
    "height_ratios": (1.2, 1.0),
    "n_offset": 0.22,                    # display offset between laws in n
    "markersize": 3.0,
    "yscale": "log",
    "ylim_e": (1.5e-3, 2.0),             # eV; floor leaves the lower-left corner to the box
    "annot_band_values": True,           # midHot / deepKE per law, empty corner of (a)
    "annot_loc_band_values": (0.03, 0.04),
    "legend_ncols": 2,
    "panel_letter_start": "c",           # continues L-1's (a) / (b) beside it
    # short legend labels for the half width (full wording in the caption)
    "legend_short": {
        "Correlated bands (calibration, condition)": "Exp. correlated bands",
        r"Experiment $\bar{E}$": "Experiment",
        "Experiment median ($n = 1$)": "Exp. median ($n = 1$)",
    },
}

FIG_L_3: dict = {
    "stem": "laws_3_ked_shapes",
    "representation": "3d",
    "n_values": (1, 2, 3, 4, 5, 6),
    "linewidth": 1.3,
    "show_mean_E_marker": False,         # four E-bar marks per panel do not read; caption
    "show_counts": False,
    "legend_ncols": 3,
}

# L-4 and L-5 at the full text width, side-by-side panels.
FIG_L_4: dict = {
    "stem": "laws_4_size_with_bare",
    "figsize": (TEXTWIDTH_IN, 3.0),
    "n_xlim_all": (-0.7, 20.7),          # (a) includes n = 0
    "bare_hatch": "////",                # n = 0: shown, not scored
    "separator": dict(color="0.6", linewidth=0.6),
    "annot_loc_bare": (0.97, 0.95),
}

FIG_L_5: dict = {
    "stem": "laws_5_ked_n012",
    "n_values": (0, 1, 2),
    "figsize": (TEXTWIDTH_IN, 3.0),
    # n = 0: the exp. 3-D columns are NaN below 0.4 eV (Abel-centre spike)
    # and end at the 6.1 eV edge of the exported speed window
    "e_xlim_by_n": {0: (0.0, 6.0), 1: (0.0, 3.5), 2: (0.0, 2.0)},
    "linewidth": 1.3,
    "show_mean_E_marker": False,         # mean / median drawn as marker rows instead
    "show_counts": False,
    "mark_y0": 1.27,                     # experiment's row; laws step down by mark_dy
    "mark_dy": 0.07,
    "mark_size": 4.0,
    # mean / median value box (x, y, ha, va in axes fraction), under the
    # curves; drawn only up to moment_values_max_rows series (experiment +
    # LQ) -- with more laws it covers the traces and the values go in the caption
    "moment_values_max_rows": 2,
    "annot_loc_by_n": {0: (0.97, 0.03, "right", "bottom"),
                       1: (0.97, 0.82, "right", "top"),
                       2: (0.97, 0.82, "right", "top")},
    "show_symmetric_line": True,         # n = 0 panel: point-charge E_C / 2 (Q = 2)
    "benchmark": dict(color=LIT_COLOR, linewidth=1.0),
    "legend_ncols": 4,
}
# Documented benchmark value [eV] (THESIS_DISCUSSION §9); the derived line
# must match it to 2 dp.
BENCHMARK_DOC_EV: dict[str, float] = {"symmetric": 2.70}

# =============================================================================
# Imports (after the settings so SHOW_FIGURE picks the backend)
# =============================================================================
import matplotlib  # noqa: E402

if not SHOW_FIGURE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

from i2_helium_md.physics.constants import E_CHARGE, coulomb_energy  # noqa: E402
from i2_helium_md.postprocess import (  # noqa: E402
    detected_ensemble_view,
    load_he_abundance_reference,
    load_ihe_ked_reference,
)
from i2_helium_md.postprocess.tier2_confirmation import read_confirmation_detection  # noqa: E402
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from scripts.post_processing.tier2_b2_thesis_figure import (  # noqa: E402
    FIG_B2_1,
    FIG_B2_2,
    FIG_B2_3,
    KEData,
    PvData,
    SizeData,
    _annot,
    _finish_axis,
    _panel_letters,
    anchor_values,
    check_anchor,
    draw_ke_mean,
    draw_ke_ratio,
    draw_pv,
    ke_data,
    ked_traces,
    merged_style,
    pv_data,
    size_data,
)
from scripts.post_processing.tier2_thesis_observable_table import scored_ke_median  # noqa: E402
from summary_sections import DETECTED_KE_N_MIN  # noqa: E402


# =============================================================================
# Data layer
# =============================================================================
@dataclass(frozen=True)
class LawData:
    """One law: legend name, style, and the B2-script data products."""
    key: str
    name: str
    color: str
    marker: str
    dash: str
    slot: int
    size: SizeData
    ke: KEData
    view: object
    read: object


def ke_moments(read, n: int) -> tuple[float, float]:
    """(mean, median) detected KE [eV] of the scored ions at shell size ``n``:
    the scorer's ``ke_by_n`` mean and the committed ``KE{n}_med`` recipe."""
    ke = read.ke_by_n(min_count=1)
    hit = np.flatnonzero(ke.n == n)
    if hit.size != 1:
        raise ValueError(f"no scored ions at n = {n}")
    return float(ke.mean_eV[hit[0]]), scored_ke_median(read, n)


def exp_ke_moments(ked_ref, n: int) -> tuple[float, float]:
    """(mean, median) KE [eV] of the experimental fragment n (``IHe_KED_reference.csv``)."""
    hit = np.flatnonzero(np.asarray(ked_ref.n) == n)
    if hit.size != 1:
        raise ValueError(f"ihe_ked reference has no unique row for n = {n}")
    return float(ked_ref.mean_KE_eV[hit[0]]), float(ked_ref.median_KE_eV[hit[0]])


def moment_anchor_values(law: LawData) -> tuple[dict[str, float], dict[str, float]]:
    """The L-4 / L-5 scalars, keyed by (g13 column, extra-table column)."""
    (m1, d1), (m2, d2), (m0, d0) = (ke_moments(law.read, n) for n in (1, 2, 0))
    g13 = {"KE1_mean": m1, "KE1_med": d1, "KE2_mean": m2, "KE2_med": d2}
    extra = {"bare_share": law.size.bare_frac, "KE0_mean_eV": m0, "KE0_median_eV": d0,
             "KE1_median_eV": d1, "KE2_median_eV": d2}
    return g13, extra


def coulomb_benchmarks_eV(cfg_path: Path) -> dict[str, float]:
    """Bare-ion KE benchmark [eV] of the production Coulomb budget.

    Pair budget E = E_coulomb_scale * e^2 / (4 pi eps0 R0) at the run's
    ``R0_GS_angstrom`` (Q = 2, point charges); ``symmetric`` = E / 2, the
    per-fragment share at equal masses.
    """
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    pair_eV = (float(cfg["E_coulomb_scale"]) * coulomb_energy(float(cfg["R0_GS_angstrom"]))
               / abs(E_CHARGE))
    return {"symmetric": pair_eV / 2.0}


def check_benchmarks(values: dict[str, float]) -> None:
    """Raise unless every derived benchmark matches ``BENCHMARK_DOC_EV`` to 2 dp."""
    bad = [f"{k}: derived {v:.4f} vs documented {BENCHMARK_DOC_EV[k]:.2f}"
           for k, v in values.items() if round(v, 2) != BENCHMARK_DOC_EV[k]]
    if bad:
        raise ValueError("benchmark mismatch:\n  " + "\n  ".join(bad))


def load_law(key: str, abundance, ked_ref) -> LawData:
    """Read one pooled container with the scorer's reader; anchor-check it."""
    pool, name, color, marker, dash, slot = LAWS[key]
    run_dir = RUNS_DIR / pool
    det = load_detection_result(run_dir / "detection.npz")
    read = read_confirmation_detection(det, label=pool)
    size, ke = size_data(read, abundance), ke_data(read, ked_ref)
    law = LawData(key=key, name=name, color=color, marker=marker, dash=dash, slot=slot,
                  size=size, ke=ke, view=detected_ensemble_view(det), read=read)
    if ANCHOR_CHECK:
        g13, extra = moment_anchor_values(law)
        check_anchor(anchor_values(size, ke) | g13, WINNERS_CSV, key, pool)
        check_anchor(extra, EXTRA_CSV, key)
        print(f"[laws] anchor OK vs {WINNERS_CSV.name} + {EXTRA_CSV.name}:{key}")
    return law


def _law_style(base: dict, law: LawData) -> dict:
    """``base`` with the simulation colour and legend name of ``law``."""
    return base | {"sim": dict(color=law.color), "sim_label": law.name}


# =============================================================================
# Draw layer
# =============================================================================
def draw_size_laws(ax, laws: list[LawData], st: dict, *, log_tail: bool) -> None:
    """Solvated f(n): experiment bars, one marker line per law (in %)."""
    ref = laws[0].size
    ax.bar(ref.n_ref_solv, 100.0 * ref.frac_ref_solv, width=st["bar_width"],
           alpha=st["exp_bar_alpha"], linewidth=0, label="Experiment", **st["exp"])
    for law in laws:
        n, f = law.size.n_sim_solv, 100.0 * law.size.frac_sim_solv
        if log_tail:
            keep = f > 0
            n, f = n[keep], f[keep]
        ax.plot(n, f, marker=law.marker, markersize=st["markersize"],
                linewidth=st["linewidth"], color=law.color, label=law.name)
    if log_tail:
        ax.set(yscale="log", xlim=(st["tail_lo"] - 0.6, st["n_xlim"][1]),
               ylim=st["tail_ylim"])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    else:
        ax.set(xlim=st["n_xlim"])
        lines = [rf"Exp.: $n_1$ {ref.ref_n1:.3f}"]
        lines += [rf"{law.name}: $n_1$ {law.size.sim_n1:.3f}, $W_1$ {law.size.w1_solvated:.3f}"
                  for law in laws]
        _annot(ax, "\n".join(lines), st, st["annot_loc"])
    ax.set(xlabel=st["xlabel"], ylabel=st["ylabel_solvated"])
    _finish_axis(ax, st, integer_x=True)


def draw_ke_overlay(ax, law: LawData, dx: float, ms: float) -> None:
    """Mean KE per n of a further law, offset by ``dx`` in n (scored + thin bins)."""
    d = law.ke
    ok, thin = d.scored_mask, ~d.scored_mask
    ax.errorbar(d.n[ok] + dx, d.mean_eV[ok], yerr=d.se_eV[ok], fmt=law.marker,
                linestyle="none", markersize=ms, capsize=1.5, color=law.color,
                label=rf"Simulation ({law.name}) $\bar{{E}}$")
    if np.any(thin):
        ax.plot(d.n[thin] + dx, d.mean_eV[thin], law.marker, markersize=ms,
                markerfacecolor="none", color=law.color)


def draw_ratio_overlay(ax, law: LawData, dx: float, ms: float) -> None:
    """Per-bin KE ratio of a further law, offset by ``dx`` in n."""
    d = law.ke
    ok = d.scored_mask & np.isfinite(d.ratio)
    thin = ~d.scored_mask & np.isfinite(d.ratio)
    ax.errorbar(d.n[ok] + dx, d.ratio[ok], yerr=d.ratio_err[ok], fmt=law.marker,
                linestyle="none", markersize=ms, capsize=1.5, color=law.color)
    if np.any(thin):
        ax.plot(d.n[thin] + dx, d.ratio[thin], law.marker, markersize=ms,
                markerfacecolor="none", color=law.color)


def draw_size_with_bare(ax, laws: list[LawData], abundance, st: dict) -> None:
    """All detected ions, n = 0 included (not renormalised), in %: experiment
    bars (n = 0 hatched), each law as a marker line over n >= 1 and an open
    marker at n = 0; bare shares annotated."""
    n_exp = np.asarray(abundance.n)
    f_exp = 100.0 * np.asarray(abundance.ion_fraction, dtype=float)
    bare = n_exp == 0
    ax.bar(n_exp[~bare], f_exp[~bare], width=st["bar_width"], alpha=st["exp_bar_alpha"],
           linewidth=0, **st["exp"])
    exp_rgba = to_rgba(st["exp"]["color"])
    ax.bar(n_exp[bare], f_exp[bare], width=st["bar_width"], linewidth=0,
           facecolor=(*exp_rgba[:3], st["exp_bar_alpha"]), edgecolor=exp_rgba,
           hatch=st["bare_hatch"])
    ax.axvline(0.5, zorder=0, **st["separator"])
    for law in laws:
        n, f = law.size.n_scored, 100.0 * law.size.frac_scored
        ax.plot(n[n >= 1], f[n >= 1], marker=law.marker, markersize=st["markersize"],
                linewidth=st["linewidth"], color=law.color)
        ax.plot(n[n <= 1], f[n <= 1], linewidth=st["linewidth"], color=law.color)
        ax.plot(n[n == 0], f[n == 0], law.marker, markersize=st["markersize"] + 1.0,
                markerfacecolor="white", color=law.color)
    lines = [r"Bare I$^+$ ($n = 0$, not scored)", f"Exp.: {f_exp[bare].sum():.1f} %"]
    lines += [f"{law.name}: {100.0 * law.size.bare_frac:.1f} %" for law in laws]
    _annot(ax, "\n".join(lines), st, st["annot_loc_bare"])
    ax.set(xlim=st["n_xlim_all"], xlabel=st["xlabel"], ylabel=st["ylabel_scored"])
    _finish_axis(ax, st, integer_x=True)


def draw_ked_moments(ax, laws: list[LawData], ked_ref, n: int, st: dict) -> None:
    """Mean (filled triangle) and median (open diamond) KE [eV] of the
    experiment and each law as marker rows above the peak-normalised traces,
    plus a mean / median value box when there are at most
    ``moment_values_max_rows`` series."""
    rows = [("Exp.", st["exp"]["color"], exp_ke_moments(ked_ref, n))]
    rows += [(law.name, law.color, ke_moments(law.read, n)) for law in laws]
    ms = st["mark_size"]
    for k, (_, color, (mean, median)) in enumerate(rows):
        y = st["mark_y0"] - k * st["mark_dy"]
        ax.plot([mean], [y], "v", markersize=ms, color=color, clip_on=False)
        ax.plot([median], [y], "D", markersize=ms * 0.85, markerfacecolor="none",
                color=color, clip_on=False)
    if len(rows) <= st["moment_values_max_rows"]:
        text = "Mean / median (eV)\n" + "\n".join(
            f"{name}: {mean:.3f} / {median:.3f}" for name, _, (mean, median) in rows)
        x, y, ha, va = st["annot_loc_by_n"][n]
        _annot(ax, text, st, (x, y), ha=ha, va=va)


def draw_benchmarks(ax, values: dict[str, float], st: dict) -> None:
    """Bare-ion Coulomb benchmark (``coulomb_benchmarks_eV``) as a vertical line."""
    if st["show_symmetric_line"]:
        ax.axvline(values["symmetric"], linestyle="-.", **st["benchmark"],
                   label=rf"Symmetric $E_C/2$ ({values['symmetric']:.2f} eV)")


def _shrink_markers(ax, ms: float) -> None:
    """Set every marker in ``ax`` to ``ms`` (the reused B2 draw layer uses 4 / 6)."""
    for line in ax.get_lines():
        if line.get_marker() not in (None, "None", ""):
            line.set_markersize(ms * 1.5 if line.get_marker() == "D" else ms)


def _short_label(label: str, st: dict) -> str:
    """Half-width legend wording: ``Simulation (X) Ē`` -> ``X``, fixed map otherwise."""
    if label.startswith("Simulation (") and label.endswith(r") $\bar{E}$"):
        return label[len("Simulation ("):-len(r") $\bar{E}$")]
    if label.startswith("Simulation, <"):
        return "Open: " + label[len("Simulation, "):]
    return st["legend_short"].get(label, label)


# =============================================================================
# Figure builders
# =============================================================================
def build_l_1(laws: list[LawData]) -> plt.Figure:
    """(a) solvated f(n), linear, over (b) the n >= tail_lo tail, log axis."""
    st = merged_style(FIG_B2_1 | FIG_L_1)
    fig, axes = plt.subplots(2, 1, figsize=st["figsize"], constrained_layout=True,
                             gridspec_kw={"height_ratios": st["height_ratios"]})
    draw_size_laws(axes[0], laws, st, log_tail=False)
    draw_size_laws(axes[1], laws, st, log_tail=True)
    axes[0].set_xlabel("")
    axes[1].set_ylabel(r"Fraction / %")
    _size_legend(fig, axes[0], st)
    _panel_letters(axes, st)
    return fig


def _size_legend(fig, ax, st: dict) -> None:
    """Figure legend from ``ax``: experiment first, laws in draw order."""
    handles, labels = ax.get_legend_handles_labels()
    first = labels.index("Experiment")
    handles = [handles[first]] + handles[:first] + handles[first + 1:]
    labels = [labels[first]] + labels[:first] + labels[first + 1:]
    ncols = len(labels) if len(labels) <= st["legend_max_cols"] else 2
    fig.legend(handles, labels, loc="outside upper center", ncols=ncols,
               frameon=st["legend_frameon"], fontsize=st["legend_fontsize"])


def build_l_2(laws: list[LawData]) -> plt.Figure:
    """(a) mean KE vs n, log axis; (b) per-bin ratio; laws offset in n."""
    st = merged_style(FIG_B2_2 | FIG_L_2)
    first, rest = laws[0], laws[1:]
    if first.slot != 0:
        raise ValueError(f"first law {first.key!r} must sit at slot 0 (drawn without offset)")
    n_max = int(np.asarray(first.ke.ref.n).max())
    fig, (ax_e, ax_r) = plt.subplots(2, 1, figsize=st["figsize"], sharex=True,
                                     constrained_layout=True,
                                     gridspec_kw={"height_ratios": st["height_ratios"]})
    st_first = _law_style(st, first)
    # the first law at slot 0 via the reused draw layer (no offset there)
    draw_ke_mean(ax_e, first.ke, st_first | {"yscale": st["yscale"]}, DETECTED_KE_N_MIN,
                 n_max, show_legend=False, show_chi2=False)
    draw_ke_ratio(ax_r, first.ke, st_first | {"annot_band_means": False})
    ms = st["markersize"]
    for ax in (ax_e, ax_r):
        _shrink_markers(ax, ms)
    for law in rest:
        dx = law.slot * st["n_offset"]
        draw_ke_overlay(ax_e, law, dx, ms)
        draw_ratio_overlay(ax_r, law, dx, ms)
    ax_e.set_xlabel("")
    if st["ylim_e"] is not None:
        ax_e.set_ylim(st["ylim_e"])
    if st["yscale"] == "log":
        ax_e.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    handles, labels = ax_e.get_legend_handles_labels()
    fig.legend(handles, [_short_label(lab, st) for lab in labels], loc="outside upper center",
               ncols=st["legend_ncols"], frameon=st["legend_frameon"],
               fontsize=st["legend_fontsize"])
    if st["annot_band_values"]:
        def fmt(band) -> str:
            return "—" if band.n_bins == 0 else f"{band.value:.3f}"
        lines = ["midHot / deepKE"]          # n bands (n 2-8 / 10-17): caption
        lines += [f"{law.name}: {fmt(law.ke.midhot)} / {fmt(law.ke.deepke)}" for law in laws]
        _annot(ax_e, "\n".join(lines), st, st["annot_loc_band_values"], ha="left", va="bottom")
    _panel_letters((ax_e, ax_r), st)
    return fig


def build_l_3(laws: list[LawData], panels: dict[str, list[PvData]]) -> plt.Figure:
    """B2-3 panels with every law's trace; reuses ``draw_pv`` for the first law."""
    st = merged_style(FIG_B2_3 | FIG_L_3)
    ncols = st["ncols"]
    nrows = -(-len(st["n_values"]) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=st["figsize"], sharex=st["sharex"],
                             sharey=st["sharey"], constrained_layout=True, squeeze=False)
    flat = axes.ravel()
    first = laws[0]
    for k, ax in enumerate(flat[:len(st["n_values"])]):
        draw_pv(ax, panels[first.key][k], _law_style(st, first))
        for law in laws[1:]:
            t = ked_traces(panels[law.key][k], st["x_axis"], st["e_min_eV"])
            if t.sim is not None:
                ax.plot(t.x_sim, t.sim, law.dash, linewidth=st["linewidth"],
                        color=law.color, label=f"Simulation ({law.name})")
        if st["outer_labels_only"]:
            if k % ncols:
                ax.set_ylabel("")
            if k + ncols < len(st["n_values"]):
                ax.set_xlabel("")
    for ax in flat[len(st["n_values"]):]:
        ax.set_axis_off()
    handles, labels = flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncols=st["legend_ncols"],
               frameon=st["legend_frameon"], fontsize=st["legend_fontsize"])
    _panel_letters(flat[:len(st["n_values"])], st)
    return fig


def build_l_4(laws: list[LawData], abundance) -> plt.Figure:
    """(a) all detected ions incl. n = 0 | (b) the solvated branch (L-1 (a))."""
    st = merged_style(FIG_B2_1 | FIG_L_1 | FIG_L_4)
    fig, axes = plt.subplots(1, 2, figsize=st["figsize"], constrained_layout=True)
    draw_size_with_bare(axes[0], laws, abundance, st)
    draw_size_laws(axes[1], laws, st, log_tail=False)
    _size_legend(fig, axes[1], st)
    _panel_letters(axes, st)
    return fig


def build_l_5(laws: list[LawData], panels: dict[str, list[PvData]], ked_ref,
              benchmarks: dict[str, float]) -> plt.Figure:
    """KED shapes n = 0, 1, 2 side by side (L-3 traces) with mean / median
    marks; the bare-ion panel carries the Coulomb benchmark."""
    st = merged_style(FIG_B2_3 | FIG_L_3 | FIG_L_5)
    fig, axes = plt.subplots(1, len(st["n_values"]), figsize=st["figsize"], sharey=True,
                             constrained_layout=True)
    first = laws[0]
    for k, (ax, n) in enumerate(zip(axes, st["n_values"])):
        draw_pv(ax, panels[first.key][k], _law_style(st, first))
        for law in laws[1:]:
            t = ked_traces(panels[law.key][k], st["x_axis"], st["e_min_eV"])
            if t.sim is not None:
                ax.plot(t.x_sim, t.sim, law.dash, linewidth=st["linewidth"],
                        color=law.color, label=f"Simulation ({law.name})")
        if n == 0:
            draw_benchmarks(ax, benchmarks, st)
        draw_ked_moments(ax, laws, ked_ref, n, st)
        if k:
            ax.set_ylabel("")
    handles, labels = axes[0].get_legend_handles_labels()
    handles += [Line2D([], [], marker="v", linestyle="none", color="0.3"),
                Line2D([], [], marker="D", linestyle="none", markerfacecolor="none",
                       color="0.3")]
    labels += ["Mean", "Median"]
    fig.legend(handles, labels, loc="outside upper center", ncols=st["legend_ncols"],
               frameon=st["legend_frameon"], fontsize=st["legend_fontsize"])
    _panel_letters(axes, st)
    return fig


# =============================================================================
# Main
# =============================================================================
def _save(fig: plt.Figure, stem: str) -> None:
    for ext in FORMATS:
        path = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(path, dpi=DPI)
        print(f"[laws] wrote {path}")


def main() -> int:
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_DIR / "IHe_KED_reference.csv")
    abundance = load_he_abundance_reference(ABUNDANCE_REF_PATH)
    needed = {key for fig_key in FIGURES for suffix in FIGURE_VARIANTS[fig_key]
              for key in VARIANTS[suffix]}
    by_key = {key: load_law(key, abundance, ked_ref) for key in LAWS if key in needed}
    rep = FIG_L_3["representation"]
    panels, panels_5 = {}, {}
    if "l_3" in FIGURES:
        panels = {key: [pv_data(law.view, IHE_KED_REFERENCE_DIR, ked_ref, n, rep)
                        for n in FIG_L_3["n_values"]] for key, law in by_key.items()}
    benchmarks: dict[str, float] = {}
    if "l_5" in FIGURES:
        panels_5 = {key: [pv_data(law.view, IHE_KED_REFERENCE_DIR, ked_ref, n, rep)
                          for n in FIG_L_5["n_values"]] for key, law in by_key.items()}
        benchmarks = coulomb_benchmarks_eV(RUNS_DIR / LAWS["B2"][0] / "cfg.json")
        check_benchmarks(benchmarks)
        print(f"[laws] benchmarks OK: " + ", ".join(f"{k} {v:.4f} eV"
                                                     for k, v in benchmarks.items()))

    builders = {
        "l_1": (FIG_L_1, lambda laws: build_l_1(laws)),
        "l_2": (FIG_L_2, lambda laws: build_l_2(laws)),
        "l_3": (FIG_L_3, lambda laws: build_l_3(laws, panels)),
        "l_4": (FIG_L_4, lambda laws: build_l_4(laws, abundance)),
        "l_5": (FIG_L_5, lambda laws: build_l_5(laws, panels_5, ked_ref, benchmarks)),
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(FIG_RC):                  # thesis_style font (Palatino)
        for fig_key in FIGURES:
            fig_style, build = builders[fig_key]
            for suffix in FIGURE_VARIANTS[fig_key]:
                fig = build([by_key[key] for key in VARIANTS[suffix]])
                _save(fig, f"{fig_style['stem']}_{suffix}")
                if not SHOW_FIGURE:
                    plt.close(fig)
    if SHOW_FIGURE:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
