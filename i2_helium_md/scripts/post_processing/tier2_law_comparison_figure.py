"""Tier-2 discussion figures: B2, B3 and h405 vs experiment (detected ensemble).

The discussion-chapter twins of the B2 thesis figures
(``tier2_b2_thesis_figure``), keyed to ``THESIS_DISCUSSION.md`` §0
(identifiability), §11 (B2 vs the re-landed h405 under identical settings)
and §15 (final observable table). Every law is read from its pooled N = 6000
thesis container (``{B2,B3,h405}_N6000_pooled``, s1-s6,
``self_bound_gate_release = decoupled``), written to
``data/figures/thesis_tier2/`` -- L-1 and L-2 at half the thesis text width
(side-by-side subfigures, equal height), L-3 at the full width -- once per
``VARIANTS`` entry: ``…_h405``
(B2 vs the re-landed h405, the discussion pair) and ``…_h405_b3`` (plus the
effective high-tail variant B3):

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

Every recipe is the B2 script's data layer (``size_data``, ``ke_data``,
``pv_data`` -- the ``tier2_confirmation`` scorers); the draw layer is reused
for the experiment and the first law, the further laws are overlaid in their
own colour / marker / dash. Nothing is re-derived.

Anchor guard: per law, the quoted scalars (``anchor_values``) are checked
against that law's decoupled row of
``data/runs/h2b_forward_model/g13_decoupled_table.csv`` to 4 dp (row label
= law, container = pool name); a mismatch raises.

Edit the ``USER SETTINGS`` block and run::

    python scripts/post_processing/tier2_law_comparison_figure.py
"""

from __future__ import annotations

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
    B3_COLOR, CONTRAST_COLOR, FIG_RC, SIM_COLOR, TEXTWIDTH_IN, THESIS_NAMES,
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
# file-stem suffix -> laws drawn; every figure is written once per variant
VARIANTS: dict[str, tuple[str, ...]] = {
    "h405": ("B2", "h405"),               # discussion: B2 vs the contrast law
    "h405_b3": ("B2", "B3", "h405"),      # with the effective high-tail variant
}
WINNERS_CSV: Path = RUNS_DIR / "h2b_forward_model" / "g13_decoupled_table.csv"
ANCHOR_CHECK: bool = True
IHE_KED_REFERENCE_DIR: Path = PROJECT_ROOT / "data" / "reference" / "ihe_ked"
ABUNDANCE_REF_PATH: Path = PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"

OUT_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_tier2"
FIGURES: tuple[str, ...] = ("l_1", "l_2", "l_3")
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

# =============================================================================
# Imports (after the settings so SHOW_FIGURE picks the backend)
# =============================================================================
import matplotlib  # noqa: E402

if not SHOW_FIGURE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

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


def load_law(key: str, abundance, ked_ref) -> LawData:
    """Read one pooled container with the scorer's reader; anchor-check it."""
    pool, name, color, marker, dash, slot = LAWS[key]
    run_dir = RUNS_DIR / pool
    det = load_detection_result(run_dir / "detection.npz")
    read = read_confirmation_detection(det, label=pool)
    size, ke = size_data(read, abundance), ke_data(read, ked_ref)
    if ANCHOR_CHECK:
        check_anchor(anchor_values(size, ke), WINNERS_CSV, key, pool)
        print(f"[laws] anchor OK vs {WINNERS_CSV.name}:{key}")
    return LawData(key=key, name=name, color=color, marker=marker, dash=dash, slot=slot,
                   size=size, ke=ke, view=detected_ensemble_view(det))


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
    handles, labels = axes[0].get_legend_handles_labels()
    first = labels.index("Experiment")           # experiment leads, laws in draw order
    handles = [handles[first]] + handles[:first] + handles[first + 1:]
    labels = [labels[first]] + labels[:first] + labels[first + 1:]
    ncols = len(labels) if len(labels) <= st["legend_max_cols"] else 2
    fig.legend(handles, labels, loc="outside upper center", ncols=ncols,
               frameon=st["legend_frameon"], fontsize=st["legend_fontsize"])
    _panel_letters(axes, st)
    return fig


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
    needed = {key for keys in VARIANTS.values() for key in keys}
    by_key = {key: load_law(key, abundance, ked_ref) for key in LAWS if key in needed}
    rep = FIG_L_3["representation"]
    panels = {}
    if "l_3" in FIGURES:
        panels = {key: [pv_data(law.view, IHE_KED_REFERENCE_DIR, ked_ref, n, rep)
                        for n in FIG_L_3["n_values"]] for key, law in by_key.items()}

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(FIG_RC):                  # thesis_style font (Palatino)
        for suffix, keys in VARIANTS.items():
            laws = [by_key[key] for key in keys]
            if "l_1" in FIGURES:
                _save(build_l_1(laws), f"{FIG_L_1['stem']}_{suffix}")
            if "l_2" in FIGURES:
                _save(build_l_2(laws), f"{FIG_L_2['stem']}_{suffix}")
            if "l_3" in FIGURES:
                _save(build_l_3(laws, panels), f"{FIG_L_3['stem']}_{suffix}")
    if SHOW_FIGURE:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
