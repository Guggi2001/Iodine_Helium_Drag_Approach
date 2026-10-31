"""Thesis figure: mean KE and size distribution of the legacy hard-sphere model, n = 0-20.

Legacy collision run ``data/runs/single_pulse_droplet_long`` (hard-sphere
Mode 3, sigma ~ v^-2, ``E_coulomb_scale`` 0.8) vs experiment, two panels on a
shared n axis:

* **(a) mean KE vs n** -- the B2-2 (a) panel (``draw_ke_mean``: correlated
  calibration / condition bands, per-point bars, n = 1 median; simulation
  with SE, open markers below ``DETECTED_KE_MIN_BIN_COUNT`` ions) on a log
  axis, here from n = 0 (bare I+) to n = 17, where the experimental KE
  reference ends (``ke_n_max``; the simulation is cut there too);
* **(b) size distribution** -- the B2-1 (b) side-by-side bars over **all**
  ions, n = 0 included, not renormalised: both sides as fractions of
  n = 0..``N_MAX`` (the experimental abundance window). The legacy run
  carries fragments above ``N_MAX``; their share of all ions is written to
  the points CSV (``frac_sim_all_ions``) for the caption.

A second figure, ``legacy_ke_split_and_size``, splits (a) as B2-2 split
does: n = 0..``split_n``-1 | n = ``split_n``..17, linear, each on its own
scale, over the same size bars (c).

The run has no Tier-2 detection stage; the ion checkpoint is read directly
through the shared mass gate: per-n KE from ``fragment_mean_kinetic_energy``
(E = m(n) v^2 / 2, SE ddof = 1), counts from ``fragment_gate_counts``. No
committed score table exists for this run, so the plotted numbers go to
``legacy_ke_and_size_points.csv`` next to the figure. The run directory is not
tracked: the figure regenerates only where the run exists.

Usage::

    python scripts/post_processing/legacy_ke_size_figure.py
"""

from __future__ import annotations

import csv
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT_DIR = Path(__file__).resolve().parent
for _p in (PROJECT_ROOT, _SCRIPT_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from scripts.post_processing.thesis_style import FIG_RC, LEGACY_COLOR, TEXTWIDTH_IN  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================
RUN_DIR: Path = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet_long"
IHE_KED_REFERENCE_DIR: Path = PROJECT_ROOT / "data" / "reference" / "ihe_ked"
ABUNDANCE_REF_PATH: Path = PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"
OUT_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_legacy"
N_MAX: int = 20                     # the abundance reference ends at n = 20
FORMATS: tuple[str, ...] = ("png",)
DPI: int = 200
SHOW_FIGURE: bool = False

FIG_LEGACY_KE_SIZE: dict = {
    # B2-2 (FIG_B2_2) + the B2-1 bar keys, with these overrides
    "stem": "legacy_ke_and_size",
    "figsize": (TEXTWIDTH_IN, 5.0),
    "height_ratios": (1.2, 1.0),
    "sim": dict(color=LEGACY_COLOR),
    "sim_label": "Hard-sphere model",
    "yscale": "log",                     # (a)
    "n_xlim": (-0.7, N_MAX + 0.7),
    "bar_width_pair": 0.4,               # (b), as B2-1 (b)
    "legend_loc": "upper right",         # (a)
    "legend_loc_size": "upper right",    # (b)
    "ylabel_size": "Fraction / %",
    # split variant (B2-2 split): (a) n = 0..split-1 | (b) n = split..17,
    # linear, own scales, over (c) the size bars
    "stem_split": "legacy_ke_split_and_size",
    "figsize_split": (TEXTWIDTH_IN, 5.0),
    "height_ratios_split": (1.15, 1.0),
    "split_n": 5,                        # = FIG_B2_2["split_n"]
    "titles_split": ("$n$ = 0–{hi}", "$n$ = {lo}–{n_max}"),
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
    fragment_gate_counts,
    fragment_mean_kinetic_energy,
    load_he_abundance_reference,
    load_ihe_ked_reference,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.post_processing.tier2_b2_thesis_figure import (  # noqa: E402
    FIG_B2_1,
    FIG_B2_2,
    KEData,
    _finish_axis,
    _legend,
    _panel_letters,
    draw_ke_mean,
    merged_style,
)
from summary_sections import DETECTED_KE_MIN_BIN_COUNT  # noqa: E402


# =============================================================================
# Data layer
# =============================================================================
@dataclass(frozen=True)
class LegacySizeData:
    """Size fractions over n = 0..n_max (both sides sum to 1 over that window).

    ``frac_sim_all_ions`` is the same count over all ions outside the droplet,
    so ``1 - frac_sim_all_ions.sum()`` is the share beyond the window.
    """
    n: np.ndarray
    count: np.ndarray
    frac_sim: np.ndarray
    frac_sim_all_ions: np.ndarray
    frac_exp: np.ndarray


def ke_n_max(ked_ref) -> int:
    """Last n of the experimental KE reference (17): the simulation KE stops there too."""
    return int(np.asarray(ked_ref.n).max())


def legacy_ke_data(ion, ked_ref) -> KEData:
    """Per-n mean KE [eV] of the mass-gated ions, n = 0..``ke_n_max`` (empty
    gates skipped), as a ``KEData`` for ``draw_ke_mean``; no scores (legacy run)."""
    rows = []
    for n in range(ke_n_max(ked_ref) + 1):
        try:
            rows.append(fragment_mean_kinetic_energy(ion, n))
        except ValueError:            # empty gate
            continue
    n_arr = np.array([r.n for r in rows], dtype=int)
    count = np.array([r.num_atoms_used for r in rows], dtype=int)
    nan = np.full(n_arr.size, np.nan)
    return KEData(
        n=n_arr, mean_eV=np.array([r.mean_KE_eV for r in rows]),
        se_eV=np.array([r.stat_err_mean_KE_eV for r in rows]), count=count,
        scored_mask=count >= DETECTED_KE_MIN_BIN_COUNT, ratio=nan, ratio_err=nan,
        ref=ked_ref, chi2_score=None, midhot=None, deepke=None,
    )


def legacy_size_data(ion, abundance, n_max: int) -> LegacySizeData:
    """Mass-gate counts per n (``fragment_gate_counts``) vs the experimental
    abundance, both as fractions of n = 0..n_max."""
    n = np.arange(n_max + 1)
    count = np.asarray(fragment_gate_counts(ion, n), dtype=int)
    num_outside = 2 * int(np.count_nonzero(ion.b_ion_outside))
    exp_of = dict(zip(np.asarray(abundance.n).tolist(),
                      np.asarray(abundance.ion_fraction, dtype=float).tolist()))
    frac_exp = np.array([exp_of.get(int(k), 0.0) for k in n])
    return LegacySizeData(n=n, count=count, frac_sim=count / count.sum(),
                          frac_sim_all_ions=count / num_outside,
                          frac_exp=frac_exp / frac_exp.sum())


def points_rows(ke: KEData, size: LegacySizeData) -> list[dict]:
    """One row per n: size fractions, simulated <E> +- SE, experimental <E> [eV]."""
    ke_of = {int(k): (m, s) for k, m, s in zip(ke.n, ke.mean_eV, ke.se_eV)}
    exp_of = dict(zip(np.asarray(ke.ref.n).tolist(), np.asarray(ke.ref.mean_KE_eV).tolist()))
    return [{"n": int(k), "count": int(c), "frac_sim": f, "frac_sim_all_ions": fa,
             "frac_exp": fe,
             "mean_E_sim_eV": ke_of.get(int(k), ("", ""))[0],
             "se_E_sim_eV": ke_of.get(int(k), ("", ""))[1],
             "mean_E_exp_eV": exp_of.get(int(k), "")}
            for k, c, f, fa, fe in zip(size.n, size.count, size.frac_sim,
                                       size.frac_sim_all_ions, size.frac_exp)]


# =============================================================================
# Draw layer + figure
# =============================================================================
def draw_size_all(ax, d: LegacySizeData, st: dict) -> None:
    """B2-1 (b) bars over all ions, n = 0 included (fractions in %)."""
    w = st["bar_width_pair"]
    ax.bar(d.n - w / 2, 100.0 * d.frac_exp, width=w, label="Experiment", **st["exp"])
    ax.bar(d.n + w / 2, 100.0 * d.frac_sim, width=w, label=st["sim_label"], **st["sim"])
    ax.set(xlim=st["n_xlim"], xlabel=st["xlabel"], ylabel=st["ylabel_size"])
    _legend(ax, st, st["legend_loc_size"])
    _finish_axis(ax, st, integer_x=True)


def build_legacy_ke_size(ke: KEData, size: LegacySizeData,
                         style: dict | None = None) -> plt.Figure:
    """(a) mean KE vs n from n = 0, log axis, over (b) the all-ion size bars."""
    st = merged_style(FIG_B2_1 | FIG_B2_2 | (FIG_LEGACY_KE_SIZE if style is None else style))
    fig, (ax_e, ax_s) = plt.subplots(2, 1, figsize=st["figsize"], sharex=True,
                                     constrained_layout=True,
                                     gridspec_kw={"height_ratios": st["height_ratios"]})
    draw_ke_mean(ax_e, ke, st, 0, ke_n_max(ke.ref), show_legend=True, show_chi2=False,
                 n_min=0)
    ax_e.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    ax_e.set_xlabel("")
    draw_size_all(ax_s, size, st)
    _panel_letters((ax_e, ax_s), st)
    return fig


def build_legacy_ke_split_size(ke: KEData, size: LegacySizeData,
                               style: dict | None = None) -> plt.Figure:
    """B2-2 split layout from n = 0: (a) n = 0..split-1 | (b) n = split..17
    (``ke_n_max``), both linear on their own scales, over (c) the all-ion
    size bars."""
    st = merged_style(FIG_B2_1 | FIG_B2_2 | (FIG_LEGACY_KE_SIZE if style is None else style))
    split, n_max = st["split_n"], ke_n_max(ke.ref)
    fig, axd = plt.subplot_mosaic([["a", "b"], ["c", "c"]], figsize=st["figsize_split"],
                                  constrained_layout=True,
                                  gridspec_kw={"height_ratios": st["height_ratios_split"]})
    st_lin = st | {"yscale": "linear"}
    draw_ke_mean(axd["a"], ke, st_lin, 0, split - 1, show_legend=False, show_chi2=False,
                 n_min=0)
    draw_ke_mean(axd["b"], ke, st_lin, split, n_max, show_legend=False, show_chi2=False,
                 n_min=0)
    axd["b"].set_ylabel("")
    t_lo, t_hi = st["titles_split"]
    axd["a"].set_title(t_lo.format(hi=split - 1), fontsize=st["title_fontsize"])
    axd["b"].set_title(t_hi.format(lo=split, n_max=n_max), fontsize=st["title_fontsize"])
    # one legend for (a) + (b) above the figure, as B2-2 split
    handles = {}
    for key in ("a", "b"):
        for h, lab in zip(*axd[key].get_legend_handles_labels()):
            handles.setdefault(lab, h)
    fig.legend(handles.values(), handles.keys(), loc="outside upper center", ncols=3,
               frameon=st["legend_frameon"], fontsize=st["legend_fontsize"])
    draw_size_all(axd["c"], size, st)
    _panel_letters((axd["a"], axd["b"], axd["c"]), st)
    return fig


def main() -> int:
    with warnings.catch_warnings():
        # schema v5 predates the E_int reservoir; the loader zero-fills it (unread here)
        warnings.filterwarnings("ignore", message="ion checkpoint predates")
        ion = RunDirectory(RUN_DIR).load_ion()
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_DIR / "IHe_KED_reference.csv")
    abundance = load_he_abundance_reference(ABUNDANCE_REF_PATH)
    ke = legacy_ke_data(ion, ked_ref)
    size = legacy_size_data(ion, abundance, N_MAX)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = FIG_LEGACY_KE_SIZE["stem"]
    with plt.rc_context(FIG_RC):                  # thesis_style font (Palatino)
        for fig_stem, build in ((stem, build_legacy_ke_size),
                                (FIG_LEGACY_KE_SIZE["stem_split"], build_legacy_ke_split_size)):
            fig = build(ke, size)
            for ext in FORMATS:
                path = OUT_DIR / f"{fig_stem}.{ext}"
                fig.savefig(path, dpi=DPI)
                print(f"[legacy] wrote {path}")
    rows = points_rows(ke, size)
    csv_path = OUT_DIR / f"{stem}_points.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  n = 0: sim {100 * size.frac_sim[0]:.2f} % / exp {100 * size.frac_exp[0]:.2f} %; "
          f"share of all ions beyond n = {N_MAX}: "
          f"{100 * (1.0 - size.frac_sim_all_ions.sum()):.2f} %")
    print(f"[legacy] wrote {csv_path}")
    if SHOW_FIGURE:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
