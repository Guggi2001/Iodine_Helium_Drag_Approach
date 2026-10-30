"""Shared style of the Tier-2 thesis figures (``data/figures/thesis_tier2/``).

One palette, font and size for every thesis figure script
(``tier2_g13_shed_timing``, ``tier2_g6_chord_map``, ``tier2_b2_thesis_figure``):

* simulation colours are samples of matplotlib's ``viridis`` (hex =
  ``viridis(x)``, x in the key; B2 blue 0.30); the experiment is one warm
  accent outside viridis, so measurement and model never share a hue; time
  and literature reference marks are black / grey;
* text and math in Palatino (the thesis font) via ``FIG_RC``, applied with
  ``plt.rc_context(FIG_RC)`` around figure creation and saving only;
* figures are drawn at the thesis text width, so the font sizes below are
  the printed sizes;
* axis labels: capitalised quantity, ``quantity / unit``; symbols in math.
"""

from __future__ import annotations

TEXTWIDTH_CM: float = 16.6                # thesis \textwidth
TEXTWIDTH_IN: float = TEXTWIDTH_CM / 2.54

VIRIDIS: dict[float, str] = {             # viridis(x) as hex
    0.00: "#440154", 0.10: "#482475", 0.30: "#355f8d", 0.33: "#31678e",
    0.50: "#21918c", 0.60: "#22a884", 0.70: "#44bf70", 0.85: "#9bd93c",
    1.00: "#fde725",
}
SIM_COLOR: str = VIRIDIS[0.30]            # production law B2 (simulation)
CONTRAST_COLOR: str = VIRIDIS[0.70]       # contrast law h405 / second simulation series
B3_COLOR: str = VIRIDIS[0.50]             # effective high-tail variant B3 (between B2 and h405)
EXP_COLOR: str = "#b2182b"                # experiment: crimson, outside viridis

# Law names as printed in the thesis figures. Keys are the repo / file-name
# labels (B2, B3, h405 stay in code, run dirs and file stems; they are
# historical); values are what a reader sees.
THESIS_NAMES: dict[str, str] = {
    "B2": "LQ",                           # production law: linear-quadratic ridge, c 1.6
    "B3": r"LQ$_\mathrm{eff}$",           # effective high-tail variant of LQ
    "h405": r"CC$_{5.4}$",                # re-landed contrast law: capped cubic, v_c 5.40 A/ps
}
LIT_COLOR: str = "0.55"                   # literature marks, shaded windows

FONT: dict[str, float] = {
    "label_fontsize": 10, "tick_fontsize": 9, "title_fontsize": 10,
    "legend_fontsize": 8.5, "annot_fontsize": 8.5,
}

FIG_RC: dict = {                          # Palatino text + math (thesis font)
    "font.family": "Palatino Linotype", "mathtext.fontset": "custom",
    "mathtext.rm": "Palatino Linotype", "mathtext.it": "Palatino Linotype:italic",
    "mathtext.bf": "Palatino Linotype:bold", "mathtext.fallback": "stix",
    "axes.linewidth": 0.8, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
}


def style_axis(ax) -> None:
    """Apply the shared tick and axis-label font sizes to ``ax``."""
    ax.tick_params(labelsize=FONT["tick_fontsize"])
    ax.xaxis.label.set_size(FONT["label_fontsize"])
    ax.yaxis.label.set_size(FONT["label_fontsize"])


def panel_letter(ax, s: str) -> None:
    """Bold ``(s)`` in the top-left corner inside ``ax``."""
    ax.text(0.01, 0.98, f"({s})", transform=ax.transAxes, va="top", ha="left",
            fontsize=FONT["label_fontsize"], fontweight="bold")
