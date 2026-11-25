"""Thesis figure: anisotropy parameter beta per fragment, LQ (B2) vs experiment.

One panel, written to ``data/figures/thesis_results_discussion/b2_5_anisotropy``:
beta of I+He_n against n for

* the production law B2 (printed as LQ), pooled N = 6000 container
  ``B2_N6000_pooled`` (``self_bound_gate_release = decoupled``), and
* the experiment, ``data/reference/ihe_ked/IHe_KED_beta.csv`` (v2, corrected
  inversion about the polarisation axis), column ``beta_mean``.

Convention (``THESIS_DISCUSSION.md`` §13.1): I(theta) ∝ 1 + beta·P2(cos theta),
theta measured from the laser polarisation (lab x in the model, image X in the
experiment); beta = 2 pure cos², 0 isotropic. Reference lines at both limits.

Model beta: 5·<P2(v_x / |v|)> over the ions of one fragment, from the terminal
detected velocities. Fragment selection is the shared mass gate
(``detected_ensemble_view`` + ``select_final_mass_gate`` at m(n)), the same
detected ensemble as every other B2 thesis figure (droplet-retained ions
excluded, suppressed ions scored at n = 0). Error bar: standard error of the
mean of 5·P2. Bins with fewer than ``MIN_COUNT`` ions are not drawn.

Experimental beta: pyabel ``rIbeta`` per speed bin, P(v)-weighted over the
central <E> signal mask. Diagnostic, **no sysErr** in the export, so drawn
without error bars. Experimental n = 0 is droplet-born I+ only (gas-phase
background subtracted); model n = 0 is every bare I+ of the detected ensemble.

Anchor guard: before drawing, model beta for n = 0-7 is compared with the
committed values of ``THESIS_DISCUSSION.md`` §13.3 (2 dp); a difference is reported.

Edit the ``USER SETTINGS`` block and run::

    python scripts/post_processing/tier2_anisotropy_figure.py
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# =============================================================================
# USER SETTINGS
# =============================================================================
RUN_DIR: Path = PROJECT_ROOT / "data" / "runs" / "B2_N6000_pooled"
BETA_REF_PATH: Path = PROJECT_ROOT / "data" / "reference" / "ihe_ked" / "IHe_KED_beta.csv"
OUT_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_results_discussion"
STEM: str = "b2_5_anisotropy"
FORMATS: tuple[str, ...] = ("png",)
DPI: int = 200
SHOW_FIGURE: bool = False

N_MAX_MODEL: int = 17          # highest n considered for the model
MIN_COUNT: int = 20            # model bins with fewer ions are not drawn
ANCHOR_CHECK: bool = True
# THESIS_DISCUSSION.md §13.3, model beta (B2), n = 0..7
ANCHOR_BETA_MODEL: dict[int, float] = {
    0: 1.88, 1: 2.04, 2: 1.97, 3: 1.98, 4: 1.98, 5: 1.83, 6: 1.93, 7: 2.04,
}

WIDTH_FRACTION: float = 1.0    # printed at WIDTH_FRACTION·\linewidth: font sizes stay the printed sizes
FIG_HEIGHT_CM: float = 6.5
Y_LIM: tuple[float, float] = (-0.25, 2.45)

from scripts.post_processing.thesis_style import (  # noqa: E402  (shared thesis style)
    EXP_COLOR, FIG_RC, FONT, LIT_COLOR, SIM_COLOR, TEXTWIDTH_IN, THESIS_NAMES, style_axis,
)

import matplotlib  # noqa: E402

if not SHOW_FIGURE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

from i2_helium_md.physics.shell_schedule import complex_mass_amu  # noqa: E402
from i2_helium_md.postprocess.detected_view import detected_ensemble_view  # noqa: E402
from i2_helium_md.postprocess.velocity_distribution import select_final_mass_gate  # noqa: E402
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from i2_helium_md.simulation.thesis_runs import require_run  # noqa: E402
from scripts.post_processing.anchors import report  # noqa: E402


# =============================================================================
# Data layer
# =============================================================================
@dataclass(frozen=True)
class BetaSeries:
    """Anisotropy parameter per fragment.

    Attributes
    ----------
    n : np.ndarray, shape (k,), int
        Number of attached He.
    beta : np.ndarray, shape (k,)
        Anisotropy parameter (dimensionless, -1 ... 2).
    beta_se : np.ndarray | None, shape (k,)
        Standard error of ``beta``; ``None`` when the source has none.
    count : np.ndarray | None, shape (k,), int
        Ions per fragment (model only).
    """

    n: np.ndarray
    beta: np.ndarray
    beta_se: np.ndarray | None = None
    count: np.ndarray | None = None


def model_beta(run_dir: Path, n_max: int, min_count: int) -> BetaSeries:
    """Model beta = 5·<P2(v_x/|v|)> per fragment of the detected ensemble.

    Parameters
    ----------
    run_dir
        Run (or pooled container) directory holding ``detection.npz``.
    n_max
        Highest He count considered.
    min_count
        Fragments with fewer gated ions are dropped (must be >= 2 for an SE).

    Returns
    -------
    BetaSeries
        n, beta, standard error and ion count for every kept fragment.

    Raises
    ------
    FileNotFoundError
        If ``detection.npz`` is missing.
    ValueError
        If ``min_count < 2`` or a gated ion has zero speed.
    """
    if min_count < 2:
        raise ValueError(f"min_count must be >= 2 for a standard error, got {min_count}.")
    path = run_dir / "detection.npz"
    if not path.is_file():
        raise FileNotFoundError(f"no detection.npz in {run_dir}")
    view = detected_ensemble_view(load_detection_result(path))
    vx, vy, vz = view.velocities_final_x, view.velocities_final_y, view.velocities_final_z
    speed_Aps = np.sqrt(vx**2 + vy**2 + vz**2)

    rows = []
    for n in range(n_max + 1):
        gate = select_final_mass_gate(view, mass_amu=complex_mass_amu(n))
        count = int(gate.sum())
        if count < min_count:
            continue
        if np.any(speed_Aps[gate] <= 0.0):
            raise ValueError(f"n = {n}: gated ion with zero speed, direction undefined.")
        cos_theta = vx[gate] / speed_Aps[gate]          # angle to the polarisation (lab x)
        five_p2 = 5.0 * 0.5 * (3.0 * cos_theta**2 - 1.0)
        rows.append((n, five_p2.mean(), five_p2.std(ddof=1) / np.sqrt(count), count))
    n, beta, se, cnt = (np.array(c) for c in zip(*rows))
    return BetaSeries(n=n.astype(int), beta=beta, beta_se=se, count=cnt.astype(int))


def experiment_beta(path: Path, column: str = "beta_mean") -> BetaSeries:
    """Experimental beta per fragment from ``IHe_KED_beta.csv``.

    Parameters
    ----------
    path
        The committed export (columns ``n, beta_mean, beta_core``).
    column
        ``beta_mean`` (P(v)-weighted over the signal mask) or ``beta_core``.

    Returns
    -------
    BetaSeries
        n and beta; no uncertainty (the export carries none).

    Raises
    ------
    FileNotFoundError, ValueError
        Missing file or missing column.
    """
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or "n" not in reader.fieldnames or column not in reader.fieldnames:
            raise ValueError(f"{path.name}: need columns 'n' and {column!r}, got {reader.fieldnames}")
        rows = [(int(r["n"]), float(r[column])) for r in reader]
    n, beta = (np.array(c) for c in zip(*rows))
    return BetaSeries(n=n.astype(int), beta=beta)


def anchor_differences(model: BetaSeries, anchor: dict[int, float]) -> list[str]:
    """Model beta values that differ from the committed §13.3 values at 2 dp."""
    lookup = dict(zip(model.n.tolist(), model.beta.tolist()))
    return [f"THESIS_DISCUSSION 13.3 beta n = {n}: {lookup.get(n)} vs {ref}"
            for n, ref in anchor.items() if n not in lookup or round(lookup[n], 2) != ref]


# =============================================================================
# Draw layer
# =============================================================================
def build_figure(model: BetaSeries, exp: BetaSeries) -> plt.Figure:
    """beta vs n, model with SE bars and experiment, cos² / isotropic limits."""
    fig, ax = plt.subplots(figsize=(WIDTH_FRACTION * TEXTWIDTH_IN, FIG_HEIGHT_CM / 2.54),
                           constrained_layout=True)

    for y, text in ((2.0, r"pure $\cos^2$"), (0.0, "isotropic")):
        ax.axhline(y, color=LIT_COLOR, ls="--", lw=0.8, zorder=0)
        ax.text(0.995, y, text, transform=ax.get_yaxis_transform(), ha="right",
                va="bottom", fontsize=FONT["annot_fontsize"], color=LIT_COLOR)

    ax.errorbar(model.n, model.beta, yerr=model.beta_se, fmt="o", ms=4.5, color=SIM_COLOR,
                capsize=2, lw=1.0, label=THESIS_NAMES["B2"])
    ax.plot(exp.n, exp.beta, "s", ms=4.5, color=EXP_COLOR, label="Experiment")

    n_hi = int(max(model.n.max(), exp.n.max()))
    ax.set_xlim(-0.6, n_hi + 0.6)
    ax.set_ylim(*Y_LIM)
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.set_xlabel(r"Number of attached He atoms $n$")
    ax.set_ylabel(r"Anisotropy parameter $\beta$")
    style_axis(ax)
    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index(THESIS_NAMES["B2"]), labels.index("Experiment")]
    ax.legend([handles[i] for i in order], [labels[i] for i in order], loc="center left",
              bbox_to_anchor=(0.0, 0.47), frameon=False, fontsize=FONT["legend_fontsize"])
    return fig


# =============================================================================
# Main
# =============================================================================
def main() -> int:
    require_run(RUN_DIR / "detection.npz", "gen_thesis_laws.py")
    model = model_beta(RUN_DIR, N_MAX_MODEL, MIN_COUNT)
    exp = experiment_beta(BETA_REF_PATH)
    if ANCHOR_CHECK:
        report(anchor_differences(model, ANCHOR_BETA_MODEL),
               "[anisotropy] anchor OK vs THESIS_DISCUSSION.md §13.3")
    for n, b, se, c in zip(model.n, model.beta, model.beta_se, model.count):
        print(f"[anisotropy] model n={n:2d}  N={c:5d}  beta={b:.3f} ± {se:.3f}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(FIG_RC):                 # thesis_style font (Palatino)
        fig = build_figure(model, exp)
        for ext in FORMATS:
            out = OUT_DIR / f"{STEM}.{ext}"
            fig.savefig(out, dpi=DPI)
            print(f"[anisotropy] wrote {out}")
    if SHOW_FIGURE:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
