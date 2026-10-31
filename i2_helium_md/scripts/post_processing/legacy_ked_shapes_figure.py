"""Thesis figure: KED shapes of the legacy hard-sphere model, n = 0-3.

The B2-3 figure (``tier2_b2_thesis_figure``, ``b2_3_ked_shapes``) for the
legacy collision run ``data/runs/single_pulse_droplet_long`` (hard-sphere
Mode 3, sigma ~ v^-2, ``E_coulomb_scale`` 0.8): experiment raw + movmean vs
simulation, 3-D P(E) on a kinetic-energy axis (the committed P(v) recipe
through the Jacobian P(E) = P(v) / (m v), ``ked_traces``), peak-normalized,
<E> markers, per-n range, here for n = 0-3 (2 x 2).

The run has no Tier-2 detection stage: the ion checkpoint (``ion.npz``) is
read directly -- ``pv_data`` takes any ``IonCheckpoint`` (the detected view
of the B2 figures is a stand-in for one). The mass gate is the shared
``select_final_mass_gate`` (final mass within 0.5 amu, outside the droplet).
Data and draw layer are imported from the B2 script; nothing is re-derived.

There is no committed score table for this run, so the quoted numbers are
written next to the figure as ``legacy_ked_shapes_points.csv`` (n, ions in
the gate, simulated and experimental <E>). The run directory is not tracked
(``ion.npz`` 528 MB): the figure regenerates only where the run exists.

Writes ``data/figures/thesis_legacy/legacy_ked_shapes.png`` (+ the CSV).

Usage::

    python scripts/post_processing/legacy_ked_shapes_figure.py
"""

from __future__ import annotations

import csv
import sys
import warnings
from pathlib import Path

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
OUT_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_legacy"
FORMATS: tuple[str, ...] = ("png",)
DPI: int = 200
SHOW_FIGURE: bool = False

FIG_LEGACY_KED: dict = {
    # B2-3 (FIG_B2_3) with these keys overridden
    "stem": "legacy_ked_shapes",
    "n_values": (0, 1, 2, 3),
    "ncols": 2,
    "figsize": (TEXTWIDTH_IN, 5.0),
    "sim": dict(color=LEGACY_COLOR),
    "sim_label": "Hard-sphere model",
    # n = 0: exp. 3-D columns NaN below 0.4 eV (Abel-centre spike), end at the
    # 6.1 eV edge of the speed window; n = 3 widened past B2-3's 1.5 eV so the
    # legacy tail (to ~1.8 eV) is not cut
    "e_xlim_by_n": {0: (0.0, 6.0), 1: (0.0, 3.5), 2: (0.0, 2.5), 3: (0.0, 2.0)},
}

# =============================================================================
# Imports (after the settings so SHOW_FIGURE picks the backend)
# =============================================================================
import matplotlib  # noqa: E402

if not SHOW_FIGURE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from i2_helium_md.postprocess import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.post_processing.tier2_b2_thesis_figure import (  # noqa: E402
    FIG_B2_3,
    PvData,
    build_b2_3,
    pv_data,
)


def legacy_style() -> dict:
    """``FIG_B2_3`` with the legacy overrides (one level, as ``merged_style``)."""
    return FIG_B2_3 | FIG_LEGACY_KED


def legacy_panels(ion, ked_ref, n_values) -> list[PvData]:
    """The B2-3 panels of an ion checkpoint (3-D representation)."""
    return [pv_data(ion, IHE_KED_REFERENCE_DIR, ked_ref, n, "3d") for n in n_values]


def points_rows(panels: list[PvData]) -> list[dict]:
    """The quoted numbers per panel: n, ions in the gate, simulated / exp. <E> [eV]."""
    return [{"n": d.n, "num_sim": d.num_sim,
             "mean_E_sim_eV": "" if d.mean_E_sim_eV is None else d.mean_E_sim_eV,
             "mean_E_exp_eV": d.mean_E_exp_eV} for d in panels]


def main() -> int:
    with warnings.catch_warnings():
        # schema v5 predates the E_int reservoir; the loader zero-fills it (unread here)
        warnings.filterwarnings("ignore", message="ion checkpoint predates")
        ion = RunDirectory(RUN_DIR).load_ion()
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_DIR / "IHe_KED_reference.csv")
    st = legacy_style()
    panels = legacy_panels(ion, ked_ref, st["n_values"])

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(FIG_RC):                  # thesis_style font (Palatino)
        fig = build_b2_3(panels, st)
        for ext in FORMATS:
            path = OUT_DIR / f"{st['stem']}.{ext}"
            fig.savefig(path, dpi=DPI)
            print(f"[legacy] wrote {path}")
    rows = points_rows(panels)
    csv_path = OUT_DIR / f"{st['stem']}_points.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for r in rows:
        print(f"  n = {r['n']}: N {r['num_sim']}, <E> sim {r['mean_E_sim_eV']:.4f} / "
              f"exp {r['mean_E_exp_eV']:.4f} eV")
    print(f"[legacy] wrote {csv_path}")
    if SHOW_FIGURE:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
