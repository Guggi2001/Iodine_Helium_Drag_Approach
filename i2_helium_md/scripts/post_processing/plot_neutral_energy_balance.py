"""Plot the neutral-stage energy-balance debug figure.

Reproduces the figure drawn at the end of
``vmi_sim_3d_neutral_propa_HeDFT_mimic.m`` (line 965) from the saved
``neutral.npz`` of a run directory. Unlike the legacy MATLAB code this
runs strictly post-hoc -- the simulation never opens a window.

Run with::

    python scripts/post_processing/plot_neutral_energy_balance.py

Output is interactive (``plt.show()``); a PNG copy is also written to
``<run>/figures/neutral_energy_balance.png``.
"""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# USER SETTINGS
# =============================================================================
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet"


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt  # noqa: E402

from i2_helium_md.postprocess import neutral_energy_totals  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def main() -> int:
    run = RunDirectory(RUN_DIR)
    neutral = run.load_neutral()
    totals = neutral_energy_totals(neutral)

    print(
        f"Loaded neutral checkpoint: N={neutral.num_molecules}, "
        f"time=[{totals.time_ps[0]:.3f}, {totals.time_ps[-1]:.3f}] ps "
        f"({totals.time_ps.size} samples)"
    )

    fig = _build_figure(totals)
    out_dir = run.root / "figures"
    out_dir.mkdir(exist_ok=True)
    fig.savefig(out_dir / "neutral_energy_balance.png", dpi=150)

    plt.show()
    return 0


def _build_figure(totals) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8.0, 4.5), constrained_layout=True)

    ax.plot(totals.time_ps, totals.E_kin_eV, label=r"$E_{kin}$", linewidth=1.5)
    ax.plot(totals.time_ps, totals.E_pot_eV, label=r"$E_{pot}$", linewidth=1.5)
    ax.plot(totals.time_ps, totals.E_dissip_eV, label=r"$E_{dissip}$", linewidth=1.5)
    ax.plot(
        totals.time_ps, totals.E_system_eV,
        label=r"$E_{system}$", color="black", linewidth=1.8,
    )

    ax.set_title("Energy balance neutral atoms")
    ax.set_xlabel("t / ps")
    ax.set_ylabel("energy / eV")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return fig


if __name__ == "__main__":
    raise SystemExit(main())
