"""Plot the ion-stage temperature-diagnostic debug figure.

Reproduces the temperature-diagnostic panel of
``vmi_sim_3d_ion_propa.m`` (lines 883-892). Two left-axis traces:
``<T'/T>`` actual and ``<T'/T>`` derived from the mass ratio. One
right-axis trace: ``<theta_lab>`` in degrees -- the lab-frame
scattering angle, *not* the COM-frame one. For heavy projectile on
light target the lab-frame cone is very narrow (max
``asin(1/rho) ~ 1.81 deg`` for I+ on He), which matches the legacy
MATLAB plot's y-axis range.

Data source: ``ion.npz`` ``temperature_diagnostic`` field (shape
``(num_steps, 3)``, columns = [T_ratio_actual, T_ratio_mass,
theta_lab_rad]). Rows where no atom collided are NaN and are masked
out before plotting.

Run with::

    python scripts/post_processing/plot_ion_temperature_diagnostic.py
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
import numpy as np  # noqa: E402

from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def main() -> int:
    run = RunDirectory(RUN_DIR)
    ion = run.load_ion()
    td = ion.temperature_diagnostic

    valid = np.isfinite(td[:, 0])
    n_valid = int(valid.sum())
    print(
        f"Loaded ion checkpoint: N={ion.num_molecules}, "
        f"{ion.time_ps.size} stored steps, {n_valid} with collisions"
    )
    if n_valid == 0:
        print(
            "WARNING: temperature_diagnostic has no valid (non-NaN) rows. "
            "Either the run had no collisions or this is a pre-v5 checkpoint "
            "(re-run the ion stage to upgrade)."
        )
        return 1

    fig = _build_figure(ion.time_ps[valid], td[valid])
    out_dir = run.root / "figures"
    out_dir.mkdir(exist_ok=True)
    fig.savefig(out_dir / "ion_temperature_diagnostic.png", dpi=150)

    plt.show()
    return 0


def _build_figure(time_ps: np.ndarray, td: np.ndarray) -> plt.Figure:
    fig, ax_left = plt.subplots(figsize=(8.5, 4.5), constrained_layout=True)
    ax_right = ax_left.twinx()

    ax_left.plot(
        time_ps, td[:, 0],
        label=r"actual $\langle T'/T\rangle$",
        linewidth=1.4,
    )
    ax_left.plot(
        time_ps, td[:, 1],
        label=r"$\langle T'/T\rangle$ from mass ratio",
        linewidth=1.4,
        linestyle="--",
    )
    ax_right.plot(
        time_ps, td[:, 2] * 180.0 / np.pi,
        label=r"$\langle\theta\rangle$",
        color="tab:red",
        linewidth=1.4,
    )

    ax_left.set_xlabel("t / ps")
    ax_left.set_ylabel(r"$\langle T'/T\rangle$")
    ax_right.set_ylabel(r"$\langle\theta\rangle$ / deg")

    lines_left, labels_left = ax_left.get_legend_handles_labels()
    lines_right, labels_right = ax_right.get_legend_handles_labels()
    # Place legend at the vertical center on the left side of the axes.
    # Using bbox_to_anchor=(0.0, 0.5) with loc="center left" anchors the
    # legend's left-center point to the axes' left-center (in axes coords).
    ax_left.legend(
        lines_left + lines_right,
        labels_left + labels_right,
        loc="center left",
        bbox_to_anchor=(0.0, 0.5),
        frameon=False,
    )

    ax_left.spines["top"].set_visible(False)
    return fig


if __name__ == "__main__":
    raise SystemExit(main())
