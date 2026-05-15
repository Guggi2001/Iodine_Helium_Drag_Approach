"""Plot the experimental VMI velocity-distribution comparison.

This figure uses the realistic experimental-condition MD run and overlays
experimental gas/droplet VMI references with simulated final-velocity
histograms. It is intentionally separate from
``plot_hedft_comparison.py``, whose velocity-vs-time plot uses the 9 A
HeDFT-comparison run.

Run with:

    python scripts/post_processing/plot_experimental_comparison.py
"""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# USER SETTINGS
# =============================================================================
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet"
VMI_HE_PATH = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_he.csv"
VMI_GAS_PATH = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_gas.csv"

# Mass selection (amu) for the simulation curves.
MASS_I_HE_AMU = 131.0      # I + 1 He
MASS_I_HE2_AMU = 135.0     # I + 2 He

# MATLAB used ``edges_velocity = 0:0.04:26`` (A/ps) and displayed
# ``xlim([0, 28])`` (A/ps). Plotting now happens in m/s, so the displayed
# range is 0..2800 m/s, but the underlying bin grid stays equivalent
# (bin width 4 m/s, max 2600 m/s = 0.04 A/ps * 100, 26 A/ps * 100).
HIST_BIN_WIDTH_APS = 0.04
HIST_EDGE_MAX_APS = 26.0
HIST_NUM_BINS = int(round(HIST_EDGE_MAX_APS / HIST_BIN_WIDTH_APS))
VELOCITY_PLOT_V_MAX_MPS = 2800.0

# MATLAB used ``movmean(h, 15)`` on the simulation histogram.
HIST_SMOOTHING_WINDOW = 15


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    VmiReference,
    compute_final_velocity_histogram,
    load_vmi_reference,
)
from i2_helium_md.postprocess._smoothing import (  # noqa: E402
    moving_mean,
    normalise_trace,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def main() -> int:
    ion = RunDirectory(RUN_DIR).load_ion()
    vmi_he = load_vmi_reference(VMI_HE_PATH)
    vmi_gas = load_vmi_reference(VMI_GAS_PATH)

    print(
        f"Loaded experimental-condition ion checkpoint: N={ion.num_molecules}, "
        f"time=[{ion.time_ps[0]:.3f}, {ion.time_ps[-1]:.3f}] ps "
        f"({ion.time_ps.size} samples)"
    )

    hist_he = compute_final_velocity_histogram(
        ion,
        mass_amu=MASS_I_HE_AMU,
        num_bins=HIST_NUM_BINS,
        v_max_Aps=HIST_EDGE_MAX_APS,
    )
    hist_he2 = compute_final_velocity_histogram(
        ion,
        mass_amu=MASS_I_HE2_AMU,
        num_bins=HIST_NUM_BINS,
        v_max_Aps=HIST_EDGE_MAX_APS,
    )

    _build_experimental_velocity_figure(
        vmi_gas=vmi_gas,
        vmi_he=vmi_he,
        sim_he=hist_he,
        sim_he2=hist_he2,
        sim_smoothing_window=HIST_SMOOTHING_WINDOW,
    )

    plt.show()
    return 0


def _build_experimental_velocity_figure(
    *,
    vmi_gas: VmiReference,
    vmi_he: VmiReference,
    sim_he,
    sim_he2,
    sim_smoothing_window: int = HIST_SMOOTHING_WINDOW,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(9.5, 4.0), constrained_layout=True)
    _draw_velocity_distribution_tile(
        ax,
        vmi_gas=vmi_gas,
        vmi_he=vmi_he,
        sim_he=sim_he,
        sim_he2=sim_he2,
        sim_smoothing_window=sim_smoothing_window,
    )
    return fig


def _draw_velocity_distribution_tile(
    ax: plt.Axes,
    *,
    vmi_gas: VmiReference,
    vmi_he: VmiReference,
    sim_he,
    sim_he2,
    sim_smoothing_window: int,
) -> None:
    """Experimental I+ gas, I+He droplet, simulation I+He, I+He2."""
    # 5-element discrete colormap that approximates colorcet('L08','N',5).
    palette = plt.colormaps["plasma"](np.linspace(0.05, 0.85, 5))
    c_gas, c_he, c_sim_he, c_sim_he2, _ = palette

    # Gas-phase background mask: > 400 m/s (was > 4 A/ps in MATLAB).
    mask_gas = vmi_gas.velocity_mps > 400.0
    max_gas = float(vmi_gas.signal_arb[mask_gas].max())
    max_he = float(vmi_he.signal_arb.max())
    sim_he_density = normalise_trace(
        moving_mean(sim_he.density, sim_smoothing_window)
    )
    sim_he2_density = normalise_trace(
        moving_mean(sim_he2.density, sim_smoothing_window)
    )

    ax.plot(
        vmi_gas.velocity_mps,
        vmi_gas.signal_arb / max_gas,
        color=c_gas,
        linewidth=2.0,
        label=r"$I_2$:$I^+$",
    )
    ax.plot(
        vmi_he.velocity_mps,
        vmi_he.signal_arb / max_he,
        linestyle=":",
        color=c_he,
        linewidth=2.0,
        label=r"$I_2 He_N$:$I^+ He$",
    )
    ax.plot(
        sim_he.bin_centers_mps,
        sim_he_density,
        linestyle="--",
        color=c_sim_he,
        linewidth=2.0,
        label=r"simulation $I^+ He$",
    )
    ax.plot(
        sim_he2.bin_centers_mps,
        sim_he2_density,
        linestyle="-.",
        color=c_sim_he2,
        linewidth=2.0,
        label=r"simulation $I^+ He_2$",
    )

    ax.set_xlim(0.0, VELOCITY_PLOT_V_MAX_MPS)
    ax.set_ylim(0.0, 1.1)
    ax.set_xlabel("v / m/s")
    ax.set_ylabel("signal / arb. units")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


if __name__ == "__main__":
    raise SystemExit(main())
