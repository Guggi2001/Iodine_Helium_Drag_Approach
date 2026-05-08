"""Recreate the minimal paper figure of post_process_single_pulse_paper_v3.

Reproduces the simulation-side panels of
``post_process_single_pulse_paper_v3.m`` that can be driven from
existing run-directory artifacts and the 1-D experimental VMI
references in ``data/reference/``:

* Tile A -- experimental gas/droplet VMI velocity references overlaid
  with simulated final-velocity histograms for I+He and I+He2
  (mirrors lines ~76-100 + line 280 of the legacy script).
* Tile B -- simulated final azimuthal phi histogram, mass-selected
  (mirrors line 314).
* Tile C -- simulated final ion mass spectrum (mirrors line 397).

Out of scope (deferred): the polar VMI image surface (line 111),
cos^2 angular anisotropy fit (Tile 2 of legacy figure 1), and the
beta(v) function (legacy figure 3). These require a 2-D polar VMI
image not present in this repository and CLAUDE.md flags them as
"full experimental VMI interpretation" which is out of default
scope.

Saves both ``compare_simulation_and_measurement.pdf`` (filename
matches the legacy MATLAB output) and a PNG copy under
``<run>/figures/``.

Run with::

    python scripts/post_processing/plot_paper_figure.py
"""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# USER SETTINGS
# =============================================================================
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet"
VMI_HE_PATH = PROJECT_ROOT / "data" / "reference" / "vmi_iplus_he.csv"
VMI_GAS_PATH = PROJECT_ROOT / "data" / "reference" / "vmi_iplus_gas.csv"

# Mass selection (amu) for the simulated curves in tiles A and B.
MASS_I_HE_AMU = 131.0      # I + 1 He
MASS_I_HE2_AMU = 135.0     # I + 2 He

# Tile A (velocity histogram) -- MATLAB used edges 0:0.04:26 and xlim 0..28.
HIST_BIN_WIDTH_APS = 0.04
HIST_EDGE_MAX_APS = 26.0
HIST_NUM_BINS = int(round(HIST_EDGE_MAX_APS / HIST_BIN_WIDTH_APS))
HIST_SMOOTHING_WINDOW = 15
VELOCITY_PLOT_V_MAX_APS = 28.0

# Tile B (phi histogram).
PHI_BIN_WIDTH_RAD = 0.05
PHI_SMOOTHING_WINDOW = 15

# Tile C (mass spectrum).
MASS_BIN_WIDTH_AMU = 1.0


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    compute_final_velocity_histogram,
    ion_energy_totals,  # noqa: F401  (kept for ad-hoc debugging)
    load_vmi_reference,
    mass_spectrum,
    phi_histogram,
)
from i2_helium_md.postprocess._smoothing import (  # noqa: E402
    moving_mean,
    normalise_trace,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def main() -> int:
    run = RunDirectory(RUN_DIR)
    ion = run.load_ion()
    vmi_he = load_vmi_reference(VMI_HE_PATH)
    vmi_gas = load_vmi_reference(VMI_GAS_PATH)

    print(
        f"Loaded ion checkpoint: N={ion.num_molecules}, "
        f"time=[{ion.time_ps[0]:.3f}, {ion.time_ps[-1]:.3f}] ps "
        f"({ion.time_ps.size} samples)"
    )

    sim_he = compute_final_velocity_histogram(
        ion, mass_amu=MASS_I_HE_AMU,
        num_bins=HIST_NUM_BINS, v_max_Aps=HIST_EDGE_MAX_APS,
    )
    sim_he2 = compute_final_velocity_histogram(
        ion, mass_amu=MASS_I_HE2_AMU,
        num_bins=HIST_NUM_BINS, v_max_Aps=HIST_EDGE_MAX_APS,
    )
    phi_he = phi_histogram(
        ion, bin_width_rad=PHI_BIN_WIDTH_RAD, mass_amu=MASS_I_HE_AMU,
    )
    masses = mass_spectrum(ion, bin_width_amu=MASS_BIN_WIDTH_AMU)

    fig = _build_figure(
        vmi_gas=vmi_gas, vmi_he=vmi_he,
        sim_he=sim_he, sim_he2=sim_he2,
        phi_he=phi_he, masses=masses,
    )

    out_dir = run.root / "figures"
    out_dir.mkdir(exist_ok=True)
    fig.savefig(out_dir / "compare_simulation_and_measurement.pdf")
    fig.savefig(out_dir / "compare_simulation_and_measurement.png", dpi=150)
    print(f"Saved figure to {out_dir / 'compare_simulation_and_measurement.pdf'}")

    plt.show()
    return 0


def _build_figure(
    *,
    vmi_gas, vmi_he, sim_he, sim_he2, phi_he, masses,
) -> plt.Figure:
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2), constrained_layout=True)
    _draw_velocity_tile(axes[0], vmi_gas, vmi_he, sim_he, sim_he2)
    _draw_phi_tile(axes[1], phi_he)
    _draw_mass_tile(axes[2], masses)
    return fig


def _draw_velocity_tile(ax, vmi_gas, vmi_he, sim_he, sim_he2) -> None:
    palette = plt.colormaps["plasma"](np.linspace(0.05, 0.85, 5))
    c_gas, c_he, c_sim_he, c_sim_he2, _ = palette

    mask_gas = vmi_gas.velocity_Aps > 4.0
    max_gas = float(vmi_gas.signal_arb[mask_gas].max())
    max_he = float(vmi_he.signal_arb.max())
    sim_he_density = normalise_trace(
        moving_mean(sim_he.density, HIST_SMOOTHING_WINDOW)
    )
    sim_he2_density = normalise_trace(
        moving_mean(sim_he2.density, HIST_SMOOTHING_WINDOW)
    )

    ax.plot(
        vmi_gas.velocity_Aps, vmi_gas.signal_arb / max_gas,
        color=c_gas, linewidth=2.0, label=r"$I_2$:$I^+$",
    )
    ax.plot(
        vmi_he.velocity_Aps, vmi_he.signal_arb / max_he,
        color=c_he, linewidth=2.0, linestyle=":",
        label=r"$I_2 He_N$:$I^+ He$",
    )
    ax.plot(
        sim_he.bin_centers_Aps, sim_he_density,
        color=c_sim_he, linewidth=2.0, linestyle="--",
        label=r"simulation $I^+ He$",
    )
    ax.plot(
        sim_he2.bin_centers_Aps, sim_he2_density,
        color=c_sim_he2, linewidth=2.0, linestyle="-.",
        label=r"simulation $I^+ He_2$",
    )

    ax.set_xlim(0.0, VELOCITY_PLOT_V_MAX_APS)
    ax.set_ylim(0.0, 1.1)
    ax.set_xlabel(r"v / $\mathrm{\AA}/\mathrm{ps}$")
    ax.set_ylabel("signal / arb. units")
    ax.set_title("(a) radial velocity")
    ax.legend(frameon=False, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _draw_phi_tile(ax, phi_he) -> None:
    if phi_he.density.sum() > 0:
        smoothed = moving_mean(phi_he.density, PHI_SMOOTHING_WINDOW)
        normalised = normalise_trace(smoothed)
    else:
        normalised = np.zeros_like(phi_he.density, dtype=float)
    ax.plot(phi_he.bin_centers_rad, normalised, linewidth=2.0)

    ax.set_xlim(0.0, 2.0 * np.pi)
    ax.set_ylim(0.0, 1.1)
    ax.set_xlabel(r"$\varphi$ / rad")
    ax.set_ylabel("signal / arb. units")
    ax.set_title(rf"(b) azimuth, m={MASS_I_HE_AMU:.0f} u, n={phi_he.num_atoms_used}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _draw_mass_tile(ax, masses) -> None:
    ax.bar(
        masses.bin_centers_amu, masses.counts,
        width=MASS_BIN_WIDTH_AMU * 0.9, edgecolor="black", linewidth=0.5,
    )
    ax.set_xlabel("m / u")
    ax.set_ylabel("count")
    ax.set_title("(c) final ion mass spectrum")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


if __name__ == "__main__":
    raise SystemExit(main())
