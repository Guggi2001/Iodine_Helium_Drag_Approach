"""Tier-0 t*-seeded clean-form test: isolate the drag FORM from the transient.

The from-onset Tier-0 run (``tier0_drag_comparison.py``) conflates two errors:
the uncalibrated pre-t* transient (which puts MD off the reference trajectory by
window entry) and the drag-form question itself. This script removes the
transient by **seeding the ion at t* with the reference state** and integrating
the drag law forward with ``run_ion_propagation`` over ``[t*, t_end]`` -- the
strictly-cleaner form isolation hinted at in DRAG_PORT_DESIGN_DECISIONS §6.4.

Data limitation (important, recorded)
-------------------------------------
The HeDFT reference CSV stores per-atom SPEED magnitudes ``|v1|, |v2|`` and the
I-I separation ``R``, but only the x,z velocity components (no y) -- so the full
3D per-atom velocity is NOT recoverable from the file. We therefore seed the
rotation-invariants the comparison actually scores: ``R(t*)``, the radial
separation rate ``dR/dt(t*)`` (split equally between the equal-mass atoms along a
constructed I-I axis), and the per-atom speeds ``|v1|, |v2|`` (leftover speed
placed transverse). The raw ``R(t)`` carries bubble-mode oscillations (the modes
CEEMDAN removes before drag extraction), so ``dR/dt`` is taken as a central
difference; its exact value barely matters because the transverse speed
dominates ``|v|``.

How to use
----------
Edit USER SETTINGS, then::

    python scripts/post_processing/tier0_tstar_seeded_comparison.py
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# USER SETTINGS
# =============================================================================
CASE = "9A"  # "9A" or "18A"
# A from-onset run for the same case, used only to read the droplet radius the
# real pipeline assigns to this droplet size (so the spatial gate matches).
ONSET_RUN_DIR = PROJECT_ROOT / "data" / "runs" / "9A_drag_tier0_N50"
SHOW_FIGURE = True


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
# This script's own directory is on sys.path when run directly, so the sibling
# harness module imports cleanly.
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

import tier0_drag_comparison as T  # noqa: E402
from i2_helium_md.physics.constants import U  # noqa: E402
from i2_helium_md.presets import (  # noqa: E402
    single_pulse_N2000_drag,
    single_pulse_N2000_18Angst_drag,
)
from i2_helium_md.simulation.checkpoint import NeutralCheckpoint  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from i2_helium_md.postprocess.hedft_loader import load_hedft_trajectory  # noqa: E402


_BUILDERS = {"9A": single_pulse_N2000_drag, "18A": single_pulse_N2000_18Angst_drag}
_POST_WINDOW_MARGIN_PS = 4.0


def _interp(t, ref_t, ref_y):
    return float(np.interp(t, ref_t, ref_y))


def seed_neutral_at_reference(hedft, t_star, *, m_eff_amu, droplet_radius):
    """Build a 1-molecule NeutralCheckpoint encoding the reference state at t*.

    Seeds the scored rotation-invariants: R(t*), dR/dt(t*) (radial, along z),
    and per-atom speeds |v1|,|v2| (leftover placed transverse, along x). See the
    module docstring for why the full 3D velocity is not recoverable.
    """
    R = _interp(t_star, hedft.time_ps, hedft.distance_A)
    v1 = _interp(t_star, hedft.time_ps, hedft.v1_magnitude_Aps)
    v2 = _interp(t_star, hedft.time_ps, hedft.v2_magnitude_Aps)
    h = 0.1
    dRdt = (
        _interp(t_star + h, hedft.time_ps, hedft.distance_A)
        - _interp(t_star - h, hedft.time_ps, hedft.distance_A)
    ) / (2.0 * h)
    vz1, vz2 = +dRdt / 2.0, -dRdt / 2.0
    vperp1 = float(np.sqrt(max(0.0, v1 * v1 - vz1 * vz1)))
    vperp2 = float(np.sqrt(max(0.0, v2 * v2 - vz2 * vz2)))
    print(f"  reference @ t*={t_star}: R={R:.3f}, dR/dt={dRdt:.3f}, "
          f"|v1|={v1:.3f} (rad{vz1:.2f},perp{vperp1:.2f}), "
          f"|v2|={v2:.3f} (rad{vz2:.2f},perp{vperp2:.2f})")

    two_N = 2
    x = np.array([0.0, 0.0])
    y = np.array([0.0, 0.0])
    z = np.array([+R / 2.0, -R / 2.0])           # I-I axis = z
    vx = np.array([+vperp1, -vperp2])            # transverse, COM-minimising signs
    vy = np.array([0.0, 0.0])
    vz = np.array([vz1, vz2])                     # radial separation rate

    def col(a):
        return a.reshape(two_N, 1)

    zeros = np.zeros((two_N, 1))
    return NeutralCheckpoint(
        num_molecules=1,
        time_ps=np.zeros(1),
        positions_x=col(x), positions_y=col(y), positions_z=col(z),
        velocities_x=col(vx), velocities_y=col(vy), velocities_z=col(vz),
        mass_kg=np.full(two_N, m_eff_amu * U),
        droplet_radii=np.full(two_N, droplet_radius),
        r0=np.zeros(1),
        E_kin_eV=zeros.copy(), E_pot_eV=zeros.copy(),
        E_initial_eV=np.zeros(1),
        E_dissip_eV=zeros.copy(), L_droplet_eV_ps=zeros.copy(),
    )


def run_case(case: str, onset_run_dir: Path):
    builder = _BUILDERS[case]
    ref = PROJECT_ROOT / "data" / "reference" / f"{case}_All_Data.csv"
    coeff = PROJECT_ROOT / "data" / "reference" / "drag" / case / "linear_and_cubic"
    hedft = load_hedft_trajectory(ref)
    t_star, t_end, meff = T.read_drag_window(coeff)

    droplet_radius = float(
        RunDirectory(onset_run_dir).load_ion().droplet_radii_angstrom[0]
    )
    print(f"[{case}] droplet_radius={droplet_radius:.3f} A, "
          f"window=[{t_star},{t_end}], m_eff={meff:.4f}")

    neutral = seed_neutral_at_reference(
        hedft, t_star, m_eff_amu=meff, droplet_radius=droplet_radius
    )

    duration = (t_end - t_star) + _POST_WINDOW_MARGIN_PS
    cfg = builder(num_molecules=1, ion_simulation_time=duration,
                  dt_ion=0.01, seed=1)
    cfg.validate()
    ion = run_ion_propagation(cfg, neutral, verbose=False)

    # Shift the ion clock so t=0 -> physical t*.
    ion_shifted = replace(ion, time_ps=ion.time_ps + t_star)

    metrics = T.score(ion_shifted, hedft, (t_star, t_end))
    T.print_summary(metrics, label=f"{case} t*-seeded")
    return ion_shifted, hedft, (t_star, t_end), metrics


def main() -> int:
    ion, hedft, window, _ = run_case(CASE, ONSET_RUN_DIR)
    if SHOW_FIGURE:
        import matplotlib.pyplot as plt
        T.build_figure(ion, hedft, window)
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
