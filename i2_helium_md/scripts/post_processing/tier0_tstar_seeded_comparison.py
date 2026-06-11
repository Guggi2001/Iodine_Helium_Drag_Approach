"""Tier-0 t*-seeded clean-form test: isolate the drag FORM from the transient.

The from-onset Tier-0 run (``tier0_drag_comparison.py``) conflates two errors:
the uncalibrated pre-t* transient (which puts MD off the reference trajectory by
window entry) and the drag-form question itself. This script removes the
transient by **seeding the ion at t* with the reference state** and integrating
the drag law forward with ``run_ion_propagation`` over ``[t*, t_end]`` -- the
strictly-cleaner form isolation hinted at in DRAG_PORT_DESIGN_DECISIONS §6.4.

Seed convention (internal-consistency, |v| placed radially)
-----------------------------------------------------------
The reference CSV now carries the full 3D per-atom velocity AND per-atom
positions, so the true I-I axis and the real radial/transverse split are
recoverable. The seed nonetheless places each atom's *full speed* ``|v_i|``
along the I-I axis (transverse zero), because the drag law was extracted with
the scalar speed ``v = |v2|`` as the radial velocity; forward-integrating gamma
from ``|v|`` is the honest Tier-0 internal-consistency check (does the BAOAB
driver reproduce the speed gamma was fit to?). This retires the previous
reconstruction, which differenced raw ``R(t)`` for the radial rate and dumped
the leftover speed into a fictitious transverse component -- inert KE the
central-force MD damped on startup, producing the spurious 9 A downcurve.

The harness additionally *reports* (does not seed) the real radial/transverse
split of each atom at t*, projected onto the true bond axis
``R_hat = (r1 - r2)/|r1 - r2|`` read from the position columns. For 9 A this
shows atom 2 is ~99.6% transverse at t* (it genuinely co-translates); for
18 A it is ~99.9% radial.

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
CASE = "18A"  # "9A" or "18A"
# A from-onset run for the same case, used only to read the droplet radius the
# real pipeline assigns to this droplet size (so the spatial gate matches).
ONSET_RUN_DIR = PROJECT_ROOT / "data" / "runs" / "18A_drag_tier0_N50"
SHOW_FIGURE = True
SHOW_POS_FIGURE = False
# Complementary diagnostic: time evolution of the three applied forces
# (drag, droplet-confining, Coulomb) projected onto the radial axis, for the
# clean extraction atom (atom 2). Reconstructed post-hoc from the stored
# trajectory -- the checkpoint does not store forces.
SHOW_FORCE_FIGURE = True


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
from i2_helium_md.simulation.ion import run_ion_propagation # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from i2_helium_md.postprocess.hedft_loader import load_hedft_trajectory  # noqa: E402


_BUILDERS = {"9A": single_pulse_N2000_drag, "18A": single_pulse_N2000_18Angst_drag}
_POST_WINDOW_MARGIN_PS = 25.0


def _interp(t, ref_t, ref_y):
    return float(np.interp(t, ref_t, ref_y))


def _report_real_radial_split(hedft, t_star):
    """Print the true radial/transverse split of each atom at t* (diagnostic).

    Reads the real 3D per-atom positions and velocities, builds the true bond
    axis ``R_hat = (r1 - r2)/|r1 - r2|``, and reports ``v . R_hat`` (radial)
    versus the transverse remainder for each atom. Also cross-checks that the
    position-derived separation matches the stored ``R_distance`` column. This
    feeds the verdict (is atom 2 genuinely near-radial?) but NOT the seed.
    """
    def vec(comps):
        return np.array([_interp(t_star, hedft.time_ps, c) for c in comps])

    r1 = vec((hedft.x1_A, hedft.y1_A, hedft.z1_A))
    r2 = vec((hedft.x2_A, hedft.y2_A, hedft.z2_A))
    v1 = vec((hedft.v1_x_Aps, hedft.v1_y_Aps, hedft.v1_z_Aps))
    v2 = vec((hedft.v2_x_Aps, hedft.v2_y_Aps, hedft.v2_z_Aps))

    r_rel = r1 - r2
    R_pos = float(np.linalg.norm(r_rel))
    R_col = _interp(t_star, hedft.time_ps, hedft.distance_A)
    # Provenance cross-check (TASK section 5): |r1 - r2| must equal R_distance.
    if not np.isclose(R_pos, R_col, rtol=1e-6, atol=1e-4):
        raise ValueError(
            f"position/separation mismatch at t*={t_star}: "
            f"|r1-r2|={R_pos:.6f} != R_distance={R_col:.6f}"
        )
    R_hat = r_rel / R_pos

    def split(v):
        rad = float(v @ R_hat)
        perp = float(np.linalg.norm(v - rad * R_hat))
        return rad, perp

    rad1, perp1 = split(v1)
    rad2, perp2 = split(v2)
    print(f"  real split @ t*={t_star} (true R_hat from positions, "
          f"|r1-r2|={R_pos:.3f}==R_distance):")
    print(f"    atom1: |v1|={np.linalg.norm(v1):.3f}  "
          f"radial={rad1:+.3f}  transverse={perp1:.3f}")
    print(f"    atom2: |v2|={np.linalg.norm(v2):.3f}  "
          f"radial={rad2:+.3f}  transverse={perp2:.3f}")


def seed_neutral_at_reference(hedft, t_star, *, m_eff_amu, droplet_radius):
    """Build a 1-molecule NeutralCheckpoint encoding the reference state at t*.

    Seeds the scored rotation-invariants: R(t*) and the per-atom speeds
    |v1|,|v2| placed *radially* along the constructed I-I axis (z), transverse
    zero. This is the internal-consistency convention -- the drag law was
    extracted with v = |v| as the radial velocity, so forward-integrating gamma
    from |v| tests whether the BAOAB driver reproduces the fitted speed. See the
    module docstring; the real radial/transverse split is reported separately by
    :func:`_report_real_radial_split` for the verdict, not used here.
    """
    R = _interp(t_star, hedft.time_ps, hedft.distance_A)
    v1 = _interp(t_star, hedft.time_ps, hedft.v1_magnitude_Aps)
    v2 = _interp(t_star, hedft.time_ps, hedft.v2_magnitude_Aps)

    _report_real_radial_split(hedft, t_star)
    print(f"  seed @ t*={t_star}: R={R:.3f}, "
          f"|v1|={v1:.3f} (radial), |v2|={v2:.3f} (radial), transverse=0")

    two_N = 2
    x = np.array([0.0, 0.0])
    y = np.array([0.0, 0.0])
    z = np.array([+R / 2.0, -R / 2.0])           # I-I axis = z
    vx = np.array([0.0, 0.0])                     # no transverse component
    vy = np.array([0.0, 0.0])
    vz = np.array([+v1, -v2])                     # full speed, outward (radial)

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


def run_case(case: str, onset_run_dir: Path, *, overrides=None):
    """Seed the ion at t* and forward-integrate the drag law over [t*, t_end].

    ``overrides`` (optional ``dict``) is forwarded to the preset builder as
    extra ``SimConfig`` field overrides -- e.g.
    ``{"binding_energy_I_ion_eV": 0.6, "potential_steepness": 8.0}`` to study the
    confining-potential influence (see ``tier0_confining_sweep.py``). The seed
    itself is parameter-independent (it encodes only the reference state and the
    fixed droplet radius), so every override starts from the same t* state.
    """
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
                  dt_ion=0.01, seed=1, **(overrides or {}))
    cfg.validate()
    ion = run_ion_propagation(cfg, neutral, verbose=False)

    # Shift the ion clock so t=0 -> physical t*.
    ion_shifted = replace(ion, time_ps=ion.time_ps + t_star)

    metrics = T.score(ion_shifted, hedft, (t_star, t_end))
    T.print_summary(metrics, label=f"{case} t*-seeded")
    return ion_shifted, hedft, (t_star, t_end), droplet_radius, cfg, metrics




def main() -> int:
    ion, hedft, window, drop_radius, cfg, _ = run_case(CASE, ONSET_RUN_DIR)
    if SHOW_FIGURE:
        import matplotlib.pyplot as plt
        T.build_figure(ion, hedft, window, SHOW_POS_FIGURE, drop_radius)
        plt.show()
    if SHOW_FORCE_FIGURE:
        import matplotlib.pyplot as plt
        # atom 2 (index 1) is the clean extraction atom; atom 1 (index 0) is its
        # partner. Both are seeded outward-radial, so the same sign conventions
        # apply (Coulomb +, droplet -, drag opposing outward motion).
        T.build_force_figure(ion, cfg, window, atom_index=1)
        T.build_force_figure(ion, cfg, window, atom_index=0)
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
