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
CASE = "9A"  # "9A" or "18A"
# A from-onset run for the same case, used only to read the droplet radius the
# real pipeline assigns to this droplet size (so the spatial gate matches).
ONSET_RUN_DIR = PROJECT_ROOT / "data" / "runs" / "9A_drag_tier0_N50"
SHOW_FIGURE = True
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
from i2_helium_md.physics.drag import drag_force  # noqa: E402
from i2_helium_md.physics.interactions import partner_interaction_ion  # noqa: E402
from i2_helium_md.physics.leapfrog import _droplet_acceleration  # noqa: E402
from i2_helium_md.presets import (  # noqa: E402
    single_pulse_N2000_drag,
    single_pulse_N2000_18Angst_drag,
)
from i2_helium_md.simulation.checkpoint import NeutralCheckpoint  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation, _drag_gate_steepness  # noqa: E402
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


def _reconstruct_radial_forces(ion, cfg, *, atom_index):
    """Reconstruct the three radial-projected forces for one atom over time.

    The ``IonCheckpoint`` stores no forces, but they are deterministic functions
    of the stored state + ``cfg``, so they are replayed here with the *same*
    physics functions the ion driver uses (CLAUDE.md rule 1 -- no duplicate
    physics). Each force is projected onto the radial unit vector ``r_hat``;
    sign convention: outward ``+``, inward ``-``.

    Returns ``(drag_radial, coulomb_radial, droplet_radial)`` arrays of shape
    ``(num_stored_steps,)`` in amu*A/ps^2, for the requested ``atom_index`` in
    the 2N layout.
    """
    px, py, pz = ion.positions_x, ion.positions_y, ion.positions_z
    vx, vy, vz = ion.velocities_x, ion.velocities_y, ion.velocities_z
    mass_kg = ion.mass_kg                       # (2N,)
    mass_amu = mass_kg / U                       # (2N,)
    droplet_radii = ion.droplet_radii_angstrom   # (2N,)
    two_N = px.shape[0]
    charge = np.ones(two_N, dtype=float)
    steepness = _drag_gate_steepness(cfg)

    n_steps = px.shape[1]
    drag_r = np.empty(n_steps)
    coul_r = np.empty(n_steps)
    drop_r = np.empty(n_steps)

    for j in range(n_steps):
        x, y, z = px[:, j], py[:, j], pz[:, j]
        ux, uy, uz = vx[:, j], vy[:, j], vz[:, j]

        r = np.sqrt(x**2 + y**2 + z**2)
        r_safe = np.where(r > 0, r, 1.0)
        rhx, rhy, rhz = x / r_safe, y / r_safe, z / r_safe

        depth = r - droplet_radii
        speed = np.sqrt(ux**2 + uy**2 + uz**2)
        s_safe = np.where(speed > 0, speed, 1.0)
        vhx, vhy, vhz = ux / s_safe, uy / s_safe, uz / s_safe

        # Drag: native force [amu*A/ps^2], magnitude form, direction -v_hat.
        f_drag = drag_force(speed, depth, cfg.drag_coefficients, steepness)
        fdx, fdy, fdz = -f_drag * vhx, -f_drag * vhy, -f_drag * vhz
        drag_radial = fdx * rhx + fdy * rhy + fdz * rhz

        # Droplet confining: accel [A/ps^2] -> force via * mass_amu.
        adx, ady, adz = _droplet_acceleration(
            x, y, z, mass_kg, droplet_radii, cfg, use_ion_binding=True,
        )
        drop_radial = (adx * rhx + ady * rhy + adz * rhz) * mass_amu

        # Coulomb partner: accel [A/ps^2] -> force via * mass_amu.
        acx, acy, acz, _ = partner_interaction_ion(
            x, y, z, mass_kg, charge, cfg,
        )
        coul_radial = (acx * rhx + acy * rhy + acz * rhz) * mass_amu

        drag_r[j] = drag_radial[atom_index]
        coul_r[j] = coul_radial[atom_index]
        drop_r[j] = drop_radial[atom_index]

    return drag_r, coul_r, drop_r


def build_force_figure(ion, cfg, window, *, atom_index=1):
    """Plot the radial-projected force evolution for the clean atom (atom 2).

    Complements :func:`tier0_drag_comparison.build_figure`: it shows *how the
    forces balance along the radial axis* over time -- drag, Coulomb, and the
    droplet-confining force, each projected onto ``r_hat`` (outward ``+``). The
    drag trace under-resolves the true dissipation for the non-radial 9 A case,
    the visual counterpart of the model-dimensionality residual.

    ``atom_index`` defaults to 1 = atom 2 (single-molecule 2N layout), the clean
    extraction atom the drag law was fit to.
    """
    import matplotlib.pyplot as plt

    drag_r, coul_r, drop_r = _reconstruct_radial_forces(
        ion, cfg, atom_index=atom_index,
    )
    t = ion.time_ps
    t_start, t_end = window

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.plot(t, coul_r, color="tab:red", label="Coulomb")
    ax.plot(t, drop_r, color="tab:blue", label="droplet (confining)")
    ax.plot(t, drag_r, color="tab:green", label="drag")
    ax.axvspan(t_start, t_end, color="tab:green", alpha=0.12,
               label="scored window")
    ax.set_xlabel("t / ps")
    ax.set_ylabel(r"$F \cdot \hat{r}$ / $\mathrm{amu}\,\mathrm{\AA}/\mathrm{ps}^2$")
    ax.set_title(
        f"{CASE} t*-seeded -- radial-projected forces, atom {atom_index + 1}"
    )
    ax.legend(frameon=False, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def main() -> int:
    ion, hedft, window, drop_radius, cfg, _ = run_case(CASE, ONSET_RUN_DIR)
    if SHOW_FIGURE:
        import matplotlib.pyplot as plt
        T.build_figure(ion, hedft, window, drop_radius)
        plt.show()
    if SHOW_FORCE_FIGURE:
        import matplotlib.pyplot as plt
        # atom 2 (index 1) is the clean extraction atom; atom 1 (index 0) is its
        # partner. Both are seeded outward-radial, so the same sign conventions
        # apply (Coulomb +, droplet -, drag opposing outward motion).
        build_force_figure(ion, cfg, window, atom_index=1)
        build_force_figure(ion, cfg, window, atom_index=0)
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
