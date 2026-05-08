"""Hard-sphere collision physics for atom-helium scattering.

This module implements **Mode 3** of the legacy MATLAB collision logic
(`vmi_sim_3d_neutral_propa_HeDFT_mimic.m`, lines ~595-820), based on
*Andreas Braun's PhD thesis section D.2*. The two public functions are:

* :func:`sample_collision_events`  -- given per-particle distance
  travelled this timestep, depth into the droplet, and current kinetic
  energy, decide stochastically which particles undergo a collision.
* :func:`apply_collision`  -- given which particles collide, sample
  the impact parameter, compute new kinetic energy and lab-frame
  scattering angle, and produce updated 3D velocity vectors.

The module is **pure physics**: it has no dependency on ``SimConfig``,
no I/O, and no plotting. The simulation driver pulls config values
and passes them in as keyword arguments. This keeps the collision
physics easily testable in isolation and reusable for both neutral
and ion stages.

Modes 1 (probability per step) and 2 (mean-free-path tracking) from
the legacy code are deliberately **NOT** ported here. They are
considered legacy code paths and will be added if needed.

Important note on the azimuthal-smearing convention
---------------------------------------------------
The scattered velocity has its azimuth in the plane perpendicular to
the incoming velocity sampled with ``COSBETA = uniform(-1, 1)``. This
is **not** a uniform azimuth -- a true uniform azimuth would use
``phi = uniform(0, 2*pi); COSBETA = cos(phi)``. We mirror the MATLAB
exactly because the MATLAB is the project reference, but the
non-uniformity may be a latent bug in the legacy code.

In practice, the **reference direction** that defines the orientation
of the perpendicular plane is itself randomized per particle per step
(via ``reference_direction = rand(...) - 0.5``), and the bias may wash
out as the perpendicular axes rotate around the incoming velocity.
A unit test in ``test_collisions.py`` verifies that the resulting
scattered velocity directions ARE isotropic in the perpendicular
plane to within statistical tolerance. If you want to compare against
a uniform-phi alternative, see the docstring of
:func:`apply_collision`.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from .constants import DENSITY_DROPLET, EV, U


class CollisionDiagnostics(NamedTuple):
    """Per-step internals exposed by ``apply_collision`` when
    ``return_diagnostics=True``.

    Attributes
    ----------
    b_collision : np.ndarray, shape (n,)
        Boolean mask: True where the atom collided this step.
    COSTHETA : np.ndarray, shape (n,)
        Centre-of-mass scattering ``cos(theta)``. ``1.0`` for non-colliders.
    COStheta_lab : np.ndarray, shape (n,)
        **Lab-frame** scattering ``cos(theta_lab)`` after the optional
        Gaussian smearing has been applied. ``1.0`` for non-colliders.
        For heavy projectile on light target (e.g. I+ on He, rho ~ 32)
        this is always close to 1, with max scattering angle ~ asin(1/rho)
        ~ 1.81 deg. The legacy MATLAB temperature-diagnostic third column
        is the mean of ``arccos`` of this quantity (``vmi_sim_3d_ion_propa.m:683``).
    rho : np.ndarray, shape (n,)
        Mass ratio ``m_atom / m_scatterer``.
    E0_eV : np.ndarray, shape (n,)
        Pre-collision kinetic energy per atom in eV.
    E1_eV : np.ndarray, shape (n,)
        Post-collision kinetic energy per atom in eV.

    The legacy MATLAB temperature diagnostic
    (``vmi_sim_3d_ion_propa.m:683``) is built from
    ``mean(E1[b]/E0[b])``, ``mean((1+rho[b]**2)/(1+rho[b])**2)`` and
    ``mean(arccos(clip(COStheta_lab[b], -1, 1)))``.
    """
    b_collision: np.ndarray
    COSTHETA: np.ndarray
    COStheta_lab: np.ndarray
    rho: np.ndarray
    E0_eV: np.ndarray
    E1_eV: np.ndarray


def temperature_diagnostic_from_collision(
    diag: CollisionDiagnostics,
) -> np.ndarray:
    """Compute the legacy MATLAB temperature diagnostic for one step.

    Mirrors the per-step accumulator at ``vmi_sim_3d_ion_propa.m:683``::

        diagnostic_values = [
            mean(E1[b]/E0[b]),
            mean((1 + rho[b]**2) / (1 + rho[b])**2),
            mean(theta_lab[b]),
        ]

    where ``theta_lab`` is the **lab-frame** scattering angle in radians
    (post-Gaussian-smearing if smearing is enabled) and the means are
    taken only over atoms with ``b_collision == True``. The MATLAB
    code's ``theta`` variable on line 561 is the lab-frame angle, not
    the COM-frame angle -- a heavy ion on light helium target has a
    very narrow lab-frame scattering cone (max ~ asin(1/rho) which is
    ~ 1.81 deg for I+ on He).

    Parameters
    ----------
    diag : CollisionDiagnostics
        Output of ``apply_collision(..., return_diagnostics=True)``.

    Returns
    -------
    np.ndarray, shape (3,)
        ``[<T'/T>_actual, <T'/T>_from_mass_ratio, <theta_lab>_rad]``.
        All ``np.nan`` if no atom collided.
    """
    b = np.asarray(diag.b_collision, dtype=bool)
    if not b.any():
        return np.full(3, np.nan, dtype=float)

    E0 = np.asarray(diag.E0_eV, dtype=float)[b]
    E1 = np.asarray(diag.E1_eV, dtype=float)[b]
    rho = np.asarray(diag.rho, dtype=float)[b]
    cos_theta_lab = np.clip(
        np.asarray(diag.COStheta_lab, dtype=float)[b], -1.0, 1.0,
    )

    # Guard E0 = 0 (would only happen below the Landau cutoff which the
    # sampler already screens out, but be defensive).
    safe_E0 = np.where(E0 > 0.0, E0, 1.0)
    t_ratio_actual = float(np.mean(np.where(E0 > 0.0, E1 / safe_E0, 1.0)))
    t_ratio_mass = float(np.mean((1.0 + rho ** 2) / (1.0 + rho) ** 2))
    theta_mean = float(np.mean(np.arccos(cos_theta_lab)))

    return np.array(
        [t_ratio_actual, t_ratio_mass, theta_mean],
        dtype=float,
    )


# ===========================================================================
# Public: velocity-dependent cross section helper
# ===========================================================================
def velocity_dependent_cross_section(
    v_angstrom_per_ps: np.ndarray,
    *,
    sigma_0_angstrom_sq: float,
    exponent: float,
) -> np.ndarray:
    """Compute per-particle cross section for the v-dependent ion model.

    The legacy MATLAB code (lines ~440-466 of
    ``vmi_sim_3d_ion_propa.m``) uses

        sigma = sigma_0 * v ** exponent

    where ``v`` is the per-particle speed in Å/ps. With the production
    setting ``exponent = -2``, slow ions get arbitrarily large cross
    sections; combined with the Landau cutoff in
    :func:`sample_collision_events` (``E0 < E_min`` blocks collisions),
    this gives the empirically-tuned behaviour of strongly-coupled
    slow ions and weakly-coupled fast ions.

    At ``v = 0`` exactly, this returns ``+inf`` (when ``exponent < 0``).
    That is mathematically reasonable: an infinitely-slow ion has an
    infinite mean collision time so any nonzero step gives p > 1 →
    "always collides". The downstream ``trial < p_scatter`` check in
    :func:`sample_collision_events` then evaluates True regardless of
    the trial value, so a ``+inf`` propagates cleanly.

    Parameters
    ----------
    v_angstrom_per_ps : np.ndarray
        Per-particle speed in Å/ps. Shape ``(n_particles,)``.
    sigma_0_angstrom_sq : float
        Reference cross section ``sigma_0`` in Å² (e.g.
        ``cfg.geometric_scattering_crosssection_Iplus``).
    exponent : float
        Power-law exponent (``cfg.sigma_ion_exponent``, typically -2).

    Returns
    -------
    np.ndarray
        Per-particle cross section in Å², same shape as ``v``.
    """
    v = np.asarray(v_angstrom_per_ps, dtype=float)
    if np.any(v < 0):
        raise ValueError("speed must be non-negative")

    # Numpy's float power handles 0**(-2) = inf with a warning.
    # We suppress the warning since inf is the intended result.
    with np.errstate(divide="ignore"):
        return sigma_0_angstrom_sq * v ** exponent


# ===========================================================================
# Public: collision-event sampling (Mode 3)
# ===========================================================================
def sample_collision_events(
    *,
    distance_travelled_angstrom: np.ndarray,
    depth_angstrom: np.ndarray,
    E0_eV: np.ndarray,
    sigma_angstrom_sq: float | np.ndarray,
    droplet_density_per_angstrom3: float = DENSITY_DROPLET,
    E_min_eV: float = 0.0,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Decide which particles collide this timestep (Mode 3).

    Mode 3 of the legacy MATLAB collision modes. Per-particle
    collision probability is

        p_scatter = distance_travelled * sigma * rho_droplet

    A collision is sampled iff:

    1. ``rand() < p_scatter``  (Poisson-thinning approximation)
    2. ``depth < 0``           (particle is inside the droplet)
    3. ``E0 >= E_min``         (kinetic energy above Landau threshold)

    All input arrays must have the same shape. The output has the
    same shape and dtype ``bool``.

    Parameters
    ----------
    distance_travelled_angstrom : np.ndarray
        Distance each particle moved this timestep, in Å.
    depth_angstrom : np.ndarray
        Distance from droplet surface in Å. Negative = inside droplet.
        Non-negative = outside; cannot collide.
    E0_eV : np.ndarray
        Current kinetic energy of each particle, in eV.
    sigma_angstrom_sq : float or np.ndarray
        Geometric scattering cross section in Å². Scalar uses the same
        value for all particles; array allows per-particle (e.g.
        velocity-dependent) cross sections.
    droplet_density_per_angstrom3 : float, optional
        Number density of helium atoms in the droplet, in atoms/Å³.
        Default ``DENSITY_DROPLET`` from ``physics.constants`` =
        ``0.8 * BULK_DENSITY_HELIUM`` ≈ 0.01752 / Å³.
    E_min_eV : float, optional
        Minimum kinetic energy below which collisions are forbidden
        (Landau cutoff). Default 0 (no threshold).
    rng : np.random.Generator, optional
        Reproducible RNG. If None, a fresh default RNG is constructed.

    Returns
    -------
    np.ndarray
        Boolean array, True where the particle collides this step.

    Raises
    ------
    ValueError
        If input arrays have inconsistent shapes or invalid values.
    """
    if rng is None:
        rng = np.random.default_rng()

    distance_travelled_angstrom = np.asarray(distance_travelled_angstrom, dtype=float)
    depth_angstrom = np.asarray(depth_angstrom, dtype=float)
    E0_eV = np.asarray(E0_eV, dtype=float)

    # Shape consistency
    if not (
        distance_travelled_angstrom.shape == depth_angstrom.shape == E0_eV.shape
    ):
        raise ValueError(
            "shape mismatch: "
            f"distance_travelled {distance_travelled_angstrom.shape}, "
            f"depth {depth_angstrom.shape}, E0 {E0_eV.shape}"
        )

    if np.any(distance_travelled_angstrom < 0):
        raise ValueError("distance_travelled_angstrom must be non-negative")

    if isinstance(sigma_angstrom_sq, np.ndarray):
        if sigma_angstrom_sq.shape != distance_travelled_angstrom.shape:
            raise ValueError(
                f"sigma shape {sigma_angstrom_sq.shape} does not match "
                f"distance shape {distance_travelled_angstrom.shape}"
            )

    # Probability of collision this step
    p_scatter = (
        distance_travelled_angstrom * sigma_angstrom_sq * droplet_density_per_angstrom3
    )
    trial = rng.uniform(0.0, 1.0, size=distance_travelled_angstrom.shape)

    inside_droplet = depth_angstrom < 0.0
    above_threshold = E0_eV >= E_min_eV

    return (trial < p_scatter) & inside_droplet & above_threshold


# ===========================================================================
# Public: apply scattering to colliders
# ===========================================================================
def apply_collision(
    *,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    masses_amu: np.ndarray,
    b_collision: np.ndarray,
    scatter_mass_amu: float,
    neutral_scatter_angle_std_deg: float = 0.0,
    rng: np.random.Generator | None = None,
    return_diagnostics: bool = False,
) -> (
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    | tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, CollisionDiagnostics]
):
    """Apply hard-sphere scattering to colliding particles.

    For each particle marked True in ``b_collision``:

    1. Sample impact parameter ``b/R`` via inverse-CDF (``b/R = sqrt(u)``).
    2. Compute centre-of-mass scattering angle ``COSTHETA = 2(b/R)² - 1``.
    3. Transform to lab-frame angle using mass ratio
       ``rho = m_particle / m_scatterer``.
    4. Optionally smear the lab angle by Gaussian noise.
    5. Build a random orthonormal basis in the plane perpendicular to
       the incoming velocity, sample a random azimuth in that plane,
       and assemble the new 3D velocity vector with the correct
       parallel/perpendicular split.

    Particles with ``b_collision = False`` retain their incoming
    velocity exactly; their reported energy loss is zero.

    Parameters
    ----------
    vx, vy, vz : np.ndarray
        Velocity components, all of shape ``(2N,)``, in Å/ps.
    masses_amu : np.ndarray
        Per-particle mass in atomic mass units. Shape ``(2N,)``.
    b_collision : np.ndarray
        Boolean mask of shape ``(2N,)``; True where the particle
        underwent a collision (typically from
        :func:`sample_collision_events`).
    scatter_mass_amu : float
        Mass of the scatterer (helium = 4.0 amu).
    neutral_scatter_angle_std_deg : float, optional
        Standard deviation in degrees of additional Gaussian smearing
        on the lab-frame scattering angle. Default 0 (no smearing).
        Applied only to colliding particles.
    rng : np.random.Generator, optional
        Reproducible RNG. If None, a fresh default RNG is constructed.
    return_diagnostics : bool, optional
        If True, return a 5-tuple appending a
        :class:`CollisionDiagnostics` namedtuple
        ``(b_collision, COSTHETA, rho, E0_eV, E1_eV)``. The
        ion-propagation driver uses this to record the legacy MATLAB
        temperature-diagnostic accumulator from
        ``vmi_sim_3d_ion_propa.m:683``. Default ``False`` keeps the
        4-tuple shape so existing neutral-side callers are unaffected.

    Returns
    -------
    vx_new, vy_new, vz_new : np.ndarray
        Updated velocity components, same shape as inputs.
    delta_E_eV : np.ndarray
        ``E0 - E1`` per particle in eV (positive = energy lost). Zero
        for non-colliding particles.
    diagnostics : CollisionDiagnostics, optional
        Only present when ``return_diagnostics=True``.

    Notes
    -----
    The scattering azimuth is sampled with ``COSBETA = uniform(-1, 1);
    SINBETA = sqrt(1 - COSBETA²)``, mirroring the legacy MATLAB code.
    This is **not** a uniform azimuth, but the orthogonal axes
    ``velocity_normal_1`` and ``velocity_normal_2`` are themselves
    rotated randomly by drawing a fresh random reference direction
    per particle per step. A unit test verifies that the resulting
    scattered velocity directions are isotropic in the perpendicular
    plane to within statistical tolerance. If you want to compare
    against a true uniform-phi alternative, replace lines forming
    ``COSBETA``/``SINBETA`` in the source with::

        phi = rng.uniform(0, 2*np.pi, size=...)
        COSBETA = np.cos(phi); SINBETA = np.sin(phi)

    Raises
    ------
    ValueError
        If input arrays have inconsistent shapes.
    """
    if rng is None:
        rng = np.random.default_rng()

    vx = np.asarray(vx, dtype=float)
    vy = np.asarray(vy, dtype=float)
    vz = np.asarray(vz, dtype=float)
    masses_amu = np.asarray(masses_amu, dtype=float)
    b_collision = np.asarray(b_collision, dtype=bool)

    if not (vx.shape == vy.shape == vz.shape == masses_amu.shape == b_collision.shape):
        raise ValueError(
            "shape mismatch: "
            f"vx {vx.shape}, vy {vy.shape}, vz {vz.shape}, "
            f"masses {masses_amu.shape}, b_collision {b_collision.shape}"
        )

    # ----- speed and unit velocity vector -----
    v_speed = np.sqrt(vx**2 + vy**2 + vz**2)
    # Avoid division by zero for stationary particles. They cannot
    # collide because depth-check or E0<E_min would have screened them,
    # but be defensive: we set unit vectors to (1,0,0) for stationary,
    # which is harmless because b_collision will be False there.
    safe = v_speed > 0.0
    v_unit_x = np.where(safe, vx / np.where(safe, v_speed, 1.0), 1.0)
    v_unit_y = np.where(safe, vy / np.where(safe, v_speed, 1.0), 0.0)
    v_unit_z = np.where(safe, vz / np.where(safe, v_speed, 1.0), 0.0)

    # ----- E0 in eV -----
    # E = (1/2) m v^2  with m in amu, v in Å/ps:
    # 1 amu = U kg
    # 1 Å/ps = 100 m/s
    # We use the project-level constants U and EV (= 1 eV in J) from
    # physics.constants so we have a single source of truth.
    A_PER_PS_TO_M_PER_S = 100.0
    E0_eV = 0.5 * (masses_amu * U) * (v_speed * A_PER_PS_TO_M_PER_S) ** 2 / EV

    # ----- Sample impact parameter b/R via inverse-CDF -----
    n = vx.shape[0]
    impact_parameter_norm = np.sqrt(rng.uniform(0.0, 1.0, size=n))   # b/R in [0,1]

    # COM-frame scattering cos(theta) from b/R:
    #   for hard sphere: cos(theta_com) = 2 (b/R)^2 - 1
    COSTHETA = 2.0 * impact_parameter_norm**2 - 1.0
    SINTHETA = np.sqrt(np.clip(1.0 - COSTHETA**2, 0.0, 1.0))

    # Non-colliders pass through unchanged in COM frame
    COSTHETA = np.where(b_collision, COSTHETA, 1.0)
    SINTHETA = np.where(b_collision, SINTHETA, 0.0)

    # ----- mass ratio and post-collision energy -----
    rho = masses_amu / scatter_mass_amu
    # E1 = E0 * (1 + 2*rho*cos + rho^2) / (1+rho)^2
    E1_eV = E0_eV * (1.0 + 2.0 * rho * COSTHETA + rho**2) / (1.0 + rho) ** 2
    delta_E_eV = E0_eV - E1_eV
    # Force exact zero where no collision (avoid roundoff noise)
    delta_E_eV = np.where(b_collision, delta_E_eV, 0.0)

    # ----- COM -> lab frame scattering angle -----
    denom = np.sqrt(np.maximum(1.0 + 2.0 * rho * COSTHETA + rho**2, 0.0))
    # Avoid division by zero: denom is zero only if rho=1 and COSTHETA=-1,
    # i.e. equal masses with head-on impact. We guard with where.
    safe_denom = denom > 0.0
    COStheta_lab = np.where(
        safe_denom, (COSTHETA + rho) / np.where(safe_denom, denom, 1.0), 1.0
    )
    SINtheta_lab = np.sqrt(np.clip(1.0 - COStheta_lab**2, 0.0, 1.0))

    # ----- optional Gaussian angular smearing (only colliders) -----
    if neutral_scatter_angle_std_deg > 0.0:
        std_rad = neutral_scatter_angle_std_deg * np.pi / 180.0
        n_coll = int(b_collision.sum())
        if n_coll > 0:
            theta_lab = np.arccos(np.clip(COStheta_lab[b_collision], -1.0, 1.0))
            theta_lab_smeared = theta_lab + rng.standard_normal(n_coll) * std_rad
            COStheta_lab[b_collision] = np.cos(theta_lab_smeared)
            SINtheta_lab[b_collision] = np.sqrt(
                np.clip(1.0 - COStheta_lab[b_collision] ** 2, 0.0, 1.0)
            )

    # ----- random orthonormal basis perpendicular to incoming velocity -----
    # 1) random reference direction (one fresh draw per particle per step)
    ref = rng.uniform(-0.5, 0.5, size=(n, 3))
    ref_norm = np.linalg.norm(ref, axis=1, keepdims=True)
    # Defensive: avoid divide-by-zero if a random reference happened to be 0
    ref_norm = np.where(ref_norm > 0, ref_norm, 1.0)
    ref = ref / ref_norm

    v_unit = np.column_stack([v_unit_x, v_unit_y, v_unit_z])

    # 2) velocity_normal_1 = normalise(v_unit x ref)
    n1 = np.cross(v_unit, ref)
    n1_norm = np.linalg.norm(n1, axis=1, keepdims=True)
    n1_norm = np.where(n1_norm > 0, n1_norm, 1.0)
    n1 = n1 / n1_norm

    # 3) velocity_normal_2 = v_unit x n1   (already unit length and orthogonal)
    n2 = np.cross(v_unit, n1)

    # ----- random azimuth (legacy MATLAB convention; see module docstring) -----
    COSBETA = (rng.uniform(0.0, 1.0, size=n) - 0.5) * 2.0   # uniform in [-1, 1]
    SINBETA = np.sqrt(np.clip(1.0 - COSBETA**2, 0.0, 1.0))

    # ----- assemble new velocity vector -----
    # v_new = sqrt(2 E1 / m) in Å/ps
    v_new_speed = np.sqrt(
        2.0 * E1_eV * EV / (masses_amu * U)
    ) / A_PER_PS_TO_M_PER_S

    v_parallel = v_new_speed * COStheta_lab
    v_perp = v_new_speed * SINtheta_lab

    new_v = (
        v_unit * v_parallel[:, None]
        + n1 * (COSBETA * v_perp)[:, None]
        + n2 * (SINBETA * v_perp)[:, None]
    )

    # Non-colliders retain their incoming velocity exactly. Note that
    # for non-colliders the construction above gives ~ v_unit * v_speed
    # (since cos_theta_lab=1, sin_theta_lab=0, and v_new_speed=v_speed
    # because rho-formula collapses for COSTHETA=1), but small roundoff
    # could accumulate; we restore exactly.
    nc = ~b_collision
    new_v[nc, 0] = vx[nc]
    new_v[nc, 1] = vy[nc]
    new_v[nc, 2] = vz[nc]

    vx_new = new_v[:, 0]
    vy_new = new_v[:, 1]
    vz_new = new_v[:, 2]

    if return_diagnostics:
        return (
            vx_new, vy_new, vz_new, delta_E_eV,
            CollisionDiagnostics(
                b_collision=b_collision,
                COSTHETA=COSTHETA,
                COStheta_lab=COStheta_lab,
                rho=rho,
                E0_eV=E0_eV,
                E1_eV=E1_eV,
            ),
        )
    return vx_new, vy_new, vz_new, delta_E_eV
