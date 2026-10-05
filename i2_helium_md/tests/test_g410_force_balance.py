"""THESIS_FINAL_GOALS G4.10: force-balance (Method A) drag cross-check.

No MD, no reference files.

1. the eV/A -> amu·A/ps^2 factor matches the leapfrog conversion;
2. a synthetic radial trajectory integrated with a known cubic drag, Coulomb
   pair force and solvation well is inverted back to its b (the end-to-end
   check of the balance, signs and units);
3. projection: velocity perpendicular to the pair axis gets no Coulomb
   drive under ``projected``, the full F_C under ``scalar``;
4. well sign (the 114c1ca amendment): outward motion -> the well term is
   positive and *lowers* the apparent drag; ``legacy_no_well`` drops it;
5. ``fit_drag_form`` recovers exact coefficients and fails loudly.
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from i2_helium_md.extraction.force_balance import (
    AMU_A_PER_PS2_PER_EV_PER_A,
    balance_geometry,
    fit_drag_form,
    force_balance_trace,
    savgol_speed_and_derivative,
)
from i2_helium_md.physics.constants import EV, EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2, U
from i2_helium_md.physics.drag import spatial_gate
from i2_helium_md.physics.potentials import droplet_force
from i2_helium_md.postprocess.hedft_loader import HedftTrajectory

M_AMU = 202.953908
B_TRUE = 2.2
E_BIND = 0.12
S_POT = 14.2
S_RHO = 3.5
R_DROP = 28.0
K_COULOMB_EV_A = 14.39964548


def _trajectory(t, x2, y2, z2, vx2, vy2, vz2) -> HedftTrajectory:
    """Symmetric pair: atom 1 mirrored through the origin."""
    speed = np.sqrt(vx2**2 + vy2**2 + vz2**2)
    dist = 2.0 * np.sqrt(x2**2 + y2**2 + z2**2)
    return HedftTrajectory(
        time_ps=t,
        v1_magnitude_Aps=speed, v2_magnitude_Aps=speed,
        v1_z_Aps=-vz2, v2_z_Aps=vz2, v1_x_Aps=-vx2, v2_x_Aps=vx2,
        v1_y_Aps=-vy2, v2_y_Aps=vy2,
        distance_A=dist,
        x1_A=-x2, x2_A=x2, y1_A=-y2, y2_A=y2, z1_A=-z2, z2_A=z2,
        droplet_radius_A=R_DROP,
        source_path=Path("synthetic"),
    )


def test_conversion_matches_leapfrog():
    # 1 eV/A on 1 amu: leapfrog gives F/m[kg] * EV*1e-4 A/ps^2.
    assert AMU_A_PER_PS2_PER_EV_PER_A == pytest.approx(EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2 / U)
    assert AMU_A_PER_PS2_PER_EV_PER_A == pytest.approx(EV * 1e-4 / U, rel=1e-15)
    assert AMU_A_PER_PS2_PER_EV_PER_A == pytest.approx(9648.533, rel=1e-6)


@pytest.fixture(scope="module")
def radial_synthetic():
    """Radial Coulomb explosion from r = 9 A with a known gated cubic drag."""

    def rhs(_t, y):
        r, v = y
        F_C = K_COULOMB_EV_A / (2.0 * r) ** 2 * AMU_A_PER_PS2_PER_EV_PER_A
        F_w = droplet_force(r - R_DROP, S_POT, E_BIND) * AMU_A_PER_PS2_PER_EV_PER_A
        F_d = B_TRUE * spatial_gate(r - R_DROP, S_RHO) * v**3
        return [v, (F_C - F_w - F_d) / M_AMU]

    t = np.arange(0.0, 8.0 + 5e-4, 1e-3)
    sol = solve_ivp(rhs, (t[0], t[-1]), [9.0, 0.0], t_eval=t, rtol=1e-11, atol=1e-12,
                    method="DOP853")
    r, v = sol.y
    zero = np.zeros_like(t)
    accel = np.array([rhs(0.0, [ri, vi])[1] for ri, vi in zip(r, v)])
    return t, r, v, accel, _trajectory(t, zero, zero, r, zero, zero, v)


def test_radial_inversion_recovers_b(radial_synthetic):
    t, r, v, accel, traj = radial_synthetic
    g = balance_geometry(traj)
    win = (t >= 1.0) & (t <= 7.0)
    kw = dict(mass_amu=M_AMU, binding_energy_eV=E_BIND, potential_steepness_A=S_POT,
              droplet_radius_A=R_DROP, gate_steepness_A=S_RHO)
    # Exact derivative: the balance itself inverts to the ODE's b (only the
    # ODE tolerance remains).
    fb_exact = force_balance_trace(g, t[win], v[win], accel[win], **kw)
    assert fb_exact.interior.sum() > 100
    np.testing.assert_allclose(fb_exact.cos_theta[fb_exact.interior], 1.0, atol=1e-12)
    fit_exact = fit_drag_form(
        fb_exact.speed_Aps[fb_exact.interior], fb_exact.F_drag_bulk[fb_exact.interior], "cubic"
    )
    assert fit_exact.b == pytest.approx(B_TRUE, rel=1e-6)
    # SG(order 2, 0.5 ps) derivative: the window's curvature bias is -0.36 %
    # here (measured: b 2.1920); 1 % bounds it, far inside the G4.10 band
    # half-width (19 %).
    speed, dvdt = savgol_speed_and_derivative(t[win], v[win], 0.5, smooth_speed=False)
    fb = force_balance_trace(g, t[win], speed, dvdt, **kw)
    fit = fit_drag_form(fb.speed_Aps[fb.interior], fb.F_drag_bulk[fb.interior], "cubic")
    assert fit.b == pytest.approx(B_TRUE, rel=1e-2)
    # Radial: scalar == projected.
    fb_s = force_balance_trace(
        g, t[win], speed, dvdt, mass_amu=M_AMU, binding_energy_eV=E_BIND,
        potential_steepness_A=S_POT, droplet_radius_A=R_DROP, gate_steepness_A=S_RHO,
        projection="scalar",
    )
    np.testing.assert_allclose(fb_s.F_drag_apparent, fb.F_drag_apparent, rtol=1e-12)


def _static_case(vx2, vz2):
    """Atom 2 at z = 10 A (inside a 28 A droplet), constant velocity."""
    t = np.arange(0.0, 2.0 + 5e-4, 1e-3)
    ones = np.ones_like(t)
    zero = np.zeros_like(t)
    traj = _trajectory(t, zero, zero, 10.0 * ones, vx2 * ones, zero, vz2 * ones)
    g = balance_geometry(traj)
    tw = t[500:1500]
    return g, tw, np.full(tw.shape, np.hypot(vx2, vz2)), np.zeros(tw.shape)


@pytest.mark.parametrize("projection", ["projected", "scalar", "legacy_no_well"])
def test_projection_perpendicular_motion(projection):
    g, tw, v, a = _static_case(vx2=3.0, vz2=0.0)
    fb = force_balance_trace(
        g, tw, v, a, mass_amu=M_AMU, binding_energy_eV=E_BIND, potential_steepness_A=S_POT,
        droplet_radius_A=R_DROP,
        gate_steepness_A=S_RHO, projection=projection,
    )
    F_C = K_COULOMB_EV_A / 20.0**2 * AMU_A_PER_PS2_PER_EV_PER_A
    np.testing.assert_allclose(fb.F_coulomb, F_C, rtol=1e-8)
    expected = 0.0 if projection == "projected" else F_C
    np.testing.assert_allclose(fb.F_coulomb_along, expected, atol=1e-9)
    # Tangential motion: the radial well does no work along v.
    np.testing.assert_allclose(fb.F_well_along, 0.0, atol=1e-9)


def test_well_sign_outward_motion_lowers_drag():
    g, tw, v, a = _static_case(vx2=0.0, vz2=3.0)
    kw = dict(mass_amu=M_AMU, binding_energy_eV=E_BIND, potential_steepness_A=S_POT,
              droplet_radius_A=R_DROP,
        gate_steepness_A=S_RHO)
    fb = force_balance_trace(g, tw, v, a, **kw)
    dU = droplet_force(10.0 - R_DROP, S_POT, E_BIND) * AMU_A_PER_PS2_PER_EV_PER_A
    assert dU > 0.0
    np.testing.assert_allclose(fb.F_well_along, dU, rtol=1e-12)
    np.testing.assert_allclose(fb.F_drag_apparent, fb.F_coulomb_along - dU, rtol=1e-12)
    legacy = force_balance_trace(g, tw, v, a, projection="legacy_no_well", **kw)
    np.testing.assert_allclose(legacy.F_drag_apparent - fb.F_drag_apparent, dU, rtol=1e-10)


def test_inertia_and_variable_mass():
    g, tw, v, _ = _static_case(vx2=0.0, vz2=3.0)
    a = np.full(tw.shape, -0.5)
    m = np.linspace(210.0, 180.0, tw.size)
    kw = dict(binding_energy_eV=E_BIND, potential_steepness_A=S_POT,
              droplet_radius_A=R_DROP,
        gate_steepness_A=S_RHO)
    fb = force_balance_trace(g, tw, v, a, mass_amu=m, **kw)
    np.testing.assert_allclose(fb.F_inertia, m * a)
    fb0 = force_balance_trace(g, tw, v, np.zeros_like(a), mass_amu=m, **kw)
    np.testing.assert_allclose(fb.F_drag_apparent - fb0.F_drag_apparent, -m * a)


def test_fit_drag_form_exact_and_validation():
    v = np.linspace(2.5, 5.0, 200)
    F = 3.0 * v + 1.5 * v**3
    lc = fit_drag_form(v, F, "linear_cubic")
    assert (lc.a, lc.b) == pytest.approx((3.0, 1.5), rel=1e-10)
    assert lc.rmse_force == pytest.approx(0.0, abs=1e-9)
    cub = fit_drag_form(v, 2.0 * v**3, "cubic")
    assert (cub.a, cub.b) == pytest.approx((0.0, 2.0), rel=1e-12)
    lin = fit_drag_form(v, 7.0 * v, "linear")
    assert (lin.a, lin.b) == pytest.approx((7.0, 0.0), rel=1e-12)
    F_nan = F.copy()
    F_nan[:10] = np.nan
    assert fit_drag_form(v, F_nan, "cubic").n_points == 190
    with pytest.raises(ValueError, match="form must be one of"):
        fit_drag_form(v, F, "quadratic")
    with pytest.raises(ValueError, match="finite samples"):
        fit_drag_form(v[:2], F[:2], "linear_cubic")


def test_trace_validation():
    g, tw, v, a = _static_case(vx2=0.0, vz2=3.0)
    kw = dict(mass_amu=M_AMU, binding_energy_eV=E_BIND, potential_steepness_A=S_POT,
              droplet_radius_A=R_DROP,
        gate_steepness_A=S_RHO)
    with pytest.raises(ValueError, match="projection must be one of"):
        force_balance_trace(g, tw, v, a, projection="cos", **kw)
    with pytest.raises(ValueError, match="outside the geometry"):
        force_balance_trace(g, tw + 5.0, v, a, **kw)
    with pytest.raises(ValueError, match="equal shape"):
        force_balance_trace(g, tw, v[:-1], a, **kw)
    t_bad = np.array([0.0, 0.1, 0.3, 0.4])
    with pytest.raises(ValueError, match="uniformly sampled"):
        savgol_speed_and_derivative(t_bad, np.ones(4), 0.2, smooth_speed=True)
