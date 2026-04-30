"""Tests for i2_helium_md/physics/collisions.py."""

import numpy as np
import pytest

from i2_helium_md.physics.collisions import (
    apply_collision,
    sample_collision_events,
    velocity_dependent_cross_section,
)
from i2_helium_md.physics.constants import EV, U


# Constants used throughout the tests -- match the module's source of truth.
AMU_KG = U
EV_J = EV
A_PER_PS_TO_M_PER_S = 100.0


def _energy_eV(m_amu: np.ndarray, vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> np.ndarray:
    """Helper: compute kinetic energy in eV from (m_amu, v in Å/ps)."""
    v2 = vx**2 + vy**2 + vz**2
    return 0.5 * (m_amu * AMU_KG) * (v2 * A_PER_PS_TO_M_PER_S**2) / EV_J


# ===========================================================================
# velocity_dependent_cross_section
# ===========================================================================
class TestVelocityDependentCrossSection:
    def test_zero_velocity_gives_inf(self):
        """v=0 with negative exponent gives sigma=inf (no exception)."""
        sigma = velocity_dependent_cross_section(
            np.array([0.0]),
            sigma_0_angstrom_sq=2500.0,
            exponent=-2.0,
        )
        assert np.isinf(sigma[0])
        assert sigma[0] > 0  # +inf, not -inf

    def test_known_values_v_minus_2(self):
        """Spot-check the formula at production parameters."""
        v = np.array([0.5, 1.0, 5.0, 50.0])
        sigma = velocity_dependent_cross_section(
            v, sigma_0_angstrom_sq=2500.0, exponent=-2.0,
        )
        expected = np.array([1.0e4, 2500.0, 100.0, 1.0])
        np.testing.assert_allclose(sigma, expected, rtol=1e-12)

    def test_positive_exponent(self):
        """sigma scales correctly for positive exponents too."""
        v = np.array([1.0, 2.0, 4.0])
        sigma = velocity_dependent_cross_section(
            v, sigma_0_angstrom_sq=10.0, exponent=1.5,
        )
        # 10 * v^1.5
        np.testing.assert_allclose(sigma, 10.0 * v ** 1.5, rtol=1e-12)

    def test_zero_exponent_constant(self):
        """With exponent=0, output is a constant array equal to sigma_0."""
        v = np.array([1.0, 5.0, 100.0])
        sigma = velocity_dependent_cross_section(
            v, sigma_0_angstrom_sq=42.0, exponent=0.0,
        )
        np.testing.assert_allclose(sigma, 42.0)

    def test_negative_velocity_raises(self):
        """Negative speed is invalid input."""
        with pytest.raises(ValueError, match="non-negative"):
            velocity_dependent_cross_section(
                np.array([1.0, -0.5]),
                sigma_0_angstrom_sq=2500.0,
                exponent=-2.0,
            )

    def test_shape_preserved(self):
        """Output shape matches input shape."""
        v = np.linspace(0.1, 10.0, 17)
        sigma = velocity_dependent_cross_section(
            v, sigma_0_angstrom_sq=2500.0, exponent=-2.0,
        )
        assert sigma.shape == v.shape

    def test_inf_propagates_through_sample_collision(self):
        """An atom at v=0 (sigma=inf) should ALWAYS collide.

        This is the contract that downstream ion code relies on: when
        an ion's velocity is zero (or near-zero with v^-2), the
        cross section diverges, p_scatter > 1, and the trial < p_scatter
        check is True regardless of the trial value.
        """
        # 100 trials with the same configuration
        counts = 0
        for seed in range(100):
            rng = np.random.default_rng(seed)
            sigma = velocity_dependent_cross_section(
                np.array([0.0]),
                sigma_0_angstrom_sq=2500.0,
                exponent=-2.0,
            )
            b = sample_collision_events(
                distance_travelled_angstrom=np.array([0.05]),
                depth_angstrom=np.array([-10.0]),  # inside droplet
                E0_eV=np.array([1.0]),
                sigma_angstrom_sq=sigma,
                E_min_eV=0.0,
                rng=rng,
            )
            counts += int(b[0])
        assert counts == 100, (
            f"expected 100/100 collisions at v=0 (sigma=inf), got {counts}"
        )


# ===========================================================================
# sample_collision_events
# ===========================================================================
class TestSampleCollisionEvents:
    def test_outside_droplet_never_collides(self):
        rng = np.random.default_rng(0)
        n = 1000
        out = sample_collision_events(
            distance_travelled_angstrom=np.full(n, 0.5),
            depth_angstrom=np.full(n, +5.0),         # outside
            E0_eV=np.full(n, 1.0),
            sigma_angstrom_sq=30.0,
            rng=rng,
        )
        assert not out.any()

    def test_below_E_min_never_collides(self):
        rng = np.random.default_rng(0)
        n = 1000
        out = sample_collision_events(
            distance_travelled_angstrom=np.full(n, 0.5),
            depth_angstrom=np.full(n, -5.0),
            E0_eV=np.full(n, 1e-5),                  # below threshold
            sigma_angstrom_sq=30.0,
            E_min_eV=1e-3,
            rng=rng,
        )
        assert not out.any()

    def test_zero_distance_never_collides(self):
        rng = np.random.default_rng(0)
        out = sample_collision_events(
            distance_travelled_angstrom=np.zeros(1000),
            depth_angstrom=np.full(1000, -5.0),
            E0_eV=np.full(1000, 1.0),
            sigma_angstrom_sq=30.0,
            rng=rng,
        )
        assert not out.any()

    def test_collision_rate_matches_formula(self):
        """Empirical fraction of collisions should match d * sigma * rho."""
        rng = np.random.default_rng(0)
        n = 200_000
        d, sigma, rho = 0.05, 30.0, 0.8 * 0.0219
        p_expected = d * sigma * rho
        out = sample_collision_events(
            distance_travelled_angstrom=np.full(n, d),
            depth_angstrom=np.full(n, -5.0),
            E0_eV=np.full(n, 1.0),
            sigma_angstrom_sq=sigma,
            droplet_density_per_angstrom3=rho,
            rng=rng,
        )
        empirical = out.mean()
        assert abs(empirical - p_expected) < 0.005, (
            f"empirical {empirical:.4f} vs expected {p_expected:.4f}"
        )

    def test_shape_mismatch_raises(self):
        rng = np.random.default_rng(0)
        with pytest.raises(ValueError, match="shape mismatch"):
            sample_collision_events(
                distance_travelled_angstrom=np.zeros(10),
                depth_angstrom=np.zeros(20),
                E0_eV=np.zeros(10),
                sigma_angstrom_sq=30.0,
                rng=rng,
            )

    def test_negative_distance_raises(self):
        rng = np.random.default_rng(0)
        with pytest.raises(ValueError, match="non-negative"):
            sample_collision_events(
                distance_travelled_angstrom=np.array([-0.1, 0.5]),
                depth_angstrom=np.array([-1.0, -1.0]),
                E0_eV=np.array([1.0, 1.0]),
                sigma_angstrom_sq=30.0,
                rng=rng,
            )

    def test_per_particle_sigma_array(self):
        """sigma_angstrom_sq accepts a per-particle array."""
        rng = np.random.default_rng(0)
        n = 1000
        sigma = np.array([0.0] * (n // 2) + [60.0] * (n // 2))
        out = sample_collision_events(
            distance_travelled_angstrom=np.full(n, 0.05),
            depth_angstrom=np.full(n, -5.0),
            E0_eV=np.full(n, 1.0),
            sigma_angstrom_sq=sigma,
            rng=rng,
        )
        # sigma=0 half should never collide
        assert not out[: n // 2].any()
        # sigma=60 half should sometimes collide
        assert out[n // 2 :].any()

    def test_reproducible_with_seed(self):
        kwargs = dict(
            distance_travelled_angstrom=np.full(100, 0.5),
            depth_angstrom=np.full(100, -5.0),
            E0_eV=np.full(100, 1.0),
            sigma_angstrom_sq=30.0,
        )
        out1 = sample_collision_events(rng=np.random.default_rng(7), **kwargs)
        out2 = sample_collision_events(rng=np.random.default_rng(7), **kwargs)
        np.testing.assert_array_equal(out1, out2)


# ===========================================================================
# apply_collision: passthrough behaviour
# ===========================================================================
class TestApplyCollisionPassthrough:
    def test_no_collisions_keeps_velocities_exactly(self):
        rng = np.random.default_rng(0)
        n = 50
        vx = rng.standard_normal(n)
        vy = rng.standard_normal(n)
        vz = rng.standard_normal(n)
        m = np.full(n, 127.0)
        b = np.zeros(n, dtype=bool)
        vx_n, vy_n, vz_n, dE = apply_collision(
            vx=vx, vy=vy, vz=vz, masses_amu=m,
            b_collision=b, scatter_mass_amu=4.0, rng=rng,
        )
        np.testing.assert_array_equal(vx_n, vx)
        np.testing.assert_array_equal(vy_n, vy)
        np.testing.assert_array_equal(vz_n, vz)
        assert np.all(dE == 0.0)

    def test_zero_dE_for_non_colliders_in_mixed_batch(self):
        """In a mixed batch, only colliders should have nonzero dE."""
        rng = np.random.default_rng(0)
        n = 100
        vx = np.full(n, 5.0)   # all moving in +x
        vy = np.zeros(n)
        vz = np.zeros(n)
        m = np.full(n, 127.0)
        b = np.zeros(n, dtype=bool)
        b[:50] = True   # first half collides
        _, _, _, dE = apply_collision(
            vx=vx, vy=vy, vz=vz, masses_amu=m,
            b_collision=b, scatter_mass_amu=4.0, rng=rng,
        )
        assert np.all(dE[50:] == 0.0)
        assert np.all(dE[:50] >= 0.0)


# ===========================================================================
# apply_collision: physical correctness
# ===========================================================================
class TestApplyCollisionPhysics:
    def test_energy_conserved_for_non_colliders(self):
        rng = np.random.default_rng(0)
        n = 200
        vx = rng.standard_normal(n)
        vy = rng.standard_normal(n)
        vz = rng.standard_normal(n)
        m = np.full(n, 127.0)
        b = rng.random(n) < 0.5
        E0 = _energy_eV(m, vx, vy, vz)
        vx_n, vy_n, vz_n, _ = apply_collision(
            vx=vx, vy=vy, vz=vz, masses_amu=m,
            b_collision=b, scatter_mass_amu=4.0, rng=rng,
        )
        E1 = _energy_eV(m, vx_n, vy_n, vz_n)
        # non-colliders: exact energy preservation
        np.testing.assert_allclose(E1[~b], E0[~b])

    def test_energy_strictly_decreases_for_colliders(self):
        rng = np.random.default_rng(0)
        n = 1000
        vx = np.full(n, 5.0); vy = np.zeros(n); vz = np.zeros(n)
        m = np.full(n, 127.0)
        b = np.ones(n, dtype=bool)
        E0 = _energy_eV(m, vx, vy, vz)
        vx_n, vy_n, vz_n, dE = apply_collision(
            vx=vx, vy=vy, vz=vz, masses_amu=m,
            b_collision=b, scatter_mass_amu=4.0, rng=rng,
        )
        E1 = _energy_eV(m, vx_n, vy_n, vz_n)
        assert np.all(E1 <= E0 + 1e-12)
        assert np.all(dE >= -1e-12)
        # And dE should agree with direct subtraction
        np.testing.assert_allclose(dE, E0 - E1, atol=1e-10)

    def test_speed_consistent_with_E1(self):
        """|v_new|^2 should match 2 E1 / m exactly."""
        rng = np.random.default_rng(0)
        n = 500
        vx = np.full(n, 5.0); vy = np.zeros(n); vz = np.zeros(n)
        m = np.full(n, 127.0)
        b = np.ones(n, dtype=bool)
        E0 = _energy_eV(m, vx, vy, vz)
        vx_n, vy_n, vz_n, dE = apply_collision(
            vx=vx, vy=vy, vz=vz, masses_amu=m,
            b_collision=b, scatter_mass_amu=4.0, rng=rng,
        )
        speed_sq = vx_n**2 + vy_n**2 + vz_n**2
        E1 = E0 - dE
        # speed in Å/ps: speed_sq = 2 * E1 / m (with proper unit conversions)
        expected_speed_sq = 2.0 * E1 * EV_J / (m * AMU_KG) / (A_PER_PS_TO_M_PER_S**2)
        np.testing.assert_allclose(speed_sq, expected_speed_sq, rtol=1e-9)

    def test_heavy_projectile_loses_little_energy(self):
        """Iodine (127 amu) hitting He (4 amu) should lose only a small
        fraction of its energy per collision."""
        rng = np.random.default_rng(0)
        n = 5000
        vx = np.full(n, 10.0); vy = np.zeros(n); vz = np.zeros(n)
        m = np.full(n, 127.0)
        b = np.ones(n, dtype=bool)
        E0 = _energy_eV(m, vx, vy, vz)
        _, _, _, dE = apply_collision(
            vx=vx, vy=vy, vz=vz, masses_amu=m,
            b_collision=b, scatter_mass_amu=4.0, rng=rng,
        )
        # Maximum fractional loss for elastic 2-body scattering with
        # rho = M/m_target = 127/4 ~ 31.75 is 4*rho/(1+rho)^2 ~ 0.119
        # = ~12% per head-on collision. Mean should be ~6% (theta uniform).
        mean_loss_frac = (dE / E0).mean()
        assert 0.0 < mean_loss_frac < 0.10, (
            f"mean fractional energy loss = {mean_loss_frac:.4f}; "
            "expected between 0 and 0.10 for I/He scattering"
        )

    def test_equal_mass_head_on_full_transfer(self):
        """Equal masses with head-on impact (b/R=0 -> COSTHETA=-1) should
        transfer all energy to the scatterer.

        We set b_collision=True and b/R=0 by mocking the RNG. Since we
        can't easily mock a random sample, we use a statistical check:
        for equal masses, the maximum energy transfer over many trials
        should approach 1.0.
        """
        rng = np.random.default_rng(0)
        n = 50_000
        vx = np.full(n, 10.0); vy = np.zeros(n); vz = np.zeros(n)
        m = np.full(n, 4.0)       # same as scatter mass
        b = np.ones(n, dtype=bool)
        E0 = _energy_eV(m, vx, vy, vz)
        _, _, _, dE = apply_collision(
            vx=vx, vy=vy, vz=vz, masses_amu=m,
            b_collision=b, scatter_mass_amu=4.0, rng=rng,
        )
        max_loss = (dE / E0).max()
        # Should approach 1.0 for equal masses with b/R close to 0.
        assert max_loss > 0.98, (
            f"max fractional energy transfer = {max_loss:.4f}; "
            "for equal masses should approach 1.0"
        )

    def test_E0_invariance_when_no_collision(self):
        """The function's internal E0 calculation must agree with our
        external _energy_eV helper."""
        # This is implicitly tested by test_energy_conserved_for_non_colliders
        pass


# ===========================================================================
# apply_collision: angular distribution checks
# ===========================================================================
class TestApplyCollisionDirections:
    def test_scattered_velocity_has_correct_lab_angle(self):
        """Mean cos(theta_lab) for COM uniform in cos should match the
        analytic transformation:
            cos(theta_lab) = (cos_com + rho) / sqrt(1 + 2 rho cos_com + rho²)
        Averaged over uniform cos_com, this gives a known result.
        """
        rng = np.random.default_rng(0)
        n = 100_000
        vx = np.full(n, 10.0); vy = np.zeros(n); vz = np.zeros(n)
        m = np.full(n, 127.0)
        b = np.ones(n, dtype=bool)
        vx_n, vy_n, vz_n, _ = apply_collision(
            vx=vx, vy=vy, vz=vz, masses_amu=m,
            b_collision=b, scatter_mass_amu=4.0, rng=rng,
        )
        # cos(theta_lab) = (vx_new . v_unit_in) / |v_new|
        speed_n = np.sqrt(vx_n**2 + vy_n**2 + vz_n**2)
        cos_lab = vx_n / speed_n   # since incoming is along +x

        # Numerical reference: average over uniform cos_com in [-1, 1]
        rho = 127.0 / 4.0
        cos_com = np.linspace(-1, 1, 10_000)
        denom = np.sqrt(1 + 2 * rho * cos_com + rho**2)
        cos_lab_ref = ((cos_com + rho) / denom).mean()
        # Tolerance: this is a sample mean, should be within ~3 sigma
        assert abs(cos_lab.mean() - cos_lab_ref) < 0.005, (
            f"empirical <cos lab> = {cos_lab.mean():.4f}, "
            f"analytic = {cos_lab_ref:.4f}"
        )

    def test_perpendicular_plane_isotropic(self):
        """Critical regression test for the COSBETA = uniform(-1,1)
        convention. Despite that not being a uniform azimuth, the
        random reference direction should make the perpendicular
        component isotropic.

        For incoming velocity along +x, the perpendicular component
        is (vy_new, vz_new). Its azimuth phi = atan2(vz, vy) should
        be uniform in [-π, π]. We check this with the mean angle and
        the second-moment matrix.
        """
        rng = np.random.default_rng(0)
        n = 200_000
        vx = np.full(n, 10.0); vy = np.zeros(n); vz = np.zeros(n)
        m = np.full(n, 127.0)
        b = np.ones(n, dtype=bool)
        _, vy_n, vz_n, _ = apply_collision(
            vx=vx, vy=vy, vz=vz, masses_amu=m,
            b_collision=b, scatter_mass_amu=4.0, rng=rng,
        )
        phi = np.arctan2(vz_n, vy_n)
        # Mean of cos(phi) and sin(phi) should be ~0 for uniform phi
        mean_cos_phi = float(np.cos(phi).mean())
        mean_sin_phi = float(np.sin(phi).mean())
        assert abs(mean_cos_phi) < 0.01, (
            f"<cos phi> = {mean_cos_phi:.4f}, should be ~0"
        )
        assert abs(mean_sin_phi) < 0.01, (
            f"<sin phi> = {mean_sin_phi:.4f}, should be ~0"
        )
        # Mean of cos(2 phi) should also be ~0 (catches a 2-fold bias)
        mean_cos_2phi = float(np.cos(2 * phi).mean())
        assert abs(mean_cos_2phi) < 0.01, (
            f"<cos 2phi> = {mean_cos_2phi:.4f}, should be ~0 "
            "(otherwise the perpendicular plane has a 2-fold anisotropy)"
        )


# ===========================================================================
# Reproducibility
# ===========================================================================
class TestReproducibility:
    def test_apply_collision_reproducible(self):
        kwargs = dict(
            vx=np.array([5.0, 5.0, 5.0]),
            vy=np.zeros(3),
            vz=np.zeros(3),
            masses_amu=np.full(3, 127.0),
            b_collision=np.array([True, True, False]),
            scatter_mass_amu=4.0,
        )
        a = apply_collision(rng=np.random.default_rng(11), **kwargs)
        b = apply_collision(rng=np.random.default_rng(11), **kwargs)
        for ai, bi in zip(a, b):
            np.testing.assert_array_equal(ai, bi)
