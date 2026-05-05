"""Tests for i2_helium_md/simulation/ion_initial_state.py."""

from dataclasses import replace

import numpy as np
import pytest

from i2_helium_md.presets import single_pulse_N2000
from i2_helium_md.simulation.checkpoint import (
    IonCheckpoint,
    _ION_SCHEMA_VERSION,
)
from i2_helium_md.simulation.ion_initial_state import build_initial_ion_state
from i2_helium_md.simulation.neutral import run_neutral_propagation


# ---------------------------------------------------------------------------
# Helper: build a "real" neutral checkpoint by running a tiny neutral
# propagation. Fixture-style so multiple tests can reuse it.
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def small_neutral_run():
    """Short neutral propagation: 5 molecules, 10 steps (= 0.1 ps)."""
    cfg = single_pulse_N2000(num_molecules=5, seed=42)
    cfg = replace(cfg, t_max_neutral=0.1, dt_neutral=0.01)
    neutral = run_neutral_propagation(cfg, verbose=False)
    return cfg, neutral


# ===========================================================================
# Basic API and shapes
# ===========================================================================
class TestApi:
    def test_returns_ion_checkpoint(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=20)
        assert isinstance(ion, IonCheckpoint)

    def test_schema_version_is_current(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=20)
        assert ion.schema_version == _ION_SCHEMA_VERSION

    def test_per_atom_traj_shapes_2N_T(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        N = cfg.num_molecules
        T = 31
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=T)
        for arr in (ion.positions_x, ion.positions_y, ion.positions_z,
                    ion.velocities_x, ion.velocities_y, ion.velocities_z,
                    ion.E_kin_eV, ion.E_pot_eV, ion.E_dissip_eV,
                    ion.relative_loss_per_ps, ion.number_of_collisions,
                    ion.mass_history_kg):
            assert arr.shape == (2 * N, T), f"got {arr.shape}"

    def test_per_atom_static_arrays(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        N = cfg.num_molecules
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=10)
        for arr in (ion.mass_kg, ion.mass_final_kg,
                    ion.droplet_radii_angstrom,
                    ion.positions_final_x, ion.positions_final_y,
                    ion.positions_final_z,
                    ion.velocities_final_x, ion.velocities_final_y,
                    ion.velocities_final_z):
            assert arr.shape == (2 * N,)

    def test_b_ion_outside_shape_N(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=10)
        assert ion.b_ion_outside.shape == (cfg.num_molecules,)
        assert ion.b_ion_outside.dtype == bool

    def test_time_ps_shape_T(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=42)
        assert ion.time_ps.shape == (42,)


# ===========================================================================
# Inheritance from neutral checkpoint
# ===========================================================================
class TestInheritance:
    def test_column_zero_is_neutral_last_column(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=10)
        # Default start_id = -1 -> last neutral column
        np.testing.assert_array_equal(ion.positions_x[:, 0], neutral.positions_x[:, -1])
        np.testing.assert_array_equal(ion.positions_y[:, 0], neutral.positions_y[:, -1])
        np.testing.assert_array_equal(ion.positions_z[:, 0], neutral.positions_z[:, -1])
        np.testing.assert_array_equal(ion.velocities_x[:, 0], neutral.velocities_x[:, -1])
        np.testing.assert_array_equal(ion.velocities_y[:, 0], neutral.velocities_y[:, -1])
        np.testing.assert_array_equal(ion.velocities_z[:, 0], neutral.velocities_z[:, -1])

    def test_droplet_radii_inherited(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=10)
        np.testing.assert_array_equal(
            ion.droplet_radii_angstrom, neutral.droplet_radii
        )

    def test_mass_kg_inherited_and_history_starts_with_it(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=10)
        np.testing.assert_array_equal(ion.mass_kg, neutral.mass_kg)
        np.testing.assert_array_equal(ion.mass_history_kg[:, 0], ion.mass_kg)
        np.testing.assert_array_equal(ion.mass_final_kg, ion.mass_kg)

    def test_explicit_start_id_zero_uses_neutral_t0(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=10, start_id=0)
        np.testing.assert_array_equal(ion.positions_x[:, 0], neutral.positions_x[:, 0])

    def test_invalid_start_id_raises(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        with pytest.raises(ValueError, match="start_id"):
            build_initial_ion_state(cfg, neutral, num_steps_ion=10, start_id=999)


# ===========================================================================
# Initial conditions: column 0 only populated, rest is zero
# ===========================================================================
class TestInitialColumnsAreEmpty:
    def test_position_columns_1plus_are_zero(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=20)
        for arr in (ion.positions_x, ion.positions_y, ion.positions_z,
                    ion.velocities_x, ion.velocities_y, ion.velocities_z,
                    ion.E_kin_eV, ion.E_pot_eV):
            assert np.all(arr[:, 1:] == 0)

    def test_dissip_n_coll_loss_all_zero(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=20)
        # E_dissip starts at 0 EVERYWHERE (including column 0).
        assert np.all(ion.E_dissip_eV == 0)
        assert np.all(ion.number_of_collisions == 0)
        assert np.all(ion.relative_loss_per_ps == 0)

    def test_b_ion_outside_starts_false(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=20)
        assert not ion.b_ion_outside.any()


def _deep_copy_neutral(neutral):
    """Return a NeutralCheckpoint with all numpy arrays deep-copied.

    ``dataclasses.replace`` only does a shallow copy: array fields are
    shared between the new and original instances. Tests that mutate
    arrays must deep-copy first or they corrupt the shared fixture.
    """
    from i2_helium_md.simulation.checkpoint import NeutralCheckpoint
    from dataclasses import fields
    kwargs = {}
    for f in fields(NeutralCheckpoint):
        v = getattr(neutral, f.name)
        if isinstance(v, np.ndarray):
            kwargs[f.name] = v.copy()
        else:
            kwargs[f.name] = v
    return NeutralCheckpoint(**kwargs)


# ===========================================================================
# Initial energies: the bug-fix tests
# ===========================================================================
class TestInitialEnergies:
    """Regression tests for the bugs we fixed vs. legacy MATLAB t=0 E_pot/E_kin."""

    def test_E_kin_uses_v_squared_not_v_to_the_fourth(self, small_neutral_run):
        """Legacy ion script line 289 has E = m * (vx²+vy²)² / 2 / eV
        which is wrong (squares v² to give v⁴). We use m*v²/2/eV."""
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=10)

        # Manual recomputation
        from i2_helium_md.physics.constants import EV
        vx = neutral.velocities_x[:, -1]
        vy = neutral.velocities_y[:, -1]
        vz = neutral.velocities_z[:, -1]
        m = neutral.mass_kg
        E_expected = 0.5 * m * (vx ** 2 + vy ** 2 + vz ** 2) * (100.0 ** 2) / EV
        np.testing.assert_allclose(ion.E_kin_eV[:, 0], E_expected, rtol=1e-12)

    def test_E_kin_includes_vz(self, small_neutral_run):
        """Legacy MATLAB omitted vz from the t=0 E_kin formula. Verify
        ours always includes it by constructing a state with vx=vy=0
        but vz != 0 -- a buggy E_kin would be zero, ours should be nonzero."""
        cfg, neutral = small_neutral_run

        # Deep-copy before mutating so we don't pollute the module-scoped fixture.
        nc_modified = _deep_copy_neutral(neutral)
        nc_modified.velocities_x[:, -1] = 0.0
        nc_modified.velocities_y[:, -1] = 0.0
        nc_modified.velocities_z[:, -1] = 3.0  # 3 Å/ps

        ion = build_initial_ion_state(cfg, nc_modified, num_steps_ion=10)
        # Expected: positive, nonzero E_kin per atom
        assert np.all(ion.E_kin_eV[:, 0] > 0)

    def test_E_pot_droplet_uses_3d_radius(self, small_neutral_run):
        """Legacy ion script line 291 used sqrt(x²+y²) instead of
        sqrt(x²+y²+z²). Verify ours uses the 3D radius."""
        cfg, neutral = small_neutral_run
        N = cfg.num_molecules

        # Deep-copy before mutating so we don't pollute the module-scoped fixture.
        nc_modified = _deep_copy_neutral(neutral)
        # Atom-1 of each molecule at z=+50 (outside droplet); atom-2 at z=-50
        # so the pair Coulomb energy is finite (not r=0 -> inf).
        nc_modified.positions_x[:, -1] = 0.0
        nc_modified.positions_y[:, -1] = 0.0
        nc_modified.positions_z[:N, -1] = 50.0
        nc_modified.positions_z[N:, -1] = -50.0
        nc_modified.velocities_x[:, -1] = 0.0
        nc_modified.velocities_y[:, -1] = 0.0
        nc_modified.velocities_z[:, -1] = 0.0

        ion = build_initial_ion_state(cfg, nc_modified, num_steps_ion=10)
        # Atom outside droplet -> droplet term ≈ binding_energy_I_ion ≈ 0.3 eV
        # Pair separation = 100 Å -> Coulomb pair = 14.4/100 = 0.144 eV
        # Half per atom = 0.072 eV
        # Total per atom: ≈ 0.372 eV
        for i in range(2 * N):
            assert 0.30 < ion.E_pot_eV[i, 0] < 0.40, (
                f"atom {i} E_pot = {ion.E_pot_eV[i, 0]}, expected ~0.37 eV"
            )

    def test_E_pot_includes_partner_coulomb(self, small_neutral_run):
        """Legacy ion script t=0 omits the partner Coulomb term entirely.
        Ours includes it. Verify by checking E_pot is dominated by Coulomb
        when atoms are deep inside the droplet at typical separation."""
        cfg, neutral = small_neutral_run
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=10)

        # Atoms are deep inside the droplet (radius ~28 Å, atoms within ~5 Å of origin).
        # Droplet potential ~ 0 there. So E_pot is dominated by Coulomb partner.
        # At R0_GS = 9 Å pair separation, E_pair = 14.4/9 = 1.6 eV.
        # Per atom (half-half split): 0.8 eV.
        # If the bug were present, E_pot would be ~0 (just droplet, near 0 inside).
        for i in range(2 * cfg.num_molecules):
            assert ion.E_pot_eV[i, 0] > 0.5, (
                f"atom {i} E_pot = {ion.E_pot_eV[i, 0]}, "
                f"expected ~0.8 eV (Coulomb half-pair)"
            )


# ===========================================================================
# Scope-check error handling
# ===========================================================================
class TestScopeChecks:
    def test_effusive_dynamics_raises(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        cfg2 = replace(cfg, effusive_dynamics=True)
        with pytest.raises(NotImplementedError, match="effusive_dynamics"):
            build_initial_ion_state(cfg2, neutral, num_steps_ion=10)

    def test_single_charge_ionization_raises(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        cfg2 = replace(cfg, single_charge_ionization_allowed=True)
        with pytest.raises(NotImplementedError,
                           match="single_charge_ionization_allowed"):
            build_initial_ion_state(cfg2, neutral, num_steps_ion=10)

    def test_additional_droplet_charges_raises(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        cfg2 = replace(cfg, additional_droplet_charges=2)
        with pytest.raises(NotImplementedError, match="additional_droplet_charges"):
            build_initial_ion_state(cfg2, neutral, num_steps_ion=10)

    def test_highly_charged_iodine_raises(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        cfg2 = replace(cfg, highly_charged_iodine=True)
        with pytest.raises(NotImplementedError, match="highly_charged_iodine"):
            build_initial_ion_state(cfg2, neutral, num_steps_ion=10)


# ===========================================================================
# Validation
# ===========================================================================
class TestValidation:
    def test_num_molecules_mismatch_raises(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        cfg2 = replace(cfg, num_molecules=cfg.num_molecules + 1)
        with pytest.raises(ValueError, match="num_molecules"):
            build_initial_ion_state(cfg2, neutral, num_steps_ion=10)

    def test_num_steps_must_be_positive(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        with pytest.raises(ValueError, match="num_steps_ion"):
            build_initial_ion_state(cfg, neutral, num_steps_ion=0)
        with pytest.raises(ValueError, match="num_steps_ion"):
            build_initial_ion_state(cfg, neutral, num_steps_ion=-1)
