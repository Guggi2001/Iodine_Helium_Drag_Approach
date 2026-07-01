"""Tests for i2_helium_md/simulation/checkpoint.py."""

import numpy as np
import pytest

from i2_helium_md import single_pulse_N2000
from i2_helium_md.simulation.checkpoint import (
    IonCheckpoint,
    NeutralCheckpoint,
    load_ion_checkpoint,
    load_neutral_checkpoint,
    save_ion_checkpoint,
    save_neutral_checkpoint,
)


# ===========================================================================
# Helper builders
# ===========================================================================
def _make_neutral_checkpoint(num_molecules: int = 5,
                              num_steps: int = 10) -> NeutralCheckpoint:
    """Construct a self-consistent NeutralCheckpoint with random arrays."""
    N = num_molecules
    T = num_steps
    rng = np.random.default_rng(0)

    return NeutralCheckpoint(
        num_molecules=N,
        time_ps=np.linspace(0, 1, T),
        positions_x=rng.standard_normal((2 * N, T)),
        positions_y=rng.standard_normal((2 * N, T)),
        positions_z=rng.standard_normal((2 * N, T)),
        velocities_x=rng.standard_normal((2 * N, T)),
        velocities_y=rng.standard_normal((2 * N, T)),
        velocities_z=rng.standard_normal((2 * N, T)),
        mass_kg=np.full(2 * N, 127 * 1.66e-27),
        droplet_radii=np.full(2 * N, 30.0),
        r0=np.full(N, 5.0),
        E_kin_eV=rng.standard_normal((2 * N, T)),
        E_pot_eV=rng.standard_normal((2 * N, T)),
        E_initial_eV=rng.standard_normal(N),
        E_dissip_eV=rng.standard_normal((2 * N, T)),
        L_droplet_eV_ps=rng.standard_normal((2 * N, T)),
    )


def _make_ion_checkpoint(num_molecules: int = 5, num_steps: int = 10) -> IonCheckpoint:
    N = num_molecules
    T = num_steps
    rng = np.random.default_rng(1)

    return IonCheckpoint(
        num_molecules=N,
        time_ps=np.linspace(0, 1, T),
        positions_x=rng.standard_normal((2 * N, T)),
        positions_y=rng.standard_normal((2 * N, T)),
        positions_z=rng.standard_normal((2 * N, T)),
        velocities_x=rng.standard_normal((2 * N, T)),
        velocities_y=rng.standard_normal((2 * N, T)),
        velocities_z=rng.standard_normal((2 * N, T)),
        positions_final_x=rng.standard_normal(2 * N),
        positions_final_y=rng.standard_normal(2 * N),
        positions_final_z=rng.standard_normal(2 * N),
        velocities_final_x=rng.standard_normal(2 * N),
        velocities_final_y=rng.standard_normal(2 * N),
        velocities_final_z=rng.standard_normal(2 * N),
        mass_kg=np.full(2 * N, 127 * 1.66e-27),
        mass_final_kg=np.full(2 * N, 127 * 1.66e-27),
        mass_history_kg=np.full((2 * N, T), 127 * 1.66e-27),
        droplet_radii_angstrom=np.full(2 * N, 30.0),
        E_kin_eV=rng.standard_normal((2 * N, T)),
        E_pot_eV=rng.standard_normal((2 * N, T)),
        E_dissip_eV=rng.standard_normal((2 * N, T)),
        E_mass_transfer_eV=rng.standard_normal((2 * N, T)),
        E_int_eV=rng.standard_normal((2 * N, T)),
        n_shell=rng.integers(14, 22, size=(2 * N, T)).astype(float),
        b_ion_outside=np.zeros(N, dtype=bool),
        relative_loss_per_ps=rng.standard_normal((2 * N, T)),
        number_of_collisions=np.zeros((2 * N, T), dtype=int),
        temperature_diagnostic=rng.standard_normal((T, 3)),
        mass_scenario="anchored_discrete",
    )


# ===========================================================================
# Round-trip tests
# ===========================================================================
class TestRoundTrip:
    def test_neutral_round_trip(self, tmp_path):
        """Save then load a neutral checkpoint -- arrays must be byte-identical."""
        ckpt = _make_neutral_checkpoint(num_molecules=5, num_steps=10)
        path = save_neutral_checkpoint(ckpt, tmp_path / "n.npz")
        assert path.exists()

        loaded = load_neutral_checkpoint(path)
        assert loaded.num_molecules == ckpt.num_molecules
        assert loaded.schema_version == ckpt.schema_version
        np.testing.assert_array_equal(loaded.positions_x, ckpt.positions_x)
        np.testing.assert_array_equal(loaded.velocities_z, ckpt.velocities_z)
        np.testing.assert_array_equal(loaded.E_kin_eV, ckpt.E_kin_eV)
        np.testing.assert_array_equal(loaded.mass_kg, ckpt.mass_kg)

    def test_ion_round_trip(self, tmp_path):
        ckpt = _make_ion_checkpoint(num_molecules=3, num_steps=20)
        path = save_ion_checkpoint(ckpt, tmp_path / "i.npz")
        loaded = load_ion_checkpoint(path)
        assert loaded.num_molecules == 3
        assert loaded.schema_version == 7
        np.testing.assert_array_equal(loaded.positions_final_x,
                                       ckpt.positions_final_x)
        np.testing.assert_array_equal(loaded.b_ion_outside, ckpt.b_ion_outside)
        # v6/v7 fields round-trip: renamed mass-transfer channel, the v7
        # E_int reservoir, n_shell, and the scalar mass_scenario metadata
        # (recovered as a Python str).
        np.testing.assert_array_equal(loaded.E_mass_transfer_eV,
                                       ckpt.E_mass_transfer_eV)
        np.testing.assert_array_equal(loaded.E_int_eV, ckpt.E_int_eV)
        np.testing.assert_array_equal(loaded.n_shell, ckpt.n_shell)
        assert loaded.mass_scenario == "anchored_discrete"
        assert isinstance(loaded.mass_scenario, str)

    def test_extension_added_automatically(self, tmp_path):
        """Saving without .npz extension should add it."""
        ckpt = _make_neutral_checkpoint()
        path = save_neutral_checkpoint(ckpt, tmp_path / "no_extension")
        assert path.suffix == ".npz"

    def test_load_without_extension(self, tmp_path):
        """load() should also tolerate a missing .npz extension."""
        ckpt = _make_neutral_checkpoint()
        save_neutral_checkpoint(ckpt, tmp_path / "foo.npz")
        loaded = load_neutral_checkpoint(tmp_path / "foo")  # no .npz
        assert loaded.num_molecules == ckpt.num_molecules


# ===========================================================================
# Validation tests
# ===========================================================================
class TestValidation:
    def test_load_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_neutral_checkpoint(tmp_path / "does_not_exist.npz")

    def test_cfg_match_passes(self, tmp_path):
        ckpt = _make_neutral_checkpoint(num_molecules=2000)
        save_neutral_checkpoint(ckpt, tmp_path / "n.npz")
        cfg = single_pulse_N2000(num_molecules=2000)
        loaded = load_neutral_checkpoint(tmp_path / "n.npz", cfg=cfg)
        assert loaded.num_molecules == 2000

    def test_cfg_mismatch_raises(self, tmp_path):
        ckpt = _make_neutral_checkpoint(num_molecules=2000)
        save_neutral_checkpoint(ckpt, tmp_path / "n.npz")
        cfg = single_pulse_N2000(num_molecules=500)
        with pytest.raises(ValueError, match="num_molecules"):
            load_neutral_checkpoint(tmp_path / "n.npz", cfg=cfg)

    def test_wrong_schema_version_raises(self, tmp_path):
        """Manually-written file with wrong version should fail loudly."""
        path = tmp_path / "bad.npz"
        np.savez_compressed(path, schema_version=99,
                            num_molecules=1, time_ps=np.zeros(1))
        with pytest.raises(ValueError, match="schema_version"):
            load_neutral_checkpoint(path)

    def test_missing_schema_version_raises(self, tmp_path):
        """An old .npz without schema_version must be rejected, not silently loaded."""
        path = tmp_path / "old.npz"
        np.savez_compressed(path, num_molecules=1, time_ps=np.zeros(1))
        with pytest.raises(ValueError, match="schema_version"):
            load_neutral_checkpoint(path)

    def test_missing_field_raises(self, tmp_path):
        """An npz that has the version but is missing required fields must fail."""
        from i2_helium_md.simulation.checkpoint import _NEUTRAL_SCHEMA_VERSION
        path = tmp_path / "incomplete.npz"
        np.savez_compressed(
            path,
            schema_version=_NEUTRAL_SCHEMA_VERSION,
            num_molecules=1,
        )
        with pytest.raises(ValueError, match="missing fields"):
            load_neutral_checkpoint(path)

    def test_validator_catches_wrong_2N_T_shape(self, tmp_path):
        """If E_kin_eV has shape (N, T) instead of (2N, T), reject with a clear error.

        Regression guard for the schema-v1 -> v2 transition.
        """
        ckpt = _make_neutral_checkpoint(num_molecules=4, num_steps=5)
        # Forcibly corrupt one field to the old (N, T) shape
        ckpt.E_kin_eV = np.zeros((ckpt.num_molecules, ckpt.time_ps.size))
        save_neutral_checkpoint(ckpt, tmp_path / "wrong.npz")
        cfg = single_pulse_N2000(num_molecules=4)
        with pytest.raises(ValueError, match="E_kin_eV"):
            load_neutral_checkpoint(tmp_path / "wrong.npz", cfg=cfg)

    def test_validator_catches_inconsistent_num_steps(self, tmp_path):
        """All trajectory arrays must agree on num_steps."""
        ckpt = _make_neutral_checkpoint(num_molecules=4, num_steps=5)
        # Make E_pot_eV's num_steps inconsistent
        ckpt.E_pot_eV = np.zeros((2 * ckpt.num_molecules, 7))
        save_neutral_checkpoint(ckpt, tmp_path / "wrong.npz")
        cfg = single_pulse_N2000(num_molecules=4)
        with pytest.raises(ValueError, match="num_steps"):
            load_neutral_checkpoint(tmp_path / "wrong.npz", cfg=cfg)


# ===========================================================================
# File hygiene
# ===========================================================================
class TestFileLayout:
    def test_creates_parent_directory(self, tmp_path):
        """Path with non-existent parent dir should be created."""
        ckpt = _make_neutral_checkpoint()
        path = save_neutral_checkpoint(ckpt, tmp_path / "deep" / "nest" / "n.npz")
        assert path.exists()

    def test_compressed_smaller_than_raw(self, tmp_path):
        """Sanity check: savez_compressed actually saves space.
        We check the file is smaller than the raw bytes of all arrays.
        """
        ckpt = _make_neutral_checkpoint(num_molecules=20, num_steps=200)
        path = save_neutral_checkpoint(ckpt, tmp_path / "big.npz")
        raw_bytes = sum(
            v.nbytes for v in ckpt.__dict__.values()
            if isinstance(v, np.ndarray)
        )
        # compressed should be clearly smaller than raw (random data still
        # compresses ~10-20% with deflate)
        assert path.stat().st_size < raw_bytes


# ===========================================================================
# Ion schema v6 (Tier-1a mass dynamics) + v5 back-compat shim
# ===========================================================================
class TestIonSchemaV6:
    def test_n_shell_wrong_shape_rejected(self, tmp_path):
        """n_shell must be (2N, num_steps) like the other trajectory arrays."""
        ckpt = _make_ion_checkpoint(num_molecules=3, num_steps=6)
        # Corrupt n_shell to the wrong leading dimension (N instead of 2N).
        ckpt.n_shell = np.zeros((ckpt.num_molecules, ckpt.time_ps.size))
        save_ion_checkpoint(ckpt, tmp_path / "wrong.npz")
        cfg = single_pulse_N2000(num_molecules=3)
        with pytest.raises(ValueError, match="n_shell"):
            load_ion_checkpoint(tmp_path / "wrong.npz", cfg=cfg)

    def test_mass_scenario_defaults_to_fixed(self):
        """A v6 checkpoint built without an explicit mass_scenario is 'fixed'."""
        ckpt = IonCheckpoint(
            **{k: v for k, v in _make_ion_checkpoint(2, 3).__dict__.items()
               if k not in ("mass_scenario", "schema_version")}
        )
        assert ckpt.mass_scenario == "fixed"

    def test_v5_backcompat_shim(self, tmp_path):
        """A legacy v5 ion .npz cascades v5->v6->v7 via the migration shim.

        The shim maps the renamed mass-transfer field, synthesizes n_shell
        from mass_history_kg, defaults mass_scenario to 'fixed' (v5->v6),
        and synthesizes an all-zero E_int_eV with a warning (v6->v7).
        """
        from i2_helium_md.physics.constants import MASS_HE_AMU, MASS_I_ION_AMU, U

        ckpt = _make_ion_checkpoint(num_molecules=2, num_steps=5)
        path = save_ion_checkpoint(ckpt, tmp_path / "legacy.npz")

        # Rewrite the .npz to mimic a v5 file: old field name, no n_shell,
        # no mass_scenario, no E_int_eV, schema_version=5.
        with np.load(path, allow_pickle=False) as z:
            data = {k: z[k] for k in z.files}
        data["E_mass_attach_defect_eV"] = data.pop("E_mass_transfer_eV")
        data.pop("n_shell")
        data.pop("mass_scenario", None)
        data.pop("E_int_eV")
        data["schema_version"] = np.asarray(5)
        np.savez_compressed(path, **data)

        # The v6->v7 arm warns about the synthesized reservoir.
        with pytest.warns(UserWarning, match="E_int"):
            loaded = load_ion_checkpoint(path)

        # Version walked all the way to v7; renamed field preserved verbatim.
        assert loaded.schema_version == 7
        np.testing.assert_array_equal(loaded.E_mass_transfer_eV,
                                       ckpt.E_mass_transfer_eV)
        # mass_scenario defaulted; n_shell synthesized from mass_history_kg
        # via the same rule the writer uses.
        assert loaded.mass_scenario == "fixed"
        assert loaded.n_shell.shape == loaded.E_mass_transfer_eV.shape
        expected_n = np.rint(
            (ckpt.mass_history_kg / U - MASS_I_ION_AMU) / MASS_HE_AMU
        )
        np.testing.assert_array_equal(loaded.n_shell, expected_n)
        # E_int reservoir synthesized as all-zero, same shape as the other
        # (2N, T) trajectory arrays (4-term-equivalent).
        assert loaded.E_int_eV.shape == loaded.E_mass_transfer_eV.shape
        np.testing.assert_array_equal(
            loaded.E_int_eV, np.zeros_like(loaded.E_mass_transfer_eV)
        )

    def test_pre_v5_still_rejected(self, tmp_path):
        """The shim only upgrades v5; older versions still fail the check."""
        ckpt = _make_ion_checkpoint(num_molecules=2, num_steps=3)
        path = save_ion_checkpoint(ckpt, tmp_path / "v4.npz")
        with np.load(path, allow_pickle=False) as z:
            data = {k: z[k] for k in z.files}
        data["schema_version"] = np.asarray(4)
        np.savez_compressed(path, **data)
        with pytest.raises(ValueError, match="schema_version"):
            load_ion_checkpoint(path)


# ===========================================================================
# Ion schema v7 (Tier-2 E_int reservoir) + v6 back-compat shim
# ===========================================================================
class TestIonSchemaV7:
    def _write_v6_file(self, tmp_path):
        """Save a current checkpoint, then strip it back to a v6 on-disk file.

        Drops ``E_int_eV`` (absent before v7) and stamps ``schema_version``
        back to 6, mimicking a Tier-1a-era ``ion.npz``.
        """
        ckpt = _make_ion_checkpoint(num_molecules=2, num_steps=5)
        path = save_ion_checkpoint(ckpt, tmp_path / "v6.npz")
        with np.load(path, allow_pickle=False) as z:
            data = {k: z[k] for k in z.files}
        data.pop("E_int_eV")
        data["schema_version"] = np.asarray(6)
        np.savez_compressed(path, **data)
        return ckpt, path

    def test_v7_round_trip_preserves_E_int(self, tmp_path):
        """A genuine v7 file round-trips E_int_eV bit-for-bit."""
        ckpt = _make_ion_checkpoint(num_molecules=3, num_steps=7)
        path = save_ion_checkpoint(ckpt, tmp_path / "v7.npz")
        loaded = load_ion_checkpoint(path)
        assert loaded.schema_version == 7
        np.testing.assert_array_equal(loaded.E_int_eV, ckpt.E_int_eV)

    def test_v6_backcompat_shim_synthesizes_zeros_and_warns(self, tmp_path):
        """A v6 ion .npz loads under v7 with an all-zero E_int_eV + a warning."""
        ckpt, path = self._write_v6_file(tmp_path)
        with pytest.warns(UserWarning, match="E_int"):
            loaded = load_ion_checkpoint(path)
        assert loaded.schema_version == 7
        # Every v6 field survives verbatim; the reservoir is synthesized zeros.
        np.testing.assert_array_equal(loaded.E_mass_transfer_eV,
                                       ckpt.E_mass_transfer_eV)
        assert loaded.E_int_eV.shape == ckpt.E_mass_transfer_eV.shape
        np.testing.assert_array_equal(
            loaded.E_int_eV, np.zeros_like(ckpt.E_mass_transfer_eV)
        )

    def test_genuine_v7_missing_E_int_raises(self, tmp_path):
        """A file stamped v7 but missing E_int_eV must fail loudly, not zero-fill.

        The migration only synthesizes zeros for pre-v7 files; a genuine v7
        that lost the field is a corrupt/incomplete write and must raise.
        """
        ckpt = _make_ion_checkpoint(num_molecules=2, num_steps=4)
        path = save_ion_checkpoint(ckpt, tmp_path / "broken_v7.npz")
        with np.load(path, allow_pickle=False) as z:
            data = {k: z[k] for k in z.files}
        data.pop("E_int_eV")  # schema_version stays 7
        np.savez_compressed(path, **data)
        with pytest.raises(ValueError, match="missing fields"):
            load_ion_checkpoint(path)

    def test_wrong_shape_E_int_rejected(self, tmp_path):
        """E_int_eV must be (2N, num_steps) like the other trajectory arrays."""
        ckpt = _make_ion_checkpoint(num_molecules=3, num_steps=6)
        # Corrupt E_int_eV to the wrong leading dimension (N instead of 2N).
        ckpt.E_int_eV = np.zeros((ckpt.num_molecules, ckpt.time_ps.size))
        save_ion_checkpoint(ckpt, tmp_path / "wrong.npz")
        cfg = single_pulse_N2000(num_molecules=3)
        with pytest.raises(ValueError, match="E_int_eV"):
            load_ion_checkpoint(tmp_path / "wrong.npz", cfg=cfg)
