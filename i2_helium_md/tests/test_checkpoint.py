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
        E_mass_attach_defect_eV=rng.standard_normal((2 * N, T)),
        b_ion_outside=np.zeros(N, dtype=bool),
        relative_loss_per_ps=rng.standard_normal((2 * N, T)),
        number_of_collisions=np.zeros((2 * N, T), dtype=int),
        temperature_diagnostic=rng.standard_normal((T, 3)),
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
        np.testing.assert_array_equal(loaded.positions_final_x,
                                       ckpt.positions_final_x)
        np.testing.assert_array_equal(loaded.b_ion_outside, ckpt.b_ion_outside)

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
