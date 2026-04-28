"""Tests for i2_helium_md/simulation/run_directory.py."""

import json
import numpy as np
import pytest

from i2_helium_md import single_pulse_N2000
from i2_helium_md.simulation.checkpoint import (
    IonCheckpoint,
    NeutralCheckpoint,
)
from i2_helium_md.simulation.run_directory import RunDirectory

# Reuse helpers from test_checkpoint
from .test_checkpoint import _make_neutral_checkpoint, _make_ion_checkpoint  # type: ignore


# ===========================================================================
# Construction and introspection
# ===========================================================================
class TestConstruction:
    def test_construction_does_not_create_dir(self, tmp_path):
        """Instantiating RunDirectory should not touch the filesystem."""
        run = RunDirectory(tmp_path / "new_run")
        assert not (tmp_path / "new_run").exists()
        assert run.exists() is False

    def test_repr_includes_path(self, tmp_path):
        run = RunDirectory(tmp_path / "x")
        s = repr(run)
        assert "RunDirectory" in s
        assert "x" in s

    def test_paths_compose_correctly(self, tmp_path):
        run = RunDirectory(tmp_path / "r")
        assert run.cfg_path.name == "cfg.json"
        assert run.neutral_path.name == "neutral.npz"
        assert run.ion_path.name == "ion.npz"
        assert run.neutral_path.parent == tmp_path / "r"

    def test_existence_flags(self, tmp_path):
        run = RunDirectory(tmp_path / "r")
        assert not run.has_cfg()
        assert not run.has_neutral()
        assert not run.has_ion()


# ===========================================================================
# SimConfig save/load
# ===========================================================================
class TestCfgRoundTrip:
    def test_save_creates_directory(self, tmp_path):
        run = RunDirectory(tmp_path / "auto_create")
        cfg = single_pulse_N2000(seed=42)
        run.save_cfg(cfg)
        assert run.exists()

    def test_round_trip(self, tmp_path):
        run = RunDirectory(tmp_path / "r")
        cfg = single_pulse_N2000(seed=42, num_molecules=500)
        run.save_cfg(cfg)
        loaded = run.load_cfg()
        assert loaded.seed == 42
        assert loaded.num_molecules == 500

    def test_load_missing_raises(self, tmp_path):
        run = RunDirectory(tmp_path / "empty")
        with pytest.raises(FileNotFoundError):
            run.load_cfg()

    def test_unknown_fields_raise(self, tmp_path):
        """A cfg.json with extra fields signals version skew -- fail loudly."""
        run = RunDirectory(tmp_path / "r")
        cfg = single_pulse_N2000(seed=1)
        run.save_cfg(cfg)
        # corrupt the file by adding an unknown field
        payload = json.loads(run.cfg_path.read_text())
        payload["mystery_field_from_the_future"] = 42
        run.cfg_path.write_text(json.dumps(payload))
        with pytest.raises(ValueError, match="unknown fields"):
            run.load_cfg()

    def test_json_is_human_readable(self, tmp_path):
        """cfg.json should be indented for legibility, not minified."""
        run = RunDirectory(tmp_path / "r")
        cfg = single_pulse_N2000(seed=42)
        run.save_cfg(cfg)
        text = run.cfg_path.read_text()
        # at least one line break per top-level field
        assert text.count("\n") > 5


# ===========================================================================
# Neutral checkpoint flow
# ===========================================================================
class TestNeutralFlow:
    def test_save_load_round_trip(self, tmp_path):
        run = RunDirectory(tmp_path / "r")
        ckpt = _make_neutral_checkpoint(num_molecules=5, num_steps=10)
        run.save_neutral(ckpt)
        assert run.has_neutral()
        loaded = run.load_neutral()
        np.testing.assert_array_equal(loaded.positions_x, ckpt.positions_x)

    def test_save_with_cfg_creates_cfg_json(self, tmp_path):
        """Passing cfg= to save_neutral should also save cfg.json."""
        run = RunDirectory(tmp_path / "r")
        cfg = single_pulse_N2000(num_molecules=5)
        ckpt = _make_neutral_checkpoint(num_molecules=5)
        run.save_neutral(ckpt, cfg=cfg)
        assert run.has_cfg()

    def test_save_with_cfg_does_not_overwrite_existing(self, tmp_path):
        """Saving twice with different cfgs should NOT overwrite the first."""
        run = RunDirectory(tmp_path / "r")
        cfg_a = single_pulse_N2000(num_molecules=5, seed=1)
        cfg_b = single_pulse_N2000(num_molecules=5, seed=99)
        ckpt = _make_neutral_checkpoint(num_molecules=5)
        run.save_neutral(ckpt, cfg=cfg_a)
        run.save_neutral(ckpt, cfg=cfg_b)
        # The first cfg should still be there
        assert run.load_cfg().seed == 1

    def test_load_uses_cfg_json_for_validation(self, tmp_path):
        """When load_neutral is called with no cfg= arg but cfg.json exists,
        validation should still happen.
        """
        run = RunDirectory(tmp_path / "r")
        cfg = single_pulse_N2000(num_molecules=5)
        ckpt = _make_neutral_checkpoint(num_molecules=5)
        run.save_neutral(ckpt, cfg=cfg)
        # passing no cfg to load -- should pull from cfg.json and validate
        loaded = run.load_neutral()
        assert loaded.num_molecules == 5

    def test_load_explicit_cfg_overrides_cfg_json(self, tmp_path):
        """If the user passes a cfg explicitly, it takes precedence over
        cfg.json (useful for ad-hoc validation against a different config).
        """
        run = RunDirectory(tmp_path / "r")
        cfg_saved = single_pulse_N2000(num_molecules=5)
        ckpt = _make_neutral_checkpoint(num_molecules=5)
        run.save_neutral(ckpt, cfg=cfg_saved)
        cfg_other = single_pulse_N2000(num_molecules=99)
        with pytest.raises(ValueError, match="num_molecules"):
            run.load_neutral(cfg=cfg_other)

    def test_load_missing_neutral_raises(self, tmp_path):
        run = RunDirectory(tmp_path / "empty")
        with pytest.raises(FileNotFoundError):
            run.load_neutral()


# ===========================================================================
# Ion checkpoint flow
# ===========================================================================
class TestIonFlow:
    def test_save_load_round_trip(self, tmp_path):
        run = RunDirectory(tmp_path / "r")
        ckpt = _make_ion_checkpoint(num_molecules=3, num_steps=8)
        run.save_ion(ckpt)
        assert run.has_ion()
        loaded = run.load_ion()
        assert loaded.num_molecules == 3
        np.testing.assert_array_equal(
            loaded.positions_final_x, ckpt.positions_final_x,
        )

    def test_independent_of_neutral(self, tmp_path):
        """Saving an ion checkpoint should not require a neutral one."""
        run = RunDirectory(tmp_path / "r")
        ckpt = _make_ion_checkpoint(num_molecules=2)
        run.save_ion(ckpt)
        assert run.has_ion()
        assert not run.has_neutral()


# ===========================================================================
# End-to-end pipeline
# ===========================================================================
class TestPipeline:
    def test_full_pipeline_two_processes(self, tmp_path):
        """Simulate two scripts using the same run dir at different times."""
        # ---- Script A: neutral stage ----
        run_a = RunDirectory(tmp_path / "p")
        cfg = single_pulse_N2000(num_molecules=4, seed=42)
        run_a.save_cfg(cfg)
        n_ckpt = _make_neutral_checkpoint(num_molecules=4)
        run_a.save_neutral(n_ckpt)

        # ---- Script B: ion stage (fresh process, only knows the path) ----
        run_b = RunDirectory(tmp_path / "p")
        cfg_b = run_b.load_cfg()
        assert cfg_b.seed == 42
        n_loaded = run_b.load_neutral()
        np.testing.assert_array_equal(n_loaded.positions_x, n_ckpt.positions_x)

        i_ckpt = _make_ion_checkpoint(num_molecules=4)
        run_b.save_ion(i_ckpt)

        # ---- Script C: postprocess (different process again) ----
        run_c = RunDirectory(tmp_path / "p")
        i_loaded = run_c.load_ion()
        np.testing.assert_array_equal(
            i_loaded.positions_final_x, i_ckpt.positions_final_x,
        )

    def test_two_runs_side_by_side(self, tmp_path):
        """Two run dirs should not interfere with each other."""
        run_x = RunDirectory(tmp_path / "x")
        run_y = RunDirectory(tmp_path / "y")
        cfg_x = single_pulse_N2000(num_molecules=4, seed=1)
        cfg_y = single_pulse_N2000(num_molecules=4, seed=2)
        run_x.save_cfg(cfg_x)
        run_y.save_cfg(cfg_y)
        assert run_x.load_cfg().seed == 1
        assert run_y.load_cfg().seed == 2
