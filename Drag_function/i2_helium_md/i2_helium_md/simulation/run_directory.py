"""A run-directory abstraction: one folder per simulation run.

Each run directory holds all artifacts produced by one end-to-end pipeline
execution -- neutral checkpoint, ion checkpoint, and a JSON dump of the
SimConfig that produced them. Convenience methods enforce the filename
conventions so user code never has to type ``"neutral.npz"`` or remember
where to put things.

Layout
------
::

    <root>/
        cfg.json        SimConfig dump (auto-saved when first artifact is written)
        neutral.npz     NeutralCheckpoint
        ion.npz         IonCheckpoint
        (figures/)      postprocess plots, optional
        (logs/)         text logs, optional

Usage
-----
::

    from i2_helium_md import single_pulse_N2000
    from i2_helium_md.simulation.run_directory import RunDirectory

    cfg = single_pulse_N2000(seed=42)
    run = RunDirectory("data/runs/my_first_run")
    run.save_cfg(cfg)

    # ... run neutral stage, get a NeutralCheckpoint ...
    run.save_neutral(neutral_ckpt)

    # later (possibly different process):
    run = RunDirectory("data/runs/my_first_run")
    cfg = run.load_cfg()
    neutral = run.load_neutral(cfg=cfg)
    # ... run ion stage ...
    run.save_ion(ion_ckpt)

Why this rather than raw paths
------------------------------
* Filenames are an enforced convention, not user invention.
* A single string (the run name) is the only thing two scripts need to
  agree on -- no risk of a typo in ``"neutral.npz"`` vs ``"Neutral.npz"``.
* The directory is **self-describing**: ``cfg.json`` next to the data
  documents exactly what produced it.
* Multiple runs sit side-by-side without name collisions.
* Future artifacts (pump.npz, probe.npz, postprocess outputs) extend the
  layout naturally without changing existing call sites.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any

from ..config import SimConfig
from .checkpoint import (
    IonCheckpoint,
    NeutralCheckpoint,
    load_ion_checkpoint,
    load_neutral_checkpoint,
    save_ion_checkpoint,
    save_neutral_checkpoint,
)


# ===========================================================================
# Filename conventions (single source of truth)
# ===========================================================================
_CFG_FILENAME: str = "cfg.json"
_NEUTRAL_FILENAME: str = "neutral.npz"
_ION_FILENAME: str = "ion.npz"


# ===========================================================================
# RunDirectory
# ===========================================================================
class RunDirectory:
    """One folder per simulation run.

    The directory is created on demand (not at construction time) so that
    instantiating a ``RunDirectory`` for inspection of an existing run is
    cheap and side-effect free. Directories are created automatically when
    the first artifact is saved.

    Parameters
    ----------
    root : str or Path
        Path to the run directory. Does not need to exist yet.
    """

    def __init__(self, root: str | Path) -> None:
        self.root: Path = Path(root)

    # -----------------------------------------------------------------
    # Identity / introspection
    # -----------------------------------------------------------------
    def __repr__(self) -> str:
        return f"RunDirectory({str(self.root)!r})"

    def exists(self) -> bool:
        """True if the directory itself has been created."""
        return self.root.is_dir()

    @property
    def cfg_path(self) -> Path:
        return self.root / _CFG_FILENAME

    @property
    def neutral_path(self) -> Path:
        return self.root / _NEUTRAL_FILENAME

    @property
    def ion_path(self) -> Path:
        return self.root / _ION_FILENAME

    def has_cfg(self) -> bool:
        return self.cfg_path.exists()

    def has_neutral(self) -> bool:
        return self.neutral_path.exists()

    def has_ion(self) -> bool:
        return self.ion_path.exists()

    # -----------------------------------------------------------------
    # SimConfig save/load (JSON)
    # -----------------------------------------------------------------
    def save_cfg(self, cfg: SimConfig) -> Path:
        """Serialise a SimConfig to ``<root>/cfg.json``.

        Returns the path written.
        """
        self.root.mkdir(parents=True, exist_ok=True)
        payload = dataclasses.asdict(cfg)
        with self.cfg_path.open("w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, sort_keys=False)
        return self.cfg_path

    def load_cfg(self) -> SimConfig:
        """Load a SimConfig from ``<root>/cfg.json``.

        Raises ``FileNotFoundError`` if the run does not have one.
        Raises ``ValueError`` if the JSON contains fields that are not in
        ``SimConfig`` (this signals a version skew between the run and the
        current code).
        """
        if not self.has_cfg():
            raise FileNotFoundError(
                f"no cfg.json in run directory {self.root}; "
                "save it with run.save_cfg(cfg) before loading."
            )
        with self.cfg_path.open("r", encoding="utf-8") as fh:
            payload: dict[str, Any] = json.load(fh)

        valid_fields = {f.name for f in dataclasses.fields(SimConfig)}
        unknown = set(payload) - valid_fields
        if unknown:
            raise ValueError(
                f"cfg.json in {self.root} has unknown fields: {sorted(unknown)}. "
                "This run was probably produced by a different version of the code."
            )
        return SimConfig(**payload)

    # -----------------------------------------------------------------
    # Neutral checkpoint
    # -----------------------------------------------------------------
    def save_neutral(
        self,
        ckpt: NeutralCheckpoint,
        cfg: SimConfig | None = None,
    ) -> Path:
        """Save a NeutralCheckpoint to ``<root>/neutral.npz``.

        Parameters
        ----------
        ckpt : NeutralCheckpoint
        cfg : SimConfig, optional
            If provided and no cfg.json exists yet, also save the cfg.
            If cfg.json already exists, this is a no-op (the existing file
            is left alone -- we never overwrite a saved config silently).
        """
        self.root.mkdir(parents=True, exist_ok=True)
        if cfg is not None and not self.has_cfg():
            self.save_cfg(cfg)
        return save_neutral_checkpoint(ckpt, self.neutral_path)

    def load_neutral(
        self,
        cfg: SimConfig | None = None,
    ) -> NeutralCheckpoint:
        """Load the NeutralCheckpoint from ``<root>/neutral.npz``.

        Parameters
        ----------
        cfg : SimConfig, optional
            If provided, validates the checkpoint shape against this cfg.
            If None and a ``cfg.json`` is present in the run directory,
            it is loaded and used for validation automatically.
        """
        if cfg is None and self.has_cfg():
            cfg = self.load_cfg()
        return load_neutral_checkpoint(self.neutral_path, cfg=cfg)

    # -----------------------------------------------------------------
    # Ion checkpoint
    # -----------------------------------------------------------------
    def save_ion(
        self,
        ckpt: IonCheckpoint,
        cfg: SimConfig | None = None,
    ) -> Path:
        """Save an IonCheckpoint to ``<root>/ion.npz``."""
        self.root.mkdir(parents=True, exist_ok=True)
        if cfg is not None and not self.has_cfg():
            self.save_cfg(cfg)
        return save_ion_checkpoint(ckpt, self.ion_path)

    def load_ion(
        self,
        cfg: SimConfig | None = None,
    ) -> IonCheckpoint:
        """Load the IonCheckpoint from ``<root>/ion.npz``."""
        if cfg is None and self.has_cfg():
            cfg = self.load_cfg()
        return load_ion_checkpoint(self.ion_path, cfg=cfg)
