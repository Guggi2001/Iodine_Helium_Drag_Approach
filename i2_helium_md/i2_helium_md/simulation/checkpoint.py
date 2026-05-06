"""Checkpoint I/O for neutral and ion simulation stages.

Replaces the legacy ``save('neutral_propagation_checkpoint', ...)`` /
``load('neutral_propagation_checkpoint')`` calls in
``vmi_sim_3d_neutral_propa_HeDFT_mimic.m`` and ``vmi_sim_3d_ion_propa.m``.

Design
------
Two dataclasses encode exactly what each stage needs to either resume or be
analyzed: positions, velocities, time axis, masses, droplet radii, plus
energy diagnostics. We deliberately do **not** persist anything that lives
in ``constants.py`` (eV, u, ...) or in ``cfg`` (binding energies, flags,
mode switches) -- those are recovered at load time from the current source
of truth, with shape validation against the checkpoint.

File format
-----------
``.npz`` (NumPy's native binary format). Pros:

- Native to NumPy -- no extra dependency.
- One file holds all named arrays.
- ~100x smaller install footprint than HDF5.
- Forward-compatible: extra fields can be added; the loader uses
  explicit field names with explicit defaults.

Schema versioning
-----------------
Each checkpoint includes a ``schema_version`` integer. The loader checks
this on read and refuses to load incompatible versions. Bump the version
whenever fields are removed or their meaning changes; additions are
backward-compatible.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from pathlib import Path
from typing import Any

import numpy as np

from ..config import SimConfig


# ===========================================================================
# Schema versions: bump when removing/renaming fields, OR when shapes change
# in a backward-incompatible way.
#
# Version history:
#   1 -- initial; energy/L_droplet diagnostics had per-molecule shape (N, T).
#   2 -- energy/L_droplet diagnostics moved to per-atom shape (2N, T)
#        to match the legacy MATLAB code and allow per-atom debugging.
#
# Ion-specific history:
#   2 -- initial ion checkpoint shape, matching neutral v2.
#   3 -- adds three fields needed for postprocess and diagnostics:
#          * droplet_radii_angstrom (2N,)   -- droplet radius per atom
#          * mass_history_kg        (2N, T) -- mass over time (changes
#                                              as helium attaches)
#          * E_dissip_eV            (2N, T) -- cumulative energy dissipated
#                                              per atom (matches neutral)
#   4 -- adds the mass-attachment kinetic-energy defect diagnostic that
#        the legacy MATLAB code tracks (vmi_sim_3d_ion_propa.m:762).
#        When 4 amu of helium attaches at the atom's current velocity,
#        the recomputed E_kin = 1/2 m_new v^2 is spuriously larger than
#        the pre-attachment E_kin by 1/2 (m_new - m_old) v^2. The defect
#        is the negative of that running sum, so that
#        E_kin + E_pot + E_dissip + E_mass_attach_defect is conserved
#        modulo Verlet drift on each side.
#          * E_mass_attach_defect_eV (2N, T) -- cumulative attach defect
# ===========================================================================
_NEUTRAL_SCHEMA_VERSION: int = 2
_ION_SCHEMA_VERSION: int = 4


# ===========================================================================
# Neutral checkpoint -- saved at end of neutral propagation
# ===========================================================================
@dataclass
class NeutralCheckpoint:
    """State at the end of neutral propagation.

    All arrays are NumPy arrays with explicit dtypes. Shapes:

    * ``num_molecules``      : scalar int
    * ``time_ps``            : (num_steps,)
    * ``positions_x``        : (2 * num_molecules, num_steps)  Angstrom
    * ``positions_y``        : (2 * num_molecules, num_steps)  Angstrom
    * ``positions_z``        : (2 * num_molecules, num_steps)  Angstrom
    * ``velocities_x``       : (2 * num_molecules, num_steps)  Angstrom/ps
    * ``velocities_y``       : (2 * num_molecules, num_steps)  Angstrom/ps
    * ``velocities_z``       : (2 * num_molecules, num_steps)  Angstrom/ps
    * ``mass_kg``            : (2 * num_molecules,)            kg
    * ``droplet_radii``      : (2 * num_molecules,)            Angstrom
    * ``r0``                 : (num_molecules,)                Angstrom (initial radial distance)
    * ``E_kin_eV``           : (2 * num_molecules, num_steps)  eV (per-atom)
    * ``E_pot_eV``           : (2 * num_molecules, num_steps)  eV (per-atom; pair energy split 50/50)
    * ``E_initial_eV``       : (num_molecules,)                eV (per-molecule total at t=0)
    * ``E_dissip_eV``        : (2 * num_molecules, num_steps)  eV (cumulative, per-atom)
    * ``L_droplet_eV_ps``    : (2 * num_molecules, num_steps)  eV*ps (droplet potential work, per-atom)

    Two-atom layout convention: indices [0, num_molecules) are the first
    atom of each molecule; indices [num_molecules, 2*num_molecules) are
    the second atom. Same convention as the rest of the codebase.

    Per-atom vs per-molecule: ``E_initial_eV`` is per-molecule because
    it represents the photon energy delivered to the molecule as a whole.
    All other energy/dissipation arrays are per-atom, matching MATLAB.
    Per-molecule values are recovered by summing atom 1 + atom 2 entries.
    """

    num_molecules: int
    time_ps: np.ndarray
    positions_x: np.ndarray
    positions_y: np.ndarray
    positions_z: np.ndarray
    velocities_x: np.ndarray
    velocities_y: np.ndarray
    velocities_z: np.ndarray
    mass_kg: np.ndarray
    droplet_radii: np.ndarray
    r0: np.ndarray
    E_kin_eV: np.ndarray
    E_pot_eV: np.ndarray
    E_initial_eV: np.ndarray
    E_dissip_eV: np.ndarray
    L_droplet_eV_ps: np.ndarray
    schema_version: int = _NEUTRAL_SCHEMA_VERSION


# ===========================================================================
# Ion checkpoint -- saved at end of ion propagation
# ===========================================================================
@dataclass
class IonCheckpoint:
    """Final state at the end of ion propagation.

    All arrays are NumPy arrays with explicit dtypes. Shapes:

    * ``num_molecules``      : scalar int
    * ``time_ps``            : (num_steps,)
    * ``positions_x``        : (2 * num_molecules, num_steps)  Angstrom
    * ``positions_y``        : (2 * num_molecules, num_steps)  Angstrom
    * ``positions_z``        : (2 * num_molecules, num_steps)  Angstrom
    * ``velocities_x``       : (2 * num_molecules, num_steps)  Angstrom/ps
    * ``velocities_y``       : (2 * num_molecules, num_steps)  Angstrom/ps
    * ``velocities_z``       : (2 * num_molecules, num_steps)  Angstrom/ps
    * ``positions_final_x``  : (2 * num_molecules,)            Angstrom (asymptotic)
    * ``positions_final_y``  : (2 * num_molecules,)            Angstrom
    * ``positions_final_z``  : (2 * num_molecules,)            Angstrom
    * ``velocities_final_x`` : (2 * num_molecules,)            Angstrom/ps (used by VMI postprocess)
    * ``velocities_final_y`` : (2 * num_molecules,)            Angstrom/ps
    * ``velocities_final_z`` : (2 * num_molecules,)            Angstrom/ps
    * ``mass_kg``            : (2 * num_molecules,)            kg (initial mass)
    * ``mass_final_kg``      : (2 * num_molecules,)            kg (after possible mass attachment)
    * ``mass_history_kg``    : (2 * num_molecules, num_steps)  kg (mass over time)
    * ``droplet_radii_angstrom``: (2 * num_molecules,)         Angstrom (per atom; same value for the two atoms of a molecule)
    * ``E_kin_eV``           : (2 * num_molecules, num_steps)  eV (per-atom)
    * ``E_pot_eV``           : (2 * num_molecules, num_steps)  eV (per-atom)
    * ``E_dissip_eV``        : (2 * num_molecules, num_steps)  eV (per-atom, cumulative)
    * ``E_mass_attach_defect_eV``: (2 * num_molecules, num_steps) eV (per-atom, cumulative)
    * ``b_ion_outside``      : (num_molecules,) bool           True if ion exited droplet
    * ``relative_loss_per_ps``: (2 * num_molecules, num_steps) 1/ps (per-atom energy loss rate)
    * ``number_of_collisions``: (2 * num_molecules, num_steps) int (cumulative, per-atom)

    Schema v4 (the current version) differs from v3 by adding
    ``E_mass_attach_defect_eV`` (the legacy MATLAB
    ``E_mass_attach_defect`` diagnostic; see the schema-version block
    above). Older files cannot be loaded; rerun the ion stage to upgrade.
    """

    num_molecules: int
    time_ps: np.ndarray
    positions_x: np.ndarray
    positions_y: np.ndarray
    positions_z: np.ndarray
    velocities_x: np.ndarray
    velocities_y: np.ndarray
    velocities_z: np.ndarray
    positions_final_x: np.ndarray
    positions_final_y: np.ndarray
    positions_final_z: np.ndarray
    velocities_final_x: np.ndarray
    velocities_final_y: np.ndarray
    velocities_final_z: np.ndarray
    mass_kg: np.ndarray
    mass_final_kg: np.ndarray
    mass_history_kg: np.ndarray
    droplet_radii_angstrom: np.ndarray
    E_kin_eV: np.ndarray
    E_pot_eV: np.ndarray
    E_dissip_eV: np.ndarray
    E_mass_attach_defect_eV: np.ndarray
    b_ion_outside: np.ndarray
    relative_loss_per_ps: np.ndarray
    number_of_collisions: np.ndarray
    schema_version: int = _ION_SCHEMA_VERSION


# ===========================================================================
# Save / load functions
# ===========================================================================
def save_neutral_checkpoint(
    checkpoint: NeutralCheckpoint,
    path: str | Path,
) -> Path:
    """Save a NeutralCheckpoint to a .npz file.

    Parameters
    ----------
    checkpoint : NeutralCheckpoint
    path : str or Path
        Output path. ``.npz`` extension added if missing.

    Returns
    -------
    Path
        Path to the file written.
    """
    return _save_checkpoint(checkpoint, path)


def load_neutral_checkpoint(
    path: str | Path,
    cfg: SimConfig | None = None,
) -> NeutralCheckpoint:
    """Load a NeutralCheckpoint from a .npz file.

    Parameters
    ----------
    path : str or Path
        Input path.
    cfg : SimConfig, optional
        If provided, validate the checkpoint shape against
        ``cfg.num_molecules`` and raise on mismatch. If None, no shape
        validation is performed against config.

    Returns
    -------
    NeutralCheckpoint
    """
    return _load_checkpoint(
        path,
        dataclass_type=NeutralCheckpoint,
        expected_version=_NEUTRAL_SCHEMA_VERSION,
        cfg=cfg,
    )


def save_ion_checkpoint(
    checkpoint: IonCheckpoint,
    path: str | Path,
) -> Path:
    """Save an IonCheckpoint to a .npz file."""
    return _save_checkpoint(checkpoint, path)


def load_ion_checkpoint(
    path: str | Path,
    cfg: SimConfig | None = None,
) -> IonCheckpoint:
    """Load an IonCheckpoint from a .npz file."""
    return _load_checkpoint(
        path,
        dataclass_type=IonCheckpoint,
        expected_version=_ION_SCHEMA_VERSION,
        cfg=cfg,
    )


# ===========================================================================
# Generic implementation (private)
# ===========================================================================
def _save_checkpoint(checkpoint: Any, path: str | Path) -> Path:
    """Write a checkpoint dataclass to a .npz file."""
    p = Path(path)
    if p.suffix != ".npz":
        p = p.with_suffix(".npz")
    p.parent.mkdir(parents=True, exist_ok=True)

    payload = {}
    for f in fields(checkpoint):
        value = getattr(checkpoint, f.name)
        # Wrap scalars as 0-d arrays so .npz round-trips cleanly
        if np.isscalar(value):
            payload[f.name] = np.asarray(value)
        else:
            payload[f.name] = np.asarray(value)

    np.savez_compressed(p, **payload)
    return p


def _load_checkpoint(
    path: str | Path,
    *,
    dataclass_type: type,
    expected_version: int,
    cfg: SimConfig | None,
) -> Any:
    """Read a .npz file into a checkpoint dataclass."""
    p = Path(path)
    if p.suffix != ".npz":
        p = p.with_suffix(".npz")
    if not p.exists():
        raise FileNotFoundError(f"checkpoint not found: {p}")

    with np.load(p, allow_pickle=False) as npz:
        # 1. Schema version check
        if "schema_version" not in npz.files:
            raise ValueError(
                f"checkpoint at {p} has no schema_version field; "
                "it was written by an older version. Re-run the simulation."
            )
        version = int(npz["schema_version"])
        if version != expected_version:
            raise ValueError(
                f"checkpoint at {p} has schema_version={version}, "
                f"this code expects {expected_version}. Re-run the simulation."
            )

        # 2. Build kwargs for dataclass construction
        expected_fields = {f.name for f in fields(dataclass_type)}
        missing = expected_fields - set(npz.files)
        if missing:
            raise ValueError(
                f"checkpoint at {p} is missing fields: {sorted(missing)}"
            )

        kwargs: dict[str, Any] = {}
        for f in fields(dataclass_type):
            arr = npz[f.name]
            # unwrap 0-d arrays for int / scalar fields
            if f.type is int or f.name == "schema_version":
                kwargs[f.name] = int(arr)
            elif f.type is float:
                kwargs[f.name] = float(arr)
            else:
                kwargs[f.name] = np.asarray(arr)

    instance = dataclass_type(**kwargs)

    # 3. Optional shape validation against cfg
    if cfg is not None:
        _validate_against_cfg(instance, cfg, path=p)

    return instance


def _validate_against_cfg(
    checkpoint: Any,
    cfg: SimConfig,
    *,
    path: Path,
) -> None:
    """Cross-check checkpoint against a SimConfig.

    Currently checks ``num_molecules`` and array shapes. Add more
    invariants here as the pipeline matures (e.g. dt consistency,
    num_steps bounds).
    """
    if checkpoint.num_molecules != cfg.num_molecules:
        raise ValueError(
            f"checkpoint at {path} has num_molecules={checkpoint.num_molecules} "
            f"but cfg has {cfg.num_molecules}. Either rerun the upstream stage "
            f"or use the matching cfg."
        )

    N = cfg.num_molecules

    # (2N,) static per-atom arrays.
    expected_2N_shape = (2 * N,)
    static_2N_fields = ("mass_kg", "droplet_radii", "droplet_radii_angstrom",
                        "mass_final_kg",
                        "positions_final_x", "positions_final_y",
                        "positions_final_z", "velocities_final_x",
                        "velocities_final_y", "velocities_final_z")
    for fname in static_2N_fields:
        if hasattr(checkpoint, fname):
            arr = getattr(checkpoint, fname)
            if arr.shape != expected_2N_shape:
                raise ValueError(
                    f"checkpoint field {fname!r} has shape {arr.shape}, "
                    f"expected {expected_2N_shape}"
                )

    # (N,) per-molecule arrays.
    expected_N_shape = (N,)
    static_N_fields = ("r0", "E_initial_eV", "b_ion_outside")
    for fname in static_N_fields:
        if hasattr(checkpoint, fname):
            arr = getattr(checkpoint, fname)
            if arr.shape != expected_N_shape:
                raise ValueError(
                    f"checkpoint field {fname!r} has shape {arr.shape}, "
                    f"expected {expected_N_shape}"
                )

    # (2N, num_steps) trajectory and per-atom diagnostic arrays.
    # We don't pin num_steps because it depends on the run, but we
    # check the leading dimension and that all such fields agree.
    trajectory_2N_T_fields = (
        "positions_x", "positions_y", "positions_z",
        "velocities_x", "velocities_y", "velocities_z",
        "E_kin_eV", "E_pot_eV", "E_dissip_eV", "L_droplet_eV_ps",
        "E_mass_attach_defect_eV",
        "relative_loss_per_ps", "number_of_collisions",
        "mass_history_kg",
    )
    num_steps_seen: int | None = None
    for fname in trajectory_2N_T_fields:
        if hasattr(checkpoint, fname):
            arr = getattr(checkpoint, fname)
            if arr.ndim != 2 or arr.shape[0] != 2 * N:
                raise ValueError(
                    f"checkpoint field {fname!r} has shape {arr.shape}, "
                    f"expected ({2*N}, num_steps)"
                )
            if num_steps_seen is None:
                num_steps_seen = arr.shape[1]
            elif arr.shape[1] != num_steps_seen:
                raise ValueError(
                    f"checkpoint field {fname!r} has num_steps={arr.shape[1]}, "
                    f"but earlier fields had num_steps={num_steps_seen}"
                )

    # time_ps must agree with the trajectory length.
    if hasattr(checkpoint, "time_ps") and num_steps_seen is not None:
        ts = checkpoint.time_ps
        if ts.shape != (num_steps_seen,):
            raise ValueError(
                f"checkpoint time_ps has shape {ts.shape}, "
                f"expected ({num_steps_seen},)"
            )
