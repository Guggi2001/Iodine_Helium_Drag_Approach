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
- Explicit field names -- one named array per dataclass field.

Schema versioning
-----------------
Each checkpoint includes a ``schema_version`` integer. The loader checks
this on read and refuses to load incompatible versions. The loader raises
on any missing field, so *any* schema change -- additions included -- needs
a version bump plus a migration arm in the shim (the v5->v6 and v6->v7
precedents; see the ion-specific history below).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Any

import numpy as np

from ..config import SimConfig
from ..physics.constants import MASS_HE_AMU, MASS_I_ION_AMU, U


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
#          * mass_history_kg        (2N, T) -- mass over time (helium
#                                              attaches OR sheds; NOT
#                                              monotone under the Tier-1a
#                                              anchored_discrete scenario)
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
#   5 -- adds the per-step legacy MATLAB temperature diagnostic from
#        vmi_sim_3d_ion_propa.m:683:
#          [<T'/T>_actual, <T'/T>_from_mass_ratio, <theta_lab>]
#        averaged over the colliding atoms in each stored step. NaN
#        in any row where no collision happened in that stored step.
#          * temperature_diagnostic (T, 3)
#   6 -- Tier-1a mass-dynamics. The mass-transfer channel now covers He
#        shedding (anchored_discrete), not only attachment, so:
#          * RENAME E_mass_attach_defect_eV -> E_mass_transfer_eV (same
#            (2N, T); positive for continuous-velocity shedding).
#          * ADD n_shell (2N, num_steps) -- per-atom integer He-shell count.
#          * ADD mass_scenario (scalar str) -- the run's mass scenario tag.
#          * DROP the mass_history_kg non-decreasing assumption (mass may
#            fall under anchored_discrete shedding).
#        Back-compat: load_ion_checkpoint migrates legacy v5 files (maps the
#        renamed field, synthesizes n_shell from mass_history_kg, defaults
#        mass_scenario='fixed') so existing v5 ion.npz run dirs still load.
#   7 -- Tier-2 internal-energy reservoir (Phase C Slice X):
#          * ADD E_int_eV (2N, T) -- per-atom internal energy [eV]; all-zero
#            under fixed / anchored_discrete, evolved by the biphasic
#            generative driver via the S1/S2/K1/K2 budget.
#        Back-compat: the migration shim is a stepwise cascade (v5->v6->v7);
#        the v6->v7 arm synthesizes an all-zero E_int_eV and emits a
#        UserWarning that the file predates the reservoir.
# ===========================================================================
_NEUTRAL_SCHEMA_VERSION: int = 2
_ION_SCHEMA_VERSION: int = 7


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
    * ``mass_history_kg``    : (2 * num_molecules, num_steps)  kg (mass over time; may rise via attachment/biphasic pickup OR fall via anchored_discrete/biphasic shedding -- not assumed monotone)
    * ``droplet_radii_angstrom``: (2 * num_molecules,)         Angstrom (per atom; same value for the two atoms of a molecule)
    * ``E_kin_eV``           : (2 * num_molecules, num_steps)  eV (per-atom)
    * ``E_pot_eV``           : (2 * num_molecules, num_steps)  eV (per-atom)
    * ``E_dissip_eV``        : (2 * num_molecules, num_steps)  eV (per-atom, cumulative)
    * ``E_mass_transfer_eV`` : (2 * num_molecules, num_steps) eV (per-atom, cumulative; mass-transfer bookkeeping, positive for Tier-1a continuous-velocity shedding; formerly ``E_mass_attach_defect_eV``)
    * ``E_int_eV``           : (2 * num_molecules, num_steps) eV (per-atom; the Tier-2 internal-energy reservoir added at schema v7. All-zero for ``fixed`` / ``anchored_discrete`` runs and for v6-origin files migrated on load; evolves only under the ``biphasic`` generative driver via the S1/S2/K1/K2 budget)
    * ``n_shell``            : (2 * num_molecules, num_steps) int-valued (per-atom He-shell count over time; constant under ``fixed``, the 21->14 staircase under ``anchored_discrete``, genuine channel-written state under ``biphasic``)
    * ``mass_scenario``      : scalar str                     (the run's mass scenario tag: ``fixed`` / ``anchored_discrete`` / ...)
    * ``b_ion_outside``      : (num_molecules,) bool           True if ion exited droplet
    * ``relative_loss_per_ps``: (2 * num_molecules, num_steps) 1/ps (per-atom energy loss rate)
    * ``number_of_collisions``: (2 * num_molecules, num_steps) int (cumulative, per-atom)
    * ``temperature_diagnostic``: (num_steps, 3)               eV-ratio, mass-ratio, rad
        Per-stored-step legacy MATLAB temperature diagnostic
        ``[<T'/T>_actual, <T'/T>_from_mass_ratio, <theta_lab>]``
        averaged over the atoms that collided in that step. ``NaN``
        in every column when no collision occurred. Mirrors the
        ``diagnostic_array`` accumulator in
        ``vmi_sim_3d_ion_propa.m:683``.

    Schema v6 differs from v5 by the Tier-1a mass-dynamics changes: it
    renames ``E_mass_attach_defect_eV`` to ``E_mass_transfer_eV`` (the
    channel now also covers shedding), adds ``n_shell`` and the
    ``mass_scenario`` tag, and drops the non-decreasing-mass assumption.
    Schema v7 (the current version) adds the Tier-2 ``E_int_eV`` internal-
    energy reservoir. Legacy v5 and v6 files are migrated on load
    (v5->v6->v7; see :func:`load_ion_checkpoint` and
    :func:`_migrate_ion_checkpoint`); a v6-origin file gets an all-zero
    synthesized ``E_int_eV`` plus a load-time warning. Pre-v5 files cannot
    be loaded.
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
    E_mass_transfer_eV: np.ndarray
    E_int_eV: np.ndarray
    n_shell: np.ndarray
    b_ion_outside: np.ndarray
    relative_loss_per_ps: np.ndarray
    number_of_collisions: np.ndarray
    temperature_diagnostic: np.ndarray
    mass_scenario: str = "fixed"
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
    """Load an IonCheckpoint from a .npz file.

    Legacy v5 and v6 checkpoints are migrated transparently on load via
    :func:`_migrate_ion_checkpoint` (a stepwise v5->v6->v7 cascade), so
    existing ``ion.npz`` run directories still load instead of failing the
    version check. The v6->v7 arm synthesizes an all-zero ``E_int_eV`` and
    emits a ``UserWarning`` that the file predates the internal-energy
    reservoir (strict-warning test/CI configs must expect or filter it);
    the v5->v6 arm is silent.
    """
    return _load_checkpoint(
        path,
        dataclass_type=IonCheckpoint,
        expected_version=_ION_SCHEMA_VERSION,
        cfg=cfg,
        migrate=_migrate_ion_checkpoint,
    )


def _migrate_ion_checkpoint(
    raw: dict[str, np.ndarray],
    version: int,
) -> tuple[dict[str, np.ndarray], int]:
    """Migrate a loaded ion checkpoint dict to the current schema.

    Migrations are applied **stepwise** so any surviving legacy file walks
    up to the current schema (v5 -> v6 -> v7); the arms are chained (each
    falls through to the next) rather than exclusive.

    **v5 -> v6** (the Tier-1a mass-dynamics bump):

    * maps ``E_mass_attach_defect_eV`` -> ``E_mass_transfer_eV`` (the
      mass-transfer channel now also covers shedding; the array is
      identical, only the name and sign-convention scope change);
    * synthesizes the absent ``n_shell`` from the already-present
      ``mass_history_kg`` via ``round((m/U - m_I+) / m_He)`` -- the same
      derivation the live driver uses, so writer and shim agree by
      construction. For a legacy fixed-mass run this is constant
      (~19 He at ``m_eff``);
    * defaults ``mass_scenario`` to ``"fixed"`` (every existing v5 run is
      a fixed / Tier-0 run).

    **v6 -> v7** (the Tier-2 internal-energy reservoir): synthesizes an
    all-zero ``E_int_eV`` (same shape as ``E_mass_transfer_eV``) because the
    file predates the reservoir, and emits a ``UserWarning`` recording that
    provenance. A file that predates ``E_int`` is 4-term-equivalent: an
    all-zero reservoir leaves the ledger closure unchanged. (The v5->v6 arm
    stays silent, as delivered; only the v6->v7 arm warns.)

    Any other version is returned unchanged (the caller's strict version
    check then raises). A genuine v7 file missing ``E_int_eV`` is therefore
    **not** repaired here -- it fails the ``missing fields`` check loudly
    rather than being silently zero-filled.
    """
    raw = dict(raw)
    if version == 5:
        if "E_mass_attach_defect_eV" in raw and "E_mass_transfer_eV" not in raw:
            raw["E_mass_transfer_eV"] = raw.pop("E_mass_attach_defect_eV")
        if "n_shell" not in raw and "mass_history_kg" in raw:
            mass_amu = np.asarray(raw["mass_history_kg"], dtype=float) / U
            raw["n_shell"] = np.rint((mass_amu - MASS_I_ION_AMU) / MASS_HE_AMU)
        if "mass_scenario" not in raw:
            raw["mass_scenario"] = np.asarray("fixed")
        raw["schema_version"] = np.asarray(6)
        version = 6
    if version == 6:
        if "E_int_eV" not in raw and "E_mass_transfer_eV" in raw:
            raw["E_int_eV"] = np.zeros_like(
                np.asarray(raw["E_mass_transfer_eV"], dtype=float)
            )
            warnings.warn(
                "ion checkpoint predates the Tier-2 E_int internal-energy "
                "reservoir (schema v6); synthesizing an all-zero E_int_eV. "
                "The energy-ledger closure for this file stays 4-term-"
                "equivalent.",
                UserWarning,
                stacklevel=2,
            )
        raw["schema_version"] = np.asarray(7)
        version = 7
    return raw, version


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
    migrate: Any = None,
) -> Any:
    """Read a .npz file into a checkpoint dataclass.

    ``migrate``, if given, is a callable ``(raw_dict, version) ->
    (raw_dict, version)`` applied before the strict version check, so an
    older on-disk schema can be upgraded in-memory to ``expected_version``
    (the ion v5->v6->v7 back-compat cascade). It receives a mutable copy of the
    loaded arrays and must return the migrated arrays plus the new version.
    """
    p = Path(path)
    if p.suffix != ".npz":
        p = p.with_suffix(".npz")
    if not p.exists():
        raise FileNotFoundError(f"checkpoint not found: {p}")

    # Materialize the .npz into a mutable dict so any migration can rename
    # / synthesize fields after the file handle is closed.
    with np.load(p, allow_pickle=False) as npz:
        if "schema_version" not in npz.files:
            raise ValueError(
                f"checkpoint at {p} has no schema_version field; "
                "it was written by an older version. Re-run the simulation."
            )
        version = int(npz["schema_version"])
        raw: dict[str, Any] = {k: npz[k] for k in npz.files}

    # 1. Migrate (if a shim is supplied) then enforce the version.
    if version != expected_version and migrate is not None:
        raw, version = migrate(raw, version)
    if version != expected_version:
        raise ValueError(
            f"checkpoint at {p} has schema_version={version}, "
            f"this code expects {expected_version}. Re-run the simulation."
        )

    # 2. Build kwargs for dataclass construction
    expected_fields = {f.name for f in fields(dataclass_type)}
    missing = expected_fields - set(raw.keys())
    if missing:
        raise ValueError(
            f"checkpoint at {p} is missing fields: {sorted(missing)}"
        )

    kwargs: dict[str, Any] = {}
    for f in fields(dataclass_type):
        arr = raw[f.name]
        # Unwrap 0-d arrays for the scalar fields, matched by name: under
        # ``from __future__ import annotations`` ``f.type`` is a string, so
        # type-based dispatch is impossible here. These are the schemas'
        # only non-array fields.
        if f.name in ("schema_version", "num_molecules"):
            kwargs[f.name] = int(arr)
        elif f.name == "mass_scenario":
            kwargs[f.name] = str(arr)
        else:
            kwargs[f.name] = np.asarray(arr)

    instance = dataclass_type(**kwargs)

    # 3. Optional shape validation against cfg
    if cfg is not None:
        _validate_against_cfg(instance, cfg, path=p)

    return instance


# ===========================================================================
# Continuation-stage boot helpers (shared by relaxation_stage / detection_stage
# -- review fix 2026-07-07: the two stages carried verbatim copies)
# ===========================================================================
def check_biphasic_seed_checkpoint(ckpt: Any, *, stage: str) -> None:
    """Seed-coherence guards for a stage that boots from a checkpoint column.

    Shared by ``run_relaxation_stage`` and ``run_detection_stage`` (rule 1 --
    one source for the guard logic). Two fail-loud checks:

    1. **Scenario:** the seed must be a ``biphasic`` checkpoint -- only that
       scenario carries the physical ``E_int``/``n_shell`` cascade state the
       continuation stages propagate (a v6-migrated file carries synthesized
       all-zero ``E_int_eV``, reading as trivially frozen).
    2. **Stride:** ``mass_history_kg[:, -1] == mass_final_kg`` -- the stage
       seeds from the last *stored* column, but a strided run whose
       allocation was exactly consumed leaves the true final state only in
       the ``*_final_*`` fields, and ``E_int_eV``/``n_shell`` have no such
       fields to recover from. ``mass_final_kg`` is written from the true
       final state, so any mass event in the dropped tail is caught here.

    Parameters
    ----------
    ckpt : IonCheckpoint
        The seed checkpoint (``ion.npz`` or ``relaxation.npz`` -- one
        contract, both are v7 ``IonCheckpoint`` s).
    stage : str
        The calling stage's name for the error message (e.g.
        ``"run_relaxation_stage"``).

    Raises
    ------
    ValueError
        On a non-biphasic scenario or a stale final stored column.
    """
    if ckpt.mass_scenario != "biphasic":
        raise ValueError(
            f"{stage} requires a biphasic seed checkpoint; got "
            f"mass_scenario={ckpt.mass_scenario!r}. The continuation stages "
            "propagate the biphasic E_int/n_shell cascade state, which other "
            "scenarios do not carry."
        )
    if not np.array_equal(
        np.asarray(ckpt.mass_history_kg)[:, -1], np.asarray(ckpt.mass_final_kg)
    ):
        raise ValueError(
            f"{stage} seeds from the seed checkpoint's final stored column, "
            "but mass_history_kg[:, -1] != mass_final_kg: the upstream run "
            "was stored with a stride that dropped the true final state "
            "(E_int_eV/n_shell have no *_final_* fields to recover from). "
            "Re-run the upstream stage with a larger max_bytes so the final "
            "state lands in the last stored column."
        )


def stage_stream_rng(seed: Any, stream_key: int) -> np.random.Generator:
    """Derive a continuation stage's private PCG64 stream (rule-1 single source).

    ``SeedSequence((seed, stream_key))`` with a fixed module-level key per
    stage (``RELAXATION_STREAM_KEY`` / ``DETECTION_STREAM_KEY``) -- the
    ion-stage stream is never touched and the stages never replay each
    other. ``seed=None`` -> a non-reproducible default generator (mirrors
    the driver convention).
    """
    if seed is None:
        return np.random.default_rng()
    return np.random.default_rng(np.random.SeedSequence((int(seed), stream_key)))


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
        "E_mass_transfer_eV", "E_int_eV", "n_shell",
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

    # temperature_diagnostic is (num_steps, 3) -- different leading dim
    # than the (2N, num_steps) trajectory arrays, so it's checked here.
    if hasattr(checkpoint, "temperature_diagnostic"):
        td = checkpoint.temperature_diagnostic
        if td.ndim != 2 or td.shape[1] != 3:
            raise ValueError(
                f"checkpoint temperature_diagnostic has shape {td.shape}, "
                f"expected (num_steps, 3)"
            )
        if num_steps_seen is not None and td.shape[0] != num_steps_seen:
            raise ValueError(
                f"checkpoint temperature_diagnostic has num_steps={td.shape[0]}, "
                f"but trajectory arrays had num_steps={num_steps_seen}"
            )
