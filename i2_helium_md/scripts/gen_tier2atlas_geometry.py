"""Atlas Axis A / geometry-correction stage G1 — the droplet-geometry grid.

Runs the `TIER2_SENSITIVITY_ATLAS_PLAN` §3.1 grid: **droplet radius × birth
law**, one fixed droplet size per cell (size *sampling* is switched off — the
radius is the controlled variable), everything else at the standing production
point `finc1v725`.

**11 cells, N = 500, one shared seed** (common random numbers, so cross-cell
differences are partly paired — plan §3.1 resolution pre-registration):

    R axis  R1 26.6 Å (N 1727, the standing point's realized mean)
            R2 34.0 Å (N 3605, the parent document's smallest droplet)
            R3 49.4 Å (N 11059, the parent ensemble's raw mean)
            R4 68.3 Å (N 29227, the parent's largest — L2/L3 only)

    law     L1 center-pin      (`single_initial_position=True`)
            L2 parent Boltzmann (`birth_position_law="boltzmann"`,
                                 well 313.2 K = the parent's DFT-fit β₂)
            L3 uniform_volume m 3 Å (the standing production law)

`R1/L3` is the grid's **internal reference cell**: the standing geometry at a
fixed R instead of a sampled one.

Three config guards are exercised deliberately (plan §3.2) and are verified by
`--dry-run` before any MD:

1. T8: the analytic prior arms are refused under `use_single_droplet_size=True`,
   so every cell runs `droplet_size_prior="legacy"` (the standing point runs the
   analytic `kornilov_lognormal` arm — this is a deliberate difference);
2. `uniform_volume` requires `single_initial_position=False` (L3);
3. `birth_position_law="boltzmann"` requires
   `initial_position_margin_angstrom == 0.0` (L1, L2).

**Scoring rule (plan §3.1, from the Axis A pre-read).** The landing is a
*mixture* property — the pooled battery's W₁ beats every single geometry bin —
so these fixed-geometry cells are compared to **each other** and to a
re-weighted mixture reconstruction, never directly against the committed
acceptance. E_bind stays at the standing 0.1168 eV in every cell (G0-3: its
R-dependence is bounded ≤ 0.009 eV over this span, D0 §9.1).

Atlas stance: instrument runs — nothing moves `finc1v725`.

Usage::

    python scripts/gen_tier2atlas_geometry.py --dry-run   # cfg + guard check
    python scripts/gen_tier2atlas_geometry.py             # all 11 cells
    python scripts/gen_tier2atlas_geometry.py r1l3 r3l2   # selected cells
"""

from __future__ import annotations

import argparse
import dataclasses
import os
from pathlib import Path
import sys
from typing import NamedTuple

for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass


# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"
VARIANT = "shared_pure_cubic"
N = 500                         # fragments per cell (plan §3.1)
SEED = 20260727                 # ONE seed shared by every cell (CRN)
ION_TIME_PS = 30.0
DT_ION_PS = 0.01

# Standing production point finc1v725 (verified against bigc1v725s1/cfg.json).
V_C_APS = 7.25
TAU_PS = 3.2
E0_EV = 0.27
P_TAIL = -1.0

BUDGET_EV = 2.70
S_EFF = 8.0
F_RET = 0.1
LAMBDA0_PER_PS = 0.9
COOLING_GATE = "density_scaled"

# Atlas plan §3.5b (2026-07-27): the anchored-radius cells leave a third to a
# half of their ions still helium-coupled at handover, which the "exclude"
# policy refuses. "exclude_all_coupled" counts that class instead of
# integrating it, keeping bound (physics) and marginal (modelling exclusion)
# decomposed. The six R <= 34 A cells have an EMPTY marginal class -- every
# violator there was provably bound -- so re-running their detection under
# this policy must reproduce them bit-for-bit (the item-10 arm oracle,
# enforced by _run_detection_only below).
RETAINED_POLICY = "exclude_all_coupled"
SHED_CONVENTION = "co_moving"
SHELL_MODEL = "density_tied"
PARTITION_LAW = "sigma_proportional"
R0_GS_ANGSTROM = 2.666
E_COULOMB_SCALE = 1.0
RELAXATION_TIME_PS = 8000.0
RELAXATION_FORCES = "coulomb"
RELAXATION_DISSIPATION = "landau_gated_drag"
V_LIMIT_M_PER_S = 58.0

# Axis A grid. Sizes are He counts; the realized radius is the bulk conversion
# R = 2.2173*N^(1/3) (D0 §15.3 — the only convention on the propagation path).
R_AXIS: dict[str, tuple[int, float]] = {     # key -> (N_He, target R [Å])
    "r1": (1727, 26.6),
    "r2": (3605, 34.0),
    "r3": (11059, 49.4),
    "r4": (29227, 68.3),
}
R_TARGET_TOL_ANGSTROM = 0.05

# Parent document's own well depth behind the Boltzmann law: β₂ = 26.99 meV
# (D0 §15.4). The config default 573.3 K stays untouched (atlas G0-2).
PARENT_WELL_K = 313.2

# Cells: full 3 × 3 on r1..r3, plus the two funded r4 cells on L2/L3.
CELL_KEYS: tuple[str, ...] = (
    "r1l1", "r1l2", "r1l3",
    "r2l1", "r2l2", "r2l3",
    "r3l1", "r3l2", "r3l3",
    "r4l2", "r4l3",
)

# Reference run for the cfg-diff guard: a standing-point battery member.
REFERENCE_RUN = "9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s1"

# Keys every cell must differ in vs the reference (N, seed, fixed-size branch).
_BASE_DIFF_KEYS = frozenset({
    "num_molecules", "seed",
    "use_single_droplet_size", "single_droplet_size", "droplet_size_prior",
    # §3.5b: the reference battery member runs "exclude"; this grid needs the
    # marginal-counting arm at the anchored radii (and runs it everywhere, so
    # the grid stays one homogeneous configuration).
    "detection_droplet_retained_policy",
})
# Additional per-law diffs vs the reference (whose law IS L3).
_LAW_DIFF_KEYS: dict[str, frozenset[str]] = {
    "l1": frozenset({"birth_position_law", "initial_position_margin_angstrom",
                     "single_initial_position"}),
    "l2": frozenset({"birth_position_law", "initial_position_margin_angstrom",
                     "binding_energy_molecule_K"}),
    "l3": frozenset(),
}

DEFAULT_CONCURRENCY = 3
SKIP_COMPLETED_RUNS = True
# §3.5b recovery: complete the (or re-verify a) cell from its stored
# relaxation.npz, running ONLY the detection stage -- zero new MD. This is how
# the five anchored-radius cells are finished under the new retained policy,
# and how the six already-scored cells are re-verified against it.
RESUME_DETECTION = True
# Safety catch for the arm oracle: with this False a detection re-run that
# does NOT reproduce an existing detection.npz refuses to write. Only set it
# True deliberately, when the result is *expected* to change.
ALLOW_DETECTION_OVERWRITE = False
# Resume semantics: complete cells are skipped untouched by SKIP_COMPLETED_RUNS
# above, so overwriting only ever rebuilds a cell that died mid-stage (the
# 2026-07-26 first launch was killed with three cells in the E2 relaxation
# stage). Set False for a strictly-once launch.
OVERWRITE_EXISTING_RUN = True


# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2_md_confirmation import rq4graded_rungs_eV  # noqa: E402
from scripts.tier0_common import run_dir_name  # noqa: E402
from scripts.tier2_common import build_biphasic_cfg, cfg_diff_vs_reference  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    read_confirmation_detection,
)
from i2_helium_md.simulation.checkpoint import load_ion_checkpoint  # noqa: E402
from i2_helium_md.simulation.detection_stage import (  # noqa: E402
    load_detection_result,
    run_detection_stage,
    save_detection_result,
)
from i2_helium_md.simulation.initial_state import (  # noqa: E402
    droplet_radius_bulk_angstrom,
)
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

DETECTION_TIME_PS = 8.53e6
_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "cfg.json", "neutral.npz", "ion.npz", "relaxation.npz", "detection.npz",
)


class GeometryCell(NamedTuple):
    """One Axis A grid cell."""

    label: str          # e.g. "r1l3"
    r_key: str          # "r1".."r4"
    law_key: str        # "l1".."l3"
    droplet_N_he: int
    target_R_angstrom: float

    @property
    def realized_R_angstrom(self) -> float:
        """Bulk-conversion radius of this cell's droplet [Å]."""
        return float(droplet_radius_bulk_angstrom(np.array([self.droplet_N_he]))[0])


def _matrix() -> tuple[GeometryCell, ...]:
    cells = []
    for key in CELL_KEYS:
        r_key, law_key = key[:2], key[2:]
        n_he, target_R = R_AXIS[r_key]
        cells.append(GeometryCell(key, r_key, law_key, n_he, target_R))
    return tuple(cells)


GEOMETRY_MATRIX = _matrix()
_SPEC_BY_LABEL = {c.label: c for c in GEOMETRY_MATRIX}


def atlas_run_dir_name(label: str) -> str:
    """Run-dir basename in the atlas namespace (no ``_tier2_`` substring)."""
    name = run_dir_name(
        CASE, VARIANT, N,
        run_tag=f"tier2atlas_conf{int(round(BUDGET_EV * 100))}_geo{label}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_cell(cell: GeometryCell) -> SimConfig:
    """Build one cell's validated config (standing point + this geometry)."""
    if cell.law_key == "l1":            # center-pin: r0 sampled then zeroed
        birth_law, margin, single_pos = "boltzmann", 0.0, True
        well_K = None                   # law inert at r0 == 0; keep the default
    elif cell.law_key == "l2":          # the parent document's thermal law
        birth_law, margin, single_pos = "boltzmann", 0.0, False
        well_K = PARENT_WELL_K
    elif cell.law_key == "l3":          # the standing production law
        birth_law, margin, single_pos = "uniform_volume", 3.0, False
        well_K = None
    else:
        raise ValueError(f"unknown law key {cell.law_key!r}")

    cfg = build_biphasic_cfg(
        CASE, VARIANT,
        num_molecules=N, ion_time_ps=ION_TIME_PS, dt_ion_ps=DT_ION_PS,
        seed=SEED, lambda0_per_ps=LAMBDA0_PER_PS,
        f_int=E0_EV / BUDGET_EV, f_ret=F_RET, tau_ps=TAU_PS,
        evap_rrk_dof=S_EFF, coulomb_available_eV=BUDGET_EV,
        relaxation_time_ps=RELAXATION_TIME_PS, relaxation_forces=RELAXATION_FORCES,
        cooling_spatial_gate=COOLING_GATE, drag_form="capped_cubic",
        drag_coefficient_overrides={"v_c": V_C_APS, "p_tail": P_TAIL},
        dissociation_ladder="tabulated", tabulated_ladder_rungs_eV=rq4graded_rungs_eV(),
        R0_GS_angstrom=R0_GS_ANGSTROM, E_coulomb_scale=E_COULOMB_SCALE,
        single_initial_position=single_pos,
        detection_time_ps=DETECTION_TIME_PS, birth_position_law=birth_law,
        initial_position_margin_angstrom=margin,
        detection_droplet_retained_policy=RETAINED_POLICY,
        evaporation_shed_convention=SHED_CONVENTION, initial_shell_model=SHELL_MODEL,
        internal_energy_partition_law=PARTITION_LAW,
        # Fixed droplet size: size sampling OFF, so the analytic prior must go
        # back to 'legacy' (T8 guard) -- the controlled-variable requirement.
        droplet_size_prior="legacy",
        use_single_droplet_size=True,
        single_droplet_size=cell.droplet_N_he,
        binding_energy_molecule_K=well_K,
    )
    cfg = dataclasses.replace(
        cfg,
        relaxation_dissipation=RELAXATION_DISSIPATION,
        v_limit_m_per_s=V_LIMIT_M_PER_S,
    )
    cfg.validate()
    return cfg


def verify_geometry(cfg: SimConfig, cell: GeometryCell) -> None:
    """Check the realized radius matches the cell's target radius."""
    realized = cell.realized_R_angstrom
    if abs(realized - cell.target_R_angstrom) > R_TARGET_TOL_ANGSTROM:
        raise AssertionError(
            f"[{cell.label}] N = {cell.droplet_N_he} He -> R = {realized:.3f} Å, "
            f"off the grid target {cell.target_R_angstrom} Å by more than "
            f"{R_TARGET_TOL_ANGSTROM} Å (bulk convention R = 2.2173*N^(1/3))."
        )
    if cfg.single_droplet_size != cell.droplet_N_he:
        raise AssertionError(f"[{cell.label}] single_droplet_size not stamped.")
    if not cfg.use_single_droplet_size or cfg.droplet_size_prior != "legacy":
        raise AssertionError(
            f"[{cell.label}] droplet size must be the fixed legacy branch; got "
            f"use_single_droplet_size={cfg.use_single_droplet_size}, "
            f"droplet_size_prior={cfg.droplet_size_prior!r}."
        )
    # E_bind stays the standing jointly-extracted value (atlas G0-3).
    stamped = cfg.drag_coefficients.effective_binding_energy_I_ion_eV
    if cfg.binding_energy_I_ion_eV != stamped:
        raise AssertionError(
            f"[{cell.label}] E_bind {cfg.binding_energy_I_ion_eV!r} != the "
            f"bundle stamp {stamped!r}; the geometry grid must not break the "
            f"§6.5.1 pairing (no unvalidated-binding hatch on this axis)."
        )


def verify_against_reference(cfg: SimConfig, cell: GeometryCell) -> list[str]:
    """Diff vs the standing-point reference; only geometry keys may differ.

    Returns the sorted diff key list (so ``--dry-run`` can print it). The
    field-set bookkeeping — including tolerating ``SimConfig`` fields added
    after the reference run was written, at their defaults only — lives in the
    shared :func:`scripts.tier2_common.cfg_diff_vs_reference`.
    """
    ref_path = PROJECT_ROOT / "data" / "runs" / REFERENCE_RUN / "cfg.json"
    diff = cfg_diff_vs_reference(cfg, ref_path, context=cell.label)
    expected = _BASE_DIFF_KEYS | _LAW_DIFF_KEYS[cell.law_key]
    if set(diff) != expected:
        raise AssertionError(
            f"[{cell.label}] cfg diff vs {REFERENCE_RUN} is {diff}; expected "
            f"exactly {sorted(expected)} (N + seed + the fixed-size branch + "
            f"this cell's birth-law keys). An unexpected key means a knob "
            f"moved off the standing point."
        )
    return diff


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


# Artifacts that make a cell resumable at the detection step: everything the
# MD produced. A cell holding these needs ZERO new MD to be completed.
_PRE_DETECTION_ARTIFACTS: tuple[str, ...] = (
    "cfg.json", "neutral.npz", "ion.npz", "relaxation.npz",
)

# Which cfg fields the detection-only path is allowed to change relative to a
# stored run's own cfg.json. Anything else means the MD on disk was produced
# under different physics and must not be reused.
_DETECTION_ONLY_DIFF_KEYS = frozenset({"detection_droplet_retained_policy"})


def _detection_results_identical(a, b) -> tuple[bool, str]:
    """Compare two ``DetectionResult`` objects field by field (bit-exact).

    Returns ``(identical, first_difference)``; ``first_difference`` is ``""``
    when they match. Used as the §3.5b item-10 arm oracle: at R <= 34 A the
    marginal class is empty, so the new policy must be a no-op there.
    """
    for field in dataclasses.fields(a):
        va, vb = getattr(a, field.name), getattr(b, field.name)
        if isinstance(va, np.ndarray) or isinstance(vb, np.ndarray):
            va, vb = np.asarray(va), np.asarray(vb)
            if va.shape != vb.shape:
                return False, f"{field.name}: shape {va.shape} != {vb.shape}"
            if va.dtype.kind in "US":
                if not np.array_equal(va, vb):
                    return False, f"{field.name}: string arrays differ"
            elif not np.array_equal(va, vb, equal_nan=True):
                bad = int(np.count_nonzero(va != vb))
                return False, f"{field.name}: {bad} element(s) differ"
        elif va != vb:
            return False, f"{field.name}: {va!r} != {vb!r}"
    return True, ""


def _run_detection_only(cell: GeometryCell, cfg: SimConfig, run_dir: Path) -> str:
    """Complete (or re-verify) a cell from its stored ``relaxation.npz``.

    The §3.5b zero-MD recovery path. The stored cfg is diffed against the
    freshly built one first: only ``detection_droplet_retained_policy`` may
    have moved, otherwise the on-disk MD was produced under different physics
    and reusing it would be silently wrong.

    When ``detection.npz`` already exists the new result is compared against
    it and the file is **only** rewritten if the two are bit-identical (so the
    rewrite is a no-op). A differing result is reported and left unwritten
    unless ``ALLOW_DETECTION_OVERWRITE`` is set -- that is the arm oracle, not
    a formality: it is what proves the new policy changed nothing at R <= 34 A.
    """
    stored_diff = cfg_diff_vs_reference(
        cfg, run_dir / "cfg.json", context=f"{cell.label} (stored)"
    )
    if not set(stored_diff) <= _DETECTION_ONLY_DIFF_KEYS:
        raise AssertionError(
            f"[{cell.label}] stored cfg differs from the rebuilt cfg in "
            f"{sorted(set(stored_diff) - _DETECTION_ONLY_DIFF_KEYS)}, not just "
            f"the detection policy — the MD on disk was produced under "
            f"different physics and must not be reused. Delete the run dir to "
            f"recompute it from scratch."
        )
    seed_ckpt = load_ion_checkpoint(run_dir / "relaxation.npz")
    print(f"[{cell.label}] detection-only (zero MD) from relaxation.npz "
          f"[cfg diff: {stored_diff or 'none'}] ...", flush=True)
    detect = run_detection_stage(seed_ckpt, cfg, save_path=None)

    det_path = run_dir / "detection.npz"
    if det_path.exists():
        previous = load_detection_result(det_path)
        identical, difference = _detection_results_identical(detect, previous)
        if identical:
            print(f"[{cell.label}] ORACLE OK: detection reproduces the stored "
                  f"result bit-for-bit under {cfg.detection_droplet_retained_policy!r}.",
                  flush=True)
        elif not ALLOW_DETECTION_OVERWRITE:
            raise AssertionError(
                f"[{cell.label}] the new detection result DIFFERS from the "
                f"stored one ({difference}) and ALLOW_DETECTION_OVERWRITE is "
                f"False, so nothing was written. At R <= 34 A this must not "
                f"happen (the marginal class is empty there); investigate "
                f"before overwriting."
            )
        else:
            print(f"[{cell.label}] OVERWRITING a differing detection result "
                  f"({difference}).", flush=True)

    RunDirectory(run_dir).save_cfg(cfg)
    save_detection_result(detect, det_path)
    read = read_confirmation_detection(detect, label=cell.label)
    print(f"[{cell.label}] done (zero MD) -> {run_dir} "
          f"(scored={read.num_scored}/{read.num_ions}, "
          f"trap_bound={read.trap_bound_frac:.3f}, "
          f"trap_marginal={read.trap_marginal_frac:.3f}, "
          f"n_mean={read.n_mean:.2f})", flush=True)
    return cell.label


def _run_one(cell: GeometryCell) -> str:
    cfg = build_cell(cell)
    verify_geometry(cfg, cell)
    verify_against_reference(cfg, cell)
    run_dir = PROJECT_ROOT / "data" / "runs" / atlas_run_dir_name(cell.label)
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir) and not RESUME_DETECTION:
            print(f"[{cell.label}] skip (complete) -> {run_dir}", flush=True)
            return cell.label
        if RESUME_DETECTION and all(
            (run_dir / name).exists() for name in _PRE_DETECTION_ARTIFACTS
        ):
            return _run_detection_only(cell, cfg, run_dir)
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(f"{run_dir} exists and is not complete.")
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{cell.label}] R={cell.realized_R_angstrom:.2f} Å "
          f"(N={cell.droplet_N_he} He), law={cfg.birth_position_law}"
          f"{'/center-pin' if cfg.single_initial_position else ''}, "
          f"margin={cfg.initial_position_margin_angstrom} Å, "
          f"well={cfg.binding_energy_molecule_K} K, seed={cfg.seed}", flush=True)
    print(f"[{cell.label}] neutral ...", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{cell.label}] ion ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    print(f"[{cell.label}] relaxation (E2) ...", flush=True)
    relax = run_relaxation_stage(ion, cfg, save_path=run.root / "relaxation.npz")
    print(f"[{cell.label}] detection ...", flush=True)
    detect = run_detection_stage(
        relax.checkpoint, cfg, save_path=run.root / "detection.npz")
    dm = detect.detected_mask
    nd = float(detect.n_detected[dm].mean()) if dm.any() else float("nan")
    print(f"[{cell.label}] done -> {run_dir} (n_detect_mean={nd:.2f}, "
          f"droplet_retained={int(np.count_nonzero(~dm))}/{dm.size})", flush=True)
    return cell.label


def _run_one_by_label(label: str) -> str:
    return _run_one(_SPEC_BY_LABEL[label])


def _select(labels: list[str]) -> list[GeometryCell]:
    if not labels:
        return list(GEOMETRY_MATRIX)
    unknown = sorted(set(labels) - set(_SPEC_BY_LABEL))
    if unknown:
        raise SystemExit(f"unknown label(s) {unknown}; expected among "
                         f"{[c.label for c in GEOMETRY_MATRIX]}")
    return [_SPEC_BY_LABEL[label] for label in labels]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cells = _select(args.labels)
    print(f"Atlas Axis A (G1) geometry grid: {len(cells)} cell(s), N={N}, "
          f"seed={SEED} (shared), R keys "
          f"{sorted(set(c.r_key for c in cells))}, laws "
          f"{sorted(set(c.law_key for c in cells))}, concurrency="
          f"{args.concurrency}{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument runs -- standing point does not move).", flush=True)

    if args.dry_run:
        for cell in cells:
            cfg = build_cell(cell)
            verify_geometry(cfg, cell)
            diff = verify_against_reference(cfg, cell)
            print(f"[{cell.label}] cfg OK R={cell.realized_R_angstrom:.2f} Å "
                  f"law={cfg.birth_position_law}"
                  f"{'+center' if cfg.single_initial_position else ''} "
                  f"margin={cfg.initial_position_margin_angstrom} "
                  f"well={cfg.binding_energy_molecule_K} K -> "
                  f"{atlas_run_dir_name(cell.label)}", flush=True)
            print(f"[{cell.label}]   diff vs standing: {diff}", flush=True)
        print("dry run complete: all cfgs pass the geometry + "
              "standing-point diff guards.", flush=True)
        return 0

    concurrency = max(1, args.concurrency)
    if concurrency == 1 or len(cells) == 1:
        for cell in cells:
            _run_one(cell)
        return 0

    from concurrent.futures import ProcessPoolExecutor, as_completed
    workers = min(concurrency, len(cells))
    print(f"launching {len(cells)} cells through a {workers}-slot pool ...",
          flush=True)
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as poolx:
        futures = {poolx.submit(_run_one_by_label, c.label): c.label for c in cells}
        for fut in as_completed(futures):
            label = futures[fut]
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[{label}] FAILED: {exc!r}", flush=True)
    if failures:
        print(f"{failures}/{len(cells)} cell(s) failed.", flush=True)
        return 1
    print(f"all {len(cells)} cells complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
