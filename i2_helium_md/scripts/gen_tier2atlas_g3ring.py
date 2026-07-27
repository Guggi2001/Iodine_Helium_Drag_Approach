"""Atlas G3 Step 3 — the MD confirmation ring at the corrected geometry.

Runs the `TIER2_SENSITIVITY_ATLAS_PLAN` §3.5d ring: **14 cells × N = 500,
one fresh shared seed**, at the G2-adopted corrected droplet geometry
(`legacy`+`raw` sampled sizes at the preset's own source conditions,
Boltzmann births at the 313.2 K parent well), confirming or refuting the
§3.5c twin basin in real MD.

Arms (plan §3.5d; per-cell knobs = (v_c, E_bind, τ, E₀), everything else
the standing finc1v725 pins):

* **A** a031/a033/a035/a037 — E₀ ladder through twin sub-basin 1
  (v_c 5.5, τ 4.8), spanning the whole n̄ bias bracket downward;
* **B** b029/b031/b033 — E₀ ladder through sub-basin 2 (v_c 6.0, τ 6.4);
* **C** c50/c65 — twin-fail edge controls (§6.6 off-needle discipline);
* **D** d030/d036 — the τ-cross cells (v_c 5.5 @ τ 6.4 near-gate;
  v_c 6.0 @ τ 4.8 twin-fail), separating τ from v_c per sub-basin;
* **E** e0482/e154 — E_bind bracket at the arm-A center, run under the
  documented §6.5.1 pairing hatch (`allow_unvalidated_binding_pairing`);
* **F** f725 — the standing point at the corrected geometry: the
  broken-landing continuity anchor (twin n̄ 17.0, n₁ 0).

**GR-P1 oracle (§1.4, enforced before any MD and by --dry-run):** the 14
frozen twin rows below must reproduce **string-exact** from the committed
`h2b_g3scan_predictions.csv`, and every cell's cfg must diff against the
standing battery reference in exactly its pre-registered key set.

Retained policy: `exclude_all_coupled` (interim per the §3.5d
adjudication; bound/marginal decomposition read by the ring scorer,
`scripts/post_processing/tier2atlas_g3ring_table.py`).

Atlas stance: instrument runs — nothing here moves `finc1v725`.

Usage::

    python scripts/gen_tier2atlas_g3ring.py --dry-run   # oracles + cfg guards
    python scripts/gen_tier2atlas_g3ring.py             # all 14 cells
    python scripts/gen_tier2atlas_g3ring.py a037 f725   # selected cells
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import os
from pathlib import Path
import sys
from typing import NamedTuple, Optional

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
N = 500                         # fragments per cell (plan §3.5d)
SEED = 20260728                 # ONE fresh seed shared by every cell (CRN)
ION_TIME_PS = 30.0
DT_ION_PS = 0.01

BUDGET_EV = 2.70
S_EFF = 8.0
F_RET = 0.1
LAMBDA0_PER_PS = 0.9
COOLING_GATE = "density_scaled"
P_TAIL = -1.0

# The G2-adopted corrected geometry (G0 spec — no new numbers): legacy
# sampled sizes read RAW from the nozzle correlation at the preset's own
# source conditions (<N> = 12794), Boltzmann births at the parent well.
SIZE_PRIOR = "legacy"
SIZE_SAMPLER_MODE = "raw"
BIRTH_LAW = "boltzmann"
BIRTH_MARGIN_ANGSTROM = 0.0
PARENT_WELL_K = 313.2           # per-run override; config default untouched

# §3.5d adjudication: exclude_all_coupled stays the INTERIM policy — the
# scorer reads the bound/marginal decomposition; the final call is G4's.
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
DETECTION_TIME_PS = 8.53e6

# Exact E_bind values behind the arm-E tags (§3.5c convention: the plan's
# display values resolved to their provenance-exact numbers).
EBIND_EV = {"eb0482": 0.0482, "eb154": 0.154}

# Reference run for the cfg-diff guard: a standing-point battery member.
REFERENCE_RUN = "9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s1"

# Keys every cell must differ in vs the reference (N + seed + the corrected
# geometry + the interim retained policy). `droplet_size_sampler_mode` is NOT
# here: the reference cfg.json predates the field, so the diff helper cannot
# verify a non-default value against it — the diff runs on a copy at the
# default and the real "raw" value is guarded by verify_corrected_geometry.
_BASE_DIFF_KEYS = frozenset({
    "num_molecules", "seed",
    "droplet_size_prior",
    "birth_position_law", "initial_position_margin_angstrom",
    "binding_energy_molecule_K",
    "detection_droplet_retained_policy",
})
_KNOB_KEYS = {
    "v_c": "drag_coefficients",
    "tau": "internal_energy_cooling_tau_ps",
    "e0": "internal_energy_partition_fraction",
}
_EBIND_KEYS = frozenset({
    "binding_energy_I_ion_eV", "allow_unvalidated_binding_pairing",
})

# Standing finc1v725 knob values (a knob equal to these leaves no cfg diff).
STANDING_V_C = 7.25
STANDING_TAU_PS = 3.2
STANDING_E0_EV = 0.27

# The committed twin-scan record the GR-P1 oracle re-reads.
TWIN_SCAN_CSV = Path("data/runs/h2b_forward_model/h2b_g3scan_predictions.csv")

DEFAULT_CONCURRENCY = 3
SKIP_COMPLETED_RUNS = True
# Resume semantics (the G1-grid precedent): complete cells are skipped
# untouched by SKIP_COMPLETED_RUNS, so overwriting only ever rebuilds a cell
# that died mid-stage — a relaunch after a crash is safe and minimal.
OVERWRITE_EXISTING_RUN = True


class RingCell(NamedTuple):
    """One §3.5d ring cell."""

    label: str
    v_c: float
    eb_tag: str                 # "eb1168" (bundle pairing) | arm-E tags
    tau_ps: float
    e0_eV: float
    role: str


RING_MATRIX: tuple[RingCell, ...] = (
    RingCell("a031", 5.5, "eb1168", 4.8, 0.31, "arm A: E0 ladder, sub-basin 1"),
    RingCell("a033", 5.5, "eb1168", 4.8, 0.33, "arm A: E0 ladder, sub-basin 1"),
    RingCell("a035", 5.5, "eb1168", 4.8, 0.35, "arm A: E0 ladder, sub-basin 1"),
    RingCell("a037", 5.5, "eb1168", 4.8, 0.37, "arm A: twin-gated top"),
    RingCell("b029", 6.0, "eb1168", 6.4, 0.29, "arm B: E0 ladder, sub-basin 2"),
    RingCell("b031", 6.0, "eb1168", 6.4, 0.31, "arm B: E0 ladder, sub-basin 2"),
    RingCell("b033", 6.0, "eb1168", 6.4, 0.33, "arm B: twin-gated top"),
    RingCell("c50", 5.0, "eb1168", 4.8, 0.34, "arm C: low-v_c edge control"),
    RingCell("c65", 6.5, "eb1168", 6.4, 0.33, "arm C: high-v_c edge control"),
    RingCell("d030", 5.5, "eb1168", 6.4, 0.30, "arm D: tau cross (sub-basin 1)"),
    RingCell("d036", 6.0, "eb1168", 4.8, 0.36, "arm D: tau cross (sub-basin 2)"),
    RingCell("e0482", 5.5, "eb0482", 4.8, 0.36, "arm E: shallow-well bracket"),
    RingCell("e154", 5.5, "eb154", 4.8, 0.36, "arm E: deep-well bracket"),
    RingCell("f725", STANDING_V_C, "eb1168", STANDING_TAU_PS, STANDING_E0_EV,
             "arm F: standing-point continuity anchor"),
)
_SPEC_BY_LABEL = {c.label: c for c in RING_MATRIX}

# The frozen twin rows (plan §3.5d table; committed h2b_g3scan_predictions.csv
# strings, compared EXACTLY): label -> (trapped_frac, suppressed_frac,
# nbar_det, n1_solv, w1_solv, gate).
TWIN_ROWS: dict[str, tuple[str, str, str, str, str, str]] = {
    "a031": ("0.0535", "0.0008", "6.717", "0.0371", "1.8778", "0"),
    "a033": ("0.0535", "0.0124", "5.93", "0.1051", "1.1757", "0"),
    "a035": ("0.0535", "0.0559", "5.211", "0.17", "0.7258", "0"),
    "a037": ("0.0535", "0.1339", "4.564", "0.2097", "0.5026", "1"),
    "b029": ("0.1294", "0.0033", "6.046", "0.0857", "1.4151", "0"),
    "b031": ("0.1294", "0.0392", "5.224", "0.1775", "0.8761", "0"),
    "b033": ("0.1294", "0.1318", "4.487", "0.2255", "0.6786", "1"),
    "c50": ("0.0111", "0.0479", "4.586", "0.1858", "0.7477", "0"),
    "c65": ("0.2132", "0.0572", "6.153", "0.1298", "1.8528", "0"),
    "d030": ("0.0535", "0.0284", "4.374", "0.2029", "0.7529", "0"),
    "d036": ("0.1294", "0.0511", "6.391", "0.1365", "1.8973", "0"),
    "e0482": ("0.0111", "0.0913", "5.356", "0.187", "1.0057", "0"),
    "e154": ("0.0761", "0.0919", "4.651", "0.1959", "0.6407", "1"),
    "f725": ("0.308", "0.0", "16.994", "0.0", "12.1053", "0"),
}


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
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "cfg.json", "neutral.npz", "ion.npz", "relaxation.npz", "detection.npz",
)


def verify_twin_preregistration() -> None:
    """GR-P1 first half: the frozen twin rows reproduce string-exact from the
    committed scan record. Fails loud — a mismatch means either the frozen
    table or the committed CSV drifted, and the ring must not launch."""
    csv_path = PROJECT_ROOT / TWIN_SCAN_CSV
    with open(csv_path, newline="") as fh:
        idx = {
            (r["v_c"], r["E_bind_tag"], r["tau_ps"], r["E0_eV"]): r
            for r in csv.DictReader(fh)
        }
    for cell in RING_MATRIX:
        # key formatting matches the scan CSV's str(float): "5.5", "7.25"
        key = (str(cell.v_c), cell.eb_tag, str(cell.tau_ps), str(cell.e0_eV))
        row = idx.get(key)
        if row is None:
            raise AssertionError(
                f"[{cell.label}] GR-P1 FAILED: no committed scan row at {key}."
            )
        got = (row["trapped_frac"], row["suppressed_frac"], row["nbar_det"],
               row["n1_solv"], row["w1_solv"], row["gate"])
        if got != TWIN_ROWS[cell.label]:
            raise AssertionError(
                f"[{cell.label}] GR-P1 FAILED: committed scan row {got} != "
                f"frozen pre-registration {TWIN_ROWS[cell.label]}."
            )
    print(f"GR-P1 twin oracle PASSED: all {len(RING_MATRIX)} frozen rows "
          "reproduce string-exact from the committed scan record.",
          flush=True)


def ring_run_dir_name(label: str) -> str:
    """Run-dir basename in the atlas namespace (no ``_tier2_`` substring)."""
    name = run_dir_name(
        CASE, VARIANT, N,
        run_tag=f"tier2atlas_conf{int(round(BUDGET_EV * 100))}_g3r{label}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_cell(cell: RingCell) -> SimConfig:
    """One ring cell's validated config: corrected geometry + this cell's
    (v_c, E_bind, tau, E0), the standing pins otherwise."""
    cfg = build_biphasic_cfg(
        CASE, VARIANT,
        num_molecules=N, ion_time_ps=ION_TIME_PS, dt_ion_ps=DT_ION_PS,
        seed=SEED, lambda0_per_ps=LAMBDA0_PER_PS,
        f_int=cell.e0_eV / BUDGET_EV, f_ret=F_RET, tau_ps=cell.tau_ps,
        evap_rrk_dof=S_EFF, coulomb_available_eV=BUDGET_EV,
        relaxation_time_ps=RELAXATION_TIME_PS,
        relaxation_forces=RELAXATION_FORCES,
        cooling_spatial_gate=COOLING_GATE, drag_form="capped_cubic",
        drag_coefficient_overrides={"v_c": cell.v_c, "p_tail": P_TAIL},
        dissociation_ladder="tabulated",
        tabulated_ladder_rungs_eV=rq4graded_rungs_eV(),
        R0_GS_angstrom=R0_GS_ANGSTROM, E_coulomb_scale=E_COULOMB_SCALE,
        single_initial_position=False,
        detection_time_ps=DETECTION_TIME_PS,
        birth_position_law=BIRTH_LAW,
        initial_position_margin_angstrom=BIRTH_MARGIN_ANGSTROM,
        detection_droplet_retained_policy=RETAINED_POLICY,
        evaporation_shed_convention=SHED_CONVENTION,
        initial_shell_model=SHELL_MODEL,
        internal_energy_partition_law=PARTITION_LAW,
        droplet_size_prior=SIZE_PRIOR,
        use_single_droplet_size=False,
        droplet_size_sampler_mode=SIZE_SAMPLER_MODE,
        binding_energy_molecule_K=PARENT_WELL_K,
    )
    extra: dict[str, object] = dict(
        relaxation_dissipation=RELAXATION_DISSIPATION,
        v_limit_m_per_s=V_LIMIT_M_PER_S,
    )
    if cell.eb_tag != "eb1168":
        # Arm E: the exit well moves off the bundle stamp — the §6.5.1
        # pairing is deliberately broken, so the documented hatch is set
        # (the §6.7 item-2 precedent; provenance stamp stays untouched).
        extra["binding_energy_I_ion_eV"] = EBIND_EV[cell.eb_tag]
        extra["allow_unvalidated_binding_pairing"] = True
    cfg = dataclasses.replace(cfg, **extra)
    cfg.validate()
    return cfg


def verify_corrected_geometry(cfg: SimConfig, cell: RingCell) -> None:
    """The G0/G2 corrected-geometry stamps, checked field by field."""
    if (cfg.droplet_size_prior != "legacy" or cfg.use_single_droplet_size
            or cfg.droplet_size_sampler_mode != "raw"):
        raise AssertionError(
            f"[{cell.label}] corrected size law broken: prior="
            f"{cfg.droplet_size_prior!r}, single={cfg.use_single_droplet_size}, "
            f"mode={cfg.droplet_size_sampler_mode!r}."
        )
    if (cfg.birth_position_law != "boltzmann"
            or cfg.initial_position_margin_angstrom != 0.0
            or cfg.binding_energy_molecule_K != PARENT_WELL_K):
        raise AssertionError(
            f"[{cell.label}] corrected birth law broken: law="
            f"{cfg.birth_position_law!r}, margin="
            f"{cfg.initial_position_margin_angstrom}, well="
            f"{cfg.binding_energy_molecule_K}."
        )
    if cfg.detection_droplet_retained_policy != RETAINED_POLICY:
        raise AssertionError(f"[{cell.label}] retained policy drifted.")
    stamped = cfg.drag_coefficients.effective_binding_energy_I_ion_eV
    if cell.eb_tag == "eb1168":
        if cfg.binding_energy_I_ion_eV != stamped:
            raise AssertionError(
                f"[{cell.label}] E_bind {cfg.binding_energy_I_ion_eV!r} != "
                f"bundle stamp {stamped!r} on a bundle-paired cell."
            )
    else:
        if cfg.binding_energy_I_ion_eV != EBIND_EV[cell.eb_tag]:
            raise AssertionError(f"[{cell.label}] arm-E well not stamped.")
        if not cfg.allow_unvalidated_binding_pairing:
            raise AssertionError(f"[{cell.label}] arm-E hatch not set.")


def expected_diff_keys(cell: RingCell) -> frozenset[str]:
    """The pre-registered cfg-diff key set for this cell vs the reference."""
    keys = set(_BASE_DIFF_KEYS)
    if cell.v_c != STANDING_V_C:
        keys.add(_KNOB_KEYS["v_c"])
    if cell.tau_ps != STANDING_TAU_PS:
        keys.add(_KNOB_KEYS["tau"])
    if cell.e0_eV != STANDING_E0_EV:
        keys.add(_KNOB_KEYS["e0"])
    if cell.eb_tag != "eb1168":
        keys |= _EBIND_KEYS
    return frozenset(keys)


def verify_against_reference(cfg: SimConfig, cell: RingCell) -> list[str]:
    """Diff vs the standing battery reference; only the pre-registered keys
    may differ (GR-P1 second half — an unexpected key = a drifted pin).

    The reference cfg.json predates ``droplet_size_sampler_mode``, and the
    shared diff helper (correctly) refuses to certify a post-reference field
    off its default. The diff therefore runs on a copy pinned to the default;
    the real ``"raw"`` value is asserted by :func:`verify_corrected_geometry`.
    """
    ref_path = PROJECT_ROOT / "data" / "runs" / REFERENCE_RUN / "cfg.json"
    cfg_for_diff = dataclasses.replace(
        cfg, droplet_size_sampler_mode="post_pickup"
    )
    diff = cfg_diff_vs_reference(cfg_for_diff, ref_path, context=cell.label)
    expected = expected_diff_keys(cell)
    if set(diff) != expected:
        raise AssertionError(
            f"[{cell.label}] cfg diff vs {REFERENCE_RUN} is {sorted(diff)}; "
            f"pre-registered exactly {sorted(expected)}."
        )
    return diff


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _run_one(cell: RingCell) -> str:
    cfg = build_cell(cell)
    verify_corrected_geometry(cfg, cell)
    verify_against_reference(cfg, cell)
    run_dir = PROJECT_ROOT / "data" / "runs" / ring_run_dir_name(cell.label)
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{cell.label}] skip (complete) -> {run_dir}", flush=True)
            return cell.label
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(f"{run_dir} exists and is not complete.")
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{cell.label}] v_c={cell.v_c} {cell.eb_tag} tau={cell.tau_ps} "
          f"E0={cell.e0_eV} ({cell.role}) seed={cfg.seed}", flush=True)
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
          f"droplet_retained={int(np.count_nonzero(~dm))}/{dm.size})",
          flush=True)
    return cell.label


def _run_one_by_label(label: str) -> str:
    return _run_one(_SPEC_BY_LABEL[label])


def _select(labels: list[str]) -> list[RingCell]:
    if not labels:
        return list(RING_MATRIX)
    unknown = sorted(set(labels) - set(_SPEC_BY_LABEL))
    if unknown:
        raise SystemExit(f"unknown label(s) {unknown}; expected among "
                         f"{[c.label for c in RING_MATRIX]}")
    return [_SPEC_BY_LABEL[label] for label in labels]


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cells = _select(args.labels)
    print(f"Atlas G3 Step 3 MD ring: {len(cells)} cell(s), N={N}, "
          f"seed={SEED} (fresh, shared), corrected geometry "
          f"({SIZE_PRIOR}+{SIZE_SAMPLER_MODE}, {BIRTH_LAW} {PARENT_WELL_K} K), "
          f"concurrency={args.concurrency}"
          f"{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument runs -- standing point does not move).", flush=True)

    verify_twin_preregistration()

    if args.dry_run:
        for cell in cells:
            cfg = build_cell(cell)
            verify_corrected_geometry(cfg, cell)
            diff = verify_against_reference(cfg, cell)
            print(f"[{cell.label}] cfg OK v_c={cell.v_c} {cell.eb_tag} "
                  f"tau={cell.tau_ps} E0={cell.e0_eV} -> "
                  f"{ring_run_dir_name(cell.label)}", flush=True)
            print(f"[{cell.label}]   diff vs standing: {sorted(diff)}",
                  flush=True)
        print("dry run complete: twin oracle + all cfg guards pass.",
              flush=True)
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
        futures = {poolx.submit(_run_one_by_label, c.label): c.label
                   for c in cells}
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
