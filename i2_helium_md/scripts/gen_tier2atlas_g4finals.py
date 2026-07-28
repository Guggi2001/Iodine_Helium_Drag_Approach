"""Atlas G4 Step 1 Block 3 — the MD finalists at the corrected geometry.

Runs the `TIER2_SENSITIVITY_ATLAS_PLAN` §3.5e Block 3: **6 cells ×
N = 1000, one fresh shared seed**, at the G2-adopted corrected droplet
geometry. MD's role has changed since the G3 ring — from breadth ("does the
basin exist", where N = 500 was right) to **discrimination** among finalists
on W₁/KE differences of 0.1–0.3, which sit inside the N = 500 single-seed
scatter of ±0.13. The budget is therefore spent on precision, not coverage.

Cells (per-cell knobs = (v_c, E_bind, τ, E₀); everything else the standing
finc1v725 pins):

* **f1/f2** — the two best CLEAN ridge cells from the Block-1 fine scan
  (`h2b_g4scan_ridge.csv`): (5.5, 4.8, 0.365) and (5.5, 5.2, 0.34). These
  are the successor-point candidates.
* **f3** — the clean **deep-KE probe**: (5.5, 4.4, 0.395) carries the
  highest twin deepKE (0.887) among the clean leaders. Block 0 measured
  deepKE as rank-untransferable (ρ +0.33), so this axis can only be
  settled in MD; f3 is that test.
* **a037 / b031** — the two G3-ring gated cells, **replicated at N = 1000**.
  Mandatory, not optional: without them the new N = 1000 numbers are not
  comparable to the ring's N = 500 numbers, and they measure the N = 500 →
  1000 sampling-noise contribution to W₁ that the §3.5e open confound
  turns on.
* **x345** — the off-ridge **fail control**, pre-registered to miss:
  (5.5, 4.8, 0.345) has twin n₁ 0.1555 (below the 0.19 band) and predicted
  MD n̄ 4.685 (above the 4.37 + 0.28 ceiling), so it must fail on n₁ low /
  n̄ high. Four E₀ steps from f1 — a deliberately tight control, following
  the c50 lesson that near-gate controls teach the most.

**G4F-P1 oracle (§1.4, enforced before any MD and by --dry-run):** the six
frozen twin rows below must reproduce **string-exact** from the committed
`h2b_g4scan_predictions.csv`, and every cell's cfg must diff against the
standing battery reference in exactly its pre-registered key set.

Retained policy: `exclude_all_coupled` (still the interim convention; the
final call is the G4 adjudication, on the ring's + these cells' measured
class sizes).

Atlas stance: instrument runs — nothing here moves `finc1v725`.

Usage::

    python scripts/gen_tier2atlas_g4finals.py --dry-run   # oracles + guards
    python scripts/gen_tier2atlas_g4finals.py             # all 6 cells
    python scripts/gen_tier2atlas_g4finals.py f1 a037     # selected cells
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
N = 1000                        # fragments per cell (plan §3.5e Block 3)
SEED = 20260729                 # ONE fresh seed shared by every cell (CRN);
#                               # distinct from the ring's 20260728 so the
#                               # replicates are an independent draw
ION_TIME_PS = 30.0
DT_ION_PS = 0.01

BUDGET_EV = 2.70
S_EFF = 8.0
F_RET = 0.1
LAMBDA0_PER_PS = 0.9
COOLING_GATE = "density_scaled"
P_TAIL = -1.0

# The G2-adopted corrected geometry (G0 spec — no new numbers).
SIZE_PRIOR = "legacy"
SIZE_SAMPLER_MODE = "raw"
BIRTH_LAW = "boltzmann"
BIRTH_MARGIN_ANGSTROM = 0.0
PARENT_WELL_K = 313.2           # per-run override; config default untouched

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

EBIND_EV = {"eb0482": 0.0482, "eb154": 0.154}

# Reference run for the cfg-diff guard: a standing-point battery member.
# It is itself N = 1000, so — unlike the G3 ring — `num_molecules` does NOT
# appear in the pre-registered diff key set.
REFERENCE_RUN = "9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s1"

_BASE_DIFF_KEYS = frozenset({
    "seed",
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

STANDING_V_C = 7.25
STANDING_TAU_PS = 3.2
STANDING_E0_EV = 0.27

# The committed fine-scan record the G4F-P1 oracle re-reads.
TWIN_SCAN_CSV = Path("data/runs/h2b_forward_model/h2b_g4scan_predictions.csv")

DEFAULT_CONCURRENCY = 3
SKIP_COMPLETED_RUNS = True
# Complete cells are skipped untouched, so a relaunch only ever rebuilds a
# cell that died mid-stage (the G1/G3 safe-relaunch precedent).
OVERWRITE_EXISTING_RUN = True


class FinalCell(NamedTuple):
    """One §3.5e Block-3 finalist."""

    label: str
    v_c: float
    eb_tag: str
    tau_ps: float
    e0_eV: float
    role: str


FINAL_MATRIX: tuple[FinalCell, ...] = (
    FinalCell("f1", 5.5, "eb1168", 4.8, 0.365, "clean ridge optimum"),
    FinalCell("f2", 5.5, "eb1168", 5.2, 0.34, "clean ridge optimum, tau arm"),
    FinalCell("f3", 5.5, "eb1168", 4.4, 0.395, "clean deep-KE probe"),
    FinalCell("a037", 5.5, "eb1168", 4.8, 0.37,
              "G3-ring bridge control (N-scaling test)"),
    FinalCell("b031", 6.0, "eb1168", 6.4, 0.31,
              "G3-ring bridge control, deepKE-favourable corner"),
    FinalCell("x345", 5.5, "eb1168", 4.8, 0.345,
              "off-ridge fail control (pre-registered to miss on n1)"),
)
_SPEC_BY_LABEL = {c.label: c for c in FINAL_MATRIX}

# The frozen twin rows (committed h2b_g4scan_predictions.csv strings, compared
# EXACTLY): label -> (trapped_frac, suppressed_frac, nbar_det, pred_md_nbar,
# n1_solv, w1_solv, midhot_geo, deepke, gate).
TWIN_ROWS: dict[str, tuple[str, ...]] = {
    "f1": ("0.0535", "0.1118", "4.717", "4.174", "0.2025", "0.5388",
           "1.1103", "0.7757", "1"),
    "f2": ("0.0535", "0.0716", "4.742", "4.194", "0.1981", "0.5378",
           "1.1085", "0.6905", "1"),
    "f3": ("0.0535", "0.1521", "4.751", "4.201", "0.1981", "0.7567",
           "1.1262", "0.8867", "1"),
    "a037": ("0.0535", "0.1339", "4.564", "4.058", "0.2097", "0.5026",
             "1.0703", "0.7552", "1"),
    "b031": ("0.1294", "0.0392", "5.224", "4.561", "0.1775", "0.8761",
             "0.5469", "0.4958", "0"),
    "x345": ("0.0535", "0.0413", "5.387", "4.685", "0.1555", "0.8157",
             "1.2889", "0.8762", "0"),
}

# Pre-registered MD acceptance (unchanged from §3.5d — bias-free).
GATE_N1 = (0.19, 0.30)
GATE_NBAR = (3.77, 4.37)


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

_TWIN_COLUMNS = ("trapped_frac", "suppressed_frac", "nbar_det",
                 "pred_md_nbar", "n1_solv", "w1_solv", "midhot_geo",
                 "deepke", "gate")


def verify_twin_preregistration() -> None:
    """G4F-P1 first half: the frozen twin rows reproduce string-exact from the
    committed fine-scan record. Fails loud — a mismatch means either the frozen
    table or the committed CSV drifted, and the finalists must not launch."""
    csv_path = PROJECT_ROOT / TWIN_SCAN_CSV
    with open(csv_path, newline="") as fh:
        idx = {
            (r["v_c"], r["E_bind_tag"], r["tau_ps"], r["E0_eV"]): r
            for r in csv.DictReader(fh)
        }
    for cell in FINAL_MATRIX:
        key = (str(cell.v_c), cell.eb_tag, str(cell.tau_ps), str(cell.e0_eV))
        row = idx.get(key)
        if row is None:
            raise AssertionError(
                f"[{cell.label}] G4F-P1 FAILED: no committed scan row at {key}."
            )
        got = tuple(row[c] for c in _TWIN_COLUMNS)
        if got != TWIN_ROWS[cell.label]:
            raise AssertionError(
                f"[{cell.label}] G4F-P1 FAILED: committed scan row {got} != "
                f"frozen pre-registration {TWIN_ROWS[cell.label]}."
            )
    print(f"G4F-P1 twin oracle PASSED: all {len(FINAL_MATRIX)} frozen rows "
          "reproduce string-exact from the committed fine-scan record.",
          flush=True)


def finals_run_dir_name(label: str) -> str:
    """Run-dir basename in the atlas namespace (no ``_tier2_`` substring)."""
    name = run_dir_name(
        CASE, VARIANT, N,
        run_tag=f"tier2atlas_conf{int(round(BUDGET_EV * 100))}_g4f{label}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_cell(cell: FinalCell) -> SimConfig:
    """One finalist's validated config: corrected geometry + this cell's
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
        extra["binding_energy_I_ion_eV"] = EBIND_EV[cell.eb_tag]
        extra["allow_unvalidated_binding_pairing"] = True
    cfg = dataclasses.replace(cfg, **extra)
    cfg.validate()
    return cfg


def verify_corrected_geometry(cfg: SimConfig, cell: FinalCell) -> None:
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
            raise AssertionError(f"[{cell.label}] off-bundle well not stamped.")
        if not cfg.allow_unvalidated_binding_pairing:
            raise AssertionError(f"[{cell.label}] pairing hatch not set.")


def expected_diff_keys(cell: FinalCell) -> frozenset[str]:
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


def verify_against_reference(cfg: SimConfig, cell: FinalCell) -> list[str]:
    """Diff vs the standing battery reference; only the pre-registered keys may
    differ (G4F-P1 second half — an unexpected key = a drifted pin).

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


def _run_one(cell: FinalCell) -> str:
    cfg = build_cell(cell)
    verify_corrected_geometry(cfg, cell)
    verify_against_reference(cfg, cell)
    run_dir = PROJECT_ROOT / "data" / "runs" / finals_run_dir_name(cell.label)
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


def _select(labels: list[str]) -> list[FinalCell]:
    if not labels:
        return list(FINAL_MATRIX)
    unknown = sorted(set(labels) - set(_SPEC_BY_LABEL))
    if unknown:
        raise SystemExit(f"unknown label(s) {unknown}; expected among "
                         f"{[c.label for c in FINAL_MATRIX]}")
    return [_SPEC_BY_LABEL[label] for label in labels]


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cells = _select(args.labels)
    print(f"Atlas G4 Block 3 MD finalists: {len(cells)} cell(s), N={N}, "
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
                  f"{finals_run_dir_name(cell.label)}", flush=True)
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
