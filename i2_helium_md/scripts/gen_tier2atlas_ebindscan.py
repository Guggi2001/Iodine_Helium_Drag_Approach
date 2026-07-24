"""Atlas §6.7 item 2 — E_bind pair-separation scan (qcc, 2 × N = 500).

The confound separator for the §6.7 item-1 battery: the lq cells' KE/supp
differences vs cubic are entangled with lq's co-extracted **shallower exit
well** (E_bind 0.0482 eV). This scan holds the lq drag law fixed and sweeps
**only the exit well** to the cubic pair's value and above, so the residual
form-vs-well attribution becomes a measurement.

Cell = qcc (`capped_linear_quadratic`, shared_lq bundle, v_c 8.8 / τ 3.4 /
E₀ 0.27), seed 20260721, N = 500 — identical to the on-disk
`…tier2atlas_conf270_qcc` cell **except** the exit well:

* **eb1168** — `binding_energy_I_ion_eV = 0.11675778353879479`
  (the *exact* standing cubic-pair value, so this cell is byte-matched on
  E_bind to `finc1v725` → a clean **form-at-matched-well** comparison:
  lq@0.1168 vs cubic finc@0.1168, same seed 20260721);
* **eb154** — `binding_energy_I_ion_eV = 0.154` (top of the §6.5 Tier-0
  extracted spread).

**These deliberately break the §6.5.1 joint pairing** — the drag bundle's
stamped `effective_binding_energy_I_ion_eV` stays 0.0482 (provenance: what
the lq coefficients were co-extracted with), while the run's well is
overridden. That is honest and requires `allow_unvalidated_binding_pairing`
(the §6.5 pre-registered exception, stamped per run). The cfg diff vs the
on-disk qcc cell is therefore **exactly {binding_energy_I_ion_eV,
allow_unvalidated_binding_pairing}** — a clean E_bind OAT.

Focus reads (§6.5): trapped, suppressed, deep-bin KE, χ²_med. Headline:
does the item-1 lq KE/supp difference survive when the well is set to the
cubic pair's 0.1168 (→ it was form) or collapse (→ it was the well)?

Atlas stance: counterfactual instrument runs — nothing moves finc1v725.

Usage::

    python scripts/gen_tier2atlas_ebindscan.py            # both cells, pool of 2
    python scripts/gen_tier2atlas_ebindscan.py --dry-run  # build + verify only
    python scripts/gen_tier2atlas_ebindscan.py eb1168     # one cell
"""

from __future__ import annotations

import argparse
import dataclasses
import json
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
VARIANT = "shared_lq"
N = 500                         # atlas cell convention (one fixed seed)
SEED = 20260721                 # the qcc / finc1v725 seed (same neutral draws)
ION_TIME_PS = 30.0
DT_ION_PS = 0.01

# qcc cell knobs (unchanged from the §6.6/§6.7 qcc cell)
V_C_APS = 8.8
TAU_PS = 3.4
E0_EV = 0.27
P_TAIL = -1.0

BUDGET_EV = 2.70
S_EFF = 8.0
F_RET = 0.1
LAMBDA0_PER_PS = 0.9
COOLING_GATE = "density_scaled"

BIRTH_LAW = "uniform_volume"
BIRTH_MARGIN_ANGSTROM = 3.0
RETAINED_POLICY = "exclude"
SHED_CONVENTION = "co_moving"
SHELL_MODEL = "density_tied"
PARTITION_LAW = "sigma_proportional"
SIZE_PRIOR = "kornilov_lognormal"
R0_GS_ANGSTROM = 2.666
E_COULOMB_SCALE = 1.0
SINGLE_INITIAL_POSITION = False
RELAXATION_TIME_PS = 8000.0
RELAXATION_FORCES = "coulomb"
RELAXATION_DISSIPATION = "landau_gated_drag"
V_LIMIT_M_PER_S = 58.0

DEFAULT_CONCURRENCY = 2
SKIP_COMPLETED_RUNS = True
OVERWRITE_EXISTING_RUN = False

# The on-disk qcc cell whose cfg.json is the field-by-field reference.
QCC_RUN_DIR_NAME = "9A_drag_shared_lq_N500_tier2atlas_conf270_qcc"
# The lq bundle's stamped effective binding (unchanged by this scan).
SHARED_LQ_E_BIND_EV = 0.048236582655347665
# The exact standing cubic-pair well (byte-matched to finc1v725).
CUBIC_PAIR_E_BIND_EV = 0.11675778353879479

# Diff vs the qcc cell: exactly the well override + its escape hatch.
ALLOWED_CFG_DIFF_KEYS = frozenset(
    {"binding_energy_I_ion_eV", "allow_unvalidated_binding_pairing"}
)
MANDATORY_CFG_DIFF_KEYS = ALLOWED_CFG_DIFF_KEYS


# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2_md_confirmation import rq4graded_rungs_eV  # noqa: E402
from scripts.tier0_common import run_dir_name  # noqa: E402
from scripts.tier2_common import build_biphasic_cfg  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

DETECTION_TIME_PS = 8.53e6

_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "cfg.json", "neutral.npz", "ion.npz", "relaxation.npz", "detection.npz",
)


class EbindSpec(NamedTuple):
    """One E_bind scan cell: a label and the overridden exit well [eV]."""

    label: str
    binding_eV: float
    role: str


EBIND_MATRIX: tuple[EbindSpec, ...] = (
    EbindSpec("eb1168", CUBIC_PAIR_E_BIND_EV,
              "well = cubic pair 0.1168 (form-at-matched-well vs finc)"),
    EbindSpec("eb154", 0.154, "well = 0.154 (top of Tier-0 spread)"),
)

_SPEC_BY_LABEL = {s.label: s for s in EBIND_MATRIX}


def atlas_run_dir_name(label: str) -> str:
    name = run_dir_name(
        CASE, VARIANT, N,
        run_tag=f"tier2atlas_conf{int(round(BUDGET_EV * 100))}_qcc_{label}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_cell(spec: EbindSpec) -> SimConfig:
    """Build one E_bind cell: the qcc cfg with only the exit well overridden."""
    cfg = build_biphasic_cfg(
        CASE,
        VARIANT,
        num_molecules=N,
        ion_time_ps=ION_TIME_PS,
        dt_ion_ps=DT_ION_PS,
        seed=SEED,
        lambda0_per_ps=LAMBDA0_PER_PS,
        f_int=E0_EV / BUDGET_EV,
        f_ret=F_RET,
        tau_ps=TAU_PS,
        evap_rrk_dof=S_EFF,
        coulomb_available_eV=BUDGET_EV,
        relaxation_time_ps=RELAXATION_TIME_PS,
        relaxation_forces=RELAXATION_FORCES,
        cooling_spatial_gate=COOLING_GATE,
        drag_form="capped_linear_quadratic",
        drag_coefficient_overrides={"v_c": V_C_APS, "p_tail": P_TAIL},
        dissociation_ladder="tabulated",
        tabulated_ladder_rungs_eV=rq4graded_rungs_eV(),
        R0_GS_angstrom=R0_GS_ANGSTROM,
        E_coulomb_scale=E_COULOMB_SCALE,
        single_initial_position=SINGLE_INITIAL_POSITION,
        detection_time_ps=DETECTION_TIME_PS,
        birth_position_law=BIRTH_LAW,
        initial_position_margin_angstrom=BIRTH_MARGIN_ANGSTROM,
        detection_droplet_retained_policy=RETAINED_POLICY,
        evaporation_shed_convention=SHED_CONVENTION,
        initial_shell_model=SHELL_MODEL,
        internal_energy_partition_law=PARTITION_LAW,
        droplet_size_prior=SIZE_PRIOR,
    )
    # Landau-gated E2 arm (finc pin, outside build_biphasic_cfg's surface) AND
    # the well override: the drag bundle's stamped effective binding stays
    # 0.0482 (provenance), only the run's climbed well changes → the §6.5.1
    # pairing is deliberately broken, so the escape hatch is set (honest).
    cfg = dataclasses.replace(
        cfg,
        relaxation_dissipation=RELAXATION_DISSIPATION,
        v_limit_m_per_s=V_LIMIT_M_PER_S,
        binding_energy_I_ion_eV=spec.binding_eV,
        allow_unvalidated_binding_pairing=True,
    )
    cfg.validate()
    return cfg


def verify_against_qcc(cfg: SimConfig, spec: EbindSpec) -> None:
    """Field-by-field diff vs the on-disk qcc cfg.json: only the well + hatch.

    The qcc cell (E_bind 0.0482, no hatch) IS the spec; an E_bind cell must
    differ from it by exactly {binding_energy_I_ion_eV,
    allow_unvalidated_binding_pairing}. Anything else = a drifted pin, aborts.
    """
    ref_path = PROJECT_ROOT / "data" / "runs" / QCC_RUN_DIR_NAME / "cfg.json"
    ref = json.loads(ref_path.read_text(encoding="utf-8"))
    mine = json.loads(json.dumps(dataclasses.asdict(cfg)))
    if set(ref) != set(mine):
        raise AssertionError(
            f"[{spec.label}] cfg field sets differ from the qcc reference: "
            f"only-in-ref={sorted(set(ref) - set(mine))}, "
            f"only-in-mine={sorted(set(mine) - set(ref))}"
        )
    if ref["seed"] != SEED or ref["num_molecules"] != N:
        raise AssertionError(
            f"[{spec.label}] qcc reference seed/N ({ref['seed']}/"
            f"{ref['num_molecules']}) != expected ({SEED}/{N})."
        )
    diff = sorted(k for k in ref if ref[k] != mine[k])
    if set(diff) != ALLOWED_CFG_DIFF_KEYS:
        raise AssertionError(
            f"[{spec.label}] cfg diff vs qcc is {diff}; the pre-registered "
            f"diff must be exactly {sorted(ALLOWED_CFG_DIFF_KEYS)} (the E_bind "
            "override + its §6.5.1 escape hatch)."
        )
    assert mine["binding_energy_I_ion_eV"] == spec.binding_eV
    assert mine["allow_unvalidated_binding_pairing"] is True
    # provenance stamp unchanged (the drag law is NOT re-paired):
    dc = mine["drag_coefficients"]
    assert abs(dc["effective_binding_energy_I_ion_eV"] - SHARED_LQ_E_BIND_EV) < 1e-15
    assert dc["form"] == "capped_linear_quadratic"
    assert dc["coefficients"]["v_c"] == V_C_APS


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _run_one(spec: EbindSpec) -> str:
    cfg = build_cell(spec)
    verify_against_qcc(cfg, spec)
    run_dir = PROJECT_ROOT / "data" / "runs" / atlas_run_dir_name(spec.label)
    label = f"{spec.label} ({spec.role})"
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{spec.label}] skip (complete) -> {run_dir}", flush=True)
            return spec.label
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(
                f"{run_dir} already exists and is not a complete run; set "
                "OVERWRITE_EXISTING_RUN=True to regenerate."
            )
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(
        f"[{label}] E_bind={cfg.binding_energy_I_ion_eV:.6f} eV "
        f"(stamp {cfg.drag_coefficients.effective_binding_energy_I_ion_eV:.4f}, "
        f"hatch={cfg.allow_unvalidated_binding_pairing}) seed={cfg.seed}",
        flush=True,
    )
    print(f"[{spec.label}] neutral ...", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{spec.label}] ion ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    print(f"[{spec.label}] relaxation (E2) ...", flush=True)
    relax = run_relaxation_stage(ion, cfg, save_path=run.root / "relaxation.npz")
    print(f"[{spec.label}] detection ...", flush=True)
    detect = run_detection_stage(
        relax.checkpoint, cfg, save_path=run.root / "detection.npz"
    )
    det_mask = detect.detected_mask
    n_detect = (
        float(detect.n_detected[det_mask].mean())
        if det_mask.any() else float("nan")
    )
    print(
        f"[{spec.label}] done -> {run_dir} (n_detect_mean={n_detect:.2f}, "
        f"droplet_retained={int(np.count_nonzero(~det_mask))}/{det_mask.size})",
        flush=True,
    )
    return spec.label


def _run_one_by_label(label: str) -> str:
    return _run_one(_SPEC_BY_LABEL[label])


def _select(labels: list[str]) -> list[EbindSpec]:
    if not labels:
        return list(EBIND_MATRIX)
    unknown = sorted(set(labels) - set(_SPEC_BY_LABEL))
    if unknown:
        raise SystemExit(
            f"unknown cell label(s) {unknown}; expected among "
            f"{[s.label for s in EBIND_MATRIX]}"
        )
    return [_SPEC_BY_LABEL[label] for label in labels]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*", help="cell labels (default: both)")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true",
                        help="build + verify every cfg, no propagation")
    args = parser.parse_args(argv)

    specs = _select(args.labels)
    print(
        f"Atlas §6.7 item 2 E_bind scan: {len(specs)} cell(s), N={N}, "
        f"seed={SEED}, qcc lq v_c {V_C_APS} / τ {TAU_PS} / E₀ {E0_EV}, "
        f"budget={BUDGET_EV} eV, concurrency={args.concurrency}"
        f"{' [DRY RUN]' if args.dry_run else ''} "
        "(instrument runs -- the standing point does not move).",
        flush=True,
    )

    if args.dry_run:
        for spec in specs:
            cfg = build_cell(spec)
            verify_against_qcc(cfg, spec)
            print(
                f"[{spec.label}] cfg OK: E_bind={cfg.binding_energy_I_ion_eV:.6f}"
                f" (stamp 0.0482, hatch on) -> {atlas_run_dir_name(spec.label)};"
                " diff vs qcc = {binding_energy_I_ion_eV, "
                "allow_unvalidated_binding_pairing}",
                flush=True,
            )
        print("dry run complete: both cfgs build + pass the qcc-diff guard.",
              flush=True)
        return 0

    concurrency = max(1, args.concurrency)
    if concurrency == 1 or len(specs) == 1:
        for spec in specs:
            _run_one(spec)
        return 0

    from concurrent.futures import ProcessPoolExecutor, as_completed

    workers = min(concurrency, len(specs))
    print(f"launching {len(specs)} cells through a {workers}-slot pool "
          "(BLAS threads pinned to 1) ...", flush=True)
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as poolx:
        futures = {poolx.submit(_run_one_by_label, s.label): s.label
                   for s in specs}
        for fut in as_completed(futures):
            label = futures[fut]
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[{label}] FAILED: {exc!r}", flush=True)
    if failures:
        print(f"{failures}/{len(specs)} cell(s) failed.", flush=True)
        return 1
    print(f"all {len(specs)} cells complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
