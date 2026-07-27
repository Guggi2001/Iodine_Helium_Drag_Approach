"""Atlas §6.7 item 2 firm-up — E_bind scan at N=1000, 3 paired seeds.

Seed-replicated escalation of the single-seed §6.7 item-2 scan
(`gen_tier2atlas_ebindscan.py`, N=500, seed 20260721). Runs the two lq
exit-well overrides **{0.1168 = cubic pair, 0.154}** at **N=1000, seeds
20260722/23/24** so each cell pairs 1:1 (matched seed + N) with the item-1
battery members:

* `eb1168@seedK` vs `bigc1v725sK` (cubic @ same well 0.1168) → **form**
  isolated at a matched well;
* `eb1168@seedK` vs `qccbigsK` (lq @ 0.0482) → **well** isolated at a
  matched form.

Cell = qcc (`capped_linear_quadratic`, shared_lq bundle, v_c 8.8 / τ 3.4 /
E₀ 0.27), only the exit well overridden; the drag bundle's stamped
`effective_binding_energy_I_ion_eV` stays 0.0482 (provenance), so the
§6.5.1 pairing is deliberately broken and `allow_unvalidated_binding_pairing`
is set (honest, per run). The cfg diff vs the paired **qccbigs{k}** member
(lq @ 0.0482, N=1000, same seed) is therefore exactly
{binding_energy_I_ion_eV, allow_unvalidated_binding_pairing} — a clean
E_bind OAT.

Atlas stance: instrument runs — nothing moves finc1v725.

Usage::

    python scripts/gen_tier2atlas_ebindseeds.py                 # 6 cells, pool of 3
    python scripts/gen_tier2atlas_ebindseeds.py --dry-run
    python scripts/gen_tier2atlas_ebindseeds.py eb1168_bigs1    # one cell
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
N = 1000                        # matched to the item-1 battery for pairing
ION_TIME_PS = 30.0
DT_ION_PS = 0.01

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

DEFAULT_CONCURRENCY = 3
SKIP_COMPLETED_RUNS = True
OVERWRITE_EXISTING_RUN = False

SHARED_LQ_E_BIND_EV = 0.048236582655347665
CUBIC_PAIR_E_BIND_EV = 0.11675778353879479

# Wells to sweep: label -> binding_energy_I_ion_eV.
WELLS: dict[str, float] = {"eb1168": CUBIC_PAIR_E_BIND_EV, "eb154": 0.154}
# Seed -> paired battery member (lq @ 0.0482, N=1000) used as the diff ref.
SEED_TO_BIGS = {
    20260722: "9A_drag_shared_lq_N1000_tier2atlas_conf270_qccbigs1",
    20260723: "9A_drag_shared_lq_N1000_tier2atlas_conf270_qccbigs2",
    20260724: "9A_drag_shared_lq_N1000_tier2atlas_conf270_qccbigs3",
}

ALLOWED_CFG_DIFF_KEYS = frozenset(
    {"binding_energy_I_ion_eV", "allow_unvalidated_binding_pairing"}
)


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

DETECTION_TIME_PS = 8.53e6
_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "cfg.json", "neutral.npz", "ion.npz", "relaxation.npz", "detection.npz",
)


class EbindCell(NamedTuple):
    label: str          # e.g. "eb1168_bigs1"
    well_label: str     # "eb1168" / "eb154"
    binding_eV: float
    seed: int
    paired_bigs: str    # the qccbigs{k} dir (lq @ 0.0482, same seed)


def _matrix() -> tuple[EbindCell, ...]:
    cells = []
    for k, (seed, bigs) in enumerate(SEED_TO_BIGS.items(), start=1):
        for well_label, binding in WELLS.items():
            cells.append(EbindCell(
                f"{well_label}_bigs{k}", well_label, binding, seed, bigs))
    return tuple(cells)


EBIND_MATRIX = _matrix()
_SPEC_BY_LABEL = {c.label: c for c in EBIND_MATRIX}


def atlas_run_dir_name(label: str) -> str:
    name = run_dir_name(
        CASE, VARIANT, N,
        run_tag=f"tier2atlas_conf{int(round(BUDGET_EV * 100))}_qcc_{label}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_cell(cell: EbindCell) -> SimConfig:
    cfg = build_biphasic_cfg(
        CASE, VARIANT,
        num_molecules=N, ion_time_ps=ION_TIME_PS, dt_ion_ps=DT_ION_PS,
        seed=cell.seed, lambda0_per_ps=LAMBDA0_PER_PS,
        f_int=E0_EV / BUDGET_EV, f_ret=F_RET, tau_ps=TAU_PS,
        evap_rrk_dof=S_EFF, coulomb_available_eV=BUDGET_EV,
        relaxation_time_ps=RELAXATION_TIME_PS, relaxation_forces=RELAXATION_FORCES,
        cooling_spatial_gate=COOLING_GATE, drag_form="capped_linear_quadratic",
        drag_coefficient_overrides={"v_c": V_C_APS, "p_tail": P_TAIL},
        dissociation_ladder="tabulated", tabulated_ladder_rungs_eV=rq4graded_rungs_eV(),
        R0_GS_angstrom=R0_GS_ANGSTROM, E_coulomb_scale=E_COULOMB_SCALE,
        single_initial_position=SINGLE_INITIAL_POSITION,
        detection_time_ps=DETECTION_TIME_PS, birth_position_law=BIRTH_LAW,
        initial_position_margin_angstrom=BIRTH_MARGIN_ANGSTROM,
        detection_droplet_retained_policy=RETAINED_POLICY,
        evaporation_shed_convention=SHED_CONVENTION, initial_shell_model=SHELL_MODEL,
        internal_energy_partition_law=PARTITION_LAW, droplet_size_prior=SIZE_PRIOR,
    )
    cfg = dataclasses.replace(
        cfg,
        relaxation_dissipation=RELAXATION_DISSIPATION,
        v_limit_m_per_s=V_LIMIT_M_PER_S,
        binding_energy_I_ion_eV=cell.binding_eV,
        allow_unvalidated_binding_pairing=True,
    )
    cfg.validate()
    return cfg


def verify_against_bigs(cfg: SimConfig, cell: EbindCell) -> None:
    """Diff vs the paired qccbigs member: only the well + hatch may differ."""
    ref_path = PROJECT_ROOT / "data" / "runs" / cell.paired_bigs / "cfg.json"
    ref = json.loads(ref_path.read_text(encoding="utf-8"))
    mine = json.loads(json.dumps(dataclasses.asdict(cfg)))
    if ref["seed"] != cell.seed or ref["num_molecules"] != N:
        raise AssertionError(
            f"[{cell.label}] paired ref seed/N ({ref['seed']}/"
            f"{ref['num_molecules']}) != expected ({cell.seed}/{N})."
        )
    diff = cfg_diff_vs_reference(cfg, ref_path, context=cell.label)
    if set(diff) != ALLOWED_CFG_DIFF_KEYS:
        raise AssertionError(
            f"[{cell.label}] cfg diff vs qccbigs is {diff}; must be exactly "
            f"{sorted(ALLOWED_CFG_DIFF_KEYS)} (well override + §6.5.1 hatch)."
        )
    assert mine["binding_energy_I_ion_eV"] == cell.binding_eV
    assert mine["allow_unvalidated_binding_pairing"] is True
    dc = mine["drag_coefficients"]
    assert abs(dc["effective_binding_energy_I_ion_eV"] - SHARED_LQ_E_BIND_EV) < 1e-15
    assert dc["coefficients"]["v_c"] == V_C_APS


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _run_one(cell: EbindCell) -> str:
    cfg = build_cell(cell)
    verify_against_bigs(cfg, cell)
    run_dir = PROJECT_ROOT / "data" / "runs" / atlas_run_dir_name(cell.label)
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{cell.label}] skip (complete) -> {run_dir}", flush=True)
            return cell.label
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(f"{run_dir} exists and is not complete.")
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{cell.label}] well={cfg.binding_energy_I_ion_eV:.6f} "
          f"(stamp 0.0482, hatch on) seed={cfg.seed}", flush=True)
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


def _select(labels: list[str]) -> list[EbindCell]:
    if not labels:
        return list(EBIND_MATRIX)
    unknown = sorted(set(labels) - set(_SPEC_BY_LABEL))
    if unknown:
        raise SystemExit(f"unknown label(s) {unknown}; expected among "
                         f"{[c.label for c in EBIND_MATRIX]}")
    return [_SPEC_BY_LABEL[label] for label in labels]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cells = _select(args.labels)
    print(f"Atlas §6.7 item-2 firm-up: {len(cells)} cell(s), N={N}, "
          f"wells {sorted(set(c.well_label for c in cells))}, "
          f"seeds {sorted(set(c.seed for c in cells))}, concurrency="
          f"{args.concurrency}{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument runs -- standing point does not move).", flush=True)

    if args.dry_run:
        for cell in cells:
            cfg = build_cell(cell)
            verify_against_bigs(cfg, cell)
            print(f"[{cell.label}] cfg OK well={cfg.binding_energy_I_ion_eV:.6f} "
                  f"seed={cell.seed} pairs {cell.paired_bigs.split('_')[-1]} -> "
                  f"{atlas_run_dir_name(cell.label)}", flush=True)
        print("dry run complete: all cfgs pass the qccbigs-diff guard.", flush=True)
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
