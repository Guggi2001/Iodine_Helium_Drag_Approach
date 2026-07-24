"""Atlas §6.7 item 1 — lq sanity battery generator (5 × N = 1000, paired-by-seed).

Escalates the §6.6 MD spot-check (3 × N = 500, single seed;
`TIER2_SENSITIVITY_ATLAS_FINDINGS.md` §6.6) to the N = 1000 battery scale
under the **paired-by-seed** design pre-registered in the findings doc
(§6.7 item 1, frozen 2026-07-24). One cell only — the §6.6 MD-merit winner
**qcc** (`capped_linear_quadratic`, `shared_lq` Tier-0 bundle: a ≈ 9.8e-05,
c = 12.792 amu/Å, jointly-extracted E_bind = 0.0482 eV — a *consistent*
§6.5.1 pairing, no binding escape hatch), v_c 8.8 / p_tail −1, τ 3.4,
E₀ 0.27, standing finc1v725 pins otherwise (leg-D uniform-volume birth,
Landau-gated E2 arm) — run at **five seeds 20260722–20260726**, the cubic
battery's exact seeds, so member ``qccbigs{k}`` pairs by seed (same neutral
draws) with ``bigc1v725s{k}``.

**Config discipline.** Every member is built through
:func:`scripts.tier2_common.build_biphasic_cfg` at the standing pins and
then **verified field-by-field against its paired cubic member's**
``cfg.json`` (matched N + seed) — the diff must be exactly
{drag_form, drag_coefficients, binding_energy_I_ion_eV,
internal_energy_cooling_tau_ps} (the form swap + the lq pair's E_bind +
τ 3.2 → 3.4). Any other diff aborts the member (battery precedent: the
paired run-dir cfg IS the spec).

**Namespace.** Run dirs carry the ``tier2atlas_conf270_<label>`` tag — the
atlas namespace (plan §1.4): no ``_tier2_`` substring (F3 campaign glob)
and no ``tier2probe`` substring (probe glob) may appear; asserted at build.

**Concurrency.** This machine is 4 cores / 17 GB (≈4 GB free at design
time); each N = 1000 member peaks ~1–1.5 GB and ~80 min (8000 ps E2
relaxation dominates). The default is a **2-slot process pool** with BLAS
threads pinned to 1 (set below, before numpy import) — safe on RAM/CPU and
a full walk-away. Bump ``--concurrency`` only after freeing RAM.

Atlas stance: counterfactual instrument runs — nothing here moves the
standing production point.

Usage::

    python scripts/gen_tier2atlas_lqbattery.py                 # all 5, pool of 2
    python scripts/gen_tier2atlas_lqbattery.py --concurrency 3 # all 5, pool of 3
    python scripts/gen_tier2atlas_lqbattery.py --dry-run       # build+verify only
    python scripts/gen_tier2atlas_lqbattery.py qccbigs1        # one member (serial)
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
from pathlib import Path
import sys
from typing import NamedTuple

# --- thread pinning (must precede numpy/BLAS import) -------------------------
# One BLAS thread per process so a pool of K workers does not oversubscribe the
# 4 physical cores (nested threading would thrash). Children inherit this env.
for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")

# UTF-8 stdout/stderr so the banner's τ/§/E₀ survive a cp1252 console or a
# piped/redirected detached log (Windows default is cp1252 when not a TTY).
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass


# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"
VARIANT = "shared_lq"           # the lq shared Tier-0 bundle (a, c, E_bind)
N = 1000                        # battery member size (paired with bigc1v725s*)
ION_TIME_PS = 30.0
DT_ION_PS = 0.01

# --- qcc cell knobs (the §6.6 MD-merit winner) -------------------------------
V_C_APS = 8.8                   # lq cap speed
TAU_PS = 3.4                    # Newton-cooling time
E0_EV = 0.27                    # f_int = E0 / budget = 0.1 (standing)
P_TAIL = -1.0                   # production-analog constant-force tail

BUDGET_EV = 2.70                # production budget
S_EFF = 8.0
F_RET = 0.1
LAMBDA0_PER_PS = 0.9
COOLING_GATE = "density_scaled"

# leg-D pins (the finc1v725 / T9 chain surface) -- identical to the spot-check.
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
RELAXATION_DISSIPATION = "landau_gated_drag"   # finc1v725 Landau-on arm
V_LIMIT_M_PER_S = 58.0                          # 0.58 A/ps

DEFAULT_CONCURRENCY = 2         # 2-slot pool (see module docstring)

SKIP_COMPLETED_RUNS = True
OVERWRITE_EXISTING_RUN = False

# Fields allowed to differ from the paired cubic member; the MANDATORY subset
# must always differ (the form swap itself). τ differs (3.2 → 3.4); seed and N
# match by construction (paired), so neither appears in the diff.
ALLOWED_CFG_DIFF_KEYS = frozenset(
    {
        "drag_form",
        "drag_coefficients",
        "binding_energy_I_ion_eV",
        "internal_energy_cooling_tau_ps",
    }
)
MANDATORY_CFG_DIFF_KEYS = frozenset(
    {"drag_form", "drag_coefficients", "binding_energy_I_ion_eV"}
)

# The lq pair's stamped effective binding (shared_lq bundle) -- asserted exact.
SHARED_LQ_E_BIND_EV = 0.048236582655347665


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

DETECTION_TIME_PS = 8.53e6      # the Sourced flight time (finc pin)

_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "cfg.json", "neutral.npz", "ion.npz", "relaxation.npz", "detection.npz",
)


class BatterySpec(NamedTuple):
    """One §6.7 battery member: a seed and its paired cubic run dir."""

    label: str
    seed: int
    paired_cubic_run_dir: str


BATTERY_MATRIX: tuple[BatterySpec, ...] = (
    BatterySpec("qccbigs1", 20260722,
                "9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s1"),
    BatterySpec("qccbigs2", 20260723,
                "9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s2"),
    BatterySpec("qccbigs3", 20260724,
                "9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s3"),
    BatterySpec("qccbigs4", 20260725,
                "9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s4"),
    BatterySpec("qccbigs5", 20260726,
                "9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s5"),
)

_SPEC_BY_LABEL = {s.label: s for s in BATTERY_MATRIX}


def atlas_run_dir_name(label: str) -> str:
    """The ``tier2atlas`` namespace dir name, with the substring locks."""
    name = run_dir_name(
        CASE, VARIANT, N,
        run_tag=f"tier2atlas_conf{int(round(BUDGET_EV * 100))}_{label}",
    )
    # Namespace locks (plan §1.4): the F3 campaign glob (*_tier2_*) and the
    # probe glob (*tier2probe*) must never match an atlas run.
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_cell(spec: BatterySpec) -> SimConfig:
    """Build one member's cfg (qcc knobs + the member's seed)."""
    cfg = build_biphasic_cfg(
        CASE,
        VARIANT,
        num_molecules=N,
        ion_time_ps=ION_TIME_PS,
        dt_ion_ps=DT_ION_PS,
        seed=spec.seed,
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
        # a and c merge from the shared_lq bundle; only the tail is supplied.
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
    # The Landau-gated E2 arm is a finc1v725 pin outside build_biphasic_cfg's
    # surface (set by the standing run's driver the same way).
    cfg = dataclasses.replace(
        cfg,
        relaxation_dissipation=RELAXATION_DISSIPATION,
        v_limit_m_per_s=V_LIMIT_M_PER_S,
    )
    cfg.validate()
    return cfg


def verify_against_paired(cfg: SimConfig, spec: BatterySpec) -> None:
    """Field-by-field diff vs the paired cubic member cfg.json.

    The paired cubic member (matched N + seed) **is** the configuration spec:
    an lq member must differ from it only by the pre-registered four fields
    (the form swap + the lq pair's E_bind + τ). An unexpected diff means a
    silently-drifted pin and aborts the member.
    """
    ref_path = (
        PROJECT_ROOT / "data" / "runs" / spec.paired_cubic_run_dir / "cfg.json"
    )
    ref = json.loads(ref_path.read_text(encoding="utf-8"))
    mine = json.loads(json.dumps(dataclasses.asdict(cfg)))
    if set(ref) != set(mine):
        raise AssertionError(
            f"[{spec.label}] cfg field sets differ from the paired cubic "
            f"member: only-in-ref={sorted(set(ref) - set(mine))}, "
            f"only-in-mine={sorted(set(mine) - set(ref))}"
        )
    # Paired-by-seed sanity: the reference must actually be this member's seed.
    if ref["seed"] != spec.seed or ref["num_molecules"] != N:
        raise AssertionError(
            f"[{spec.label}] paired cubic member seed/N "
            f"({ref['seed']}/{ref['num_molecules']}) != expected "
            f"({spec.seed}/{N}) -- pairing broken."
        )
    diff = sorted(k for k in ref if ref[k] != mine[k])
    if not (MANDATORY_CFG_DIFF_KEYS <= set(diff) <= ALLOWED_CFG_DIFF_KEYS):
        raise AssertionError(
            f"[{spec.label}] cfg diff vs paired cubic member is {diff}; the "
            f"pre-registered diff must contain "
            f"{sorted(MANDATORY_CFG_DIFF_KEYS)} and stay within "
            f"{sorted(ALLOWED_CFG_DIFF_KEYS)}"
        )
    # The diff *values* are the pre-registered knobs themselves.
    dc = mine["drag_coefficients"]
    assert dc["form"] == "capped_linear_quadratic"
    assert dc["coefficients"]["v_c"] == V_C_APS
    assert dc["coefficients"]["p_tail"] == P_TAIL
    assert abs(mine["binding_energy_I_ion_eV"] - SHARED_LQ_E_BIND_EV) < 1e-15
    assert mine["internal_energy_cooling_tau_ps"] == TAU_PS


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _run_one(spec: BatterySpec) -> str:
    """Write one member run dir: neutral -> ion -> relaxation -> detection."""
    cfg = build_cell(spec)
    verify_against_paired(cfg, spec)
    run_dir = PROJECT_ROOT / "data" / "runs" / atlas_run_dir_name(spec.label)
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
    coeffs = cfg.drag_coefficients.coefficients
    drag_desc = f"{cfg.drag_form} " + " ".join(
        f"{k}={float(v):g}" for k, v in sorted(coeffs.items())
    )
    print(
        f"[{spec.label}] seed={cfg.seed} drag=({drag_desc}) "
        f"E_bind={cfg.binding_energy_I_ion_eV:.6f} "
        f"tau={cfg.internal_energy_cooling_tau_ps:.2f} ps "
        f"f_int={cfg.internal_energy_partition_fraction:.6f}",
        flush=True,
    )
    print(f"[{spec.label}] neutral propagation ...", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{spec.label}] ion propagation ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    print(f"[{spec.label}] relaxation stage (E2) ...", flush=True)
    relax = run_relaxation_stage(ion, cfg, save_path=run.root / "relaxation.npz")
    print(f"[{spec.label}] detection stage ...", flush=True)
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
    """Pool entry point (picklable): resolve the label and run it."""
    return _run_one(_SPEC_BY_LABEL[label])


def _select(labels: list[str]) -> list[BatterySpec]:
    if not labels:
        return list(BATTERY_MATRIX)
    unknown = sorted(set(labels) - set(_SPEC_BY_LABEL))
    if unknown:
        raise SystemExit(
            f"unknown member label(s) {unknown}; expected among "
            f"{[s.label for s in BATTERY_MATRIX]}"
        )
    return [_SPEC_BY_LABEL[label] for label in labels]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*", help="member labels (default: all)")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY,
                        help=f"parallel members (default {DEFAULT_CONCURRENCY})")
    parser.add_argument("--dry-run", action="store_true",
                        help="build + verify every cfg, no propagation")
    args = parser.parse_args(argv)

    specs = _select(args.labels)
    print(
        f"Atlas §6.7 lq sanity battery: {len(specs)} member(s), N={N}, "
        f"cell=qcc (lq v_c {V_C_APS} / τ {TAU_PS} / E₀ {E0_EV}), "
        f"budget={BUDGET_EV} eV, concurrency={args.concurrency}"
        f"{' [DRY RUN]' if args.dry_run else ''} "
        "(counterfactual instrument runs -- the standing point does not move).",
        flush=True,
    )

    if args.dry_run:
        for spec in specs:
            cfg = build_cell(spec)
            verify_against_paired(cfg, spec)
            print(
                f"[{spec.label}] cfg OK (seed={spec.seed}) -> "
                f"{atlas_run_dir_name(spec.label)}; diff vs paired cubic "
                "member = {drag_form, drag_coefficients, binding, τ}",
                flush=True,
            )
        print("dry run complete: all cfgs build + pass the paired-diff guard.",
              flush=True)
        return 0

    concurrency = max(1, args.concurrency)
    if concurrency == 1 or len(specs) == 1:
        for spec in specs:
            _run_one(spec)
        return 0

    # Bounded process pool: at most `concurrency` members propagate at once.
    from concurrent.futures import ProcessPoolExecutor, as_completed

    workers = min(concurrency, len(specs))
    print(f"launching {len(specs)} members through a {workers}-slot pool "
          "(BLAS threads pinned to 1) ...", flush=True)
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_run_one_by_label, s.label): s.label
                   for s in specs}
        for fut in as_completed(futures):
            label = futures[fut]
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001 - report, keep others going
                failures += 1
                print(f"[{label}] FAILED: {exc!r}", flush=True)
    if failures:
        print(f"{failures}/{len(specs)} member(s) failed.", flush=True)
        return 1
    print(f"all {len(specs)} members complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
