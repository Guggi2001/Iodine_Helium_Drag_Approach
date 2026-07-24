"""Atlas §6.6 MD spot-check generator — three capped-lq counterfactual cells.

Runs the pre-registered MD spot-check of the twin outcome (b)
(`TIER2_SENSITIVITY_ATLAS_FINDINGS.md` §6.6; user-selected cases A/B/C,
2026-07-24): the linear+quadratic drag form with its **own** Tier-0 shared
artifacts (bundle ``shared_lq``: a ≈ 0, c = 12.792 amu/Å, jointly-extracted
E_bind = 0.0482 eV — a *consistent* §6.5.1 pairing, no binding escape
hatch) under the production-analog constant-force tail
(``capped_linear_quadratic``, p_tail = −1), at the standing finc1v725
configuration otherwise:

* **A (qca)** — v_c 9.0 / τ 3.5: the twin-best four-way pass (the
  (b)-confirmation cell).
* **B (qcb)** — v_c 9.0 / τ 3.2: the off-needle negative control (twin
  predicts a clear W₁ failure; tests needle-geometry transfer).
* **C (qcc)** — v_c 8.8 / τ 3.4: the second pass corner (window extent).

E₀ = 0.27 eV → f_int = 0.1 for all three (unchanged vs the standing point).

**Config discipline:** every cell is built through
:func:`scripts.tier2_common.build_biphasic_cfg` at the standing pins (leg-D
levers incl. Landau-gated E2 dissipation) and then **verified field-by-field
against the on-disk finc1v725 ``cfg.json``** — the diff must be exactly
{drag_coefficients, binding_energy_I_ion_eV, internal_energy_cooling_tau_ps}
(the battery precedent: the run-dir cfg IS the spec). Any other diff aborts.

**Namespace:** run dirs carry the ``tier2atlas_conf270_<label>`` tag — the
atlas namespace (plan §1.4): no ``_tier2_`` substring (F3 campaign glob) and
no ``tier2probe`` substring (probe glob) may appear; asserted at build time.

Atlas stance: these are counterfactual instrument runs — nothing here moves
the standing production point.

Usage::

    python scripts/gen_tier2atlas_spotcheck.py            # all three cells
    python scripts/gen_tier2atlas_spotcheck.py qca        # one cell
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
import sys
from typing import NamedTuple, Optional


# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"
VARIANT = "shared_lq"           # the lq shared Tier-0 bundle (a, c, E_bind)
N = 500                         # atlas cell convention (one fixed seed)
SEED = 20260721                 # the finc1v725 seed (same neutral draws)
ION_TIME_PS = 30.0
DT_ION_PS = 0.01

BUDGET_EV = 2.70                # production budget (f_int = E0 / budget)
E0_EV = 0.27                    # all three cells (f_int = 0.1, standing)
S_EFF = 8.0
F_RET = 0.1
LAMBDA0_PER_PS = 0.9
COOLING_GATE = "density_scaled"

# leg-D pins (the finc1v725 / T9 chain surface)
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

P_TAIL = -1.0                   # the production-analog constant-force tail

SKIP_COMPLETED_RUNS = True
OVERWRITE_EXISTING_RUN = False

# The standing run whose cfg.json is the field-by-field reference.
FINC_RUN_DIR_NAME = "9A_drag_shared_pure_cubic_N500_tier2probe_conf270_finc1v725"
# Fields allowed to differ from the finc1v725 reference (tau is absent from
# cell B's diff -- its 3.2 ps equals the standing value); the MANDATORY
# subset must always differ (the form swap itself).
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


class SpotCheckSpec(NamedTuple):
    """One §6.6 spot-check cell (knob values only — pins above)."""

    label: str
    v_c_aps: float
    tau_ps: float
    role: str


SPOTCHECK_MATRIX: tuple[SpotCheckSpec, ...] = (
    SpotCheckSpec("qca", 9.0, 3.5, "(b)-confirmation: twin-best pass"),
    SpotCheckSpec("qcb", 9.0, 3.2, "off-needle negative control"),
    SpotCheckSpec("qcc", 8.8, 3.4, "second pass corner (window extent)"),
)


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


def build_cell(spec: SpotCheckSpec) -> SimConfig:
    """Build one cell's cfg (standing pins + the cell's (v_c, tau))."""
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
        tau_ps=spec.tau_ps,
        evap_rrk_dof=S_EFF,
        coulomb_available_eV=BUDGET_EV,
        relaxation_time_ps=RELAXATION_TIME_PS,
        relaxation_forces=RELAXATION_FORCES,
        cooling_spatial_gate=COOLING_GATE,
        drag_form="capped_linear_quadratic",
        # a and c merge from the shared_lq bundle; only the tail is supplied.
        drag_coefficient_overrides={"v_c": spec.v_c_aps, "p_tail": P_TAIL},
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


def verify_against_standing(cfg: SimConfig, spec: SpotCheckSpec) -> None:
    """Field-by-field diff vs the finc1v725 cfg.json: exact allowed set only.

    The standing run's cfg **is** the configuration spec (battery precedent);
    an unexpected diff means a silently-drifted pin and aborts the cell.
    """
    ref_path = (
        PROJECT_ROOT / "data" / "runs" / FINC_RUN_DIR_NAME / "cfg.json"
    )
    ref = json.loads(ref_path.read_text(encoding="utf-8"))
    mine = json.loads(json.dumps(dataclasses.asdict(cfg)))
    if set(ref) != set(mine):
        raise AssertionError(
            f"cfg field sets differ from the standing reference: "
            f"only-in-ref={sorted(set(ref) - set(mine))}, "
            f"only-in-mine={sorted(set(mine) - set(ref))}"
        )
    diff = sorted(k for k in ref if ref[k] != mine[k])
    if not (MANDATORY_CFG_DIFF_KEYS <= set(diff) <= ALLOWED_CFG_DIFF_KEYS):
        raise AssertionError(
            f"[{spec.label}] cfg diff vs finc1v725 is {diff}; the "
            f"pre-registered diff must contain "
            f"{sorted(MANDATORY_CFG_DIFF_KEYS)} and stay within "
            f"{sorted(ALLOWED_CFG_DIFF_KEYS)}"
        )
    # The diff *values* are the pre-registered knobs themselves.
    dc = mine["drag_coefficients"]
    assert dc["form"] == "capped_linear_quadratic"
    assert dc["coefficients"]["v_c"] == spec.v_c_aps
    assert dc["coefficients"]["p_tail"] == P_TAIL
    assert abs(mine["binding_energy_I_ion_eV"] - 0.048236582655347665) < 1e-15
    assert mine["internal_energy_cooling_tau_ps"] == spec.tau_ps


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _run_one(spec: SpotCheckSpec) -> None:
    """Write one spot-check run dir: neutral -> ion -> relaxation -> detection."""
    cfg = build_cell(spec)
    verify_against_standing(cfg, spec)
    run_dir = PROJECT_ROOT / "data" / "runs" / atlas_run_dir_name(spec.label)
    label = f"{spec.label} ({spec.role})"
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{label}] skip (complete) -> {run_dir}", flush=True)
            return
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
        f"[{label}] drag=({drag_desc}) E_bind={cfg.binding_energy_I_ion_eV:.6f} "
        f"tau={cfg.internal_energy_cooling_tau_ps:.2f} ps "
        f"f_int={cfg.internal_energy_partition_fraction:.6f} seed={cfg.seed}",
        flush=True,
    )
    print(f"[{label}] neutral propagation ...", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{label}] ion propagation ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    print(f"[{label}] relaxation stage (E2) ...", flush=True)
    relax = run_relaxation_stage(ion, cfg, save_path=run.root / "relaxation.npz")
    print(f"[{label}] detection stage ...", flush=True)
    detect = run_detection_stage(
        relax.checkpoint, cfg, save_path=run.root / "detection.npz"
    )
    det_mask = detect.detected_mask
    n_detect = (
        float(detect.n_detected[det_mask].mean())
        if det_mask.any() else float("nan")
    )
    print(
        f"[{label}] done -> {run_dir} (n_detect_mean={n_detect:.2f}, "
        f"droplet_retained={int(np.count_nonzero(~det_mask))}/{det_mask.size})",
        flush=True,
    )


def main() -> int:
    wanted = set(sys.argv[1:])
    specs = [
        s for s in SPOTCHECK_MATRIX if not wanted or s.label in wanted
    ]
    if wanted and len(specs) != len(wanted):
        unknown = sorted(wanted - {s.label for s in SPOTCHECK_MATRIX})
        raise SystemExit(f"unknown cell label(s) {unknown}; "
                         f"expected among {[s.label for s in SPOTCHECK_MATRIX]}")
    print(
        f"Atlas §6.6 MD spot-check: {len(specs)} cell(s), N={N}, "
        f"seed={SEED}, budget={BUDGET_EV} eV (counterfactual instrument runs "
        "-- the standing point does not move).",
        flush=True,
    )
    for spec in specs:
        _run_one(spec)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
