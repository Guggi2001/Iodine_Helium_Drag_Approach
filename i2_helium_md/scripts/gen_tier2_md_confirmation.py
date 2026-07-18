"""Generate the Slice-T3 MD-confirmation pilot matrix C1-C4 (plan §I.10).

The Addendum-I Step-2 MD confirmation: the joint (v_c, tau) drag closure found
by the zero-MD Steps 0-1c (calibrated ``capped_cubic`` tails on the rq4graded
ladder, §I.8/§I.9) is confirmed -- or broken -- by real biphasic MD at
**production kinematics** (R0_GS = 2.666 A, E_coulomb_scale = 1.0 -> the
2.70 eV/fragment channel), with the real RRK cascade replacing the 1D frozen
fate map. Four configs (S2-D3: 2 targets + 2 controls), N = 50 pilot each
(S2-D2; the winner re-runs at N = 500 in Slice T4):

======  ======================  =====  ==========  =======  ==========================
config  drag                    tau    ladder      E0 [eV]  role
======  ======================  =====  ==========  =======  ==========================
C1      p = -1, v_c = 7.5       3.8    rq4graded   0.25     full-house target
C2      p = 0,  v_c = 6.5       4.0    rq4graded   0.24     form discrimination
C3      linear_cubic (current)  6.55   flat form_u 0.25     baseline control (S2-P1)
C4      p = -1, v_c = 7.5       4.4    floor1      0.23     bounded-physics claim (I51)
======  ======================  =====  ==========  =======  ==========================

Common pins (§I.10): biphasic; s_eff = 8; ``cooling_spatial_gate =
"density_scaled"``; ``single_initial_position = False`` at the fixed N = 2000
droplet (S2-D4; margin 0 -- the 1D scan used 3-6 A birth margins, recorded
caveat); f_int = E0/2.70 stamped exactly; f_ret = 0.1; lambda_0 = 0.9/ps;
kappa/picture ride the config defaults (physics-live only in C3 -- under a
tabulated ladder they are D_0/Sigma-dead but stay reportable); 30 ps ion
window on the Tier-0 dt; E2 relaxation at the 1000 ps cap; detection stage at
the Sourced t_detect = 8.53 us (CALIBRATION_MAP row 24).

The rq4graded / floor1 rung tables are **constructed here, generator-side**,
from ``d0_of_n`` x the taper multipliers -- no taper physics enters the
package (the Slice-T2 boundary):

* ``rq4graded`` -- the RQ4-graded diagnostic taper **2.2 : 1.5 : 1.3** on
  rungs 1-3 of the flat mixture bottom (findings I51; outside the bounded
  set, reported as a prediction FOR the external RQ4 calculation);
* ``floor1`` -- the knob-free X2/3Pi transition-at-n=1 floor variant:
  rung 1 = **13.3 meV** absolute ([IHe05]; ~1.44x the 9.2 meV flat bottom),
  rest Form-U verbatim (plan L3).

Run dirs live in the probe namespace under the C-label tag
(``tier2probe_conf270_c1`` ...): knobs are read from the authoritative
``cfg.json``, never parsed from the tag (the F3 convention; the two-decimal
f_int encoding would alias C1 with C2 at production). Scoring is Slice T4's
job (ihe_ked run-summary layer + solvated W1/n1/ratio) -- this script runs
the pilots and reports nothing. Pre-registered predictions S2-P1..S2-P4 are
frozen in the plan. **Nothing here discharges F5.**

Usage::

    python scripts/gen_tier2_md_confirmation.py
"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Mapping, NamedTuple, Optional


# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"                     # fixed N = 2000 droplet (S2-D4; §6.6 defense)
VARIANT = "shared_pure_cubic"   # the locked Tier-0 bundle (supplies b)
N = 50                          # S2-D2 pilot scale; winner -> N = 500 (T4)

ION_TIME_PS = 30.0              # probe convention (Tier-0 dt bundle)
DT_ION_PS = 0.01
SEED = 20260604                 # one shared pilot seed (the bridge seed)

# --- production kinematics (§I.10; first in-repo off-9 A MD) ----------------
BUDGET_EV = 2.70                # guard below: this generator is production-only
R0_GS_ANGSTROM = 2.666          # real I2 bond length -> 14.4/2.666/2 = 2.70 eV
E_COULOMB_SCALE = 1.0           # stamped explicitly (droplet preset uses 0.8!)
SINGLE_INITIAL_POSITION = False  # off-center births, margin 0 (recorded caveat)

# --- T9 oracle-chain leg (plan §I.11; leg A'' ACTIVE 2026-07-17) --------------
# "a"         = the delivered Slice-T3 configuration (boltzmann births —
#               byte-identical run-dir names and cfgs; the regression anchor);
# "aprime"    = + birth_position_law="uniform_volume" at the Step-1c
#               full-house margin 3 A (exactly one lever flipped vs leg A;
#               run dirs carry the "ap" config-label prefix — apc1..apc4 —
#               same conf namespace);
# "aprime_cm" = leg A'' — A' + evaporation_shed_convention="co_moving" (the
#               OQ-J working-convention adjudication, 2026-07-17; exactly one
#               lever flipped vs leg A'; run dirs carry the "apcm" prefix —
#               apcmc1..apcmc4). Validates the enum end-to-end and re-baselines
#               the KE axis on the twin's co-moving basis before leg B (T5).
# "b"         = leg B — A'' + initial_shell_model="density_tied" (the Slice-T5
#               dressing arm, 2026-07-18; exactly one lever flipped vs A'';
#               run dirs carry the "b" prefix — bc1..bc4). The dressed A/B
#               against the certified apcm baseline and the twin's legb
#               re-score (pre-registered histograms + per-bin mean KE on the
#               co-moving basis).
LEG = "b"
BIRTH_MARGIN_ANGSTROM = 3.0     # the Step-1c full-house cells' margin (§4k)
# The position axis creates near-barrier transients: marginal E > 0 ions
# (apc3 ion 9: escape margin +1.1 meV) fly *conservative* scattering orbits
# through the droplet in E2 -- the relaxation "coulomb" mode is zero-gamma,
# so no capture is possible and the orbit decouples on its own clock
# (~4000 ps at v_inf ~ 0.32 A/ps for ion 9). The cap must outlast that
# clock. Numerical adequacy, not a physics knob (detection re-reads at
# 8.53 us); the missing-drag-in-E2 question is recorded as an OQ in the log.
# 8000 ps: ion 9 was outbound at d ~ +10 A at 5000 ps and needs d ~ +73 A
# (exposure underflow) at ~0.3 A/ps -- 8000 carries ~10x slack on that leg.
APRIME_RELAXATION_TIME_PS = 8000.0

# --- common mechanism pins ---------------------------------------------------
S_EFF = 8.0                     # the floors-reached kinetics arm (Wave-8 I27)
COOLING_GATE = "density_scaled"  # B.1(1): the intended primary production arm
F_RET = 0.1                     # prior (not identifiable at 9 A only)
LAMBDA0_PER_PS = 0.9            # pickup live (bridge central value)

# --- post-ion stages ---------------------------------------------------------
RELAXATION_TIME_PS = 1000.0     # finite cap (gated arm never freezes; Wave 4)
RELAXATION_FORCES = "coulomb"   # the two I+ fragments still repel

# --- safety ------------------------------------------------------------------
SKIP_COMPLETED_RUNS = True      # resume: skip a run dir with all five artifacts
OVERWRITE_EXISTING_RUN = False  # else refuse to clobber a partial/existing dir


# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.tier2_common import (  # noqa: E402
    EXPERIMENTAL_RELAXATION_TIME_PS,
    PRODUCTION_BUDGET_EV,
    build_biphasic_cfg,
    tier2_confirmation_run_dir_name,
)
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.physics.constants import N_STAR  # noqa: E402
from i2_helium_md.physics.dissociation_ladder import d0_of_n  # noqa: E402
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


# The detector read is the same Sourced 8.53 us flight time the relaxation
# constant records (CALIBRATION_MAP row 24) -- one number, one source.
DETECTION_TIME_PS = EXPERIMENTAL_RELAXATION_TIME_PS

# A conf run dir is "complete" (safe to skip on resume) only with all five.
_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "cfg.json",
    "neutral.npz",
    "ion.npz",
    "relaxation.npz",
    "detection.npz",
)


# =============================================================================
# Ladder tables (generator-side; the Slice-T2 boundary)
# =============================================================================

# Base-table pins: the delivered mixture kappa=1 Form-U ladder (the config
# defaults, stated explicitly so the tables are self-describing).
LADDER_PICTURE_PIN = "statistical_mixture"
LADDER_KAPPA_PIN = 1.0
# Table height: N_STAR + 11 = 32, the Form-U cache height (T2 convention);
# out-of-table occupancy above fails loudly at the lookup.
LADDER_TABLE_HEIGHT = N_STAR + 11

# rq4graded: the RQ4-graded diagnostic taper on rungs 1-3 (findings I51).
RQ4GRADED_TAPER_MULTIPLIERS = (2.2, 1.5, 1.3)
# floor1: the knob-free [IHe05] X2/3Pi first-rung depth [eV] (plan L3).
FLOOR1_RUNG1_EV = 0.0133


def form_u_rungs_eV() -> tuple[float, ...]:
    """The delivered Form-U ladder read through ``d0_of_n(1..32)`` [eV].

    Feeding this table back through ``dissociation_ladder="tabulated"`` is
    bit-identical to Form-U (the Slice-T2 equivalence oracle) -- it is the
    base every variant edits.
    """
    return tuple(
        float(d0_of_n(n, picture=LADDER_PICTURE_PIN, kappa=LADDER_KAPPA_PIN))
        for n in range(1, LADDER_TABLE_HEIGHT + 1)
    )


def rq4graded_rungs_eV() -> tuple[float, ...]:
    """The rq4graded table [eV]: rungs 1-3 x (2.2, 1.5, 1.3), rest Form-U."""
    rungs = list(form_u_rungs_eV())
    for i, mult in enumerate(RQ4GRADED_TAPER_MULTIPLIERS):
        rungs[i] *= mult
    return tuple(rungs)


def floor1_rungs_eV() -> tuple[float, ...]:
    """The floor1 table [eV]: rung 1 = 13.3 meV absolute, rest Form-U."""
    rungs = list(form_u_rungs_eV())
    rungs[0] = FLOOR1_RUNG1_EV
    return tuple(rungs)


def ladder_selector_and_table(
    ladder_key: str,
) -> tuple[Optional[str], Optional[tuple[float, ...]]]:
    """Map a matrix ladder key to ``(dissociation_ladder, rungs)`` kwargs.

    ``"form_u"`` returns ``(None, None)`` -- ride the config default, byte-inert
    (the C3 control). The tabulated keys return their generator-built table.
    """
    if ladder_key == "form_u":
        return None, None
    if ladder_key == "rq4graded":
        return "tabulated", rq4graded_rungs_eV()
    if ladder_key == "floor1":
        return "tabulated", floor1_rungs_eV()
    raise ValueError(
        f"unknown ladder key {ladder_key!r}; expected one of "
        "['form_u', 'rq4graded', 'floor1'] (a new variant is a plan "
        "amendment, not a config value -- §I.10)."
    )


# =============================================================================
# The frozen C1-C4 matrix (plan §I.10 Slice T3)
# =============================================================================


class ConfirmationSpec(NamedTuple):
    """One row of the §I.10 pilot matrix (knob values only -- pins above)."""

    label: str
    drag_form: Optional[str]                                # None = current law
    drag_coefficient_overrides: Optional[Mapping[str, float]]
    tau_ps: float
    ladder_key: str
    e0_eV: float                                            # E_int(0) [eV]
    role: str


CONFIRMATION_MATRIX: tuple[ConfirmationSpec, ...] = (
    ConfirmationSpec(
        "c1", "capped_cubic", {"v_c": 7.5, "p_tail": -1.0},
        3.8, "rq4graded", 0.25, "full-house target",
    ),
    ConfirmationSpec(
        "c2", "capped_cubic", {"v_c": 6.5, "p_tail": 0.0},
        4.0, "rq4graded", 0.24, "form discrimination",
    ),
    ConfirmationSpec(
        "c3", None, None,
        6.55, "form_u", 0.25, "baseline control",
    ),
    ConfirmationSpec(
        "c4", "capped_cubic", {"v_c": 7.5, "p_tail": -1.0},
        4.4, "floor1", 0.23, "bounded-physics claim (I51)",
    ),
)


def _run_is_complete(run_dir: Path) -> bool:
    """True when the run dir holds every conf artifact (resume guard)."""
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def build_confirmation(project_root: Path) -> list[tuple[str, SimConfig, Path]]:
    """Build ``(label, cfg, run_dir)`` for every C-config (no propagation).

    Refuses any budget other than the 2.70 eV production stamp: the Step-2
    confirmation exists only at the production channel (the calibrated
    (v_c, tau) closure was derived for production kinematics; the 0.80 eV
    program is the delivered probe's job).
    """
    if not np.isclose(BUDGET_EV, PRODUCTION_BUDGET_EV):
        raise ValueError(
            f"BUDGET_EV={BUDGET_EV} -- the MD-confirmation matrix runs only at "
            f"the {PRODUCTION_BUDGET_EV} eV production budget (§I.10; the "
            "0.80 eV validation program is the staircase probe's namespace)."
        )
    if LEG == "a":
        birth_law: Optional[str] = None       # ride the boltzmann default
        birth_margin: Optional[float] = None  # byte-inert (the T3 leg exactly)
        retained_policy: Optional[str] = None  # delivered loud guard
        relaxation_ps = RELAXATION_TIME_PS    # the T3 cap exactly
        shed_convention: Optional[str] = None  # ride the cold default
        shell_model: Optional[str] = None      # ride the `full` default
        label_prefix = ""
    elif LEG in ("aprime", "aprime_cm", "b"):
        birth_law = "uniform_volume"
        birth_margin = BIRTH_MARGIN_ANGSTROM
        # The position axis populates the V0-2 droplet-retained class (well-
        # trapped sub-barrier fragments; twin predicts 3-7 %), which can never
        # decouple -- the leg runs under the exclude policy so those ions are
        # classified, not fatal (bookkeeping convention, not a physics lever).
        retained_policy = "exclude"
        relaxation_ps = APRIME_RELAXATION_TIME_PS
        # Leg A'': exactly one lever vs A' -- the OQ-J co-moving shed
        # convention (adjudicated working convention, 2026-07-17). None on A'
        # keeps the delivered cold default byte-identically. Leg B rides A''.
        shed_convention = None if LEG == "aprime" else "co_moving"
        # Leg B: exactly one lever vs A'' -- the Slice-T5 dressing arm. None
        # elsewhere keeps the delivered `full` default byte-identically.
        shell_model = "density_tied" if LEG == "b" else None
        # tag charset is [a-z0-9] -> "apc1" (A') / "apcmc1" (A'') / "bc1" (B)
        label_prefix = {"aprime": "ap", "aprime_cm": "apcm", "b": "b"}[LEG]
    else:
        raise ValueError(
            f"unknown LEG {LEG!r}; expected 'a', 'aprime', 'aprime_cm', or "
            "'b' (a new oracle-chain leg is a plan amendment, not a config "
            "value -- §I.11)."
        )
    scheduled: list[tuple[str, SimConfig, Path]] = []
    for spec in CONFIRMATION_MATRIX:
        selector, rungs = ladder_selector_and_table(spec.ladder_key)
        cfg = build_biphasic_cfg(
            CASE,
            VARIANT,
            num_molecules=N,
            ion_time_ps=ION_TIME_PS,
            dt_ion_ps=DT_ION_PS,
            seed=SEED,
            lambda0_per_ps=LAMBDA0_PER_PS,
            # E_int(0) enters as the exact partition quotient (never rounded;
            # a two-decimal f_int would alias C1 with C2).
            f_int=spec.e0_eV / BUDGET_EV,
            f_ret=F_RET,
            tau_ps=spec.tau_ps,
            evap_rrk_dof=S_EFF,
            coulomb_available_eV=BUDGET_EV,
            relaxation_time_ps=relaxation_ps,
            relaxation_forces=RELAXATION_FORCES,
            cooling_spatial_gate=COOLING_GATE,
            drag_form=spec.drag_form,
            drag_coefficient_overrides=spec.drag_coefficient_overrides,
            dissociation_ladder=selector,
            tabulated_ladder_rungs_eV=rungs,
            R0_GS_angstrom=R0_GS_ANGSTROM,
            E_coulomb_scale=E_COULOMB_SCALE,
            single_initial_position=SINGLE_INITIAL_POSITION,
            detection_time_ps=DETECTION_TIME_PS,
            birth_position_law=birth_law,
            initial_position_margin_angstrom=birth_margin,
            detection_droplet_retained_policy=retained_policy,
            evaporation_shed_convention=shed_convention,
            initial_shell_model=shell_model,
        )
        run_dir = project_root / "data" / "runs" / tier2_confirmation_run_dir_name(
            CASE,
            VARIANT,
            N,
            config_label=f"{label_prefix}{spec.label}",
            budget_eV=BUDGET_EV,
        )
        label = f"{CASE} {VARIANT} N={N} conf[{LEG}] {spec.role}: {spec.label}"
        scheduled.append((label, cfg, run_dir))
    return scheduled


def _run_one(label: str, cfg, run_dir: Path) -> None:
    """Write one conf run dir: neutral -> ion -> relaxation (E2) -> detection."""
    if run_dir.exists():
        # Resume: a fully-written run dir is skipped so a re-run recovers only
        # the points that failed. A *partial* dir (crash mid-run) is not
        # complete, so it still trips the overwrite guard.
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{label}] skip (complete) -> {run_dir}")
            return
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(
                f"{run_dir} already exists and is not a complete run; set "
                "OVERWRITE_EXISTING_RUN=True to regenerate (conf runs are "
                "scored by the generating code version -- regeneration is the "
                "recovery path, not migration)."
            )
    cfg.validate()
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)

    coeffs = cfg.drag_coefficients.coefficients
    drag_desc = f"{cfg.drag_form} " + " ".join(
        f"{k}={float(v):g}" for k, v in sorted(coeffs.items())
    )
    print(
        f"[{label}] budget={cfg.coulomb_available_eV:.2f} eV "
        f"R0={cfg.R0_GS_angstrom:g} A drag=({drag_desc}) "
        f"ladder={cfg.dissociation_ladder} "
        f"tau={cfg.internal_energy_cooling_tau_ps:.2f} ps "
        f"f_int={cfg.internal_energy_partition_fraction:.6f} "
        f"s_eff={cfg.evap_rrk_dof:g} gate={cfg.cooling_spatial_gate}"
    )
    print(f"[{label}] neutral propagation ...")
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)

    print(f"[{label}] ion propagation ...")
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)

    print(f"[{label}] relaxation stage (E2) ...")
    relax = run_relaxation_stage(ion, cfg, save_path=run.root / "relaxation.npz")

    print(f"[{label}] detection stage (t_detect={cfg.detection_time_ps:g} ps) ...")
    detect = run_detection_stage(
        relax.checkpoint, cfg, save_path=run.root / "detection.npz"
    )

    n_relax = float(relax.terminal_n.mean())
    # droplet_retained ions carry their in-droplet handover state verbatim --
    # they must never enter a detected aggregate (review fix 2026-07-18).
    det_mask = detect.detected_mask
    n_retained = int(np.count_nonzero(~det_mask))
    n_detect = (
        float(detect.n_detected[det_mask].mean())
        if det_mask.any() else float("nan")
    )
    print(
        f"[{label}] done -> {run_dir} "
        f"(ion {ion.time_ps.size} steps; relaxed {relax.time_relaxed_ps:.2f} ps, "
        f"n_relax_mean={n_relax:.2f}, n_detect_mean={n_detect:.2f}, "
        f"droplet_retained={n_retained}/{det_mask.size})"
    )


def main() -> int:
    """Generate every C-config of the §I.10 pilot matrix (USER SETTINGS)."""
    scheduled = build_confirmation(PROJECT_ROOT)
    print(
        f"Tier-2 MD-confirmation pilot: {len(scheduled)} config(s), "
        f"budget={BUDGET_EV} eV, N={N} (S2-P1..P4 pre-registered -- scored by "
        "Slice T4, reported not auto-adjudicated)."
    )
    for label, cfg, run_dir in scheduled:
        _run_one(label, cfg, run_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
