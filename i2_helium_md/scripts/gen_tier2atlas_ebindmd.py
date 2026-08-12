"""Atlas §6.5 Step 3 — the E_bind MD confirmation cell (capped arm).

Purpose
-------
The §6.5 Step-2 twin scan measured the E_bind **KE transfer function**
``dKE1/dE_bind = -0.5421`` on the capped (h405) arm and ``-0.5486`` on the
linear arm. Only the *linear* arm has an MD cross-check (the committed
free-form ring's CRN well pairs: -0.5268 / -0.5062, mean **-0.5165**), so
the twin's measured ~6 % over-steepness is known on one drag form only.
This script closes that: **one** MD cell, the shallow provenance well
``E_bind = 0.0482`` at the h405 pins, CRN-paired against the *committed*
``linrh405p`` run (same seed 20260731, same N = 500, same everything
else). No second cell is generated — the partner already exists.

Why a single override is enough
-------------------------------
``U(r) = E_bind * (1 - rho_hat(r))``: the droplet exit well and the He
density gate are the SAME erf at the SAME width (14.2 A). Births sit deep
(``U(birth) = 0``), so for a fixed trajectory E_bind is a purely additive,
velocity-independent per-fragment exit toll. Everything the scan measures
is the *deviation* from that ledger, which is why the pair must differ in
the well and nothing else — enforced below by a cfg diff against the
partner run that must come back as exactly two keys.

Provenance posture (the §6.5.1 joint-pairing exception)
-------------------------------------------------------
E_bind is **Derived** — jointly extracted with the drag coefficients by
Tier-0 Method B. This cell overrides the well while **holding the drag
bundle stamp** (0.1168), which deliberately breaks that pairing and so
runs under the documented ``allow_unvalidated_binding_pairing`` hatch —
the §6.7 item-2 precedent ("well overridden, drag stamp held -> honest
under the hatch"). The result is a **sensitivity read, never a candidate
point**; §9.1 bounds the physically defensible variation at <= 0.009 eV.

Cost: 1 x N = 500, ~1 MD cell. Nothing here moves any standing point.

Usage::

    python scripts/gen_tier2atlas_ebindmd.py --dry-run   # guards only
    python scripts/gen_tier2atlas_ebindmd.py
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path
from typing import Optional

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2atlas_linring import (  # noqa: E402
    EBIND_EV,
    RING_MATRIX,
    SEED,
    N,
    build_cell,
)
from scripts.tier0_common import run_dir_name  # noqa: E402
from scripts.tier2_common import cfg_diff_vs_reference  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================

LABEL = "ebmdh405s"          # h405 pins, shallow well
WELL_TAG = "eb0482"
BUDGET_EV = 2.70             # namespace tag only (matches the ring)
CASE, VARIANT = "9A", "shared_pure_cubic"   # run_dir_name prepends "drag_"

# The committed CRN partner: the free-form ring's h405 incumbent cell.
PARTNER_RUN = "9A_drag_shared_pure_cubic_N500_tier2atlas_conf270_linrh405p"

# The ONLY two keys this cell may differ from its partner in.
EXPECTED_DIFF_KEYS = frozenset({
    "binding_energy_I_ion_eV",
    "allow_unvalidated_binding_pairing",
})

_REQUIRED_ARTIFACTS = (
    "cfg.json", "neutral.npz", "ion.npz", "relaxation.npz", "detection.npz",
)

# --- Pre-registered predictions (FROZEN 2026-08-10, before the run) ---------
# Read by the scorer; recorded here so the freeze lives with the generator.
MDP1_SLOPE_BAND = (-0.60, -0.44)   # the refund exists in MD
MDP2_SLOPE_BAND = (-0.55, -0.47)   # twin over-steepness transfers (~0.510)
MDP3_KE2_MAX_GAP = 0.10            # rigid translation across bins
MDP4_TRAP_MAX = 0.025              # trap(eb0482) given lever 0.6-1.1/eV
MDP5_GRADING_BAND = (0.05, 0.25)   # MD (mid - n1) exceeds the twin's +0.020
MDP6_KE1_CEILING = 0.75            # the (A)-ceiling band still not reached
# The Step-2 twin forecast this cell tests (committed atlas_ebind_twin.csv).
TWIN_SLOPE_CAPPED = -0.5421
TWIN_KE1_EB0482_CAPPED = 0.6608


def h405_spec():
    """The ring's h405 partner spec — the pins this cell inherits verbatim."""
    for cell in RING_MATRIX:
        if cell.label == "h405p":
            return cell
    raise AssertionError("no h405p cell in the committed RING_MATRIX")


def build_shallow_well_cell() -> SimConfig:
    """The h405 config with the well moved to 0.0482 eV and nothing else.

    Built by calling the ring's own ``build_cell`` on the h405p spec, so
    every pin is inherited rather than re-typed; the two overrides are
    applied afterwards and are exactly the keys ``EXPECTED_DIFF_KEYS``
    names. The drag bundle stamp is deliberately NOT touched — holding it
    is what makes this an honest, hatch-declared pairing exception.
    """
    cfg = build_cell(h405_spec())
    cfg = dataclasses.replace(
        cfg,
        binding_energy_I_ion_eV=EBIND_EV[WELL_TAG],
        allow_unvalidated_binding_pairing=True,
    )
    cfg.validate()
    return cfg


def verify_partner_pairing(cfg: SimConfig) -> list[str]:
    """The CRN guard: diff vs the committed partner must be exactly the well.

    This is the whole experiment's validity condition — a third key in the
    diff means the pair is not a common-random-numbers pair and the
    measured slope would carry an uncontrolled second difference.

    It also **subsumes** the ring's ``verify_corrected_geometry`` stamp
    check: the partner is itself a geometry-verified committed cell, so a
    two-key diff proves this cell inherits that geometry exactly. (The ring
    helper cannot be reused directly — on the capped branch it asserts
    ``well == bundle stamp``, which is precisely the pairing this cell
    breaks on purpose.)
    """
    ref_path = PROJECT_ROOT / "data" / "runs" / PARTNER_RUN / "cfg.json"
    if not ref_path.exists():
        raise AssertionError(
            f"CRN partner {PARTNER_RUN} is missing — the committed h405 ring "
            "cell must exist before this cell may run")
    diff = cfg_diff_vs_reference(cfg, ref_path, context=LABEL)
    if set(diff) != EXPECTED_DIFF_KEYS:
        raise AssertionError(
            f"[{LABEL}] cfg diff vs the CRN partner is {sorted(diff)}; "
            f"pre-registered exactly {sorted(EXPECTED_DIFF_KEYS)}."
        )
    if cfg.binding_energy_I_ion_eV != EBIND_EV[WELL_TAG]:
        raise AssertionError(f"[{LABEL}] well not stamped.")
    stamped = cfg.drag_coefficients.effective_binding_energy_I_ion_eV
    if stamped != EBIND_EV["eb1168"]:
        raise AssertionError(
            f"[{LABEL}] the drag bundle stamp moved ({stamped!r}); this cell "
            "must hold it so the pairing exception stays the declared one.")
    if not cfg.allow_unvalidated_binding_pairing:
        raise AssertionError(f"[{LABEL}] the §6.5.1 hatch is not set.")
    print(f"[{LABEL}] CRN guard PASSED: diff vs {PARTNER_RUN} is exactly "
          f"{sorted(diff)}", flush=True)
    return sorted(diff)


def md_run_dir_name() -> str:
    name = run_dir_name(
        CASE, VARIANT, N,
        run_tag=f"tier2atlas_conf{int(round(BUDGET_EV * 100))}_{LABEL}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def run_cell(cfg: SimConfig) -> Path:
    run_dir = PROJECT_ROOT / "data" / "runs" / md_run_dir_name()
    if run_dir.exists() and all(
            (run_dir / n).exists() for n in _REQUIRED_ARTIFACTS):
        print(f"[{LABEL}] skip (complete) -> {run_dir}", flush=True)
        return run_dir
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{LABEL}] capped_cubic v_c={h405_spec().v_c} {WELL_TAG} "
          f"(E_bind={cfg.binding_energy_I_ion_eV}) tau={h405_spec().tau_ps} "
          f"E0={h405_spec().e0_eV} seed={cfg.seed} N={N}", flush=True)
    print(f"[{LABEL}] neutral ...", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{LABEL}] ion ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    print(f"[{LABEL}] relaxation (E2) ...", flush=True)
    relax = run_relaxation_stage(ion, cfg,
                                 save_path=run.root / "relaxation.npz")
    print(f"[{LABEL}] detection ...", flush=True)
    detect = run_detection_stage(relax.checkpoint, cfg,
                                 save_path=run.root / "detection.npz")
    dm = detect.detected_mask
    nd = float(detect.n_detected[dm].mean()) if dm.any() else float("nan")
    print(f"[{LABEL}] done -> {run_dir} (n_detect_mean={nd:.2f}, "
          f"droplet_retained={int(np.count_nonzero(~dm))}/{dm.size})",
          flush=True)
    return run_dir


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    print(f"Atlas §6.5 Step 3: 1 cell, N={N}, seed={SEED}, well "
          f"{WELL_TAG}={EBIND_EV[WELL_TAG]} eV, CRN partner {PARTNER_RUN}"
          f"{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument run -- standing point does not move).", flush=True)

    cfg = build_shallow_well_cell()
    verify_partner_pairing(cfg)
    print(f"[{LABEL}] twin forecast under test: slope "
          f"{TWIN_SLOPE_CAPPED:+.4f} eV/eV, KE1({WELL_TAG}) "
          f"{TWIN_KE1_EB0482_CAPPED:.4f}", flush=True)

    if args.dry_run:
        print(f"dry run complete: all guards pass -> {md_run_dir_name()}",
              flush=True)
        return 0
    run_cell(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
