"""Shared Tier-2 run wiring (Phase D bridge; Phase F extends this module).

Tier 2 deliberately reuses the locked Tier-0 drag bundle plumbing, exactly as
Tier 1a did: :func:`build_biphasic_cfg` builds a Tier-0 drag config and swaps
only the mass scenario plus the generative-knob values. Seeded at Phase D
(Slice Z) with the single pinned representative priored point; the Phase-F
campaign (F1) extends this module with the run-matrix conventions instead of
re-implementing it.

Pinned representative priored point (plan §0/§4, 2026-07-02 — chosen, not
fitted; Phase F owns the real values):

* λ₀ = 0.9 /ps (central of the 0.7–1.1 band),
* f_int = 0.5 (inside the derived 0.80 eV window [0.21–0.24, ~0.6]; lands the
  closed-form t× in the GAH25 5–6.5 ps prior),
* f_ret = 0.1 (prior small, nonzero to exercise S1),
* τ = the config default 6.55 ps (the geometric mid √(2.6·16.5) delivered at
  Slice K; the plan's "6.5" is rounded shorthand for the same mid),
* κ = 1.0 and picture = ``statistical_mixture`` (the config defaults),
* ν = 2.42 /ps and s = 3n−3 (Sourced/Derived — no choice, config defaults).
"""

from __future__ import annotations

from dataclasses import replace
from typing import Mapping, Optional

from i2_helium_md.config import SimConfig
from i2_helium_md.physics.shell_schedule import complex_mass_amu

from scripts.tier0_common import build_drag_cfg, run_dir_name


TIER2_BRIDGE_TAG = "tier2_bridge_biphasic"

# The pinned representative priored point (demonstration choice, not a fit).
BRIDGE_LAMBDA0_PER_PS = 0.9
BRIDGE_F_INT = 0.5
BRIDGE_F_RET = 0.1


def build_biphasic_cfg(
    case: str,
    variant: str,
    *,
    num_molecules: int,
    ion_time_ps: float,
    dt_ion_ps: float,
    seed: int,
    lambda0_per_ps: float = BRIDGE_LAMBDA0_PER_PS,
    f_int: float = BRIDGE_F_INT,
    f_ret: float = BRIDGE_F_RET,
    coeff_overrides: Optional[Mapping[str, float]] = None,
    e_bind_override: Optional[float] = None,
) -> SimConfig:
    """Build a Tier-2 ``biphasic`` config from a Tier-0 drag config.

    The drag coefficients, effective binding, geometry, timestep, and seed all
    come from :func:`scripts.tier0_common.build_drag_cfg` (the run parameters
    inherit Tier-1a exactly; plan §0). The Tier-2-specific changes are the
    mass-scenario switch, the generative-knob values (pinned priored point by
    default; Phase F sweeps them by argument), the scenario-stamped 0.80 eV
    budget, and the §6.6 escape hatch for the intentional constant-``m_eff``
    coefficient pairing (the same structural §6.5 trip as Tier 1a).

    τ, κ, picture, ν, s, p, and the cap/rate-form selectors stay at their
    config defaults (module docstring). ``mass_initial_amu`` is set to the
    n = 21 complex mass for metadata consistency; the run-time source of truth
    remains ``simulation/ion_initial_state.py`` (the biphasic branch).
    """
    fixed_cfg = build_drag_cfg(
        case,
        variant,
        num_molecules=num_molecules,
        ion_time_ps=ion_time_ps,
        dt_ion_ps=dt_ion_ps,
        seed=seed,
        coeff_overrides=coeff_overrides,
        e_bind_override=e_bind_override,
    )
    cfg = replace(
        fixed_cfg,
        mass_scenario="biphasic",
        coulomb_available_eV=0.80,
        pickup_rate_coefficient=float(lambda0_per_ps),
        internal_energy_partition_fraction=float(f_int),
        internal_energy_retained_fraction=float(f_ret),
        allow_inconsistent_mass_pairing=True,
        mass_initial_amu=complex_mass_amu(21),
    )
    cfg.validate()
    return cfg


def tier2_bridge_run_dir_name(case: str, variant: str, n: int) -> str:
    """Return the Tier-0-style run directory basename for the bridge run."""
    return run_dir_name(case, variant, n, run_tag=TIER2_BRIDGE_TAG)


__all__ = [
    "BRIDGE_F_INT",
    "BRIDGE_F_RET",
    "BRIDGE_LAMBDA0_PER_PS",
    "TIER2_BRIDGE_TAG",
    "build_biphasic_cfg",
    "tier2_bridge_run_dir_name",
]
