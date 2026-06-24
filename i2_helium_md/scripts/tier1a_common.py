"""Shared Tier-1a run wiring for the anchored mass-dynamics RMSE table.

Tier 1a deliberately reuses the locked Tier-0 drag bundle plumbing. This module
only adds the anchored-discrete config mutation and the run-tag convention that
the Tier-1a generator and scorer share. Anchored tags include ``continuous`` so
regenerated physical Tier-1a runs are not confused with stale cold-shed runs.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Mapping, Optional

from i2_helium_md.config import SimConfig
from i2_helium_md.physics.shell_schedule import complex_mass_amu

from scripts.tier0_common import M_EFF_AMU, build_drag_cfg, run_dir_name


TIER1A_FIXED_TAG = "tier1a_fixed"
TIER1A_T_STAR_VALUES_PS: tuple[float, ...] = (0.5, 5.0, 9.0)


def build_anchored_cfg(
    case: str,
    variant: str,
    *,
    t_star_ps: float,
    num_molecules: int,
    ion_time_ps: float,
    dt_ion_ps: float,
    seed: int,
    coeff_overrides: Optional[Mapping[str, float]] = None,
    e_bind_override: Optional[float] = None,
) -> SimConfig:
    """Build a Tier-1a ``anchored_discrete`` config from a Tier-0 drag config.

    The underlying drag coefficients, effective binding, geometry, timestep, and
    seed all come from :func:`scripts.tier0_common.build_drag_cfg`. The only
    Tier-1a-specific changes are the mass-scenario switch, schedule onset, the
    provenance stamp, and the section-6.6 escape hatch for the intentional
    constant-``m_eff`` coefficient pairing.

    ``mass_initial_amu`` is set to the n=21 complex mass for metadata
    consistency. The run-time source of truth remains
    ``simulation/ion_initial_state.py``, which owns the anchored initial-mass
    override.
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
        mass_scenario="anchored_discrete",
        t_star_ps=float(t_star_ps),
        anchor_mode="time",
        coulomb_available_eV=0.80,
        allow_inconsistent_mass_pairing=True,
        mass_initial_amu=complex_mass_amu(21),
    )
    cfg.validate()
    return cfg


def tier1a_run_tag(scenario: str, t_star_ps: float | None) -> str:
    """Return the run-tag suffix shared by the generator and scorer."""
    if scenario == "fixed":
        if t_star_ps is not None:
            raise ValueError("fixed Tier-1a run does not take t_star_ps")
        return TIER1A_FIXED_TAG
    if scenario == "anchored_discrete":
        if t_star_ps is None:
            raise ValueError("anchored_discrete Tier-1a run requires t_star_ps")
        return f"tier1a_anchored_continuous_t{float(t_star_ps):.1f}"
    raise ValueError(
        "scenario must be 'fixed' or 'anchored_discrete', "
        f"got {scenario!r}"
    )


def tier1a_run_dir_name(
    case: str,
    variant: str,
    n: int,
    scenario: str,
    t_star_ps: float | None,
) -> str:
    """Return the Tier-0-style run directory basename for a Tier-1a run."""
    return run_dir_name(
        case,
        variant,
        n,
        run_tag=tier1a_run_tag(scenario, t_star_ps),
    )


__all__ = [
    "M_EFF_AMU",
    "TIER1A_FIXED_TAG",
    "TIER1A_T_STAR_VALUES_PS",
    "build_anchored_cfg",
    "tier1a_run_dir_name",
    "tier1a_run_tag",
]
