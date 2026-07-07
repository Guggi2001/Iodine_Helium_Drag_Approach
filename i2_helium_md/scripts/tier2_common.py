"""Shared Tier-2 run wiring (Phase D bridge + Phase F campaign harness).

Tier 2 deliberately reuses the locked Tier-0 drag bundle plumbing, exactly as
Tier 1a did: :func:`build_biphasic_cfg` builds a Tier-0 drag config and swaps
only the mass scenario plus the generative-knob values. Seeded at Phase D
(Slice Z) with the single pinned representative priored point; **Phase F (Slice
F1)** extends the same builder with the campaign knobs (picture, kappa, tau,
scenario-stamped budget, the E2 relaxation stage) and adds the knob-encoding
run-tag / run-dir helpers the run-matrix generator and scoreboard share.

Back-compat contract (Phase F extends, does not repurpose): the Phase-D scripts
``gen_tier2_bridge_run.py`` / ``tier2_bridge_report.py`` import this module. The
campaign knobs are **None-sentinel pass-through** kwargs (``None`` -> ride the
config default, do not inject), so a call with only the bridge kwargs produces
the byte-identical bridge config. The delivered bridge defaults
(``lambda0_per_ps``/``f_int``/``f_ret``) are **not** changed.

Pinned representative priored point (plan §0/§4, 2026-07-02 -- chosen, not
fitted; the Phase-F campaign sweeps the real values):

* lambda_0 = 0.9 /ps (central of the 0.7-1.1 band),
* f_int = 0.5 (inside the derived 0.80 eV window; lands t_x in the GAH25 prior),
* f_ret = 0.1 (prior small, nonzero to exercise S1),
* tau = the config default 6.55 ps (geometric mid sqrt(2.6*16.5)),
* kappa = 1.0 and picture = ``statistical_mixture`` (the config defaults),
* nu = 2.42 /ps and s = 3n-3 (Sourced/Derived -- no choice, config defaults).

Campaign scope note (user, 2026-07-03): the Tier-2 calibration campaign runs the
**9 A case only** -- there is no reference shell evolution for the 18 A droplet,
so the 9/18 A density-contrast route to ``f_ret`` (plan F2/F4) is unavailable and
``f_ret`` is not identifiable from the size distribution alone. Recorded in
``drag_migration_log_tier2.md`` (Phase F entry).
"""

from __future__ import annotations

from dataclasses import replace
import typing
from typing import Mapping, Optional

from i2_helium_md.config import LadderElectronicPicture, SimConfig
from i2_helium_md.physics.shell_schedule import complex_mass_amu

from scripts.tier0_common import build_drag_cfg, run_dir_name


TIER2_BRIDGE_TAG = "tier2_bridge_biphasic"

# The pinned representative priored point (demonstration choice, not a fit).
BRIDGE_LAMBDA0_PER_PS = 0.9
BRIDGE_F_INT = 0.5
BRIDGE_F_RET = 0.1

# Experimental relaxation time of I+ in this setup (user, 2026-07-03): 8530 ns.
# The E2 relaxation stage propagates the energy-gated cascade to this timescale
# before the terminal size distribution is read (it exits early once every
# fragment is frozen, so the large cap is only a ceiling, not the executed step
# count). Stored in ps -- the config field ``relaxation_time_ps`` is in ps.
EXPERIMENTAL_RELAXATION_TIME_PS = 8530.0 * 1.0e3  # 8530 ns -> 8.53e6 ps

# Scenario-stamped Coulomb budgets (§6.5 pairing guard): 0.80 eV validation
# (d = 9 A) -> 2.70 eV production.
VALIDATION_BUDGET_EV = 0.80
PRODUCTION_BUDGET_EV = 2.70

# Filesystem-safe electronic-picture abbreviations for the run tag. The keys are
# kept in lockstep with the config ``LadderElectronicPicture`` Literal (single
# source of truth): a picture added to the enum without a tag fails loudly at
# import here, rather than silently blocking run-dir naming at campaign time.
_PICTURE_TAGS: dict[str, str] = {
    "statistical_mixture": "mix",
    "x2_only": "x2",
    "cooling_relaxed": "cool",
}
_PICTURE_ENUM: frozenset[str] = frozenset(typing.get_args(LadderElectronicPicture))
if frozenset(_PICTURE_TAGS) != _PICTURE_ENUM:
    raise RuntimeError(
        "_PICTURE_TAGS is out of sync with config.LadderElectronicPicture: "
        f"tags={sorted(_PICTURE_TAGS)} vs enum={sorted(_PICTURE_ENUM)}. "
        "Add the missing abbreviation(s) so the run tag stays well-defined."
    )


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
    picture: Optional[str] = None,
    kappa: Optional[float] = None,
    tau_ps: Optional[float] = None,
    evap_rrk_dof: Optional[float] = None,
    coulomb_available_eV: float = VALIDATION_BUDGET_EV,
    relaxation_time_ps: Optional[float] = None,
    relaxation_forces: Optional[str] = None,
    cooling_spatial_gate: Optional[str] = None,
    coeff_overrides: Optional[Mapping[str, float]] = None,
    e_bind_override: Optional[float] = None,
) -> SimConfig:
    """Build a Tier-2 ``biphasic`` config from a Tier-0 drag config.

    The drag coefficients, effective binding, geometry, timestep, and seed all
    come from :func:`scripts.tier0_common.build_drag_cfg` (the run parameters
    inherit Tier-1a exactly; plan §0). The Tier-2-specific changes are the
    mass-scenario switch, the generative-knob values (pinned priored point by
    default; the campaign sweeps them by argument), the scenario-stamped budget,
    and the §6.6 escape hatch for the intentional constant-``m_eff`` coefficient
    pairing (the same structural §6.5 trip as Tier 1a).

    Campaign knobs are **None-sentinel pass-through**: ``picture``, ``kappa``,
    ``tau_ps`` default to ``None`` and are only injected when set, so they ride
    the config defaults otherwise (keeping the Phase-D bridge byte-reproducible
    and not duplicating the config defaults in the script layer).

    Parameters
    ----------
    case, variant : str
        Droplet geometry (``"9A"``/``"18A"``) and drag-bundle key (Tier-0).
    num_molecules, ion_time_ps, dt_ion_ps, seed :
        Ensemble size, ion-stage duration/timestep [ps], and RNG seed.
    lambda0_per_ps, f_int, f_ret : float
        Pickup rate lambda_0 [1/ps] and the E_int partition/retained fractions
        (default = the pinned bridge point).
    picture : str or None
        Ladder electronic picture (``statistical_mixture``/``x2_only``/
        ``cooling_relaxed``). ``None`` -> ride the config default.
    kappa : float or None
        Ladder steepness (Free knob). ``None`` -> ride the config default.
    tau_ps : float or None
        Newton-cooling time [ps] (Bounded). ``None`` -> ride the config default.
    evap_rrk_dof : float or None
        Constant RRK effective-dof override ``s`` (the ``TIER2_STAIRCASE_PROBE_PLAN``
        Addendum-A falsification lever; guarded ``s >= 1`` at config-load, physics-live
        at the driver's evaporation step). ``None`` -> ride the config default (the
        per-n ``s = 3n-3`` convention). Diagnostic, not a production knob.
    coulomb_available_eV : float
        Scenario-stamped per-ion Coulomb budget [eV]; default 0.80 (validation),
        2.70 for production. Guards the budget<->drag/shell pairing (§6.5).
    relaxation_time_ps : float or None
        When set, enables the Phase-E E2 post-ejection relaxation stage
        (``relaxation_stage_enabled=True``) and stamps this cap [ps]. ``None`` ->
        stage stays disabled (the Phase-D bridge path; default scope unchanged).
    relaxation_forces : str or None
        Relaxation translation arm (``"coulomb"``/``"free_flight"``). ``None`` ->
        ride the config default (``coulomb`` -- the two I+ fragments still repel).
    cooling_spatial_gate : str or None
        K2 cooling spatial-gate arm (``"none"``/``"density_scaled"``). ``None`` ->
        ride the config default (``"none"``, ungated bath dissipation). The
        ``"density_scaled"`` arm attenuates the Newton-cooling drain by the local
        He density (cooling off outside the bubble) -- the total-strip capability
        lever probed by the staircase-probe total-strip A/B grid.
    coeff_overrides, e_bind_override :
        Passed through to :func:`build_drag_cfg`.

    Returns
    -------
    SimConfig
        A validated ``biphasic`` config.

    Notes
    -----
    tau/kappa/picture/nu/s and the cap/rate-form selectors stay at their config
    defaults unless overridden here. ``mass_initial_amu`` is set to the n = 21
    complex mass for metadata consistency; the run-time source of truth remains
    ``simulation/ion_initial_state.py`` (the biphasic branch).
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

    overrides: dict[str, object] = dict(
        mass_scenario="biphasic",
        coulomb_available_eV=float(coulomb_available_eV),
        pickup_rate_coefficient=float(lambda0_per_ps),
        internal_energy_partition_fraction=float(f_int),
        internal_energy_retained_fraction=float(f_ret),
        allow_inconsistent_mass_pairing=True,
        mass_initial_amu=complex_mass_amu(21),
    )
    if picture is not None:
        overrides["ladder_electronic_picture"] = picture
    if kappa is not None:
        overrides["ladder_steepness"] = float(kappa)
    if tau_ps is not None:
        overrides["internal_energy_cooling_tau_ps"] = float(tau_ps)
    if evap_rrk_dof is not None:
        overrides["evap_rrk_dof"] = float(evap_rrk_dof)
    if relaxation_time_ps is not None:
        overrides["relaxation_stage_enabled"] = True
        overrides["relaxation_time_ps"] = float(relaxation_time_ps)
    if relaxation_forces is not None:
        overrides["relaxation_forces"] = relaxation_forces
    if cooling_spatial_gate is not None:
        overrides["cooling_spatial_gate"] = cooling_spatial_gate

    cfg = replace(fixed_cfg, **overrides)
    cfg.validate()
    return cfg


def tier2_bridge_run_dir_name(case: str, variant: str, n: int) -> str:
    """Return the Tier-0-style run directory basename for the bridge run."""
    return run_dir_name(case, variant, n, run_tag=TIER2_BRIDGE_TAG)


def _budget_tag(budget_eV: float) -> str:
    """Encode a scenario budget [eV] as a two-decimal ``bNNN`` filesystem tag."""
    return f"b{int(round(float(budget_eV) * 100)):03d}"


def tier2_run_tag(
    *,
    picture: str,
    kappa: float,
    lambda0_per_ps: float,
    f_int: float,
    f_ret: float,
    tau_ps: float,
    budget_eV: float,
    evap_rrk_dof: Optional[float] = None,
    cooling_spatial_gate: Optional[str] = None,
    total_strip: bool = False,
) -> str:
    """Return the campaign run-tag encoding one grid point.

    Encodes **every knob that can vary across the campaign** so the run-dir name
    uniquely identifies the point (the F3 scoreboard attributes each terminal
    size distribution by dir name -- a collision either aborts a run or
    mis-attributes a score). Covered: budget, electronic picture, kappa, lambda_0
    (pickup rate), f_int, f_ret, tau, and (since the 2026-07-06 s_eff promotion)
    the constant RRK effective-dof override. Precision: kappa/lambda_0/f_int/
    f_ret and tau all to **two decimals** -- tau needs it for the Stage-2
    [2.6, 16.5] ps sweep (one decimal would alias the 6.55 default with a 6.5
    sweep value). The case and N are carried by :func:`tier2_run_dir_name`; the
    optional ``total_strip`` suffix marks the §6.11 secondary regime-axis variant.

    ``evap_rrk_dof`` (the promoted Bounded s_eff knob) appends ``_sNN.NN`` when
    set, before any ``_totalstrip`` suffix; ``None`` (the per-n ``s = 3n-3``
    classical arm) appends **nothing**, so pre-promotion campaign tags stay
    byte-identical. Mirrors :func:`tier2_probe_run_tag` (parity-locked by test).

    ``cooling_spatial_gate`` mirrors the probe encoder: ``"density_scaled"`` appends
    ``_cgds`` (after the s_eff suffix, before ``_totalstrip``); ``None`` / ``"none"``
    append **nothing**, so campaigns that do not sweep the gate stay byte-identical.
    Encoded here (not left probe-only) so a campaign that *does* pass the gate cannot
    collide two distinct-physics runs onto one run-dir name.

    Raises
    ------
    ValueError
        On an unrecognised ``picture`` (mirrors ``check_ladder_config``).
    """
    if picture not in _PICTURE_TAGS:
        raise ValueError(
            f"picture must be one of {sorted(_PICTURE_TAGS)}; got {picture!r}"
        )
    tag = (
        f"tier2_{_budget_tag(budget_eV)}_{_PICTURE_TAGS[picture]}"
        f"_k{float(kappa):.2f}_l{float(lambda0_per_ps):.2f}"
        f"_fi{float(f_int):.2f}_fr{float(f_ret):.2f}_tau{float(tau_ps):.2f}"
    )
    if evap_rrk_dof is not None:
        tag += f"_s{float(evap_rrk_dof):.2f}"
    if cooling_spatial_gate == "density_scaled":
        tag += "_cgds"
    if total_strip:
        tag += "_totalstrip"
    return tag


def tier2_run_dir_name(
    case: str,
    variant: str,
    n: int,
    *,
    picture: str,
    kappa: float,
    lambda0_per_ps: float,
    f_int: float,
    f_ret: float,
    tau_ps: float,
    budget_eV: float,
    evap_rrk_dof: Optional[float] = None,
    cooling_spatial_gate: Optional[str] = None,
    total_strip: bool = False,
) -> str:
    """Return the Tier-0-style run directory basename for a campaign grid point."""
    return run_dir_name(
        case,
        variant,
        n,
        run_tag=tier2_run_tag(
            picture=picture,
            kappa=kappa,
            lambda0_per_ps=lambda0_per_ps,
            f_int=f_int,
            f_ret=f_ret,
            tau_ps=tau_ps,
            budget_eV=budget_eV,
            evap_rrk_dof=evap_rrk_dof,
            cooling_spatial_gate=cooling_spatial_gate,
            total_strip=total_strip,
        ),
    )


def tier2_probe_run_tag(
    *,
    picture: str,
    kappa: float,
    lambda0_per_ps: float,
    f_int: float,
    f_ret: float,
    tau_ps: float,
    budget_eV: float,
    evap_rrk_dof: Optional[float] = None,
    cooling_spatial_gate: Optional[str] = None,
) -> str:
    """Return the staircase-capability-probe run-tag for one grid point.

    The pre-F5 probe (``TIER2_STAIRCASE_PROBE_PLAN.md``) lives in its **own
    run-dir namespace**: the tag starts ``tier2probe_`` (no ``_tier2_``
    substring anywhere in the resulting dir name), so the F3 campaign glob
    ``*_tier2_*`` can never sweep a probe run into the scoreboard, and the
    probe glob ``*_tier2probe_*`` can never match a campaign run. Knob
    encoding and precision mirror :func:`tier2_run_tag` (two decimals
    everywhere; same collision rationale). No ``total_strip`` variant — the
    probe is a 0.80 eV validation-side sweep only.

    ``evap_rrk_dof`` (the Addendum-A RRK-dof mini-probe lever) appends
    ``_sNN.NN`` when set; ``None`` (the per-n ``s = 3n-3`` default) appends
    **nothing**, so the delivered 45-run probe tags stay byte-identical.

    ``cooling_spatial_gate`` (the total-strip A/B lever) appends ``_cgds``
    for ``"density_scaled"`` and **nothing** for ``None``/``"none"`` (the
    default ungated arm), so pre-arm probe tags stay byte-identical. This
    dimension is **probe-only** and is intentionally NOT mirrored on the
    campaign encoder :func:`tier2_run_tag` -- the parity-lock test exercises
    only the base knobs + ``evap_rrk_dof``, so the two encoders stay locked on
    that shared body while the probe carries this extra axis until promotion.

    Raises
    ------
    ValueError
        On an unrecognised ``picture`` (mirrors ``check_ladder_config``).
    """
    if picture not in _PICTURE_TAGS:
        raise ValueError(
            f"picture must be one of {sorted(_PICTURE_TAGS)}; got {picture!r}"
        )
    tag = (
        f"tier2probe_{_budget_tag(budget_eV)}_{_PICTURE_TAGS[picture]}"
        f"_k{float(kappa):.2f}_l{float(lambda0_per_ps):.2f}"
        f"_fi{float(f_int):.2f}_fr{float(f_ret):.2f}_tau{float(tau_ps):.2f}"
    )
    if evap_rrk_dof is not None:
        tag += f"_s{float(evap_rrk_dof):.2f}"
    if cooling_spatial_gate == "density_scaled":
        tag += "_cgds"
    return tag


def tier2_probe_run_dir_name(
    case: str,
    variant: str,
    n: int,
    *,
    picture: str,
    kappa: float,
    lambda0_per_ps: float,
    f_int: float,
    f_ret: float,
    tau_ps: float,
    budget_eV: float,
    evap_rrk_dof: Optional[float] = None,
    cooling_spatial_gate: Optional[str] = None,
) -> str:
    """Return the Tier-0-style run directory basename for a probe grid point."""
    return run_dir_name(
        case,
        variant,
        n,
        run_tag=tier2_probe_run_tag(
            picture=picture,
            kappa=kappa,
            lambda0_per_ps=lambda0_per_ps,
            f_int=f_int,
            f_ret=f_ret,
            tau_ps=tau_ps,
            budget_eV=budget_eV,
            evap_rrk_dof=evap_rrk_dof,
            cooling_spatial_gate=cooling_spatial_gate,
        ),
    )


__all__ = [
    "BRIDGE_F_INT",
    "BRIDGE_F_RET",
    "BRIDGE_LAMBDA0_PER_PS",
    "EXPERIMENTAL_RELAXATION_TIME_PS",
    "PRODUCTION_BUDGET_EV",
    "TIER2_BRIDGE_TAG",
    "VALIDATION_BUDGET_EV",
    "build_biphasic_cfg",
    "tier2_bridge_run_dir_name",
    "tier2_probe_run_dir_name",
    "tier2_probe_run_tag",
    "tier2_run_dir_name",
    "tier2_run_tag",
]
