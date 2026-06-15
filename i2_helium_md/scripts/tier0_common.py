"""Shared Tier-0 run wiring for Method-B form/parameter sweeps.

Both ``scripts/gen_tier0_runs.py`` (generates a run) and
``scripts/post_processing/tier0_drag_comparison.py`` (scores it) drive off a
``CASE`` + ``VARIANT`` selection. This module is the single source for:

- the **named catalog** of committed Method-B (``trajectory_matching``) bundles,
- resolving a catalog key to its on-disk bundle directory (with the
  geometry guard: per-case keys substitute ``CASE`` so a 9 A bundle can never
  land on an 18 A run),
- :func:`build_drag_cfg`, which loads any bundle, applies an optional inline
  hand-tuning override, and wires the resulting ``{form, coefficients, E_bind}``
  into the *base* preset (mirroring ``single_pulse_N2000_drag`` but
  form-parameterized -- the production drag presets are frozen to
  ``linear_cubic`` / ``shared_pure_cubic`` and must not be used here),
- the run-directory naming both scripts agree on, and
- the **per-case** scoring-window source (a shared bundle stamps the 9 A onset
  ``t_start=2.67`` for both cases, so the 18 A scoring window must come from the
  per-case bundle, decoupled from the variant being run).

Method-A bundles (``linear_and_cubic`` / ``power`` / ``quadratic`` /
``linear_and_quadratic``) are intentionally **not** in the catalog: they carry
no jointly-validated effective binding and would need
``allow_unvalidated_binding_pairing``. Scope here is Method B.

Import from either script with ``from scripts.tier0_common import ...`` after
inserting ``PROJECT_ROOT`` on ``sys.path`` (``scripts/`` is an implicit
namespace package, no ``__init__.py`` needed).
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Mapping, Optional

from i2_helium_md import (
    REFERENCE_DRAG_ROOT,
    load_drag_coefficients,
    single_pulse_N2000,
    single_pulse_N2000_18Angst,
)
from i2_helium_md.config import SimConfig
from i2_helium_md.physics.drag import _REQUIRED_COEFF_KEYS

# The single reference mass every Method-B bundle stamps (= presets._DRAG_M_EFF_AMU).
# The loader compares this exactly (float-eps) against each bundle's meff_amu.
M_EFF_AMU = 202.953908

# Valid droplet-geometry cases -> the non-drag base preset for each. The *drag*
# presets are form-frozen (linear_cubic / shared_pure_cubic), so this module
# wires the drag surface onto the base presets directly instead.
BASE_PRESETS = {
    "9A": single_pulse_N2000,
    "18A": single_pulse_N2000_18Angst,
}

# Named catalog: variant key -> bundle dir relative to data/reference/drag.
# ``{case}`` is substituted with the run's CASE for per-case keys, which makes
# cross-geometry application impossible by construction (the geometry guard).
# Shared keys carry no ``{case}`` and run on either geometry.
CATALOG: dict[str, str] = {
    # --- shared joint-refit bundles (run on either 9A or 18A geometry) ---
    "shared_pure_cubic": "shared/trajectory_matching/shared_pure_cubic",
    "shared_3param": "shared/trajectory_matching/shared_3param",
    "shared_pl": "shared/trajectory_matching/pl_shared_3param",
    "shared_lq": "shared/trajectory_matching/lq_shared_3param",
    "shared_lq_pq": "shared/trajectory_matching/lq_shared_pure_quadratic",
    # --- per-case single-curve bundles (resolve against the run's CASE) ---
    "percase_linear_cubic": "{case}/trajectory_matching",
    "percase_pl": "{case}/trajectory_matching/pl_shared_3param",
    "percase_lq": "{case}/trajectory_matching/lq_shared_3param",
    "percase_lq_pq": "{case}/trajectory_matching/lq_shared_pure_quadratic",
}


def _validate_case(case: str) -> None:
    if case not in BASE_PRESETS:
        raise ValueError(
            f"CASE must be one of {sorted(BASE_PRESETS)}, got {case!r}"
        )


def resolve_bundle_dir(case: str, variant: str) -> Path:
    """Resolve a (``case``, ``variant``) selection to its bundle directory.

    Parameters
    ----------
    case : str
        Droplet geometry, ``"9A"`` or ``"18A"``.
    variant : str
        A key of :data:`CATALOG`.

    Returns
    -------
    Path
        The absolute directory holding ``fit_parameters.json``. Per-case keys
        have ``{case}`` substituted, so the returned path can only ever point at
        the run's own geometry (the geometry guard).

    Raises
    ------
    ValueError
        If ``case`` or ``variant`` is unknown.
    """
    _validate_case(case)
    if variant not in CATALOG:
        raise ValueError(
            f"VARIANT must be one of {sorted(CATALOG)}, got {variant!r}"
        )
    rel = CATALOG[variant].format(case=case)
    return REFERENCE_DRAG_ROOT / Path(rel)


def build_drag_cfg(
    case: str,
    variant: str,
    *,
    num_molecules: int,
    ion_time_ps: float,
    dt_ion_ps: float,
    seed: int,
    coeff_overrides: Optional[Mapping[str, float]] = None,
    e_bind_override: Optional[float] = None,
) -> SimConfig:
    """Build a deterministic Tier-0 drag ``SimConfig`` from a committed bundle.

    Loads the bundle for (``case``, ``variant``), optionally applies an inline
    hand-tuning override, and wires the resulting ``{form, coefficients,
    E_bind}`` onto the case's base preset under the Tier-0 envelope
    (``mass_scenario="fixed"`` at ``m_eff``, no escape hatch). With no overrides
    the produced cfg is identical to the committed-bundle run.

    Parameters
    ----------
    case : str
        Droplet geometry, ``"9A"`` or ``"18A"`` (selects the base preset).
    variant : str
        A key of :data:`CATALOG`.
    num_molecules : int
        Ensemble size ``N``.
    ion_time_ps, dt_ion_ps : float
        Ion-stage duration and timestep [ps].
    seed : int
        RNG seed.
    coeff_overrides : Mapping[str, float] or None
        Optional per-coefficient overrides. Keys **must** be a subset of the
        bundle form's coefficient keys (``linear_cubic`` {a, b},
        ``linear_quadratic`` {a, c}, ``power_law`` {C, n}); an unknown key is a
        hard error. Default ``None`` = use the committed coefficients exactly.
    e_bind_override : float or None
        Optional effective-binding override [eV]. When given it is wired into
        **both** ``binding_energy_I_ion_eV`` and the bundle's stamped
        ``effective_binding_energy_I_ion_eV`` so the section-6.5.1 pairing
        identity stays exact and ``validate()`` passes with no escape hatch.
        Default ``None`` = use the bundle's stamped binding.

    Returns
    -------
    SimConfig
        A validated drag cfg dispatching to the BAOAB drag path.

    Raises
    ------
    ValueError
        On unknown case/variant, an override key not in the form's coefficient
        set, or any downstream bundle/config validation failure.
    """
    coeff_dir = resolve_bundle_dir(case, variant)
    coeffs = load_drag_coefficients(coeff_dir, expected_m_eff_amu=M_EFF_AMU)

    new_coefficients = dict(coeffs.coefficients)
    if coeff_overrides:
        allowed = set(_REQUIRED_COEFF_KEYS[coeffs.form])
        unknown = sorted(set(coeff_overrides) - allowed)
        if unknown:
            raise ValueError(
                f"coeff_overrides keys {unknown} are not coefficients of form "
                f"{coeffs.form!r}; allowed keys are {sorted(allowed)}"
            )
        new_coefficients.update({k: float(v) for k, v in coeff_overrides.items()})

    e_bind = (
        float(e_bind_override)
        if e_bind_override is not None
        else coeffs.effective_binding_energy_I_ion_eV
    )

    coeffs = replace(
        coeffs,
        coefficients=new_coefficients,
        effective_binding_energy_I_ion_eV=e_bind,
    )

    base_preset = BASE_PRESETS[case]
    cfg = base_preset(
        num_molecules=num_molecules,
        ion_simulation_time=ion_time_ps,
        dt_ion=dt_ion_ps,
        seed=seed,
        drag_form=coeffs.form,
        drag_coefficients=coeffs,
        mass_scenario="fixed",
        m_eff_amu=M_EFF_AMU,
        mass_initial_amu=M_EFF_AMU,
        binding_energy_I_ion_eV=e_bind,
    )
    return cfg


def run_dir_name(case: str, variant: str, n: int, run_tag: str = "") -> str:
    """Build the agreed run-directory basename under ``data/runs``.

    Both scripts derive the run dir from this, so the comparison finds exactly
    what the generator wrote. ``run_tag`` (default ``""``) is the marker for a
    hand-tuned run; setting it keeps a tuned run from silently overwriting the
    committed-variant run dir.
    """
    suffix = f"_{run_tag}" if run_tag else ""
    return f"{case}_drag_{variant}_N{n}{suffix}"


def window_source_dir(case: str) -> Path:
    """Directory whose ``fit_parameters.json`` carries the per-case scoring window.

    The scoring window is a **case** property, not a property of the variant
    being run: a shared bundle stamps the 9 A onset ``t_start=2.67`` for both
    cases, so the 18 A window must be read from the per-case bundle
    (9 A ``[2.67, 14.081]``, 18 A ``[4.54, 14.768]``). ``meff_amu`` is identical
    everywhere, so reading it here too is fine.
    """
    _validate_case(case)
    return REFERENCE_DRAG_ROOT / case / "trajectory_matching"
