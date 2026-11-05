"""Shared Tier-0 run wiring for the thesis drag laws.

Both ``scripts/gen_tier0_runs.py`` (generates a run) and
``scripts/post_processing/tier0_drag_comparison.py`` (scores it) use:

- :func:`build_tier0_cfg`, the Tier-0 config of a thesis law (``LQ``, ``LQ_eff``,
  ``CC`` with its cap, or an m = 21 calibration-ridge point);
- :func:`build_drag_cfg`, which loads a bundle, applies optional coefficient /
  binding / mass overrides, and wires them into the *base* preset (mirroring
  ``single_pulse_N2000_drag`` but form-parameterized -- the drag presets are
  frozen to ``linear_cubic`` / ``CC`` and must not be used here);
- the run-directory naming both scripts agree on;
- the per-case scoring windows (:data:`TIER0_WINDOWS_PS`; the window belongs to
  the droplet case, not to the law being run).

Import with ``from scripts.tier0_common import ...`` after inserting
``PROJECT_ROOT`` on ``sys.path`` (``scripts/`` is an implicit namespace package,
no ``__init__.py`` needed).
"""

from __future__ import annotations

from dataclasses import replace
import json
from typing import Mapping, Optional

from i2_helium_md import (
    DRAG_BUNDLE_CATALOG,
    load_capped_cc,
    load_drag_bundle,
    single_pulse_N2000,
    single_pulse_N2000_18Angst,
)
from i2_helium_md.config import HISTORICAL_DRAG_BUNDLE_IDS, SimConfig
from i2_helium_md.physics.drag import _REQUIRED_COEFF_KEYS
from i2_helium_md.presets import CALIBRATION_RIDGE_ROOT

# The single reference mass every Method-B bundle stamps (= presets._DRAG_M_EFF_AMU).
# The loader compares this exactly (float-eps) against each bundle's meff_amu.
M_EFF_AMU = 202.953908

# Valid droplet-geometry cases -> the non-drag base preset for each. The *drag*
# presets are form-frozen (linear_cubic / CC), so this module
# wires the drag surface onto the base presets directly instead.
BASE_PRESETS = {
    "9A": single_pulse_N2000,
    "18A": single_pulse_N2000_18Angst,
}

# In-window Tier-0 scoring range per droplet case [ps] (METHOD_B extraction
# windows: from the onset to the end of the smoothed reference trace).
TIER0_WINDOWS_PS: dict[str, tuple[float, float]] = {
    "9A": (2.67, 14.081),
    "18A": (4.54, 14.768),
}


def _validate_case(case: str) -> None:
    if case not in BASE_PRESETS:
        raise ValueError(
            f"CASE must be one of {sorted(BASE_PRESETS)}, got {case!r}"
        )


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
    m_eff_amu: Optional[float] = None,
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
        A key of ``DRAG_BUNDLE_CATALOG``.
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
    m_eff_amu : float or None
        Optional fixed ion mass [amu] (THESIS_FINAL_GOALS G11.2: the Tier-0
        extraction mass as a parameter). The cfg's ``m_eff_amu`` /
        ``mass_initial_amu`` and the run cfg's ``extraction_mass_amu`` stamp
        are all set to it (the §6.5 fixed-mass guard needs the pair equal).
        The bundle file is not touched: its ``meff_amu`` stays the mass it was
        fitted at (:data:`M_EFF_AMU`), so a run at another mass is a
        cross-mass evaluation of that law, recorded in the run's cfg.json.
        Default ``None`` = :data:`M_EFF_AMU`.

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
    _validate_case(case)
    if variant not in DRAG_BUNDLE_CATALOG:
        raise ValueError(
            f"VARIANT must be one of {sorted(DRAG_BUNDLE_CATALOG)}, got {variant!r}"
        )
    coeffs = load_drag_bundle(variant, expected_m_eff_amu=M_EFF_AMU)

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

    mass_amu = M_EFF_AMU if m_eff_amu is None else float(m_eff_amu)
    if not mass_amu > 0.0:
        raise ValueError(f"m_eff_amu must be positive, got {m_eff_amu!r}")

    coeffs = replace(
        coeffs,
        coefficients=new_coefficients,
        effective_binding_energy_I_ion_eV=e_bind,
        extraction_mass_amu=mass_amu,
    )

    # The selector (G4.3 change 6): the bundle id rides on the cfg, and a bundle
    # that stamps the gate width it was extracted under brings that gate with it
    # (one VARIANT switches bundle + gate together; the §6.5.2 guard
    # cross-checks). Unstamped (pre-G4) bundles leave the gate at its default.
    gate_fields: dict[str, object] = {}
    if coeffs.extraction_gate_steepness_angstrom is not None:
        gate_fields = dict(
            drag_spatial_gate="erf_independent",
            drag_gate_steepness=coeffs.extraction_gate_steepness_angstrom,
        )

    base_preset = BASE_PRESETS[case]
    cfg = base_preset(
        num_molecules=num_molecules,
        ion_simulation_time=ion_time_ps,
        dt_ion=dt_ion_ps,
        seed=seed,
        drag_form=coeffs.form,
        drag_coefficients=coeffs,
        drag_bundle=variant,
        mass_scenario="fixed",
        m_eff_amu=mass_amu,
        mass_initial_amu=mass_amu,
        binding_energy_I_ion_eV=e_bind,
        **gate_fields,
    )
    return cfg


def build_tier0_cfg(
    case: str,
    law: str,
    *,
    num_molecules: int,
    ion_time_ps: float,
    dt_ion_ps: float,
    seed: int,
    m_eff_amu: Optional[float] = None,
    e_bind_override: Optional[float] = None,
    ridge_c: Optional[float] = None,
) -> SimConfig:
    """The validated Tier-0 config of a thesis law.

    ``CC`` runs capped at the 3.5 A gate, as in the thesis. Its bundle was
    extracted at the legacy 14.2 A gate, so that pairing needs
    ``allow_inconsistent_gate_pairing``.

    Parameters
    ----------
    case : str
        Droplet geometry, ``"9A"`` or ``"18A"``.
    law : str
        ``"LQ"``, ``"LQ_eff"`` or ``"CC"``.
    num_molecules : int
        Ensemble size ``N``.
    ion_time_ps, dt_ion_ps : float
        Ion-stage duration and timestep [ps].
    seed : int
        RNG seed.
    m_eff_amu : float or None
        Fixed ion mass [amu]; see :func:`build_drag_cfg`.
    e_bind_override : float or None
        Effective binding [eV]; default the law's own.
    ridge_c : float or None
        With ``law="LQ"``: take a and E_bind from the ``calibration_ridge`` point
        at this c [amu/A] (the m = 21 fit, not B2). ``e_bind_override`` still wins.

    Returns
    -------
    SimConfig
        A validated drag cfg dispatching to the BAOAB drag path.

    Raises
    ------
    ValueError
        On an unknown case or law, or ``ridge_c`` with a law other than ``LQ``.
    FileNotFoundError
        If the ridge has no point at ``ridge_c``.
    """
    coeff_overrides = None
    if ridge_c is not None:
        if law != "LQ":
            raise ValueError(f"ridge_c needs law='LQ', got {law!r}")
        label = f"c{ridge_c:.3f}".replace(".", "p")
        path = CALIBRATION_RIDGE_ROOT / "ridge" / f"{label}.json"
        point = json.loads(path.read_text(encoding="utf-8"))
        coeff_overrides = {"a": point["a_amu_per_ps"], "c": point["c_amu_per_A"]}
        if e_bind_override is None:
            e_bind_override = point["e_bind_eV"]
    cfg = build_drag_cfg(
        case, law, num_molecules=num_molecules, ion_time_ps=ion_time_ps,
        dt_ion_ps=dt_ion_ps, seed=seed, coeff_overrides=coeff_overrides,
        e_bind_override=e_bind_override, m_eff_amu=m_eff_amu,
    )
    if law == "CC":
        capped = load_capped_cc()
        cfg = replace(
            cfg, drag_form=capped.form,
            drag_coefficients=replace(cfg.drag_coefficients, form=capped.form,
                                      coefficients=capped.coefficients),
            drag_spatial_gate="erf_independent", drag_gate_steepness=3.5,
            allow_inconsistent_gate_pairing=True,
        )
    cfg.validate()
    return cfg


def run_dir_name(case: str, variant: str, n: int, run_tag: str = "") -> str:
    """Build the agreed run-directory basename under ``data/runs``.

    Both scripts derive the run dir from this, so the comparison finds exactly
    what the generator wrote. ``run_tag`` (default ``""``) is the marker for a
    hand-tuned run; setting it keeps a tuned run from silently overwriting the
    committed-variant run dir. A renamed bundle keeps its historical id here, so
    the name still finds the finished runs on disk.
    """
    suffix = f"_{run_tag}" if run_tag else ""
    stem = HISTORICAL_DRAG_BUNDLE_IDS.get(variant, variant)
    return f"{case}_drag_{stem}_N{n}{suffix}"
