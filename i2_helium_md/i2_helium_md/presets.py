"""Preset ``SimConfig`` builders.

Each preset corresponds to one of the old ``inputfiles_*/*.m`` scripts.
Start from a preset and override whichever fields you need:

    >>> cfg = single_pulse_N2000(num_molecules=500, seed=123)
    >>> cfg = single_pulse_N2000_18Angst(seed=123)
    >>> cfg = single_pulse_droplet_distribution(seed=123)
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

from .config import SimConfig
from .physics.drag import DragCoefficients, LINEAR_CUBIC

# Anchor for the frozen drag-coefficient layout (data/reference/drag/<case>/...).
# ``parents[1]`` is the project root: this file is i2_helium_md/i2_helium_md/
# presets.py, so parents[0] = the package dir, parents[1] = the repo root that
# holds ``data/``. Defined once here so the directory is not a scattered literal.
REFERENCE_DRAG_ROOT = Path(__file__).resolve().parents[1] / "data" / "reference" / "drag"

# Required keys in a linear_cubic fit_parameters.json (a/b plus the 1-sigma
# errors and the stamped extraction mass). Errors are read to assert presence
# (provenance completeness) even though Slice 3 does not yet consume them.
_FIT_PARAM_REQUIRED_KEYS = ("a", "b", "a_err", "b_err", "meff_amu")


def load_drag_coefficients(
    coeff_dir: Path, *, expected_m_eff_amu: float
) -> DragCoefficients:
    """Load a ``linear_cubic`` :class:`DragCoefficients` bundle from disk.

    Thin, content-validating, provenance-enforcing loader living in the
    presets/config layer so ``physics/`` stays I/O-free (CLAUDE.md rule 6). The
    case (9 A vs 18 A) is the droplet geometry the *caller* (a preset) already
    encodes; this loader is case-agnostic and takes only a directory.

    Reads ``coeff_dir/fit_parameters.json`` and returns a validated bundle whose
    ``extraction_mass_amu`` is stamped **from the JSON** (single source of
    truth), so the bundle's provenance can never silently diverge from the file
    it came from.

    Parameters
    ----------
    coeff_dir : Path
        Directory containing ``fit_parameters.json`` (e.g.
        ``REFERENCE_DRAG_ROOT / "9A" / "linear_and_cubic"``).
    expected_m_eff_amu : float
        The reference mass the *caller* (preset) expects this case to carry.
        Compared **exactly** (float-eps) against the JSON ``meff_amu`` -- a
        plumbing/provenance identity, DISTINCT from the guard's ~8 amu physics
        band (config.py ``_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU``).

    Returns
    -------
    DragCoefficients
        A ``linear_cubic`` bundle with ``extraction_mass_model="constant"`` and
        ``extraction_mass_amu`` taken from the JSON.

    Raises
    ------
    FileNotFoundError
        If ``fit_parameters.json`` is absent (the message names the path).
    ValueError
        If the JSON is malformed, missing required keys, or its ``meff_amu``
        disagrees with ``expected_m_eff_amu``.
    """
    json_path = coeff_dir / "fit_parameters.json"
    if not json_path.is_file():
        raise FileNotFoundError(
            f"drag coefficient file not found: {json_path}"
        )
    try:
        with open(json_path, "r", encoding="utf-8") as fh:
            raw = json.load(fh)
    except (json.JSONDecodeError, OSError) as exc:
        raise ValueError(f"malformed drag coefficient file {json_path}: {exc}")

    missing = [k for k in _FIT_PARAM_REQUIRED_KEYS if k not in raw]
    if missing:
        raise ValueError(
            f"drag coefficient file {json_path} missing required keys "
            f"{missing}; need {_FIT_PARAM_REQUIRED_KEYS}"
        )

    json_m_eff = float(raw["meff_amu"])
    # Exact-match provenance identity (NOT the physics band): the same number
    # should be flowing two ways. A real difference is a wiring/provenance bug.
    if abs(json_m_eff - float(expected_m_eff_amu)) > 1e-6:
        raise ValueError(
            f"drag coefficient provenance mismatch in {json_path}: JSON "
            f"meff_amu={json_m_eff} != expected_m_eff_amu={expected_m_eff_amu}"
        )

    return DragCoefficients(
        form=LINEAR_CUBIC,
        coefficients={"a": float(raw["a"]), "b": float(raw["b"])},
        extraction_mass_model="constant",
        extraction_mass_amu=json_m_eff,  # stamped from JSON, single source.
    )


def single_pulse_N2000(**overrides) -> SimConfig:
    """Reproduces ``inputfiles_dft_comparison/single_pulse_N2000.m``.

    This is the canonical preset for He-DFT comparison at R0 = 9 A.

    Parameters
    ----------
    **overrides
        Any ``SimConfig`` field to override from the preset default.
    """
    cfg = SimConfig(
        # --- from single_pulse_N2000.m ---
        R0_GS_angstrom=9.0,
        E_coulomb_scale=1.0,
        single_initial_position=True,
        custom_DFT_start=False,
        deltaR0_angstrom=0.0,
        T_particles_K=0.4,
        sigma_dependent_on_v=True,
        single_pulse=True,
        partner_interaction=True,
        additional_droplet_charges=0,
        highly_charged_iodine=False,
        num_molecules=2000,
        effusive_dynamics=False,
        hard_sphere_collision_mode=3,
        scattering_probability=0.004,
        geometric_scattering_crosssection_I=30.0,
        scatter_mass_neutral_amu=4.0,
        scatter_mass_ion_amu=4.0,
        geometric_scattering_crosssection_Iplus=2500.0,
        binding_energy_I_ion_eV=0.3,
        neutral_scatter_angle_std_deg=0.0,
        ion_scatter_angle_std_deg=0.0,
        mass_attach_probability=0.09,
        single_charge_ionization_allowed=False,
        use_single_droplet_size=True,
        single_droplet_size=2000,
        p_source_mbar=40.0,
        T_source_K=14.0,
        # --- from run_simulation.m ---
        Xdip_active=True,
        debug=False,
        num_neutral_export_timesteps=40,
        v_limit_m_per_s=40.0,
        sigma_ion_exponent=-2.0,
        lambda_pump_nm=630.0,
        E_diss_eV=1.556,
    )
    return replace(cfg, **overrides)


def single_pulse_N2000_18Angst(**overrides) -> SimConfig:
    """Reproduces ``inputfiles_dft_comparison/single_pulse_N2000_18Angst.m``.

    This is the fixed-droplet 18 A He-DFT comparison input file. It shares
    most settings with :func:`single_pulse_N2000`, but the MATLAB file has
    active assignments for the larger initial I-I distance, smaller ensemble,
    weaker I+ hard-sphere cross section, weaker ion binding, and lower helium
    attachment probability.

    Parameters
    ----------
    **overrides
        Any ``SimConfig`` field to override from the preset default.
    """
    cfg = single_pulse_N2000(
        # --- active differences from single_pulse_N2000_18Angst.m ---
        R0_GS_angstrom=18.0,
        num_molecules=2000,
        geometric_scattering_crosssection_Iplus= 1600.0,
        binding_energy_I_ion_eV=0.05,
        mass_attach_probability=0.005,
    )
    return replace(cfg, **overrides)


def single_pulse_droplet_distribution(**overrides) -> SimConfig:
    """Reproduces ``inputfiles_dft_comparison/single_pulse_droplet_distribution.m``.

    This preset keeps the same single-pulse ion/neutral physics as
    :func:`single_pulse_N2000`, but switches from the fixed 2000-atom
    droplet to the source-condition droplet-size sampler. It also uses
    the ground-state I2 equilibrium distance and a larger ensemble,
    matching the active assignments in the MATLAB input file.

    Parameters
    ----------
    **overrides
        Any ``SimConfig`` field to override from the preset default.
    """
    cfg = single_pulse_N2000(
        # --- active differences from single_pulse_droplet_distribution.m ---
        R0_GS_angstrom=2.666,
        E_coulomb_scale=0.8,
        single_initial_position=False,
        num_molecules=8000,
        use_single_droplet_size=False,
    )
    return replace(cfg, **overrides)


# Drag-law reference mass (~19 He, §2.2). Both extracted cases (9 A and 18 A)
# were fit at this value; it is the exact provenance the loader cross-checks.
_DRAG_M_EFF_AMU = 202.953908


def single_pulse_N2000_drag(**overrides) -> SimConfig:
    """:func:`single_pulse_N2000` (9 A) wired with the linear_cubic drag law.

    Loads and validates the frozen 9 A ``linear_cubic`` coefficients and exposes
    them via ``drag_coefficients`` + the mass surface. **No behavioral change in
    Slice 3:** the hard-sphere collision path still runs (Slice 4 swaps it); the
    coefficients are loaded and config-validated but not yet consumed by a
    stepper. The hard-sphere fields are inherited unchanged from
    :func:`single_pulse_N2000`.

    Parameters
    ----------
    **overrides
        Any ``SimConfig`` field to override from the preset default.
    """
    coeffs = load_drag_coefficients(
        REFERENCE_DRAG_ROOT / "9A" / "linear_and_cubic",
        expected_m_eff_amu=_DRAG_M_EFF_AMU,
    )
    cfg = single_pulse_N2000(
        drag_form="linear_cubic",
        drag_coefficients=coeffs,
        mass_scenario="fixed",
        m_eff_amu=_DRAG_M_EFF_AMU,
        mass_initial_amu=_DRAG_M_EFF_AMU,
        binding_energy_I_ion_eV = 0.23,
    )
    return replace(cfg, **overrides)


def single_pulse_N2000_18Angst_drag(**overrides) -> SimConfig:
    """:func:`single_pulse_N2000_18Angst` (18 A) wired with the linear_cubic drag law.

    The 18 A analog of :func:`single_pulse_N2000_drag`: builds the ``18A`` case
    path, loads and validates the frozen 18 A ``linear_cubic`` coefficients. Same
    Slice 3 caveat -- coefficients are loaded and validated, the collision path
    still runs until Slice 4.

    Parameters
    ----------
    **overrides
        Any ``SimConfig`` field to override from the preset default.
    """
    coeffs = load_drag_coefficients(
        REFERENCE_DRAG_ROOT / "18A" / "linear_and_cubic",
        expected_m_eff_amu=_DRAG_M_EFF_AMU,
    )
    cfg = single_pulse_N2000_18Angst(
        drag_form="linear_cubic",
        drag_coefficients=coeffs,
        mass_scenario="fixed",
        m_eff_amu=_DRAG_M_EFF_AMU,
        mass_initial_amu=_DRAG_M_EFF_AMU,
    )
    return replace(cfg, **overrides)
