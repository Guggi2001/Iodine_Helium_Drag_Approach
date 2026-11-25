"""Preset ``SimConfig`` builders.

Each preset corresponds to one of the old ``inputfiles_*/*.m`` scripts.
Start from a preset and override whichever fields you need:

    >>> cfg = single_pulse_N2000(num_molecules=500, seed=123)
    >>> cfg = single_pulse_N2000_18Angst(seed=123)
    >>> cfg = single_pulse_droplet_distribution(seed=123)
"""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import fields, replace
from functools import partial
from pathlib import Path

from .config import SimConfig
from .physics.dissociation_ladder import rq4graded_rungs_eV
from .physics.drag import (
    DragCoefficients,
    CAPPED_CUBIC,
    LINEAR_CUBIC,
    REALIZED_FORMS,
    _REQUIRED_COEFF_KEYS,
)
from .physics.shell_schedule import complex_mass_amu

# Anchor for the frozen drag-coefficient layout (data/reference/drag/<bundle>/...).
# ``parents[1]`` is the project root: this file is i2_helium_md/i2_helium_md/
# presets.py, so parents[0] = the package dir, parents[1] = the repo root that
# holds ``data/``. Defined once here so the directory is not a scattered literal.
REFERENCE_DRAG_ROOT = Path(__file__).resolve().parents[1] / "data" / "reference" / "drag"

# Anchor for the TDDFT references: one folder per droplet case, holding the raw
# trajectory export ``<case>_All_Data.csv`` and its smoothed speed traces.
REFERENCE_TDDFT_ROOT = Path(__file__).resolve().parents[1] / "data" / "reference" / "TDDFT"
TDDFT_CASES = ("9A", "18A")


def tddft_reference_path(case: str, filename: str | None = None) -> Path:
    """Path of a TDDFT reference file of one droplet case.

    Parameters
    ----------
    case : str
        Droplet case, one of ``TDDFT_CASES``.
    filename : str, optional
        File inside the case folder. ``None`` selects the raw trajectory
        export ``<case>_All_Data.csv``.

    Returns
    -------
    Path
        ``REFERENCE_TDDFT_ROOT / case / filename``; existence is not checked,
        the loaders do that.

    Raises
    ------
    ValueError
        If ``case`` is not one of ``TDDFT_CASES``.
    """
    if case not in TDDFT_CASES:
        raise ValueError(f"unknown TDDFT case {case!r}; expected one of {TDDFT_CASES}")
    return REFERENCE_TDDFT_ROOT / case / (filename or f"{case}_All_Data.csv")

# Required keys in a fit_parameters.json, per realised form: the form's
# coefficients (from the drag module's single-source key table) plus their
# 1-sigma/sensitivity errors, and the stamped extraction mass. Errors are read
# to assert presence (provenance completeness) even where not yet consumed.
_FIT_PARAM_REQUIRED_KEYS = {
    form: keys + tuple(f"{k}_err" for k in keys) + ("meff_amu",)
    for form, keys in _REQUIRED_COEFF_KEYS.items()
    if form in REALIZED_FORMS
}

# Additional keys a Method-B (trajectory_matching) fit_parameters.json must
# carry (METHOD_B doc §6): the §6.5.1 coupled binding, the mass model it ran
# under, the calibration window, and the reference it was fit against.
_FIT_PARAM_METHOD_B_KEYS = (
    "extraction_mass_model",
    "effective_binding_energy_I_ion_eV",
    "t_start",
    "t_end",
    "reference_file",
)


def load_drag_coefficients(
    coeff_dir: Path, *, expected_m_eff_amu: float
) -> DragCoefficients:
    """Load a realised-form :class:`DragCoefficients` bundle from disk.

    Thin, content-validating, provenance-enforcing loader living in the
    presets/config layer so ``physics/`` stays I/O-free (CLAUDE.md rule 6). The
    case (9 A vs 18 A) is the droplet geometry the *caller* (a preset) already
    encodes; this loader is case-agnostic and takes only a directory.

    Reads ``coeff_dir/fit_parameters.json`` and returns a validated bundle whose
    ``extraction_mass_amu`` is stamped **from the JSON** (single source of
    truth), so the bundle's provenance can never silently diverge from the file
    it came from. The form is read from the JSON's optional ``form`` key; a
    file without one is a legacy ``linear_cubic`` bundle and loads exactly as
    before (loader contract unchanged). Realised forms and their per-form
    coefficient keys come from the drag module's single-source tables
    (``linear_cubic`` {a, b}, ``linear_quadratic`` {a, c}, ``power_law``
    {C, n}; METHOD_B §10.3).

    Parameters
    ----------
    coeff_dir : Path
        Directory containing ``fit_parameters.json`` (e.g.
        ``REFERENCE_DRAG_ROOT / "LQ"``).
    expected_m_eff_amu : float
        The reference mass the *caller* (preset) expects this case to carry.
        Compared **exactly** (float-eps) against the JSON ``meff_amu`` -- a
        plumbing/provenance identity, DISTINCT from the guard's ~8 amu physics
        band (config.py ``_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU``).

    Optional keys ``drag_gate_steepness_angstrom`` and ``bundle_id`` (written
    by post-G4 bundles) are read into the bundle's provenance stamps; absent,
    both are ``None`` (legacy).

    Returns
    -------
    DragCoefficients
        A bundle of the JSON's declared form with ``extraction_mass_amu``
        taken from the JSON. Legacy (Method-A) files yield
        ``extraction_method="force_balance"`` with
        ``extraction_mass_model="constant"`` and no stamped binding; files
        declaring ``extraction_method="trajectory_matching"`` (Method B) must
        also carry the joint-calibration provenance keys
        (``_FIT_PARAM_METHOD_B_KEYS``) and yield a bundle stamped with
        ``effective_binding_energy_I_ion_eV``.

    Raises
    ------
    FileNotFoundError
        If ``fit_parameters.json`` is absent (the message names the path).
    ValueError
        If the JSON is malformed, declares an unrealised form, is missing
        required keys (per form and per method), carries an unknown
        ``extraction_method``, or its ``meff_amu`` disagrees with
        ``expected_m_eff_amu``.
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

    # Form discrimination: absent key = legacy linear_cubic (the pre-§10
    # bundles never stamped a form; their contract is unchanged).
    form = str(raw.get("form", LINEAR_CUBIC))
    if form not in REALIZED_FORMS:
        raise ValueError(
            f"drag coefficient file {json_path} declares form {form!r}; "
            f"loadable forms are {REALIZED_FORMS} ('threshold' is reserved)"
        )
    required = _FIT_PARAM_REQUIRED_KEYS[form]

    missing = [k for k in required if k not in raw]
    if missing:
        raise ValueError(
            f"drag coefficient file {json_path} (form {form!r}) missing "
            f"required keys {missing}; need {required}"
        )

    json_m_eff = float(raw["meff_amu"])
    # Exact-match provenance identity (NOT the physics band): the same number
    # should be flowing two ways. A real difference is a wiring/provenance bug.
    if abs(json_m_eff - float(expected_m_eff_amu)) > 1e-6:
        raise ValueError(
            f"drag coefficient provenance mismatch in {json_path}: JSON "
            f"meff_amu={json_m_eff} != expected_m_eff_amu={expected_m_eff_amu}"
        )

    coefficients = {k: float(raw[k]) for k in _REQUIRED_COEFF_KEYS[form]}

    # Method discrimination: a legacy (Method-A) file carries no
    # ``extraction_method`` key and loads exactly as before -- it is NOT
    # re-stamped (Method A stays the frozen independent cross-reference; "no
    # jointly-validated binding" is honestly represented as None and handled
    # by the §6.5.1 guard). A Method-B file must declare itself and carry the
    # full joint-calibration provenance.
    # Optional G4 stamps (THESIS_FINAL_GOALS G4.3 changes 4/6): absent in every
    # pre-G4 file, which therefore loads exactly as before (None = legacy; the
    # §6.5.2 guard reads a None gate stamp as the 14.2 A legacy width).
    gate_stamp = raw.get("drag_gate_steepness_angstrom")
    stamps = {
        "extraction_gate_steepness_angstrom": (
            None if gate_stamp is None else float(gate_stamp)
        ),
        "bundle_id": raw.get("bundle_id"),
    }

    extraction_method = str(raw.get("extraction_method", "force_balance"))
    if extraction_method == "trajectory_matching":
        missing_b = [k for k in _FIT_PARAM_METHOD_B_KEYS if k not in raw]
        if missing_b:
            raise ValueError(
                f"trajectory_matching drag coefficient file {json_path} "
                f"missing required Method-B keys {missing_b}; need "
                f"{_FIT_PARAM_METHOD_B_KEYS}"
            )
        return DragCoefficients(
            form=form,
            coefficients=coefficients,
            extraction_mass_model=str(raw["extraction_mass_model"]),
            extraction_mass_amu=json_m_eff,
            extraction_method="trajectory_matching",
            effective_binding_energy_I_ion_eV=float(
                raw["effective_binding_energy_I_ion_eV"]
            ),
            **stamps,
        )
    if extraction_method != "force_balance":
        raise ValueError(
            f"drag coefficient file {json_path} carries unknown "
            f"extraction_method {extraction_method!r}; expected "
            f"'force_balance' (or absent, legacy) or 'trajectory_matching'"
        )

    return DragCoefficients(
        form=form,
        coefficients=coefficients,
        extraction_mass_model="constant",
        extraction_mass_amu=json_m_eff,  # stamped from JSON, single source.
        **stamps,
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

# Named drag-bundle catalog (THESIS_FINAL_GOALS G4.3 change 6): bundle id ->
# directory relative to ``REFERENCE_DRAG_ROOT``. The id is what
# ``SimConfig.drag_bundle`` selects (key set == ``config._DRAG_BUNDLES``, enforced
# by a test). The three thesis laws: LQ = B2 (production), LQ_eff = B3 (effective
# variant), CC = h405 (contrast law; the file holds the uncapped cubic b, see
# :func:`load_capped_cc`).
DRAG_BUNDLE_CATALOG: dict[str, str] = {
    "LQ": "LQ",
    "LQ_eff": "LQ_eff",
    "CC": "CC",
}

# The Tier-0 ridge at m = I+ + 21 He and its 9 A exit bound.
CALIBRATION_RIDGE_ROOT = REFERENCE_DRAG_ROOT / "calibration_ridge"

def resolve_drag_bundle_dir(bundle_id: str) -> Path:
    """Resolve a catalog id to its bundle directory.

    Parameters
    ----------
    bundle_id : str
        A key of :data:`DRAG_BUNDLE_CATALOG`.

    Returns
    -------
    Path
        Absolute directory expected to hold ``fit_parameters.json``.

    Raises
    ------
    ValueError
        On an unknown ``bundle_id``.
    """
    if bundle_id not in DRAG_BUNDLE_CATALOG:
        raise ValueError(
            f"drag bundle must be one of {sorted(DRAG_BUNDLE_CATALOG)}, "
            f"got {bundle_id!r}"
        )
    return REFERENCE_DRAG_ROOT / DRAG_BUNDLE_CATALOG[bundle_id]


def load_drag_bundle(bundle_id: str, *, expected_m_eff_amu: float) -> DragCoefficients:
    """Load a catalog drag bundle by id (the ``SimConfig.drag_bundle`` selector).

    Thin wrapper over :func:`resolve_drag_bundle_dir` +
    :func:`load_drag_coefficients` that also enforces the id identity: a bundle
    whose JSON stamps a ``bundle_id`` must stamp exactly the requested id
    (a misplaced or copied bundle fails loudly). The unstamped CC bundle keeps
    ``bundle_id=None``; the config guard then requires ``cfg.drag_bundle`` to be
    an unstamped id. Every bundle applies to either droplet geometry.

    Raises
    ------
    ValueError
        On an unknown id, a bundle-id mismatch, or any loader refusal.
    FileNotFoundError
        If the bundle directory has no ``fit_parameters.json``.
    """
    coeffs = load_drag_coefficients(
        resolve_drag_bundle_dir(bundle_id),
        expected_m_eff_amu=expected_m_eff_amu,
    )
    if coeffs.bundle_id is not None and coeffs.bundle_id != bundle_id:
        raise ValueError(
            f"drag bundle {bundle_id!r} resolves to a file stamped "
            f"bundle_id={coeffs.bundle_id!r}"
        )
    return coeffs


def load_capped_cc(*, expected_m_eff_amu: float = _DRAG_M_EFF_AMU) -> DragCoefficients:
    """Load the thesis law CC: the CC bundle's cubic ``b`` with its stamped cap.

    :func:`load_drag_bundle` returns the file's uncapped ``linear_cubic`` law,
    which the frozen drag presets use; this adds the ``cap`` block (``v_c``
    [A/ps], ``p_tail``) to give the ``capped_cubic`` law of the thesis runs.
    Thesis runs pair it with the 3.5 A gate under
    ``allow_inconsistent_gate_pairing``, because the bundle was extracted at the
    legacy 14.2 A gate.

    Parameters
    ----------
    expected_m_eff_amu : float
        Provenance mass the loader cross-checks [amu].

    Returns
    -------
    DragCoefficients
        Form ``capped_cubic`` with coefficients ``{b, v_c, p_tail}``; binding and
        provenance as stamped in the bundle.

    Raises
    ------
    KeyError
        If the bundle file has no ``cap`` block.
    """
    bundle = load_drag_bundle("CC", expected_m_eff_amu=expected_m_eff_amu)
    path = resolve_drag_bundle_dir("CC") / "fit_parameters.json"
    cap = json.loads(path.read_text(encoding="utf-8"))["cap"]
    return replace(
        bundle,
        form=CAPPED_CUBIC,
        coefficients={"b": float(bundle.coefficients["b"]),
                      "v_c": float(cap["v_c"]), "p_tail": float(cap["p_tail"])},
    )


# Validated preset drag bundle (METHOD_B §9.7 verdict + §10 decision,
# 2026-06-12): the cross-case shared-form joint refit passed every
# pre-registered band, and the pure-cubic variant (a = 0 exactly,
# gamma = g*b*v**2) is the chosen production candidate. One size-independent
# bundle serves both droplet cases; the spatial gate handles the surface.
# The frozen drag presets stay on this law.
_SHARED_DRAG_BUNDLE_ID = "CC"


def single_pulse_N2000_drag(**overrides) -> SimConfig:
    """:func:`single_pulse_N2000` (9 A) wired with the validated drag law.

    Loads the shared Method-B ``CC`` bundle, uncapped (trajectory-matching
    joint refit over both HeDFT cases, METHOD_B §9.7) and exposes it via
    ``drag_coefficients`` + the mass surface; the ion stage dispatches to the
    BAOAB drag path. The jointly calibrated effective binding stamped in the
    bundle is wired into ``binding_energy_I_ion_eV`` directly from the loaded
    coefficients, so the §6.5.1 drag<->binding pairing identity holds by
    construction and the config validates silently (no escape hatch). The
    hard-sphere fields are inherited unchanged from :func:`single_pulse_N2000`
    but the collision path is not entered.

    Parameters
    ----------
    **overrides
        Any ``SimConfig`` field to override from the preset default.
    """
    coeffs = load_drag_bundle(_SHARED_DRAG_BUNDLE_ID, expected_m_eff_amu=_DRAG_M_EFF_AMU)
    cfg = single_pulse_N2000(
        drag_form="linear_cubic",
        drag_coefficients=coeffs,
        mass_scenario="fixed",
        m_eff_amu=_DRAG_M_EFF_AMU,
        mass_initial_amu=_DRAG_M_EFF_AMU,
        binding_energy_I_ion_eV=coeffs.effective_binding_energy_I_ion_eV,
    )
    return replace(cfg, **overrides)


def single_pulse_N2000_18Angst_drag(**overrides) -> SimConfig:
    """:func:`single_pulse_N2000_18Angst` (18 A) wired with the validated drag law.

    The 18 A analog of :func:`single_pulse_N2000_drag`: same shared
    ``CC`` Method-B bundle (the law is size-independent by the
    §9.7 verdict), same stamped-binding wiring, on the 18 A base preset.

    Parameters
    ----------
    **overrides
        Any ``SimConfig`` field to override from the preset default.
    """
    coeffs = load_drag_bundle(_SHARED_DRAG_BUNDLE_ID, expected_m_eff_amu=_DRAG_M_EFF_AMU)
    cfg = single_pulse_N2000_18Angst(
        drag_form="linear_cubic",
        drag_coefficients=coeffs,
        mass_scenario="fixed",
        m_eff_amu=_DRAG_M_EFF_AMU,
        mass_initial_amu=_DRAG_M_EFF_AMU,
        binding_energy_I_ion_eV=coeffs.effective_binding_energy_I_ion_eV,
    )
    return replace(cfg, **overrides)


# =============================================================================
# Thesis production laws: B2 ("LQ" in the thesis), B3, h405
# =============================================================================

#: Battery members -> RNG seeds; the three laws share them.
THESIS_MEMBER_SEEDS: dict[str, int] = {
    "s1": 20260730, "s2": 20260731, "s3": 20260732,
    "s4": 20260733, "s5": 20260734, "s6": 20260729,
}

# Coulomb energy available to one ion [eV]: the normalisation of the birth
# internal energy E_int(0) = internal_energy_partition_fraction * this value.
_THESIS_COULOMB_AVAILABLE_EV = 2.7


def _thesis_battery_base(member: str) -> SimConfig:
    """Everything B2, B3 and h405 share: source, mass model, cooling, detection."""
    if member not in THESIS_MEMBER_SEEDS:
        raise KeyError(f"unknown member {member!r}; expected one of {list(THESIS_MEMBER_SEEDS)}")
    return single_pulse_N2000(
        seed=THESIS_MEMBER_SEEDS[member],
        num_molecules=1000,
        ion_simulation_time=500.005,          # [ps]; 50001 steps, the last at t_h = 500 ps
        # Source and geometry: I2 placed thermally in post-pickup droplets.
        R0_GS_angstrom=2.666,
        single_initial_position=False,
        use_single_droplet_size=False,
        droplet_size_sampler_mode="post_pickup",
        binding_energy_molecule_K=313.2,
        # Drag acts where helium is: 3.5 A density gate.
        drag_spatial_gate="erf_independent",
        drag_gate_steepness=3.5,
        # Mass model: 21 He at birth, pickup + RRK evaporation. The drag
        # coefficients were extracted at constant mass; B2 passes Tier 0
        # unchanged at 21 He, hence the accepted mass pairing.
        mass_scenario="biphasic",
        mass_initial_amu=complex_mass_amu(21),
        allow_inconsistent_mass_pairing=True,
        initial_shell_model="density_tied",
        pickup_rate_coefficient=0.9,          # [1/ps]
        evap_rrk_dof=8.0,
        evaporation_shed_convention="co_moving",
        dissociation_ladder="tabulated",
        tabulated_ladder_rungs_eV=rq4graded_rungs_eV(),
        # Ion internal energy: birth partition and cooling to the helium.
        coulomb_available_eV=_THESIS_COULOMB_AVAILABLE_EV,
        internal_energy_partition_law="sigma_proportional",
        internal_energy_retained_fraction=0.1,
        cooling_spatial_gate="density_scaled",
        # Detection after 8.53 us; the self-bound gate is lifted outside the droplet.
        detection_stage_enabled=True,
        detection_time_ps=8530000.0,
        detection_droplet_retained_policy="exclude_all_coupled",
        detection_coulomb_closure="partner_aware",
        self_bound_gate_release="decoupled",
        # Stored by the thesis runs but without effect here (E2 stage off,
        # hard-sphere path not entered).
        relaxation_time_ps=8000.0,
        relaxation_dissipation="landau_gated_drag",
        v_limit_m_per_s=58.0,
    )


def _with_thesis_law(
    base: SimConfig, coeffs: DragCoefficients, bundle: str, *,
    E0_eV: float, tau_ps: float, **extra,
) -> SimConfig:
    """``base`` under one drag law, its birth energy ``E0_eV`` and cooling time."""
    return replace(
        base,
        drag_form=coeffs.form,
        drag_coefficients=coeffs,
        drag_bundle=bundle,
        binding_energy_I_ion_eV=coeffs.effective_binding_energy_I_ion_eV,
        internal_energy_partition_fraction=E0_eV / _THESIS_COULOMB_AVAILABLE_EV,
        internal_energy_cooling_tau_ps=tau_ps,
        **extra,
    )


def drag_production_B2(member: str = "s1", **overrides) -> SimConfig:
    """B2, the thesis production law ("LQ"): the TDDFT ridge point c = 1.6 amu/A.

    Parameters
    ----------
    member : str
        Battery member, a key of :data:`THESIS_MEMBER_SEEDS` (sets the seed).
    **overrides
        Any ``SimConfig`` field to override.

    Raises
    ------
    KeyError
        On an unknown member.
    """
    coeffs = load_drag_bundle("LQ", expected_m_eff_amu=_DRAG_M_EFF_AMU)
    cfg = _with_thesis_law(_thesis_battery_base(member), coeffs, "LQ",
                           E0_eV=0.3297, tau_ps=6.4)
    return replace(cfg, **overrides)


def drag_effective_B3(member: str = "s1", **overrides) -> SimConfig:
    """B3, the effective high-tail variant (a = 19.8 amu/ps, 48.2 meV well).

    Not a TDDFT (Tier-0) law. Arguments as :func:`drag_production_B2`.
    """
    coeffs = load_drag_bundle("LQ_eff", expected_m_eff_amu=_DRAG_M_EFF_AMU)
    cfg = _with_thesis_law(_thesis_battery_base(member), coeffs, "LQ_eff",
                           E0_eV=0.3389, tau_ps=6.4)
    return replace(cfg, **overrides)


def drag_contrast_h405(member: str = "s1", **overrides) -> SimConfig:
    """h405, the capped-cubic contrast law, re-landed at the 3.5 A gate.

    The bundle was extracted at the 14.2 A gate, hence the accepted gate
    pairing. Arguments as :func:`drag_production_B2`.
    """
    coeffs = load_capped_cc()
    cfg = _with_thesis_law(_thesis_battery_base(member), coeffs, "CC",
                           E0_eV=0.424, tau_ps=4.4,
                           allow_inconsistent_gate_pairing=True)
    return replace(cfg, **overrides)


#: Thesis law name -> preset.
THESIS_PRESETS: dict[str, Callable[..., SimConfig]] = {
    "B2": drag_production_B2,
    "B3": drag_effective_B3,
    "h405": drag_contrast_h405,
}


#: Seeded hard-sphere runs behind the legacy figures: run directory -> builder.
#: The runs of the submitted figures were unseeded, hence the new names.
LEGACY_THESIS_RUNS: dict[str, Callable[[], SimConfig]] = {
    "legacy_9A_N2000": partial(single_pulse_N2000, num_molecules=2000, seed=20261125),
    "legacy_18A_N2000": partial(
        single_pulse_N2000_18Angst, num_molecules=2000, seed=20261126),
    "legacy_droplet_N20000": partial(
        single_pulse_droplet_distribution, num_molecules=20000, seed=20261127),
}


def same_model(a: SimConfig, b: SimConfig) -> list[str]:
    """Names of the fields in which two configs differ; empty = the same model.

    The labels ``drag_bundle`` and ``drag_coefficients.bundle_id`` are ignored:
    the coefficients are stored inline, so the labels carry no physics (the
    stored B3 runs are labelled ``LQ``, finished runs keep historical ids).

    Parameters
    ----------
    a, b : SimConfig
        Configs to compare, built or loaded from a ``cfg.json``.

    Returns
    -------
    list of str
        Differing field names in declaration order; values compare exactly.
    """
    def unlabelled(cfg: SimConfig) -> SimConfig:
        coeffs = cfg.drag_coefficients
        if coeffs is not None:
            coeffs = replace(coeffs, bundle_id=None)
        return replace(cfg, drag_bundle="CC", drag_coefficients=coeffs)

    def value(cfg: SimConfig, name: str):
        v = getattr(cfg, name)
        return tuple(v) if isinstance(v, list) else v    # a run loads lists as tuples

    a, b = unlabelled(a), unlabelled(b)
    return [f.name for f in fields(SimConfig) if value(a, f.name) != value(b, f.name)]
