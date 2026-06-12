"""Simulation configuration.

This replaces the ~36 MATLAB ``global`` variables scattered across
``run_simulation.m``, ``physical_constants.m`` and the various ``inputfiles_*/*.m``
preset scripts with a single strongly-typed dataclass.

Rule of thumb: every physical parameter lives here. Nothing else in the
code should have a tunable magic number.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Literal, Optional

from .physics.constants import EV, K_B
from .physics.drag import (
    DragCoefficients,
    LINEAR_CUBIC,
    LINEAR_QUADRATIC,
    POWER_LAW,
    THRESHOLD,
)


# ---------------------------------------------------------------------------
# Enumerations for readability
# ---------------------------------------------------------------------------
CollisionMode = Literal[1, 2, 3]
# 1: constant scattering probability per timestep
# 2: scatter after traveling one mean free path
# 3: scatter w/ probability sigma * dR * rho_droplet  (default, "sigma mode")

# ---------------------------------------------------------------------------
# Drag-model enum aliases (Slice 3, DRAG_PORT_DESIGN_DECISIONS.md §1-§6).
#
# Named module-scope ``Literal`` aliases mirror the existing ``CollisionMode``
# house style: a single source of truth for each enum's *members* (greppable,
# all members declared) without introducing ``enum.Enum``. The full member sets
# are declared now -- even members unreachable at Tier 0 -- because the §6.5
# guard's refusal logic references the non-``fixed`` mass scenarios, so they
# must exist as types to be referenced. Field defaults (below) pick each enum's
# *inert* member, not the design's "primary" (a default config must do nothing
# surprising).
# ---------------------------------------------------------------------------
DragForm = Literal["linear_cubic", "linear_quadratic", "threshold", "power_law"]
DragSpatialGate = Literal["density_proportional", "erf_tied", "erf_independent", "sharp"]
MassScenario = Literal["fixed", "scenario_A_accretion", "scenario_B_stripping", "biphasic"]
NoiseForm = Literal["none", "multiplicative_local_fdt", "empirical_residual"]
NoiseCalibration = Literal["hard_sphere_variance", "tddft_residual", "strict_fdt_bath"]
NoiseGeometry = Literal["longitudinal", "isotropic", "anisotropic"]
NoiseLowVBehavior = Literal["vanish", "blend_to_isotropic"]
MassRateForm = Literal["density_only", "sweeping", "dwell_time"]
ValidationHistogramMetric = Literal["wasserstein", "chi2", "ks"]

# Mass<->coefficient consistency band (§6.5/§6.6). A *physical* statement -- the
# drag curve is mass-insensitive within ~1-2 He -- NOT a user knob. 8.0 amu is
# the 2-He edge (2 x 4.0026), the looser, safer-against-false-refuse choice.
# DISTINCT from the loader's exact-match provenance check (presets.py): that one
# is a plumbing identity (same number flowing two ways), this is a physics band.
_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU = 8.0

# The non-``fixed`` mass scenarios (those requiring time-resolved coefficients).
_EVOLVING_MASS_SCENARIOS = (
    "scenario_A_accretion",
    "scenario_B_stripping",
    "biphasic",
)

# Recognised drag-form members (reuses the drag.py tags -- no duplicate string
# literals). The guard rejects anything outside this set: ``Literal`` is not
# enforced at runtime, so this is the recovery for the typo-catching given up by
# choosing named ``Literal`` aliases over ``enum.Enum``.
_KNOWN_DRAG_FORMS = (LINEAR_CUBIC, LINEAR_QUADRATIC, THRESHOLD, POWER_LAW)


@dataclass
class SimConfig:
    """All tunable parameters for a single-pulse I2-in-He-droplet simulation.

    Defaults reproduce ``inputfiles_dft_comparison/single_pulse_N2000.m``
    combined with the constants set in ``run_simulation.m``.
    """

    # ------------------------------------------------------------------
    # Reproducibility
    # ------------------------------------------------------------------
    seed: Optional[int] = None          # None -> fresh randomness each run

    # ------------------------------------------------------------------
    # Simulation mode flags
    # ------------------------------------------------------------------
    single_pulse: bool = True
    effusive_dynamics: bool = False
    debug: bool = False                 # was global DEBUG

    # ------------------------------------------------------------------
    # Time grid
    # ------------------------------------------------------------------
    t_max_neutral: float = 200.0        # ps  (hardcoded in neutral script)
    dt_neutral: float = 0.01            # ps

    ion_simulation_time: float = 20.0   # ps  (from ion script)
    dt_ion: float = 0.01                # ps

    # number of timesteps stored from neutral run for export
    num_neutral_export_timesteps: int = 40

    # ------------------------------------------------------------------
    # Laser / pump
    # ------------------------------------------------------------------
    lambda_pump_nm: float = 630.0       # nm
    fwhm_lambda_nm: float = 33.0        # nm
    E_diss_eV: float = 1.556            # dissociation energy of I2 X state

    # ------------------------------------------------------------------
    # Molecule ensemble
    # ------------------------------------------------------------------
    num_molecules: int = 2000
    R0_GS_angstrom: float = 9.0         # ground-state equilibrium distance
    deltaR0_angstrom: float = 0.0       # width of initial R distribution
    T_particles_K: float = 0.4          # translational temperature in droplet
    single_initial_position: bool = True    # all at droplet center
    partner_interaction: bool = True    # include I-I X potential

    # ------------------------------------------------------------------
    # Droplet parameters
    # ------------------------------------------------------------------
    use_single_droplet_size: bool = True
    single_droplet_size: int = 2000     # number of He atoms per droplet
    p_source_mbar: float = 40.0         # nozzle pressure (only for size dist)
    T_source_K: float = 14.0            # nozzle temperature

    # droplet solvation potential
    potential_steepness: float = 14.2                # atoms
    potential_steepness_molecule: float = 14.3324    # from DFT fit
    binding_energy_I_atom_K: float = 318.43          # K -> converted below
    binding_energy_molecule_K: float = 573.3         # K -> converted below

    # Xdip: additional Gaussian dip in X ground-state potential
    Xdip_active: bool = True

    # ------------------------------------------------------------------
    # Landau / minimum velocity
    # ------------------------------------------------------------------
    v_limit_m_per_s: float = 40.0       # Landau velocity, m/s

    # ------------------------------------------------------------------
    # Hard-sphere collisions
    # ------------------------------------------------------------------
    hard_sphere_collision_mode: CollisionMode = 3
    scattering_probability: float = 0.004              # mode 1 only
    geometric_scattering_crosssection_I: float = 30.0  # A^2, neutral I
    geometric_scattering_crosssection_Iplus: float = 2500.0  # A^2, ion I+
    scatter_mass_neutral_amu: float = 4.0              # He mass
    scatter_mass_ion_amu: float = 4.0
    sigma_dependent_on_v: bool = True   # v-dependent cross section for ions
    sigma_ion_exponent: float = -2.0    # sigma ~ v^sigma_ion_exponent

    neutral_scatter_angle_std_deg: float = 0.0
    ion_scatter_angle_std_deg: float = 0.0

    # ------------------------------------------------------------------
    # Ion-specific
    # ------------------------------------------------------------------
    binding_energy_I_ion_eV: float = 0.3
    mass_attach_probability: float = 0.09
    single_charge_ionization_allowed: bool = False
    additional_droplet_charges: int = 0
    highly_charged_iodine: bool = False
    E_coulomb_scale: float = 1.0        # scaling factor for Coulomb potential

    # ------------------------------------------------------------------
    # HeDFT "custom start" mimic
    # ------------------------------------------------------------------
    custom_DFT_start: bool = False

    # ------------------------------------------------------------------
    # Drag model (Slice 3 -- declarative + validation surface only; no
    # behavioral change. The collision path still runs until Slice 4 swaps
    # it. All ~18 fields declared now with INERT defaults; most are read by
    # the guard or Slice 4. See CLAUDE.md "Slice 3 declared-field exception"
    # for the field -> activating-slice table.
    # ------------------------------------------------------------------
    # -- Tier-0-live (read by the guard now and/or the Slice 4 stepper) --
    drag_form: DragForm = "linear_cubic"                 # guard + Slice 4 build
    drag_coefficients: Optional[DragCoefficients] = None  # guard + Slice 4 gamma
    drag_spatial_gate: DragSpatialGate = "density_proportional"  # erf-tied (§5.5)
    drag_gate_steepness: float = 14.2     # A; = potential_steepness (Slice 4 gate)
    mass_scenario: MassScenario = "fixed"                # inert member (guard)
    m_eff_amu: float = 202.953908         # amu; drag-law reference mass (§2.2)
    mass_initial_amu: float = 202.953908  # amu; = m_eff_amu under `fixed` (A-ii)
    allow_inconsistent_mass_pairing: bool = False        # refuse -> warn downgrade
    allow_unvalidated_binding_pairing: bool = False      # §6.5.1 refuse -> warn downgrade
    drag_low_v_floor: float = 0.0         # A/ps; inert for linear_cubic (power_law n<0 only)

    # -- Deferred (declared now, no Tier-0 reader; activated later) --
    noise_form: NoiseForm = "none"                       # Slice >=4 / Tier 3
    noise_calibration: NoiseCalibration = "hard_sphere_variance"   # Tier 3
    noise_geometry: NoiseGeometry = "longitudinal"       # Tier 3
    noise_low_v_behavior: NoiseLowVBehavior = "vanish"   # Tier 3 (anisotropic only)
    mass_rate_form: MassRateForm = "density_only"        # Tier 1 (scenario != fixed)
    mass_rate_coefficient: float = 0.0    # kappa0/eta0; Tier 1
    mass_relaxation_tau_ps: float = 0.0   # ps; biphasic only; Tier 1
    helium_density_profile: Optional[object] = None      # future G4 density profile
    validation_histogram_metric: ValidationHistogramMetric = "wasserstein"  # Tier 2

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    output_dir: str = "results"
    data_dir: str = "data/reference"

    # ==================================================================
    # Derived quantities (do not set by hand)
    # ==================================================================
    @property
    def num_timesteps_neutral(self) -> int:
        """Number of leapfrog steps in the neutral propagation stage."""
        import math
        return math.ceil(self.t_max_neutral / self.dt_neutral)

    @property
    def v_limit_angstrom_per_ps(self) -> float:
        """Landau velocity converted to A/ps (1 A/ps = 100 m/s)."""
        return self.v_limit_m_per_s / 100.0

    @property
    def binding_energy_I_atom_eV(self) -> float:
        """binding_energy_I_atom_K converted to eV (MATLAB: 318.43 * k_B / eV)."""
        return self.binding_energy_I_atom_K * K_B / EV

    @property
    def binding_energy_molecule_meV(self) -> float:
        """binding_energy_molecule_K converted to meV (MATLAB: 573.3 * k_B / eV * 1000)."""
        return self.binding_energy_molecule_K * K_B / EV * 1000.0

    @property
    def E_min_eV(self) -> float:
        """Minimum allowed kinetic energy (Landau cutoff), eV.

        MATLAB: E_min = (127*u) * v_limit^2 / 2 / eV, with v_limit in A/ps.
        """
        from .physics.constants import U, EV, MASS_I_AMU
        v = self.v_limit_angstrom_per_ps * 100.0          # back to m/s
        return (MASS_I_AMU * U) * v ** 2 / 2.0 / EV

    # ==================================================================
    # Validation
    # ==================================================================
    def validate(self) -> None:
        """Sanity checks — fail fast rather than produce garbage."""
        if self.E_min_eV > self.binding_energy_I_atom_eV:
            # matches the MATLAB warning "all neutrals will escape!"
            import warnings
            warnings.warn(
                f"E_min ({self.E_min_eV*1000:.2f} meV) > binding energy "
                f"({self.binding_energy_I_atom_eV*1000:.2f} meV): all neutrals "
                "will escape the droplet.",
                RuntimeWarning,
            )
        if self.num_molecules <= 0:
            raise ValueError("num_molecules must be positive")
        if self.dt_neutral <= 0 or self.dt_ion <= 0:
            raise ValueError("timesteps must be positive")
        if self.hard_sphere_collision_mode not in (1, 2, 3):
            raise ValueError("hard_sphere_collision_mode must be 1, 2, or 3")
        check_drag_config(self)


# ---------------------------------------------------------------------------
# Drag config-load guard (Slice 3, DRAG_PORT_DESIGN_DECISIONS.md §3.3 + §6.5)
# ---------------------------------------------------------------------------
def check_drag_config(cfg: "SimConfig") -> None:
    """Validate the drag surface of ``cfg`` at config-load (§3.3 + §6.5).

    A separate, unit-testable function called from :meth:`SimConfig.validate`.

    First it **always** rejects an unrecognised ``cfg.drag_form`` (the runtime
    recovery for the typo-catching given up by named ``Literal`` aliases vs.
    ``enum.Enum``; mirrors the ``mass_scenario`` reject arm). This runs before
    the ``drag_coefficients is None`` early return, so a typo is caught even on a
    non-drag config; the default ``"linear_cubic"`` passes unchanged.

    The remaining two checks no-op when ``cfg.drag_coefficients is None`` (the
    inert default -- a default config must do nothing surprising). When
    coefficients are present, ``coeffs.form`` must equal ``cfg.drag_form`` (the
    Slice 4 stepper builds from ``cfg.drag_form`` while consuming ``coeffs``, so
    they must name the same form), then:

    1. **Per-form dissipativity (§3.3, relaxed per METHOD_B §9.5)** -- the
       only live, non-vacuous refusal on Tier-0-reachable input. For
       ``linear_cubic`` requires ``a >= 0`` (``a = 0`` is the pure-cubic
       variant: ``gamma = g*(a + b*v^2) >= 0`` still vanishes only at
       ``v = 0`` where no energy can be added, so it is strictly dissipative
       -- but only while ``b > 0``, so ``a = 0`` with ``b <= 0`` is refused)
       and asserts there is no real turnover speed ``v_dagger = sqrt(-a/b)``
       (true whenever ``b > 0``, which both extracted cases satisfy; the
       max-trajectory-speed ceiling that a ``b < 0`` re-extraction would need
       is unsourced and recorded as a §7 open item, not invented here). The
       reserved forms' branches are written for completeness but unreachable --
       ``physics/drag.py`` raises ``NotImplementedError`` for them upstream.
    2. **Mass <-> coefficient consistency (§6.5)** -- ``fixed`` is
       self-consistent only with ``constant``-mass coefficients whose
       ``extraction_mass_amu`` matches ``m_eff_amu`` within
       :data:`_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU` (~2 He); the evolving
       scenarios require ``time_resolved`` coefficients. An inconsistent pairing
       is a hard error unless ``allow_inconsistent_mass_pairing`` downgrades it
       to a warning.
    3. **Drag <-> binding consistency (§6.5.1)** -- the coefficients and the
       effective droplet binding are a jointly-calibrated coupled pair. A
       bundle stamped with ``effective_binding_energy_I_ion_eV`` requires
       ``cfg.binding_energy_I_ion_eV`` to equal the stamp exactly (1e-9 eV --
       a wiring identity, not a physics band); an unstamped (legacy /
       Method-A) bundle is refused outright because its pairing was never
       jointly validated. Either refusal downgrades to a ``RuntimeWarning``
       under ``allow_unvalidated_binding_pairing`` (used deliberately by the
       Method-B fit loop, which by construction evaluates not-yet-validated
       pairings).

    Parameters
    ----------
    cfg : SimConfig
        The configuration to validate.

    Raises
    ------
    ValueError
        On an unrecognised ``drag_form``, a ``drag_form``/coefficient-form
        mismatch, a non-dissipative form, or an inconsistent mass<->coefficient
        pairing when ``allow_inconsistent_mass_pairing`` is ``False``.
    """
    # Typo-recovery arm: Literal is not runtime-enforced. Runs unconditionally
    # so a typo'd drag_form fails loudly even on a non-drag config.
    if cfg.drag_form not in _KNOWN_DRAG_FORMS:
        raise ValueError(
            f"unknown drag_form {cfg.drag_form!r}; expected one of "
            f"{_KNOWN_DRAG_FORMS}"
        )

    coeffs = cfg.drag_coefficients
    if coeffs is None:
        return  # inert default: nothing further to validate.

    # The stepper (Slice 4) builds from cfg.drag_form while consuming coeffs, so
    # the two must name the same form.
    if coeffs.form != cfg.drag_form:
        raise ValueError(
            f"drag_form {cfg.drag_form!r} does not match coefficient form "
            f"{coeffs.form!r}"
        )

    # --- 1. Per-form dissipativity (§3.3) ---
    form = cfg.drag_form  # == coeffs.form (cross-checked above); §2.1 "reads drag_form"
    c = coeffs.coefficients
    if form == LINEAR_CUBIC:
        a = float(c["a"])
        b = float(c["b"])
        if a < 0.0:
            raise ValueError(
                f"linear_cubic drag requires a >= 0 (dissipative at low v; "
                f"a = 0 is the METHOD_B §9.5 pure-cubic variant), got a={a!r}"
            )
        if a == 0.0 and not (b > 0.0):
            raise ValueError(
                f"linear_cubic with a == 0 (pure-cubic variant, METHOD_B "
                f"§9.5) requires b > 0, got b={b!r}"
            )
        # Turnover v_dagger = sqrt(-a/b): real only if b < 0. b > 0 => no real
        # turnover => vacuously dissipative everywhere (assert-and-skip, §3.2).
        # A b < 0 re-extraction would need a max-trajectory-speed ceiling to
        # bound v_dagger against; that ceiling is unsourced (§7 open item) and
        # is deliberately NOT invented here.
        assert b > 0.0 or a > 0.0  # documents the intent (see refusals above)
    elif form == LINEAR_QUADRATIC:  # unreachable: NotImplemented upstream
        if not (float(c["a"]) > 0.0 and float(c["c"]) >= 0.0):
            raise ValueError("linear_quadratic drag requires a > 0, c >= 0")
    elif form == THRESHOLD:  # unreachable: NotImplemented upstream
        if not (float(c["F_sat"]) > 0.0 and float(c["v0"]) > 0.0):
            raise ValueError("threshold drag requires F_sat > 0, v0 > 0")
    elif form == POWER_LAW:  # unreachable: NotImplemented upstream
        if not (float(c["gamma"]) > 0.0):
            raise ValueError("power_law drag requires gamma > 0")
    else:  # pragma: no cover -- membership already enforced above
        raise ValueError(f"unknown drag form {form!r}")

    # --- 2. Mass <-> coefficient consistency (§6.5) ---
    if cfg.mass_scenario == "fixed":
        consistent = (
            coeffs.extraction_mass_model == "constant"
            and abs(coeffs.extraction_mass_amu - cfg.m_eff_amu)
            <= _MASS_COEFFICIENT_CONSISTENCY_TOL_AMU
        )
        detail = (
            f"mass_scenario='fixed' requires constant-mass coefficients within "
            f"{_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU} amu of m_eff_amu="
            f"{cfg.m_eff_amu}; got extraction_mass_model="
            f"{coeffs.extraction_mass_model!r}, extraction_mass_amu="
            f"{coeffs.extraction_mass_amu}"
        )
    elif cfg.mass_scenario in _EVOLVING_MASS_SCENARIOS:
        consistent = coeffs.extraction_mass_model == "time_resolved"
        detail = (
            f"mass_scenario={cfg.mass_scenario!r} requires time_resolved "
            f"coefficients; got extraction_mass_model="
            f"{coeffs.extraction_mass_model!r}"
        )
    else:  # pragma: no cover -- Literal type forbids other members
        raise ValueError(f"unknown mass_scenario {cfg.mass_scenario!r}")

    if not consistent:
        msg = f"inconsistent mass<->coefficient pairing: {detail}"
        if cfg.allow_inconsistent_mass_pairing:
            warnings.warn(msg, RuntimeWarning)
        else:
            raise ValueError(
                msg + " (set allow_inconsistent_mass_pairing=True to override)"
            )

    # --- 3. Drag <-> binding consistency (§6.5.1) ---
    # The drag coefficients and the effective droplet binding are a
    # jointly-calibrated COUPLED PAIR (TIER0_FINDINGS "Correct drag traps the
    # ions"): an in-window-correct drag delivers sub-barrier surface KE, so the
    # binding the run uses must be the one the coefficients were jointly
    # validated with. Enforced like the §6.5 mass arm.
    stamped_binding = coeffs.effective_binding_energy_I_ion_eV
    if stamped_binding is None:
        binding_msg = (
            "drag<->binding pairing not jointly validated (§6.5.1): the "
            "coefficient bundle carries no effective_binding_energy_I_ion_eV "
            f"stamp (extraction_method={coeffs.extraction_method!r}), so "
            f"binding_energy_I_ion_eV={cfg.binding_energy_I_ion_eV} was never "
            "validated against these coefficients"
        )
    elif abs(cfg.binding_energy_I_ion_eV - stamped_binding) > 1e-9:
        # Exact pairing: presets copy the value from the bundle stamp, so a
        # real difference is a wiring bug, not a physics band.
        binding_msg = (
            "drag<->binding pairing mismatch (§6.5.1): "
            f"cfg.binding_energy_I_ion_eV={cfg.binding_energy_I_ion_eV} != "
            f"stamped effective_binding_energy_I_ion_eV={stamped_binding}"
        )
    else:
        binding_msg = None

    if binding_msg is not None:
        if cfg.allow_unvalidated_binding_pairing:
            warnings.warn(binding_msg, RuntimeWarning)
        else:
            raise ValueError(
                binding_msg
                + " (set allow_unvalidated_binding_pairing=True to override)"
            )
