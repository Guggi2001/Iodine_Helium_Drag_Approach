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

import numpy as np

from .physics.constants import EV, K_B, S_ABS_EV
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
MassScenario = Literal["fixed", "biphasic", "anchored_discrete"]
AnchorMode = Literal["time", "onset_strip"]   # Tier-1a schedule anchor; radial cross-check deferred (plan §6)
NoiseForm = Literal["none", "multiplicative_local_fdt", "empirical_residual"]
NoiseCalibration = Literal["hard_sphere_variance", "tddft_residual", "strict_fdt_bath"]
NoiseGeometry = Literal["longitudinal", "isotropic", "anisotropic"]
NoiseLowVBehavior = Literal["vanish", "blend_to_isotropic"]
PickupRateForm = Literal["density_only", "sweeping", "dwell_time"]
ValidationHistogramMetric = Literal["wasserstein", "chi2", "ks"]

# Tier-2 Phase-B pickup channel (Slice P). ``PickupOccupancyCap`` selects the
# Langmuir shell-saturation factor ``(1 - n/n*)_+^p`` (``langmuir``, default) vs the
# density-only limit (``none``, cap inert). ``HeCaptureVelocity`` selects the incoming
# He velocity in the momentum-conserving capture reset: ``at_rest`` (default, ``u_He=0``)
# vs ``thermal`` (a Tier-3 rule-2 arm -- declared-but-unbuilt; the pickup step raises
# ``NotImplementedError`` if it is selected). Full member sets are declared now so the
# enum-reject guard and the deferred arms both have types to reference.
PickupOccupancyCap = Literal["langmuir", "none"]
HeCaptureVelocity = Literal["at_rest", "thermal"]

# Tier-2 Phase-B helium-density gate (Slice rho). Selects the ``rho_He/rho_bulk``
# profile that gates pickup: ``erf_complement`` (default) reuses the drag
# erf-complement surface via ``physics.helium_density.rho_he_ratio``; ``tabulated``
# is the declared sourced-profile (baseline/TDDFT) fallback whose machinery is built
# but whose data array is deferred (rule-2; CALIBRATION_MAP row 8). rho stays G2 --
# this is the surface-density occupancy gate, NOT a G2->G4 drag-gate promotion.
HeliumDensityProfile = Literal["erf_complement", "tabulated"]

# Tier-2 Phase-A dissociation ladder (Slice L). ``LadderElectronicPicture`` must
# list exactly the picture keys of ``physics.dissociation_ladder._FIRST_RUNG_EV``
# (the single source of truth for D_0(1); enforced by a parity test). The
# ``cooling_relaxed`` arm ships a provisional-average stub (rule-2; pinned at
# Phase F). ``LadderForm`` selects the Form-U sigmoid vs the tabulated fallback.
LadderForm = Literal["form_u", "tabulated"]
LadderElectronicPicture = Literal["statistical_mixture", "x2_only", "cooling_relaxed"]

# Mass<->coefficient consistency band (§6.5/§6.6). A *physical* statement -- the
# drag curve is mass-insensitive within ~1-2 He -- NOT a user knob. 8.0 amu is
# the 2-He edge (2 x 4.0026), the looser, safer-against-false-refuse choice.
# DISTINCT from the loader's exact-match provenance check (presets.py): that one
# is a plumbing identity (same number flowing two ways), this is a physics band.
_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU = 8.0

# The non-``fixed`` mass scenarios (those requiring time-resolved coefficients).
# A/B retired 2026-06-23/24 (superseded by `biphasic`); `anchored_discrete` is the
# Tier-1a anchored kinematic scenario (runs against the constant-m_eff coefficients
# under the §6.6 mid-window defence, i.e. allow_inconsistent_mass_pairing=True).
_EVOLVING_MASS_SCENARIOS = (
    "biphasic",
    "anchored_discrete",
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

    # -- Tier-1a anchored kinematic validation (Slice I*) --
    # Read only when mass_scenario == "anchored_discrete" (the driver builds the
    # He-shell schedule from t_star_ps); inert under `fixed`. The schedule sheds
    # one He per event 21 -> 14 with continuous velocity; cold-shed is diagnostic.
    t_star_ps: float = 5.0                # ps; schedule onset (n=21 held t<=t*; sweep {0.5,5,9})
    anchor_mode: AnchorMode = "time"      # time-anchored (radial cross-check deferred, plan §6)
    anchor_n_final: int = 14              # endpoint for onset_strip stress only
    coulomb_available_eV: float = 0.80    # eV; provenance stamp (d=9A); NO hard refuse (plan §8)

    # -- Tier-2 Phase-A dissociation ladder (Slice L) --
    # The Form-U ladder D_0(n) + its cumulative gate threshold Sigma(n). Read by
    # physics/dissociation_ladder.py (kappa, picture, selector) and by
    # check_ladder_config below (picture / selector reject + D_0(1) > D_floor).
    # All three fields are LIVE at Slice L (removed from the rule-2 exception
    # table). ladder_steepness kappa stays a Free knob (Phase-F calibration).
    dissociation_ladder: LadderForm = "form_u"                 # Form-U sigmoid vs tabulated fallback
    ladder_steepness: float = 1.0                              # kappa [per-unit-n]; Free (Phase-F fit)
    ladder_electronic_picture: LadderElectronicPicture = "statistical_mixture"  # sets D_0(1)

    # -- Tier-2 Phase-A solvation cooling (Slice K) --
    # Newton cooling of E_solv.struct toward the occupancy-resolved asymptote
    # E_inf(N) = -|S(N)|, with the pair/electrostriction binding split. Read by
    # physics/solvation_cooling.py and by check_solvation_cooling_config below.
    # |S| is Sourced (defaults to the constants.py anchor, single source); tau is
    # Bounded -- the [2.6,16.5] ps R8 band is a soft prior (NOT guard-enforced),
    # default = its geometric mid (Phase-F sweeps tau).
    solv_struct_asymptote_eV: float = S_ABS_EV     # |S|; Sourced (MASS K2), eV-primary 0.308
    internal_energy_cooling_tau_ps: float = 6.55   # tau [ps]; Bounded, geometric-mid of [2.6,16.5]

    # -- Tier-2 Phase-A internal-energy budget (Slice U) --
    # S2 onset partition f_int and S1 pickup retained fraction f_ret. Both Bounded
    # [0,1] (CALIBRATION rows 14/13). DECLARED-BUT-UNREAD in Phase A (rule-2): the
    # E_int budget *module* (physics/internal_energy_budget.py) takes them as kwargs;
    # the Phase-C generative driver reads these config fields. Defaults stay None --
    # the committed values are Phase-F calibration outputs (the f_int floor is
    # scenario-split, 0.21-0.35 @ 0.80 eV vs 0.065-0.10 @ 2.70 eV, so no single
    # default serves both budgets; f_ret is fit from the 9/18 A density contrast).
    # check_internal_energy_budget_config below bounds them to [0,1] WHEN set (the
    # soft ~0.2 f_int ceiling and the self-unbound floor are advisory, NOT enforced).
    internal_energy_partition_fraction: Optional[float] = None  # f_int; Bounded [0,1], Phase-F
    internal_energy_retained_fraction: Optional[float] = None   # f_ret; Bounded [0,1], Phase-F

    # -- Tier-2 Phase-B helium density gate (Slice rho) --
    # The rho_He/rho_bulk occupancy profile that gates pickup. Read by
    # check_helium_density_config below (enum reject arm) -- LIVE at Slice rho
    # (repurposed from the former Optional[object] G4 placeholder; NOT on the
    # rule-2 exception table). rho stays G2: this is the surface-density gate, not
    # a G2->G4 drag-gate promotion. Steepness is NOT a field here -- the Phase-C
    # driver passes _drag_gate_steepness(cfg) so density and drag share one surface.
    helium_density_profile: HeliumDensityProfile = "erf_complement"  # Slice rho (G2 gate)

    # -- Tier-2 Phase-B pickup channel (Slice P) --
    # The Poisson He-capture gain channel. ``pickup_rate_form`` /
    # ``pickup_occupancy_cap`` / ``he_capture_velocity`` are read at config-load by
    # check_pickup_config (enum reject) -- LIVE at Slice P. ``pickup_rate_coefficient``
    # (lambda_0, ps^-1) and ``pickup_occupancy_exponent`` (p) are DECLARED-BUT-UNREAD in
    # Phase B (rule-2): the pickup *module* (physics/pickup.py) takes them as kwargs; the
    # Phase-C generative driver reads these config fields. ``pickup_rate_coefficient``
    # defaults to 0.0 (inert -- no pickup until a scenario sets it); lambda_0 is a
    # Sourced+Bounded prior (~0.7-1.1 ps^-1, GAH25) pinned at Phase F, not here.
    # ``pickup_occupancy_exponent`` p is held fixed at 1.0 (NOT p=kappa -- the A12 tie is
    # physically inverse; freed at Phase F only if the size-dist first-shell edge demands
    # it). Renamed atomically from the retired Tier-1 mass_rate_* fields (no alias);
    # mass_relaxation_tau_ps was dead (superseded by internal_energy_cooling_tau_ps).
    pickup_rate_form: PickupRateForm = "density_only"        # Slice P (LIVE, enum guard)
    pickup_rate_coefficient: float = 0.0    # lambda_0 [ps^-1]; Bounded; Phase-C reader
    pickup_occupancy_cap: PickupOccupancyCap = "langmuir"    # Slice P (LIVE, enum guard)
    pickup_occupancy_exponent: float = 1.0  # p; fixed (not p=kappa); Phase-C reader
    he_capture_velocity: HeCaptureVelocity = "at_rest"       # Slice P declared / Phase-C read

    # -- Deferred (declared now, no Tier-0 reader; activated later) --
    noise_form: NoiseForm = "none"                       # Slice >=4 / Tier 3
    noise_calibration: NoiseCalibration = "hard_sphere_variance"   # Tier 3
    noise_geometry: NoiseGeometry = "longitudinal"       # Tier 3
    noise_low_v_behavior: NoiseLowVBehavior = "vanish"   # Tier 3 (anisotropic only)
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

    @property
    def electrostriction_binding_eV(self) -> float:
        """Full-shell electrostriction binding ``E_elec(n*)`` [eV] (derived).

        Surfaced for inspection only (Slice K): the non-positive collective binding
        marginal at full occupancy, ``-(|S| - Sigma(n*))``, computed from the
        configured picture/kappa/|S|. Not a stored or tunable field -- it is N-,
        picture- and kappa-dependent, hence a derived property.
        """
        from .physics.constants import N_STAR
        from .physics.solvation_cooling import e_electrostriction_eV
        return e_electrostriction_eV(
            N_STAR,
            picture=self.ladder_electronic_picture,
            kappa=self.ladder_steepness,
            s_abs_eV=self.solv_struct_asymptote_eV,
        )

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
        if self.anchor_mode not in ("time", "onset_strip"):
            raise ValueError(
                f"unknown anchor_mode {self.anchor_mode!r}; expected one of "
                "('time', 'onset_strip')"
            )
        if self.anchor_mode == "onset_strip":
            if (
                isinstance(self.anchor_n_final, (bool, np.bool_))
                or not isinstance(self.anchor_n_final, (int, np.integer))
            ):
                raise ValueError(
                    f"anchor_n_final must be an integer; got {self.anchor_n_final!r}."
                )
            if not (0 <= int(self.anchor_n_final) < 21):
                raise ValueError(
                    f"anchor_n_final must be in [0, 20] for onset_strip; "
                    f"got {self.anchor_n_final!r}."
                )
        check_drag_config(self)
        check_ladder_config(self)
        check_solvation_cooling_config(self)
        check_internal_energy_budget_config(self)
        check_helium_density_config(self)
        check_pickup_config(self)


# ---------------------------------------------------------------------------
# Tier-2 Phase-A dissociation-ladder config-load guard (Slice L)
# ---------------------------------------------------------------------------
# Known ladder selector forms (mirrors _KNOWN_DRAG_FORMS): Literal is not
# runtime-enforced, so this is the typo-recovery set the guard rejects against.
_KNOWN_LADDER_FORMS = ("form_u", "tabulated")


def check_ladder_config(cfg: "SimConfig") -> None:
    """Validate the dissociation-ladder surface of ``cfg`` at config-load (Slice L).

    A separate, unit-testable guard called from :meth:`SimConfig.validate`. It is a
    load-time fail-loud check (no silent clamp):

    1. **Selector reject** -- ``cfg.dissociation_ladder`` must be a known Form-U /
       tabulated selector (mirrors the ``drag_form`` typo-recovery arm).
    2. **Steepness positivity** -- ``cfg.ladder_steepness`` (kappa) must be > 0. A
       non-positive steepness inverts/flattens the Form-U cliff and drives the rung
       normalisation ``(1 - sigma(1))`` to zero (inf rungs), so it is refused at
       load (CLAUDE.md principle 4). This is the field's live read at Slice L.
    3. **Picture reject** -- ``cfg.ladder_electronic_picture`` must be a recognised
       electronic picture (resolved via ``dissociation_ladder.first_rung_d0_eV``,
       which raises ``ValueError`` on an unknown key; mirrors the ``mass_scenario``
       reject arm).
    4. **D_0(1) > D_floor** -- the resolved first rung must exceed the bulk-He floor
       for the configured picture. Defensive: the three sourced pictures all
       satisfy it, but a future picture / floor edit that violated it would be a
       physics error, caught loudly here rather than producing a non-dissipative
       inverted ladder.

    Raises
    ------
    ValueError
        On an unrecognised ladder selector, a non-positive ``ladder_steepness``, an
        unrecognised electronic picture, or a resolved ``D_0(1)`` that does not
        exceed ``D_floor``.
    """
    # Local import avoids a module-load cycle (dissociation_ladder imports only
    # constants; config is imported widely) and keeps the guard self-contained.
    from .physics.dissociation_ladder import D_FLOOR_EV, first_rung_d0_eV

    if cfg.dissociation_ladder not in _KNOWN_LADDER_FORMS:
        raise ValueError(
            f"unknown dissociation_ladder {cfg.dissociation_ladder!r}; expected "
            f"one of {_KNOWN_LADDER_FORMS}"
        )

    if not (cfg.ladder_steepness > 0.0):
        raise ValueError(
            f"ladder_steepness (kappa) must be > 0 (a non-positive steepness "
            f"inverts/flattens the Form-U cliff and diverges the rung "
            f"normalisation); got {cfg.ladder_steepness!r}"
        )

    # Picture reject (raises ValueError on an unknown key) + resolved first rung.
    d0_1_eV = first_rung_d0_eV(cfg.ladder_electronic_picture)

    if not (d0_1_eV > D_FLOOR_EV):
        raise ValueError(
            f"ladder first rung D_0(1)={d0_1_eV!r} eV must exceed the bulk-He "
            f"floor D_floor={D_FLOOR_EV!r} eV (picture="
            f"{cfg.ladder_electronic_picture!r})"
        )


# ---------------------------------------------------------------------------
# Tier-2 Phase-B helium-density config-load guard (Slice rho)
# ---------------------------------------------------------------------------
# Known density-profile selectors (mirrors _KNOWN_LADDER_FORMS): Literal is not
# runtime-enforced, so this is the typo-recovery set the guard rejects against.
_KNOWN_HELIUM_DENSITY_PROFILES = ("erf_complement", "tabulated")


def check_helium_density_config(cfg: "SimConfig") -> None:
    """Validate the helium-density-gate surface of ``cfg`` at config-load (Slice rho).

    A separate, unit-testable guard called from :meth:`SimConfig.validate`. It is a
    load-time fail-loud enum reject arm (no silent clamp): ``cfg.helium_density_profile``
    must be a known ``erf_complement`` / ``tabulated`` selector (mirrors the
    ``dissociation_ladder`` / ``drag_form`` typo-recovery arms). This is the field's
    live read at Slice rho -- it repurposes the former ``Optional[object]`` G4
    placeholder into the born-live G2 surface-density selector.

    The steepness is **not** a helium-density field: the Phase-C driver passes
    ``_drag_gate_steepness(cfg)`` so the density and drag gates share one surface
    (CALIBRATION_MAP row 5); there is nothing steepness-related to guard here.

    Raises
    ------
    ValueError
        On an unrecognised ``helium_density_profile`` selector.
    """
    if cfg.helium_density_profile not in _KNOWN_HELIUM_DENSITY_PROFILES:
        raise ValueError(
            f"unknown helium_density_profile {cfg.helium_density_profile!r}; "
            f"expected one of {_KNOWN_HELIUM_DENSITY_PROFILES}"
        )


# ---------------------------------------------------------------------------
# Tier-2 Phase-B pickup-channel config-load guard (Slice P)
# ---------------------------------------------------------------------------
# Known selectors (mirrors _KNOWN_LADDER_FORMS / _KNOWN_HELIUM_DENSITY_PROFILES):
# Literal is not runtime-enforced, so these are the typo-recovery sets the guard
# rejects against. NOTE the sets include the valid-but-unbuilt rule-2 arms
# (sweeping/dwell_time, thermal): a config carrying them round-trips (the
# declared-but-unread contract); the NotImplementedError fires LAZILY at point-of-use
# in physics/pickup.py, never at config-load. The guard rejects only genuine typos.
_KNOWN_PICKUP_RATE_FORMS = ("density_only", "sweeping", "dwell_time")
_KNOWN_PICKUP_OCCUPANCY_CAPS = ("langmuir", "none")
_KNOWN_HE_CAPTURE_VELOCITIES = ("at_rest", "thermal")


def check_pickup_config(cfg: "SimConfig") -> None:
    """Validate the pickup-channel surface of ``cfg`` at config-load (Slice P).

    A separate, unit-testable guard called from :meth:`SimConfig.validate`. It is a
    load-time fail-loud enum reject arm (no silent clamp) over the three pickup
    selectors -- ``pickup_rate_form``, ``pickup_occupancy_cap``, ``he_capture_velocity``
    -- mirroring the ``dissociation_ladder`` / ``helium_density_profile`` typo-recovery
    arms. These are the fields' live reads at Slice P.

    The guard accepts the **valid-but-unbuilt** rule-2 arms (``sweeping`` / ``dwell_time``
    / ``thermal``): they are legal enum members, so a config carrying them round-trips;
    the ``NotImplementedError`` fires lazily inside :mod:`physics.pickup` when the arm is
    actually selected (the declared-but-unread contract). Only genuine typos are refused
    here. ``pickup_rate_coefficient`` (lambda_0) and ``pickup_occupancy_exponent`` (p)
    carry no load-time bound -- they are Phase-C driver reads, priored/pinned at Phase F.

    Raises
    ------
    ValueError
        On an unrecognised ``pickup_rate_form``, ``pickup_occupancy_cap``, or
        ``he_capture_velocity`` selector.
    """
    if cfg.pickup_rate_form not in _KNOWN_PICKUP_RATE_FORMS:
        raise ValueError(
            f"unknown pickup_rate_form {cfg.pickup_rate_form!r}; "
            f"expected one of {_KNOWN_PICKUP_RATE_FORMS}"
        )
    if cfg.pickup_occupancy_cap not in _KNOWN_PICKUP_OCCUPANCY_CAPS:
        raise ValueError(
            f"unknown pickup_occupancy_cap {cfg.pickup_occupancy_cap!r}; "
            f"expected one of {_KNOWN_PICKUP_OCCUPANCY_CAPS}"
        )
    if cfg.he_capture_velocity not in _KNOWN_HE_CAPTURE_VELOCITIES:
        raise ValueError(
            f"unknown he_capture_velocity {cfg.he_capture_velocity!r}; "
            f"expected one of {_KNOWN_HE_CAPTURE_VELOCITIES}"
        )


# ---------------------------------------------------------------------------
# Tier-2 Phase-A solvation-cooling config-load guard (Slice K)
# ---------------------------------------------------------------------------
def check_solvation_cooling_config(cfg: "SimConfig") -> None:
    """Validate the solvation-cooling surface of ``cfg`` at config-load (Slice K).

    A separate, unit-testable guard called from :meth:`SimConfig.validate`. Both are
    load-time fail-loud checks (no silent clamp):

    1. **tau positivity** -- ``cfg.internal_energy_cooling_tau_ps`` must be > 0. A
       non-positive relaxation time is unphysical and inverts the Newton cooling.
       The cooling-time band [2.6, 16.5] ps is a *soft* R8 prior and is deliberately
       NOT enforced (decision 2026-06-30) -- Phase F may probe its edges.
    2. **|S| positivity** -- ``cfg.solv_struct_asymptote_eV`` must be > 0. A
       non-positive collective asymptote inverts the binding split (the
       electrostriction marginal ``-(|S(N)| - Sigma(N))`` would flip sign).

    Raises
    ------
    ValueError
        On a non-positive ``internal_energy_cooling_tau_ps`` or
        ``solv_struct_asymptote_eV``.
    """
    if not (cfg.internal_energy_cooling_tau_ps > 0.0):
        raise ValueError(
            f"internal_energy_cooling_tau_ps (tau) must be > 0 (a non-positive "
            f"relaxation time is unphysical and inverts the cooling); got "
            f"{cfg.internal_energy_cooling_tau_ps!r}"
        )

    if not (cfg.solv_struct_asymptote_eV > 0.0):
        raise ValueError(
            f"solv_struct_asymptote_eV (|S|) must be > 0 (a non-positive collective "
            f"asymptote inverts the binding split); got "
            f"{cfg.solv_struct_asymptote_eV!r}"
        )


# ---------------------------------------------------------------------------
# Tier-2 Phase-A internal-energy-budget config-load guard (Slice U)
# ---------------------------------------------------------------------------
def check_internal_energy_budget_config(cfg: "SimConfig") -> None:
    """Validate the E_int-budget surface of ``cfg`` at config-load (Slice U).

    A separate, unit-testable guard called from :meth:`SimConfig.validate`. The two
    partition fractions are **declared-but-unread** in Phase A (rule-2): the budget
    *module* takes them as kwargs; the Phase-C driver reads the config. Their defaults
    are ``None`` (committed values are Phase-F calibration outputs), so the guard is a
    **bounded-when-set** load-time fail-loud check:

    * ``None`` -> no-op (the field is not yet pinned).
    * set -> must satisfy ``0 <= f <= 1`` (the hard cap; MASS S2 / CALIBRATION rows
      14 (f_int) and 13 (f_ret)).

    The soft ~0.2 ``f_int`` ceiling and the self-unbound floor ``Sigma(n*)/E_avail``
    are **advisory, not constraints** (MASS S2) and are deliberately NOT enforced --
    Tier 2 may probe sub-floor / above-ceiling values.

    Raises
    ------
    ValueError
        On an ``internal_energy_partition_fraction`` or
        ``internal_energy_retained_fraction`` outside ``[0, 1]`` when set.
    """
    for name in (
        "internal_energy_partition_fraction",
        "internal_energy_retained_fraction",
    ):
        value = getattr(cfg, name)
        if value is None:
            continue
        if not (0.0 <= value <= 1.0):
            raise ValueError(
                f"{name} must be in [0, 1] when set (a partition fraction; the "
                f"~0.2 ceiling and self-unbound floor are advisory, not enforced); "
                f"got {value!r}"
            )


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
    elif form == LINEAR_QUADRATIC:
        a = float(c["a"])
        cc = float(c["c"])
        # §10.3 dissipativity arm (mirrors the §9.5 linear_cubic relaxation):
        # a >= 0, c >= 0, a + c > 0. gamma = g*(a + c*v) >= 0 with no
        # turnover; pure-quadratic (METHOD_B §10, the Method-A n~+2
        # hypothesis) is the a == 0, c > 0 corner.
        if a < 0.0 or cc < 0.0:
            raise ValueError(
                f"linear_quadratic drag requires a >= 0 and c >= 0 "
                f"(dissipative; a = 0 is the pure-quadratic variant, "
                f"METHOD_B §10.3), got a={a!r}, c={cc!r}"
            )
        if not (a + cc > 0.0):
            raise ValueError(
                f"linear_quadratic drag requires a + c > 0 (a zero-drag law "
                f"is not a drag law), got a={a!r}, c={cc!r}"
            )
    elif form == THRESHOLD:  # unreachable: NotImplemented upstream
        if not (float(c["F_sat"]) > 0.0 and float(c["v0"]) > 0.0):
            raise ValueError("threshold drag requires F_sat > 0, v0 > 0")
    elif form == POWER_LAW:
        C = float(c["C"])
        n = float(c["n"])
        # §10.3: C > 0 keeps the law strictly dissipative for v > 0; n >= 1
        # keeps gamma = g*C*v**(n-1) finite at v = 0 (an n < 1 law would
        # diverge at rest and activate the §3.8 drag_low_v_floor obligation
        # -- the floor stays inert, refused here instead).
        if not (C > 0.0):
            raise ValueError(
                f"power_law drag requires C > 0 (dissipative), got C={C!r}"
            )
        if not (n >= 1.0):
            raise ValueError(
                f"power_law drag requires n >= 1 (gamma finite at v = 0; "
                f"METHOD_B §10.3 -- n < 1 would need the inert §3.8 "
                f"low-velocity floor), got n={n!r}"
            )
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
