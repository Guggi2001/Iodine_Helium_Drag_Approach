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

from .physics.constants import EV, K_B, NU_EVAP_PER_PS, S_ABS_EV
from .physics.drag import (
    CAPPED_CUBIC,
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
DragForm = Literal[
    "linear_cubic", "linear_quadratic", "threshold", "power_law", "capped_cubic"
]
DragSpatialGate = Literal["density_proportional", "erf_tied", "erf_independent", "sharp"]
MassScenario = Literal["fixed", "biphasic", "anchored_discrete"]
AnchorMode = Literal["time", "onset_strip"]   # Tier-1a schedule anchor; radial cross-check deferred (plan §6)
NoiseForm = Literal["none", "multiplicative_local_fdt", "empirical_residual"]
NoiseCalibration = Literal["hard_sphere_variance", "tddft_residual", "strict_fdt_bath"]
NoiseGeometry = Literal["longitudinal", "isotropic", "anisotropic"]
NoiseLowVBehavior = Literal["vanish", "blend_to_isotropic"]
PickupRateForm = Literal["density_only", "sweeping", "dwell_time"]
ValidationHistogramMetric = Literal["wasserstein", "chi2", "ks"]
RelaxationForces = Literal["coulomb", "free_flight"]   # Tier-2 Phase-E relaxation translation arm
RelaxationDissipation = Literal["zero_gamma", "landau_gated_drag"]  # Tier-2 §I.11.2 item 2 arm (c)

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

# Tier-2 drag-port cooling spatial gate. Selects whether the K2 Newton-cooling drain
# (``physics.solvation_cooling.newton_cool_step``) is attenuated by the local He
# density: ``none`` (default) is the locked ungated bath dissipation (``-E_int/tau``
# everywhere, byte-identical to the pre-arm code); ``density_scaled`` scales the drain
# by ``rho_He/rho_bulk`` from the SAME erf-complement surface drag+pickup share
# (``tau_eff = tau/rho_ratio``), so cooling switches off once the complex is ejected
# into vacuum (no droplet bath to radiate into). A first-class interchangeable model
# arm (DRAG_PORT_DESIGN_DECISIONS); reuses ``drag_gate_steepness(cfg)`` -- no new
# steepness knob. A residual out-of-bubble floor ``rho_min`` is a deferred OQ.
CoolingSpatialGate = Literal["none", "density_scaled"]

# Tier-2 Phase-A dissociation ladder (Slice L). ``LadderElectronicPicture`` must
# list exactly the picture keys of ``physics.dissociation_ladder._FIRST_RUNG_EV``
# (the single source of truth for D_0(1); enforced by a parity test). The
# ``cooling_relaxed`` arm ships a provisional-average stub (rule-2; pinned at
# Phase F). ``LadderForm`` selects the Form-U sigmoid vs the tabulated fallback.
LadderForm = Literal["form_u", "tabulated"]
LadderElectronicPicture = Literal["statistical_mixture", "x2_only", "cooling_relaxed"]

# T9 leg A′ (Tier-2 plan §I.11): what the detection stage does with an ion
# that is energetically bound in the droplet well at handover (the V0-2
# droplet-retained class). ``refuse`` keeps the delivered loud guard;
# ``exclude`` classifies bound ions as ``droplet_retained`` and excludes
# them from the free-flight event loop (they can never decouple).
DetectionDropletRetainedPolicy = Literal["refuse", "exclude"]

# Slice T7 (Tier-2 plan §I.11): birth-position law for the molecule centre.
# ``boltzmann`` is the delivered thermal sampler (byte-inert default);
# ``uniform_volume`` is the 1D twin's L1 ensemble law — p(r) ∝ r² on
# [0, R − margin], hard surface margin, no smoothing. A twin-parity /
# capability arm, not a physical claim (plan §I.11.0 NB 2026-07-16).
BirthPositionLaw = Literal["boltzmann", "uniform_volume"]

# OQ-J adjudication (2026-07-17; Tier-2 plan §I.11.2 item 1, findings §4n/I59):
# the velocity/momentum convention of the generative evaporation channel.
# ``cold`` (byte-inert default) leaves the shed He at rest in the LAB frame —
# the complex keeps its full momentum, ``v -> m/(m - m_He)*v``, and the ledger
# books the (negative) reduced-mass KE injection; the delivered behaviour,
# retained as the diagnostic bound arm (a directed backward kick, exact only
# for a complex at rest — the [Nat23] regime). ``co_moving`` is the physical
# zeroth order for thermal evaporation from a moving complex: the He leaves
# co-moving, the velocity is unchanged, and the ledger books the (positive)
# ``+0.5*m_He*|v|^2`` the He carries away (Tier-1a's continuous-velocity
# path). Read by evaporation_step(_components) via biphasic_step (ion + E2
# stages) and by the detection-stage event loop — one field, three stages.
# RQ10 (co-moving + thermal recoil ε, coupled to RQ2) stays open.
EvaporationShedConvention = Literal["cold", "co_moving"]

# Slice T5 (Tier-2 plan §I.11; revives H.3b): the t = 0 shell-dressing law of
# the biphasic seed. ``full`` (byte-inert default) keeps the delivered
# 21-for-all convention — every ion born with the complete n* = 21 first
# shell. ``density_tied`` dresses each ion by the local He availability at
# its birth position, ``n_0i = round(n* * rho_He/rho_bulk(d_birth,i))``,
# through the SAME erf-complement surface the drag/pickup/cooling gates share
# (``helium_density.rho_he_ratio`` at ``drag_gate_steepness(cfg)``) — zero
# new free parameters. Read by the biphasic column-0 seed only (per-ion
# ``n_shell(0)`` / mass / ``e_bind_pair`` E_pot fold; the E_int(0) onset is
# NOT coupled here — that is Slice T6's p-law); biphasic-only by guard
# (check_initial_shell_config). A first-order occupancy statement, not
# shell-restructuring dynamics — under-dressed ions may re-fill via the live
# pickup channel (a reported twin-divergence candidate, plan §T5 boundary).
InitialShellModel = Literal["full", "density_tied"]

# Slice T6 (Tier-2 plan §I.11; the D2 p-law): how the S2 Coulomb onset E_int(0)
# couples to the T5 initial-shell dressing. ``constant`` (byte-inert default,
# p = 0) keeps the delivered per-ion onset E_int(0) = f_int * E_avail (no
# coupling). ``sigma_proportional`` (p = 1) scales the onset by the
# ladder-resolved ratio (Sigma(n0)/Sigma(n*))^1, so under-dressed births
# (n0 < n*) get less onset (the §4j direction — p = 0 over-suppresses). The
# factor is supplied by ``internal_energy_budget.sigma_partition_factor``;
# structurally inert at n0 = n* (ratio 1) regardless of law, and read by the
# biphasic column-0 seed only (the E_int(0) onset). Biphasic-only by guard
# (check_internal_energy_partition_config); a parameter-free arm, not a knob.
InternalEnergyPartitionLaw = Literal["constant", "sigma_proportional"]

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
_KNOWN_DRAG_FORMS = (
    LINEAR_CUBIC,
    LINEAR_QUADRATIC,
    THRESHOLD,
    POWER_LAW,
    CAPPED_CUBIC,
)


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
    # Slice T7: molecule-centre birth law + hard surface margin [Å].
    # The margin is read only under ``uniform_volume`` (guard-refused
    # otherwise — no silent carry); firm band {3, 4.67, 6} Å (H.2b D5).
    birth_position_law: BirthPositionLaw = "boltzmann"
    initial_position_margin_angstrom: float = 0.0
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
    # Tabulated-ladder data path (Slice T2, §I.10): the 1-indexed per-rung D_0
    # table [eV] (rungs[i-1] = D_0(i)) consumed only under
    # dissociation_ladder="tabulated" -- dissociation_ladder.resolve_ladder
    # builds the TabulatedLadder the three stages inject into every ladder
    # consumer. None (default) is inert. Validation (single-sourced in
    # resolve_ladder, called from check_ladder_config): 'tabulated' without a
    # table, or a table under 'form_u', is refused; floor >= N_STAR positive
    # finite entries (Sigma(n*) must be table-covered). Under 'tabulated',
    # kappa/picture are D_0/Sigma-dead (the table replaces the Form-U
    # parametrisation) but stay guard-checked and reportable. Rung tables are
    # built generator-side (rq4graded / floor1 -- no taper physics in the
    # package; Slice T3).
    tabulated_ladder_rungs_eV: Optional[tuple[float, ...]] = None

    # -- Tier-2 Phase-A solvation cooling (Slice K) --
    # Newton cooling of E_solv.struct toward the occupancy-resolved asymptote
    # E_inf(N) = -|S(N)|, with the pair/electrostriction binding split. Read by
    # physics/solvation_cooling.py and by check_solvation_cooling_config below.
    # |S| is Sourced (defaults to the constants.py anchor, single source); tau is
    # Bounded -- the [2.6,16.5] ps R8 band is a soft prior (NOT guard-enforced),
    # default = its geometric mid (Phase-F sweeps tau).
    solv_struct_asymptote_eV: float = S_ABS_EV     # |S|; Sourced (MASS K2), eV-primary 0.308
    internal_energy_cooling_tau_ps: float = 6.55   # tau [ps]; Bounded, geometric-mid of [2.6,16.5]
    # Cooling spatial gate arm (see CoolingSpatialGate above). Default "none" =
    # locked ungated bath dissipation (byte-identical); "density_scaled" attenuates
    # the K2 drain by rho_He/rho_bulk. Read by biphasic_step; guard-checked below.
    cooling_spatial_gate: CoolingSpatialGate = "none"   # inert member (arm default)

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
    # driver passes drag_gate_steepness(cfg) so density and drag share one surface.
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

    # -- Tier-2 Phase-B evaporation channel (Slice Q) --
    # The energy-gated, RRK-rate-limited He-loss channel. ``evap_rate_prefactor_per_ps``
    # (nu) defaults to the sourced/pinned NU_EVAP_PER_PS constants anchor (single source,
    # the Slice-K |S| precedent); Sourced, never tuned. ``evap_rrk_dof`` (s) is an OPTIONAL
    # effective-scalar override of the per-n effective_dof(n): None -> per-n mode count
    # (n=2->4, n>=3->3n-3), a set float -> a fixed s applied to the n>=2 bracket only (n=1
    # stays the direct k=nu branch), GUARDED s>=1 WHEN SET by check_evaporation_config
    # (s<1 diverges the rate). ``gate_onset_override_eV`` forces a fixed self-bound gate
    # threshold in place of the parameter-free Sigma(n): None -> Derived Sigma(n)
    # (production); a float is a DIAGNOSTIC / falsification lever that can never silently
    # enter production -- check_evaporation_config REFUSES it unless allow_gate_onset_override
    # is explicitly True (a refuse->warn provenance guard mirroring
    # allow_unvalidated_binding_pairing; R10 diagnostic-lever). All three are read by the
    # Phase-C generative driver / the evaporation module (declared-but-unread in Phase B).
    evap_rate_prefactor_per_ps: float = NU_EVAP_PER_PS  # nu [ps^-1]; Sourced (=constants anchor)
    evap_rrk_dof: Optional[float] = None          # s override; None=per-n; guarded s>=1 when set
    gate_onset_override_eV: Optional[float] = None      # None=parameter-free Sigma(n); else diagnostic
    allow_gate_onset_override: bool = False             # provenance guard for the override above
    # Shed velocity/momentum convention (OQ-J arm; see EvaporationShedConvention
    # above). "cold" = delivered momentum-conserving reset (byte-inert default);
    # "co_moving" = the physical co-moving He (T5+/leg-B legs stamp it explicitly).
    evaporation_shed_convention: EvaporationShedConvention = "cold"

    # -- Tier-2 Slice T5 initial-shell dressing (plan §I.11; H.3b revived) --
    # "full" = delivered 21-for-all (byte-inert default); "density_tied" =
    # per-ion n_0 = round(n* * rho_hat(d_birth)) through the shared
    # erf-complement surface (see InitialShellModel above). Biphasic-only
    # (check_initial_shell_config); no preset or generator selects the arm —
    # the T5+/leg-B legs stamp it explicitly (the T7 precedent).
    initial_shell_model: InitialShellModel = "full"

    # -- Tier-2 Slice T6 E_int(0)-dressing coupling (plan §I.11; the D2 p-law) --
    # "constant" = delivered onset E_int(0) = f_int * E_avail (byte-inert
    # default, p = 0); "sigma_proportional" = onset * (Sigma(n0)/Sigma(n*))^1
    # (see InternalEnergyPartitionLaw above). Structurally inert without the T5
    # dressing (ratio 1 at n0 = n*); biphasic-only
    # (check_internal_energy_partition_config); no preset or generator selects
    # the arm -- the T9 legs stamp it explicitly (the T5/T7 precedent).
    internal_energy_partition_law: InternalEnergyPartitionLaw = "constant"

    # -- Deferred (declared now, no Tier-0 reader; activated later) --
    noise_form: NoiseForm = "none"                       # Slice >=4 / Tier 3
    noise_calibration: NoiseCalibration = "hard_sphere_variance"   # Tier 3
    noise_geometry: NoiseGeometry = "longitudinal"       # Tier 3
    noise_low_v_behavior: NoiseLowVBehavior = "vanish"   # Tier 3 (anisotropic only)
    validation_histogram_metric: ValidationHistogramMetric = "wasserstein"  # Tier 2

    # -- Tier-2 Phase-E post-ejection relaxation stage (Slice E2; opt-in) --
    # The R5 mitigation: propagate the biphasic mass subsystem past the 20 ps ion
    # stage to the experimental timescale (pickup off via lambda_0=0, drag off via
    # gamma=0) so the terminal size distribution is read at matched time rather than
    # at the truncated sim-end upper bound. Default off keeps the default simulation
    # scope unchanged (forbidden-list compliant); check_relaxation_config gates them
    # and they are read by simulation/relaxation_stage.run_relaxation_stage.
    relaxation_stage_enabled: bool = False               # opt-in; default off = default scope
    relaxation_time_ps: Optional[float] = None           # required-when-enabled (no sourced t_exp)
    relaxation_dt_ps: Optional[float] = None             # None -> dt_ion; guard nu*dt <= 0.1
    relaxation_forces: RelaxationForces = "coulomb"      # translation arm (coulomb default; free_flight)
    # Dissipation arm (Tier-2 plan §I.11.2 item 2, arm (c)). "zero_gamma" (default,
    # byte-inert) is the delivered conservative closure; "landau_gated_drag" turns
    # the locked pure-cubic drag on ABOVE the Landau cutoff (cfg.v_limit) and keeps
    # sub-Landau motion frictionless -- captures the marginal droplet-retained ions
    # the zero-gamma convention leaves as artificial long-lived resonances. Read by
    # simulation/relaxation_stage; only the coulomb translation arm consumes it.
    relaxation_dissipation: RelaxationDissipation = "zero_gamma"

    # -- Tier-2 Slice DS detection-time continuation stage (opt-in) --
    # Event-driven (Gillespie) continuation of the post-ejection evaporation
    # cascade to the detector arrival time -- an *exact solver* of the delivered
    # mechanism under the P1-P3 handover guard (TIER2_DETECTION_STAGE_DESIGN.md
    # §2), not a model choice, so the enabled flag alone satisfies the
    # every-model-choice-behind-an-enum rule (design §3.4). The stage seeds from
    # relaxation.npz when E2 ran, or directly from ion.npz on the skip path
    # (design §1 item 5) -- the relaxation stage is NOT required.
    # ``detection_time_ps`` is **Sourced calibration data with no baked-in
    # default** (8.53e6 ps = 8.53 us TOF flight time, CALIBRATION_MAP row 24);
    # required-when-enabled. Both read by check_detection_config below and by
    # simulation/detection_stage.run_detection_stage.
    detection_stage_enabled: bool = False                # opt-in; default off = default scope
    detection_time_ps: Optional[float] = None            # Sourced (8.53e6); required-when-enabled
    # Droplet-retained (well-trapped) ions at handover (T9 leg A′, 2026-07-16):
    # "refuse" (default) = the delivered loud P1–P3 guard; "exclude" = classify
    # energetically bound violators `droplet_retained` per the V0-2 scoring
    # convention (excluded from the event loop and the IHe_n read; unbound
    # violators still refuse loudly).
    detection_droplet_retained_policy: DetectionDropletRetainedPolicy = "refuse"

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
        configured ladder (Form-U picture/kappa, or the tabulated rung table --
        review fix 2026-07-16: under ``dissociation_ladder='tabulated'`` the
        Form-U parametrisation is Sigma-dead and must not be reported). Not a
        stored or tunable field -- it is N- and ladder-dependent, hence a
        derived property.
        """
        from .physics.constants import N_STAR
        from .physics.dissociation_ladder import resolve_ladder
        from .physics.solvation_cooling import e_electrostriction_eV
        return e_electrostriction_eV(
            N_STAR,
            picture=self.ladder_electronic_picture,
            kappa=self.ladder_steepness,
            s_abs_eV=self.solv_struct_asymptote_eV,
            ladder=resolve_ladder(
                self.dissociation_ladder, self.tabulated_ladder_rungs_eV
            ),
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
        check_birth_position_config(self)
        check_drag_config(self)
        check_ladder_config(self)
        check_solvation_cooling_config(self)
        check_cooling_spatial_gate_config(self)
        check_internal_energy_budget_config(self)
        check_helium_density_config(self)
        check_pickup_config(self)
        check_evaporation_config(self)
        check_initial_shell_config(self)
        check_internal_energy_partition_config(self)
        check_biphasic_config(self)
        check_relaxation_config(self)
        check_detection_config(self)


# ---------------------------------------------------------------------------
# Shared string-enum typo-recovery reject arm
# ---------------------------------------------------------------------------
def _reject_unknown_enum(value: object, known: tuple, *, field: str) -> None:
    """Fail-loud on an unrecognised string-enum selector (CLAUDE.md principle 4).

    ``Literal`` type hints are not runtime-enforced, so every interchangeable
    string-enum ``SimConfig`` field (``drag_form``, ``dissociation_ladder``,
    ``helium_density_profile``, ``cooling_spatial_gate``, the pickup selectors, ...)
    is typo-guarded at config-load against its ``_KNOWN_*`` tuple. This is the single
    shared reject arm; each caller passes its own ``field`` name so the message stays
    field-specific (byte-identical to the per-field copies it replaces).

    Raises
    ------
    ValueError
        When ``value`` is not in ``known``.
    """
    if value not in known:
        raise ValueError(f"unknown {field} {value!r}; expected one of {known}")


# ---------------------------------------------------------------------------
# Slice T7 birth-position-law config-load guard (Tier-2 plan §I.11)
# ---------------------------------------------------------------------------
_KNOWN_BIRTH_POSITION_LAWS = ("boltzmann", "uniform_volume")


def check_birth_position_config(cfg: "SimConfig") -> None:
    """Validate the birth-position surface of ``cfg`` at config-load (Slice T7).

    Rules
    -----
    1. ``birth_position_law`` must be a known selector (typo guard).
    2. ``initial_position_margin_angstrom`` must be finite and >= 0.
    3. The margin belongs to the ``uniform_volume`` law only: a non-zero
       margin under ``boltzmann`` is refused loudly (no silent carry —
       CLAUDE.md rule 2 convention).
    4. ``uniform_volume`` under ``single_initial_position=True`` is refused
       (review fix 2026-07-18): ``build_initial_state`` zeroes ``r0`` after
       sampling, so the advertised law would be silently inert on positions
       while still shifting the RNG draw stream (the uniform draw count
       differs from the Boltzmann rejection sampler's) — the run would
       reproduce neither the boltzmann leg nor the uniform_volume one.

    The margin-vs-droplet-radius bound is enforced at sampling time
    (:func:`i2_helium_md.sampling.radial_positions.sample_radial_positions`),
    where the per-molecule radii are known.
    """
    _reject_unknown_enum(
        cfg.birth_position_law,
        _KNOWN_BIRTH_POSITION_LAWS,
        field="birth_position_law",
    )
    margin = cfg.initial_position_margin_angstrom
    if not np.isfinite(margin) or margin < 0.0:
        raise ValueError(
            "initial_position_margin_angstrom must be finite and >= 0; "
            f"got {margin!r}"
        )
    if cfg.birth_position_law == "boltzmann" and margin != 0.0:
        raise ValueError(
            "initial_position_margin_angstrom is defined for "
            "birth_position_law='uniform_volume' only; got "
            f"margin={margin} under 'boltzmann' (set the margin to 0.0 "
            "or select the uniform_volume law)"
        )
    if cfg.birth_position_law == "uniform_volume" and cfg.single_initial_position:
        raise ValueError(
            "birth_position_law='uniform_volume' requires "
            "single_initial_position=False: build_initial_state zeroes r0 "
            "under single_initial_position=True, so the uniform_volume law "
            "would be silently inert on positions while still shifting the "
            "RNG draw stream (set single_initial_position=False or keep the "
            "boltzmann default)"
        )


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
    5. **Tabulated data path** (Slice T2, §I.10) -- the
       ``(dissociation_ladder, tabulated_ladder_rungs_eV)`` pairing resolves
       through ``dissociation_ladder.resolve_ladder`` (single source; the stages
       call the same resolver at point-of-use): ``'tabulated'`` requires a table
       with >= ``N_STAR`` positive finite rungs, and a table under ``'form_u'``
       is refused (the off-diagonal stale-intent hazard).

    Raises
    ------
    ValueError
        On an unrecognised ladder selector, a non-positive ``ladder_steepness``, an
        unrecognised electronic picture, a resolved ``D_0(1)`` that does not
        exceed ``D_floor``, or an invalid tabulated-ladder pairing/table
        (see :func:`~i2_helium_md.physics.dissociation_ladder.resolve_ladder`).
    """
    # Local import avoids a module-load cycle (dissociation_ladder imports only
    # constants; config is imported widely) and keeps the guard self-contained.
    from .physics.dissociation_ladder import (
        D_FLOOR_EV,
        first_rung_d0_eV,
        resolve_ladder,
    )

    _reject_unknown_enum(
        cfg.dissociation_ladder, _KNOWN_LADDER_FORMS, field="dissociation_ladder"
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

    # Tabulated data path (arm 5): pairing + table validation, single-sourced in
    # the resolver the stages also call at point-of-use.
    resolve_ladder(cfg.dissociation_ladder, cfg.tabulated_ladder_rungs_eV)


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
    ``drag_gate_steepness(cfg)`` so the density and drag gates share one surface
    (CALIBRATION_MAP row 5); there is nothing steepness-related to guard here.

    Raises
    ------
    ValueError
        On an unrecognised ``helium_density_profile`` selector.
    """
    _reject_unknown_enum(
        cfg.helium_density_profile,
        _KNOWN_HELIUM_DENSITY_PROFILES,
        field="helium_density_profile",
    )


# ---------------------------------------------------------------------------
# Tier-2 drag-port cooling-spatial-gate config-load guard
# ---------------------------------------------------------------------------
# Known cooling-gate selectors (mirrors _KNOWN_HELIUM_DENSITY_PROFILES): Literal is
# not runtime-enforced, so this is the typo-recovery set the guard rejects against.
_KNOWN_COOLING_SPATIAL_GATES = ("none", "density_scaled")


def check_cooling_spatial_gate_config(cfg: "SimConfig") -> None:
    """Validate the K2 cooling-spatial-gate arm of ``cfg`` at config-load.

    A separate, unit-testable guard called from :meth:`SimConfig.validate`. It is a
    load-time fail-loud enum reject arm (no silent clamp): ``cfg.cooling_spatial_gate``
    must be a known ``none`` / ``density_scaled`` selector (mirrors the
    ``helium_density_profile`` / ``drag_form`` typo-recovery arms). This is the field's
    live read at config-load; the physics read is in
    :func:`~i2_helium_md.simulation.ion_propagation_step.biphasic_step`.

    The gate steepness is **not** a cooling field: the ``density_scaled`` arm reuses
    ``drag_gate_steepness(cfg)`` so the cooling, density, and drag gates share one
    surface (no new steepness knob). There is nothing steepness-related to guard here.

    Raises
    ------
    ValueError
        On an unrecognised ``cooling_spatial_gate`` selector.
    """
    _reject_unknown_enum(
        cfg.cooling_spatial_gate,
        _KNOWN_COOLING_SPATIAL_GATES,
        field="cooling_spatial_gate",
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
    _reject_unknown_enum(
        cfg.pickup_rate_form, _KNOWN_PICKUP_RATE_FORMS, field="pickup_rate_form"
    )
    _reject_unknown_enum(
        cfg.pickup_occupancy_cap,
        _KNOWN_PICKUP_OCCUPANCY_CAPS,
        field="pickup_occupancy_cap",
    )
    _reject_unknown_enum(
        cfg.he_capture_velocity,
        _KNOWN_HE_CAPTURE_VELOCITIES,
        field="he_capture_velocity",
    )


# ---------------------------------------------------------------------------
# Tier-2 Phase-B evaporation-channel config-load guard (Slice Q)
# ---------------------------------------------------------------------------
def check_evaporation_config(cfg: "SimConfig") -> None:
    """Validate the evaporation-channel surface of ``cfg`` at config-load (Slice Q).

    A separate, unit-testable guard called from :meth:`SimConfig.validate`. Two
    load-time fail-loud checks (no silent clamp):

    1. **RRK dof ``s >= 1`` (when set).** ``cfg.evap_rrk_dof`` is an *optional*
       effective-scalar override of the per-``n`` mode count; ``None`` is the production
       path (no check). When set it must satisfy ``s >= 1`` -- ``s < 1`` makes the RRK
       bracket exponent ``s - 1 < 0`` diverge the rate as ``E_int -> D_0`` (the
       divergent-rate regime). The nu *value* prior and the parameter-free ``Sigma(n)``
       gate carry no bound here (nu is Sourced/pinned, the gate is Derived); the nu
       *sign* guard lives in :func:`check_biphasic_config` plus the
       ``physics.evaporation.rrk_rate`` module defense (Phase-B post-review fix
       2026-07-02).

    2. **Gate-onset override provenance refuse.** ``cfg.gate_onset_override_eV`` forces a
       fixed self-bound threshold in place of the parameter-free ``Sigma(n)``. It is a
       diagnostic / falsification lever (R10), **not** a production knob, so a non-``None``
       value is **refused** unless ``cfg.allow_gate_onset_override`` is explicitly ``True``
       -- and even then the refusal only downgrades to a ``RuntimeWarning`` (mirroring the
       §6.5.1 ``allow_unvalidated_binding_pairing`` refuse->warn arm), so a forced gate can
       never *silently* enter a production / Tier-2-lock run.

    3. **Shed-convention enum typo guard (OQ-J arm, 2026-07-17).**
       ``cfg.evaporation_shed_convention`` selects the evaporation channel's
       velocity/momentum convention (``cold`` = delivered momentum-conserving
       reset, byte-inert default; ``co_moving`` = the physical co-moving He).
       An unknown value is refused here (earliest), mirroring the other enum
       surfaces (``_reject_unknown_enum``).

    Raises
    ------
    ValueError
        On ``evap_rrk_dof < 1`` when set, a non-``None`` ``gate_onset_override_eV``
        without ``allow_gate_onset_override``, or an unknown
        ``evaporation_shed_convention``.
    """
    _reject_unknown_enum(
        cfg.evaporation_shed_convention,
        ("cold", "co_moving"),
        field="evaporation_shed_convention",
    )

    if cfg.evap_rrk_dof is not None and not (cfg.evap_rrk_dof >= 1.0):
        raise ValueError(
            f"evap_rrk_dof (s) must be >= 1 when set (s < 1 diverges the RRK rate as "
            f"E_int -> D_0); got {cfg.evap_rrk_dof!r}. Leave it None for the per-n "
            f"effective_dof(n) default."
        )

    if cfg.gate_onset_override_eV is not None:
        msg = (
            f"gate_onset_override_eV={cfg.gate_onset_override_eV!r} forces a fixed "
            f"self-bound gate threshold in place of the parameter-free Sigma(n); it is a "
            f"diagnostic / falsification lever (R10), not a production knob"
        )
        if cfg.allow_gate_onset_override:
            warnings.warn(msg, RuntimeWarning)
        else:
            raise ValueError(
                msg + " (set allow_gate_onset_override=True to override)"
            )


# ---------------------------------------------------------------------------
# Tier-2 Slice T5 initial-shell-model config-load guard (plan §I.11)
# ---------------------------------------------------------------------------
_KNOWN_INITIAL_SHELL_MODELS = ("full", "density_tied")


def check_initial_shell_config(cfg: "SimConfig") -> None:
    """Validate the initial-shell dressing surface of ``cfg`` at config-load (T5).

    Rules
    -----
    1. ``initial_shell_model`` must be a known selector (typo guard).
    2. ``density_tied`` is a biphasic-seed law: the dressing is read only by
       the biphasic column-0 seed (per-ion ``n_0`` / mass / E_pot binding
       fold), so under any other ``mass_scenario`` the advertised arm would
       be silently inert — refused loudly instead (the T7/F3
       no-silent-inert convention).

    ``density_tied`` under center-pinned births (``single_initial_position=
    True``) is deliberately **legal**: the arm is then *physically* inert
    (``rho_hat(center) ≈ 1`` → ``n_0 = 21`` exactly — the H.3b wiring
    oracle), and unlike the T7 ``uniform_volume`` case there is no RNG
    draw-stream shift, so nothing is silently broken.

    Raises
    ------
    ValueError
        On an unknown ``initial_shell_model``, or ``density_tied`` outside
        ``mass_scenario='biphasic'``.
    """
    _reject_unknown_enum(
        cfg.initial_shell_model,
        _KNOWN_INITIAL_SHELL_MODELS,
        field="initial_shell_model",
    )
    if (
        cfg.initial_shell_model == "density_tied"
        and cfg.mass_scenario != "biphasic"
    ):
        raise ValueError(
            "initial_shell_model='density_tied' requires mass_scenario="
            "'biphasic': the dressing law is read only by the biphasic "
            "column-0 seed (per-ion n_0 / mass / E_pot binding fold), so "
            f"under mass_scenario={cfg.mass_scenario!r} it would be silently "
            "inert (keep the 'full' default there)"
        )


# ---------------------------------------------------------------------------
# Tier-2 Slice T6 E_int(0)-dressing p-law config-load guard (plan §I.11)
# ---------------------------------------------------------------------------
_KNOWN_INTERNAL_ENERGY_PARTITION_LAWS = ("constant", "sigma_proportional")


def check_internal_energy_partition_config(cfg: "SimConfig") -> None:
    """Validate the E_int(0)-dressing coupling surface of ``cfg`` at config-load (T6).

    Rules
    -----
    1. ``internal_energy_partition_law`` must be a known selector (typo guard).
    2. ``sigma_proportional`` couples the S2 onset to the t0 shell, and that
       onset is deposited only by the biphasic column-0 seed, so under any other
       ``mass_scenario`` the advertised arm would be silently inert -- refused
       loudly instead (the T5/T7 no-silent-inert convention). ``constant`` (the
       byte-inert default) is legal everywhere.

    Raises
    ------
    ValueError
        On an unknown ``internal_energy_partition_law``, or
        ``sigma_proportional`` outside ``mass_scenario='biphasic'``.
    """
    _reject_unknown_enum(
        cfg.internal_energy_partition_law,
        _KNOWN_INTERNAL_ENERGY_PARTITION_LAWS,
        field="internal_energy_partition_law",
    )
    if (
        cfg.internal_energy_partition_law == "sigma_proportional"
        and cfg.mass_scenario != "biphasic"
    ):
        raise ValueError(
            "internal_energy_partition_law='sigma_proportional' requires "
            "mass_scenario='biphasic': the p-law couples the S2 onset E_int(0) "
            "to the t0 shell, which is deposited only by the biphasic column-0 "
            f"seed, so under mass_scenario={cfg.mass_scenario!r} it would be "
            "silently inert (keep the 'constant' default there)"
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
# Biphasic generative-driver config-load guard (Tier-2 Phase-C Slice G)
# ---------------------------------------------------------------------------
def check_biphasic_config(cfg: "SimConfig") -> None:
    """Validate the ``biphasic`` generative-driver surface at config-load (Slice G).

    Fires **only** when ``cfg.mass_scenario == "biphasic"`` (the Tier-2 production
    scenario); a no-op for ``fixed`` / ``anchored_discrete``. Four checks:

    1. **Required partition fractions (fail-loud).** ``internal_energy_partition_
       fraction`` (f_int, the S2 onset ``E_int(0)=f_int*E_avail``) and
       ``internal_energy_retained_fraction`` (f_ret, the S1 pickup heat
       ``+f_ret*D_0``) are ``None``-defaulted (their committed values are Phase-F
       calibration outputs). Under ``biphasic`` the generative driver *cannot*
       compute the onset or S1 heat without them, so a ``None`` is refused here at
       the earliest point (CLAUDE.md principle 4), not deep in the loop. Their
       ``[0, 1]`` range is still enforced by
       :func:`check_internal_energy_budget_config`; this guard only enforces
       *presence*.
    2. **Required drag bundle (fail-loud; Slice-G review fix).** ``biphasic`` is a
       *drag-path* scenario: ``run_ion_propagation`` dispatches on
       ``drag_coefficients is not None``, so a biphasic config without a bundle
       would silently fall onto the hard-sphere collision path -- where nothing
       evolves ``E_int`` or re-applies the E_pot binding fold and the seeded
       column-0 physics corrupts the ledger. Refused here (earliest) and again at
       the driver (``ion._check_scope_ion_driver``) for non-validated configs.
    3. **Non-negative channel rates (fail-loud; Slice-G review fix + Phase-B
       post-review fix 2026-07-02).** A negative ``pickup_rate_coefficient``
       (lambda_0) gives ``P_attach < 0`` and a negative
       ``evap_rate_prefactor_per_ps`` (nu) gives ``P_shed < 0``, so the channel
       silently never fires (draws live in ``[0, 1)``) -- unphysical and
       indistinguishable from "channel on" without this refusal. The Slice-P/Q
       decisions that lambda_0 / nu carry no load-time bound cover the *value
       priors* (lambda_0 pinned at Phase F, nu Sourced 2.42), not the sign;
       ``physics.pickup.lambda_attach`` and ``physics.evaporation.rrk_rate`` carry
       the mirroring module-level defenses.
    4. **Rate-zero inert warnings (advisory, not refuses).** ``pickup_rate_coefficient
       == 0.0`` (the default) -> the Poisson pickup channel is structurally inert
       (``P_attach = 0``). That is a *legitimate* biphasic run -- the
       evaporation-only cascade from the seeded ``n_0`` shell (close cousin of the
       Phase-E relaxation stage) -- so it warns rather than raises (Slice-G decision
       #5, user 2026-07-01). ``evap_rate_prefactor_per_ps == 0.0`` mirrors it
       symmetrically (Slice-G re-review fix 2026-07-02): the evaporation channel is
       structurally inert (``k = 0``, ``P_shed = 0``), the run is pickup-only growth
       -- and since nu defaults to the Sourced 2.42 an explicit 0 is a deliberate
       diagnostic, warned so it can never pass silently. The genuinely-undefined
       knobs (f_int/f_ret) fail loud; the merely-inert ones are advisory.

    Raises
    ------
    ValueError
        When ``mass_scenario == "biphasic"`` and ``internal_energy_partition_
        fraction`` or ``internal_energy_retained_fraction`` is ``None``, or
        ``drag_coefficients`` is ``None``, or ``pickup_rate_coefficient < 0``, or
        ``evap_rate_prefactor_per_ps < 0``.
    """
    if cfg.mass_scenario != "biphasic":
        return

    for name in (
        "internal_energy_partition_fraction",
        "internal_energy_retained_fraction",
    ):
        if getattr(cfg, name) is None:
            raise ValueError(
                f"mass_scenario='biphasic' requires {name} to be set (not None): "
                "the generative driver needs both partition fractions (f_int for the "
                "S2 onset, f_ret for the S1 pickup heat) to compute the E_int budget. "
                "Their committed values are a Phase-F calibration output; set them "
                "explicitly to run a biphasic simulation."
            )

    if cfg.drag_coefficients is None:
        raise ValueError(
            "mass_scenario='biphasic' requires a drag_coefficients bundle: the "
            "generative biphasic mechanism runs on the BAOAB drag path "
            "(run_ion_propagation dispatches on drag_coefficients is not None), "
            "and without a bundle the run would silently fall onto the hard-sphere "
            "collision path, where E_int is never evolved and the E_pot binding "
            "fold is lost after column 0."
        )

    if cfg.pickup_rate_coefficient < 0.0:
        raise ValueError(
            "pickup_rate_coefficient (lambda_0) must be >= 0: a negative rate "
            "gives P_attach = 1 - exp(+|lambda_0|*dt) < 0, so the pickup channel "
            f"silently never fires; got {cfg.pickup_rate_coefficient!r}."
        )

    if cfg.evap_rate_prefactor_per_ps < 0.0:
        raise ValueError(
            "evap_rate_prefactor_per_ps (nu) must be >= 0: a negative RRK "
            "prefactor gives k < 0 -> P_shed < 0, so the evaporation channel "
            "silently never fires (the lambda_0 sign class; Phase-B post-review "
            f"fix 2026-07-02); got {cfg.evap_rate_prefactor_per_ps!r}."
        )

    if cfg.pickup_rate_coefficient == 0.0:
        warnings.warn(
            "mass_scenario='biphasic' with pickup_rate_coefficient=0.0: the Poisson "
            "He-pickup channel is structurally inert (P_attach=0), so this run is an "
            "evaporation-only cascade from the seeded n_0 shell. Set "
            "pickup_rate_coefficient>0 to enable pickup.",
            UserWarning,
        )

    if cfg.evap_rate_prefactor_per_ps == 0.0:
        warnings.warn(
            "mass_scenario='biphasic' with evap_rate_prefactor_per_ps=0.0: the "
            "energy-gated RRK evaporation channel is structurally inert (k=0, "
            "P_shed=0), so this run is pickup-only growth from the seeded n_0 "
            "shell. nu defaults to the Sourced NU_EVAP_PER_PS=2.42; an explicit 0 "
            "is a diagnostic configuration.",
            UserWarning,
        )


# Step-budget ceiling for the relaxation stage under a cooling spatial gate that can
# void the freeze early-exit (see check_relaxation_config). ~1e7 steps is a heavy but
# completable run; the 8.53e6 ps experimental cap at dt_ion=0.01 would be ~8.5e8 steps
# (a multi-hour runaway with no guaranteed freeze). Raise deliberately if a gated run
# genuinely needs a longer *bounded* relaxation.
_DENSITY_SCALED_MAX_RELAX_STEPS = 10_000_000


def check_relaxation_config(cfg: "SimConfig") -> None:
    """Validate the Tier-2 Phase-E relaxation-stage surface (Slice E2).

    Fires **only** when ``cfg.relaxation_stage_enabled`` is True (the
    ``check_biphasic_config`` pattern); a no-op otherwise, so the default config
    (stage disabled) is unaffected and the default simulation scope is unchanged.
    When enabled it requires:

    1. ``mass_scenario == "biphasic"`` -- the relaxation stage composes the
       biphasic ``biphasic_step`` under a lambda_0 = 0 view; no other scenario
       carries the ``E_int`` reservoir / channel machinery it propagates.
    2. ``relaxation_time_ps`` set and > 0 -- there is **no sourced** experimental
       flight time (MASS doc §R5 gives only "hundreds of ps"); the concrete value
       is a Phase-F campaign choice, so it is required-when-enabled with no default
       (the freeze early-exit bounds the cost of generous values).
    3. ``nu * dt_relax <= 0.1`` (``dt_relax`` = ``relaxation_dt_ps`` or ``dt_ion``)
       -- since ``k <= nu`` this bounds the one-event-per-step bias at
       ``(k*dt)^2/2 <~ 0.5%`` (MASS authorizes no larger step). ``dt_relax`` must
       be positive.
    4. ``relaxation_forces in {"coulomb", "free_flight"}`` -- the translation-arm
       enum reject guard.
    5. Under ``cooling_spatial_gate == "density_scaled"`` the freeze early-exit that
       makes check 2's "generous values are cheap" true is **void** (an ejected
       fragment stops cooling and may never freeze), so ``relaxation_time_ps /
       dt_relax`` is bounded to ``_DENSITY_SCALED_MAX_RELAX_STEPS`` -- otherwise the
       loop would silently run to the full (e.g. experimental) cap.

    Raises
    ------
    ValueError
        On any failed check above (fail-loud; CLAUDE.md principle 4).
    """
    if not cfg.relaxation_stage_enabled:
        return

    if cfg.mass_scenario != "biphasic":
        raise ValueError(
            "relaxation_stage_enabled=True requires mass_scenario='biphasic': the "
            "post-ejection relaxation stage propagates the biphasic E_int reservoir "
            "and pickup/evaporation channels (under lambda_0=0), which no other "
            f"mass scenario carries; got mass_scenario={cfg.mass_scenario!r}."
        )

    if cfg.relaxation_time_ps is None or cfg.relaxation_time_ps <= 0.0:
        raise ValueError(
            "relaxation_stage_enabled=True requires relaxation_time_ps to be set "
            "and > 0: MASS §R5 gives no sourced experimental flight time (only "
            "'hundreds of ps'), so the value is a required Phase-F campaign choice; "
            f"got relaxation_time_ps={cfg.relaxation_time_ps!r}."
        )

    dt_relax = cfg.dt_ion if cfg.relaxation_dt_ps is None else cfg.relaxation_dt_ps
    if dt_relax <= 0.0:
        raise ValueError(
            f"relaxation_dt_ps must be > 0 (defaults to dt_ion); got {dt_relax!r}."
        )
    nu_dt = cfg.evap_rate_prefactor_per_ps * dt_relax
    if nu_dt > 0.1:
        raise ValueError(
            f"nu*dt_relax must be <= 0.1 (one-event-per-step bias bound "
            f"(k*dt)^2/2 <~ 0.5%, since k <= nu); got nu={cfg.evap_rate_prefactor_per_ps} "
            f"* dt_relax={dt_relax} = {nu_dt}. Reduce relaxation_dt_ps."
        )

    # Under cooling_spatial_gate="density_scaled" the K2 freeze early-exit is NOT
    # guaranteed: a fragment ejected into vacuum (rho_He -> 0) stops cooling, so E_int
    # is held and can sit above D_0(n) indefinitely (evaporation's RRK rate -> 0 near
    # threshold). The generous-relaxation_time_ps rationale (the freeze early-exit
    # bounds the cost of large values) is therefore void, so bound the step budget
    # explicitly here to turn a silent multi-hour runaway into a load-time error.
    if cfg.cooling_spatial_gate == "density_scaled":
        n_relax_steps = cfg.relaxation_time_ps / dt_relax
        if n_relax_steps > _DENSITY_SCALED_MAX_RELAX_STEPS:
            raise ValueError(
                "cooling_spatial_gate='density_scaled' voids the relaxation freeze "
                "early-exit (an ejected fragment stops cooling and may never freeze), "
                f"so relaxation_time_ps/dt_relax = {cfg.relaxation_time_ps}/{dt_relax} "
                f"= {n_relax_steps:.3g} steps must not exceed "
                f"{_DENSITY_SCALED_MAX_RELAX_STEPS:.3g}. Set a finite, modest "
                "relaxation_time_ps (the total-strip probe uses 1000 ps) instead of "
                "the experimental flight cap."
            )

    if cfg.relaxation_forces not in ("coulomb", "free_flight"):
        raise ValueError(
            "relaxation_forces must be 'coulomb' (default) or 'free_flight'; got "
            f"{cfg.relaxation_forces!r}."
        )

    # 6. Dissipation arm (§I.11.2 item 2, arm (c)): typo-reject, then enforce the
    #    no-silent-inert pairing -- only the coulomb translation has a drag O-step,
    #    so landau_gated_drag under free_flight would be read by nothing (the
    #    T5/T7 sigma_proportional precedent: refuse rather than silently ignore).
    _reject_unknown_enum(
        cfg.relaxation_dissipation, _KNOWN_RELAXATION_DISSIPATIONS,
        field="relaxation_dissipation",
    )
    if (cfg.relaxation_dissipation == "landau_gated_drag"
            and cfg.relaxation_forces != "coulomb"):
        raise ValueError(
            "relaxation_dissipation='landau_gated_drag' requires "
            "relaxation_forces='coulomb': the free_flight arm is ballistic (no "
            "drag O-step), so the dissipation coefficient would be silently "
            f"inert; got relaxation_forces={cfg.relaxation_forces!r}."
        )


_KNOWN_RELAXATION_DISSIPATIONS = ("zero_gamma", "landau_gated_drag")


_KNOWN_DETECTION_RETAINED_POLICIES = ("refuse", "exclude")


def check_detection_config(cfg: "SimConfig") -> None:
    """Validate the Tier-2 Slice-DS detection-stage surface at config-load.

    Fires **only** when ``cfg.detection_stage_enabled`` is True (the
    ``check_relaxation_config`` pattern); a no-op otherwise, so the default
    config (stage disabled) is unaffected and existing ``cfg.json`` files that
    predate the fields load with the defaults (the Wave-7 back-compat
    criterion, probe plan Addendum C.2). When enabled it requires:

    1. ``mass_scenario == "biphasic"`` -- the detection stage continues the
       biphasic evaporation cascade (``E_int`` reservoir + RRK channel), which
       no other scenario carries.
    2. ``detection_time_ps`` set and > 0 -- **Sourced calibration data with no
       baked-in default** (8.53e6 ps TOF flight time, CALIBRATION_MAP row 24);
       the constant is data, not code, so it must be supplied explicitly.
    3. ``detection_time_ps`` beyond the earliest possible handover on the
       absolute time axis: ``ion_simulation_time``, plus ``relaxation_time_ps``
       when the relaxation stage is enabled (the stage seeds from the seed
       checkpoint's *final* column, so a detection time inside the seed window
       can never be reached forward). The stage re-checks against the actual
       handover time ``t_h`` at entry (the config cannot know realized times).
    4. ``evap_rate_prefactor_per_ps > 0`` -- with ``nu == 0`` the RRK rate is
       identically zero, so every in-band ion would sit at ``k = 0`` *without*
       being in a permanent state: the frozen/suppressed/time-exhausted
       taxonomy (design §2.3) becomes unsound and the detector read is
       meaningless. The biphasic guard merely *warns* at ``nu == 0`` (a legal
       pickup-only diagnostic run); the detection stage refuses it.

    The relaxation stage is **not** required (design §1 item 5, the skip
    path): with ``relaxation_stage_enabled=False`` the stage seeds directly
    from ``ion.npz`` and the P1-P3 handover guard is the sole defense.

    Raises
    ------
    ValueError
        On any failed check above (fail-loud; CLAUDE.md principle 4).
    """
    _reject_unknown_enum(
        cfg.detection_droplet_retained_policy,
        _KNOWN_DETECTION_RETAINED_POLICIES,
        field="detection_droplet_retained_policy",
    )
    if not cfg.detection_stage_enabled:
        return

    if cfg.mass_scenario != "biphasic":
        raise ValueError(
            "detection_stage_enabled=True requires mass_scenario='biphasic': the "
            "detection stage continues the biphasic E_int/RRK evaporation "
            "cascade, which no other mass scenario carries; got "
            f"mass_scenario={cfg.mass_scenario!r}."
        )

    if cfg.detection_time_ps is None or cfg.detection_time_ps <= 0.0:
        raise ValueError(
            "detection_stage_enabled=True requires detection_time_ps to be set "
            "and > 0: the detector arrival time is Sourced calibration data "
            "(8.53e6 ps TOF flight time, CALIBRATION_MAP row 24) with no "
            f"baked-in default; got detection_time_ps={cfg.detection_time_ps!r}."
        )

    earliest_handover_ps = float(cfg.ion_simulation_time)
    if cfg.relaxation_stage_enabled and cfg.relaxation_time_ps is not None:
        earliest_handover_ps += float(cfg.relaxation_time_ps)
    if cfg.detection_time_ps <= earliest_handover_ps:
        raise ValueError(
            f"detection_time_ps={cfg.detection_time_ps!r} must lie beyond the "
            f"seed window end on the absolute time axis (>= {earliest_handover_ps} "
            "ps = ion_simulation_time"
            + (" + relaxation_time_ps" if cfg.relaxation_stage_enabled else "")
            + "): the stage continues the cascade forward from the handover "
            "state and cannot read a time inside the seed window."
        )

    if not (cfg.evap_rate_prefactor_per_ps > 0.0):
        raise ValueError(
            "detection_stage_enabled=True requires evap_rate_prefactor_per_ps "
            "(nu) > 0: with nu = 0 the RRK rate is identically zero, so in-band "
            "ions sit at k = 0 without being frozen or suppressed -- the "
            "detection stage's permanent-state taxonomy (frozen / suppressed / "
            "time_exhausted) is unsound and the detector read meaningless; got "
            f"{cfg.evap_rate_prefactor_per_ps!r}."
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
       ``capped_cubic`` arm (Tier-2 §I.10 Slice T1) requires ``b > 0``,
       ``v_c > 0``, and ``p_tail`` in the Step-1c-adjudicated set ``{0, -1}``.
       The reserved forms' branches are written for completeness but
       unreachable -- ``physics/drag.py`` raises ``NotImplementedError`` for
       them upstream.
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
    _reject_unknown_enum(cfg.drag_form, _KNOWN_DRAG_FORMS, field="drag_form")

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
    elif form == CAPPED_CUBIC:
        # Tier-2 Addendum I §I.10 Slice T1: locked pure cubic in-band, tail
        # gamma = b*v_c^2*(v/v_c)^p_tail above the cap. gamma >= 0 everywhere
        # with no turnover for b > 0, so dissipativity reduces to b > 0 plus a
        # positive cap speed.
        b = float(c["b"])
        v_c = float(c["v_c"])
        p_tail = float(c["p_tail"])
        if not (b > 0.0):
            raise ValueError(
                f"capped_cubic drag requires b > 0 (the locked in-band pure "
                f"cubic is dissipative only for b > 0), got b={b!r}"
            )
        if not (v_c > 0.0):
            raise ValueError(
                f"capped_cubic drag requires v_c > 0 (the cap is a speed; "
                f"v_c = inf is the pure-cubic byte-identity limit), got "
                f"v_c={v_c!r}"
            )
        if p_tail not in (0.0, -1.0):
            raise ValueError(
                f"capped_cubic drag requires p_tail in {{0, -1}} (the "
                f"Step-1c-surviving tail set, TIER2_STAIRCASE_PROBE_PLAN "
                f"§I.10; any other exponent needs a fresh adjudication), "
                f"got p_tail={p_tail!r}"
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
