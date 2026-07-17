"""Detection-time continuation stage (Tier-2 Slice DS -- the I11 resolution).

Continues the post-ejection evaporation cascade from the seed checkpoint's
final state to the detector arrival time ``cfg.detection_time_ps`` (Sourced:
8.53e6 ps TOF flight, CALIBRATION_MAP row 24) with an **event-driven
(Gillespie) solver**: exponential waiting times per shed, ≤ ~n_h events per
ion, *any* detection time reachable. Design: ``TIER2_DETECTION_STAGE_DESIGN.md``.

Exactness (P1-P3, design §2.1)
------------------------------
The stage is an **exact solver of the delivered mechanism** -- not new physics
-- valid iff at the handover state (the seed's final column, time ``t_h``):

* **P1** pickup off / **P2** drag off: the stage never calls those channels,
  but omission is not proof of absence -- both are density-gated, so the
  guard verifies they are physically dead **by position** for *every* ion
  (frozen included: a pickup re-heats ``E_int`` by ``f_ret*D_0`` per capture
  and can unfreeze the cascade): ``rho * max(lambda0, 1/tau) * (t_detect -
  t_h) <= EPS_DRAIN`` with ``rho`` the local erfc ``rho_He/rho_bulk``;
* **P3** cooling negligible over the whole remaining flight, for every ion
  **not at the energetic floor** -- live in-band AND suppressed ions alike,
  because suppression is *cooling-reversible* (K2 drains ``E_int`` below
  ``Sigma(n)`` and re-opens the gate): ``rho_cool * (t_detect - t_h)/tau <=
  EPS_DRAIN`` where ``rho_cool`` is the local density under
  ``cooling_spatial_gate == "density_scaled"`` and **1.0 under ``"none"``**
  (ungated cooling acts everywhere, so a non-frozen ungated ion always
  fails -- ungated runs must freeze in E2 and arrive here as per-ion
  no-ops). Energetically frozen ions are exempt from this bound only
  (cooling cannot unfreeze; the exposure bound above still applies to them).

Under the guard-verified P1-P3, ``E_int`` is constant between sheds, so the
locked mechanism reduces *exactly* to a Markov jump chain at the delivered
rate (:func:`~i2_helium_md.physics.evaporation.rrk_rate` verbatim, all
gating encoded). Between fires the gate status cannot change; across a fire
the self-bound margin ``G = E_int - Sigma(n)`` is shed-invariant, so a
suppressed or frozen ion is suppressed/frozen **forever** (permanent states,
§2.3) -- a claim that is true precisely because the guard has verified the
cooling and exposure channels are dead.

Event loop (per ion, exact; design §2.2)
----------------------------------------
``k = rrk_rate(E_int, n)``; ``k == 0`` -> permanent state (reason recorded);
else draw ``dt = -ln(u)/k``; past ``t_detect`` -> ``time_exhausted``; else
fire: ``n -> n-1``, the **cold-shed velocity/mass/KE-defect reset**
(:func:`~i2_helium_md.physics.mass_jump.cold_shed`: ``v -> m/(m-m_He) * v`` --
NB the design doc's early "v unchanged" phrasing was corrected as-built; the
kick never feeds back into the (E_int, n) jump chain, which reads neither
``v`` nor positions), the K1 drain ``E_int -= D_0(n_pre)``
(:func:`~i2_helium_md.physics.internal_energy_budget.dE_int_shed_eV`), and
the ``e_bind_pair`` fold ``E_pot += D_0(n_pre)``. Termination is
unconditional: ``n`` strictly decreases per fire.

There is **no timestep**, hence no ``nu*dt <= 0.1`` guard: exponential
waiting times have zero discretization bias. The E2 per-step Bernoulli chain
``P = 1 - exp(-k*dt)`` at constant ``k`` converges in distribution to this
exact process as ``dt -> 0`` -- the fixed-dt path is the biased approximation
of *this* stage, not the reverse.

Ledger (design §1 fork 4)
-------------------------
Every fire books its K1 drain, e_bind fold, and cold-shed mass-transfer
defect per event, so the 5-term invariant ``E_kin + E_pot + E_dissip +
E_mass_transfer + E_int`` closes across the stage boundary: per fire
``dE_int = -D_0`` cancels the fold ``dE_pot = +D_0``, and the cold-shed
``dE_kin`` cancels the booked ``dE_mass_transfer`` exactly (the reduced-mass
form). ``E_dissip`` never moves (P3: zero cooling drain). ``E_kin`` at
detection is recomputed from the terminal ``(m, v)`` by the existing idiom.

Seeding (design §1 item 5, the skip path)
-----------------------------------------
``seed_ckpt`` is ``relaxation.npz`` when E2 ran, or ``ion.npz`` directly when
the run skips E2 (``relaxation_stage_enabled=False``) -- both are v7
``IonCheckpoint`` s, one loading contract. The E2 seed-coherence guards apply
identically (biphasic scenario; the ``mass_history_kg[:, -1] ==
mass_final_kg`` stride guard); the P1-P3 guard is the sole physics defense on
the skip path.

RNG contract
------------
A fresh stage-private stream via ``SeedSequence((cfg.seed,
DETECTION_STREAM_KEY))`` (the ``RELAXATION_STREAM_KEY`` pattern); no existing
stream is touched. Draws are per-ion sequential in ion-major order, one
scalar uniform per waiting time; the draw count is data-dependent, which is
acceptable because the stream is stage-private (no cross-stage draw-order
contract exists).

Artifact
--------
``detection.npz`` -- a small dedicated file (NOT an ``IonCheckpoint``; schema
v7 untouched): per-ion terminal arrays + the ragged per-ion event records as
flat arrays with ``(2N+1,)`` offsets. ``allow_pickle=False``; shape-validated
on load.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Optional

import numpy as np

from ..config import SimConfig, check_detection_config
from ..physics.constants import MASS_HE_AMU, MASS_I_ION_AMU, U
from ..physics.dissociation_ladder import d0_of_n, resolve_ladder
from ..physics.evaporation import gate_margin_eV, rrk_rate
from ..physics.helium_density import rho_he_ratio
from ..physics.internal_energy_budget import dE_int_shed_eV
from ..physics.mass_jump import cold_shed
from ..physics.potentials import droplet_potential
from .checkpoint import check_biphasic_seed_checkpoint, stage_stream_rng
from .ion import drag_gate_steepness
from .ion_propagation_step import (
    _amu_ang2_ps2_to_eV,
    _depth,
    _E_kin_eV,
    ion_state_from_checkpoint_column,
)

# Fixed module-level key that spawns the detection stage's own PCG64 stream via
# SeedSequence((cfg.seed, DETECTION_STREAM_KEY)) -- the RELAXATION_STREAM_KEY
# pattern; the ion-stage and relaxation-stage streams are never touched.
DETECTION_STREAM_KEY: int = 0xD5_2026

#: P3 handover-guard bound on the neglected K2 drain fraction
#: ``rho_factor * (t_detect - t_h)/tau`` (design §2.1) -- a *defensive bound*,
#: not a knob: gated ejected ions have erfc-suppressed rho/rho_bulk that
#: underflows to 0, so the guard passes with enormous margin or fails loudly
#: on a genuinely wrong input (a live in-bubble or ungated ion). Dimensionless;
#: equals the fractional E_int error committed by treating cooling as zero.
EPS_DRAIN: float = 1e-6

#: Legal per-ion permanent/terminal state reasons (design §2.3).
STATE_REASONS: tuple[str, ...] = ("frozen", "suppressed", "time_exhausted")

# detection.npz schema (this artifact's own version counter; it is NOT an
# IonCheckpoint and does not participate in the v5->v6->v7 cascade).
_DETECTION_SCHEMA_VERSION: int = 1

# Non-array scalar fields of the artifact (unwrapped from 0-d arrays on load).
_SCALAR_FIELDS: tuple[str, ...] = (
    "num_molecules", "t_handover_ps", "detection_time_ps",
)


@dataclass(frozen=True)
class DetectionResult:
    """Result of one detection-stage run (the ``detection.npz`` payload).

    Per-ion arrays have shape ``(2N,)`` (both iodine fragments of every
    molecule are independent I+He_n detections, the E1 convention). The ragged
    per-ion event records are stored as flat arrays with ``event_offsets``:
    ion ``i``'s events occupy ``slice(event_offsets[i], event_offsets[i+1])``.

    Attributes
    ----------
    num_molecules : int
        ``N`` (the ensemble carries ``2N`` ions).
    t_handover_ps : float
        Absolute handover time ``t_h`` (the seed checkpoint's final column
        time; continues the ion-stage/E2 axis).
    detection_time_ps : float
        Absolute detector arrival time ``t_detect`` (Sourced).
    n_detected : np.ndarray, shape (2N,)
        Terminal He occupancy at ``t_detect`` (int-valued float, the v7
        convention).
    E_int_detected_eV, mass_detected_kg : np.ndarray, shape (2N,)
        Terminal internal energy [eV] and complex mass [kg].
    vx_detected, vy_detected, vz_detected : np.ndarray, shape (2N,)
        Terminal velocity components [A/ps] (cold-shed kicked per fire,
        constant between fires; positions are never propagated -- design §2.4).
    E_kin_detected_eV, E_pot_detected_eV : np.ndarray, shape (2N,)
        Terminal kinetic energy (recomputed from ``(m, v)`` by the existing
        idiom) and potential energy (held MD potential + the ``e_bind_pair``
        fold at the terminal ``n`` -- the free-flight convention) [eV].
    E_dissip_detected_eV, E_mass_transfer_detected_eV : np.ndarray, shape (2N,)
        Terminal cumulative ledger channels [eV]: ``E_dissip`` is carried
        through **unchanged** (P3: zero cooling drain); ``E_mass_transfer``
        accumulates the per-fire cold-shed defects. Together with the three
        arrays above these make the 5-term closure checkable across the
        stage boundary without the seed checkpoint.
    state_reason : np.ndarray, shape (2N,), unicode
        Per-ion terminal reason: ``"frozen"`` (energetic floor -- converged
        terminal), ``"suppressed"`` (ejected net-self-unbound complex, rides
        to the detector at its handover ``n``; OQ-B untouched), or
        ``"time_exhausted"`` (cascade live at the detector -- the physical
        in-flight snapshot).
    event_offsets : np.ndarray, shape (2N+1,), int
        Ragged-record offsets into the flat event arrays; ``event_offsets[0]
        == 0`` and ``event_offsets[-1] == num_events``.
    event_time_ps : np.ndarray, shape (num_events,)
        Absolute shed times [ps] (enables the free detection-time
        sensitivity re-read, design §3.3 -- valid only *up to* the generated
        ``t_detect``; events beyond it were never sampled).
    event_pre_shed_n : np.ndarray, shape (num_events,)
        Pre-shed rung ``n`` of each fire.
    event_dE_int_eV : np.ndarray, shape (num_events,)
        K1 drain per fire (``-D_0(n_pre)`` < 0) [eV].
    event_dE_bind_fold_eV : np.ndarray, shape (num_events,)
        ``e_bind_pair`` E_pot fold per fire (``+D_0(n_pre)`` > 0) [eV].
    event_dE_mass_transfer_eV : np.ndarray, shape (num_events,)
        Cold-shed KE defect per fire (< 0) [eV].
    """

    num_molecules: int
    t_handover_ps: float
    detection_time_ps: float
    n_detected: np.ndarray
    E_int_detected_eV: np.ndarray
    mass_detected_kg: np.ndarray
    vx_detected: np.ndarray
    vy_detected: np.ndarray
    vz_detected: np.ndarray
    E_kin_detected_eV: np.ndarray
    E_pot_detected_eV: np.ndarray
    E_dissip_detected_eV: np.ndarray
    E_mass_transfer_detected_eV: np.ndarray
    state_reason: np.ndarray
    event_offsets: np.ndarray
    event_time_ps: np.ndarray
    event_pre_shed_n: np.ndarray
    event_dE_int_eV: np.ndarray
    event_dE_bind_fold_eV: np.ndarray
    event_dE_mass_transfer_eV: np.ndarray

    @property
    def terminal_n(self) -> np.ndarray:
        """The detected terminal shell count -- the E1 input hook.

        ``compute_terminal_shell_distribution`` duck-types on a
        ``terminal_n``-bearing input; the presence of :attr:`state_reason`
        disambiguates the ``"detected"`` provenance from E2's ``"relaxed"``.
        """
        return self.n_detected

    def events_for_ion(self, ion_id: int) -> slice:
        """The flat-array slice holding ion ``ion_id``'s event records."""
        return slice(
            int(self.event_offsets[ion_id]), int(self.event_offsets[ion_id + 1])
        )

    def reason_fractions(self) -> dict[str, float]:
        """Per-reason ensemble fractions (the report's state-reason columns)."""
        reasons = np.asarray(self.state_reason)
        return {
            reason: float(np.mean(reasons == reason)) for reason in STATE_REASONS
        }


def _permanent_reason(E_int_eV: float, n: int, *, picture: str, kappa: float,
                      gate_onset_eV, ladder=None) -> str:
    """Classify a ``k == 0`` ion into its permanent state (design §2.3).

    Called only when :func:`rrk_rate` returned exactly 0, so with ``nu > 0``
    (enforced at config-load) the lanes are: ``n <= 0`` or ``E_int`` at/below
    the top rung -> ``"frozen"`` (the energetic floor); ``n >= 2`` at/above
    the self-bound threshold -> ``"suppressed"``; and the **float-underflow
    band** -- ``n >= 2`` strictly in-band but with ``E_int`` within ~1e-6
    relative of ``D_0(n)``, where the RRK bracket ``(1 - D_0/E_int)^(s-1)``
    underflows to exactly 0 (s-1 up to 59) -> ``"frozen"`` (review fix
    2026-07-07: the expected waiting time exceeds any physical flight time
    by hundreds of orders of magnitude, so this IS the energetic floor;
    E2's strict freeze mask ``E_int < D_0`` hands such ions over un-frozen,
    making the lane reachable -- the earlier defensive raise crashed the
    whole run on it). Anything else contradicts the delivered gating and
    raises.
    """
    if n <= 0:
        return "frozen"
    d0 = float(d0_of_n(n, picture=picture, kappa=kappa, ladder=ladder))
    if E_int_eV <= d0:
        return "frozen"
    if n >= 2 and float(gate_margin_eV(
        E_int_eV, n, picture=picture, kappa=kappa, gate_onset_eV=gate_onset_eV,
        ladder=ladder,
    )) >= 0.0:
        return "suppressed"
    if n >= 2 and E_int_eV > d0:
        # In-band with k == 0: the RRK-bracket underflow band just above the
        # top rung -- physically frozen (see docstring).
        return "frozen"
    raise RuntimeError(
        f"detection stage: rrk_rate returned 0 for an unclassifiable state "
        f"(n={n}, E_int={E_int_eV!r} eV) -- inconsistent with the delivered "
        "gating; this indicates a rate/classification drift bug."
    )


def run_detection_stage(
    seed_ckpt,
    cfg: SimConfig,
    *,
    rng: Optional[np.random.Generator] = None,
    save_path: Optional[str | Path] = None,
) -> DetectionResult:
    """Continue the evaporation cascade to ``cfg.detection_time_ps``. Exact.

    Parameters
    ----------
    seed_ckpt : IonCheckpoint
        A finished ``biphasic`` v7 checkpoint: ``relaxation.npz`` when E2 ran,
        or ``ion.npz`` on the skip path (design §1 item 5). The stage seeds
        from its **final** stored column.
    cfg : SimConfig
        The run config. Must have ``detection_stage_enabled=True`` and pass
        :func:`~i2_helium_md.config.check_detection_config`. The rate surface
        (``nu``, picture, kappa, ``evap_rrk_dof``, ``gate_onset_override_eV``)
        is read verbatim -- the stage adds no rate physics.
    rng : np.random.Generator, optional
        The stage's generator. Default ``None`` -> a fresh stream seeded
        ``SeedSequence((cfg.seed, DETECTION_STREAM_KEY))`` (or
        non-reproducible when ``cfg.seed is None``).
    save_path : str or Path, optional
        If given, the result is written there via :func:`save_detection_result`
        (the pipeline passes ``<run_dir>/detection.npz``).

    Returns
    -------
    DetectionResult

    Raises
    ------
    ValueError
        If the stage is disabled, any ``check_detection_config`` arm fails,
        the seed checkpoint is not ``biphasic``, the seed's final stored
        column is stale (stride guard), ``detection_time_ps`` does not lie
        beyond the realized handover time ``t_h``, or the P1-P3 handover
        guard fails (the per-ion violator list, with the remedy).
    NotImplementedError
        If the config selects the declared-but-unbuilt ``tabulated``
        density-profile arm (point-of-use refusal, the biphasic_step
        precedent). The tabulated *ladder* is wired at Slice T2 (§I.10):
        ``resolve_ladder`` injects it; selecting it without a rung table is
        a loud ``ValueError`` instead.
    """
    if not cfg.detection_stage_enabled:
        raise ValueError(
            "run_detection_stage requires cfg.detection_stage_enabled=True "
            "(the stage is opt-in; enable it and set detection_time_ps)."
        )
    check_detection_config(cfg)   # fires the enabled-path checks (fail-loud)

    # Ladder resolution (Slice T2, §I.10): rrk_rate / d0_of_n /
    # dE_int_shed_eV below thread the resolved ladder -- None rides Form-U
    # (byte-identical), a TabulatedLadder overrides it; 'tabulated' without a
    # table fails loudly in resolve_ladder (the pre-T2 lazy refusal moved to
    # the data path). The density-profile arm stays a point-of-use refusal
    # (the biphasic_step precedent -- review fix 2026-07-07): selecting it
    # would otherwise silently run erfc physics.
    ladder = resolve_ladder(cfg.dissociation_ladder, cfg.tabulated_ladder_rungs_eV)
    if cfg.helium_density_profile != "erf_complement":
        raise NotImplementedError(
            f"helium_density_profile={cfg.helium_density_profile!r} is a "
            "declared-but-unbuilt rule-2 arm; the P1-P3 handover guard "
            "composes only the analytic 'erf_complement' density "
            "(physics/helium_density.py)."
        )

    # Seed-coherence guards (shared with E2 -- one source for both
    # continuation stages), applied to whichever checkpoint the pipeline
    # hands over: relaxation.npz or, on the skip path, ion.npz.
    check_biphasic_seed_checkpoint(seed_ckpt, stage="run_detection_stage")

    seed = ion_state_from_checkpoint_column(seed_ckpt, -1)
    t_h = float(seed.time_ps)
    t_detect = float(cfg.detection_time_ps)
    if t_detect <= t_h:
        raise ValueError(
            f"detection_time_ps={t_detect} must lie beyond the realized "
            f"handover time t_h={t_h} ps on the absolute axis (the seed "
            "checkpoint's final column time); the config-load bound uses the "
            "nominal window ends and cannot see realized times."
        )

    num_molecules = int(seed_ckpt.num_molecules)
    two_n = 2 * num_molecules
    picture = cfg.ladder_electronic_picture
    kappa = float(cfg.ladder_steepness)
    nu = float(cfg.evap_rate_prefactor_per_ps)
    s_override = cfg.evap_rrk_dof
    gate_onset = cfg.gate_onset_override_eV
    tau_ps = float(cfg.internal_energy_cooling_tau_ps)

    n0 = np.rint(np.asarray(seed.n_shell, dtype=float)).astype(int)
    E_int0 = np.asarray(seed.E_int_eV, dtype=float)

    # ---- P1-P3 handover guard (design §2.1; review-hardened 2026-07-07) --
    # Two per-ion dimensionless bounds over the remaining flight
    # dt_flight = t_detect - t_h:
    #
    # * **Cooling drain (P3)** -- applies to every ion NOT at the energetic
    #   floor, i.e. live in-band AND suppressed ions alike: suppression is
    #   **cooling-reversible** (K2 drains E_int below Sigma(n) and re-opens
    #   the gate -- biphasic_step runs cooling before the gate for exactly
    #   this reason), so "suppressed forever" holds only where cooling is
    #   dead. rho_cool mirrors biphasic_step's cool_rho: the local erfc
    #   density under density_scaled, exactly 1.0 under "none" (ungated
    #   cooling acts everywhere -- a non-frozen ungated ion can never hand
    #   over; it must freeze in E2).
    # * **Helium exposure (P1/P2)** -- applies to EVERY ion, frozen
    #   included: pickup (lambda ~ lambda0*rho) adds He and re-heats E_int by
    #   f_ret*D_0 per capture (which can unfreeze the cascade), and drag
    #   (gamma ~ rho) moves v. Both channels are merely *omitted* here, so
    #   the guard must verify they are physically dead by position:
    #   rho * max(lambda0, 1/tau) * dt_flight <= EPS_DRAIN. The erfc density
    #   underflows to exactly 0 for properly ejected ions, so this passes
    #   with enormous margin or fails loudly on an in-/near-bubble ion --
    #   the skip path's actual defense.
    d0_h = np.asarray(
        d0_of_n(np.maximum(n0, 1), picture=picture, kappa=kappa, ladder=ladder),
        dtype=float,
    )
    frozen_h = (n0 <= 0) | (E_int0 < d0_h)
    depth = _depth(seed.x, seed.y, seed.z,
                   np.asarray(seed_ckpt.droplet_radii_angstrom, dtype=float))
    rho_density = np.asarray(
        rho_he_ratio(depth, steepness=drag_gate_steepness(cfg)), dtype=float
    )
    rho_cool = (
        rho_density if cfg.cooling_spatial_gate == "density_scaled"
        else np.ones(two_n, dtype=float)
    )
    dt_flight = t_detect - t_h
    cooling_fraction = rho_cool * dt_flight / tau_ps
    exposure_rate = max(float(cfg.pickup_rate_coefficient), 1.0 / tau_ps)
    exposure_fraction = rho_density * exposure_rate * dt_flight
    cooling_violation = (~frozen_h) & (cooling_fraction > EPS_DRAIN)
    exposure_violation = exposure_fraction > EPS_DRAIN
    violators = cooling_violation | exposure_violation
    # T9 leg A' (V0-2 convention, 2026-07-16): under the "exclude" policy an
    # energetically *bound* violator (KE + U < E_bind: turned around by the
    # droplet well -- it can never decouple at any relaxation length) is
    # classified `droplet_retained` and excluded from the free-flight event
    # loop instead of failing the run. Unbound violators (slow escapers that
    # a longer relaxation WOULD decouple) still refuse loudly, under either
    # policy -- the trapped class is physics, an undecoupled escaper is a
    # configuration error.
    droplet_retained = np.zeros(two_n, dtype=bool)
    if np.any(violators) and cfg.detection_droplet_retained_policy == "exclude":
        bound = _conservatively_bound(seed, seed_ckpt, cfg)
        droplet_retained = violators & bound
        violators = violators & ~bound
    if np.any(violators):
        idx = np.flatnonzero(violators)
        badness = np.maximum(
            np.where(frozen_h, 0.0, cooling_fraction), exposure_fraction
        )
        worst = idx[np.argmax(badness[idx])]
        raise ValueError(
            "detection-stage P1-P3 handover guard failed: "
            f"{idx.size}/{two_n} ions are not decoupled from helium over the "
            f"remaining flight ({int(np.sum(cooling_violation))} with a "
            "non-negligible neglected K2 cooling drain [non-frozen ions: "
            "suppression is cooling-reversible], "
            f"{int(np.sum(exposure_violation))} with live pickup/drag "
            f"exposure; bounds vs EPS_DRAIN = {EPS_DRAIN:g}). Worst ion "
            f"{int(worst)}: rho={rho_density[worst]:.3e}, "
            f"cooling_fraction={cooling_fraction[worst]:.3e}, "
            f"exposure_fraction={exposure_fraction[worst]:.3e}. Violating "
            f"ion ids: {idx.tolist()}. Remedy: extend relaxation_time_ps (or "
            "enable the relaxation stage on the skip path) so ions decouple "
            "before handover; an ungated (cooling_spatial_gate='none') run "
            "must freeze in E2 -- its non-frozen ions can never hand over. "
            "A genuinely well-trapped (energetically bound) ion never "
            "decouples: select detection_droplet_retained_policy='exclude' "
            "to classify it droplet_retained per the V0-2 convention."
        )

    # ---- RNG (stage-private stream; shared derivation helper) -----------
    if rng is None:
        rng = stage_stream_rng(cfg.seed, DETECTION_STREAM_KEY)

    # ---- Per-ion event loop (ion-major order; design §2.2) ---------------
    n_out = n0.astype(float).copy()
    E_int_out = E_int0.copy()
    m_amu_out = np.asarray(seed.mass_kg, dtype=float) / U
    vx_out = np.asarray(seed.vx, dtype=float).copy()
    vy_out = np.asarray(seed.vy, dtype=float).copy()
    vz_out = np.asarray(seed.vz, dtype=float).copy()
    E_pot_out = np.asarray(seed.E_pot_eV, dtype=float).copy()
    E_mt_out = np.asarray(seed.E_mass_transfer_eV, dtype=float).copy()
    E_dissip_out = np.asarray(seed.E_dissip_eV, dtype=float).copy()  # unchanged
    reasons: list[str] = []

    ev_offsets = np.zeros(two_n + 1, dtype=int)
    ev_time: list[float] = []
    ev_n_pre: list[int] = []
    ev_dE_int: list[float] = []
    ev_fold: list[float] = []
    ev_dE_mt: list[float] = []

    for i in range(two_n):
        if droplet_retained[i]:
            # V0-2 convention: droplet-retained ions carry their handover
            # state verbatim (no events) and are excluded from the IHe_n
            # read by their state reason.
            reasons.append("droplet_retained")
            ev_offsets[i + 1] = len(ev_time)
            continue
        t = t_h
        n_i = int(n_out[i])
        E_i = float(E_int_out[i])
        m_i = float(m_amu_out[i])
        v_i = np.array([vx_out[i], vy_out[i], vz_out[i]], dtype=float)

        while True:
            k = float(rrk_rate(
                E_i, n_i, nu=nu, picture=picture, kappa=kappa,
                evap_rrk_dof=s_override, gate_onset_eV=gate_onset,
                ladder=ladder,
            ))
            if k == 0.0:
                reasons.append(_permanent_reason(
                    E_i, n_i, picture=picture, kappa=kappa,
                    gate_onset_eV=gate_onset, ladder=ladder,
                ))
                break
            u = float(rng.random())
            # rng.random() lives in [0, 1): u == 0 (prob ~2^-53) maps to an
            # infinite waiting time -- statistically exact (P(u=0) = 0 in the
            # continuum) and handled without a log(0) warning.
            dt_wait = math.inf if u <= 0.0 else -math.log(u) / k
            if t + dt_wait > t_detect:
                reasons.append("time_exhausted")
                break
            t = t + dt_wait

            # Fire: compose the delivered primitives verbatim (design §2.6).
            n_pre = n_i
            shed = cold_shed(v_i, m_i, m_he_amu=MASS_HE_AMU)
            dE_int = float(
                dE_int_shed_eV(n_pre, picture=picture, kappa=kappa, ladder=ladder)
            )
            # E_pot fold: exactly +D_0(n_pre) = -dE_int, so the K1/fold pair
            # cancels bit-exactly in the 5-term sum (review fix 2026-07-07:
            # recomputing the fold as an e_bind_pair difference left ~1e-17 eV
            # per-fire residues and two redundant ladder lookups).
            fold = -dE_int
            dE_mt_eV = float(_amu_ang2_ps2_to_eV(shed.dE_mass_transfer))

            ev_time.append(t)
            ev_n_pre.append(n_pre)
            ev_dE_int.append(dE_int)
            ev_fold.append(fold)
            ev_dE_mt.append(dE_mt_eV)

            n_i = n_pre - 1
            m_i = float(shed.m_plus_amu)
            v_i = np.asarray(shed.v_plus, dtype=float)
            E_i = E_i + dE_int
            E_pot_out[i] = E_pot_out[i] + fold
            E_mt_out[i] = E_mt_out[i] + dE_mt_eV

        ev_offsets[i + 1] = len(ev_time)
        n_out[i] = float(n_i)
        E_int_out[i] = E_i
        m_amu_out[i] = m_i
        vx_out[i], vy_out[i], vz_out[i] = v_i

    # m <-> n lockstep (the biphasic_step Q3 guard, cheap structural catch).
    if not np.allclose(
        m_amu_out, MASS_I_ION_AMU + n_out * MASS_HE_AMU, rtol=0.0, atol=1e-6,
    ):
        raise AssertionError(
            "detection stage m<->n consistency drift: mass and n disagree "
            "beyond 1e-6 amu after the event loop (cold_shed is the single "
            "mass source; the counters must stay in lockstep)."
        )

    mass_out_kg = m_amu_out * U
    result = DetectionResult(
        num_molecules=num_molecules,
        t_handover_ps=t_h,
        detection_time_ps=t_detect,
        n_detected=n_out,
        E_int_detected_eV=E_int_out,
        mass_detected_kg=mass_out_kg,
        vx_detected=vx_out,
        vy_detected=vy_out,
        vz_detected=vz_out,
        E_kin_detected_eV=_E_kin_eV(mass_out_kg, vx_out, vy_out, vz_out),
        E_pot_detected_eV=E_pot_out,
        E_dissip_detected_eV=E_dissip_out,
        E_mass_transfer_detected_eV=E_mt_out,
        state_reason=np.asarray(reasons),
        event_offsets=ev_offsets,
        event_time_ps=np.asarray(ev_time, dtype=float),
        event_pre_shed_n=np.asarray(ev_n_pre, dtype=float),
        event_dE_int_eV=np.asarray(ev_dE_int, dtype=float),
        event_dE_bind_fold_eV=np.asarray(ev_fold, dtype=float),
        event_dE_mass_transfer_eV=np.asarray(ev_dE_mt, dtype=float),
    )

    if save_path is not None:
        save_detection_result(result, save_path)

    return result


def _conservatively_bound(seed, seed_ckpt, cfg: SimConfig) -> np.ndarray:
    """Per-ion mask: cannot reach infinity under the E2 conservative dynamics.

    The droplet-retained criterion for the ``exclude`` policy (V0-2
    convention, barrier-corrected 2026-07-16). The E2 relaxation "coulomb"
    mode is **zero-gamma** (conservative BAOAB: droplet well + residual
    Coulomb, no drag), so an ion escapes iff its total energy clears the
    **effective potential** along the outward path:

        E_tot = E_kin + U(depth)   must reach every
        V_eff(r') = U(r' - R) + L^2 / (2 m r'^2)   for r' >= r,

    with L = m |r x v| the (conserved) angular momentum about the droplet
    center. The radial-only criterion (E_kin + U < E_bind) is the L = 0
    special case; a near-threshold ion on a mostly-tangential orbit can sit
    *above* the radial threshold yet *below* its centrifugal barrier — a
    quasi-bound resonance that never decouples (the leg-A' apc3 "ion 9":
    +1.1 meV radial margin, ~5 A of outward drift in 3000 ps). The residual
    pair Coulomb (repulsive, partner >> 100 A away) only aids escape, so
    ignoring it keeps the criterion conservative in the correct direction:
    no genuinely escaping ion is ever classified retained... the reverse —
    Coulomb could push a marginally-retained ion out — is bounded by the
    partner distance and is negligible at the guard's operating point
    (< 0.1 meV beyond ~10^4 A).

    Returns a boolean (2N,) mask. Units: masses kg -> amu; energies eV via
    the amu*A^2/ps^2 conversion (strict dimensional bookkeeping).
    """
    x = np.asarray(seed.x, dtype=float)
    y = np.asarray(seed.y, dtype=float)
    z = np.asarray(seed.z, dtype=float)
    vx = np.asarray(seed.vx, dtype=float)
    vy = np.asarray(seed.vy, dtype=float)
    vz = np.asarray(seed.vz, dtype=float)
    m_amu = np.asarray(seed.mass_kg, dtype=float) / U
    radii = np.asarray(seed_ckpt.droplet_radii_angstrom, dtype=float)

    r = np.sqrt(x**2 + y**2 + z**2)
    r_safe = np.maximum(r, 1e-12)
    E_kin = _E_kin_eV(
        np.asarray(seed.mass_kg, dtype=float), vx, vy, vz
    )
    U_here = droplet_potential(
        r - radii,
        steepness=cfg.potential_steepness,
        binding_energy=cfg.binding_energy_I_ion_eV,
    )
    E_tot = E_kin + U_here                                   # eV

    # |r x v|^2 per ion [A^4/ps^2]; L^2/(2 m r'^2) = m*h2/(2 r'^2) in
    # amu*A^2/ps^2 -> eV.
    hx = y * vz - z * vy
    hy = z * vx - x * vz
    hz = x * vy - y * vx
    h2 = hx**2 + hy**2 + hz**2

    # Outward-path barrier on a grid r' in [r, r + 12*steepness] (U is
    # asymptotic well beyond ~7 widths; the centrifugal term only falls).
    n_grid = 400
    span = 12.0 * float(cfg.potential_steepness)
    frac = np.linspace(0.0, 1.0, n_grid)[None, :]            # (1, G)
    r_grid = r_safe[:, None] + span * frac                   # (2N, G)
    U_grid = droplet_potential(
        r_grid - radii[:, None],
        steepness=cfg.potential_steepness,
        binding_energy=cfg.binding_energy_I_ion_eV,
    )
    E_cf_grid = _amu_ang2_ps2_to_eV(
        0.5 * m_amu[:, None] * h2[:, None] / r_grid**2
    )
    barrier = np.maximum(
        (U_grid + E_cf_grid).max(axis=1),
        cfg.binding_energy_I_ion_eV,      # the r' -> inf asymptote
    )
    return E_tot < barrier


# ===========================================================================
# detection.npz I/O (dedicated small artifact -- NOT an IonCheckpoint)
# ===========================================================================
def save_detection_result(result: DetectionResult, path: str | Path) -> Path:
    """Write a :class:`DetectionResult` to ``.npz`` (``detection.npz``).

    Scalars are wrapped as 0-d arrays; the unicode ``state_reason`` array
    round-trips under ``allow_pickle=False``. Returns the path written.
    """
    p = Path(path)
    if p.suffix != ".npz":
        p = p.with_suffix(".npz")
    p.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        f.name: np.asarray(getattr(result, f.name)) for f in fields(result)
    }
    payload["detection_schema_version"] = np.asarray(_DETECTION_SCHEMA_VERSION)
    np.savez_compressed(p, **payload)
    return p


def load_detection_result(path: str | Path) -> DetectionResult:
    """Load a :class:`DetectionResult` from ``detection.npz`` (fail-loud).

    ``allow_pickle=False``; validates the artifact's own schema version, field
    completeness, and the ragged-record shape invariants (offsets ``(2N+1,)``
    starting at 0, monotone non-decreasing, ending at the flat event length;
    all flat event arrays equal-length; per-ion arrays ``(2N,)``).
    """
    p = Path(path)
    if p.suffix != ".npz":
        p = p.with_suffix(".npz")
    if not p.exists():
        raise FileNotFoundError(f"detection artifact not found: {p}")

    with np.load(p, allow_pickle=False) as npz:
        raw = {k: npz[k] for k in npz.files}

    if "detection_schema_version" not in raw:
        raise ValueError(
            f"detection artifact at {p} has no detection_schema_version field; "
            "it is not a detection.npz file."
        )
    version = int(raw.pop("detection_schema_version"))
    if version != _DETECTION_SCHEMA_VERSION:
        raise ValueError(
            f"detection artifact at {p} has detection_schema_version="
            f"{version}, this code expects {_DETECTION_SCHEMA_VERSION}. "
            "Re-run the detection stage."
        )

    expected = {f.name for f in fields(DetectionResult)}
    missing = expected - set(raw)
    if missing:
        raise ValueError(
            f"detection artifact at {p} is missing fields: {sorted(missing)}"
        )
    unknown = set(raw) - expected
    if unknown:
        raise ValueError(
            f"detection artifact at {p} has unknown fields: {sorted(unknown)}. "
            "This file was probably produced by a different version of the code."
        )

    kwargs: dict = {}
    for f in fields(DetectionResult):
        arr = raw[f.name]
        if f.name == "num_molecules":
            kwargs[f.name] = int(arr)
        elif f.name in _SCALAR_FIELDS:
            kwargs[f.name] = float(arr)
        else:
            kwargs[f.name] = np.asarray(arr)

    two_n = 2 * kwargs["num_molecules"]
    per_ion = (
        "n_detected", "E_int_detected_eV", "mass_detected_kg",
        "vx_detected", "vy_detected", "vz_detected",
        "E_kin_detected_eV", "E_pot_detected_eV", "E_dissip_detected_eV",
        "E_mass_transfer_detected_eV", "state_reason",
    )
    for name in per_ion:
        if kwargs[name].shape != (two_n,):
            raise ValueError(
                f"detection artifact field {name!r} has shape "
                f"{kwargs[name].shape}, expected ({two_n},)"
            )
    offsets = kwargs["event_offsets"]
    if offsets.shape != (two_n + 1,):
        raise ValueError(
            f"detection artifact event_offsets has shape {offsets.shape}, "
            f"expected ({two_n + 1},)"
        )
    num_events = int(offsets[-1])
    if int(offsets[0]) != 0 or np.any(np.diff(offsets) < 0):
        raise ValueError(
            "detection artifact event_offsets must start at 0 and be monotone "
            "non-decreasing."
        )
    for name in ("event_time_ps", "event_pre_shed_n", "event_dE_int_eV",
                 "event_dE_bind_fold_eV", "event_dE_mass_transfer_eV"):
        if kwargs[name].shape != (num_events,):
            raise ValueError(
                f"detection artifact field {name!r} has shape "
                f"{kwargs[name].shape}, expected ({num_events},) per "
                "event_offsets[-1]."
            )

    bad_reasons = set(np.asarray(kwargs["state_reason"]).tolist()) - set(
        STATE_REASONS
    )
    if bad_reasons:
        raise ValueError(
            f"detection artifact state_reason holds unknown reasons "
            f"{sorted(bad_reasons)}; expected members of {STATE_REASONS}."
        )

    return DetectionResult(**kwargs)
