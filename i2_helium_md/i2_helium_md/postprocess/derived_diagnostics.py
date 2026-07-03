"""Derived diagnostics that *explain* a finished ``biphasic`` run (Tier-2 Phase E).

Post-hoc reconstructions read from a finished ``biphasic`` run's v7 ion
checkpoint. The three thin helpers were delivered at Phase D (Slice Z, the
generative-vs-anchored bridge report, `TIER2_PHASE_D_IMPLEMENTATION_PLAN.md`
§2.2/§3); Phase-E Slice E5 composes them into the full regime-determination
diagnostic :func:`reconstruct_diagnostics` (module renamed ``bridge_diagnostics``
→ ``derived_diagnostics`` at E5a for its concern, Quality Principle 6).

Thin helpers (delivered Phase D):

* :func:`mean_shell_count` — the emergent mean ``n(t)`` (average of the
  ``n_shell`` array over the 2N ion rows), overlaid on the anchored
  21→19→14 staircase family.
* :func:`crossing_time_ps` — the per-ion gate-crossing time
  ``t× = min{t : E_int(t) < Σ(n(t))}`` [ps] (MASS §6.11; literally GAH25's
  ``t₀``), with ``Σ(n)`` from the Phase-A ladder cumulative. NaN where the
  gate never opens in-window.
* :func:`regime_parameter` — the regime order parameter
  ``Π(t) = λ(n)·f_ret·τ`` (dimensionless; ``Π > 1`` shedding persists,
  ``Π < 1`` freeze), with ``ρ_He/ρ_bulk`` **re-derived** from the checkpoint
  positions + per-atom droplet radii + the same resolved gate steepness the
  run used (the §2.2 reconstruction contract).

Composed diagnostic (Phase E, Slice E5b):

* :func:`reconstruct_diagnostics` — assembles ``t×`` + an ensemble summary,
  ``Π(t)``, the shell-retaining ↔ total-strip regime label, total-strip
  reachability (``E_∞(0) → 0``), and advisory sanity flags into a
  :class:`Diagnostics`. Field names are pinned (the Phase-F F3/F4 scoreboard
  binds to ``t_cross_ps`` / ``Pi_t`` / ``regime_label`` /
  ``total_strip_reachable`` / ``sanity_flags``).

These are diagnostics of a *finished run* — pure reads, no RNG, no integrator
state, no mutation. Zero schema cost: everything is reconstructed from the v7
arrays.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import SimConfig
from ..physics.evaporation import is_self_bound
from ..physics.helium_density import rho_he_ratio
from ..physics.pickup import lambda_attach
from ..physics.solvation_cooling import e_infinity_eV
from ..simulation.ion import drag_gate_steepness


def mean_shell_count(n_shell) -> np.ndarray:
    """Emergent mean He-shell count over the run's ions, per stored step.

    Parameters
    ----------
    n_shell : array_like, shape (2N, T)
        Per-atom integer-valued shell-count trajectory (the v7
        ``IonCheckpoint.n_shell`` array). Both CE fragments are ions, so the
        mean runs over all 2N rows (no row selection).

    Returns
    -------
    np.ndarray, shape (T,)
        ``mean(n_shell, axis=0)`` as float [dimensionless count].

    Raises
    ------
    ValueError
        If ``n_shell`` is not 2-D.
    """
    n = np.asarray(n_shell, dtype=float)
    if n.ndim != 2:
        raise ValueError(
            f"n_shell must be 2-D (2N, T); got shape {n.shape}"
        )
    return n.mean(axis=0)


def crossing_time_ps(
    E_int_eV,
    n_shell,
    time_ps,
    *,
    picture: str,
    kappa: float,
    gate_onset_eV: float | None = None,
) -> np.ndarray:
    """Per-ion gate-crossing time ``t× = min{t : E_int(t) < Σ(n(t))}`` [ps].

    MASS §6.11: ``t×`` is the instant the internal-energy reservoir drains
    below the cumulative dissociation threshold ``Σ(n) = Σ_{i≤n} D₀(i)`` and
    the evaporation gate opens — literally GAH25's ``t₀``. Reconstructed
    post-hoc from the stored v7 arrays; the threshold follows the *stored*
    ``n(t)``, not a constant. The gate is the driver's own
    :func:`~i2_helium_md.physics.evaporation.is_self_bound` (shared surface,
    not a comparison copy): strict (``E_int == Σ`` is still self-bound), and
    resolving ``gate_onset_eV`` exactly as the run did — a
    ``gate_onset_override_eV`` run must pass its override here or the
    reconstruction would silently misrepresent it.

    In the bridge run all ions must agree (pre-crossing dynamics are fully
    deterministic: evaporation suppressed, pickup dead at ``n = n*``, pure K2
    cooling) — per-ion disagreement is a wiring flag, not physics (plan
    §2.2). The reconstruction lands within one stored dt above the continuous
    crossing (sampling discreteness).

    Parameters
    ----------
    E_int_eV : array_like, shape (2N, T)
        Internal-energy reservoir trajectory [eV] (v7 ``E_int_eV``).
    n_shell : array_like, shape (2N, T)
        Shell-count trajectory (integer-valued; v7 ``n_shell``).
    time_ps : array_like, shape (T,)
        Stored time axis [ps].
    picture : str
        Ladder electronic picture (``cfg.ladder_electronic_picture``).
    kappa : float
        Ladder steepness κ (``cfg.ladder_steepness``).
    gate_onset_eV : float, optional
        The run's ``cfg.gate_onset_override_eV``. ``None`` (default) → the
        parameter-free ladder cumulative ``Σ(n)``; a set float → that fixed
        diagnostic threshold, mirroring the driver's resolution.

    Returns
    -------
    np.ndarray, shape (2N,)
        First time [ps] with ``E_int < Σ(n)`` per ion; ``NaN`` where the gate
        never opens within the stored window (an edge guard only — the §6.11
        Lyapunov argument guarantees the gate always opens eventually).

    Raises
    ------
    ValueError
        If the array shapes are inconsistent.
    """
    E = np.asarray(E_int_eV, dtype=float)
    n = np.asarray(n_shell)
    t = np.asarray(time_ps, dtype=float)
    if E.ndim != 2 or n.shape != E.shape:
        raise ValueError(
            f"E_int_eV and n_shell must share one (2N, T) shape; got "
            f"{E.shape} and {n.shape}"
        )
    if t.ndim != 1 or t.shape[0] != E.shape[1]:
        raise ValueError(
            f"time_ps must have shape (T,) = ({E.shape[1]},); got {t.shape}"
        )
    if t.shape[0] > 1 and not np.all(np.diff(t) > 0.0):
        raise ValueError(
            "time_ps must be strictly increasing (a monotone stored time axis); "
            f"got {t}"
        )

    below = is_self_bound(
        E, n, picture=picture, kappa=kappa, gate_onset_eV=gate_onset_eV
    )                                         # (2N, T) bool: E_int < threshold
    first_idx = below.argmax(axis=1)          # 0 when a row is all-False ...
    ever_below = below.any(axis=1)            # ... so mask those rows to NaN
    return np.where(ever_below, t[first_idx], np.nan)


def regime_parameter(ckpt, cfg: SimConfig) -> np.ndarray:
    """Regime order parameter ``Π(t) = λ(n(t))·f_ret·τ`` per ion [dimensionless].

    MASS §6.11: ``Π > 1`` → pickup-sustained shedding persists (strip side);
    ``Π < 1`` → freeze side (the shell self-binds). At the 9 Å / 0.80 eV
    priored central point the expectation is the freeze side throughout, with
    ``Π = 0`` at gate-open (``n = n*`` zeroes the Langmuir cap) and ``Π → 0``
    at droplet exit (``ρ_He → 0``) — plan §2.2 corrected oracle.

    Reconstruction contract (§2.2): ``ρ_He/ρ_bulk`` is **re-derived** from the
    checkpoint positions and per-atom droplet radii through the same
    erf-complement gate and the same resolved steepness the run used
    (``simulation.ion.drag_gate_steepness`` — shared surface, not a copy),
    then fed to the Phase-B ``lambda_attach``. The stored columns sample the
    post-step positions, so this is the stored-grid diagnostic of the in-step
    rate, exact on the stored samples.

    Parameters
    ----------
    ckpt : IonCheckpoint (or any object with the fields read here)
        Needs ``positions_x/y/z`` (2N, T) [Å], ``droplet_radii_angstrom``
        (2N,) [Å], ``n_shell`` (2N, T).
    cfg : SimConfig
        The run's config: supplies λ₀, p, cap, rate form, f_ret, τ, and the
        gate-steepness resolution. ``helium_density_profile`` must be the
        built ``'erf_complement'`` arm.

    Returns
    -------
    np.ndarray, shape (2N, T)
        ``Π`` per ion per stored step (report as mean ± envelope over ions).

    Raises
    ------
    ValueError
        If ``f_ret`` is unset (a non-biphasic config cannot be
        reconstructed; ``f_ret = 0`` is in-band and gives ``Π ≡ 0``), or if
        ``τ`` is unset/non-positive (unphysical — mirrors the
        ``check_solvation_cooling_config`` load guard for hand-built
        configs that bypassed ``validate()``).
    NotImplementedError
        If ``cfg.helium_density_profile`` selects the declared-but-unbuilt
        ``'tabulated'`` arm — the erf re-derivation would silently
        misrepresent such a run (the driver's point-of-use contract).
    """
    if cfg.helium_density_profile != "erf_complement":
        raise NotImplementedError(
            f"helium_density_profile={cfg.helium_density_profile!r}: only the "
            "analytic 'erf_complement' gate can be re-derived from positions; "
            "a tabulated-profile run cannot be reconstructed this way."
        )
    f_ret = cfg.internal_energy_retained_fraction
    if f_ret is None:
        raise ValueError(
            "internal_energy_retained_fraction (f_ret) is unset; Π = "
            "λ·f_ret·τ needs the run's biphasic budget knobs."
        )
    tau_ps = cfg.internal_energy_cooling_tau_ps
    if tau_ps is None or not (tau_ps > 0.0):
        raise ValueError(
            f"internal_energy_cooling_tau_ps (tau) must be > 0 (a non-positive "
            f"relaxation time is unphysical); got {tau_ps!r}"
        )

    x = np.asarray(ckpt.positions_x, dtype=float)
    y = np.asarray(ckpt.positions_y, dtype=float)
    z = np.asarray(ckpt.positions_z, dtype=float)
    droplet_radii = np.asarray(ckpt.droplet_radii_angstrom, dtype=float)
    n_shell = np.asarray(ckpt.n_shell)

    if x.ndim != 2 or y.shape != x.shape or z.shape != x.shape or n_shell.shape != x.shape:
        raise ValueError(
            "positions_x/y/z and n_shell must share one (2N, T) shape; got "
            f"x={x.shape}, y={y.shape}, z={z.shape}, n_shell={n_shell.shape}"
        )
    if droplet_radii.ndim != 1 or droplet_radii.shape[0] != x.shape[0]:
        raise ValueError(
            f"droplet_radii_angstrom must have shape (2N,) = ({x.shape[0]},); "
            f"got {droplet_radii.shape}"
        )

    # Depth into the droplet r - R_droplet [Å] (negative inside), the same
    # geometry the driver's _depth lift computes each step.
    r = np.sqrt(x * x + y * y + z * z)
    depth = r - droplet_radii[:, np.newaxis]

    rho = rho_he_ratio(depth, steepness=drag_gate_steepness(cfg))
    lam = lambda_attach(
        rho, n_shell,
        lambda0=cfg.pickup_rate_coefficient,
        p=cfg.pickup_occupancy_exponent,
        cap=cfg.pickup_occupancy_cap,
        pickup_rate_form=cfg.pickup_rate_form,
    )
    return lam * f_ret * tau_ps


# ===========================================================================
# Composed regime-determination diagnostic (Phase-E Slice E5b)
# ===========================================================================
# MASS §6.11 sanity band on t× [ps]: GAH25 t0 ≈ 5–6.5 ps is a ±factor-2 prior,
# so a factor-10 miss (< 1 ps or > 15 ps) is the genuine flag, not the prior.
_T_CROSS_BAND_PS: tuple[float, float] = (1.0, 15.0)
# Regime convention (a documented reporting rule, NOT a physics knob): total-strip
# if the ensemble never self-binds (majority no finite t×) or is stripped to the
# core (median terminal n ≤ this).
_TOTAL_STRIP_MEDIAN_N: float = 1.0


@dataclass(frozen=True)
class TCrossSummary:
    """Ensemble summary of the per-ion gate-crossing time ``t×`` [ps].

    Attributes
    ----------
    median_ps : float
        Median of the finite per-ion ``t×`` (``NaN`` if none crossed).
    spread_ps : float
        Peak-to-peak of the finite per-ion ``t×`` (``NaN`` if none crossed).
    all_agree : bool
        True iff every ion crossed and all crossings fall within one stored dt --
        the wiring check for the bridge, whose pre-crossing dynamics are fully
        deterministic (evaporation suppressed, pickup dead at ``n = n*``, pure K2).
    """

    median_ps: float
    spread_ps: float
    all_agree: bool


@dataclass(frozen=True)
class Diagnostics:
    """Post-hoc diagnostics that explain a finished biphasic run (MASS §6.11).

    Reconstructed at **zero schema cost** from the v7 arrays. Field names are
    pinned: the Phase-F F3/F4 scoreboard binds to ``t_cross_ps`` / ``Pi_t`` /
    ``regime_label`` / ``total_strip_reachable`` / ``sanity_flags``.

    Attributes
    ----------
    t_cross_ps : np.ndarray, shape (2N,)
        Per-ion gate-crossing time [ps]; ``NaN`` where the gate never opens.
    t_cross_summary : TCrossSummary
        Ensemble median / spread / all-agree flag of ``t_cross_ps``.
    Pi_t : np.ndarray, shape (2N, T)
        Regime order parameter ``Π(t) = λ(n)·f_ret·τ`` per ion per stored step.
    regime_label : str
        ``"shell_retaining"`` (self-binds early → moderate terminal n) or
        ``"total_strip"`` (never self-binds in window, or stripped to the core).
    total_strip_reachable : bool
        Whether ``E_∞(N) → 0`` as ``N → 0`` under the run's picture/κ (OQ6; a
        Slice-K structural consistency confirmation, not an independent anchor).
    sanity_flags : tuple[str, ...]
        Advisory flags (empty when all sane): a ``t×`` out of the [1, 15] ps band,
        or ``Π > 1`` where a freeze-side point is expected -- units/wiring smells,
        not physics verdicts.
    """

    t_cross_ps: np.ndarray
    t_cross_summary: TCrossSummary
    Pi_t: np.ndarray
    regime_label: str
    total_strip_reachable: bool
    sanity_flags: tuple[str, ...]


def reconstruct_diagnostics(
    ckpt,
    cfg: SimConfig,
    *,
    freeze_side_expected: bool = True,
) -> Diagnostics:
    """Reconstruct the regime-determination diagnostics from a finished v7 run.

    Composes the three delivered helpers (:func:`crossing_time_ps`,
    :func:`regime_parameter`) with the ensemble summary, the regime label, the
    total-strip reachability, and the sanity flags. Re-derives nothing; zero
    schema cost.

    Parameters
    ----------
    ckpt : IonCheckpoint (or any object with the fields read here)
        A finished ``biphasic`` run: needs ``E_int_eV`` (2N, T), ``n_shell``
        (2N, T), ``time_ps`` (T,), and (for ``Π``) ``positions_x/y/z`` +
        ``droplet_radii_angstrom``.
    cfg : SimConfig
        The run's config (picture, κ, the biphasic budget knobs, the gate-onset
        override). Must be a reconstructable biphasic config -- ``regime_parameter``
        raises if ``f_ret``/``τ`` are unset or the density profile is not
        ``'erf_complement'``.
    freeze_side_expected
        Whether the run sits at a freeze-side condition (e.g. the 9 Å /
        0.80 eV pinned point), where ``Π > 1`` is a units/wiring smell and
        earns an advisory flag. Pass ``False`` at shedding-persists conditions
        (e.g. 2.70 eV production), where ``Π > 1`` is legitimate physics
        (MASS §6.11) -- the freeze-side reading must not be carried there.

    Returns
    -------
    Diagnostics
        The frozen bundle of ``t×`` + summary, ``Π(t)``, regime label,
        reachability, and sanity flags.

    Raises
    ------
    ValueError
        Empty ensemble (no atoms or no stored steps), or a fractional
        ``n_shell`` entry (violates the v7 int-valued contract -> upstream
        corruption; the same fail-loud convention as the E1 extractor).
    """
    picture = cfg.ladder_electronic_picture
    kappa = cfg.ladder_steepness

    E_int = np.asarray(ckpt.E_int_eV, dtype=float)
    n_shell = np.asarray(ckpt.n_shell, dtype=float)
    time_ps = np.asarray(ckpt.time_ps, dtype=float)

    if E_int.ndim != 2 or E_int.size == 0:
        raise ValueError(
            f"empty ensemble: E_int_eV has shape {E_int.shape}, nothing to "
            f"reconstruct diagnostics from."
        )
    # v7 stores n_shell as an int-valued float array; a fractional entry is
    # upstream corruption. Fail loud (E1's convention on the same field)
    # instead of letting np.rint feed the regime label rounded garbage.
    residual = n_shell - np.rint(n_shell)
    if np.any(residual != 0.0):
        raise ValueError(
            "n_shell must be integer-valued (v7 contract); found fractional "
            "entries -- upstream corruption."
        )

    t_cross = crossing_time_ps(
        E_int, n_shell, time_ps, picture=picture, kappa=kappa,
        gate_onset_eV=cfg.gate_onset_override_eV,
    )
    Pi_t = regime_parameter(ckpt, cfg)

    # Ensemble t× summary.
    finite = np.isfinite(t_cross)
    if finite.any():
        median_ps = float(np.median(t_cross[finite]))
        spread_ps = float(np.ptp(t_cross[finite]))
    else:
        median_ps = float("nan")
        spread_ps = float("nan")
    dt = float(np.min(np.diff(time_ps))) if time_ps.size > 1 else 0.0
    all_agree = bool(finite.all() and spread_ps <= dt + 1e-12)
    summary = TCrossSummary(median_ps=median_ps, spread_ps=spread_ps, all_agree=all_agree)

    # Regime label (documented reporting rule; not a config knob).
    terminal_n = np.rint(n_shell[:, -1]).astype(int)
    median_terminal_n = float(np.median(terminal_n))
    n_no_cross = int((~finite).sum())
    majority_no_cross = n_no_cross * 2 > t_cross.size
    if majority_no_cross or median_terminal_n <= _TOTAL_STRIP_MEDIAN_N:
        regime_label = "total_strip"
    else:
        regime_label = "shell_retaining"

    # Total-strip reachability: E_∞(N) → 0 as N → 0 (occupancy-resolved, OQ6).
    e_inf_0 = float(
        e_infinity_eV(0, picture=picture, kappa=kappa, s_abs_eV=cfg.solv_struct_asymptote_eV)
    )
    total_strip_reachable = bool(np.isclose(e_inf_0, 0.0, atol=1e-9))

    # Advisory sanity flags (§6.11) -- units/wiring smells, not physics verdicts.
    flags: list[str] = []
    lo, hi = _T_CROSS_BAND_PS
    if finite.any():
        tc = t_cross[finite]
        if np.any(tc < lo) or np.any(tc > hi):
            flags.append(
                f"t_cross out of the [{lo}, {hi}] ps sanity band "
                f"(min {float(tc.min()):.3g}, max {float(tc.max()):.3g} ps)"
            )
    # The Pi > 1 smell is condition-specific: it flags only where the caller
    # declared a freeze-side expectation (at a shedding-persists condition,
    # e.g. 2.70 eV production, Pi > 1 is legitimate physics, not a smell).
    pi_max = float(np.nanmax(Pi_t)) if Pi_t.size else 0.0
    if freeze_side_expected and pi_max > 1.0:
        flags.append(
            f"Pi exceeds 1 (max {pi_max:.3g}): shedding-persists regime -- a "
            "units/wiring flag where a freeze-side point is expected "
            "(e.g. 9A / 0.80 eV)"
        )

    return Diagnostics(
        t_cross_ps=t_cross,
        t_cross_summary=summary,
        Pi_t=Pi_t,
        regime_label=regime_label,
        total_strip_reachable=total_strip_reachable,
        sanity_flags=tuple(flags),
    )
