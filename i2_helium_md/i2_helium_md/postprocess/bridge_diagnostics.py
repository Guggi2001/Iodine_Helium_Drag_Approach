"""Phase-D bridge reconstruction diagnostics (Slice Z — the thin versions).

Post-hoc reconstructions read from a finished ``biphasic`` run's v7 ion
checkpoint for the generative-vs-anchored bridge report
(`TIER2_PHASE_D_IMPLEMENTATION_PLAN.md` §2.2/§3):

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

These are diagnostics of a *finished run* — pure reads, no RNG, no
integrator state, no mutation. Phase-E Slice D2 generalizes this module into
the full regime-determination diagnostic; keep additions thin and reusable.
"""

from __future__ import annotations

import numpy as np

from ..config import SimConfig
from ..physics.dissociation_ladder import ladder_cumsum
from ..physics.helium_density import rho_he_ratio
from ..physics.pickup import lambda_attach
from ..simulation.ion import _drag_gate_steepness


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
) -> np.ndarray:
    """Per-ion gate-crossing time ``t× = min{t : E_int(t) < Σ(n(t))}`` [ps].

    MASS §6.11: ``t×`` is the instant the internal-energy reservoir drains
    below the cumulative dissociation threshold ``Σ(n) = Σ_{i≤n} D₀(i)`` and
    the evaporation gate opens — literally GAH25's ``t₀``. Reconstructed
    post-hoc from the stored v7 arrays; the threshold follows the *stored*
    ``n(t)``, not a constant. The crossing is strict (``E_int == Σ`` is still
    self-bound), matching the driver's gate.

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

    sigma = ladder_cumsum(n, picture=picture, kappa=kappa)   # (2N, T) [eV]
    below = E < sigma
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
    (``simulation.ion._drag_gate_steepness`` — shared surface, not a copy),
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
        If ``f_ret`` or ``τ`` is unset/non-positive (a non-biphasic config
        cannot be reconstructed).
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

    x = np.asarray(ckpt.positions_x, dtype=float)
    y = np.asarray(ckpt.positions_y, dtype=float)
    z = np.asarray(ckpt.positions_z, dtype=float)
    droplet_radii = np.asarray(ckpt.droplet_radii_angstrom, dtype=float)
    n_shell = np.asarray(ckpt.n_shell)

    # Depth into the droplet r - R_droplet [Å] (negative inside), the same
    # geometry the driver's _depth lift computes each step.
    r = np.sqrt(x * x + y * y + z * z)
    depth = r - droplet_radii[:, np.newaxis]

    rho = rho_he_ratio(depth, steepness=_drag_gate_steepness(cfg))
    lam = lambda_attach(
        rho, n_shell,
        lambda0=cfg.pickup_rate_coefficient,
        p=cfg.pickup_occupancy_exponent,
        cap=cfg.pickup_occupancy_cap,
        pickup_rate_form=cfg.pickup_rate_form,
    )
    return lam * f_ret * tau_ps
