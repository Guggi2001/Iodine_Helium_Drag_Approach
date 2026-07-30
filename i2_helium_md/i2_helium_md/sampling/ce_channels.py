"""CE channel-mixture sampler (Tier-2 atlas §3.5l (C) design — (B) source).

Per-molecule Coulomb-explosion channel draw for the (C) probe
(``TIER2_CE_CHANNEL_EXIT_STRIP_DESIGN.md`` §3.1). Each I₂ pair fires one of
a small set of CE channels, independently of the droplet geometry::

    c_m ~ Categorical(w_single, w_Q2, w_Q3),      Σ w_c = 1
    E_m ~ TruncNormal(E_c, σ_c;  E_m > 0)         [eV per I⁺ fragment]

with channel means (per I⁺, from the point-charge reference at
R₀ = 2.666 Å times the shared Hatherly fraction f)::

    E_single = cfg.ce_single_ker_eV            (not Coulombic — no f)
    E_Q2     = f · E_REF_PER_ION_EV            (= 2.16 eV at f = 0.8)
    E_Q3     = f · 2 · E_REF_PER_ION_EV        (= 4.32 eV at f = 0.8)

**Wiring (the emulation route, OQ-A):** the pair Coulomb drive is scaled per
molecule with the dimensionless ``s_m = E_m / E_REF_PER_ION_EV`` at unit
charges — kinematically exact for the scored I⁺ (the force enters only
through the ``q₁q₂`` product). In Q3-labeled molecules one fragment is
tagged ``CE_CHANNEL_Q3_PARTNER`` (the emulated I²⁺) and excluded from every
I⁺-scored observable (OQ-I, non-optional in v1); its propagation is kept —
the pair kinematics need it.

RNG contract (design §5): the draw lives on a **dedicated named stream**
seeded ``SeedSequence((cfg.seed, CE_CHANNEL_STREAM_KEY))`` — appended after
all existing streams; channels-off consumes nothing, so every pre-existing
stream is byte-identical (the s(n) S1–S4 off-mode precedent). The frozen
draw order is (1) one uniform per molecule for the channel, (2) one standard
normal per molecule for ``E_m`` (rejection-resampled per-molecule until
``E_m > 0``), (3) one uniform per molecule for the Q3 partner pick.

Units: energies eV; ``s_m``, weights, f dimensionless.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Dedicated CE channel-draw stream key (``SeedSequence((seed, key))`` via
#: ``simulation.checkpoint.stage_stream_rng``); distinct from the relaxation
#: (0xE2_2026) and detection (0xD5_2026) stage keys.
CE_CHANNEL_STREAM_KEY: int = 0xCE1_2026

#: Point-charge per-I⁺ Coulomb reference E_ref = k_C / (2·R₀) [eV] at
#: R₀ = 2.666 Å (design §3.1). The Coulomb constant literal matches the
#: single force-path source in ``physics/interactions.py`` (14.39964548
#: eV·Å = e²/(4πε₀)); the per-ion share is half the pair energy.
CE_R0_ANGSTROM: float = 2.666
E_REF_PER_ION_EV: float = 14.39964548 / (2.0 * CE_R0_ANGSTROM)

#: Per-ion channel codes stored in the checkpoint v8 ``ce_channel`` field.
#: ``CE_CHANNEL_NONE`` marks a run without channel sampling (the v7→v8
#: migration sentinel and the ``ce_channel_mode='off'`` fill value).
CE_CHANNEL_NONE: int = -1
CE_CHANNEL_SINGLE: int = 0
CE_CHANNEL_Q2: int = 1
CE_CHANNEL_Q3: int = 2
CE_CHANNEL_Q3_PARTNER: int = 3

#: Scored-channel labels by code (the partner shares the Q3 physics but is
#: excluded from I⁺ scoring).
CE_CHANNEL_LABELS: dict[int, str] = {
    CE_CHANNEL_NONE: "none",
    CE_CHANNEL_SINGLE: "single",
    CE_CHANNEL_Q2: "q2",
    CE_CHANNEL_Q3: "q3",
    CE_CHANNEL_Q3_PARTNER: "q3_partner",
}


@dataclass(frozen=True)
class CeChannelDraw:
    """One ensemble's CE channel assignment (2N per-ion layout).

    Attributes
    ----------
    channel : np.ndarray, shape (2N,), int
        Per-ion channel code (``CE_CHANNEL_*``). Both fragments of a
        molecule share the molecule's channel; in a Q3 molecule exactly one
        fragment carries ``CE_CHANNEL_Q3_PARTNER`` (uniform pick).
    E_m_eV : np.ndarray, shape (2N,), float
        Per-ion sampled channel KER stamp [eV per I⁺] (> 0), identical for
        the two fragments of a molecule.
    pair_scale : np.ndarray, shape (N,), float
        Per-molecule dimensionless Coulomb-drive scale
        ``s_m = E_m / E_REF_PER_ION_EV`` (> 0). Multiplies the pair Coulomb
        on top of the global ``cfg.E_coulomb_scale``.
    """

    channel: np.ndarray
    E_m_eV: np.ndarray
    pair_scale: np.ndarray


def ce_channel_means_eV(
    *, fraction_f: float, single_ker_eV: float,
) -> tuple[float, float, float]:
    """Per-I⁺ channel means ``(E_single, E_Q2, E_Q3)`` [eV] (design §3.1).

    ``E_Q2 = f·E_ref`` and ``E_Q3 = 2f·E_ref`` from the shared
    channel-independent fraction f (Hatherly); the single channel's mean is
    the Bounded ``ce_single_ker_eV`` (not Coulombic — no f).
    """
    return (
        float(single_ker_eV),
        float(fraction_f) * E_REF_PER_ION_EV,
        2.0 * float(fraction_f) * E_REF_PER_ION_EV,
    )


def sample_ce_channels(
    num_molecules: int,
    rng: np.random.Generator,
    *,
    weights: tuple[float, float, float],
    fraction_f: float,
    sigma_eV: tuple[float, float, float],
    single_ker_eV: float,
) -> CeChannelDraw:
    """Draw the per-molecule CE channel mixture (frozen three-draw order).

    Parameters
    ----------
    num_molecules : int
        N molecules (2N ions in the standard layout).
    rng : np.random.Generator
        The **dedicated** CE channel stream (``CE_CHANNEL_STREAM_KEY``);
        never the ion-stage stream.
    weights : (w_single, w_Q2, w_Q3)
        Categorical channel weights, non-negative, summing to 1.
    fraction_f : float
        Shared Coulomb fraction f ∈ (0, 1] (prior 0.80 [0.65, 0.90]).
    sigma_eV : (σ_single, σ_Q2, σ_Q3)
        Per-channel TruncNormal widths [eV], ≥ 0 (P2-measured 0.42/0.31/0.55).
    single_ker_eV : float
        Single-channel mean KER [eV per I⁺] (> 0; P2-measured 0.53).

    Returns
    -------
    CeChannelDraw
        Channel codes, per-ion E_m stamps, per-molecule pair scales.

    Raises
    ------
    ValueError
        On invalid weights / fraction / widths / mean (mirrors the
        config-load guard so a bypassing caller still fails loudly).
    """
    w = np.asarray(weights, dtype=float)
    if w.shape != (3,) or np.any(w < 0) or not np.isclose(w.sum(), 1.0):
        raise ValueError(
            f"ce channel weights must be 3 non-negative values summing to 1, "
            f"got {weights!r}"
        )
    if not (0.0 < fraction_f <= 1.0):
        raise ValueError(f"ce fraction f must be in (0, 1], got {fraction_f!r}")
    sig = np.asarray(sigma_eV, dtype=float)
    if sig.shape != (3,) or np.any(~np.isfinite(sig)) or np.any(sig < 0):
        raise ValueError(
            f"ce channel sigmas must be 3 finite values >= 0, got {sigma_eV!r}"
        )
    if not (single_ker_eV > 0):
        raise ValueError(
            f"ce single-channel KER must be > 0 eV, got {single_ker_eV!r}"
        )

    means = np.asarray(
        ce_channel_means_eV(fraction_f=fraction_f, single_ker_eV=single_ker_eV)
    )

    # Draw 1: channel categorical (one uniform per molecule).
    u_channel = rng.random(num_molecules)
    edges = np.cumsum(w)
    channel_mol = np.searchsorted(edges, u_channel, side="right")
    channel_mol = np.minimum(channel_mol, 2)          # float-edge safety

    # Draw 2: E_m ~ TruncNormal(E_c, sigma_c; E_m > 0) via per-molecule
    # rejection resampling (data-dependent draw count — acceptable on a
    # dedicated stream; the detection-stage precedent).
    z = rng.standard_normal(num_molecules)
    E_mol = means[channel_mol] + sig[channel_mol] * z
    bad = E_mol <= 0.0
    while np.any(bad):
        z_new = rng.standard_normal(int(np.count_nonzero(bad)))
        E_mol[bad] = means[channel_mol[bad]] + sig[channel_mol[bad]] * z_new
        bad = E_mol <= 0.0

    # Draw 3: Q3 partner pick (one uniform per molecule, consumed for every
    # molecule so the draw count is composition-independent).
    u_partner = rng.random(num_molecules)

    # Assemble the (2N,) per-ion layout: [atom-1 block | atom-2 block].
    channel = np.concatenate([channel_mol, channel_mol]).astype(int)
    is_q3 = channel_mol == CE_CHANNEL_Q3
    partner_is_atom1 = is_q3 & (u_partner < 0.5)
    partner_is_atom2 = is_q3 & ~(u_partner < 0.5)
    channel[:num_molecules][partner_is_atom1] = CE_CHANNEL_Q3_PARTNER
    channel[num_molecules:][partner_is_atom2] = CE_CHANNEL_Q3_PARTNER

    E_m = np.concatenate([E_mol, E_mol])
    return CeChannelDraw(
        channel=channel,
        E_m_eV=E_m,
        pair_scale=E_mol / E_REF_PER_ION_EV,
    )


def ce_pair_scale_from_checkpoint(ckpt) -> np.ndarray | None:
    """Per-molecule Coulomb pair scale from a v8 ``IonCheckpoint``.

    Returns the ``(N,)`` array ``s_m = ce_E_m_eV / E_REF_PER_ION_EV`` when
    the run carries a sampled channel assignment, or ``None`` when every ion
    is ``CE_CHANNEL_NONE`` (channels off / migrated pre-v8 file) — the
    ``None`` return keeps the off path structurally bit-identical (no scale
    array is ever built or applied).

    Raises
    ------
    ValueError
        On a mixed assignment (some ions labeled, some not) or a labeled ion
        with a non-finite / non-positive ``ce_E_m_eV`` stamp — both indicate
        a corrupted checkpoint.
    """
    channel = np.asarray(ckpt.ce_channel, dtype=int)
    if np.all(channel == CE_CHANNEL_NONE):
        return None
    if np.any(channel == CE_CHANNEL_NONE):
        raise ValueError(
            "checkpoint carries a mixed ce_channel assignment (some ions "
            "labeled, some CE_CHANNEL_NONE) — the channel draw is "
            "per-molecule and all-or-nothing."
        )
    E_m = np.asarray(ckpt.ce_E_m_eV, dtype=float)
    if np.any(~np.isfinite(E_m)) or np.any(E_m <= 0.0):
        raise ValueError(
            "checkpoint ce_E_m_eV must be finite and > 0 on a "
            "channel-sampled run."
        )
    n_mol = channel.size // 2
    if not np.array_equal(E_m[:n_mol], E_m[n_mol:]):
        raise ValueError(
            "checkpoint ce_E_m_eV differs between the two fragments of a "
            "molecule — the per-molecule stamp must be shared."
        )
    return E_m[:n_mol] / E_REF_PER_ION_EV


def ce_scored_ion_mask(ce_channel) -> np.ndarray:
    """Boolean mask of ions that enter I⁺-scored observables (OQ-I).

    ``True`` for every ion except the ``CE_CHANNEL_Q3_PARTNER`` tags — the
    emulated I²⁺ partners are excluded from the histogram, KED and fate
    reads (they are detected on a different m/q row in the experiment).
    All-``True`` for a channels-off run (no partner tags exist).
    """
    return np.asarray(ce_channel, dtype=int) != CE_CHANNEL_Q3_PARTNER
