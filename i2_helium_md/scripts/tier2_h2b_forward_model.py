"""The 1D chord-model forward model ("the twin") - zero MD.

Recovered VERBATIM from the H.2b scratchpad driver ``h2b_feasibility.py``
(analytic feasibility pass, Addendum H design freeze D1-D7; executed
2026-07-11; both surviving scratchpad copies byte-identical) and committed
under Tier-2 plan SI.11 V0-3 (2026-07-16) so every T5-T9 oracle comparison
is reproducible in-repo. V0-3 scope additions on top of the recovered file
(everything else is the recovered code unchanged):

* repo-relative path bootstrap (was a hardcoded absolute path),
* outputs to ``data/runs/h2b_forward_model/`` (gitignored; was the
  scratchpad directory),
* a **birth-law input** (``BIRTH_LAW``): ``"uniform_volume"`` (the native
  L1 law: r^2 in droplet volume + hard margin) or ``"boltzmann"`` (the
  repo's realized sampler, ``sampling/radial_positions.py`` - center-pinned
  at T = 0.4 K), so each T9 leg re-scores at the realized MD law,
* the ``birthlaw`` stage - the V0-1 quantification, re-issued from this
  committed script as its first reproducible output.

Forward model: position (L1) x dressing (L2) x ladder (L3) x droplet prior
-> per-fragment chord dynamics (1D two-body ODE: Coulomb + solvation well +
gated pure-cubic drag; ballistic bracket = drag off) -> cooling exposure K
-> F.2b fate map (E_ej = E0_i * exp(-K); suppressed iff E_ej > Sigma(n_eject);
else exact eps=0 ladder descent; detected ~ energetic floor, s_eff=8
convention) -> solvated histogram + (n, mean-KE) curve -> thresholds T1-T5.

Pins (probe pins, production kinematics): tau = 6.55 ps, kappa = 1, mixture
picture, R0_sep = 2.666 A, b = 2.5153508541760052 amu*ps/A^2 (shared_pure_cubic
bundle), E_bind_ion = 0.11675778353879479 eV (bundle stamp), erf steepness
14.2 A (SimConfig default; drag/pickup/cooling share it), R(N) bulk-density
convention, complex mass = 126.90 + 4.0026*n amu (Tier-1a shell reference).

Recorded wiring oracles (``oracles`` stage; plan SI.11 V0-3): production
center-pin K = 0.74460, 9 A kinematics K = 0.89767, Sigma(21) = 0.18783720 eV.

Usage:  python tier2_h2b_forward_model.py oracles | levers | scan | report |
        w12pred | birthlaw | legaprime | legb | legc | legd | repilots1 |
        repilots2 | g3landmarks | g3scan
Outputs: CSVs + text summaries in OUT (see USER SETTINGS).

``g3landmarks`` (atlas §3.5 G3 Step 1, 2026-07-27): the twin landmark
re-issue at the G2-adopted corrected geometry -- S6 wiring oracle, the 11
Axis-A G1 cells vs their MD rows (the long-chord twin-authority
measurement), center-pin landmarks per grid radius, and the
corrected-ensemble row vs the D2b §4.3 forecast.

``g3scan`` (atlas §3.5c G3 Step 2, 2026-07-27): the nested Route A/B twin
scan at the corrected geometry -- landmark + S6 oracles, the
zero-integration pre-scan (analytic Route-A kill criterion), then the chord
surface (v_c x E_bind, one integration each, npz-cached) scored over the
free surface (tau x E0) against the §3.5c pre-registered gate.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from i2_helium_md.physics.constants import (  # noqa: E402
    D0_1_X2_WAVENUMBER,
    EV_PER_WAVENUMBER,
    MASS_HE_AMU,
    MASS_I_ION_AMU,
    droplet_radius_bulk_angstrom,
)
from i2_helium_md.physics.dissociation_ladder import d0_of_n  # noqa: E402
from i2_helium_md.physics.helium_density import rho_he_ratio  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.sampling.radial_positions import (  # noqa: E402
    sample_radial_positions,
)

# --------------------------------------------------------------------------
# USER SETTINGS
# --------------------------------------------------------------------------
# Output directory for CSVs / cached fragment tables (gitignored).
OUT = REPO / "data" / "runs" / "h2b_forward_model"

# Birth-position law for the master sample (V0-3 scope addition):
#   "uniform_volume" - the native L1 law (r^2 in droplet volume; hard margin
#                      applied downstream via margin_weight) - the law the
#                      Step-1b/1c closure was scored on;
#   "boltzmann"      - the repo's realized sampler (r^2 * exp(-U/kT) at
#                      T = 0.4 K; center-pinned) for T9 leg-A re-scores.
#                      Margins are meaningless under this law: only
#                      MARGINS_A == (0.0,) is accepted (guard below).
BIRTH_LAW = "uniform_volume"

# --------------------------------------------------------------------------
# Pinned constants (probe pins; provenance in module docstring)
# --------------------------------------------------------------------------
_BUNDLE = json.loads(
    (
        REPO
        / "data/reference/drag/shared/trajectory_matching/shared_pure_cubic/fit_parameters.json"
    ).read_text()
)
B_DRAG = float(_BUNDLE["b"])  # amu*ps/A^2 (pure cubic: F = rho_hat*b*v^3)
E_BIND_ION_EV = float(_BUNDLE["effective_binding_energy_I_ion_eV"])
assert abs(B_DRAG - 2.5153508541760052) < 1e-12
assert abs(E_BIND_ION_EV - 0.11675778353879479) < 1e-12

_CFG = SimConfig()
STEEP_A = float(_CFG.potential_steepness)  # 14.2 A; == drag_gate_steepness
assert STEEP_A == float(_CFG.drag_gate_steepness) == 14.2
KE_COUL_EVA = float(_CFG.E_coulomb_scale) * 14.39964548  # eV*A (q1=q2=1)

TAU_PS = 6.55  # GAH25 cooling time (probe pin)
KAPPA = 1.0
PICTURE = "statistical_mixture"
R0_SEP_PROD_A = 2.666  # production kinematics (E_avail = 5.401 eV pair)
R0_SEP_9A_A = 9.0  # validation kinematics (oracle only)
N_STAR = 21
X2_RUNG_EV = D0_1_X2_WAVENUMBER * EV_PER_WAVENUMBER  # 13.25 meV

EV_TO_AMU_A2_PS2 = 9648.53320892  # 1 eV in amu*A^2/ps^2  (EV*1e-4/U)
AMU_A2_PS2_TO_EV = 1.0 / EV_TO_AMU_A2_PS2

DT_PS = 0.01
T_END_PS = 150.0

E0_GRID = np.round(np.arange(0.20, 0.5001, 0.01), 3)  # sharp E0 scan (D1)
P_GRID = (0, 1)  # E0-dressing coupling brackets (D2)
# L3: floor variant + knobbed slid family; 'rq4graded' is a DIAGNOSTIC beyond
# the bounded set (the RQ4 pre-registered target ratios 2.2:1.5:1.3 on the
# first rungs) - reported separately as a prediction FOR the external
# calculation, never as a bounded lever.
LADDERS = ("flat", "floor1", "slid2", "slid3", "rq4graded")
MARGINS_A = (3.0, 4.67, 6.0)  # L1 firm band (D5)
PRIORS = ("kor0625", "d040", "d080", "pickup")  # D4
BRACKETS = ("current_law", "ballistic")  # D3
N_MOLECULES = 20000
SEED = 20260711

# Thresholds (D6)
T1_W1_MAX = 0.5
T2_LO, T2_HI = 1.75, 2.6
T3_LO, T3_HI = 0.26, 0.36
T4_FACTOR = 1.5
T5_HARD, T5_FLAG = 0.435, 0.15


def complex_mass_amu(n):
    return MASS_I_ION_AMU + np.asarray(n) * MASS_HE_AMU


def ladder_rungs(variant: str) -> np.ndarray:
    """1-indexed rung array D0(n), n = 1..N_STAR [eV]."""
    r = np.atleast_1d(
        d0_of_n(np.arange(1, N_STAR + 1), picture=PICTURE, kappa=KAPPA)
    ).astype(float)
    if variant == "flat":
        pass
    elif variant == "floor1":
        r[0] = X2_RUNG_EV
    elif variant == "slid2":
        r[:2] = X2_RUNG_EV
    elif variant == "slid3":
        r[:3] = X2_RUNG_EV
    elif variant == "rq4graded":
        r[0] *= 2.2
        r[1] *= 1.5
        r[2] *= 1.3
    else:
        raise ValueError(variant)
    return r


def sigma_cum(rungs: np.ndarray) -> np.ndarray:
    """Sigma[k] = sum of first k rungs, k = 0..N_STAR."""
    return np.concatenate(([0.0], np.cumsum(rungs)))


# --------------------------------------------------------------------------
# Chord dynamics: vectorised 1D two-body integrator (Heun / RK2)
# --------------------------------------------------------------------------
def drag_gamma_tail_amu_per_ps(v_abs, v_c, p_tail):
    """Capped-cubic friction coefficient gamma(v) [amu/ps] (Slice-T1 form).

    gamma = b*v^2 for v <= v_c; b*v_c^2*(v/v_c)^p_tail above the cap. Used
    only when the capped tail is selected (``integrate_pairs(v_c=...)``);
    the pure-cubic default keeps its original arithmetic verbatim.
    """
    v_abs = np.asarray(v_abs, dtype=float)
    gam_in = B_DRAG * v_abs**2
    ratio = np.where(v_abs > v_c, v_abs / v_c, 1.0)  # >= 1 where used
    gam_tail = B_DRAG * v_c**2 * ratio**p_tail
    return np.where(v_abs > v_c, gam_tail, gam_in)


def integrate_pairs(
    r0,
    mu,
    R_drop,
    mass_amu,
    r0_sep=R0_SEP_PROD_A,
    drag_on=True,
    dt=DT_PS,
    t_end=T_END_PS,
    track_profile=False,
    v_c=None,
    p_tail=None,
    e_bind_ev=None,
):
    """Integrate fragment pairs born at radius r0 (cos(pos,axis) = mu).

    Fragment A moves along +axis (cosine mu), B along -axis (cosine -mu).
    Forces per fragment: Coulomb (pair, along axis), droplet solvation
    (projected radial erf well, depth E_bind), gated cubic drag rho_hat*b*v^3.

    ``v_c``/``p_tail`` (leg-A' extension, 2026-07-16) select the Slice-T1
    ``capped_cubic`` tail: gamma = b*v^2 in-band, b*v_c^2*(v/v_c)^p_tail
    above the cap. ``v_c=None`` (default) is the pure-cubic law with its
    original arithmetic — byte-inert.

    ``e_bind_ev`` (g3scan extension, 2026-07-27) overrides the solvation-well
    depth [eV] for the E_bind chord axis (plan §3.5c; the §6.5/§6.7
    joint-pairing exception — off-bundle wells against the bundle-b drag
    law). ``None`` (default) keeps the bundle stamp ``E_BIND_ION_EV`` with
    identical arithmetic — byte-inert.

    Returns dict with per-fragment (2, M) arrays: K (cooling exposure),
    v_inf (asymptotic speed, A/ps, residual-Coulomb-corrected), t_exit (ps),
    v_peak; optional center-pin profile samples.
    """
    if v_c is not None and p_tail is None:
        raise ValueError("v_c requires p_tail (the capped-cubic tail exponent)")
    e_bind = E_BIND_ION_EV if e_bind_ev is None else float(e_bind_ev)
    r0 = np.atleast_1d(np.asarray(r0, dtype=float))
    mu = np.atleast_1d(np.asarray(mu, dtype=float))
    R_drop = np.atleast_1d(np.asarray(R_drop, dtype=float))
    mass_amu = np.atleast_1d(np.asarray(mass_amu, dtype=float))
    M = r0.size

    # state: displacement s and speed v per fragment (A, B)
    s = np.zeros((2, M))
    v = np.zeros((2, M))
    K = np.zeros((2, M))
    t_exit = np.full((2, M), np.nan)
    v_peak = np.zeros((2, M))
    mu2 = np.stack([mu, -mu])  # per-fragment cosine
    prof = {"t": [], "v": [], "depth": [], "rsep": []} if track_profile else None

    def geom(s_arr):
        a = 0.5 * r0_sep + s_arr  # distance from birth point, per fragment
        r = np.sqrt(r0**2 + a**2 + 2.0 * r0 * mu2 * a)
        dr_ds = (r0 * mu2 + a) / np.maximum(r, 1e-12)
        return r, dr_ds

    def accel(s_arr, v_arr):
        r, dr_ds = geom(s_arr)
        depth = r - R_drop
        rho = rho_he_ratio(depth, steepness=STEEP_A)
        r_sep = r0_sep + s_arr[0] + s_arr[1]
        F_c = KE_COUL_EVA / r_sep**2  # eV/A, outward along axis, both
        dUdr = (
            e_bind / (STEEP_A * np.sqrt(np.pi)) * np.exp(-((depth / STEEP_A) ** 2))
        )  # eV/A
        F_solv = -dUdr * dr_ds  # eV/A along motion direction
        F_ev = F_c + F_solv
        F = F_ev * EV_TO_AMU_A2_PS2
        if drag_on:
            if v_c is None:
                F = F - rho * B_DRAG * np.abs(v_arr) ** 2 * v_arr
            else:
                gam = drag_gamma_tail_amu_per_ps(np.abs(v_arr), v_c, p_tail)
                F = F - rho * gam * v_arr
        return F / mass_amu, rho, depth

    n_steps = int(round(t_end / dt))
    for i in range(n_steps):
        a1, rho1, depth = accel(s, v)
        # exposure (trapezoid via midpoint of rho over step ~ rho1 is fine at dt)
        K += rho1 * (dt / TAU_PS)
        s_pred = s + v * dt
        v_pred = v + a1 * dt
        a2, _, _ = accel(s_pred, v_pred)
        s = s + v * dt + 0.5 * a1 * dt**2
        v = v + 0.5 * (a1 + a2) * dt
        v_peak = np.maximum(v_peak, v)
        newly_out = np.isnan(t_exit) & (depth >= 0.0)
        t_exit[newly_out] = i * dt
        if track_profile and (i % 5 == 0):
            r, _ = geom(s)
            prof["t"].append(i * dt)
            prof["v"].append(v[0, 0])
            prof["depth"].append(r[0, 0] - R_drop[0])
            prof["rsep"].append(r0_sep + s[0, 0] + s[1, 0])

    # residual Coulomb -> asymptotic speed (half per fragment; equal-split
    # approximation, residual ~ 10 meV pair at t_end, error << 1 meV)
    r_sep_end = r0_sep + s[0] + s[1]
    resid_ev = KE_COUL_EVA / r_sep_end
    v_inf = np.sqrt(v**2 + 0.5 * resid_ev * EV_TO_AMU_A2_PS2 * 2.0 / mass_amu)
    r_end, _ = geom(s)
    depth_end = r_end - R_drop
    # trapped: dissipated below the solvation barrier -> parked in/near the
    # droplet at t_end (an escaping fragment at >= 1 A/ps has depth >> 2*steep)
    trapped = depth_end < 2.0 * STEEP_A
    out = {
        "K": K, "v_inf": v_inf, "t_exit": t_exit, "v_peak": v_peak,
        "depth_end": depth_end, "trapped": trapped,
    }
    if track_profile:
        out["profile"] = {k: np.asarray(vv) for k, vv in prof.items()}
    return out


# --------------------------------------------------------------------------
# Fate map (F.2b closed form, unified: E_ej = E0_i e^-K; searchsorted descent)
# --------------------------------------------------------------------------
def fate_map(n_eject, K, E0, p, sig):
    """Return (n_det, suppressed_mask) per fragment.

    E0_i = E0 * (Sigma(n_eject)/Sigma(21))**p ; E_ej = E0_i * exp(-K);
    suppressed iff E_ej > Sigma(n_eject) (self-unbound at ejection -> bare,
    RQ3 spec-b); opened: descend to smallest k with Sigma(k) > Sigma(n0)-E_ej.
    Born-bare (n_eject = 0) counted bare. The descent formula gives n_det = 0
    exactly for the suppressed branch too (same bookkeeping).
    """
    sig_n0 = sig[n_eject]
    E0_i = E0 * (sig_n0 / sig[N_STAR]) ** p if p else np.full(n_eject.shape, E0)
    E_ej = E0_i * np.exp(-K)
    suppressed = (E_ej > sig_n0) & (n_eject >= 1)
    target = sig_n0 - E_ej
    n_det = np.searchsorted(sig, target, side="right")
    n_det = np.minimum(n_det, n_eject)  # born-bare stays 0; no ties above n0
    return n_det.astype(int), suppressed


# --------------------------------------------------------------------------
# Experimental targets
# --------------------------------------------------------------------------
def load_experiment():
    rows = []
    with open(REPO / "data/reference/integrated_i_he_abundance.csv") as fh:
        for row in csv.DictReader(fh):
            rows.append((int(row["n"]), float(row["ionPercent"])))
    n_arr = np.array([r[0] for r in rows])
    pct = np.array([r[1] for r in rows])
    n_max = int(n_arr.max())
    h = np.zeros(n_max + 1)
    h[n_arr] = pct
    solv = h[1:] / h[1:].sum()  # renormalised n >= 1
    return solv, n_max


# user-supplied mean-KE table (design freeze D7; eV; n>=13 one-sided < 0.1)
KE_EXP_EV = {
    0: 2.9, 1: 0.974, 2: 0.545, 3: 0.390, 4: 0.341, 5: 0.289, 6: 0.248,
    7: 0.201, 8: 0.176, 9: 0.154, 10: 0.138, 11: 0.122, 12: 0.105,
}
KE_BAND_N13 = 0.1


def kinetic_energy_eV(mass_amu, v_aps):
    return 0.5 * mass_amu * v_aps**2 * AMU_A2_PS2_TO_EV


def wasserstein_bins(h_model, h_exp):
    n = max(h_model.size, h_exp.size)
    a = np.zeros(n)
    b = np.zeros(n)
    a[: h_model.size] = h_model
    b[: h_exp.size] = h_exp
    return float(np.abs(np.cumsum(a - b)).sum())


# --------------------------------------------------------------------------
# Sampling (importance): one master sample reweighted per prior x margin
# --------------------------------------------------------------------------
N_LO, N_HI = 250.0, 16000.0


def draw_master(rng, m, birth_law=None):
    """Master molecule sample: droplet size (uniform-in-lnN proposal), birth
    radius per the selected law (V0-3 birth-law input), axis cosine.

    ``uniform_volume``: x = U^(1/3) -> r^2 proposal on [0, R]; margins are
    applied downstream as importance weights (``margin_weight``).
    ``boltzmann``: the repo sampler's realized law (one rejection loop per
    unique droplet radius - slow under a continuous droplet prior; intended
    for fixed-N leg-A re-scores).
    """
    law = BIRTH_LAW if birth_law is None else birth_law
    lnN = rng.uniform(np.log(N_LO), np.log(N_HI), m)
    N = np.exp(lnN)
    R = droplet_radius_bulk_angstrom(N)
    if law == "uniform_volume":
        x = rng.uniform(0.0, 1.0, m) ** (1.0 / 3.0)  # r/R proposal ~ r^2 on [0, R]
    elif law == "boltzmann":
        r0 = sample_radial_positions(_CFG, R, rng=rng)
        x = r0 / R
    else:
        raise ValueError(f"unknown birth law {law!r}")
    mu = rng.uniform(-1.0, 1.0, m)
    return {"N": N, "R": R, "x": x, "mu": mu}


def prior_weight(prior, N):
    lnN = np.log(N)
    if prior == "kor0625":
        d = 0.625
    elif prior == "d040":
        d = 0.40
    elif prior == "d080":
        d = 0.80
    elif prior == "pickup":
        d = 0.625
    else:
        raise ValueError(prior)
    mu_ln = np.log(2000.0) - 0.5 * d**2  # mean <N> = 2000 convention
    w = np.exp(-0.5 * ((lnN - mu_ln) / d) ** 2)  # / proposal (uniform in lnN)
    if prior == "pickup":
        w = w * N ** (2.0 / 3.0)
    return w


def margin_weight(margin, x, R, birth_law=None):
    """Importance weight of the hard surface margin for a ``draw_master``
    sample. ``birth_law`` resolves exactly like ``draw_master``'s (None ->
    the module-global ``BIRTH_LAW``): a caller that drew with an explicit
    law must weight with the same one, or the 1/cap^3 uniform-volume
    weights would silently reweight a Boltzmann sample (review fix
    2026-07-18 -- the two halves of the importance scheme previously could
    not be kept in sync per-call)."""
    law = BIRTH_LAW if birth_law is None else birth_law
    if law == "boltzmann":
        # Margins belong to the uniform_volume law; under the realized
        # Boltzmann law the sample already carries its own radial density.
        if margin != 0.0:
            raise ValueError(
                "margins are defined for BIRTH_LAW='uniform_volume' only; "
                f"got margin={margin} under 'boltzmann' (set MARGINS_A = (0.0,))"
            )
        return np.ones_like(np.asarray(x, dtype=float))
    cap = 1.0 - margin / R
    ok = (x <= cap) & (cap > 0)
    with np.errstate(divide="ignore"):
        w = np.where(ok, 1.0 / cap**3, 0.0)
    return w


# --------------------------------------------------------------------------
# Stages
# --------------------------------------------------------------------------
def stage_oracles():
    print("=== H.2b wiring oracles ===")
    rungs = ladder_rungs("flat")
    sig = sigma_cum(rungs)
    print(f"O-L  Sigma(21) mixture k=1 = {sig[N_STAR]:.8f} eV  (target 0.18783720)")

    rho_c = rho_he_ratio(-droplet_radius_bulk_angstrom(2000.0), steepness=STEEP_A)
    print(f"O-L2 rho_hat(center, N=2000) = {rho_c:.4f} -> n_eject = "
          f"{int(np.rint(N_STAR * rho_c))}  (target 21)")

    m21 = complex_mass_amu(N_STAR)
    R2000 = droplet_radius_bulk_angstrom(2000.0)
    print(f"     R(2000) = {R2000:.3f} A ; m(21) = {m21:.4f} amu")

    # O2: center-pin production, current law
    res = integrate_pairs(
        np.array([0.0]), np.array([1.0]), np.array([R2000]), np.array([m21]),
        r0_sep=R0_SEP_PROD_A, drag_on=True, track_profile=True,
    )
    ke_det = kinetic_energy_eV(m21, res["v_inf"][0, 0])
    print(f"O2   production center-pin (current law):")
    print(f"     K       = {res['K'][0, 0]:.5f}   (Wave-11: 0.74603)")
    print(f"     expo*tau= {res['K'][0, 0] * TAU_PS:.4f} ps (Wave-11: 4.886)")
    print(f"     t_exit  = {res['t_exit'][0, 0]:.2f} ps  (Wave-11: 4.62)")
    print(f"     v_peak  = {res['v_peak'][0, 0]:.2f} A/ps (Wave-11: 10.51)")
    print(f"     v_inf   = {res['v_inf'][0, 0]:.2f} A/ps (Wave-11 detected: 4.11)")
    print(f"     KE_det(n=21 mass) = {ke_det:.4f} eV (Wave-11: ~0.18 at bare mass"
          f" -> {kinetic_energy_eV(complex_mass_amu(0), res['v_inf'][0, 0]):.4f})")

    # O3: 9 A kinematics center-pin
    res9 = integrate_pairs(
        np.array([0.0]), np.array([1.0]), np.array([R2000]), np.array([m21]),
        r0_sep=R0_SEP_9A_A, drag_on=True,
    )
    print(f"O3   9A center-pin K = {res9['K'][0, 0]:.5f}  (Wave-10/11: 0.89754)")

    # O-ball: ballistic production center-pin vs analytic
    resb = integrate_pairs(
        np.array([0.0]), np.array([1.0]), np.array([R2000]), np.array([m21]),
        r0_sep=R0_SEP_PROD_A, drag_on=False,
    )
    e_frag = 0.5 * KE_COUL_EVA / R0_SEP_PROD_A - E_BIND_ION_EV
    v_ana = np.sqrt(2.0 * e_frag * EV_TO_AMU_A2_PS2 / m21)
    print(f"O-B  ballistic v_inf = {resb['v_inf'][0, 0]:.4f} vs analytic "
          f"{v_ana:.4f} A/ps (rel {abs(resb['v_inf'][0, 0] / v_ana - 1):.2e})")

    # O4: fate map at center pin, flat ladder, measured K
    Kc = res["K"][0, 0]
    Estar = sig[N_STAR] * np.exp(Kc)
    n_det, sup = fate_map(
        np.array([N_STAR]), np.array([Kc]), 0.28, 0, sig
    )
    print(f"O4   E*(K) = {Estar:.4f} eV (Wave-11: 0.396); E0=0.28 -> n_det = "
          f"{n_det[0]} (Wave-11 closed form: 6; MD companion 5.76)")
    n_det2, sup2 = fate_map(np.array([N_STAR]), np.array([Kc]), 0.52, 0, sig)
    print(f"     E0=0.52 -> suppressed={bool(sup2[0])}, n_det={n_det2[0]} "
          f"(Wave-11 primary: suppressed -> bare)")

    # dt sensitivity spot check
    res_h = integrate_pairs(
        np.array([0.0]), np.array([1.0]), np.array([R2000]), np.array([m21]),
        r0_sep=R0_SEP_PROD_A, drag_on=True, dt=0.005,
    )
    print(f"O-dt K(dt=0.005) = {res_h['K'][0, 0]:.5f} (delta "
          f"{abs(res_h['K'][0, 0] - res['K'][0, 0]):.1e})")


def build_fragment_table(force_rebuild=False):
    """Master sample + both-bracket dynamics; cached to npz (per birth law)."""
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / f"h2b_fragments_{BIRTH_LAW}.npz"
    if cache.exists() and not force_rebuild:
        return dict(np.load(cache))
    rng = np.random.default_rng(SEED)
    ms = draw_master(rng, N_MOLECULES)
    r0 = ms["x"] * ms["R"]
    depth_birth = r0 - ms["R"]
    rho_b = rho_he_ratio(depth_birth, steepness=STEEP_A)
    n_eject = np.clip(np.rint(N_STAR * rho_b), 0, N_STAR).astype(int)
    mass = complex_mass_amu(n_eject)

    out = {
        "N": ms["N"], "R": ms["R"], "x": ms["x"], "mu": ms["mu"],
        "depth_birth": depth_birth, "n_eject": n_eject,
    }
    for tag, drag in (("cl", True), ("ball", False)):
        res = integrate_pairs(r0, ms["mu"], ms["R"], mass, drag_on=drag)
        out[f"K_{tag}"] = res["K"]
        out[f"vinf_{tag}"] = res["v_inf"]
        out[f"texit_{tag}"] = res["t_exit"]
        out[f"trapped_{tag}"] = res["trapped"]
    np.savez_compressed(cache, **out)
    return out


def stage_levers(tab):
    """Lever-interaction map (D5): reachable n_eject / chord-K per margin."""
    rows = []
    for margin in MARGINS_A:
        for prior in PRIORS:
            w_mol = prior_weight(prior, tab["N"]) * margin_weight(
                margin, tab["x"], tab["R"]
            )
            w = np.concatenate([w_mol, w_mol])
            ne = np.concatenate([tab["n_eject"], tab["n_eject"]])
            K = tab["K_cl"].reshape(-1)
            sel = w > 0
            q = lambda a, p: float(
                np.quantile(a[sel], p, weights=w[sel], method="inverted_cdf")
            )
            rows.append({
                "margin_A": margin, "prior": prior,
                "n_eject_min": int(ne[sel].min()),
                "n_eject_q05": q(ne.astype(float), 0.05),
                "n_eject_q50": q(ne.astype(float), 0.50),
                "K_q05": round(q(K, 0.05), 4),
                "K_q50": round(q(K, 0.50), 4),
                "K_q95": round(q(K, 0.95), 4),
                "frac_neject_le_10": round(
                    float(w[sel & (ne <= 10)].sum() / w[sel].sum()), 4),
            })
    with open(OUT / "h2b_lever_interaction_map.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wtr.writeheader()
        wtr.writerows(rows)
    # n_eject(depth) reference curve at N=2000
    R2000 = droplet_radius_bulk_angstrom(2000.0)
    print("depth[A] rho_hat n_eject Sigma(n_eject)[eV]  (N=2000)")
    sig = sigma_cum(ladder_rungs("flat"))
    for dpt in (-R2000, -20, -14, -10, -8, -6, -4.67, -4, -3):
        rr = rho_he_ratio(dpt, steepness=STEEP_A)
        ne = int(np.clip(np.rint(N_STAR * rr), 0, N_STAR))
        print(f"{dpt:8.2f} {rr:7.4f} {ne:7d} {sig[ne]:.4f}")
    print(f"lever map -> {OUT / 'h2b_lever_interaction_map.csv'}")


BRACKET_KEY = {"current_law": "cl", "ballistic": "ball"}


def score_cell(tab, sig, E0, p, w_frag, bracket, solv_exp):
    key = BRACKET_KEY.get(bracket, bracket)
    K = tab[f"K_{key}"].reshape(-1)
    trapped = tab[f"trapped_{key}"].reshape(-1).astype(bool)
    ne = np.concatenate([tab["n_eject"], tab["n_eject"]])
    n_det, sup = fate_map(ne, K, E0, p, sig)
    wt = w_frag / w_frag.sum()
    trapped_frac = float(wt[trapped].sum())
    w_det = np.where(trapped, 0.0, wt)  # trapped ions never reach detection
    det_tot = w_det.sum()
    if det_tot < 1e-12:
        return None
    bare = float(w_det[n_det == 0].sum() / det_tot)  # among detected fragments
    hist = np.bincount(n_det, weights=w_det, minlength=N_STAR + 1)
    solv = hist[1:]
    tot = solv.sum()
    if tot < 1e-12:
        return None
    solv = solv / tot
    W1 = wasserstein_bins(solv, solv_exp)
    n1, n2 = solv[0], solv[1]
    return {
        "n_det": n_det, "bare": bare, "solv": solv, "W1": W1,
        "n1": n1, "n2": n2, "ratio": n1 / n2 if n2 > 0 else np.inf,
        "trapped_frac": trapped_frac, "w_det": w_det,
    }


def ke_curve(tab, n_det, w_frag, bracket):
    v = tab[f"vinf_{bracket}"].reshape(-1)
    ke = kinetic_energy_eV(complex_mass_amu(n_det), v)
    means = {}
    for n in range(0, 22):
        m = n_det == n
        wsum = w_frag[m].sum()
        means[n] = float((w_frag[m] * ke[m]).sum() / wsum) if wsum > 0 else np.nan
    return means


def t4_check(tab, n_det, w_frag):
    """KE envelope + single-lambda feasibility on current-law binning."""
    ke_cl = ke_curve(tab, n_det, w_frag, "cl")
    ke_ball = ke_curve(tab, n_det, w_frag, "ball")
    lam_grid = np.linspace(0.0, 1.0, 101)
    best_lam, best_err = None, np.inf
    for lam in lam_grid:
        errs = []
        for n in range(1, 13):
            lo, hi = ke_cl.get(n, np.nan), ke_ball.get(n, np.nan)
            if not np.isfinite(lo) or not np.isfinite(hi):
                continue
            model = (1 - lam) * lo + lam * hi
            errs.append(abs(np.log(model / KE_EXP_EV[n])))
        if errs and max(errs) < best_err:
            best_err, best_lam = max(errs), lam
    enveloped = all(
        (ke_cl.get(n, np.inf) <= KE_EXP_EV[n] * T4_FACTOR)
        and (ke_ball.get(n, 0.0) >= KE_EXP_EV[n] / T4_FACTOR)
        for n in range(1, 13)
        if np.isfinite(ke_cl.get(n, np.nan)) and np.isfinite(ke_ball.get(n, np.nan))
    )
    ok = best_err < np.log(T4_FACTOR) if best_lam is not None else False
    return ok, enveloped, best_lam, np.exp(best_err) if np.isfinite(best_err) else np.nan, ke_cl, ke_ball


def stage_scan(tab):
    solv_exp, _ = load_experiment()
    rows = []
    n_frag_w = {}
    for prior in PRIORS:
        for margin in MARGINS_A:
            w_mol = prior_weight(prior, tab["N"]) * margin_weight(
                margin, tab["x"], tab["R"]
            )
            n_frag_w[(prior, margin)] = np.concatenate([w_mol, w_mol])
    ladd = {v: sigma_cum(ladder_rungs(v)) for v in LADDERS}
    for bracket in BRACKETS:
        for prior in PRIORS:
            for margin in MARGINS_A:
                w = n_frag_w[(prior, margin)]
                for lv in LADDERS:
                    sig = ladd[lv]
                    for p in P_GRID:
                        for E0 in E0_GRID:
                            sc = score_cell(tab, sig, E0, p, w, bracket, solv_exp)
                            if sc is None:
                                continue
                            t1 = sc["W1"] <= T1_W1_MAX
                            t2 = T2_LO <= sc["ratio"] <= T2_HI
                            t3 = T3_LO <= sc["n1"] <= T3_HI
                            t5_hard = sc["bare"] <= T5_HARD
                            row = {
                                "bracket": bracket, "prior": prior,
                                "margin_A": margin, "ladder": lv, "p": p,
                                "E0_eV": E0, "W1_bins": round(sc["W1"], 4),
                                "n1": round(sc["n1"], 4), "n2": round(sc["n2"], 4),
                                "n1_over_n2": round(sc["ratio"], 3),
                                "bare_frac": round(sc["bare"], 4),
                                "trapped_frac": round(sc["trapped_frac"], 4),
                                "T1": int(t1), "T2": int(t2), "T3": int(t3),
                                "T5_hard_ok": int(t5_hard),
                                "T5_flag": int(sc["bare"] > T5_FLAG),
                                "T123": int(t1 and t2 and t3),
                            }
                            rows.append(row)
    with open(OUT / "h2b_scan.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wtr.writeheader()
        wtr.writerows(rows)
    print(f"scan rows: {len(rows)} -> {OUT / 'h2b_scan.csv'}")

    # summary: best rows
    import pandas as pd  # local convenience if available

    df = pd.DataFrame(rows)
    print("\n=== best W1 per (bracket, ladder) [any prior/margin/p/E0] ===")
    idx = df.groupby(["bracket", "ladder"])["W1_bins"].idxmin()
    cols = ["bracket", "ladder", "prior", "margin_A", "p", "E0_eV", "W1_bins",
            "n1", "n1_over_n2", "bare_frac", "trapped_frac", "T123"]
    print(df.loc[idx, cols].to_string(index=False))
    print("\n=== T1&T2&T3 landing cells ===")
    landed = df[df["T123"] == 1]
    print(f"count = {len(landed)}")
    if len(landed):
        print(landed.sort_values("W1_bins").head(30)[cols].to_string(index=False))
    return df


def stage_report(tab):
    """T4 + feeder map on the best landed (or best-W1) current-law cells."""
    import pandas as pd

    solv_exp, _ = load_experiment()
    df = pd.read_csv(OUT / "h2b_scan.csv")
    dcl = df[df.bracket == "current_law"]
    landed = dcl[dcl.T123 == 1]
    pick = (landed if len(landed) else dcl).sort_values("W1_bins").head(8)
    print("=== T4 (KE feasibility) + detail on top current-law cells ===")
    frows = []
    for _, r in pick.iterrows():
        w = prior_weight(r.prior, tab["N"]) * margin_weight(
            r.margin_A, tab["x"], tab["R"]
        )
        w = np.concatenate([w, w])
        sig = sigma_cum(ladder_rungs(r.ladder))
        sc = score_cell(tab, sig, r.E0_eV, int(r.p), w, "cl", solv_exp)
        w = sc["w_det"]  # trapped fragments excluded from every detected read
        ok, env, lam, worst, ke_cl, ke_ball = t4_check(tab, sc["n_det"], w)
        print(f"\n[{r.ladder} p={int(r.p)} {r.prior} m={r.margin_A} "
              f"E0={r.E0_eV:.2f}]  W1={r.W1_bins:.3f} n1={r.n1:.3f} "
              f"ratio={r.n1_over_n2:.2f} bare={r.bare_frac:.3f} "
              f"trapped={r.trapped_frac:.3f}")
        print(f"  T4: single-lambda ok={ok} (lambda={lam}, worst x{worst:.2f}); "
              f"envelope={env}")
        print("  n : exp_KE  model_cl  model_ball")
        for n in range(1, 13):
            print(f"  {n:2d}: {KE_EXP_EV[n]:.3f}   {ke_cl.get(n, np.nan):.3f}"
                  f"     {ke_ball.get(n, np.nan):.3f}")
        hi = [n for n in range(13, 22) if np.isfinite(ke_cl.get(n, np.nan))]
        band_ok = all(
            (1 - lam) * ke_cl[n] + lam * ke_ball.get(n, 0.0) < KE_BAND_N13
            for n in hi if np.isfinite(ke_ball.get(n, np.nan))
        ) if lam is not None else False
        print(f"  n>=13 band (<0.1 eV at lambda): ok={band_ok}")
        if hi:
            print(f"  n>=13 (cl): " + ", ".join(
                f"{n}:{ke_cl[n]:.3f}" for n in hi))
        print("  solvated bins 1..8 model vs exp: "
              + " ".join(f"{sc['solv'][k]:.3f}/{solv_exp[k]:.3f}"
                         for k in range(8)))
        # feeder map: per-bin birth-depth decomposition
        ne = np.concatenate([tab["n_eject"], tab["n_eject"]])
        depth = np.concatenate([tab["depth_birth"], tab["depth_birth"]])
        wt = w / w.sum()
        for n in (0, 1, 2, 3):
            m = sc["n_det"] == n
            if wt[m].sum() < 1e-9:
                continue
            frows.append({
                "cell": f"{r.ladder}|p{int(r.p)}|{r.prior}|m{r.margin_A}|E0{r.E0_eV:.2f}",
                "bin_n": n, "weight": round(float(wt[m].sum()), 4),
                "mean_birth_depth_A": round(float(
                    (wt[m] * depth[m]).sum() / wt[m].sum()), 2),
                "mean_n_eject": round(float(
                    (wt[m] * ne[m]).sum() / wt[m].sum()), 2),
                "mean_K": round(float(
                    (wt[m] * tab["K_cl"].reshape(-1)[m]).sum() / wt[m].sum()), 3),
            })
    if frows:
        with open(OUT / "h2b_feeder_map.csv", "w", newline="") as fh:
            wtr = csv.DictWriter(fh, fieldnames=list(frows[0]))
            wtr.writeheader()
            wtr.writerows(frows)
        print(f"\nfeeder map -> {OUT / 'h2b_feeder_map.csv'}")


def stage_w12pred():
    """Forward-model predictions for the W12 MD leg as designed: DELIVERED
    physics (full dressing n_eject=21, flat mixture ladder, current law),
    fixed N = 2000, position sampling on; margin band {0, 4.67} A pending
    the W12 sampler audit. N_mc chords; fate at E0 in {0.28, 0.38} (p moot
    at full dressing)."""
    rng = np.random.default_rng(SEED + 1)
    m = 20000
    R2000 = float(droplet_radius_bulk_angstrom(2000.0))
    sig = sigma_cum(ladder_rungs("flat"))
    for margin in (0.0, 4.67):
        cap = 1.0 - margin / R2000
        r0 = (rng.uniform(0.0, 1.0, m) ** (1.0 / 3.0)) * cap * R2000
        mu = rng.uniform(-1.0, 1.0, m)
        mass = np.full(m, complex_mass_amu(N_STAR))
        print(f"\n=== W12 predictions (N=2000, full dressing, flat ladder, "
              f"margin={margin} A) ===")
        res = {}
        for tag, drag in (("cl", True), ("ball", False)):
            res[tag] = integrate_pairs(r0, mu, np.full(m, R2000), mass,
                                       drag_on=drag)
        K = res["cl"]["K"].reshape(-1)
        trapped = res["cl"]["trapped"].reshape(-1)
        v = res["cl"]["v_inf"].reshape(-1)
        print(f"trapped fraction (cl, 150 ps read): {trapped.mean():.3f}")
        print(f"K quantiles 5/25/50/75/95: "
              + "/".join(f"{np.quantile(K, q):.3f}"
                         for q in (0.05, 0.25, 0.5, 0.75, 0.95)))
        kA, kB = res["cl"]["K"][0], res["cl"]["K"][1]
        vA, vB = res["cl"]["v_inf"][0], res["cl"]["v_inf"][1]
        print(f"pair anti-correlation: corr(K_A,K_B) = "
              f"{np.corrcoef(kA, kB)[0, 1]:.3f}; corr(v_A,v_B) = "
              f"{np.corrcoef(vA, vB)[0, 1]:.3f}")
        ok = ~trapped
        print(f"detected speed range (cl, 5-95%): "
              f"{np.quantile(v[ok], 0.05):.2f}-{np.quantile(v[ok], 0.95):.2f} "
              f"A/ps (center-pin ref 4.13)")
        ne = np.full(2 * m, N_STAR)
        for E0 in (0.28, 0.38):
            n_det, sup = fate_map(ne, K, E0, 0, sig)
            w = np.where(trapped, 0.0, 1.0)
            w /= w.sum()
            supf = float(w[sup].sum())
            hist = np.bincount(n_det, weights=w, minlength=22)
            nbar = float((np.arange(22) * hist).sum() / hist.sum())
            top = ", ".join(f"n{k}:{hist[k]:.3f}" for k in range(0, 8)
                            if hist[k] > 0.005)
            print(f"E0={E0:.2f}: suppressed(bare) = {supf:.3f}; "
                  f"nbar_det = {nbar:.2f}; bins: {top}")
            ke = kinetic_energy_eV(complex_mass_amu(n_det), v)
            for k in (1, 2, 3):
                mm = (n_det == k) & ok
                if mm.sum() > 20:
                    print(f"    mean KE(n={k}) = "
                          f"{(w[mm] * ke[mm]).sum() / w[mm].sum():.3f} eV")


def stage_birthlaw():
    """V0-1 quantification, re-issued from the committed twin (plan SI.11 V0).

    Characterizes the realized birth-radius law at fixed N = 2000
    (R = 27.936 A): the repo Boltzmann sampler (analytic density on a fine
    grid, cross-checked by an empirical draw through
    ``sampling/radial_positions.py``) vs the native uniform_volume law at
    margin 0 and the three firm margins. Writes h2b_birth_law_quantiles.csv.
    Recorded V0-1 numbers (adjudication 2026-07-16): Boltzmann median
    1.37 A, 95% 2.75 A, 99% 3.41 A, mode 1.16 A - center-pinned.
    """
    from i2_helium_md.physics.constants import EV, K_B
    from i2_helium_md.physics.potentials import droplet_potential

    OUT.mkdir(parents=True, exist_ok=True)
    R2000 = float(droplet_radius_bulk_angstrom(2000.0))
    kT_eV = K_B / EV * float(_CFG.T_particles_K)
    binding_eV = float(_CFG.binding_energy_molecule_meV) / 1000.0
    dr = 0.01
    r = np.arange(dr / 2, 2.0 * R2000, dr)
    U = droplet_potential(
        r - R2000,
        steepness=float(_CFG.potential_steepness_molecule),
        binding_energy=binding_eV,
    )
    p = r**2 * np.exp(-U / kT_eV)
    p /= p.sum() * dr
    cdf = np.cumsum(p) * dr

    def q_ana(f):
        return float(r[np.searchsorted(cdf, f)])

    print("=== V0-1 birth-law quantification (N = 2000) ===")
    print(f"R(2000) = {R2000:.3f} A ; kT = {kT_eV * 1000:.4f} meV ; "
          f"E_b(molecule) = {binding_eV * 1000:.2f} meV ; "
          f"U(center)/kT = {U[0] / kT_eV:.2f}")
    qs = (0.05, 0.25, 0.50, 0.75, 0.95, 0.99)
    rows = [{
        "law": "boltzmann_analytic", "margin_A": 0.0,
        **{f"q{int(100 * f):02d}_A": round(q_ana(f), 3) for f in qs},
        "mode_A": round(float(r[np.argmax(p)]), 3),
    }]
    print("boltzmann (analytic):  "
          + "  ".join(f"q{int(100 * f)}={q_ana(f):.2f}" for f in qs)
          + f"  mode={r[np.argmax(p)]:.2f} A")

    # empirical cross-check through the actual repo sampler
    rng = np.random.default_rng(SEED)
    r_emp = sample_radial_positions(_CFG, np.full(4000, R2000), rng=rng)
    rows.append({
        "law": "boltzmann_sampled_n4000", "margin_A": 0.0,
        **{f"q{int(100 * f):02d}_A": round(float(np.quantile(r_emp, f)), 3)
           for f in qs},
        "mode_A": np.nan,
    })
    print("boltzmann (sampled):   "
          + "  ".join(f"q{int(100 * f)}={np.quantile(r_emp, f):.2f}" for f in qs))

    # uniform_volume comparators: quantile of r^2 law on [0, R - m]
    for margin in (0.0,) + MARGINS_A:
        cap = R2000 - margin
        rows.append({
            "law": "uniform_volume", "margin_A": margin,
            **{f"q{int(100 * f):02d}_A": round(cap * f ** (1.0 / 3.0), 3)
               for f in qs},
            "mode_A": round(cap, 3),
        })
        print(f"uniform_volume m={margin:4.2f}: "
              + "  ".join(f"q{int(100 * f)}={cap * f ** (1 / 3):.2f}"
                          for f in qs))
    path = OUT / "h2b_birth_law_quantiles.csv"
    with open(path, "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wtr.writeheader()
        wtr.writerows(rows)
    print(f"quantiles -> {path}")


# --------------------------------------------------------------------------
# T9 leg A' twin re-score (plan SI.11; executed 2026-07-16 leg trigger)
# --------------------------------------------------------------------------
LEG_APRIME_MARGIN_A = 3.0  # the Step-1c full-house cells' margin (findings S4k)

# The SI.10 confirmation matrix knobs (label, v_c, p_tail, tau_ps, ladder, E0):
# v_c=None = the current pure-cubic law (C3 control).
APRIME_CONFIGS = (
    ("c1", 7.5, -1.0, 3.8, "rq4graded", 0.25),
    ("c2", 6.5, 0.0, 4.0, "rq4graded", 0.24),
    ("c3", None, None, 6.55, "flat", 0.25),
    ("c4", 7.5, -1.0, 4.4, "floor1", 0.23),
)


def stage_legaprime(m=20000):
    """Twin re-score at the T9 leg-A / leg-A' MD configurations exactly.

    Leg-A' configuration (one lever flipped vs the delivered T3 dirs):
    fixed N = 2000 droplet (delta prior), birth law uniform_volume with the
    Step-1c full-house margin 3 A, **undressed** shell (n_eject = 21 for
    every fragment — the MD ANCHOR_N_START convention; T5 not flipped),
    E0-coupling p = 0 (T6 not flipped), per-C-config (drag tail, tau,
    ladder, E0) from the SI.10 matrix. Leg-A rows (center-pinned delta —
    the twin image of the delivered T3 dirs) are the oracle anchor.

    Ladder tables are imported from the MD generator so the fate map uses
    the MD's *exact* rung tables (incl. floor1's 13.3 meV absolute rung —
    the twin's own `ladder_rungs("floor1")` uses 13.25 meV and is NOT used
    here). tau enters as the exact K rescale K(tau) = K(6.55)*(6.55/tau)
    (SI.8 mechanics — tau appears only in the exposure bookkeeping).

    Writes h2b_leg_aprime_predictions.csv (histograms + classes) and
    h2b_leg_aprime_ke.csv (per-bin mean detected KE) — the pre-registered
    prediction record for the leg-A' MD pilots.
    """
    from scripts.gen_tier2_md_confirmation import (
        floor1_rungs_eV,
        form_u_rungs_eV,
        rq4graded_rungs_eV,
    )

    OUT.mkdir(parents=True, exist_ok=True)
    tables = {
        "flat": form_u_rungs_eV(),
        "rq4graded": rq4graded_rungs_eV(),
        "floor1": floor1_rungs_eV(),
    }
    R2000 = float(droplet_radius_bulk_angstrom(2000.0))
    m21 = float(complex_mass_amu(N_STAR))

    rng = np.random.default_rng(SEED)
    r0_ap = (R2000 - LEG_APRIME_MARGIN_A) * np.cbrt(rng.uniform(0.0, 1.0, m))
    mu_ap = rng.uniform(-1.0, 1.0, m)
    legs = {
        "a": (np.zeros(1), np.ones(1)),          # center-pinned delta (T3 twin)
        "aprime": (r0_ap, mu_ap),                # uniform_volume, margin 3 A
    }

    hist_rows, ke_rows = [], []
    for label, v_c, p_tail, tau, ladder_key, e0 in APRIME_CONFIGS:
        sig = sigma_cum(np.asarray(tables[ladder_key][:N_STAR], dtype=float))
        for leg, (r0, mu) in legs.items():
            M = r0.size
            res = integrate_pairs(
                r0, mu, np.full(M, R2000), np.full(M, m21),
                r0_sep=R0_SEP_PROD_A, drag_on=True, v_c=v_c, p_tail=p_tail,
            )
            K = res["K"].reshape(-1) * (TAU_PS / tau)  # exact tau rescale
            trapped = res["trapped"].reshape(-1).astype(bool)
            v = res["v_inf"].reshape(-1)
            ne = np.full(2 * M, N_STAR)
            n_det, sup = fate_map(ne, K, e0, 0, sig)
            w = np.where(trapped, 0.0, 1.0)
            w_tot = w.sum()
            trapped_frac = float(trapped.mean())
            sup_frac = float(w[sup].sum() / w_tot)
            hist = np.bincount(n_det, weights=w, minlength=N_STAR + 1) / w_tot
            nbar = float((np.arange(N_STAR + 1) * hist).sum())
            row = {
                "leg": leg, "config": label, "ladder": ladder_key,
                "tau_ps": tau, "E0_eV": e0,
                "v_c": "" if v_c is None else v_c,
                "p_tail": "" if p_tail is None else p_tail,
                "trapped_frac": round(trapped_frac, 4),
                "suppressed_frac": round(sup_frac, 4),
                "bare_frac": round(float(hist[0]), 4),
                "nbar_det": round(nbar, 3),
                "K_q05": round(float(np.quantile(K, 0.05)), 4),
                "K_q50": round(float(np.quantile(K, 0.50)), 4),
                "K_q95": round(float(np.quantile(K, 0.95)), 4),
            }
            row.update({f"h{k}": round(float(hist[k]), 4)
                        for k in range(N_STAR + 1)})
            hist_rows.append(row)
            ke = kinetic_energy_eV(complex_mass_amu(n_det), v)
            for k in range(0, N_STAR + 1):
                mask = (n_det == k) & ~trapped
                if mask.sum() >= (1 if leg == "a" else 20):
                    ke_rows.append({
                        "leg": leg, "config": label, "n": k,
                        "weight": round(float(w[mask].sum() / w_tot), 4),
                        "mean_KE_eV": round(
                            float((w[mask] * ke[mask]).sum() / w[mask].sum()), 4
                        ),
                    })
            top = ", ".join(
                f"n{k}:{hist[k]:.3f}" for k in range(N_STAR + 1) if hist[k] > 0.02
            )
            print(f"[{leg:6s} {label}] trapped={trapped_frac:.3f} "
                  f"supp={sup_frac:.3f} nbar={nbar:.2f}  {top}")

    with open(OUT / "h2b_leg_aprime_predictions.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(hist_rows[0]))
        wtr.writeheader()
        wtr.writerows(hist_rows)
    with open(OUT / "h2b_leg_aprime_ke.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(ke_rows[0]))
        wtr.writeheader()
        wtr.writerows(ke_rows)
    print(f"predictions -> {OUT / 'h2b_leg_aprime_predictions.csv'}")
    print(f"KE table    -> {OUT / 'h2b_leg_aprime_ke.csv'}")


# --------------------------------------------------------------------------
# T9 leg B twin re-score (plan SI.11 Slice T5; dressed configuration)
# --------------------------------------------------------------------------
def stage_legb(m=20000):
    """Twin re-score at the T9 leg-B MD configuration exactly.

    Leg-B configuration (one lever flipped vs leg A''): the Slice-T5
    initial-shell dressing ``density_tied`` — per-fragment
    n_eject = clip(rint(N_STAR * rho_hat(r_birth - R2000)), 0, N_STAR)
    through the shared erf-complement surface (STEEP_A == the MD's
    drag_gate_steepness), chords integrated at the dressed birth mass
    complex_mass_amu(n_eject) (the H.2b L2 convention, build_fragment_table
    verbatim). Everything else rides the leg-A' twin: fixed N = 2000 delta
    prior, uniform_volume births at margin 3 A, E0-coupling p = 0 (T6 not
    flipped), per-C-config (drag tail, tau, ladder, E0) from the SI.10
    matrix via the MD generator's exact rung tables. The MD's co-moving
    shed convention (leg A'') is the twin's native KE bookkeeping, so the
    KE table is directly comparable — every KE value below is on the
    **co-moving basis** (the item-3 stated-basis amendment).

    The ``b_undressed`` rows re-run the leg-A' ensemble (same seed, same
    draws, n_eject = 21) as the in-stage wiring oracle — they must equal
    stage_legaprime's ``aprime`` rows exactly.

    Twin-divergence channels *listed* for the MD A/B (S2c-P4): pickup
    re-filling after under-dressed birth (the twin has no live Langmuir
    channel — the MD may land above the twin at small n_eject); the MD
    dresses per ATOM at ion-t0 positions (±R0/2 chord offset + neutral
    drift) while the twin dresses per molecule center at birth; no trapped
    dynamics in the twin beyond the 150 ps chord read.

    Writes h2b_leg_b_predictions.csv (histograms + classes + n_eject
    quantiles) and h2b_leg_b_ke.csv (per-bin mean detected KE) — the
    pre-registered prediction record for the leg-B MD pilots.
    """
    from scripts.gen_tier2_md_confirmation import (
        floor1_rungs_eV,
        form_u_rungs_eV,
        rq4graded_rungs_eV,
    )

    OUT.mkdir(parents=True, exist_ok=True)
    tables = {
        "flat": form_u_rungs_eV(),
        "rq4graded": rq4graded_rungs_eV(),
        "floor1": floor1_rungs_eV(),
    }
    R2000 = float(droplet_radius_bulk_angstrom(2000.0))
    m21 = float(complex_mass_amu(N_STAR))

    # Identical draws to stage_legaprime (same SEED, same order) so the
    # undressed anchor is bit-comparable and the dressed leg differs only
    # by the dressing law — never by sampling noise.
    rng = np.random.default_rng(SEED)
    r0_ap = (R2000 - LEG_APRIME_MARGIN_A) * np.cbrt(rng.uniform(0.0, 1.0, m))
    mu_ap = rng.uniform(-1.0, 1.0, m)

    # Slice-T5 dressing at the molecule birth center (twin convention;
    # fragments of a pair share n_eject — mass-symmetric pairs).
    rho_b = rho_he_ratio(r0_ap - R2000, steepness=STEEP_A)
    ne_mol = np.clip(np.rint(N_STAR * rho_b), 0, N_STAR).astype(int)

    legs = {
        "b_undressed": (np.full(m, N_STAR), np.full(m, m21)),
        "b": (ne_mol, complex_mass_amu(ne_mol).astype(float)),
    }

    hist_rows, ke_rows = [], []
    for label, v_c, p_tail, tau, ladder_key, e0 in APRIME_CONFIGS:
        sig = sigma_cum(np.asarray(tables[ladder_key][:N_STAR], dtype=float))
        for leg, (ne_m, mass_m) in legs.items():
            res = integrate_pairs(
                r0_ap, mu_ap, np.full(m, R2000), mass_m,
                r0_sep=R0_SEP_PROD_A, drag_on=True, v_c=v_c, p_tail=p_tail,
            )
            K = res["K"].reshape(-1) * (TAU_PS / tau)  # exact tau rescale
            trapped = res["trapped"].reshape(-1).astype(bool)
            v = res["v_inf"].reshape(-1)
            ne = np.concatenate([ne_m, ne_m])  # (2, M).reshape(-1) order
            n_det, sup = fate_map(ne, K, e0, 0, sig)
            w = np.where(trapped, 0.0, 1.0)
            w_tot = w.sum()
            if w_tot == 0.0:
                raise ValueError(
                    f"[{leg} {label}] every fragment is trapped (w_tot == 0): the "
                    "histogram / fraction normalisation would be 0/0 and write NaN "
                    "into the pre-registered prediction CSV. A config that traps the "
                    "whole ensemble cannot produce a prediction record."
                )
            trapped_frac = float(trapped.mean())
            sup_frac = float(w[sup].sum() / w_tot)
            hist = np.bincount(n_det, weights=w, minlength=N_STAR + 1) / w_tot
            nbar = float((np.arange(N_STAR + 1) * hist).sum())
            nef = ne.astype(float)
            row = {
                "leg": leg, "config": label, "ladder": ladder_key,
                "tau_ps": tau, "E0_eV": e0,
                "v_c": "" if v_c is None else v_c,
                "p_tail": "" if p_tail is None else p_tail,
                "trapped_frac": round(trapped_frac, 4),
                "suppressed_frac": round(sup_frac, 4),
                "bare_frac": round(float(hist[0]), 4),
                "nbar_det": round(nbar, 3),
                "K_q05": round(float(np.quantile(K, 0.05)), 4),
                "K_q50": round(float(np.quantile(K, 0.50)), 4),
                "K_q95": round(float(np.quantile(K, 0.95)), 4),
                "n_eject_q05": round(float(np.quantile(nef, 0.05)), 2),
                "n_eject_q50": round(float(np.quantile(nef, 0.50)), 2),
                "n_eject_mean": round(float(nef.mean()), 3),
            }
            row.update({f"h{k}": round(float(hist[k]), 4)
                        for k in range(N_STAR + 1)})
            hist_rows.append(row)
            ke = kinetic_energy_eV(complex_mass_amu(n_det), v)
            for k in range(0, N_STAR + 1):
                mask = (n_det == k) & ~trapped
                if mask.sum() >= 20:
                    ke_rows.append({
                        "leg": leg, "config": label, "n": k,
                        "weight": round(float(w[mask].sum() / w_tot), 4),
                        "mean_KE_eV": round(
                            float((w[mask] * ke[mask]).sum() / w[mask].sum()), 4
                        ),
                    })
            top = ", ".join(
                f"n{k}:{hist[k]:.3f}" for k in range(N_STAR + 1) if hist[k] > 0.02
            )
            print(f"[{leg:11s} {label}] trapped={trapped_frac:.3f} "
                  f"supp={sup_frac:.3f} nbar={nbar:.2f} "
                  f"ne_mean={nef.mean():.2f}  {top}")

    with open(OUT / "h2b_leg_b_predictions.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(hist_rows[0]))
        wtr.writeheader()
        wtr.writerows(hist_rows)
    with open(OUT / "h2b_leg_b_ke.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(ke_rows[0]))
        wtr.writeheader()
        wtr.writerows(ke_rows)
    print(f"predictions -> {OUT / 'h2b_leg_b_predictions.csv'}")
    print(f"KE table    -> {OUT / 'h2b_leg_b_ke.csv'}")


# --------------------------------------------------------------------------
# T9 leg C twin re-score (plan SI.11 Slice T6; E_int(0)-dressing p-law)
# --------------------------------------------------------------------------
def stage_legc(m=20000):
    """Twin re-score at the T9 leg-C MD configuration exactly.

    Leg-C configuration (one lever flipped vs leg B): the Slice-T6
    E_int(0)-dressing coupling ``sigma_proportional`` (p = 1) — the S2 onset
    is regularised to the dressed shell,

        E0_i = E0 * (Sigma(n_eject)/Sigma(21))**1 ,

    so under-dressed births (small n_eject near the surface) get proportionally
    less onset. Everything else rides leg B verbatim: the Slice-T5
    ``density_tied`` dressing (per-fragment n_eject via the shared
    erf-complement surface, chords at the dressed mass), fixed N = 2000 delta
    prior, ``uniform_volume`` births at margin 3 A, the co-moving KE basis
    (leg A''), and the per-C-config knobs (drag tail, tau, ladder, E0) from the
    SI.10 matrix. The MD flips exactly the same one lever
    (``internal_energy_partition_law = "sigma_proportional"``, T6) on the
    certified leg-B ``bc`` baseline — the twin's p-law (fate_map arg 4) and the
    MD's (``internal_energy_budget.sigma_partition_factor``) are the same
    ``(Sigma(n0)/Sigma(n*))**p`` closed form, so the A/B is valid.

    The ``c_p0`` rows re-run the leg-B **dressed** ensemble (same seed, same
    draws, p = 0) as the in-stage wiring oracle — they must equal stage_legb's
    ``b`` rows exactly.

    S2c-P2 mechanism (pre-registered direction): p = 1 lowers the onset
    (ratio <= 1), so it de-suppresses the over-suppressed under-dressed births
    (p = 0 over-suppresses — the SS4j direction), shifting weight off the
    suppressed/bare side into the shallow solvated bins (suppressed_frac down,
    nbar_det up vs p = 0).

    Twin-divergence channels *listed* for the MD A/B (S2c-P4): as leg B (pickup
    re-filling after under-dressed birth; per-atom vs per-molecule-center
    dressing; no trapped dynamics beyond the 150 ps chord read). Experimental
    scoring is deferred to T9 on the **solvated branch** (bare renormalised out;
    post-leg-B decision).

    Writes h2b_leg_c_predictions.csv (histograms + classes + n_eject quantiles)
    and h2b_leg_c_ke.csv (per-bin mean detected KE, co-moving basis) — the
    pre-registered prediction record for the leg-C MD pilots.
    """
    from scripts.gen_tier2_md_confirmation import (
        floor1_rungs_eV,
        form_u_rungs_eV,
        rq4graded_rungs_eV,
    )

    OUT.mkdir(parents=True, exist_ok=True)
    tables = {
        "flat": form_u_rungs_eV(),
        "rq4graded": rq4graded_rungs_eV(),
        "floor1": floor1_rungs_eV(),
    }
    R2000 = float(droplet_radius_bulk_angstrom(2000.0))

    # Identical draws to stage_legb (same SEED, same order) so the dressed
    # p = 0 oracle is bit-comparable and leg C differs from leg B only by the
    # p-law exponent — never by sampling noise.
    rng = np.random.default_rng(SEED)
    r0_ap = (R2000 - LEG_APRIME_MARGIN_A) * np.cbrt(rng.uniform(0.0, 1.0, m))
    mu_ap = rng.uniform(-1.0, 1.0, m)

    # Slice-T5 dressing at the molecule birth center (leg-B convention verbatim).
    rho_b = rho_he_ratio(r0_ap - R2000, steepness=STEEP_A)
    ne_mol = np.clip(np.rint(N_STAR * rho_b), 0, N_STAR).astype(int)
    mass_mol = complex_mass_amu(ne_mol).astype(float)

    # Both legs share the dressed shell/mass; only the onset coupling p differs.
    legs = {
        "c_p0": 0,   # dressed, p = 0 == leg-B `b` (the wiring oracle)
        "c": 1,      # dressed, p = 1 (leg C)
    }

    hist_rows, ke_rows = [], []
    for label, v_c, p_tail, tau, ladder_key, e0 in APRIME_CONFIGS:
        sig = sigma_cum(np.asarray(tables[ladder_key][:N_STAR], dtype=float))
        for leg, p_law in legs.items():
            res = integrate_pairs(
                r0_ap, mu_ap, np.full(m, R2000), mass_mol,
                r0_sep=R0_SEP_PROD_A, drag_on=True, v_c=v_c, p_tail=p_tail,
            )
            K = res["K"].reshape(-1) * (TAU_PS / tau)  # exact tau rescale
            trapped = res["trapped"].reshape(-1).astype(bool)
            v = res["v_inf"].reshape(-1)
            ne = np.concatenate([ne_mol, ne_mol])  # (2, M).reshape(-1) order
            n_det, sup = fate_map(ne, K, e0, p_law, sig)
            w = np.where(trapped, 0.0, 1.0)
            w_tot = w.sum()
            if w_tot == 0.0:
                raise ValueError(
                    f"[{leg} {label}] every fragment is trapped (w_tot == 0): the "
                    "histogram / fraction normalisation would be 0/0 and write NaN "
                    "into the pre-registered prediction CSV. A config that traps the "
                    "whole ensemble cannot produce a prediction record."
                )
            trapped_frac = float(trapped.mean())
            sup_frac = float(w[sup].sum() / w_tot)
            hist = np.bincount(n_det, weights=w, minlength=N_STAR + 1) / w_tot
            nbar = float((np.arange(N_STAR + 1) * hist).sum())
            nef = ne.astype(float)
            row = {
                "leg": leg, "config": label, "ladder": ladder_key,
                "tau_ps": tau, "E0_eV": e0,
                "v_c": "" if v_c is None else v_c,
                "p_tail": "" if p_tail is None else p_tail,
                "trapped_frac": round(trapped_frac, 4),
                "suppressed_frac": round(sup_frac, 4),
                "bare_frac": round(float(hist[0]), 4),
                "nbar_det": round(nbar, 3),
                "K_q05": round(float(np.quantile(K, 0.05)), 4),
                "K_q50": round(float(np.quantile(K, 0.50)), 4),
                "K_q95": round(float(np.quantile(K, 0.95)), 4),
                "n_eject_q05": round(float(np.quantile(nef, 0.05)), 2),
                "n_eject_q50": round(float(np.quantile(nef, 0.50)), 2),
                "n_eject_mean": round(float(nef.mean()), 3),
            }
            row.update({f"h{k}": round(float(hist[k]), 4)
                        for k in range(N_STAR + 1)})
            hist_rows.append(row)
            ke = kinetic_energy_eV(complex_mass_amu(n_det), v)
            for k in range(0, N_STAR + 1):
                mask = (n_det == k) & ~trapped
                if mask.sum() >= 20:
                    ke_rows.append({
                        "leg": leg, "config": label, "n": k,
                        "weight": round(float(w[mask].sum() / w_tot), 4),
                        "mean_KE_eV": round(
                            float((w[mask] * ke[mask]).sum() / w[mask].sum()), 4
                        ),
                    })
            top = ", ".join(
                f"n{k}:{hist[k]:.3f}" for k in range(N_STAR + 1) if hist[k] > 0.02
            )
            print(f"[{leg:5s} {label}] trapped={trapped_frac:.3f} "
                  f"supp={sup_frac:.3f} nbar={nbar:.2f} "
                  f"ne_mean={nef.mean():.2f}  {top}")

    with open(OUT / "h2b_leg_c_predictions.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(hist_rows[0]))
        wtr.writeheader()
        wtr.writerows(hist_rows)
    with open(OUT / "h2b_leg_c_ke.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(ke_rows[0]))
        wtr.writeheader()
        wtr.writerows(ke_rows)
    print(f"predictions -> {OUT / 'h2b_leg_c_predictions.csv'}")
    print(f"KE table    -> {OUT / 'h2b_leg_c_ke.csv'}")


def stage_legd(m=20000):
    """Twin re-score at the T9 leg-D MD configuration exactly.

    Leg-D configuration (one lever flipped vs leg C): the Slice-T8 droplet
    prior ``kornilov_lognormal`` — per-molecule N drawn exactly from the
    truncated ln-normal (mu_ln = ln(2000) − 0.625²/2 on [N_LO, N_HI]; the
    package sampler ``sample_droplet_sizes_analytic``, T8-D2/D3, so the twin
    and the MD share the prior by construction) with R(N) the bulk-density
    convention. Everything else rides leg C verbatim: ``uniform_volume``
    births at the absolute margin 3 A *per droplet* ([0, R_i − 3]), the
    Slice-T5 ``density_tied`` dressing at the per-molecule surface
    (rho_hat(r_birth − R_i)), the Slice-T6 p = 1 onset coupling, the
    co-moving KE basis, and the per-C-config knobs (drag tail, tau, ladder,
    E0) from the SI.10 matrix. The MD flips exactly the same one lever
    (``droplet_size_prior = "kornilov_lognormal"``, T8) on the certified
    leg-C ``cc`` baseline.

    The ``d_delta`` rows re-run the leg-C dressed p = 1 ensemble (same seed,
    same u/mu draws, delta N = 2000) as the in-stage wiring oracle — they
    must equal stage_legc's ``c`` rows exactly. The ``d`` rows scale the SAME
    u/mu draws to the per-molecule radii (maximum draw coupling: d vs
    d_delta differs only through the droplet axis, never sampling noise).

    Twin-divergence channels *listed* for the MD A/B (S2c-P4): as leg C
    (pickup re-filling after under-dressed birth; per-atom vs
    per-molecule-center dressing; no trapped dynamics beyond the 150 ps
    chord read), plus the droplet-axis-specific pair: the MD well/gate
    depths follow R_i dynamically during the cascade (the twin's chord is
    frozen at birth geometry), and small-droplet chords truncate earlier.
    Experimental scoring is deferred to T9 on the **solvated branch**.

    Writes h2b_leg_d_predictions.csv (histograms + classes + n_eject and N
    quantiles) and h2b_leg_d_ke.csv (per-bin mean detected KE, co-moving
    basis) — the pre-registered prediction record for the leg-D MD pilots.
    """
    from dataclasses import replace as _replace

    from i2_helium_md.sampling.droplet_sizes import sample_droplet_sizes_analytic

    from scripts.gen_tier2_md_confirmation import (
        floor1_rungs_eV,
        form_u_rungs_eV,
        rq4graded_rungs_eV,
    )

    OUT.mkdir(parents=True, exist_ok=True)
    tables = {
        "flat": form_u_rungs_eV(),
        "rq4graded": rq4graded_rungs_eV(),
        "floor1": floor1_rungs_eV(),
    }
    R2000 = float(droplet_radius_bulk_angstrom(2000.0))

    # Identical u/mu draws to stage_legc (same SEED, same order) so the
    # d_delta oracle is bit-comparable; the prior draw comes AFTER the legc
    # stream so it cannot disturb it.
    rng = np.random.default_rng(SEED)
    u_ap = rng.uniform(0.0, 1.0, m)
    mu_ap = rng.uniform(-1.0, 1.0, m)
    cfg_d = _replace(
        _CFG,
        droplet_size_prior="kornilov_lognormal",
        use_single_droplet_size=False,
        num_molecules=m,
    )
    N_d = sample_droplet_sizes_analytic(cfg_d, rng=rng)
    R_d = droplet_radius_bulk_angstrom(N_d)

    # Per-leg geometry: births r^2 on [0, R − 3] at each leg's radii, and the
    # T5 dressing at each leg's surface (leg-C convention, per-molecule R).
    legs = {}
    for leg, R_leg, N_leg in (
        ("d_delta", np.full(m, R2000), np.full(m, 2000.0)),
        ("d", np.asarray(R_d, dtype=float), np.asarray(N_d, dtype=float)),
    ):
        r0 = (R_leg - LEG_APRIME_MARGIN_A) * np.cbrt(u_ap)
        rho = rho_he_ratio(r0 - R_leg, steepness=STEEP_A)
        ne_mol = np.clip(np.rint(N_STAR * rho), 0, N_STAR).astype(int)
        legs[leg] = {
            "R": R_leg, "N": N_leg, "r0": r0, "ne_mol": ne_mol,
            "mass_mol": complex_mass_amu(ne_mol).astype(float),
        }

    hist_rows, ke_rows = [], []
    for label, v_c, p_tail, tau, ladder_key, e0 in APRIME_CONFIGS:
        sig = sigma_cum(np.asarray(tables[ladder_key][:N_STAR], dtype=float))
        for leg, g in legs.items():
            res = integrate_pairs(
                g["r0"], mu_ap, g["R"], g["mass_mol"],
                r0_sep=R0_SEP_PROD_A, drag_on=True, v_c=v_c, p_tail=p_tail,
            )
            K = res["K"].reshape(-1) * (TAU_PS / tau)  # exact tau rescale
            trapped = res["trapped"].reshape(-1).astype(bool)
            v = res["v_inf"].reshape(-1)
            ne = np.concatenate([g["ne_mol"], g["ne_mol"]])  # (2, M) order
            n_det, sup = fate_map(ne, K, e0, 1, sig)  # p = 1 (leg-C carry)
            w = np.where(trapped, 0.0, 1.0)
            w_tot = w.sum()
            if w_tot == 0.0:
                raise ValueError(
                    f"[{leg} {label}] every fragment is trapped (w_tot == 0): the "
                    "histogram / fraction normalisation would be 0/0 and write NaN "
                    "into the pre-registered prediction CSV. A config that traps the "
                    "whole ensemble cannot produce a prediction record."
                )
            trapped_frac = float(trapped.mean())
            sup_frac = float(w[sup].sum() / w_tot)
            hist = np.bincount(n_det, weights=w, minlength=N_STAR + 1) / w_tot
            nbar = float((np.arange(N_STAR + 1) * hist).sum())
            nef = ne.astype(float)
            row = {
                "leg": leg, "config": label, "ladder": ladder_key,
                "tau_ps": tau, "E0_eV": e0,
                "v_c": "" if v_c is None else v_c,
                "p_tail": "" if p_tail is None else p_tail,
                "trapped_frac": round(trapped_frac, 4),
                "suppressed_frac": round(sup_frac, 4),
                "bare_frac": round(float(hist[0]), 4),
                "nbar_det": round(nbar, 3),
                "K_q05": round(float(np.quantile(K, 0.05)), 4),
                "K_q50": round(float(np.quantile(K, 0.50)), 4),
                "K_q95": round(float(np.quantile(K, 0.95)), 4),
                "n_eject_q05": round(float(np.quantile(nef, 0.05)), 2),
                "n_eject_q50": round(float(np.quantile(nef, 0.50)), 2),
                "n_eject_mean": round(float(nef.mean()), 3),
                "N_q05": round(float(np.quantile(g["N"], 0.05)), 1),
                "N_q50": round(float(np.quantile(g["N"], 0.50)), 1),
                "N_q95": round(float(np.quantile(g["N"], 0.95)), 1),
            }
            row.update({f"h{k}": round(float(hist[k]), 4)
                        for k in range(N_STAR + 1)})
            hist_rows.append(row)
            ke = kinetic_energy_eV(complex_mass_amu(n_det), v)
            for k in range(0, N_STAR + 1):
                mask = (n_det == k) & ~trapped
                if mask.sum() >= 20:
                    ke_rows.append({
                        "leg": leg, "config": label, "n": k,
                        "weight": round(float(w[mask].sum() / w_tot), 4),
                        "mean_KE_eV": round(
                            float((w[mask] * ke[mask]).sum() / w[mask].sum()), 4
                        ),
                    })
            top = ", ".join(
                f"n{k}:{hist[k]:.3f}" for k in range(N_STAR + 1) if hist[k] > 0.02
            )
            print(f"[{leg:7s} {label}] trapped={trapped_frac:.3f} "
                  f"supp={sup_frac:.3f} nbar={nbar:.2f} "
                  f"ne_mean={nef.mean():.2f}  {top}")

    with open(OUT / "h2b_leg_d_predictions.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(hist_rows[0]))
        wtr.writeheader()
        wtr.writerows(hist_rows)
    with open(OUT / "h2b_leg_d_ke.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(ke_rows[0]))
        wtr.writeheader()
        wtr.writerows(ke_rows)
    print(f"predictions -> {OUT / 'h2b_leg_d_predictions.csv'}")
    print(f"KE table    -> {OUT / 'h2b_leg_d_ke.csv'}")


# --------------------------------------------------------------------------
# SI.11.4 re-pilot Stage-1 pre-registration (v_c brackets at the leg-D config)
# --------------------------------------------------------------------------
# RP-D1/D2 (adjudicated 2026-07-20): core arms c1 (rq4graded, p = -1) and c4
# (floor1, p = -1) plus the c2 (p = 0) spot leg; C3 dropped (control
# discharged). Brackets: 5 v_c values spanning +-1.0 A/ps about the Step-1c
# centers (c1/c4: 7.5; c2: 6.5) -- the SI.11.4.2 working proposal, frozen at
# this trigger.
REPILOT_S1_VC_BRACKETS = {
    "c1": (6.5, 7.0, 7.5, 8.0, 8.5),
    "c2": (5.5, 6.0, 6.5, 7.0, 7.5),
    "c4": (6.5, 7.0, 7.5, 8.0, 8.5),
}


def repilot_s1_cell_label(config, v_c):
    """``s1c1v75``-style cell label (v_c x 10, two digits -- run-tag-safe)."""
    return f"s1{config}v{round(v_c * 10):02d}"


def stage_repilot1(m=20000):
    """Twin re-score over the Stage-1 v_c brackets at the leg-D configuration.

    The SI.11.4.5 pre-registration record for the Stage-1 KE-pin sweep:
    every cell rides the full leg-D configuration (kornilov_lognormal
    delta = 0.625 about <N> = 2000; uniform_volume births at margin 3 A per
    droplet; density_tied dressing; p = 1 onset coupling; co-moving KE
    basis; per-config tau / ladder / E0 from the SI.10 matrix) with ONLY
    v_c moving through the bracket. Draw discipline is stage_legd's
    verbatim (same SEED, same u/mu order, prior draw after): cells differ
    only through the drag tail, never through sampling noise, and the
    chord integration is cached per unique (v_c, p_tail) -- c1 and c4
    share it by construction, which IS the registered KE-degeneracy claim.

    In-stage wiring oracle (S1-P1): the bracket-center cells (the C-matrix
    v_c values) must reproduce the on-disk ``h2b_leg_d_predictions.csv`` /
    ``h2b_leg_d_ke.csv`` ``d`` rows exactly on every shared column (same
    draws, same computation -- run at the same m as the committed record).

    Registered-authority note (I69/RP-D5): the histogram/suppression
    columns of these rows carry ordering/direction authority only; the KE
    columns carry the x1.10 quantitative band.

    Writes h2b_repilot_s1_predictions.csv / h2b_repilot_s1_ke.csv
    (leg = "s1", config = the ``s1c1v75``-style cell label).
    """
    from dataclasses import replace as _replace

    from i2_helium_md.sampling.droplet_sizes import sample_droplet_sizes_analytic

    from scripts.gen_tier2_md_confirmation import (
        floor1_rungs_eV,
        form_u_rungs_eV,
        rq4graded_rungs_eV,
    )

    OUT.mkdir(parents=True, exist_ok=True)
    tables = {
        "flat": form_u_rungs_eV(),
        "rq4graded": rq4graded_rungs_eV(),
        "floor1": floor1_rungs_eV(),
    }

    # stage_legd's draw discipline, verbatim (u/mu first, prior after).
    rng = np.random.default_rng(SEED)
    u_ap = rng.uniform(0.0, 1.0, m)
    mu_ap = rng.uniform(-1.0, 1.0, m)
    cfg_d = _replace(
        _CFG,
        droplet_size_prior="kornilov_lognormal",
        use_single_droplet_size=False,
        num_molecules=m,
    )
    N_d = np.asarray(sample_droplet_sizes_analytic(cfg_d, rng=rng), dtype=float)
    R_d = np.asarray(droplet_radius_bulk_angstrom(N_d), dtype=float)

    r0 = (R_d - LEG_APRIME_MARGIN_A) * np.cbrt(u_ap)
    rho = rho_he_ratio(r0 - R_d, steepness=STEEP_A)
    ne_mol = np.clip(np.rint(N_STAR * rho), 0, N_STAR).astype(int)
    mass_mol = complex_mass_amu(ne_mol).astype(float)
    ne = np.concatenate([ne_mol, ne_mol])  # (2, M) order
    nef = ne.astype(float)

    configs = {c[0]: c for c in APRIME_CONFIGS}
    chord_cache = {}
    hist_rows, ke_rows = [], []
    for config, bracket in REPILOT_S1_VC_BRACKETS.items():
        _, v_c_center, p_tail, tau, ladder_key, e0 = configs[config]
        sig = sigma_cum(np.asarray(tables[ladder_key][:N_STAR], dtype=float))
        for v_c in bracket:
            key = (v_c, p_tail)
            if key not in chord_cache:
                chord_cache[key] = integrate_pairs(
                    r0, mu_ap, R_d, mass_mol,
                    r0_sep=R0_SEP_PROD_A, drag_on=True,
                    v_c=v_c, p_tail=p_tail,
                )
            res = chord_cache[key]
            K = res["K"].reshape(-1) * (TAU_PS / tau)  # exact tau rescale
            trapped = res["trapped"].reshape(-1).astype(bool)
            v = res["v_inf"].reshape(-1)
            n_det, sup = fate_map(ne, K, e0, 1, sig)  # p = 1 (leg-C carry)
            w = np.where(trapped, 0.0, 1.0)
            w_tot = w.sum()
            if w_tot == 0.0:
                raise ValueError(
                    f"[s1 {config} v_c={v_c}] every fragment is trapped "
                    "(w_tot == 0) -- no prediction record possible."
                )
            cell = repilot_s1_cell_label(config, v_c)
            trapped_frac = float(trapped.mean())
            sup_frac = float(w[sup].sum() / w_tot)
            hist = np.bincount(n_det, weights=w, minlength=N_STAR + 1) / w_tot
            nbar = float((np.arange(N_STAR + 1) * hist).sum())
            row = {
                "leg": "s1", "config": cell, "ladder": ladder_key,
                "tau_ps": tau, "E0_eV": e0,
                "v_c": v_c, "p_tail": p_tail,
                "trapped_frac": round(trapped_frac, 4),
                "suppressed_frac": round(sup_frac, 4),
                "bare_frac": round(float(hist[0]), 4),
                "nbar_det": round(nbar, 3),
                "K_q05": round(float(np.quantile(K, 0.05)), 4),
                "K_q50": round(float(np.quantile(K, 0.50)), 4),
                "K_q95": round(float(np.quantile(K, 0.95)), 4),
                "n_eject_q05": round(float(np.quantile(nef, 0.05)), 2),
                "n_eject_q50": round(float(np.quantile(nef, 0.50)), 2),
                "n_eject_mean": round(float(nef.mean()), 3),
                "N_q05": round(float(np.quantile(N_d, 0.05)), 1),
                "N_q50": round(float(np.quantile(N_d, 0.50)), 1),
                "N_q95": round(float(np.quantile(N_d, 0.95)), 1),
            }
            row.update({f"h{k}": round(float(hist[k]), 4)
                        for k in range(N_STAR + 1)})
            hist_rows.append(row)
            ke = kinetic_energy_eV(complex_mass_amu(n_det), v)
            for k in range(0, N_STAR + 1):
                mask = (n_det == k) & ~trapped
                if mask.sum() >= 20:
                    ke_rows.append({
                        "leg": "s1", "config": cell, "n": k,
                        "weight": round(float(w[mask].sum() / w_tot), 4),
                        "mean_KE_eV": round(
                            float((w[mask] * ke[mask]).sum() / w[mask].sum()), 4
                        ),
                    })
            top = ", ".join(
                f"n{k}:{hist[k]:.3f}" for k in range(N_STAR + 1) if hist[k] > 0.02
            )
            print(f"[s1 {cell}] trapped={trapped_frac:.3f} "
                  f"supp={sup_frac:.3f} nbar={nbar:.2f}  {top}")

    # ---- S1-P1 in-stage wiring oracle: bracket centers == the leg-D `d` rows
    legd_hist_path = OUT / "h2b_leg_d_predictions.csv"
    legd_ke_path = OUT / "h2b_leg_d_ke.csv"
    if not (legd_hist_path.exists() and legd_ke_path.exists()):
        raise FileNotFoundError(
            "leg-D prediction CSVs not found next to the Stage-1 output -- "
            "the S1-P1 wiring oracle needs the committed leg-D record (run "
            "stage_legd at the same m first)."
        )
    with open(legd_hist_path, newline="") as fh:
        legd_hist = {r["config"]: r for r in csv.DictReader(fh)
                     if r["leg"] == "d"}
    with open(legd_ke_path, newline="") as fh:
        legd_ke = [r for r in csv.DictReader(fh) if r["leg"] == "d"]
    for config, bracket in REPILOT_S1_VC_BRACKETS.items():
        v_c_center = configs[config][1]
        cell = repilot_s1_cell_label(config, v_c_center)
        mine = next(r for r in hist_rows if r["config"] == cell)
        ref = legd_hist[config]
        for key, val in mine.items():
            if key in ("leg", "config"):
                continue
            if str(val) != ref[key]:
                raise AssertionError(
                    f"S1-P1 oracle FAILED: cell {cell} vs leg-D `d` {config} "
                    f"differ at {key}: {val!r} != {ref[key]!r}"
                )
        mine_ke = [(r["n"], str(r["weight"]), str(r["mean_KE_eV"]))
                   for r in ke_rows if r["config"] == cell]
        ref_ke = [(int(r["n"]), r["weight"], r["mean_KE_eV"])
                  for r in legd_ke if r["config"] == config]
        if mine_ke != ref_ke:
            raise AssertionError(
                f"S1-P1 oracle FAILED: KE rows of {cell} differ from the "
                f"leg-D `d` {config} rows."
            )
    print("S1-P1 wiring oracle PASSED: all three bracket-center cells "
          "reproduce the leg-D `d` rows exactly (predictions + KE).")

    with open(OUT / "h2b_repilot_s1_predictions.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(hist_rows[0]))
        wtr.writeheader()
        wtr.writerows(hist_rows)
    with open(OUT / "h2b_repilot_s1_ke.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(ke_rows[0]))
        wtr.writeheader()
        wtr.writerows(ke_rows)
    print(f"predictions -> {OUT / 'h2b_repilot_s1_predictions.csv'}")
    print(f"KE table    -> {OUT / 'h2b_repilot_s1_ke.csv'}")


# --------------------------------------------------------------------------
# SI.11.4 re-pilot Stage-2 pre-registration ((tau, E0) grid, two-candidate
# v_c carry -- user adjudication (b), 2026-07-20)
# --------------------------------------------------------------------------
# Core arms c1/c4 at v_c in {6.5, 8.5} (the anchor/ratio cell + the chi2
# edge cell -- Stage 1 measured the two-sidedness, I72); tau x E0 grid
# centered on the C-matrix values (grid centers are the on-disk Stage-1
# cells). Spot leg: c2 at v_c = 6.5 along a grid diagonal.
REPILOT_S2_TAU_GRID = (3.2, 3.8, 4.4, 5.0)
REPILOT_S2_E0_GRID = (0.21, 0.23, 0.25, 0.27)
REPILOT_S2_VC_CARRY = (6.5, 8.5)
REPILOT_S2_SPOT = (  # (config, v_c, tau_ps, E0_eV)
    ("c2", 6.5, 3.2, 0.21),
    ("c2", 6.5, 4.0, 0.24),   # the on-disk s1c2v65 point (oracle)
    ("c2", 6.5, 5.0, 0.27),
)


def repilot_s2_cell_label(config, v_c, tau, e0):
    """``s2c1v65t38e25``-style cell label (all knobs x10/x100, tag-safe)."""
    return (f"s2{config}v{round(v_c * 10):02d}"
            f"t{round(tau * 10):02d}e{round(e0 * 100):02d}")


def repilot_s2_cells():
    """The frozen Stage-2 cell list: (config, v_c, tau_ps, E0_eV)."""
    cells = [
        (config, v_c, tau, e0)
        for config in ("c1", "c4")
        for v_c in REPILOT_S2_VC_CARRY
        for tau in REPILOT_S2_TAU_GRID
        for e0 in REPILOT_S2_E0_GRID
    ]
    cells.extend(REPILOT_S2_SPOT)
    return cells


def stage_repilot2(m=20000):
    """Twin re-score over the Stage-2 (tau, E0) grid at the carried v_c.

    The SI.11.4.5 pre-registration record for the Stage-2 histogram sweep
    under adjudication (b): every cell rides the full leg-D configuration
    with (v_c, tau, E0) set per cell; the drag chord is cached per
    (v_c, p_tail) (tau and E0 enter only the exposure bookkeeping and the
    fate map -- SI.8 mechanics), so all 35 cells of one drag tail share
    one integration.

    In-stage wiring oracle (S2s-P1): the four grid-center cells (the
    C-matrix (tau, E0) at each carried v_c) plus the spot center must
    reproduce the committed ``h2b_repilot_s1_*`` rows exactly on every
    shared column (same draws, same computation, same m).

    Registered-authority note (I73 refinement of I69/RP-D5): the v65
    cells' histogram rows carry *measured near-quantitative* authority
    (Stage 1: dnbar ~ 0, W1 0.21-0.26 at v65); the v85 rows carry the
    measured -0.6..-0.8 He channel-(d) softening; KE columns keep the
    x1.10 band.

    Writes h2b_repilot_s2_predictions.csv / h2b_repilot_s2_ke.csv
    (leg = "s2", config = the ``s2c1v65t38e25``-style cell label).
    """
    from dataclasses import replace as _replace

    from i2_helium_md.sampling.droplet_sizes import sample_droplet_sizes_analytic

    from scripts.gen_tier2_md_confirmation import (
        floor1_rungs_eV,
        form_u_rungs_eV,
        rq4graded_rungs_eV,
    )

    OUT.mkdir(parents=True, exist_ok=True)
    tables = {
        "flat": form_u_rungs_eV(),
        "rq4graded": rq4graded_rungs_eV(),
        "floor1": floor1_rungs_eV(),
    }

    # stage_legd / stage_repilot1 draw discipline, verbatim.
    rng = np.random.default_rng(SEED)
    u_ap = rng.uniform(0.0, 1.0, m)
    mu_ap = rng.uniform(-1.0, 1.0, m)
    cfg_d = _replace(
        _CFG,
        droplet_size_prior="kornilov_lognormal",
        use_single_droplet_size=False,
        num_molecules=m,
    )
    N_d = np.asarray(sample_droplet_sizes_analytic(cfg_d, rng=rng), dtype=float)
    R_d = np.asarray(droplet_radius_bulk_angstrom(N_d), dtype=float)

    r0 = (R_d - LEG_APRIME_MARGIN_A) * np.cbrt(u_ap)
    rho = rho_he_ratio(r0 - R_d, steepness=STEEP_A)
    ne_mol = np.clip(np.rint(N_STAR * rho), 0, N_STAR).astype(int)
    mass_mol = complex_mass_amu(ne_mol).astype(float)
    ne = np.concatenate([ne_mol, ne_mol])
    nef = ne.astype(float)

    configs = {c[0]: c for c in APRIME_CONFIGS}
    chord_cache = {}
    hist_rows, ke_rows = [], []
    for config, v_c, tau, e0 in repilot_s2_cells():
        _, _, p_tail, _, ladder_key, _ = configs[config]
        sig = sigma_cum(np.asarray(tables[ladder_key][:N_STAR], dtype=float))
        key = (v_c, p_tail)
        if key not in chord_cache:
            chord_cache[key] = integrate_pairs(
                r0, mu_ap, R_d, mass_mol,
                r0_sep=R0_SEP_PROD_A, drag_on=True,
                v_c=v_c, p_tail=p_tail,
            )
        res = chord_cache[key]
        K = res["K"].reshape(-1) * (TAU_PS / tau)
        trapped = res["trapped"].reshape(-1).astype(bool)
        v = res["v_inf"].reshape(-1)
        n_det, sup = fate_map(ne, K, e0, 1, sig)  # p = 1 (leg-C carry)
        w = np.where(trapped, 0.0, 1.0)
        w_tot = w.sum()
        if w_tot == 0.0:
            raise ValueError(
                f"[s2 {config} v_c={v_c} tau={tau} E0={e0}] every fragment "
                "trapped -- no prediction record possible."
            )
        cell = repilot_s2_cell_label(config, v_c, tau, e0)
        trapped_frac = float(trapped.mean())
        sup_frac = float(w[sup].sum() / w_tot)
        hist = np.bincount(n_det, weights=w, minlength=N_STAR + 1) / w_tot
        nbar = float((np.arange(N_STAR + 1) * hist).sum())
        row = {
            "leg": "s2", "config": cell, "ladder": ladder_key,
            "tau_ps": tau, "E0_eV": e0,
            "v_c": v_c, "p_tail": p_tail,
            "trapped_frac": round(trapped_frac, 4),
            "suppressed_frac": round(sup_frac, 4),
            "bare_frac": round(float(hist[0]), 4),
            "nbar_det": round(nbar, 3),
            "K_q05": round(float(np.quantile(K, 0.05)), 4),
            "K_q50": round(float(np.quantile(K, 0.50)), 4),
            "K_q95": round(float(np.quantile(K, 0.95)), 4),
            "n_eject_q05": round(float(np.quantile(nef, 0.05)), 2),
            "n_eject_q50": round(float(np.quantile(nef, 0.50)), 2),
            "n_eject_mean": round(float(nef.mean()), 3),
            "N_q05": round(float(np.quantile(N_d, 0.05)), 1),
            "N_q50": round(float(np.quantile(N_d, 0.50)), 1),
            "N_q95": round(float(np.quantile(N_d, 0.95)), 1),
        }
        row.update({f"h{k}": round(float(hist[k]), 4)
                    for k in range(N_STAR + 1)})
        hist_rows.append(row)
        ke = kinetic_energy_eV(complex_mass_amu(n_det), v)
        for k in range(0, N_STAR + 1):
            mask = (n_det == k) & ~trapped
            if mask.sum() >= 20:
                ke_rows.append({
                    "leg": "s2", "config": cell, "n": k,
                    "weight": round(float(w[mask].sum() / w_tot), 4),
                    "mean_KE_eV": round(
                        float((w[mask] * ke[mask]).sum() / w[mask].sum()), 4
                    ),
                })
        print(f"[s2 {cell}] trapped={trapped_frac:.3f} supp={sup_frac:.3f} "
              f"nbar={nbar:.2f}")

    # ---- S2s-P1 in-stage wiring oracle: grid centers == the Stage-1 rows
    s1_hist_path = OUT / "h2b_repilot_s1_predictions.csv"
    s1_ke_path = OUT / "h2b_repilot_s1_ke.csv"
    if not (s1_hist_path.exists() and s1_ke_path.exists()):
        raise FileNotFoundError(
            "Stage-1 prediction CSVs not found -- the S2s-P1 oracle needs "
            "the committed Stage-1 record (run stage_repilot1 at the same "
            "m first)."
        )
    with open(s1_hist_path, newline="") as fh:
        s1_hist = {r["config"]: r for r in csv.DictReader(fh)}
    with open(s1_ke_path, newline="") as fh:
        s1_ke_rows = [r for r in csv.DictReader(fh)]
    center_map = {
        repilot_s2_cell_label("c1", 6.5, 3.8, 0.25): "s1c1v65",
        repilot_s2_cell_label("c1", 8.5, 3.8, 0.25): "s1c1v85",
        repilot_s2_cell_label("c4", 6.5, 4.4, 0.23): "s1c4v65",
        repilot_s2_cell_label("c4", 8.5, 4.4, 0.23): "s1c4v85",
        repilot_s2_cell_label("c2", 6.5, 4.0, 0.24): "s1c2v65",
    }
    for s2_cell, s1_cell in center_map.items():
        mine = next(r for r in hist_rows if r["config"] == s2_cell)
        ref = s1_hist[s1_cell]
        for key, val in mine.items():
            if key in ("leg", "config"):
                continue
            if str(val) != ref[key]:
                raise AssertionError(
                    f"S2s-P1 oracle FAILED: {s2_cell} vs {s1_cell} differ "
                    f"at {key}: {val!r} != {ref[key]!r}"
                )
        mine_ke = [(r["n"], str(r["weight"]), str(r["mean_KE_eV"]))
                   for r in ke_rows if r["config"] == s2_cell]
        ref_ke = [(int(r["n"]), r["weight"], r["mean_KE_eV"])
                  for r in s1_ke_rows if r["config"] == s1_cell]
        if mine_ke != ref_ke:
            raise AssertionError(
                f"S2s-P1 oracle FAILED: KE rows of {s2_cell} differ from "
                f"{s1_cell}."
            )
    print("S2s-P1 wiring oracle PASSED: all five grid-center cells "
          "reproduce the Stage-1 rows exactly (predictions + KE).")

    with open(OUT / "h2b_repilot_s2_predictions.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(hist_rows[0]))
        wtr.writeheader()
        wtr.writerows(hist_rows)
    with open(OUT / "h2b_repilot_s2_ke.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(ke_rows[0]))
        wtr.writeheader()
        wtr.writerows(ke_rows)
    print(f"predictions -> {OUT / 'h2b_repilot_s2_predictions.csv'}")
    print(f"KE table    -> {OUT / 'h2b_repilot_s2_ke.csv'}")


# --------------------------------------------------------------------------
# G3 Step 1 — twin landmark re-issue at the corrected geometry
# (TIER2_SENSITIVITY_ATLAS_PLAN.md §3.5 G3, first concrete step; executed
# under its own `[PROCEED TO IMPLEMENTATION]`, 2026-07-27)
# --------------------------------------------------------------------------
# The standing production point finc1v725 (drag_migration_log_tier2.md):
# rq4graded ladder, capped-cubic v_c 7.25 / p_tail -1, tau 3.2 ps, E0 0.27 eV.
G3_STANDING = ("rq4graded", 7.25, -1.0, 3.2, 0.27)

# Parent-consistent Boltzmann well depth (G0 decision 2: per-run override,
# the config default stays 573.3 K).
G3_CORRECTED_WELL_K = 313.2

# One seed for every G3 draw — the Axis-A grid seed family (CRN across cells:
# each Part-B cell re-seeds the SAME value, so cross-cell differences are
# partly paired, mirroring the MD grid's shared-seed discipline).
G3_SEED = 20260727

# Part-A oracle cells: the three committed S6 finalist rows
# (h2b_s6_final_{predictions,ke}.csv; probe plan I.11.5.5).
G3_S6_CELLS = (
    ("rq4graded", 7.25, -1.0, 3.2, 0.27),
    ("rq4graded", 7.5, -1.0, 3.2, 0.27),
    ("floor1", 7.0, -1.0, 3.8, 0.25),
)

# Part-B cells: the Axis-A G1 grid (gen_tier2atlas_geometry.py; R via the
# bulk convention from the cell's fixed N_he). Laws: l1 center-pin (r0 = 0),
# l2 parent Boltzmann at 313.2 K, l3 uniform_volume margin 3 A.
G3_GRID_CELLS = (
    ("r1l1", 1727, "l1"), ("r1l2", 1727, "l2"), ("r1l3", 1727, "l3"),
    ("r2l1", 3605, "l1"), ("r2l2", 3605, "l2"), ("r2l3", 3605, "l3"),
    ("r3l1", 11059, "l1"), ("r3l2", 11059, "l2"), ("r3l3", 11059, "l3"),
    ("r4l2", 29227, "l2"), ("r4l3", 29227, "l3"),
)

# The MD-scored G1 rows this stage measures the twin against
# (TIER2_PARAMETER_INFLUENCE.md §14.1; committed scorer: midHot geometric
# n2-8, deepKE arithmetic n10-17; trap = t_b + t_m, i.e. the twin's single
# trapped class is compared against the MD bound+marginal TOTAL; NaN = the
# MD band was empty). Columns: trap, supp, nbar_det, n1_solv, W1_solv,
# midhot_geo, deepke.
G3_MD_G1_ROWS = {
    "r1l1": (0.0, 0.0, 8.39, 0.0, 4.98, 1.04, 1.83),
    "r1l2": (0.0, 0.0, 8.35, 0.0, 4.76, 1.25, 1.58),
    "r1l3": (0.059, 0.167, 4.46, 0.235, 0.81, 1.02, 0.51),
    "r2l1": (0.0, 0.0, 12.58, 0.0, 7.96, np.nan, 1.49),
    "r2l2": (0.0, 0.0, 12.30, 0.0, 7.53, 1.42, 1.38),
    "r2l3": (0.196, 0.140, 5.33, 0.184, 1.40, 1.11, 0.75),
    "r3l1": (0.538, 0.0, 18.35, 0.0, 13.46, np.nan, 1.94),
    "r3l2": (0.476, 0.0, 15.12, 0.0, 10.24, np.nan, 2.19),
    "r3l3": (0.404, 0.123, 6.30, 0.151, 2.31, 1.31, 1.59),
    "r4l2": (0.800, 0.0, 15.28, 0.0, 10.39, np.nan, 1.93),
    "r4l3": (0.570, 0.126, 6.71, 0.128, 2.83, 1.33, 1.97),
}

G3_MD_OBS_KEYS = ("trap", "supp", "nbar", "n1_solv", "w1", "midhot_geo",
                  "deepke")


def g3_ref_mean_ke():
    """Committed per-n reference mean KE (IHe_KED_reference.csv, n = 0..17)."""
    from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference

    ref = load_ihe_ked_reference(
        REPO / "data/reference/ihe_ked/IHe_KED_reference.csv"
    )
    return {int(n): float(k) for n, k in zip(ref.n, ref.mean_KE_eV)}


def g3_score(n_det, sup, trapped, v, ref_ke, solv_exp, min_bin_count=20):
    """Score one twin cell with the S6 columns + the committed KE bands.

    Returns ``(obs, ke_bins)``. ``obs`` holds unrounded values:
    ``trapped_frac`` / ``sup_frac`` / ``hist`` (bare+solvated, non-trapped
    normalised) / ``nbar`` / ``n1_solv`` / ``ratio`` / ``w1`` /
    ``midhot_arith`` (the S6-frozen arithmetic n2-8 convention) /
    ``midhot_geo`` + ``deepke`` (+ their bin counts; the committed
    ``tier2_confirmation`` conventions: geometric n2-8 / arithmetic n10-17
    of the per-bin sim/ref-mean ratios) / ``n1_ke``. ``ke_bins`` lists
    ``(n, weight, mean_KE_eV)`` for non-trapped bins with at least
    ``min_bin_count`` fragments (the S6 KE-row convention). Bands read only
    bins present in ``ke_bins``; empty band -> NaN with count 0.
    """
    w = np.where(trapped, 0.0, 1.0)
    w_tot = w.sum()
    if w_tot == 0.0:
        raise ValueError("every fragment is trapped -- no scoreable ensemble")
    trapped_frac = float(trapped.mean())
    sup_frac = float(w[sup].sum() / w_tot)
    hist = np.bincount(n_det, weights=w, minlength=N_STAR + 1) / w_tot
    nbar = float((np.arange(N_STAR + 1) * hist).sum())
    solv_tot = hist[1:].sum()
    if solv_tot <= 0.0:
        raise ValueError("no solvated mass -- the solvated columns are 0/0")
    solv = hist[1:] / solv_tot
    n1_solv = float(solv[0])
    ratio = float(hist[1] / hist[2]) if hist[2] > 0 else np.nan
    w1 = wasserstein_bins(solv, solv_exp)
    ke = kinetic_energy_eV(complex_mass_amu(n_det), v)
    ke_mean, ke_bins = {}, []
    for k in range(0, N_STAR + 1):
        mask = (n_det == k) & ~trapped
        if mask.sum() >= min_bin_count:
            mk = float((w[mask] * ke[mask]).sum() / w[mask].sum())
            ke_mean[k] = mk
            ke_bins.append((k, float(w[mask].sum() / w_tot), mk))

    def band(lo, hi, aggregate):
        ratios = [ke_mean[n] / ref_ke[n] for n in range(lo, hi + 1)
                  if n in ke_mean and n in ref_ke]
        if not ratios:
            return np.nan, 0
        if aggregate == "geometric":
            return float(np.exp(np.mean(np.log(ratios)))), len(ratios)
        return float(np.mean(ratios)), len(ratios)

    midhot_arith, _ = band(2, 8, "arithmetic")
    midhot_geo, midhot_bins = band(2, 8, "geometric")
    deepke, deepke_bins = band(10, 17, "arithmetic")
    return {
        "trapped_frac": trapped_frac, "sup_frac": sup_frac, "hist": hist,
        "nbar": nbar, "n1_solv": n1_solv, "ratio": ratio, "w1": w1,
        "midhot_arith": midhot_arith, "midhot_geo": midhot_geo,
        "midhot_bins": midhot_bins, "deepke": deepke,
        "deepke_bins": deepke_bins, "n1_ke": ke_mean.get(1, np.nan),
    }, ke_bins


def _g3_md_rung_tables():
    """The MD generator's exact rung tables (the leg-A'/S6 convention)."""
    from scripts.gen_tier2_md_confirmation import (
        floor1_rungs_eV,
        rq4graded_rungs_eV,
    )

    return {
        "rq4graded": sigma_cum(
            np.asarray(rq4graded_rungs_eV()[:N_STAR], dtype=float)
        ),
        "floor1": sigma_cum(
            np.asarray(floor1_rungs_eV()[:N_STAR], dtype=float)
        ),
    }


def _g3_s6_oracle(m, ref_ke, solv_exp):
    """Part A: re-derive the three committed h2b_s6_final rows bit-for-bit.

    Draw discipline is stage_legd's `d` leg verbatim (SEED 20260711; u/mu
    first, the analytic kornilov prior after), density_tied dressing,
    p = 1 onset coupling, exact tau rescale from TAU_PS.
    """
    from dataclasses import replace as _replace

    from i2_helium_md.sampling.droplet_sizes import sample_droplet_sizes_analytic

    sig_tables = _g3_md_rung_tables()

    rng = np.random.default_rng(SEED)
    u_ap = rng.uniform(0.0, 1.0, m)
    mu_ap = rng.uniform(-1.0, 1.0, m)
    cfg_d = _replace(
        _CFG,
        droplet_size_prior="kornilov_lognormal",
        use_single_droplet_size=False,
        num_molecules=m,
    )
    N_d = np.asarray(sample_droplet_sizes_analytic(cfg_d, rng=rng), dtype=float)
    R_d = np.asarray(droplet_radius_bulk_angstrom(N_d), dtype=float)
    r0 = (R_d - LEG_APRIME_MARGIN_A) * np.cbrt(u_ap)
    rho = rho_he_ratio(r0 - R_d, steepness=STEEP_A)
    ne_mol = np.clip(np.rint(N_STAR * rho), 0, N_STAR).astype(int)
    mass_mol = complex_mass_amu(ne_mol).astype(float)
    ne = np.concatenate([ne_mol, ne_mol])

    with open(OUT / "h2b_s6_final_predictions.csv", newline="") as fh:
        ref_rows = {
            (r["ladder"], r["v_c"], r["tau_ps"], r["E0_eV"]): r
            for r in csv.DictReader(fh)
        }
    with open(OUT / "h2b_s6_final_ke.csv", newline="") as fh:
        ref_ke_rows = list(csv.DictReader(fh))

    chord_cache = {}
    for ladder, v_c, p_tail, tau, e0 in G3_S6_CELLS:
        key = (v_c, p_tail)
        if key not in chord_cache:
            chord_cache[key] = integrate_pairs(
                r0, mu_ap, R_d, mass_mol,
                r0_sep=R0_SEP_PROD_A, drag_on=True, v_c=v_c, p_tail=p_tail,
            )
        res = chord_cache[key]
        K = res["K"].reshape(-1) * (TAU_PS / tau)
        n_det, sup = fate_map(ne, K, e0, 1, sig_tables[ladder])
        obs, ke_bins = g3_score(
            n_det, sup, res["trapped"].reshape(-1).astype(bool),
            res["v_inf"].reshape(-1), ref_ke, solv_exp,
        )
        row = {
            "ladder": ladder, "v_c": v_c, "tau_ps": tau, "E0_eV": e0,
            "p_tail": p_tail, "m": m,
            "trapped_frac": round(obs["trapped_frac"], 4),
            "suppressed_frac": round(obs["sup_frac"], 4),
            "nbar_det": round(obs["nbar"], 3),
            "n1_solv": round(obs["n1_solv"], 4),
            "ratio_n1_n2": round(obs["ratio"], 3),
            "w1_solv": round(obs["w1"], 4),
            "midhot_n2_8": round(obs["midhot_arith"], 4),
            "n1_ke_eV": round(obs["n1_ke"], 4),
        }
        row.update({f"h{k}": round(float(obs["hist"][k]), 4)
                    for k in range(N_STAR + 1)})
        cell_id = (ladder, str(v_c), str(tau), str(e0))
        ref = ref_rows[cell_id]
        for col, val in row.items():
            if str(val) != ref[col]:
                raise AssertionError(
                    f"G3-P1 oracle FAILED: s6 cell {cell_id} differs at "
                    f"{col}: {val!r} != {ref[col]!r}"
                )
        mine_ke = [(n, str(round(wgt, 4)), str(round(mke, 4)))
                   for n, wgt, mke in ke_bins]
        want_ke = [(int(r["n"]), r["weight"], r["mean_KE_eV"])
                   for r in ref_ke_rows
                   if (r["ladder"], r["v_c"], r["tau_ps"], r["E0_eV"])
                   == cell_id]
        if mine_ke != want_ke:
            raise AssertionError(
                f"G3-P1 oracle FAILED: KE rows of s6 cell {cell_id} differ "
                "from the committed record."
            )
        print(f"[g3 oracle] s6 {ladder} v{v_c} t{tau} e{e0}: bit-exact")
    print("G3-P1 wiring oracle PASSED: all three committed h2b_s6_final "
          "rows reproduced string-identically (predictions + KE).")
    # The standing twin row (v7.25) is Part C's baseline for the
    # corrected-vs-standing delta.
    lad, v_c, p_tail, tau, e0 = G3_STANDING
    res = chord_cache[(v_c, p_tail)]
    K = res["K"].reshape(-1) * (TAU_PS / tau)
    n_det, sup = fate_map(ne, K, e0, 1, sig_tables[lad])
    obs, _ = g3_score(
        n_det, sup, res["trapped"].reshape(-1).astype(bool),
        res["v_inf"].reshape(-1), ref_ke, solv_exp,
    )
    return obs


def _g3_grid_cell_draw(law, R, m, rng):
    """Per-cell birth draw for Part B. Draw order (documented, per cell,
    fresh rng at G3_SEED): positions first, axis cosines after.

    l1: the exact center-pin -- ONE chord (r0 = 0; the twin is deterministic
    there, so m collapses to 1 and the row has zero mechanism spread).
    l2: the repo Boltzmann sampler at the parent well 313.2 K.
    l3: uniform_volume, hard margin 3 A (exact inverse-CDF, the L3 law).
    """
    from dataclasses import replace as _replace

    if law == "l1":
        r0 = np.zeros(1)
    elif law == "l2":
        cfg_l2 = _replace(_CFG, binding_energy_molecule_K=G3_CORRECTED_WELL_K)
        r0 = sample_radial_positions(cfg_l2, np.full(m, R), rng=rng)
    elif law == "l3":
        r0 = (R - 3.0) * np.cbrt(rng.uniform(0.0, 1.0, m))
    else:
        raise ValueError(f"unknown grid law {law!r}")
    mu = rng.uniform(-1.0, 1.0, r0.size)
    return r0, mu


def _g3_grid_twin(m, ref_ke, solv_exp):
    """Part B: the 11 G1 geometry cells through the twin at the standing
    point; every twin observable is printed next to its MD value and the
    delta -- the re-issued twin-authority measurement at long chords."""
    lad, v_c, p_tail, tau, e0 = G3_STANDING
    sig = _g3_md_rung_tables()[lad]

    hist_rows, ke_rows = [], []
    for label, n_he, law in G3_GRID_CELLS:
        R = float(droplet_radius_bulk_angstrom(float(n_he)))
        rng = np.random.default_rng(G3_SEED)  # CRN: same seed every cell
        r0, mu = _g3_grid_cell_draw(law, R, m, rng)
        M = r0.size
        depth = R - r0
        rho = rho_he_ratio(r0 - R, steepness=STEEP_A)
        ne_mol = np.clip(np.rint(N_STAR * rho), 0, N_STAR).astype(int)
        res = integrate_pairs(
            r0, mu, np.full(M, R), complex_mass_amu(ne_mol).astype(float),
            r0_sep=R0_SEP_PROD_A, drag_on=True, v_c=v_c, p_tail=p_tail,
        )
        K = res["K"].reshape(-1) * (TAU_PS / tau)
        ne = np.concatenate([ne_mol, ne_mol])
        n_det, sup = fate_map(ne, K, e0, 1, sig)
        obs, ke_bins = g3_score(
            n_det, sup, res["trapped"].reshape(-1).astype(bool),
            res["v_inf"].reshape(-1), ref_ke, solv_exp,
            min_bin_count=(1 if law == "l1" else 20),
        )
        md = dict(zip(G3_MD_OBS_KEYS, G3_MD_G1_ROWS[label]))
        twin = {
            "trap": obs["trapped_frac"], "supp": obs["sup_frac"],
            "nbar": obs["nbar"], "n1_solv": obs["n1_solv"], "w1": obs["w1"],
            "midhot_geo": obs["midhot_geo"], "deepke": obs["deepke"],
        }
        row = {
            "cell": label, "law": law, "N_he": n_he, "R_A": round(R, 2),
            "m_chords": M,
            "depth_mean_A": round(float(depth.mean()), 2),
            "depth_sd_A": round(float(depth.std()), 2),
            "n_eject_mean": round(float(ne.mean()), 2),
            "K_q05": round(float(np.quantile(K, 0.05)), 4),
            "K_q50": round(float(np.quantile(K, 0.50)), 4),
            "K_q95": round(float(np.quantile(K, 0.95)), 4),
            "midhot_bins": obs["midhot_bins"],
            "deepke_bins": obs["deepke_bins"],
        }
        for key in G3_MD_OBS_KEYS:
            row[f"twin_{key}"] = round(twin[key], 4)
            row[f"md_{key}"] = md[key]
            row[f"d_{key}"] = round(twin[key] - md[key], 4)
        row.update({f"h{k}": round(float(obs["hist"][k]), 4)
                    for k in range(N_STAR + 1)})
        hist_rows.append(row)
        ke_rows.extend(
            {"cell": label, "n": n, "weight": round(wgt, 4),
             "mean_KE_eV": round(mke, 4)}
            for n, wgt, mke in ke_bins
        )
        print(f"[g3 grid {label}] twin trap={twin['trap']:.3f}"
              f" (MD {md['trap']:.3f})  nbar={twin['nbar']:.2f}"
              f" (MD {md['nbar']:.2f})  W1={twin['w1']:.2f}"
              f" (MD {md['w1']:.2f})  deepKE={twin['deepke']:.2f}"
              f" (MD {md['deepke']:.2f})")

    with open(OUT / "h2b_g3_grid_twin.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(hist_rows[0]))
        wtr.writeheader()
        wtr.writerows(hist_rows)
    with open(OUT / "h2b_g3_grid_twin_ke.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(ke_rows[0]))
        wtr.writeheader()
        wtr.writerows(ke_rows)
    print(f"grid twin  -> {OUT / 'h2b_g3_grid_twin.csv'}")
    print(f"grid KE    -> {OUT / 'h2b_g3_grid_twin_ke.csv'}")
    return hist_rows


def _g3_corrected_master(m, force_rebuild=False):
    """Part C master draw: the G2-adopted corrected geometry.

    ``legacy`` + ``raw`` sizes (nozzle correlation at the preset's own
    p = 40 mbar / T = 14 K => <N> = 12794), Boltzmann births at the
    313.2 K parent well. Draw order (documented, seed G3_SEED): N raw
    ln-normal -> Boltzmann birth radii -> axis cosines. Cached to npz
    (the Boltzmann sampler loops per unique radius -- minutes at m = 20000).
    """
    from dataclasses import replace as _replace

    from i2_helium_md.sampling.droplet_sizes import sample_droplet_sizes

    cache = OUT / "h2b_g3_corrected_master.npz"
    if cache.exists() and not force_rebuild:
        return dict(np.load(cache))
    cfg_c = _replace(
        _CFG,
        num_molecules=m,
        use_single_droplet_size=False,
        binding_energy_molecule_K=G3_CORRECTED_WELL_K,
    )
    rng = np.random.default_rng(G3_SEED)
    N = np.asarray(sample_droplet_sizes(cfg_c, mode="raw", rng=rng),
                   dtype=float)
    R = np.asarray(droplet_radius_bulk_angstrom(N), dtype=float)
    r0 = np.asarray(sample_radial_positions(cfg_c, R, rng=rng), dtype=float)
    mu = rng.uniform(-1.0, 1.0, m)
    out = {"N": N, "R": R, "r0": r0, "mu": mu}
    np.savez_compressed(cache, **out)
    return out


def _g3_corrected_ensemble(m, ref_ke, solv_exp):
    """The corrected-ensemble computation at the standing point (Step-1
    Part C (ii), factored value-identically for reuse): committed master
    -> standing chord -> fate map -> S6 scoring.

    Returns a dict with ``row`` (the exact ``h2b_g3_corrected_row.csv``
    schema), ``ke_bins``, ``obs``, and the raw per-fragment arrays
    (``res``/``K655``/``ne``/``trapped`` etc.) so stage_g3scan's oracle and
    zero-integration pre-scan read the same standing chord (§1.4).
    """
    lad, v_c, p_tail, tau, e0 = G3_STANDING
    sig = _g3_md_rung_tables()[lad]
    ms = _g3_corrected_master(m)
    N, R, r0, mu = ms["N"], ms["R"], ms["r0"], ms["mu"]
    depth = R - r0
    rho = rho_he_ratio(r0 - R, steepness=STEEP_A)
    ne_mol = np.clip(np.rint(N_STAR * rho), 0, N_STAR).astype(int)
    res = integrate_pairs(
        r0, mu, R, complex_mass_amu(ne_mol).astype(float),
        r0_sep=R0_SEP_PROD_A, drag_on=True, v_c=v_c, p_tail=p_tail,
    )
    K = res["K"].reshape(-1) * (TAU_PS / tau)
    ne = np.concatenate([ne_mol, ne_mol])
    n_det, sup = fate_map(ne, K, e0, 1, sig)
    obs, ke_bins = g3_score(
        n_det, sup, res["trapped"].reshape(-1).astype(bool),
        res["v_inf"].reshape(-1), ref_ke, solv_exp,
    )
    t_exit = res["t_exit"].reshape(-1)
    row = {
        "geometry": "corrected", "ladder": lad, "v_c": v_c, "tau_ps": tau,
        "E0_eV": e0, "p_tail": p_tail, "m": m,
        "N_q05": round(float(np.quantile(N, 0.05)), 1),
        "N_q50": round(float(np.quantile(N, 0.50)), 1),
        "N_q95": round(float(np.quantile(N, 0.95)), 1),
        "N_mean": round(float(N.mean()), 1),
        "R_q05": round(float(np.quantile(R, 0.05)), 2),
        "R_q50": round(float(np.quantile(R, 0.50)), 2),
        "R_q95": round(float(np.quantile(R, 0.95)), 2),
        "depth_q05": round(float(np.quantile(depth, 0.05)), 2),
        "depth_q50": round(float(np.quantile(depth, 0.50)), 2),
        "depth_q95": round(float(np.quantile(depth, 0.95)), 2),
        "n_eject_mean": round(float(ne.mean()), 3),
        "K_q05": round(float(np.quantile(K, 0.05)), 4),
        "K_q50": round(float(np.quantile(K, 0.50)), 4),
        "K_q95": round(float(np.quantile(K, 0.95)), 4),
        "t_exit_q50": round(float(np.nanquantile(t_exit, 0.50)), 2),
        "t_exit_q95": round(float(np.nanquantile(t_exit, 0.95)), 2),
        "trapped_frac": round(obs["trapped_frac"], 4),
        "suppressed_frac": round(obs["sup_frac"], 4),
        "nbar_det": round(obs["nbar"], 3),
        "n1_solv": round(obs["n1_solv"], 4),
        "ratio_n1_n2": round(obs["ratio"], 3),
        "w1_solv": round(obs["w1"], 4),
        "midhot_n2_8": round(obs["midhot_arith"], 4),
        "midhot_geo": round(obs["midhot_geo"], 4),
        "midhot_bins": obs["midhot_bins"],
        "deepke": round(obs["deepke"], 4),
        "deepke_bins": obs["deepke_bins"],
        "n1_ke_eV": round(obs["n1_ke"], 4),
    }
    row.update({f"h{k}": round(float(obs["hist"][k]), 4)
                for k in range(N_STAR + 1)})
    return {
        "row": row, "ke_bins": ke_bins, "obs": obs, "sig": sig,
        "N": N, "R": R, "r0": r0, "mu": mu, "depth": depth,
        "ne_mol": ne_mol, "ne": ne, "res": res,
        "K655": res["K"].reshape(-1),
        "trapped": res["trapped"].reshape(-1).astype(bool),
        "v_inf": res["v_inf"].reshape(-1),
    }


def _g3_corrected(m, ref_ke, solv_exp, standing_obs):
    """Part C: landmark re-issue at the corrected geometry.

    (i) center-pin landmark rows per grid radius + the R(2000) continuity
    anchor (must reproduce the recorded K = 0.74460); (ii) the corrected-
    ensemble row at the standing point, printed against the standing twin
    row (Part A) and the D2b §4.3 corr x L2 re-weighting forecast.
    """
    lad, v_c, p_tail, tau, e0 = G3_STANDING
    sig = _g3_md_rung_tables()[lad]

    # ---- (i) center-pin landmarks (single chord, n_eject = 21 dressing).
    # Two laws per pin: the PURE-CUBIC column is the continuity family (the
    # recorded K = 0.74460 landmark was measured under it -- stage_oracles
    # O2 runs v_c=None), the capped-tail columns are the operative standing
    # law (finc1v725) being re-issued.
    pin_rows = []
    for tag, n_he in (("R2000_anchor", 2000.0), ("r1", 1727.0),
                      ("r2", 3605.0), ("r3", 11059.0),
                      ("corr_meanN", 12794.0), ("r4", 29227.0)):
        R = float(droplet_radius_bulk_angstrom(n_he))
        args = (
            np.array([0.0]), np.array([1.0]), np.array([R]),
            np.array([float(complex_mass_amu(N_STAR))]),
        )
        res = integrate_pairs(
            *args, r0_sep=R0_SEP_PROD_A, drag_on=True, v_c=v_c, p_tail=p_tail,
        )
        res_pc = integrate_pairs(*args, r0_sep=R0_SEP_PROD_A, drag_on=True)
        K655 = float(res["K"][0, 0])
        K655_pc = float(res_pc["K"][0, 0])
        n_det, _ = fate_map(
            np.array([N_STAR]), np.array([K655 * (TAU_PS / tau)]), e0, 1, sig
        )
        pin_rows.append({
            "tag": tag, "N_he": int(n_he), "R_A": round(R, 2),
            "K_tau655_cubic": round(K655_pc, 5),
            "exposure_ps": round(K655 * TAU_PS, 4),
            "K_tau655": round(K655, 5),
            "K_tau32": round(K655 * (TAU_PS / tau), 5),
            "t_exit_ps": round(float(res["t_exit"][0, 0]), 2),
            "v_peak_Aps": round(float(res["v_peak"][0, 0]), 2),
            "v_inf_Aps": round(float(res["v_inf"][0, 0]), 2),
            "n_det": int(n_det[0]),
            "KE_det_eV": round(float(kinetic_energy_eV(
                complex_mass_amu(n_det[0]), res["v_inf"][0, 0])), 4),
        })
        print(f"[g3 pin {tag:12s}] R={R:6.2f}  K_cubic={K655_pc:.5f}  "
              f"K_capped={K655:.5f}  t_exit={res['t_exit'][0, 0]:.2f}  "
              f"v_inf={res['v_inf'][0, 0]:.2f}  n_det={n_det[0]}")
    anchor = next(r for r in pin_rows if r["tag"] == "R2000_anchor")
    if abs(anchor["K_tau655_cubic"] - 0.74460) > 1.5e-5:
        raise AssertionError(
            f"G3-P2 continuity anchor FAILED: pure-cubic K(R2000, tau 6.55) "
            f"= {anchor['K_tau655_cubic']} != recorded 0.74460"
        )
    print("G3-P2 continuity anchor PASSED: the recorded pure-cubic "
          "center-pin landmark K = 0.74460 is reproduced.")

    # ---- (ii) the corrected-ensemble row at the standing point
    ens = _g3_corrected_ensemble(m, ref_ke, solv_exp)
    row, ke_bins, obs, N = ens["row"], ens["ke_bins"], ens["obs"], ens["N"]

    with open(OUT / "h2b_g3_corrected_landmarks.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(pin_rows[0]))
        wtr.writeheader()
        wtr.writerows(pin_rows)
    with open(OUT / "h2b_g3_corrected_row.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(row))
        wtr.writeheader()
        wtr.writerows([row])
    with open(OUT / "h2b_g3_corrected_ke.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=["n", "weight", "mean_KE_eV"])
        wtr.writeheader()
        wtr.writerows(
            {"n": n, "weight": round(wgt, 4), "mean_KE_eV": round(mke, 4)}
            for n, wgt, mke in ke_bins
        )

    print(f"\n[g3 corrected ensemble] drawn <N> = {N.mean():.0f} "
          "(nozzle correlation 12794)")
    print("               standing(twin)  corrected(twin)   D2b corr x L2 "
          "forecast")
    fc = {"trap": "0.31-0.42", "supp": "0", "nbar": "13.9-14.7",
          "n1_solv": "0", "w1": "9.0-9.8", "midhot_geo": "1.256",
          "deepke": "1.80-1.90"}
    twin_now = {
        "trap": obs["trapped_frac"], "supp": obs["sup_frac"],
        "nbar": obs["nbar"], "n1_solv": obs["n1_solv"], "w1": obs["w1"],
        "midhot_geo": obs["midhot_geo"], "deepke": obs["deepke"],
    }
    std = {
        "trap": standing_obs["trapped_frac"], "supp": standing_obs["sup_frac"],
        "nbar": standing_obs["nbar"], "n1_solv": standing_obs["n1_solv"],
        "w1": standing_obs["w1"], "midhot_geo": standing_obs["midhot_geo"],
        "deepke": standing_obs["deepke"],
    }
    for key in G3_MD_OBS_KEYS:
        print(f"  {key:11s} {std[key]:14.4f}  {twin_now[key]:15.4f}   "
              f"{fc[key]}")
    print(f"landmarks  -> {OUT / 'h2b_g3_corrected_landmarks.csv'}")
    print(f"ensemble   -> {OUT / 'h2b_g3_corrected_row.csv'}")


def stage_g3landmarks(m=20000):
    """G3 Step 1: twin landmark re-issue at the corrected geometry.

    Part A -- G3-P1 wiring oracle: the three committed ``h2b_s6_final``
    rows re-derived string-identically before any new number is read
    (atlas §1.4 oracle discipline). Part B -- the 11 Axis-A G1 cells
    through the twin at the standing point, each printed against its MD
    row: the re-issued twin<->MD transfer measurement at long chords
    (channel-(d) authority box). Part C -- center-pin landmark rows per
    grid radius (with the K = 0.74460 continuity anchor) and the
    corrected-ensemble twin row vs the D2b §4.3 forecast. Zero MD.
    """
    OUT.mkdir(parents=True, exist_ok=True)
    ref_ke = g3_ref_mean_ke()
    solv_exp, _ = load_experiment()
    print("=== G3 Step 1: Part A (S6 wiring oracle) ===")
    standing_obs = _g3_s6_oracle(m, ref_ke, solv_exp)
    print("\n=== G3 Step 1: Part B (G1 grid twin re-issue) ===")
    _g3_grid_twin(m, ref_ke, solv_exp)
    print("\n=== G3 Step 1: Part C (corrected-geometry landmarks) ===")
    _g3_corrected(m, ref_ke, solv_exp, standing_obs)


# --------------------------------------------------------------------------
# G3 Step 2 — the Route A/B twin scan at the corrected geometry
# (plan §3.5c, designed + user-approved 2026-07-27; executed under its own
# `[PROCEED TO IMPLEMENTATION]`)
# --------------------------------------------------------------------------
# Chord surface (Route B): every (v_c, E_bind) pair is one integration of the
# committed corrected master (G3_SEED, m = 20000). The 5.0 floor respects the
# TDDFT band top (4.95 A/ps); 3.5/4.25 are outside-Tier-0-authority
# DIAGNOSTIC arms (they override the calibrated band — if only these land,
# the finding points at the collaborator ask, not at a parameter point).
G3SCAN_VC_TIER0 = (5.0, 5.5, 6.0, 6.5, 7.25, 8.0, 9.0, 10.0)
G3SCAN_VC_DIAG = (3.5, 4.25)
# E_bind: the Tier-0 extracted spread (findings stage-2a form table; the
# §6.5/§6.7 joint-pairing exception — off-bundle wells run against the
# bundle-b chord law, stamped per row). Exact extracted values, tagged with
# the §6.7 item-2 labels: 0.0482 = lq co-extracted, bundle stamp = standing,
# 0.154 = 9 A per-case cubic.
G3SCAN_EBIND = (
    ("eb0482", 0.0482),
    ("eb1168", E_BIND_ION_EV),
    ("eb154", 0.154),
)
# Free surface (Route A): exact tau rescale + fate-map re-run — pure
# post-processing of a stored chord (the Step-1 structural fact).
G3SCAN_TAU_PS = (2.4, 3.2, 4.8, 6.4, 9.6, 12.8)
G3SCAN_TAU_SOURCED = 6.55  # a landing needing tau >> this is flagged
G3SCAN_E0_GRID = np.round(np.arange(0.17, 0.5201, 0.01), 3)  # 36 values
# Pre-registered hard gate (twin level, item-3 authority box applied):
# n1_solv in [0.19, 0.30] AND nbar in [4.4, 7.1] (target 4.07 + the
# residence-conditional bias bracket [+0.3, +3]). W1/KE/supp non-gating.
G3SCAN_GATE_N1SOLV = (0.19, 0.30)
G3SCAN_GATE_NBAR = (4.4, 7.1)
# Pre-scan (§3.5c block 1): exposure scalings X -> f*X of the stored
# standing chord; the analytic Route-A kill criterion reads f = 0.25 (the
# maximal reduction the chord axes could plausibly deliver).
G3SCAN_PRESCAN_F = np.round(np.arange(0.05, 1.0001, 0.05), 2)
G3SCAN_PRESCAN_FMIN = 0.25


def g3scan_gate(n1_solv, nbar):
    """The §3.5c pre-registered hard gate; NaN observables never gate."""
    return bool(
        G3SCAN_GATE_N1SOLV[0] <= n1_solv <= G3SCAN_GATE_N1SOLV[1]
        and G3SCAN_GATE_NBAR[0] <= nbar <= G3SCAN_GATE_NBAR[1]
    )


def _g3scan_light(n_det, trapped):
    """Light score for the pre-scan: (n1_solv, nbar, frac_le8) over the
    non-trapped read (no KE, no W1). n1_solv is NaN if nothing solvates."""
    w = np.where(trapped, 0.0, 1.0)
    hist = np.bincount(n_det, weights=w, minlength=N_STAR + 1) / w.sum()
    nbar = float((np.arange(N_STAR + 1) * hist).sum())
    solv_tot = hist[1:].sum()
    n1 = float(hist[1] / solv_tot) if solv_tot > 0 else np.nan
    frac_le8 = float(hist[: 9].sum())
    return n1, nbar, frac_le8


def _g3scan_prescan(K655, trapped, ne, sig):
    """§3.5c block 1 — the required-exposure-reduction map (zero
    integration). For every (tau, E0) on the free surface, scale the stored
    standing chord's exposure X -> f*X and record where the hard gate holds;
    the analytic Route-A kill criterion fires iff no f >= 0.25 produces
    n1_solv inside the gate band at any (tau, E0).

    Returns (rows, verdict) — verdict holds the kill flag + counts.
    """
    rows = []
    n1_ok_any = 0
    gate_ok_any = 0
    for tau in G3SCAN_TAU_PS:
        for e0 in G3SCAN_E0_GRID:
            n1_by_f, nbar_by_f, gate_by_f = {}, {}, {}
            for f in G3SCAN_PRESCAN_F:
                K = K655 * (f * (TAU_PS / tau))
                n_det, _ = fate_map(ne, K, e0, 1, sig)
                n1, nbar, frac_le8 = _g3scan_light(n_det, trapped)
                n1_by_f[f] = n1
                nbar_by_f[f] = nbar
                gate_by_f[f] = g3scan_gate(n1, nbar)
                if f == 1.0:
                    strip_le8_f100 = frac_le8
            f_gate = [f for f, ok in gate_by_f.items() if ok]
            f_admiss = [f for f in G3SCAN_PRESCAN_F
                        if f >= G3SCAN_PRESCAN_FMIN - 1e-9]
            n1_ok = any(
                G3SCAN_GATE_N1SOLV[0] <= n1_by_f[f] <= G3SCAN_GATE_N1SOLV[1]
                for f in f_admiss
            )
            gate_ok = any(gate_by_f[f] for f in f_admiss)
            n1_ok_any += int(n1_ok)
            gate_ok_any += int(gate_ok)
            rows.append({
                "tau_ps": tau,
                "tau_flag": int(tau > G3SCAN_TAU_SOURCED),
                "E0_eV": e0,
                "n1_solv_f100": round(n1_by_f[1.0], 4),
                "nbar_f100": round(nbar_by_f[1.0], 3),
                "strip_le8_f100": round(strip_le8_f100, 4),
                "gate_f100": int(gate_by_f[1.0]),
                "f_gate_lo": min(f_gate) if f_gate else np.nan,
                "f_gate_hi": max(f_gate) if f_gate else np.nan,
                "n1_solv_f025": round(n1_by_f[0.25], 4),
                "nbar_f025": round(nbar_by_f[0.25], 3),
                "gate_f025": int(gate_by_f[0.25]),
                "n1_ok_fge025": int(n1_ok),
                "gate_ok_fge025": int(gate_ok),
            })
    verdict = {
        "route_a_killed": n1_ok_any == 0,
        "cells_n1_ok_fge025": n1_ok_any,
        "cells_gate_ok_fge025": gate_ok_any,
        "cells_total": len(rows),
    }
    return rows, verdict


def _g3scan_chord_tag(v_c, eb_tag):
    return f"vc{str(v_c).replace('.', 'p')}_{eb_tag}"


def _g3scan_chord_family(v_c, eb_tag, e_bind_ev, m, ens, force_rebuild=False):
    """One chord-surface integration of the corrected master, npz-cached.

    The standing family (v_c 7.25, bundle well) is served from the oracle's
    in-memory chord (``ens``) — never re-integrated, never cached separately.
    Cache staleness is guarded by the stamped (v_c, e_bind, m, seed).
    """
    lad, vc_std, p_tail, _, _ = G3_STANDING
    if v_c == vc_std and eb_tag == "eb1168":
        return {"K": ens["res"]["K"], "v_inf": ens["res"]["v_inf"],
                "t_exit": ens["res"]["t_exit"],
                "trapped": ens["res"]["trapped"]}
    cache = OUT / f"h2b_g3scan_chord_{_g3scan_chord_tag(v_c, eb_tag)}.npz"
    if cache.exists() and not force_rebuild:
        dat = dict(np.load(cache))
        stamps = (float(dat["v_c"]), float(dat["e_bind_ev"]),
                  int(dat["m"]), int(dat["seed"]))
        if stamps != (float(v_c), float(e_bind_ev), int(m), G3_SEED):
            raise AssertionError(
                f"stale g3scan chord cache {cache.name}: stamps {stamps}"
            )
        return dat
    ms = _g3_corrected_master(m)
    rho = rho_he_ratio(ms["r0"] - ms["R"], steepness=STEEP_A)
    ne_mol = np.clip(np.rint(N_STAR * rho), 0, N_STAR).astype(int)
    res = integrate_pairs(
        ms["r0"], ms["mu"], ms["R"], complex_mass_amu(ne_mol).astype(float),
        r0_sep=R0_SEP_PROD_A, drag_on=True, v_c=v_c, p_tail=p_tail,
        e_bind_ev=e_bind_ev,
    )
    out = {"K": res["K"], "v_inf": res["v_inf"], "t_exit": res["t_exit"],
           "trapped": res["trapped"],
           "v_c": float(v_c), "e_bind_ev": float(e_bind_ev),
           "m": float(m), "seed": float(G3_SEED)}
    np.savez_compressed(cache, **out)
    return out


def stage_g3scan(m=20000, force_rebuild=False):
    """G3 Step 2: the nested Route A/B twin scan at the corrected geometry
    (plan §3.5c; zero MD).

    Block order (frozen): 0. oracles — the S6 machinery oracle (G3-P1) and
    the committed ``h2b_g3_corrected_row{,_ke}.csv`` re-derived bit-exact at
    the standing cell (the Step-1 landmark-continuity convention); 1. the
    zero-integration pre-scan with the analytic Route-A kill criterion;
    2. the chord families (v_c x E_bind, npz-cached), each scored over the
    full free surface (tau x E0) as it lands.

    Fixed: the committed corrected master (G3_SEED, m = 20000), rq4graded
    ladder, p = 1, p_tail = -1, detection = the non-trapped chord read.
    Gate/failure criterion: pre-registered in §3.5c (hard gate n1_solv AND
    nbar; W1/KE/supp non-gating; trap reported as a floor).
    """
    import time as _time

    OUT.mkdir(parents=True, exist_ok=True)
    ref_ke = g3_ref_mean_ke()
    solv_exp, _ = load_experiment()
    lad, _, p_tail, _, _ = G3_STANDING

    # ---- Block 0a: the S6 machinery oracle (standing geometry, bit-exact).
    print("=== G3 Step 2: Block 0 (oracles) ===")
    _g3_s6_oracle(m, ref_ke, solv_exp)

    # ---- Block 0b: the corrected-row landmark oracle (bit-exact).
    ens = _g3_corrected_ensemble(m, ref_ke, solv_exp)
    with open(OUT / "h2b_g3_corrected_row.csv", newline="") as fh:
        ref_rows = list(csv.DictReader(fh))
    if len(ref_rows) != 1:
        raise AssertionError("h2b_g3_corrected_row.csv must hold one row")
    ref_row = ref_rows[0]
    if set(ref_row) != set(map(str, ens["row"])):
        raise AssertionError(
            "g3scan corrected-row oracle FAILED: column set differs"
        )
    for col, val in ens["row"].items():
        if str(val) != ref_row[col]:
            raise AssertionError(
                f"g3scan corrected-row oracle FAILED at {col}: "
                f"{val!r} != {ref_row[col]!r}"
            )
    with open(OUT / "h2b_g3_corrected_ke.csv", newline="") as fh:
        ref_ke_rows = [(int(r["n"]), r["weight"], r["mean_KE_eV"])
                       for r in csv.DictReader(fh)]
    mine_ke = [(n, str(round(w, 4)), str(round(k, 4)))
               for n, w, k in ens["ke_bins"]]
    if mine_ke != ref_ke_rows:
        raise AssertionError(
            "g3scan corrected-row oracle FAILED: KE rows differ from the "
            "committed h2b_g3_corrected_ke.csv"
        )
    print("g3scan landmark oracle PASSED: h2b_g3_corrected_row.csv + _ke.csv "
          "re-derived string-identically at the standing cell.")

    sig = ens["sig"]
    ne = ens["ne"]
    R_frag = np.concatenate([ens["R"], ens["R"]])

    # ---- Block 1: the zero-integration pre-scan (Route-A kill criterion).
    print("\n=== G3 Step 2: Block 1 (zero-integration pre-scan) ===")
    pre_rows, pre_verdict = _g3scan_prescan(
        ens["K655"], ens["trapped"], ne, sig
    )
    with open(OUT / "h2b_g3scan_prescan.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(pre_rows[0]))
        wtr.writeheader()
        wtr.writerows(pre_rows)
    print(f"pre-scan   -> {OUT / 'h2b_g3scan_prescan.csv'}")
    print(f"[pre-scan] (tau, E0) cells with n1_solv in gate band at some "
          f"f >= {G3SCAN_PRESCAN_FMIN}: {pre_verdict['cells_n1_ok_fge025']}"
          f"/{pre_verdict['cells_total']}"
          f" (full hard gate: {pre_verdict['cells_gate_ok_fge025']})")
    if pre_verdict["route_a_killed"]:
        print("[pre-scan] ANALYTIC ROUTE-A KILL CRITERION FIRED: X -> 0.25*X "
              "cannot produce n1_solv in the gate band at any (tau, E0). "
              "Route A is dead at this surface (§3.5c block 1); the chord "
              "scan proceeds — Route B carries the remaining freedom.")
    else:
        print("[pre-scan] Route-A kill criterion NOT fired: exposure "
              "reduction inside the credible chord-axis range can reach the "
              "n1_solv band.")

    # ---- Block 2: chord families + the full nested scoring.
    print("\n=== G3 Step 2: Block 2 (chord families + nested scoring) ===")
    eb_by_tag = dict(G3SCAN_EBIND)
    chord_rows, scan_rows = [], []
    gated_hist_rows, gated_ke_rows = [], []
    n_gated = 0
    t0 = _time.time()
    for v_c in G3SCAN_VC_TIER0 + G3SCAN_VC_DIAG:
        authority = "tier0" if v_c in G3SCAN_VC_TIER0 else "diagnostic"
        for eb_tag, e_bind in G3SCAN_EBIND:
            fam = _g3scan_chord_family(
                v_c, eb_tag, e_bind, m, ens, force_rebuild=force_rebuild
            )
            K655 = np.asarray(fam["K"]).reshape(-1)
            trapped = np.asarray(fam["trapped"]).reshape(-1).astype(bool)
            v_inf = np.asarray(fam["v_inf"]).reshape(-1)
            t_exit = np.asarray(fam["t_exit"]).reshape(-1)
            det = ~trapped
            det_yield = float(det.mean())
            if det.sum() == 0:
                raise AssertionError(
                    f"chord family {_g3scan_chord_tag(v_c, eb_tag)}: every "
                    "fragment trapped — the non-trapped read is empty"
                )
            chord_rows.append({
                "v_c": v_c, "p_tail": p_tail, "E_bind_tag": eb_tag,
                "E_bind_eV": e_bind, "authority": authority,
                "ebind_exception": int(eb_tag != "eb1168"), "m": m,
                "trapped_frac": round(float(trapped.mean()), 4),
                "det_yield": round(det_yield, 4),
                "K655_q05": round(float(np.quantile(K655, 0.05)), 4),
                "K655_q50": round(float(np.quantile(K655, 0.50)), 4),
                "K655_q95": round(float(np.quantile(K655, 0.95)), 4),
                "v_inf_det_q50": round(float(np.quantile(v_inf[det], 0.50)), 3),
                "t_exit_q50": round(float(np.nanquantile(t_exit, 0.50)), 2),
                "R_det_q05": round(float(np.quantile(R_frag[det], 0.05)), 2),
                "R_det_q50": round(float(np.quantile(R_frag[det], 0.50)), 2),
                "R_det_q95": round(float(np.quantile(R_frag[det], 0.95)), 2),
                "R_src_q05": round(float(np.quantile(R_frag, 0.05)), 2),
                "R_src_q50": round(float(np.quantile(R_frag, 0.50)), 2),
                "R_src_q95": round(float(np.quantile(R_frag, 0.95)), 2),
            })
            fam_gated = 0
            for tau in G3SCAN_TAU_PS:
                K = K655 * (TAU_PS / tau)
                for e0 in G3SCAN_E0_GRID:
                    n_det, sup = fate_map(ne, K, e0, 1, sig)
                    row = {
                        "leg": "g3scan", "ladder": lad, "v_c": v_c,
                        "p_tail": p_tail, "E_bind_tag": eb_tag,
                        "E_bind_eV": e_bind, "authority": authority,
                        "ebind_exception": int(eb_tag != "eb1168"),
                        "tau_ps": tau,
                        "tau_flag": int(tau > G3SCAN_TAU_SOURCED),
                        "E0_eV": e0, "m": m,
                        "det_yield": round(det_yield, 4),
                    }
                    try:
                        obs, ke_bins = g3_score(
                            n_det, sup, trapped, v_inf, ref_ke, solv_exp
                        )
                    except ValueError:
                        # a legitimate scan outcome (e.g. nothing solvates):
                        # record the cell as un-scoreable, never gated
                        row.update({k: np.nan for k in (
                            "trapped_frac", "suppressed_frac", "nbar_det",
                            "n1_solv", "ratio_n1_n2", "w1_solv",
                            "midhot_arith", "midhot_geo", "deepke",
                            "n1_ke_eV")})
                        row.update({"midhot_bins": 0, "deepke_bins": 0,
                                    "gate_n1": 0, "gate_nbar": 0, "gate": 0})
                        scan_rows.append(row)
                        continue
                    g_n1 = (G3SCAN_GATE_N1SOLV[0] <= obs["n1_solv"]
                            <= G3SCAN_GATE_N1SOLV[1])
                    g_nb = (G3SCAN_GATE_NBAR[0] <= obs["nbar"]
                            <= G3SCAN_GATE_NBAR[1])
                    gate = bool(g_n1 and g_nb)
                    row.update({
                        "trapped_frac": round(obs["trapped_frac"], 4),
                        "suppressed_frac": round(obs["sup_frac"], 4),
                        "nbar_det": round(obs["nbar"], 3),
                        "n1_solv": round(obs["n1_solv"], 4),
                        "ratio_n1_n2": round(obs["ratio"], 3),
                        "w1_solv": round(obs["w1"], 4),
                        "midhot_arith": round(obs["midhot_arith"], 4),
                        "midhot_geo": round(obs["midhot_geo"], 4),
                        "midhot_bins": obs["midhot_bins"],
                        "deepke": round(obs["deepke"], 4),
                        "deepke_bins": obs["deepke_bins"],
                        "n1_ke_eV": round(obs["n1_ke"], 4),
                        "gate_n1": int(g_n1), "gate_nbar": int(g_nb),
                        "gate": int(gate),
                    })
                    scan_rows.append(row)
                    if gate:
                        fam_gated += 1
                        n_gated += 1
                        hrow = {k: row[k] for k in (
                            "v_c", "E_bind_tag", "tau_ps", "E0_eV")}
                        hrow.update({f"h{k}": round(float(obs["hist"][k]), 4)
                                     for k in range(N_STAR + 1)})
                        gated_hist_rows.append(hrow)
                        gated_ke_rows.extend(
                            {"v_c": v_c, "E_bind_tag": eb_tag, "tau_ps": tau,
                             "E0_eV": e0, "n": n, "weight": round(w, 4),
                             "mean_KE_eV": round(k_, 4)}
                            for n, w, k_ in ke_bins
                        )
            print(f"[g3scan {_g3scan_chord_tag(v_c, eb_tag):14s}] "
                  f"({authority:10s}) trap={trapped.mean():.3f} "
                  f"K655_q50={np.quantile(K655, 0.50):.3f} "
                  f"gated={fam_gated}/216  "
                  f"[{_time.time() - t0:6.0f} s]")

    with open(OUT / "h2b_g3scan_chords.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(chord_rows[0]))
        wtr.writeheader()
        wtr.writerows(chord_rows)
    with open(OUT / "h2b_g3scan_predictions.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(scan_rows[0]))
        wtr.writeheader()
        wtr.writerows(scan_rows)
    gated_hist_fields = (["v_c", "E_bind_tag", "tau_ps", "E0_eV"]
                         + [f"h{k}" for k in range(N_STAR + 1)])
    with open(OUT / "h2b_g3scan_gated_predictions.csv", "w",
              newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=gated_hist_fields)
        wtr.writeheader()
        wtr.writerows(gated_hist_rows)
    with open(OUT / "h2b_g3scan_gated_ke.csv", "w", newline="") as fh:
        wtr = csv.DictWriter(
            fh, fieldnames=["v_c", "E_bind_tag", "tau_ps", "E0_eV", "n",
                            "weight", "mean_KE_eV"]
        )
        wtr.writeheader()
        wtr.writerows(gated_ke_rows)
    print(f"\nchords     -> {OUT / 'h2b_g3scan_chords.csv'}")
    print(f"scan       -> {OUT / 'h2b_g3scan_predictions.csv'}")
    print(f"gated hist -> {OUT / 'h2b_g3scan_gated_predictions.csv'}")
    print(f"gated KE   -> {OUT / 'h2b_g3scan_gated_ke.csv'}")

    # ---- Verdict block (pre-registered gate / failure criterion).
    gated = [r for r in scan_rows if r["gate"] == 1]
    print(f"\n=== G3 Step 2 verdict: {len(gated)}/{len(scan_rows)} cells "
          "inside the hard gate ===")
    if gated:
        by_class = {}
        for r in gated:
            key = (r["authority"], "tau_flagged" if r["tau_flag"]
                   else "tau_sourced")
            by_class[key] = by_class.get(key, 0) + 1
        for key, cnt in sorted(by_class.items()):
            print(f"  {key[0]:10s} / {key[1]:11s}: {cnt}")
        best = sorted(gated, key=lambda r: r["w1_solv"])[:12]
        print("  gated cells by W1 (reported, NOT gating):")
        for r in best:
            print(f"    vc={r['v_c']:<5} {r['E_bind_tag']:6s} "
                  f"tau={r['tau_ps']:<4} E0={r['E0_eV']:<5} "
                  f"n1={r['n1_solv']:.3f} nbar={r['nbar_det']:.2f} "
                  f"W1={r['w1_solv']:.3f} trap={r['trapped_frac']:.3f} "
                  f"supp={r['suppressed_frac']:.3f} "
                  f"deepKE={r['deepke']:.2f} [{r['authority']}"
                  f"{', TAU-FLAG' if r['tau_flag'] else ''}"
                  f"{', EBIND-EXC' if r['ebind_exception'] else ''}]")
    else:
        print("  PRE-REGISTERED FAILURE CRITERION FIRED: no cell inside the "
              "gate anywhere on the grid INCLUDING the diagnostic arms.")
        print("  => Route A and Route B fail at the (v_c, tau, E0, E_bind) "
              "surface — the corrected-geometry failure localizes to the")
        print("     mechanism (pickup re-filling, per-shed eps, ladder "
              "shape — the knobs the twin does not carry), and the standing")
        print("     point reads as an effective model of the detected "
              "subset (§3.5c).")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "oracles"
    OUT.mkdir(parents=True, exist_ok=True)
    if mode == "oracles":
        stage_oracles()
    elif mode == "w12pred":
        stage_w12pred()
    elif mode == "birthlaw":
        stage_birthlaw()
    elif mode == "legaprime":
        stage_legaprime()
    elif mode == "legb":
        stage_legb()
    elif mode == "legc":
        stage_legc()
    elif mode == "legd":
        stage_legd()
    elif mode == "repilots1":
        stage_repilot1()
    elif mode == "repilots2":
        stage_repilot2()
    elif mode == "g3landmarks":
        stage_g3landmarks()
    elif mode == "g3scan":
        stage_g3scan()
    elif mode == "levers":
        tab = build_fragment_table()
        stage_levers(tab)
    elif mode == "scan":
        tab = build_fragment_table()
        stage_scan(tab)
    elif mode == "report":
        tab = build_fragment_table()
        stage_report(tab)
    else:
        raise SystemExit(
            f"unknown mode {mode!r} "
            "(oracles | levers | scan | report | w12pred | birthlaw | "
            "legaprime | legb | legc | legd | repilots1 | repilots2 | "
            "g3landmarks | g3scan)"
        )


if __name__ == "__main__":
    main()
