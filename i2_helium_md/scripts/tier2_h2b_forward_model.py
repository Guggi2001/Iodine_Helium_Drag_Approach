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
        w12pred | birthlaw | legaprime | legb | legc
Outputs: CSVs + text summaries in OUT (see USER SETTINGS).
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
):
    """Integrate fragment pairs born at radius r0 (cos(pos,axis) = mu).

    Fragment A moves along +axis (cosine mu), B along -axis (cosine -mu).
    Forces per fragment: Coulomb (pair, along axis), droplet solvation
    (projected radial erf well, depth E_bind), gated cubic drag rho_hat*b*v^3.

    ``v_c``/``p_tail`` (leg-A' extension, 2026-07-16) select the Slice-T1
    ``capped_cubic`` tail: gamma = b*v^2 in-band, b*v_c^2*(v/v_c)^p_tail
    above the cap. ``v_c=None`` (default) is the pure-cubic law with its
    original arithmetic — byte-inert.

    Returns dict with per-fragment (2, M) arrays: K (cooling exposure),
    v_inf (asymptotic speed, A/ps, residual-Coulomb-corrected), t_exit (ps),
    v_peak; optional center-pin profile samples.
    """
    if v_c is not None and p_tail is None:
        raise ValueError("v_c requires p_tail (the capped-cubic tail exponent)")
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
            E_BIND_ION_EV / (STEEP_A * np.sqrt(np.pi)) * np.exp(-((depth / STEEP_A) ** 2))
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
            "legaprime | legb | legc)"
        )


if __name__ == "__main__":
    main()
