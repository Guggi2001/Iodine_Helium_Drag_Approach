"""Slice-T8 analytic droplet-prior sampler (plan §I.11.3, T8-D2/T8-D3).

The MD arm draws *exactly* from the closed-form family (inverse-CDF on the
truncated ln-normal — no importance weighting on the MD side):

    kornilov_lognormal:        lnN ~ Normal(mu_ln, delta),
                               mu_ln = ln(<N>) − delta²/2
    pickup_weighted_lognormal: lnN ~ Normal(mu_ln + (2/3)·delta², delta)
                               (N^(2/3)·LogNormal(mu_ln) ∝ LogNormal(mu_ln +
                               (2/3)delta²) — the exact ln-normal identity)

both truncated to the twin's proposal window [250, 16000] (rule-1: the two
windows may not drift — pinned against the twin module here).

Tolerances: uniformity/mean checks are sample-size-based (3-sigma bands at
n = 20000, stated per assert); the truncated-mass pins are hand-computed
normal-CDF values with a half-last-digit allowance.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.config import (
    DROPLET_PRIOR_N_HI,
    DROPLET_PRIOR_N_LO,
    SimConfig,
)
from i2_helium_md.sampling.droplet_sizes import (
    analytic_prior_truncated_mass,
    sample_droplet_sizes_analytic,
)
from i2_helium_md.simulation import initial_state
from i2_helium_md.simulation.initial_state import build_initial_state

TWIN_SCRIPT = (
    Path(__file__).resolve().parents[1] / "scripts" / "tier2_h2b_forward_model.py"
)


def _cfg(prior="kornilov_lognormal", **kw):
    kw.setdefault("use_single_droplet_size", False)
    kw.setdefault("num_molecules", 20000)
    kw.setdefault("seed", 20260720)
    return SimConfig(droplet_size_prior=prior, **kw)


def test_truncation_window_matches_the_twin():
    # Rule-1 pin: the package window and the twin's uniform-in-lnN proposal
    # support are the same two numbers (T8-D2).
    spec = importlib.util.spec_from_file_location(
        "tier2_h2b_forward_model_window_pin", TWIN_SCRIPT
    )
    twin = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(twin)
    assert DROPLET_PRIOR_N_LO == twin.N_LO == 250.0
    assert DROPLET_PRIOR_N_HI == twin.N_HI == 16000.0


def test_samples_stay_inside_the_window():
    N = sample_droplet_sizes_analytic(_cfg())
    assert N.shape == (20000,)
    assert N.min() >= DROPLET_PRIOR_N_LO
    assert N.max() <= DROPLET_PRIOR_N_HI


def test_kornilov_is_the_truncated_lnnormal_exactly():
    # Invert the sampler's own transform analytically: with
    # z = (lnN − mu_ln)/delta, u = (Phi(z) − Phi(a)) / (Phi(b) − Phi(a))
    # must be U(0, 1) if and only if lnN is the truncated Normal(mu_ln, delta).
    from scipy.special import ndtr

    cfg = _cfg()
    N = sample_droplet_sizes_analytic(cfg)
    mu_ln = np.log(2000.0) - 0.5 * 0.625**2
    z = (np.log(N) - mu_ln) / 0.625
    a = (np.log(250.0) - mu_ln) / 0.625
    b = (np.log(16000.0) - mu_ln) / 0.625
    u = (ndtr(z) - ndtr(a)) / (ndtr(b) - ndtr(a))
    assert u.min() >= 0.0 and u.max() <= 1.0
    # mean of U(0,1) = 0.5, sigma = 1/sqrt(12); 3-sigma band at n = 20000
    # is 3 * 0.2887 / sqrt(20000) = 0.0061.
    assert abs(u.mean() - 0.5) < 0.0061
    # var of U(0,1) = 1/12; 3-sigma band on the sample var at n = 20000
    # (sigma_var ~ sqrt(1/180)/sqrt(n)) is 0.0016.
    assert abs(u.var() - 1.0 / 12.0) < 0.0016


def test_pickup_weight_is_the_exact_lnnormal_shift():
    # E[ln N]_pickup − E[ln N]_kornilov = (2/3)·delta² = 0.2604 up to the
    # (small, opposite-sign) truncation biases; 3-sigma statistical band on
    # the difference of means at n = 20000 is 3·delta·sqrt(2/n) = 0.019,
    # truncation bias <= 0.004 -> tolerance 0.023.
    N_k = sample_droplet_sizes_analytic(_cfg(seed=1))
    N_p = sample_droplet_sizes_analytic(_cfg("pickup_weighted_lognormal", seed=2))
    shift = float(np.log(N_p).mean() - np.log(N_k).mean())
    assert abs(shift - (2.0 / 3.0) * 0.625**2) < 0.023


def test_truncated_mass_closed_form_pins():
    # Hand-computed: delta = 0.625, <N> = 2000 -> mu_ln = 7.4056;
    # P(lnN < ln250) = Phi(-3.015) = 0.00129, P(lnN > ln16000) =
    # 1 - Phi(3.640) = 0.00014 -> 0.00142 total.
    assert abs(analytic_prior_truncated_mass(_cfg()) - 0.00142) < 1e-4
    # delta = 0.80 (the widest D4 family member): mu_ln = 7.2809;
    # Phi(-2.199) + (1 - Phi(2.999)) = 0.01393 + 0.00135 = 0.01528.
    cfg80 = _cfg(droplet_prior_delta=0.80)
    assert abs(analytic_prior_truncated_mass(cfg80) - 0.01528) < 3e-4


def test_heavy_truncation_warns():
    # No-silent-caps: a family placing > 5 % of its mass outside the window
    # must say so loudly. mean 15000, delta 0.8 puts ~31 % above N_HI.
    cfg = _cfg(droplet_prior_mean_N=15000.0, droplet_prior_delta=0.80)
    with pytest.warns(RuntimeWarning, match="truncat"):
        sample_droplet_sizes_analytic(cfg)


def test_initial_state_dispatch_uses_the_analytic_arm():
    from i2_helium_md.physics.constants import droplet_radius_bulk_angstrom

    cfg = _cfg(num_molecules=16)
    ckpt = build_initial_state(cfg, num_steps=1)
    radii = np.asarray(ckpt.droplet_radii, dtype=float)
    r_lo = float(droplet_radius_bulk_angstrom(DROPLET_PRIOR_N_LO))
    r_hi = float(droplet_radius_bulk_angstrom(DROPLET_PRIOR_N_HI))
    # per-atom layout (2N): each molecule's two atoms share its radius
    assert radii.shape == (32,)
    assert radii.min() >= r_lo and radii.max() <= r_hi
    # A continuous prior at n = 16 collides with probability ~0: the radii
    # must actually vary (the fixed-N path would give one unique value).
    assert np.unique(radii).size == 16


def test_default_path_never_touches_the_analytic_sampler(monkeypatch):
    # Byte-inert default (T8-D1): under ``legacy`` the dispatch is literally
    # the pre-T8 branch pair; the analytic sampler must not even be called.
    from i2_helium_md.physics.constants import droplet_radius_bulk_angstrom

    def _boom(*a, **k):
        raise AssertionError("analytic sampler called on the legacy path")

    monkeypatch.setattr(initial_state, "sample_droplet_sizes_analytic", _boom)
    ckpt = build_initial_state(SimConfig(num_molecules=4), num_steps=1)
    radii = np.asarray(ckpt.droplet_radii, dtype=float)
    assert np.allclose(radii, float(droplet_radius_bulk_angstrom(2000.0)))
