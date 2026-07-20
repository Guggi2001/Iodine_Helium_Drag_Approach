"""Tests for the committed 1D chord-model twin (V0-3, plan SI.11).

The script is the verbatim recovery of the H.2b scratchpad driver plus the
V0-3 scope additions (repo-relative bootstrap, birth-law input, ``birthlaw``
stage). These tests pin the recorded wiring oracles and the V0-3 additions;
the physics itself was validated against MD landmarks at execution time
(findings SS4j; oracles below reproduce the recorded values).

Tolerances: the wiring-oracle K values were recorded at 5 decimals in plan
SI.11 V0-3, so the pins allow 1.5e-5 (half-ulp of the recording plus float
jitter); Sigma(21) was recorded at 8 decimals -> 1e-8. The Boltzmann-draw
quantile band is +/-0.25 A on a 500-sample median of a distribution whose
recorded median is 1.35-1.37 A (generous sample-size band, not a physics
tolerance).
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

SCRIPT = (
    Path(__file__).resolve().parents[1] / "scripts" / "tier2_h2b_forward_model.py"
)


@pytest.fixture(scope="module")
def twin():
    spec = importlib.util.spec_from_file_location(
        "tier2_h2b_forward_model_under_test", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# Recorded wiring oracles (plan SI.11 V0-3)
# ---------------------------------------------------------------------------
def test_sigma21_pin(twin):
    sig = twin.sigma_cum(twin.ladder_rungs("flat"))
    assert abs(sig[twin.N_STAR] - 0.18783720) < 1e-8


def test_production_and_9a_center_pin_K(twin):
    from i2_helium_md.physics.constants import droplet_radius_bulk_angstrom

    R2000 = float(droplet_radius_bulk_angstrom(2000.0))
    m21 = float(twin.complex_mass_amu(twin.N_STAR))
    args = (
        np.array([0.0]), np.array([1.0]), np.array([R2000]), np.array([m21]),
    )
    res = twin.integrate_pairs(*args, r0_sep=twin.R0_SEP_PROD_A, drag_on=True)
    assert abs(res["K"][0, 0] - 0.74460) < 1.5e-5
    res9 = twin.integrate_pairs(*args, r0_sep=twin.R0_SEP_9A_A, drag_on=True)
    assert abs(res9["K"][0, 0] - 0.89767) < 1.5e-5


def test_fate_map_recorded_o4(twin):
    """O4: measured production K -> n_det = 6 at E0 = 0.28; suppressed at 0.52."""
    sig = twin.sigma_cum(twin.ladder_rungs("flat"))
    K = np.array([0.74460])
    ne = np.array([twin.N_STAR])
    n_det, sup = twin.fate_map(ne, K, 0.28, 0, sig)
    assert n_det[0] == 6 and not sup[0]
    n_det2, sup2 = twin.fate_map(ne, K, 0.52, 0, sig)
    assert sup2[0] and n_det2[0] == 0


def test_ladder_variants(twin):
    flat = twin.ladder_rungs("flat")
    graded = twin.ladder_rungs("rq4graded")
    assert np.allclose(graded[:3], flat[:3] * [2.2, 1.5, 1.3])
    assert np.array_equal(graded[3:], flat[3:])
    floor1 = twin.ladder_rungs("floor1")
    assert abs(floor1[0] - twin.X2_RUNG_EV) < 1e-12
    assert np.array_equal(floor1[1:], flat[1:])
    with pytest.raises(ValueError):
        twin.ladder_rungs("nope")


# ---------------------------------------------------------------------------
# V0-3 birth-law input
# ---------------------------------------------------------------------------
def test_uniform_volume_master_is_r_squared(twin):
    """x = U^(1/3) -> x^3 uniform on [0, 1] (the r^2 volume law)."""
    rng = np.random.default_rng(1234)
    ms = twin.draw_master(rng, 20000, birth_law="uniform_volume")
    x3 = ms["x"] ** 3
    for f in (0.25, 0.5, 0.75):
        assert abs(np.quantile(x3, f) - f) < 0.02  # n=20000 sampling band


def test_margin_weight_uniform_volume(twin):
    R = np.full(4, 28.0)
    x = np.array([0.0, 0.5, 0.892, 0.95])  # cap at m=3: 1 - 3/28 = 0.8929
    w = twin.margin_weight(3.0, x, R)
    cap = 1.0 - 3.0 / 28.0
    assert np.allclose(w[:3], 1.0 / cap**3)
    assert w[3] == 0.0


def test_boltzmann_draw_center_pinned(twin, monkeypatch):
    """Fixed N = 2000: the realized Boltzmann law is center-pinned
    (recorded V0-1 median 1.35-1.37 A)."""
    monkeypatch.setattr(twin, "N_LO", 2000.0)
    monkeypatch.setattr(twin, "N_HI", 2000.0)
    rng = np.random.default_rng(twin.SEED)
    ms = twin.draw_master(rng, 500, birth_law="boltzmann")
    r0 = ms["x"] * ms["R"]
    assert abs(np.median(r0) - 1.36) < 0.25
    assert np.quantile(r0, 0.99) < 5.0  # nothing anywhere near the surface


def test_boltzmann_margin_guard(twin, monkeypatch):
    monkeypatch.setattr(twin, "BIRTH_LAW", "boltzmann")
    x = np.array([0.1])
    R = np.array([28.0])
    assert np.array_equal(twin.margin_weight(0.0, x, R), np.ones(1))
    with pytest.raises(ValueError, match="uniform_volume"):
        twin.margin_weight(3.0, x, R)


def test_unknown_birth_law_raises(twin):
    with pytest.raises(ValueError, match="birth law"):
        twin.draw_master(np.random.default_rng(0), 10, birth_law="nope")


def test_stage_birthlaw_smoke(twin, monkeypatch, tmp_path):
    """The V0-1 re-issue runs and writes its CSV (no production dirs touched)."""
    monkeypatch.setattr(twin, "OUT", tmp_path)
    twin.stage_birthlaw()
    out = tmp_path / "h2b_birth_law_quantiles.csv"
    assert out.exists()
    text = out.read_text()
    assert "boltzmann_analytic" in text and "uniform_volume" in text


# ---------------------------------------------------------------------------
# Leg-A' extension: capped tail + re-score stage
# ---------------------------------------------------------------------------
def test_drag_gamma_tail_hand_values(twin):
    """gamma = b*v^2 in-band; b*v_c^2*(v/v_c)^p above the cap (Slice-T1 form)."""
    b = twin.B_DRAG
    v = np.array([3.0, 10.0])
    g0 = twin.drag_gamma_tail_amu_per_ps(v, v_c=5.0, p_tail=0.0)
    assert np.allclose(g0, [b * 9.0, b * 25.0])          # tail p=0: constant gamma
    gm1 = twin.drag_gamma_tail_amu_per_ps(v, v_c=5.0, p_tail=-1.0)
    assert np.allclose(gm1, [b * 9.0, b * 25.0 * 0.5])   # tail p=-1: gamma ~ 1/v
    assert twin.drag_gamma_tail_amu_per_ps(np.array([0.0]), 5.0, -1.0)[0] == 0.0


def test_capped_tail_in_band_matches_pure_cubic(twin):
    """v_c above the peak speed: capped result ~ the pure-cubic result."""
    from i2_helium_md.physics.constants import droplet_radius_bulk_angstrom

    R2000 = float(droplet_radius_bulk_angstrom(2000.0))
    m21 = float(twin.complex_mass_amu(twin.N_STAR))
    args = (
        np.array([0.0]), np.array([1.0]), np.array([R2000]), np.array([m21]),
    )
    ref = twin.integrate_pairs(*args, drag_on=True)
    capped = twin.integrate_pairs(*args, drag_on=True, v_c=1e6, p_tail=0.0)
    # different multiplication grouping -> allclose, not array_equal
    assert np.allclose(capped["K"], ref["K"], rtol=1e-10)
    assert np.allclose(capped["v_inf"], ref["v_inf"], rtol=1e-10)


def test_capped_tail_lifts_speed(twin):
    """A saturated tail (p=-1) drags less at high v -> larger v_inf, smaller K."""
    from i2_helium_md.physics.constants import droplet_radius_bulk_angstrom

    R2000 = float(droplet_radius_bulk_angstrom(2000.0))
    m21 = float(twin.complex_mass_amu(twin.N_STAR))
    args = (
        np.array([0.0]), np.array([1.0]), np.array([R2000]), np.array([m21]),
    )
    ref = twin.integrate_pairs(*args, drag_on=True)
    tail = twin.integrate_pairs(*args, drag_on=True, v_c=7.5, p_tail=-1.0)
    assert tail["v_inf"][0, 0] > ref["v_inf"][0, 0]
    assert tail["K"][0, 0] < ref["K"][0, 0]


def test_v_c_without_p_tail_raises(twin):
    with pytest.raises(ValueError, match="p_tail"):
        twin.integrate_pairs(
            np.array([0.0]), np.array([1.0]), np.array([28.0]),
            np.array([210.0]), v_c=7.5,
        )


def test_stage_legaprime_smoke(twin, monkeypatch, tmp_path):
    """The leg-A' re-score runs end-to-end at small m and writes both CSVs."""
    monkeypatch.setattr(twin, "OUT", tmp_path)
    twin.stage_legaprime(m=200)
    hist = tmp_path / "h2b_leg_aprime_predictions.csv"
    ke = tmp_path / "h2b_leg_aprime_ke.csv"
    assert hist.exists() and ke.exists()
    text = hist.read_text()
    # both legs x all four configs present
    assert text.count("aprime,") == 4 and text.count("\na,") == 4


# ---------------------------------------------------------------------------
# Leg-B extension (Slice T5 dressed re-score)
# ---------------------------------------------------------------------------
def test_stage_legb_smoke_and_anchor(twin, monkeypatch, tmp_path):
    """The leg-B re-score runs at small m, writes both CSVs, and its
    undressed anchor rows reproduce stage_legaprime's aprime rows exactly
    (same seed, same draws — the in-stage wiring oracle)."""
    import csv as _csv

    monkeypatch.setattr(twin, "OUT", tmp_path)
    twin.stage_legaprime(m=200)
    twin.stage_legb(m=200)
    hist = tmp_path / "h2b_leg_b_predictions.csv"
    ke = tmp_path / "h2b_leg_b_ke.csv"
    assert hist.exists() and ke.exists()

    def rows(path, leg):
        with open(path, newline="") as fh:
            return [r for r in _csv.DictReader(fh) if r["leg"] == leg]

    anchor = rows(hist, "b_undressed")
    ap = rows(tmp_path / "h2b_leg_aprime_predictions.csv", "aprime")
    assert len(anchor) == 4 and len(ap) == 4
    for ra, rb in zip(ap, anchor):
        for key in ra:
            if key == "leg":
                continue
            assert ra[key] == rb[key], (key, ra[key], rb[key])
    # the dressed leg exists and genuinely moves: n_eject < 21 fragments
    dressed = rows(hist, "b")
    assert len(dressed) == 4
    for r in dressed:
        assert float(r["n_eject_q05"]) < twin.N_STAR


# ---------------------------------------------------------------------------
# Leg-C extension (Slice T6 E_int(0)-dressing p-law re-score)
# ---------------------------------------------------------------------------
def test_stage_legc_smoke_and_oracle(twin, monkeypatch, tmp_path):
    """The leg-C re-score (one lever flipped vs leg B: the T6 p = 1 onset
    coupling) runs at small m, writes both CSVs, and its ``c_p0`` rows (dressed,
    p = 0) reproduce stage_legb's ``b`` rows exactly (the in-stage wiring
    oracle). p = 1 can only lower the onset (ratio <= 1), so it regularizes the
    over-suppression: suppressed_frac falls and nbar_det rises vs p = 0."""
    import csv as _csv

    monkeypatch.setattr(twin, "OUT", tmp_path)
    twin.stage_legb(m=200)
    twin.stage_legc(m=200)
    hist = tmp_path / "h2b_leg_c_predictions.csv"
    ke = tmp_path / "h2b_leg_c_ke.csv"
    assert hist.exists() and ke.exists()

    def rows(path, leg):
        with open(path, newline="") as fh:
            return [r for r in _csv.DictReader(fh) if r["leg"] == leg]

    # wiring oracle: c_p0 (dressed, p = 0) == leg-B `b` (dressed, p = 0) exactly
    anchor = rows(tmp_path / "h2b_leg_b_predictions.csv", "b")
    c_p0 = rows(hist, "c_p0")
    assert len(anchor) == 4 and len(c_p0) == 4
    for rb, rc in zip(anchor, c_p0):
        for key in rb:
            if key == "leg":
                continue
            assert rb[key] == rc[key], (key, rb[key], rc[key])

    # the p = 1 leg genuinely moves in the pre-registered direction
    dressed = {r["config"]: r for r in rows(hist, "c")}
    base = {r["config"]: r for r in c_p0}
    assert len(dressed) == 4
    for cfg in base:
        assert float(dressed[cfg]["suppressed_frac"]) <= float(
            base[cfg]["suppressed_frac"]
        ) + 1e-9
        assert float(dressed[cfg]["nbar_det"]) >= float(base[cfg]["nbar_det"]) - 1e-9
    # at least one config strictly de-suppresses (the mechanism is live)
    assert any(
        float(dressed[cfg]["suppressed_frac"]) < float(base[cfg]["suppressed_frac"])
        for cfg in base
    )
