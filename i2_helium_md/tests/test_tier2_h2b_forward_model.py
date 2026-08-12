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


# ---------------------------------------------------------------------------
# Leg-D extension (Slice T8 droplet-prior re-score)
# ---------------------------------------------------------------------------
def test_stage_legd_smoke_and_oracle(twin, monkeypatch, tmp_path):
    """The leg-D re-score (one lever flipped vs leg C: the T8
    ``kornilov_lognormal`` droplet prior, delta = 0.625 about <N> = 2000) runs
    at small m, writes both CSVs, and its ``d_delta`` rows (prior off — delta
    N = 2000) reproduce stage_legc's ``c`` rows exactly (the in-stage wiring
    oracle: same seed, same u/mu draws — leg D differs from leg C only through
    the droplet axis, never through sampling noise). The ``d`` rows genuinely
    spread the droplet axis (T8-D4)."""
    import csv as _csv

    monkeypatch.setattr(twin, "OUT", tmp_path)
    twin.stage_legc(m=200)
    twin.stage_legd(m=200)
    hist = tmp_path / "h2b_leg_d_predictions.csv"
    ke = tmp_path / "h2b_leg_d_ke.csv"
    assert hist.exists() and ke.exists()

    def rows(path, leg):
        with open(path, newline="") as fh:
            return [r for r in _csv.DictReader(fh) if r["leg"] == leg]

    # wiring oracle: d_delta (dressed, p = 1, delta prior) == leg-C `c` exactly.
    # Compared over the leg-C schema (the d rows add the N_q* columns).
    anchor = rows(tmp_path / "h2b_leg_c_predictions.csv", "c")
    d_delta = rows(hist, "d_delta")
    assert len(anchor) == 4 and len(d_delta) == 4
    for rc, rd in zip(anchor, d_delta):
        for key in rc:
            if key == "leg":
                continue
            assert rc[key] == rd[key], (key, rc[key], rd[key])

    # the kornilov leg genuinely samples the droplet axis
    d = rows(hist, "d")
    assert len(d) == 4
    for r in d:
        assert float(r["N_q05"]) < 2000.0 < float(r["N_q95"])
        hist_sum = sum(float(r[f"h{k}"]) for k in range(twin.N_STAR + 1))
        # 22 columns rounded at 4 decimals -> 2e-3 accumulation allowance
        assert abs(hist_sum - 1.0) < 2e-3
    # and the delta rows pin the N quantiles at exactly 2000
    for r in d_delta:
        assert float(r["N_q05"]) == 2000.0 and float(r["N_q95"]) == 2000.0


def test_stage_repilot1_smoke_and_oracle(twin, monkeypatch, tmp_path):
    """The Stage-1 v_c-bracket re-score (SI.11.4 pre-registration) runs at
    small m, writes both CSVs, spans exactly the three adjudicated brackets
    (c1/c4: 6.5..8.5; c2: 5.5..7.5 -- C3 dropped, RP-D2), and its in-stage
    S1-P1 wiring oracle passes: the bracket-center cells reproduce the leg-D
    ``d`` rows exactly (same SEED/draw order -- cells differ only through the
    drag tail, never through sampling noise)."""
    import csv as _csv

    monkeypatch.setattr(twin, "OUT", tmp_path)
    twin.stage_legd(m=200)  # provides the in-stage oracle's reference record
    twin.stage_repilot1(m=200)

    hist = tmp_path / "h2b_repilot_s1_predictions.csv"
    ke = tmp_path / "h2b_repilot_s1_ke.csv"
    assert hist.exists() and ke.exists()

    with open(hist, newline="") as fh:
        rows = [r for r in _csv.DictReader(fh) if r["leg"] == "s1"]
    assert len(rows) == 15
    labels = [r["config"] for r in rows]
    assert len(set(labels)) == 15
    # bracket coverage + cell-label convention (v_c x 10)
    for config, bracket in twin.REPILOT_S1_VC_BRACKETS.items():
        for v_c in bracket:
            cell = twin.repilot_s1_cell_label(config, v_c)
            row = next(r for r in rows if r["config"] == cell)
            assert float(row["v_c"]) == v_c
    assert "s1c1v75" in labels and "s1c2v55" in labels
    # C3 never rides (RP-D2)
    assert not any("c3" in lab for lab in labels)

    # c1/c4 share the chord integration (same drag tail): at equal v_c the
    # K quantiles must be identical -- the registered KE-degeneracy claim
    for v_c in twin.REPILOT_S1_VC_BRACKETS["c1"]:
        r1 = next(r for r in rows
                  if r["config"] == twin.repilot_s1_cell_label("c1", v_c))
        r4 = next(r for r in rows
                  if r["config"] == twin.repilot_s1_cell_label("c4", v_c))
        # tau rescale differs (3.8 vs 4.4), so compare tau-unscaled K; the
        # CSV stores K_q50 rounded at 4 decimals, so the tau re-multiply
        # carries up to ~tau * 5e-5 of pure round-off -- tolerance sized to
        # that, far below any physical K difference between drag tails
        k1 = float(r1["K_q50"]) * float(r1["tau_ps"])
        k4 = float(r4["K_q50"]) * float(r4["tau_ps"])
        assert abs(k1 - k4) < 5e-4


def test_stage_repilot1_requires_legd_record(twin, monkeypatch, tmp_path):
    """Without the leg-D prediction record on disk the S1-P1 oracle cannot
    run -- the stage must fail loud, not silently pre-register unverified
    rows."""
    monkeypatch.setattr(twin, "OUT", tmp_path)
    with pytest.raises(FileNotFoundError, match="leg-D"):
        twin.stage_repilot1(m=50)


def test_stage_repilot2_smoke_and_oracle(twin, monkeypatch, tmp_path):
    """The Stage-2 (tau, E0)-grid re-score (adjudication (b): two-candidate
    v_c carry) runs at small m, writes both CSVs, spans exactly the frozen
    67-cell grid (c1/c4 x {6.5, 8.5} x 4tau x 4E0 + the 3-cell c2 spot
    diagonal), and its in-stage S2s-P1 oracle passes: the five grid-center
    cells reproduce the Stage-1 rows exactly."""
    import csv as _csv

    monkeypatch.setattr(twin, "OUT", tmp_path)
    twin.stage_legd(m=200)      # leg-D record (repilot1's oracle input)
    twin.stage_repilot1(m=200)  # Stage-1 record (repilot2's oracle input)
    twin.stage_repilot2(m=200)

    hist = tmp_path / "h2b_repilot_s2_predictions.csv"
    ke = tmp_path / "h2b_repilot_s2_ke.csv"
    assert hist.exists() and ke.exists()

    with open(hist, newline="") as fh:
        rows = [r for r in _csv.DictReader(fh) if r["leg"] == "s2"]
    assert len(rows) == 67
    labels = [r["config"] for r in rows]
    assert len(set(labels)) == 67
    assert "s2c1v65t38e25" in labels and "s2c4v85t50e27" in labels
    assert "s2c2v65t40e24" in labels
    # tau/E0 budget-direction sanity on one v65 column: at fixed tau,
    # higher E0 must strip deeper (nbar down) -- the fate-map mechanics
    col = sorted(
        (float(r["E0_eV"]), float(r["nbar_det"]))
        for r in rows
        if r["config"].startswith("s2c1v65t38")
    )
    assert all(a[1] >= b[1] for a, b in zip(col, col[1:]))


def test_stage_repilot2_requires_stage1_record(twin, monkeypatch, tmp_path):
    """Without the Stage-1 record the S2s-P1 oracle cannot run -- fail loud."""
    monkeypatch.setattr(twin, "OUT", tmp_path)
    with pytest.raises(FileNotFoundError, match="Stage-1"):
        twin.stage_repilot2(m=50)


# ---------------------------------------------------------------------------
# G3 Step 1 -- twin landmark re-issue (atlas plan §3.5 G3)
# ---------------------------------------------------------------------------
def test_g3_md_reference_table_matches_grid(twin):
    """The embedded MD table covers exactly the 11 G1 cells, with the D0
    §14.1 spot values (r1l3 = the standing-geometry reference cell; r4l2 =
    the widest-marginal cell)."""
    labels = [c[0] for c in twin.G3_GRID_CELLS]
    assert len(labels) == 11 and len(set(labels)) == 11
    assert set(labels) == set(twin.G3_MD_G1_ROWS)
    r1l3 = dict(zip(twin.G3_MD_OBS_KEYS, twin.G3_MD_G1_ROWS["r1l3"]))
    assert r1l3["trap"] == 0.059 and r1l3["w1"] == 0.81
    assert r1l3["deepke"] == 0.51
    r4l2 = dict(zip(twin.G3_MD_OBS_KEYS, twin.G3_MD_G1_ROWS["r4l2"]))
    assert r4l2["trap"] == 0.800 and np.isnan(r4l2["midhot_geo"])
    # every anchored-law cell (l1/l2) has supp exactly 0 -- the confirmed
    # G1.1 prediction the twin rows are read against
    for label, _, law in twin.G3_GRID_CELLS:
        if law in ("l1", "l2"):
            assert twin.G3_MD_G1_ROWS[label][1] == 0.0


def test_g3_standing_cell_is_committed_s6_row(twin):
    """G3_STANDING must be one of the committed h2b_s6_final oracle rows."""
    import csv as _csv

    lad, v_c, p_tail, tau, e0 = twin.G3_STANDING
    with open(twin.OUT / "h2b_s6_final_predictions.csv", newline="") as fh:
        keys = {
            (r["ladder"], r["v_c"], r["tau_ps"], r["E0_eV"], r["p_tail"])
            for r in _csv.DictReader(fh)
        }
    assert (lad, str(v_c), str(tau), str(e0), str(p_tail)) in keys
    assert (lad, v_c, p_tail, tau, e0) in twin.G3_S6_CELLS


def test_g3_score_synthetic_conventions(twin):
    """Hand-checkable synthetic ensemble: trapped exclusion, solvated
    normalisation, both midHot conventions, deep-band membership, and the
    min_bin_count skip."""
    ref_ke = {1: 1.0, 2: 0.5, 8: 0.2, 10: 0.1, 17: 0.05}
    solv_exp = np.array([0.5, 0.5])
    m1 = float(twin.complex_mass_amu(1))
    m2 = float(twin.complex_mass_amu(2))
    m8 = float(twin.complex_mass_amu(8))
    m10 = float(twin.complex_mass_amu(10))

    def v_for(mass, ke_eV):
        return float(np.sqrt(2.0 * ke_eV * twin.EV_TO_AMU_A2_PS2 / mass))

    n_det = np.array([1, 1, 2, 8, 10, 5])
    trapped = np.array([False, False, False, False, False, True])
    sup = np.array([False, False, False, False, False, False])
    v = np.array([
        v_for(m1, 1.2), v_for(m1, 0.8),      # n=1 mean KE 1.0 eV
        v_for(m2, 0.75),                     # n=2 ratio 1.5
        v_for(m8, 0.4),                      # n=8 ratio 2.0
        v_for(m10, 0.05),                    # n=10 ratio 0.5
        99.0,                                # trapped -- must be ignored
    ])
    obs, ke_bins = twin.g3_score(
        n_det, sup, trapped, v, ref_ke, solv_exp, min_bin_count=1
    )
    assert abs(obs["trapped_frac"] - 1.0 / 6.0) < 1e-12  # over ALL fragments
    assert abs(obs["nbar"] - (1 + 1 + 2 + 8 + 10) / 5.0) < 1e-12
    assert abs(obs["n1_solv"] - 0.4) < 1e-12   # 2 of 5 detected, all solvated
    assert abs(obs["ratio"] - 2.0) < 1e-12
    # band n2-8: ratios {1.5, 2.0} -- arithmetic vs geometric must differ
    assert abs(obs["midhot_arith"] - (1.5 + 2.0) / 2.0) < 1e-9
    assert abs(obs["midhot_geo"] - np.sqrt(1.5 * 2.0)) < 1e-9
    assert obs["midhot_bins"] == 2
    assert abs(obs["deepke"] - 0.5) < 1e-9 and obs["deepke_bins"] == 1
    assert abs(obs["n1_ke"] - 1.0) < 1e-9
    assert [b[0] for b in ke_bins] == [1, 2, 8, 10]
    # min_bin_count=2 drops the singleton bins -> deep band empties to NaN
    obs2, ke_bins2 = twin.g3_score(
        n_det, sup, trapped, v, ref_ke, solv_exp, min_bin_count=2
    )
    assert [b[0] for b in ke_bins2] == [1]
    assert np.isnan(obs2["deepke"]) and obs2["deepke_bins"] == 0


def test_g3_grid_cell_draw_laws(twin):
    """l1 collapses to one center chord; l3 is r^2 on [0, R-3]; l2 (parent
    Boltzmann 313.2 K) reproduces the pre-registered exposure-table depth
    (mean depth 29.3 A at R = 34)."""
    R = 34.0
    rng = np.random.default_rng(twin.G3_SEED)
    r0, mu = twin._g3_grid_cell_draw("l1", R, 5000, rng)
    assert r0.shape == (1,) and r0[0] == 0.0 and abs(mu[0]) <= 1.0

    rng = np.random.default_rng(twin.G3_SEED)
    r0, _ = twin._g3_grid_cell_draw("l3", R, 20000, rng)
    x3 = (r0 / (R - 3.0)) ** 3
    for f in (0.25, 0.5, 0.75):
        assert abs(np.quantile(x3, f) - f) < 0.02

    rng = np.random.default_rng(twin.G3_SEED)
    r0, _ = twin._g3_grid_cell_draw("l2", R, 4000, rng)
    depth = R - r0
    assert abs(depth.mean() - 29.3) < 0.5     # plan §3.1 exposure table

    with pytest.raises(ValueError, match="unknown grid law"):
        twin._g3_grid_cell_draw("l9", R, 10, rng)


# ---------------------------------------------------------------------------
# G3 Step 2 -- the Route A/B twin scan (atlas plan §3.5c)
# ---------------------------------------------------------------------------
def _center_pin_args(twin):
    from i2_helium_md.physics.constants import droplet_radius_bulk_angstrom

    R2000 = float(droplet_radius_bulk_angstrom(2000.0))
    m21 = float(twin.complex_mass_amu(twin.N_STAR))
    return (
        np.array([0.0]), np.array([1.0]), np.array([R2000]), np.array([m21]),
    )


def test_e_bind_default_is_byte_inert(twin):
    """e_bind_ev=None and the explicit bundle stamp are the identical
    arithmetic path -- array_equal, not allclose (the g3scan extension must
    not perturb any recorded landmark)."""
    args = _center_pin_args(twin)
    ref = twin.integrate_pairs(*args, drag_on=True, v_c=7.25, p_tail=-1.0)
    explicit = twin.integrate_pairs(
        *args, drag_on=True, v_c=7.25, p_tail=-1.0,
        e_bind_ev=twin.E_BIND_ION_EV,
    )
    for key in ("K", "v_inf", "t_exit", "v_peak"):
        assert np.array_equal(ref[key], explicit[key]), key


def test_e_bind_lever_direction(twin):
    """A shallower well decelerates the exit less: smaller K, larger v_inf
    (the §6.7 item-2 lever direction); deeper well the other way."""
    args = _center_pin_args(twin)
    mid = twin.integrate_pairs(*args, drag_on=True, v_c=7.25, p_tail=-1.0)
    shallow = twin.integrate_pairs(
        *args, drag_on=True, v_c=7.25, p_tail=-1.0, e_bind_ev=0.0482
    )
    deep = twin.integrate_pairs(
        *args, drag_on=True, v_c=7.25, p_tail=-1.0, e_bind_ev=0.154
    )
    assert shallow["K"][0, 0] < mid["K"][0, 0] < deep["K"][0, 0]
    assert shallow["v_inf"][0, 0] > mid["v_inf"][0, 0] > deep["v_inf"][0, 0]


def test_g3scan_grid_matches_frozen_spec(twin):
    """The §3.5c frozen grids: chord surface 10 x 3 (30 integrations), free
    surface 6 x 36, nested total 6480; the v_c floor is the TDDFT band top;
    the diagnostic arms and the standing cell are present; gate bounds."""
    assert twin.G3SCAN_VC_TIER0 == (5.0, 5.5, 6.0, 6.5, 7.25, 8.0, 9.0, 10.0)
    assert min(twin.G3SCAN_VC_TIER0) == 5.0   # TDDFT band top 4.95 respected
    assert twin.G3SCAN_VC_DIAG == (3.5, 4.25)
    tags = [t for t, _ in twin.G3SCAN_EBIND]
    vals = dict(twin.G3SCAN_EBIND)
    assert tags == ["eb0482", "eb1168", "eb154"]  # §6.7 item-2 continuity
    assert vals["eb1168"] == twin.E_BIND_ION_EV   # standing = bundle stamp
    assert vals["eb0482"] == 0.0482 and vals["eb154"] == 0.154
    assert twin.G3SCAN_TAU_PS == (2.4, 3.2, 4.8, 6.4, 9.6, 12.8)
    assert twin.G3SCAN_TAU_SOURCED == 6.55
    e0 = twin.G3SCAN_E0_GRID
    assert e0[0] == 0.17 and e0[-1] == 0.52 and e0.size == 36
    assert np.allclose(np.diff(e0), 0.01)
    n_chord = len(twin.G3SCAN_VC_TIER0 + twin.G3SCAN_VC_DIAG) * len(
        twin.G3SCAN_EBIND
    )
    assert n_chord == 30
    assert n_chord * len(twin.G3SCAN_TAU_PS) * e0.size == 6480
    assert twin.G3SCAN_GATE_N1SOLV == (0.19, 0.30)
    assert twin.G3SCAN_GATE_NBAR == (4.4, 7.1)
    # the standing cell is on the grid (landmark continuity)
    assert twin.G3_STANDING[1] in twin.G3SCAN_VC_TIER0
    # pre-scan: f grid carries the 0.25 kill point and full exposure
    assert 0.25 in twin.G3SCAN_PRESCAN_F and 1.0 in twin.G3SCAN_PRESCAN_F
    assert twin.G3SCAN_PRESCAN_FMIN == 0.25


def test_g3scan_gate_boundaries(twin):
    """Hard gate = n1_solv AND nbar clauses, closed intervals; NaN never
    gates (an un-scoreable cell must not pass)."""
    assert twin.g3scan_gate(0.24, 5.0)
    assert twin.g3scan_gate(0.19, 4.4) and twin.g3scan_gate(0.30, 7.1)
    assert not twin.g3scan_gate(0.18, 5.0)
    assert not twin.g3scan_gate(0.31, 5.0)
    assert not twin.g3scan_gate(0.24, 4.3)
    assert not twin.g3scan_gate(0.24, 7.2)
    assert not twin.g3scan_gate(np.nan, 5.0)
    assert not twin.g3scan_gate(0.24, np.nan)


def test_g3scan_prescan_synthetic(twin):
    """Hand-built ensemble on the flat ladder: exposure scalings f sweep the
    fate map from no-strip to deep-strip, the recorded gate interval brackets
    exactly the f values whose light score lands both clauses, and the
    f = 1 / f = 0.25 columns match direct fate-map evaluations."""
    sig = twin.sigma_cum(twin.ladder_rungs("flat"))
    rng = np.random.default_rng(7)
    M = 4000
    ne = np.full(M, twin.N_STAR)
    K655 = rng.uniform(1.0, 3.0, M)
    trapped = np.zeros(M, dtype=bool)
    trapped[:200] = True  # 5 % trapped -- excluded from the light score
    rows, verdict = twin._g3scan_prescan(K655, trapped, ne, sig)
    assert len(rows) == len(twin.G3SCAN_TAU_PS) * twin.G3SCAN_E0_GRID.size
    assert verdict["cells_total"] == len(rows)
    # verdict counts are consistent with the per-row flags
    assert verdict["cells_n1_ok_fge025"] == sum(r["n1_ok_fge025"]
                                                for r in rows)
    assert verdict["cells_gate_ok_fge025"] == sum(r["gate_ok_fge025"]
                                                  for r in rows)
    assert verdict["route_a_killed"] == (verdict["cells_n1_ok_fge025"] == 0)
    # spot-check one row's f = 1.0 column against a direct evaluation
    r = next(row for row in rows
             if row["tau_ps"] == 3.2 and row["E0_eV"] == 0.27)
    K = K655 * (twin.TAU_PS / 3.2)
    n_det, _ = twin.fate_map(ne, K, 0.27, 1, sig)
    n1, nbar, _ = twin._g3scan_light(n_det, trapped)
    assert r["n1_solv_f100"] == round(n1, 4)
    assert r["nbar_f100"] == round(nbar, 3)
    assert r["gate_f100"] == int(twin.g3scan_gate(n1, nbar))
    # and the f = 0.25 column
    n_det, _ = twin.fate_map(ne, 0.25 * K, 0.27, 1, sig)
    n1_q, nbar_q, _ = twin._g3scan_light(n_det, trapped)
    assert r["n1_solv_f025"] == round(n1_q, 4)
    assert r["gate_f025"] == int(twin.g3scan_gate(n1_q, nbar_q))
    # gate interval, when present, is inside the f grid and ordered
    for row in rows:
        lo, hi = row["f_gate_lo"], row["f_gate_hi"]
        if not (np.isnan(lo) or np.isnan(hi)):
            assert 0.05 <= lo <= hi <= 1.0


def test_g3scan_light_score_conventions(twin):
    """Trapped exclusion + solvated normalisation + the strip diagnostic."""
    n_det = np.array([0, 1, 1, 3, 9, 21])
    trapped = np.array([False, False, False, False, False, True])
    n1, nbar, frac_le8 = twin._g3scan_light(n_det, trapped)
    assert abs(nbar - (0 + 1 + 1 + 3 + 9) / 5.0) < 1e-12
    assert abs(n1 - 2.0 / 4.0) < 1e-12         # 4 solvated, 2 at n=1
    assert abs(frac_le8 - 4.0 / 5.0) < 1e-12   # n<=8 of the detected read
    # nothing solvates -> NaN n1 (never gated), nbar still defined
    n1b, nbarb, _ = twin._g3scan_light(np.zeros(4, dtype=int),
                                       np.zeros(4, dtype=bool))
    assert np.isnan(n1b) and nbarb == 0.0


def test_g3scan_chord_tag_and_cache_guard(twin, tmp_path, monkeypatch):
    """Cache filenames are unambiguous per (v_c, E_bind); a stale stamp
    fails loud instead of silently serving another family's chord."""
    assert twin._g3scan_chord_tag(7.25, "eb1168") == "vc7p25_eb1168"
    assert twin._g3scan_chord_tag(5.0, "eb0482") == "vc5p0_eb0482"
    monkeypatch.setattr(twin, "OUT", tmp_path)
    bad = {"K": np.zeros((2, 1)), "v_inf": np.zeros((2, 1)),
           "t_exit": np.zeros((2, 1)), "trapped": np.zeros((2, 1)),
           "v_c": 9.0, "e_bind_ev": 0.154, "m": 20000.0,
           "seed": float(twin.G3_SEED)}
    np.savez_compressed(tmp_path / "h2b_g3scan_chord_vc9p0_eb154.npz", **bad)
    with pytest.raises(AssertionError, match="stale g3scan chord cache"):
        twin._g3scan_chord_family(9.0, "eb154", 0.154, 999, ens=None)


# ---------------------------------------------------------------------------
# G4 Step 1 / Block 1 — the fine ridge scan (atlas plan §3.5e)
# ---------------------------------------------------------------------------


def test_g4_grids_contain_the_g3_grids_as_exact_sublattices(twin):
    """G4-P1 rests on this: every Step-2 grid point survives in the fine grid.

    If a coarse point were lost the bit-exact reproduction oracle could not be
    evaluated there, and the fine scan would silently stop being a refinement.
    """
    assert set(twin.G4SCAN_VC) >= {5.0, 5.5, 6.0, 6.5}
    assert set(twin.G4SCAN_TAU_PS) >= {4.8, 6.4}
    fine = {round(float(x), 3) for x in twin.G4SCAN_E0_GRID}
    coarse = {round(float(x), 3) for x in twin.G3SCAN_E0_GRID
              if twin.G4SCAN_E0_GRID[0] <= float(x) <= twin.G4SCAN_E0_GRID[-1]}
    assert coarse <= fine
    assert twin.G4SCAN_EBIND == twin.G3SCAN_EBIND


def test_g4_vc_grid_is_quarter_stepped_across_the_basin(twin):
    assert sorted(float(v) for v in twin.G4SCAN_VC) == [
        5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5
    ]


def test_g4_tau_grid_is_0p4_stepped_and_reaches_past_the_sourced_value(twin):
    tau = sorted(float(t) for t in twin.G4SCAN_TAU_PS)
    assert tau == [4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 6.4, 6.8]
    assert max(tau) > twin.G3SCAN_TAU_SOURCED   # the flagged cell exists


def test_g4_e0_grid_is_half_stepped(twin):
    e0 = [float(x) for x in twin.G4SCAN_E0_GRID]
    assert e0[0] == 0.26 and e0[-1] == 0.42
    assert all(abs((b - a) - 0.005) < 1e-9 for a, b in zip(e0, e0[1:]))


def test_g4_gate_retires_the_crude_nbar_bracket(twin):
    """n1 stays the MD band (transfer is MD-grade); nbar moves to the MD band
    itself, because Block 0's bias model predicts MD nbar to +/- 0.3 He."""
    assert twin.G4SCAN_GATE_N1SOLV == (0.19, 0.30)
    assert twin.G4SCAN_GATE_NBAR == (3.77, 4.37)
    assert twin.G4SCAN_GATE_NBAR != twin.G3SCAN_GATE_NBAR


def test_g4_ranking_axes_exclude_the_unlicensed_deep_ke(twin):
    """Block 0 measured deepKE rho = +0.33 (< 0.7): the twin may not rank it."""
    assert "deepke" not in twin.G4SCAN_RANK_AXES
    assert set(twin.G4SCAN_RANK_AXES) == {"w1", "midhot"}


def test_g4_pareto_front_keeps_only_nondominated_cells(twin):
    """(W1, |ln midHot|) front: a cell dominated on both axes is dropped."""
    cells = [
        {"w1_solv": 1.0, "midhot_geo": 1.0},    # best midHot
        {"w1_solv": 0.5, "midhot_geo": 1.5},    # best W1
        {"w1_solv": 1.2, "midhot_geo": 1.6},    # dominated by both
    ]
    front = twin._g4_pareto_front(cells)
    assert front == [True, True, False]


def test_g4_ridge_connectivity_detects_a_gap(twin):
    """G4-P2 is a 4-neighbour path test on the (v_c, tau) lattice."""
    vc, tau = twin.G4SCAN_VC, twin.G4SCAN_TAU_PS
    connected = [(5.5, 4.8), (5.75, 4.8), (5.75, 5.2), (5.75, 5.6),
                 (5.75, 6.0), (5.75, 6.4), (6.0, 6.4)]
    assert twin._g4_ridge_connected(connected, (5.5, 4.8), (6.0, 6.4), vc, tau)
    assert not twin._g4_ridge_connected(
        [(5.5, 4.8), (6.0, 6.4)], (5.5, 4.8), (6.0, 6.4), vc, tau
    )


# ---------------------------------------------------------------------------
# Atlas §6.5 Step 2 — the E_bind twin scan (designed 2026-08-10)
# ---------------------------------------------------------------------------
def _ebind_rows(twin, arm="H", slope=-0.55, mid_slope=-0.70, trap_slope=0.85,
                zero_ke1=None):
    """Synthetic scan rows on an exact line, in the stage's row schema."""
    rows = []
    for tag, e_bind, in_prov in twin.EBINDSCAN_WELLS:
        ke1 = 0.90 + slope * e_bind
        if tag == "eb0" and zero_ke1 is not None:
            ke1 = zero_ke1
        rows.append({
            "arm": arm, "cell": "synthetic", "form": "capped",
            "tau_ps": 4.4, "E0_eV": 0.405,
            "E_bind_tag": tag, "E_bind_eV": e_bind, "in_provenance": in_prov,
            "n1_ke_eV": ke1,
            "ke2_eV": 0.80 + slope * e_bind,
            "ke_mid_geo_eV": 0.60 + mid_slope * e_bind,
            "ke_deep_eV": 0.20 + mid_slope * e_bind,
            "trapped_frac": 0.02 + trap_slope * e_bind,
        })
    return rows


def test_ebindscan_well_grid_is_the_tier0_provenance_spread(twin):
    """Plan §6.5: five co-extracted wells + two labelled diagnostics."""
    tags = [w[0] for w in twin.EBINDSCAN_WELLS]
    vals = [w[1] for w in twin.EBINDSCAN_WELLS]
    assert len(set(tags)) == len(tags) == 7
    assert vals == sorted(vals)
    prov = [w[1] for w in twin.EBINDSCAN_WELLS if w[2] == 1]
    diag = [w[1] for w in twin.EBINDSCAN_WELLS if w[2] == 0]
    assert prov == [0.0482, 0.071, 0.113, twin.E_BIND_ION_EV, 0.154]
    assert diag == [0.0, 0.2168]
    # the bundle stamp is the standing well, carried exactly (no re-typing)
    assert dict((t, v) for t, v, _ in twin.EBINDSCAN_WELLS)["eb1168"] == (
        twin.E_BIND_ION_EV
    )


def test_ebindscan_arms_pin_h405_and_the_md_measured_lin_chord(twin):
    arms = {a[0]: a for a in twin.EBINDSCAN_ARMS}
    assert arms["H"] == ("H", "h405", "capped", 5.5, 4.4, 0.405)
    assert arms["L"] == ("L", "lr1", "lin", 27.5, 4.8, 0.35)


def test_ebindscan_prediction_bands_are_frozen(twin):
    """EB-P1..P6 thresholds frozen 2026-08-10 before any number was read."""
    assert twin.EBINDSCAN_P1_MAX == 0.56
    assert twin.EBINDSCAN_P2_MAX_RESID_EV == 0.005
    assert twin.EBINDSCAN_P3_BAND == (0.10, 0.20)
    assert twin.EBINDSCAN_P4_BAND == (0.6, 1.1)
    assert twin.EBINDSCAN_P5_KE1 == 0.75


def test_ebindscan_ols_recovers_an_exact_line(twin):
    x = [0.0482, 0.071, 0.113, 0.1168, 0.154]
    y = [0.9 - 0.55 * xi for xi in x]
    slope, icept, resid, n = twin._ebindscan_ols(x, y)
    assert abs(slope + 0.55) < 1e-12
    assert abs(icept - 0.9) < 1e-12
    assert resid < 1e-12
    assert n == 5


def test_ebindscan_ols_rejects_a_two_point_fit(twin):
    with pytest.raises(ValueError, match=">= 3 finite"):
        twin._ebindscan_ols([0.0482, 0.1168], [0.9, 0.85])


def test_ebindscan_band_ev_aggregations(twin):
    ke_bins = [(1, 0.2, 1.0), (2, 0.1, 0.8), (4, 0.1, 0.2),
               (10, 0.1, 0.3), (12, 0.1, 0.1)]
    mid, mid_n = twin._ebindscan_band_ev(ke_bins, 2, 8, "geometric")
    deep, deep_n = twin._ebindscan_band_ev(ke_bins, 10, 17, "arithmetic")
    assert (mid_n, deep_n) == (2, 2)
    assert abs(mid - np.sqrt(0.8 * 0.2)) < 1e-12
    assert abs(deep - 0.2) < 1e-12
    empty, empty_n = twin._ebindscan_band_ev(ke_bins, 15, 17, "arithmetic")
    assert empty_n == 0 and np.isnan(empty)


def test_ebindscan_md_ring_slopes_read_the_committed_artifact(twin):
    """EB-P1's reference number is read from the ring CSV, never hardcoded."""
    slopes = twin._ebindscan_md_ring_slopes()
    assert set(slopes) == {"lr1_lr2", "lr3_lr4", "mean"}
    # both CRN pairs measured a partial refund: strictly between the rigid
    # -1 prediction and no response at all
    for key in ("lr1_lr2", "lr3_lr4"):
        assert -1.0 < slopes[key] < 0.0
    assert abs(slopes["mean"] + 0.52) < 0.03


def test_ebindscan_summary_passes_its_predictions_on_the_designed_case(twin):
    md_ring = {"mean": -0.52}
    rows = _ebind_rows(twin, slope=-0.50, mid_slope=-0.65, trap_slope=0.85,
                       zero_ke1=0.95)
    s = twin._ebindscan_arm_summary("H", rows, md_ring)
    assert s["EB_P1"] == "PASS"          # |slope| 0.50 < 0.56
    assert s["EB_P2"] == "PASS"          # exact line
    assert s["EB_P3"] == "PASS"          # 0.65 - 0.50 = 0.15 in [0.10, 0.20]
    assert s["EB_P4"] == "PASS"          # 0.85/eV in [0.6, 1.1]
    assert s["EB_P6"] == "PASS"          # 0.95 above the 0.90 extrapolation
    assert abs(s["EB_P3_delta"] - 0.15) < 1e-9
    assert s["n_fit"] == 5


def test_ebindscan_summary_fires_each_falsifier(twin):
    md_ring = {"mean": -0.52}
    # EB-P1/P3: a rigid (unrefunded) response
    rigid = twin._ebindscan_arm_summary(
        "H", _ebind_rows(twin, slope=-1.0, mid_slope=-1.0), md_ring)
    assert rigid["EB_P1"] == "FAIL" and rigid["EB_P3"] == "FAIL"
    # EB-P4: trap lever outside the D0 §9 band
    flat_trap = twin._ebindscan_arm_summary(
        "H", _ebind_rows(twin, trap_slope=0.05), md_ring)
    assert flat_trap["EB_P4"] == "FAIL"
    # EB-P6: E_bind = 0 sitting below the provenance extrapolation
    below = twin._ebindscan_arm_summary(
        "H", _ebind_rows(twin, zero_ke1=0.80), md_ring)
    assert below["EB_P6"] == "FAIL"
    # EB-P5 is arm-H only and fires when a provenance well clears the band
    high = _ebind_rows(twin, slope=-0.50)
    for r in high:
        r["n1_ke_eV"] += 0.20
    assert twin._ebindscan_arm_summary("H", high, md_ring)["EB_P5"] == "FAIL"
    high_l = _ebind_rows(twin, arm="L", slope=-0.50)
    for r in high_l:
        r["n1_ke_eV"] += 0.20
    assert twin._ebindscan_arm_summary("L", high_l, md_ring)["EB_P5"] == ""


def test_ebindscan_summary_flags_a_nonlinear_response(twin):
    """EB-P2 is the linearity envelope: a bent response must fail it."""
    md_ring = {"mean": -0.52}
    rows = _ebind_rows(twin, slope=-0.50)
    for r in rows:
        r["n1_ke_eV"] += 10.0 * r["E_bind_eV"] ** 2  # bend ~ 0.014 eV resid
    s = twin._ebindscan_arm_summary("H", rows, md_ring)
    assert s["EB_P2"] == "FAIL"
    assert s["resid_max_KE1_eV"] > twin.EBINDSCAN_P2_MAX_RESID_EV


class TestLinTauRefinement:
    """Plan §6.7: the τ refinement stage must not disturb the committed
    sweep, and its anchors must reproduce it."""

    def test_default_tau_grid_is_unchanged(self):
        from scripts.tier2_h2b_forward_model import (
            G3SCAN_TAU_PS, LINTAU_GRID, LINTAU_ANCHOR_TAUS,
        )
        # The committed grid must not have moved — every earlier stage and
        # every committed CSV depends on it.
        assert G3SCAN_TAU_PS == (2.4, 3.2, 4.8, 6.4, 9.6, 12.8)
        # The refinement fills the 3.2 -> 4.8 gap and carries anchors.
        assert LINTAU_ANCHOR_TAUS == (4.8, 6.4)
        assert set(LINTAU_ANCHOR_TAUS) <= set(LINTAU_GRID)
        assert set(LINTAU_ANCHOR_TAUS) <= set(G3SCAN_TAU_PS)
        new = sorted(set(LINTAU_GRID) - set(G3SCAN_TAU_PS))
        assert new == [3.6, 4.0, 4.4, 5.2, 5.6]
        # every refined value sits inside the committed grid's span
        assert min(LINTAU_GRID) >= min(G3SCAN_TAU_PS)
        assert max(LINTAU_GRID) <= max(G3SCAN_TAU_PS)

    def test_scan_signature_defaults_preserve_committed_behaviour(self):
        import inspect

        from scripts.tier2_h2b_forward_model import _linsweep_scan

        sig = inspect.signature(_linsweep_scan)
        # Both new knobs must default to the committed behaviour, or the
        # linscan/linqscan CSVs would silently change schema/values.
        assert sig.parameters["tau_grid"].default is None
        assert sig.parameters["with_tail"].default is False

    def test_anchor_oracle_rejects_a_drifted_row(self):
        from scripts.tier2_h2b_forward_model import _lintau_anchor_oracle

        # A row claiming to be at an anchor τ but carrying a wrong value
        # must fail loudly — this oracle is the licence for the whole scan.
        bad = [{
            "a": "35.0", "E_bind_tag": "eb0482", "tau_ps": "4.8",
            "E0_eV": "0.37", "trapped_frac": 0.0128, "suppressed_frac": 0.1242,
            "nbar_det": 4.625, "n1_solv": 0.2008, "w1_solv": 0.9999,
            "n1_ke_eV": 0.7027, "ke2_eV": 0.6322, "deepke": 1.8634,
            "midhot_arith": 1.5300, "gate": 1,
        }]
        with pytest.raises(AssertionError, match="lintau anchor"):
            _lintau_anchor_oracle(bad)
