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
