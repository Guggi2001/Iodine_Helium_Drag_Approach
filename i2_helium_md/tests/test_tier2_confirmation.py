"""Tests for the Tier-2 confirmation-matrix scorer (§I.11.4 Stage 0).

Covers ``i2_helium_md/postprocess/tier2_confirmation.py`` — the repo-ization
of the T9 leg-scoring conventions (findings §4m–§4r) plus the experimental
scoring layer (T4-absorbed):

- detection-read conventions: ``droplet_retained`` excluded from every read,
  ``suppressed`` scored in bin 0, histogram on integer support 0..n_max;
- twin CSV loaders (``h2b_leg_*_predictions.csv`` / ``_ke.csv`` contract);
- W1 between fraction vectors (reuses ``wasserstein_integer_support`` — no
  second metric implementation);
- the solvated n >= 1 renormalization (RQ8: bare is not a model target);
- the committed ihe_ked error model: per-point sqrt(statErr^2 + sysErr^2)
  scatter plus the calib/condition fractional bands as two *correlated*
  whole-curve shifts, profiled as unit-normal nuisances (closed-form 2x2).

All detection inputs are tiny synthetic ``DetectionResult`` objects — no
production checkpoints, no figures (testing rules).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.postprocess.abundance_loader import HeAbundanceReference
from i2_helium_md.postprocess.ihe_ked import IHeKedReference
from i2_helium_md.postprocess.tier2_confirmation import (
    read_confirmation_detection,
    load_twin_prediction,
    load_twin_ke_curve,
    wasserstein_between,
    solvated_renormalized,
    score_histogram_vs_reference,
    score_ke_curve_vs_reference,
)
from i2_helium_md.simulation.detection_stage import DetectionResult


# ---------------------------------------------------------------------------
# synthetic-input builders
# ---------------------------------------------------------------------------

def make_detection(
    n_detected,
    state_reason,
    ke_eV,
) -> DetectionResult:
    """Minimal synthetic ``DetectionResult`` (no events, zero ledger fields)."""
    n_arr = np.asarray(n_detected, dtype=float)
    m = n_arr.size
    zeros = np.zeros(m, dtype=float)
    return DetectionResult(
        num_molecules=max(m // 2, 1),
        t_handover_ps=30.0,
        detection_time_ps=8.53e6,
        n_detected=n_arr,
        E_int_detected_eV=zeros.copy(),
        mass_detected_kg=zeros.copy(),
        vx_detected=zeros.copy(),
        vy_detected=zeros.copy(),
        vz_detected=zeros.copy(),
        E_kin_detected_eV=np.asarray(ke_eV, dtype=float),
        E_pot_detected_eV=zeros.copy(),
        E_dissip_detected_eV=zeros.copy(),
        E_mass_transfer_detected_eV=zeros.copy(),
        state_reason=np.asarray(state_reason, dtype="<U16"),
        event_offsets=np.zeros(m + 1, dtype=np.int64),
        event_time_ps=np.zeros(0, dtype=float),
        event_pre_shed_n=np.zeros(0, dtype=float),
        event_dE_int_eV=np.zeros(0, dtype=float),
        event_dE_bind_fold_eV=np.zeros(0, dtype=float),
        event_dE_mass_transfer_eV=np.zeros(0, dtype=float),
    )


def make_abundance_reference(n, fraction) -> HeAbundanceReference:
    """Tiny in-memory abundance reference (fractions must sum to 1)."""
    n = np.asarray(n, dtype=int)
    fraction = np.asarray(fraction, dtype=float)
    return HeAbundanceReference(
        n=n,
        label=np.asarray([f"I^+He_{k}" for k in n], dtype="<U16"),
        mass_center_u=127.0 + 4.0026 * n.astype(float),
        ion_counts=fraction.copy(),
        ion_fraction=fraction,
        source_path=Path("synthetic"),
    )


def make_ked_reference(
    n,
    mean_KE_eV,
    *,
    stat_err=0.01,
    sys_err=0.0,
    calib_frac=0.04,
    condition_frac=0.06,
) -> IHeKedReference:
    """Tiny in-memory ihe_ked reference with uniform error columns."""
    n = np.asarray(n, dtype=int)
    mean = np.asarray(mean_KE_eV, dtype=float)
    ones = np.ones(n.size, dtype=float)
    return IHeKedReference(
        n=n,
        label=np.asarray([f"I+He_{k}" for k in n], dtype="<U16"),
        mass_center_u=127.0 + 4.0026 * n.astype(float),
        N_counts=1000.0 * ones,
        N_eff=800.0 * ones,
        mean_KE_eV=mean,
        mode_KE_eV=mean.copy(),
        median_KE_eV=mean.copy(),
        sigma_KE_eV=0.5 * ones,
        stat_err_mean_KE_eV=stat_err * ones,
        sys_err_mean_KE_eV=sys_err * ones,
        calib_syst_frac=calib_frac * ones,
        condition_syst_frac=condition_frac * ones,
        bg_off_shift_eV=0.0 * ones,
        dominant_error=np.asarray(["calib"] * n.size, dtype="<U16"),
        noise_limited=np.zeros(n.size, dtype=int),
        source_path=Path("synthetic"),
    )


_TWIN_PRED_HEADER = (
    "leg,config,ladder,tau_ps,E0_eV,v_c,p_tail,trapped_frac,suppressed_frac,"
    "bare_frac,nbar_det,K_q05,K_q50,K_q95,n_eject_q05,n_eject_q50,"
    "n_eject_mean,N_q05,N_q50,N_q95,"
    + ",".join(f"h{k}" for k in range(22))
)


def write_twin_prediction_csv(path: Path, rows: list[dict]) -> Path:
    lines = [_TWIN_PRED_HEADER]
    for row in rows:
        h = row["h"]
        lines.append(
            f"{row['leg']},{row['config']},rq4graded,3.8,0.25,7.5,-1.0,"
            f"{row.get('trapped_frac', 0.0)},{row.get('suppressed_frac', h[0])},"
            f"{row.get('bare_frac', h[0])},{row['nbar_det']},"
            "0.1,0.5,25.0,13.0,17.0,16.7,2000.0,2000.0,2000.0,"
            + ",".join(f"{v}" for v in h)
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_twin_ke_csv(path: Path, rows: list[tuple[str, str, int, float, float]]) -> Path:
    lines = ["leg,config,n,weight,mean_KE_eV"]
    for leg, config, n, weight, ke in rows:
        lines.append(f"{leg},{config},{n},{weight},{ke}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# read_confirmation_detection — the §4r conventions
# ---------------------------------------------------------------------------

class TestReadConfirmationDetection:
    def test_droplet_retained_excluded_from_every_read(self):
        det = make_detection(
            n_detected=[2.0, 3.0, 21.0],
            state_reason=["frozen", "frozen", "droplet_retained"],
            ke_eV=[1.0, 2.0, 50.0],
        )
        read = read_confirmation_detection(det, label="t")
        assert read.num_ions == 3
        assert read.num_scored == 2
        assert read.trapped_frac == pytest.approx(1.0 / 3.0)
        # the retained ion's verbatim handover state leaks into no aggregate
        assert read.n_mean == pytest.approx(2.5)
        assert float(read.fraction.sum()) == pytest.approx(1.0)
        assert read.fraction[21] == 0.0

    def test_suppressed_scores_in_bin_zero(self):
        det = make_detection(
            n_detected=[21.0, 4.0],
            state_reason=["suppressed", "frozen"],
            ke_eV=[3.0, 1.0],
        )
        read = read_confirmation_detection(det, label="t")
        assert read.suppressed_frac == pytest.approx(0.5)
        assert read.fraction[0] == pytest.approx(0.5)
        assert read.fraction[21] == 0.0
        assert read.n_mean == pytest.approx(2.0)  # (0 + 4) / 2

    def test_n1_fraction_and_histogram_support(self):
        det = make_detection(
            n_detected=[1.0, 1.0, 2.0, 5.0],
            state_reason=["frozen"] * 4,
            ke_eV=[1.0] * 4,
        )
        read = read_confirmation_detection(det, label="t")
        assert read.n_values.tolist() == list(range(22))
        assert read.n1_frac == pytest.approx(0.5)

    def test_all_retained_raises(self):
        det = make_detection(
            n_detected=[21.0], state_reason=["droplet_retained"], ke_eV=[9.0]
        )
        with pytest.raises(ValueError, match="droplet_retained"):
            read_confirmation_detection(det, label="t")

    def test_out_of_range_n_raises(self):
        det = make_detection(
            n_detected=[25.0], state_reason=["frozen"], ke_eV=[1.0]
        )
        with pytest.raises(ValueError, match="n_max"):
            read_confirmation_detection(det, label="t")

    def test_non_integer_n_raises(self):
        det = make_detection(
            n_detected=[2.5], state_reason=["frozen"], ke_eV=[1.0]
        )
        with pytest.raises(ValueError, match="int"):
            read_confirmation_detection(det, label="t")

    def test_ke_by_n_mean_and_se(self):
        det = make_detection(
            n_detected=[1.0, 1.0, 2.0],
            state_reason=["frozen"] * 3,
            ke_eV=[1.0, 3.0, 5.0],
        )
        read = read_confirmation_detection(det, label="t")
        ke = read.ke_by_n()
        assert ke.n.tolist() == [1, 2]
        assert ke.mean_eV[0] == pytest.approx(2.0)
        # ddof=1 sample SE: std([1,3], ddof=1)/sqrt(2) = sqrt(2)/sqrt(2) = 1
        assert ke.se_eV[0] == pytest.approx(1.0)
        # single-member bin: SE defined as 0 (the fragment_mean_kinetic_energy
        # precedent)
        assert ke.se_eV[1] == 0.0
        assert ke.count.tolist() == [2, 1]

    def test_suppressed_ke_lands_in_bin_zero(self):
        det = make_detection(
            n_detected=[21.0, 21.0],
            state_reason=["suppressed", "suppressed"],
            ke_eV=[2.0, 4.0],
        )
        read = read_confirmation_detection(det, label="t")
        ke = read.ke_by_n()
        assert ke.n.tolist() == [0]
        assert ke.mean_eV[0] == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# twin CSV loaders
# ---------------------------------------------------------------------------

class TestTwinLoaders:
    def test_prediction_row_round_trip(self, tmp_path):
        h = [0.5, 0.25, 0.25] + [0.0] * 19
        csv_path = write_twin_prediction_csv(
            tmp_path / "pred.csv",
            [
                {"leg": "d", "config": "c1", "h": h, "nbar_det": 0.75,
                 "trapped_frac": 0.1, "suppressed_frac": 0.5},
                {"leg": "d_delta", "config": "c1", "h": h, "nbar_det": 0.75},
            ],
        )
        row = load_twin_prediction(csv_path, leg="d", config="c1")
        assert row.leg == "d"
        assert row.config == "c1"
        assert row.h.shape == (22,)
        assert row.h[0] == pytest.approx(0.5)
        assert float(row.h.sum()) == pytest.approx(1.0)
        assert row.nbar_det == pytest.approx(0.75)
        assert row.trapped_frac == pytest.approx(0.1)
        assert row.suppressed_frac == pytest.approx(0.5)

    def test_prediction_rounding_renormalized(self, tmp_path):
        # 4-decimal CSV rounding: sum 0.9997 is accepted and renormalized
        h = [0.4999, 0.4998] + [0.0] * 20
        csv_path = write_twin_prediction_csv(
            tmp_path / "pred.csv",
            [{"leg": "d", "config": "c1", "h": h, "nbar_det": 0.5}],
        )
        row = load_twin_prediction(csv_path, leg="d", config="c1")
        assert float(row.h.sum()) == pytest.approx(1.0, abs=1e-12)

    def test_prediction_missing_config_raises(self, tmp_path):
        csv_path = write_twin_prediction_csv(
            tmp_path / "pred.csv",
            [{"leg": "d", "config": "c1", "h": [1.0] + [0.0] * 21,
              "nbar_det": 0.0}],
        )
        with pytest.raises(ValueError, match="c9"):
            load_twin_prediction(csv_path, leg="d", config="c9")

    def test_prediction_grossly_unnormalized_raises(self, tmp_path):
        h = [0.5, 0.4] + [0.0] * 20  # sum 0.9 — beyond rounding tolerance
        csv_path = write_twin_prediction_csv(
            tmp_path / "pred.csv",
            [{"leg": "d", "config": "c1", "h": h, "nbar_det": 0.5}],
        )
        with pytest.raises(ValueError, match="sum"):
            load_twin_prediction(csv_path, leg="d", config="c1")

    def test_ke_curve_loader(self, tmp_path):
        csv_path = write_twin_ke_csv(
            tmp_path / "ke.csv",
            [
                ("d", "c1", 0, 0.5, 1.34),
                ("d", "c1", 1, 0.5, 0.99),
                ("d", "c2", 1, 1.0, 0.90),
            ],
        )
        curve = load_twin_ke_curve(csv_path, leg="d", config="c1")
        assert curve.n.tolist() == [0, 1]
        assert curve.mean_ke_for(1) == pytest.approx(0.99)
        assert curve.mean_ke_for(7) is None

    def test_ke_curve_missing_selection_raises(self, tmp_path):
        csv_path = write_twin_ke_csv(
            tmp_path / "ke.csv", [("d", "c1", 0, 1.0, 1.0)]
        )
        with pytest.raises(ValueError, match="c2"):
            load_twin_ke_curve(csv_path, leg="d", config="c2")


# ---------------------------------------------------------------------------
# W1 adapter + solvated renormalization
# ---------------------------------------------------------------------------

class TestWassersteinBetween:
    def test_identity_is_zero(self):
        n = np.arange(4)
        f = np.array([0.25, 0.25, 0.25, 0.25])
        assert wasserstein_between(n, f, n, f) == 0.0

    def test_one_bin_shift_is_mass(self):
        # moving 0.5 mass one rung: W1 = 0.5
        n = np.arange(3)
        a = np.array([0.5, 0.5, 0.0])
        b = np.array([0.5, 0.0, 0.5])
        assert wasserstein_between(n, a, n, b) == pytest.approx(0.5)

    def test_matches_manual_cdf_sum(self):
        n = np.arange(5)
        a = np.array([0.1, 0.2, 0.3, 0.2, 0.2])
        b = np.array([0.3, 0.3, 0.2, 0.1, 0.1])
        manual = float(np.abs(np.cumsum(a - b)).sum())
        assert wasserstein_between(n, a, n, b) == pytest.approx(manual)


class TestSolvatedRenormalized:
    def test_drops_bare_and_renormalizes(self):
        n = np.arange(4)
        f = np.array([0.5, 0.25, 0.15, 0.10])
        n_out, f_out = solvated_renormalized(n, f)
        assert n_out.tolist() == [1, 2, 3]
        assert f_out[0] == pytest.approx(0.5)
        assert float(f_out.sum()) == pytest.approx(1.0)

    def test_zero_solvated_mass_raises(self):
        n = np.arange(3)
        f = np.array([1.0, 0.0, 0.0])
        with pytest.raises(ValueError, match="solvated"):
            solvated_renormalized(n, f)


# ---------------------------------------------------------------------------
# histogram score vs the abundance reference
# ---------------------------------------------------------------------------

class TestScoreHistogramVsReference:
    def test_hand_oracle(self):
        det = make_detection(
            n_detected=[0.0, 1.0, 1.0, 2.0],
            state_reason=["frozen"] * 4,
            ke_eV=[1.0] * 4,
        )
        read = read_confirmation_detection(det, label="t")
        # reference: bare 0.5, n1 0.3, n2 0.2 -> solvated (0.6, 0.4)
        ref = make_abundance_reference([0, 1, 2], [0.5, 0.3, 0.2])
        score = score_histogram_vs_reference(read, ref)
        # sim solvated: n1 2/3, n2 1/3
        assert score.sim_n1 == pytest.approx(2.0 / 3.0)
        assert score.sim_n2 == pytest.approx(1.0 / 3.0)
        assert score.ref_n1 == pytest.approx(0.6)
        assert score.ref_n2 == pytest.approx(0.4)
        assert score.ratio_n1_over_n2 == pytest.approx(2.0)
        # W1 on support {1,2}: |2/3-0.6| = 1/15
        assert score.w1_solvated == pytest.approx(1.0 / 15.0)

    def test_real_reference_solvated_targets(self):
        # the H.2b targets: n>=1 renormalization of the committed abundance
        # reference gives n1 ~ 0.310, n2 ~ 0.142 (the 31 % / 14.2 % targets)
        from i2_helium_md.postprocess.abundance_loader import (
            load_he_abundance_reference,
        )
        path = (
            Path(__file__).resolve().parents[1]
            / "data" / "reference" / "integrated_i_he_abundance.csv"
        )
        if not path.exists():
            pytest.skip("committed abundance reference not present")
        ref = load_he_abundance_reference(path)
        n_out, f_out = solvated_renormalized(ref.n, ref.ion_fraction)
        assert f_out[n_out.tolist().index(1)] == pytest.approx(0.310, abs=0.002)
        assert f_out[n_out.tolist().index(2)] == pytest.approx(0.142, abs=0.002)


# ---------------------------------------------------------------------------
# KE-curve score — the committed error model
# ---------------------------------------------------------------------------

class TestScoreKeCurveVsReference:
    def test_perfect_match_scores_zero(self):
        det = make_detection(
            n_detected=[1.0, 2.0],
            state_reason=["frozen"] * 2,
            ke_eV=[1.3, 0.7],
        )
        read = read_confirmation_detection(det, label="t")
        ref = make_ked_reference([0, 1, 2], [3.7, 1.3, 0.7])
        score = score_ke_curve_vs_reference(read, ref)
        assert score.n.tolist() == [1, 2]  # bare excluded by n_min=1
        assert score.chi2_profiled == pytest.approx(0.0, abs=1e-24)
        assert score.a_calib == pytest.approx(0.0, abs=1e-12)
        assert score.b_condition == pytest.approx(0.0, abs=1e-12)

    def test_bands_off_reduces_to_plain_chi2(self):
        det = make_detection(
            n_detected=[1.0, 2.0],
            state_reason=["frozen"] * 2,
            ke_eV=[1.4, 0.8],
        )
        read = read_confirmation_detection(det, label="t")
        ref = make_ked_reference(
            [1, 2], [1.3, 0.7],
            stat_err=0.1, calib_frac=0.0, condition_frac=0.0,
        )
        score = score_ke_curve_vs_reference(read, ref, include_sim_se=False)
        expected = (0.1 / 0.1) ** 2 + (0.1 / 0.1) ** 2
        assert score.chi2_profiled == pytest.approx(expected)
        assert score.chi2_unprofiled == pytest.approx(expected)
        assert score.a_calib == pytest.approx(0.0, abs=1e-12)

    def test_single_band_hand_oracle(self):
        # one point, one active band (condition off): the profiled optimum is
        # a = w*u*rho / (w*u^2 + 1), chi2 = w*(rho - a*u)^2 + a^2 with
        # u = calib_frac * r, w = 1/sigma^2.
        det = make_detection(
            n_detected=[1.0], state_reason=["frozen"], ke_eV=[1.4]
        )
        read = read_confirmation_detection(det, label="t")
        sigma = 0.05
        r, s, c = 1.3, 1.4, 0.04
        ref = make_ked_reference(
            [1], [r], stat_err=sigma, calib_frac=c, condition_frac=0.0
        )
        score = score_ke_curve_vs_reference(read, ref, include_sim_se=False)
        w = 1.0 / sigma**2
        u = c * r
        rho = s - r
        a_expected = w * u * rho / (w * u * u + 1.0)
        chi2_expected = w * (rho - a_expected * u) ** 2 + a_expected**2
        assert score.a_calib == pytest.approx(a_expected)
        assert score.chi2_profiled == pytest.approx(chi2_expected)
        # profiling can only lower chi2
        assert score.chi2_profiled < score.chi2_unprofiled

    def test_coherent_shift_absorbed_by_bands(self):
        # sim = ref * 1.04 exactly: a 4 % coherent lift is one calib-band
        # sigma, so the profiled chi2 must be far below the unprofiled one
        det = make_detection(
            n_detected=[1.0, 2.0, 3.0],
            state_reason=["frozen"] * 3,
            ke_eV=[1.3 * 1.04, 0.7 * 1.04, 0.5 * 1.04],
        )
        read = read_confirmation_detection(det, label="t")
        ref = make_ked_reference([1, 2, 3], [1.3, 0.7, 0.5], stat_err=0.005)
        score = score_ke_curve_vs_reference(read, ref, include_sim_se=False)
        assert score.chi2_profiled < 0.15 * score.chi2_unprofiled
        # the pull should land mostly on the (larger-lever) shift pair with
        # total shift ~ +4 %: a*0.04 + b*0.06 ~ 0.04
        total_shift = score.a_calib * 0.04 + score.b_condition * 0.06
        assert total_shift == pytest.approx(0.04, abs=0.01)

    def test_within_1p25_fallback_counter(self):
        det = make_detection(
            n_detected=[1.0, 2.0],
            state_reason=["frozen"] * 2,
            ke_eV=[1.3 * 1.2, 0.7 * 1.3],
        )
        read = read_confirmation_detection(det, label="t")
        ref = make_ked_reference([1, 2], [1.3, 0.7])
        score = score_ke_curve_vs_reference(read, ref)
        assert score.n_points == 2
        assert score.n_within_1p25 == 1

    def test_min_count_filters_thin_bins(self):
        det = make_detection(
            n_detected=[1.0, 1.0, 2.0],
            state_reason=["frozen"] * 3,
            ke_eV=[1.3, 1.3, 0.7],
        )
        read = read_confirmation_detection(det, label="t")
        ref = make_ked_reference([1, 2], [1.3, 0.7])
        score = score_ke_curve_vs_reference(read, ref, min_count=2)
        assert score.n.tolist() == [1]

    def test_no_scorable_bins_raises(self):
        det = make_detection(
            n_detected=[20.0], state_reason=["frozen"], ke_eV=[0.1]
        )
        read = read_confirmation_detection(det, label="t")
        ref = make_ked_reference([1, 2], [1.3, 0.7])
        with pytest.raises(ValueError, match="scorable"):
            score_ke_curve_vs_reference(read, ref)

    def test_sim_se_widens_sigma(self):
        det = make_detection(
            n_detected=[1.0, 1.0],
            state_reason=["frozen"] * 2,
            ke_eV=[1.0, 2.0],  # spread -> nonzero SE
        )
        read = read_confirmation_detection(det, label="t")
        ref = make_ked_reference(
            [1], [1.3], stat_err=0.1, calib_frac=0.0, condition_frac=0.0
        )
        with_se = score_ke_curve_vs_reference(read, ref, include_sim_se=True)
        without = score_ke_curve_vs_reference(read, ref, include_sim_se=False)
        assert with_se.chi2_profiled < without.chi2_profiled
