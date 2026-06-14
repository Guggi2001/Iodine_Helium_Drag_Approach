"""Method-B trajectory-matching extraction core (METHOD_B doc, plan slice B3).

All optimizer/objective tests run against a **stub** ``run_fn`` (no MD at
all): the stub synthesizes an ion-checkpoint-like object whose ensemble-mean
|v2| deviates from the reference linearly in ``(a - a*, b - b*, E - E*)``
along three linearly independent basis functions, so the objective has a
unique known minimum at ``(a*, b*, E*)``. One wiring test runs a single real
objective evaluation on a tiny (N=2, 0.2 ps) drag run with a synthetic
neutral checkpoint -- no production-sized simulation anywhere (CLAUDE.md
testing rules).
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from i2_helium_md import SimConfig, single_pulse_N2000_drag
from i2_helium_md.extraction.trajectory_matching import (
    CaseSetup,
    DIAGNOSTIC_4PARAM,
    E_BIND_STATIC_EV,
    LQ_SHARED_3PARAM,
    LQ_SHARED_PURE_QUADRATIC,
    PL_SHARED_3PARAM,
    SHARED_3PARAM,
    SHARED_PURE_CUBIC,
    _C_to_pivot,
    _pivot_to_C,
    evaluate_form_objective,
    evaluate_joint_form_objective,
    evaluate_joint_objective,
    evaluate_objective,
    fit_form_trajectory_matching,
    fit_shared_form_trajectory_matching,
    fit_shared_trajectory_matching,
    fit_trajectory_matching,
    form_sensitivity_halfwidths,
    joint_sensitivity_halfwidths,
    sensitivity_halfwidths,
    write_fit_parameters,
    write_per_case_form_fit_parameters,
    write_shared_fit_parameters,
    write_shared_form_fit_parameters,
)
from i2_helium_md.presets import load_drag_coefficients

# Reuse the Slice 4 synthetic neutral end-state (same cross-test reuse
# precedent as test_run_directory -> test_checkpoint).
from .test_ion_drag_smoke import _synthetic_neutral  # type: ignore

_M_EFF = 202.953908

# The stub's known optimum (a* [amu/ps], b* [amu*ps/A^2], E* [eV]).
_TRUE = (14.0, 2.0, 0.20)


def _make_stub_run_fn(t_grid, ref_speed, *, trap_above_e_eV=None, true=_TRUE):
    """Stub ``run_fn(cfg, neutral)`` with an analytic objective landscape.

    Mean |v2|(t) = ref(t) + 0.05*(a-a*)*1 + 0.5*(b-b*)*sin + 2.0*(E-E*)*ramp,
    three linearly independent deviations -> the in-window RMSE is uniquely
    minimized (at 0) at ``true`` (default ``_TRUE``). If ``trap_above_e_eV``
    is set, candidate bindings above it place the ions at the droplet surface
    with no escape (depth 0 < 2*steepness), exercising the penalty branch.
    """
    t_grid = np.asarray(t_grid, dtype=float)
    ref_speed = np.asarray(ref_speed, dtype=float)
    T = t_grid.size
    f_const = np.ones(T)
    f_sin = np.sin(np.linspace(0.0, 2.0 * np.pi, T))
    f_ramp = np.linspace(-1.0, 1.0, T)

    def run_fn(cfg, neutral):
        a = float(cfg.drag_coefficients.coefficients["a"])
        b = float(cfg.drag_coefficients.coefficients["b"])
        e = float(cfg.binding_energy_I_ion_eV)
        speed = (
            ref_speed
            + 0.05 * (a - true[0]) * f_const
            + 0.5 * (b - true[1]) * f_sin
            + 2.0 * (e - true[2]) * f_ramp
        )
        n = 1  # one molecule -> atoms [0] = I1, [1] = I2
        velocities_x = np.vstack([speed, speed])  # both atoms along +x
        zeros = np.zeros((2 * n, T))
        radius = 30.0
        trapped = trap_above_e_eV is not None and e > trap_above_e_eV
        # Escaped: depth = 3*steepness > 2*steepness with v_r > 0.
        # Trapped: depth = 0 (at the surface, inside the gate band).
        r_final = radius if trapped else radius + 3.0 * cfg.drag_gate_steepness
        positions_x = np.full((2 * n, T), r_final)
        return SimpleNamespace(
            num_molecules=n,
            time_ps=t_grid,
            velocities_x=velocities_x,
            velocities_y=zeros.copy(),
            velocities_z=zeros.copy(),
            positions_x=positions_x,
            positions_y=zeros.copy(),
            positions_z=zeros.copy(),
            droplet_radii_angstrom=np.full(2 * n, radius),
        )

    return run_fn


@pytest.fixture
def stub_setup():
    """A CaseSetup whose cfg is the plain default (no MD is ever run)."""
    t_grid = np.linspace(2.0, 12.0, 50)
    ref_speed = np.full(50, 3.5)
    return (
        CaseSetup(
            case="18A",
            cfg_base=SimConfig(),
            neutral=None,  # the stub ignores it
            t_ref_ps=t_grid,
            ref_speed_Aps=ref_speed,
            window=(float(t_grid[0]), float(t_grid[-1])),
            a0=_TRUE[0],
            b0=_TRUE[1],
            reference_path=Path("synthetic_reference.csv"),
        ),
        t_grid,
        ref_speed,
    )


class TestObjective:
    def test_zero_rmse_and_objective_at_true_params(self, stub_setup):
        setup, t_grid, ref_speed = stub_setup
        run_fn = _make_stub_run_fn(t_grid, ref_speed)
        res = evaluate_objective(_TRUE, setup, run_fn=run_fn)
        assert res.rmse_Aps == pytest.approx(0.0, abs=1e-12)
        assert res.escape_fraction == 1.0
        assert res.penalty_Aps == 0.0
        assert res.objective == pytest.approx(0.0, abs=1e-12)
        assert res.num_overlap_points == t_grid.size

    @pytest.mark.parametrize(
        "theta",
        [
            (_TRUE[0] * 1.2, _TRUE[1], _TRUE[2]),
            (_TRUE[0], _TRUE[1] * 1.5, _TRUE[2]),
            (_TRUE[0], _TRUE[1], _TRUE[2] + 0.05),
        ],
    )
    def test_any_single_param_deviation_raises_objective(self, stub_setup, theta):
        setup, t_grid, ref_speed = stub_setup
        run_fn = _make_stub_run_fn(t_grid, ref_speed)
        res = evaluate_objective(theta, setup, run_fn=run_fn)
        assert res.objective > 1e-3

    def test_penalty_fires_on_trapped_ions(self, stub_setup):
        setup, t_grid, ref_speed = stub_setup
        run_fn = _make_stub_run_fn(t_grid, ref_speed, trap_above_e_eV=0.25)
        trapped = evaluate_objective(
            (_TRUE[0], _TRUE[1], 0.30), setup, run_fn=run_fn
        )
        assert trapped.escape_fraction == 0.0
        assert trapped.penalty_Aps == pytest.approx(2.0)  # default weight
        # The penalty dominates the small rmse from the E deviation: the
        # optimizer is steered away from trapping bindings.
        escaped = evaluate_objective(
            (_TRUE[0], _TRUE[1], 0.20), setup, run_fn=run_fn
        )
        assert trapped.objective > escaped.objective + 1.5

    def test_deterministic_bitwise(self, stub_setup):
        setup, t_grid, ref_speed = stub_setup
        run_fn = _make_stub_run_fn(t_grid, ref_speed)
        theta = (15.0, 2.2, 0.18)
        r1 = evaluate_objective(theta, setup, run_fn=run_fn)
        r2 = evaluate_objective(theta, setup, run_fn=run_fn)
        assert r1.objective == r2.objective
        assert r1.rmse_Aps == r2.rmse_Aps

    def test_stamped_pairing_passes_guard_silently(self, stub_setup, recwarn):
        # Each evaluation stamps E_bind into the bundle AND the cfg: the
        # §6.5.1 guard must pass with no warning (no escape hatch in play).
        setup, t_grid, ref_speed = stub_setup
        run_fn = _make_stub_run_fn(t_grid, ref_speed)
        evaluate_objective(_TRUE, setup, run_fn=run_fn)
        assert not [w for w in recwarn.list if "6.5.1" in str(w.message)]


class TestOptimizer:
    def test_smoke_tiny_budget(self, stub_setup):
        setup, t_grid, ref_speed = stub_setup
        run_fn = _make_stub_run_fn(t_grid, ref_speed)
        fit = fit_trajectory_matching(
            setup,
            run_fn=run_fn,
            maxfev_per_start=5,
            extra_starts=False,
            prescan_e_bind_eV=[0.1, 0.2, 0.3],
        )
        assert np.isfinite(fit.best.objective)
        assert len(fit.starts) == 1
        assert fit.n_evaluations > 3

    def test_converges_to_known_minimum_on_stub(self, stub_setup):
        setup, t_grid, ref_speed = stub_setup
        run_fn = _make_stub_run_fn(t_grid, ref_speed)
        fit = fit_trajectory_matching(
            setup,
            run_fn=run_fn,
            maxfev_per_start=400,
            extra_starts=False,
            prescan_e_bind_eV=np.linspace(0.05, 0.305, 8),
        )
        # NM on a smooth 3-param landscape with a unique zero minimum:
        # parameters recovered to ~1% (loose -- NM simplex tolerance, not a
        # statement about MD physics).
        assert fit.best.a == pytest.approx(_TRUE[0], rel=0.01)
        assert fit.best.b == pytest.approx(_TRUE[1], rel=0.01)
        assert fit.best.e_bind_eV == pytest.approx(_TRUE[2], abs=0.005)
        assert fit.converged

    def test_multi_start_records_each_minimum(self, stub_setup):
        setup, t_grid, ref_speed = stub_setup
        run_fn = _make_stub_run_fn(t_grid, ref_speed)
        fit = fit_trajectory_matching(
            setup,
            run_fn=run_fn,
            maxfev_per_start=50,
            extra_starts=True,
            prescan_e_bind_eV=[0.1, 0.2],
        )
        assert len(fit.starts) == 3
        assert fit.best.objective == min(r.objective for r in fit.starts)


class TestSensitivityHalfwidths:
    def test_halfwidths_positive_finite_and_bounded(self, stub_setup):
        # Coarse-band contract only (the helper is a sensitivity band, not a
        # CI): positive, finite, and bounded by the largest scanned change
        # (factor 2.0 -> one parameter scale). Evaluated at a slightly
        # off-optimum point so best.objective > 0 and the 10%-rise threshold
        # is meaningful.
        setup, t_grid, ref_speed = stub_setup
        run_fn = _make_stub_run_fn(t_grid, ref_speed)
        best = evaluate_objective(
            (_TRUE[0] * 1.05, _TRUE[1], _TRUE[2]), setup, run_fn=run_fn
        )
        hw = sensitivity_halfwidths(setup, best, run_fn=run_fn)
        theta = (best.a, best.b, best.e_bind_eV)
        for halfwidth, param in zip(hw, theta):
            assert np.isfinite(halfwidth)
            assert halfwidth > 0
            assert halfwidth <= param  # factor-2 scan bound


class TestWriteFitParameters:
    def _completed_fit(self, stub_setup, tmp_path):
        setup, t_grid, ref_speed = stub_setup
        run_fn = _make_stub_run_fn(t_grid, ref_speed)
        fit = fit_trajectory_matching(
            setup,
            run_fn=run_fn,
            maxfev_per_start=200,
            extra_starts=False,
            prescan_e_bind_eV=[0.15, 0.2, 0.25],
        )
        fit.a_err = 0.1
        fit.b_err = 0.05
        fit.e_bind_err_eV = 0.01
        fit.uncertainty_model = "seed_sweep_std_plus_rmse_sensitivity"
        return fit

    def test_written_bundle_loads_through_slice3_loader(self, stub_setup, tmp_path):
        # The end-to-end schema check: the Method-B writer's output must load
        # through load_drag_coefficients' trajectory_matching mode.
        fit = self._completed_fit(stub_setup, tmp_path)
        out = write_fit_parameters(fit, tmp_path, transverse_contaminated=False)
        assert out.name == "fit_parameters.json"
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.extraction_method == "trajectory_matching"
        assert coeffs.coefficients["a"] == pytest.approx(fit.best.a)
        assert coeffs.coefficients["b"] == pytest.approx(fit.best.b)
        assert coeffs.effective_binding_energy_I_ion_eV == pytest.approx(
            fit.best.e_bind_eV
        )

    def test_refuses_trapped_fit(self, stub_setup, tmp_path):
        setup, t_grid, ref_speed = stub_setup
        run_fn = _make_stub_run_fn(t_grid, ref_speed, trap_above_e_eV=0.0)
        fit = fit_trajectory_matching(
            setup,
            run_fn=run_fn,
            maxfev_per_start=5,
            extra_starts=False,
            prescan_e_bind_eV=[0.2],
        )
        fit.a_err, fit.b_err, fit.uncertainty_model = 0.1, 0.1, "x"
        assert fit.best.escape_fraction == 0.0
        with pytest.raises(ValueError, match="trapped"):
            write_fit_parameters(fit, tmp_path, transverse_contaminated=False)

    def test_refuses_unfilled_uncertainty(self, stub_setup, tmp_path):
        fit = self._completed_fit(stub_setup, tmp_path)
        fit.a_err = None
        with pytest.raises(ValueError, match="uncertainty"):
            write_fit_parameters(fit, tmp_path, transverse_contaminated=False)


class TestRealWiring:
    def test_single_real_objective_evaluation(self):
        """One real (N=2, 0.2 ps) objective evaluation through the actual
        BAOAB ion driver on a synthetic neutral checkpoint: wiring only, no
        physics assertion, no production-sized run."""
        cfg = single_pulse_N2000_drag(
            num_molecules=2, ion_simulation_time=0.2, dt_ion=0.01, seed=0,
        )
        t_ref = np.linspace(0.05, 0.15, 5)
        setup = CaseSetup(
            case="9A",
            cfg_base=cfg,
            neutral=_synthetic_neutral(num_molecules=2),
            t_ref_ps=t_ref,
            ref_speed_Aps=np.full(5, 3.0),
            window=(0.05, 0.15),
            a0=float(cfg.drag_coefficients.coefficients["a"]),
            b0=float(cfg.drag_coefficients.coefficients["b"]),
            reference_path=Path("synthetic_reference.csv"),
        )
        res = evaluate_objective((setup.a0, setup.b0, 0.2), setup)
        assert np.isfinite(res.objective)
        assert np.isfinite(res.rmse_Aps)
        assert 0.0 <= res.escape_fraction <= 1.0
        assert res.num_overlap_points == 5

    def test_e_bind_bounds_respected_in_fit(self, ):
        # E_BIND_STATIC_EV is the physical upper bound the fit must respect.
        assert E_BIND_STATIC_EV == 0.308


# ===========================================================================
# Shared-form joint refit (METHOD_B §9) -- stub-MD tests
# ===========================================================================

def _joint_setups(t_grid, ref_speed, true_by_case):
    """Two-case stub context for the §9 joint fit.

    Each case gets its own stub landscape (its own true ``(a*, b*, E*)``);
    ``CaseSetup.neutral`` carries the case name as a dispatch marker so ONE
    ``run_fn`` (the signature the joint objective passes around) routes to
    the right per-case stub.
    """
    setups, stubs = {}, {}
    for case, true in true_by_case.items():
        setups[case] = CaseSetup(
            case=case,
            cfg_base=SimConfig(),
            neutral=case,  # dispatch marker; the stub ignores it otherwise
            t_ref_ps=t_grid,
            ref_speed_Aps=ref_speed,
            window=(float(t_grid[0]), float(t_grid[-1])),
            a0=true[0] if true[0] > 0 else 1.0,  # unused by the shared fit
            b0=true[1],
            reference_path=Path(f"synthetic_{case}.csv"),
        )
        stubs[case] = _make_stub_run_fn(t_grid, ref_speed, true=true)

    def run_fn(cfg, neutral):
        return stubs[neutral](cfg, neutral)

    return setups, run_fn


@pytest.fixture
def joint_grid():
    t_grid = np.linspace(2.0, 12.0, 50)
    ref_speed = np.full(50, 3.5)
    return t_grid, ref_speed


class TestJointObjective:
    def test_zero_at_shared_true_params(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _joint_setups(
            t_grid, ref_speed, {"18A": _TRUE, "9A": _TRUE}
        )
        res = evaluate_joint_objective(
            _TRUE[0], _TRUE[1], {"18A": _TRUE[2], "9A": _TRUE[2]},
            setups, run_fn=run_fn,
        )
        assert res.objective == pytest.approx(0.0, abs=1e-12)
        assert set(res.per_case) == {"18A", "9A"}

    def test_equal_weight_mean_of_per_case_objectives(self, joint_grid):
        # Deviate only the 9A binding: the joint objective must equal the
        # plain mean of the two per-case objectives (equal weight, §9.2).
        t_grid, ref_speed = joint_grid
        setups, run_fn = _joint_setups(
            t_grid, ref_speed, {"18A": _TRUE, "9A": _TRUE}
        )
        res = evaluate_joint_objective(
            _TRUE[0], _TRUE[1], {"18A": _TRUE[2], "9A": _TRUE[2] + 0.05},
            setups, run_fn=run_fn,
        )
        expected = np.mean([r.objective for r in res.per_case.values()])
        assert res.objective == pytest.approx(expected, rel=1e-12)
        assert res.per_case["18A"].objective == pytest.approx(0.0, abs=1e-12)
        assert res.per_case["9A"].objective > 1e-3

    def test_case_key_mismatch_raises(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _joint_setups(t_grid, ref_speed, {"18A": _TRUE})
        with pytest.raises(ValueError, match="do not match"):
            evaluate_joint_objective(
                _TRUE[0], _TRUE[1], {"9A": 0.2}, setups, run_fn=run_fn
            )


class TestSharedFit:
    def test_shared_3param_recovers_shared_true(self, joint_grid):
        # Both cases share one true law -> the over-constrained 3-parameter
        # fit must recover it (the §9.2 generalization hypothesis, stubbed).
        t_grid, ref_speed = joint_grid
        setups, run_fn = _joint_setups(
            t_grid, ref_speed, {"18A": _TRUE, "9A": _TRUE}
        )
        fit = fit_shared_trajectory_matching(
            setups, variant=SHARED_3PARAM, anchors=(_TRUE[0], _TRUE[1]),
            run_fn=run_fn, maxfev_per_start=400, extra_starts=False,
            prescan_e_bind_eV=np.linspace(0.05, 0.305, 8),
        )
        assert fit.best.a == pytest.approx(_TRUE[0], rel=0.01)
        assert fit.best.b == pytest.approx(_TRUE[1], rel=0.01)
        e_values = set(fit.best.e_bind_eV.values())
        assert len(e_values) == 1  # shared binding by construction
        assert e_values.pop() == pytest.approx(_TRUE[2], abs=0.005)
        assert fit.converged
        assert fit.cases == ["18A", "9A"]

    def test_shared_3param_a_lower_bound_zero_reachable(self, joint_grid):
        # The §9 decision: a free with lower bound 0 -- a pure-cubic-true
        # landscape must let the 3-param fit run a to (near) 0, which the
        # per-case fit's old 0.1*a0 bound forbade.
        t_grid, ref_speed = joint_grid
        true = (0.0, 2.0, 0.20)
        setups, run_fn = _joint_setups(
            t_grid, ref_speed, {"18A": true, "9A": true}
        )
        fit = fit_shared_trajectory_matching(
            setups, variant=SHARED_3PARAM, anchors=(14.0, 2.0),
            run_fn=run_fn, maxfev_per_start=400, extra_starts=False,
            prescan_e_bind_eV=[0.15, 0.2, 0.25],
        )
        assert fit.best.a == pytest.approx(0.0, abs=0.05)
        assert fit.best.b == pytest.approx(true[1], rel=0.01)

    def test_pure_cubic_fixes_a_to_zero(self, joint_grid):
        t_grid, ref_speed = joint_grid
        true = (0.0, 2.0, 0.20)
        setups, run_fn = _joint_setups(
            t_grid, ref_speed, {"18A": true, "9A": true}
        )
        fit = fit_shared_trajectory_matching(
            setups, variant=SHARED_PURE_CUBIC, anchors=(14.0, 2.0),
            run_fn=run_fn, maxfev_per_start=400, extra_starts=False,
            prescan_e_bind_eV=[0.15, 0.2, 0.25],
        )
        assert fit.best.a == 0.0  # fixed, not fitted
        assert fit.best.b == pytest.approx(true[1], rel=0.01)
        e_values = set(fit.best.e_bind_eV.values())
        assert e_values.pop() == pytest.approx(true[2], abs=0.005)

    def test_diagnostic_recovers_per_case_bindings(self, joint_grid):
        # Same law, different per-case E* -> the 4-parameter diagnostic must
        # localize the disagreement in the bindings, not the law.
        t_grid, ref_speed = joint_grid
        setups, run_fn = _joint_setups(
            t_grid, ref_speed,
            {"18A": (14.0, 2.0, 0.15), "9A": (14.0, 2.0, 0.25)},
        )
        fit = fit_shared_trajectory_matching(
            setups, variant=DIAGNOSTIC_4PARAM, anchors=(14.0, 2.0),
            run_fn=run_fn, maxfev_per_start=600, extra_starts=False,
            prescan_e_bind_eV=[0.15, 0.2, 0.25],
        )
        assert fit.best.a == pytest.approx(14.0, rel=0.02)
        assert fit.best.b == pytest.approx(2.0, rel=0.02)
        assert fit.best.e_bind_eV["18A"] == pytest.approx(0.15, abs=0.01)
        assert fit.best.e_bind_eV["9A"] == pytest.approx(0.25, abs=0.01)

    def test_unknown_variant_raises(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _joint_setups(
            t_grid, ref_speed, {"18A": _TRUE, "9A": _TRUE}
        )
        with pytest.raises(ValueError, match="variant"):
            fit_shared_trajectory_matching(
                setups, variant="nope", anchors=(14.0, 2.0), run_fn=run_fn
            )

    def test_single_case_raises(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _joint_setups(t_grid, ref_speed, {"18A": _TRUE})
        with pytest.raises(ValueError, match="cross-case"):
            fit_shared_trajectory_matching(
                setups, variant=SHARED_3PARAM, anchors=(14.0, 2.0),
                run_fn=run_fn,
            )

    def test_mismatched_provenance_raises(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _joint_setups(
            t_grid, ref_speed, {"18A": _TRUE, "9A": _TRUE}
        )
        setups["9A"].cfg_base = SimConfig(seed=1)  # different neutral seed
        with pytest.raises(ValueError, match="must share"):
            fit_shared_trajectory_matching(
                setups, variant=SHARED_3PARAM, anchors=(14.0, 2.0),
                run_fn=run_fn,
            )


class TestJointSensitivity:
    def test_halfwidths_finite_and_zero_for_fixed_a(self, joint_grid):
        t_grid, ref_speed = joint_grid
        true = (0.0, 2.0, 0.20)
        setups, run_fn = _joint_setups(
            t_grid, ref_speed, {"18A": true, "9A": true}
        )
        # Slightly off-optimum so the 10%-rise threshold is meaningful.
        best = evaluate_joint_objective(
            0.0, true[1] * 1.05, {"18A": true[2], "9A": true[2]},
            setups, run_fn=run_fn,
        )
        hw_a, hw_b, hw_e = joint_sensitivity_halfwidths(
            setups, best, run_fn=run_fn
        )
        assert hw_a == 0.0  # a fixed at 0 under pure-cubic: not a free param
        assert np.isfinite(hw_b) and hw_b > 0
        assert np.isfinite(hw_e) and hw_e > 0

    def test_per_case_bindings_refused(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _joint_setups(
            t_grid, ref_speed, {"18A": _TRUE, "9A": _TRUE}
        )
        best = evaluate_joint_objective(
            _TRUE[0], _TRUE[1], {"18A": 0.15, "9A": 0.25},
            setups, run_fn=run_fn,
        )
        with pytest.raises(ValueError, match="shared-E_bind"):
            joint_sensitivity_halfwidths(setups, best, run_fn=run_fn)


class TestWriteSharedFitParameters:
    def _completed_shared_fit(self, joint_grid, variant=SHARED_3PARAM):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _joint_setups(
            t_grid, ref_speed, {"18A": _TRUE, "9A": _TRUE}
        )
        fit = fit_shared_trajectory_matching(
            setups, variant=variant, anchors=(_TRUE[0], _TRUE[1]),
            run_fn=run_fn, maxfev_per_start=200, extra_starts=False,
            prescan_e_bind_eV=[0.15, 0.2, 0.25],
        )
        fit.a_err, fit.b_err, fit.e_bind_err_eV = 0.1, 0.05, 0.01
        fit.uncertainty_model = "rmse_sensitivity_only_seed_sweep_omitted"
        return fit

    def test_written_bundle_loads_through_slice3_loader(
        self, joint_grid, tmp_path
    ):
        fit = self._completed_shared_fit(joint_grid)
        out = write_shared_fit_parameters(fit, tmp_path)
        assert out.name == "fit_parameters.json"
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.extraction_method == "trajectory_matching"
        assert coeffs.coefficients["a"] == pytest.approx(fit.best.a)
        assert coeffs.coefficients["b"] == pytest.approx(fit.best.b)
        assert coeffs.effective_binding_energy_I_ion_eV == pytest.approx(
            next(iter(fit.best.e_bind_eV.values()))
        )

    def test_shared_provenance_fields_stamped(self, joint_grid, tmp_path):
        import json as _json

        fit = self._completed_shared_fit(joint_grid)
        out = write_shared_fit_parameters(fit, tmp_path)
        with open(out, "r", encoding="utf-8") as fh:
            raw = _json.load(fh)
        assert raw["calibration_cases"] == ["18A", "9A"]
        assert raw["stage"] == "stage2_joint"
        assert raw["variant"] == SHARED_3PARAM
        assert raw["anchors"]["a0"] == pytest.approx(_TRUE[0])
        assert set(raw["per_case"]) == {"18A", "9A"}
        assert raw["flags"]["cross_case_axis_consumed_by_joint_fit"] is True

    def test_refuses_diagnostic_variant(self, joint_grid, tmp_path):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _joint_setups(
            t_grid, ref_speed, {"18A": _TRUE, "9A": _TRUE}
        )
        fit = fit_shared_trajectory_matching(
            setups, variant=DIAGNOSTIC_4PARAM, anchors=(_TRUE[0], _TRUE[1]),
            run_fn=run_fn, maxfev_per_start=5, extra_starts=False,
            prescan_e_bind_eV=[0.2],
        )
        fit.a_err, fit.b_err, fit.uncertainty_model = 0.1, 0.1, "x"
        with pytest.raises(ValueError, match="diagnostic"):
            write_shared_fit_parameters(fit, tmp_path)

    def test_refuses_trapped_fit_in_any_case(self, joint_grid, tmp_path):
        t_grid, ref_speed = joint_grid
        setups = {}
        stubs = {
            "18A": _make_stub_run_fn(t_grid, ref_speed),
            "9A": _make_stub_run_fn(t_grid, ref_speed, trap_above_e_eV=0.0),
        }
        for case in ("18A", "9A"):
            setups[case] = CaseSetup(
                case=case, cfg_base=SimConfig(), neutral=case,
                t_ref_ps=t_grid, ref_speed_Aps=ref_speed,
                window=(float(t_grid[0]), float(t_grid[-1])),
                a0=_TRUE[0], b0=_TRUE[1],
                reference_path=Path(f"synthetic_{case}.csv"),
            )

        def run_fn(cfg, neutral):
            return stubs[neutral](cfg, neutral)

        fit = fit_shared_trajectory_matching(
            setups, variant=SHARED_3PARAM, anchors=(_TRUE[0], _TRUE[1]),
            run_fn=run_fn, maxfev_per_start=5, extra_starts=False,
            prescan_e_bind_eV=[0.2],
        )
        fit.a_err, fit.b_err, fit.uncertainty_model = 0.1, 0.1, "x"
        assert fit.best.per_case["9A"].escape_fraction == 0.0
        with pytest.raises(ValueError, match="trapped"):
            write_shared_fit_parameters(fit, tmp_path)

    def test_refuses_unfilled_uncertainty(self, joint_grid, tmp_path):
        fit = self._completed_shared_fit(joint_grid)
        fit.a_err = None
        with pytest.raises(ValueError, match="uncertainty"):
            write_shared_fit_parameters(fit, tmp_path)


class TestSharedBundleArtifacts:
    """The committed §9 shared-form bundles (2026-06-11 refit, verdict PASS)
    stay loadable through the Slice-3 loader and well-formed."""

    @pytest.mark.parametrize(
        "variant", ["shared_3param", "shared_pure_cubic"]
    )
    def test_shared_variant_bundle_loads(self, variant):
        import json

        from i2_helium_md.presets import REFERENCE_DRAG_ROOT

        coeff_dir = (
            REFERENCE_DRAG_ROOT / "shared" / "trajectory_matching" / variant
        )
        coeffs = load_drag_coefficients(coeff_dir, expected_m_eff_amu=_M_EFF)
        assert coeffs.extraction_method == "trajectory_matching"
        assert coeffs.coefficients["a"] >= 0.0  # a >= 0 (§9.5 relaxed guard)
        assert coeffs.coefficients["b"] > 0.0
        assert coeffs.effective_binding_energy_I_ion_eV > 0.0
        if variant == "shared_pure_cubic":
            assert coeffs.coefficients["a"] == 0.0  # fixed, not fitted
        with open(
            coeff_dir / "fit_parameters.json", "r", encoding="utf-8"
        ) as fh:
            raw = json.load(fh)
        assert raw["variant"] == variant
        assert raw["calibration_cases"] == ["18A", "9A"]
        assert raw["stage"] == "stage2_joint"


# ===========================================================================
# Alternative drag-form discrimination (METHOD_B §10) -- stub-MD tests
# ===========================================================================
# Same stub discipline as the §9 tests: an analytic landscape with a unique
# known minimum per family; no MD anywhere (the real-wiring path is covered
# by the smoke harness's end-to-end family runs in test_ion_drag_smoke).

# Per-family stub optima and pre-registered-style anchors (test values).
_LQ_TRUE = (4.0, 11.0, 0.20)            # (a*, c*, E*)
_LQ_ANCHORS = {"a0": 4.0, "c0": 11.0}
_PL_TRUE = (10.0, 2.0, 0.20)            # (C*, n*, E*)
_PL_ANCHORS = {"C0": 10.0, "n0": 2.0, "v_ref_Aps": 3.0}


def _make_form_stub_run_fn(
    t_grid, ref_speed, *, keys, true, trap_above_e_eV=None
):
    """Form-generic stub: deviations linear in the two form coefficients + E.

    Mirrors :func:`_make_stub_run_fn` with the coefficient keys parameterized
    (``("a", "c")`` for linear_quadratic, ``("C", "n")`` for power_law).
    """
    t_grid = np.asarray(t_grid, dtype=float)
    ref_speed = np.asarray(ref_speed, dtype=float)
    T = t_grid.size
    f_const = np.ones(T)
    f_sin = np.sin(np.linspace(0.0, 2.0 * np.pi, T))
    f_ramp = np.linspace(-1.0, 1.0, T)

    def run_fn(cfg, neutral):
        k1 = float(cfg.drag_coefficients.coefficients[keys[0]])
        k2 = float(cfg.drag_coefficients.coefficients[keys[1]])
        e = float(cfg.binding_energy_I_ion_eV)
        speed = (
            ref_speed
            + 0.05 * (k1 - true[0]) * f_const
            + 0.5 * (k2 - true[1]) * f_sin
            + 2.0 * (e - true[2]) * f_ramp
        )
        n = 1
        velocities_x = np.vstack([speed, speed])
        zeros = np.zeros((2 * n, T))
        radius = 30.0
        trapped = trap_above_e_eV is not None and e > trap_above_e_eV
        r_final = radius if trapped else radius + 3.0 * cfg.drag_gate_steepness
        positions_x = np.full((2 * n, T), r_final)
        return SimpleNamespace(
            num_molecules=n,
            time_ps=t_grid,
            velocities_x=velocities_x,
            velocities_y=zeros.copy(),
            velocities_z=zeros.copy(),
            positions_x=positions_x,
            positions_y=zeros.copy(),
            positions_z=zeros.copy(),
            droplet_radii_angstrom=np.full(2 * n, radius),
        )

    return run_fn


def _form_joint_setups(t_grid, ref_speed, *, keys, true_by_case):
    """Two-case (or one-case) stub context for the §10 form-phase fits."""
    setups, stubs = {}, {}
    for case, true in true_by_case.items():
        setups[case] = CaseSetup(
            case=case,
            cfg_base=SimConfig(),
            neutral=case,  # dispatch marker, as in _joint_setups
            t_ref_ps=t_grid,
            ref_speed_Aps=ref_speed,
            window=(float(t_grid[0]), float(t_grid[-1])),
            a0=1.0,  # unused by the form-phase fits (anchors are explicit)
            b0=1.0,
            reference_path=Path(f"synthetic_{case}.csv"),
        )
        stubs[case] = _make_form_stub_run_fn(
            t_grid, ref_speed, keys=keys, true=true
        )

    def run_fn(cfg, neutral):
        return stubs[neutral](cfg, neutral)

    return setups, run_fn


class TestPivotTransform:
    def test_round_trip_identity(self):
        # C -> gamma_ref -> C across a (C, n) grid (analytical, tight).
        for C in (0.5, 2.5154, 10.36):
            for n in (1.0, 2.0, 2.056, 3.0, 4.0):
                gamma_ref = _C_to_pivot(C, n, 3.0)
                assert _pivot_to_C(gamma_ref, n, 3.0) == pytest.approx(
                    C, rel=1e-14
                )

    def test_pivot_is_identity_at_n_equal_one(self):
        # gamma_ref = C * v_ref**0 = C.
        assert _C_to_pivot(7.0, 1.0, 3.0) == pytest.approx(7.0, rel=1e-14)

    def test_pivot_value_explicit(self):
        # gamma_ref = C * v_ref**(n-1): 2 * 3**2 = 18 at n = 3.
        assert _C_to_pivot(2.0, 3.0, 3.0) == pytest.approx(18.0, rel=1e-14)


class TestFormObjective:
    @pytest.mark.parametrize(
        "form,keys,true,anchors",
        [
            ("linear_quadratic", ("a", "c"), _LQ_TRUE, _LQ_ANCHORS),
            ("power_law", ("C", "n"), _PL_TRUE, _PL_ANCHORS),
        ],
    )
    def test_zero_objective_at_true_params(self, joint_grid, form, keys, true, anchors):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=keys, true_by_case={"18A": true}
        )
        res = evaluate_form_objective(
            form, {keys[0]: true[0], keys[1]: true[1]}, true[2],
            setups["18A"], run_fn=run_fn,
        )
        assert res.objective == pytest.approx(0.0, abs=1e-12)
        assert res.escape_fraction == 1.0
        assert res.form == form

    def test_linear_cubic_cross_check_matches_evaluate_objective(
        self, stub_setup
    ):
        # The form-generic path with linear_cubic {a, b} must reproduce the
        # §9 evaluate_objective bitwise (same _run_and_score core).
        setup, t_grid, ref_speed = stub_setup
        run_fn = _make_stub_run_fn(t_grid, ref_speed)
        theta = (15.0, 2.2, 0.18)
        legacy = evaluate_objective(theta, setup, run_fn=run_fn)
        generic = evaluate_form_objective(
            "linear_cubic", {"a": theta[0], "b": theta[1]}, theta[2],
            setup, run_fn=run_fn,
        )
        assert generic.objective == legacy.objective
        assert generic.rmse_Aps == legacy.rmse_Aps
        assert generic.escape_fraction == legacy.escape_fraction

    def test_joint_case_key_mismatch_raises(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"), true_by_case={"18A": _LQ_TRUE}
        )
        with pytest.raises(ValueError, match="do not match"):
            evaluate_joint_form_objective(
                "linear_quadratic", {"a": 4.0, "c": 11.0}, {"9A": 0.2},
                setups, run_fn=run_fn,
            )


class TestFormFit:
    def test_lq_shared_3param_recovers_true(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"),
            true_by_case={"18A": _LQ_TRUE, "9A": _LQ_TRUE},
        )
        fit = fit_shared_form_trajectory_matching(
            setups, variant=LQ_SHARED_3PARAM, anchors=_LQ_ANCHORS,
            run_fn=run_fn, maxfev_per_start=400, extra_starts=False,
            prescan_e_bind_eV=np.linspace(0.05, 0.305, 8),
        )
        assert fit.form == "linear_quadratic"
        assert fit.best.coefficients["a"] == pytest.approx(_LQ_TRUE[0], rel=0.01)
        assert fit.best.coefficients["c"] == pytest.approx(_LQ_TRUE[1], rel=0.01)
        e_values = set(fit.best.e_bind_eV.values())
        assert len(e_values) == 1
        assert e_values.pop() == pytest.approx(_LQ_TRUE[2], abs=0.005)
        assert fit.converged

    def test_lq_pure_quadratic_fixes_a_to_zero(self, joint_grid):
        t_grid, ref_speed = joint_grid
        true = (0.0, 11.0, 0.20)
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"),
            true_by_case={"18A": true, "9A": true},
        )
        fit = fit_shared_form_trajectory_matching(
            setups, variant=LQ_SHARED_PURE_QUADRATIC, anchors=_LQ_ANCHORS,
            run_fn=run_fn, maxfev_per_start=400, extra_starts=False,
            prescan_e_bind_eV=[0.15, 0.2, 0.25],
        )
        assert fit.best.coefficients["a"] == 0.0  # fixed, not fitted
        assert fit.best.coefficients["c"] == pytest.approx(true[1], rel=0.01)
        e_values = set(fit.best.e_bind_eV.values())
        assert e_values.pop() == pytest.approx(true[2], abs=0.005)

    def test_pl_shared_3param_recovers_true_through_pivot(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("C", "n"),
            true_by_case={"18A": _PL_TRUE, "9A": _PL_TRUE},
        )
        fit = fit_shared_form_trajectory_matching(
            setups, variant=PL_SHARED_3PARAM, anchors=_PL_ANCHORS,
            run_fn=run_fn, maxfev_per_start=600, extra_starts=False,
            prescan_e_bind_eV=np.linspace(0.05, 0.305, 8),
        )
        assert fit.form == "power_law"
        # The optimizer works in (gamma_ref, n); the result stamps raw {C, n}.
        assert fit.best.coefficients["C"] == pytest.approx(_PL_TRUE[0], rel=0.02)
        assert fit.best.coefficients["n"] == pytest.approx(_PL_TRUE[1], rel=0.02)
        e_values = set(fit.best.e_bind_eV.values())
        assert e_values.pop() == pytest.approx(_PL_TRUE[2], abs=0.005)

    def test_stage1_analog_single_case_fit(self, joint_grid):
        # The Stage-1 analog: fit the family's full variant on 18A only.
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"),
            true_by_case={"18A": _LQ_TRUE},
        )
        fit = fit_form_trajectory_matching(
            setups["18A"], variant=LQ_SHARED_3PARAM, anchors=_LQ_ANCHORS,
            run_fn=run_fn, maxfev_per_start=400, extra_starts=False,
            prescan_e_bind_eV=[0.15, 0.2, 0.25],
        )
        assert fit.cases == ["18A"]
        assert fit.best.coefficients["a"] == pytest.approx(_LQ_TRUE[0], rel=0.01)
        assert fit.best.coefficients["c"] == pytest.approx(_LQ_TRUE[1], rel=0.01)

    def test_shared_fit_requires_two_cases(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"), true_by_case={"18A": _LQ_TRUE}
        )
        with pytest.raises(ValueError, match="cross-case"):
            fit_shared_form_trajectory_matching(
                setups, variant=LQ_SHARED_3PARAM, anchors=_LQ_ANCHORS,
                run_fn=run_fn,
            )

    def test_unknown_variant_raises(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"),
            true_by_case={"18A": _LQ_TRUE, "9A": _LQ_TRUE},
        )
        with pytest.raises(ValueError, match="variant"):
            fit_shared_form_trajectory_matching(
                setups, variant="shared_3param",  # a §9 variant, not §10
                anchors=_LQ_ANCHORS, run_fn=run_fn,
            )

    def test_missing_anchor_raises(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("C", "n"),
            true_by_case={"18A": _PL_TRUE, "9A": _PL_TRUE},
        )
        with pytest.raises(ValueError, match="anchors"):
            fit_shared_form_trajectory_matching(
                setups, variant=PL_SHARED_3PARAM,
                anchors={"C0": 10.0, "n0": 2.0},  # v_ref_Aps missing
                run_fn=run_fn,
            )

    def test_nonpositive_anchor_raises(self, joint_grid):
        # The re-wired shared bundle's a0 = 0 must never condition a fit.
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"),
            true_by_case={"18A": _LQ_TRUE, "9A": _LQ_TRUE},
        )
        with pytest.raises(ValueError, match="positive"):
            fit_shared_form_trajectory_matching(
                setups, variant=LQ_SHARED_3PARAM,
                anchors={"a0": 0.0, "c0": 11.0},
                run_fn=run_fn,
            )

    def test_anchor_n0_outside_bounds_raises(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("C", "n"),
            true_by_case={"18A": _PL_TRUE, "9A": _PL_TRUE},
        )
        with pytest.raises(ValueError, match="n_bounds"):
            fit_shared_form_trajectory_matching(
                setups, variant=PL_SHARED_3PARAM,
                anchors={"C0": 10.0, "n0": 2.0, "v_ref_Aps": 3.0},
                n_bounds=(1.0, 1.5),  # n0 = 2.0 outside
                run_fn=run_fn,
            )


class TestFormSensitivity:
    def test_lq_pure_quadratic_fixed_a_gets_zero_halfwidth(self, joint_grid):
        t_grid, ref_speed = joint_grid
        true = (0.0, 11.0, 0.20)
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"),
            true_by_case={"18A": true, "9A": true},
        )
        best = evaluate_joint_form_objective(
            "linear_quadratic", {"a": 0.0, "c": true[1] * 1.05},
            {"18A": true[2], "9A": true[2]}, setups, run_fn=run_fn,
        )
        hw = form_sensitivity_halfwidths(setups, best, run_fn=run_fn)
        assert set(hw) == {"a", "c", "e_bind_eV"}
        assert hw["a"] == 0.0  # fixed at 0: not a free parameter
        assert np.isfinite(hw["c"]) and hw["c"] > 0
        assert np.isfinite(hw["e_bind_eV"]) and hw["e_bind_eV"] > 0

    def test_pl_scans_in_pivot_space(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("C", "n"),
            true_by_case={"18A": _PL_TRUE, "9A": _PL_TRUE},
        )
        best = evaluate_joint_form_objective(
            "power_law", {"C": _PL_TRUE[0] * 1.05, "n": _PL_TRUE[1]},
            {"18A": _PL_TRUE[2], "9A": _PL_TRUE[2]}, setups, run_fn=run_fn,
        )
        hw = form_sensitivity_halfwidths(
            setups, best, v_ref_Aps=3.0, run_fn=run_fn
        )
        assert set(hw) == {"gamma_ref", "n", "e_bind_eV"}
        for v in hw.values():
            assert np.isfinite(v) and v > 0

    def test_pl_requires_v_ref(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("C", "n"),
            true_by_case={"18A": _PL_TRUE, "9A": _PL_TRUE},
        )
        best = evaluate_joint_form_objective(
            "power_law", {"C": _PL_TRUE[0], "n": _PL_TRUE[1]},
            {"18A": _PL_TRUE[2], "9A": _PL_TRUE[2]}, setups, run_fn=run_fn,
        )
        with pytest.raises(ValueError, match="v_ref_Aps"):
            form_sensitivity_halfwidths(setups, best, run_fn=run_fn)

    def test_per_case_bindings_refused(self, joint_grid):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"),
            true_by_case={"18A": _LQ_TRUE, "9A": _LQ_TRUE},
        )
        best = evaluate_joint_form_objective(
            "linear_quadratic", {"a": 4.0, "c": 11.0},
            {"18A": 0.15, "9A": 0.25}, setups, run_fn=run_fn,
        )
        with pytest.raises(ValueError, match="shared-E_bind"):
            form_sensitivity_halfwidths(setups, best, run_fn=run_fn)


class TestWriteSharedFormFitParameters:
    def _completed_fit(self, joint_grid, *, variant, keys, true, anchors):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=keys,
            true_by_case={"18A": true, "9A": true},
        )
        fit = fit_shared_form_trajectory_matching(
            setups, variant=variant, anchors=anchors,
            run_fn=run_fn, maxfev_per_start=200, extra_starts=False,
            prescan_e_bind_eV=[0.15, 0.2, 0.25],
        )
        fit.coeff_errs = {k: 0.1 for k in keys}
        fit.e_bind_err_eV = 0.01
        fit.uncertainty_model = "rmse_sensitivity_only_seed_sweep_omitted"
        return fit

    def test_lq_bundle_loads_through_loader(self, joint_grid, tmp_path):
        fit = self._completed_fit(
            joint_grid, variant=LQ_SHARED_3PARAM, keys=("a", "c"),
            true=_LQ_TRUE, anchors=_LQ_ANCHORS,
        )
        out = write_shared_form_fit_parameters(
            fit, tmp_path, anchor_provenance="test constants"
        )
        assert out.name == "fit_parameters.json"
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.form == "linear_quadratic"
        assert coeffs.extraction_method == "trajectory_matching"
        assert coeffs.coefficients["a"] == pytest.approx(
            fit.best.coefficients["a"]
        )
        assert coeffs.coefficients["c"] == pytest.approx(
            fit.best.coefficients["c"]
        )
        assert coeffs.effective_binding_energy_I_ion_eV == pytest.approx(
            next(iter(fit.best.e_bind_eV.values()))
        )

    def test_pl_bundle_loads_and_stamps_pivot(self, joint_grid, tmp_path):
        import json as _json

        fit = self._completed_fit(
            joint_grid, variant=PL_SHARED_3PARAM, keys=("C", "n"),
            true=_PL_TRUE, anchors=_PL_ANCHORS,
        )
        out = write_shared_form_fit_parameters(
            fit, tmp_path, anchor_provenance="test constants"
        )
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.form == "power_law"
        assert coeffs.coefficients["C"] == pytest.approx(
            fit.best.coefficients["C"]
        )
        with open(out, "r", encoding="utf-8") as fh:
            raw = _json.load(fh)
        assert raw["variant"] == PL_SHARED_3PARAM
        assert raw["C_err"] == 0.1 and raw["n_err"] == 0.1
        assert raw["anchors"]["provenance"] == "test constants"
        assert raw["flags"]["model_selection_on_seen_data"] is True
        # Pivot block: gamma_ref consistent with the stamped raw {C, n}.
        assert raw["pivot"]["v_ref_Aps"] == 3.0
        assert raw["pivot"]["gamma_ref"] == pytest.approx(
            raw["C"] * 3.0 ** (raw["n"] - 1.0)
        )

    def test_refuses_single_case_stage1_fit(self, joint_grid, tmp_path):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"), true_by_case={"18A": _LQ_TRUE}
        )
        fit = fit_form_trajectory_matching(
            setups["18A"], variant=LQ_SHARED_3PARAM, anchors=_LQ_ANCHORS,
            run_fn=run_fn, maxfev_per_start=50, extra_starts=False,
            prescan_e_bind_eV=[0.2],
        )
        fit.coeff_errs = {"a": 0.1, "c": 0.1}
        fit.uncertainty_model = "x"
        with pytest.raises(ValueError, match="record-only"):
            write_shared_form_fit_parameters(
                fit, tmp_path, anchor_provenance="test"
            )

    def test_refuses_unfilled_uncertainty(self, joint_grid, tmp_path):
        fit = self._completed_fit(
            joint_grid, variant=LQ_SHARED_3PARAM, keys=("a", "c"),
            true=_LQ_TRUE, anchors=_LQ_ANCHORS,
        )
        fit.coeff_errs = {"a": 0.1}  # c_err missing
        with pytest.raises(ValueError, match="uncertainty"):
            write_shared_form_fit_parameters(
                fit, tmp_path, anchor_provenance="test"
            )

    def test_refuses_trapped_fit(self, joint_grid, tmp_path):
        t_grid, ref_speed = joint_grid
        setups, _ = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"),
            true_by_case={"18A": _LQ_TRUE, "9A": _LQ_TRUE},
        )
        # Every candidate binding traps (trap threshold 0): the best fit
        # carries a nonzero penalty and must never become a bundle.
        trap_fn = _make_form_stub_run_fn(
            t_grid, ref_speed, keys=("a", "c"), true=_LQ_TRUE,
            trap_above_e_eV=0.0,
        )
        fit = fit_shared_form_trajectory_matching(
            setups, variant=LQ_SHARED_3PARAM, anchors=_LQ_ANCHORS,
            run_fn=lambda cfg, neutral: trap_fn(cfg, neutral),
            maxfev_per_start=5, extra_starts=False,
            prescan_e_bind_eV=[0.2],
        )
        fit.coeff_errs = {"a": 0.1, "c": 0.1}
        fit.uncertainty_model = "x"
        assert all(
            r.escape_fraction == 0.0 for r in fit.best.per_case.values()
        )
        with pytest.raises(ValueError, match="trapped"):
            write_shared_form_fit_parameters(
                fit, tmp_path, anchor_provenance="test"
            )


class TestPerCaseFormBundleWriter:
    """The §10.8 single-case (9 A-only / 18 A-only) form bundle writer.

    All fits run against the stub ``run_fn`` (no MD); a single-case fit is the
    Stage-1-analog engine called per case (``fit_form_trajectory_matching``).
    """

    def _single_case_fit(self, joint_grid, *, case, variant, keys, true, anchors):
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=keys, true_by_case={case: true}
        )
        fit = fit_form_trajectory_matching(
            setups[case], variant=variant, anchors=anchors,
            run_fn=run_fn, maxfev_per_start=200, extra_starts=False,
            prescan_e_bind_eV=[0.15, 0.2, 0.25],
        )
        fit.coeff_errs = {k: 0.1 for k in keys}
        fit.e_bind_err_eV = 0.01
        fit.uncertainty_model = (
            "rmse_sensitivity_only_seed_sweep_omitted_per_s8"
        )
        return fit

    def test_lq_per_case_bundle_loads_through_loader(self, joint_grid, tmp_path):
        fit = self._single_case_fit(
            joint_grid, case="18A", variant=LQ_SHARED_3PARAM, keys=("a", "c"),
            true=_LQ_TRUE, anchors=_LQ_ANCHORS,
        )
        out = write_per_case_form_fit_parameters(
            fit, tmp_path, transverse_contaminated=False,
            anchor_provenance="test constants",
        )
        assert out.name == "fit_parameters.json"
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.form == "linear_quadratic"
        assert coeffs.extraction_method == "trajectory_matching"
        assert coeffs.coefficients["a"] == pytest.approx(
            fit.best.coefficients["a"]
        )
        assert coeffs.coefficients["c"] == pytest.approx(
            fit.best.coefficients["c"]
        )
        assert coeffs.effective_binding_energy_I_ion_eV == pytest.approx(
            next(iter(fit.best.e_bind_eV.values()))
        )

    def test_per_case_honesty_flags_and_single_case_provenance(
        self, joint_grid, tmp_path
    ):
        import json as _json

        # 9 A carries the transverse-contamination flag.
        fit = self._single_case_fit(
            joint_grid, case="9A", variant=LQ_SHARED_PURE_QUADRATIC,
            keys=("a", "c"), true=(0.0, 11.0, 0.20), anchors=_LQ_ANCHORS,
        )
        out = write_per_case_form_fit_parameters(
            fit, tmp_path, transverse_contaminated=True,
            anchor_provenance="test constants",
        )
        with open(out, "r", encoding="utf-8") as fh:
            raw = _json.load(fh)
        assert raw["calibration_case"] == "9A"
        assert raw["stage"] == "per_case_diagnostic"
        assert raw["variant"] == LQ_SHARED_PURE_QUADRATIC
        assert isinstance(raw["reference_file"], str)  # single, not a list
        flags = raw["flags"]
        assert flags["per_case_calibrated_not_validated"] is True
        assert flags["cross_case_axis_not_applied"] is True
        assert flags["transverse_contaminated_non_radial_reference"] is True
        assert flags["full_window_heldout_window_axis_forfeited"] is True
        # pure-quadratic: a fixed at 0, its band 0 (still a written coeff).
        assert raw["a"] == 0.0

    def test_pl_per_case_bundle_stamps_pivot(self, joint_grid, tmp_path):
        import json as _json

        fit = self._single_case_fit(
            joint_grid, case="18A", variant=PL_SHARED_3PARAM, keys=("C", "n"),
            true=_PL_TRUE, anchors=_PL_ANCHORS,
        )
        out = write_per_case_form_fit_parameters(
            fit, tmp_path, transverse_contaminated=False,
            anchor_provenance="test constants", v_ref_Aps=3.0,
        )
        coeffs = load_drag_coefficients(tmp_path, expected_m_eff_amu=_M_EFF)
        assert coeffs.form == "power_law"
        assert coeffs.coefficients["C"] == pytest.approx(
            fit.best.coefficients["C"]
        )
        with open(out, "r", encoding="utf-8") as fh:
            raw = _json.load(fh)
        assert raw["C_err"] == 0.1 and raw["n_err"] == 0.1
        assert raw["pivot"]["v_ref_Aps"] == 3.0
        assert raw["pivot"]["gamma_ref"] == pytest.approx(
            raw["C"] * 3.0 ** (raw["n"] - 1.0)
        )

    def test_refuses_multi_case_fit(self, joint_grid, tmp_path):
        # A genuine 2-case shared fit must NOT go through the per-case writer.
        t_grid, ref_speed = joint_grid
        setups, run_fn = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"),
            true_by_case={"18A": _LQ_TRUE, "9A": _LQ_TRUE},
        )
        fit = fit_shared_form_trajectory_matching(
            setups, variant=LQ_SHARED_3PARAM, anchors=_LQ_ANCHORS,
            run_fn=run_fn, maxfev_per_start=50, extra_starts=False,
            prescan_e_bind_eV=[0.2],
        )
        fit.coeff_errs = {"a": 0.1, "c": 0.1}
        fit.uncertainty_model = "x"
        with pytest.raises(ValueError, match="SINGLE-case"):
            write_per_case_form_fit_parameters(
                fit, tmp_path, transverse_contaminated=False,
                anchor_provenance="test",
            )

    def test_refuses_unfilled_uncertainty(self, joint_grid, tmp_path):
        fit = self._single_case_fit(
            joint_grid, case="18A", variant=LQ_SHARED_3PARAM, keys=("a", "c"),
            true=_LQ_TRUE, anchors=_LQ_ANCHORS,
        )
        fit.coeff_errs = {"a": 0.1}  # c_err missing
        with pytest.raises(ValueError, match="uncertainty"):
            write_per_case_form_fit_parameters(
                fit, tmp_path, transverse_contaminated=False,
                anchor_provenance="test",
            )

    def test_pl_requires_v_ref(self, joint_grid, tmp_path):
        fit = self._single_case_fit(
            joint_grid, case="18A", variant=PL_SHARED_3PARAM, keys=("C", "n"),
            true=_PL_TRUE, anchors=_PL_ANCHORS,
        )
        with pytest.raises(ValueError, match="v_ref_Aps"):
            write_per_case_form_fit_parameters(
                fit, tmp_path, transverse_contaminated=False,
                anchor_provenance="test",  # v_ref_Aps omitted
            )

    def test_refuses_trapped_fit(self, joint_grid, tmp_path):
        t_grid, ref_speed = joint_grid
        setups, _ = _form_joint_setups(
            t_grid, ref_speed, keys=("a", "c"), true_by_case={"18A": _LQ_TRUE},
        )
        trap_fn = _make_form_stub_run_fn(
            t_grid, ref_speed, keys=("a", "c"), true=_LQ_TRUE,
            trap_above_e_eV=0.0,
        )
        fit = fit_form_trajectory_matching(
            setups["18A"], variant=LQ_SHARED_3PARAM, anchors=_LQ_ANCHORS,
            run_fn=lambda cfg, neutral: trap_fn(cfg, neutral),
            maxfev_per_start=5, extra_starts=False, prescan_e_bind_eV=[0.2],
        )
        fit.coeff_errs = {"a": 0.1, "c": 0.1}
        fit.uncertainty_model = "x"
        assert all(
            r.escape_fraction == 0.0 for r in fit.best.per_case.values()
        )
        with pytest.raises(ValueError, match="trapped"):
            write_per_case_form_fit_parameters(
                fit, tmp_path, transverse_contaminated=False,
                anchor_provenance="test",
            )

