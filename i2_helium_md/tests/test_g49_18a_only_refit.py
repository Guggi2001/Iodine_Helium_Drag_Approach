"""THESIS_FINAL_GOALS G4.9: 18 A-only refit surface (9 A held out as non-radial).

No MD: stub landscapes as in ``test_extraction_trajectory_matching``.

1. ``fit_shared_trajectory_matching`` refuses one case by default and accepts it
   with ``allow_single_case=True``; the single-case objective is that case's
   objective and the fit recovers the stub truth;
2. the writer's flags record a single-case 18 A fit honestly (no cross-case axis
   consumed, no non-radial reference in the objective);
3. the new catalog ids resolve to the 18 A folder, on either geometry;
4. the driver's pure helpers: the G4.9 bands and the step-1 residual shares.
"""

import json

import numpy as np
import pytest

from i2_helium_md import SimConfig
from i2_helium_md.extraction.trajectory_matching import (
    SHARED_CUBIC,
    evaluate_objective,
    fit_shared_trajectory_matching,
    write_shared_fit_parameters,
)
from i2_helium_md.presets import REFERENCE_DRAG_ROOT, resolve_drag_bundle_dir

from .test_extraction_trajectory_matching import (  # type: ignore
    _TRUE,
    _joint_setups,
)

_TRUE_CUBIC = (0.0, 2.0, 0.20)


@pytest.fixture
def grid():
    t_grid = np.linspace(2.0, 12.0, 50)
    return t_grid, np.full(50, 3.5)


class TestSingleCaseFit:
    def test_single_case_refused_by_default(self, grid):
        t_grid, ref = grid
        setups, run_fn = _joint_setups(t_grid, ref, {"18A": _TRUE})
        with pytest.raises(ValueError, match="allow_single_case"):
            fit_shared_trajectory_matching(
                setups, variant=SHARED_CUBIC, anchors=(14.0, 2.0), run_fn=run_fn,
            )

    def test_empty_refused_even_with_flag(self):
        with pytest.raises(ValueError, match="cross-case"):
            fit_shared_trajectory_matching(
                {}, variant=SHARED_CUBIC, anchors=(14.0, 2.0), allow_single_case=True,
            )

    def test_single_case_recovers_truth_and_objective_is_the_case_objective(self, grid):
        t_grid, ref = grid
        setups, run_fn = _joint_setups(t_grid, ref, {"18A": _TRUE_CUBIC})
        fit = fit_shared_trajectory_matching(
            setups, variant=SHARED_CUBIC, anchors=(14.0, 2.0), run_fn=run_fn,
            maxfev_per_start=400, extra_starts=False,
            prescan_e_bind_eV=[0.15, 0.2, 0.25], allow_single_case=True,
        )
        assert fit.cases == ["18A"]
        assert fit.best.a == 0.0
        assert fit.best.b == pytest.approx(_TRUE_CUBIC[1], rel=0.01)
        assert fit.best.e_bind_eV["18A"] == pytest.approx(_TRUE_CUBIC[2], abs=0.005)
        # Mean over one case == that case's objective, bit for bit.
        direct = evaluate_objective(
            (fit.best.a, fit.best.b, fit.best.e_bind_eV["18A"]), setups["18A"],
            run_fn=run_fn,
        )
        assert fit.best.objective == direct.objective


class TestWriterFlags:
    def test_single_case_18a_flags(self, grid, tmp_path):
        t_grid, ref = grid
        setups, stub = _joint_setups(t_grid, ref, {"18A": _TRUE_CUBIC})
        setups["18A"].cfg_base = SimConfig(
            drag_spatial_gate="erf_independent", drag_gate_steepness=3.5)

        def run_fn(cfg, neutral):
            # Move the stub's parked ions past 2 x potential width (escape).
            ion = stub(cfg, neutral)
            ion.positions_x = ion.positions_x + 3.0 * cfg.potential_steepness
            return ion

        fit = fit_shared_trajectory_matching(
            setups, variant=SHARED_CUBIC, anchors=(14.0, 2.0), run_fn=run_fn,
            maxfev_per_start=5, extra_starts=False, prescan_e_bind_eV=[0.2],
            allow_single_case=True,
        )
        fit.a_err, fit.b_err, fit.e_bind_err_eV, fit.uncertainty_model = 0, 0.1, 0.01, "x"
        out = write_shared_fit_parameters(
            fit, tmp_path, stage="single_case_18A", bundle_id="srho3p5A18A_cubic",
        )
        raw = json.loads(out.read_text(encoding="utf-8"))
        assert raw["calibration_cases"] == ["18A"]
        assert raw["stage"] == "single_case_18A"
        assert raw["flags"]["transverse_contaminated_non_radial_reference"] is False
        assert raw["flags"]["cross_case_axis_consumed_by_joint_fit"] is False


class TestCatalog:
    @pytest.mark.parametrize("family", ["linear", "cubic", "linear_cubic"])
    @pytest.mark.parametrize("case", ["9A", "18A"])
    def test_ids_resolve_to_the_18a_folder_on_either_geometry(self, family, case):
        path = resolve_drag_bundle_dir(f"srho3p5A18A_{family}", case)
        assert path == REFERENCE_DRAG_ROOT / "18A" / "trajectory_matching_srho3p5A" / family

    def test_driver_bundle_id_and_export_folder(self):
        from scripts.extraction.method_b_18A_only_refit import bundle_id
        from scripts.tier0_common import reference_export_dir_name

        assert bundle_id("cubic") == "srho3p5A18A_cubic"
        # Distinct from the G4 export folder (tier0_srho3p5A): no overwrite.
        assert reference_export_dir_name("tier0", "srho3p5A18A_cubic") == "tier0_srho3p5A18A"


class TestDriverHelpers:
    def test_bands(self):
        from scripts.extraction.method_b_18A_only_refit import band_checks_18a

        assert band_checks_18a({"rmse_Aps": 0.19, "escape_fraction": 1.0})["all_pass"]
        assert not band_checks_18a({"rmse_Aps": 0.191, "escape_fraction": 1.0})["all_pass"]
        assert not band_checks_18a({"rmse_Aps": 0.1, "escape_fraction": 0.99})["all_pass"]

    def test_residual_shares_locate_the_extra_residual(self):
        from scripts.extraction.method_b_18A_only_refit import residual_shares

        t = np.linspace(2.0, 14.0, 13)               # 2..14 ps, 1 ps steps
        frac = np.where(t >= 8.0, 0.8, 0.1)          # transverse from 8 ps
        r14 = np.full(t.size, 0.1)
        r35 = np.where(t >= 8.0, 0.3, 0.1)           # extra residual only there
        s = residual_shares(t, r14, r35, frac, fraction_min=0.5, t_split_ps=6.0)
        assert s["transverse"]["share_of_extra_sq_residual"] == pytest.approx(1.0)
        assert s["transverse"]["share_of_samples"] == pytest.approx(7 / 13)
        assert s["late"]["share_of_extra_sq_residual"] == pytest.approx(1.0)
        assert s["late"]["share_of_samples"] == pytest.approx(9 / 13)

    def test_residual_shares_refuse_no_extra_residual_and_bad_shapes(self):
        from scripts.extraction.method_b_18A_only_refit import residual_shares

        t = np.linspace(0, 1, 5)
        with pytest.raises(ValueError, match="no extra residual"):
            residual_shares(t, t, t, t, fraction_min=0.5, t_split_ps=0.5)
        with pytest.raises(ValueError, match="same shape"):
            residual_shares(t, t[:3], t, t, fraction_min=0.5, t_split_ps=0.5)
