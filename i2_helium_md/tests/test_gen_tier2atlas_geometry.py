"""Unit coverage for scripts/gen_tier2atlas_geometry.py (atlas Axis A / G1).

Config-build + guard wiring only; no MD. The standing-point diff test needs the
battery reference run dir (``data/runs`` uncommitted) and skips when absent.

The three deliberately-exercised config guards (plan §3.2) are covered here so
they are verified at generator-test time rather than discovered at runtime:
T8 analytic-vs-fixed-size, ``uniform_volume`` vs ``single_initial_position``,
and ``boltzmann`` vs a non-zero birth margin.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from i2_helium_md.simulation.detection_stage import DetectionResult
from scripts.gen_tier2atlas_geometry import (
    CELL_KEYS,
    GEOMETRY_MATRIX,
    N,
    PARENT_WELL_K,
    PROJECT_ROOT,
    R_AXIS,
    R_TARGET_TOL_ANGSTROM,
    REFERENCE_RUN,
    RETAINED_POLICY,
    SEED,
    V_C_APS,
    _BASE_DIFF_KEYS,
    _DETECTION_ONLY_DIFF_KEYS,
    _detection_results_identical,
    atlas_run_dir_name,
    build_cell,
    verify_against_reference,
    verify_geometry,
)


class TestMatrix:
    def test_eleven_cells_three_by_three_plus_two_large_R(self):
        assert len(GEOMETRY_MATRIX) == 11
        labels = [c.label for c in GEOMETRY_MATRIX]
        assert labels == list(CELL_KEYS)
        # full factorial on r1..r3, L2/L3 only at r4 (the funded large-R pair)
        for r in ("r1", "r2", "r3"):
            assert {c.law_key for c in GEOMETRY_MATRIX if c.r_key == r} == {
                "l1", "l2", "l3"}
        assert {c.law_key for c in GEOMETRY_MATRIX if c.r_key == "r4"} == {
            "l2", "l3"}

    def test_one_shared_seed_and_cell_N(self):
        assert {c.label for c in GEOMETRY_MATRIX} == set(CELL_KEYS)
        cfgs = [build_cell(c) for c in GEOMETRY_MATRIX]
        assert {c.seed for c in cfgs} == {SEED}      # common random numbers
        assert {c.num_molecules for c in cfgs} == {N}

    @pytest.mark.parametrize("r_key,n_he,target_R", [
        (k, v[0], v[1]) for k, v in R_AXIS.items()
    ])
    def test_bulk_conversion_hits_each_target_radius(self, r_key, n_he, target_R):
        # R = 2.2173*N^(1/3) (D0 §15.3) -- the grid's He counts are pinned to
        # the plan's radii, so a convention change here must fail loudly.
        cell = next(c for c in GEOMETRY_MATRIX if c.r_key == r_key)
        assert cell.droplet_N_he == n_he
        assert abs(cell.realized_R_angstrom - target_R) <= R_TARGET_TOL_ANGSTROM

    def test_namespace_locks(self):
        for cell in GEOMETRY_MATRIX:
            name = atlas_run_dir_name(cell.label)
            assert "_tier2atlas_" in name
            assert "_tier2_" not in name and "tier2probe" not in name
            assert f"_N{N}_" in name
        assert len({atlas_run_dir_name(c.label) for c in GEOMETRY_MATRIX}) == 11


class TestCellBuild:
    @pytest.mark.parametrize("cell", GEOMETRY_MATRIX, ids=lambda c: c.label)
    def test_builds_validates_and_pins_the_geometry(self, cell):
        cfg = build_cell(cell)
        cfg.validate()
        verify_geometry(cfg, cell)
        # fixed size: sampling off, legacy prior (the T8 pairing)
        assert cfg.use_single_droplet_size is True
        assert cfg.single_droplet_size == cell.droplet_N_he
        assert cfg.droplet_size_prior == "legacy"
        # standing production drag point rides unchanged
        assert cfg.drag_form == "capped_cubic"
        assert cfg.drag_coefficients.coefficients["v_c"] == V_C_APS
        assert cfg.internal_energy_cooling_tau_ps == 3.2

    def test_law_wiring_per_column(self):
        by_label = {c.label: build_cell(c) for c in GEOMETRY_MATRIX}
        # L1 center-pin: r0 sampled then zeroed; margin must be 0 under boltzmann
        l1 = by_label["r1l1"]
        assert l1.single_initial_position is True
        assert l1.birth_position_law == "boltzmann"
        assert l1.initial_position_margin_angstrom == 0.0
        assert l1.binding_energy_molecule_K == 573.3   # inert at r0 == 0
        # L2 parent law: the parent's own DFT-fit well, per-run override
        l2 = by_label["r1l2"]
        assert l2.single_initial_position is False
        assert l2.birth_position_law == "boltzmann"
        assert l2.initial_position_margin_angstrom == 0.0
        assert l2.binding_energy_molecule_K == PARENT_WELL_K == 313.2
        # L3 standing law
        l3 = by_label["r1l3"]
        assert l3.single_initial_position is False
        assert l3.birth_position_law == "uniform_volume"
        assert l3.initial_position_margin_angstrom == 3.0
        assert l3.binding_energy_molecule_K == 573.3

    def test_e_bind_stays_paired_with_the_drag_bundle(self):
        # Atlas G0-3: the geometry axis does NOT touch the well, so the §6.5.1
        # pairing must hold on every cell (no unvalidated-binding hatch here).
        for cell in GEOMETRY_MATRIX:
            cfg = build_cell(cell)
            assert cfg.binding_energy_I_ion_eV == pytest.approx(
                cfg.drag_coefficients.effective_binding_energy_I_ion_eV)
            assert cfg.allow_unvalidated_binding_pairing is False

    def test_only_geometry_differs_across_cells(self):
        # Everything outside the geometry keys must be identical cell-to-cell.
        geometry_keys = {
            "single_droplet_size", "birth_position_law",
            "initial_position_margin_angstrom", "single_initial_position",
            "binding_energy_molecule_K",
        }
        base = dataclasses.asdict(build_cell(GEOMETRY_MATRIX[0]))
        for cell in GEOMETRY_MATRIX[1:]:
            other = dataclasses.asdict(build_cell(cell))
            diff = {k for k in base if base[k] != other[k]}
            assert diff <= geometry_keys, f"{cell.label} moved {diff - geometry_keys}"


class TestExercisedGuards:
    """The three guards the grid trips on purpose must be live (plan §3.2)."""

    def test_analytic_prior_refused_under_fixed_size(self):
        cfg = build_cell(GEOMETRY_MATRIX[0])
        bad = dataclasses.replace(cfg, droplet_size_prior="kornilov_lognormal")
        with pytest.raises(ValueError, match="use_single_droplet_size"):
            bad.validate()

    def test_uniform_volume_refused_with_center_pin(self):
        l3 = build_cell(next(c for c in GEOMETRY_MATRIX if c.law_key == "l3"))
        bad = dataclasses.replace(l3, single_initial_position=True)
        with pytest.raises(ValueError, match="single_initial_position"):
            bad.validate()

    def test_boltzmann_refused_with_nonzero_margin(self):
        l2 = build_cell(next(c for c in GEOMETRY_MATRIX if c.law_key == "l2"))
        bad = dataclasses.replace(l2, initial_position_margin_angstrom=3.0)
        with pytest.raises(ValueError, match="initial_position_margin_angstrom"):
            bad.validate()

    def test_off_target_droplet_size_is_caught(self):
        cell = GEOMETRY_MATRIX[0]
        cfg = build_cell(cell)
        wrong = cell._replace(droplet_N_he=cell.droplet_N_he * 2)
        with pytest.raises(AssertionError, match="off the grid target"):
            verify_geometry(cfg, wrong)


def _reference_present() -> bool:
    return (PROJECT_ROOT / "data" / "runs" / REFERENCE_RUN / "cfg.json").exists()


@pytest.mark.skipif(not _reference_present(),
                    reason="standing battery run dir absent (data/runs uncommitted)")
class TestStandingPointDiff:
    @pytest.mark.parametrize("cell", GEOMETRY_MATRIX, ids=lambda c: c.label)
    def test_diff_vs_standing_point_is_geometry_only(self, cell):
        diff = verify_against_reference(build_cell(cell), cell)
        assert "num_molecules" in diff and "seed" in diff
        assert "use_single_droplet_size" in diff
        # the drag/mass/ladder/cooling knobs must NOT appear
        for forbidden in ("drag_form", "internal_energy_cooling_tau_ps",
                          "internal_energy_partition_fraction",
                          "binding_energy_I_ion_eV", "v_limit_m_per_s",
                          "cooling_spatial_gate", "tabulated_ladder_rungs_eV"):
            assert forbidden not in diff

    def test_a_moved_knob_is_reported(self):
        cell = GEOMETRY_MATRIX[0]
        tampered = dataclasses.replace(
            build_cell(cell), internal_energy_cooling_tau_ps=4.8)
        with pytest.raises(AssertionError, match="cfg diff"):
            verify_against_reference(tampered, cell)


class TestRetainedPolicyWiring:
    """Plan §3.5b: the grid runs the marginal-counting arm on every cell, so
    the grid stays one homogeneous configuration and the R <= 34 A cells double
    as the arm's bit-for-bit oracle."""

    def test_every_cell_runs_exclude_all_coupled(self):
        assert RETAINED_POLICY == "exclude_all_coupled"
        for cell in GEOMETRY_MATRIX:
            cfg = build_cell(cell)
            assert cfg.detection_droplet_retained_policy == "exclude_all_coupled"

    def test_policy_is_a_declared_diff_key(self):
        # It differs from the standing battery member (which runs "exclude"),
        # so it must be pre-registered or verify_against_reference would fail.
        assert "detection_droplet_retained_policy" in _BASE_DIFF_KEYS

    def test_detection_only_may_change_the_policy_and_nothing_else(self):
        assert _DETECTION_ONLY_DIFF_KEYS == {"detection_droplet_retained_policy"}


class TestDetectionResultComparison:
    """The oracle's comparison helper: it must catch a single changed element,
    not just a changed shape (it is the only thing standing between a silent
    overwrite and the item-10 bit-for-bit claim)."""

    def _result(self, **overrides):
        two_n = 4
        base = dict(
            num_molecules=2, t_handover_ps=20.0, detection_time_ps=8.53e6,
            n_detected=np.zeros(two_n), E_int_detected_eV=np.zeros(two_n),
            mass_detected_kg=np.ones(two_n), vx_detected=np.zeros(two_n),
            vy_detected=np.zeros(two_n), vz_detected=np.zeros(two_n),
            E_kin_detected_eV=np.zeros(two_n),
            E_pot_detected_eV=np.zeros(two_n),
            E_dissip_detected_eV=np.zeros(two_n),
            E_mass_transfer_detected_eV=np.zeros(two_n),
            state_reason=np.array(["frozen"] * two_n),
            event_offsets=np.zeros(two_n + 1, dtype=int),
            event_time_ps=np.zeros(0), event_pre_shed_n=np.zeros(0),
            event_dE_int_eV=np.zeros(0), event_dE_bind_fold_eV=np.zeros(0),
            event_dE_mass_transfer_eV=np.zeros(0),
        )
        base.update(overrides)
        return DetectionResult(**base)

    def test_identical_results_compare_equal(self):
        identical, diff = _detection_results_identical(
            self._result(), self._result()
        )
        assert identical and diff == ""

    def test_one_changed_float_is_caught(self):
        other = self._result()
        other.E_kin_detected_eV[2] = 1e-9
        identical, diff = _detection_results_identical(self._result(), other)
        assert not identical
        assert "E_kin_detected_eV" in diff

    def test_one_changed_reason_is_caught(self):
        other = self._result(
            state_reason=np.array(
                ["frozen", "droplet_retained_marginal", "frozen", "frozen"]
            )
        )
        identical, diff = _detection_results_identical(self._result(), other)
        assert not identical
        assert "state_reason" in diff

    def test_a_changed_scalar_is_caught(self):
        identical, diff = _detection_results_identical(
            self._result(), self._result(t_handover_ps=21.0)
        )
        assert not identical
        assert "t_handover_ps" in diff
