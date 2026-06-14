"""Deterministic mechanical-correctness smoke run for the Tier-0 drag branch.

Slice 4 §8 "thorough debug" artifact: a tiny fixed-seed drag run on a *synthetic*
neutral checkpoint (decoupled from the neutral stage). Proves the machine runs
*correctly* -- finite trajectories, tight energy closure, monotone positive
dissipation, a v5 checkpoint round-trip, and scope-guard rejection. It asserts
**nothing** about matching TDDFT; that is the separate Tier-0 analysis task.
"""

from dataclasses import replace

import numpy as np
import pytest

from i2_helium_md.physics.constants import U
from i2_helium_md.physics.drag import (
    DragCoefficients,
    LINEAR_QUADRATIC,
    POWER_LAW,
    THRESHOLD,
)
from i2_helium_md.presets import single_pulse_N2000_drag
from i2_helium_md.simulation.checkpoint import (
    NeutralCheckpoint,
    load_ion_checkpoint,
    save_ion_checkpoint,
)
from i2_helium_md.simulation.ion import run_ion_propagation


def _synthetic_neutral(num_molecules=2, mass_amu=202.953908, droplet_radius=30.0):
    """Tiny hand-built neutral end-state: ion pairs inside a droplet, set to repel.

    Synthetic (not a neutral run), per the Slice 4 smoke-harness decision. Mass
    is set to ``m_eff`` so the fixed-mass drag law is self-consistent (§6.5). The
    2N layout puts atom 1 of molecule ``i`` at index ``i`` and atom 2 at index
    ``N+i``; each molecule's two atoms sit ~2.6 A apart so the Coulomb repulsion
    drives an explosion the gated drag then damps. All atoms are well inside the
    droplet (depth < 0 -> gate ~ 1, drag active).
    """
    N = num_molecules
    two_N = 2 * N
    x = np.zeros(two_N)
    y = np.zeros(two_N)
    z = np.zeros(two_N)
    sep = 1.3
    for i in range(N):
        axis = i % 3
        offset = 4.0 * i  # keep molecules apart but inside the droplet
        if axis == 0:
            x[i], x[N + i] = sep, -sep
        elif axis == 1:
            y[i], y[N + i] = sep, -sep
        else:
            z[i], z[N + i] = sep, -sep
        z[i] += offset
        z[N + i] += offset

    def col(a):
        return a.reshape(two_N, 1)

    zeros_2Nk = np.zeros((two_N, 1))
    return NeutralCheckpoint(
        num_molecules=N,
        time_ps=np.zeros(1),
        positions_x=col(x), positions_y=col(y), positions_z=col(z),
        velocities_x=zeros_2Nk.copy(),
        velocities_y=zeros_2Nk.copy(),
        velocities_z=zeros_2Nk.copy(),
        mass_kg=np.full(two_N, mass_amu * U),
        droplet_radii=np.full(two_N, droplet_radius),
        r0=np.zeros(N),
        E_kin_eV=zeros_2Nk.copy(),
        E_pot_eV=zeros_2Nk.copy(),
        E_initial_eV=np.zeros(N),
        E_dissip_eV=zeros_2Nk.copy(),
        L_droplet_eV_ps=zeros_2Nk.copy(),
    )


@pytest.fixture
def drag_cfg():
    """Tiny deterministic Tier-0 drag config (20 internal steps, fixed seed)."""
    return single_pulse_N2000_drag(
        num_molecules=2, ion_simulation_time=0.2, dt_ion=0.01, seed=0,
    )


@pytest.fixture
def neutral():
    return _synthetic_neutral(num_molecules=2)


def test_dispatches_to_drag_branch(drag_cfg, neutral):
    """drag_coefficients present -> BAOAB path -> no collisions ever recorded."""
    ck = run_ion_propagation(drag_cfg, neutral)
    assert np.all(ck.number_of_collisions == 0)


def test_finite_trajectories(drag_cfg, neutral):
    ck = run_ion_propagation(drag_cfg, neutral)
    for name in (
        "positions_x", "positions_y", "positions_z",
        "velocities_x", "velocities_y", "velocities_z",
        "E_kin_eV", "E_pot_eV", "E_dissip_eV",
    ):
        assert np.all(np.isfinite(getattr(ck, name))), name


def test_energy_closes_tight(drag_cfg, neutral):
    """E_kin + E_pot + E_dissip conserved to baseline conservative-Verlet level.

    Dissipation is booked analytically (exact), so the only drift is the
    conservative Verlet drift on E_kin + E_pot.
    """
    ck = run_ion_propagation(drag_cfg, neutral)
    E0 = (ck.E_kin_eV[:, 0] + ck.E_pot_eV[:, 0] + ck.E_dissip_eV[:, 0]).sum()
    E1 = (ck.E_kin_eV[:, -1] + ck.E_pot_eV[:, -1] + ck.E_dissip_eV[:, -1]).sum()
    rel = abs(E1 - E0) / abs(E0)
    assert rel < 1e-2, f"energy drift {rel * 100:.4f}% exceeds Verlet tolerance"


def test_dissipation_monotone_and_positive(drag_cfg, neutral):
    ck = run_ion_propagation(drag_cfg, neutral)
    assert np.all(np.diff(ck.E_dissip_eV, axis=1) >= 0)
    assert ck.E_dissip_eV[:, -1].sum() > 0.0


def test_checkpoint_v5_round_trip(tmp_path, drag_cfg, neutral):
    ck = run_ion_propagation(drag_cfg, neutral)
    path = save_ion_checkpoint(ck, tmp_path / "ion_drag.npz")
    loaded = load_ion_checkpoint(path)
    assert loaded.schema_version == 5
    np.testing.assert_array_equal(loaded.E_dissip_eV, ck.E_dissip_eV)
    # Tier-0 fills survive the round-trip:
    assert np.all(loaded.E_mass_attach_defect_eV == 0.0)
    assert np.all(loaded.number_of_collisions == 0)
    assert np.all(np.isnan(loaded.temperature_diagnostic))


def test_scope_guard_rejects_active_noise(drag_cfg, neutral):
    bad = replace(drag_cfg, noise_form="multiplicative_local_fdt")
    with pytest.raises(NotImplementedError, match="noise_form"):
        run_ion_propagation(bad, neutral)


def test_scope_guard_rejects_mass_scenario(drag_cfg, neutral):
    bad = replace(
        drag_cfg,
        mass_scenario="scenario_A_accretion",
        allow_inconsistent_mass_pairing=True,  # bypass the §6.5 config guard
    )
    with pytest.raises(NotImplementedError, match="mass_scenario"):
        run_ion_propagation(bad, neutral)


def test_scope_guard_rejects_wrong_realized_mass(drag_cfg):
    """Slice 4 fix Change B: the realized-mass trip-wire fires on a ~127 amu
    integration mass (76 amu below m_eff) even though the config is consistent.

    This is the direct coverage that makes the smoke harness's m_eff sidestep no
    longer a silent blind spot.
    """
    from i2_helium_md.simulation.ion_propagation_step import _check_drag_scope

    mass_iodine_kg = np.full(4, 126.9 * U)  # bare I+, not m_eff
    with pytest.raises(NotImplementedError, match="m_eff"):
        _check_drag_scope(drag_cfg, mass_iodine_kg)


def test_scope_guard_accepts_m_eff_mass(drag_cfg):
    """The same guard passes when the realized mass is m_eff (happy path)."""
    from i2_helium_md.simulation.ion_propagation_step import _check_drag_scope

    mass_meff_kg = np.full(4, drag_cfg.m_eff_amu * U)
    _check_drag_scope(drag_cfg, mass_meff_kg)  # must not raise


# ---------------------------------------------------------------------------
# METHOD_B §10 form phase: the realised families run through the same driver
# ---------------------------------------------------------------------------
def _form_cfg(drag_cfg, form, coefficients):
    """The drag preset re-pointed at a form-phase coefficient bundle.

    Stamps the config's own binding (exact §6.5.1 pairing) and m_eff (exact
    §6.5 consistency), mirroring how the extraction objective builds configs.
    """
    coeffs = DragCoefficients(
        form=form,
        coefficients=coefficients,
        extraction_mass_model="constant",
        extraction_mass_amu=drag_cfg.m_eff_amu,
        extraction_method="trajectory_matching",
        effective_binding_energy_I_ion_eV=drag_cfg.binding_energy_I_ion_eV,
    )
    return replace(drag_cfg, drag_form=form, drag_coefficients=coeffs)


_FORM_PHASE_PARAMS = [
    (LINEAR_QUADRATIC, {"a": 4.0, "c": 11.0}),
    (LINEAR_QUADRATIC, {"a": 0.0, "c": 11.0}),   # pure-quadratic variant
    (POWER_LAW, {"C": 10.36, "n": 2.056}),
]


@pytest.mark.parametrize("form,coefficients", _FORM_PHASE_PARAMS)
def test_scope_guard_accepts_realised_forms(drag_cfg, form, coefficients):
    from i2_helium_md.simulation.ion_propagation_step import _check_drag_scope

    cfg = _form_cfg(drag_cfg, form, coefficients)
    cfg.validate()  # config-level guards (form agreement, §10.3 arms) pass
    mass_meff_kg = np.full(4, cfg.m_eff_amu * U)
    _check_drag_scope(cfg, mass_meff_kg)  # must not raise


def test_scope_guard_rejects_reserved_threshold(drag_cfg):
    from i2_helium_md.simulation.ion_propagation_step import _check_drag_scope

    cfg = _form_cfg(drag_cfg, THRESHOLD, {"F_sat": 1.0, "v0": 1.0})
    mass_meff_kg = np.full(4, cfg.m_eff_amu * U)
    with pytest.raises(NotImplementedError, match="drag_form"):
        _check_drag_scope(cfg, mass_meff_kg)


@pytest.mark.parametrize("form,coefficients", _FORM_PHASE_PARAMS)
def test_form_phase_families_run_end_to_end(drag_cfg, neutral, form, coefficients):
    """Tiny deterministic run per family: finite, dissipative, collision-free."""
    cfg = _form_cfg(drag_cfg, form, coefficients)
    ck = run_ion_propagation(cfg, neutral)
    for name in (
        "positions_x", "positions_y", "positions_z",
        "velocities_x", "velocities_y", "velocities_z",
        "E_kin_eV", "E_pot_eV", "E_dissip_eV",
    ):
        assert np.all(np.isfinite(getattr(ck, name))), name
    assert np.all(np.diff(ck.E_dissip_eV, axis=1) >= 0)
    assert ck.E_dissip_eV[:, -1].sum() > 0.0
    assert np.all(ck.number_of_collisions == 0)
