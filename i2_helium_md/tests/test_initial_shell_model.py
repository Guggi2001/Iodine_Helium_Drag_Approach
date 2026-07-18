"""Slice T5 (Tier-2 plan §I.11; revives H.3b): the initial-shell dressing arm.

Config enum ``initial_shell_model ∈ {"full" (byte-inert default),
"density_tied"}``. Under ``density_tied`` the biphasic column-0 seed dresses
each ion by the local He availability at its birth position,

    n_0i = round(n* · rho_He/rho_bulk(d_birth,i)),

through the SAME erf-complement surface the drag/pickup/cooling gates share
(``helium_density.rho_he_ratio`` at ``drag_gate_steepness(cfg)``) — zero new
free parameters. Per-ion consequences at the seed: ``n_shell(0)``,
``mass_kg(0) = complex_mass_amu(n_0i)``, and the t0 ``e_bind_pair(n_0i)``
E_pot fold. The E_int(0) onset is NOT coupled to the dressing here — that is
Slice T6's p-law.

H.3b wiring oracles covered below:
* the ``full`` arm keeps the delivered 21-for-all values exactly;
* center-pinned births + ``density_tied`` give n_0 = 21 exactly (the arm is
  inert without the position axis);
* n_0 is monotone non-increasing in birth radius.
"""

from __future__ import annotations

import json
import math
from dataclasses import fields, replace

import numpy as np
import pytest

from i2_helium_md.config import SimConfig, check_initial_shell_config
from i2_helium_md.physics.constants import U
from i2_helium_md.physics.dissociation_ladder import resolve_ladder
from i2_helium_md.physics.shell_schedule import ANCHOR_N_START, complex_mass_amu
from i2_helium_md.physics.solvation_cooling import e_bind_pair_eV
from i2_helium_md.presets import single_pulse_N2000, single_pulse_N2000_drag
from i2_helium_md.simulation.checkpoint import NeutralCheckpoint
from i2_helium_md.simulation.ion import drag_gate_steepness
from i2_helium_md.simulation.ion_initial_state import build_initial_ion_state
from i2_helium_md.simulation.neutral import run_neutral_propagation
from i2_helium_md.simulation.run_directory import RunDirectory


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def small_neutral_run():
    """Short neutral propagation: 5 molecules, 10 steps (= 0.1 ps)."""
    cfg = single_pulse_N2000(num_molecules=5, seed=42)
    cfg = replace(cfg, t_max_neutral=0.1, dt_neutral=0.01)
    neutral = run_neutral_propagation(cfg, verbose=False)
    return cfg, neutral


def _biphasic_cfg(neutral_cfg, **extra):
    """A valid Tier-2 biphasic config on the Tier-0 drag bundle.

    ``pickup_rate_coefficient=0.9`` (the bridge lambda_0) keeps
    ``check_biphasic_config`` warning-free so ``validate()`` runs clean.
    """
    drag_cfg = single_pulse_N2000_drag(
        num_molecules=neutral_cfg.num_molecules,
        t_max_neutral=neutral_cfg.t_max_neutral,
        dt_neutral=neutral_cfg.dt_neutral,
    )
    return replace(
        drag_cfg,
        mass_scenario="biphasic",
        internal_energy_partition_fraction=0.3,
        internal_energy_retained_fraction=0.2,
        pickup_rate_coefficient=0.9,
        allow_inconsistent_mass_pairing=True,  # §6.5 structural trip (§6.6 defence)
        **extra,
    )


def _deep_copy_neutral(neutral):
    """NeutralCheckpoint with all numpy arrays deep-copied.

    ``dataclasses.replace`` shares array references; tests that mutate
    positions must copy first (the test_ion_initial_state idiom).
    """
    kwargs = {}
    for f in fields(NeutralCheckpoint):
        v = getattr(neutral, f.name)
        kwargs[f.name] = v.copy() if isinstance(v, np.ndarray) else v
    return NeutralCheckpoint(**kwargs)


def _neutral_with_radii(neutral, radii_angstrom):
    """Deep-copied checkpoint with last-column atom radii set along +x."""
    ckpt = _deep_copy_neutral(neutral)
    r = np.asarray(radii_angstrom, dtype=float)
    ckpt.positions_x[:, -1] = r
    ckpt.positions_y[:, -1] = 0.0
    ckpt.positions_z[:, -1] = 0.0
    return ckpt


def _expected_n0(cfg, ckpt):
    """Independent hand oracle: n_0 = rint(n* · 0.5·(1 − erf(depth/steepness)))."""
    r = np.sqrt(
        ckpt.positions_x[:, -1] ** 2
        + ckpt.positions_y[:, -1] ** 2
        + ckpt.positions_z[:, -1] ** 2
    )
    depth = r - ckpt.droplet_radii
    steepness = drag_gate_steepness(cfg)
    rho = np.array([0.5 * (1.0 - math.erf(d / steepness)) for d in depth])
    return np.rint(ANCHOR_N_START * rho)


def _fold_eV(cfg, n):
    """The e_bind_pair fold at the cfg's resolved ladder (scalar or array n)."""
    return e_bind_pair_eV(
        n,
        picture=cfg.ladder_electronic_picture,
        kappa=cfg.ladder_steepness,
        ladder=resolve_ladder(
            cfg.dissociation_ladder, cfg.tabulated_ladder_rungs_eV
        ),
    )


# ===========================================================================
# Config surface
# ===========================================================================
class TestConfigSurface:
    def test_default_is_full(self):
        assert SimConfig().initial_shell_model == "full"

    def test_default_config_still_validates(self):
        SimConfig().validate()

    def test_unknown_selector_refused(self, small_neutral_run):
        cfg, _ = small_neutral_run
        bad = _biphasic_cfg(cfg, initial_shell_model="densty_tied")
        with pytest.raises(ValueError, match="initial_shell_model"):
            bad.validate()

    def test_density_tied_requires_biphasic(self):
        # mass_scenario defaults to "fixed": the dressing law would be
        # silently inert there (the seed only reads it under biphasic) —
        # refused loudly instead (the T7/F3 no-silent-inert convention).
        cfg = SimConfig(initial_shell_model="density_tied")
        with pytest.raises(ValueError, match="biphasic"):
            check_initial_shell_config(cfg)

    def test_density_tied_under_biphasic_validates(self, small_neutral_run):
        cfg, _ = small_neutral_run
        _biphasic_cfg(cfg, initial_shell_model="density_tied").validate()

    def test_pre_slice_cfg_json_loads_with_default(self, tmp_path):
        # Slice-DS/T7 precedent: a cfg.json written before the field existed
        # must load with the byte-inert default.
        run = RunDirectory(tmp_path / "9A_drag_tier2probe_backcompat_shell")
        run.save_cfg(SimConfig())
        payload = json.loads(run.cfg_path.read_text(encoding="utf-8"))
        assert payload.pop("initial_shell_model") == "full"
        run.cfg_path.write_text(json.dumps(payload), encoding="utf-8")
        assert run.load_cfg().initial_shell_model == "full"

    def test_cfg_json_roundtrip(self, tmp_path, small_neutral_run):
        cfg, _ = small_neutral_run
        run = RunDirectory(tmp_path / "9A_drag_tier2probe_shell_roundtrip")
        run.save_cfg(_biphasic_cfg(cfg, initial_shell_model="density_tied"))
        assert run.load_cfg().initial_shell_model == "density_tied"


# ===========================================================================
# The dressed biphasic seed
# ===========================================================================
class TestDressedSeed:
    def test_full_arm_values_unchanged(self, small_neutral_run):
        # The delivered 21-for-all convention: uniform n*=21 mass, scalar
        # e_bind_pair(21) fold on top of the undressed (fixed-scenario,
        # fold-free) E_pot — same cfg geometry, E_pot is mass-independent.
        cfg, neutral = small_neutral_run
        bip = _biphasic_cfg(cfg)  # initial_shell_model rides the default
        ion_full = build_initial_ion_state(bip, neutral, num_steps_ion=5)
        np.testing.assert_allclose(
            ion_full.mass_kg, complex_mass_amu(ANCHOR_N_START) * U
        )
        np.testing.assert_array_equal(
            ion_full.n_shell[:, 0], float(ANCHOR_N_START)
        )
        drag_cfg = replace(
            bip,
            mass_scenario="fixed",
            internal_energy_partition_fraction=None,
            internal_energy_retained_fraction=None,
            pickup_rate_coefficient=0.0,
            allow_inconsistent_mass_pairing=False,
        )
        ion_fixed = build_initial_ion_state(drag_cfg, neutral, num_steps_ion=5)
        np.testing.assert_allclose(
            ion_full.E_pot_eV[:, 0] - ion_fixed.E_pot_eV[:, 0],
            _fold_eV(bip, ANCHOR_N_START),
            rtol=0, atol=1e-12,
        )

    def test_density_tied_per_ion_hand_oracle(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        two_N = 2 * cfg.num_molecules
        ckpt = _neutral_with_radii(neutral, np.linspace(0.0, 26.0, two_N))
        bip = _biphasic_cfg(cfg, initial_shell_model="density_tied")
        ion = build_initial_ion_state(bip, ckpt, num_steps_ion=5)
        n0 = _expected_n0(bip, ckpt)
        # The dressing genuinely moves: surface births are under-dressed.
        assert n0.min() < ANCHOR_N_START
        np.testing.assert_allclose(ion.mass_kg, complex_mass_amu(n0) * U)
        np.testing.assert_array_equal(ion.n_shell[:, 0], n0)
        np.testing.assert_allclose(ion.mass_history_kg[:, 0], ion.mass_kg)
        np.testing.assert_allclose(ion.mass_final_kg, ion.mass_kg)

    def test_density_tied_fold_is_per_ion(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        two_N = 2 * cfg.num_molecules
        ckpt = _neutral_with_radii(neutral, np.linspace(0.0, 26.0, two_N))
        bip_tied = _biphasic_cfg(cfg, initial_shell_model="density_tied")
        bip_full = _biphasic_cfg(cfg)
        ion_tied = build_initial_ion_state(bip_tied, ckpt, num_steps_ion=5)
        ion_full = build_initial_ion_state(bip_full, ckpt, num_steps_ion=5)
        n0 = _expected_n0(bip_tied, ckpt)
        # E_pot is mass-independent, so the two arms differ exactly by the
        # per-ion vs scalar binding fold.
        np.testing.assert_allclose(
            ion_tied.E_pot_eV[:, 0] - ion_full.E_pot_eV[:, 0],
            _fold_eV(bip_tied, n0) - _fold_eV(bip_full, ANCHOR_N_START),
            rtol=0, atol=1e-12,
        )

    def test_E_kin_rides_the_dressed_mass(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        two_N = 2 * cfg.num_molecules
        ckpt = _neutral_with_radii(neutral, np.linspace(0.0, 26.0, two_N))
        bip_tied = _biphasic_cfg(cfg, initial_shell_model="density_tied")
        ion_tied = build_initial_ion_state(bip_tied, ckpt, num_steps_ion=5)
        ion_full = build_initial_ion_state(
            _biphasic_cfg(cfg), ckpt, num_steps_ion=5
        )
        n0 = _expected_n0(bip_tied, ckpt)
        # Same velocities, per-ion mass: E_kin scales as m(n_0)/m(21).
        np.testing.assert_allclose(
            ion_tied.E_kin_eV[:, 0] * complex_mass_amu(ANCHOR_N_START),
            ion_full.E_kin_eV[:, 0] * complex_mass_amu(n0),
            rtol=1e-12,
        )

    def test_center_pinned_births_are_inert(self, small_neutral_run):
        # H.3b oracle: rho_hat(center) ≈ 1 -> n_0 = 21 exactly; the dressed
        # arm is inert without the position axis. Atom 1 at +4.5 Å, atom 2 at
        # -4.5 Å (the center-pinned ±R0/2 molecule geometry — partners must
        # not coincide or the pair Coulomb diverges), droplet R ≈ 27.9 Å.
        cfg, neutral = small_neutral_run
        N = cfg.num_molecules
        ckpt = _neutral_with_radii(
            neutral, np.concatenate([np.full(N, 4.5), np.full(N, -4.5)])
        )
        bip_tied = _biphasic_cfg(cfg, initial_shell_model="density_tied")
        ion_tied = build_initial_ion_state(bip_tied, ckpt, num_steps_ion=5)
        ion_full = build_initial_ion_state(
            _biphasic_cfg(cfg), ckpt, num_steps_ion=5
        )
        np.testing.assert_array_equal(
            ion_tied.n_shell[:, 0], float(ANCHOR_N_START)
        )
        np.testing.assert_array_equal(ion_tied.mass_kg, ion_full.mass_kg)
        np.testing.assert_array_equal(
            ion_tied.E_pot_eV[:, 0], ion_full.E_pot_eV[:, 0]
        )

    def test_n0_monotone_nonincreasing_in_birth_radius(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        two_N = 2 * cfg.num_molecules
        ckpt = _neutral_with_radii(neutral, np.linspace(0.0, 26.0, two_N))
        bip = _biphasic_cfg(cfg, initial_shell_model="density_tied")
        ion = build_initial_ion_state(bip, ckpt, num_steps_ion=5)
        n0 = ion.n_shell[:, 0]  # radii are sorted ascending by construction
        assert np.all(np.diff(n0) <= 0.0)

    def test_E_int_onset_not_coupled_to_dressing(self, small_neutral_run):
        # T5/T6 boundary: the S2 onset stays the constant f_int * E_avail per
        # ion — the Sigma-ratio coupling is Slice T6's p-law, not this arm.
        cfg, neutral = small_neutral_run
        two_N = 2 * cfg.num_molecules
        ckpt = _neutral_with_radii(neutral, np.linspace(0.0, 26.0, two_N))
        bip = _biphasic_cfg(cfg, initial_shell_model="density_tied")
        ion = build_initial_ion_state(bip, ckpt, num_steps_ion=5)
        expected = (
            bip.internal_energy_partition_fraction * bip.coulomb_available_eV
        )
        np.testing.assert_allclose(ion.E_int_eV[:, 0], expected)


# ===========================================================================
# build_biphasic_cfg pass-through (scripts/tier2_common.py)
# ===========================================================================
class TestBuildBiphasicCfgPassThrough:
    @staticmethod
    def _build(**kwargs):
        from scripts.tier2_common import build_biphasic_cfg

        return build_biphasic_cfg(
            "9A",
            "shared_pure_cubic",
            num_molecules=2,
            ion_time_ps=0.02,
            dt_ion_ps=0.01,
            seed=123,
            **kwargs,
        )

    def test_none_rides_the_default(self):
        assert self._build().initial_shell_model == "full"

    def test_kwarg_stamps_the_arm(self):
        cfg = self._build(initial_shell_model="density_tied")
        assert cfg.initial_shell_model == "density_tied"
