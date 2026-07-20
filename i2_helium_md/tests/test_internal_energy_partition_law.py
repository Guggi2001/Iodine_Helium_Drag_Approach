"""Slice T6 (Tier-2 plan §I.11; the D2 p-law): the E_int(0)-dressing coupling.

Config enum ``internal_energy_partition_law ∈ {"constant" (byte-inert default,
p = 0), "sigma_proportional" (p = 1)}``. Under ``sigma_proportional`` the S2
Coulomb onset rides the T5 initial-shell dressing through the ladder-resolved
ratio:

    E_int(0)_i = f_int · E_avail · (Σ(n_0i) / Σ(n*))^p ,

with Σ = ``dissociation_ladder.ladder_cumsum`` at the cfg picture/kappa/ladder
and n* = ``N_STAR`` = 21. Structurally inert at n_0 = n* (ratio = 1) regardless
of law, so the arm is a no-op without the T5 position/dressing axis; the default
``constant`` arm is byte-inert (factor ≡ 1.0). Under-dressed births (n_0 < n*)
get *less* onset — the §4j p = 1 direction (p = 0 over-suppresses).

Boundary: T6 touches the E_int(0) onset **only**. The per-ion t0 mass, n_shell,
and the E_pot binding fold are Slice T5's and are unchanged by the law.
"""

from __future__ import annotations

import json
from dataclasses import fields, replace

import numpy as np
import pytest

from i2_helium_md.config import (
    SimConfig,
    check_internal_energy_partition_config,
)
from i2_helium_md.physics.constants import N_STAR, U
from i2_helium_md.physics.dissociation_ladder import ladder_cumsum, resolve_ladder
from i2_helium_md.physics.internal_energy_budget import sigma_partition_factor
from i2_helium_md.physics.shell_schedule import ANCHOR_N_START, complex_mass_amu
from i2_helium_md.presets import single_pulse_N2000, single_pulse_N2000_drag
from i2_helium_md.simulation.checkpoint import NeutralCheckpoint
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
    """A valid Tier-2 biphasic config on the Tier-0 drag bundle (T5 idiom)."""
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
    """NeutralCheckpoint with all numpy arrays deep-copied (T5 idiom)."""
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


def _sigma(cfg, n):
    """Σ(n) at the cfg's resolved ladder (scalar or array n) — the p-law ratio
    numerator/denominator; the independent oracle."""
    return ladder_cumsum(
        n,
        picture=cfg.ladder_electronic_picture,
        kappa=cfg.ladder_steepness,
        ladder=resolve_ladder(
            cfg.dissociation_ladder, cfg.tabulated_ladder_rungs_eV
        ),
    )


# ===========================================================================
# The pure p-law factor
# ===========================================================================
class TestSigmaPartitionFactor:
    _KW = dict(kappa=1.0)  # Form-U default picture; ladder=None (byte-inert path)

    def test_constant_is_unity_scalar(self):
        f = sigma_partition_factor(19, law="constant", **self._KW)
        assert f == 1.0
        assert np.ndim(f) == 0

    def test_constant_is_unity_for_array_and_ignores_n0(self):
        # p = 0: the factor is the multiplicative identity regardless of n0,
        # so the delivered constant onset is preserved byte-for-byte.
        f = sigma_partition_factor(
            np.array([21, 14, 5, 0]), law="constant", **self._KW
        )
        assert f == 1.0
        assert np.ndim(f) == 0

    def test_sigma_proportional_scalar_hand_oracle(self):
        # Σ(n0)/Σ(n*), the ladder-resolved ratio.
        n0 = 14
        got = sigma_partition_factor(n0, law="sigma_proportional", **self._KW)
        expected = ladder_cumsum(n0, **self._KW) / ladder_cumsum(N_STAR, **self._KW)
        assert got == pytest.approx(expected, rel=0, abs=1e-15)
        assert np.ndim(got) == 0

    def test_sigma_proportional_inert_at_nstar(self):
        # ratio = Σ(n*)/Σ(n*) = 1 exactly — inert without the T5 dressing axis.
        f = sigma_partition_factor(N_STAR, law="sigma_proportional", **self._KW)
        assert f == 1.0

    def test_sigma_proportional_array_hand_oracle(self):
        n0 = np.array([21, 19, 14, 7, 2])
        got = sigma_partition_factor(n0, law="sigma_proportional", **self._KW)
        expected = ladder_cumsum(n0, **self._KW) / ladder_cumsum(N_STAR, **self._KW)
        np.testing.assert_allclose(got, expected, rtol=0, atol=1e-15)
        assert isinstance(got, np.ndarray)

    def test_sigma_proportional_under_dressed_gives_less_onset(self):
        # p = 1: n0 < n* -> factor < 1 (the §4j "p = 0 over-suppresses"
        # direction); strictly increasing in n0.
        n0 = np.array([2, 7, 14, 19, 21])
        f = sigma_partition_factor(n0, law="sigma_proportional", **self._KW)
        assert np.all(f[:-1] < f[1:])          # strictly increasing in n0
        assert np.all(f[:-1] < 1.0)            # under-dressed -> < 1
        assert f[-1] == 1.0                    # n0 = n* -> 1

    def test_unknown_law_refused(self):
        with pytest.raises(ValueError, match="internal_energy_partition_law"):
            sigma_partition_factor(19, law="sigma_propotional", **self._KW)


# ===========================================================================
# Config surface
# ===========================================================================
class TestConfigSurface:
    def test_default_is_constant(self):
        assert SimConfig().internal_energy_partition_law == "constant"

    def test_default_config_still_validates(self):
        SimConfig().validate()

    def test_unknown_selector_refused(self, small_neutral_run):
        cfg, _ = small_neutral_run
        bad = _biphasic_cfg(cfg, internal_energy_partition_law="sigma_propotional")
        with pytest.raises(ValueError, match="internal_energy_partition_law"):
            bad.validate()

    def test_sigma_proportional_requires_biphasic(self):
        # mass_scenario defaults to "fixed": the p-law would be silently inert
        # there (the onset seed only runs under biphasic) — refused loudly
        # instead (the T5/T7 no-silent-inert convention).
        cfg = SimConfig(internal_energy_partition_law="sigma_proportional")
        with pytest.raises(ValueError, match="biphasic"):
            check_internal_energy_partition_config(cfg)

    def test_sigma_proportional_under_biphasic_validates(self, small_neutral_run):
        cfg, _ = small_neutral_run
        _biphasic_cfg(
            cfg, internal_energy_partition_law="sigma_proportional"
        ).validate()

    def test_pre_slice_cfg_json_loads_with_default(self, tmp_path):
        # Slice-DS/T5/T7 precedent: a cfg.json written before the field existed
        # must load with the byte-inert default.
        run = RunDirectory(tmp_path / "9A_drag_tier2probe_backcompat_plaw")
        run.save_cfg(SimConfig())
        payload = json.loads(run.cfg_path.read_text(encoding="utf-8"))
        assert payload.pop("internal_energy_partition_law") == "constant"
        run.cfg_path.write_text(json.dumps(payload), encoding="utf-8")
        assert run.load_cfg().internal_energy_partition_law == "constant"

    def test_cfg_json_roundtrip(self, tmp_path, small_neutral_run):
        cfg, _ = small_neutral_run
        run = RunDirectory(tmp_path / "9A_drag_tier2probe_plaw_roundtrip")
        run.save_cfg(
            _biphasic_cfg(cfg, internal_energy_partition_law="sigma_proportional")
        )
        assert (
            run.load_cfg().internal_energy_partition_law == "sigma_proportional"
        )


# ===========================================================================
# The coupled onset seed (build_initial_ion_state)
# ===========================================================================
class TestCoupledOnsetSeed:
    def test_constant_arm_onset_is_delivered_value(self, small_neutral_run):
        # Byte-inert default: E_int(0) = f_int * E_avail per ion, even when the
        # T5 dressing is live (the p = 0 boundary preserved).
        cfg, neutral = small_neutral_run
        two_N = 2 * cfg.num_molecules
        ckpt = _neutral_with_radii(neutral, np.linspace(0.0, 26.0, two_N))
        bip = _biphasic_cfg(cfg, initial_shell_model="density_tied")
        ion = build_initial_ion_state(bip, ckpt, num_steps_ion=5)
        expected = (
            bip.internal_energy_partition_fraction * bip.coulomb_available_eV
        )
        np.testing.assert_allclose(ion.E_int_eV[:, 0], expected)

    def test_sigma_proportional_full_arm_is_inert(self, small_neutral_run):
        # n0 = 21 for all ions under `full` -> ratio 1 -> onset unchanged.
        cfg, neutral = small_neutral_run
        bip = _biphasic_cfg(cfg, internal_energy_partition_law="sigma_proportional")
        ion = build_initial_ion_state(bip, neutral, num_steps_ion=5)
        expected = (
            bip.internal_energy_partition_fraction * bip.coulomb_available_eV
        )
        np.testing.assert_allclose(ion.E_int_eV[:, 0], expected)

    def test_sigma_proportional_dressed_per_ion_hand_oracle(self, small_neutral_run):
        cfg, neutral = small_neutral_run
        two_N = 2 * cfg.num_molecules
        ckpt = _neutral_with_radii(neutral, np.linspace(0.0, 26.0, two_N))
        bip = _biphasic_cfg(
            cfg,
            initial_shell_model="density_tied",
            internal_energy_partition_law="sigma_proportional",
        )
        ion = build_initial_ion_state(bip, ckpt, num_steps_ion=5)
        n0 = ion.n_shell[:, 0]  # the T5 dressed shell (already covered by T5 tests)
        onset0 = (
            bip.internal_energy_partition_fraction * bip.coulomb_available_eV
        )
        expected = onset0 * _sigma(bip, n0) / _sigma(bip, N_STAR)
        np.testing.assert_allclose(ion.E_int_eV[:, 0], expected, rtol=0, atol=1e-15)
        # The dressing genuinely moves: under-dressed births get less onset.
        assert n0.min() < N_STAR
        assert ion.E_int_eV[:, 0].min() < onset0

    def test_law_touches_only_E_int_not_mass_or_fold(self, small_neutral_run):
        # T5/T6 boundary the other way: flipping the p-law leaves the per-ion
        # mass, n_shell, and the E_pot binding fold identical — only E_int moves.
        cfg, neutral = small_neutral_run
        two_N = 2 * cfg.num_molecules
        ckpt = _neutral_with_radii(neutral, np.linspace(0.0, 26.0, two_N))
        ion_const = build_initial_ion_state(
            _biphasic_cfg(cfg, initial_shell_model="density_tied"),
            ckpt, num_steps_ion=5,
        )
        ion_prop = build_initial_ion_state(
            _biphasic_cfg(
                cfg,
                initial_shell_model="density_tied",
                internal_energy_partition_law="sigma_proportional",
            ),
            ckpt, num_steps_ion=5,
        )
        np.testing.assert_array_equal(ion_prop.mass_kg, ion_const.mass_kg)
        np.testing.assert_array_equal(
            ion_prop.n_shell[:, 0], ion_const.n_shell[:, 0]
        )
        np.testing.assert_array_equal(
            ion_prop.E_pot_eV[:, 0], ion_const.E_pot_eV[:, 0]
        )
        # E_int strictly differs on the under-dressed ions.
        assert np.any(ion_prop.E_int_eV[:, 0] < ion_const.E_int_eV[:, 0])


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
        assert self._build().internal_energy_partition_law == "constant"

    def test_kwarg_stamps_the_arm(self):
        cfg = self._build(internal_energy_partition_law="sigma_proportional")
        assert cfg.internal_energy_partition_law == "sigma_proportional"
