"""Slice T7 (Tier-2 plan SI.11): the birth-position-law arm.

``birth_position_law = "uniform_volume"`` realizes the 1D twin's L1 ensemble
law -- p(r) proportional to r^2 on [0, R - margin], zero outside (hard
margin, no smoothing; user decision 2026-07-16) -- as an interchangeable
arm. ``"boltzmann"`` (the delivered ``sample_radial_positions`` rejection
sampler) stays the byte-inert default; its draw path must not change (RNG
draw-order rule), which the existing sampler suites regression-guard.

Tolerances: distribution checks use n = 20000 fixed-seed draws; on
u = (r / cap)^3 (analytically Uniform(0,1) under the r^2 law) the quantile
sampling band at that n is well under +/-0.02, which is the assert width.
Identity and guard checks are exact.
"""

from __future__ import annotations

import dataclasses
import json

import numpy as np
import pytest

from i2_helium_md.config import SimConfig, check_birth_position_config
from i2_helium_md.sampling.radial_positions import sample_radial_positions
from i2_helium_md.simulation.run_directory import RunDirectory


def _uniform_cfg(margin: float = 0.0, **kw) -> SimConfig:
    # uniform_volume is refused under single_initial_position=True (the
    # SimConfig default zeroes r0, leaving the law silently inert -- review
    # fix 2026-07-18), so the fixture selects the position-live arm.
    kw.setdefault("single_initial_position", False)
    return SimConfig(
        birth_position_law="uniform_volume",
        initial_position_margin_angstrom=margin,
        **kw,
    )


# ---------------------------------------------------------------------------
# Config surface: defaults, guards, back-compat
# ---------------------------------------------------------------------------
class TestConfigSurface:
    def test_defaults_are_byte_inert(self):
        cfg = SimConfig()
        assert cfg.birth_position_law == "boltzmann"
        assert cfg.initial_position_margin_angstrom == 0.0
        cfg.validate()  # default surface passes untouched

    def test_unknown_law_rejected(self):
        cfg = SimConfig(birth_position_law="nope")
        with pytest.raises(ValueError, match="birth_position_law"):
            cfg.validate()

    def test_margin_under_boltzmann_refused(self):
        # No silent carry: the margin is uniform_volume-only (plan T7).
        cfg = SimConfig(initial_position_margin_angstrom=3.0)
        with pytest.raises(ValueError, match="uniform_volume"):
            cfg.validate()

    def test_uniform_volume_under_single_position_refused(self):
        # review fix 2026-07-18: under single_initial_position=True the
        # sampled r0 is zeroed, so the advertised law would be silently
        # inert on positions while still shifting the RNG draw stream --
        # the combination is refused loudly at config-load.
        cfg = _uniform_cfg(single_initial_position=True)
        with pytest.raises(ValueError, match="single_initial_position"):
            cfg.validate()

    @pytest.mark.parametrize("bad", [-1.0, float("nan"), float("inf")])
    def test_bad_margin_rejected(self, bad):
        cfg = _uniform_cfg(margin=bad)
        with pytest.raises(ValueError, match="margin"):
            cfg.validate()

    def test_firm_band_margins_pass(self):
        for margin in (0.0, 3.0, 4.67, 6.0):
            check_birth_position_config(_uniform_cfg(margin=margin))

    def test_pre_t7_cfg_json_loads_with_defaults(self, tmp_path):
        # Slice-DS precedent: a cfg.json written before the two T7 fields
        # existed must load with the byte-inert defaults.
        run = RunDirectory(tmp_path / "9A_drag_tier2probe_backcompat_t7")
        run.save_cfg(SimConfig())
        payload = json.loads(run.cfg_path.read_text(encoding="utf-8"))
        assert payload.pop("birth_position_law") == "boltzmann"
        assert payload.pop("initial_position_margin_angstrom") == 0.0
        run.cfg_path.write_text(json.dumps(payload), encoding="utf-8")

        cfg = run.load_cfg()
        assert cfg.birth_position_law == "boltzmann"
        assert cfg.initial_position_margin_angstrom == 0.0

    def test_cfg_json_roundtrip(self, tmp_path):
        run = RunDirectory(tmp_path / "9A_drag_tier2probe_t7_roundtrip")
        cfg = _uniform_cfg(margin=4.67)
        run.save_cfg(cfg)
        loaded = run.load_cfg()
        assert loaded.birth_position_law == "uniform_volume"
        assert loaded.initial_position_margin_angstrom == 4.67


# ---------------------------------------------------------------------------
# The uniform_volume law itself
# ---------------------------------------------------------------------------
class TestUniformVolumeLaw:
    R = 27.936  # R(N=2000), the probe convention

    def _draw(self, margin: float, n: int = 20000, seed: int = 7):
        cfg = _uniform_cfg(margin=margin, seed=seed)
        radii = np.full(n, self.R)
        return sample_radial_positions(
            cfg, radii, rng=np.random.default_rng(seed)
        )

    def test_r_squared_law_margin_zero(self):
        r = self._draw(margin=0.0)
        u = (r / self.R) ** 3  # analytically Uniform(0, 1)
        for f in (0.05, 0.25, 0.5, 0.75, 0.95):
            assert abs(np.quantile(u, f) - f) < 0.02
        # margin 0 reaches the surface (plan T7 oracle)
        assert r.max() > 0.995 * self.R
        assert r.max() <= self.R

    def test_hard_margin_support(self):
        margin = 4.67
        cap = self.R - margin
        r = self._draw(margin=margin)
        assert r.max() <= cap
        assert r.max() > 0.995 * cap  # hard edge, not a soft rolloff
        u = (r / cap) ** 3
        for f in (0.25, 0.5, 0.75):
            assert abs(np.quantile(u, f) - f) < 0.02

    def test_median_matches_twin_comparator(self):
        # The committed twin's uniform_volume comparator row at m = 4.67:
        # median = (R - m) * 0.5^(1/3) = 18.47 A (V0-1 quantiles CSV).
        r = self._draw(margin=4.67)
        assert abs(np.median(r) - 18.47) < 0.15  # n=20000 sampling band

    def test_margin_at_least_radius_raises(self):
        with pytest.raises(ValueError, match="margin"):
            self._draw(margin=self.R)

    def test_per_droplet_radii_respected(self):
        # Mixed droplet sizes: each molecule's support is its own droplet's.
        cfg = _uniform_cfg(margin=3.0, seed=11)
        radii = np.array([20.0] * 5000 + [30.0] * 5000)
        r = sample_radial_positions(cfg, radii, rng=np.random.default_rng(11))
        assert r[:5000].max() <= 17.0
        assert r[5000:].max() <= 27.0
        assert r[5000:].max() > 17.0  # the larger droplet uses its range

    def test_differs_from_boltzmann(self):
        # Liveness: the arm actually changes the law (boltzmann is
        # center-pinned at T = 0.4 K; uniform_volume median ~ 0.8 R).
        r_uni = self._draw(margin=0.0, n=2000)
        cfg_b = SimConfig(seed=7)
        r_bol = sample_radial_positions(
            cfg_b, np.full(2000, self.R), rng=np.random.default_rng(7)
        )
        assert np.median(r_bol) < 3.0
        assert np.median(r_uni) > 20.0


# ---------------------------------------------------------------------------
# Integration: the neutral initial state consumes the arm
# ---------------------------------------------------------------------------
class TestInitialStateIntegration:
    def test_build_initial_state_uniform_volume(self):
        from i2_helium_md.simulation.initial_state import build_initial_state

        margin = 4.67
        cfg = dataclasses.replace(
            _uniform_cfg(margin=margin),
            num_molecules=200,
            single_initial_position=False,
            seed=3,
        )
        ckpt = build_initial_state(cfg, num_steps=2)
        x = ckpt.positions_x[:, 0]
        y = ckpt.positions_y[:, 0]
        z = ckpt.positions_z[:, 0]
        # molecule centre radius: the two atoms sit symmetrically about the
        # centre, so the pair midpoint radius is the sampled r0
        N = cfg.num_molecules
        cx = (x[:N] + x[N:]) / 2
        cy = (y[:N] + y[N:]) / 2
        cz = (z[:N] + z[N:]) / 2
        r0 = np.sqrt(cx**2 + cy**2 + cz**2)
        R = 27.936  # R(N=2000) bulk convention
        assert r0.max() <= (R - margin) + 1e-9
        assert np.median(r0) > 10.0  # not center-pinned
