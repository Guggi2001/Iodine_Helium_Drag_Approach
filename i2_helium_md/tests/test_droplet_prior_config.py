"""Config surface for the Slice-T8 droplet-size prior (plan §I.11.3, T8-D1).

The prior selection is a byte-inert enum ``droplet_size_prior``: ``legacy``
(default — the boolean-driven dispatch exactly as delivered: fixed N under
``use_single_droplet_size=True``, else the ported pickup-cell MC) vs the
analytic D4 family ``kornilov_lognormal`` / ``pickup_weighted_lognormal``
(T8-D2: mu_ln = ln<N> − delta²/2 about <N> = 2000, delta = 0.625 Kornilov;
pickup = × N^(2/3) un-re-centered), drawn exactly on the twin's truncation
window (T8-D3).

Guarded at config-load like every other string-enum surface
(``_reject_unknown_enum``), with the no-silent-inert rule in *both*
directions:

* analytic arms are refused under ``use_single_droplet_size=True`` (the
  fixed-size branch never reads the family — the T5/T6/E2
  ``_require_pairing`` precedent);
* non-default family parameters are refused under ``legacy`` (they would be
  silently ignored — the T7 margin-under-boltzmann precedent).
"""

import pytest

from i2_helium_md.config import (
    DROPLET_PRIOR_N_HI,
    DROPLET_PRIOR_N_LO,
    SimConfig,
    check_droplet_prior_config,
)


def test_default_is_legacy_with_d4_family_params():
    cfg = SimConfig()
    assert cfg.droplet_size_prior == "legacy"
    assert cfg.droplet_prior_mean_N == 2000.0
    assert cfg.droplet_prior_delta == 0.625


def test_field_defaults_back_missing_json_keys():
    # An old cfg.json predating the fields omits the keys; RunDirectory.load_cfg
    # reconstructs via SimConfig(**payload), so the dataclass defaults are the
    # back-compat guarantee (the relaxation_dissipation precedent).
    fields = SimConfig.__dataclass_fields__
    assert fields["droplet_size_prior"].default == "legacy"
    assert fields["droplet_prior_mean_N"].default == 2000.0
    assert fields["droplet_prior_delta"].default == 0.625


def test_unknown_prior_rejected():
    cfg = SimConfig(droplet_size_prior="bogus")
    with pytest.raises(ValueError, match="droplet_size_prior"):
        check_droplet_prior_config(cfg)


def test_kornilov_refused_under_single_droplet_size():
    # No-silent-inert: the fixed-size branch never samples, so the advertised
    # family would be dead config.
    cfg = SimConfig(droplet_size_prior="kornilov_lognormal")
    with pytest.raises(ValueError, match="use_single_droplet_size"):
        check_droplet_prior_config(cfg)


def test_pickup_weighted_refused_under_single_droplet_size():
    cfg = SimConfig(droplet_size_prior="pickup_weighted_lognormal")
    with pytest.raises(ValueError, match="use_single_droplet_size"):
        check_droplet_prior_config(cfg)


def test_analytic_arms_accepted_with_sampled_sizes():
    for arm in ("kornilov_lognormal", "pickup_weighted_lognormal"):
        cfg = SimConfig(droplet_size_prior=arm, use_single_droplet_size=False)
        check_droplet_prior_config(cfg)  # no raise


def test_non_default_delta_refused_under_legacy():
    # Silently-ignored parameter refusal (the T7 margin precedent): under
    # ``legacy`` the family parameters are never read.
    cfg = SimConfig(droplet_prior_delta=0.40)
    with pytest.raises(ValueError, match="droplet_prior_delta"):
        check_droplet_prior_config(cfg)


def test_non_default_mean_refused_under_legacy():
    cfg = SimConfig(droplet_prior_mean_N=3000.0)
    with pytest.raises(ValueError, match="droplet_prior_mean_N"):
        check_droplet_prior_config(cfg)


def test_nonpositive_delta_refused_under_analytic_arm():
    cfg = SimConfig(
        droplet_size_prior="kornilov_lognormal",
        use_single_droplet_size=False,
        droplet_prior_delta=0.0,
    )
    with pytest.raises(ValueError, match="droplet_prior_delta"):
        check_droplet_prior_config(cfg)


class TestSamplerMode:
    """Atlas G0-1: ``droplet_size_sampler_mode`` (``raw`` / ``post_pickup``).

    Before the field the mode was hardcoded ``"post_pickup"`` at the single
    call site, so the default is the byte-inert value and ``raw`` — the parent
    document's own ensemble (D0 §15.6) — becomes reachable. Guarded in the
    rule-3 direction: refused wherever the legacy *sampled* branch does not
    run, since nothing would read it there.
    """

    def test_default_is_post_pickup_the_byte_inert_value(self):
        assert SimConfig().droplet_size_sampler_mode == "post_pickup"
        fields = SimConfig.__dataclass_fields__
        assert fields["droplet_size_sampler_mode"].default == "post_pickup"

    def test_unknown_mode_rejected(self):
        cfg = SimConfig(droplet_size_sampler_mode="bogus",
                        use_single_droplet_size=False)
        with pytest.raises(ValueError, match="droplet_size_sampler_mode"):
            check_droplet_prior_config(cfg)

    def test_raw_accepted_on_the_legacy_sampled_branch(self):
        cfg = SimConfig(droplet_size_sampler_mode="raw",
                        use_single_droplet_size=False)
        check_droplet_prior_config(cfg)  # no raise

    def test_raw_refused_under_fixed_size(self):
        # The fixed-size branch never samples, so the mode would be inert.
        cfg = SimConfig(droplet_size_sampler_mode="raw",
                        use_single_droplet_size=True)
        with pytest.raises(ValueError, match="droplet_size_sampler_mode"):
            check_droplet_prior_config(cfg)

    def test_raw_refused_under_an_analytic_prior_arm(self):
        # The analytic arms bypass sample_droplet_sizes entirely.
        cfg = SimConfig(droplet_size_sampler_mode="raw",
                        droplet_size_prior="kornilov_lognormal",
                        use_single_droplet_size=False)
        with pytest.raises(ValueError, match="droplet_size_sampler_mode"):
            check_droplet_prior_config(cfg)

    def test_default_mode_stays_legal_everywhere(self):
        # The default must never trip the guard, or every existing cfg breaks.
        for prior, single in (("legacy", True), ("legacy", False),
                              ("kornilov_lognormal", False)):
            check_droplet_prior_config(
                SimConfig(droplet_size_prior=prior,
                          use_single_droplet_size=single)
            )


def test_mean_outside_truncation_window_refused():
    # The truncation window [250, 16000] is the twin's proposal support
    # (T8-D2); a mean outside it would silently discard most of the prior
    # mass (no-silent-caps).
    assert DROPLET_PRIOR_N_LO == 250.0 and DROPLET_PRIOR_N_HI == 16000.0
    for bad_mean in (100.0, 20000.0):
        cfg = SimConfig(
            droplet_size_prior="kornilov_lognormal",
            use_single_droplet_size=False,
            droplet_prior_mean_N=bad_mean,
        )
        with pytest.raises(ValueError, match="droplet_prior_mean_N"):
            check_droplet_prior_config(cfg)
