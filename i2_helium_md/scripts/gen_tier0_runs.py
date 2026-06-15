"""Generate a Tier-0 drag run for any Method-B form/parameter set.

Produces a full-duration drag run at a chosen ``N`` for one droplet geometry
(9 A or 18 A) wired with a committed Method-B (``trajectory_matching``) bundle
selected from the named catalog in :mod:`scripts.tier0_common`, with an optional
inline hand-tuning override. Not a committed artifact -- ``data/runs`` is
gitignored; the small mean-series CSV exported by ``tier0_drag_comparison.py`` is
the committed regression reference.

The drag law is chosen by ``VARIANT`` (a catalog key), independent of the
production drag presets (which are frozen to ``linear_cubic`` /
``shared_pure_cubic``). With the defaults below the run is identical to
``single_pulse_N2000_18Angst_drag``.

Usage::

    # edit USER SETTINGS, then:
    python scripts/gen_tier0_runs.py
"""

from __future__ import annotations

from pathlib import Path
import sys

# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"      # droplet geometry: "9A" or "18A"

# Drag law: a key of scripts.tier0_common.CATALOG. Shared keys run on either
# geometry; per-case keys (percase_*) resolve against CASE.
#   shared_pure_cubic  shared_3param  shared_pl  shared_lq  shared_lq_pq
#   percase_linear_cubic  percase_pl  percase_lq  percase_lq_pq
VARIANT = "percase_linear_cubic"

N = 50            # ensemble size, e.g. 50 or 2000

ION_TIME_PS = 30.0
DT_ION_PS = 0.01
SEED = 20260604

# Optional inline hand-tuning (Method-B). Leave empty / None for the committed
# bundle exactly. Override keys must be coefficients of the chosen form
# (linear_cubic {a,b}, linear_quadratic {a,c}, power_law {C,n}).
COEFF_OVERRIDES: dict[str, float] = {}
E_BIND_OVERRIDE = None  # eV, or None to use the bundle's stamped binding

# Marker appended to the run-dir name for a tuned run, so it never overwrites the
# committed-variant run. Set the SAME value in tier0_drag_comparison.py to score
# it. Leave "" for an untuned run.
RUN_TAG = ""

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tier0_common import build_drag_cfg, run_dir_name  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

# =============================================================================
# MAIN SCRIPT
# =============================================================================

cfg = build_drag_cfg(
    CASE,
    VARIANT,
    num_molecules=N,
    ion_time_ps=ION_TIME_PS,
    dt_ion_ps=DT_ION_PS,
    seed=SEED,
    coeff_overrides=COEFF_OVERRIDES,
    e_bind_override=E_BIND_OVERRIDE,
)
cfg.validate()

tuned = bool(COEFF_OVERRIDES) or E_BIND_OVERRIDE is not None
run_dir = PROJECT_ROOT / "data" / "runs" / run_dir_name(CASE, VARIANT, N, RUN_TAG)
run = RunDirectory(run_dir)
run.save_cfg(cfg)

print(
    f"[{CASE} {VARIANT} N={N}] form={cfg.drag_form} "
    f"coeffs={dict(cfg.drag_coefficients.coefficients)} "
    f"E_bind={cfg.binding_energy_I_ion_eV:.4f} eV"
    + ("  (HAND-TUNED)" if tuned else "")
)

print(f"[{CASE} {VARIANT} N={N}] neutral propagation ...")
neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)

print(f"[{CASE} {VARIANT} N={N}] ion propagation ...")
ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)

print(
    f"[{CASE} {VARIANT} N={N}] done -> {run_dir}  "
    f"(neutral {neutral.time_ps.size} steps, ion {ion.time_ps.size} steps)"
)
