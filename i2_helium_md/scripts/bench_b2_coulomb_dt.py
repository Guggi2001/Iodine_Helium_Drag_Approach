"""Deterministic B2 Coulomb-explosion timestep-refinement benchmark (audit F1).

One I+-I+ pair, born at rest at R0 = 2.666 A at depth -35 A inside a
53.47 A droplet (the B2 s1 median geometry), propagated 40 ps with the
**production** conservative force (``make_ion_accel_fn``: FD Coulomb +
erf droplet well) and the **production** BAOAB stepper
(``make_ion_baoab_step``) under the B2 drag law (bundle
``srho3p5A18A_linear_quadratic``, a 18.007, c 1.6, 3.5 A gate) or with
gamma = 0. Fixed fragment mass n = 21 (``complex_mass_amu(21)``), no
pickup / evaporation / cooling events, noise off (T_eff = 0). Pure
deterministic numerics -- NOT the stochastic B2 model.

Geometries: ``radial`` (axis along the droplet radius: one ion outward,
one inward through the droplet) and ``tangential`` (axis perpendicular;
both ions equivalent). Each with drag off / on. dt series 0.02, 0.01,
0.005, 0.0025, 0.00125 ps; the finest step is a numerical reference,
not the exact solution.

The cfg is the B2 production cfg of member s1, built by the generator
(``gen_tier2_detfix_postpickup.build_lq_member``) -- only its force,
drag and gate fields are read.

Usage::

    python scripts/bench_b2_coulomb_dt.py                 # write the CSV
    python scripts/bench_b2_coulomb_dt.py --check         # rerun, compare to committed CSV
    python scripts/bench_b2_coulomb_dt.py --session-original
        # recovery mode: the 2026-10-20 session check's mass (210.95858 amu,
        # I atomic mass - m_e + 21 He CODATA) and eV constant 9648.533212;
        # prints, writes nothing (used to confirm the recovered script)

Output: ``data/runs/b2_coulomb_dt_benchmark/b2_coulomb_dt_benchmark.csv``
(force-add; ``data/runs/`` is gitignored). Note:
``docs/drag_port/Tier2/B2_COULOMB_DT_BENCHMARK.md``.
"""

from __future__ import annotations

import argparse
import csv
import os
from functools import partial
from pathlib import Path
import sys
import warnings

for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from i2_helium_md.config import drag_gate_steepness  # noqa: E402
from i2_helium_md.physics.baoab import make_ion_baoab_step  # noqa: E402
from i2_helium_md.physics.constants import U  # noqa: E402
from i2_helium_md.physics.drag import drag_gamma  # noqa: E402
from i2_helium_md.physics.leapfrog import make_ion_accel_fn  # noqa: E402
from i2_helium_md.physics.shell_schedule import complex_mass_amu  # noqa: E402
from i2_helium_md.simulation.ion_propagation_step import _amu_ang2_ps2_to_eV  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================
R_DROPLET_ANGSTROM = 53.47
BIRTH_DEPTH_ANGSTROM = -35.0     # session used -34.99; see the note
R0_ANGSTROM = 2.666
N_SHELL = 21
T_END_PS = 40.0
DT_SERIES_PS = (0.02, 0.01, 0.005, 0.0025, 0.00125)
CSV_PATH = PROJECT_ROOT / "data" / "runs" / "b2_coulomb_dt_benchmark" / "b2_coulomb_dt_benchmark.csv"
# Recovery mode (the original session script's literals).
SESSION_MASS_AMU = 126.90447 - 0.000548579909 + 21 * 4.002602
SESSION_EV_IN_AMU_A2_PS2 = 9648.533212
SESSION_BIRTH_DEPTH_ANGSTROM = -34.99

COLUMNS = ("geometry", "drag", "ion", "dt_ps", "n_steps", "speed_A_per_ps",
           "KE_eV", "dKE_vs_finest_eV", "dKE_rel_vs_finest",
           "ratio_vs_finest", "increment_eV", "increment_ratio", "depth_final_A", "E_dissip_pair_eV",
           "E_pair_coulomb_final_eV")


def b2_cfg():
    """The committed B2 production cfg (member s1) via its generator."""
    from scripts.gen_tier2_detfix_postpickup import (
        B2_C, B2_E0_EV, B2_TAU_PS, build_lq_member,
    )
    with warnings.catch_warnings():
        # The biphasic<->constant-coefficient pairing warning is the known,
        # allowed production pairing (allow_inconsistent_mass_pairing=True).
        warnings.simplefilter("ignore", RuntimeWarning)
        cfg = build_lq_member("s1", B2_E0_EV, B2_TAU_PS, B2_C)
    coeffs = cfg.drag_coefficients.coefficients
    if (cfg.drag_form != "linear_quadratic" or coeffs["c"] != 1.6
            or cfg.drag_gate_steepness != 3.5 or cfg.noise_form != "none"):
        raise AssertionError(f"not the B2 law: {cfg.drag_form} {coeffs} "
                             f"gate {cfg.drag_gate_steepness} noise {cfg.noise_form}")
    return cfg


def run_pair(cfg, dt, *, drag, axis, m_amu, birth_depth, to_eV):
    """Propagate one pair for T_END_PS at step dt. Returns per-ion arrays
    (speed [A/ps], KE [eV], final depth [A]) and pair scalars
    (E_dissip [eV], final pair Coulomb E_pot [eV], step count)."""
    centre = np.array([R_DROPLET_ANGSTROM + birth_depth, 0.0, 0.0])
    u = np.array([1.0, 0.0, 0.0]) if axis == "radial" else np.array([0.0, 1.0, 0.0])
    p1, p2 = centre + 0.5 * R0_ANGSTROM * u, centre - 0.5 * R0_ANGSTROM * u
    pos = (np.array([p1[0], p2[0]]), np.array([p1[1], p2[1]]), np.array([p1[2], p2[2]]))
    vel = (np.zeros(2), np.zeros(2), np.zeros(2))
    radii = np.full(2, R_DROPLET_ANGSTROM)
    acc = make_ion_accel_fn(cfg, np.full(2, m_amu * U), radii, np.ones(2))
    if drag:
        gamma = partial(drag_gamma, coeffs=cfg.drag_coefficients,
                        steepness=drag_gate_steepness(cfg))
    else:
        def gamma(speed, depth):
            return np.zeros_like(speed)
    step = make_ion_baoab_step(np.full(2, m_amu), radii, acc, gamma, T_eff=0.0)
    n_steps = int(round(T_END_PS / dt))
    e_dissip = np.zeros(2)
    e_pot = np.zeros(1)
    for _ in range(n_steps):
        pos, vel, e_pot, de = step(pos, vel, dt)
        e_dissip += de
    speed = np.sqrt(vel[0] ** 2 + vel[1] ** 2 + vel[2] ** 2)
    ke = to_eV(0.5 * m_amu * speed ** 2)
    depth = np.sqrt(pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2) - R_DROPLET_ANGSTROM
    return speed, ke, depth, float(to_eV(e_dissip).sum()), float(e_pot[0]), n_steps


def benchmark(*, m_amu, birth_depth, to_eV):
    """All rows (dicts keyed by COLUMNS), deterministic."""
    cfg = b2_cfg()
    rows = []
    for axis in ("radial", "tangential"):
        ions = ("outward", "inward") if axis == "radial" else ("tangential_1", "tangential_2")
        for drag in (False, True):
            res = {dt: run_pair(cfg, dt, drag=drag, axis=axis, m_amu=m_amu,
                                birth_depth=birth_depth, to_eV=to_eV)
                   for dt in DT_SERIES_PS}
            fine = res[min(DT_SERIES_PS)][1]
            for k, ion in enumerate(ions):
                diff = {dt: float(res[dt][1][k] - fine[k]) for dt in DT_SERIES_PS}
                # Reference-free increment KE(dt/2) - KE(dt) and its ratio to the
                # next one (inc(dt) / inc(dt/2)); blank where dt/2 or dt/4 is absent.
                inc = {dt: float(res[dt / 2][1][k] - res[dt][1][k])
                       for dt in DT_SERIES_PS if dt / 2 in res}
                for dt in DT_SERIES_PS:
                    half = dt / 2
                    ratio = ("" if half not in diff or diff[half] == 0.0
                             else diff[dt] / diff[half])
                    inc_ratio = ("" if dt not in inc or half not in inc
                                 else inc[dt] / inc[half])
                    sp, ke, dep, ed, ep, ns = res[dt]
                    rows.append(dict(
                        geometry=axis, drag="B2" if drag else "off", ion=ion,
                        dt_ps=dt, n_steps=ns, speed_A_per_ps=float(sp[k]),
                        KE_eV=float(ke[k]), dKE_vs_finest_eV=diff[dt],
                        dKE_rel_vs_finest=diff[dt] / float(fine[k]),
                        ratio_vs_finest=ratio, increment_eV=inc.get(dt, ""),
                        increment_ratio=inc_ratio, depth_final_A=float(dep[k]),
                        E_dissip_pair_eV=ed, E_pair_coulomb_final_eV=ep,
                    ))
    return rows


def _fmt(v):
    return repr(v) if isinstance(v, float) else str(v)


def write_csv(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([_fmt(r[c]) for c in COLUMNS])


def check_csv(rows, path):
    """Exact (string-repr) comparison of a fresh run against the CSV."""
    with open(path, newline="", encoding="utf-8") as fh:
        ref = list(csv.reader(fh))
    new = [list(COLUMNS)] + [[_fmt(r[c]) for c in COLUMNS] for r in rows]
    if ref != new:
        bad = [i for i, (a, b) in enumerate(zip(ref, new)) if a != b]
        raise SystemExit(f"CHECK FAILED: {len(bad)} differing rows "
                         f"(first {bad[:5]}); row counts {len(ref)} / {len(new)}")
    print(f"CHECK PASSED: {len(rows)} rows reproduce {path} exactly.")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="rerun and compare to the CSV")
    ap.add_argument("--session-original", action="store_true",
                    help="recovery mode with the original session literals; print only")
    args = ap.parse_args(argv)
    if args.session_original:
        rows = benchmark(m_amu=SESSION_MASS_AMU,
                         birth_depth=SESSION_BIRTH_DEPTH_ANGSTROM,
                         to_eV=lambda e: e / SESSION_EV_IN_AMU_A2_PS2)
        for r in rows:
            print(r["geometry"], r["drag"], r["ion"], r["dt_ps"],
                  f"speed={r['speed_A_per_ps']:.8f}", f"KE={r['KE_eV']:.8f}",
                  f"dKE={r['dKE_vs_finest_eV']:.8e}")
        return 0
    rows = benchmark(m_amu=complex_mass_amu(N_SHELL),
                     birth_depth=BIRTH_DEPTH_ANGSTROM, to_eV=_amu_ang2_ps2_to_eV)
    if args.check:
        check_csv(rows, CSV_PATH)
    else:
        write_csv(rows, CSV_PATH)
        print(f"wrote {len(rows)} rows -> {CSV_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
