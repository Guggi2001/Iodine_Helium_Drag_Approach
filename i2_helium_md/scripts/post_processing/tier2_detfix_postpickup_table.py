"""Droplet-size sampler A/B at h405 — scorer (raw production vs post_pickup).

Scores the ``gen_tier2_detfix_postpickup.py`` members against the production
battery (``gen_tier2_detfix_battery.py``). The two batteries share seeds and
pipeline and differ in exactly ``droplet_size_sampler_mode`` (raw → post_pickup),
but they are **not CRN** (the pickup Monte Carlo consumes the shared RNG stream
first), so the comparison is pooled battery vs pooled battery: s1–s5, N = 5000
each (r6 excluded — its post_pickup run tripped the D4 safeguard; see the
generator's ``POOL_MEMBERS``).

Pure scorer — runs no MD and writes only the CSV below. Observable vector =
the detfix battery scorer's ``_row`` (``observable_columns`` + ``lowke_columns``
+ ``md_gate`` + ``S``), extended by the realized droplet-size ensemble (read
from each member's ``ion.npz`` radii; N = (R / 2.2173)³, the propagation-path
bulk convention, D0 §15.3).

Order of operations (oracle first; failure stops before a new number):

1. **Production anchor** — production s1–s5 re-pooled here equal the
   committed ``pooled_new`` row (N = 5000) of ``detfix_battery_table.csv`` on
   :data:`ANCHOR_COLS` to 4 decimals.
2. **Arm identity** — each A/B member's cfg differs from its production
   member's in exactly ``droplet_size_sampler_mode``, ``raw`` → ``post_pickup``.

Output (committed; every column, memory "scorer CSVs carry full vector"):
``data/runs/h2b_forward_model/detfix_postpickup_table.csv`` — one row per
member per arm, ``pooled_raw5000`` / ``pooled_pp5000``, and
``delta_pp_minus_raw``; and ``detfix_postpickup_chord.csv`` — the
chord -> (n, KE) map per arm (emission-axis chord bins, D0 §15.8).

Invocation::

    python scripts/post_processing/tier2_detfix_postpickup_table.py
"""

from __future__ import annotations

import dataclasses
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    pool_confirmation_reads,
    read_confirmation_detection,
)
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.gen_tier2_detfix_battery import (  # noqa: E402
    member_run_dir_name as production_run_dir_name,
)
from scripts.gen_tier2_detfix_postpickup import (  # noqa: E402
    POOL_MEMBERS,
    member_run_dir_name as postpickup_run_dir_name,
)
from scripts.post_processing.tier2_detfix_battery_table import (  # noqa: E402
    OracleFailure,
    TABLE_CSV as DETFIX_TABLE_CSV,
    _check_4dp,
    _csv_row,
    _row,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    RUNS_ROOT,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

OUT_CSV = RUNS_ROOT / "h2b_forward_model" / "detfix_postpickup_table.csv"
CHORD_CSV = RUNS_ROOT / "h2b_forward_model" / "detfix_postpickup_chord.csv"
#: Emission-axis chord bins [Å] for the chord -> (n, KE) map (D0 §15.8).
CHORD_EDGES_ANGSTROM = (0, 25, 30, 33, 35, 37, 39, 41, 44, 48, 55, 65, 80, 200)
ANCHOR_COLS = ("nbar_det", "n1_solv", "w1_solv", "midHot", "deepKE",
               "chi2_med", "S")
#: Propagation-path N -> R convention (bulk density, D0 §15.3).
R_PER_N_CUBE_ROOT = 2.2173

# ---------------------------------------------------------------------------


def size_columns(run_dirs: list[Path]) -> dict[str, float]:
    """Realized droplet-size ensemble of the given runs' molecules.

    Reads ``droplet_radii_angstrom`` (one radius per molecule) from each
    ``ion.npz``. Returns mean / median N [He atoms] and mean / q05 / q95 R [Å].
    """
    radii = []
    for run_dir in run_dirs:
        with np.load(run_dir / "ion.npz", allow_pickle=False) as z:
            radii.append(np.asarray(z["droplet_radii_angstrom"], dtype=float))
    R = np.concatenate(radii)
    N = (R / R_PER_N_CUBE_ROOT) ** 3
    q05, q95 = np.quantile(R, [0.05, 0.95])
    return {"N_mean": float(N.mean()), "N_median": float(np.median(N)),
            "R_mean": float(R.mean()), "R_q05": float(q05), "R_q95": float(q95)}


def emission_chords(run_dir: Path) -> np.ndarray:
    """Straight-line path [Å] from each fragment's birth point to the droplet
    surface along the Coulomb-explosion axis (partner -> fragment).

    Uses the first stored ion-stage column; fragment ``i`` pairs with
    ``i ± num_molecules``. Chord = -p·ê + sqrt((p·ê)² - |p|² + R²) for birth
    position p, unit axis ê, droplet radius R. Ignores deflection and the
    soft erf edge of the density gate. Shape (2 * num_molecules,).
    """
    with np.load(run_dir / "ion.npz", allow_pickle=False) as z:
        R = np.asarray(z["droplet_radii_angstrom"], dtype=float)
        P = np.stack([z[k][:, 0] for k in ("positions_x", "positions_y",
                                           "positions_z")], axis=1)
    n_mol = R.size // 2
    partner = np.r_[np.arange(n_mol, 2 * n_mol), np.arange(n_mol)]
    axis = P - P[partner]
    axis /= np.linalg.norm(axis, axis=1, keepdims=True)
    pe = np.sum(P * axis, axis=1)
    return -pe + np.sqrt(pe ** 2 - np.sum(P * P, axis=1) + R ** 2)


def chord_rows(arm: str, run_dirs: list[Path]) -> list[dict[str, Any]]:
    """Chord-binned outcome map of one arm (pooled members): per bin the ion
    count, retained fraction, and over scored ions the n = 1 fraction, mean n
    and mean detected KE [eV]."""
    chords, n_det, scored, ke = [], [], [], []
    for run_dir in run_dirs:
        det = load_detection_result(run_dir / "detection.npz")
        chords.append(emission_chords(run_dir))
        n_det.append(det.n_detected)
        scored.append(det.detected_mask)
        ke.append(det.E_kin_detected_eV)
    c, n, ok, k = map(np.concatenate, (chords, n_det, scored, ke))
    rows = []
    edges = CHORD_EDGES_ANGSTROM
    for lo, hi in zip(edges[:-1], edges[1:]):
        in_bin = (c >= lo) & (c < hi)
        s = in_bin & ok
        rows.append({
            "arm": arm, "chord_lo_A": lo, "chord_hi_A": hi,
            "ions": int(in_bin.sum()),
            "retained_frac": float(1 - ok[in_bin].mean()) if in_bin.any() else "",
            "n1_frac": float((n[s] == 1).mean()) if s.any() else "",
            "mean_n": float(n[s].mean()) if s.any() else "",
            "mean_KE_eV": float(np.nanmean(k[s])) if s.any() else "",
        })
    return rows


def _arm_member_row(label, member, run_dir, ab, ked) -> tuple[dict[str, Any], Any]:
    cfg = RunDirectory(run_dir).load_cfg()
    read = read_confirmation_detection(load_detection_result(run_dir / "detection.npz"),
                                       label=label)
    row = _row(label, member, cfg.droplet_size_sampler_mode, cfg.seed, read, ab, ked)
    row.update(size_columns([run_dir]))
    return row, read


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    arms = {"raw": production_run_dir_name, "pp": postpickup_run_dir_name}
    missing = [arms[a](m) for a in arms for m in POOL_MEMBERS
               if not (RUNS_ROOT / arms[a](m) / "detection.npz").exists()]
    if missing:
        print(f"not yet run: {missing} — stopping before any number.")
        return 1

    # Oracle 2: arm identity, member by member.
    for member in POOL_MEMBERS:
        cfg_raw = RunDirectory(RUNS_ROOT / arms["raw"](member)).load_cfg()
        cfg_pp = RunDirectory(RUNS_ROOT / arms["pp"](member)).load_cfg()
        diff = sorted(f.name for f in dataclasses.fields(cfg_raw)
                      if getattr(cfg_raw, f.name) != getattr(cfg_pp, f.name))
        if (diff != ["droplet_size_sampler_mode"]
                or (cfg_raw.droplet_size_sampler_mode, cfg_pp.droplet_size_sampler_mode)
                != ("raw", "post_pickup")):
            raise OracleFailure(f"[{member}] arm cfg diff {diff} / modes "
                                f"{cfg_raw.droplet_size_sampler_mode!r} → "
                                f"{cfg_pp.droplet_size_sampler_mode!r}")

    rows: list[dict[str, Any]] = []
    pooled: dict[str, dict[str, Any]] = {}
    for arm, name_fn in arms.items():
        reads = []
        for member in POOL_MEMBERS:
            row, read = _arm_member_row(f"{member}_{arm}", member,
                                        RUNS_ROOT / name_fn(member), ab, ked)
            rows.append(row)
            reads.append(read)
        mode = "raw" if arm == "raw" else "post_pickup"
        pooled[arm] = _row(f"pooled_{arm}5000", "pooled5000", mode, "-",
                           pool_confirmation_reads(reads, label=f"pooled_{arm}"),
                           ab, ked)
        pooled[arm].update(size_columns([RUNS_ROOT / name_fn(m) for m in POOL_MEMBERS]))
        if arm == "raw":
            # Oracle 1: the production reference reproduces before any new number.
            _check_4dp("production pooled_new (s1-s5) anchor", pooled[arm],
                       _csv_row(DETFIX_TABLE_CSV, "pooled_new"), ANCHOR_COLS)
            print("production anchor PASSED (pooled_new, s1-s5, 4 dp)", flush=True)
        rows.append(pooled[arm])

    delta: dict[str, Any] = {"label": "delta_pp_minus_raw", "member": "pooled5000",
                             "pipeline": "post_pickup - raw", "seed": "-"}
    for key, value in pooled["pp"].items():
        if key in delta:
            continue
        try:
            delta[key] = float(value) - float(pooled["raw"][key])
        except (TypeError, ValueError):
            delta[key] = ""
    rows.append(delta)
    keys = list(rows[0].keys())
    rows = [{k: r.get(k, "") for k in keys} for r in rows]

    show = ("label", "pipeline", "num_scored", "trap", "trap_bound", "trap_marg",
            "supp", "nbar_det", "n1_solv", "w1_solv", "midHot", "deepKE",
            "chi2_med", "KE1_mean", "KE2_mean", "S", "N_mean", "R_mean")
    show = tuple(c for c in show if c in keys)
    print("\n=== droplet-size sampler A/B at h405 — observable vector ===")
    print(format_table([{c: r[c] for c in show} for r in rows]))
    print(f"\nCSV -> {write_rows_csv(OUT_CSV, rows)}")

    chord = []
    for arm, name_fn in arms.items():
        chord += chord_rows(arm, [RUNS_ROOT / name_fn(m) for m in POOL_MEMBERS])
    print("\n=== chord -> (n, KE) map, s1–s5 pooled per arm ===")
    print(format_table(chord))
    print(f"\nCSV -> {write_rows_csv(CHORD_CSV, chord)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
