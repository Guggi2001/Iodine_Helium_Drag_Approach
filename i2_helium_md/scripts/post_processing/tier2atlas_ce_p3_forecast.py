"""Atlas (C) pre-step P3 — weight-composition forecast table (zero MD).

The design-doc §9 P3: the mixture arithmetic over the measured §3.5k budget
rows, computed from committed artifacts only and frozen as the registered
placement expectation for the (C) probe
(`TIER2_CE_CHANNEL_EXIT_STRIP_DESIGN.md`; findings "(C) pre-steps P1 + P2").

Inputs (ALL committed — the gate-on-committed-artifacts rule):

- ``atlas_budget_probe.csv``   — the §3.5k measured rows (h405 / bud226k /
  bud411k), the per-channel proxies.
- ``atlas_ce_strip_p1.csv``    — the P1 grid; its analytic Poisson-binomial
  head (P(0), P(1), P(2)) at the pinned box is the suppressed-conversion
  kernel.
- P2 measured widths (sigma_Q2 0.311, sigma_Q3 0.553 eV; committed findings
  values, re-derivable by ``tier2atlas_ce_gas_widths.py``).

Conventions (each stated once here, carried into every table):

- **Channel proxies**: Q2 <- bud226k (budget 2.26 vs channel 2.16), Q3 <-
  bud411k (4.11 vs 4.32). KE1 is budget-corrected along the measured
  kinematic slope S_k = 0.389 eV/eV; n1/supp/trap are NOT corrected (their
  budget slopes are unmeasured; the 4-5 % proxy bias is a stated risk).
- **Single channel** (E_single 0.53 eV, P2): no measured row exists. The
  §3.5k line extrapolates KE1 < 0 at this budget — these ions sit below the
  strip velocity gate and produce no n = 1 at reference KE. Registered
  assumption: n1_single in [0, 0.05] (band ends carried through), supp = 0,
  conversions = 0.
- **Conversion kernel**: suppressed fraction s_c converts through the P1
  analytic head at the twin speed (P0 = 1 above v_strip, so the kernel
  depends on (j0, w_j) only): P(0) -> bare, P(1) -> n = 1, P(2) -> n = 2,
  remainder -> small-n solvated (>= 3). Evaluated at the four pinned-box
  corners (j0 in {1.5, 2}, w_j in {0.5, 1.0}).
- **Converted KE**: KE_conv,c = KE1_c * (0.774 / 0.637) — the h405-measured
  converted/existing n = 1 ratio (P1 KE1_conv vs the committed KE1); bare
  conversions carry KE_bare,c = KE_conv,c * m_I/m(1). Approximation, stated:
  it assumes the conv/existing speed ratio is budget-independent.
- **Widths**: per-channel n = 1 KE dispersion sigma_c(KE1) =
  sqrt((S_k * sigma_c)^2 + needle^2), needle = 0.04 (the §3.5k measured
  per-cell KE1 SD); components treated as Gaussians for the blend moments,
  the mode scan, and the band weights.
- **Strip toll NOT applied** (P1/§3.5j limitation carried): forecast KE1 is
  toll-high by up to ~0.2 eV at full strip; the MD probe measures it.

Oracles (hard-fail before any new number):

- O1: the budget-CSV proxy rows reproduce the committed §3.5k findings
  values (n1/KE1/supp/trap to 3 decimals).
- O2: the P1-CSV kernel cells reproduce this script's own analytic
  Poisson-binomial head (self-consistency of the committed artifact).

Pure arithmetic — runs nothing, mutates nothing.

Invocation:

    python scripts/post_processing/tier2atlas_ce_p3_forecast.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Any

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.physics.constants import MASS_HE_AMU, MASS_I_ION_AMU  # noqa: E402
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS — the frozen anchors (design §3.5 table, P1/P2-updated)
# ---------------------------------------------------------------------------

RUNS_ROOT = PROJECT_ROOT / "data" / "runs"
BUDGET_CSV = RUNS_ROOT / "h2b_forward_model" / "atlas_budget_probe.csv"
P1_CSV = RUNS_ROOT / "h2b_forward_model" / "atlas_ce_strip_p1.csv"
SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_ce_p3_forecast.csv"

# OQ-C weight vector (adjudicated procedure; OQ-K closed 2026-07-29).
WEIGHTS = {"single": 0.30, "Q2": 0.50, "Q3": 0.20}

# Channel source energies [eV per I+] (design §3.1, f = 0.8 paper anchor).
E_CHANNEL = {"Q2": 2.16, "Q3": 4.32}
PROXY_ROW = {"Q2": "bud226k", "Q3": "bud411k"}

S_K = 0.389                     # measured KER->KE1 slope [eV/eV] (§3.5k)
SIGMA_KER = {"Q2": 0.311, "Q3": 0.553}   # P2 measured [eV per I+]
NEEDLE_SD = 0.04                # §3.5k per-cell KE1 SD scale
KE_CONV_RATIO = 0.774 / 0.637   # P1/h405 converted-vs-existing n = 1 KE
BARE_MASS_RATIO = MASS_I_ION_AMU / (MASS_I_ION_AMU + MASS_HE_AMU)

N1_SINGLE_BAND = (0.0, 0.05)    # registered single-channel assumption

# The pinned P1 box corners (a = 2; kernel is a-independent above v_strip).
KERNEL_CELLS = [(2.0, 1.5, 0.5), (2.0, 1.5, 1.0),
                (2.0, 2.0, 0.5), (2.0, 2.0, 1.0)]

# O1 oracle: committed §3.5k findings values for the proxy rows.
O1_ROWS = {
    "h405":    {"n1_solv": 0.208, "KE1_mean": 0.637, "supp": 0.183,
                "trap": 0.087},
    "bud226k": {"n1_solv": 0.126, "KE1_mean": 0.486, "supp": 0.073,
                "trap": 0.179},
    "bud411k": {"n1_solv": 0.329, "KE1_mean": 1.200, "supp": 0.546,
                "trap": 0.001},
}
O1_TOL = 0.0015

# Experimental n = 1 KED reference row (committed IHe_KED_reference.csv,
# §3.5j Read 1) + the bare-row band reads (§3.5l read (i)) — display targets.
REF_N1_KED = {"mean": 1.302, "mode": 0.891, "sigma": 0.697,
              "above_1p15": 0.487}
REF_N1_SOLV = 0.31
REF_BARE_BELOW_1EV = 0.019

# ---------------------------------------------------------------------------


def _read_csv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _logistic_head(j0: float, wj: float, n_rungs: int = 21,
                   kmax: int = 2) -> np.ndarray:
    """P(exactly k survivors), k = 0..kmax, at P0 = 1 (the P1 kernel)."""
    j = np.arange(1, n_rungs + 1)
    surv = 1.0 - 1.0 / (1.0 + np.exp(-(j - j0) / wj))
    pmf = np.zeros(kmax + 1)
    pmf[0] = 1.0
    for s in surv:
        carry = pmf[:-1] * s
        pmf *= 1.0 - s
        pmf[1:] += carry
    return pmf


def _blend_stats(comps: list[tuple[float, float, float]]) -> dict[str, float]:
    """Weight/mu/sigma Gaussian components -> blend moments, mode, band weight."""
    w = np.array([c[0] for c in comps])
    mu = np.array([c[1] for c in comps])
    sg = np.array([c[2] for c in comps])
    if w.sum() <= 0:
        return {"weight": 0.0, "mean": float("nan"), "sd": float("nan"),
                "mode": float("nan"), "above_1p15": float("nan")}
    p = w / w.sum()
    mean = float((p * mu).sum())
    var = float((p * (sg ** 2 + mu ** 2)).sum() - mean ** 2)
    grid = np.linspace(0.0, 4.0, 4001)
    dens = np.zeros_like(grid)
    for pi, mi, si in zip(p, mu, sg):
        dens += pi / si * np.exp(-0.5 * ((grid - mi) / si) ** 2)
    from math import erf, sqrt
    above = float(sum(
        pi * 0.5 * (1.0 - erf((1.15 - mi) / (si * sqrt(2.0))))
        for pi, mi, si in zip(p, mu, sg)
    ))
    return {"weight": float(w.sum()), "mean": mean,
            "sd": float(np.sqrt(max(var, 0.0))),
            "mode": float(grid[int(np.argmax(dens))]), "above_1p15": above}


def main() -> None:
    budget_rows = {r["label"]: r for r in _read_csv(BUDGET_CSV)}
    p1_rows = {(float(r["a"]), float(r["j0"]), float(r["w_j"])): r
               for r in _read_csv(P1_CSV)}

    # ---- O1: proxy rows reproduce the committed §3.5k values.
    bad = []
    for label, expect in O1_ROWS.items():
        row = budget_rows[label]
        for col, ref in expect.items():
            got = float(row[col])
            if abs(got - ref) > O1_TOL:
                bad.append(f"{label}.{col}: {got:.4f} vs committed {ref:.4f}")
    if bad:
        print("*** O1 PROXY-ROW ORACLE FAILED — do not read the forecast ***")
        for line in bad:
            print(f"    {line}")
        return
    print("O1 proxy-row oracle OK (§3.5k committed values reproduced).")

    # ---- O2: the P1 kernel cells match this script's analytic head.
    bad = []
    for cell in KERNEL_CELLS:
        row = p1_rows[cell]
        pmf = _logistic_head(cell[1], cell[2])
        for k, col in enumerate(("P_strip0", "P_surv1", "P_surv2")):
            if abs(float(row[col]) - pmf[k]) > 1e-9:
                bad.append(f"{cell} {col}: csv {float(row[col]):.6f} vs "
                           f"analytic {pmf[k]:.6f}")
    if bad:
        print("*** O2 KERNEL ORACLE FAILED — do not read the forecast ***")
        for line in bad:
            print(f"    {line}")
        return
    print("O2 kernel oracle OK (P1 CSV head = analytic Poisson-binomial).\n")

    # ---- Per-channel input table.
    chan: dict[str, dict[str, float]] = {}
    for c in ("Q2", "Q3"):
        row = budget_rows[PROXY_ROW[c]]
        ke1 = float(row["KE1_mean"]) + S_K * (E_CHANNEL[c]
                                              - float(row["budget_eV"]))
        chan[c] = {
            "n1": float(row["n1_solv"]), "supp": float(row["supp"]),
            "trap": float(row["trap"]), "nbar": float(row["nbar_det"]),
            "KE1": ke1,
            "sd": float(np.hypot(S_K * SIGMA_KER[c], NEEDLE_SD)),
        }
    print("=== Per-channel inputs (proxy rows, S_k-corrected KE1) ===")
    for c in ("Q2", "Q3"):
        d = chan[c]
        print(f"  {c} (w {WEIGHTS[c]:.2f}, E {E_CHANNEL[c]:.2f} eV <- "
              f"{PROXY_ROW[c]}): n1 {d['n1']:.3f}, KE1 {d['KE1']:.3f} "
              f"+/- {d['sd']:.3f}, supp {d['supp']:.3f}, trap {d['trap']:.3f}")
    print(f"  single (w {WEIGHTS['single']:.2f}, E 0.53 eV): no measured row "
          f"— n1 in {N1_SINGLE_BAND}, supp 0, no conversions (registered "
          "assumption).\n")

    # ---- B-only (strip-off): the §3.5k complementarity forecast.
    rows_out: list[dict[str, Any]] = []
    for n1s in N1_SINGLE_BAND:
        comps = [(WEIGHTS["Q2"] * chan["Q2"]["n1"], chan["Q2"]["KE1"],
                  chan["Q2"]["sd"]),
                 (WEIGHTS["Q3"] * chan["Q3"]["n1"], chan["Q3"]["KE1"],
                  chan["Q3"]["sd"])]
        if n1s > 0:
            comps.append((WEIGHTS["single"] * n1s, 0.10, 0.05))
        st = _blend_stats(comps)
        supp_mix = sum(WEIGHTS[c] * chan[c]["supp"] for c in ("Q2", "Q3"))
        rows_out.append({
            "variant": "B-only", "kernel": "-", "n1_single": n1s,
            "n1": st["weight"], "KE1_mean": st["mean"], "KE1_sd": st["sd"],
            "KE1_mode": st["mode"], "above_1p15": st["above_1p15"],
            "supp": supp_mix, "bare": 0.0, "slow_bare": 0.0,
            "n2_conv_extra": 0.0,
            "q3_share_n1": (WEIGHTS["Q3"] * chan["Q3"]["n1"]) / st["weight"]
            if st["weight"] else float("nan"),
        })

    # ---- C-full (strip-on) at the four kernel corners x single-band ends.
    for cell in KERNEL_CELLS:
        pmf = _logistic_head(cell[1], cell[2])
        for n1s in N1_SINGLE_BAND:
            comps = []
            bare_w = 0.0
            slow_bare_w = 0.0
            n1_tot = 0.0
            n2_extra = 0.0
            for c in ("Q2", "Q3"):
                d = chan[c]
                w = WEIGHTS[c]
                ke_conv = d["KE1"] * KE_CONV_RATIO
                sd_conv = d["sd"] * KE_CONV_RATIO
                comps.append((w * d["n1"], d["KE1"], d["sd"]))
                comps.append((w * d["supp"] * pmf[1], ke_conv, sd_conv))
                n1_tot += w * (d["n1"] + d["supp"] * pmf[1])
                n2_extra += w * d["supp"] * pmf[2]
                kb = ke_conv * BARE_MASS_RATIO
                bw = w * d["supp"] * pmf[0]
                bare_w += bw
                if kb < 1.0:
                    slow_bare_w += bw
            if n1s > 0:
                comps.append((WEIGHTS["single"] * n1s, 0.10, 0.05))
                n1_tot += WEIGHTS["single"] * n1s
            st = _blend_stats(comps)
            q3_n1 = (WEIGHTS["Q3"] * (chan["Q3"]["n1"]
                                      + chan["Q3"]["supp"] * pmf[1]))
            rows_out.append({
                "variant": "C-full",
                "kernel": f"j0 {cell[1]:g} w {cell[2]:g}",
                "n1_single": n1s,
                "n1": n1_tot, "KE1_mean": st["mean"], "KE1_sd": st["sd"],
                "KE1_mode": st["mode"], "above_1p15": st["above_1p15"],
                "supp": 0.0, "bare": bare_w, "slow_bare": slow_bare_w,
                "n2_conv_extra": n2_extra,
                "q3_share_n1": q3_n1 / n1_tot if n1_tot else float("nan"),
            })

    print("=== P3 forecast (weights single/Q2/Q3 = "
          f"{WEIGHTS['single']}/{WEIGHTS['Q2']}/{WEIGHTS['Q3']}) ===")
    print(f"{'variant':>7} {'kernel':>12} {'n1_s':>5} | {'n1':>6} "
          f"{'KE1':>6} {'sd':>6} {'mode':>6} {'>1.15':>6} | {'supp':>5} "
          f"{'bare':>6} {'slowB':>6} {'Q3@n1':>6}")
    for r in rows_out:
        print(f"{r['variant']:>7} {r['kernel']:>12} {r['n1_single']:>5.2f} | "
              f"{r['n1']:>6.3f} {r['KE1_mean']:>6.3f} {r['KE1_sd']:>6.3f} "
              f"{r['KE1_mode']:>6.3f} {r['above_1p15']:>6.3f} | "
              f"{r['supp']:>5.3f} {r['bare']:>6.3f} {r['slow_bare']:>6.3f} "
              f"{r['q3_share_n1']:>6.3f}")
    print(f"\nreference n = 1 KED: mean {REF_N1_KED['mean']}, mode "
          f"{REF_N1_KED['mode']}, sigma {REF_N1_KED['sigma']}, above-1.15 "
          f"{REF_N1_KED['above_1p15']}; reference n1_solv {REF_N1_SOLV}; "
          f"experimental bare < 1 eV: {REF_BARE_BELOW_1EV}.")
    print("Caveats carried: toll not applied (KE forecast toll-high); "
          "proxy budgets 4-5 % off the channel means (KE1 corrected via "
          "S_k, occupancies not); KE_conv budget-independence assumed; "
          "single-channel row assumed, not measured.")

    if SAVE_CSV_PATH is not None:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows_out)}")


if __name__ == "__main__":
    main()
