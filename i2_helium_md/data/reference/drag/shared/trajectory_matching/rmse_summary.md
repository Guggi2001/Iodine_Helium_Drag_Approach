# Drag-form fit results — RMSE summary

One-glance comparison of every shared-form trajectory-matching fit (Method-B
§9 + §10), pulled directly from the `fit_parameters.json` / `stage1_*.json`
artifacts in this directory. **Generated 2026-06-14.**

- All RMSE and objective values are in **Å/ps**; the objective is the
  equal-weight mean of the two per-case in-window `|v2|` RMSEs (plus an escape
  penalty, which is 0 for every row below — all fits fully escape).
- **Escape = yes** means `escape_fraction = 1.0` (all 2N ions escaped the
  droplet); no fit below trapped.
- Parameter units: `a` [amu/ps], `b` [amu·ps/Å²], `c` [amu/Å],
  `C` [amu·Å^(1−n)·ps^(n−2)], `n` dimensionless, `E_bind` [eV].
- Gating bands (§9.4, reused by §10): 18 Å RMSE ≤ 0.19, 9 Å RMSE ≤ 0.45,
  escape = 1.0 (both cases).
- Incumbent reference objective (`shared_pure_cubic`): **0.1292649398514104**.

## Stage-2 — shared joint fits (both cases in the objective; gating)

| Variant (try) | Form | Parameters | E_bind | 18 Å RMSE | 9 Å RMSE | Escape (18 Å / 9 Å) | Objective (mean) | §9.4 bands | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| `shared_pure_cubic` *(§9, PRODUCTION)* | linear_cubic | a = 0, b = 2.5154 | 0.1168 | 0.1344 | 0.1241 | yes / yes | **0.129265** | PASS | incumbent / production candidate |
| `shared_3param` *(§9)* | linear_cubic | a = 0.00025, b = 2.5159 | 0.1168 | 0.1345 | 0.1240 | yes / yes | 0.129265 | PASS | `a → 0`; `T_a0`-equivalent to incumbent |
| `pl_shared_3param` *(§10)* | power_law | C = 2.8352, n = 2.9269 | 0.1126 | 0.1305 | 0.1278 | yes / yes | **0.129171** | PASS | **EQUIVALENT** (Δ = −0.0001); n̂ ≈ 3 confirms pure-cubic |
| `lq_shared_3param` *(§10)* | linear_quadratic | a = 0.0001, c = 12.7922 | 0.0482 | 0.1308 | 0.1955 | yes / yes | 0.163152 | PASS | **WORSE** (Δ = +0.0339); rejected-by-objective |
| `lq_shared_pure_quadratic` *(§10)* | linear_quadratic | a = 0, c = 12.7973 | 0.0483 | 0.1313 | 0.1950 | yes / yes | 0.163150 | PASS | **WORSE** (Δ = +0.0339); rejected-by-objective |

`Δ = variant objective − incumbent objective`. Every variant passes the §9.4
bands; the **objective** (not the bands) is what separates the forms.

## Stage-1 analog — 18 Å-only fit, predict held-out 9 Å (recorded, NON-gating)

Each family is fit on 18 Å alone, then the untouched 9 Å trajectory is scored
as a prediction (band ≤ 0.45 Å/ps, for comparability with pure-cubic's 0.2685).

| Family | Fit parameters (18 Å only) | E_bind | 9 Å predicted RMSE | Escape | Band ≤ 0.45 |
|---|---|---|---|---|---|
| pure-cubic *(§9, reused §8 18 Å bundle)* | a = 1.456, b = 3.316 | 0.0709 | **0.2685** | yes | PASS |
| `power_law` *(§10)* | C = 1.0529, n = 4.000 | 0.0825 | 0.4110 | yes | PASS |
| `linear_quadratic` *(§10)* | a = 16.6254, c = 6.1313 | 0.0532 | 0.6987 | yes | **FAIL** |

## Reading the table

- **`power_law` (free n) recovers the pure-cubic exponent**: the joint fit
  lands at **n̂ = 2.927 ≈ 3** with the *lowest* objective of all (0.129171,
  marginally under the incumbent but inside the ±0.005 equivalence band), and
  the held-out 9 Å prediction (0.411) passes. This is **not** the Method-A
  power-law exponent (n ≈ 2.06).
- **Forced `v²` (`linear_quadratic`) is measurably worse**: both lq variants
  collapse to the pure-quadratic corner (`a → 0`) and sit at objective 0.1632
  (Δ = +0.034, well outside the equivalence band), and the lq 9 Å prediction
  (0.699) fails its band.
- **Conclusion:** the trajectory objective *does* discriminate the drag
  exponent and it picks **n = 3**. No alternative form beats `shared_pure_cubic`
  → incumbent confirmed, no escalation, presets unchanged. (Full record:
  METHOD_B §10.7; `drag_migration_log.md`.)
