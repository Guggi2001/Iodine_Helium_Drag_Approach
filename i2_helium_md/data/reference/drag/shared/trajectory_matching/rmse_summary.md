# Drag-form fit results — RMSE summary

One-glance comparison of every trajectory-matching fit (Method-B §9 + §10):
the shared-form joint fits and Stage-1 analogs below, plus the §10.8 per-case
single-curve diagnostic fits in the final section. Pulled directly from the
`fit_parameters.json` / `stage1_*.json` artifacts in this directory and in
`data/reference/drag/<case>/trajectory_matching/<variant>/`. **Generated
2026-06-14.**

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

## Per-case single-curve fits (§8 + §10.8 — DIAGNOSTIC, never preset-wired)

Each row is a fit to **one curve only** — the complete per-case diagonal across
all three realized forms: `linear_cubic` (§8) at
`data/reference/drag/<case>/trajectory_matching/fit_parameters.json`, and the
`linear_quadratic` / `power_law` variants (§10.8) at
`.../trajectory_matching/<variant>/fit_parameters.json`. These are
**calibrated-not-validated**: each fit trivially matches its own curve, so the
RMSE below is the *self-fit* on the very curve being fitted (**not** held-out, not
comparable to the Stage-1 prediction RMSEs above). **9 Å is the genuinely
non-radial / transverse-contaminated case.** The production law is the §9 *shared*
`shared_pure_cubic` joint fit — **no per-case row here is a production candidate**
(the §8 per-case `linear_cubic` fits motivated the shared refit and were then
superseded by it; the §10.8 fits are diagnostic only).

Half-width conventions differ by provenance: the `linear_cubic` (§8) bands are
**seed-sweep std + RMSE sensitivity**; the §10.8 form bands are
**sensitivity-only** (seed sweep omitted per §8). For `power_law`, `C_err` is the
fixed-`n` partial band `C_err = γ_ref_err / v_ref**(n−1)` (not a marginal
uncertainty).

| Case | Variant | Form | Parameters (± half-width) | E_bind | Self RMSE | Escape | Diagnostic note |
|---|---|---|---|---|---|---|---|
| 18 Å | `linear_cubic` *(§8)* | linear_cubic | a = 1.46 ± 1.46, b = 3.32 ± 0.16 | 0.0709 ± 0.0068 | 0.0968 | yes | `a_err = a` exactly → `a` consistent with 0: the §8 weak-`a` ridge that motivated the shared refit |
| 18 Å | `lq_shared_3param` *(§10.8)* | linear_quadratic | a = 16.63 ± 1.59, c = 6.13 ± 0.59 | 0.0532 ± 0.0071 | 0.0994 | yes | `a` *not* collapsed — but degenerate with the `a=0` corner (next row, equal RMSE): the §10.4.1 weak-`a` ridge |
| 18 Å | `lq_shared_pure_quadratic` *(§10.8)* | linear_quadratic | a ≡ 0, c = 11.27 ± 0.55 | 0.0598 ± 0.0080 | 0.0993 | yes | same RMSE as 3-param → `a` carries no in-curve information |
| 18 Å | `pl_shared_3param` *(§10.8)* | power_law | C = 1.05 ± 0.10, **n = 4.000 ± 1.333** | 0.0825 ± 0.0079 | 0.0929 | yes | **exponent unidentified** — railed to the n = 4 bound; 18 Å's narrow speed range has no exponent leverage |
| 9 Å | `linear_cubic` *(§8)* | linear_cubic | a = 7.50 ± 0.72, b = 1.75 ± 0.03 | 0.1543 ± 0.0031 | 0.0409 | yes | `a` resolved away from 0 here (unlike 18 Å) — but jointly the §9 shared fit still drives `a → 0` |
| 9 Å | `lq_shared_3param` *(§10.8)* | linear_quadratic | **a → 0**, c = 10.32 ± 0.10 | 0.1303 ± 0.0044 | 0.0446 | yes | **collapses to pure-quadratic** independently |
| 9 Å | `lq_shared_pure_quadratic` *(§10.8)* | linear_quadratic | a ≡ 0, c = 10.32 ± 0.10 | 0.1303 ± 0.0064 | 0.0446 | yes | identical to the 3-param fit (already at `a = 0`) |
| 9 Å | `pl_shared_3param` *(§10.8)* | power_law | C = 3.59 ± 0.04, **n = 2.651 ± 0.026** | 0.1566 ± 0.0031 | 0.0406 | yes | **exponent identified tightly** — 9 Å's broad speed range pins `n` within the single curve |

*The 18 Å rows are the same fits as the 18 Å-only Stage-1 analogs above (the §8
`linear_cubic` row = the "pure-cubic, reused §8 18 Å bundle" Stage-1 row; the two
§10.8 rows = the lq / pl Stage-1 analogs). The RMSE here is their 18 Å self-fit;
the Stage-1 table reports their held-out 9 Å prediction instead.*

**Reading the per-case table.** Per-case exponent leverage is **case-asymmetric**,
refining the §10.4.1 "exponent unidentified within a single curve" expectation:
9 Å alone pins **n̂ = 2.651 ± 0.026** and collapses `a → 0`, while 18 Å alone
**rails to n = 4** (half-width 1.333) and its `a` is the flat weak-`a` ridge. The
shared joint fit's **n̂ ≈ 2.927 ≈ 3** is the *cross-case* combination of the two
speed scales — neither single case lands at 3 on its own (9 Å biased low and
contaminated; 18 Å unidentified). The §8 `linear_cubic` rows show the same
per-case-vs-cross-case pattern in the *linear* term: 18 Å's `a` is consistent
with 0 (the original weak-`a` finding, `a_err = a`) while 9 Å's is resolved
(`a = 7.50`), yet the §9 shared fit collapses `a → 0` jointly. Diagnostic only;
the incumbent is unaffected. (Full record: METHOD_B §8 / §10.8;
`drag_migration_log.md`.)
