# Phase 2 pre-registration — discriminating test of the model-implied tail-risk shape

> **Editorial amendment (2026-09-09).** The original frozen document is retained in the
> private working repository at git blob `c415a6f1024de8a0a074c52aba19f4f801f1da8a`
> (SHA-256 of the original file content:
> `b56ec870dfdc002e82081e9405b117836aa5c3486f2a2ecf3801b4af56025f97`). Two lines (original
> lines 12 and 67) were reworded on 2026-08-24 to remove internal phrasing. No element
> of the data specification, candidate models, metrics, or decision rule was changed.
> The original is available from the author on request.

**Status: PRE-REGISTERED. Written and committed BEFORE the fitting script
(`phase2_discriminating_test.py`) was created or run.** The success criterion is an
HONEST verdict, not a positive one. No specification search, no dataset search: the data,
predictor, outcome, candidate models, metrics, and decision rule below are fixed in advance,
and the result is reported whatever it is (including "cannot distinguish" or "model shape
disfavoured").

## Question
The existing empirical evidence is a window-level rank correlation ρ(ρ̄, tail risk) ≈ 0.5–0.65.
A central methodological critique: this cannot separate the **specific** model-implied dependence
of tail-risk on coupling from a **generic monotone** "more coupling → more tail risk" null.
This test asks: **does the S&P 500 tail-risk data prefer the model-implied functional SHAPE
of P(tail | J) over flexible monotone alternatives, or can it not distinguish them?**

## Data (fixed)
- Panel: the committed S&P 500 sector-ETF daily-return panel (`close_prices.csv`), 9 sectors
  pre-2015-10-08 + XLRE from 2015, 2004–2025, exactly as built by `tail_risk_frequency.py`
  (`window_metrics`, 60-day rolling windows).
- Per window: ρ̄ = mean off-diagonal pairwise correlation; **predictor** x = J = ρ̄/(1−ρ̄);
  **outcome** y = 1[max drawdown in window > 5%] (primary). Secondary, reported but not
  decision-bearing: thresholds 3% and 7%.
- Overlap is acknowledged: ~5,414 raw windows but effective ~47–90 independent (integrated
  autocorrelation 60–115, per §S10). Therefore the **primary discriminator is out-of-sample
  (held-out era), not in-sample p-values**, and all in-sample AIC/LR numbers are reported with
  the explicit effective-N caveat (raw-N likelihoods overstate evidence).

## Candidate models (fixed; all binomial/Bernoulli for y on the same windows)
- **M0 — null:** intercept only (constant P).
- **M1 — MODEL-IMPLIED SHAPE:** the predictor is the model's own collapse curve
  P_model(J) = P(collapse | J, μ=100) from
  `simulation/results/ai_coupling_overlay/combined_sweep_summary.csv` (the Suzuki-type
  saturating shape). Fit `logit(P_emp) = a + b·logit(P_model(J))` — a 2-parameter affine
  recalibration **on the logit of the model's own shape**, so M1 is forced to follow the
  model's specific J-curvature. (Equivalent analytical form, reported as a cross-check:
  probit on the Suzuki selection index s(J) = −(c/ξ)√(2/λ(J)), c=tanh(h₀/T),
  λ(J)=−1+(J/T)sech²(h₀/T), committed constants T=2, ξ=0.5, h₀=0.4.)
- **M2 — generic logistic in J:** `logit(P) = a + b·J`. 2 params. Generic monotone S-curve,
  NOT the model's curvature.
- **M3 — generic logistic in log J (≈ in ρ̄):** `logit(P) = a + b·log J`. 2 params.
- **M4 — isotonic (nonparametric monotone):** PAV monotone fit of P(y) on J. Used as the
  **upper bound** on monotone log-likelihood (how much any monotone shape could gain);
  effective df reported via the number of distinct levels.
- **M5 — flexible monotone spline:** logistic with a natural cubic spline basis in J (4 df).
  Generic flexible shape; M2/M3 are (approximately) nested in it.

## Metrics (fixed)
1. **In-sample:** maximized log-likelihood, AIC, BIC for M0–M3, M5; isotonic log-likelihood
   for M4. (Raw-N; reported with effective-N caveat.)
2. **Nested LR test:** M2 and M3 vs M5 (spline) — is a generic 1-parameter logistic shape
   rejected in favour of a flexible monotone shape? (Tests whether *any* simple monotone form
   is missing real curvature, independent of the model.)
3. **Out-of-sample (PRIMARY):** era split — train on pre-2015-10-08 windows, test on
   2015–2025 windows; also the reverse split. Report mean held-out Bernoulli log-loss (nats)
   for M1, M2, M3, M5. Lower is better. Each model's parameters are estimated on train only.

## DECISION RULE (fixed, decided in advance)
Let ΔAIC = AIC(M1) − min(AIC(M2), AIC(M3)) and ΔOOS = mean-OOS-logloss(M1) −
min(mean-OOS-logloss over {M2, M3}) (averaged over the two era splits).

- **DISCRIMINATING CONFIRMATION (model shape genuinely preferred):**
  ΔAIC ≤ −2 (M1 lower by ≥2) **AND** ΔOOS ≤ −0.002 nats (M1 better OOS) **AND** M1 is not
  rejected against M5 (spline) by the nested LR at the effective-N-adjusted level.
- **CANNOT DISTINGUISH (data lack resolution to prefer the specific shape):**
  |ΔAIC| < 2 **AND** |ΔOOS| < 0.002 nats — M1 and the generic monotone alternatives fit
  comparably. (This converts that critique into a precise, stated limitation.)
- **MODEL SHAPE DISFAVOURED:** min(AIC(M2),AIC(M3)) or AIC(M5) lower than AIC(M1) by ≥2
  **AND** the corresponding model better OOS by ≥0.002 nats — the data prefer a shape the
  model does not predict; the specific feature (e.g. curvature/saturation) is reported.
- If the three criteria conflict (e.g. AIC and OOS disagree), the **OOS verdict is decisive**
  (it best respects the window overlap), and the conflict is reported explicitly.

## Integrity commitments
- Run once, with the seed and code committed. Report all of M0–M5, both era splits, and all
  three thresholds — no dropping of unfavourable rows.
- A null ("cannot distinguish") or a "disfavoured" result is a fully acceptable, reportable
  outcome and will be reported plainly; it strengthens the honest-limitation framing and is
  NOT to be rescued by trying alternative predictors, links, or window definitions.
- The prior finding that the model's predicted *convexity* is not supported
  (`convexity_test.py`) is noted in advance; this test does not re-open or re-tune that.
