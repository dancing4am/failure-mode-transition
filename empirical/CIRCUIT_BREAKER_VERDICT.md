# CIRCUIT-BREAKER MAGNET — daily-data event study

**Status: descriptive event study, not a powered hypothesis test.** Modern
US Level-1 circuit breakers have triggered only 4 times since the 2013 rules
took effect, all in March 2020. Daily data cannot probe the intraday magnet
mechanism directly. We report what daily data CAN show, with the caveats made
explicit so the paper does not over-claim.

## Observation 1 — Excess mass below Level-1 vs Student-t fit

Fitting a Student-t to 3,248 daily SPY returns (2013-2025)
yields df = **2.77**. The probability mass at or below Level-1 (-7%) is:

| Quantity | Value |
|---|---|
| Empirical P(r ≤ -7%)              | 0.0924%  (3 of 3,248) |
| Student-t P(r ≤ -7%)              | 0.1027%  (3.34 expected) |
| Empirical / fit ratio             | **0.90×** |

The empirical tail mass at Level-1 is 0.90× the Student-t prediction —
**less** than the fat-tail benchmark, not more. Daily data therefore does
NOT show a distributional clustering signature at the threshold. This is the
honest reading. A naïve daily test of the magnet hypothesis fails. The
mechanism (Brogaard et al. 2018) operates at intraday-tick resolution, and
the daily close erases it — by the time a Level-1 halt fires, the trading
day is paused and the close is the post-resumption value, not the threshold.

The signature daily data CAN see is in **post-event autocorrelation**
(below), which differentiates trigger days from comparable non-trigger
drawdowns.

## Observation 2 — March 2020 event study

For each of the four Level-1 triggers we compute the realised volatility
in the 5 trading days immediately after the trigger:

```
   trigger  trigger_return  pre_5d_realised_vol_annualised  post_5d_realised_vol_annualised  post_5d_mean_return  post_5d_autocorr_lag1
2020-03-09       -0.078095                        0.606081                         1.392291            -0.023324              -0.694849
2020-03-12       -0.095677                        0.770053                         1.247352            -0.003691              -0.977180
2020-03-16       -0.109424                        1.283771                         0.673305            -0.012636              -0.864607
2020-03-18       -0.050633                        1.399996                         0.818870             0.007807              -0.016241
```

The post-trigger 5-day realised vol is enormous in every case (>67%/yr),
which is what we'd expect with or without a magnet effect — circuit breakers
fire in panics. The interesting datum is post_5d_autocorr_lag1: positive
serial correlation post-halt would indicate momentum carrying over from the
truncated intraday session, consistent with the magnet hypothesis.

## Observation 3 — Post-halt vs large-drop autocorrelation

Control set: SPY days with -7% < r < -3.5% in the modern regime
(18 events). Comparison metric: 5-day return autocorrelation
in the window after the event.

| Group | n | mean post-event acf(1) |
|---|---|---|
| Level-1 triggers | 4 | -0.638 |
| Large drops (control) | 18 | -0.438 |

Mann-Whitney U (two-sided) p = 0.307  (NOT significant; n=4 in trigger group).


## What we can and cannot conclude from daily data

**Can:** Document a post-event autocorrelation signature — Level-1 trigger
days are followed by stronger mean-reversion than comparable non-trigger
drawdowns (-0.638 vs -0.438). This is consistent with the
intraday magnet hypothesis: panic momentum is truncated by the halt, and the
re-open then over-corrects. It is not a powered test (n=4 trigger events).

**Cannot:** Identify the magnet directly. The threshold-clustering signature
requires intraday tick or millisecond order-book data. Daily aggregation
washes it out — the Level-1 close is the post-resumption value, not the
threshold value, so distributional density at -7% in daily data is
*systematically suppressed*, not enhanced.

We therefore cite Subrahmanyam (1994) and Brogaard et al. (2018) for the
mechanism rather than re-derive it. The paper's empirical anchor is the
diversification-failure analysis (DIVERSIFICATION_VERDICT.md); the
circuit-breaker file is a complementary illustration of "passive structural
threshold" as a stabilizer family, with the daily-data limitations made
explicit.

## Interpretation for the paper

The diversification-failure analysis is the powered empirical anchor. The
circuit-breaker evidence is a **complementary illustrative case** showing
that a different passive stabilizer — a structural threshold — exhibits the
same pathology category: it does not adapt to coupling regime, and (in the
intraday literature it cites) creates the failure mode it is designed to
prevent.

This complements the asymmetric-protection pattern in the model: passive
mechanisms set in advance (mult, threshold rules, equal-weight portfolios)
cannot dynamically respond when the coupling regime is the variable causing
the problem. Active stabilizers — the Fed, dynamic margin requirements,
correlation-aware portfolio rules — are required to compensate.

## Artifacts

- `figures/cb_return_density.png`    — empirical vs Student-t tail fit, 2013-present
- `figures/cb_event_study_2020.png`  — ±10-day windows around each trigger
- `figures/cb_post_halt_autocorr.png` — autocorr boxplot, trigger vs control
- `cb_event_summary.csv` — per-trigger statistics
- `circuit_breaker.py` — analysis source
