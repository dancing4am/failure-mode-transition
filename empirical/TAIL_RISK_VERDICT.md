# Tail-risk frequency vs coupling

## Headline (MODERATE)

The empirical tail-risk frequency curve is correlated with the model's predicted P(collapse | J) curve (Spearman = 0.81 at 5% drawdown), but the correspondence is partial. The directional finding is robust; the quantitative shape match is partial.

## Per-bin empirical results

| J bin mid | ρ̄ mid | n windows | P(DD>3%) | P(DD>5%) | P(DD>7%) | P(DD>10%) | model P(coll) μ=100 |
|---|---|---|---|---|---|---|---|
| 0.46 | 0.32 | 542 | 0.18 | 0.01 | 0.01 | 0.00 | 0.00 |
| 0.71 | 0.42 | 541 | 0.69 | 0.31 | 0.19 | 0.07 | 0.00 |
| 0.90 | 0.47 | 541 | 0.78 | 0.32 | 0.19 | 0.02 | 0.00 |
| 1.14 | 0.53 | 542 | 0.77 | 0.38 | 0.20 | 0.02 | 0.00 |
| 1.43 | 0.59 | 541 | 0.92 | 0.53 | 0.24 | 0.07 | 0.00 |
| 1.70 | 0.63 | 541 | 0.90 | 0.67 | 0.36 | 0.16 | 0.00 |
| 2.00 | 0.67 | 542 | 0.92 | 0.64 | 0.49 | 0.28 | 0.00 |
| 2.37 | 0.70 | 541 | 0.94 | 0.71 | 0.34 | 0.18 | 0.04 |
| 2.98 | 0.75 | 541 | 0.93 | 0.79 | 0.49 | 0.28 | 0.14 |
| 4.90 | 0.83 | 542 | 1.00 | 0.94 | 0.86 | 0.61 | 0.27 |

## Spearman rank correlation per threshold

| threshold | Spearman(model, empirical) | p | log-log slope (empirical) | log-log slope (model) |
|---|---|---|---|---|
| 3% | +0.798 | 0.00568 | +0.53 | +2.51 |
| 5% | +0.813 | 0.00426 | +1.45 | +2.51 |
| 7% | +0.693 | 0.0262 | +1.48 | +2.51 |
| 10% | +0.753 | 0.0119 | +1.99 | +2.51 |

## Method

5,414 sixty-day rolling windows of daily log-returns for SPY +
9 sector ETFs (10 from 2015 with XLRE), 2004–2025. For each
window: (i) ρ̄ = mean off-diagonal pairwise correlation; 
(ii) equal-weight portfolio daily return; (iii) max drawdown
in the window (peak-to-trough). Tail event = max drawdown
exceeds the threshold (3 / 5 / 7 / 10 %). Map ρ̄ → J = ρ̄/(1-ρ̄).
Bin into 10 J deciles; compute empirical P(tail event |
J bin). Compare to the minimal-model P(collapse | J, μ = 100)
curve (extended sweep, J ∈ [0.5, 10], 100 seeds, dt = 0.05).

## Recommendation

Add to Results as a moderate-strength mechanism test, with honest discussion of where the shapes match and where they diverge.

## Key sentence for the paper

> Across 5,414 sixty-day rolling windows of S&P 500 sector data binned into 10 coupling deciles, the empirical frequency of tail events (max drawdown > 5%) tracks the minimal-model P(collapse | J) curve (Spearman ρ = +0.81). This tests the model's mechanism — branch-selection probability under coupling — rather than the variance-algebra identity that the diversification-benefit metric satisfies by construction.
