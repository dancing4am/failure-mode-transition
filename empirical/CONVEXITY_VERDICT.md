# Convexity / nonlinearity test of P(tail | J)

## Headline (NEGATIVE — but informative)

The empirical tail-risk curve P(max DD > 5% | J) over 10 coupling deciles is **strongly nonlinear, but concave** (saturating toward the P=1 ceiling) rather than convex (accelerating). Quadratic fit gives J² coefficient *c* = **−0.057** (95% bootstrap CI [−0.204, −0.034]); AIC(linear) − AIC(quadratic) = +11.46 (quadratic strongly preferred over linear, and adj R² rises 0.73 → 0.92). Piecewise slope: low-half slope = 0.47, high-half slope = 0.10 — the slope **flattens** by 5× at high coupling, not steepens. Power-law log-log slope = 1.45 (minimal-model log-log slope ≈ 2.51 over the same J range). Mean finite-difference second derivative = −0.61.

The model's *predicted* P(collapse | J) curve is convex over the tested J range because P(collapse) is still rising toward its high-J asymptote and has not saturated. The empirical P(tail event) is bounded by 1.0 by construction and saturates earlier, so its second derivative becomes negative once P(tail) approaches the ceiling. The convexity prediction is not a useful defense for the exponent mismatch — the empirical curve is concave on the same J axis. This is the honest negative result.

## Model fits

| model | parameters | RSS | adj R² | AIC |
|---|---|---|---|---|
| linear | a=0.186, b=0.185 | 0.142 | 0.733 | -36.56 |
| quadratic | c·J² + b·J + a, c=-0.057 | 0.037 | 0.919 | -48.02 |
| power law | log y = 1.45 log J + -1.58 (n=10) | — | — | — |

## Method

Three regressions on the 10 decile-binned data points (J_bin_mid, P_tail_5pct) from the tail-risk analysis: linear (P = a + bJ), quadratic (P = a + bJ + cJ²), and power-law (log P = a + b log J on positive points). AIC uses the Gaussian-likelihood approximation; quadratic is preferred over linear if AIC(linear) - AIC(quadratic) > 2. Piecewise slopes are computed by linear regression on the lower-half (5 deciles) and upper-half (5 deciles) separately. Bootstrap CI on the quadratic coefficient *c* uses 1000 resamples of the 10 deciles.

## Key sentence for the paper

> The empirical tail-risk curve P(DD > 5% | J) is strongly nonlinear (quadratic fit improves AIC by 11.5 over linear; adj R² rises from 0.73 to 0.92), but concave rather than convex: the high-J half-slope is 0.20 of the low-J half-slope, indicating saturation as P(tail) approaches the [0, 1] ceiling. The model's predicted convex P(collapse) curve has not saturated in the tested J range. The exponent mismatch (model log-log slope 2.51 vs. empirical 1.45) reflects this saturation — not a failure of monotonicity or coupling-dependence, both of which are robust.

## Recommendation

**Do NOT integrate as a convexity-defense for the exponent mismatch** — the empirical curve is concave on the tested J axis, so the prediction "the model captures accelerating failure that linear models miss" is not supported. **Do consider integrating as an honest characterization** of the exponent mismatch: the empirical curve saturates because P_tail is ceiling-bounded, while the model is still in its rising regime; this is consistent with the model's *asymptotic* h/J² scaling but not with its intermediate-coupling acceleration. One sentence in the Limitations exponent paragraph would suffice.
