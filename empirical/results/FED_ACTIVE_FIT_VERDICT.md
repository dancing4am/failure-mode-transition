# Fed active-stabilizer interaction test

**Events:** 12 FOMC cuts >= 25 bp, 2004-2025

## Headline statistics

- Mean cut size: 0.29 percentage points
- Mean J at time of cut: 0.489
- Mean 90-day post-cut SPY return: +0.058

## Regression: ret_90d ~ J × |Δr|

Active-stabilizer form predicts: cut effect *grows* with coupling.
Test the interaction term's sign and significance.

| term | coefficient | std err | t | p |
|---|---|---|---|---|
| const | -0.0202 | 0.1010 | -0.20 | 0.847 |
| J_proxy | +0.3943 | 0.1866 | +2.11 | 0.068 |
| cut_size_abs | -0.7874 | 0.4554 | -1.73 | 0.122 |
| interaction | +0.6582 | 0.5898 | +1.12 | 0.297 |

*R² = 0.713, adj R² = 0.605, N = 12*

## Interpretation

- The interaction term is positive (coefficient = +0.6582) but not significant at p < 0.10 (p = 0.297). The point estimate is consistent with the active-stabilizer hypothesis but the small sample (N = 12) limits statistical power.
