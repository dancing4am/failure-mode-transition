# Eq (2) sensitivity: SUMMARY

Setup: h0 = 2.0, xi = 0.5, mult = 100, J in {2.5, 3.0, 3.5, 4.0, 4.5, 5.0}, n_seeds = 2000, seeds and bootstrap resamples identical across variants (paired comparisons). Exponent re-measured at target P = 0.3.

## P(collapse | J) at h0 = 2.0 by variant

| variant | J=2.5 | J=3.0 | J=3.5 | J=4.0 | J=4.5 | J=5.0 |
|---|---|---|---|---|---|---|
| baseline | 0.058 | 0.141 | 0.202 | 0.253 | 0.280 | 0.308 |
| const-30 | 0.018 | 0.093 | 0.165 | 0.224 | 0.258 | 0.292 |
| const+30 | 0.100 | 0.184 | 0.231 | 0.274 | 0.299 | 0.324 |
| coeff-30 | 0.057 | 0.139 | 0.200 | 0.251 | 0.278 | 0.308 |
| coeff+30 | 0.061 | 0.142 | 0.202 | 0.253 | 0.280 | 0.309 |
| linear | 0.000 | 0.000 | 0.000 | 0.000 | 0.004 | 0.086 |

Maximum |delta P| vs baseline across all variants and J: **0.275** (linear at J=4.5 (delta = -0.275)).

## Exponent p (target 0.30, xi = 0.5) by variant

| variant | p | 95% CI | R^2 | delta p vs baseline (paired 95% CI) |
|---------|---|--------|-----|--------------------------------------|
| baseline | 3.5463 | [3.33424, 3.74634] | 0.881062 | 0 [0, 0] |
| const-30 | 5.18862 | [4.98078, 5.39312] | 0.878347 | 1.64233 [1.52571, 1.77891] |
| const+30 | 2.54586 | [2.28189, 2.76651] | 0.909682 | -1.00044 [-1.1264, -0.915063] |
| coeff-30 | 3.76853 | [3.48894, 4.07311] | 0.867326 | 0.222232 [0.126464, 0.381382] |
| coeff+30 | 3.42747 | [3.1419, 3.60143] | 0.88295 | -0.118824 [-0.226436, -0.0955765] |
| linear | NA | [NA, NA] | NA | NA [NA, NA] |

Range of p across the six variants: **[2.546, 5.189]** (spread 2.643).

## Range over which conclusions hold

- +/-30% perturbations of the consumption constant and coefficient, and a structurally different purely-linear consumption law, shift P(collapse) at the manuscript's operating point (h0=2, xi=0.5) by at most 0.275 in absolute probability.
- The stabilizer-scaling exponent stays within [2.546, 5.189] across all six consumption specifications.

## Monotonicity / NA notes

- No monotonicity violations beyond Monte-Carlo noise.
- NA: linear, J=2.5, target=0.3 unreachable (P range [0.000,0.000])
- NA: linear, J=3.0, target=0.3 unreachable (P range [0.000,0.000])
- NA: linear, J=3.5, target=0.3 unreachable (P range [0.000,0.000])
- NA: linear, J=4.0, target=0.3 unreachable (P range [0.000,0.000])
- NA: linear, J=4.5, target=0.3 unreachable (P range [0.000,0.072])
