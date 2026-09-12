# Free-J_c fit h0*(J) = A (J - J_c)^p': SUMMARY

**VERDICT: J_c departs from 2 — drop (J-T) language; report the empirical J-dependence / fitted J_c as a measured quantity.**

- 2 of 8 comparable cells (n >= 5 J points; the n=4 cell xi=0.65/target=0.10 excluded) have 2 inside the 95% bootstrap CI of J_c.
- J_c point estimates (comparable cells) range from 2.155 to 2.435.

Method: profile grid search over J_c in [-1.0, 2.45] (step 0.005; J_c <= min(J) - 0.05); at each J_c, OLS of log(h0*) on log(J - J_c); J_c chosen to minimize SSE. Bootstrap: 1000 path-level resamples reusing the amplitude-requirement fit's stored per-path collapse matrices, same resample across the h0 axis within each (J, xi) file and shared across J within a xi block (seed 0xB007 + xi_idx, identical to h0_scaling.py); inversion of h0* mirrors h0_scaling.invert_h0_star exactly (imported). Replicates: unreachable targets drop that J from the replicate's fit; replicates with < 4 valid J points are dropped. 95% percentile CIs.

## Free-J_c fits (comparable cells, n = 6 J points)

| xi | target | n | A | p' | 95% CI p' | J_c | 95% CI J_c | R^2 | CI contains 2? |
|----|--------|---|---|----|-----------|-----|------------|-----|----|
| 0.35 | 0.1 | 6 | 1.896 | 0.634855 | [0.52397, 0.841878] | 2.17 | [1.83475, 2.305] | 0.999724 | YES |
| 0.35 | 0.2 | 6 | 1.3277 | 0.584859 | [0.502807, 0.728531] | 2.37 | [2.20987, 2.42512] | 0.998782 | NO |
| 0.35 | 0.3 | 6 | 0.805292 | 0.604782 | [0.530162, 0.785846] | 2.435 | [2.345, 2.45] | 0.999123 | NO |
| 0.5 | 0.1 | 6 | 3.29136 | 0.516508 | [0.405591, 0.639458] | 2.33 | [2.185, 2.42] | 0.995525 | NO |
| 0.5 | 0.2 | 6 | 1.86014 | 0.647201 | [0.520025, 0.802135] | 2.37 | [2.23, 2.45] | 0.998376 | NO |
| 0.5 | 0.3 | 6 | 1.01189 | 0.801296 | [0.663395, 0.955514] | 2.405 | [2.34, 2.45] | 0.998352 | NO |
| 0.65 | 0.2 | 6 | 2.08589 | 0.704075 | [0.565844, 0.958655] | 2.345 | [2.09987, 2.435] | 0.99681 | NO |
| 0.65 | 0.3 | 6 | 0.693475 | 1.23368 | [0.920659, 1.55536] | 2.155 | [1.91987, 2.33] | 0.98697 | YES |

## Flagged cell (n=4, not comparable — excluded from the count)

| xi | target | n | A | p' | 95% CI p' | J_c | 95% CI J_c | R^2 | CI contains 2? |
|----|--------|---|---|----|-----------|-----|------------|-----|----|
| 0.65 | 0.1 | 4 | 2.5685 | 0.960278 | [0.57483, 2.7785] | 1.79 | [-1, 2.24] | 0.992209 | YES |

## Pattern across cells

- xi=0.35: target 0.1: J_c=2.17, target 0.2: J_c=2.37, target 0.3: J_c=2.435
- xi=0.5: target 0.1: J_c=2.33, target 0.2: J_c=2.37, target 0.3: J_c=2.405
- xi=0.65: target 0.1: J_c=1.79 (n=4, flagged), target 0.2: J_c=2.345, target 0.3: J_c=2.155
- target=0.1: xi 0.35: J_c=2.17, xi 0.5: J_c=2.33, xi 0.65: J_c=1.79 (n=4, flagged)
- target=0.2: xi 0.35: J_c=2.37, xi 0.5: J_c=2.37, xi 0.65: J_c=2.345
- target=0.3: xi 0.35: J_c=2.435, xi 0.5: J_c=2.405, xi 0.65: J_c=2.155

## Bootstrap pathologies

| xi | target | reps used | reps dropped (<4 valid J) | reps with J_c at grid boundary | unreachable-target counts per J (of 1000 reps) |
|----|--------|-----------|----------------------------|-------------------------------|------------------------|
| 0.35 | 0.1 | 1000 | 0 | 0 | none |
| 0.35 | 0.2 | 1000 | 0 | 2 | none |
| 0.35 | 0.3 | 1000 | 0 | 261 | none |
| 0.5 | 0.1 | 1000 | 0 | 0 | J=5.0: 32 |
| 0.5 | 0.2 | 1000 | 0 | 31 | none |
| 0.5 | 0.3 | 1000 | 0 | 45 | none |
| 0.65 | 0.1 | 998 | 2 | 54 | J=4.0: 2; J=4.5: 986; J=5.0: 1000 |
| 0.65 | 0.2 | 1000 | 0 | 6 | none |
| 0.65 | 0.3 | 1000 | 0 | 0 | none |

No point-estimate J_c pinned at a grid boundary.
