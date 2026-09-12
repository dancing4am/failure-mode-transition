# Horizon stability check — minimal model

Cells: J in (2.5, 3.5, 5.0), mult = 100; dt = 0.05; long horizon = 80000 steps (t = 4000).

## Part 1 — exact reproduction check (n_seeds = 100, standard seeds)

Same cell seed, n_steps = 80,000; outcomes restricted to collapse_step < 20,000 must exactly match the archive.

- J = 2.5, mult = 100 (seed = 4057458): archived n_collapsed = 5, reproduced (step < 20,000) = 5; per-seed collapsed flags match: True; per-seed collapse_step match: True -> EXACT MATCH
- J = 3.5, mult = 100 (seed = 6057464): archived n_collapsed = 21, reproduced (step < 20,000) = 21; per-seed collapsed flags match: True; per-seed collapse_step match: True -> EXACT MATCH
- J = 5.0, mult = 100 (seed = 9057473): archived n_collapsed = 28, reproduced (step < 20,000) = 28; per-seed collapsed flags match: True; per-seed collapse_step match: True -> EXACT MATCH

**Part 1 verdict: EXACT per-seed reproduction of the archive.**

## Part 2 — high-power horizon scan (n_seeds = 4000, seed = standard + 7777777)

### J = 2.5, mult = 100 (seed = 11835235)

| horizon t | n_collapsed | P(collapse) | 95% CI |
|---|---|---|---|
| 250 | 227 | 0.0568 | [0.0498, 0.0640] |
| 500 | 227 | 0.0568 | [0.0495, 0.0640] |
| 1000 | 227 | 0.0568 | [0.0500, 0.0640] |
| 2000 | 227 | 0.0568 | [0.0498, 0.0640] |
| 4000 | 227 | 0.0568 | [0.0498, 0.0640] |

### J = 3.5, mult = 100 (seed = 13835241)

| horizon t | n_collapsed | P(collapse) | 95% CI |
|---|---|---|---|
| 250 | 839 | 0.2097 | [0.1973, 0.2228] |
| 500 | 839 | 0.2097 | [0.1975, 0.2225] |
| 1000 | 839 | 0.2097 | [0.1975, 0.2223] |
| 2000 | 839 | 0.2097 | [0.1970, 0.2225] |
| 4000 | 839 | 0.2097 | [0.1968, 0.2225] |

### J = 5.0, mult = 100 (seed = 16835250)

| horizon t | n_collapsed | P(collapse) | 95% CI |
|---|---|---|---|
| 250 | 1178 | 0.2945 | [0.2802, 0.3088] |
| 500 | 1178 | 0.2945 | [0.2800, 0.3088] |
| 1000 | 1178 | 0.2945 | [0.2802, 0.3085] |
| 2000 | 1178 | 0.2945 | [0.2800, 0.3090] |
| 4000 | 1178 | 0.2945 | [0.2805, 0.3090] |

### Drift verdict

- J = 2.5: P@1000 = 0.0568 (CI [0.0500, 0.0640]); P@2000 = 0.0568; P@4000 = 0.0568 -> no drift
- J = 3.5: P@1000 = 0.2097 (CI [0.1975, 0.2223]); P@2000 = 0.2097; P@4000 = 0.2097 -> no drift
- J = 5.0: P@1000 = 0.2945 (CI [0.2802, 0.3085]); P@2000 = 0.2945; P@4000 = 0.2945 -> no drift

**VERDICT: NO DRIFT — P(collapse) at t = 2000 and t = 4000 lies inside the t = 1000 bootstrap CI for every cell.**

## Part 3 — collapse-time distribution (n = 4000, full t = 4000 run)

| J | n_collapsed | median t | p90 t | p99 t | max t | max/250 | max/500 | max/1000 | max/2000 | max/4000 | after 250 | after 500 | after 1000 | after 2000 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2.5 | 227 | 16.5 | 23.3 | 29.8 | 34.6 | 0.138 | 0.069 | 0.035 | 0.017 | 0.009 | 0 | 0 | 0 | 0 |
| 3.5 | 839 | 15.3 | 19.5 | 26.3 | 31.3 | 0.125 | 0.063 | 0.031 | 0.016 | 0.008 | 0 | 0 | 0 | 0 |
| 5.0 | 1178 | 14.6 | 17.0 | 20.7 | 28.8 | 0.115 | 0.058 | 0.029 | 0.014 | 0.007 | 0 | 0 | 0 | 0 |

Max collapse time saturates in ABSOLUTE units: the latest collapse across all three cells occurs at t = 34.6, and no collapses occur after t = 2000 even though the horizon extends to t = 4000. The max does NOT track the horizon.

## Part 4 — context: where do very late collapses live in the archive?

Archived rows with collapse_step >= 18,000 (of 3573 collapses total):

| J | mult | seed_idx | collapse_step | t (physical) |
|---|---|---|---|---|
| 5.0 | 40 | 15 | 18870 | 943.5 |

Latest archived collapse overall: step 18870 (t = 943.5, 94.3% of the t = 1000 horizon) in cell J = 5.0, mult = 40.
Latest archived collapse restricted to mult = 100 cells: step 671 (t = 33.6).

