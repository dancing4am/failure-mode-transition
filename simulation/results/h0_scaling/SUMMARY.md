# Stabilizer-scaling exponent h0*(J): SUMMARY

Model: minimal Curie-Weiss + wealth (local generalized engine, bit-identical to `minimal_model.run_cell` at defaults; see regression check in h0_scaling.py). mult=100, n_seeds=2000 per (J, xi) with matched seeds across the h0 axis (seed base 0xB05CA1E + j_idx*1_000_003 + xi_idx*9973). Bootstrap: 1000 path-level resamples, 95% percentile CIs.

**Scope note:** the log(h0*) = a + p log(J) fit is restricted to the six supercritical couplings J in {2.5, ..., 5.0}, i.e. J > T = 2, because the Curie-Weiss bifurcation at J = T sits inside the manuscript's sweep; a pure power law in J cannot hold across the critical point. The secondary descriptor log(h0*) = a' + p' log(J - 2) uses the distance to that bifurcation.

## Exponent fits

| xi | target | p (log J) | 95% CI | R^2 | p' (log(J-2)) | 95% CI | R^2 |
|----|--------|-----------|--------|-----|---------------|--------|-----|
| 0.35 | 0.1 | 1.90388 | [1.73693, 2.05384] | 0.957543 | 0.758647 | [0.688805, 0.818898] | 0.99734 |
| 0.35 | 0.2 | 2.38874 | [2.12216, 2.58929] | 0.909429 | 0.96621 | [0.851062, 1.05124] | 0.976025 |
| 0.35 | 0.3 | 2.94662 | [2.66612, 3.26695] | 0.862503 | 1.20765 | [1.09038, 1.33967] | 0.950346 |
| 0.5 | 0.1 | 1.94557 | [1.77905, 2.12953] | 0.919439 | 0.784244 | [0.715405, 0.854071] | 0.979981 |
| 0.5 | 0.2 | 2.63729 | [2.35729, 2.93706] | 0.904893 | 1.06843 | [0.952707, 1.19918] | 0.974229 |
| 0.5 | 0.3 | 3.5463 | [3.33424, 3.74634] | 0.881062 | 1.44677 | [1.36345, 1.52234] | 0.96193 |
| 0.65 | 0.1 | 2.31755 | [1.89265, 2.64547] | 0.981788 | 0.783978 | [0.633775, 0.893935] | 0.990898 |
| 0.65 | 0.2 | 2.71959 | [2.423, 3.04392] | 0.909919 | 1.10063 | [0.970639, 1.23701] | 0.97761 |
| 0.65 | 0.3 | 3.62474 | [3.32895, 3.86763] | 0.940975 | 1.44779 | [1.33582, 1.54393] | 0.984735 |

## Verdict

- Measured exponent p ranges from 1.904 to 3.625 across the 9 (target, xi) combinations.
- p = 1 lies OUTSIDE the 95% CI for every combination: proportional scaling h0* ~ J is NOT supported; the measured exponent is as tabled.
- Stability: spread of p across xi in {0.35, 0.5, 0.65} and targets {0.10, 0.20, 0.30} is 1.721 (noticeable variation; see table).

## Monotonicity of P(collapse) in h0

No violations beyond Monte-Carlo noise detected (paired matched-seed McNemar-style check, z > 2 threshold). P(collapse) is monotone non-increasing in h0 up to noise.

## Unreachable targets (recorded NA, no extrapolation)

- J=4.5, xi=0.65, target=0.1: unreachable (P range [0.115,0.991] on h0 in [0.00,6.00])
- J=5.0, xi=0.65, target=0.1: unreachable (P range [0.130,0.991] on h0 in [0.00,6.00])
