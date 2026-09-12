# Reduction validity — matched σ_eff comparison (hand-compiled from the CSVs in this directory)

Script: `simulation/scripts/reduction_validity.py`. Data: `reduction_validity.csv`.
Cell: J = 5.0, μ = 100. σ_eff = sqrt(σ_ind²/N + σ_com²); total noise budget σ_ind² + σ_com² = 0.25.
Scalar runs use `minimal_model.run_cell` at ξ = σ_eff. Bootstrap 95% CIs (path resampling).

## Part 1 — idiosyncratic endpoint (r = 0), matched ξ_eff = 0.5/√N

| N | ABM P(collapse) [95% CI] (n=50) | Scalar at ξ_eff [95% CI] (n=1000) |
|---|---|---|
| 1 | 0.24 [0.12, 0.36] | 0.287 [0.259, 0.316] |
| 5 | 0.10 [0.02, 0.20] | 0.122 [0.102, 0.143] |
| 25 | 0.02 [0.00, 0.06] | 0.005 [0.001, 0.010] |
| 200 | 0.00 [0.00, 0.00] | 0.000 |

## Part 2 — mixed noise at N = 200, r = σ_com²/0.25

| r | ABM P(collapse) [95% CI] (n=100) | Scalar at matched ξ_eff [95% CI] (n=1000) |
|---|---|---|
| 0.00 | 0.00 | 0.000 |
| 0.25 | 0.17 [0.10, 0.24] | 0.152 [0.130, 0.175] |
| 0.50 | 0.18 [0.11, 0.26] | 0.212 [0.187, 0.238] |
| 0.75 | 0.32 [0.23, 0.41] | 0.257 [0.230, 0.284] |
| 1.00 | 0.23 [0.15, 0.31] | 0.314 [0.285, 0.344] |

## Conclusion

The one-dimensional reduction with σ_eff = sqrt(σ_ind²/N + σ_com²) reproduces the
ABM collapse statistics across the full range of noise compositions, within sampling
error (ABM n = 50–100, binomial SE ≈ 0.045–0.065). This documents the validity range
of the reduction. It is not evidence that N or the network does work: σ_eff is an
identity stating that independent noise averages out.
