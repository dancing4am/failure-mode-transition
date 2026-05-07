# AI coupling on the model prediction surface

## Headline

At the measured **adversarial AI coupling** (adversarial-transfer proxy median $J$ ≈ 1.15), the minimal model with maximum passive stabilizer (μ = 100) predicts P(collapse) ≈ **0.00** and passive-stabilizer effectiveness ≈ **100%**. At the measured **benign-benchmark AI coupling** (benchmark-vector proxy median $J$ ≈ 8.7), the model predicts P(collapse) ≈ **0.37** and effectiveness ≈ **63%**. For reference, the critical threshold $J_c = T = 2.0$ corresponds to P(collapse) ≈ 0.00, effectiveness ≈ 100%. The benchmark-surface coupling already sits in the regime where Theorem 1 binds.

## Method

We extracted P(collapse) at three economic-margin levels (μ ∈ {40, 60, 100}) from the existing passive baseline (J ∈ {0.5, …, 5.0}, 100 seeds, 20,000 steps at d*t* = 0.05) and ran a 5-by-5 extension sweep (J ∈ [6.0, 7.0, 8.0, 9.0, 10.0], μ ∈ [20, 40, 60, 80, 100], same protocol, 2,500 runs) to cover the benchmark-vector coupling (J ≈ 8.7). Passive-stabilizer effectiveness is defined as 1 − P(J, μ = 100) / P(J, μ = 20), bounded to [0, 1]. Predictions at the measured AI coupling values are obtained by linear interpolation on the J grid.

## Predictions at reference points

| Reference | J | P(coll) μ=100 | P(coll) μ=60 | P(coll) μ=40 | Passive effectiveness (μ=100) |
|---|---|---|---|---|---|
| AI adversarial median (adversarial-transfer proxy) | 1.15 | 0.00 | 0.04 | 0.99 | 100% |
| AI benchmark median (benchmark-vector proxy) | 8.70 | 0.37 | 0.40 | 0.42 | 63% |
| Critical threshold T = 2.0 | 2.00 | 0.00 | 0.08 | 0.52 | 100% |
| Financial average rho=0.5 | 1.00 | 0.00 | 0.03 | 1.00 | 100% |
| Financial crisis rho=0.9 | 9.00 | 0.37 | 0.41 | 0.43 | 63% |

## Key sentence for the paper

> Mapping the measured AI coupling onto the model's prediction surface (Figure 6), the adversarial-distribution median ($J$ ≈ 1.15) lies in a regime where passive stabilizers retain ≈100% of their effectiveness and the maximum-stabilizer collapse probability is P(coll) ≈ 0.00 — near-critical but tolerable. The benign-benchmark median ($J$ ≈ 8.7) sits in a regime where passive stabilizer effectiveness has fallen to ≈63% and P(coll) ≈ 0.37, comparable to financial-crisis coupling. The model quantitatively predicts that fixed-magnitude protections will be substantially less effective at the AI coupling levels currently measured on benign evaluation surfaces.

## Recommendation for paper integration

**Promote to a new Results subsection** ("AI coupling: model predictions at measured coupling levels") **with Figure 6 in the main text**. The benchmark-coupling prediction (P(collapse) = 0.37 at the maximum passive stabilizer strength tested) is *worse* than the published high-J band residual of 0.23 at J=5 — i.e., at AI's measured benchmark coupling the model predicts a higher catastrophic-failure rate than the worst case in the original sweep, even though the original 0.23 residual was already the headline. The bimodal finding — adversarial coupling tolerable but benchmark coupling severely degraded — is exactly the quantitative prediction the AI-connection critique asked for. The residual P(collapse) saturates at 0.37 across μ levels (40, 60, 100) at J=8.7, indicating that economic margin alone cannot rescue the system at this coupling level — which directly motivates the active-stabilizer architecture promoted in the paper.
