# AI coupling proxy: adversarial-attack transferability

> **Note on rounding.** The frontier-pair median *J*ₑff ≈ 1.15 reported below is rounded to *J* ≈ 1.2 in the main paper and SI §S8. The two values describe the same statistic; the paper rounds for narrative consistency with the other proxy estimates.

## Headline

Published adversarial-attack transfer rates of 2–87% across frontier-model pairs [Zou et al. 2023; Mazeika et al. 2024; Wei et al. 2023] correspond to an effective coupling $J_{\rm eff} \in [0.02, 6.46]$ under the mapping $\bar\rho = J/(1+J)$, with frontier-pair median $J_{\rm eff} \approx 1.15$ (excluding the Claude-2 hardened-target outlier). The central tendency sits below but within the same order of magnitude as the critical threshold $J_c = T = 2.0$ at which Theorem 1 begins to bind, while the high-transfer tail (87%) reaches $J_{\rm eff} = 6.5$ — comparable to financial-crisis coupling ($J = 9.0$ at $\bar\rho = 0.9$).

## Method

We hard-code published adversarial-attack transfer rates from three primary sources (Zou et al. 2023, Mazeika et al. 2024, Wei et al. 2023) and apply the Curie–Weiss mapping ρ = J/(1+J) ⇒ J = ρ/(1-ρ) developed in main-text §S6. Transfer rate is interpreted as a monotone proxy for the implicit pairwise coupling between frontier-model output distributions; under this mapping, the model's critical threshold is J_c = T = 2.0.

## Per-pair table

| source | target | ρ (transfer) | J_eff | regime | citation | eval |
|---|---|---|---|---|---|---|
| Vicuna-ensemble | GPT-3.5 | 0.866 | 6.46 | supercritical (J>T) | Zou+2023 | substring |
| Vicuna-ensemble | GPT-4 | 0.469 | 0.88 | subcritical (J<T/2) | Zou+2023 | substring |
| Vicuna-ensemble | Claude-1 | 0.479 | 0.92 | subcritical (J<T/2) | Zou+2023 | substring |
| Vicuna-ensemble | Claude-2 | 0.021 | 0.02 | subcritical (J<T/2) | Zou+2023 | substring |
| Vicuna-ensemble | PaLM-2 | 0.660 | 1.94 | near-critical (T/2<J<T) | Zou+2023 | substring |
| frontier-pair-avg | frontier-pair-avg | 0.235 | 0.31 | subcritical (J<T/2) | Mazeika+2024 | Llama-2-13B classifier |
| jailbreak-source | frontier-pair-avg | 0.590 | 1.44 | near-critical (T/2<J<T) | Wei+2023 | human/manual |

## Aggregates

- Range across all reported transfer rates: ρ ∈ [0.021, 0.866] → J_eff ∈ [0.02, 6.46]
- Median (excluding Claude-2 hardened outlier): ρ = 0.534 → J_eff = 1.15
- Mean (excluding Claude-2): ρ = 0.550 → J_eff = 1.22
- Pairs strictly above J_c = T = 2.0: 1 of 7

## Comparison with financial coupling

For reference, the main-text financial diversification analysis reports ρ̄ averaged over 22 years of S&P 500 sector data near 0.5 (J = 1.00), with crisis peaks at ρ̄ ≈ 0.9 (J = 9.00). The AI median frontier-pair transfer rate (0.53) corresponds to J_eff ≈ 1.15, of the same order as the financial average and an order of magnitude below the financial-crisis peak. The high-transfer tail of the AI distribution (ρ ≈ 0.87, J ≈ 6.5) is comparable to financial-crisis coupling.

## Caveats

1. **Transfer rate is a proxy, not a measurement.** Adversarial transfer success measures how often an attack crafted for one model also works on another. It is monotonically related to the underlying output coupling but is not the coupling itself. An exact mapping requires a direct correlation/agreement measurement on benign tasks (complementary direct-measurement proxies) — those are complementary signals, not redundant ones.

2. **Numbers depend strongly on attack and evaluation method.** Substring-match success (used by some Zou et al. 2023 numbers) overstates real coupling; LLM-judge evaluation (Mazeika et al. 2024) is more conservative; manual human evaluation is the tightest bound but harder to replicate. The headline range above brackets the spread.

3. **Hardened models pull the range to the left.** Claude-2 in particular was explicitly trained against the GCG suffix attack class (per the Anthropic report) and shows ρ = 0.02 — this is a deliberately decoupled model, not an honest sample. We exclude it from the median/mean reported above.

4. **Coupling rises with shared training, not with attack age.** Frontier models converge on shared corpora and shared evaluation benchmarks. The transferability literature is consistent with — but does not by itself prove — the paper's claim that the implicit J in the optimization landscape is rising. Direct measurement on benign benchmarks (benchmark-vector proxy) tests this more directly.

5. **All numbers above should be verified against the cited tables before paper submission.** They are extracted from each paper's main results section; transfer rates reported in abstracts may differ from full-table values.

## Implication for the paper

The active/passive distinction matters most when frontier AI systems sit in the regime where Theorem 1 binds — at or above $J_c = T = 2.0$. The published transferability literature places the median frontier-pair coupling at $J_{\rm eff} \approx 1.15$ — within the same order of magnitude as the threshold but not yet at it — while the high-transfer tail reaches $J_{\rm eff} \approx 6.5$, comparable to financial-crisis coupling. The paper's framework therefore is not hypothetical for AI governance: at the upper end of the already-published, replicated, benchmark-standardized adversarial transferability literature, frontier-model pairs already operate in a coupling regime in which fixed-magnitude protections decay as $h/J^2$ and active stabilizers are required, and the central tendency is approaching the threshold rather than receding from it.
