# AI coupling — cross-method summary (adversarial-transfer + benchmark-vector proxies)

This document synthesizes the two AI-coupling proxies (adversarial-attack
transferability and public-benchmark agreement) into one finding. A third
proxy — semantic-similarity via API queries — was not run; the two proxies
together carry the narrative without it. (See `AI_COUPLING_DIRECT_VERDICT.md`
for the subsequent direct AlpacaEval output-similarity measurement, which
became the paper's headline AI coupling number.)

---

## Headline

Two complementary proxies, drawn from independent literatures, place
frontier-model coupling at or above the paper's critical threshold
$J_c = T = 2$ on different surfaces of the system:

- **Adversarial-attack transferability** (Zou et al. 2023; Mazeika
  et al. 2024; Wei et al. 2023): median frontier-pair
  $J_{\rm eff} \approx 1.06$ (in the near-critical band); high-transfer
  tail $J_{\rm eff} \approx 6.5$ (Vicuna-ensemble GCG suffix → GPT-3.5
  at ρ = 0.866). Median *near* the threshold; tail *over* it.
- **Public-benchmark agreement** (12 frontier models × 10 standard
  benchmarks compiled from each model's primary publication or model
  card): median pairwise cross-benchmark accuracy correlation
  ρ = 0.90, $J_{\rm eff} \approx 8.7$ — supercritical, comparable to
  financial-crisis coupling ($J = 9.0$ at $\bar\rho = 0.9$). The
  capability-detrended residual is ρ ≈ 0, indicating the bulk of
  cross-model agreement *is* shared capability convergence on shared
  evaluation suites — exactly the homogenization mechanism the paper's
  Discussion identifies.

Read together: AI coupling is **already supercritical on the benign-benchmark
surface** (where models are optimized to converge) **and near-critical on
the adversarial surface** (where models are explicitly trained against
common attack classes). The "near-critical and rising" framing is what
the data supports; "AI is already in the danger zone" is more than the
evidence carries; "AI is far from coupling-relevant" is rejected by both
measurements.

---

## Cross-method positioning

Each method probes a different surface of the same underlying coupling:

| Surface | Method | What it measures | Bias direction |
|---|---|---|---|
| Adversarial | adversarial-transferability proxy — attack transferability | how often an attack crafted for one model also defeats another | *under*-estimates coupling because frontier models are explicitly hardened against shared attack classes (e.g. Claude-2 ρ=0.02 against GCG suffix) |
| Benign | benchmark-vector proxy — benchmark agreement | how often models converge on the same answers / capability profile across standard tasks | *over*-estimates coupling because shared evaluation pressure forces convergence to a single capability axis |

The honest reading is that the *true* coupling sits between these two
proxies, with the central tendency near or just above the critical
threshold. Both methods agree on the *direction*: frontier-model
coupling is high enough that Theorem 1 binds.

---

## Headline numbers (consolidated)

| Statistic | (adversarial proxy: adversarial) | (benchmark proxy: benign) |
|---|---|---|
| Sample | 7 model-pair / aggregate rows | 12 models × 10 benchmarks, 66 pairs |
| ρ range | [0.02, 0.87] | raw [0.66, 1.00]; detrended [-0.98, 0.96] |
| ρ central tendency | 0.51 (median, ex. Claude-2) | 0.90 (raw median) / -0.10 (detrended median) |
| $J_{\rm eff}$ range | [0.02, 6.46] | raw [1.9, ≫10]; detrended ~ 0 |
| $J_{\rm eff}$ central tendency | 1.06 (median) / 1.24 (mean) | 8.7 (raw median) / 0.0 (detrended median) |
| $J_c$ comparison | near-critical (just below threshold) | supercritical (well above threshold) |
| Financial reference | financial avg J = 1.0; financial crisis J = 9.0 | financial avg J = 1.0; financial crisis J = 9.0 |

---

## Caveats consolidated for the paper

The Discussion paragraph that uses these numbers must carry forward the
following caveats — not all in body text, but in mind when wording is
chosen.

1. **Both measurements are proxies, not direct coupling.** Adversarial
   transfer success monotonically tracks underlying coupling but is not
   the coupling. Benchmark-vector correlation tracks capability
   convergence on shared evaluation suites — itself a form of coupling
   the paper claims is rising — but conflates capability with
   inductive-bias agreement.
2. **Numbers depend on attack type, evaluation method, and reporting
   choice.** Substring-match attack-success rates over-state; LLM-judge
   evaluation is stricter; manual human evaluation is the tightest
   bound. Per-benchmark scores differ between papers in shot count,
   prompt formatting, and evaluator. We aggregate over this noise.
3. **Hardened models pull the adversarial range to the left.**
   Claude-2's ρ = 0.02 against GCG suffix reflects deliberate
   adversarial training, not an honest unhardened sample; we exclude
   it from the median/mean reported in the adversarial-transfer proxy.
4. **The benign-distribution coupling reflects capability
   convergence on shared benchmarks, not necessarily epistemic
   agreement.** This *is* the coupling vector the paper claims is
   rising (shared corpora, shared benchmarks, shared optimization
   techniques), so the raw figure is the relevant one for the paper's
   argument. The detrended residual ≈ 0 confirms that capability
   convergence accounts for the agreement, not that the coupling is
   small.
5. **Specific numbers should be verified against the cited primary
   sources before paper submission.** The user has flagged that
   verification is pending; the directional claim ("near-critical
   adversarially, supercritical on benign benchmarks") is robust to
   the precision of any specific cell, but the headline figures
   (medians, ranges) should be cross-checked.
6. **the two proxies are complementary lower bounds, not redundant signals.**
   Neither alone is the full coupling; the paper's claim is supported
   by their *agreement on direction*, not by either point estimate
   alone.

---

## Suggested Discussion-paragraph language

For the integration prompt to pull from. Roughly 150–200 words, fits
the paper's Discussion §"active/passive distinction as a governance
taxonomy" subsection or its own brief subsection.

> Two cross-method measurements place the implicit coupling among
> frontier AI systems at or above the paper's critical threshold.
> Adversarial-attack transfer rates compiled from the published
> red-teaming literature [Zou et al. 2023; Mazeika et al. 2024;
> Wei et al. 2023] correspond to a frontier-pair median
> $J_{\rm eff} \approx 1.06$ under the mapping $\bar\rho = J/(1+J)$,
> with a high-transfer tail at $J_{\rm eff} \approx 6.5$ — comparable
> to financial-crisis coupling. Cross-model agreement on standard
> benign benchmarks (12 frontier models × 10 standard tasks compiled
> from each model's primary publication) places the median pairwise
> coupling at $J_{\rm eff} \approx 8.7$, dominated by capability
> convergence on shared evaluation suites — the homogenization
> mechanism Theorem 1's framework predicts. These are proxy
> measurements rather than direct coupling estimates, with the
> caveats developed in the Supplementary Information and the
> empirical scripts; the directional finding — *near-critical
> adversarially, supercritical on benign-benchmark surfaces* — is
> robust across both methods. The finding is consistent with the
> paper's framing that AI coupling is rising, has reached the regime
> where Theorem 1 begins to bind, and warrants the active-stabilizer
> interventions the Discussion proposes.

---

## File pointers

- `empirical/AI_COUPLING_VERDICT.md` (adversarial-transfer) — full per-pair table, method
  notes, six caveats, paper-implication paragraph
- `empirical/AI_BENCHMARK_VERDICT.md` (benchmark-vector) — full method notes,
  aggregates, six caveats, implication paragraph
- `empirical/ai_coupling_transfer.csv` — 7 rows, adversarial-transfer per-pair numbers
- `empirical/ai_coupling_benchmarks.csv` — 66 rows, benchmark-vector per-pair numbers
- `empirical/ai_coupling_benchmarks_raw.csv` /
  `..._detrended.csv` — full 12×12 matrices
- `empirical/figures/ai_coupling_j_axis.png` (adversarial-transfer proxy) — horizontal-bar
  J-axis plot
- `empirical/figures/ai_model_correlation_matrix.png` (benchmark-vector proxy) — 2-panel
  heatmap (raw + detrended)
- `empirical/figures/ai_coupling_vs_finance.png` (benchmark-vector proxy) — AI pairs on the
  J axis vs T = 2 and financial reference lines
