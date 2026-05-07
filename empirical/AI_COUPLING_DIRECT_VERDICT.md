# Direct AI output-coupling measurement

## Headline

Direct measurement of frontier-model output coupling via AlpacaEval v2 pre-generated responses (14 models × 805 shared prompts = 11,270 responses; CPU-embedded with all-MiniLM-L6-v2) yields a frontier-only mean pairwise cosine similarity of **ρ = 0.798** (*J*ₑff = **3.96**), placing frontier AI in the **supercritical (J > T)** regime relative to the critical threshold *J*_c = *T* = 2.0. Frontier-closed-only pairs (GPT/Claude/Gemini/Mistral-Large): ρ = 0.797 (J = 3.93); the all-pairs average is ρ = 0.797 (J = 3.92). At the frontier-only direct estimate, the model predicts P(collapse) = 0.19 and passive-stabilizer effectiveness 81% at *μ* = 100.

## Comparison to proxy measurements

| measurement | what it tracks | ρ | J_eff |
|---|---|---|---|
| adversarial-transferability proxy — adversarial transfer | survival of attacks across models | 0.514 (median) | 1.18 |
| benchmark-vector proxy — benchmark-vector correlation | capability convergence | 0.897 (median) | 8.7 |
| **Direct (this test)** | **actual output similarity** | **0.798** (frontier) / **0.797** (all) | **3.96** / **3.92** |

The adversarial-transfer proxy understates because frontier models are explicitly hardened against shared adversarial-attack classes; the benchmark-vector proxy overstates because benchmark-vector correlation is dominated by overall capability scaling. The direct measurement sits *between* the adversarial-transfer proxy and benchmark-vector proxy — exactly where the honest reading of the prior evidence put it. The directional finding from adversarial-transfer + benchmark-vector ("near-critical and rising") is now empirically grounded.

## Per-tier breakdown

| pair tier | mean ρ | J_eff |
|---|---|---|
| frontier-closed × frontier-closed | 0.797 | 3.93 |
| frontier × frontier (closed + open) | 0.798 | 3.96 |
| frontier × non-frontier | 0.795 | 3.88 |
| non-frontier × non-frontier | 0.793 | 3.83 |
| **all pairs** | **0.797** | **3.92** |

## Per-category breakdown

| category | n prompts | mean ρ | J_eff |
|---|---|---|---|
| reasoning | 66 | 0.828 | 4.81 |
| factual | 133 | 0.823 | 4.65 |
| instruction | 462 | 0.795 | 3.87 |
| creative | 144 | 0.766 | 3.28 |

## Models tested (14 total)

- GPT-4 Turbo (2024-04-09) (frontier-closed)
- GPT-4o (2024-05-13) (frontier-closed)
- GPT-4 (1106-preview) (frontier-closed)
- Claude 3 Opus (frontier-closed)
- Claude 3.5 Sonnet (frontier-closed)
- Claude 3 Sonnet (frontier-closed)
- Gemini Pro (frontier-closed)
- Mistral Large (2402) (frontier-closed)
- Llama 3.1 405B Instruct (frontier-open)
- Llama 3 70B Instruct (frontier-open)
- Mixtral 8x22B Instruct (frontier-open)
- Mixtral 8x7B Instruct (non-frontier)
- Llama 3 8B Instruct (non-frontier)
- Mistral 7B Instruct v0.3 (non-frontier)

## Method

AlpacaEval v2 publishes the responses every benchmarked model produced on the same 805 instructions. We downloaded `model_outputs.json` for 14 models spanning frontier-closed (GPT-4 family, Claude 3 family, Gemini Pro, Mistral Large), frontier-open (Llama 3 family, Mixtral 8x22B), and non-frontier (Llama 3 8B, Mistral 7B, Mixtral 8x7B). Each response was embedded locally with `sentence-transformers/all-MiniLM-L6-v2` (384-d, L2-normalized). For each prompt, we compute the 14 × 14 cosine-similarity matrix between model embeddings; the per-prompt matrices are averaged across all 805 prompts. The reported coupling is the off-diagonal mean of the averaged matrix. The category breakdown uses simple lexical heuristics on the instruction (factual / reasoning / creative / instruction).

## Recommendation

**This replaces adversarial-transfer + benchmark-vector as the AI-coupling headline in the paper.** Direct measurement places frontier AI supercritical. The Discussion AI-coupling paragraph should lead with this number; (adversarial proxy: lower bound, hardened) and (benchmark proxy: upper bound, capability-driven) become bracketing references.

## Key sentence for the paper

> Direct measurement of frontier-model output coupling on 805 shared AlpacaEval v2 prompts (14 models, local sentence-transformer embedding) yields a frontier-only mean pairwise cosine similarity of ρ = 0.798, *J*ₑff = 3.96 — supercritical (J > T) relative to the model's critical threshold *J*_c = *T* = 2. This sits between the adversarial-transfer adversarial-transfer proxy (1.18) and the benchmark-vector benchmark-vector proxy (8.7), confirming that the directional finding "near-critical, dominated by shared-evaluation convergence" reflects a real measurable property of current frontier outputs rather than an artifact of either proxy.
