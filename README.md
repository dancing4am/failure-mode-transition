# A Structural Failure-Mode Transition in Coupled Systems

Code and data for *"Why Fixed Protections Fail Under Rising Coordination: A Structural Failure-Mode Transition in Coupled Systems"* (Chowon Jung, submitted to *Nature Communications*).

## What's here

- **`manuscript/`** — paper, supplementary information, cover letter, arXiv metadata, and bibliography.
- **`simulation/`** — Python scripts and outputs for the minimal Curie–Weiss model, the L1–L3 ablation cascade, the six-topology agent-based extension, the dynamic-coupling sweep, the active-stabilizer experiment, the sensitivity analyses, and the alternative-model-class sweep.
- **`empirical/`** — Python scripts and outputs for the financial (S&P 500), cryptocurrency, AI-coupling (AlpacaEval), and political-science (V-Dem) analyses, plus per-analysis result writeups.
- **`REPRODUCIBILITY.md`** — full manifest mapping every figure and headline number to the script that produced it.

## Headline finding

As coupling rises in a coordinated system, the *type* of failure transitions from **fragmentation** (uncoordinated disorder) to **rigidity** (synchronized lock-in on the wrong branch). In 22 years of S&P 500 sector data the empirical type-composition curve correlates with the model's prediction at Spearman ρ = +0.93. The framework also predicts the same regime in cryptocurrency markets (64% of windows supercritical vs ~25% for equities) and locates frontier AI in the regime where passive protections decay.

The underlying scaling law: a fixed stabilizer field *h* opposing coupling *J* loses effectiveness as *h*/*J*² with an irreducible residual at high coupling. The pattern survives across alternative bistable model classes (Curie–Weiss, voter, Kuramoto, cubic Landau), six network topologies, and two empirical domains.

## Reproducing the results

See `REPRODUCIBILITY.md` for the full script-by-script manifest. Quick start:

```bash
# Simulation results (no external data needed)
python simulation/scripts/run_sweep.py
python simulation/scripts/network_abm.py
python simulation/scripts/sensitivity_sweep.py
python simulation/scripts/alternative_class_sweep.py

# Empirical results (see REPRODUCIBILITY.md for one-time data setup)
python empirical/scripts/diversification_failure.py
python empirical/scripts/tail_risk_frequency.py
python empirical/scripts/rigidity_fragmentation.py
python empirical/scripts/crypto_extension.py
python empirical/scripts/ai_coupling_direct.py
```

Python 3.13 with numpy, pandas, matplotlib, scipy, statsmodels, yfinance, sentence-transformers, and networkx.

## License

Code: MIT (see `LICENSE`). Data redistributed in `empirical/data/` is subject to upstream licenses (V-Dem CC-BY 4.0; AlpacaEval Apache 2.0). Financial and FRED data are not redistributed; download instructions are in `REPRODUCIBILITY.md`.

## Citation

```
Jung, C. Why Fixed Protections Fail Under Rising Coordination:
A Structural Failure-Mode Transition in Coupled Systems.
Manuscript submitted to Nature Communications (2026).
```

## Contact

Issues, reproductions, or alternative analyses → open an issue in this repository.
