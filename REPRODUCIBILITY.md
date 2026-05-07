# Reproducibility

**Manuscript:** *"Why Fixed Protections Fail Under Rising Coordination: A Structural Failure-Mode Transition in Coupled Systems"* (submitted to *Nature Communications*).

**Author:** Chowon Jung.

This file is the public-facing reproducibility manifest for the empirical and computational results in the manuscript and Supplementary Information. Every analysis is described here with the script that produced it.

## Repository layout

```
.
├── manuscript/
│   ├── paper.md
│   ├── supplementary_information.md
│   ├── cover_letter.md
│   ├── arxiv_metadata.md
│   └── references.bib
├── simulation/
│   ├── scripts/         (21 Python scripts: model + sweeps + figures)
│   └── results/         (CSVs and figure PNGs by analysis)
├── empirical/
│   ├── scripts/         (15 Python scripts: financial + AI + crypto + politics)
│   ├── data/            (AlpacaEval JSONs and OWID-mirrored V-Dem CSVs)
│   ├── results/         (CSVs)
│   ├── figures/         (PNG figures)
│   └── *_VERDICT.md     (per-analysis result writeups)
├── README.md
├── REPRODUCIBILITY.md   (this file)
└── LICENSE
```

## One-time data setup (read before running the empirical scripts)

Most empirical scripts read from a local cache directory at `~/experiments/feasibility/data/`. To reproduce, populate that directory with the following CSVs first:

- `close_prices.csv` — daily adjusted-close prices for SPY plus the 10 sector ETFs (XLB, XLC, XLE, XLF, XLI, XLK, XLP, XLRE, XLU, XLV, XLY) and `^VIX`, with columns `Date,SPY,XLB,...` from 2004-01-02 onwards. Source: Yahoo Finance via `yfinance`.
- `FEDFUNDS.csv`, `NFCI.csv`, `BAA10Y.csv` — FRED series with `DATE` in column 1 and the value in column 2. Source: FRED via `fredapi`.

A 30-line setup script that downloads and writes these CSVs is straightforward (`yfinance.download(tickers, start='2004-01-02')` and `fredapi.Fred(api_key=...).get_series(...)`); we do not redistribute the data because Yahoo and FRED license terms vary.

`empirical/scripts/crypto_extension.py` and `empirical/scripts/ai_coupling_direct.py` are self-contained: the first downloads crypto prices on first run and caches locally; the second reads the included AlpacaEval JSON files in `empirical/data/alpacaeval_outputs/` (these are redistributed under the AlpacaEval Apache 2.0 license).

`simulation/scripts/*.py` need no external data — they generate their own simulation outputs from random seeds.

## How to reproduce the headline numbers

### Main-text empirical results (Figures 7, 8, S2, S3, S6)

| Result | Script | Output |
|---|---|---|
| Diversification benefit (Equation 3, Figure S2) | `empirical/scripts/diversification_failure.py` | `empirical/diversification_summary.csv`, `empirical/figures/diversification_failure_*.png` |
| Tail-risk frequency by coupling decile (Figure 7) | `empirical/scripts/tail_risk_frequency.py` | `empirical/results/tail_risk_by_coupling.csv` |
| VaR exceedance by decile (Figure 7) | `empirical/scripts/var_exceedance.py` | `empirical/results/var_exceedance_by_coupling.csv` |
| Failure-mode transition (Figure 8) | `empirical/scripts/rigidity_fragmentation.py` | `empirical/results/rigidity_fragmentation_by_coupling.csv` |
| Out-of-sample temporal split | `empirical/scripts/out_of_sample_test.py` | `empirical/results/oos_prediction.csv` |
| AI direct coupling, AlpacaEval (Figure S3) | `empirical/scripts/ai_coupling_direct.py` | `empirical/results/ai_coupling_direct.csv`, `empirical/figures/ai_coupling_direct_heatmap.png` |
| Crypto cross-domain extension (§S5.4, Figure S6) | `empirical/scripts/crypto_extension.py` | `empirical/results/crypto_*.csv`, `empirical/results/CRYPTO_VERDICT.md` |
| Fed active-stabilizer interaction test (SI §S10) | `empirical/scripts/fed_active_stabilizer_fit.py` | `empirical/results/FED_ACTIVE_FIT_VERDICT.md` |

### Main-text simulation results (Figures 1–6)

| Result | Script | Output |
|---|---|---|
| Minimal model (Figure 1, 9,000 runs) | `simulation/scripts/run_sweep.py` | `simulation/results/minimal_model/*.csv` |
| L0 mean-field validation (Figure S1) | `simulation/scripts/level0_pure_cw.py` + `check_l0.py` | `simulation/results/ablation/level0_*.csv` |
| Ablation cascade L1→L2→L3 (Figure 2) | `simulation/scripts/level2_with_beliefs.py`, `level3_with_credit.py` | `simulation/results/ablation/*.csv` |
| Network ABM 6 topologies (Figure 3) | `simulation/scripts/network_abm.py` | `simulation/results/network_abm/*.csv` |
| Network ABM 100-seed re-run (SI §S9) | `simulation/scripts/network_abm_100seed.py` | `simulation/results/network_abm/headline_100seed.csv` |
| Dynamic-J / ramp speed (Figure 4) | `simulation/scripts/dynamic_j.py`, `ramp_speed_sweep.py` | `simulation/results/dynamic_j/*.csv` |
| Active vs passive stabilizer (Figure 5) | `simulation/scripts/active_stabilizer.py` | `simulation/results/active_stabilizer/*.csv` |
| AI overlay (Figure 6) | `simulation/scripts/ai_coupling_overlay.py` | `simulation/results/ai_coupling_overlay/ai_coupling_overlay_figure.png` |
| Convergence check (Methods, SI §S4) | `simulation/scripts/convergence_check.py` | `simulation/results/minimal_model/convergence_check.csv` |
| Sensitivity sweep (SI §S5.1, Figure S4) | `simulation/scripts/sensitivity_sweep.py` + `sensitivity_figure.py` | `simulation/results/sensitivity/*.csv`, `figure_s4_sensitivity_sweep.png` |
| OAT sensitivity (SI §S5.2) | `simulation/scripts/oat_sensitivity.py` | `simulation/results/sensitivity/oat_sensitivity.csv` |
| Alternative-class sweep (SI §S5.3, Figure S5) | `simulation/scripts/alternative_class_sweep.py` + `alternative_class_figure.py` | `simulation/results/alternative_class/*.csv`, `figure_s5_alternative_class.png` |

## Pre-registered analysis design

The empirical analyses in the manuscript are pre-registered in this repository in the sense that:

- The temporal out-of-sample split (2004–2014 train, 2015–2025 test) is the only split reported in the main text.
- The 60-day rolling window, the 5% drawdown threshold, the 0.8/0.3 rigidity/fragmentation correlation thresholds, and the J = ρ̄/(1 − ρ̄) mapping are all documented in the manuscript Methods and SI §S10 before any output is computed.
- The 9,000-run sweep grid (J ∈ {0.5, …, 5.0}, μ ∈ {20, …, 100}, 100 seeds) and dt = 0.05 integration step are pre-specified in the script defaults; all reported results derive from the same configuration.

If reviewers request alternative analyses (different temporal splits, alternative sector weights, instrumental-variable specifications), these would be run with the existing pipeline and added to this manifest as appendices marked **post-review**.

## Software environment

Python 3.13 with the following pinned-by-import dependencies:

```
numpy >= 1.26
pandas >= 2.0
matplotlib >= 3.8
scipy >= 1.11
statsmodels >= 0.14
yfinance >= 0.2  (for crypto download)
sentence-transformers >= 2.7  (for AI direct coupling; CPU-compatible)
networkx >= 3.0  (for network ABM topologies)
```

## Data sources

| Dataset | Source | Citation |
|---|---|---|
| Daily SPY + 9 sector ETFs (2004–2025) | Yahoo Finance via `yfinance` | — |
| FRED macroeconomic series (FEDFUNDS, NFCI, BAA10Y) | FRED | Federal Reserve Bank of St. Louis |
| 10-token cryptocurrency prices (2020–2025) | Yahoo Finance | — |
| AlpacaEval v2 model outputs (14 frontier models × 805 prompts) | tatsu-lab/alpaca_eval | Li et al. 2023 |
| V-Dem democracy indices | OWID mirror of V-Dem | Coppedge et al. (V-Dem) |

## License

Code: MIT (see LICENSE). Data: as licensed by upstream providers (FRED public domain, V-Dem CC-BY 4.0, AlpacaEval Apache 2.0).
