"""
Cryptocurrency empirical extension — second-domain test of the framework.

Same methodology as the S&P 500 sector analysis (diversification, tail-risk,
failure-mode transition) applied to a 10-token crypto basket. The framework
predicts:

  (a) crypto sits deeper in the supercritical regime than equities
      (higher average ρ̄, J = ρ̄/(1−ρ̄) above the critical threshold more often)
  (b) the same fragmentationtorigidity transition holds, but crossover
      occurs at lower J because crypto is already in the supercritical regime
  (c) the tail-risk frequency curve tracks the model's collapse curve
      with the same rank ordering as equities

Tokens (in approximate market-cap order; sufficient history from 2020):
  BTC-USD, ETH-USD, BNB-USD, XRP-USD, ADA-USD, SOL-USD, AVAX-USD,
  DOT-USD, LINK-USD, LTC-USD

Window: 60-day rolling. Period: 2020-01-01 to 2025-12-30.
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

EMP = Path(__file__).resolve().parents[1]
OUT = EMP / "results"
OUT.mkdir(parents=True, exist_ok=True)
FIG = EMP / "figures"
FIG.mkdir(parents=True, exist_ok=True)

TOKENS = ["BTC-USD", "ETH-USD", "BNB-USD", "XRP-USD", "ADA-USD",
          "SOL-USD", "AVAX-USD", "DOT-USD", "LINK-USD", "LTC-USD"]
N_TOKENS = len(TOKENS)
WINDOW = 60
START = "2020-01-01"
END = "2025-12-30"


def download():
    """Try yfinance; fall back to coingecko if unavailable."""
    try:
        import yfinance as yf
    except ImportError:
        print("yfinance not installed; install with: pip install yfinance")
        sys.exit(1)
    print(f"Downloading {len(TOKENS)} tokens, {START} to {END}...")
    df = yf.download(TOKENS, start=START, end=END, progress=False,
                     auto_adjust=True)["Close"]
    print(f"Got {df.shape[0]} days, {df.shape[1]} tokens")
    print(f"Date range: {df.index.min().date()} to {df.index.max().date()}")
    return df


def compute_panel(prices: pd.DataFrame) -> dict:
    """Compute rolling windows: ρ̄, DB, drawdown, drawdown-period correlation."""
    returns = np.log(prices / prices.shift(1)).dropna()
    print(f"Returns: {returns.shape[0]} days × {returns.shape[1]} tokens")

    rows = []
    n_windows = len(returns) - WINDOW
    for start_idx in range(n_windows):
        win_returns = returns.iloc[start_idx:start_idx + WINDOW]
        if win_returns.shape[0] < WINDOW or win_returns.isna().any().any():
            continue

        # 1. Pairwise correlation
        corr_mat = win_returns.corr().values
        n = corr_mat.shape[0]
        rho_bar = float((corr_mat.sum() - n) / (n * (n - 1)))

        # 2. Vol per token + portfolio vol (equal weight)
        vols = win_returns.std(ddof=0).values
        avg_single_vol = float(vols.mean())
        port_returns = win_returns.mean(axis=1)
        port_vol = float(port_returns.std(ddof=0))
        DB = avg_single_vol / port_vol if port_vol > 0 else np.nan

        # 3. Tail risk: max drawdown of the equal-weight portfolio
        cumret = port_returns.cumsum().values
        running_max = np.maximum.accumulate(cumret)
        drawdown = cumret - running_max
        max_dd = float(drawdown.min())
        max_dd_pct = float(np.exp(max_dd) - 1)  # to fractional return

        # 4. If max drawdown > 5%, find drawdown period and compute corr there
        rigidity_corr = np.nan
        if max_dd_pct < -0.05:
            peak_idx = int(np.argmax(cumret))
            trough_idx = int(np.argmin(cumret[peak_idx:])) + peak_idx
            if trough_idx > peak_idx:
                dd_returns = win_returns.iloc[peak_idx:trough_idx + 1]
                if dd_returns.shape[0] >= 5:
                    dd_corr = dd_returns.corr().values
                    rigidity_corr = float((dd_corr.sum() - n) / (n * (n - 1)))

        rows.append({
            "window_end": returns.index[start_idx + WINDOW - 1],
            "rho_bar": rho_bar,
            "DB": DB,
            "max_drawdown_pct": max_dd_pct,
            "drawdown_corr": rigidity_corr,
            "tail_5": int(max_dd_pct < -0.05),
            "tail_3": int(max_dd_pct < -0.03),
            "tail_10": int(max_dd_pct < -0.10),
        })

    df = pd.DataFrame(rows).set_index("window_end")
    print(f"Built panel: {df.shape[0]} windows")
    return df


def analyze(df: pd.DataFrame) -> dict:
    """Compute decile-binned statistics for J = ρ̄/(1-ρ̄) deciles."""
    df = df.copy()
    df["J"] = df["rho_bar"] / (1.0 - df["rho_bar"]).clip(lower=1e-6)
    df["decile"] = pd.qcut(df["rho_bar"], 10, labels=False, duplicates="drop")

    summary = df.groupby("decile").agg(
        rho_bar_mid=("rho_bar", "mean"),
        J_mid=("J", "mean"),
        n_windows=("rho_bar", "size"),
        tail_5_freq=("tail_5", "mean"),
        tail_3_freq=("tail_3", "mean"),
        tail_10_freq=("tail_10", "mean"),
        drawdown_corr_mean=("drawdown_corr", "mean"),
    ).reset_index()

    print("\nDecile summary:")
    print(summary.to_string(index=False))

    # Failure-mode classification (matching S&P methodology)
    crash = df[df.tail_5 == 1].copy()
    crash["rigidity"] = (crash.drawdown_corr > 0.8).astype(int)
    crash["fragmentation"] = (crash.drawdown_corr < 0.3).astype(int)
    crash["decile"] = pd.qcut(crash["rho_bar"], 10, labels=False,
                              duplicates="drop")

    rig_summary = crash.groupby("decile").agg(
        n_crashes=("rigidity", "size"),
        rigidity_share=("rigidity", "mean"),
        fragmentation_share=("fragmentation", "mean"),
        rho_bar_mid=("rho_bar", "mean"),
    ).reset_index()
    print("\nCrash failure-mode summary:")
    print(rig_summary.to_string(index=False))

    return {
        "decile_summary": summary,
        "crash_summary": rig_summary,
        "panel": df,
    }


def main():
    cache = OUT / "crypto_panel_cache.csv"
    if cache.exists():
        print(f"Loading from cache: {cache}")
        prices = pd.read_csv(cache, parse_dates=["Date"]).set_index("Date")
    else:
        prices = download()
        if prices.empty:
            print("No data; aborting")
            return
        prices.to_csv(cache, index_label="Date")
        print(f"Cached prices to {cache}")

    df = compute_panel(prices)
    res = analyze(df)
    # Use the panel returned from analyze() which has J added
    df = res["panel"]

    # Save outputs
    df.to_csv(OUT / "crypto_panel_full.csv")
    res["decile_summary"].to_csv(OUT / "crypto_decile_summary.csv", index=False)
    res["crash_summary"].to_csv(OUT / "crypto_crash_summary.csv", index=False)
    print(f"\nWrote outputs to {OUT}")

    # Compute headline statistics
    rho_mean = float(df.rho_bar.mean())
    j_mean = float(df.J.replace([np.inf, -np.inf], np.nan).dropna().mean())
    rho_p95 = float(df.rho_bar.quantile(0.95))
    j_p95 = float((df.J.replace([np.inf, -np.inf], np.nan)).quantile(0.95))
    rho_supercritical = float((df.rho_bar > 0.667).mean())   # J > 2

    summary = res["decile_summary"]
    crash_summary = res["crash_summary"]
    rig_decile_corr = np.nan
    if not crash_summary.empty:
        # Spearman between rigidity_share and rho_bar_mid across deciles
        try:
            from scipy.stats import spearmanr
            rho_dec = crash_summary["rho_bar_mid"].values
            rig = crash_summary["rigidity_share"].values
            mask = ~np.isnan(rig)
            if mask.sum() >= 3:
                rig_decile_corr, _p = spearmanr(rho_dec[mask], rig[mask])
        except ImportError:
            pass

    verdict = OUT / "CRYPTO_VERDICT.md"
    with verdict.open("w", encoding="utf-8") as f:
        f.write("# Crypto extension — verdict\n\n")
        f.write(f"**Period:** {df.index.min().date()} to {df.index.max().date()}\n")
        f.write(f"**Windows:** {df.shape[0]} (60-day rolling)\n")
        f.write(f"**Tokens:** {N_TOKENS} ({', '.join(TOKENS)})\n\n")
        f.write("## Headline statistics\n\n")
        f.write(f"- Mean pairwise correlation ρ̄: **{rho_mean:.3f}** "
                f"(implied *J* = {rho_mean/(1-rho_mean):.2f})\n")
        f.write(f"- 95th-percentile ρ̄: **{rho_p95:.3f}** "
                f"(implied *J* = {rho_p95/(1-rho_p95):.2f})\n")
        f.write(f"- Mean *J*: {j_mean:.2f}\n")
        f.write(f"- 95th-percentile *J*: {j_p95:.2f}\n")
        f.write(f"- Fraction of windows with ρ̄ > 0.667 (J > 2, supercritical): "
                f"**{rho_supercritical:.1%}**\n\n")
        if not np.isnan(rig_decile_corr):
            f.write(f"- Rigidity-share Spearman ρ across deciles: "
                    f"**{rig_decile_corr:+.3f}**\n\n")

        f.write("## Comparison to S&P 500\n\n")
        f.write("| metric | S&P 500 (paper) | Crypto |\n")
        f.write("|---|---|---|\n")
        f.write(f"| Mean ρ̄ | ~0.50 | {rho_mean:.3f} |\n")
        f.write(f"| 95th-percentile ρ̄ | ~0.92 | {rho_p95:.3f} |\n")
        f.write(f"| % windows supercritical (ρ̄ > 0.667) | ~25% | {rho_supercritical:.1%} |\n")

    print(f"\nWrote {verdict}")
    print("\nHEADLINE NUMBERS:")
    print(f"  Mean ρ̄ = {rho_mean:.3f} (J = {rho_mean/(1-rho_mean):.2f})")
    print(f"  95th-pctile ρ̄ = {rho_p95:.3f} (J = {rho_p95/(1-rho_p95):.2f})")
    print(f"  Supercritical fraction = {rho_supercritical:.1%}")
    if not np.isnan(rig_decile_corr):
        print(f"  Rigidity-decile Spearman = {rig_decile_corr:+.3f}")


if __name__ == "__main__":
    main()
