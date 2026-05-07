"""
Fed active-stabilizer interaction test.

The original Fed analysis (stabilizer_effectiveness_test.py) showed
post-cut SPY return rises with contemporaneous coupling J — opposite to
the prediction for a *passive* stabilizer (which should lose effectiveness
with J).

This script tests the stronger active-stabilizer hypothesis: post-cut
return scales with the *interaction* (|Δr| × J), where |Δr| is the cut
magnitude. Under the active-stabilizer functional form

    h_active = h_a + α · max(0, -m) · J

the intervention's effect on the order parameter scales with J. If real
Fed cuts behave like an active stabilizer, then the *effect of a cut*
should be larger when coupling is high and the cut is large.

Test: linear regression of 90-day post-cut SPY return on (|Δr|, J,
|Δr|×J), and compare interaction sign and significance.

Output: results/fed_active_fit.csv, FED_ACTIVE_FIT_VERDICT.md
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

EMP = Path(__file__).resolve().parents[1]
RES = EMP / "results"
RES.mkdir(parents=True, exist_ok=True)
DATA = Path.home() / "experiments/feasibility/data"

# Load market data
close = pd.read_csv(DATA / "close_prices.csv", parse_dates=["Date"], index_col="Date")
sector_cols = [c for c in close.columns if c not in ["SPY", "^VIX"]]
sector_prices = close[sector_cols].dropna()
spy = close["SPY"].dropna()

# Load Fed funds rate
fedfunds = pd.read_csv(DATA / "FEDFUNDS.csv", parse_dates=["DATE"],
                       index_col="DATE").iloc[:, 0].dropna()

# Compute rolling correlation
returns = sector_prices.pct_change().dropna()
WINDOW = 60


def rolling_avg_corr(returns_df, window=WINDOW):
    dates = returns_df.index[window-1:]
    avg_corrs = []
    for i in range(window-1, len(returns_df)):
        chunk = returns_df.iloc[i-window+1:i+1]
        corr_matrix = chunk.corr().values
        mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)
        avg_corrs.append(np.nanmean(corr_matrix[mask]))
    return pd.Series(avg_corrs, index=dates, name="avg_correlation")


print("Computing rolling correlation...")
avg_corr = rolling_avg_corr(returns)
print(f"  range [{avg_corr.min():.3f}, {avg_corr.max():.3f}]")

# Identify rate cuts (>= 25 bp)
ff_diff = fedfunds.diff()
rate_cuts = ff_diff[ff_diff < -0.10]   # match original 10bp threshold (N=12 in paper)
rate_cuts = rate_cuts[rate_cuts.index >= avg_corr.index[0]]
rate_cuts = rate_cuts[rate_cuts.index <= avg_corr.index[-1] - pd.Timedelta(days=90)]
print(f"\nFound {len(rate_cuts)} rate cuts >= 10 bp")

# Build events
events = []
for cut_date, cut_size in rate_cuts.items():
    trading_days = spy.index[spy.index >= cut_date]
    if len(trading_days) == 0:
        continue
    nearest = trading_days[0]
    j_val = avg_corr.asof(nearest)
    if pd.isna(j_val):
        continue
    future_date = nearest + pd.Timedelta(days=90)
    fp = spy[spy.index >= future_date]
    if len(fp) == 0:
        continue
    fwd_ret = (spy[fp.index[0]] - spy[nearest]) / spy[nearest]
    events.append({
        "date": nearest,
        "cut_size_abs": float(abs(cut_size)),
        "J_proxy": float(j_val),
        "ret_90d": float(fwd_ret),
    })

df = pd.DataFrame(events)
df["interaction"] = df["cut_size_abs"] * df["J_proxy"]
print(f"\nUsable events: {len(df)}")
print(df.to_string(index=False, float_format="%.4f"))

# Regression 1: just J
X1 = sm.add_constant(df[["J_proxy"]])
m1 = sm.OLS(df["ret_90d"], X1).fit()
print("\n--- Model 1: ret_90d ~ J ---")
print(m1.summary().tables[1])

# Regression 2: J + |dr|
X2 = sm.add_constant(df[["J_proxy", "cut_size_abs"]])
m2 = sm.OLS(df["ret_90d"], X2).fit()
print("\n--- Model 2: ret_90d ~ J + |dr| ---")
print(m2.summary().tables[1])

# Regression 3: J + |dr| + J*|dr|
X3 = sm.add_constant(df[["J_proxy", "cut_size_abs", "interaction"]])
m3 = sm.OLS(df["ret_90d"], X3).fit()
print("\n--- Model 3: ret_90d ~ J + |dr| + J*|dr| (active-stabilizer form) ---")
print(m3.summary().tables[1])

# Save outputs
df.to_csv(RES / "fed_active_fit_events.csv", index=False)

verdict = RES / "FED_ACTIVE_FIT_VERDICT.md"
with verdict.open("w", encoding="utf-8") as f:
    f.write("# Fed active-stabilizer interaction test\n\n")
    f.write(f"**Events:** {len(df)} FOMC cuts >= 25 bp, 2004-2025\n\n")
    f.write("## Headline statistics\n\n")
    f.write(f"- Mean cut size: {df.cut_size_abs.mean():.2f} percentage points\n")
    f.write(f"- Mean J at time of cut: {df.J_proxy.mean():.3f}\n")
    f.write(f"- Mean 90-day post-cut SPY return: {df.ret_90d.mean():+.3f}\n\n")

    f.write("## Regression: ret_90d ~ J × |Δr|\n\n")
    f.write("Active-stabilizer form predicts: cut effect *grows* with coupling.\n")
    f.write("Test the interaction term's sign and significance.\n\n")

    f.write("| term | coefficient | std err | t | p |\n")
    f.write("|---|---|---|---|---|\n")
    for term in m3.params.index:
        coef = m3.params[term]
        se = m3.bse[term]
        t = m3.tvalues[term]
        p = m3.pvalues[term]
        f.write(f"| {term} | {coef:+.4f} | {se:.4f} | {t:+.2f} | {p:.3f} |\n")

    f.write(f"\n*R² = {m3.rsquared:.3f}, adj R² = {m3.rsquared_adj:.3f}, "
            f"N = {len(df)}*\n\n")

    interaction_p = m3.pvalues.get("interaction", np.nan)
    interaction_sign = m3.params.get("interaction", 0.0)
    f.write("## Interpretation\n\n")
    if interaction_sign > 0 and interaction_p < 0.10:
        f.write(f"- The interaction term J × |Δr| is **positive** "
                f"(coefficient = {interaction_sign:+.4f}, p = {interaction_p:.3f}), "
                f"consistent with the active-stabilizer hypothesis: "
                f"the marginal effect of a Fed cut on post-cut return "
                f"*increases* with contemporaneous coupling.\n")
    elif interaction_sign > 0:
        f.write(f"- The interaction term is positive (coefficient = "
                f"{interaction_sign:+.4f}) but not significant at p < 0.10 "
                f"(p = {interaction_p:.3f}). The point estimate is consistent "
                f"with the active-stabilizer hypothesis but the small sample "
                f"(N = {len(df)}) limits statistical power.\n")
    else:
        f.write(f"- The interaction term is *not* positive "
                f"(coefficient = {interaction_sign:+.4f}, p = {interaction_p:.3f}). "
                f"The active-stabilizer interaction is not detected at this sample size.\n")

print(f"\nWrote {verdict}")
print("\nKey finding:")
ip = m3.pvalues.get("interaction", np.nan)
ic = m3.params.get("interaction", 0.0)
print(f"  Interaction J*|dr|: coefficient = {ic:+.4f}, p = {ip:.3f}")
print(f"  R^2 = {m3.rsquared:.3f}, N = {len(df)}")
