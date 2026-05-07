"""
Stabilizer Effectiveness × Coupling Level
==========================================================================

MODEL PREDICTION:
  When coupling (J) is LOW:
    → Fed rate cuts (stabilizer) should IMPROVE market outcomes (recovery, reduced drawdown)
  When coupling (J) is HIGH:
    → Fed rate cuts should have DIMINISHED or ZERO effect
  This ASYMMETRY is the core testable prediction. The test maps directly
  to the model's asymmetric protection finding (passive stabilizer
  effectiveness decays with coupling).

METHODOLOGY:
  1. Identify all Fed rate cut events (FEDFUNDS decreases) since 2004
  2. For each cut: measure the coupling level (J) at the time
  3. For each cut: measure the market response (SPY return over 30/60/90 days after)
  4. Test: does market response to rate cuts WEAKEN as J increases?
  
  Additional: split by stress level (NFCI) to test the viability-margin interaction

DATA: All from experiments/feasibility/data/ (already downloaded)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# ══════════════════════════════════════════════════════════════════════════════
# PATHS
# ══════════════════════════════════════════════════════════════════════════════

FEAS_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                         "..", "..", "experiments", "feasibility", "data")
FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ══════════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════

print("Loading data...")

# Market prices
close = pd.read_csv(os.path.join(FEAS_DATA, "close_prices.csv"), parse_dates=["Date"], index_col="Date")
sector_cols = [c for c in close.columns if c not in ["SPY", "^VIX"]]
sector_prices = close[sector_cols].dropna()
spy = close["SPY"].dropna()

# FRED data
def load_fred(name):
    path = os.path.join(FEAS_DATA, f"{name}.csv")
    if os.path.exists(path):
        df = pd.read_csv(path, parse_dates=["DATE"], index_col="DATE")
        return df.iloc[:, 0].dropna()
    return None

fedfunds = load_fred("FEDFUNDS")
nfci = load_fred("NFCI")
baa_spread = load_fred("BAA10Y")

# ══════════════════════════════════════════════════════════════════════════════
# 2. COMPUTE J PROXY: Rolling Pairwise Correlation
# ══════════════════════════════════════════════════════════════════════════════

print("Computing rolling correlations...")
returns = sector_prices.pct_change().dropna()
WINDOW = 60

def rolling_avg_correlation(returns_df, window=WINDOW):
    dates = returns_df.index[window-1:]
    avg_corrs = []
    for i in range(window-1, len(returns_df)):
        chunk = returns_df.iloc[i-window+1:i+1]
        corr_matrix = chunk.corr().values
        mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)
        pairwise = corr_matrix[mask]
        avg_corrs.append(np.nanmean(pairwise))
    return pd.Series(avg_corrs, index=dates, name="avg_correlation")

avg_corr = rolling_avg_correlation(returns)
print(f"  Correlation: {len(avg_corr)} days, range [{avg_corr.min():.3f}, {avg_corr.max():.3f}]")

# ══════════════════════════════════════════════════════════════════════════════
# 3. IDENTIFY FED RATE CUT EVENTS
# ══════════════════════════════════════════════════════════════════════════════

print("\nIdentifying Fed rate cuts...")

if fedfunds is None:
    print("ERROR: No Fed Funds data")
    sys.exit(1)

# Monthly data — find months where rate decreased
ff_diff = fedfunds.diff()
rate_cuts = ff_diff[ff_diff < -0.10]  # cuts of at least 10bp

# Filter to our data range
rate_cuts = rate_cuts[rate_cuts.index >= avg_corr.index[0]]
rate_cuts = rate_cuts[rate_cuts.index <= avg_corr.index[-1] - pd.Timedelta(days=90)]  # need 90d forward

print(f"  Found {len(rate_cuts)} rate cut events (≥10bp)")
for date, cut in rate_cuts.items():
    print(f"    {date.strftime('%Y-%m-%d')}: {cut:+.2f}%")

# ══════════════════════════════════════════════════════════════════════════════
# 4. MEASURE STABILIZER EFFECTIVENESS
# ══════════════════════════════════════════════════════════════════════════════

print("\nMeasuring stabilizer effectiveness...")

FORWARD_WINDOWS = [30, 60, 90]  # days after rate cut

events = []
for cut_date, cut_size in rate_cuts.items():
    # Find nearest trading day
    trading_days = spy.index[spy.index >= cut_date]
    if len(trading_days) == 0:
        continue
    nearest = trading_days[0]
    
    # J at time of cut
    j_val = avg_corr.asof(nearest)
    if pd.isna(j_val):
        continue
    
    # NFCI stress at time of cut
    stress_val = nfci.asof(nearest) if nfci is not None else np.nan
    
    # Forward returns after cut
    forward = {}
    for w in FORWARD_WINDOWS:
        future_date = nearest + pd.Timedelta(days=w)
        future_prices = spy[spy.index >= future_date]
        if len(future_prices) > 0:
            future_nearest = future_prices.index[0]
            fwd_ret = (spy[future_nearest] - spy[nearest]) / spy[nearest]
            forward[f"ret_{w}d"] = fwd_ret
        else:
            forward[f"ret_{w}d"] = np.nan
    
    # SPY return in 30 days BEFORE cut (to control for pre-existing trend)
    past_date = nearest - pd.Timedelta(days=30)
    past_prices = spy[spy.index <= nearest]
    past_prices_start = past_prices[past_prices.index >= past_date]
    if len(past_prices_start) > 0:
        pre_ret = (spy[nearest] - past_prices_start.iloc[0]) / past_prices_start.iloc[0]
    else:
        pre_ret = np.nan
    
    events.append({
        "date": nearest,
        "cut_size": cut_size,
        "J": j_val,
        "stress": stress_val,
        "pre_30d_ret": pre_ret,
        **forward,
    })

df = pd.DataFrame(events)
print(f"\n{'='*80}")
print(f"RATE CUT EVENTS WITH COUPLING LEVEL ({len(df)} events)")
print(f"{'='*80}")
print(df.to_string(index=False, float_format="%.4f"))

# ══════════════════════════════════════════════════════════════════════════════
# 5. THE ASYMMETRY TEST
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print("ASYMMETRY TEST: Does stabilizer effectiveness decline with coupling?")
print(f"{'='*80}")

for window in FORWARD_WINDOWS:
    col = f"ret_{window}d"
    valid = df.dropna(subset=["J", col])
    if len(valid) < 4:
        print(f"\n  {window}-day: Too few events ({len(valid)})")
        continue
    
    r, p = stats.pearsonr(valid["J"], valid[col])
    # Also Spearman (more robust to outliers)
    rho, p_spearman = stats.spearmanr(valid["J"], valid[col])
    
    print(f"\n  {window}-day forward return after rate cut:")
    print(f"    Pearson:  r = {r:+.3f}, p = {p:.3f}")
    print(f"    Spearman: ρ = {rho:+.3f}, p = {p_spearman:.3f}")
    
    if r < -0.15:
        print(f"    → CONSISTENT WITH MODEL: Higher J → weaker stabilizer effect")
    elif r > 0.15:
        print(f"    → OPPOSITE OF MODEL: Higher J → stronger stabilizer effect")
    else:
        print(f"    → NO CLEAR PATTERN")

# Split by stress level for interaction test
if nfci is not None and len(df.dropna(subset=["stress"])) >= 6:
    print(f"\n{'='*80}")
    print("INTERACTION TEST: Is the J effect stronger when economy is stressed?")
    print(f"{'='*80}")
    
    stress_median = df["stress"].median()
    for window in [60]:
        col = f"ret_{window}d"
        
        high_stress = df[(df["stress"] >= stress_median)].dropna(subset=["J", col])
        low_stress = df[(df["stress"] < stress_median)].dropna(subset=["J", col])
        
        print(f"\n  {window}-day, stress median = {stress_median:.3f}")
        
        if len(high_stress) >= 3:
            r_hi, p_hi = stats.pearsonr(high_stress["J"], high_stress[col])
            print(f"  HIGH stress (N={len(high_stress)}): r = {r_hi:+.3f}, p = {p_hi:.3f}")
        
        if len(low_stress) >= 3:
            r_lo, p_lo = stats.pearsonr(low_stress["J"], low_stress[col])
            print(f"  LOW stress  (N={len(low_stress)}): r = {r_lo:+.3f}, p = {p_lo:.3f}")
        
        if len(high_stress) >= 3 and len(low_stress) >= 3:
            print(f"\n  Model prediction: r_high_stress should be MORE NEGATIVE than r_low_stress")
            print(f"  Result: r_high={r_hi:+.3f} vs r_low={r_lo:+.3f}")
            if r_hi < r_lo:
                print(f"  → CONSISTENT: Stabilizers less effective in stressed + high-J conditions")
            else:
                print(f"  → NOT CONSISTENT")

# ══════════════════════════════════════════════════════════════════════════════
# 6. ALSO TEST: Unconditional correlation increase over time
#    (Evidence that J is actually increasing — needed for the AI narrative)
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print("SECULAR COUPLING TREND: Is J increasing over time?")
print(f"{'='*80}")

# Annual average correlation
annual_corr = avg_corr.resample("YE").mean()
years = np.array([(d - annual_corr.index[0]).days / 365.25 for d in annual_corr.index])
if len(annual_corr) >= 3:
    slope, intercept, r, p, se = stats.linregress(years, annual_corr.values)
    print(f"  Annual avg correlation trend: slope={slope:+.4f}/year, r={r:.3f}, p={p:.3f}")
    if slope > 0 and p < 0.1:
        print(f"  → CONSISTENT: Market coupling is increasing over time")
    else:
        print(f"  → NO SIGNIFICANT INCREASE detected in this window")

# ══════════════════════════════════════════════════════════════════════════════
# 7. FIGURES
# ══════════════════════════════════════════════════════════════════════════════

print(f"\nGenerating figures...")

# Figure 1: Stabilizer effectiveness vs J
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for i, window in enumerate(FORWARD_WINDOWS):
    ax = axes[i]
    col = f"ret_{window}d"
    valid = df.dropna(subset=["J", col])
    
    if len(valid) > 0:
        # Color by stress
        if "stress" in valid.columns:
            scatter = ax.scatter(valid["J"], valid[col]*100, 
                               c=valid["stress"], cmap="RdYlGn_r", 
                               s=np.abs(valid["cut_size"])*500, 
                               edgecolors="black", linewidth=0.5,
                               alpha=0.8)
            plt.colorbar(scatter, ax=ax, label="Financial Stress (NFCI)")
        else:
            ax.scatter(valid["J"], valid[col]*100, s=80, edgecolors="black", linewidth=0.5)
        
        # Regression line
        if len(valid) >= 3:
            z = np.polyfit(valid["J"], valid[col]*100, 1)
            x_line = np.linspace(valid["J"].min(), valid["J"].max(), 50)
            ax.plot(x_line, np.polyval(z, x_line), "r--", alpha=0.7, linewidth=2)
            r, p = stats.pearsonr(valid["J"], valid[col])
            ax.set_title(f"{window}-Day Return After Rate Cut\nr={r:.2f}, p={p:.2f}")
        
        # Label events
        for _, row in valid.iterrows():
            ax.annotate(row["date"].strftime("%Y-%m"), 
                       (row["J"], row[col]*100),
                       fontsize=6, alpha=0.6, ha="center", va="bottom")
    
    ax.set_xlabel("Coupling (J) at Rate Cut")
    ax.set_ylabel(f"{window}-Day Forward Return (%)")
    ax.axhline(0, color="gray", ls=":", alpha=0.5)

plt.suptitle("Does Higher Coupling Weaken Stabilizer (Rate Cut) Effectiveness?", 
             fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "stabilizer_effectiveness.png"), dpi=150, bbox_inches="tight")
print(f"  Saved stabilizer_effectiveness.png")

# Figure 2: Coupling trend over time
fig2, ax2 = plt.subplots(figsize=(12, 5))
ax2.plot(avg_corr.index, avg_corr.values, color="steelblue", linewidth=0.6, alpha=0.7)
# Add trend line
yearly = avg_corr.resample("YE").mean()
ax2.plot(yearly.index, yearly.values, "ko-", markersize=5, linewidth=2, label="Annual Average")
# Mark rate cuts
for _, row in df.iterrows():
    ax2.axvline(row["date"], color="red", alpha=0.3, linewidth=1)
ax2.set_ylabel("Average Pairwise Sector Correlation (J proxy)")
ax2.set_title("Coupling Strength Over Time with Fed Rate Cut Events (red lines)")
ax2.legend()
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "coupling_trend.png"), dpi=150, bbox_inches="tight")
print(f"  Saved coupling_trend.png")

# ══════════════════════════════════════════════════════════════════════════════
# 8. VERDICT
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print("STABILIZER EFFECTIVENESS VERDICT")
print(f"{'='*80}")
print("""
GO signals (any 2 of 3):
  1. J vs forward return: NEGATIVE correlation (r < -0.2) for at least one window
  2. Interaction: J effect is STRONGER under high stress
  3. Secular coupling trend: J is increasing over time

STOP signals:
  1. No pattern in J vs forward returns
  2. Interaction effect absent or reversed
  3. Coupling not increasing (undermines the AI narrative)

DECISION: _______________
""")
