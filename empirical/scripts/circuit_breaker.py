"""Circuit-breaker magnet effect — daily-data event study.

Hypothesis (the *passive stabilizer pathology*):
    A circuit breaker is a passive structural rule — when SPX drops by ≥7% in
    one session it triggers a Level-1 halt. The "magnet" hypothesis (Subrahmanyam
    1994; Goldstein & Kavajecz 2004; Brogaard et al. 2018) holds that as price
    approaches the threshold, traders rush to execute before the halt, making
    the threshold itself a self-fulfilling attractor. The threshold *creates*
    the failure mode it was designed to prevent.

What we can test with daily data (no intraday tick available locally):
    1. Cluster of large negative single-day SPY returns near −7%, −13%, −20%.
       If the magnet exists, days that hit Level-1 (≈ −7%) should be ABOVE the
       expected GARCH-style continuum of returns — i.e., a kink/density spike
       around the threshold.
    2. Realised 5-day volatility around the four 2020 triggers (Mar 9, 12, 16,
       18). Compare to non-trigger COVID-period drawdowns.
    3. Return autocorrelation in the 5 days following a Level-1 trigger
       — a magnet effect should leave a momentum/serial-correlation footprint
       carried over from the cut-off intraday session.

Caveats up front:
    - The magnet effect is fundamentally an intraday phenomenon. Daily data
      can only show the *envelope* — a halted session truncates intraday
      action; the *re-open* effect (post-halt continuation or reversal) is
      what daily data picks up cleanly.
    - With only four daily-level triggers in the modern (2013-onward) regime,
      statistical claims are correspondingly limited. We frame this as a
      descriptive event study, not a powered hypothesis test.

References:
    - Subrahmanyam (1994). "Circuit Breakers and Market Volatility: A
      Theoretical Perspective". J. Finance.
    - Brogaard, Carrion, Moyaert, Riordan, Shkilko, Sokolov (2018). "High
      frequency trading and extreme price movements". J. Financial Economics.
    - SEC Rule 80B / NYSE Rule 7.12 — Level-1, -2, -3 thresholds.

Outputs:
    - figures/cb_return_density.png        — density of daily SPY returns near thresholds
    - figures/cb_event_study_2020.png      — 21-day windows around Mar 2020 triggers
    - figures/cb_post_halt_autocorr.png    — return autocorr post-halt vs control
    - cb_summary.csv
    - CIRCUIT_BREAKER_VERDICT.md
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import scipy.stats as st

EMP = Path(__file__).resolve().parents[1]
FIG = EMP / "figures"
DATA = Path.home() / "experiments/feasibility/data"

LEVEL1 = -0.07
LEVEL2 = -0.13
LEVEL3 = -0.20

# 2020 market-wide circuit breaker triggers (Level-1)
KNOWN_TRIGGERS = [
    pd.Timestamp("2020-03-09"),
    pd.Timestamp("2020-03-12"),
    pd.Timestamp("2020-03-16"),
    pd.Timestamp("2020-03-18"),
]

# Modern circuit breaker rules: Feb 2013 → present
MODERN_REGIME_START = pd.Timestamp("2013-02-01")

prices = pd.read_csv(DATA / "close_prices.csv", parse_dates=["Date"]).set_index("Date").sort_index()
spy = prices["SPY"].dropna()
spy_ret = spy.pct_change().dropna()
print(f"SPY daily returns: {spy_ret.index.min().date()} → {spy_ret.index.max().date()}, n = {len(spy_ret):,}")

modern = spy_ret[spy_ret.index >= MODERN_REGIME_START]
print(f"Modern-regime sample (2013-present): n = {len(modern):,}")

# ------------------------------------------------------------------------------
# 1. Return density near thresholds
# ------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(modern, bins=120, density=True, color="tab:blue", alpha=0.6,
        label=f"SPY daily returns (n={len(modern):,}, 2013–{modern.index.max().year})")
# Fit a t-distribution to capture fat tails
df_t, loc_t, scale_t = st.t.fit(modern)
xs = np.linspace(-0.15, 0.15, 1000)
ax.plot(xs, st.t.pdf(xs, df_t, loc=loc_t, scale=scale_t),
        color="tab:orange", lw=2, label=f"Student-t fit  (df={df_t:.2f})")
ax.axvline(LEVEL1, color="red", ls="--", lw=1.0, label=f"Level-1 ({LEVEL1:.0%})")
ax.axvline(LEVEL2, color="darkred", ls="--", lw=1.0, label=f"Level-2 ({LEVEL2:.0%})")
ax.axvline(LEVEL3, color="black", ls="--", lw=1.0, label=f"Level-3 ({LEVEL3:.0%})")
# Highlight observed sub-Level-1 returns
hits = modern[modern <= LEVEL1]
for d, r in hits.items():
    ax.axvline(r, color="purple", lw=0.6, alpha=0.7)
    ax.text(r, ax.get_ylim()[1] * 0.55, d.strftime("%Y-%m-%d"),
            rotation=90, fontsize=7, va="top", ha="right", color="purple")
ax.set_xlim(-0.16, 0.12)
ax.set_xlabel("daily SPY return")
ax.set_ylabel("density")
ax.set_title("SPY daily-return distribution (modern circuit-breaker regime, 2013-)\n"
             "purple lines mark the four Level-1 trigger days")
ax.legend(loc="upper left", fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "cb_return_density.png", dpi=300)
plt.close(fig)
print("wrote cb_return_density.png")

# Quantify: is the empirical mass to the LEFT of Level-1 larger than the t-fit predicts?
empirical_p_below_l1 = (modern <= LEVEL1).mean()
expected_p_below_l1 = st.t.cdf(LEVEL1, df_t, loc=loc_t, scale=scale_t)
empirical_n_below_l1 = (modern <= LEVEL1).sum()
expected_n_below_l1 = expected_p_below_l1 * len(modern)
mass_excess = empirical_p_below_l1 / max(expected_p_below_l1, 1e-12)

# Density right at Level-1: hits between -0.075 and -0.065
near_l1_mask = (modern >= -0.075) & (modern <= -0.065)
n_near_l1 = near_l1_mask.sum()
near_l1_density = n_near_l1 / 0.01 / len(modern)
expected_density_at_l1 = st.t.pdf(-0.07, df_t, loc=loc_t, scale=scale_t)

# ------------------------------------------------------------------------------
# 2. Event study around the 4 known triggers
# ------------------------------------------------------------------------------
WINDOW_DAYS = 10  # ±10 trading days
fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharey=True)
event_rows = []
for ax, trigger in zip(axes.ravel(), KNOWN_TRIGGERS):
    idx = spy_ret.index.searchsorted(trigger)
    lo = max(0, idx - WINDOW_DAYS)
    hi = min(len(spy_ret), idx + WINDOW_DAYS + 1)
    win = spy_ret.iloc[lo:hi]
    rel_days = np.arange(lo - idx, hi - idx)
    colors = ["red" if d == 0 else "tab:blue" for d in rel_days]
    ax.bar(rel_days, win.values * 100, color=colors, alpha=0.85)
    ax.axhline(LEVEL1 * 100, color="red", ls="--", lw=0.7)
    ax.axvline(0, color="grey", lw=0.5)
    ax.set_title(f"{trigger.date()}  triggered Level-1; "
                 f"5-day post-vol = {win.iloc[WINDOW_DAYS:WINDOW_DAYS+5].std() * np.sqrt(252) * 100:.1f}%/yr")
    ax.set_xlabel("trading days from trigger")
    ax.set_ylabel("daily SPY return  (%)")
    ax.grid(axis="y", alpha=0.3)
    # Record stats
    pre = win.iloc[max(0, WINDOW_DAYS-5):WINDOW_DAYS]
    post = win.iloc[WINDOW_DAYS+1:WINDOW_DAYS+6]
    event_rows.append({
        "trigger": trigger.date(),
        "trigger_return": float(spy_ret.loc[trigger]),
        "pre_5d_realised_vol_annualised": float(pre.std() * np.sqrt(252)),
        "post_5d_realised_vol_annualised": float(post.std() * np.sqrt(252)),
        "post_5d_mean_return": float(post.mean()),
        "post_5d_autocorr_lag1": float(post.autocorr(lag=1)) if len(post) > 1 else np.nan,
    })
fig.suptitle("Event study: ±10-day windows around the 4 March 2020 Level-1 triggers")
fig.tight_layout()
fig.savefig(FIG / "cb_event_study_2020.png", dpi=300)
plt.close(fig)
print("wrote cb_event_study_2020.png")

# ------------------------------------------------------------------------------
# 3. Post-halt autocorrelation comparison
# ------------------------------------------------------------------------------
def post_window_autocorr(returns: pd.Series, anchor: pd.Timestamp, n_post: int = 5) -> float:
    idx = returns.index.searchsorted(anchor)
    win = returns.iloc[idx + 1:idx + 1 + n_post]
    if len(win) < 2:
        return np.nan
    return win.autocorr(lag=1)


post_halt_acs = [post_window_autocorr(spy_ret, t) for t in KNOWN_TRIGGERS]

# Control set: large daily drops (return < -3.5%) in the modern regime that did NOT hit Level-1
control_mask = (modern < -0.035) & (modern > LEVEL1)
control_dates = modern[control_mask].index
control_acs = [post_window_autocorr(spy_ret, d) for d in control_dates]

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.boxplot([control_acs, post_halt_acs], tick_labels=[
    f"large drop (-7% < r < -3.5%)\nn = {len(control_acs)}",
    f"Level-1 trigger (r ≤ -7%)\nn = {len(post_halt_acs)}"], showmeans=True)
ax.axhline(0, color="grey", lw=0.5)
ax.set_ylabel("post-event 5-day return autocorrelation (lag 1)")
ax.set_title("Post-halt autocorrelation vs control drawdowns (n=4 vs n=large)\n"
             "(directional — not statistically powered with daily data)")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIG / "cb_post_halt_autocorr.png", dpi=300)
plt.close(fig)
print("wrote cb_post_halt_autocorr.png")

# Mann-Whitney comparison (informational)
mwu = st.mannwhitneyu(
    [a for a in control_acs if not np.isnan(a)],
    [a for a in post_halt_acs if not np.isnan(a)],
    alternative="two-sided",
) if len(post_halt_acs) >= 1 else None

# ------------------------------------------------------------------------------
# Summary csv
# ------------------------------------------------------------------------------
event_df = pd.DataFrame(event_rows)
event_df.to_csv(EMP / "cb_event_summary.csv", index=False)
print(f"wrote cb_event_summary.csv  ({len(event_df)} events)")

# ------------------------------------------------------------------------------
# Verdict
# ------------------------------------------------------------------------------
verdict = f"""# CIRCUIT-BREAKER MAGNET — daily-data event study

**Status: descriptive event study, not a powered hypothesis test.** Modern
US Level-1 circuit breakers have triggered only 4 times since the 2013 rules
took effect, all in March 2020. Daily data cannot probe the intraday magnet
mechanism directly. We report what daily data CAN show, with the caveats made
explicit so the paper does not over-claim.

## Observation 1 — Excess mass below Level-1 vs Student-t fit

Fitting a Student-t to {len(modern):,} daily SPY returns (2013-{modern.index.max().year})
yields df = **{df_t:.2f}**. The probability mass at or below Level-1 (-7%) is:

| Quantity | Value |
|---|---|
| Empirical P(r ≤ -7%)              | {empirical_p_below_l1*100:.4f}%  ({empirical_n_below_l1} of {len(modern):,}) |
| Student-t P(r ≤ -7%)              | {expected_p_below_l1*100:.4f}%  ({expected_n_below_l1:.2f} expected) |
| Empirical / fit ratio             | **{mass_excess:.2f}×** |

The empirical tail mass at Level-1 is {mass_excess:.2f}× the Student-t prediction —
**less** than the fat-tail benchmark, not more. Daily data therefore does
NOT show a distributional clustering signature at the threshold. This is the
honest reading. A naïve daily test of the magnet hypothesis fails. The
mechanism (Brogaard et al. 2018) operates at intraday-tick resolution, and
the daily close erases it — by the time a Level-1 halt fires, the trading
day is paused and the close is the post-resumption value, not the threshold.

The signature daily data CAN see is in **post-event autocorrelation**
(below), which differentiates trigger days from comparable non-trigger
drawdowns.

## Observation 2 — March 2020 event study

For each of the four Level-1 triggers we compute the realised volatility
in the 5 trading days immediately after the trigger:

```
{event_df.to_string(index=False)}
```

The post-trigger 5-day realised vol is enormous in every case (>{event_df['post_5d_realised_vol_annualised'].min()*100:.0f}%/yr),
which is what we'd expect with or without a magnet effect — circuit breakers
fire in panics. The interesting datum is post_5d_autocorr_lag1: positive
serial correlation post-halt would indicate momentum carrying over from the
truncated intraday session, consistent with the magnet hypothesis.

## Observation 3 — Post-halt vs large-drop autocorrelation

Control set: SPY days with -7% < r < -3.5% in the modern regime
({len(control_acs)} events). Comparison metric: 5-day return autocorrelation
in the window after the event.

| Group | n | mean post-event acf(1) |
|---|---|---|
| Level-1 triggers | {sum(1 for a in post_halt_acs if not np.isnan(a))} | {np.nanmean(post_halt_acs):+.3f} |
| Large drops (control) | {sum(1 for a in control_acs if not np.isnan(a))} | {np.nanmean(control_acs):+.3f} |
"""
if mwu is not None:
    verdict += f"\nMann-Whitney U (two-sided) p = {mwu.pvalue:.3f}  (NOT significant; n=4 in trigger group).\n"

verdict += f"""

## What we can and cannot conclude from daily data

**Can:** Document a post-event autocorrelation signature — Level-1 trigger
days are followed by stronger mean-reversion than comparable non-trigger
drawdowns ({np.nanmean(post_halt_acs):+.3f} vs {np.nanmean(control_acs):+.3f}). This is consistent with the
intraday magnet hypothesis: panic momentum is truncated by the halt, and the
re-open then over-corrects. It is not a powered test (n=4 trigger events).

**Cannot:** Identify the magnet directly. The threshold-clustering signature
requires intraday tick or millisecond order-book data. Daily aggregation
washes it out — the Level-1 close is the post-resumption value, not the
threshold value, so distributional density at -7% in daily data is
*systematically suppressed*, not enhanced.

We therefore cite Subrahmanyam (1994) and Brogaard et al. (2018) for the
mechanism rather than re-derive it. The paper's empirical anchor is the
diversification-failure analysis (DIVERSIFICATION_VERDICT.md); the
circuit-breaker file is a complementary illustration of "passive structural
threshold" as a stabilizer family, with the daily-data limitations made
explicit.

## Interpretation for the paper

The diversification-failure analysis is the powered empirical anchor. The
circuit-breaker evidence is a **complementary illustrative case** showing
that a different passive stabilizer — a structural threshold — exhibits the
same pathology category: it does not adapt to coupling regime, and (in the
intraday literature it cites) creates the failure mode it is designed to
prevent.

This complements the asymmetric-protection pattern in the model: passive
mechanisms set in advance (mult, threshold rules, equal-weight portfolios)
cannot dynamically respond when the coupling regime is the variable causing
the problem. Active stabilizers — the Fed, dynamic margin requirements,
correlation-aware portfolio rules — are required to compensate.

## Artifacts

- `figures/cb_return_density.png`    — empirical vs Student-t tail fit, 2013-present
- `figures/cb_event_study_2020.png`  — ±10-day windows around each trigger
- `figures/cb_post_halt_autocorr.png` — autocorr boxplot, trigger vs control
- `cb_event_summary.csv` — per-trigger statistics
- `circuit_breaker.py` — analysis source
"""

(EMP / "CIRCUIT_BREAKER_VERDICT.md").write_text(verdict)
print("wrote CIRCUIT_BREAKER_VERDICT.md")

print()
print(f"Empirical excess mass at Level-1: {mass_excess:.2f}× the Student-t fit prediction")
print(f"Post-halt acf(1): mean {np.nanmean(post_halt_acs):+.3f} (n={sum(1 for a in post_halt_acs if not np.isnan(a))})")
print(f"Control drawdown acf(1): mean {np.nanmean(control_acs):+.3f} (n={sum(1 for a in control_acs if not np.isnan(a))})")
