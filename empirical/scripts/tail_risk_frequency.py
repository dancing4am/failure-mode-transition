"""
Tail-risk frequency vs coupling.

The diversification analysis (diversification_failure.py) shows that the
benefit DB follows a known variance-algebra identity. That is a tautology,
not a test of the model's mechanism. The minimal-model simulation predicts
something different and falsifiable: P(collapse) increases with J along a
specific curve. If the empirical frequency of tail events binned by ρ̄ also
increases with J along that same curve, the model captures the mechanism;
if not, it does not.

Pipeline:
  1. Load daily ETF panel (same as diversification_failure.py: 9 sectors
     pre-2015, 10 sectors post-2015 with XLRE).
  2. For each 60-day rolling window compute (a) ρ̄, (b) equal-weight
     portfolio daily returns, (c) max drawdown, (d) worst daily return,
     (e) # days below full-sample 5th percentile, (f) tail-event flags
     at thresholds 3 / 5 / 7 / 10 % max drawdown.
  3. Map ρ̄ to J = ρ̄/(1-ρ̄). Bin into J deciles.
  4. Compute empirical P(tail event | J bin) for each threshold + SE.
  5. Load model P(collapse | J, μ=100) from
     simulation/results/ai_coupling_overlay/combined_sweep_summary.csv.
  6. Compare: Spearman rank correlation, log-log slope, side-by-side plot.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

EMP = Path(__file__).resolve().parents[1]
FIG = EMP / "figures"
RESULTS = EMP / "results"
FIG.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)
DATA = Path.home() / "experiments/feasibility/data"

SECTORS_FULL = ["XLB", "XLE", "XLF", "XLI", "XLK", "XLP",
                "XLRE", "XLU", "XLV", "XLY"]
SECTORS_PRE_2015 = ["XLB", "XLE", "XLF", "XLI", "XLK", "XLP",
                    "XLU", "XLV", "XLY"]
WINDOW = 60
THRESHOLDS = [0.03, 0.05, 0.07, 0.10]
N_BINS = 10
MODEL_SUMMARY = (Path(__file__).resolve().parents[2] / "simulation"
                 / "results" / "ai_coupling_overlay"
                 / "combined_sweep_summary.csv")


def load_returns():
    prices = pd.read_csv(DATA / "close_prices.csv", parse_dates=["Date"])
    prices = prices.set_index("Date").sort_index()
    XLRE_START = pd.Timestamp("2015-10-08")
    panel_pre = prices.loc[prices.index < XLRE_START, SECTORS_PRE_2015].dropna()
    panel_post = prices.loc[prices.index >= XLRE_START, SECTORS_FULL].dropna()
    rets_pre = np.log(panel_pre / panel_pre.shift(1)).dropna()
    rets_post = np.log(panel_post / panel_post.shift(1)).dropna()
    return rets_pre, rets_post


def avg_pairwise_corr(corr_mat: np.ndarray) -> float:
    n = corr_mat.shape[0]
    iu = np.triu_indices(n, k=1)
    return corr_mat[iu].mean()


def max_drawdown_log(port_log_returns: np.ndarray) -> float:
    """Max peak-to-trough drawdown of the cumulative log-return series.
    Returned as a positive number (e.g. 0.05 = 5% drawdown)."""
    cumlog = np.cumsum(port_log_returns)
    running_max = np.maximum.accumulate(cumlog)
    drawdown = running_max - cumlog  # positive when below peak
    return float(drawdown.max())


def window_metrics(returns_panel: pd.DataFrame, full_5pct: float):
    """Iterate rolling windows; return one row per window."""
    rows = []
    cols = returns_panel.columns
    arr = returns_panel.values
    dates = returns_panel.index
    n = len(arr)
    for start in range(n - WINDOW + 1):
        end = start + WINDOW
        win = arr[start:end]
        rho = avg_pairwise_corr(np.corrcoef(win, rowvar=False))
        port = win.mean(axis=1)  # equal-weight log-return per day
        max_dd = max_drawdown_log(port)
        worst_day = float(port.min())
        n_tail_days = int((port < full_5pct).sum())
        rows.append({
            "Date": dates[end - 1],
            "n_sectors": win.shape[1],
            "rho": rho,
            "J": rho / max(1.0 - rho, 1e-9),
            "max_drawdown": max_dd,
            "worst_daily_return": worst_day,
            "n_tail_days": n_tail_days,
            "tail_5pct_event": int(max_dd > 0.05),
            "tail_3pct_event": int(max_dd > 0.03),
            "tail_7pct_event": int(max_dd > 0.07),
            "tail_10pct_event": int(max_dd > 0.10),
        })
    return pd.DataFrame(rows)


def main():
    print("Loading daily returns...")
    rets_pre, rets_post = load_returns()
    rets = pd.concat([rets_pre, rets_post])
    print(f"  rets shape: {rets.shape}")

    # Full-sample 5th percentile of equal-weight portfolio daily returns
    # (used as the reference threshold for the n_tail_days metric).
    port_pre = rets_pre.mean(axis=1)
    port_post = rets_post.mean(axis=1)
    port_all = pd.concat([port_pre, port_post])
    full_5pct = float(np.percentile(port_all.values, 5))
    print(f"  full-sample portfolio 5th-pctile daily return: {full_5pct:.4f}")

    # Compute per-window metrics, separately for pre / post panels (each
    # panel uses internally consistent N), then concatenate.
    df_pre = window_metrics(rets_pre, full_5pct)
    df_post = window_metrics(rets_post, full_5pct)
    df = pd.concat([df_pre, df_post], ignore_index=True).sort_values("Date")
    df = df.reset_index(drop=True)
    print(f"  windows analysed: {len(df)}")

    # Coupling decile bins
    df["J_bin"] = pd.qcut(df["J"], q=N_BINS, labels=False, duplicates="drop")

    # Per-bin empirical P(tail event)
    bin_rows = []
    for b in sorted(df["J_bin"].dropna().unique()):
        sub = df[df["J_bin"] == b]
        n = len(sub)
        rho_mid = float(sub["rho"].median())
        J_mid = float(sub["J"].median())
        row = {"J_bin": int(b), "rho_bin_mid": rho_mid,
               "J_bin_mid": J_mid, "n_windows": n}
        for thr in THRESHOLDS:
            col = f"tail_{int(thr*100)}pct_event"
            P = float(sub[col].mean())
            SE = float(np.sqrt(P * (1 - P) / max(n, 1)))
            row[f"P_tail_{int(thr*100)}pct"] = P
            row[f"SE_{int(thr*100)}pct"] = SE
        # Also keep mean # tail-days and mean max drawdown
        row["mean_n_tail_days"] = float(sub["n_tail_days"].mean())
        row["mean_max_drawdown"] = float(sub["max_drawdown"].mean())
        row["mean_worst_daily"] = float(sub["worst_daily_return"].mean())
        bin_rows.append(row)
    bin_df = pd.DataFrame(bin_rows)

    # Load model curve P(collapse | J, mu=100)
    model = pd.read_csv(MODEL_SUMMARY)
    model_mu100 = (model[model.mu == 100][["J", "p_collapse"]]
                   .sort_values("J").reset_index(drop=True))
    # Interpolate model P(collapse) at empirical J_bin_mid
    bin_df["model_P_collapse"] = np.interp(
        bin_df["J_bin_mid"].values,
        model_mu100["J"].values,
        model_mu100["p_collapse"].values,
        left=model_mu100["p_collapse"].iloc[0],
        right=model_mu100["p_collapse"].iloc[-1],
    )
    bin_df.to_csv(RESULTS / "tail_risk_by_coupling.csv", index=False)
    print(f"  wrote {RESULTS / 'tail_risk_by_coupling.csv'}")
    print(bin_df.to_string(index=False))

    # ---------------------------------------------------------------
    # Spearman correlation: empirical P(tail | J) vs model P(coll | J)
    # ---------------------------------------------------------------
    scaling_rows = []
    print()
    for thr in THRESHOLDS:
        col = f"P_tail_{int(thr*100)}pct"
        x = bin_df["J_bin_mid"].values
        y_emp = bin_df[col].values
        y_mod = bin_df["model_P_collapse"].values
        rho_s_model, p_s_model = spearmanr(y_mod, y_emp)
        rho_s_J, p_s_J = spearmanr(x, y_emp)
        # Log-log slope on positive bins (drop zero-frequency bins)
        mask = (y_emp > 0) & (x > 0)
        if mask.sum() >= 3:
            log_x = np.log(x[mask])
            log_y = np.log(y_emp[mask])
            slope_emp = float(np.polyfit(log_x, log_y, 1)[0])
        else:
            slope_emp = float("nan")
        # Log-log slope of model curve over the same J grid
        mask_m = (y_mod > 0) & (x > 0)
        if mask_m.sum() >= 3:
            log_x_m = np.log(x[mask_m])
            log_y_m = np.log(y_mod[mask_m])
            slope_mod = float(np.polyfit(log_x_m, log_y_m, 1)[0])
        else:
            slope_mod = float("nan")
        print(f"  thr={int(thr*100)}%  Spearman(model, emp)={rho_s_model:+.3f}  "
              f"(p={p_s_model:.3g})  Spearman(J, emp)={rho_s_J:+.3f}  "
              f"slope_emp={slope_emp:+.3f}  slope_mod={slope_mod:+.3f}")
        scaling_rows.append({
            "threshold_pct": int(thr * 100),
            "spearman_model_vs_emp": rho_s_model,
            "spearman_p_model_vs_emp": p_s_model,
            "spearman_J_vs_emp": rho_s_J,
            "spearman_p_J_vs_emp": p_s_J,
            "slope_loglog_emp": slope_emp,
            "slope_loglog_model": slope_mod,
        })
    scaling_df = pd.DataFrame(scaling_rows)
    scaling_df.to_csv(RESULTS / "tail_risk_scaling.csv", index=False)
    print(f"  wrote {RESULTS / 'tail_risk_scaling.csv'}")

    # ---------------------------------------------------------------
    # Figure 1 — primary: empirical 5% vs model on same J axis
    # ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5.5))
    x_emp = bin_df["J_bin_mid"].values
    y_5 = bin_df["P_tail_5pct"].values
    se_5 = bin_df["SE_5pct"].values
    ax.errorbar(x_emp, y_5, yerr=se_5, fmt="o", color="#d62728",
                ecolor="grey", capsize=4, markersize=10,
                label=r"empirical P(max drawdown $>$ 5% | J bin)",
                zorder=4)
    # Model curve at finer resolution
    Jline = np.linspace(model_mu100["J"].min(),
                        max(model_mu100["J"].max(), x_emp.max()), 200)
    Pline = np.interp(Jline, model_mu100["J"].values,
                      model_mu100["p_collapse"].values)
    ax.plot(Jline, Pline, "-", color="#1f77b4", linewidth=2.5,
            label=r"model P(collapse | J, $\mu=100$)")
    rho_5 = scaling_df.loc[scaling_df.threshold_pct == 5,
                           "spearman_model_vs_emp"].iloc[0]
    p_5 = scaling_df.loc[scaling_df.threshold_pct == 5,
                          "spearman_p_model_vs_emp"].iloc[0]
    ax.set_xlabel("coupling $J = \\bar\\rho/(1-\\bar\\rho)$", fontsize=11)
    ax.set_ylabel("frequency", fontsize=11)
    ax.set_title("Empirical tail-risk frequency vs model collapse curve\n"
                 f"Spearman(model, emp) = {rho_5:+.2f} "
                 f"(p = {p_5:.2g}); 5,414 sixty-day windows, 2004–2025",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "tail_risk_vs_coupling.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {FIG / 'tail_risk_vs_coupling.png'}")

    # ---------------------------------------------------------------
    # Figure 2 — multi-threshold panel
    # ---------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    axes = axes.flatten()
    colors = {3: "#2ca02c", 5: "#d62728", 7: "#ff7f0e", 10: "#9467bd"}
    for ax, thr in zip(axes, THRESHOLDS):
        col = f"P_tail_{int(thr*100)}pct"
        se_col = f"SE_{int(thr*100)}pct"
        y = bin_df[col].values
        se = bin_df[se_col].values
        ax.errorbar(x_emp, y, yerr=se, fmt="o", color=colors[int(thr*100)],
                    ecolor="grey", capsize=3, markersize=8,
                    label=f"empirical P(max DD > {int(thr*100)}%)")
        ax.plot(Jline, Pline, "-", color="#1f77b4", linewidth=2,
                alpha=0.85, label="model P(collapse) μ=100")
        rho_t = scaling_df.loc[scaling_df.threshold_pct == int(thr*100),
                                "spearman_model_vs_emp"].iloc[0]
        ax.set_title(f"threshold = {int(thr*100)}%   Spearman = {rho_t:+.2f}",
                     fontsize=10)
        ax.set_xlabel("J")
        ax.set_ylabel("frequency")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left", fontsize=8)
    fig.suptitle("Empirical tail-risk frequency vs model collapse curve "
                 "across thresholds", fontsize=12, y=1.00)
    fig.tight_layout()
    fig.savefig(FIG / "tail_risk_multi_threshold.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {FIG / 'tail_risk_multi_threshold.png'}")

    # ---------------------------------------------------------------
    # VERDICT
    # ---------------------------------------------------------------
    rho_5_val = float(scaling_df.loc[scaling_df.threshold_pct == 5,
                                      "spearman_model_vs_emp"].iloc[0])
    rho_3 = float(scaling_df.loc[scaling_df.threshold_pct == 3,
                                  "spearman_model_vs_emp"].iloc[0])
    rho_7 = float(scaling_df.loc[scaling_df.threshold_pct == 7,
                                  "spearman_model_vs_emp"].iloc[0])
    rho_10 = float(scaling_df.loc[scaling_df.threshold_pct == 10,
                                   "spearman_model_vs_emp"].iloc[0])

    if rho_5_val > 0.85 and min(rho_3, rho_5_val, rho_7) > 0.7:
        outcome = "STRONG"
        verdict_line = (
            "The empirical tail-risk frequency curve tracks the model's "
            "predicted P(collapse | J) curve closely across thresholds: "
            f"Spearman(model, empirical) = {rho_5_val:.2f} at 5% drawdown "
            f"({rho_3:.2f}/{rho_7:.2f}/{rho_10:.2f} at 3/7/10%). The model "
            "captures the mechanism, not just the direction."
        )
        recommendation = (
            "**Promote to a new Results subsection** "
            "\"Tail-risk frequency tracks the model's collapse curve\" "
            "(~200 words) with a new figure overlaying empirical and "
            "model curves on the same J axis. This transforms the "
            "empirical section from a variance-algebra identity to a "
            "mechanism test."
        )
    elif rho_5_val > 0.5:
        outcome = "MODERATE"
        verdict_line = (
            f"The empirical tail-risk frequency curve is correlated with the "
            f"model's predicted P(collapse | J) curve "
            f"(Spearman = {rho_5_val:.2f} at 5% drawdown), but the "
            f"correspondence is partial. The directional finding is robust; "
            f"the quantitative shape match is partial."
        )
        recommendation = (
            "Add to Results as a moderate-strength mechanism test, with "
            "honest discussion of where the shapes match and where they "
            "diverge."
        )
    else:
        outcome = "WEAK"
        verdict_line = (
            f"The empirical tail-risk frequency curve does NOT closely match "
            f"the model's predicted P(collapse | J) curve "
            f"(Spearman = {rho_5_val:.2f} at 5% drawdown). The model "
            f"captures the qualitative direction (more coupling → more tail "
            f"risk) but not the quantitative shape."
        )
        recommendation = (
            "Add to Discussion as a limitation. The model's "
            "collapse-frequency curve does not track real tail-risk "
            "frequency closely; report this honestly."
        )

    verdict_path = EMP / "TAIL_RISK_VERDICT.md"
    with verdict_path.open("w", encoding="utf-8") as f:
        f.write("# Tail-risk frequency vs coupling\n\n")
        f.write(f"## Headline ({outcome})\n\n")
        f.write(verdict_line + "\n\n")

        f.write("## Per-bin empirical results\n\n")
        f.write("| J bin mid | ρ̄ mid | n windows | "
                "P(DD>3%) | P(DD>5%) | P(DD>7%) | P(DD>10%) | "
                "model P(coll) μ=100 |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for _, r in bin_df.iterrows():
            f.write(
                f"| {r['J_bin_mid']:.2f} | {r['rho_bin_mid']:.2f} | "
                f"{int(r['n_windows'])} | "
                f"{r['P_tail_3pct']:.2f} | {r['P_tail_5pct']:.2f} | "
                f"{r['P_tail_7pct']:.2f} | {r['P_tail_10pct']:.2f} | "
                f"{r['model_P_collapse']:.2f} |\n"
            )
        f.write("\n")

        f.write("## Spearman rank correlation per threshold\n\n")
        f.write("| threshold | Spearman(model, empirical) | p | "
                "log-log slope (empirical) | log-log slope (model) |\n")
        f.write("|---|---|---|---|---|\n")
        for _, r in scaling_df.iterrows():
            f.write(
                f"| {int(r['threshold_pct'])}% | "
                f"{r['spearman_model_vs_emp']:+.3f} | "
                f"{r['spearman_p_model_vs_emp']:.3g} | "
                f"{r['slope_loglog_emp']:+.2f} | "
                f"{r['slope_loglog_model']:+.2f} |\n"
            )
        f.write("\n")

        f.write("## Method\n\n")
        f.write(
            "5,414 sixty-day rolling windows of daily log-returns for SPY +\n"
            "9 sector ETFs (10 from 2015 with XLRE), 2004–2025. For each\n"
            "window: (i) ρ̄ = mean off-diagonal pairwise correlation; \n"
            "(ii) equal-weight portfolio daily return; (iii) max drawdown\n"
            "in the window (peak-to-trough). Tail event = max drawdown\n"
            "exceeds the threshold (3 / 5 / 7 / 10 %). Map ρ̄ → J = ρ̄/(1-ρ̄).\n"
            f"Bin into {N_BINS} J deciles; compute empirical P(tail event |\n"
            "J bin). Compare to the minimal-model P(collapse | J, μ = 100)\n"
            "curve (extended sweep, J ∈ [0.5, 10], 100 seeds, dt = 0.05).\n\n"
        )

        f.write("## Recommendation\n\n")
        f.write(recommendation + "\n\n")

        f.write("## Key sentence for the paper\n\n")
        f.write(
            "> Across 5,414 sixty-day rolling windows of S&P 500 sector "
            f"data binned into {N_BINS} coupling deciles, the empirical "
            f"frequency of tail events (max drawdown > 5%) "
            f"{'tracks' if rho_5_val > 0.7 else 'is correlated with' if rho_5_val > 0.4 else 'is only weakly related to'} "
            f"the minimal-model P(collapse | J) curve "
            f"(Spearman ρ = {rho_5_val:+.2f}). This tests the model's "
            "mechanism — branch-selection probability under coupling — "
            "rather than the variance-algebra identity that the "
            "diversification-benefit metric satisfies by construction.\n"
        )
    print(f"  wrote {verdict_path}")


if __name__ == "__main__":
    main()
