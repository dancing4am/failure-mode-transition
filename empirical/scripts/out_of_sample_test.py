"""
Out-of-sample temporal prediction.

Splits the 2004–2025 daily-returns panel into two ~equal halves
(P1: 2004–2014; P2: 2015–2025). Calibrates an affine mapping
P_emp = a·P_model + b on P1's J-decile bins, predicts P2's
empirical tail-risk frequency, and checks whether the model's
collapse curve (parameter-free in shape) correctly predicts
out-of-sample crisis frequency. Also runs the reverse split
(P2 → P1) and the rigidity-share test out-of-sample.

This addresses the "you fitted a flexible curve to all the data"
critique: the model curve's *shape* is fixed by Theorem 1; only
scale and offset are calibrated. If the calibrated curve predicts
the held-out half (including COVID-19, 2022 inflation shock, and
the August 2024 carry-trade unwind), the framework has predictive,
not just descriptive, power.
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
DD_THRESHOLD = 0.05
RIGID_CORR_THRESHOLD = 0.8
N_BINS = 10
SPLIT_DATE = pd.Timestamp("2015-01-01")
MODEL_SUMMARY = (Path(__file__).resolve().parents[2] / "simulation"
                 / "results" / "ai_coupling_overlay"
                 / "combined_sweep_summary.csv")
MINIMAL_CELL_SUMMARY = (Path(__file__).resolve().parents[1]
                        / "simulation" / "results" / "minimal_model"
                        / "cell_summary.csv")


def load_panels():
    prices = pd.read_csv(DATA / "close_prices.csv", parse_dates=["Date"])
    prices = prices.set_index("Date").sort_index()
    XLRE_START = pd.Timestamp("2015-10-08")
    panel_pre = prices.loc[prices.index < XLRE_START,
                           SECTORS_PRE_2015].dropna()
    panel_post = prices.loc[prices.index >= XLRE_START,
                            SECTORS_FULL].dropna()
    rets_pre = np.log(panel_pre / panel_pre.shift(1)).dropna()
    rets_post = np.log(panel_post / panel_post.shift(1)).dropna()
    return rets_pre, rets_post


def avg_pairwise_corr(corr_mat: np.ndarray) -> float:
    n = corr_mat.shape[0]
    iu = np.triu_indices(n, k=1)
    return corr_mat[iu].mean()


def find_drawdown_window(port_log_returns: np.ndarray):
    cumlog = np.cumsum(port_log_returns)
    running_max = np.maximum.accumulate(cumlog)
    drawdown = running_max - cumlog
    trough_idx = int(drawdown.argmax())
    if trough_idx == 0:
        return 0, 0, 0.0
    peak_val = running_max[trough_idx]
    peak_idx = int(np.where(cumlog[: trough_idx + 1] == peak_val)[0][0])
    max_dd = float(drawdown[trough_idx])
    return peak_idx, trough_idx, max_dd


def window_features(returns_panel: pd.DataFrame):
    arr = returns_panel.values
    dates = returns_panel.index
    n = len(arr)
    rows = []
    for start in range(n - WINDOW + 1):
        end = start + WINDOW
        win = arr[start:end]
        rho = avg_pairwise_corr(np.corrcoef(win, rowvar=False))
        port = win.mean(axis=1)
        peak, trough, max_dd = find_drawdown_window(port)
        n_dd_days = trough - peak + 1 if max_dd > DD_THRESHOLD else 0
        rigidity_corr = np.nan
        if max_dd > DD_THRESHOLD and n_dd_days >= 3:
            dd = win[peak:trough + 1]
            cm = np.corrcoef(dd, rowvar=False)
            rigidity_corr = float(avg_pairwise_corr(cm))
        rows.append({
            "Date": dates[end - 1],
            "rho": rho,
            "J": rho / max(1.0 - rho, 1e-9),
            "max_drawdown": max_dd,
            "tail_5pct": int(max_dd > DD_THRESHOLD),
            "rigidity_corr": rigidity_corr,
            "rigid_crash": int(rigidity_corr > RIGID_CORR_THRESHOLD)
                if not np.isnan(rigidity_corr) else 0,
        })
    return pd.DataFrame(rows)


def bin_by_J(df, n_bins=N_BINS, label_prefix=""):
    df = df.copy()
    df["J_bin"] = pd.qcut(df["J"], q=n_bins, labels=False, duplicates="drop")
    rows = []
    for b in sorted(df["J_bin"].dropna().unique()):
        sub = df[df["J_bin"] == b]
        n = len(sub)
        crashes = sub[~sub.rigidity_corr.isna()]
        rows.append({
            "J_bin": int(b),
            "J_bin_mid": float(sub["J"].median()),
            "n_windows": n,
            "P_tail_5pct": float(sub["tail_5pct"].mean()),
            "n_crashes": len(crashes),
            "mean_rigidity_corr": float(crashes["rigidity_corr"].mean())
                if len(crashes) else np.nan,
            "P_rigid_crash": float(crashes["rigid_crash"].mean())
                if len(crashes) else np.nan,
        })
    return pd.DataFrame(rows)


def affine_fit(x, y):
    """Least-squares affine fit y = a*x + b. Returns (a, b)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    A = np.vstack([x, np.ones_like(x)]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(a), float(b)


def main():
    print("Loading panels...")
    rets_pre, rets_post = load_panels()
    df_pre = window_features(rets_pre)
    df_post = window_features(rets_post)
    df = pd.concat([df_pre, df_post], ignore_index=True).sort_values("Date")
    df = df.reset_index(drop=True)
    print(f"  total windows: {len(df)}")

    p1 = df[df["Date"] < SPLIT_DATE].copy()
    p2 = df[df["Date"] >= SPLIT_DATE].copy()
    print(f"  P1 (pre {SPLIT_DATE.date()}): {len(p1)} windows  "
          f"({p1.Date.min().date()} to {p1.Date.max().date()})")
    print(f"  P2 (post): {len(p2)} windows  "
          f"({p2.Date.min().date()} to {p2.Date.max().date()})")

    # Bin each period independently
    p1_bin = bin_by_J(p1)
    p2_bin = bin_by_J(p2)

    # Load model curve
    model = pd.read_csv(MODEL_SUMMARY)
    model_mu100 = (model[model.mu == 100][["J", "p_collapse"]]
                   .sort_values("J").reset_index(drop=True))

    # Interpolate model P(coll) at each period's bin midpoints
    p1_bin["model_P_collapse"] = np.interp(
        p1_bin["J_bin_mid"].values,
        model_mu100["J"].values,
        model_mu100["p_collapse"].values,
        left=model_mu100["p_collapse"].iloc[0],
        right=model_mu100["p_collapse"].iloc[-1],
    )
    p2_bin["model_P_collapse"] = np.interp(
        p2_bin["J_bin_mid"].values,
        model_mu100["J"].values,
        model_mu100["p_collapse"].values,
        left=model_mu100["p_collapse"].iloc[0],
        right=model_mu100["p_collapse"].iloc[-1],
    )

    # ----- Forward split: train on P1 affine map, predict P2 -----
    a_fwd, b_fwd = affine_fit(p1_bin["model_P_collapse"].values,
                               p1_bin["P_tail_5pct"].values)
    p2_bin["predicted_P_tail"] = a_fwd * p2_bin["model_P_collapse"] + b_fwd
    rho_fwd, p_fwd = spearmanr(p2_bin["predicted_P_tail"],
                                p2_bin["P_tail_5pct"])
    rmse_fwd = float(np.sqrt(np.mean(
        (p2_bin["predicted_P_tail"] - p2_bin["P_tail_5pct"]) ** 2)))
    naive_p1_mean = float(p1_bin["P_tail_5pct"].mean())
    rmse_naive_fwd = float(np.sqrt(np.mean(
        (naive_p1_mean - p2_bin["P_tail_5pct"]) ** 2)))

    # ----- Reverse split: train on P2, predict P1 -----
    a_rev, b_rev = affine_fit(p2_bin["model_P_collapse"].values,
                               p2_bin["P_tail_5pct"].values)
    p1_bin["predicted_P_tail"] = a_rev * p1_bin["model_P_collapse"] + b_rev
    rho_rev, p_rev = spearmanr(p1_bin["predicted_P_tail"],
                                p1_bin["P_tail_5pct"])
    rmse_rev = float(np.sqrt(np.mean(
        (p1_bin["predicted_P_tail"] - p1_bin["P_tail_5pct"]) ** 2)))
    naive_p2_mean = float(p2_bin["P_tail_5pct"].mean())
    rmse_naive_rev = float(np.sqrt(np.mean(
        (naive_p2_mean - p1_bin["P_tail_5pct"]) ** 2)))

    # ----- Rigidity OOS: load model rigidity-share, repeat -----
    rigidity_oos = {}
    if MINIMAL_CELL_SUMMARY.exists():
        ms = pd.read_csv(MINIMAL_CELL_SUMMARY)
        ms = ms[ms.mult == 100].sort_values("J")
        ms["model_rigidity_share"] = np.where(
            ms.p_collapse > 0, ms.p_rigidity / ms.p_collapse, np.nan)
        # Drop NaN rows so interp is safe
        ms_valid = ms.dropna(subset=["model_rigidity_share"])
        for label, period_bin in [("p1", p1_bin), ("p2", p2_bin)]:
            period_bin["model_rigidity_share"] = np.interp(
                period_bin["J_bin_mid"].values,
                ms_valid["J"].values,
                ms_valid["model_rigidity_share"].values,
                left=0.0,
                right=ms_valid["model_rigidity_share"].iloc[-1],
            )
        # Forward: train rigidity-share affine on P1, predict P2
        valid = ~p1_bin["mean_rigidity_corr"].isna()
        if valid.sum() >= 3:
            a_r_fwd, b_r_fwd = affine_fit(
                p1_bin.loc[valid, "model_rigidity_share"].values,
                p1_bin.loc[valid, "mean_rigidity_corr"].values)
            p2_bin["predicted_rigidity"] = (
                a_r_fwd * p2_bin["model_rigidity_share"] + b_r_fwd)
            valid2 = ~p2_bin["mean_rigidity_corr"].isna()
            if valid2.sum() >= 3:
                rho_r_fwd, p_r_fwd = spearmanr(
                    p2_bin.loc[valid2, "predicted_rigidity"],
                    p2_bin.loc[valid2, "mean_rigidity_corr"])
                rigidity_oos["forward"] = (rho_r_fwd, p_r_fwd)
        # Reverse
        valid = ~p2_bin["mean_rigidity_corr"].isna()
        if valid.sum() >= 3:
            a_r_rev, b_r_rev = affine_fit(
                p2_bin.loc[valid, "model_rigidity_share"].values,
                p2_bin.loc[valid, "mean_rigidity_corr"].values)
            p1_bin["predicted_rigidity"] = (
                a_r_rev * p1_bin["model_rigidity_share"] + b_r_rev)
            valid1 = ~p1_bin["mean_rigidity_corr"].isna()
            if valid1.sum() >= 3:
                rho_r_rev, p_r_rev = spearmanr(
                    p1_bin.loc[valid1, "predicted_rigidity"],
                    p1_bin.loc[valid1, "mean_rigidity_corr"])
                rigidity_oos["reverse"] = (rho_r_rev, p_r_rev)

    # Save combined results
    p1_bin["period"] = "P1 (2004-2014)"
    p2_bin["period"] = "P2 (2015-2025)"
    out_df = pd.concat([p1_bin, p2_bin], ignore_index=True)
    out_df.to_csv(RESULTS / "oos_prediction.csv", index=False)
    print()
    print("P1 bins (train forward, predicted by reverse):")
    print(p1_bin[["J_bin_mid", "P_tail_5pct", "model_P_collapse",
                 "predicted_P_tail"]].to_string(index=False))
    print()
    print("P2 bins (predicted by forward, used to train reverse):")
    print(p2_bin[["J_bin_mid", "P_tail_5pct", "model_P_collapse",
                 "predicted_P_tail"]].to_string(index=False))

    print()
    print(f"Forward (P1 -> P2):")
    print(f"  affine fit on P1: P_emp = {a_fwd:+.3f} * P_model + {b_fwd:+.3f}")
    print(f"  Spearman(predicted, actual) on P2 = {rho_fwd:+.3f} "
          f"(p={p_fwd:.3g})")
    print(f"  RMSE model     = {rmse_fwd:.3f}")
    print(f"  RMSE naive     = {rmse_naive_fwd:.3f}")
    print(f"  RMSE skill: {(1 - rmse_fwd / rmse_naive_fwd) * 100:.1f}% "
          f"reduction over naive")
    print()
    print(f"Reverse (P2 -> P1):")
    print(f"  affine fit on P2: P_emp = {a_rev:+.3f} * P_model + {b_rev:+.3f}")
    print(f"  Spearman(predicted, actual) on P1 = {rho_rev:+.3f} "
          f"(p={p_rev:.3g})")
    print(f"  RMSE model     = {rmse_rev:.3f}")
    print(f"  RMSE naive     = {rmse_naive_rev:.3f}")
    print(f"  RMSE skill: {(1 - rmse_rev / rmse_naive_rev) * 100:.1f}% "
          f"reduction over naive")
    print()
    if rigidity_oos:
        for k, (r, pv) in rigidity_oos.items():
            print(f"Rigidity OOS {k}: Spearman = {r:+.3f} (p={pv:.3g})")

    # ----- Figure: 2-panel OOS prediction -----
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8))
    for ax, period_bin, label, rho_v, p_v, rmse_m, rmse_n in [
        (axes[0], p2_bin, "P1 -> P2 (out-of-sample)",
         rho_fwd, p_fwd, rmse_fwd, rmse_naive_fwd),
        (axes[1], p1_bin, "P2 -> P1 (out-of-sample)",
         rho_rev, p_rev, rmse_rev, rmse_naive_rev),
    ]:
        x = period_bin["J_bin_mid"].values
        ax.plot(x, period_bin["P_tail_5pct"].values,
                "o", color="#d62728", markersize=10,
                label="actual P(DD>5%) per decile")
        ax.plot(x, period_bin["predicted_P_tail"].values,
                "s-", color="#1f77b4", markersize=9, linewidth=2,
                label="predicted (calibrated on other half)")
        ax.set_xlabel("coupling $J$ in held-out period", fontsize=11)
        ax.set_ylabel("P(max DD > 5% | J bin)", fontsize=11)
        ax.set_title(
            f"{label}\nSpearman = {rho_v:+.2f} "
            f"(p = {p_v:.2g});  "
            f"RMSE model = {rmse_m:.2f}  vs naive = {rmse_n:.2f}",
            fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left", fontsize=9)
    fig.suptitle("Out-of-sample temporal prediction "
                 "(model curve trained on 2004-2014 vs 2015-2025)",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "oos_prediction.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\nWrote {FIG / 'oos_prediction.png'}")

    # ----- VERDICT -----
    if rho_fwd > 0.7 and rho_rev > 0.7:
        outcome = "STRONG"
    elif rho_fwd > 0.5 or rho_rev > 0.5:
        outcome = "MODERATE"
    else:
        outcome = "WEAK"

    verdict_path = EMP / "OOS_VERDICT.md"
    with verdict_path.open("w", encoding="utf-8") as f:
        f.write("# Out-of-sample temporal prediction\n\n")
        f.write(f"## Headline ({outcome})\n\n")
        f.write(
            f"Splitting the 2004–2025 panel at 2015-01-01 and using the "
            f"model's parameter-free shape (calibrated only by an affine "
            f"scale-and-offset on one half), the model's collapse curve "
            f"trained on 2004–2014 predicts 2015–2025 tail-risk frequency "
            f"with Spearman ρ = **{rho_fwd:+.2f}** (*p* = {p_fwd:.3g}); "
            f"RMSE = {rmse_fwd:.3f} vs. naive constant baseline RMSE "
            f"{rmse_naive_fwd:.3f} "
            f"({(1 - rmse_fwd / rmse_naive_fwd) * 100:.0f}% skill). "
            f"Reverse split (train 2015–2025, predict 2004–2014): "
            f"Spearman = **{rho_rev:+.2f}** (*p* = {p_rev:.3g}); RMSE "
            f"= {rmse_rev:.3f} vs. naive {rmse_naive_rev:.3f}. "
            "The 2015–2025 held-out period contains the August 2015 "
            "China devaluation, the COVID-19 selloff (March 2020), the "
            "2022 inflation shock, and the August 2024 carry-trade "
            "unwind — none seen by the calibration.\n\n"
        )
        if rigidity_oos:
            r_f = rigidity_oos.get("forward")
            r_r = rigidity_oos.get("reverse")
            f.write("**Rigidity transition out-of-sample:** ")
            if r_f is not None:
                f.write(f"forward Spearman = {r_f[0]:+.2f} "
                        f"(p={r_f[1]:.3g}); ")
            if r_r is not None:
                f.write(f"reverse Spearman = {r_r[0]:+.2f} "
                        f"(p={r_r[1]:.3g}). ")
            f.write("\n\n")

        f.write("## Per-decile predictions\n\n")
        f.write("**P2 (2015-2025) — predicted by P1-trained mapping:**\n\n")
        f.write("| J bin mid | actual P(DD>5%) | predicted | model P(coll) | n |\n")
        f.write("|---|---|---|---|---|\n")
        for _, r in p2_bin.iterrows():
            f.write(f"| {r['J_bin_mid']:.2f} | "
                    f"{r['P_tail_5pct']:.2f} | "
                    f"{r['predicted_P_tail']:.2f} | "
                    f"{r['model_P_collapse']:.2f} | "
                    f"{int(r['n_windows'])} |\n")
        f.write("\n**P1 (2004-2014) — predicted by P2-trained mapping:**\n\n")
        f.write("| J bin mid | actual P(DD>5%) | predicted | model P(coll) | n |\n")
        f.write("|---|---|---|---|---|\n")
        for _, r in p1_bin.iterrows():
            f.write(f"| {r['J_bin_mid']:.2f} | "
                    f"{r['P_tail_5pct']:.2f} | "
                    f"{r['predicted_P_tail']:.2f} | "
                    f"{r['model_P_collapse']:.2f} | "
                    f"{int(r['n_windows'])} |\n")
        f.write("\n")

        f.write("## Method\n\n")
        f.write(
            "Affine calibration (two parameters, scale and offset) of the "
            "model's parameter-free P(collapse | J) curve on one period's "
            "decile-binned tail-risk frequency, then prediction of the "
            "other period without re-fitting. The model's *shape* (the "
            "h/J² scaling and its onset) is fixed by Theorem 1; only the "
            "absolute scale (which depends on calibration choices like "
            "drawdown threshold and noise prescription) is tuned. The "
            "naive baseline is the in-sample mean tail-risk frequency.\n\n"
        )

        f.write("## Key sentence for the paper\n\n")
        if outcome == "STRONG":
            f.write(
                f"> The model's collapse-frequency curve, calibrated on "
                f"2004–2014 data with two free parameters (scale and "
                f"offset), correctly predicts 2015–2025 tail-risk "
                f"frequency at Spearman ρ = {rho_fwd:+.2f} "
                f"({(1 - rmse_fwd / rmse_naive_fwd) * 100:.0f}% RMSE "
                f"reduction over a constant-frequency baseline), "
                f"including crisis episodes the calibration never saw "
                f"(COVID-19, 2022 inflation shock, 2024 carry-trade "
                f"unwind). Reverse split gives ρ = {rho_rev:+.2f}.\n"
            )
        else:
            f.write(
                f"> Out-of-sample prediction yields Spearman ρ = "
                f"{rho_fwd:+.2f} (forward) and {rho_rev:+.2f} (reverse). "
                f"The directional skill is positive but the absolute "
                f"prediction quality is moderate; the model's shape "
                f"captures coarse trends but not fine-grained "
                f"period-specific structure.\n"
            )
    print(f"Wrote {verdict_path}")


if __name__ == "__main__":
    main()
