"""
Convexity / nonlinearity test of empirical tail-risk vs J.

The model predicts that P(collapse) is convex in J at intermediate
coupling — accelerating failure as coupling rises. A linear correlation-
risk model predicts constant slope. If the empirical tail-risk curve is
convex, the model captures the acceleration that linear models miss.
This also helps defend the exponent mismatch: even if the absolute
slope is steeper in the model than in the data, both curves *accelerate*
in J, which is itself a non-trivial prediction.

Inputs:
  - empirical/results/tail_risk_by_coupling.csv (10 deciles, P_tail_5pct)

Outputs:
  - empirical/results/convexity_test.csv  (model-fit comparison)
  - empirical/figures/convexity_comparison.png  (data + 3 fits)
  - empirical/CONVEXITY_VERDICT.md
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EMP = Path(__file__).resolve().parents[1]
FIG = EMP / "figures"
RESULTS = EMP / "results"
TAIL_CSV = RESULTS / "tail_risk_by_coupling.csv"
N_BOOTS = 1000
SEED = 0xCAFE


def fit_linear(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    rss = float(((y - yhat) ** 2).sum())
    return coef, yhat, rss, len(coef)


def fit_quadratic(x, y):
    A = np.vstack([x ** 2, x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    rss = float(((y - yhat) ** 2).sum())
    return coef, yhat, rss, len(coef)


def fit_powerlaw(x, y, eps=1e-6):
    # log y = a + b log x  on positive y values only
    mask = (y > eps) & (x > 0)
    if mask.sum() < 3:
        return None
    lx = np.log(x[mask])
    ly = np.log(y[mask])
    A = np.vstack([lx, np.ones_like(lx)]).T
    coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
    return float(coef[0]), float(coef[1]), int(mask.sum())


def aic(rss, n, k):
    """Simple Gaussian-approximation AIC for least-squares fits."""
    if rss <= 0:
        return -np.inf
    return 2 * k + n * np.log(rss / n)


def adjusted_r2(rss, y, k):
    n = len(y)
    tss = float(((y - y.mean()) ** 2).sum())
    if tss <= 0:
        return float("nan")
    r2 = 1.0 - rss / tss
    if n - k - 1 <= 0:
        return r2
    return 1.0 - (1.0 - r2) * (n - 1) / (n - k - 1)


def main():
    df = pd.read_csv(TAIL_CSV).sort_values("J_bin_mid").reset_index(drop=True)
    x = df["J_bin_mid"].to_numpy(dtype=float)
    y = df["P_tail_5pct"].to_numpy(dtype=float)
    n = len(x)
    print(f"Loaded {n} bins from {TAIL_CSV.name}")

    # Linear fit
    coef_lin, yhat_lin, rss_lin, k_lin = fit_linear(x, y)
    aic_lin = aic(rss_lin, n, k_lin + 1)  # +1 for sigma
    r2adj_lin = adjusted_r2(rss_lin, y, k_lin)

    # Quadratic fit
    coef_q, yhat_q, rss_q, k_q = fit_quadratic(x, y)
    aic_q = aic(rss_q, n, k_q + 1)
    r2adj_q = adjusted_r2(rss_q, y, k_q)
    aic_diff = aic_lin - aic_q  # positive means quadratic is better

    # Power-law fit
    pl = fit_powerlaw(x, y)
    if pl is not None:
        slope_pl, intercept_pl, n_pl = pl
    else:
        slope_pl = float("nan")
        intercept_pl = float("nan")
        n_pl = 0

    # Piecewise slope analysis: low half (deciles 0-4) vs high half (5-9)
    half = n // 2
    coef_low, _, _, _ = fit_linear(x[:half], y[:half])
    coef_high, _, _, _ = fit_linear(x[half:], y[half:])
    slope_low = float(coef_low[0])
    slope_high = float(coef_high[0])
    slope_ratio = slope_high / slope_low if slope_low != 0 else float("nan")

    # Second-derivative estimate via central finite differences
    d2 = np.full(n, np.nan)
    for i in range(1, n - 1):
        h1 = x[i] - x[i - 1]
        h2 = x[i + 1] - x[i]
        denom = 0.5 * (h1 + h2)
        if denom > 0:
            d2[i] = ((y[i + 1] - y[i]) / h2 - (y[i] - y[i - 1]) / h1) / denom
    d2_mean = float(np.nanmean(d2))
    d2_pos_frac = float(np.mean(d2[~np.isnan(d2)] > 0))

    # Bootstrap confidence interval on the quadratic coefficient c
    rng = np.random.default_rng(SEED)
    boot_c = np.empty(N_BOOTS)
    for i in range(N_BOOTS):
        idx = rng.integers(0, n, n)
        try:
            cb, *_ = fit_quadratic(x[idx], y[idx])
            boot_c[i] = cb[0]
        except Exception:
            boot_c[i] = np.nan
    boot_c = boot_c[~np.isnan(boot_c)]
    ci_lo = float(np.percentile(boot_c, 2.5))
    ci_hi = float(np.percentile(boot_c, 97.5))
    c_quad = float(coef_q[0])

    # Save results
    rows = [
        {"model": "linear", "params": f"a={coef_lin[0]:.3f}, b={coef_lin[1]:.3f}",
         "rss": rss_lin, "aic": aic_lin, "adj_R2": r2adj_lin, "k": k_lin},
        {"model": "quadratic",
         "params": f"a={coef_q[0]:.3f} (J^2), b={coef_q[1]:.3f} (J), c={coef_q[2]:.3f}",
         "rss": rss_q, "aic": aic_q, "adj_R2": r2adj_q, "k": k_q},
        {"model": "power law",
         "params": (f"slope (log-log) = {slope_pl:.3f}, "
                    f"intercept = {intercept_pl:.3f}, n_used = {n_pl}"),
         "rss": float("nan"), "aic": float("nan"),
         "adj_R2": float("nan"), "k": 2},
    ]
    out_df = pd.DataFrame(rows)
    out_df.to_csv(RESULTS / "convexity_test.csv", index=False)
    print()
    print(out_df.to_string(index=False))

    print()
    print(f"AIC(linear) - AIC(quadratic) = {aic_diff:+.2f}  "
          f"(positive = quadratic better)")
    print(f"adj R² linear / quadratic = {r2adj_lin:.3f} / {r2adj_q:.3f}")
    print(f"Quadratic coefficient c (J^2): {c_quad:+.4f}  "
          f"95% CI [{ci_lo:+.4f}, {ci_hi:+.4f}]  -> "
          f"{'CONVEX' if c_quad > 0 and ci_lo > 0 else 'inconclusive'}")
    print(f"Piecewise slopes: low half = {slope_low:.3f}, "
          f"high half = {slope_high:.3f}, ratio = {slope_ratio:.2f}")
    print(f"Second-derivative mean (finite differences) = {d2_mean:+.4f}, "
          f"fraction positive = {d2_pos_frac:.2f}")
    print(f"Power-law log-log slope = {slope_pl:.3f}  "
          f"(model log-log slope ~ 2.51)")

    # Figure: data + 3 fits
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.scatter(x, y, color="#d62728", s=85, label="empirical P(DD>5% | J)",
               zorder=4, edgecolor="black", linewidth=0.5)
    xline = np.linspace(x.min(), x.max(), 200)
    yhat_lin_line = coef_lin[0] * xline + coef_lin[1]
    yhat_q_line = coef_q[0] * xline ** 2 + coef_q[1] * xline + coef_q[2]
    ax.plot(xline, yhat_lin_line, "--", color="#1f77b4", linewidth=2,
            label=f"linear: slope = {coef_lin[0]:.3f}, "
                  f"adj $R^2$ = {r2adj_lin:.2f}")
    ax.plot(xline, yhat_q_line, "-", color="#2ca02c", linewidth=2,
            label=f"quadratic: $c$(J²) = {c_quad:+.3f} "
                  f"[{ci_lo:+.3f}, {ci_hi:+.3f}], adj $R^2$ = "
                  f"{r2adj_q:.2f}")
    if pl is not None:
        # Power-law on log axes only — overlay equivalent on linear axes
        yhat_pl_line = np.exp(intercept_pl) * xline ** slope_pl
        ax.plot(xline, yhat_pl_line, ":", color="#9467bd", linewidth=2,
                label=f"power law: slope = {slope_pl:.2f} "
                      f"(model ≈ 2.51)")
    ax.set_xlabel("coupling $J = \\bar\\rho/(1-\\bar\\rho)$",
                  fontsize=11)
    ax.set_ylabel("P(max drawdown > 5% | J bin)", fontsize=11)
    ax.set_title("Empirical tail-risk vs coupling: testing for convexity\n"
                 f"AIC(linear) - AIC(quadratic) = {aic_diff:+.2f}; "
                 f"piecewise slope ratio (high/low) = {slope_ratio:.2f}",
                 fontsize=11)
    ax.set_ylim(-0.03, 1.05)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "convexity_comparison.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    print(f"\nWrote {FIG / 'convexity_comparison.png'}")

    # ----- VERDICT -----
    convex = (c_quad > 0) and (ci_lo > 0)
    accelerating = slope_ratio > 1.5
    if convex and accelerating and aic_diff > 2:
        outcome = "STRONG"
    elif convex or accelerating:
        outcome = "MODERATE"
    else:
        outcome = "WEAK"

    verdict_path = EMP / "CONVEXITY_VERDICT.md"
    with verdict_path.open("w", encoding="utf-8") as f:
        f.write("# Convexity / nonlinearity test of P(tail | J)\n\n")
        f.write(f"## Headline ({outcome})\n\n")
        f.write(
            f"The empirical tail-risk curve P(max DD > 5% | J) over "
            f"{n} coupling deciles is "
            f"**{'convex' if convex else 'not clearly convex'}** in J. "
            f"Quadratic-fit J² coefficient *c* = {c_quad:+.3f} "
            f"(95% bootstrap CI [{ci_lo:+.3f}, {ci_hi:+.3f}]); "
            f"AIC(linear) − AIC(quadratic) = {aic_diff:+.2f} "
            f"({'quadratic preferred' if aic_diff > 2 else 'preferences mixed' if aic_diff > 0 else 'linear preferred'}). "
            f"Piecewise slope ratio (high-J half ÷ low-J half) "
            f"= **{slope_ratio:.2f}** "
            f"(linear-model expectation = 1.0; "
            f"{'accelerating' if slope_ratio > 1.5 else 'roughly linear'}). "
            f"Power-law log-log slope = {slope_pl:.2f} "
            f"(minimal-model log-log slope ≈ 2.51). "
            f"Mean finite-difference second derivative = {d2_mean:+.4f} "
            f"(positive → convex); "
            f"{int(d2_pos_frac * 100)}% of internal points have positive "
            f"second derivative.\n\n"
        )

        f.write("## Model fits\n\n")
        f.write("| model | parameters | RSS | adj R² | AIC |\n")
        f.write("|---|---|---|---|---|\n")
        f.write(f"| linear | a={coef_lin[0]:.3f}, b={coef_lin[1]:.3f} | "
                f"{rss_lin:.3f} | {r2adj_lin:.3f} | "
                f"{aic_lin:.2f} |\n")
        f.write(f"| quadratic | "
                f"c·J² + b·J + a, c={c_quad:+.3f} | "
                f"{rss_q:.3f} | {r2adj_q:.3f} | "
                f"{aic_q:.2f} |\n")
        f.write(f"| power law | "
                f"log y = {slope_pl:.2f} log J + {intercept_pl:.2f} "
                f"(n={n_pl}) | — | — | — |\n\n")

        f.write("## Method\n\n")
        f.write(
            "Three regressions on the 10 decile-binned data points "
            "(J_bin_mid, P_tail_5pct) from the tail-risk analysis: "
            "linear (P = a + bJ), quadratic (P = a + bJ + cJ²), and "
            "power-law (log P = a + b log J on positive points). AIC "
            "uses the Gaussian-likelihood approximation; quadratic is "
            "preferred over linear if AIC(linear) - AIC(quadratic) > 2. "
            "Piecewise slopes are computed by linear regression on the "
            "lower-half (5 deciles) and upper-half (5 deciles) "
            "separately. Bootstrap CI on the quadratic coefficient *c* "
            f"uses {N_BOOTS} resamples of the 10 deciles.\n\n"
        )

        f.write("## Key sentence for the paper\n\n")
        if outcome == "STRONG":
            f.write(
                f"> Empirical tail-risk frequency is convex in coupling: "
                f"the J²-term coefficient is {c_quad:+.3f} (95% CI "
                f"[{ci_lo:+.3f}, {ci_hi:+.3f}]; quadratic AIC improvement "
                f"of {aic_diff:.1f} over linear), and the regression "
                f"slope in the high-coupling half is {slope_ratio:.1f}× "
                f"the slope in the low-coupling half. Risk does not just "
                f"increase with coupling — it accelerates, consistent "
                f"with the model's predicted convex P(collapse | J) "
                f"surface and inconsistent with linear "
                f"correlation-risk models.\n"
            )
        elif outcome == "MODERATE":
            f.write(
                f"> The empirical tail-risk curve shows accelerating "
                f"behavior at high coupling (high-J slope = "
                f"{slope_ratio:.1f}× low-J slope; quadratic-fit J² "
                f"coefficient {c_quad:+.3f}, 95% CI [{ci_lo:+.3f}, "
                f"{ci_hi:+.3f}]). The convexity is directionally "
                f"consistent with the model's predicted convex "
                f"P(collapse | J) surface, with the caveat that the "
                f"signal is partial across the alternative tests.\n"
            )
        else:
            f.write(
                f"> The empirical tail-risk curve is approximately "
                f"linear in J (high/low slope ratio = {slope_ratio:.1f}; "
                f"quadratic improvement over linear is small). The "
                f"model's predicted convexity is not clearly visible at "
                f"this level of binning.\n"
            )

        f.write("\n## Recommendation\n\n")
        if outcome == "STRONG":
            f.write(
                "Add one sentence to the tail-risk Results subsection or "
                "as a separate Convexity subsection. Strengthens the "
                "exponent-mismatch defense: even though the model's "
                "log-log slope (2.51) differs from the empirical slope "
                f"({slope_pl:.2f}), the empirical curve is still convex "
                "in J — accelerating failure that linear correlation-risk "
                "models cannot reproduce.\n"
            )
        elif outcome == "MODERATE":
            f.write(
                "Add one cautious sentence to the exponent-mismatch "
                "Limitations paragraph noting that the empirical curve "
                "shows accelerating behavior, even if the absolute "
                "exponent differs from the model's.\n"
            )
        else:
            f.write(
                "Skip integration. The convexity signal is not strong "
                "enough to support the \"accelerating failure\" framing.\n"
            )
    print(f"Wrote {verdict_path}")


if __name__ == "__main__":
    main()
