"""
Generate Figure 8: heuristic AI coupling proxy on the model's
prediction surface. Marks the AlpacaEval output-similarity proxy
(J_eff = 3.92, rho = 0.797) as the primary reference, with the
adversarial-transfer proxy (J = 1.18) and benchmark-vector proxy
(J = 8.7) as secondary background reference lines.

Reads combined_sweep_summary.csv produced by the extended J-sweep.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = (ROOT / "results" / "ai_coupling_overlay"
           / "combined_sweep_summary.csv")
OUT = ROOT / "results" / "ai_coupling_overlay"

AI_DIRECT_J = 3.92         # direct measurement (AlpacaEval, 14 models; rho = 0.797 → J = 0.797/0.203 = 3.92)
AI_DIRECT_RHO = 0.797      # mean off-diagonal cosine similarity (precise; 0.80 rounded)
AI_ADV_J = 1.18            # adversarial-transfer proxy
AI_BEN_J = 8.7             # benchmark-vector proxy
J_CRITICAL = 2.0
FIN_AVG_J = 1.0
FIN_CRISIS_J = 9.0


def curve(df, mu):
    sub = df[df.mu == mu].sort_values("J")
    return sub.J.values, sub.p_collapse.values


def interp(xs, ys, q):
    return float(np.interp(q, xs, ys))


def main():
    df = pd.read_csv(SUMMARY)
    J100, P100 = curve(df, 100)
    J60, P60 = curve(df, 60)
    J40, P40 = curve(df, 40)
    J20, P20 = curve(df, 20)

    eff = np.zeros_like(P100)
    for i in range(len(eff)):
        if P20[i] <= 0:
            eff[i] = 1.0 if P100[i] <= 0 else 0.0
        else:
            eff[i] = max(0.0, 1.0 - P100[i] / P20[i])

    p_direct = interp(J100, P100, AI_DIRECT_J)
    eff_direct = interp(J100, eff, AI_DIRECT_J)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8))

    def reference_lines(ax, ymax_label_y):
        # Primary: AlpacaEval output-similarity proxy (heuristic, not calibrated)
        ax.axvline(AI_DIRECT_J, color="#d62728", linestyle="--",
                   linewidth=2.5,
                   label=f"AlpacaEval proxy (heuristic): $J$ = {AI_DIRECT_J} "
                         f"($\\rho$ = {AI_DIRECT_RHO})")
        # Secondary proxies
        ax.axvline(AI_ADV_J, color="#2ca02c", linestyle=":",
                   linewidth=1.5, alpha=0.75,
                   label=f"Adversarial-transfer proxy: $J$ = {AI_ADV_J}")
        ax.axvline(AI_BEN_J, color="#9467bd", linestyle=":",
                   linewidth=1.5, alpha=0.75,
                   label=f"Benchmark-vector proxy: $J$ = {AI_BEN_J}")
        # Critical threshold
        ax.axvline(J_CRITICAL, color="black", linestyle="--",
                   linewidth=1.5,
                   label=f"$J_c = T$ = {J_CRITICAL}")
        # Financial reference
        ax.axvline(FIN_AVG_J, color="grey", linestyle=":",
                   linewidth=1, alpha=0.5,
                   label=f"$J_{{\\rm fin,avg}}$ = {FIN_AVG_J}")
        ax.axvline(FIN_CRISIS_J, color="grey", linestyle="-.",
                   linewidth=1, alpha=0.5,
                   label=f"$J_{{\\rm fin,crisis}}$ = {FIN_CRISIS_J}")
        # Shaded regimes
        ax.axvspan(J_CRITICAL, 11.0, alpha=0.06, color="red")
        ax.axvspan(0, J_CRITICAL / 2, alpha=0.06, color="green")

    # Panel A: P(collapse) vs J
    axA = axes[0]
    axA.plot(J40, P40, marker="o", color="#d62728", linestyle="--",
             linewidth=2, label=r"$\mu$ = 40 (unviable)")
    axA.plot(J60, P60, marker="s", color="#ff7f0e", linestyle="--",
             linewidth=2, label=r"$\mu$ = 60 (resilient)")
    axA.plot(J100, P100, marker="^", color="#1f77b4", linestyle="-",
             linewidth=2.5, label=r"$\mu$ = 100 (max stabilizer)")
    reference_lines(axA, 1.05)
    axA.scatter([AI_DIRECT_J], [p_direct], color="#d62728", s=180,
                zorder=6, edgecolor="black", linewidth=1.5)
    axA.annotate(f"P(coll) = {p_direct:.2f}",
                 xy=(AI_DIRECT_J, p_direct), xytext=(10, 12),
                 textcoords="offset points",
                 fontsize=10, color="#d62728", fontweight="bold")
    axA.set_xlabel("coupling $J$", fontsize=12)
    axA.set_ylabel("P(collapse)", fontsize=12)
    axA.set_xlim(0, 10.5)
    axA.set_ylim(-0.03, 1.05)
    axA.grid(alpha=0.3)
    axA.set_title("A.  Model P(collapse) at the AlpacaEval coupling proxy",
                  fontsize=11)
    axA.legend(loc="upper right", fontsize=8, ncol=1)

    # Panel B: passive effectiveness
    axB = axes[1]
    axB.plot(J100, eff, marker="^", color="#1f77b4", linestyle="-",
             linewidth=2.5,
             label=r"empirical effectiveness ($\mu$=100 vs $\mu$=20)")
    h_J2 = 1.0 / (J100 ** 2)
    h_J2 = h_J2 * (eff.max() / h_J2.max())
    axB.plot(J100, h_J2, color="#9467bd", linestyle=":",
             linewidth=2, alpha=0.85,
             label=r"theoretical $h/J^2$ envelope (normalized)")
    reference_lines(axB, 1.05)
    axB.scatter([AI_DIRECT_J], [eff_direct], color="#d62728", s=180,
                zorder=6, edgecolor="black", linewidth=1.5)
    axB.annotate(f"{eff_direct * 100:.0f}%",
                 xy=(AI_DIRECT_J, eff_direct), xytext=(10, 12),
                 textcoords="offset points",
                 fontsize=10, color="#d62728", fontweight="bold")
    axB.set_xlabel("coupling $J$", fontsize=12)
    axB.set_ylabel("passive stabilizer effectiveness", fontsize=11)
    axB.set_xlim(0, 10.5)
    axB.set_ylim(-0.03, 1.05)
    axB.grid(alpha=0.3)
    axB.set_title("B.  Passive effectiveness at the AlpacaEval coupling proxy",
                  fontsize=11)
    axB.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        f"Figure 8.  Heuristic AI coupling proxy on the model's "
        f"prediction surface (AlpacaEval, 14 models × 805 prompts; "
        f"$J_{{\\rm eff}}$ = {AI_DIRECT_J}, $\\rho$ = {AI_DIRECT_RHO})",
        fontsize=12, y=1.02
    )
    fig.tight_layout()

    out = OUT / "ai_coupling_overlay_figure.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")
    print(f"  At J = {AI_DIRECT_J}: P(collapse) = {p_direct:.2f}, "
          f"effectiveness = {eff_direct * 100:.0f}%")


if __name__ == "__main__":
    main()
