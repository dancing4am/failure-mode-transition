"""
Generate SI Figure S4 — sensitivity sweep variant comparison.

Reads the four CSV outputs from sensitivity_sweep.py and produces a
2×2 figure (one panel per variant) showing P(collapse) heatmap over
the (J, mu) sub-grid. A bottom strip shows rigidity-share for each
variant at mu=100.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SENS = ROOT / "results" / "sensitivity"
OUT = ROOT / "results" / "sensitivity" / "figure_s4_sensitivity_sweep.png"

VARIANTS = [
    ("baseline", "Baseline (multiplicative noise, linear $h(W)$)"),
    ("additive", "Additive noise (linear $h$)"),
    ("log_h",    "Log $h$: $h(W) = c\\log(1 + W/100)$"),
    ("const_h",  "Constant $h = 0.4$"),
]


def load(variant: str) -> pd.DataFrame:
    return pd.read_csv(SENS / f"{variant}_summary.csv")


def heatmap(df: pd.DataFrame, value_col: str) -> tuple[np.ndarray, list, list]:
    pivot = df.pivot(index="mu", columns="J", values=value_col)
    return pivot.values, list(pivot.columns), list(pivot.index)


def main():
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.5),
                             gridspec_kw={"height_ratios": [3, 1.2]})

    # Top row: P(collapse) heatmaps
    for i, (variant, title) in enumerate(VARIANTS):
        df = load(variant)
        Z, j_axis, mu_axis = heatmap(df, "p_collapse")
        ax = axes[0, i]
        im = ax.imshow(Z, origin="lower", aspect="auto", cmap="Reds",
                       vmin=0.0, vmax=0.8,
                       extent=[min(j_axis), max(j_axis),
                               min(mu_axis), max(mu_axis)])
        ax.set_xticks(j_axis)
        ax.set_yticks(mu_axis)
        ax.set_xlabel("coupling $J$")
        if i == 0:
            ax.set_ylabel("economic margin $\\mu$")
        ax.set_title(title, fontsize=10)
        # annotate values
        for jj, J in enumerate(j_axis):
            for mm, mu in enumerate(mu_axis):
                v = Z[mm, jj]
                ax.text(J, mu, f"{v:.2f}", ha="center", va="center",
                        color="white" if v > 0.4 else "black", fontsize=9)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="P(collapse)")

    # Bottom row: rigidity-share bar at mu=100 across J
    for i, (variant, _title) in enumerate(VARIANTS):
        df = load(variant)
        sub = df[df.mu == 100].sort_values("J")
        ax = axes[1, i]
        ax.bar(sub.J.values, sub.rigidity_share.values, width=0.35,
               color="#d62728", alpha=0.85)
        ax.set_xticks(sub.J.values)
        ax.set_xlabel("coupling $J$")
        if i == 0:
            ax.set_ylabel("rigidity share\n(at $\\mu$=100)")
        ax.set_ylim(0, 1.05)
        ax.axhline(0.5, linestyle=":", color="grey", alpha=0.5)
        ax.grid(alpha=0.3, axis="y")

    fig.suptitle(
        "Figure S4. Sensitivity of high-$J$ residual to noise prescription "
        "and $h(W)$ functional form (4,800 runs, 100 seeds per cell).",
        fontsize=11, y=1.00
    )
    fig.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
