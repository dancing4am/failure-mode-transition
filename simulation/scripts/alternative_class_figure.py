"""
Generate SI Figure S5 — alternative-model-class comparison.

Reads the four CSVs produced by alternative_class_sweep.py and produces
a 2x4 figure: P(collapse) heatmaps (top) + rigidity-share bar charts (bottom).
"""

from __future__ import annotations

def _txt_color(cmap_name, v):
    """Annotation colour by cell-colour luminance (not by value threshold)."""
    import matplotlib.pyplot as _plt
    r, g, b, _ = plt.get_cmap(cmap_name)(v)
    lin = [(c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4)
           for c in (r, g, b)]
    L = 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    # pick whichever of black/white text maximizes WCAG contrast; the
    # worst case of this rule is sqrt(21) ~ 4.58:1, above the 4.5:1 floor
    return "black" if (L + 0.05) / 0.05 >= (1.05) / (L + 0.05) else "white"


from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "results" / "alternative_class"
OUT = DIR / "figure_s5_alternative_class.png"

VARIANTS = [
    ("curie_weiss", "Curie–Weiss\n$-m + \\tanh((Jm+h)/T)$"),
    ("voter_like",  "Voter-like\n$-m + \\mathrm{sign}(\\cdot)\\min(1,|\\cdot|)$"),
    ("kuramoto_1d", "Kuramoto-1D\n$-\\sin(m\\pi/2) + \\tanh(\\cdot)$"),
    ("cubic",       "Cubic Landau\n$-(m - m^3/3) + \\tanh(\\cdot)$"),
]


def main():
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.5),
                             gridspec_kw={"height_ratios": [3, 1.2]})

    for i, (variant, title) in enumerate(VARIANTS):
        df = pd.read_csv(DIR / f"{variant}_summary.csv")
        pivot = df.pivot(index="mu", columns="J", values="p_collapse")
        Z = pivot.values
        j_axis = list(pivot.columns)
        mu_axis = list(pivot.index)

        ax = axes[0, i]
        # half-cell padding: centers cells on the tick values and keeps the
        # edge-column annotations (incl. J = 5.0) inside the axes
        dj = (max(j_axis) - min(j_axis)) / max(len(j_axis) - 1, 1) / 2
        dm = (max(mu_axis) - min(mu_axis)) / max(len(mu_axis) - 1, 1) / 2
        im = ax.imshow(Z, origin="lower", aspect="auto", cmap="Reds",
                       vmin=0.0, vmax=0.8,
                       extent=[min(j_axis) - dj, max(j_axis) + dj,
                               min(mu_axis) - dm, max(mu_axis) + dm])
        ax.set_xticks(j_axis)
        ax.set_yticks(mu_axis)
        ax.set_xlabel("coupling $J$")
        if i == 0:
            ax.set_ylabel("economic margin $\\mu$")
        ax.set_title(title, fontsize=10)
        for jj, J in enumerate(j_axis):
            for mm, mu in enumerate(mu_axis):
                v = Z[mm, jj]
                ax.text(J, mu, f"{v:.2f}", ha="center", va="center",
                        color=_txt_color("Reds", v), fontsize=9)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="P(collapse)")

    for i, (variant, _t) in enumerate(VARIANTS):
        df = pd.read_csv(DIR / f"{variant}_summary.csv")
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
        "The high-$J$ residual asymmetry is not Curie–Weiss-specific. "
        "Four bistable mean-field SDEs with field bias all show the same pattern "
        "(4,800 runs, 100 seeds per cell).",
        fontsize=10, y=1.00
    )
    fig.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
