"""Render manuscript Figure 1 (static asymmetry in the minimal model) from the
deposited sweep summary `simulation/results/minimal_model/cell_summary.csv`.

Replaces the archived rendering with an identical-content re-render whose
lower-row annotations are print-legible and whose in-cell text colour is chosen
by the cell colour's luminance rather than a fixed value threshold (a value
threshold puts white text on the pale mid-scale colours of a value-thresholded map).

Layout matches the archived composite: top panel P(collapse) over the (J, mu)
grid; lower row P(rigidity collapse) and P(fragmentation collapse).
Outputs: manuscript/arxiv/figures/fig1.png (300 dpi) and fig1.pdf (vector).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT.parent / "manuscript" / "arxiv" / "figures"
SUMMARY = ROOT / "results" / "minimal_model" / "cell_summary.csv"


def _txt_color(cmap_name, v):
    r, g, b, _ = plt.get_cmap(cmap_name)(v)
    lin = [(c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4)
           for c in (r, g, b)]
    L = 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    # pick whichever of black/white text maximizes WCAG contrast; the
    # worst case of this rule is sqrt(21) ~ 4.58:1, above the 4.5:1 floor
    return "black" if (L + 0.05) / 0.05 >= (1.05) / (L + 0.05) else "white"


def _heat(ax, piv, cmap, ann_fs):
    im = ax.imshow(piv.values, aspect="auto", cmap=cmap, vmin=0.0, vmax=1.0)
    ax.set_xticks(range(len(piv.columns)))
    ax.set_xticklabels(piv.columns, fontsize=14)
    ax.set_yticks(range(len(piv.index)))
    ax.set_yticklabels([f"{j:.1f}" for j in piv.index], fontsize=14)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=ann_fs, color=_txt_color(cmap, v))
    return im


def main():
    df = pd.read_csv(SUMMARY)

    def piv(field):
        p = df.pivot(index="J", columns="mult", values=field)
        return p.reindex(index=sorted(p.index, reverse=True),
                         columns=sorted(p.columns))

    fig = plt.figure(figsize=(15, 17.25))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.75, 1.0],
                          hspace=0.28, wspace=0.25,
                          left=0.07, right=0.90, top=0.95, bottom=0.06)

    # --- top: P(collapse) over the full grid -----------------------------
    ax0 = fig.add_subplot(gs[0, :])
    im0 = _heat(ax0, piv("p_collapse"), "cividis", ann_fs=14)
    ax0.set_title("Minimal model: P(collapse) over 100 seeds × 20,000 steps "
                  "(dt=0.05)", fontsize=17)
    ax0.set_xlabel("income multiplier  (passive stabilizer)", fontsize=15)
    ax0.set_ylabel("coupling J", fontsize=15)
    cb0 = fig.colorbar(im0, ax=ax0, fraction=0.04, pad=0.02)
    cb0.set_label("P(collapse)", fontsize=14)
    cb0.ax.tick_params(labelsize=13)

    # --- lower row: collapse-type shares ---------------------------------
    fig.text(0.485, 0.385,
             "Collapse type share: rigidity dominates at high J; "
             "fragmentation at low J", ha="center", fontsize=15)
    for col, (field, title) in enumerate(
            [("p_rigidity", "P(rigidity collapse)\n|m|>0.9 at collapse"),
             ("p_fragmentation",
              "P(fragmentation collapse)\n|m|<0.3 at collapse")]):
        ax = fig.add_subplot(gs[1, col])
        im = _heat(ax, piv(field), "Reds", ann_fs=13)
        ax.set_title(title, fontsize=15)
        ax.set_xlabel("income multiplier", fontsize=14)
        if col == 0:
            ax.set_ylabel("coupling J", fontsize=14)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=12)

    FIGDIR.mkdir(parents=True, exist_ok=True)
    for out in (FIGDIR / "fig1.png", FIGDIR / "fig1.pdf"):
        fig.savefig(out, dpi=300)
        print(f"Wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
