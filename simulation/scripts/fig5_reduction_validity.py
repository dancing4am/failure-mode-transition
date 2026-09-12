"""
Generate main-text Figure 4 — reduction-validity comparison.

Reads simulation/results/reduction_validity/reduction_validity.csv and
produces a two-panel figure:

  Panel A: Part 1 N-scan (idiosyncratic noise endpoint). ABM P(collapse)
           vs N (log x-axis) against the one-dimensional (scalar) model
           run at the matched effective noise sigma_eff = xi / sqrt(N).
  Panel B: Part 2 mixed noise at N = 200. P(collapse) vs the common-noise
           share r = sigma_com^2 / xi^2, ABM vs scalar model at matched
           sigma_eff.

Both series carry 95% CI error bars taken directly from the CSV
(ci_lo / ci_hi columns). No figure title is baked in (caption lives in
the manuscript); panels carry only bold lowercase letters a/b (SR
convention) plus a small in-panel regime annotation.

Output: manuscript/arxiv/figures/fig5.png (300 dpi) + fig5.pdf (vector).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "results" / "reduction_validity" / "reduction_validity.csv"
OUT = ROOT.parent / "manuscript" / "arxiv" / "figures" / "fig5.png"

# Repo-standard two-series pairing (cf. network_abm_figures.py,
# sensitivity_figure.py): matplotlib default blue for the ABM,
# default red for the reduced model. Identity is additionally encoded
# by marker shape and linestyle, never by color alone.
STYLE = {
    "abm": dict(color="#1f77b4", marker="o", linestyle="-"),
    "scalar": dict(color="#d62728", marker="s", linestyle="--"),
}


def yerr(sub: pd.DataFrame):
    """Asymmetric error-bar half-widths from absolute CI bounds."""
    lo = (sub.P_collapse - sub.ci_lo).clip(lower=0.0)
    hi = (sub.ci_hi - sub.P_collapse).clip(lower=0.0)
    return [lo.values, hi.values]


def draw_series(ax, sub: pd.DataFrame, x, source: str, label: str):
    st = STYLE[source]
    ax.errorbar(
        x, sub.P_collapse.values, yerr=yerr(sub),
        color=st["color"], marker=st["marker"], linestyle=st["linestyle"],
        linewidth=2, markersize=8, capsize=4, elinewidth=1.4,
        markeredgecolor="white", markeredgewidth=0.8, label=label,
    )


def main():
    df = pd.read_csv(CSV)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 5.4))

    # ------------------------------------------------------------------
    # Panel A — Part 1: N-scan, idiosyncratic-noise endpoint
    # ------------------------------------------------------------------
    p1 = df[df.part == "part1_nscan"]
    for source, label in [
        ("abm", "network ABM ($n=50$ seeds)"),
        ("scalar", "one-dimensional model, matched "
                   r"$\sigma_{\mathrm{eff}}$ ($n=1000$)"),
    ]:
        sub = p1[p1.source == source].sort_values("N")
        # Slight multiplicative x-offset so the two series' error bars
        # do not overprint on the log axis.
        offset = 0.94 if source == "abm" else 1.06
        draw_series(ax_a, sub, sub.N.values * offset, source, label)

    ax_a.set_xscale("log")
    ax_a.set_xticks([1, 5, 25, 200])
    ax_a.set_xticklabels(["1", "5", "25", "200"])
    ax_a.minorticks_off()
    ax_a.set_xlabel("number of agents $N$", fontsize=12)
    ax_a.set_ylabel("P(collapse)", fontsize=12)
    ax_a.set_ylim(-0.02, 0.45)
    ax_a.grid(alpha=0.3)
    ax_a.legend(loc="upper right", fontsize=10, framealpha=0.9)
    ax_a.text(0.03, 0.05,
              "idiosyncratic noise:  "
              r"$\sigma_{\mathrm{eff}} = \xi/\sqrt{N}$",
              transform=ax_a.transAxes, fontsize=11, va="bottom")
    ax_a.text(-0.10, 1.02, "a", transform=ax_a.transAxes,
              fontsize=18, fontweight="bold", va="bottom")

    # ------------------------------------------------------------------
    # Panel B — Part 2: mixed noise at N = 200
    # ------------------------------------------------------------------
    p2 = df[df.part == "part2_mixed"]
    for source, label in [
        ("abm", "network ABM, $N=200$ ($n=100$ seeds)"),
        ("scalar", "one-dimensional model, matched "
                   r"$\sigma_{\mathrm{eff}}$ ($n=1000$)"),
    ]:
        sub = p2[p2.source == source].sort_values("r")
        offset = -0.008 if source == "abm" else 0.008
        draw_series(ax_b, sub, sub.r.values + offset, source, label)

    ax_b.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_b.set_xlabel(r"common-noise share $r = \sigma_{\mathrm{com}}^2/\xi^2$",
                    fontsize=12)
    ax_b.set_ylabel("P(collapse)", fontsize=12)
    ax_b.set_xlim(-0.07, 1.07)
    ax_b.set_ylim(-0.02, 0.45)
    ax_b.grid(alpha=0.3)
    ax_b.legend(loc="upper left", fontsize=10, framealpha=0.9)
    ax_b.text(0.97, 0.05,
              "mixed noise at $N=200$:  "
              r"$\sigma_{\mathrm{eff}}^2 = \sigma_{\mathrm{com}}^2"
              r" + \sigma_{\mathrm{idio}}^2/N$",
              transform=ax_b.transAxes, fontsize=11,
              va="bottom", ha="right")
    ax_b.text(-0.10, 1.02, "b", transform=ax_b.transAxes,
              fontsize=18, fontweight="bold", va="bottom")

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    for out in (OUT, OUT.with_suffix(".pdf")):
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"Wrote {out}")
    plt.close(fig)

    # Echo the plotted values for verification against the CSV.
    for part, ax_name in [("part1_nscan", "a"), ("part2_mixed", "b")]:
        print(f"\nPanel {ax_name} ({part}):")
        sub = df[df.part == part].sort_values(["source", "N", "r"])
        for _, row in sub.iterrows():
            xcol = "N" if part == "part1_nscan" else "r"
            print(f"  {row.source:6s} {xcol}={row[xcol]:<6g} "
                  f"P={row.P_collapse:.3f} "
                  f"CI=[{row.ci_lo:.3f}, {row.ci_hi:.3f}] n={row.n}")


if __name__ == "__main__":
    main()
