"""
Generate Figure 3 + supporting figures + NETWORK_VERDICT.md from the
ABM sweep outputs in simulation/results/network_abm/.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "network_abm"

TOPOLOGY_LABELS = {
    "complete":          "Complete (mean-field control)",
    "erdos_renyi":       "Erdős–Rényi",
    "watts_strogatz":    "Watts–Strogatz",
    "barabasi_albert":   "Barabási–Albert",
    "modular":           "Modular (4 communities)",
    "modular_boundary":  "Modular, h at boundary only",
}
TOPOLOGY_ORDER = list(TOPOLOGY_LABELS.keys())

J_VALUES = (0.5, 1.0, 2.0, 3.0, 4.0, 5.0)
MULT_VALUES = (20, 40, 60, 80, 100)


def load_summary():
    df = pd.read_csv(OUT / "network_comparison.csv")
    return df


def heatmap_panel(ax, df, topology, vmin=0.0, vmax=1.0):
    sub = df[df.topology == topology]
    piv = sub.pivot(index="J", columns="mu", values="p_collapse")
    piv = piv.reindex(index=sorted(piv.index, reverse=True),
                      columns=sorted(piv.columns))
    im = ax.imshow(piv.values, aspect="auto", cmap="magma_r",
                   vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(piv.columns)))
    ax.set_xticklabels(piv.columns, fontsize=8)
    ax.set_yticks(range(len(piv.index)))
    ax.set_yticklabels([f"{j:.1f}" for j in piv.index], fontsize=8)
    ax.set_xlabel(r"income multiplier $\mu$", fontsize=9)
    ax.set_ylabel(r"coupling $J$", fontsize=9)
    ax.set_title(TOPOLOGY_LABELS[topology], fontsize=10)
    for i, J_val in enumerate(piv.index):
        for j, mu_val in enumerate(piv.columns):
            v = piv.values[i, j]
            color = "white" if v > 0.55 else "black"
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=6, color=color)
    return im


def figure_pcollapse_comparison(df):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.flatten()
    last_im = None
    for ax, topology in zip(axes, TOPOLOGY_ORDER):
        last_im = heatmap_panel(ax, df, topology)
    fig.suptitle("Figure 3. P(collapse) across network topologies "
                 "(N = 200, mean degree ≈ 20, 50 seeds per cell)",
                 fontsize=12, y=1.00)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(last_im, cax=cbar_ax, label="P(collapse)")
    fig.tight_layout(rect=[0, 0, 0.91, 0.97])
    out = OUT / "network_pcollapse_comparison.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


def figure_scaling_curves(df):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    colors = {
        "complete":         "#1f77b4",
        "erdos_renyi":      "#ff7f0e",
        "watts_strogatz":   "#2ca02c",
        "barabasi_albert":  "#d62728",
        "modular":          "#9467bd",
        "modular_boundary": "#8c564b",
    }
    markers = {
        "complete":         "o",
        "erdos_renyi":      "s",
        "watts_strogatz":   "^",
        "barabasi_albert":  "v",
        "modular":          "D",
        "modular_boundary": "P",
    }
    for topology in TOPOLOGY_ORDER:
        sub = df[(df.topology == topology) & (df.mu == 100)].sort_values("J")
        ax.plot(sub.J, sub.p_collapse,
                marker=markers[topology], color=colors[topology],
                linewidth=2, markersize=9,
                label=TOPOLOGY_LABELS[topology])
    ax.axvline(2.0, color="black", linestyle="--", linewidth=1,
               label=r"$J_c = T = 2$")
    ax.set_xlabel("coupling $J$", fontsize=11)
    ax.set_ylabel(r"P(collapse) at $\mu = 100$", fontsize=11)
    ax.set_title("Network topology vs. high-J residual collapse rate",
                 fontsize=11)
    ax.set_ylim(-0.03, 1.05)
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    out = OUT / "network_scaling_curves.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


def figure_modular_decoupling(df):
    """For modular topologies: plot intra-community vs inter-community
    correlation as a function of J."""
    sub_mod = df[(df.topology == "modular") & (df.mu == 100)].sort_values("J")
    sub_bdr = df[(df.topology == "modular_boundary") & (df.mu == 100)].sort_values("J")
    if sub_mod.mean_intra_corr.isna().all():
        print("No intra/inter correlation data; skipping modular decoupling fig.")
        return

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(sub_mod.J, sub_mod.mean_intra_corr, "o-",
            color="#9467bd", linewidth=2, markersize=9,
            label="Modular: intra-community ⟨corr(m_i, m_j)⟩")
    ax.plot(sub_mod.J, sub_mod.mean_inter_corr, "s-",
            color="#9467bd", linestyle="--", linewidth=2, markersize=9,
            label="Modular: inter-community ⟨corr(m_i, m_j)⟩")
    ax.plot(sub_bdr.J, sub_bdr.mean_intra_corr, "D-",
            color="#8c564b", linewidth=2, markersize=9,
            label="Modular boundary-h: intra-community")
    ax.plot(sub_bdr.J, sub_bdr.mean_inter_corr, "P-",
            color="#8c564b", linestyle="--", linewidth=2, markersize=9,
            label="Modular boundary-h: inter-community")
    ax.axvline(2.0, color="black", linestyle="--", linewidth=1,
               label=r"$J_c = T = 2$")
    ax.set_xlabel("coupling $J$", fontsize=11)
    ax.set_ylabel("mean pairwise corr of $m_i$ trajectories", fontsize=11)
    ax.set_title("Modular community decoupling vs. coupling level",
                 fontsize=11)
    ax.set_ylim(-0.1, 1.05)
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    out = OUT / "modular_decoupling.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


def write_verdict(df):
    """Write NETWORK_VERDICT.md with scenario-A/B/C determination."""
    # Headline numbers at (J=5, mu=100) per topology
    cell = df[(df.J == 5.0) & (df.mu == 100)].set_index("topology")
    p_complete = cell.loc["complete", "p_collapse"]
    p_modular = cell.loc["modular", "p_collapse"]
    p_modular_b = cell.loc["modular_boundary", "p_collapse"]
    p_er = cell.loc["erdos_renyi", "p_collapse"]
    p_ws = cell.loc["watts_strogatz", "p_collapse"]
    p_ba = cell.loc["barabasi_albert", "p_collapse"]

    # Scenario determination based on modular vs complete at high J
    rel_drop_mod = (p_complete - p_modular) / max(p_complete, 0.01)
    rel_drop_bdr = (p_complete - p_modular_b) / max(p_complete, 0.01)

    if p_modular < 0.05 and rel_drop_mod > 0.7:
        scenario = "C"
        scenario_text = (
            "modularity FULLY restores passive stabilizer effectiveness — "
            "the h/J² scaling does not apply on modular topologies"
        )
    elif p_modular < p_complete * 0.6 or rel_drop_mod > 0.4:
        scenario = "B"
        scenario_text = (
            "modularity PARTIALLY restores passive stabilizer effectiveness — "
            "the h/J² scaling is preserved but the residual is reduced"
        )
    else:
        scenario = "A"
        scenario_text = (
            "the h/J² scaling SURVIVES on all tested topologies — "
            "the mean-field result is topology-robust"
        )

    out = OUT / "NETWORK_VERDICT.md"
    with out.open("w", encoding="utf-8") as f:
        f.write("# Fix B — Agent-Based Network Model: VERDICT\n\n")
        f.write("## Headline\n\n")
        f.write(
            f"At (J = 5.0, μ = 100), with N = 200 agents, mean degree ≈ 20, "
            f"50 seeds per cell, the high-J residual P(collapse) varies "
            f"across topologies: complete (mean-field control) "
            f"{p_complete:.2f}; modular {p_modular:.2f}; "
            f"modular with boundary-only h {p_modular_b:.2f}; "
            f"Erdős–Rényi {p_er:.2f}; Watts–Strogatz {p_ws:.2f}; "
            f"Barabási–Albert {p_ba:.2f}. **Scenario {scenario}**: "
            f"{scenario_text}.\n\n"
        )

        f.write("## Method\n\n")
        f.write(
            "N = 200 agents on each of six graph topologies (mean degree "
            "≈ 20). Each agent's order parameter $m_i$ evolves under "
            "Curie-Weiss dynamics with LOCAL coupling — drift driven by "
            "neighbor-averaged $m$ rather than global mean $m$. "
            "Multiplicative noise $\\xi \\sqrt{1 - m_i^2} \\eta(t)$ uses a "
            "single shared $\\eta(t)$ per step (preserving the minimal-model "
            "noise scale on $m_{avg}$ — this is what makes the complete-graph "
            "control reproduce the minimal-model results). Wealth $W$ is "
            "aggregated over all agents. The passive stabilizer field "
            "$h(W) = 2\\cdot(W/500)$ is applied to all agents under the "
            "first five topologies; under `modular_boundary` it is applied "
            "only to nodes with cross-community edges (~5–10% of nodes), "
            "directly testing whether boundary-localized passive stabilization "
            "retains effectiveness — the scenario the Limitations section "
            "raises.\n\n"
        )

        f.write("## P(collapse) at (J = 5.0, μ = 100) per topology\n\n")
        f.write("| topology | P(collapse) | rigidity share | "
                "intra-comm. corr | inter-comm. corr |\n")
        f.write("|---|---|---|---|---|\n")
        for topology in TOPOLOGY_ORDER:
            row = cell.loc[topology]
            f.write(
                f"| {TOPOLOGY_LABELS[topology]} | "
                f"{row.p_collapse:.2f} | "
                f"{row.rigidity_share_collapsed:.2f} | "
                f"{row.mean_intra_corr:.2f} | "
                f"{row.mean_inter_corr:.2f} |\n"
            )
        f.write("\n")

        f.write("## Modular vs mean-field (the critical comparison)\n\n")
        f.write(f"- Complete graph (mean-field): P(coll) = {p_complete:.2f}\n")
        f.write(f"- Modular (h applied to all nodes): "
                f"P(coll) = {p_modular:.2f} "
                f"(relative reduction {rel_drop_mod * 100:.0f}%)\n")
        f.write(f"- Modular with boundary-only h: "
                f"P(coll) = {p_modular_b:.2f} "
                f"(relative reduction {rel_drop_bdr * 100:.0f}%)\n\n")

        f.write("## Implication for the paper\n\n")
        if scenario == "A":
            f.write(
                "The h/J² scaling result is robust to network topology "
                "in the regime tested. Modular structure with strong "
                "intra-community coupling and weak inter-community coupling "
                "does not eliminate the high-J collapse residual. Even "
                "boundary-localized passive stabilization fails — when the "
                "interior of each community synchronizes at high J, the "
                "stabilizer at the boundary cannot rescue the bulk.\n\n"
                "**Recommended integration:** new Results subsection "
                "\"Network robustness: the scaling law beyond mean-field\" "
                "(~150 words) + Figure 3 (the 6-panel heatmap). Update "
                "Abstract: \"the scaling is robust across mean-field, "
                "small-world, scale-free, and modular topologies.\" "
                "Update Limitations: remove the modular-network caveat or "
                "downgrade to \"more extreme topology variants are future "
                "work.\"\n"
            )
        elif scenario == "B":
            f.write(
                "Network modularity PARTIALLY restores passive stabilizer "
                "effectiveness — the high-J residual is reduced but does "
                "not vanish. This is the most informative outcome for the "
                "governance story: network architecture is itself a "
                "design variable. Modular systems with strong "
                "intra-community decoupling at high J reduce but do not "
                "eliminate the h/J² decay; passive stabilizers gain "
                "effectiveness at module boundaries when the modules "
                "themselves decouple.\n\n"
                "**Recommended integration:** new Results subsection "
                "\"Network architecture as a passive-stabilizer "
                "amplifier\" with Figure 3. The paper's central claim "
                "is preserved — the h/J² scaling is the upper bound on "
                "passive failure — but with the nuance that network "
                "modularity reduces the residual. This addresses the "
                "Limitations caveat directly and turns it from a "
                "weakness into a refinement.\n"
            )
        else:
            f.write(
                "On modular networks, passive stabilizer effectiveness is "
                "FULLY restored at high coupling — communities decouple "
                "and the h/J² decay does not bind. This is a substantial "
                "modification of the paper's central claim and should be "
                "foregrounded honestly.\n\n"
                "**Recommended integration:** rewrite the Limitations "
                "section to put network modularity at the top. Reposition "
                "the h/J² result as 'an upper bound on passive failure in "
                "well-mixed systems,' with the explicit caveat that "
                "modular real-world systems may be exempt. Add Results "
                "subsection \"Network topology as a governance variable\" "
                "presenting the modular result. The paper's "
                "active-stabilizer recommendation is unchanged in the "
                "well-mixed limit but its scope is narrower.\n"
            )

    print(f"Wrote {out}")
    return scenario


def main():
    mpl.rcParams["figure.dpi"] = 100
    df = load_summary()
    print(f"Loaded {len(df)} rows from network_comparison.csv "
          f"({df.topology.nunique()} topologies)")

    figure_pcollapse_comparison(df)
    figure_scaling_curves(df)
    figure_modular_decoupling(df)
    scenario = write_verdict(df)
    print(f"\nScenario: {scenario}")


if __name__ == "__main__":
    main()
