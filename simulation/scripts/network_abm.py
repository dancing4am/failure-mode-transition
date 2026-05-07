"""
Fix B — Agent-Based Network Model.

N agents on a graph G = (V, E) evolve coupled Curie-Weiss-like dynamics
with LOCAL (neighbor-averaged) coupling rather than mean-field. The
shared-economy wealth W is aggregated over agents.

Tests whether Theorem 1's h/J^2 scaling survives across network
topologies — directly addressing the Limitations note that "structural
deviations from mean-field (strong heterogeneity, network modularity)
might in principle restore some protection."

Topologies (N = 200, mean degree k ~ 20):
  1. complete            — control, must reproduce minimal-model results
  2. erdos_renyi         — G(N, p) with p = k/(N-1)
  3. watts_strogatz      — small-world, k=20, beta=0.1
  4. barabasi_albert     — scale-free, m_attach=10
  5. modular             — 4 communities of 50, p_in=0.35, p_out=0.01
  6. modular_boundary    — modular, but h applied only to boundary nodes

Sweep: J in {0.5, 1.0, 2.0, 3.0, 4.0, 5.0} × mu in {20, 40, 60, 80, 100}
× 50 seeds × 20,000 steps at dt=0.05.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import time
import argparse

import numpy as np
import networkx as nx
import pandas as pd
import scipy.sparse as sp

T_TEMP = 2.0
XI = 0.5
DT = 0.05
N_STEPS = 20_000
INITIAL_M = 0.0
INITIAL_WEALTH = 100.0
COLLAPSE_WEALTH = 10.0
COLLAPSE_STREAK = 200
RIGIDITY_M = 0.9
FRAGMENTATION_M = 0.3
EPS = 1e-6
SEED_BASE = 0xC0DE_BEEF

N_AGENTS = 200
MEAN_DEGREE = 20

J_VALUES = (0.5, 1.0, 2.0, 3.0, 4.0, 5.0)
MULT_VALUES = (20, 40, 60, 80, 100)
N_SEEDS = 50

TOPOLOGIES = [
    "complete",
    "erdos_renyi",
    "watts_strogatz",
    "barabasi_albert",
    "modular",
    "modular_boundary",
]

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "network_abm"
OUT.mkdir(parents=True, exist_ok=True)


# -------------------------------------------------------------------------
# Graph construction
# -------------------------------------------------------------------------
def make_graph(topology, N=N_AGENTS, k=MEAN_DEGREE, seed=42):
    if topology == "complete":
        return nx.complete_graph(N)
    if topology == "erdos_renyi":
        p = k / (N - 1)
        return nx.erdos_renyi_graph(N, p, seed=seed)
    if topology == "watts_strogatz":
        return nx.watts_strogatz_graph(N, k, 0.1, seed=seed)
    if topology == "barabasi_albert":
        return nx.barabasi_albert_graph(N, k // 2, seed=seed)
    if topology in ("modular", "modular_boundary"):
        # 4 communities of 50; p_in = 0.35, p_out = 0.01
        sizes = [50, 50, 50, 50]
        p_in, p_out = 0.35, 0.01
        p_matrix = [[p_in if i == j else p_out for j in range(4)]
                    for i in range(4)]
        return nx.stochastic_block_model(sizes, p_matrix, seed=seed)
    raise ValueError(f"Unknown topology: {topology}")


def graph_communities(G):
    """Return a node->community-index dict from the SBM block attribute."""
    if "block" in next(iter(G.nodes(data=True)))[1]:
        return {n: G.nodes[n]["block"] for n in G.nodes}
    # fallback: each node its own community
    return {n: n for n in G.nodes}


def boundary_mask_for_modular(G, N=N_AGENTS):
    """Return a length-N boolean array marking nodes with cross-community
    neighbors (boundary nodes)."""
    comm = graph_communities(G)
    mask = np.zeros(N, dtype=bool)
    for i in G.nodes:
        c_i = comm[i]
        for j in G.neighbors(i):
            if comm[j] != c_i:
                mask[i] = True
                break
    return mask


def graph_summary(G, topology):
    degs = dict(G.degree())
    return {
        "topology": topology,
        "N": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
        "mean_degree": float(np.mean(list(degs.values()))),
        "min_degree": int(min(degs.values())),
        "max_degree": int(max(degs.values())),
    }


# -------------------------------------------------------------------------
# Simulation core
# -------------------------------------------------------------------------
@dataclass
class RunResult:
    collapsed: bool
    collapse_type: str  # 'none', 'rigidity', 'fragmentation', 'mixed'
    collapse_step: int
    final_m_avg: float
    final_W: float
    final_m_std: float  # cross-agent dispersion (coherence measure)
    intra_corr: float   # mean within-community pairwise corr of m trace
    inter_corr: float   # mean across-community pairwise corr of m trace


def run_abm(adj_norm_sparse, communities, J, mu, seed,
            n_steps=N_STEPS, boundary_only_h=False,
            boundary_mask=None,
            record_traces=False, n_traces=20):
    """Run one ABM trajectory. adj_norm_sparse is row-normalized adjacency
    in CSR form. communities is a length-N int array of community labels.
    Returns RunResult."""
    rng = np.random.default_rng(seed)
    N = adj_norm_sparse.shape[0]

    m = np.full(N, INITIAL_M)
    W = INITIAL_WEALTH
    streak = 0
    collapsed = False
    collapse_step = -1
    collapse_m_avg = 0.0
    sqrt_dt = np.sqrt(DT)

    if boundary_only_h:
        h_mask = boundary_mask.astype(np.float64)
    else:
        h_mask = np.ones(N)

    # For correlation measurement, record m at every step over a sub-sample
    if record_traces:
        sample_idx = rng.choice(N, size=min(n_traces, N), replace=False)
        sample_idx.sort()
        traces = np.empty((n_steps, len(sample_idx)))
    else:
        sample_idx = None
        traces = None

    for step in range(n_steps):
        m_local = adj_norm_sparse @ m  # neighbor-averaged
        h_field = (W / 500.0) * 2.0 * h_mask
        drift = -m + np.tanh((J * m_local + h_field) / T_TEMP)
        # Shared noise eta(t) across agents at each step. This preserves
        # the minimal-model noise scale on m_avg (independent per-agent
        # noise would average to xi/sqrt(N), which would suppress noise-
        # driven branch escape and break the complete-graph control).
        eta = rng.standard_normal()
        noise = XI * np.sqrt(np.maximum(1.0 - m * m, 0.0)) * eta * sqrt_dt
        m = np.clip(m + drift * DT + noise, -1 + EPS, 1 - EPS)

        m_avg = m.mean()
        employment = 0.5 + 0.5 * m_avg
        W = max(0.0, W + (mu * employment - 30.0 - 10.0 * (W / 500.0)) * DT)

        if W < COLLAPSE_WEALTH:
            streak += 1
            if streak >= COLLAPSE_STREAK and not collapsed:
                collapsed = True
                collapse_step = step
                collapse_m_avg = m_avg
        else:
            streak = 0

        if record_traces:
            traces[step] = m[sample_idx]

    final_m_avg = float(m.mean())
    final_W = float(W)
    final_m_std = float(m.std())

    if collapsed:
        if abs(collapse_m_avg) > RIGIDITY_M:
            ctype = "rigidity"
        elif abs(collapse_m_avg) < FRAGMENTATION_M:
            ctype = "fragmentation"
        else:
            ctype = "mixed"
    else:
        ctype = "none"

    intra_corr = np.nan
    inter_corr = np.nan
    if record_traces:
        # Pairwise correlation within / across communities on the sampled
        # nodes. We use the second half of the trace to avoid the burn-in.
        burn = n_steps // 2
        x = traces[burn:].T  # (n_traces, T/2)
        # Pearson correlation matrix
        x_centered = x - x.mean(axis=1, keepdims=True)
        x_std = x.std(axis=1, keepdims=True)
        x_std[x_std == 0] = 1.0
        x_norm = x_centered / x_std
        corr_mat = (x_norm @ x_norm.T) / x.shape[1]
        # Tag pairs as intra/inter
        s_comm = communities[sample_idx]
        n = len(sample_idx)
        intra_pairs, inter_pairs = [], []
        for i in range(n):
            for j in range(i + 1, n):
                if s_comm[i] == s_comm[j]:
                    intra_pairs.append(corr_mat[i, j])
                else:
                    inter_pairs.append(corr_mat[i, j])
        if intra_pairs:
            intra_corr = float(np.mean(intra_pairs))
        if inter_pairs:
            inter_corr = float(np.mean(inter_pairs))

    return RunResult(
        collapsed=collapsed,
        collapse_type=ctype,
        collapse_step=collapse_step,
        final_m_avg=final_m_avg,
        final_W=final_W,
        final_m_std=final_m_std,
        intra_corr=intra_corr,
        inter_corr=inter_corr,
    )


def normalized_adjacency(G):
    """Return row-normalized adjacency as a sparse CSR matrix."""
    A = nx.adjacency_matrix(G).astype(np.float64)
    deg = np.array(A.sum(axis=1)).ravel()
    deg[deg == 0] = 1.0
    inv_deg = sp.diags(1.0 / deg)
    return (inv_deg @ A).tocsr()


# -------------------------------------------------------------------------
# Sweep
# -------------------------------------------------------------------------
def sweep_topology(topology, j_values=J_VALUES, mu_values=MULT_VALUES,
                   n_seeds=N_SEEDS, record_traces_for_modular=True):
    G = make_graph(topology)
    summary = graph_summary(G, topology)
    print(f"\n[{topology}] graph: N={summary['N']}, "
          f"edges={summary['n_edges']}, mean_deg={summary['mean_degree']:.1f}, "
          f"min_deg={summary['min_degree']}, max_deg={summary['max_degree']}")

    adj_norm = normalized_adjacency(G)
    comm_dict = graph_communities(G)
    communities = np.array([comm_dict[i] for i in G.nodes])

    boundary_only_h = (topology == "modular_boundary")
    boundary_mask = boundary_mask_for_modular(G) if "modular" in topology else None

    record_traces = (topology in ("modular", "modular_boundary"))

    rows = []
    t0 = time.time()
    n_total = len(j_values) * len(mu_values) * n_seeds
    n_done = 0
    last_print = t0
    for j_idx, J in enumerate(j_values):
        for m_idx, mu in enumerate(mu_values):
            for s in range(n_seeds):
                seed = (SEED_BASE
                        + hash(topology) % 1_000_000
                        + j_idx * 1_000_003
                        + m_idx * 1009
                        + s * 13)
                # only record traces for a few seeds per cell to save mem
                rt = record_traces and s < 3
                r = run_abm(adj_norm, communities, J, mu, seed,
                            boundary_only_h=boundary_only_h,
                            boundary_mask=boundary_mask,
                            record_traces=rt, n_traces=20)
                rows.append({
                    "topology": topology,
                    "J": J, "mu": mu, "seed_idx": s,
                    "collapsed": int(r.collapsed),
                    "collapse_type": r.collapse_type,
                    "collapse_step": int(r.collapse_step),
                    "final_m_avg": r.final_m_avg,
                    "final_W": r.final_W,
                    "final_m_std": r.final_m_std,
                    "intra_corr": r.intra_corr,
                    "inter_corr": r.inter_corr,
                })
                n_done += 1
                now = time.time()
                if now - last_print > 30:
                    elapsed = now - t0
                    rate = n_done / elapsed
                    eta = (n_total - n_done) / max(rate, 0.01)
                    print(f"  [{topology}] {n_done}/{n_total} "
                          f"({elapsed:.0f}s elapsed, ETA {eta:.0f}s)")
                    last_print = now
    elapsed = time.time() - t0
    print(f"  [{topology}] done in {elapsed:.0f}s")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"{topology}_results.csv", index=False)

    # cell summary
    summary_rows = []
    for J in j_values:
        for mu in mu_values:
            sub = df[(df.J == J) & (df.mu == mu)]
            n_coll = int(sub.collapsed.sum())
            summary_rows.append({
                "topology": topology,
                "J": J, "mu": mu,
                "n_seeds": len(sub),
                "p_collapse": float(sub.collapsed.mean()),
                "p_rigidity": float((sub.collapse_type == "rigidity").mean()),
                "p_fragmentation": float((sub.collapse_type == "fragmentation").mean()),
                "rigidity_share_collapsed": (
                    float((sub.collapse_type == "rigidity").sum() / max(n_coll, 1))
                ),
                "mean_intra_corr": float(sub.intra_corr.mean()),
                "mean_inter_corr": float(sub.inter_corr.mean()),
            })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT / f"{topology}_summary.csv", index=False)
    return summary_df, summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--topologies", nargs="*", default=TOPOLOGIES)
    parser.add_argument("--seeds", type=int, default=N_SEEDS)
    args = parser.parse_args()

    all_summaries = []
    graph_summaries = []
    for topology in args.topologies:
        s_df, g_summary = sweep_topology(topology, n_seeds=args.seeds)
        all_summaries.append(s_df)
        graph_summaries.append(g_summary)

    pd.DataFrame(graph_summaries).to_csv(OUT / "graph_summaries.csv", index=False)
    if all_summaries:
        all_df = pd.concat(all_summaries, ignore_index=True)
        all_df.to_csv(OUT / "network_comparison.csv", index=False)
        print(f"\nWrote {OUT / 'network_comparison.csv'} ({len(all_df)} rows)")

        # Print headline numbers at (J=5, mu=100)
        print("\nP(collapse) at (J=5.0, mu=100):")
        for topology in args.topologies:
            sub = all_df[(all_df.topology == topology)
                         & (all_df.J == 5.0) & (all_df.mu == 100)]
            if len(sub):
                row = sub.iloc[0]
                print(f"  {topology:25s}  P_coll = {row.p_collapse:.2f}  "
                      f"rigidity_share = {row.rigidity_share_collapsed:.2f}")


if __name__ == "__main__":
    main()
