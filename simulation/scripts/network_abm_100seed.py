"""
Targeted ABM 100-seed re-run for the headline cell (J=5, mu=100) on
all six topologies. Tightens the SE on the most-quoted cell so the
ABM SE matches the scalar SDE baseline (which used 100 seeds).

Reuses graph construction and run_abm from network_abm.py. Only the
single (J, mu) cell is rerun, but at 100 seeds.

Output: results/network_abm/headline_100seed.csv
"""

from __future__ import annotations

from pathlib import Path
import time

import numpy as np
import pandas as pd

from network_abm import (
    make_graph, graph_communities, boundary_mask_for_modular,
    normalized_adjacency, run_abm, TOPOLOGIES, N_AGENTS,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "network_abm" / "headline_100seed.csv"

J = 5.0
MU = 100
N_SEEDS = 100
SEED_BASE = 0xC0DE_BEEF


def main():
    rows = []
    for topology in TOPOLOGIES:
        t0 = time.time()
        G = make_graph(topology, seed=42)
        comm = graph_communities(G)
        comm_array = np.array([comm[i] for i in range(N_AGENTS)])
        boundary_only = (topology == "modular_boundary")
        boundary_mask = (boundary_mask_for_modular(G)
                         if boundary_only else None)
        adj_norm = normalized_adjacency(G)

        n_collapsed = 0
        rigidity = 0
        fragmentation = 0
        mixed = 0

        for s in range(N_SEEDS):
            seed = SEED_BASE + s * 7919
            r = run_abm(adj_norm, comm_array, J, MU, seed,
                        boundary_only_h=boundary_only,
                        boundary_mask=boundary_mask)
            if r.collapsed:
                n_collapsed += 1
                if r.collapse_type == "rigidity":
                    rigidity += 1
                elif r.collapse_type == "fragmentation":
                    fragmentation += 1
                elif r.collapse_type == "mixed":
                    mixed += 1

        p = n_collapsed / N_SEEDS
        se = float(np.sqrt(p * (1 - p) / N_SEEDS))
        rigidity_share = rigidity / max(n_collapsed, 1)
        elapsed = time.time() - t0
        print(f"{topology:20s} P(coll) = {p:.2f} (SE = {se:.3f})  "
              f"rigidity = {rigidity_share:.2f}  [{elapsed:.0f}s]")

        rows.append({
            "topology": topology,
            "J": J,
            "mu": MU,
            "n_seeds": N_SEEDS,
            "p_collapse": p,
            "SE": se,
            "rigidity_share": rigidity_share,
            "n_collapsed": n_collapsed,
            "n_rigidity": rigidity,
            "n_fragmentation": fragmentation,
            "n_mixed": mixed,
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
