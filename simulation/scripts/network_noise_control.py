# -*- coding: utf-8 -*-
"""
Common-mode vs idiosyncratic noise control for the network ABM (SI §S9).

Re-runs the headline cell (J=5, mu=100) across all six topologies under (a) the
committed SHARED common-mode shock eta(t) and (b) INDEPENDENT per-agent noise
eta_i(t), plus the modular intra/inter-community correlation under each mode and
an N-scan of the passive collapse rate under independent noise. Faithful to the
committed dynamics (imports network_abm); deterministic seeds. Writes
results/network_abm/noise_control.csv.

Run: py -3.13 network_noise_control.py   (requires networkx, as network_abm.py does)
"""
import sys
import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import network_abm as nab  # noqa: E402

J, MU, NSEEDS = 5.0, 100, 50
OUT = HERE.parents[0] / "results" / "network_abm" / "noise_control.csv"


def run_traj(adj, J, mu, seed, mode, bmask=None, rec_nodes=None):
    rng = np.random.default_rng(seed)
    N = adj.shape[0]
    m = np.full(N, nab.INITIAL_M); W = nab.INITIAL_WEALTH; streak = 0; collapsed = False
    sdt = np.sqrt(nab.DT)
    hmask = np.ones(N) if bmask is None else bmask.astype(float)
    rec = np.empty((nab.N_STEPS, len(rec_nodes))) if rec_nodes is not None else None
    for step in range(nab.N_STEPS):
        m_local = adj @ m
        drift = -m + np.tanh((J * m_local + (W / 500.0) * 2.0 * hmask) / nab.T_TEMP)
        eta = rng.standard_normal() if mode == "shared" else rng.standard_normal(N)
        m = np.clip(m + drift * nab.DT + nab.XI * np.sqrt(np.maximum(1 - m * m, 0.0)) * eta * sdt,
                    -1 + nab.EPS, 1 - nab.EPS)
        mavg = m.mean()
        W = max(0.0, W + (mu * (0.5 + 0.5 * mavg) - 30.0 - 10.0 * (W / 500.0)) * nab.DT)
        if W < nab.COLLAPSE_WEALTH:
            streak += 1
            collapsed = collapsed or (streak >= nab.COLLAPSE_STREAK)
        else:
            streak = 0
        if rec is not None:
            rec[step] = m[rec_nodes]
    return collapsed, rec


def intra_inter(rec, comm_s):
    x = rec[nab.N_STEPS // 2:].T
    xc = x - x.mean(1, keepdims=True); xs = x.std(1, keepdims=True); xs[xs == 0] = 1.0
    C = (xc / xs) @ (xc / xs).T / x.shape[1]
    intra, inter = [], []
    for i in range(len(comm_s)):
        for j in range(i + 1, len(comm_s)):
            (intra if comm_s[i] == comm_s[j] else inter).append(C[i, j])
    return (float(np.mean(intra)) if intra else np.nan,
            float(np.mean(inter)) if inter else np.nan)


def main():
    rows = []
    print(f"=== Network noise control: J={J}, mu={MU}, {NSEEDS} seeds ===")
    print(f"{'topology':18s} {'P_shared':>10s} {'P_indep':>10s}")
    for ti, topo in enumerate(nab.TOPOLOGIES):
        G = nab.make_graph(topo); adj = nab.normalized_adjacency(G)
        bmask = nab.boundary_mask_for_modular(G) if topo == "modular_boundary" else None
        out = {}
        for mode in ("shared", "independent"):
            nc = sum(int(run_traj(adj, J, MU,
                                  700_000 + ti * 101 + s * 13 + (0 if mode == "shared" else 5_000_000),
                                  mode, bmask=bmask)[0]) for s in range(NSEEDS))
            out[mode] = nc / NSEEDS
        rows.append({"kind": "p_collapse", "topology": topo,
                     "shared": out["shared"], "independent": out["independent"]})
        print(f"{topo:18s} {out['shared']:>10.3f} {out['independent']:>10.3f}")

    # modular intra/inter correlation under each mode
    G = nab.make_graph("modular"); adj = nab.normalized_adjacency(G)
    cdict = nab.graph_communities(G); comm = np.array([cdict[i] for i in G.nodes])
    samp = np.sort(np.random.default_rng(12345).choice(adj.shape[0], 20, replace=False))
    for mode in ("shared", "independent"):
        ia, it = [], []
        for s in range(5):
            _, rec = run_traj(adj, J, MU, 800_000 + s * 13 + (0 if mode == "shared" else 5_000_000),
                              mode, rec_nodes=samp)
            a, b = intra_inter(rec, comm[samp]); ia.append(a); it.append(b)
        rows.append({"kind": f"modular_corr_{mode}", "topology": "modular",
                     "shared": float(np.nanmean(ia)), "independent": float(np.nanmean(it))})
        print(f"modular {mode:11s}: intra={np.nanmean(ia):+.3f} inter={np.nanmean(it):+.3f}")

    # N-scan: passive, independent per-agent noise, mean-field coupling (J*m_avg).
    # N=1 recovers the scalar model; larger N shows the xi/sqrt(N) suppression.
    def run_mf(N, seed):
        rng = np.random.default_rng(seed)
        m = np.full(N, nab.INITIAL_M); W = nab.INITIAL_WEALTH; streak = 0; collapsed = False
        sdt = np.sqrt(nab.DT)
        for _ in range(nab.N_STEPS):
            mavg = m.mean()
            drift = -m + np.tanh((J * mavg + (W / 500.0) * 2.0) / nab.T_TEMP)
            eta = rng.standard_normal(N)
            m = np.clip(m + drift * nab.DT + nab.XI * np.sqrt(np.maximum(1 - m * m, 0.0)) * eta * sdt,
                        -1 + nab.EPS, 1 - nab.EPS)
            mavg = m.mean()
            W = max(0.0, W + (MU * (0.5 + 0.5 * mavg) - 30.0 - 10.0 * (W / 500.0)) * nab.DT)
            if W < nab.COLLAPSE_WEALTH:
                streak += 1; collapsed = collapsed or (streak >= nab.COLLAPSE_STREAK)
            else:
                streak = 0
        return collapsed
    print("N-scan (passive, independent, mean-field coupling):")
    for N in (1, 5, 25, 200):
        pN = sum(int(run_mf(N, 6_000_000 + N + s * 13)) for s in range(NSEEDS)) / NSEEDS
        rows.append({"kind": "Nscan_indep", "topology": f"N={N}", "shared": np.nan, "independent": pN})
        print(f"  N={N:>3d}: {pN:.3f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["kind", "topology", "shared", "independent"])
        w.writeheader(); w.writerows(rows)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
