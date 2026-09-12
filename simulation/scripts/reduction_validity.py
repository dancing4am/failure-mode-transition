# -*- coding: utf-8 -*-
"""
Validity range of the one-dimensional reduction (D0 supporting computation).

The manuscript's scalar model is the N-agent ABM reduced to one dimension.
The reduction's noise identity is

    sigma_eff = sqrt(sigma_ind^2 / N + sigma_com^2),

which is nothing more than the statement that independent per-agent noise
averages out in the population mean while common-mode noise does not.
This script documents the RANGE OF NOISE COMPOSITIONS over which the 1-D
reduction with sigma_eff reproduces the ABM collapse statistics.

Endpoints:
  * Common-mode endpoint (r = 1, all noise shared): EXACT. With uniform
    initialization, uniform h, row-stochastic adjacency, and a single shared
    per-step shock, every agent receives identical increments, so the N-agent
    system IS the scalar SDE (bit-identical trajectories; see
    results/network_abm/headline_100seed.csv where 5 of 6 topologies give
    literally identical counts 33/30/0/3). No rerun needed; cited only.
  * Idiosyncratic endpoint (r = 0, all noise independent): the ABM N-scan of
    network_noise_control.py (J=5, mu=100, xi=0.5, N in {1,5,25,200}) gave
    P(collapse) = 0.24 / 0.10 / 0.02 / 0.00. Part 1 below runs the MATCHED
    scalar model at xi_eff = 0.5/sqrt(N) that was never run.
  * Mixed noise (Part 2): ABM at N=200 (complete graph, committed dynamics)
    with noise split sigma_com^2 + sigma_ind^2 = 0.25 at mixing ratio
    r = sigma_com^2 / 0.25 in {0, 0.25, 0.5, 0.75, 1}, vs the scalar model
    at matched xi_eff.

Everything at the headline cell (J = 5.0, mu = 100), dt = 0.05, 20,000 steps,
passive stabilizer h = 2*(W/500), collapse = W < 10 for 200 consecutive steps.
Deterministic seeds throughout. Results are written incrementally to
results/reduction_validity/reduction_validity.csv.

Run: py -3.13 reduction_validity.py
"""
from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import minimal_model as mm  # noqa: E402
import network_abm as nab   # noqa: E402

# ---------------------------------------------------------------------------
# Configuration (mirrors network_noise_control.py's N-scan cell)
# ---------------------------------------------------------------------------
J, MU = 5.0, 100
XI_TOTAL = 0.5                      # total noise scale: sigma_ind^2+sigma_com^2 = 0.25
N_SCAN = (1, 5, 25, 200)            # N values of the committed ABM N-scan
NSEEDS_ABM_NSCAN = 50               # committed N-scan seed count
NSEEDS_ABM_MIXED = 100              # Part 2 ABM seeds per mixing ratio
NSEEDS_SCALAR = 1000                # matched scalar seeds
RATIOS = (0.0, 0.25, 0.5, 0.75, 1.0)
N_MIXED = 200
BOOT_B = 10_000
BOOT_SEED = 20260819

OUT_DIR = HERE.parents[0] / "results" / "reduction_validity"
OUT_CSV = OUT_DIR / "reduction_validity.csv"
FIELDS = ["source", "N", "r", "xi_eff", "P_collapse", "ci_lo", "ci_hi", "n", "part"]

MAX_PART2_SECONDS = 3 * 3600        # run reduced scope if estimate exceeds this


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def bootstrap_ci(binary: np.ndarray, B: int = BOOT_B, seed: int = BOOT_SEED):
    """Percentile bootstrap 95% CI for a proportion from a 0/1 array."""
    rng = np.random.default_rng(seed)
    n = len(binary)
    idx = rng.integers(0, n, size=(B, n))
    means = binary[idx].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def append_rows(rows):
    new_file = not OUT_CSV.exists()
    with open(OUT_CSV, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new_file:
            w.writeheader()
        w.writerows(rows)


def make_row(source, N, r, xi_eff, binary, part):
    lo, hi = bootstrap_ci(binary)
    return {
        "source": source, "N": N, "r": r, "xi_eff": round(xi_eff, 6),
        "P_collapse": float(binary.mean()), "ci_lo": lo, "ci_hi": hi,
        "n": len(binary), "part": part,
    }


# ---------------------------------------------------------------------------
# Part 1a — replicate the committed ABM N-scan (independent per-agent noise)
# to recover the per-seed binary outcomes for bootstrap CIs.
# EXACT replication of network_noise_control.run_mf (same dynamics, same
# deterministic seeds 6_000_000 + N + s*13).
# ---------------------------------------------------------------------------
def run_mf_indep(N, seed):
    """One trajectory: mean-field coupling J*m_avg, passive uniform h,
    INDEPENDENT per-agent multiplicative noise at xi = XI_TOTAL."""
    rng = np.random.default_rng(seed)
    m = np.full(N, nab.INITIAL_M)
    W = nab.INITIAL_WEALTH
    streak = 0
    collapsed = False
    sdt = np.sqrt(nab.DT)
    for _ in range(nab.N_STEPS):
        mavg = m.mean()
        drift = -m + np.tanh((J * mavg + (W / 500.0) * 2.0) / nab.T_TEMP)
        eta = rng.standard_normal(N)
        m = np.clip(m + drift * nab.DT
                    + nab.XI * np.sqrt(np.maximum(1 - m * m, 0.0)) * eta * sdt,
                    -1 + nab.EPS, 1 - nab.EPS)
        mavg = m.mean()
        W = max(0.0, W + (MU * (0.5 + 0.5 * mavg) - 30.0 - 10.0 * (W / 500.0)) * nab.DT)
        if W < nab.COLLAPSE_WEALTH:
            streak += 1
            collapsed = collapsed or (streak >= nab.COLLAPSE_STREAK)
        else:
            streak = 0
    return collapsed


def part1():
    print("=== Part 1: idiosyncratic endpoint (r = 0), N-scan ===")
    ref = {1: 0.24, 5: 0.10, 25: 0.02, 200: 0.00}  # committed noise_control.csv values
    for N in N_SCAN:
        t0 = time.time()
        binary = np.array([run_mf_indep(N, 6_000_000 + N + s * 13)
                           for s in range(NSEEDS_ABM_NSCAN)], dtype=float)
        p = binary.mean()
        flag = "OK" if abs(p - ref[N]) < 1e-12 else f"MISMATCH vs committed {ref[N]}"
        print(f"  ABM  N={N:>3d} indep-noise: P={p:.3f} "
              f"({time.time()-t0:.0f}s) [{flag}]")
        xi_eff = XI_TOTAL / np.sqrt(N)
        append_rows([make_row("abm", N, 0.0, xi_eff, binary, "part1_nscan")])

        # Matched scalar run at xi_eff = 0.5/sqrt(N), same (J, mult) cell.
        t0 = time.time()
        res = mm.run_cell(J, MU, n_seeds=NSEEDS_SCALAR, xi=xi_eff,
                          seed=0xA11CE + N)
        sb = res.collapsed.astype(float)
        print(f"  SCALAR xi_eff={xi_eff:.4f} (=0.5/sqrt({N})): "
              f"P={sb.mean():.3f} ({time.time()-t0:.0f}s, n={NSEEDS_SCALAR})")
        append_rows([make_row("scalar", N, 0.0, xi_eff, sb, "part1_nscan")])


# ---------------------------------------------------------------------------
# Part 2 — mixed-noise sweep at N = 200 on the committed complete graph.
# Noise per agent: sqrt(1-m_i^2) * (sigma_ind*eta_i + sigma_com*eta_shared),
# with sigma_com^2 + sigma_ind^2 = XI_TOTAL^2 = 0.25 and
# r = sigma_com^2 / (sigma_com^2 + sigma_ind^2).
# Coupling is the committed local (neighbor-averaged) coupling of
# network_abm.py: J * (A_rownorm @ m). h is uniform (passive).
# Vectorized across seeds.
# ---------------------------------------------------------------------------
def run_abm_mixed_batch(adj_dense, r, n_seeds, seed, n_steps=nab.N_STEPS):
    sigma_com = XI_TOTAL * np.sqrt(r)
    sigma_ind = XI_TOTAL * np.sqrt(1.0 - r)
    rng = np.random.default_rng(seed)
    S, N = n_seeds, adj_dense.shape[0]
    m = np.full((S, N), nab.INITIAL_M)
    W = np.full(S, nab.INITIAL_WEALTH)
    streak = np.zeros(S, dtype=np.int32)
    collapsed = np.zeros(S, dtype=bool)
    sdt = np.sqrt(nab.DT)
    adjT = adj_dense.T.copy()
    for _ in range(n_steps):
        m_local = m @ adjT                      # == (adj @ m[s]) per seed
        h = (W / 500.0) * 2.0                   # uniform passive field, per seed
        drift = -m + np.tanh((J * m_local + h[:, None]) / nab.T_TEMP)
        eta = sigma_ind * rng.standard_normal((S, N)) \
            + sigma_com * rng.standard_normal((S, 1))
        m = np.clip(m + drift * nab.DT
                    + np.sqrt(np.maximum(1 - m * m, 0.0)) * eta * sdt,
                    -1 + nab.EPS, 1 - nab.EPS)
        mavg = m.mean(axis=1)
        W = np.maximum(0.0, W + (MU * (0.5 + 0.5 * mavg) - 30.0
                                 - 10.0 * (W / 500.0)) * nab.DT)
        low = W < nab.COLLAPSE_WEALTH
        streak = np.where(low, streak + 1, 0)
        collapsed |= streak >= nab.COLLAPSE_STREAK
    return collapsed.astype(float)


def part2():
    print("=== Part 2: mixed-noise sweep at N = 200 (complete graph) ===")
    G = nab.make_graph("complete", N=N_MIXED)
    adj = np.asarray(nab.normalized_adjacency(G).todense())

    # --- timing estimate: one small pilot cell, extrapolated ---
    pilot_seeds, pilot_steps = 10, 1000
    t0 = time.time()
    run_abm_mixed_batch(adj, 0.5, pilot_seeds, seed=1, n_steps=pilot_steps)
    pilot = time.time() - t0
    per_full_cell = pilot * (NSEEDS_ABM_MIXED / pilot_seeds) * (nab.N_STEPS / pilot_steps)
    est_total = per_full_cell * len(RATIOS)
    print(f"  pilot: {pilot:.1f}s for {pilot_seeds} seeds x {pilot_steps} steps "
          f"-> est {per_full_cell/60:.1f} min/cell, {est_total/60:.1f} min total")
    if est_total > MAX_PART2_SECONDS:
        print("  ESTIMATE EXCEEDS BUDGET — this is already the reduced scope "
              "(N=200 only, 5 ratios, 100 seeds); proceeding anyway is not "
              "possible within budget. Aborting Part 2.")
        return False

    for r in RATIOS:
        sigma_com2 = XI_TOTAL**2 * r
        sigma_ind2 = XI_TOTAL**2 * (1.0 - r)
        xi_eff = float(np.sqrt(sigma_ind2 / N_MIXED + sigma_com2))

        t0 = time.time()
        ab = run_abm_mixed_batch(adj, r, NSEEDS_ABM_MIXED,
                                 seed=9_000_000 + int(round(r * 100)))
        print(f"  ABM    r={r:.2f}  xi_eff={xi_eff:.4f}  P={ab.mean():.3f} "
              f"({time.time()-t0:.0f}s, n={NSEEDS_ABM_MIXED})")
        append_rows([make_row("abm", N_MIXED, r, xi_eff, ab, "part2_mixed")])

        t0 = time.time()
        res = mm.run_cell(J, MU, n_seeds=NSEEDS_SCALAR, xi=xi_eff,
                          seed=0xB0B0 + int(round(r * 100)))
        sb = res.collapsed.astype(float)
        print(f"  SCALAR r={r:.2f}  xi_eff={xi_eff:.4f}  P={sb.mean():.3f} "
              f"({time.time()-t0:.0f}s, n={NSEEDS_SCALAR})")
        append_rows([make_row("scalar", N_MIXED, r, xi_eff, sb, "part2_mixed")])
    return True


# ---------------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if OUT_CSV.exists():
        OUT_CSV.unlink()  # fresh deterministic rebuild
    t_start = time.time()
    print(f"Reduction-validity runs: J={J}, mu={MU}, dt={nab.DT}, "
          f"steps={nab.N_STEPS}, total noise xi={XI_TOTAL}")
    print("Common-mode endpoint (r=1): EXACT by construction — cited, not rerun "
          "(headline_100seed.csv: 5/6 topologies bit-identical to scalar SDE).")
    part1()
    ran_part2 = part2()
    print(f"\nDone in {(time.time()-t_start)/60:.1f} min. "
          f"Part 2 {'completed' if ran_part2 else 'SKIPPED'}. Wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
