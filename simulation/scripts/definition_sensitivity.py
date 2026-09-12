"""Collapse-definition sensitivity at the headline cell.

One set of n=2000 matched-seed trajectories at the headline cell
(J=5, mu=100, passive stabilizer), base integrator exactly as in
minimal_model.py (EM, dt=0.05, 20,000 steps, T=2, xi=0.5, clip to
[-1+1e-6, 1-1e-6], h = 2*(W/500), wealth floored at 0).

All collapse definitions are evaluated ONLINE on the same paths in a single
pass:

  * 3x3 definition grid: W-threshold {5, 10, 20} x duration {100, 200, 400}
    consecutive steps  ->  P(collapse) for each of the 9 definitions.
  * Type gates for the baseline definition (threshold 10, duration 200):
    rigidity/fragmentation/mixed shares under |m| gate pairs
    (rigidity-threshold, fragmentation-threshold):
    (0.8,0.3), (0.9,0.3), (0.95,0.3), (0.9,0.2), (0.9,0.4).

Outputs -> simulation/results/definition_sensitivity/
  definition_grid.csv, type_gates.csv, SUMMARY.md
"""

from __future__ import annotations

import os
import time

import numpy as np

# Base-integrator constants (must match minimal_model.py)
T_TEMP = 2.0
XI = 0.5
DT = 0.05
N_STEPS = 20_000
EPS = 1e-6
INITIAL_M = 0.0
INITIAL_WEALTH = 100.0

J = 5.0
MULT = 100.0
N_TRAJ = 2000
SEED = 0xC0DE + 9 * 1_000_003 + 8 * 1009   # headline cell (j_idx=9, m_idx=8)

W_THRESHOLDS = (5.0, 10.0, 20.0)
DURATIONS = (100, 200, 400)
BASELINE = (10.0, 200)
GATE_PAIRS = ((0.8, 0.3), (0.9, 0.3), (0.95, 0.3), (0.9, 0.2), (0.9, 0.4))

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.normpath(os.path.join(BASE, "..", "results", "definition_sensitivity"))


def run():
    rng = np.random.default_rng(SEED)
    n = N_TRAJ
    m = np.full(n, INITIAL_M, dtype=np.float64)
    wealth = np.full(n, INITIAL_WEALTH, dtype=np.float64)
    sqrt_dt = np.sqrt(DT)

    nthr = len(W_THRESHOLDS)
    streaks = np.zeros((nthr, n), dtype=np.int64)

    combos = [(wt, du) for wt in W_THRESHOLDS for du in DURATIONS]
    collapsed = {cb: np.zeros(n, dtype=bool) for cb in combos}
    collapse_step = {cb: np.full(n, -1, dtype=np.int64) for cb in combos}
    collapse_absm = {cb: np.full(n, np.nan) for cb in combos}

    for t in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = (wealth / 500.0) * 2.0
        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)

        employment = 0.5 + 0.5 * m
        income = MULT * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        abs_m = None
        for k, wt in enumerate(W_THRESHOLDS):
            is_low = wealth < wt
            streaks[k] = np.where(is_low, streaks[k] + 1, 0)
            for du in DURATIONS:
                cb = (wt, du)
                new = (streaks[k] >= du) & (~collapsed[cb])
                if new.any():
                    if abs_m is None:
                        abs_m = np.abs(m)
                    collapsed[cb] |= new
                    collapse_step[cb][new] = t
                    collapse_absm[cb][new] = abs_m[new]

    return collapsed, collapse_step, collapse_absm


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    t0 = time.time()
    print(f"E7: headline cell J={J:g}, mu={MULT:g}, n={N_TRAJ}, seed={SEED} ...")
    collapsed, collapse_step, collapse_absm = run()
    print(f"  run done ({time.time()-t0:.0f}s)")

    # ---- 9-cell definition grid -------------------------------------------
    grid_rows = []
    for wt in W_THRESHOLDS:
        for du in DURATIONS:
            cb = (wt, du)
            col = collapsed[cb]
            p = col.mean()
            se = np.sqrt(p * (1 - p) / N_TRAJ)
            med = np.median(collapse_step[cb][col]) if col.any() else np.nan
            grid_rows.append((wt, du, int(col.sum()), p, 1.96 * se, med))

    with open(os.path.join(OUTDIR, "definition_grid.csv"), "w", encoding="utf-8") as fh:
        fh.write("w_threshold,duration_steps,n_collapsed,p_collapse,"
                 "ci95_halfwidth,collapse_step_median\n")
        for wt, du, nc, p, ci, med in grid_rows:
            fh.write(f"{wt:g},{du},{nc},{p:.5f},{ci:.5f},{med:.1f}\n")

    # ---- type gates on the baseline definition ----------------------------
    cb0 = BASELINE
    col0 = collapsed[cb0]
    absm0 = collapse_absm[cb0][col0]
    nc0 = int(col0.sum())

    gate_rows = []
    for rg, fg in GATE_PAIRS:
        rig = (absm0 > rg).mean()
        frag = (absm0 < fg).mean()
        mixed = 1.0 - rig - frag
        gate_rows.append((rg, fg, rig, frag, mixed))

    with open(os.path.join(OUTDIR, "type_gates.csv"), "w", encoding="utf-8") as fh:
        fh.write("rigidity_gate,fragmentation_gate,rigidity_share,"
                 "fragmentation_share,mixed_share,n_collapsed\n")
        for rg, fg, r_, f_, mx in gate_rows:
            fh.write(f"{rg},{fg},{r_:.5f},{f_:.5f},{mx:.5f},{nc0}\n")

    # ---- SUMMARY ----------------------------------------------------------
    L = []
    A = L.append
    A("# Collapse-Definition Sensitivity")
    A("")
    A(f"Headline cell J={J:g}, mu={MULT:g}, passive stabilizer; base EM "
      f"integrator (dt={DT}, {N_STEPS} steps, T={T_TEMP:g}, xi={XI:g}); "
      f"n={N_TRAJ} trajectories, matched seed {SEED} "
      "(0xC0DE + 9*1000003 + 8*1009). All definitions evaluated online on the "
      "SAME paths in one pass (independent streak counters per W-threshold).")
    A("")
    A("## P(collapse) under 9 definitions (W-threshold x duration)")
    A("")
    A("Collapse = W < threshold for `duration` consecutive steps "
      "(duration 200 steps = 10 physical time units at dt=0.05).")
    A("")
    hdr = " | ".join(f"dur {du}" for du in DURATIONS)
    A(f"| W-thr \\\\ duration | {hdr} |")
    A("|---|" + "---|" * len(DURATIONS))
    for wt in W_THRESHOLDS:
        cells = []
        for du in DURATIONS:
            p = collapsed[(wt, du)].mean()
            se = 1.96 * np.sqrt(p * (1 - p) / N_TRAJ)
            cells.append(f"{p:.4f} +/- {se:.4f}")
        A(f"| W<{wt:g} | " + " | ".join(cells) + " |")
    A("")
    base_p = collapsed[cb0].mean()
    ps = [collapsed[cb].mean() for cb in collapsed]
    A(f"Baseline (W<10, 200 steps): P = {base_p:.4f} "
      f"(n_collapsed = {nc0}). Range across all 9 definitions: "
      f"{min(ps):.4f} - {max(ps):.4f}.")
    A("")
    A("## Type shares under |m| gate pairs (baseline definition, W<10 x 200)")
    A("")
    A("Gate pair = (rigidity threshold r, fragmentation threshold f): "
      "rigidity if |m| > r at collapse, fragmentation if |m| < f, else mixed. "
      f"Shares among the {nc0} collapsed trajectories.")
    A("")
    A("| gates (r, f) | rigidity | fragmentation | mixed |")
    A("|---|---|---|---|")
    for rg, fg, r_, f_, mx in gate_rows:
        tag = " (baseline)" if (rg, fg) == (0.9, 0.3) else ""
        A(f"| ({rg:g}, {fg:g}){tag} | {r_:.4f} | {f_:.4f} | {mx:.4f} |")
    A("")
    A(f"Runtime: {time.time()-t0:.0f}s. Files: definition_grid.csv, "
      "type_gates.csv.")

    with open(os.path.join(OUTDIR, "SUMMARY.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")

    print(f"Done in {time.time()-t0:.0f}s -> {OUTDIR}")


if __name__ == "__main__":
    main()
