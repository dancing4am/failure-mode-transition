"""Delivered mean of the realized field over selection-window sizes, per arm.

Contribution 2 needs the DELIVERED MEAN exactly as Methods defines it — the
realized field h(W(t), m(t)) time-averaged over an interval — not a
convention-pinned proxy. This script measures it for all four static-sweep
arms (passive; stress at alpha=2.0 and at the informative alpha=0.05; coupling;
quadratic) at J in {4.0, 4.5, 5.0}, mu=100, n=2000 matched seeds, over windows
t <= 1, 5, 14, 50 and the full trajectory, and reports each arm's P(collapse)
on the same seeds.

Field forms (Methods "Simulation protocols"):
  passive   : h = 2*(W/500)
  stress    : h = 2*(W/500) + alpha*max(0,-m)*J          (alpha = 2.0 or 0.05)
  coupling  : h = 2*(W/500)*(1 + alpha*J)                (alpha = 1.0)
  quadratic : h = 2*(W/500)*(1 + alpha*(J/5)^2)          (alpha = 1.0)

Constants identical to minimal_model.py (dt=0.05, 20000 steps, T=2, xi=0.5,
m0=0, W0=100, collapse = W<10 for 200 consecutive steps, mu=100).

Matched seeds: one seed per J, shared across all arms at that J, so the
per-step noise draws are identical (paired). Output:
  simulation/results/magnitude_matched/delivered_mean_windows.csv
"""
from __future__ import annotations

from pathlib import Path
import time

import numpy as np
import pandas as pd

T_TEMP, XI, DT = 2.0, 0.5, 0.05
N_STEPS = 20_000
INITIAL_M, INITIAL_WEALTH = 0.0, 100.0
COLLAPSE_WEALTH, COLLAPSE_STREAK = 10.0, 200
EPS = 1e-6
MU = 100.0
N_SEEDS = 2000
SEED_BASE = 0xDE_11A7           # distinct base for this measurement
CELLS_J = (4.0, 4.5, 5.0)
WINDOWS = (1.0, 5.0, 14.0, 50.0)     # physical-time window edges
ARMS = (("passive", 0.0), ("stress_a2.0", 2.0), ("stress_a0.05", 0.05),
        ("coupling", 1.0), ("quadratic", 1.0))
OUT = Path(__file__).resolve().parents[1] / "results" / "magnitude_matched"


def field(arm: str, alpha: float, wealth: np.ndarray, m_c: np.ndarray, J: float):
    base = 2.0 * (wealth / 500.0)
    if arm == "passive":
        return base
    if arm.startswith("stress"):
        return base + alpha * np.maximum(0.0, -m_c) * J
    if arm == "coupling":
        return base * (1.0 + alpha * J)
    if arm == "quadratic":
        return base * (1.0 + alpha * (J / 5.0) ** 2)
    raise ValueError(arm)


def run_arm(arm: str, alpha: float, J: float, seed: int, n_seeds: int = N_SEEDS):
    """Returns (P(collapse), {window_edge: delivered_mean}, full_mean)."""
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M)
    wealth = np.full(n_seeds, INITIAL_WEALTH)
    collapsed = np.zeros(n_seeds, dtype=bool)
    streak = np.zeros(n_seeds, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)
    win_steps = {w: int(round(w / DT)) for w in WINDOWS}
    win_sum = {w: 0.0 for w in WINDOWS}
    win_n = {w: 0 for w in WINDOWS}
    full_sum = 0.0
    for t in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = field(arm, alpha, wealth, m_c, J)
        hsum = float(h.sum())
        full_sum += hsum
        for w, ns in win_steps.items():
            if t < ns:
                win_sum[w] += hsum
                win_n[w] += n_seeds
        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)
        income = MU * (0.5 + 0.5 * m)
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)
        streak = np.where(wealth < COLLAPSE_WEALTH, streak + 1, 0)
        collapsed |= streak >= COLLAPSE_STREAK
    means = {w: (win_sum[w] / win_n[w] if win_n[w] else float("nan")) for w in WINDOWS}
    full_mean = full_sum / (N_STEPS * n_seeds)
    return float(collapsed.mean()), means, full_mean


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for J in CELLS_J:
        seed = SEED_BASE + int(round(J * 10)) * 1_000_003
        for arm, alpha in ARMS:
            t0 = time.time()
            p, means, full_mean = run_arm(arm, alpha, J, seed)
            row = dict(J=J, mu=MU, arm=arm, alpha=alpha, p_collapse=p)
            for w in WINDOWS:
                row[f"h_mean_t<={w:g}"] = means[w]
            row["h_mean_full"] = full_mean
            rows.append(row)
            print(f"J={J} {arm:13s} a={alpha:<4g} P={p:.4f}  "
                  f"h[t<=1]={means[1.0]:.3f} h[t<=5]={means[5.0]:.3f} "
                  f"h[t<=14]={means[14.0]:.3f} h[t<=50]={means[50.0]:.3f} "
                  f"h_full={full_mean:.3f}  ({time.time()-t0:.0f}s)")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "delivered_mean_windows.csv", index=False)
    print(f"\nWrote {OUT / 'delivered_mean_windows.csv'}")
    # headline-cell summary for the main text (J=5)
    print("\n=== J=5, mu=100 delivered means by window ===")
    h5 = df[df.J == 5.0]
    print(h5[["arm", "alpha", "p_collapse", "h_mean_t<=1", "h_mean_t<=5",
              "h_mean_t<=14", "h_mean_t<=50", "h_mean_full"]].to_string(index=False))


if __name__ == "__main__":
    main()
