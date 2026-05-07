"""
Active stabilizer simulation.

Tests whether stabilizers whose *intervention magnitude* scales with stress
or coupling can eliminate the high-J residual that bounded passive
stabilizers cannot reach (per Theorem 1 of the paper).

Four conditions, same (J, mult) grid and matched seeds as the minimal model:

1. passive            h = 2*(W/500)                               (baseline)
2. active_stress      h = base + alpha*max(0,-m)*J                (Fed analogue)
3. active_coupling    h = base*(1 + alpha*J)                      (h ~ J → eta ~ 1/J)
4. active_quadratic   h = base*(1 + alpha*(J/5)**2)               (h ~ J^2 → eta ~ const)

Plus an alpha sweep on the best-performing condition at the high-J band.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import time

import numpy as np
import pandas as pd

# Constants identical to minimal_model.py
J_VALUES = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)
MULT_VALUES = (20, 30, 40, 50, 60, 70, 80, 90, 100)
HIGH_J_BAND = (4.0, 4.5, 5.0)
ALPHA_SWEEP_VALUES = (0.1, 0.25, 0.5, 1.0, 2.0, 5.0)
N_SEEDS = 100
N_STEPS = 20_000

T_TEMP = 2.0
XI = 0.5
DT = 0.05
INITIAL_M = 0.0
INITIAL_WEALTH = 100.0
COLLAPSE_WEALTH = 10.0
COLLAPSE_STREAK = 200
RIGIDITY_M = 0.9
FRAGMENTATION_M = 0.3
EPS = 1e-6

SEED_BASE = 0xC0DE


# Stabilizer functions — vectorized over seeds. m and wealth are arrays.
def h_passive(W, m, J, alpha):
    return (W / 500.0) * 2.0


def h_active_stress(W, m, J, alpha):
    base = (W / 500.0) * 2.0
    stress = np.maximum(0.0, -m)
    return base + alpha * stress * J


def h_active_coupling(W, m, J, alpha):
    base = (W / 500.0) * 2.0
    return base * (1.0 + alpha * J)


def h_active_quadratic(W, m, J, alpha):
    base = (W / 500.0) * 2.0
    return base * (1.0 + alpha * (J / 5.0) ** 2)


CONDITIONS = {
    "passive":           (h_passive,           0.0),
    "active_stress":     (h_active_stress,     2.0),
    "active_coupling":   (h_active_coupling,   1.0),
    "active_quadratic":  (h_active_quadratic,  1.0),
}


@dataclass
class CellResult:
    J: float
    mult: float
    condition: str
    alpha: float
    collapsed: np.ndarray
    collapse_step: np.ndarray
    collapse_type: np.ndarray
    final_m: np.ndarray
    final_wealth: np.ndarray


def run_cell(J, mult, h_fn, alpha, n_seeds=N_SEEDS, n_steps=N_STEPS, seed=0):
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M, dtype=np.float64)
    wealth = np.full(n_seeds, INITIAL_WEALTH, dtype=np.float64)
    collapsed = np.zeros(n_seeds, dtype=bool)
    collapse_step = np.full(n_seeds, -1, dtype=np.int32)
    collapse_type = np.zeros(n_seeds, dtype=np.int8)
    streak = np.zeros(n_seeds, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)

    for t in range(n_steps):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = h_fn(wealth, m_c, J, alpha)
        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)

        employment = 0.5 + 0.5 * m
        income = mult * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        is_low = wealth < COLLAPSE_WEALTH
        streak = np.where(is_low, streak + 1, 0)
        new_collapse = (streak >= COLLAPSE_STREAK) & (~collapsed)
        if new_collapse.any():
            collapsed |= new_collapse
            collapse_step[new_collapse] = t
            abs_m = np.abs(m)
            rigid = new_collapse & (abs_m > RIGIDITY_M)
            frag = new_collapse & (abs_m < FRAGMENTATION_M)
            mixed = new_collapse & ~(rigid | frag)
            collapse_type[rigid] = 1
            collapse_type[frag] = 2
            collapse_type[mixed] = 3

    return CellResult(J, mult, "", alpha, collapsed, collapse_step,
                      collapse_type, m, wealth)


def sweep_condition(condition_name, h_fn, alpha, j_values=J_VALUES,
                    mult_values=MULT_VALUES, seed_base=SEED_BASE):
    cells = []
    j_list = list(j_values)
    mult_list = list(mult_values)
    for j_idx, J in enumerate(j_list):
        for m_idx, mult in enumerate(mult_list):
            seed = seed_base + j_idx * 1_000_003 + m_idx * 1009
            res = run_cell(J, mult, h_fn, alpha, seed=seed)
            res.condition = condition_name
            res.alpha = alpha
            cells.append(res)
    return cells


def cells_to_long_df(cells):
    rows = []
    for c in cells:
        for s in range(len(c.collapsed)):
            rows.append({
                "J": c.J,
                "mu": c.mult,
                "seed": s,
                "collapsed": int(c.collapsed[s]),
                "collapse_type": int(c.collapse_type[s]),
                "collapse_step": int(c.collapse_step[s]),
                "final_m": float(c.final_m[s]),
                "final_W": float(c.final_wealth[s]),
                "condition": c.condition,
                "alpha": float(c.alpha),
            })
    return pd.DataFrame(rows)


def cell_summary_df(cells):
    rows = []
    for c in cells:
        n = len(c.collapsed)
        rows.append({
            "J": c.J,
            "mu": c.mult,
            "condition": c.condition,
            "alpha": c.alpha,
            "p_collapse": c.collapsed.mean(),
            "p_rigidity": (c.collapse_type == 1).mean(),
            "p_fragmentation": (c.collapse_type == 2).mean(),
            "p_mixed": (c.collapse_type == 3).mean(),
            "n_collapsed": int(c.collapsed.sum()),
            "rigidity_share_of_collapsed": (
                float((c.collapse_type == 1).sum() / max(1, c.collapsed.sum()))
            ),
            "mean_final_m": float(c.final_m.mean()),
            "mean_final_wealth": float(c.final_wealth.mean()),
        })
    return pd.DataFrame(rows)


def main():
    OUT = Path(__file__).resolve().parents[1] / "results" / "active_stabilizer"
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Output dir: {OUT}")

    # ---- Main 4-condition sweep on full grid -----------------------------
    summaries = {}
    for cond, (h_fn, alpha) in CONDITIONS.items():
        t0 = time.time()
        cells = sweep_condition(cond, h_fn, alpha)
        elapsed = time.time() - t0
        long_df = cells_to_long_df(cells)
        summary = cell_summary_df(cells)
        long_df.to_csv(OUT / f"{cond}_results.csv", index=False)
        summary.to_csv(OUT / f"{cond}_summary.csv", index=False)
        summaries[cond] = summary

        # Headline numbers
        high_J_band_at_max_mu = summary[
            summary.J.isin(HIGH_J_BAND) & (summary.mu == 100)
        ]
        print(f"\n[{cond}, alpha={alpha}] elapsed={elapsed:.1f}s")
        print(f"  P(coll) at high-J band, mu=100: "
              f"{high_J_band_at_max_mu.p_collapse.mean():.3f}")
        print(f"  Per-cell J=4.0/4.5/5.0 at mu=100: "
              f"{high_J_band_at_max_mu.sort_values('J').p_collapse.values}")

    # ---- Alpha sweep on the most promising condition ---------------------
    # Run alpha sweep for stress, coupling, and quadratic at high-J band
    # full mu grid, so we can build alpha_threshold.png comparisons.
    alpha_rows = []
    for cond_name, h_fn in [
        ("active_stress", h_active_stress),
        ("active_coupling", h_active_coupling),
        ("active_quadratic", h_active_quadratic),
    ]:
        for alpha in ALPHA_SWEEP_VALUES:
            t0 = time.time()
            cells = sweep_condition(
                cond_name, h_fn, alpha,
                j_values=HIGH_J_BAND, mult_values=MULT_VALUES,
            )
            elapsed = time.time() - t0
            summary = cell_summary_df(cells)
            print(f"  alpha-sweep {cond_name} alpha={alpha} "
                  f"P(coll@J>=4,mu=100)={summary[(summary.mu==100)].p_collapse.mean():.3f} "
                  f"({elapsed:.1f}s)")
            for _, row in summary.iterrows():
                alpha_rows.append(row.to_dict())

    alpha_df = pd.DataFrame(alpha_rows)
    alpha_df.to_csv(OUT / "alpha_sweep_results.csv", index=False)
    print(f"\nWrote {OUT / 'alpha_sweep_results.csv'} ({len(alpha_df)} rows)")


if __name__ == "__main__":
    main()
