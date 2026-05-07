"""Convergence check for the minimal model integrator.

The Results section claims that the asymmetry survives at smaller dt.
This script verifies. Re-runs two diagnostic cells at dt = 0.1, 0.05,
0.025, holding the physical horizon (dt * n_steps = 1000) fixed.

Cells:
    A. J = 5.0, mult = 100  — residual rigidity regime (the key claim:
       the 25% / 85% rigidity number is not an integration artifact).
    B. J = 0.5, mult = 100  — fragmentation-cured corner (sanity check
       that the easy direction is also stable under refinement).

Output: results/minimal_model/convergence_check.csv  +  console summary.
"""

from __future__ import annotations

from pathlib import Path
import time

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "results" / "minimal_model"
OUT.mkdir(parents=True, exist_ok=True)

# Parameters mirror minimal_model.py exactly.
T = 2.0
XI = 0.5
EPS = 1e-6
INITIAL_M = 0.0
INITIAL_WEALTH = 100.0
COLLAPSE_WEALTH = 10.0
COLLAPSE_STREAK_PHYS_TIME = 10.0  # physical time units (= 100 steps at dt=0.1)
RIGIDITY_M = 0.9
FRAGMENTATION_M = 0.3

PHYSICAL_HORIZON = 1000.0   # 10,000 steps at dt = 0.1
N_SEEDS = 100

CELLS = [
    {"name": "A_residual_rigidity", "J": 5.0, "mult": 100.0},
    {"name": "B_fragmentation_cured", "J": 0.5, "mult": 100.0},
]
DT_VALUES = (0.1, 0.05, 0.025)


def run_cell(J: float, mult: float, dt: float, n_seeds: int, seed: int) -> dict:
    n_steps = int(round(PHYSICAL_HORIZON / dt))
    streak_threshold = int(round(COLLAPSE_STREAK_PHYS_TIME / dt))
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M, dtype=np.float64)
    wealth = np.full(n_seeds, INITIAL_WEALTH, dtype=np.float64)
    collapsed = np.zeros(n_seeds, dtype=bool)
    collapse_type = np.zeros(n_seeds, dtype=np.int8)
    wealth_low_streak = np.zeros(n_seeds, dtype=np.int32)
    sqrt_dt = np.sqrt(dt)

    for _ in range(n_steps):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = (wealth / 500.0) * 2.0
        f_drift = -m_c + np.tanh((J * m_c + h) / T)
        g = XI * np.sqrt(1 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f_drift * dt + g * dW, -1 + EPS, 1 - EPS)

        employment = 0.5 + 0.5 * m
        income = mult * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * dt)

        is_low = wealth < COLLAPSE_WEALTH
        wealth_low_streak = np.where(is_low, wealth_low_streak + 1, 0)
        new_collapse = (wealth_low_streak >= streak_threshold) & (~collapsed)
        if new_collapse.any():
            collapsed |= new_collapse
            abs_m = np.abs(m)
            rigid = new_collapse & (abs_m > RIGIDITY_M)
            frag = new_collapse & (abs_m < FRAGMENTATION_M)
            mixed = new_collapse & ~(rigid | frag)
            collapse_type[rigid] = 1
            collapse_type[frag] = 2
            collapse_type[mixed] = 3

    n_coll = int(collapsed.sum())
    return {
        "p_collapse": collapsed.mean(),
        "n_collapsed": n_coll,
        "p_rigidity_of_total": (collapse_type == 1).mean(),
        "p_rigidity_of_collapsed": (collapse_type == 1).sum() / n_coll if n_coll else float("nan"),
        "p_fragmentation_of_collapsed": (collapse_type == 2).sum() / n_coll if n_coll else float("nan"),
        "n_steps": n_steps,
    }


def main():
    rows = []
    t0 = time.time()
    for cell in CELLS:
        for dt in DT_VALUES:
            seed = 1009 * int(cell["J"] * 10) + int(cell["mult"]) + int(round(1 / dt))
            res = run_cell(cell["J"], cell["mult"], dt, N_SEEDS, seed)
            rows.append({
                "cell": cell["name"],
                "J": cell["J"],
                "mult": cell["mult"],
                "dt": dt,
                "n_steps": res["n_steps"],
                "p_collapse": res["p_collapse"],
                "n_collapsed": res["n_collapsed"],
                "p_rigidity_of_collapsed": res["p_rigidity_of_collapsed"],
                "p_fragmentation_of_collapsed": res["p_fragmentation_of_collapsed"],
            })
            print(f"{cell['name']:>30s}  dt={dt:.3f}  n_steps={res['n_steps']:>6}  "
                  f"P(coll)={res['p_collapse']:.3f}  rig_share="
                  f"{res['p_rigidity_of_collapsed']:.3f}")
    elapsed = time.time() - t0
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "convergence_check.csv", index=False)
    print(f"\nDone in {elapsed:.1f}s — wrote convergence_check.csv")

    # --- Verdict ---------------------------------------------------------
    print("\nConvergence verdict:")
    for cell in CELLS:
        sub = df[df["cell"] == cell["name"]].sort_values("dt", ascending=False)
        ps = sub["p_collapse"].values
        rs = sub["p_rigidity_of_collapsed"].values
        p_spread = float(ps.max() - ps.min())
        # 100 seeds → SE on a 0.25 collapse rate is ~0.043
        sampling_se = 2 * np.sqrt(0.25 * 0.75 / N_SEEDS)
        verdict = "CONVERGED" if p_spread <= sampling_se else "NOT CONVERGED"
        print(f"  {cell['name']}: P(coll) range = [{ps.min():.3f}, {ps.max():.3f}], "
              f"spread = {p_spread:.3f}, 2·SE bound = {sampling_se:.3f} -> {verdict}")
        if not np.isnan(rs).all():
            r_valid = rs[~np.isnan(rs)]
            r_spread = float(r_valid.max() - r_valid.min()) if len(r_valid) > 1 else 0.0
            print(f"    rigidity share range = [{r_valid.min():.3f}, {r_valid.max():.3f}], "
                  f"spread = {r_spread:.3f}")


if __name__ == "__main__":
    main()
