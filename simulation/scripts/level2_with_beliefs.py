"""Level 2 — Minimal model + replicator-mutator on K=6 belief types.

Adds belief-population dynamics on top of L1. Each simulation tracks
fractions f_k (k = 0..K-1) evolving by replicator-mutator with fitness
biased by the current coordination state m. Belief diversity (epistemic
commons) attenuates the passive stabilizer field:

    h_eff = h_passive * EC

where EC = entropy(f) / log(K) ∈ [0, 1].

Mechanism: at high coupling J, m saturates → fitness selects one belief
type → EC → 0 → h_eff → 0 → passive stabilizer is gone. Should AMPLIFY
the rigidity asymmetry compared to L1.

Same (J, mult) sweep as L1 for direct comparison.
"""

from pathlib import Path
import time

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "results" / "ablation"
OUT.mkdir(parents=True, exist_ok=True)

J_VALUES = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)
MULT_VALUES = (20, 30, 40, 50, 60, 70, 80, 90, 100)
N_SEEDS = 100
N_STEPS = 20_000

T = 2.0
XI = 0.5
DT = 0.05
EPS = 1e-6
INITIAL_M = 0.0
INITIAL_WEALTH = 100.0
COLLAPSE_WEALTH = 10.0
COLLAPSE_STREAK = 200
RIGIDITY_M = 0.9
FRAGMENTATION_M = 0.3

# Belief layer
K = 6                       # number of belief types
GAMMA = 1.0                 # m → fitness coupling strength
MU = 0.01                   # mutator rate
BELIEF_AXIS = (np.arange(K) - (K - 1) / 2) / ((K - 1) / 2)  # spans -1..1
LOG_K = np.log(K)


def run_cell_l2(J: float, mult: float, n_seeds: int = N_SEEDS,
                n_steps: int = N_STEPS, seed: int = 0):
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M, dtype=np.float64)
    wealth = np.full(n_seeds, INITIAL_WEALTH, dtype=np.float64)
    f = np.full((n_seeds, K), 1.0 / K, dtype=np.float64)  # uniform start

    collapsed = np.zeros(n_seeds, dtype=bool)
    collapse_step = np.full(n_seeds, -1, dtype=np.int32)
    collapse_type = np.zeros(n_seeds, dtype=np.int8)
    wealth_low_streak = np.zeros(n_seeds, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)

    for t in range(n_steps):
        # --- Belief replicator-mutator (uses current m) ---
        # fitness w_k(m) = exp(GAMMA * belief_axis_k * m)
        w = np.exp(GAMMA * BELIEF_AXIS[None, :] * m[:, None])  # (seeds, K)
        wf = f * w
        mean_w = wf.sum(axis=1, keepdims=True)
        f = wf / np.maximum(mean_w, 1e-30)
        f = (1 - MU) * f + MU / K
        f /= f.sum(axis=1, keepdims=True)  # numeric safety

        # Epistemic commons (entropy-based diversity)
        # EC = H(f) / log(K),  H = -sum f log f
        f_safe = np.maximum(f, 1e-30)
        H = -(f_safe * np.log(f_safe)).sum(axis=1)
        EC = H / LOG_K  # (seeds,) ∈ [0, 1]

        # --- Curie-Weiss SDE with EC-modulated stabilizer field ---
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h_passive = (wealth / 500.0) * 2.0
        h_eff = h_passive * EC
        f_drift = -m_c + np.tanh((J * m_c + h_eff) / T)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f_drift * DT + g * dW, -1 + EPS, 1 - EPS)

        # --- Economic feedback (unchanged from L1) ---
        employment = 0.5 + 0.5 * m
        income = mult * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        # --- Collapse detection ---
        is_low = wealth < COLLAPSE_WEALTH
        wealth_low_streak = np.where(is_low, wealth_low_streak + 1, 0)
        new_collapse = (wealth_low_streak >= COLLAPSE_STREAK) & (~collapsed)
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

    return {
        "collapsed": collapsed,
        "collapse_step": collapse_step,
        "collapse_type": collapse_type,
        "final_m": m,
        "final_wealth": wealth,
        "final_EC": EC,
    }


def main():
    rows = []
    t0 = time.time()
    for j_idx, J in enumerate(J_VALUES):
        for m_idx, mult in enumerate(MULT_VALUES):
            seed = 0xBEEF + j_idx * 1_000_003 + m_idx * 1009
            res = run_cell_l2(J, mult, seed=seed)
            for s in range(N_SEEDS):
                rows.append({
                    "level": "L2",
                    "J": J,
                    "mult": mult,
                    "seed_idx": s,
                    "collapsed": int(res["collapsed"][s]),
                    "collapse_step": int(res["collapse_step"][s]),
                    "collapse_type": int(res["collapse_type"][s]),
                    "final_m": float(res["final_m"][s]),
                    "final_wealth": float(res["final_wealth"][s]),
                    "final_EC": float(res["final_EC"][s]),
                })
    elapsed = time.time() - t0
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "level2_raw_results.csv", index=False)
    summary = (
        df.groupby(["J", "mult"])
          .agg(
              p_collapse=("collapsed", "mean"),
              p_rigidity=("collapse_type", lambda s: (s == 1).mean()),
              p_fragmentation=("collapse_type", lambda s: (s == 2).mean()),
              p_mixed=("collapse_type", lambda s: (s == 3).mean()),
              mean_final_m=("final_m", "mean"),
              mean_final_wealth=("final_wealth", "mean"),
              mean_final_EC=("final_EC", "mean"),
          )
          .reset_index()
    )
    summary.to_csv(OUT / "level2_cell_summary.csv", index=False)

    print(f"L2 sweep done in {elapsed:.1f}s — {len(df):,} rows, {len(summary)} cells")
    print("\nL2: P(collapse) by cell:")
    print(summary.pivot(index="J", columns="mult", values="p_collapse").round(2).to_string())
    print("\nL2: mean final EC by cell  (low EC = monoculture):")
    print(summary.pivot(index="J", columns="mult", values="mean_final_EC").round(2).to_string())


if __name__ == "__main__":
    main()
