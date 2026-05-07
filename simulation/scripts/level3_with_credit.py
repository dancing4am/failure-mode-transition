"""Level 3 — L2 + Minsky-Keen credit cycle.

Adds aggregate debt and a 4-phase machine (Expansion → Euphoria → MinskyMoment
→ Deleveraging → Expansion) on top of the L2 model. Debt accumulates with
optimism (m > 0); Minsky moments wipe wealth and spike noise. Same (J, mult)
sweep as L1/L2 for direct comparison.

Phases:
  0 Expansion       — baseline income; debt grows when m > 0
  1 Euphoria        — income × 1.1; debt accelerates
  2 MinskyMoment    — wealth × 0.7 (one-shot), D reset, income × 0.5, noise spike
  3 Deleveraging    — income × 0.85; debt paid down

Transitions (vectorised by boolean masks):
  0 → 1 if debt/wealth > 0.5 AND m > 0.3
  1 → 2 if debt/wealth > 1.0
  2 → 3 after MINSKY_DURATION = 100 steps in MinskyMoment (5 physical
        time units at dt = 0.05)
  3 → 0 if debt/wealth < 0.1
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

# Belief layer (from L2)
K = 6
GAMMA = 1.0
MU = 0.01
BELIEF_AXIS = (np.arange(K) - (K - 1) / 2) / ((K - 1) / 2)
LOG_K = np.log(K)

# Credit layer (Minsky-Keen)
ALPHA_DEBT = 0.5            # debt growth coefficient * m_pos * wealth
BETA_DEBT = 0.05            # debt amortisation rate
RATIO_TO_EUPHORIA = 0.5
M_TO_EUPHORIA = 0.3
RATIO_TO_MINSKY = 1.0
MINSKY_DURATION = 100      # = 5 physical time units at dt=0.05 (was 50 at dt=0.1)
RATIO_TO_EXPANSION = 0.1
MINSKY_WEALTH_HAIRCUT = 0.7
INCOME_FACTOR = np.array([1.0, 1.10, 0.50, 0.85])  # by phase
NOISE_FACTOR = np.array([1.0, 1.0, 1.5, 1.0])      # noise spike during MinskyMoment


def run_cell_l3(J: float, mult: float, n_seeds: int = N_SEEDS,
                n_steps: int = N_STEPS, seed: int = 0):
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M, dtype=np.float64)
    wealth = np.full(n_seeds, INITIAL_WEALTH, dtype=np.float64)
    f = np.full((n_seeds, K), 1.0 / K, dtype=np.float64)

    D = np.zeros(n_seeds, dtype=np.float64)
    phase = np.zeros(n_seeds, dtype=np.int8)
    phase_timer = np.zeros(n_seeds, dtype=np.int32)
    minsky_count = np.zeros(n_seeds, dtype=np.int32)

    collapsed = np.zeros(n_seeds, dtype=bool)
    collapse_step = np.full(n_seeds, -1, dtype=np.int32)
    collapse_type = np.zeros(n_seeds, dtype=np.int8)
    wealth_low_streak = np.zeros(n_seeds, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)

    for t in range(n_steps):
        # --- Belief layer ---
        w = np.exp(GAMMA * BELIEF_AXIS[None, :] * m[:, None])
        wf = f * w
        f = wf / np.maximum(wf.sum(axis=1, keepdims=True), 1e-30)
        f = (1 - MU) * f + MU / K
        f /= f.sum(axis=1, keepdims=True)
        f_safe = np.maximum(f, 1e-30)
        EC = -(f_safe * np.log(f_safe)).sum(axis=1) / LOG_K

        # --- Curie-Weiss SDE (Minsky moment spikes noise) ---
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h_passive = (wealth / 500.0) * 2.0
        h_eff = h_passive * EC
        f_drift = -m_c + np.tanh((J * m_c + h_eff) / T)
        g = XI * NOISE_FACTOR[phase] * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f_drift * DT + g * dW, -1 + EPS, 1 - EPS)

        # --- Credit / Minsky layer ---
        D = np.maximum(0.0, D + (ALPHA_DEBT * np.maximum(0.0, m) * wealth - BETA_DEBT * D) * DT)
        debt_ratio = D / np.maximum(wealth, 1.0)

        to_euphoria = (phase == 0) & (debt_ratio > RATIO_TO_EUPHORIA) & (m > M_TO_EUPHORIA)
        to_minsky = (phase == 1) & (debt_ratio > RATIO_TO_MINSKY)
        to_delev = (phase == 2) & (phase_timer > MINSKY_DURATION)
        to_expand = (phase == 3) & (debt_ratio < RATIO_TO_EXPANSION)

        # Apply one-shot Minsky wealth haircut + debt wipe
        wealth = np.where(to_minsky, wealth * MINSKY_WEALTH_HAIRCUT, wealth)
        D = np.where(to_minsky, 0.0, D)
        minsky_count = minsky_count + to_minsky.astype(np.int32)

        new_phase = phase.copy()
        new_phase = np.where(to_euphoria, 1, new_phase)
        new_phase = np.where(to_minsky, 2, new_phase)
        new_phase = np.where(to_delev, 3, new_phase)
        new_phase = np.where(to_expand, 0, new_phase)
        phase_changed = (new_phase != phase)
        phase = new_phase
        phase_timer = np.where(phase_changed, 0, phase_timer + 1)

        # --- Economic feedback (phase-modulated income) ---
        employment = 0.5 + 0.5 * m
        income = mult * employment * INCOME_FACTOR[phase]
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
        "minsky_count": minsky_count,
    }


def main():
    rows = []
    t0 = time.time()
    for j_idx, J in enumerate(J_VALUES):
        for m_idx, mult in enumerate(MULT_VALUES):
            seed = 0xCAFE + j_idx * 1_000_003 + m_idx * 1009
            res = run_cell_l3(J, mult, seed=seed)
            for s in range(N_SEEDS):
                rows.append({
                    "level": "L3",
                    "J": J,
                    "mult": mult,
                    "seed_idx": s,
                    "collapsed": int(res["collapsed"][s]),
                    "collapse_step": int(res["collapse_step"][s]),
                    "collapse_type": int(res["collapse_type"][s]),
                    "final_m": float(res["final_m"][s]),
                    "final_wealth": float(res["final_wealth"][s]),
                    "final_EC": float(res["final_EC"][s]),
                    "minsky_count": int(res["minsky_count"][s]),
                })
    elapsed = time.time() - t0
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "level3_raw_results.csv", index=False)
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
              mean_minsky_count=("minsky_count", "mean"),
          )
          .reset_index()
    )
    summary.to_csv(OUT / "level3_cell_summary.csv", index=False)

    print(f"L3 sweep done in {elapsed:.1f}s — {len(df):,} rows, {len(summary)} cells")
    print("\nL3: P(collapse) by cell:")
    print(summary.pivot(index="J", columns="mult", values="p_collapse").round(2).to_string())
    print("\nL3: P(rigidity) by cell:")
    print(summary.pivot(index="J", columns="mult", values="p_rigidity").round(2).to_string())
    print("\nL3: mean Minsky moments per run:")
    print(summary.pivot(index="J", columns="mult", values="mean_minsky_count").round(1).to_string())


if __name__ == "__main__":
    main()
