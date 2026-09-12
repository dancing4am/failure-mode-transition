"""Third arm for the informative magnitude-matched control: pure passive (alpha=0).

The informative control (magnitude_matched_informative.py) compares, at the
informative operating point (b=2, alpha=0.05, mu=100, J in {4, 4.5, 5}), the
stress form against a pure constant matched to the stress form's
full-trajectory delivered mean. This script adds the THIRD arm needed to
decide whether the stress arm's rectified term max(0,-m)*J does anything at
alpha=0.05: the pure passive field h = 2*(W/500) (alpha=0), on the SAME
matched seeds and same n=2000.

It reports, per J in {4.0, 4.5, 5.0}:
  P(passive, alpha=0), P(stress, alpha=0.05), P(constant),
  paired differences stress-vs-passive and stress-vs-constant with 95% CIs
  (10,000-rep paired bootstrap), and the stress-term activation at alpha=0.05:
  the share of DECISION-WINDOW (t <= 50) seed-steps with max(0,-m) > 0 and the
  mean added field alpha*max(0,-m)*J it contributes over those active steps.

Reuses the integrator and pairing of magnitude_matched_informative.py verbatim
(imported), so the stress arm here is bit-identical to informative_control.csv.

Output: simulation/results/magnitude_matched/informative_passive_arm.csv
"""
from __future__ import annotations

from pathlib import Path
import time

import numpy as np
import pandas as pd

from magnitude_matched_informative import (
    T_TEMP, XI, DT, N_STEPS, INITIAL_M, INITIAL_WEALTH, COLLAPSE_WEALTH,
    COLLAPSE_STREAK, EPS, MU, N_SEEDS, N_BOOT, BOOT_SEED,
    run_stress, run_constant, paired_boot_diff_ci, _seed_for_J,
)

B_STAR, ALPHA_STAR = 2.0, 0.05      # the informative operating point
CELLS_J = (5.0, 4.0, 4.5)
DECISION_WINDOW_STEPS = int(round(50.0 / DT))    # t <= 50 -> first 1000 steps

OUT = Path(__file__).resolve().parents[1] / "results" / "magnitude_matched"


def run_passive(J: float, seed: int, n_seeds: int = N_SEEDS):
    """Pure passive field h = 2*(W/500) (alpha=0). Same RNG call pattern as
    run_stress/run_constant so the seed gives identical per-step noise draws."""
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M)
    wealth = np.full(n_seeds, INITIAL_WEALTH)
    collapsed = np.zeros(n_seeds, dtype=bool)
    streak = np.zeros(n_seeds, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)
    for _ in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = 2.0 * (wealth / 500.0)
        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)
        income = MU * (0.5 + 0.5 * m)
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)
        streak = np.where(wealth < COLLAPSE_WEALTH, streak + 1, 0)
        collapsed |= streak >= COLLAPSE_STREAK
    return collapsed


def stress_activation(J: float, b: float, alpha: float, seed: int,
                      n_seeds: int = N_SEEDS, window_steps: int = DECISION_WINDOW_STEPS):
    """Measure the rectified stress term over the decision window (t <= 50).

    Returns (active_frac, mean_added_active) where active_frac is the share of
    (seed, step) pairs in the window with max(0,-m) > 0, and mean_added_active
    is the mean of alpha*max(0,-m)*J over those active pairs (0 if none).
    """
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M)
    wealth = np.full(n_seeds, INITIAL_WEALTH)
    streak = np.zeros(n_seeds, dtype=np.int32)
    collapsed = np.zeros(n_seeds, dtype=bool)
    sqrt_dt = np.sqrt(DT)
    active_pairs = 0
    total_pairs = 0
    added_sum = 0.0
    for t in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        stress = np.maximum(0.0, -m_c)
        h = b * (wealth / 500.0) + alpha * stress * J
        if t < window_steps:
            act = stress > 0.0
            active_pairs += int(act.sum())
            total_pairs += n_seeds
            added_sum += float((alpha * stress * J)[act].sum())
        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)
        income = MU * (0.5 + 0.5 * m)
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)
        streak = np.where(wealth < COLLAPSE_WEALTH, streak + 1, 0)
        collapsed |= streak >= COLLAPSE_STREAK
    active_frac = active_pairs / total_pairs if total_pairs else 0.0
    mean_added_active = added_sum / active_pairs if active_pairs else 0.0
    return active_frac, mean_added_active


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    boot_rng = np.random.default_rng(BOOT_SEED)
    rows = []
    for J in CELLS_J:
        seed = _seed_for_J(J)
        t0 = time.time()
        coll_stress, h_bar = run_stress(J, B_STAR, ALPHA_STAR, seed)
        coll_passive = run_passive(J, seed)
        coll_const = run_constant(J, h_bar, seed)
        p_s = float(coll_stress.mean())
        p_p = float(coll_passive.mean())
        p_c = float(coll_const.mean())
        d_sp = p_s - p_p
        lo_sp, hi_sp = paired_boot_diff_ci(coll_stress, coll_passive, boot_rng)
        d_sc = p_s - p_c
        lo_sc, hi_sc = paired_boot_diff_ci(coll_stress, coll_const, boot_rng)
        act_frac, mean_added = stress_activation(J, B_STAR, ALPHA_STAR, seed)
        sp_ci_contains_0 = (lo_sp <= 0.0 <= hi_sp)
        rows.append(dict(
            J=J, mu=MU, b=B_STAR, alpha=ALPHA_STAR, h_bar_full=h_bar,
            p_passive=p_p, p_stress=p_s, p_const=p_c,
            diff_stress_minus_passive=d_sp, dsp_lo=lo_sp, dsp_hi=hi_sp,
            diff_stress_minus_const=d_sc, dsc_lo=lo_sc, dsc_hi=hi_sc,
            stress_active_frac_window=act_frac,
            mean_added_field_active=mean_added,
            stress_vs_passive_CI_contains_0=sp_ci_contains_0, n=N_SEEDS,
        ))
        print(f"J={J:.1f}  P(passive a=0)={p_p:.4f}  P(stress a=0.05)={p_s:.4f}  "
              f"P(const)={p_c:.4f}")
        print(f"   stress-passive diff={d_sp:+.4f} CI=[{lo_sp:+.4f},{hi_sp:+.4f}] "
              f"(contains 0: {sp_ci_contains_0})")
        print(f"   stress-const   diff={d_sc:+.4f} CI=[{lo_sc:+.4f},{hi_sc:+.4f}]")
        print(f"   stress-term active on {act_frac*100:.3f}% of decision-window steps; "
              f"mean added field there = {mean_added:.4f}  ({time.time()-t0:.0f}s)")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "informative_passive_arm.csv", index=False)
    print(f"\nWrote {OUT / 'informative_passive_arm.csv'}")


if __name__ == "__main__":
    main()
