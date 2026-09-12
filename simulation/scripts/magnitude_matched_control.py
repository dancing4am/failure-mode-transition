"""
Magnitude-matched constant-field controls.

Question: within static sweeps (J constant per run), does the *form* of the
stabilizer response matter beyond the delivered *magnitude* of the field?

1a. Constancy tabulation: delivered field per variant at the wealth
    equilibrium W*(mu=100, m=+1) and at W0=100, per J.
    Wealth ODE: dW/dt = mu*(1+m)/2 - (30 + 10*(W/500)).
    At mu=100, m=+1: income = 100, so 100 = 30 + W*/50 -> W* = 3500.

1b. Stress variant (alpha=2.0, mu=100) rerun at J in {4,4.5,5,7,10,15,20},
    n=2000 matched seeds (seed = 0xE1_0001 + j_idx*1_000_003), recording the
    realized field h in the decision window t <= 50 (steps <= 1000).
    Then constant-field controls with the same seeds:
      (i)  h = per-J realized window-mean  (constant, no wealth feedback)
      (ii) h = per-J realized full-run mean
    Compare P(collapse) stress vs matched-constant, bootstrap 95% CIs.

Outputs to simulation/results/magnitude_matched/:
    constancy_table.csv, stress_field_stats.csv, collapse_comparison.csv,
    SUMMARY.md
"""

from __future__ import annotations

from pathlib import Path
import time

import numpy as np
import pandas as pd

# Integrator constants (identical to minimal_model.py / active_stabilizer.py)
T_TEMP = 2.0
XI = 0.5
DT = 0.05
N_STEPS = 20_000
INITIAL_M = 0.0
INITIAL_WEALTH = 100.0
COLLAPSE_WEALTH = 10.0
COLLAPSE_STREAK = 200
EPS = 1e-6

MU = 100.0
ALPHA_STRESS = 2.0

J_VALUES_1A = (3.0, 4.0, 4.5, 5.0, 7.0, 10.0, 15.0, 20.0)
J_VALUES_1B = (4.0, 4.5, 5.0, 7.0, 10.0, 15.0, 20.0)

N_SEEDS = 2000
SEED_BASE = 0xE1_0001            # seed base: 0xE1_0001 + j_idx*1_000_003
WINDOW_STEPS = 1000              # t <= 50  <=>  step index <= 1000 (inclusive)

N_BOOT = 10_000
BOOT_SEED = 0xB007

OUT = Path(__file__).resolve().parents[1] / "results" / "magnitude_matched"


# ---------------------------------------------------------------------------
# 1a — constancy tabulation (analytic, no simulation)
# ---------------------------------------------------------------------------

def wealth_equilibrium(mu: float, m: float) -> float:
    """Solve mu*(1+m)/2 = 30 + 10*(W/500) for W."""
    income = mu * (1.0 + m) / 2.0
    return (income - 30.0) * 50.0


def delivered_field(variant: str, W: float, J: float) -> float:
    base = 2.0 * (W / 500.0)
    if variant == "passive":
        return base
    if variant == "coupling_a1":
        return base * (1.0 + 1.0 * J)
    if variant == "quadratic_a1":
        return base * (1.0 + 1.0 * (J / 5.0) ** 2)
    if variant == "fixed_2.4":
        return 2.4
    if variant == "fixed_5.0":
        return 5.0
    raise ValueError(variant)


def build_constancy_table() -> pd.DataFrame:
    W_star = wealth_equilibrium(MU, +1.0)   # 3500
    rows = []
    for variant in ("passive", "coupling_a1", "quadratic_a1",
                    "fixed_2.4", "fixed_5.0"):
        for J in J_VALUES_1A:
            rows.append({
                "variant": variant,
                "J": J,
                "h_at_Wstar_3500": delivered_field(variant, W_star, J),
                "h_at_W0_100": delivered_field(variant, INITIAL_WEALTH, J),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 1b — simulation
# ---------------------------------------------------------------------------

def run_stress(J: float, seed: int, n_seeds: int = N_SEEDS):
    """Stress variant; record realized-field stats in window and full run."""
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M)
    wealth = np.full(n_seeds, INITIAL_WEALTH)
    collapsed = np.zeros(n_seeds, dtype=bool)
    streak = np.zeros(n_seeds, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)

    win_sum = 0.0
    win_n = 0
    win_max = -np.inf
    win_active = 0            # (seed,step) cells with stress term active (m<0)
    full_sum = 0.0
    full_n = 0
    full_active = 0

    for t in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        base = 2.0 * (wealth / 500.0)
        stress = np.maximum(0.0, -m_c)
        h = base + ALPHA_STRESS * stress * J

        full_sum += h.sum()
        full_n += n_seeds
        full_active += int((stress > 0).sum())
        if t <= WINDOW_STEPS:
            win_sum += h.sum()
            win_n += n_seeds
            win_max = max(win_max, float(h.max()))
            win_active += int((stress > 0).sum())

        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)

        employment = 0.5 + 0.5 * m
        income = MU * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        is_low = wealth < COLLAPSE_WEALTH
        streak = np.where(is_low, streak + 1, 0)
        collapsed |= streak >= COLLAPSE_STREAK

    stats = {
        "window_mean_h": win_sum / win_n,
        "window_max_h": win_max,
        "window_frac_stress_active": win_active / win_n,
        "full_mean_h": full_sum / full_n,
        "full_frac_stress_active": full_active / full_n,
    }
    return collapsed, stats


def run_constant(J: float, h_const: float, seed: int, n_seeds: int = N_SEEDS):
    """Constant field h_const: no wealth feedback, no stress term."""
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M)
    wealth = np.full(n_seeds, INITIAL_WEALTH)
    collapsed = np.zeros(n_seeds, dtype=bool)
    streak = np.zeros(n_seeds, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)

    for t in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        f = -m_c + np.tanh((J * m_c + h_const) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)

        employment = 0.5 + 0.5 * m
        income = MU * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        is_low = wealth < COLLAPSE_WEALTH
        streak = np.where(is_low, streak + 1, 0)
        collapsed |= streak >= COLLAPSE_STREAK

    return collapsed


def boot_ci(x: np.ndarray, rng: np.random.Generator, n_boot: int = N_BOOT):
    """Percentile 95% CI on the mean of a boolean array."""
    n = len(x)
    idx = rng.integers(0, n, size=(n_boot, n))
    means = x[idx].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def paired_boot_diff_ci(a: np.ndarray, b: np.ndarray,
                        rng: np.random.Generator, n_boot: int = N_BOOT):
    """Paired (matched-seed) percentile 95% CI on mean(a) - mean(b)."""
    n = len(a)
    d = a.astype(float) - b.astype(float)
    idx = rng.integers(0, n, size=(n_boot, n))
    diffs = d[idx].mean(axis=1)
    return float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    # ---- 1a -------------------------------------------------------------
    constancy = build_constancy_table()
    constancy.to_csv(OUT / "constancy_table.csv", index=False)
    print("1a constancy table written.")
    print(constancy.to_string(index=False))

    # ---- 1b -------------------------------------------------------------
    boot_rng = np.random.default_rng(BOOT_SEED)
    field_rows = []
    comp_rows = []
    for j_idx, J in enumerate(J_VALUES_1B):
        seed = SEED_BASE + j_idx * 1_000_003
        t0 = time.time()
        coll_stress, stats = run_stress(J, seed)
        h_win = stats["window_mean_h"]
        h_full = stats["full_mean_h"]
        coll_cwin = run_constant(J, h_win, seed)
        coll_cfull = run_constant(J, h_full, seed)
        elapsed = time.time() - t0

        field_rows.append({"J": J, "seed": seed, **stats})

        p_s = coll_stress.mean()
        p_w = coll_cwin.mean()
        p_f = coll_cfull.mean()
        ci_s = boot_ci(coll_stress, boot_rng)
        ci_w = boot_ci(coll_cwin, boot_rng)
        ci_f = boot_ci(coll_cfull, boot_rng)
        dw_ci = paired_boot_diff_ci(coll_stress, coll_cwin, boot_rng)
        df_ci = paired_boot_diff_ci(coll_stress, coll_cfull, boot_rng)
        comp_rows.append({
            "J": J,
            "h_const_window": h_win,
            "h_const_full": h_full,
            "p_stress": p_s, "p_stress_lo": ci_s[0], "p_stress_hi": ci_s[1],
            "p_const_window": p_w,
            "p_const_window_lo": ci_w[0], "p_const_window_hi": ci_w[1],
            "p_const_full": p_f,
            "p_const_full_lo": ci_f[0], "p_const_full_hi": ci_f[1],
            "diff_stress_minus_window": p_s - p_w,
            "diff_window_lo": dw_ci[0], "diff_window_hi": dw_ci[1],
            "diff_stress_minus_full": p_s - p_f,
            "diff_full_lo": df_ci[0], "diff_full_hi": df_ci[1],
        })
        print(f"J={J:5.1f}  h_win={h_win:.4f}  h_full={h_full:.4f}  "
              f"P(stress)={p_s:.4f}  P(cwin)={p_w:.4f}  P(cfull)={p_f:.4f}  "
              f"frac_active_win={stats['window_frac_stress_active']:.4f}  "
              f"({elapsed:.1f}s)")

    field_df = pd.DataFrame(field_rows)
    comp_df = pd.DataFrame(comp_rows)
    field_df.to_csv(OUT / "stress_field_stats.csv", index=False)
    comp_df.to_csv(OUT / "collapse_comparison.csv", index=False)
    print(f"\nWrote CSVs to {OUT}")
    return constancy, field_df, comp_df


if __name__ == "__main__":
    main()
