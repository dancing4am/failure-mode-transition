"""
Sensitivity sweep — robustness check on noise prescription and h(W)
functional form. Runs three variants of the minimal model and compares
the high-J residual collapse rate against the canonical L1 baseline.

Variants:
  baseline   — multiplicative noise xi*sqrt(1-m^2), h(W) = 2*(W/500)
  additive   — additive noise xi (no m-dependence), h(W) unchanged
  log_h      — multiplicative noise unchanged, h(W) = c_log * log(1 + W/100)
               with c_log chosen so h matches baseline at W=100 (initial)
  const_h    — multiplicative noise unchanged, h(W) = c_const (constant)
               c_const = 2*(100/500) = 0.4 (= h at initial wealth)

For tractability we sweep the high-J band only (J in {3.5, 4, 4.5, 5}) at
the protective-margin band (mu in {60, 80, 100}), 100 seeds per cell.
12 cells x 4 variants x 100 seeds = 4,800 runs. ~3 minutes on CPU.

Output:
  results/sensitivity/{baseline,additive,log_h,const_h}_summary.csv
  results/sensitivity/comparison_summary.md
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "sensitivity"
OUT.mkdir(parents=True, exist_ok=True)

# -- Sweep grid (high-J band; protective margins) --
J_VALUES = (3.5, 4.0, 4.5, 5.0)
MULT_VALUES = (60, 80, 100)
N_SEEDS = 100
N_STEPS = 20_000

# -- Fixed model constants (mirror minimal_model.py) --
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

# -- h(W) variants --
def h_linear(wealth: np.ndarray) -> np.ndarray:
    return 2.0 * (wealth / 500.0)

def h_log(wealth: np.ndarray) -> np.ndarray:
    # Calibrate so that h(W=100) matches linear: 2*(100/500) = 0.4
    # h_log(W) = c * log(1 + W/100); at W=100, h = c*log(2) = 0.4 => c = 0.577
    c_log = 0.4 / np.log(2.0)
    return c_log * np.log1p(wealth / 100.0)

def h_const(wealth: np.ndarray) -> np.ndarray:
    return np.full_like(wealth, 0.4)  # constant at h(W=100) value

H_FUNCS = {
    "baseline": h_linear,
    "additive": h_linear,
    "log_h":    h_log,
    "const_h":  h_const,
}


@dataclass
class CellResult:
    variant: str
    J: float
    mult: float
    p_collapse: float
    rigidity_share: float
    fragmentation_share: float
    se: float


def run_cell(variant: str, J: float, mult: float, h_func: Callable, seed: int) -> CellResult:
    rng = np.random.default_rng(seed)
    m = np.full(N_SEEDS, INITIAL_M, dtype=np.float64)
    wealth = np.full(N_SEEDS, INITIAL_WEALTH, dtype=np.float64)
    collapsed = np.zeros(N_SEEDS, dtype=bool)
    collapse_type = np.zeros(N_SEEDS, dtype=np.int8)
    streak = np.zeros(N_SEEDS, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)

    for _t in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = h_func(wealth)
        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        if variant == "additive":
            g = np.full_like(m_c, XI)        # additive: g = xi (constant)
        else:
            g = XI * np.sqrt(1.0 - m_c * m_c)  # multiplicative
        dW = sqrt_dt * rng.standard_normal(N_SEEDS)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)

        employment = 0.5 + 0.5 * m
        income = mult * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        is_low = wealth < COLLAPSE_WEALTH
        streak = np.where(is_low, streak + 1, 0)
        new_coll = (streak >= COLLAPSE_STREAK) & (~collapsed)
        if new_coll.any():
            collapsed |= new_coll
            abs_m = np.abs(m)
            collapse_type[new_coll & (abs_m > RIGIDITY_M)] = 1
            collapse_type[new_coll & (abs_m < FRAGMENTATION_M)] = 2
            collapse_type[new_coll & (collapse_type == 0)] = 3

    n_coll = int(collapsed.sum())
    p = n_coll / N_SEEDS
    se = float(np.sqrt(p * (1 - p) / N_SEEDS))
    rigidity = float((collapse_type == 1).sum()) / max(n_coll, 1)
    fragmentation = float((collapse_type == 2).sum()) / max(n_coll, 1)
    return CellResult(variant, J, mult, p, rigidity, fragmentation, se)


def main():
    results: list[CellResult] = []
    for variant in ("baseline", "additive", "log_h", "const_h"):
        h_func = H_FUNCS[variant]
        for j_idx, J in enumerate(J_VALUES):
            for m_idx, mult in enumerate(MULT_VALUES):
                seed = (hash(variant) & 0xFFFF) + j_idx * 1_000_003 + m_idx * 1009
                cell = run_cell(variant, J, mult, h_func, seed)
                results.append(cell)
                print(f"{variant:10s} J={J} mu={mult:3d}  P(coll)={cell.p_collapse:.2f} "
                      f"rigid={cell.rigidity_share:.2f}  SE={cell.se:.3f}")

    # CSVs per variant
    for variant in ("baseline", "additive", "log_h", "const_h"):
        rows = [r for r in results if r.variant == variant]
        out_csv = OUT / f"{variant}_summary.csv"
        with out_csv.open("w", encoding="utf-8") as f:
            f.write("J,mu,p_collapse,rigidity_share,fragmentation_share,SE\n")
            for r in rows:
                f.write(f"{r.J},{r.mult},{r.p_collapse:.4f},"
                        f"{r.rigidity_share:.4f},{r.fragmentation_share:.4f},"
                        f"{r.se:.4f}\n")
        print(f"Wrote {out_csv}")

    # Comparison summary at the headline cell (J=5, mu=100)
    headline = {r.variant: r for r in results if r.J == 5.0 and r.mult == 100}
    band_avg = {}
    for variant in ("baseline", "additive", "log_h", "const_h"):
        cells = [r for r in results if r.variant == variant
                 and r.J >= 4.0 and r.mult == 100]
        band_avg[variant] = float(np.mean([c.p_collapse for c in cells]))

    md = OUT / "comparison_summary.md"
    with md.open("w", encoding="utf-8") as f:
        f.write("# Sensitivity sweep — variant comparison\n\n")
        f.write("Headline cell (J=5, mu=100):\n\n")
        f.write("| variant | P(collapse) | rigidity share | SE |\n")
        f.write("|---|---|---|---|\n")
        for v, r in headline.items():
            f.write(f"| {v} | {r.p_collapse:.2f} | {r.rigidity_share:.2f} | {r.se:.3f} |\n")
        f.write("\nHigh-J band average (J in {4, 4.5, 5}, mu=100):\n\n")
        f.write("| variant | band-avg P(collapse) |\n")
        f.write("|---|---|\n")
        for v, p in band_avg.items():
            f.write(f"| {v} | {p:.2f} |\n")
        f.write("\nReports verify whether the qualitative residual-collapse finding\n")
        f.write("survives changes to the noise prescription and the h(W) functional form.\n")

    print(f"\nWrote {md}")
    print("\nHeadline (J=5, mu=100):")
    for v, r in headline.items():
        print(f"  {v:10s}  P(coll) = {r.p_collapse:.2f}  rigidity = {r.rigidity_share:.2f}")
    print("\nHigh-J band averages:")
    for v, p in band_avg.items():
        print(f"  {v:10s}  P(coll) = {p:.2f}")


if __name__ == "__main__":
    main()
