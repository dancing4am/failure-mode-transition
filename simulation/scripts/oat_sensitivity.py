"""
One-at-a-time (OAT) sensitivity sweep over the six fixed minimal-model
parameters. Faster substitute for full Sobol (which would require
~10 hours). Each parameter is varied at ±50% from baseline; headline
cell (J=5, mu=100) reported with 100 seeds. 12 perturbation cells +
1 baseline = 13 cells × 100 seeds = 1,300 runs (~1 minute CPU).

Output: results/sensitivity/oat_sensitivity.csv,
        results/sensitivity/oat_sensitivity_summary.md
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "sensitivity"
OUT.mkdir(parents=True, exist_ok=True)

# -- Baseline parameters --
@dataclass
class Params:
    T: float = 2.0
    xi: float = 0.5
    h_coef: float = 2.0          # in h(W) = h_coef * (W/500)
    cons_baseline: float = 30.0
    cons_slope: float = 10.0     # in consumption(W) = baseline + slope * (W/500)
    W0: float = 100.0

J_HEAD = 5.0
MU_HEAD = 100
N_SEEDS = 100
N_STEPS = 20_000
DT = 0.05
INITIAL_M = 0.0
COLLAPSE_WEALTH = 10.0
COLLAPSE_STREAK = 200
RIGIDITY_M = 0.9
FRAGMENTATION_M = 0.3
EPS = 1e-6


def run_cell(p: Params, J: float, mu: float, n_seeds: int, seed: int):
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M, dtype=np.float64)
    wealth = np.full(n_seeds, p.W0, dtype=np.float64)
    collapsed = np.zeros(n_seeds, dtype=bool)
    rigidity = np.zeros(n_seeds, dtype=bool)
    streak = np.zeros(n_seeds, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)

    for _ in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = p.h_coef * (wealth / 500.0)
        f = -m_c + np.tanh((J * m_c + h) / p.T)
        g = p.xi * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)

        employment = 0.5 + 0.5 * m
        income = mu * employment
        consumption = p.cons_baseline + p.cons_slope * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        is_low = wealth < COLLAPSE_WEALTH
        streak = np.where(is_low, streak + 1, 0)
        new_coll = (streak >= COLLAPSE_STREAK) & (~collapsed)
        if new_coll.any():
            collapsed |= new_coll
            rigidity |= (new_coll & (np.abs(m) > RIGIDITY_M))

    n_coll = int(collapsed.sum())
    p_coll = n_coll / n_seeds
    se = float(np.sqrt(p_coll * (1 - p_coll) / n_seeds))
    rig_share = float(rigidity.sum()) / max(n_coll, 1)
    return p_coll, se, rig_share


def main():
    base = Params()
    rows = []

    # Baseline cell first
    t0 = time.time()
    p_b, se_b, rig_b = run_cell(base, J_HEAD, MU_HEAD, N_SEEDS, seed=42)
    elapsed = time.time() - t0
    print(f"baseline                                  P(coll) = {p_b:.2f}  "
          f"SE = {se_b:.3f}  rigidity = {rig_b:.2f}  [{elapsed:.0f}s]")
    rows.append({
        "perturbation": "baseline", "param": "—", "value": "—",
        "p_collapse": p_b, "SE": se_b, "rigidity_share": rig_b,
    })

    # OAT perturbations: each parameter at -50% and +50%
    perturbations = [
        ("T", base.T * 0.5, base.T * 1.5),
        ("xi", base.xi * 0.5, base.xi * 1.5),
        ("h_coef", base.h_coef * 0.5, base.h_coef * 1.5),
        ("cons_baseline", base.cons_baseline * 0.5, base.cons_baseline * 1.5),
        ("cons_slope", base.cons_slope * 0.5, base.cons_slope * 1.5),
        ("W0", base.W0 * 0.5, base.W0 * 1.5),
    ]
    for name, low, high in perturbations:
        for label, val in (("-50%", low), ("+50%", high)):
            kwargs = {name: val}
            p_perturbed = replace(base, **kwargs)
            t0 = time.time()
            seed = (hash(name + label) & 0xFFFF) + 100003
            p_coll, se, rig = run_cell(p_perturbed, J_HEAD, MU_HEAD,
                                        N_SEEDS, seed=seed)
            elapsed = time.time() - t0
            print(f"{name:14s} {label:5s}  ({val:.4g})  "
                  f"P(coll) = {p_coll:.2f}  SE = {se:.3f}  "
                  f"rigidity = {rig:.2f}  [{elapsed:.0f}s]")
            rows.append({
                "perturbation": label, "param": name, "value": val,
                "p_collapse": p_coll, "SE": se, "rigidity_share": rig,
            })

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "oat_sensitivity.csv", index=False)
    print(f"\nWrote {OUT / 'oat_sensitivity.csv'}")

    # Summary md: how far does each perturbation move the residual?
    md = OUT / "oat_sensitivity_summary.md"
    with md.open("w", encoding="utf-8") as f:
        f.write("# OAT sensitivity — single-cell summary\n\n")
        f.write("Headline cell (*J*=5, *μ*=100), 100 seeds. "
                f"Baseline P(collapse) = {p_b:.2f} (SE = {se_b:.3f}).\n\n")
        f.write("| param | -50% P(coll) | +50% P(coll) | range | "
                "rigidity stable? |\n")
        f.write("|---|---|---|---|---|\n")
        for name in ("T", "xi", "h_coef", "cons_baseline", "cons_slope", "W0"):
            sub = df[df.param == name]
            low = sub[sub.perturbation == "-50%"].iloc[0]
            high = sub[sub.perturbation == "+50%"].iloc[0]
            rng = abs(high.p_collapse - low.p_collapse)
            rig_stable = "yes" if abs(high.rigidity_share - rig_b) < 0.2 \
                         and abs(low.rigidity_share - rig_b) < 0.2 else "shifts"
            f.write(f"| {name} | {low.p_collapse:.2f} | {high.p_collapse:.2f} "
                    f"| {rng:.2f} | {rig_stable} |\n")
        f.write("\nRange = |P(coll, +50%) − P(coll, -50%)|.\n")

    print(f"Wrote {md}")


if __name__ == "__main__":
    main()
