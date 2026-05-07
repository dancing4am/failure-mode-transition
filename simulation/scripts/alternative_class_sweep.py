"""
Alternative-model-class sweep — addresses the Limitations open question
"Whether equivalent scaling emerges in Kuramoto, voter, or compartmental
dynamics is open."

We test three non-Curie-Weiss order-parameter dynamics for the same
qualitative passive-stabilizer asymmetry the main text reports. Same
(J, mu) sweep grid, same wealth feedback, same collapse criteria as the
minimal model — only the order-parameter dynamics differ.

Variants:
  curie_weiss  — baseline:  -m + tanh((Jm + h)/T)              (paper main text)
  voter_like   — -m + sign((Jm + h)/T) · min(1, |(Jm+h)/T|)    (piecewise-linear sign)
  kuramoto_1d  — -sin(m*pi/2) + tanh((Jm + h)/T)               (sine-driven restoring)
  cubic        — -m + tanh((Jm + h)/T) but f(m) = m - m^3/3    (cubic landscape)

The voter_like and kuramoto_1d variants test "is the asymmetry tanh-
specific?". cubic tests "is it Curie-Weiss free-energy specific?". If
the high-J residual persists across all four with comparable rigidity
dominance, the result is generic to bistable mean-field SDEs with a
field bias, not to the Curie-Weiss form alone.

Compute: 4 variants × 12 cells (high-J band: J∈{3.5,4,4.5,5} × μ∈{60,80,100})
× 100 seeds = 4,800 runs. ~3 minutes CPU.

Output:
  results/alternative_class/{curie_weiss,voter_like,kuramoto_1d,cubic}_summary.csv
  results/alternative_class/comparison_summary.md
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "alternative_class"
OUT.mkdir(parents=True, exist_ok=True)

J_VALUES = (3.5, 4.0, 4.5, 5.0)
MULT_VALUES = (60, 80, 100)
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


def f_curie_weiss(m, J, h, T):
    """Standard Curie-Weiss mean-field self-consistency drift."""
    return -m + np.tanh((J * m + h) / T)


def f_voter_like(m, J, h, T):
    """Piecewise-linear sign function — 'voter-style' bistable update.

    Replaces tanh with sign(x)·min(1, |x|), which has the same saturating
    behavior at large argument but is not differentiable at zero.
    Tests whether the asymmetry depends on the smooth tanh form."""
    arg = (J * m + h) / T
    return -m + np.sign(arg) * np.minimum(1.0, np.abs(arg))


def f_kuramoto_1d(m, J, h, T):
    """Sine-driven mean-field — Kuramoto-inspired restoring force.

    Replaces -m + tanh(x) with -sin(m·π/2) + tanh(x). The first term
    has the same +/- bistable structure but as a sine wave; the second
    term retains the field-driven mean-field bias."""
    return -np.sin(m * np.pi / 2) + np.tanh((J * m + h) / T)


def f_cubic(m, J, h, T):
    """Cubic Landau-style restoring force.

    Replaces the tanh saturation with a cubic potential (Landau theory
    of phase transitions): f(m) = m - m^3/3, plus the same field
    contribution. Tests whether the asymmetry depends on Curie-Weiss
    free-energy form versus generic Landau bistability."""
    return -(m - m**3 / 3.0) + np.tanh((J * m + h) / T)


VARIANTS = {
    "curie_weiss": f_curie_weiss,
    "voter_like":  f_voter_like,
    "kuramoto_1d": f_kuramoto_1d,
    "cubic":       f_cubic,
}


@dataclass
class CellResult:
    variant: str
    J: float
    mult: float
    p_collapse: float
    rigidity_share: float
    se: float


def run_cell(variant: str, drift: Callable, J: float, mu: float, seed: int) -> CellResult:
    rng = np.random.default_rng(seed)
    m = np.full(N_SEEDS, INITIAL_M, dtype=np.float64)
    wealth = np.full(N_SEEDS, INITIAL_WEALTH, dtype=np.float64)
    collapsed = np.zeros(N_SEEDS, dtype=bool)
    rigidity = np.zeros(N_SEEDS, dtype=bool)
    streak = np.zeros(N_SEEDS, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)

    for _ in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = 2.0 * (wealth / 500.0)
        f = drift(m_c, J, h, T_TEMP)
        g = XI * np.sqrt(np.maximum(0.0, 1.0 - m_c * m_c))
        dW = sqrt_dt * rng.standard_normal(N_SEEDS)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)

        employment = 0.5 + 0.5 * m
        income = mu * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        is_low = wealth < COLLAPSE_WEALTH
        streak = np.where(is_low, streak + 1, 0)
        new_coll = (streak >= COLLAPSE_STREAK) & (~collapsed)
        if new_coll.any():
            collapsed |= new_coll
            rigidity |= (new_coll & (np.abs(m) > RIGIDITY_M))

    n_coll = int(collapsed.sum())
    p = n_coll / N_SEEDS
    se = float(np.sqrt(p * (1 - p) / N_SEEDS))
    rig_share = float(rigidity.sum()) / max(n_coll, 1)
    return CellResult(variant, J, mu, p, rig_share, se)


def main():
    results: list[CellResult] = []
    for variant, drift in VARIANTS.items():
        for j_idx, J in enumerate(J_VALUES):
            for m_idx, mu in enumerate(MULT_VALUES):
                seed = (hash(variant) & 0xFFFF) + j_idx * 1_000_003 + m_idx * 1009
                cell = run_cell(variant, drift, J, mu, seed)
                results.append(cell)
                print(f"{variant:12s} J={J} mu={mu:3d}  "
                      f"P(coll)={cell.p_collapse:.2f}  "
                      f"rigidity={cell.rigidity_share:.2f}  "
                      f"SE={cell.se:.3f}")

    for variant in VARIANTS:
        rows = [r for r in results if r.variant == variant]
        out_csv = OUT / f"{variant}_summary.csv"
        with out_csv.open("w", encoding="utf-8") as f:
            f.write("J,mu,p_collapse,rigidity_share,SE\n")
            for r in rows:
                f.write(f"{r.J},{r.mult},{r.p_collapse:.4f},"
                        f"{r.rigidity_share:.4f},{r.se:.4f}\n")
        print(f"Wrote {out_csv}")

    # Comparison summary
    headline = {r.variant: r for r in results
                if r.J == 5.0 and r.mult == 100}
    band_avg = {}
    band_rig = {}
    for variant in VARIANTS:
        cells = [r for r in results if r.variant == variant
                 and r.J >= 4.0 and r.mult == 100]
        band_avg[variant] = float(np.mean([c.p_collapse for c in cells]))
        band_rig[variant] = float(np.mean([c.rigidity_share for c in cells]))

    md = OUT / "comparison_summary.md"
    with md.open("w", encoding="utf-8") as f:
        f.write("# Alternative-model-class sweep — variant comparison\n\n")
        f.write("Headline cell (J=5, mu=100):\n\n")
        f.write("| variant | P(collapse) | rigidity share | SE |\n")
        f.write("|---|---|---|---|\n")
        for v, r in headline.items():
            f.write(f"| {v} | {r.p_collapse:.2f} | {r.rigidity_share:.2f} "
                    f"| {r.se:.3f} |\n")
        f.write("\nHigh-J band (J in {4,4.5,5}, mu=100):\n\n")
        f.write("| variant | band-avg P(collapse) | band-avg rigidity share |\n")
        f.write("|---|---|---|\n")
        for v in VARIANTS:
            f.write(f"| {v} | {band_avg[v]:.2f} | {band_rig[v]:.2f} |\n")
        f.write("\nIf all four variants show comparable residual P(collapse)\n")
        f.write("and rigidity dominance, the asymmetry is generic to bistable\n")
        f.write("mean-field SDEs with field bias, not Curie-Weiss-specific.\n")

    print(f"\nWrote {md}")
    print("\nHeadline (J=5, mu=100):")
    for v, r in headline.items():
        print(f"  {v:12s}  P(coll)={r.p_collapse:.2f}  rigidity={r.rigidity_share:.2f}")
    print("\nHigh-J band averages (J >= 4, mu=100):")
    for v in VARIANTS:
        print(f"  {v:12s}  P(coll)={band_avg[v]:.2f}  rigidity={band_rig[v]:.2f}")


if __name__ == "__main__":
    main()
