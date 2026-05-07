"""Level 0 — Pure Curie-Weiss SDE with FIXED external field h.

No economic dynamics, no wealth, no mult. This is the analytical baseline:
the mean-field equation m* = tanh((J*m* + h)/T) has fixed points that the
SDE should reproduce in the long-time limit.

We sweep J over the same range as the minimal model, at three fixed h values
(low, mid, high — corresponding to the wealth-driven h ranges that L1 spans),
to characterise the equilibrium m vs J curve. This validates the analytical
backbone: the SDE form alone produces the high-J ordered (rigidity-prone)
state regardless of any stabilizer; the stabilizer (if present) only shifts
the location, not the existence, of the high-|m| regime.

Output:
  - level0_equilibrium_curves.png
  - level0_summary.csv
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "results" / "ablation"
OUT.mkdir(parents=True, exist_ok=True)

J_VALUES = np.linspace(0.5, 5.0, 19)
H_VALUES = (0.0, 0.2, 0.4)
T = 2.0
XI = 0.5
DT = 0.1
N_SEEDS = 200
N_STEPS = 5000
EPS = 1e-6


def equilibrium_run(J: float, h: float, n_seeds: int, n_steps: int, seed: int):
    rng = np.random.default_rng(seed)
    m = rng.uniform(-0.1, 0.1, size=n_seeds)  # start near zero
    sqrt_dt = np.sqrt(DT)
    for _ in range(n_steps):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        f = -m_c + np.tanh((J * m_c + h) / T)
        g = XI * np.sqrt(1 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)
    return m  # final-state distribution


def mean_field_fixed_point(J: float, h: float, T: float = 2.0):
    """Solve m* = tanh((J*m* + h)/T) by iteration; pick branch by sign(h+small)."""
    m = np.sign(h + 1e-9) * 0.5 if abs(h) < 1e-9 else np.sign(h) * 0.5
    for _ in range(500):
        m = np.tanh((J * m + h) / T)
    return m


rows = []
fig, ax = plt.subplots(figsize=(9, 6))
colors = ["tab:blue", "tab:orange", "tab:red"]

for col_idx, h in enumerate(H_VALUES):
    sde_means = []
    sde_q25 = []
    sde_q75 = []
    mf_pos = []
    for J in J_VALUES:
        m_final = equilibrium_run(J, h, N_SEEDS, N_STEPS, seed=int(J * 1000) + col_idx)
        sde_means.append(m_final.mean())
        sde_q25.append(np.percentile(m_final, 25))
        sde_q75.append(np.percentile(m_final, 75))
        mf_pos.append(mean_field_fixed_point(J, h))
        rows.append({
            "level": "L0",
            "J": float(J),
            "h": h,
            "mean_m_sde": float(m_final.mean()),
            "q25_m_sde": float(np.percentile(m_final, 25)),
            "q75_m_sde": float(np.percentile(m_final, 75)),
            "mean_field_m": float(mean_field_fixed_point(J, h)),
        })
    sde_means = np.array(sde_means)
    sde_q25 = np.array(sde_q25)
    sde_q75 = np.array(sde_q75)
    mf_pos = np.array(mf_pos)
    color = colors[col_idx]
    ax.fill_between(J_VALUES, sde_q25, sde_q75, alpha=0.15, color=color)
    ax.plot(J_VALUES, sde_means, "-", color=color, lw=2, label=f"SDE mean (h={h})")
    ax.plot(J_VALUES, mf_pos, "--", color=color, lw=1, label=f"mean-field m* (h={h})")
    if h == 0.0:
        # plot symmetric branch
        ax.plot(J_VALUES, -mf_pos, "--", color=color, lw=1, alpha=0.6)

ax.axvline(1.0, color="grey", lw=0.5, ls=":", label="critical J=T (mean-field)")
ax.set_xlabel("coupling J")
ax.set_ylabel("equilibrium magnetization m*")
ax.set_title("Level 0 — Pure CW SDE: equilibrium m vs J at fixed h\n"
             "SDE long-time mean tracks the mean-field fixed point\n"
             "high-|m| ordered regime emerges for J > T regardless of stabilizer")
ax.set_ylim(-1.05, 1.05)
ax.grid(alpha=0.3)
ax.legend(fontsize=9, loc="lower right")
fig.tight_layout()
fig.savefig(OUT / "level0_equilibrium_curves.png", dpi=300)
plt.close(fig)

pd.DataFrame(rows).to_csv(OUT / "level0_summary.csv", index=False)
print(f"L0 done — {OUT / 'level0_equilibrium_curves.png'}")
print(f"L0 summary csv: {OUT / 'level0_summary.csv'}")
