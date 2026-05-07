"""
Minimal Model — Curie-Weiss SDE coupled to economic feedback.

Two coupled equations only. No Hawkes, no Lotka-Volterra, no mortality, no
demographics, no ecology. Implements the spec in MINIMAL_MODEL_TASK.md.

dm/dt = -m + tanh((J*m + h)/T)        + xi*sqrt(1-m^2)*dW   (Curie-Weiss SDE)
d(wealth)/dt = mult * employment(m) - consumption(wealth)    (economic feedback)

with h = (wealth/500) * 2.0 — the passive stabilizer field.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

# ---------------------------------------------------------------------------
# Sweep configuration (matches MINIMAL_MODEL_TASK.md)
# ---------------------------------------------------------------------------

J_VALUES = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)
MULT_VALUES = (20, 30, 40, 50, 60, 70, 80, 90, 100)
N_SEEDS = 100
N_STEPS = 20_000   # dt = 0.05 -> 1,000 physical time units (matches dt=0.1, N=10,000)

# Fixed model constants
T_TEMP = 2.0
XI = 0.5
DT = 0.05          # dt=0.1 was unconverged by ~10 pp at high J; dt=0.05 saturates
INITIAL_M = 0.0
INITIAL_WEALTH = 100.0

# Collapse thresholds (per task spec; scaled to preserve 10 physical time units)
COLLAPSE_WEALTH = 10.0          # wealth < this for...
COLLAPSE_STREAK = 200           # ...this many consecutive steps  (= 10 physical units at dt=0.05)
RIGIDITY_M = 0.9                # |m| > this at collapse → rigidity
FRAGMENTATION_M = 0.3           # |m| < this at collapse → fragmentation

EPS = 1e-6  # epsilon-clamp to keep sqrt(1-m^2) safe


@dataclass
class CellResult:
    J: float
    mult: float
    collapsed: np.ndarray       # bool[n_seeds]
    collapse_step: np.ndarray   # int32[n_seeds]; -1 if never
    collapse_type: np.ndarray   # int8[n_seeds]; 0=none, 1=rigidity, 2=fragmentation, 3=mixed
    final_m: np.ndarray         # float64[n_seeds]
    final_wealth: np.ndarray    # float64[n_seeds]


def run_cell(
    J: float,
    mult: float,
    n_seeds: int = N_SEEDS,
    n_steps: int = N_STEPS,
    T: float = T_TEMP,
    xi: float = XI,
    dt: float = DT,
    seed: int = 0,
    record_trajectories: bool = False,
    record_every: int = 100,
) -> CellResult | tuple[CellResult, np.ndarray, np.ndarray]:
    """Vectorize n_seeds independent trajectories under one (J, mult) cell."""
    rng = np.random.default_rng(seed)

    m = np.full(n_seeds, INITIAL_M, dtype=np.float64)
    wealth = np.full(n_seeds, INITIAL_WEALTH, dtype=np.float64)

    collapsed = np.zeros(n_seeds, dtype=bool)
    collapse_step = np.full(n_seeds, -1, dtype=np.int32)
    collapse_type = np.zeros(n_seeds, dtype=np.int8)
    wealth_low_streak = np.zeros(n_seeds, dtype=np.int32)

    sqrt_dt = np.sqrt(dt)

    if record_trajectories:
        n_snaps = n_steps // record_every + 1
        m_traj = np.empty((n_snaps, n_seeds), dtype=np.float64)
        w_traj = np.empty((n_snaps, n_seeds), dtype=np.float64)
        m_traj[0] = m
        w_traj[0] = wealth
        snap_idx = 1

    for t in range(n_steps):
        # --- Curie-Weiss SDE ---
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = (wealth / 500.0) * 2.0
        f = -m_c + np.tanh((J * m_c + h) / T)
        g = xi * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f * dt + g * dW, -1 + EPS, 1 - EPS)

        # --- Economic feedback ---
        employment = 0.5 + 0.5 * m
        income = mult * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * dt)

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

        if record_trajectories and (t + 1) % record_every == 0:
            m_traj[snap_idx] = m
            w_traj[snap_idx] = wealth
            snap_idx += 1

    result = CellResult(
        J=J,
        mult=mult,
        collapsed=collapsed,
        collapse_step=collapse_step,
        collapse_type=collapse_type,
        final_m=m,
        final_wealth=wealth,
    )

    if record_trajectories:
        return result, m_traj[:snap_idx], w_traj[:snap_idx]
    return result


def sweep(
    j_values: Iterable[float] = J_VALUES,
    mult_values: Iterable[int] = MULT_VALUES,
    n_seeds: int = N_SEEDS,
    n_steps: int = N_STEPS,
    seed_base: int = 0xC0DE,
) -> list[CellResult]:
    """Run the full (J, mult) sweep. Each cell gets a deterministic unique base seed."""
    cells: list[CellResult] = []
    j_list = list(j_values)
    mult_list = list(mult_values)
    for j_idx, J in enumerate(j_list):
        for m_idx, mult in enumerate(mult_list):
            seed = seed_base + j_idx * 1_000_003 + m_idx * 1009
            cells.append(run_cell(J, mult, n_seeds=n_seeds, n_steps=n_steps, seed=seed))
    return cells
