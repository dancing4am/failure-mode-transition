"""Critical-ramp-speed sweep on the minimal model.

Question: between the gradual (15k-step) and sudden (~100-step) extremes
already characterised in dynamic_j.py, at what ramp duration does collapse
risk cross 50% as a function of economic margin (`mult`)?

Holds total simulation horizon fixed at N_STEPS = 20,000 so that even the
slowest ramp (15k) has 5,000 post-ramp settling steps for collapse to
express. Each ramp is linear from J_START to J_END over `ramp_steps`,
then J holds at J_END for the remainder.

Same model parameters as dynamic_j.py (Curie-Weiss SDE coupled to
economic feedback, h = (wealth/500)*2.0).

Outputs (under results/dynamic_j_ramp_speed/):
    raw_results.csv            — one row per (ramp_steps, mult, seed)
    summary.csv                — one row per (ramp_steps, mult)
    threshold_by_mult.csv      — derived: ramp_steps_at_50pct vs mult
    collapse_vs_ramp_speed.png — Figure A: P(collapse) vs ramp_steps, by mult
    collapse_types_vs_ramp.png — Figure B: collapse-type composition
    threshold_vs_mult.png      — Figure C: 50% threshold vs mult
    VERDICT.md                 — headline numbers + governance framing
"""

from __future__ import annotations

from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Project-relative output dir (works on Windows or WSL)
OUT = Path(__file__).resolve().parents[1] / "results" / "dynamic_j_ramp_speed"
OUT.mkdir(parents=True, exist_ok=True)

# --- Sweep grid ----------------------------------------------------------
# Ramp durations doubled to preserve physical time at dt=0.05 (was {500..15000} at dt=0.1).
RAMP_STEPS_VALUES = (1_000, 2_000, 4_000, 10_000, 20_000, 30_000)
MULT_VALUES = (40, 50, 60, 65, 70, 75, 80, 90, 100)
N_SEEDS = 100
N_STEPS = 40_000  # 2,000 physical time units; slowest ramp (1,500) + 500 settle

# --- Model constants (must match dynamic_j.py) ---------------------------
J_START, J_END = 0.5, 5.0
T = 2.0
XI = 0.5
DT = 0.05
EPS = 1e-6
INITIAL_M = 0.0
INITIAL_WEALTH = 100.0
COLLAPSE_WEALTH = 10.0
COLLAPSE_STREAK = 200       # 10 physical time units at dt=0.05
RIGIDITY_M = 0.9
FRAGMENTATION_M = 0.3


def schedule_linear_hold(step: int, ramp_steps: int) -> float:
    """Ramp linearly J_START -> J_END over ramp_steps, then hold at J_END."""
    if step >= ramp_steps:
        return J_END
    return J_START + (J_END - J_START) * step / ramp_steps


def run_cell(mult: float, ramp_steps: int, n_seeds: int, n_steps: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M, dtype=np.float64)
    wealth = np.full(n_seeds, INITIAL_WEALTH, dtype=np.float64)
    collapsed = np.zeros(n_seeds, dtype=bool)
    collapse_step = np.full(n_seeds, -1, dtype=np.int32)
    collapse_J = np.full(n_seeds, np.nan, dtype=np.float64)
    collapse_type = np.zeros(n_seeds, dtype=np.int8)
    wealth_low_streak = np.zeros(n_seeds, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)

    for t in range(n_steps):
        J_t = schedule_linear_hold(t, ramp_steps)

        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = (wealth / 500.0) * 2.0
        f_drift = -m_c + np.tanh((J_t * m_c + h) / T)
        g = XI * np.sqrt(1 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f_drift * DT + g * dW, -1 + EPS, 1 - EPS)

        employment = 0.5 + 0.5 * m
        income = mult * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        is_low = wealth < COLLAPSE_WEALTH
        wealth_low_streak = np.where(is_low, wealth_low_streak + 1, 0)
        new_collapse = (wealth_low_streak >= COLLAPSE_STREAK) & (~collapsed)
        if new_collapse.any():
            collapsed |= new_collapse
            collapse_step[new_collapse] = t
            collapse_J[new_collapse] = J_t
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
        "collapse_J": collapse_J,
        "collapse_type": collapse_type,
        "final_m": m,
        "final_wealth": wealth,
    }


def interp_threshold(ramps: np.ndarray, p_coll: np.ndarray, level: float = 0.5) -> float:
    """Linear-interpolate (in log10(ramp_steps), p) the ramp duration at which
    p_collapse crosses `level`. Returns nan if curve never crosses."""
    order = np.argsort(ramps)
    r = ramps[order].astype(float)
    p = p_coll[order].astype(float)
    if not (p.min() <= level <= p.max()):
        return float("nan")
    # find the first crossing
    log_r = np.log10(r)
    for i in range(len(p) - 1):
        a, b = p[i], p[i + 1]
        if (a - level) * (b - level) <= 0 and a != b:
            frac = (level - a) / (b - a)
            return float(10 ** (log_r[i] + frac * (log_r[i + 1] - log_r[i])))
    return float("nan")


def main():
    rows = []
    t0 = time.time()
    n_cells = len(RAMP_STEPS_VALUES) * len(MULT_VALUES)
    cell_idx = 0
    for r_idx, ramp_steps in enumerate(RAMP_STEPS_VALUES):
        for m_idx, mult in enumerate(MULT_VALUES):
            cell_idx += 1
            seed = 7919 + r_idx * 1009 + m_idx * 31
            res = run_cell(mult, ramp_steps, N_SEEDS, N_STEPS, seed)
            for s in range(N_SEEDS):
                rows.append({
                    "ramp_steps": ramp_steps,
                    "mult": mult,
                    "seed_idx": s,
                    "collapsed": int(res["collapsed"][s]),
                    "collapse_step": int(res["collapse_step"][s]),
                    "collapse_J": float(res["collapse_J"][s]) if res["collapsed"][s] else np.nan,
                    "collapse_type": int(res["collapse_type"][s]),
                    "final_m": float(res["final_m"][s]),
                    "final_wealth": float(res["final_wealth"][s]),
                })
            print(f"[{cell_idx:>2}/{n_cells}] ramp={ramp_steps:>5} mult={mult:>3} "
                  f"-> p_collapse={res['collapsed'].mean():.2f}")
    elapsed = time.time() - t0
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "raw_results.csv", index=False)
    print(f"\nSweep done in {elapsed:.1f}s — {len(df):,} runs")

    # --- Summary table ------------------------------------------------------
    summary_rows = []
    for ramp_steps in RAMP_STEPS_VALUES:
        for mult in MULT_VALUES:
            cell = df[(df["ramp_steps"] == ramp_steps) & (df["mult"] == mult)]
            mask = cell["collapsed"] == 1
            n = len(cell)
            n_coll = int(mask.sum())
            ctypes = cell.loc[mask, "collapse_type"]
            collapse_Js = cell.loc[mask, "collapse_J"]
            summary_rows.append({
                "ramp_steps": ramp_steps,
                "mult": mult,
                "n": n,
                "n_collapsed": n_coll,
                "p_collapsed": n_coll / n,
                "median_J_at_collapse": float(collapse_Js.median()) if n_coll else np.nan,
                "n_rigidity": int((ctypes == 1).sum()),
                "n_fragmentation": int((ctypes == 2).sum()),
                "n_mixed": int((ctypes == 3).sum()),
                "p_rigidity_of_collapsed": (ctypes == 1).sum() / n_coll if n_coll else np.nan,
                "p_fragmentation_of_collapsed": (ctypes == 2).sum() / n_coll if n_coll else np.nan,
            })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT / "summary.csv", index=False)

    # --- Derived: 50% threshold ramp_steps per mult ------------------------
    threshold_rows = []
    for mult in MULT_VALUES:
        sub = summary_df[summary_df["mult"] == mult].sort_values("ramp_steps")
        ramps = sub["ramp_steps"].values
        p_coll = sub["p_collapsed"].values
        t50 = interp_threshold(ramps, p_coll, 0.5)
        t25 = interp_threshold(ramps, p_coll, 0.25)
        t10 = interp_threshold(ramps, p_coll, 0.10)
        threshold_rows.append({
            "mult": mult,
            "p_at_500": float(sub.iloc[0]["p_collapsed"]),
            "p_at_15000": float(sub.iloc[-1]["p_collapsed"]),
            "ramp_steps_at_50pct": t50,
            "ramp_steps_at_25pct": t25,
            "ramp_steps_at_10pct": t10,
        })
    thresh_df = pd.DataFrame(threshold_rows)
    thresh_df.to_csv(OUT / "threshold_by_mult.csv", index=False)

    # --- Figure A: P(collapse) vs ramp_steps, one curve per mult -----------
    fig, ax = plt.subplots(figsize=(10, 6))
    cmap = plt.get_cmap("viridis")
    for k, mult in enumerate(MULT_VALUES):
        sub = summary_df[summary_df["mult"] == mult].sort_values("ramp_steps")
        ax.plot(sub["ramp_steps"], sub["p_collapsed"],
                marker="o", lw=2, ms=7,
                color=cmap(k / (len(MULT_VALUES) - 1)),
                label=f"mult={mult}")
    ax.axhline(0.5, color="grey", ls="--", lw=0.7, label="50% line")
    ax.set_xscale("log")
    ax.set_xlabel("Ramp duration (steps, log scale)")
    ax.set_ylabel("P(collapse)")
    ax.set_title("Critical ramp speed: P(collapse) vs deployment rate\n"
                 "slower ramps protect; the threshold depends on economic margin")
    ax.set_ylim(-0.02, 1.05)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=9, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "collapse_vs_ramp_speed.png", dpi=300)
    plt.close(fig)

    # --- Figure B: collapse-type composition vs ramp_steps -----------------
    fig, axes = plt.subplots(1, len(MULT_VALUES), figsize=(2.4 * len(MULT_VALUES), 4),
                              sharey=True)
    width = 0.7
    for ax, mult in zip(axes, MULT_VALUES):
        sub = summary_df[summary_df["mult"] == mult].sort_values("ramp_steps")
        x = np.arange(len(sub))
        rig = sub["n_rigidity"].values
        frag = sub["n_fragmentation"].values
        mix = sub["n_mixed"].values
        total = (rig + frag + mix).astype(float)
        total_safe = np.where(total == 0, 1, total)
        rig_p = rig / total_safe
        frag_p = frag / total_safe
        mix_p = mix / total_safe
        ax.bar(x, frag_p, width=width, color="tab:blue", label="fragmentation")
        ax.bar(x, mix_p, width=width, bottom=frag_p, color="tab:gray", label="mixed")
        ax.bar(x, rig_p, width=width, bottom=frag_p + mix_p, color="tab:red", label="rigidity")
        ax.set_xticks(x)
        ax.set_xticklabels([str(int(r)) for r in sub["ramp_steps"]],
                           rotation=45, ha="right", fontsize=8)
        ax.set_title(f"mult={mult}", fontsize=9)
        ax.set_ylim(0, 1.05)
        for xi_, t_n in zip(x, total):
            ax.text(xi_, 1.02, f"n={int(t_n)}", ha="center", fontsize=6, color="grey")
    axes[0].set_ylabel("collapse-type fraction (of collapsed seeds)")
    axes[-1].legend(fontsize=8, loc="lower right")
    fig.suptitle("Collapse-type composition by ramp duration (per mult)\n"
                 "rigidity dominates short ramps; fragmentation appears at slow ramps + low mult")
    fig.tight_layout()
    fig.savefig(OUT / "collapse_types_vs_ramp.png", dpi=300)
    plt.close(fig)

    # --- Figure C: 50% threshold vs mult -----------------------------------
    fig, ax = plt.subplots(figsize=(9, 6))
    valid_50 = thresh_df.dropna(subset=["ramp_steps_at_50pct"])
    valid_25 = thresh_df.dropna(subset=["ramp_steps_at_25pct"])
    valid_10 = thresh_df.dropna(subset=["ramp_steps_at_10pct"])
    if len(valid_50):
        ax.plot(valid_50["mult"], valid_50["ramp_steps_at_50pct"],
                "o-", lw=2, ms=8, color="tab:red", label="50% collapse threshold")
    if len(valid_25):
        ax.plot(valid_25["mult"], valid_25["ramp_steps_at_25pct"],
                "s--", lw=2, ms=7, color="tab:orange", label="25% collapse threshold")
    if len(valid_10):
        ax.plot(valid_10["mult"], valid_10["ramp_steps_at_10pct"],
                "^:", lw=2, ms=7, color="tab:green", label="10% collapse threshold")
    ax.set_yscale("log")
    ax.set_xlabel("income multiplier (economic margin)")
    ax.set_ylabel("ramp duration at threshold (steps, log scale)")
    ax.set_title("Critical ramp speed vs economic margin\n"
                 "richer economies tolerate faster AI deployment")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "threshold_vs_mult.png", dpi=300)
    plt.close(fig)

    # --- Console digest ----------------------------------------------------
    print("\nP(collapse) by (ramp_steps, mult):")
    pivot_p = summary_df.pivot(index="ramp_steps", columns="mult",
                                values="p_collapsed").round(2)
    print(pivot_p.to_string())
    print("\nRigidity share of collapses by (ramp_steps, mult):")
    pivot_r = summary_df.pivot(index="ramp_steps", columns="mult",
                                values="p_rigidity_of_collapsed").round(2)
    print(pivot_r.to_string())
    print("\n50% collapse threshold (ramp_steps) by mult:")
    print(thresh_df[["mult", "ramp_steps_at_50pct",
                     "ramp_steps_at_25pct", "ramp_steps_at_10pct"]].to_string(index=False))


if __name__ == "__main__":
    main()
