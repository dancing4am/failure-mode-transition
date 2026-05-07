"""Dynamic J(t) — singularity-trajectory experiment on the minimal model.

J ramps over time to model AI-driven coupling tightening. The question is:
which economies (parameterised by income multiplier `mult`) survive to the
highest J before collapsing, and through which failure mode?

Schedules:
  linear:        J(t) = J0 + (J1 - J0) * t / T
  exponential:   J(t) = J0 * (J1/J0)^(t/T)         (S-curve adoption proxy)

Sweep:
  mult ∈ [40, 50, 60, 65, 70, 75, 80, 90, 100]
  100 seeds per mult per schedule
  5000 steps per run

For each run we record:
  - first-collapse step (-1 if never collapsed)
  - J at that step
  - collapse type at that step (rigidity / fragmentation / mixed)

Outputs:
  results/dynamic_j/raw_results.csv
  results/dynamic_j/survival_curves.png
  results/dynamic_j/collapse_types_by_J.png
  results/dynamic_j/critical_j_vs_mult.png
  results/dynamic_j/VERDICT.md
"""

from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "results" / "dynamic_j"
OUT.mkdir(parents=True, exist_ok=True)

MULT_VALUES = (40, 50, 60, 65, 70, 75, 80, 90, 100)
N_SEEDS = 100
# 30,000 steps at dt=0.05 = 1,500 physical time units (matches dt=0.1 × 15,000).
N_STEPS = 30_000
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


def schedule_linear(step: int, n_steps: int) -> float:
    """Gradual: J(t) = J0 + (J1-J0)·t/T. Adiabatic AI deployment."""
    return J_START + (J_END - J_START) * step / n_steps


def schedule_exponential(step: int, n_steps: int) -> float:
    """S-curve adoption: J(t) = J0·(J1/J0)^(t/T). Rapid late."""
    frac = step / n_steps
    return J_START * (J_END / J_START) ** frac


# Sudden = step function: J=J_START until t=t_jump, then J=J_END.
# Models a discontinuous AI capability gain (singularity-style onset).
SUDDEN_JUMP_STEP = 200      # 10 physical time units at dt=0.05 (was 100 at dt=0.1)


def schedule_sudden(step: int, n_steps: int) -> float:
    return J_END if step >= SUDDEN_JUMP_STEP else J_START


SCHEDULES = {
    "linear": schedule_linear,
    "exponential": schedule_exponential,
    "sudden": schedule_sudden,
}


def run_cell_dynamic(
    mult: float, schedule_name: str,
    n_seeds: int = N_SEEDS, n_steps: int = N_STEPS, seed: int = 0,
) -> dict:
    schedule = SCHEDULES[schedule_name]
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
        J_t = schedule(t, n_steps)

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


def main():
    rows = []
    t0 = time.time()
    for sched_name in SCHEDULES:
        for m_idx, mult in enumerate(MULT_VALUES):
            seed = (hash(sched_name) & 0xFFFF) + m_idx * 1009
            res = run_cell_dynamic(mult, sched_name, seed=seed)
            for s in range(N_SEEDS):
                rows.append({
                    "schedule": sched_name,
                    "mult": mult,
                    "seed_idx": s,
                    "collapsed": int(res["collapsed"][s]),
                    "collapse_step": int(res["collapse_step"][s]),
                    "collapse_J": float(res["collapse_J"][s]) if res["collapsed"][s] else np.nan,
                    "collapse_type": int(res["collapse_type"][s]),
                    "final_m": float(res["final_m"][s]),
                    "final_wealth": float(res["final_wealth"][s]),
                })
    elapsed = time.time() - t0
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "raw_results.csv", index=False)
    print(f"Sweep done in {elapsed:.1f}s — {len(df):,} runs")

    # ---- Figure A: survival curves -----------------------------------------------
    # For each schedule and each mult, plot fraction surviving vs J(t).
    # Build a grid of J values and compute frac surviving by counting collapse_J > J_grid.
    J_grid = np.linspace(J_START, J_END, 200)

    fig, axes = plt.subplots(1, len(SCHEDULES), figsize=(7 * len(SCHEDULES), 6), sharey=True)
    cmap = plt.get_cmap("viridis")
    for ax, sched in zip(axes, SCHEDULES):
        sub = df[df["schedule"] == sched]
        for k, mult in enumerate(MULT_VALUES):
            cell = sub[sub["mult"] == mult]
            n = len(cell)
            # frac surviving at J = (# seeds with collapse_J > J or never collapsed) / n
            collapse_Js = cell["collapse_J"].values
            never = ~cell["collapsed"].astype(bool).values
            survival = np.array([
                (np.sum(collapse_Js > J) + never.sum()) / n
                for J in J_grid
            ])
            ax.plot(J_grid, survival, color=cmap(k / (len(MULT_VALUES) - 1)),
                    lw=2, label=f"mult={mult}")
        ax.set_xlabel("J(t)  (coupling at the moment of collapse)")
        ax.set_title(f"Schedule: {sched}")
        ax.set_ylim(-0.02, 1.05)
        ax.grid(alpha=0.3)
        ax.legend(title="income multiplier", fontsize=9, loc="lower left", ncol=2)
    axes[0].set_ylabel("fraction surviving")
    fig.suptitle("Singularity trajectory: survival curves under dynamic J(t)\n"
                 "richer economies survive to higher coupling before collapsing")
    fig.tight_layout()
    fig.savefig(OUT / "survival_curves.png", dpi=300)
    plt.close(fig)
    print("wrote survival_curves.png")

    # ---- Figure B: collapse type by J at collapse --------------------------------
    fig, axes = plt.subplots(1, len(SCHEDULES), figsize=(7 * len(SCHEDULES), 6), sharey=True)
    bins = np.linspace(J_START, J_END, 11)
    bin_centres = 0.5 * (bins[:-1] + bins[1:])
    width = (bins[1] - bins[0]) * 0.85
    for ax, sched in zip(axes, SCHEDULES):
        sub = df[(df["schedule"] == sched) & (df["collapsed"] == 1)]
        rig = np.histogram(sub.loc[sub["collapse_type"] == 1, "collapse_J"], bins=bins)[0]
        frag = np.histogram(sub.loc[sub["collapse_type"] == 2, "collapse_J"], bins=bins)[0]
        mixed = np.histogram(sub.loc[sub["collapse_type"] == 3, "collapse_J"], bins=bins)[0]
        total = rig + frag + mixed
        # Normalise to proportions
        total_safe = np.where(total == 0, 1, total)
        rig_p = rig / total_safe
        frag_p = frag / total_safe
        mixed_p = mixed / total_safe
        # Stacked bars
        ax.bar(bin_centres, frag_p, width=width, label="fragmentation",
               color="tab:blue", alpha=0.85)
        ax.bar(bin_centres, mixed_p, width=width, bottom=frag_p, label="mixed",
               color="tab:gray", alpha=0.7)
        ax.bar(bin_centres, rig_p, width=width, bottom=frag_p + mixed_p,
               label="rigidity", color="tab:red", alpha=0.85)
        # Annotate total count per bin
        for c, t_n in zip(bin_centres, total):
            ax.text(c, 1.02, f"n={int(t_n)}", ha="center", fontsize=7, color="grey")
        ax.set_xlabel("J at collapse")
        ax.set_title(f"Schedule: {sched}")
        ax.set_ylim(0, 1.12)
        ax.legend(loc="lower left", fontsize=9)
    axes[0].set_ylabel("collapse-type fraction")
    fig.suptitle("Collapse type vs J at collapse: fragmentation dominates at low J,\n"
                 "rigidity dominates at high J — the asymmetry, dynamically resolved")
    fig.tight_layout()
    fig.savefig(OUT / "collapse_types_by_J.png", dpi=300)
    plt.close(fig)
    print("wrote collapse_types_by_J.png")

    # ---- Figure C: critical J* vs mult --------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 6))
    markers = {"linear": "o", "exponential": "s", "sudden": "^"}
    for sched in SCHEDULES:
        sub = df[df["schedule"] == sched]
        rows_c = []
        for mult in MULT_VALUES:
            cell = sub[sub["mult"] == mult]
            collapsed_J = cell.loc[cell["collapsed"] == 1, "collapse_J"].dropna().values
            n_collapsed = len(collapsed_J)
            n_total = len(cell)
            if n_collapsed > 0:
                med = np.median(collapsed_J)
                q25 = np.percentile(collapsed_J, 25)
                q75 = np.percentile(collapsed_J, 75)
            else:
                med = q25 = q75 = np.nan
            rows_c.append({
                "schedule": sched, "mult": mult,
                "median_J_at_collapse": med,
                "q25": q25, "q75": q75,
                "frac_collapsed": n_collapsed / n_total,
            })
        cdf = pd.DataFrame(rows_c)
        cdf.to_csv(OUT / f"critical_j_summary_{sched}.csv", index=False)
        valid = cdf.dropna(subset=["median_J_at_collapse"])
        ax.errorbar(
            valid["mult"], valid["median_J_at_collapse"],
            yerr=[valid["median_J_at_collapse"] - valid["q25"],
                  valid["q75"] - valid["median_J_at_collapse"]],
            marker=markers[sched], lw=2, ms=8, capsize=4,
            label=f"{sched} ramp",
        )
    ax.axhline(J_START, color="grey", ls=":", lw=0.7, label=f"J_start={J_START}")
    ax.axhline(J_END, color="grey", ls="--", lw=0.7, label=f"J_end={J_END}")
    ax.set_xlabel("income multiplier (economic margin)")
    ax.set_ylabel("median J at collapse  (IQR shown)")
    ax.set_title("Critical J* vs economic margin\n"
                 "richer economies tolerate more coupling before collapse")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "critical_j_vs_mult.png", dpi=300)
    plt.close(fig)
    print("wrote critical_j_vs_mult.png")

    # ---- Verdict numbers ---------------------------------------------------------
    summary_rows = []
    for sched in SCHEDULES:
        sub = df[df["schedule"] == sched]
        for mult in MULT_VALUES:
            cell = sub[sub["mult"] == mult]
            collapsed_mask = cell["collapsed"] == 1
            n = len(cell)
            n_coll = int(collapsed_mask.sum())
            collapse_Js = cell.loc[collapsed_mask, "collapse_J"]
            ctypes = cell.loc[collapsed_mask, "collapse_type"]
            summary_rows.append({
                "schedule": sched, "mult": mult,
                "n": n, "n_collapsed": n_coll,
                "p_collapsed_total": n_coll / n,
                "median_J_at_collapse": float(collapse_Js.median()) if n_coll else np.nan,
                "n_rigidity": int((ctypes == 1).sum()),
                "n_fragmentation": int((ctypes == 2).sum()),
                "n_mixed": int((ctypes == 3).sum()),
                "p_rigidity_of_collapsed": ((ctypes == 1).sum() / n_coll) if n_coll else np.nan,
                "p_fragmentation_of_collapsed": ((ctypes == 2).sum() / n_coll) if n_coll else np.nan,
            })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT / "summary.csv", index=False)
    print(f"wrote summary.csv  ({len(summary_df)} rows)")

    # Pretty digest
    print("\nMedian J at collapse by (schedule, mult):")
    pivot_med = summary_df.pivot(index="schedule", columns="mult",
                                 values="median_J_at_collapse").round(2)
    print(pivot_med.to_string())
    print("\nFraction collapsed (schedule, mult):")
    pivot_p = summary_df.pivot(index="schedule", columns="mult",
                                values="p_collapsed_total").round(2)
    print(pivot_p.to_string())
    print("\nRigidity share of collapsed by (schedule, mult):")
    pivot_rig = summary_df.pivot(index="schedule", columns="mult",
                                  values="p_rigidity_of_collapsed").round(2)
    print(pivot_rig.to_string())


if __name__ == "__main__":
    main()
