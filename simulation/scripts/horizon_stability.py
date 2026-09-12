"""Horizon (t_max) stability check for the minimal-model collapse probability.

Verifies that the cumulative collapse fraction P(collapse) is stable as the
simulation horizon grows from t = 250 to t = 4000 (dt = 0.05, so up to
80,000 steps), for the three headline cells J in {2.5, 3.5, 5.0} at mult = 100.

Because run_cell draws exactly ONE rng.standard_normal(n_seeds) per step,
re-running the SAME seed with a LARGER n_steps reproduces the shorter run's
trajectory prefix exactly, and the collapse flag is permanent (first passage).
Part 1 exploits this to check exact agreement with the archived
results/minimal_model/raw_results.csv (t = 1000 horizon, 20,000 steps).

Parts:
  1. Exact reproduction check vs archived raw_results.csv (n_seeds = 100,
     standard cell seeds, outcomes restricted to collapse_step < 20,000).
  2. High-power run (n_seeds = 4000, seed = standard cell seed + 7_777_777,
     the horizon-stability offset), cumulative P(collapse) at
     t in {250, 500, 1000, 2000, 4000} with 10,000-rep bootstrap 95% CIs.
  3. Collapse-time distribution: median / p90 / p99 / max collapse time in
     physical units; whether the max saturates in absolute time or tracks
     the horizon; counts of collapses after t = 250 / 500 / 1000 / 2000.
  4. Context: locate the very late collapse near step 18,870 in the archived
     raw_results.csv (which J, mult cell it belongs to).

Outputs (results/horizon_stability/):
  horizon_stability.csv, collapse_time_stats.csv, SUMMARY.md
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Mirror the sibling scripts (run_sweep.py etc.): minimal_model lives in the
# same directory; make the import robust to any working directory.
SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from minimal_model import run_cell, J_VALUES, MULT_VALUES, DT  # noqa: E402

OUT = SCRIPTS_DIR.parents[0] / "results" / "horizon_stability"
OUT.mkdir(parents=True, exist_ok=True)

RAW_CSV = SCRIPTS_DIR.parents[0] / "results" / "minimal_model" / "raw_results.csv"

SEED_BASE = 0xC0DE
HS_SEED_OFFSET = 7_777_777          # horizon-stability offset (documented)

J_TEST = (2.5, 3.5, 5.0)
MULT_TEST = 100
HORIZONS_T = (250, 500, 1000, 2000, 4000)   # physical time
N_STEPS_LONG = 80_000                        # t = 4000 at dt = 0.05
N_SEEDS_REPRO = 100
N_SEEDS_POWER = 4000
N_BOOT = 10_000
BOOT_SEED = 20260819

ARCHIVE_STEPS = 20_000                       # archived horizon (t = 1000)


def cell_seed(J: float, mult: int) -> int:
    j_idx = J_VALUES.index(J)
    m_idx = MULT_VALUES.index(mult)
    return SEED_BASE + j_idx * 1_000_003 + m_idx * 1009


def bootstrap_ci(indicator: np.ndarray, n_boot: int, rng: np.random.Generator,
                 chunk: int = 1000) -> tuple[float, float]:
    """Percentile 95% CI for the mean of a 0/1 path indicator, resampling paths."""
    n = indicator.size
    means = np.empty(n_boot, dtype=np.float64)
    done = 0
    while done < n_boot:
        k = min(chunk, n_boot - done)
        idx = rng.integers(0, n, size=(k, n))
        means[done:done + k] = indicator[idx].mean(axis=1)
        done += k
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def main() -> None:
    lines: list[str] = []

    def log(msg: str = "") -> None:
        print(msg)
        lines.append(msg)

    log("# Horizon stability check — minimal model")
    log("")
    log(f"Cells: J in {J_TEST}, mult = {MULT_TEST}; dt = {DT}; "
        f"long horizon = {N_STEPS_LONG} steps (t = {N_STEPS_LONG * DT:.0f}).")
    log("")

    raw = pd.read_csv(RAW_CSV)

    # ------------------------------------------------------------------
    # Part 1 — exact reproduction check vs archived raw_results.csv
    # ------------------------------------------------------------------
    log("## Part 1 — exact reproduction check (n_seeds = 100, standard seeds)")
    log("")
    log("Same cell seed, n_steps = 80,000; outcomes restricted to "
        "collapse_step < 20,000 must exactly match the archive.")
    log("")

    all_exact = True
    for J in J_TEST:
        seed = cell_seed(J, MULT_TEST)
        res = run_cell(J, MULT_TEST, n_seeds=N_SEEDS_REPRO,
                       n_steps=N_STEPS_LONG, seed=seed)
        within = (res.collapse_step >= 0) & (res.collapse_step < ARCHIVE_STEPS)
        new_flag = within.astype(int)
        new_step = np.where(within, res.collapse_step, -1)

        arc = raw[(raw["J"] == J) & (raw["mult"] == MULT_TEST)].sort_values("seed_idx")
        assert len(arc) == N_SEEDS_REPRO, f"archive cell J={J} has {len(arc)} rows"
        arc_flag = arc["collapsed"].to_numpy()
        arc_step = arc["collapse_step"].to_numpy()

        flag_match = bool(np.array_equal(new_flag, arc_flag))
        step_match = bool(np.array_equal(new_step, arc_step))
        exact = flag_match and step_match
        all_exact &= exact
        log(f"- J = {J}, mult = {MULT_TEST} (seed = {seed}): "
            f"archived n_collapsed = {int(arc_flag.sum())}, "
            f"reproduced (step < 20,000) = {int(new_flag.sum())}; "
            f"per-seed collapsed flags match: {flag_match}; "
            f"per-seed collapse_step match: {step_match} -> "
            f"{'EXACT MATCH' if exact else 'MISMATCH'}")
    log("")
    log(f"**Part 1 verdict: {'EXACT per-seed reproduction of the archive' if all_exact else 'MISMATCH — prefix reproduction FAILED'}.**")
    log("")

    # ------------------------------------------------------------------
    # Part 2 — high-power run + cumulative P(collapse) per horizon
    # ------------------------------------------------------------------
    log("## Part 2 — high-power horizon scan "
        f"(n_seeds = {N_SEEDS_POWER}, seed = standard + {HS_SEED_OFFSET})")
    log("")

    boot_rng = np.random.default_rng(BOOT_SEED)
    horizon_rows = []
    stats_rows = []
    power_results = {}

    for J in J_TEST:
        seed = cell_seed(J, MULT_TEST) + HS_SEED_OFFSET
        t0 = time.time()
        res = run_cell(J, MULT_TEST, n_seeds=N_SEEDS_POWER,
                       n_steps=N_STEPS_LONG, seed=seed)
        print(f"  [run] J={J} n={N_SEEDS_POWER} x {N_STEPS_LONG} steps "
              f"in {time.time() - t0:.1f}s")
        power_results[J] = res
        cs = res.collapse_step

        log(f"### J = {J}, mult = {MULT_TEST} (seed = {seed})")
        log("")
        log("| horizon t | n_collapsed | P(collapse) | 95% CI |")
        log("|---|---|---|---|")
        for t_h in HORIZONS_T:
            thr = int(round(t_h / DT))
            ind = ((cs >= 0) & (cs < thr)).astype(np.float64)
            p = float(ind.mean())
            lo, hi = bootstrap_ci(ind, N_BOOT, boot_rng)
            horizon_rows.append({
                "J": J, "mult": MULT_TEST, "n_seeds": N_SEEDS_POWER,
                "horizon_t": t_h, "P_collapse": p,
                "ci_lo": lo, "ci_hi": hi, "n_collapsed": int(ind.sum()),
            })
            log(f"| {t_h} | {int(ind.sum())} | {p:.4f} | [{lo:.4f}, {hi:.4f}] |")
        log("")

    hs_df = pd.DataFrame(horizon_rows)
    hs_df.to_csv(OUT / "horizon_stability.csv", index=False)

    # Drift verdict: t=2000 / t=4000 values vs the t=1000 bootstrap CI
    log("### Drift verdict")
    log("")
    any_drift = False
    for J in J_TEST:
        sub = hs_df[hs_df["J"] == J].set_index("horizon_t")
        p1000 = sub.loc[1000, "P_collapse"]
        lo1000, hi1000 = sub.loc[1000, "ci_lo"], sub.loc[1000, "ci_hi"]
        drift = False
        for t_h in (2000, 4000):
            p = sub.loc[t_h, "P_collapse"]
            if p < lo1000 or p > hi1000:
                drift = True
        any_drift |= drift
        log(f"- J = {J}: P@1000 = {p1000:.4f} (CI [{lo1000:.4f}, {hi1000:.4f}]); "
            f"P@2000 = {sub.loc[2000, 'P_collapse']:.4f}; "
            f"P@4000 = {sub.loc[4000, 'P_collapse']:.4f} -> "
            f"{'DRIFT' if drift else 'no drift'}")
    log("")
    if any_drift:
        log("**VERDICT: DRIFT DETECTED — P(collapse) at t = 2000/4000 falls "
            "outside the t = 1000 bootstrap CI. STOP-level finding.**")
    else:
        log("**VERDICT: NO DRIFT — P(collapse) at t = 2000 and t = 4000 lies "
            "inside the t = 1000 bootstrap CI for every cell.**")
    log("")

    # ------------------------------------------------------------------
    # Part 3 — collapse-time distribution
    # ------------------------------------------------------------------
    log("## Part 3 — collapse-time distribution (n = 4000, full t = 4000 run)")
    log("")
    log("| J | n_collapsed | median t | p90 t | p99 t | max t | max/250 | max/500 | max/1000 | max/2000 | max/4000 | after 250 | after 500 | after 1000 | after 2000 |")
    log("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for J in J_TEST:
        cs = power_results[J].collapse_step
        ct = cs[cs >= 0] * DT   # physical collapse times
        med = float(np.median(ct))
        p90 = float(np.percentile(ct, 90))
        p99 = float(np.percentile(ct, 99))
        mx = float(ct.max())
        after = {t_h: int((ct >= t_h).sum()) for t_h in (250, 500, 1000, 2000)}
        stats_rows.append({
            "J": J, "mult": MULT_TEST, "n_seeds": N_SEEDS_POWER,
            "n_collapsed": int(ct.size),
            "median_t": med, "p90_t": p90, "p99_t": p99, "max_t": mx,
            "max_frac_of_250": mx / 250, "max_frac_of_500": mx / 500,
            "max_frac_of_1000": mx / 1000, "max_frac_of_2000": mx / 2000,
            "max_frac_of_4000": mx / 4000,
            "n_after_250": after[250], "n_after_500": after[500],
            "n_after_1000": after[1000], "n_after_2000": after[2000],
        })
        log(f"| {J} | {ct.size} | {med:.1f} | {p90:.1f} | {p99:.1f} | {mx:.1f} "
            f"| {mx / 250:.3f} | {mx / 500:.3f} | {mx / 1000:.3f} "
            f"| {mx / 2000:.3f} | {mx / 4000:.3f} "
            f"| {after[250]} | {after[500]} | {after[1000]} | {after[2000]} |")
    log("")

    st_df = pd.DataFrame(stats_rows)
    st_df.to_csv(OUT / "collapse_time_stats.csv", index=False)

    max_overall = st_df["max_t"].max()
    if (st_df["n_after_2000"] == 0).all():
        sat = (f"Max collapse time saturates in ABSOLUTE units: the latest "
               f"collapse across all three cells occurs at t = {max_overall:.1f}, "
               f"and no collapses occur after t = 2000 even though the horizon "
               f"extends to t = 4000. The max does NOT track the horizon.")
    else:
        sat = (f"Max collapse time may track the horizon: collapses still occur "
               f"after t = 2000 (latest at t = {max_overall:.1f} of 4000).")
    log(sat)
    log("")

    # ------------------------------------------------------------------
    # Part 4 — locate the archived very-late collapse near step 18,870
    # ------------------------------------------------------------------
    log("## Part 4 — context: where do very late collapses live in the archive?")
    log("")
    late = raw[raw["collapse_step"] >= 18_000].sort_values(
        "collapse_step", ascending=False)
    log(f"Archived rows with collapse_step >= 18,000 "
        f"(of {int(raw['collapsed'].sum())} collapses total):")
    log("")
    if len(late):
        log("| J | mult | seed_idx | collapse_step | t (physical) |")
        log("|---|---|---|---|---|")
        for _, r in late.iterrows():
            log(f"| {r['J']} | {int(r['mult'])} | {int(r['seed_idx'])} "
                f"| {int(r['collapse_step'])} | {r['collapse_step'] * DT:.1f} |")
    else:
        log("(none)")
    log("")
    arc_max = raw.loc[raw["collapse_step"].idxmax()]
    log(f"Latest archived collapse overall: step {int(arc_max['collapse_step'])} "
        f"(t = {arc_max['collapse_step'] * DT:.1f}, "
        f"{arc_max['collapse_step'] / ARCHIVE_STEPS:.1%} of the t = 1000 horizon) "
        f"in cell J = {arc_max['J']}, mult = {int(arc_max['mult'])}.")
    mu100_max = raw[raw["mult"] == MULT_TEST]["collapse_step"].max()
    log(f"Latest archived collapse restricted to mult = 100 cells: "
        f"step {int(mu100_max)} (t = {mu100_max * DT:.1f}).")
    log("")

    (OUT / "SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nWrote {OUT / 'horizon_stability.csv'}")
    print(f"Wrote {OUT / 'collapse_time_stats.csv'}")
    print(f"Wrote {OUT / 'SUMMARY.md'}")


if __name__ == "__main__":
    main()
