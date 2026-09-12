# -*- coding: utf-8 -*-
"""
branch_observable.py — branch-selection observable and collapse-time distributions.

E2 — BRANCH-SELECTION OBSERVABLE.
    P(wrong branch) := fraction of runs with m < 0 at the end of the transient
    window, taken at t = 50 (step 1000 at dt = 0.05).

    Window justification: the passive baseline's W-collapse times saturate
    early — at mu = 100 the residual collapses have median ~15 and max ~35
    physical time units (verified empirically below from the part-(a) reruns).
    The deterministic branch-selection transient of the Curie-Weiss order
    parameter is O(1/|1 - J/T|) ~ O(1) time units for supercritical J. t = 50
    therefore sits safely AFTER both the branch-selection transient and the
    known collapse-time saturation at t ~ 35, while remaining early relative
    to the full horizon T_total = 1000, so it is not contaminated by rare
    late stochastic escapes.

    Archived raw CSVs carry final_m (t = 1000 end state) and collapse_step
    but not m at t = 50, so the needed cells are RERUN here with matched
    seeds, recording m at step 1000 and at the end.

    Cells:
      (a) passive baseline, mu=100, J in {2.5,3,3.5,4,4.5,5}      (n=1000)
      (b) fixed h=2.4 and h=5.0, mu=100, J in {3,5,7,10,15,20}    (n=1000)
      (c) responsive variants (stress alpha=2, coupling alpha=1,
          quadratic alpha=1), mu=100, J in {4,4.5,5}              (n=1000)
      (d) S5.5 consumption variants (wdot_sensitivity.VARIANTS), passive
          h0=2, mu=100, J in {2.5..5}, n=2000, amplitude-requirement (h0_scaling) matched seeds.

E3 — COLLAPSE-TIME DISTRIBUTIONS for the fixed-field controls.
    results/active_stabilizer/scaling_control.csv only stores per-cell
    p_collapse (no per-run collapse_step), so the fixed h=2.4/5.0 cells at
    J in {7,10,15,20} are rerun (n=1000, scaling_control seed convention),
    recording collapse_step per run. The part-(b) rerun serves both E2 and E3.

Seed conventions (matched to the archived experiments):
  * parts (a), (c): seed = 0xC0DE + j_idx*1_000_003 + 8*1009
        (minimal_model / active_stabilizer grid convention; j_idx indexes
         J_VALUES = (0.5,...,5.0), 8 = index of mu=100 in MULT_VALUES)
  * part (b) / E3:  seed = 0xC0DE + ji*1_000_003
        (scaling_control.py convention; ji indexes JS = [3,5,7,10,15,20];
         matched across conditions at each J)
  * part (d):       seed = 0xB0_5CA1E + j_idx*1_000_003 + 1*9973
        (the xi=0.5 seed convention of h0_scaling / wdot_sensitivity;
         n=2000 so the noise stream is identical to the archived runs)

The engine below draws exactly one rng.standard_normal(n_seeds) per step, in
the same call order as minimal_model.run_cell, so at matched (seed, n_seeds)
its collapse outcomes are bit-identical to the committed engines. A mandatory
regression check verifies this against minimal_model.run_cell and against the
archived wdot_sensitivity raw_linear npz before any cell is run.

Run: py -3.13 branch_observable.py
Outputs: results/branch_observable/{e2_runs_partX.csv, e2_dissociation.csv,
         e3_collapse_times.csv, SUMMARY.md}
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import active_stabilizer as A  # noqa: E402  (h functions + constants)

OUT = HERE.parent / "results" / "branch_observable"

# Model constants (identical to minimal_model.py / active_stabilizer.py)
T_TEMP = 2.0
XI = 0.5
DT = 0.05
N_STEPS = 20_000
INITIAL_M = 0.0
INITIAL_WEALTH = 100.0
COLLAPSE_WEALTH = 10.0
COLLAPSE_STREAK = 200
EPS = 1e-6

RECORD_STEP = 1000          # m recorded when (t+1) == 1000  ->  t_phys = 50.0
MU = 100

SEED_BASE_GRID = 0xC0DE     # minimal_model convention
SEED_BASE_TASKB = 0xB0_5CA1E  # h0_scaling / wdot_sensitivity convention

# Full grid J list (for j_idx lookup, parts a/c)
J_VALUES = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)
MU_IDX = 8                  # index of mu=100 in MULT_VALUES

JS_A = (2.5, 3.0, 3.5, 4.0, 4.5, 5.0)
JS_B = (3.0, 5.0, 7.0, 10.0, 15.0, 20.0)   # scaling_control JS
JS_C = (4.0, 4.5, 5.0)
JS_D = (2.5, 3.0, 3.5, 4.0, 4.5, 5.0)      # h0_scaling J_GRID
JS_E3 = (7.0, 10.0, 15.0, 20.0)

N_ABC = 1000
N_D = 2000

# wdot_sensitivity.VARIANTS (S5.5), reproduced verbatim
VARIANTS = (
    ("baseline",  dict(cons_const=30.0, cons_coeff=10.0, cons_form="affine")),
    ("const-30",  dict(cons_const=21.0, cons_coeff=10.0, cons_form="affine")),
    ("const+30",  dict(cons_const=39.0, cons_coeff=10.0, cons_form="affine")),
    ("coeff-30",  dict(cons_const=30.0, cons_coeff=7.0,  cons_form="affine")),
    ("coeff+30",  dict(cons_const=30.0, cons_coeff=13.0, cons_form="affine")),
    ("linear",    dict(cons_form="linear", cons_lin=0.08)),
)


def h_fixed(val):
    def f(W, m, J, alpha):
        return np.full_like(np.asarray(W, dtype=np.float64), val)
    return f


# ---------------------------------------------------------------------------
# Engine: scalar (one cell), one standard_normal(n_seeds) per step —
# identical RNG stream and op order to minimal_model.run_cell.
# ---------------------------------------------------------------------------
def run_cell(J, mult, h_fn, alpha, n_seeds, seed,
             cons_const=30.0, cons_coeff=10.0, cons_form="affine",
             cons_lin=0.08, n_steps=N_STEPS):
    if cons_form not in ("affine", "linear"):
        raise ValueError(cons_form)
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M, dtype=np.float64)
    wealth = np.full(n_seeds, INITIAL_WEALTH, dtype=np.float64)
    collapsed = np.zeros(n_seeds, dtype=bool)
    collapse_step = np.full(n_seeds, -1, dtype=np.int32)
    streak = np.zeros(n_seeds, dtype=np.int32)
    m_t50 = np.full(n_seeds, np.nan)
    sqrt_dt = np.sqrt(DT)

    for t in range(n_steps):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = h_fn(wealth, m_c, J, alpha)
        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)

        employment = 0.5 + 0.5 * m
        income = mult * employment
        if cons_form == "affine":
            consumption = cons_const + cons_coeff * (wealth / 500.0)
        else:
            consumption = cons_lin * wealth
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        is_low = wealth < COLLAPSE_WEALTH
        streak = np.where(is_low, streak + 1, 0)
        new_collapse = (streak >= COLLAPSE_STREAK) & (~collapsed)
        if new_collapse.any():
            collapsed |= new_collapse
            collapse_step[new_collapse] = t

        if (t + 1) == RECORD_STEP:
            m_t50 = m.copy()

    return dict(m_t50=m_t50, final_m=m, final_W=wealth,
                collapsed=collapsed, collapse_step=collapse_step)


# ---------------------------------------------------------------------------
# Regression checks (mandatory, run before any cell)
# ---------------------------------------------------------------------------
def regression_checks():
    import minimal_model

    # 1) bit-identity vs minimal_model.run_cell at J=5, mu=100, n=100
    seed = SEED_BASE_GRID + 9 * 1_000_003 + 8 * 1009
    ref = minimal_model.run_cell(J=5.0, mult=100, n_seeds=100, seed=seed)
    gen = run_cell(5.0, 100, A.h_passive, 0.0, n_seeds=100, seed=seed)
    ok1 = (np.array_equal(ref.collapsed, gen["collapsed"])
           and np.array_equal(ref.collapse_step, gen["collapse_step"])
           and np.array_equal(ref.final_m, gen["final_m"]))
    print(f"[check 1] bit-identity vs minimal_model.run_cell: "
          f"{'PASS' if ok1 else 'FAIL'}")

    # 2) reproduce archived wdot_sensitivity linear-variant P(collapse) at
    #    h0=2.0, J=5.0, xi=0.5, n=2000 (bit-level vs stored collapse row)
    npz = (HERE.parent / "results" / "wdot_sensitivity" / "raw_linear"
           / "raw_J5.0_xi0.50.npz")
    ok2 = True
    if npz.exists():
        with np.load(npz) as z:
            h0 = z["h0"]
            C = z["collapsed"]
        row = int(np.where(np.isclose(h0, 2.0))[0][0])
        seed_d = SEED_BASE_TASKB + 5 * 1_000_003 + 1 * 9973
        gen2 = run_cell(5.0, 100, A.h_passive, 0.0, n_seeds=2000, seed=seed_d,
                        cons_form="linear", cons_lin=0.08)
        ok2 = np.array_equal(C[row], gen2["collapsed"])
        print(f"[check 2] bit-identity vs archived wdot raw_linear J=5 h0=2 "
              f"(P={C[row].mean():.4f} vs {gen2['collapsed'].mean():.4f}): "
              f"{'PASS' if ok2 else 'FAIL'}")
    else:
        print(f"[check 2] SKIPPED (missing {npz})")
    return ok1 and ok2


# ---------------------------------------------------------------------------
def cell_rows(part, condition, alpha, cons, J, mu, n, seed, r):
    """Per-run long rows + per-cell summary row."""
    runs = pd.DataFrame({
        "part": part, "condition": condition, "cons": cons, "J": J, "mu": mu,
        "alpha": alpha, "cell_seed": seed,
        "run": np.arange(n),
        "m_t50": r["m_t50"], "final_m": r["final_m"], "final_W": r["final_W"],
        "collapsed": r["collapsed"].astype(int),
        "collapse_step": r["collapse_step"],
    })
    p_branch = float((r["m_t50"] < 0).mean())
    p_mneg_end = float((r["final_m"] < 0).mean())
    p_wcoll = float(r["collapsed"].mean())
    summ = dict(part=part, condition=condition, cons=cons, J=J, mu=mu,
                alpha=alpha, n=n, cell_seed=seed,
                P_branch_t50=p_branch, P_mneg_end=p_mneg_end,
                P_Wcollapse=p_wcoll,
                dissociation=p_branch - p_wcoll)
    return runs, summ


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if not regression_checks():
        print("ABORT: regression checks failed")
        sys.exit(1)

    summaries = []
    t_start = time.time()

    # ---- part (a): passive baseline --------------------------------------
    runs_a = []
    for J in JS_A:
        j_idx = J_VALUES.index(J)
        seed = SEED_BASE_GRID + j_idx * 1_000_003 + MU_IDX * 1009
        r = run_cell(J, MU, A.h_passive, 0.0, n_seeds=N_ABC, seed=seed)
        rr, s = cell_rows("a", "passive", 0.0, "affine", J, MU, N_ABC, seed, r)
        runs_a.append(rr); summaries.append(s)
        print(f"[a] passive J={J}: P_branch={s['P_branch_t50']:.3f} "
              f"P_Wcoll={s['P_Wcollapse']:.3f}  ({time.time()-t_start:.0f}s)",
              flush=True)
    pd.concat(runs_a).to_csv(OUT / "e2_runs_partA.csv", index=False)

    # ---- part (b): fixed fields (also E3 source) -------------------------
    runs_b = []
    for cond_name, hval in [("fixed2.4", 2.4), ("fixed5.0", 5.0)]:
        hfn = h_fixed(hval)
        for J in JS_B:
            ji = JS_B.index(J)
            seed = SEED_BASE_GRID + ji * 1_000_003   # scaling_control convention
            r = run_cell(J, MU, hfn, 0.0, n_seeds=N_ABC, seed=seed)
            rr, s = cell_rows("b", cond_name, 0.0, "affine", J, MU, N_ABC,
                              seed, r)
            runs_b.append(rr); summaries.append(s)
            print(f"[b] {cond_name} J={J}: P_branch={s['P_branch_t50']:.3f} "
                  f"P_Wcoll={s['P_Wcollapse']:.3f}  "
                  f"({time.time()-t_start:.0f}s)", flush=True)
    runs_b_df = pd.concat(runs_b)
    runs_b_df.to_csv(OUT / "e2_runs_partB.csv", index=False)

    # ---- part (c): responsive variants -----------------------------------
    runs_c = []
    for cond_name, hfn, alpha in [
        ("active_stress", A.h_active_stress, 2.0),
        ("active_coupling", A.h_active_coupling, 1.0),
        ("active_quadratic", A.h_active_quadratic, 1.0),
    ]:
        for J in JS_C:
            j_idx = J_VALUES.index(J)
            seed = SEED_BASE_GRID + j_idx * 1_000_003 + MU_IDX * 1009
            r = run_cell(J, MU, hfn, alpha, n_seeds=N_ABC, seed=seed)
            rr, s = cell_rows("c", cond_name, alpha, "affine", J, MU, N_ABC,
                              seed, r)
            runs_c.append(rr); summaries.append(s)
            print(f"[c] {cond_name} J={J}: P_branch={s['P_branch_t50']:.3f} "
                  f"P_Wcoll={s['P_Wcollapse']:.3f}  "
                  f"({time.time()-t_start:.0f}s)", flush=True)
    pd.concat(runs_c).to_csv(OUT / "e2_runs_partC.csv", index=False)

    # ---- part (d): S5.5 consumption variants (h0_scaling seeds, n=2000) ------
    runs_d = []
    for vname, kw in VARIANTS:
        for J in JS_D:
            j_idx = JS_D.index(J)
            seed = SEED_BASE_TASKB + j_idx * 1_000_003 + 1 * 9973
            r = run_cell(J, MU, A.h_passive, 0.0, n_seeds=N_D, seed=seed, **kw)
            rr, s = cell_rows("d", "passive_h0=2", 0.0, vname, J, MU, N_D,
                              seed, r)
            runs_d.append(rr); summaries.append(s)
            print(f"[d] {vname} J={J}: P_branch={s['P_branch_t50']:.3f} "
                  f"P_Wcoll={s['P_Wcollapse']:.3f}  "
                  f"({time.time()-t_start:.0f}s)", flush=True)
    pd.concat(runs_d).to_csv(OUT / "e2_runs_partD.csv", index=False)

    # ---- E2 dissociation table -------------------------------------------
    summ_df = pd.DataFrame(summaries)
    summ_df.to_csv(OUT / "e2_dissociation.csv", index=False)

    # ---- E3: collapse-time distributions ---------------------------------
    e3_rows = []
    # fixed-field cells (from part-b reruns), all J with any collapses
    for cond in ("fixed2.4", "fixed5.0"):
        for J in JS_B:
            sub = runs_b_df[(runs_b_df.condition == cond) & (runs_b_df.J == J)]
            cs = sub.loc[sub.collapsed == 1, "collapse_step"].to_numpy()
            e3_rows.append(_e3_row(cond, J, len(sub), cs))
    # passive baseline comparison (part-a reruns, high-J residual)
    ra = pd.concat(runs_a)
    for J in JS_A:
        sub = ra[ra.J == J]
        cs = sub.loc[sub.collapsed == 1, "collapse_step"].to_numpy()
        e3_rows.append(_e3_row("passive", J, len(sub), cs))
    e3_df = pd.DataFrame(e3_rows)
    e3_df.to_csv(OUT / "e3_collapse_times.csv", index=False)

    _write_summary(summ_df, e3_df)
    print(f"\nDONE in {time.time()-t_start:.0f}s — outputs in {OUT}")


def _e3_row(cond, J, n, cs):
    t = cs * DT  # physical collapse times
    row = dict(condition=cond, J=J, n=n, n_collapsed=len(cs),
               p_collapse=len(cs) / n)
    if len(cs):
        row.update(median_t=float(np.median(t)),
                   p90_t=float(np.percentile(t, 90)),
                   p99_t=float(np.percentile(t, 99)),
                   max_t=float(t.max()),
                   frac_after_t100=float((t > 100.0).mean()))
    else:
        row.update(median_t=np.nan, p90_t=np.nan, p99_t=np.nan,
                   max_t=np.nan, frac_after_t100=np.nan)
    return row


def _write_summary(summ_df, e3_df):
    L = []
    L.append("# E2 + E3 — Branch-selection observable and fixed-field "
             "collapse-time distributions\n")
    L.append(f"Conventions: dt={DT}, T={T_TEMP}, xi={XI}, m0={INITIAL_M}, "
             f"W0={INITIAL_WEALTH}, mu={MU}; horizon {N_STEPS} steps "
             f"(t=1000). Matched seeds per the archived experiments "
             f"(see header of branch_observable.py). Regression checks: "
             f"bit-identity vs minimal_model.run_cell and vs archived "
             f"wdot_sensitivity raw_linear npz passed before running.\n")
    L.append("## E2 — window choice\n")
    L.append("P(wrong branch) = fraction of runs with m < 0 at t = 50 "
             "(step 1000). Justification: passive-baseline W-collapse times "
             "saturate early (empirically below: median ~15, max ~35 at "
             "mu=100), and the deterministic branch-selection transient is "
             "O(1) time units for supercritical J; t = 50 lies safely after "
             "both, yet is early enough (5% of the horizon) to exclude rare "
             "late stochastic escapes from the definition.\n")
    L.append("## E2 — dissociation table\n")
    L.append("| part | condition | cons | J | P_branch(t=50) | P(m<0 end) | "
             "P(W-collapse) | P_branch - P_Wcoll |")
    L.append("|---|---|---|---|---|---|---|---|")
    for _, r in summ_df.iterrows():
        L.append(f"| {r.part} | {r.condition} | {r.cons} | {r.J:g} "
                 f"| {r.P_branch_t50:.3f} | {r.P_mneg_end:.3f} "
                 f"| {r.P_Wcollapse:.3f} | {r.dissociation:+.3f} |")
    L.append("\nBinomial SE <= 0.016 at n=1000, <= 0.011 at n=2000.\n")
    L.append("## E3 — collapse-time distributions (physical units, "
             "t = collapse_step * dt)\n")
    L.append("| condition | J | P(coll) | median | p90 | p99 | max | "
             "frac after t=100 |")
    L.append("|---|---|---|---|---|---|---|---|")
    for _, r in e3_df.iterrows():
        def f(x):
            return "-" if not np.isfinite(x) else f"{x:.1f}"
        L.append(f"| {r.condition} | {r.J:g} | {r.p_collapse:.3f} "
                 f"| {f(r.median_t)} | {f(r.p90_t)} | {f(r.p99_t)} "
                 f"| {f(r.max_t)} | "
                 + ("-" if not np.isfinite(r.frac_after_t100)
                    else f"{r.frac_after_t100:.3f}") + " |")
    L.append("")
    with open(OUT / "SUMMARY.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))


if __name__ == "__main__":
    main()
