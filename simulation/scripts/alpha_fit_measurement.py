"""
Alpha-fit measurement of the collapse-observable scaling exponent.

Implements simulation/prereg/PREREG_alpha_fit_2026-05-30.md (frozen
2026-05-30, commit 5544672). NO toy: the passive-stabilizer condition is
imported verbatim from active_stabilizer.py (run_cell + h_passive), which
is the manuscript Methods pipeline.

Fit target: passive-stabilizer effectiveness E(J) = 1 - P(J,mu_high)/
P(J,mu_low) (Figure 8 Panel B), plus raw gap dP(J) = P(J,mu_low) -
P(J,mu_high). Form E(J) ~ (J - T)^(-alpha), T = J_c = 2.

Usage:
    python alpha_fit_measurement.py gate          # validation gate only
    python alpha_fit_measurement.py full [nseeds]  # gate, then full fit
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from active_stabilizer import run_cell, h_passive, T_TEMP, N_STEPS, SEED_BASE

# ---- Frozen pre-registered design ------------------------------------
J_GRID = (2.25, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0)
SUBWINDOW = (6.0, 7.0, 8.0, 9.0, 10.0, 12.0)          # high-J, J >= 6
MU_VALUES = (20, 40, 60, 100)                          # 20 = saturating control
MU_HIGH = 100
MU_LOW_PRIMARY = 40
MU_LOW_ALT = (20, 60)
T = T_TEMP                                             # J_c = T = 2.0
N_SEEDS_DEFAULT = 800                                  # floor >= 200; escalate if SE(alpha) > 0.05
N_BOOT = 4000
GATE_J, GATE_MU, GATE_TARGET, GATE_TOL = 5.0, 100, 0.30, 0.03


def cell_seed(j_idx: int, m_idx: int) -> int:
    # Same matched-seed scheme as active_stabilizer.sweep_condition.
    return SEED_BASE + j_idx * 1_000_003 + m_idx * 1009


def p_collapse(J: float, mu: int, j_idx: int, m_idx: int, n_seeds: int):
    res = run_cell(float(J), mu, h_passive, 0.0, n_seeds=n_seeds,
                   n_steps=N_STEPS, seed=cell_seed(j_idx, m_idx))
    return res.collapsed.astype(np.float64)            # per-seed 0/1 array


def run_gate(n_seeds: int):
    ji = J_GRID.index(GATE_J)
    mi = MU_VALUES.index(GATE_MU)
    arr = p_collapse(GATE_J, GATE_MU, ji, mi, n_seeds)
    p = arr.mean()
    lo, hi = GATE_TARGET - GATE_TOL, GATE_TARGET + GATE_TOL
    ok = lo <= p <= hi
    print(f"[GATE] P(J={GATE_J}, mu={GATE_MU}) = {p:.4f}  "
          f"(n={n_seeds}; target [{lo:.2f},{hi:.2f}]) -> "
          f"{'PASS' if ok else 'FAIL'}")
    return p, ok, arr


def ols_alpha(J_arr, E_arr):
    """alpha from log E = c - alpha*log(J-T). Returns (alpha, r2, n_used)."""
    x = np.log(np.asarray(J_arr) - T)
    E = np.asarray(E_arr, dtype=np.float64)
    mask = np.isfinite(E) & (E > 0)
    if mask.sum() < 3:
        return np.nan, np.nan, int(mask.sum())
    xv, yv = x[mask], np.log(E[mask])
    slope, intercept = np.polyfit(xv, yv, 1)
    yhat = slope * xv + intercept
    ss_res = float(((yv - yhat) ** 2).sum())
    ss_tot = float(((yv - yv.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return -slope, r2, int(mask.sum())


def effectiveness(P, mu_low, J_sub):
    return np.array([
        (1.0 - P[(J, MU_HIGH)].mean() / P[(J, mu_low)].mean())
        if P[(J, mu_low)].mean() > 0 else np.nan
        for J in J_sub
    ])


def raw_gap(P, mu_low, J_sub):
    return np.array([P[(J, mu_low)].mean() - P[(J, MU_HIGH)].mean()
                     for J in J_sub])


def bootstrap_alpha(P, mu_low, J_sub, target="E", n_boot=N_BOOT, seed=0xA1FA):
    rng = np.random.default_rng(seed)
    alphas = []
    for _ in range(n_boot):
        Pb = {}
        for J in J_sub:
            for mu in (mu_low, MU_HIGH):
                a = P[(J, mu)]
                idx = rng.integers(0, a.size, a.size)
                Pb[(J, mu)] = a[idx].mean()
        if target == "E":
            vals = np.array([
                (1.0 - Pb[(J, MU_HIGH)] / Pb[(J, mu_low)])
                if Pb[(J, mu_low)] > 0 else np.nan for J in J_sub])
        else:
            vals = np.array([Pb[(J, mu_low)] - Pb[(J, MU_HIGH)] for J in J_sub])
        a, _, _ = ols_alpha(J_sub, vals)
        if np.isfinite(a):
            alphas.append(a)
    alphas = np.array(alphas)
    if alphas.size == 0:
        return dict(mean=np.nan, se=np.nan, lo=np.nan, hi=np.nan, n=0)
    return dict(mean=float(alphas.mean()), se=float(alphas.std(ddof=1)),
                lo=float(np.percentile(alphas, 2.5)),
                hi=float(np.percentile(alphas, 97.5)), n=int(alphas.size))


def local_slopes(J_sub, vals):
    x = np.log(np.asarray(J_sub) - T)
    y = np.asarray(vals, dtype=np.float64)
    out = []
    for i in range(len(J_sub) - 1):
        if np.isfinite(y[i]) and np.isfinite(y[i + 1]) and y[i] > 0 and y[i + 1] > 0:
            a = -(np.log(y[i + 1]) - np.log(y[i])) / (x[i + 1] - x[i])
        else:
            a = np.nan
        out.append((J_sub[i], J_sub[i + 1], a))
    return out


def classify(full_ci, sub_ci, r2_full, r2_sub):
    """Pre-registered decision rule (PREREG §5)."""
    def contains(ci, v):
        return ci["lo"] <= v <= ci["hi"]

    notes = []
    # alpha = 2 excluded a priori; flag if a CI nonetheless contains it.
    if contains(full_ci, 2.0) or contains(sub_ci, 2.0):
        notes.append("UNEXPECTED: a CI contains 2 (excluded a priori) -> surface")
    half = contains(full_ci, 0.5) and not contains(full_ci, 1.0)
    one = contains(full_ci, 1.0) and not contains(full_ci, 0.5)
    half_sub = contains(sub_ci, 0.5) and not contains(sub_ci, 1.0)
    one_sub = contains(sub_ci, 1.0) and not contains(sub_ci, 0.5)
    drift = ((half and one_sub) or (one and half_sub)
             or (not np.isclose(full_ci["mean"], sub_ci["mean"], atol=0.25)
                 and min(r2_full, r2_sub) > 0.8))
    if (r2_full is not None and r2_full < 0.7) and (r2_sub is not None and r2_sub < 0.7):
        return "NO-CLEAN-EXPONENT (low R^2 both windows) -> qualitative reframe", notes
    if drift:
        return "PRE-ASYMPTOTIC DRIFT -> resolve to J>=6 subwindow value", notes
    if half:
        return "CONSISTENT-WITH-1/2 -> resolve held cohort to h/sqrt(J)", notes
    if one:
        return "CONSISTENT-WITH-1 -> resolve held cohort to h/J", notes
    return "AMBIGUOUS (CI spans 1/2 and 1) -> report CI, no auto-resolution", notes


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "gate"
    n_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else N_SEEDS_DEFAULT

    if mode == "gate":
        run_gate(n_seeds)
        return

    # full: gate first, halt-and-surface on failure (no tuning).
    t0 = time.time()
    p_gate, ok, _ = run_gate(n_seeds)
    if not ok:
        print("[HALT] Validation gate failed. Per pre-reg discipline: "
              "NOT tuning, NOT proceeding to the alpha fit. Surfacing.")
        sys.exit(2)

    OUT = Path(__file__).resolve().parents[1] / "results" / "alpha_fit"
    OUT.mkdir(parents=True, exist_ok=True)

    # Full grid x margins.
    P = {}
    rows = []
    for ji, J in enumerate(J_GRID):
        for mi, mu in enumerate(MU_VALUES):
            arr = p_collapse(J, mu, ji, mi, n_seeds)
            P[(J, mu)] = arr
            rows.append(dict(J=J, mu=mu, p_collapse=arr.mean(), n=arr.size))
    grid_df = pd.DataFrame(rows)
    grid_df.to_csv(OUT / "p_collapse_grid.csv", index=False)
    print(f"\n[GRID] P(collapse) by (J, mu), n={n_seeds}/cell "
          f"({time.time()-t0:.1f}s)")
    print(grid_df.pivot(index="J", columns="mu", values="p_collapse")
          .round(4).to_string())

    pairings = [("primary 40/100", MU_LOW_PRIMARY),
                ("control 20/100 (SATURATING)", 20),
                ("alt 60/100", 60)]

    summary_rows = []
    for label, mu_low in pairings:
        E_full = effectiveness(P, mu_low, J_GRID)
        E_sub = effectiveness(P, mu_low, SUBWINDOW)
        dP_full = raw_gap(P, mu_low, J_GRID)
        a_full, r2_full, n_full = ols_alpha(J_GRID, E_full)
        a_sub, r2_sub, n_sub = ols_alpha(SUBWINDOW, E_sub)
        adp_full, r2dp_full, _ = ols_alpha(J_GRID, dP_full)
        bf = bootstrap_alpha(P, mu_low, J_GRID, "E")
        bs = bootstrap_alpha(P, mu_low, SUBWINDOW, "E")
        print(f"\n=== {label} ===")
        print(f"  E(J) full  : alpha={a_full:.3f} R^2={r2_full:.3f} "
              f"(n_used={n_full}); boot mean={bf['mean']:.3f} "
              f"SE={bf['se']:.3f} 95%CI[{bf['lo']:.3f},{bf['hi']:.3f}]")
        print(f"  E(J) J>=6  : alpha={a_sub:.3f} R^2={r2_sub:.3f} "
              f"(n_used={n_sub}); boot mean={bs['mean']:.3f} "
              f"SE={bs['se']:.3f} 95%CI[{bs['lo']:.3f},{bs['hi']:.3f}]")
        print(f"  dP(J) full : alpha={adp_full:.3f} R^2={r2dp_full:.3f}")
        print(f"  E(J) values: "
              + ", ".join(f"{J}:{e:.3f}" for J, e in zip(J_GRID, E_full)))
        if label.startswith("primary"):
            ls = local_slopes(J_GRID, E_full)
            print("  model-free local slopes (J_i->J_{i+1} : alpha):")
            print("    " + ", ".join(
                f"{a}->{b}:{('%.2f' % s) if np.isfinite(s) else 'nan'}"
                for a, b, s in ls))
            cls, notes = classify(bf, bs, r2_full, r2_sub)
            print(f"  [DECISION] {cls}")
            for nt in notes:
                print(f"  [FLAG] {nt}")
            if bf["se"] > 0.05 or bs["se"] > 0.05:
                print(f"  [NOTE] bootstrap SE(alpha) full={bf['se']:.3f} "
                      f"sub={bs['se']:.3f}; target <=0.05 -> "
                      f"escalate seeds if either exceeds 0.05.")
        summary_rows.append(dict(
            pairing=label, mu_low=mu_low,
            alpha_E_full=a_full, r2_E_full=r2_full, n_E_full=n_full,
            boot_mean_full=bf["mean"], boot_se_full=bf["se"],
            ci_lo_full=bf["lo"], ci_hi_full=bf["hi"],
            alpha_E_sub=a_sub, r2_E_sub=r2_sub, n_E_sub=n_sub,
            boot_mean_sub=bs["mean"], boot_se_sub=bs["se"],
            ci_lo_sub=bs["lo"], ci_hi_sub=bs["hi"],
            alpha_dP_full=adp_full, r2_dP_full=r2dp_full))
    pd.DataFrame(summary_rows).to_csv(OUT / "alpha_fit_summary.csv", index=False)
    print(f"\nWrote {OUT/'p_collapse_grid.csv'} and "
          f"{OUT/'alpha_fit_summary.csv'} ({time.time()-t0:.1f}s total)")


if __name__ == "__main__":
    main()
