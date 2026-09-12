"""F0 discretization-convergence check for the late collapse channel.

Reruns the paired clipped-EM + arcsin-Lamperti comparison of
late_channel_check.py at HALVED step size dt = 0.025 (n_steps = 40,000,
same physical horizon t_max = 1000) and compares against the archived
dt = 0.05 results (late_channel_summary.csv, n=2000/cell).

Settings (identical to late_channel_check except dt and n):
  * cells (h, J) in row-major order [(2.4,15), (2.4,20), (5.0,15), (5.0,20)],
    mu=100, constant h, T=2, xi=0.5, m0=0, W0=100;
  * paired dW per (seed, step) fed to BOTH integrators;
  * n = 1000 trajectories per cell;
  * seed_i = 0xF0_0001 + i * 1_000_003 (FRESH seeds — paths are NOT shared
    across dt, so the dt comparison is cross-sample, not paired);
  * collapse = W < 10 for 400 consecutive steps at dt=0.025.  The archived
    run used 200 steps at dt=0.05; the streak length is scaled by 1/dt
    (200 * 0.05/dt) so the collapse definition stays "W < 10 sustained for
    10 PHYSICAL time units" at every dt.

Cross-sample dt=0.025 vs dt=0.05 comparison, per scheme per cell:
  * P(collapse): two-sample binomial diff, unpooled 95% Wald CI, pooled z;
  * late fraction (t > 100 among collapsed): same, with n = n_collapsed;
  * timing medians/p99 reported descriptively.

Verdict: CONVERGED iff every dt difference (per scheme, per cell, for P
and late fraction) has 0 inside its two-sample 95% CI.

Outputs:
  simulation/results/boundary_integrity/late_channel/dt_convergence.csv
  appended section in simulation/results/boundary_integrity/late_channel/SUMMARY_late_channel.md
"""

from __future__ import annotations

import csv
import os
import time

import numpy as np

from late_channel_check import (
    run_cell_paired, CELLS, MU, LATE_T, OUTDIR, SUMDIR,
)
from boundary_integrity import COLLAPSE_STREAK, DT as DT_REF, N_STEPS as N_STEPS_REF

DT_FINE = 0.025
N_STEPS_FINE = 40_000                     # same horizon: 0.025 * 40000 = 1000
STREAK_FINE = int(round(COLLAPSE_STREAK * DT_REF / DT_FINE))   # = 400
N_TRAJ_FINE = 1000
SEED_BASE_F0 = 0xF0_0001                  # = 15728641

ARCHIVE_CSV = os.path.join(OUTDIR, "late_channel_summary.csv")
OUT_CSV = os.path.join(OUTDIR, "dt_convergence.csv")
SUMMARY_MD = os.path.join(SUMDIR, "SUMMARY_late_channel.md")

ARM_KEY = {"clipped_EM": "em", "lamperti": "lam"}


def load_archive():
    """Archived dt=0.05 per-cell x per-scheme counts and timing stats."""
    ref = {}
    with open(ARCHIVE_CSV, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            key = (float(row["h"]), float(row["J"]), row["integrator"])
            ref[key] = dict(
                n=int(row["n"]),
                n_collapsed=int(row["n_collapsed"]),
                n_late=int(row["n_late_t100"]),
                p=float(row["p_collapse"]),
                late_frac=float(row["late_frac"]) if row["late_frac"] else np.nan,
                median_t=float(row["median_t"]) if row["median_t"] else np.nan,
                p99_t=float(row["p99_t"]) if row["p99_t"] else np.nan,
            )
    return ref


def arm_stats_dt(collapsed, steps, dt):
    """P(collapse), late count/fraction, timing stats at arbitrary dt."""
    col = collapsed
    t = steps[col].astype(float) * dt
    n_col = int(col.sum())
    n_late = int((t > LATE_T).sum())
    if n_col:
        med, p90, p99, tmax = (float(np.median(t)),
                               float(np.percentile(t, 90)),
                               float(np.percentile(t, 99)), float(t.max()))
    else:
        med = p90 = p99 = tmax = np.nan
    return dict(p=float(col.mean()), n=len(col), n_collapsed=n_col,
                n_late=n_late,
                late_frac=(n_late / n_col if n_col else np.nan),
                median_t=med, p90_t=p90, p99_t=p99, max_t=tmax)


def two_sample(k1, n1, k2, n2, z_crit=1.96):
    """diff = p1 - p2 with unpooled 95% Wald CI and pooled z."""
    if n1 == 0 or n2 == 0:
        return (np.nan,) * 5
    p1, p2 = k1 / n1, k2 / n2
    diff = p1 - p2
    se_u = np.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    pp = (k1 + k2) / (n1 + n2)
    se_p = np.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
    z = diff / se_p if se_p > 0 else np.nan
    return diff, diff - z_crit * se_u, diff + z_crit * se_u, z, se_u


def fmt(x, d=3):
    return "-" if (x is None or (isinstance(x, float) and np.isnan(x))) \
        else f"{x:.{d}f}"


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    ref = load_archive()
    t0 = time.time()

    rows = []
    drift_cells = []
    for i, (h, J) in enumerate(CELLS):
        seed = SEED_BASE_F0 + i * 1_000_003
        print(f"cell {i}: h={h} J={J:g} seed={seed} dt={DT_FINE} "
              f"steps={N_STEPS_FINE} streak={STREAK_FINE} ...", flush=True)
        r = run_cell_paired(h, J, N_TRAJ_FINE, seed,
                            n_steps=N_STEPS_FINE, dt=DT_FINE,
                            collapse_streak=STREAK_FINE)
        print(f"  done ({time.time() - t0:.0f}s)", flush=True)

        for arm_name, key in ARM_KEY.items():
            a = r[key]
            st = arm_stats_dt(a["collapsed"], a["collapse_step"], DT_FINE)
            rf = ref[(h, J, arm_name)]

            dP, dP_lo, dP_hi, zP, _ = two_sample(
                st["n_collapsed"], st["n"], rf["n_collapsed"], rf["n"])
            dLF, dLF_lo, dLF_hi, zLF, _ = two_sample(
                st["n_late"], st["n_collapsed"], rf["n_late"],
                rf["n_collapsed"])

            ok_P = (not np.isnan(dP_lo)) and dP_lo <= 0.0 <= dP_hi
            ok_LF = (not np.isnan(dLF_lo)) and dLF_lo <= 0.0 <= dLF_hi
            if not (ok_P and ok_LF):
                what = []
                if not ok_P:
                    what.append(f"P {dP:+.4f}")
                if not ok_LF:
                    what.append(f"late-frac {dLF:+.4f}")
                drift_cells.append((h, J, arm_name, ", ".join(what)))

            diag = (dict(clip_step_rate=a["clip_step_rate"])
                    if key == "em" else
                    dict(cap_rate=a["cap_rate"], fold_rate=a["fold_rate"]))
            rows.append(dict(
                cell=f"h{h}_J{J:g}", h=h, J=J, mu=MU, scheme=arm_name,
                seed_fine=seed,
                dt_fine=DT_FINE, n_steps_fine=N_STEPS_FINE,
                streak_fine=STREAK_FINE, n_fine=st["n"],
                p_fine=st["p"], n_collapsed_fine=st["n_collapsed"],
                n_late_fine=st["n_late"], late_frac_fine=st["late_frac"],
                median_t_fine=st["median_t"], p90_t_fine=st["p90_t"],
                p99_t_fine=st["p99_t"], max_t_fine=st["max_t"],
                dt_ref=DT_REF, n_ref=rf["n"], p_ref=rf["p"],
                n_collapsed_ref=rf["n_collapsed"], n_late_ref=rf["n_late"],
                late_frac_ref=rf["late_frac"], median_t_ref=rf["median_t"],
                p99_t_ref=rf["p99_t"],
                dP=dP, dP_ci_lo=dP_lo, dP_ci_hi=dP_hi, z_P=zP,
                dLateFrac=dLF, dLateFrac_ci_lo=dLF_lo,
                dLateFrac_ci_hi=dLF_hi, z_LateFrac=zLF,
                P_within_CI=int(ok_P), LateFrac_within_CI=int(ok_LF),
                **diag))

    converged = not drift_cells
    verdict = "CONVERGED" if converged else "NOT CONVERGED"

    # -------------------------------------------------------------- CSV
    cols = list(dict.fromkeys(k for row in rows for k in row))
    with open(OUT_CSV, "w", encoding="utf-8") as fh:
        fh.write(",".join(cols) + "\n")
        for row in rows:
            fh.write(",".join(
                "" if (k not in row or row[k] is None or
                       (isinstance(row[k], float) and np.isnan(row[k])))
                else (f"{row[k]:.6g}" if isinstance(row[k], float)
                      else str(row[k]))
                for k in cols) + "\n")

    # ------------------------------------------------- SUMMARY.md append
    L = []
    A = L.append
    A("")
    A("---")
    A("")
    A("# Discretization convergence of the late channel (dt 0.05 -> 0.025)")
    A("")
    A("**Question.** Do the late-channel results above survive halving the")
    A("step size?  Rerun of the paired EM + Lamperti comparison at dt=0.025,")
    A(f"n_steps={N_STEPS_FINE} (same physical horizon t_max=1000), n={N_TRAJ_FINE}")
    A("per cell, paired dW per (seed, step), same four cells, mu=100,")
    A("constant h, T=2, xi=0.5, m0=0, W0=100.  The collapse streak is scaled")
    A(f"with 1/dt — {STREAK_FINE} consecutive steps at dt=0.025 (= 200 at")
    A("dt=0.05) — so collapse remains 'W<10 sustained for 10 physical time")
    A(f"units' at every dt.  Seeds 0xF00001 + i*1000003 (fresh; paths are NOT")
    A("shared across dt, so the dt comparison is CROSS-SAMPLE: two-sample")
    A("binomial z / unpooled 95% Wald CI on the differences).  Script:")
    A("simulation/scripts/late_channel_dt_convergence.py; data:")
    A("simulation/results/boundary_integrity/late_channel/dt_convergence.csv.")
    A("")
    A("## Per-cell comparison (scheme x dt)")
    A("")
    A("| h | J | scheme | dt | n | P(collapse) | late frac | median t | p99 t |")
    A("|---|---|---|---|---|---|---|---|---|")
    for row in rows:
        A(f"| {row['h']} | {row['J']:g} | {row['scheme']} | 0.05 "
          f"| {row['n_ref']} | {row['p_ref']:.4f} "
          f"| {fmt(row['late_frac_ref'])} | {fmt(row['median_t_ref'], 1)} "
          f"| {fmt(row['p99_t_ref'], 1)} |")
        A(f"| {row['h']} | {row['J']:g} | {row['scheme']} | 0.025 "
          f"| {row['n_fine']} | {row['p_fine']:.4f} "
          f"| {fmt(row['late_frac_fine'])} | {fmt(row['median_t_fine'], 1)} "
          f"| {fmt(row['p99_t_fine'], 1)} |")
    A("")
    A("## dt differences (dt=0.025 minus dt=0.05, two-sample 95% CI)")
    A("")
    A("| h | J | scheme | dP [95% CI] | z(P) | dLateFrac [95% CI] | z(LF) | within CI? |")
    A("|---|---|---|---|---|---|---|---|")
    for row in rows:
        oks = ("P yes" if row["P_within_CI"] else "P NO") + ", " + \
              ("LF yes" if row["LateFrac_within_CI"] else "LF NO")
        A(f"| {row['h']} | {row['J']:g} | {row['scheme']} "
          f"| {row['dP']:+.4f} [{row['dP_ci_lo']:+.4f}, {row['dP_ci_hi']:+.4f}] "
          f"| {row['z_P']:+.2f} "
          f"| {row['dLateFrac']:+.4f} [{row['dLateFrac_ci_lo']:+.4f}, "
          f"{row['dLateFrac_ci_hi']:+.4f}] | {row['z_LateFrac']:+.2f} "
          f"| {oks} |")
    A("")
    A("## Verdict")
    A("")
    if converged:
        A(f"**{verdict}.** All dt=0.025 vs dt=0.05 differences (both schemes,")
        A("all four cells, P(collapse) and late fraction) lie within their")
        A("two-sample 95% CIs; the late channel and its magnitude are not a")
        A("dt=0.05 discretization artifact at this resolution.")
    else:
        A(f"**{verdict}.** Cells outside the two-sample 95% CI:")
        A("")
        for h, J, arm, what in drift_cells:
            A(f"- (h={h}, J={J:g}, {arm}): {what} (dt=0.025 minus dt=0.05)")
    A("")
    A(f"Runtime: {time.time() - t0:.0f}s.")
    A("")

    with open(SUMMARY_MD, "a", encoding="utf-8") as fh:
        fh.write("\n".join(L))

    print(f"Verdict: {verdict}")
    if drift_cells:
        for h, J, arm, what in drift_cells:
            print(f"  drift: h={h} J={J:g} {arm}: {what}")
    print(f"Done in {time.time() - t0:.0f}s")
    print(f"  {OUT_CSV}")
    print(f"  {SUMMARY_MD} (appended)")


if __name__ == "__main__":
    main()
