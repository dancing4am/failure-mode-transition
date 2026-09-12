"""Ordering test for candidate delivered-amplitude summaries of a stabilizer.

Contribution 2 asks whether any single delivered-amplitude statistic orders the
five static-sweep arms (passive; stress alpha=2.0 and 0.05; coupling; quadratic)
the same way P(collapse) does. This script measures, for each arm at
J in {4.0, 4.5, 5.0}, mu=100, n=2,000 matched seeds:
  (i)   full-trajectory delivered mean
  (ii)  selection-window delivered mean, t <= 1
  (iii) delivered mean conditional on m < 0, over t <= 1
  (iv)  (iii) restricted to runs that eventually collapse
  (v)   (ii) / lambda, lambda = -1 + (J/T) sech^2(h_eff/T) with h_eff = (ii)
        (amplitude relative to the origin-instability rate; lambda<=0 -> stable)
and P(collapse). It then reports, per statistic, whether it orders the five arms
consistently with P(collapse) at every J (PASS/FAIL, with the first inversion).

Field forms identical to delivered_mean_windows.py / Methods. Post-hoc
diagnostic (not part of the two-arm pre-registered plan).

Output: simulation/results/magnitude_matched/delivered_mean_ordering.csv
"""
from __future__ import annotations

from pathlib import Path
import time

import numpy as np
import pandas as pd

T_TEMP, XI, DT = 2.0, 0.5, 0.05
N_STEPS = 20_000
INITIAL_M, INITIAL_WEALTH = 0.0, 100.0
COLLAPSE_WEALTH, COLLAPSE_STREAK = 10.0, 200
EPS = 1e-6
MU = 100.0
N_SEEDS = 2000
SEED_BASE = 0xDE_11A7          # same base as delivered_mean_windows for pairing
CELLS_J = (4.0, 4.5, 5.0)
SELWIN_STEPS = int(round(1.0 / DT))     # t <= 1 -> first 20 steps
ARMS = (("passive", 0.0), ("stress_a2.0", 2.0), ("stress_a0.05", 0.05),
        ("coupling", 1.0), ("quadratic", 1.0))
OUT = Path(__file__).resolve().parents[1] / "results" / "magnitude_matched"


def field(arm, alpha, wealth, m_c, J):
    base = 2.0 * (wealth / 500.0)
    if arm == "passive":
        return base
    if arm.startswith("stress"):
        return base + alpha * np.maximum(0.0, -m_c) * J
    if arm == "coupling":
        return base * (1.0 + alpha * J)
    if arm == "quadratic":
        return base * (1.0 + alpha * (J / 5.0) ** 2)
    raise ValueError(arm)


def run_arm(arm, alpha, J, seed, n=N_SEEDS):
    rng = np.random.default_rng(seed)
    m = np.full(n, INITIAL_M)
    wealth = np.full(n, INITIAL_WEALTH)
    collapsed = np.zeros(n, dtype=bool)
    streak = np.zeros(n, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)
    full_sum = 0.0
    win_sum = 0.0
    win_n = 0
    neg_sum = np.zeros(n)          # sum of h where m<0 in t<=1, per seed
    neg_cnt = np.zeros(n, dtype=np.int64)
    for t in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = field(arm, alpha, wealth, m_c, J)
        full_sum += float(h.sum())
        if t < SELWIN_STEPS:
            win_sum += float(h.sum())
            win_n += n
            neg = m_c < 0.0
            neg_sum += np.where(neg, h, 0.0)
            neg_cnt += neg
        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)
        income = MU * (0.5 + 0.5 * m)
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)
        streak = np.where(wealth < COLLAPSE_WEALTH, streak + 1, 0)
        collapsed |= streak >= COLLAPSE_STREAK
    p = float(collapsed.mean())
    full_mean = full_sum / (N_STEPS * n)
    win_mean = win_sum / win_n
    cond_all = neg_sum.sum() / neg_cnt.sum() if neg_cnt.sum() else float("nan")
    csum = neg_sum[collapsed].sum(); ccnt = neg_cnt[collapsed].sum()
    cond_coll = csum / ccnt if ccnt else float("nan")
    lam = -1.0 + (J / T_TEMP) / np.cosh(win_mean / T_TEMP) ** 2
    rate_rel = (win_mean / lam) if lam > 0 else float("nan")
    return dict(p_collapse=p, full_mean=full_mean, win_mean=win_mean,
                cond_m_neg=cond_all, cond_m_neg_collapsers=cond_coll,
                lam_at_hwin=lam, amp_over_lam=rate_rel,
                in_domain=bool(lam > 0))


def orders_consistently(df, stat, in_domain_only=False):
    """Does `stat`, sorted ascending, give P(collapse) non-increasing at every J?
    Return (bool, first_inversion_str). With in_domain_only, restrict to arms
    whose origin is linearly unstable (lambda>0) at the measured window field."""
    for J in CELLS_J:
        d = df[df.J == J]
        if in_domain_only:
            d = d[d.in_domain]
        d = d.sort_values(stat)
        ps = d.p_collapse.values
        arms = d.arm.values
        for i in range(len(ps) - 1):
            if ps[i + 1] > ps[i] + 1e-9:
                return False, (f"J={J}: {arms[i]}({d[stat].values[i]:.3f}->P={ps[i]:.3f}) "
                               f"< {arms[i+1]}({d[stat].values[i+1]:.3f}->P={ps[i+1]:.3f})")
    return True, ""


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for J in CELLS_J:
        seed = SEED_BASE + int(round(J * 10)) * 1_000_003
        for arm, alpha in ARMS:
            t0 = time.time()
            r = run_arm(arm, alpha, J, seed)
            r.update(J=J, arm=arm, alpha=alpha)
            rows.append(r)
            print(f"J={J} {arm:13s} P={r['p_collapse']:.4f}  full={r['full_mean']:.3f} "
                  f"win<=1={r['win_mean']:.3f}  cond(m<0)={r['cond_m_neg']:.3f} "
                  f"cond|coll={r['cond_m_neg_collapsers']:.3f}  lam={r['lam_at_hwin']:.3f} "
                  f"amp/lam={r['amp_over_lam']:.3f}  ({time.time()-t0:.0f}s)")
    df = pd.DataFrame(rows)
    cols = ["J", "arm", "alpha", "p_collapse", "full_mean", "win_mean",
            "cond_m_neg", "cond_m_neg_collapsers", "lam_at_hwin", "amp_over_lam",
            "in_domain"]
    df[cols].to_csv(OUT / "delivered_mean_ordering.csv", index=False)
    print(f"\nWrote {OUT / 'delivered_mean_ordering.csv'}")
    stats = [("full_mean", "(i) full-trajectory mean"),
             ("win_mean", "(ii) selection-window mean (t<=1)"),
             ("cond_m_neg", "(iii) mean | m<0, t<=1"),
             ("cond_m_neg_collapsers", "(iv) mean | m<0, collapsers, t<=1"),
             ("amp_over_lam", "(v) window mean / lambda")]
    print("\n=== SECONDARY ordering test (ALL 5 arms; out-of-domain coupling included) ===")
    for stat, label in stats:
        ok, inv = orders_consistently(df, stat)
        print(f"  {label:38s}: {'PASS' if ok else 'FAIL'}" + ("" if ok else f"  first inversion: {inv}"))
    print("\n=== PRIMARY ordering test (in-domain arms only; lambda>0 at window field) ===")
    in_dom = sorted(df[df.in_domain].arm.unique())
    print(f"  in-domain arms: {in_dom}")
    for stat, label in stats:
        ok, inv = orders_consistently(df, stat, in_domain_only=True)
        print(f"  {label:38s}: {'PASS' if ok else 'FAIL'}" + ("" if ok else f"  first inversion: {inv}"))


if __name__ == "__main__":
    main()
