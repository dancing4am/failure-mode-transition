"""Mechanism test for the late collapse channel: activated (Kramers) escape?

Pre-registered in:
  simulation/prereg/PRESPEC_late_channel_mechanism_2026-08-22.md

Fixed-field cells (constant h, wealth ledger driven by m, collapse = W<10 for
200 consecutive steps), boundary-preserving Lamperti integrator (arcsin
transform, exact reflecting fold, trust-region drift cap), reused verbatim from
boundary_integrity.  Base settings: dt=0.05, 20000 steps, T=2, xi=0.5, mu=100,
W0=100.

L1  initial-condition dependence  (m0=0 vs m0=+0.9, same seeds)
L2  waiting-time memorylessness   (CV, KS vs exponential, Q-Q data)
L3  barrier scaling               (log rate vs deterministic + effective barrier)

Outputs -> simulation/results/late_channel_mechanism/
"""
from __future__ import annotations

import os
import time
import numpy as np
from scipy import stats, integrate, optimize

from boundary_integrity import (
    T_TEMP, XI, DT, N_STEPS, INITIAL_WEALTH, COLLAPSE_WEALTH, COLLAPSE_STREAK,
    DRIFT_CAP,
)

MU = 100.0
N_TRAJ = 2000
HS = (2.4, 5.0)
JS = (7.0, 10.0, 15.0, 20.0)
CELLS = [(h, J) for h in HS for J in JS]          # row-major, 8 cells
SEED_BASE = 0xB0_0001
EARLY_T = 50.0
LATE_T = 250.0
T_MAX = N_STEPS * DT                               # 1000
OUTDIR = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "results", "late_channel_mechanism"))


# ----------------------------------------------------------------------------
# Deterministic drift, fixed points, barriers
# ----------------------------------------------------------------------------
def drift(m, J, h):
    return -m + np.tanh((J * m + h) / T_TEMP)


def dfdm(m, J, h):
    return -1.0 + (J / T_TEMP) / np.cosh((J * m + h) / T_TEMP) ** 2


def fixed_points(J, h):
    """Return sorted roots of drift on (-1,1) with stability flags.

    For large J the stable wells sit within ~1e-7 of +-1, so the grid is
    geometrically refined toward both boundaries to bracket them.
    """
    edges_hi = 1.0 - np.logspace(-1.0, -13.0, 400)     # 0.9 ... 1-1e-13
    edges_lo = -1.0 + np.logspace(-1.0, -13.0, 400)
    grid = np.unique(np.concatenate([
        np.linspace(-0.9, 0.9, 180001), edges_hi, edges_lo]))
    fg = drift(grid, J, h)
    sign = np.sign(fg)
    idx = np.where(np.diff(sign) != 0)[0]
    roots = []
    for i in idx:
        r = optimize.brentq(drift, grid[i], grid[i + 1], args=(J, h))
        roots.append(r)
    roots = sorted(roots)
    return [(r, "stable" if dfdm(r, J, h) < 0 else "unstable") for r in roots]


def barriers(J, h):
    """Deltas from high-m well up to the saddle just below it.

    Returns (m_high, m_saddle, dU, G) where
      dU = U(saddle) - U(well),  U(m) = int (m - tanh((Jm+h)/T)) dm   [dU>0]
      G  = Ghat(well) - Ghat(saddle),  Ghat(m) = int f/(1-m^2) dm     [G>0]
    Effective escape exponent = (2/xi^2) * G  (Ito multiplicative noise).
    """
    fps = fixed_points(J, h)
    stable = [r for r, s in fps if s == "stable"]
    unstable = [r for r, s in fps if s == "unstable"]
    m_high = max(stable)
    # saddle just below the high well (if monostable, no saddle -> None)
    below = [u for u in unstable if u < m_high]
    if not below:
        return m_high, None, np.nan, np.nan
    m_s = max(below)
    dU, _ = integrate.quad(lambda m: -drift(m, J, h), m_high, m_s)   # U(s)-U(high)
    G, _ = integrate.quad(lambda m: drift(m, J, h) / (1.0 - m * m), m_s, m_high)
    return m_high, m_s, dU, G


# ----------------------------------------------------------------------------
# Boundary-preserving (Lamperti) integrator, constant h, wealth ledger
# ----------------------------------------------------------------------------
def lamperti_run(h, J, n, seed, m0, m_saddle, xi=XI, n_steps=N_STEPS, dt=DT,
                 collapse_streak=COLLAPSE_STREAK):
    """One Lamperti ensemble at constant h.  Records, per trajectory:
       collapsed, collapse_step, escape_step (first m<saddle after entering the
       well), and alive_at_late (uncollapsed at t=LATE_T).
    """
    rng = np.random.default_rng(seed)
    sqrt_dt = np.sqrt(dt)
    u_bound = np.pi / (2.0 * xi)
    late_step = int(round(LATE_T / dt))

    u = np.full(n, np.arcsin(np.clip(m0, -1 + 1e-9, 1 - 1e-9)) / xi)
    w = np.full(n, INITIAL_WEALTH)
    col = np.zeros(n, dtype=bool)
    cstep = np.full(n, -1, dtype=np.int64)
    stk = np.zeros(n, dtype=np.int64)
    esc = np.full(n, -1, dtype=np.int64)
    entered = np.zeros(n, dtype=bool)          # has reached above the saddle
    cap_events = fold_events = 0

    enter_gate = m_saddle if m_saddle is not None else 0.0
    for t in range(n_steps):
        dW = sqrt_dt * rng.standard_normal(n)
        s = np.sin(xi * u)                     # = m
        c = np.cos(xi * u)                     # = sqrt(1-m^2) on |u|<=u_bound
        f = -s + np.tanh((J * s + h) / T_TEMP)
        drift_dt = (f + 0.5 * xi * xi * s) / (xi * np.maximum(c, 1e-12)) * dt
        over = np.abs(drift_dt) > DRIFT_CAP
        cap_events += int(over.sum())
        drift_dt = np.clip(drift_dt, -DRIFT_CAP, DRIFT_CAP)
        u = u + drift_dt + dW
        out = (u > u_bound) | (u < -u_bound)
        fold_events += int(out.sum())
        if out.any():
            width = 2.0 * u_bound
            z = np.mod(u + u_bound, 2.0 * width)
            z = np.where(z > width, 2.0 * width - z, z)
            u = z - u_bound

        m = np.sin(xi * u)
        # escape bookkeeping (only meaningful when a saddle exists)
        if m_saddle is not None:
            entered |= (m > enter_gate)
            new_esc = entered & (m < m_saddle) & (esc < 0)
            if new_esc.any():
                esc[new_esc] = t

        income = MU * (0.5 + 0.5 * m)
        consumption = 30.0 + 10.0 * (w / 500.0)
        w = np.maximum(0.0, w + (income - consumption) * dt)
        is_low = w < COLLAPSE_WEALTH
        stk = np.where(is_low, stk + 1, 0)
        new = (stk >= collapse_streak) & (~col)
        if new.any():
            col |= new
            cstep[new] = t

    alive_at_late = ~(col & (cstep < late_step))     # uncollapsed at t=LATE_T
    return dict(collapsed=col, cstep=cstep, esc=esc, alive_at_late=alive_at_late,
                cap_rate=cap_events / (n * n_steps),
                fold_rate=fold_events / (n * n_steps))


# ----------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------
def late_hazard(res, n):
    """MLE constant hazard in the window (LATE_T, T_MAX] among at-risk-at-LATE_T.

    Returns (rate, n_late, n_at_risk, person_time).
    """
    late_step = int(round(LATE_T / DT))
    col, cstep = res["collapsed"], res["cstep"]
    at_risk = res["alive_at_late"]
    nar = int(at_risk.sum())
    late = col & (cstep > late_step)
    n_late = int(late.sum())
    # person-time after LATE_T for at-risk trajectories
    ct = np.where(col, cstep.astype(float) * DT, T_MAX)
    ct = np.minimum(ct, T_MAX)
    exposure = np.where(at_risk, np.maximum(ct - LATE_T, 0.0), 0.0).sum()
    rate = n_late / exposure if exposure > 0 else np.nan
    return rate, n_late, nar, exposure


def timing(res):
    col, cstep = res["collapsed"], res["cstep"]
    t = cstep[col].astype(float) * DT
    n_col = int(col.sum())
    n_early = int((t <= EARLY_T).sum())
    n_late = int((t > LATE_T).sum())
    med = float(np.median(t)) if n_col else np.nan
    p99 = float(np.percentile(t, 99)) if n_col else np.nan
    return dict(p=col.mean(), n_col=n_col, n_early=n_early, n_late=n_late,
                median_t=med, p99_t=p99)


def fmt(x, d=4):
    return "nan" if (isinstance(x, float) and np.isnan(x)) else f"{x:.{d}f}"


# ----------------------------------------------------------------------------
def main():
    os.makedirs(OUTDIR, exist_ok=True)
    t0 = time.time()
    log = []
    def P(s=""):
        print(s, flush=True); log.append(s)

    # ---- barriers per cell ----
    bar = {}
    P("# Late-channel mechanism test\n")
    P("## Fixed points and barriers (deterministic drift)\n")
    P("| h | J | m_high | m_saddle | dU (drift) | G (geom) | 2G/xi^2 |")
    P("|---|---|---|---|---|---|---|")
    for (h, J) in CELLS:
        m_high, m_s, dU, G = barriers(J, h)
        bar[(h, J)] = (m_high, m_s, dU, G)
        eff = (2.0 / XI**2) * G if not np.isnan(G) else np.nan
        P(f"| {h} | {J:g} | {fmt(m_high,4)} | {fmt(m_s,4) if m_s is not None else 'none'} "
          f"| {fmt(dU,4)} | {fmt(G,4)} | {fmt(eff,3)} |")
    P("")

    # ---- L1: IC dependence ----
    P("## L1 — initial-condition dependence (Lamperti, n=2000, same seeds)\n")
    P("| h | J | IC | P(coll) | n_col | n_early(t<=50) | n_late(t>250) | median t | p99 t | late rate r | n at-risk |")
    P("|---|---|---|---|---|---|---|---|---|---|---|")
    L1 = {}
    for i, (h, J) in enumerate(CELLS):
        seed = SEED_BASE + i * 1_000_003
        m_s = bar[(h, J)][1]
        for ic, m0 in (("a_m0=0", 0.0), ("b_m0=0.9", 0.9)):
            res = lamperti_run(h, J, N_TRAJ, seed, m0, m_s)
            tm = timing(res)
            r, n_late, nar, expo = late_hazard(res, N_TRAJ)
            L1[(h, J, ic)] = dict(res=res, tm=tm, rate=r, n_late=n_late,
                                  nar=nar, expo=expo)
            P(f"| {h} | {J:g} | {ic} | {tm['p']:.4f} | {tm['n_col']} "
              f"| {tm['n_early']} | {tm['n_late']} | {fmt(tm['median_t'],1)} "
              f"| {fmt(tm['p99_t'],1)} | {fmt(r,6)} | {nar} |")
        P("")

    # pooled early removal and per-cell / pooled rate ratios
    tot_early_a = sum(L1[(h, J, "a_m0=0")]["tm"]["n_early"] for (h, J) in CELLS)
    tot_early_b = sum(L1[(h, J, "b_m0=0.9")]["tm"]["n_early"] for (h, J) in CELLS)
    late_a = sum(L1[(h, J, "a_m0=0")]["n_late"] for (h, J) in CELLS)
    late_b = sum(L1[(h, J, "b_m0=0.9")]["n_late"] for (h, J) in CELLS)
    expo_a = sum(L1[(h, J, "a_m0=0")]["expo"] for (h, J) in CELLS)
    expo_b = sum(L1[(h, J, "b_m0=0.9")]["expo"] for (h, J) in CELLS)
    pooled_ra = late_a / expo_a if expo_a else np.nan
    pooled_rb = late_b / expo_b if expo_b else np.nan
    P("### L1 pooled\n")
    P(f"- early collapses: IC(a)={tot_early_a}, IC(b)={tot_early_b}, "
      f"reduction = {1 - tot_early_b / max(tot_early_a,1):.3f}")
    P(f"- pooled late rate: r(a)={fmt(pooled_ra,6)}, r(b)={fmt(pooled_rb,6)}, "
      f"ratio r(b)/r(a) = {fmt(pooled_rb/pooled_ra,3) if pooled_ra else 'nan'}")
    P("- per-cell late-rate ratio (cells with >=10 late in both arms):")
    ratios = []
    for (h, J) in CELLS:
        ra = L1[(h, J, "a_m0=0")]["rate"]; rb = L1[(h, J, "b_m0=0.9")]["rate"]
        na = L1[(h, J, "a_m0=0")]["n_late"]; nb = L1[(h, J, "b_m0=0.9")]["n_late"]
        if na >= 10 and nb >= 10 and ra > 0:
            ratios.append(rb / ra)
            P(f"    (h={h}, J={J:g}): r(b)/r(a) = {rb/ra:.3f} (n_late a={na}, b={nb})")
    P("")

    # L1 verdict
    early_removed = tot_early_b <= 0.20 * tot_early_a
    pooled_ratio = pooled_rb / pooled_ra if pooled_ra else np.nan
    rate_preserved = (not np.isnan(pooled_ratio)) and (0.5 <= pooled_ratio <= 2.0)
    rate_vanished = (not np.isnan(pooled_ratio)) and (pooled_ratio < 0.25)
    if early_removed and rate_preserved:
        L1_verdict = "SUPPORTED"
    elif rate_vanished:
        L1_verdict = "REFUTED"
    else:
        L1_verdict = "INCONCLUSIVE"
    P(f"**L1 verdict: {L1_verdict}** "
      f"(early_removed={early_removed}, pooled ratio={fmt(pooled_ratio,3)})\n")

    # ---- L2: waiting-time memorylessness (IC a) ----
    P("## L2 — waiting-time distribution of late collapses (IC a, m0=0)\n")
    P("| h | J | series | n | mean | CV | KS stat | KS p | consistent? |")
    P("|---|---|---|---|---|---|---|---|---|")
    l2_rows = []
    qq = {}
    for (h, J) in CELLS:
        res = L1[(h, J, "a_m0=0")]["res"]
        col, cstep, esc = res["collapsed"], res["cstep"], res["esc"]
        late_step = int(round(LATE_T / DT))
        # collapse-time waiting series (shifted by observed min)
        ctimes = cstep[col & (cstep > late_step)].astype(float) * DT
        # m-escape waiting series (cleaner Poisson)
        etimes = esc[(esc > late_step)].astype(float) * DT
        for label, arr in (("collapse", ctimes), ("escape", etimes)):
            n = len(arr)
            if n < 20:
                P(f"| {h} | {J:g} | {label} | {n} | - | - | - | - | (too few) |")
                continue
            w = arr - arr.min()                 # shift so support starts at 0
            mean = w.mean(); cv = w.std(ddof=1) / mean if mean > 0 else np.nan
            ks, p = stats.kstest(w, "expon", args=(0.0, mean))
            ok = (0.75 <= cv <= 1.25) and (p > 0.05)
            P(f"| {h} | {J:g} | {label} | {n} | {mean:.1f} | {cv:.3f} "
              f"| {ks:.3f} | {p:.3f} | {'yes' if ok else 'no'} |")
            l2_rows.append((h, J, label, n, mean, cv, ks, p, ok))
            if label == "escape":
                qq[(h, J)] = np.sort(w)
    P("")

    # ---- L3: barrier scaling regression (IC a late hazard) ----
    P("## L3 — barrier scaling of the late hazard (IC a)\n")
    reg = []
    for (h, J) in CELLS:
        r = L1[(h, J, "a_m0=0")]["rate"]; n_late = L1[(h, J, "a_m0=0")]["n_late"]
        m_high, m_s, dU, G = bar[(h, J)]
        if n_late >= 10 and r and r > 0 and not np.isnan(dU):
            reg.append((h, J, dU, G, r, np.log(r), n_late))
    P("Cells entering the regression (n_late>=10):")
    P("| h | J | dU | G | late rate r | ln r | n_late |")
    P("|---|---|---|---|---|---|---|")
    for (h, J, dU, G, r, lr, nl) in reg:
        P(f"| {h} | {J:g} | {dU:.4f} | {G:.4f} | {r:.6f} | {lr:.3f} | {nl} |")
    P("")
    def linfit(xs, ys):
        xs = np.asarray(xs); ys = np.asarray(ys)
        b, a = np.polyfit(xs, ys, 1)
        yhat = a + b * xs
        ss_res = ((ys - yhat) ** 2).sum()
        ss_tot = ((ys - ys.mean()) ** 2).sum()
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        return b, a, r2
    l3 = {}
    if len(reg) >= 3:
        dUs = [x[2] for x in reg]; Gs = [x[3] for x in reg]; lrs = [x[5] for x in reg]
        bU, aU, r2U = linfit(dUs, lrs)
        bG, aG, r2G = linfit(Gs, lrs)
        pred_slope = -2.0 / XI**2
        slope_ok = (-11.2 <= bG <= -4.8)
        l3 = dict(bU=bU, r2U=r2U, bG=bG, r2G=r2G, pred=pred_slope, slope_ok=slope_ok)
        P(f"- ln r vs dU (drift potential): slope={bU:.3f}, R^2={r2U:.3f}")
        P(f"- ln r vs G (geometric barrier): slope={bG:.3f}, R^2={r2G:.3f}; "
          f"Kramers predicts -2/xi^2 = {pred_slope:.2f}; within +-40%? {slope_ok}")
    else:
        P(f"- too few cells with >=10 late events ({len(reg)}) for a regression.")
    P("")

    # ---- noise sweep at one cell (h=2.4, J=15) ----
    P("## L3b — noise sweep at (h=2.4, J=15)\n")
    P("| xi | late rate r | n_late | ln r | 2G/xi^2 |")
    P("|---|---|---|---|---|")
    hs_cell = (2.4, 15.0); m_s_cell = bar[hs_cell][1]; G_cell = bar[hs_cell][3]
    seed_ns = SEED_BASE + CELLS.index(hs_cell) * 1_000_003
    ns_rows = []
    for xi in (0.35, 0.5, 0.65):
        res = lamperti_run(hs_cell[0], hs_cell[1], N_TRAJ, seed_ns, 0.0, m_s_cell, xi=xi)
        r, n_late, nar, expo = late_hazard(res, N_TRAJ)
        eff = (2.0 / xi**2) * G_cell
        lr = np.log(r) if (r and r > 0) else np.nan
        ns_rows.append((xi, r, n_late, lr, eff))
        P(f"| {xi} | {fmt(r,6)} | {n_late} | {fmt(lr,3)} | {eff:.3f} |")
    rates_ns = [x[1] for x in ns_rows]
    monotone = all(rates_ns[k] < rates_ns[k+1] for k in range(len(rates_ns)-1)
                   if not (np.isnan(rates_ns[k]) or np.isnan(rates_ns[k+1])))
    P(f"\n- rate monotonically increasing with xi? {monotone}\n")

    # ---- write CSVs ----
    with open(os.path.join(OUTDIR, "L1_ic_dependence.csv"), "w") as fh:
        fh.write("h,J,ic,p_collapse,n_col,n_early,n_late,median_t,p99_t,late_rate,n_at_risk,exposure\n")
        for (h, J) in CELLS:
            for ic in ("a_m0=0", "b_m0=0.9"):
                d = L1[(h, J, ic)]; tm = d["tm"]
                fh.write(f"{h},{J:g},{ic},{tm['p']:.6f},{tm['n_col']},{tm['n_early']},"
                         f"{tm['n_late']},{fmt(tm['median_t'],3)},{fmt(tm['p99_t'],3)},"
                         f"{fmt(d['rate'],8)},{d['nar']},{d['expo']:.3f}\n")
    with open(os.path.join(OUTDIR, "L2_waiting_times.csv"), "w") as fh:
        fh.write("h,J,series,n,mean,cv,ks_stat,ks_p,consistent\n")
        for row in l2_rows:
            fh.write(",".join(str(x) if not isinstance(x, float) else f"{x:.6g}"
                              for x in row) + "\n")
    with open(os.path.join(OUTDIR, "L3_barriers.csv"), "w") as fh:
        fh.write("h,J,m_high,m_saddle,dU,G,eff_2G_over_xi2,late_rate,ln_rate,n_late\n")
        for (h, J) in CELLS:
            m_high, m_s, dU, G = bar[(h, J)]
            r = L1[(h, J, "a_m0=0")]["rate"]; nl = L1[(h, J, "a_m0=0")]["n_late"]
            eff = (2.0 / XI**2) * G if not np.isnan(G) else np.nan
            lr = np.log(r) if (r and r > 0) else np.nan
            fh.write(f"{h},{J:g},{fmt(m_high,6)},{fmt(m_s,6) if m_s is not None else 'nan'},"
                     f"{fmt(dU,6)},{fmt(G,6)},{fmt(eff,6)},{fmt(r,8)},{fmt(lr,6)},{nl}\n")
    with open(os.path.join(OUTDIR, "L3b_noise_sweep.csv"), "w") as fh:
        fh.write("xi,late_rate,n_late,ln_rate,eff_2G_over_xi2\n")
        for (xi, r, nl, lr, eff) in ns_rows:
            fh.write(f"{xi},{fmt(r,8)},{nl},{fmt(lr,6)},{eff:.6f}\n")
    for (h, J), w in qq.items():
        q_theory = -np.log(1 - (np.arange(1, len(w) + 1) - 0.5) / len(w))
        with open(os.path.join(OUTDIR, f"L2_qq_escape_h{h}_J{J:g}.csv"), "w") as fh:
            fh.write("theoretical_exp_quantile,observed_shifted_waiting\n")
            for a, b in zip(q_theory * w.mean(), w):
                fh.write(f"{a:.6f},{b:.6f}\n")

    with open(os.path.join(OUTDIR, "SUMMARY.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(log) + "\n")
    P(f"Runtime: {time.time() - t0:.0f}s  -> {OUTDIR}")


if __name__ == "__main__":
    main()
