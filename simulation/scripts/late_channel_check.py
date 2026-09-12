"""Integrity check: is the fixed-field LATE collapse channel
(collapse t > 100 at constant h in {2.4, 5.0}, J in {15, 20}, mu=100; see
simulation/results/branch_observable/e3_collapse_times.csv) an artifact of
the clipped Euler-Maruyama integrator, or does it survive under a
boundary-preserving scheme?

Method: paired-noise comparison, per cell, n=2000:
  * clipped EM in m — exact replica of branch_observable.run_cell dynamics
    with fixed h (clip m to [-1+EPS, 1-EPS] before and after the update),
    but WITH the wealth/collapse block (income mu*(1+m)/2, consumption
    30 + 10*(W/500), wealth floored at 0; collapse = W < 10 for 200
    consecutive steps);
  * boundary-preserving Lamperti — the CORRECT constant-diffusion arcsin
    transform already implemented in boundary_integrity.lamperti_run,
    reused verbatim (formulas and constants imported from that module):
      u = arcsin(m)/xi,
      du = [(f + xi^2 m / 2) / (xi sqrt(1-m^2))] dt + dW,   m = sin(xi u),
      exact reflecting fold at |u| = pi/(2 xi)  (sin(pi-x) = sin(x)),
      trust-region cap |drift*dt| <= 0.5 (engagements counted),
    here with CONSTANT h (no wealth feedback into h) and the same
    wealth/collapse block as above for collapse detection.
  * PAIRED: one dW ~ sqrt(dt) N(0,1) draw per (seed, step), fed to BOTH
    integrators, so per-trajectory differences are integrator-only.

Cells and seeds: (h, J) in {2.4, 5.0} x {15, 20}, mu=100, cell index i in
row-major order [(2.4,15), (2.4,20), (5.0,15), (5.0,20)],
seed_i = 0xB0_0001 + i * 1_000_003.
(A fixed hexadecimal seed base is used for reproducibility.)
NOTE: e3_collapse_times.csv used n=1000 and seed base 0xC0DE + ji*1_000_003
(ji = index in JS_B), so the EM arm here is a fresh-seed, larger-n
reproduction — agreement is expected only within binomial error.

Base integrator settings (identical to minimal_model / boundary_integrity):
dt=0.05, 20000 steps (t_max=1000), T=2, xi=0.5, m0=0, W0=100.

Outputs -> simulation/results/boundary_integrity/late_channel/
  late_channel_paired.csv   (one row per trajectory per cell)
  late_channel_summary.csv  (one row per cell x integrator + paired diffs)
Verdict text -> simulation/results/boundary_integrity/late_channel/SUMMARY_late_channel.md
"""

from __future__ import annotations

import os
import time

import numpy as np

# Reuse the boundary_integrity implementation constants verbatim.
from boundary_integrity import (
    T_TEMP, XI, DT, N_STEPS, EPS,
    INITIAL_M, INITIAL_WEALTH, COLLAPSE_WEALTH, COLLAPSE_STREAK,
    SQRT_DT, U_BOUND, DRIFT_CAP,
)

MU = 100.0
N_TRAJ = 2000
CELLS = [(2.4, 15.0), (2.4, 20.0), (5.0, 15.0), (5.0, 20.0)]  # (h, J)
SEED_BASE = 0xB0_0001            # = 11534337
LATE_T = 100.0                   # "late" collapse threshold in physical time

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.normpath(os.path.join(
    BASE, "..", "results", "boundary_integrity", "late_channel"))
SUMDIR = os.path.normpath(os.path.join(
    BASE, "..", "results", "boundary_integrity", "late_channel"))

# e3_collapse_times.csv reference values (n=1000, seed base 0xC0DE) for the
# sanity check of the EM arm.
E3_REF = {  # (h, J): (n, p_collapse, frac_after_t100, n_collapsed)
    (2.4, 15.0): (1000, 0.178, 0.2808988764044944, 178),
    (2.4, 20.0): (1000, 0.269, 0.1970260223048327, 269),
    (5.0, 15.0): (1000, 0.035, 0.22857142857142856, 35),
    (5.0, 20.0): (1000, 0.084, 0.15476190476190477, 84),
}


def run_cell_paired(h, J, n, seed, n_steps=N_STEPS, dt=DT,
                    collapse_streak=COLLAPSE_STREAK):
    """Advance clipped-EM and Lamperti trajectories on the SAME dW stream.

    EM block replicates branch_observable.run_cell (fixed h) exactly;
    Lamperti block replicates boundary_integrity.lamperti_run exactly,
    except h is constant instead of 2*(W/500).  Both arms carry their own
    wealth/collapse ledger (identical update rule).

    dt / collapse_streak default to the module constants (dt=0.05, 200
    steps = 10 physical time units below W<10).  Callers changing dt MUST
    rescale collapse_streak by 1/dt to preserve the 10-unit streak.
    """
    rng = np.random.default_rng(seed)
    sqrt_dt = np.sqrt(dt)

    # --- EM state ---
    m_em = np.full(n, INITIAL_M, dtype=np.float64)
    w_em = np.full(n, INITIAL_WEALTH, dtype=np.float64)
    col_em = np.zeros(n, dtype=bool)
    step_em = np.full(n, -1, dtype=np.int64)
    stk_em = np.zeros(n, dtype=np.int64)
    clip_steps = 0

    # --- Lamperti state ---
    u = np.full(n, np.arcsin(INITIAL_M) / XI)   # = 0
    w_lm = np.full(n, INITIAL_WEALTH, dtype=np.float64)
    col_lm = np.zeros(n, dtype=bool)
    step_lm = np.full(n, -1, dtype=np.int64)
    stk_lm = np.zeros(n, dtype=np.int64)
    cap_events = 0
    fold_events = 0

    for t in range(n_steps):
        dW = sqrt_dt * rng.standard_normal(n)   # ONE draw, fed to both arms

        # ---------------- clipped EM (branch_observable.run_cell, fixed h)
        m_c = np.clip(m_em, -1 + EPS, 1 - EPS)
        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        prop = m_c + f * dt + g * dW
        clip_steps += int(((prop > 1 - EPS) | (prop < -1 + EPS)).sum())
        m_em = np.clip(prop, -1 + EPS, 1 - EPS)

        income = MU * (0.5 + 0.5 * m_em)
        consumption = 30.0 + 10.0 * (w_em / 500.0)
        w_em = np.maximum(0.0, w_em + (income - consumption) * dt)
        is_low = w_em < COLLAPSE_WEALTH
        stk_em = np.where(is_low, stk_em + 1, 0)
        new = (stk_em >= collapse_streak) & (~col_em)
        if new.any():
            col_em |= new
            step_em[new] = t

        # ---------------- Lamperti (boundary_integrity.lamperti_run, fixed h)
        s = np.sin(XI * u)                 # = m, exactly in [-1, 1]
        c = np.cos(XI * u)                 # = sqrt(1-m^2) >= 0 on |u|<=U_BOUND
        f = -s + np.tanh((J * s + h) / T_TEMP)
        drift_dt = (f + 0.5 * XI * XI * s) / (XI * np.maximum(c, 1e-12)) * dt
        over = np.abs(drift_dt) > DRIFT_CAP
        cap_events += int(over.sum())
        drift_dt = np.clip(drift_dt, -DRIFT_CAP, DRIFT_CAP)
        u = u + drift_dt + dW              # SAME dW as the EM arm

        out = (u > U_BOUND) | (u < -U_BOUND)
        fold_events += int(out.sum())
        if out.any():
            # reflect u into [-U_BOUND, U_BOUND]; sin(xi*u) is invariant
            width = 2.0 * U_BOUND
            z = np.mod(u + U_BOUND, 2.0 * width)
            z = np.where(z > width, 2.0 * width - z, z)
            u = z - U_BOUND

        m_lm = np.sin(XI * u)
        income = MU * (0.5 + 0.5 * m_lm)
        consumption = 30.0 + 10.0 * (w_lm / 500.0)
        w_lm = np.maximum(0.0, w_lm + (income - consumption) * dt)
        is_low = w_lm < COLLAPSE_WEALTH
        stk_lm = np.where(is_low, stk_lm + 1, 0)
        new = (stk_lm >= collapse_streak) & (~col_lm)
        if new.any():
            col_lm |= new
            step_lm[new] = t

    return dict(
        em=dict(collapsed=col_em, collapse_step=step_em,
                clip_step_rate=clip_steps / (n * n_steps)),
        lam=dict(collapsed=col_lm, collapse_step=step_lm,
                 cap_rate=cap_events / (n * n_steps),
                 fold_rate=fold_events / (n * n_steps)),
    )


def paired_ci(d, z=1.96):
    d = np.asarray(d, dtype=float)
    mean = d.mean()
    se = d.std(ddof=1) / np.sqrt(len(d))
    return mean, mean - z * se, mean + z * se


def timing_stats(steps):
    """median/p90/p99/max of physical collapse times; NaNs if empty."""
    if len(steps) == 0:
        return (np.nan,) * 4
    t = np.asarray(steps, dtype=float) * DT
    return (float(np.median(t)), float(np.percentile(t, 90)),
            float(np.percentile(t, 99)), float(t.max()))


def arm_stats(res):
    col = res["collapsed"]
    steps = res["collapse_step"][col]
    t = steps.astype(float) * DT
    n_col = int(col.sum())
    n_late = int((t > LATE_T).sum())
    late_frac = n_late / n_col if n_col else np.nan
    med, p90, p99, tmax = timing_stats(steps)
    return dict(p=col.mean(), n_collapsed=n_col, n_late=n_late,
                late_frac=late_frac, median_t=med, p90_t=p90,
                p99_t=p99, max_t=tmax)


def two_prop_z(k1, n1, k2, n2):
    """Pooled two-proportion z (for the e3 sanity check)."""
    p1, p2 = k1 / n1, k2 / n2
    pp = (k1 + k2) / (n1 + n2)
    se = np.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
    return (p1 - p2) / se if se > 0 else np.nan


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(SUMDIR, exist_ok=True)
    t0 = time.time()

    paired_path = os.path.join(OUTDIR, "late_channel_paired.csv")
    summary_path = os.path.join(OUTDIR, "late_channel_summary.csv")

    cell_results = []
    with open(paired_path, "w", encoding="utf-8") as fh:
        fh.write("cell,h,J,mu,seed,traj,em_collapsed,em_step,em_t,"
                 "lam_collapsed,lam_step,lam_t\n")
        for i, (h, J) in enumerate(CELLS):
            seed = SEED_BASE + i * 1_000_003
            print(f"cell {i}: h={h} J={J:g} seed={seed} ...", flush=True)
            r = run_cell_paired(h, J, N_TRAJ, seed)
            cell_results.append((h, J, seed, r))
            em, lm = r["em"], r["lam"]
            for k in range(N_TRAJ):
                es, ls = em["collapse_step"][k], lm["collapse_step"][k]
                fh.write(f"h{h}_J{J:g},{h},{J:g},{MU:g},{seed},{k},"
                         f"{int(em['collapsed'][k])},{es},"
                         f"{es * DT if es >= 0 else ''},"
                         f"{int(lm['collapsed'][k])},{ls},"
                         f"{ls * DT if ls >= 0 else ''}\n")
            print(f"  done ({time.time() - t0:.0f}s)", flush=True)

    rows = []
    verdicts = {}
    for (h, J, seed, r) in cell_results:
        em, lm = r["em"], r["lam"]
        se_, le_ = arm_stats(em), arm_stats(lm)

        d_col = em["collapsed"].astype(float) - lm["collapsed"].astype(float)
        dc, dc_lo, dc_hi = paired_ci(d_col)
        em_late = em["collapsed"] & (em["collapse_step"] * DT > LATE_T)
        lm_late = lm["collapsed"] & (lm["collapse_step"] * DT > LATE_T)
        d_late = em_late.astype(float) - lm_late.astype(float)
        dl, dl_lo, dl_hi = paired_ci(d_late)

        # e3 sanity (EM arm vs published n=1000 run, fresh seeds)
        ref_n, ref_p, ref_lf, ref_k = E3_REF[(h, J)]
        z_p = two_prop_z(se_["n_collapsed"], N_TRAJ, ref_k, ref_n)
        z_late = two_prop_z(se_["n_late"], se_["n_collapsed"],
                            round(ref_lf * ref_k), ref_k) \
            if se_["n_collapsed"] else np.nan

        # verdict per cell on the late channel
        if le_["n_late"] == 0:
            v = "evaporates"
        elif le_["n_late"] < 0.5 * se_["n_late"]:
            v = "shrinks materially"
        else:
            v = "survives"
        verdicts[(h, J)] = (v, se_, le_, (dc, dc_lo, dc_hi),
                            (dl, dl_lo, dl_hi), z_p, z_late)

        for arm, st, extra in (
            ("clipped_EM", se_, dict(clip_step_rate=em["clip_step_rate"])),
            ("lamperti", le_, dict(cap_rate=lm["cap_rate"],
                                   fold_rate=lm["fold_rate"])),
        ):
            rows.append(dict(
                cell=f"h{h}_J{J:g}", h=h, J=J, mu=MU, n=N_TRAJ, seed=seed,
                integrator=arm, p_collapse=st["p"],
                n_collapsed=st["n_collapsed"], n_late_t100=st["n_late"],
                late_frac=st["late_frac"], median_t=st["median_t"],
                p90_t=st["p90_t"], p99_t=st["p99_t"], max_t=st["max_t"],
                paired_dP_EMminusLam=dc, paired_dP_ci_lo=dc_lo,
                paired_dP_ci_hi=dc_hi, paired_dPlate_EMminusLam=dl,
                paired_dPlate_ci_lo=dl_lo, paired_dPlate_ci_hi=dl_hi,
                e3_sanity_z_pcollapse=z_p, e3_sanity_z_latefrac=z_late,
                verdict=v, **extra))

    # column union across arms (EM has clip_step_rate; Lamperti cap/fold)
    cols = list(dict.fromkeys(list(rows[0].keys()) + list(rows[1].keys())))
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write(",".join(cols) + "\n")
        for row in rows:
            fh.write(",".join(
                "" if (k not in row or row[k] is None or
                       (isinstance(row[k], float) and np.isnan(row[k])))
                else (f"{row[k]:.6g}" if isinstance(row[k], float)
                      else str(row[k]))
                for k in cols) + "\n")

    # ------------------------------------------------------------- SUMMARY.md
    L = []
    A = L.append
    A("# Late collapse channel under a boundary-preserving integrator")
    A("")
    A("**Question.** e3_collapse_times.csv (clipped EM, n=1000) shows that at")
    A("constant field h in {2.4, 5.0}, J in {15, 20}, mu=100, 8-28% of wealth")
    A("collapses occur after t=100 (p99 up to ~985). Is this late channel a")
    A("clipped-EM boundary artifact?")
    A("")
    A("**Method.** Paired-noise comparison, n=2000 per cell: clipped EM (exact")
    A("replica of branch_observable.run_cell with fixed h + wealth/collapse")
    A("ledger) vs the constant-diffusion Lamperti scheme of")
    A("boundary_integrity.lamperti_run (u = arcsin(m)/xi, exact reflecting")
    A(f"fold at |u| = pi/(2 xi), trust cap |drift dt| <= {DRIFT_CAP}), with")
    A("constant h and the identical wealth/collapse ledger. One dW draw per")
    A("(seed, step) feeds BOTH integrators. dt=0.05, 20000 steps, T=2,")
    A("xi=0.5, m0=0, W0=100; collapse = W<10 for 200 consecutive steps;")
    A("late = collapse time t > 100. Seeds 0xB00001 + i*1000003, cell order")
    A("(2.4,15), (2.4,20), (5.0,15), (5.0,20). Script:")
    A("simulation/scripts/late_channel_check.py; data:")
    A("simulation/results/boundary_integrity/late_channel/.")
    A("")
    A("## Per-cell results (n=2000, paired)")
    A("")
    A("| h | J | arm | P(collapse) | n coll | n late (t>100) | late frac | median t | p90 | p99 | max |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for (h, J, seed, r) in cell_results:
        v, se_, le_, dci, dlci, z_p, z_late = verdicts[(h, J)]
        for arm, st in (("EM", se_), ("Lamperti", le_)):
            def fmt(x, d=1):
                return "-" if (isinstance(x, float) and np.isnan(x)) \
                    else f"{x:.{d}f}"
            A(f"| {h} | {J:g} | {arm} | {st['p']:.4f} | {st['n_collapsed']} "
              f"| {st['n_late']} | {fmt(st['late_frac'], 3)} "
              f"| {fmt(st['median_t'])} | {fmt(st['p90_t'])} "
              f"| {fmt(st['p99_t'])} | {fmt(st['max_t'])} |")
    A("")
    A("## Paired differences (EM minus Lamperti, 95% CI)")
    A("")
    A("| h | J | dP(collapse) | dP(late collapse) | verdict on late channel |")
    A("|---|---|---|---|---|")
    for (h, J, seed, r) in cell_results:
        v, se_, le_, (dc, dcl, dch), (dl, dll, dlh), z_p, z_late = \
            verdicts[(h, J)]
        A(f"| {h} | {J:g} | {dc:+.4f} [{dcl:+.4f}, {dch:+.4f}] "
          f"| {dl:+.4f} [{dll:+.4f}, {dlh:+.4f}] | **{v}** |")
    A("")
    A("## Sanity vs e3_collapse_times.csv (EM arm)")
    A("")
    A("Different n (2000 vs 1000) and different seed base (0xB00001 vs")
    A("0xC0DE-derived), so only binomial-level agreement is expected.")
    A("Two-proportion z-scores (this run vs e3):")
    A("")
    A("| h | J | e3 P | this EM P | z(P) | e3 late frac | this EM late frac | z(late) |")
    A("|---|---|---|---|---|---|---|---|")
    for (h, J, seed, r) in cell_results:
        v, se_, le_, dci, dlci, z_p, z_late = verdicts[(h, J)]
        ref_n, ref_p, ref_lf, ref_k = E3_REF[(h, J)]
        A(f"| {h} | {J:g} | {ref_p:.3f} | {se_['p']:.3f} | {z_p:+.2f} "
          f"| {ref_lf:.3f} | {se_['late_frac']:.3f} | {z_late:+.2f} |")
    A("")
    A("## Integrator diagnostics")
    A("")
    A("| h | J | EM clip-step rate | Lamperti cap rate | Lamperti fold rate |")
    A("|---|---|---|---|---|")
    for (h, J, seed, r) in cell_results:
        em, lm = r["em"], r["lam"]
        A(f"| {h} | {J:g} | {em['clip_step_rate']:.3e} "
          f"| {lm['cap_rate']:.3e} | {lm['fold_rate']:.3e} |")
    A("")
    A("## Verdict")
    A("")
    for (h, J, seed, r) in cell_results:
        v, se_, le_, dci, dlci, z_p, z_late = verdicts[(h, J)]
        A(f"- (h={h}, J={J:g}): late channel **{v}** under the")
        A(f"  boundary-preserving scheme (EM {se_['n_late']} vs Lamperti "
          f"{le_['n_late']} late collapses of n=2000).")
    A("")
    A(f"Runtime: {time.time() - t0:.0f}s.")

    with open(os.path.join(SUMDIR, "SUMMARY_late_channel.md"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")

    print(f"Done in {time.time() - t0:.0f}s")
    print(f"  {paired_path}")
    print(f"  {summary_path}")
    print(f"  {os.path.join(SUMDIR, 'SUMMARY_late_channel.md')}")


if __name__ == "__main__":
    main()
