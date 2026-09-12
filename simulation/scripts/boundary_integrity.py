"""Boundary/clipping integrity audit for the minimal model.

Three parts, all at the base-integrator settings of minimal_model.py
(Euler-Maruyama, dt=0.05, 20,000 steps, T=2, xi=0.5, clip m to [-1+1e-6, 1-1e-6],
h = 2*(W/500), wealth floored at 0; collapse = W<10 for 200 consecutive steps;
types by |m| at collapse):

(i)  Clip-hit instrumentation: fraction of STEPS clipped and fraction of
     TRAJECTORIES ever clipped, split upper/lower, at the headline cell
     (J=5, mu=100, passive stabilizer) and at fixed-field cells
     (constant h in {2.4, 5.0}, J in {5, 10, 15, 20}).  n=1000 per cell.

(ii) Boundary-preserving comparison at the headline cell, n=2000, paired
     noise: clipped-EM in m vs a Lamperti-transformed integrator.
     NOTE ON THE TRANSFORM: with dm = f dt + xi*sqrt(1-m^2) dW,
       u = artanh(m):  du = [(f + xi^2 m)/(1-m^2)] dt + [xi/sqrt(1-m^2)] dW
     (Ito correction +xi^2 m/(1-m^2); the diffusion is xi*cosh(u), NOT
     constant xi -- the naïve final algebra step drops a factor).
     Because m=+/-1 is an ATTAINABLE boundary here (see part iii), the
     artanh coordinate reaches +/-inf in finite time and the artanh scheme
     explodes numerically (demonstrated empirically below).
     The correct constant-diffusion Lamperti transform for sigma=xi*sqrt(1-m^2)
     is u = arcsin(m)/xi:
       du = [(f + xi^2 m/2) / (xi*sqrt(1-m^2))] dt + dW,   m = sin(xi*u),
     which preserves m in [-1,1] with no clipping; boundary crossings in u are
     folded back by exact reflection (sin(pi-x)=sin(x)), i.e. the
     instantaneously-reflecting behaviour appropriate to a regular boundary.
     A trust-region cap |drift*dt| <= 0.5 tames the integrable drift
     singularity at the boundary (engagement counted and reported).

(iii) Analytic Feller boundary classification at m=+/-1 (written into
      SUMMARY.md with the calculation).

Outputs -> simulation/results/boundary_integrity/
  clip_stats.csv, lamperti_paired.csv, lamperti_summary.csv, SUMMARY.md
"""

from __future__ import annotations

import os
import time

import numpy as np

# ---------------------------------------------------------------------------
# Base-integrator constants (must match minimal_model.py)
# ---------------------------------------------------------------------------
T_TEMP = 2.0
XI = 0.5
DT = 0.05
N_STEPS = 20_000
EPS = 1e-6
INITIAL_M = 0.0
INITIAL_WEALTH = 100.0
COLLAPSE_WEALTH = 10.0
COLLAPSE_STREAK = 200
RIGIDITY_M = 0.9
FRAGMENTATION_M = 0.3

SEED_BASE = 0xC0DE
# Headline cell: J=5 (j_idx=9 in the 10-value grid), mu=100 (m_idx=8).
SEED_HEADLINE = SEED_BASE + 9 * 1_000_003 + 8 * 1009  # = 9,057,473
# Fixed-field cells: J in {5,10,15,20} -> j_idx 9,10,11,12 (natural grid
# extension), m_idx=8 (mu=100 slot).  Same seed for both h values at a given J.
FIXED_J = (5.0, 10.0, 15.0, 20.0)
FIXED_J_IDX = (9, 10, 11, 12)
FIXED_H = (2.4, 5.0)

HEAD_J = 5.0
HEAD_MULT = 100.0

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.normpath(os.path.join(BASE, "..", "results", "boundary_integrity"))

SQRT_DT = np.sqrt(DT)
U_BOUND = np.pi / (2.0 * XI)   # |u| bound for m = sin(xi*u); = pi for xi=0.5
DRIFT_CAP = 0.5                # trust-region cap on |drift*dt| in u-coordinates


# ---------------------------------------------------------------------------
# (i) Clipped EM with clip-hit instrumentation
# ---------------------------------------------------------------------------
def em_run(J, mult, n, seed, fixed_h=None, n_steps=N_STEPS):
    """Exact replica of minimal_model.run_cell dynamics + clip counting.

    fixed_h: if not None, h is held constant and the wealth/collapse block is
    skipped (m-dynamics are then decoupled from wealth).
    Returns dict of clip stats and (for passive cells) collapse arrays.
    """
    rng = np.random.default_rng(seed)
    m = np.full(n, INITIAL_M, dtype=np.float64)
    wealth = np.full(n, INITIAL_WEALTH, dtype=np.float64)

    up_steps = 0
    dn_steps = 0
    up_traj = np.zeros(n, dtype=bool)
    dn_traj = np.zeros(n, dtype=bool)

    collapsed = np.zeros(n, dtype=bool)
    collapse_step = np.full(n, -1, dtype=np.int64)
    collapse_absm = np.full(n, np.nan)
    streak = np.zeros(n, dtype=np.int64)

    for t in range(n_steps):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = fixed_h if fixed_h is not None else (wealth / 500.0) * 2.0
        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = SQRT_DT * rng.standard_normal(n)
        prop = m_c + f * DT + g * dW

        up = prop > 1 - EPS
        dn = prop < -1 + EPS
        up_steps += int(up.sum())
        dn_steps += int(dn.sum())
        up_traj |= up
        dn_traj |= dn

        m = np.clip(prop, -1 + EPS, 1 - EPS)

        if fixed_h is None:
            employment = 0.5 + 0.5 * m
            income = mult * employment
            consumption = 30.0 + 10.0 * (wealth / 500.0)
            wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

            is_low = wealth < COLLAPSE_WEALTH
            streak = np.where(is_low, streak + 1, 0)
            new = (streak >= COLLAPSE_STREAK) & (~collapsed)
            if new.any():
                collapsed |= new
                collapse_step[new] = t
                collapse_absm[new] = np.abs(m[new])

    return dict(
        n=n, n_steps=n_steps,
        frac_steps_upper=up_steps / (n * n_steps),
        frac_steps_lower=dn_steps / (n * n_steps),
        frac_traj_upper=up_traj.mean(),
        frac_traj_lower=dn_traj.mean(),
        collapsed=collapsed, collapse_step=collapse_step,
        collapse_absm=collapse_absm,
    )


# ---------------------------------------------------------------------------
# (ii) Boundary-preserving Lamperti (arcsin) integrator, reflecting fold
# ---------------------------------------------------------------------------
def lamperti_run(J, mult, n, seed, n_steps=N_STEPS):
    rng = np.random.default_rng(seed)
    u = np.full(n, np.arcsin(INITIAL_M) / XI)   # = 0
    wealth = np.full(n, INITIAL_WEALTH, dtype=np.float64)

    collapsed = np.zeros(n, dtype=bool)
    collapse_step = np.full(n, -1, dtype=np.int64)
    collapse_absm = np.full(n, np.nan)
    streak = np.zeros(n, dtype=np.int64)

    cap_events = 0
    fold_events = 0

    for t in range(n_steps):
        s = np.sin(XI * u)                 # = m, exactly in [-1, 1]
        c = np.cos(XI * u)                 # = sqrt(1-m^2) >= 0 on |u|<=U_BOUND
        h = (wealth / 500.0) * 2.0
        f = -s + np.tanh((J * s + h) / T_TEMP)
        drift_dt = (f + 0.5 * XI * XI * s) / (XI * np.maximum(c, 1e-12)) * DT
        over = np.abs(drift_dt) > DRIFT_CAP
        cap_events += int(over.sum())
        drift_dt = np.clip(drift_dt, -DRIFT_CAP, DRIFT_CAP)

        dW = SQRT_DT * rng.standard_normal(n)   # identical draws to em_run
        u = u + drift_dt + dW

        out = (u > U_BOUND) | (u < -U_BOUND)
        fold_events += int(out.sum())
        if out.any():
            # reflect u into [-U_BOUND, U_BOUND]; sin(xi*u) is invariant
            width = 2.0 * U_BOUND
            z = np.mod(u + U_BOUND, 2.0 * width)
            z = np.where(z > width, 2.0 * width - z, z)
            u = z - U_BOUND

        m = np.sin(XI * u)
        employment = 0.5 + 0.5 * m
        income = mult * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        is_low = wealth < COLLAPSE_WEALTH
        streak = np.where(is_low, streak + 1, 0)
        new = (streak >= COLLAPSE_STREAK) & (~collapsed)
        if new.any():
            collapsed |= new
            collapse_step[new] = t
            collapse_absm[new] = np.abs(m[new])

    return dict(
        collapsed=collapsed, collapse_step=collapse_step,
        collapse_absm=collapse_absm,
        cap_rate=cap_events / (n * n_steps),
        fold_rate=fold_events / (n * n_steps),
    )


# ---------------------------------------------------------------------------
# artanh-scheme explosion demonstration (small n)
# ---------------------------------------------------------------------------
def artanh_run(J, mult, n, seed, n_steps=N_STEPS, u_explode=50.0):
    """Correct-Ito artanh scheme: du = [(f+xi^2 m)/(1-m^2)]dt + [xi/sqrt(1-m^2)]dW.
    Trajectories are frozen (marked exploded) once |u| > u_explode
    (|m| within ~1e-43 of 1)."""
    rng = np.random.default_rng(seed)
    u = np.zeros(n)
    wealth = np.full(n, INITIAL_WEALTH, dtype=np.float64)
    exploded = np.zeros(n, dtype=bool)
    explode_step = np.full(n, -1, dtype=np.int64)

    for t in range(n_steps):
        m = np.tanh(u)
        sech2 = 1.0 - m * m                 # fine for |u| <= 50 (>= 1e-43... underflow-safe)
        sech2 = np.maximum(sech2, 1e-300)
        h = (wealth / 500.0) * 2.0
        f = -m + np.tanh((J * m + h) / T_TEMP)
        drift = (f + XI * XI * m) / sech2
        diff = XI / np.sqrt(sech2)
        dW = SQRT_DT * rng.standard_normal(n)
        u_new = u + drift * DT + diff * dW
        u = np.where(exploded, u, u_new)

        newly = (~exploded) & (np.abs(u) > u_explode)
        exploded |= newly
        explode_step[newly] = t
        u = np.clip(u, -u_explode, u_explode)

        employment = 0.5 + 0.5 * np.tanh(u)
        income = mult * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

    return dict(exploded=exploded, explode_step=explode_step)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def classify(absm, rig=RIGIDITY_M, frag=FRAGMENTATION_M):
    """0=none, 1=rigidity, 2=fragmentation, 3=mixed (NaN -> 0)."""
    out = np.zeros(absm.shape, dtype=np.int8)
    ok = np.isfinite(absm)
    out[ok & (absm > rig)] = 1
    out[ok & (absm < frag)] = 2
    out[ok & (out == 0)] = 3
    return out


def paired_ci(d, z=1.96):
    d = np.asarray(d, dtype=float)
    mean = d.mean()
    se = d.std(ddof=1) / np.sqrt(len(d))
    return mean, mean - z * se, mean + z * se


def feller_numbers(J, h):
    """Inward drift c and CIR/Bessel-comparison quantities at m=+1 and m=-1."""
    c_up = 1.0 - np.tanh((J + h) / T_TEMP)     # drift magnitude toward interior at m=+1
    c_dn = 1.0 + np.tanh((-J + h) / T_TEMP)    # at m=-1
    out = {}
    for name, cc in (("+1", c_up), ("-1", c_dn)):
        delta = 2.0 * cc / XI**2               # Feller dimension 4a/sigma_loc^2, sigma_loc^2=2 xi^2
        cls = "entrance (unattainable)" if cc >= XI**2 else "regular (attainable, reflecting; not absorbing)"
        out[name] = (cc, delta, cls)
    return out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    t0 = time.time()

    # ------------------------------------------------------------------ (i)
    print("E6(i): clip-hit instrumentation ...")
    clip_rows = []
    r = em_run(HEAD_J, HEAD_MULT, n=1000, seed=SEED_HEADLINE)
    clip_rows.append(("headline_passive", HEAD_J, "passive", "", HEAD_MULT, 1000,
                      r["frac_steps_upper"], r["frac_steps_lower"],
                      r["frac_traj_upper"], r["frac_traj_lower"],
                      r["collapsed"].mean()))
    print(f"  headline done ({time.time()-t0:.0f}s)")
    for J, j_idx in zip(FIXED_J, FIXED_J_IDX):
        seed = SEED_BASE + j_idx * 1_000_003 + 8 * 1009
        for hfix in FIXED_H:
            r = em_run(J, HEAD_MULT, n=1000, seed=seed, fixed_h=hfix)
            clip_rows.append((f"fixed_h{hfix}_J{J:g}", J, "fixed", hfix, HEAD_MULT, 1000,
                              r["frac_steps_upper"], r["frac_steps_lower"],
                              r["frac_traj_upper"], r["frac_traj_lower"], ""))
            print(f"  J={J:g} h={hfix} done ({time.time()-t0:.0f}s)")

    with open(os.path.join(OUTDIR, "clip_stats.csv"), "w", encoding="utf-8") as fh:
        fh.write("cell,J,h_mode,h_value,mult,n,frac_steps_clipped_upper,"
                 "frac_steps_clipped_lower,frac_traj_ever_upper,frac_traj_ever_lower,"
                 "p_collapse\n")
        for row in clip_rows:
            fh.write(",".join(str(x) for x in row) + "\n")

    # ------------------------------------------------------------------ (ii)
    print("E6(ii): paired clipped-EM vs Lamperti, n=2000 ...")
    N2 = 2000
    em = em_run(HEAD_J, HEAD_MULT, n=N2, seed=SEED_HEADLINE)
    lam = lamperti_run(HEAD_J, HEAD_MULT, n=N2, seed=SEED_HEADLINE)
    print(f"  runs done ({time.time()-t0:.0f}s)")

    em_type = classify(em["collapse_absm"])
    lam_type = classify(lam["collapse_absm"])

    with open(os.path.join(OUTDIR, "lamperti_paired.csv"), "w", encoding="utf-8") as fh:
        fh.write("traj,em_collapsed,em_step,em_type,lam_collapsed,lam_step,lam_type\n")
        for i in range(N2):
            fh.write(f"{i},{int(em['collapsed'][i])},{em['collapse_step'][i]},"
                     f"{em_type[i]},{int(lam['collapsed'][i])},"
                     f"{lam['collapse_step'][i]},{lam_type[i]}\n")

    def stats(res, types):
        col = res["collapsed"]
        p = col.mean()
        rig_sh = (types == 1)[col].mean() if col.any() else np.nan
        frag_sh = (types == 2)[col].mean() if col.any() else np.nan
        mix_sh = (types == 3)[col].mean() if col.any() else np.nan
        steps = res["collapse_step"][col]
        med = np.median(steps) if col.any() else np.nan
        p99 = np.percentile(steps, 99) if col.any() else np.nan
        return p, rig_sh, frag_sh, mix_sh, med, p99

    em_s = stats(em, em_type)
    lam_s = stats(lam, lam_type)

    d_col = em["collapsed"].astype(float) - lam["collapsed"].astype(float)
    dc_mean, dc_lo, dc_hi = paired_ci(d_col)
    d_rig = (em_type == 1).astype(float) - (lam_type == 1).astype(float)
    dr_mean, dr_lo, dr_hi = paired_ci(d_rig)
    both = em["collapsed"] & lam["collapsed"]
    d_step = (em["collapse_step"][both] - lam["collapse_step"][both]).astype(float)
    ds_mean, ds_lo, ds_hi = paired_ci(d_step) if both.sum() > 1 else (np.nan,) * 3
    n_em_only = int((em["collapsed"] & ~lam["collapsed"]).sum())
    n_lam_only = int((lam["collapsed"] & ~em["collapsed"]).sum())

    with open(os.path.join(OUTDIR, "lamperti_summary.csv"), "w", encoding="utf-8") as fh:
        fh.write("metric,clipped_EM,lamperti,paired_diff,ci_lo,ci_hi\n")
        fh.write(f"P_collapse,{em_s[0]:.5f},{lam_s[0]:.5f},{dc_mean:.5f},{dc_lo:.5f},{dc_hi:.5f}\n")
        fh.write(f"rigidity_share_of_collapsed,{em_s[1]:.5f},{lam_s[1]:.5f},,,\n")
        fh.write(f"fragmentation_share_of_collapsed,{em_s[2]:.5f},{lam_s[2]:.5f},,,\n")
        fh.write(f"mixed_share_of_collapsed,{em_s[3]:.5f},{lam_s[3]:.5f},,,\n")
        fh.write(f"P_rigidity_collapse_uncond,{(em_type==1).mean():.5f},{(lam_type==1).mean():.5f},"
                 f"{dr_mean:.5f},{dr_lo:.5f},{dr_hi:.5f}\n")
        fh.write(f"collapse_step_median,{em_s[4]:.1f},{lam_s[4]:.1f},,,\n")
        fh.write(f"collapse_step_p99,{em_s[5]:.1f},{lam_s[5]:.1f},,,\n")
        fh.write(f"collapse_step_paired_mean_diff_bothcollapsed,,,{ds_mean:.2f},{ds_lo:.2f},{ds_hi:.2f}\n")
        fh.write(f"n_discordant_EM_only,{n_em_only},,,,\n")
        fh.write(f"n_discordant_Lamperti_only,{n_lam_only},,,,\n")
        fh.write(f"lamperti_drift_cap_rate,,{lam['cap_rate']:.3e},,,\n")
        fh.write(f"lamperti_fold_rate,,{lam['fold_rate']:.3e},,,\n")

    # artanh explosion demo
    print("E6(ii): artanh-scheme explosion demo, n=200 ...")
    at = artanh_run(HEAD_J, HEAD_MULT, n=200, seed=SEED_HEADLINE)
    at_frac = at["exploded"].mean()
    at_med = (np.median(at["explode_step"][at["exploded"]])
              if at["exploded"].any() else np.nan)
    print(f"  artanh: {at_frac:.1%} exploded, median step {at_med}")

    # ------------------------------------------------------------------ (iii)
    fell = {h: feller_numbers(HEAD_J, h) for h in (0.4, 2.4, 5.0)}

    # ------------------------------------------------------------------ SUMMARY
    lines = []
    A = lines.append
    A("# Boundary/Clipping Integrity")
    A("")
    A(f"Base integrator: EM, dt={DT}, {N_STEPS} steps, T={T_TEMP}, xi={XI}, "
      f"clip to [-1+1e-6, 1-1e-6]; collapse = W<{COLLAPSE_WEALTH:g} for "
      f"{COLLAPSE_STREAK} steps; headline cell J={HEAD_J:g}, mu={HEAD_MULT:g} "
      f"(passive), seed {SEED_HEADLINE} (0xC0DE + 9*1000003 + 8*1009).")
    A("Fixed-field cells use seeds 0xC0DE + j_idx*1000003 + 8*1009 with "
      "j_idx = 9,10,11,12 for J = 5,10,15,20 (natural grid extension); the same "
      "seed is shared by both h values at a given J. With h fixed the m-SDE "
      "decouples from wealth, so mu is irrelevant to clip statistics.")
    A("")
    A("## (i) Clip-hit statistics (n=1000 per cell)")
    A("")
    A("A 'clip hit' = the proposed EM update m + f dt + g dW falls outside "
      "[-1+1e-6, 1-1e-6] and is projected back. Fractions of steps are over "
      "all n * 20000 trajectory-steps.")
    A("")
    A("| cell | J | h | steps clipped upper | steps clipped lower | traj ever upper | traj ever lower |")
    A("|---|---|---|---|---|---|---|")
    for row in clip_rows:
        hlabel = "passive W-coupled" if row[2] == "passive" else f"{row[3]:g} (fixed)"
        A(f"| {row[0]} | {row[1]:g} | {hlabel} | {row[6]:.4f} | {row[7]:.4f} "
          f"| {row[8]:.3f} | {row[9]:.3f} |")
    A("")
    A(f"Headline-cell P(collapse) at n=1000: {clip_rows[0][10]:.3f}.")
    A("")
    A("## (ii) Clipped-EM vs boundary-preserving Lamperti (headline cell, "
      "n=2000, paired noise)")
    A("")
    A("Transform check (Ito, dm = f dt + xi sqrt(1-m^2) dW):")
    A("")
    A("- u = artanh(m):  u' = 1/(1-m^2), u'' = 2m/(1-m^2)^2  =>")
    A("  du = [(f + xi^2 m)/(1-m^2)] dt + [xi/sqrt(1-m^2)] dW.")
    A("  The Ito correction is +xi^2 m/(1-m^2) as stated in the task, but the")
    A("  diffusion is xi*cosh(u), NOT constant xi (the naïve final step")
    A("  multiplied by a spurious sqrt(1-m^2)). Since m=+/-1 is attainable")
    A("  (part iii), u = artanh(m) reaches +/-inf in finite time and the scheme")
    A("  explodes: empirically, with n=200 paired trajectories, "
      f"{at_frac:.1%} exploded past |u|>50 (median explosion step "
      f"{at_med:.0f}). The artanh route is unusable here.")
    A("- The correct constant-diffusion (Lamperti) transform is u = arcsin(m)/xi:")
    A("  du = [(f + xi^2 m/2)/(xi sqrt(1-m^2))] dt + dW,  m = sin(xi u).")
    A("  m stays in [-1,1] identically (no clipping); u is folded back into")
    A("  [-pi/(2 xi), pi/(2 xi)] by exact reflection (sin is invariant), which")
    A("  realises the instantaneously-reflecting behaviour of the regular")
    A("  boundary. The integrable drift singularity at the boundary is tamed by")
    A(f"  a trust-region cap |drift dt| <= {DRIFT_CAP} (engaged on "
      f"{lam['cap_rate']:.2e} of trajectory-steps; folds on "
      f"{lam['fold_rate']:.2e}).")
    A("  Identical per-step dW draws as the EM run (same generator, same call")
    A("  sequence) give a paired comparison.")
    A("")
    A("| metric | clipped EM | Lamperti | paired diff [95% CI] |")
    A("|---|---|---|---|")
    A(f"| P(collapse) | {em_s[0]:.4f} | {lam_s[0]:.4f} | "
      f"{dc_mean:+.4f} [{dc_lo:+.4f}, {dc_hi:+.4f}] |")
    A(f"| rigidity share of collapsed | {em_s[1]:.4f} | {lam_s[1]:.4f} |  |")
    A(f"| P(rigidity collapse), uncond. | {(em_type==1).mean():.4f} | "
      f"{(lam_type==1).mean():.4f} | {dr_mean:+.4f} [{dr_lo:+.4f}, {dr_hi:+.4f}] |")
    A(f"| collapse step median | {em_s[4]:.0f} | {lam_s[4]:.0f} | "
      f"{ds_mean:+.1f} [{ds_lo:+.1f}, {ds_hi:+.1f}] (both-collapsed) |")
    A(f"| collapse step p99 | {em_s[5]:.0f} | {lam_s[5]:.0f} |  |")
    A("")
    A(f"Discordant pairs: EM-only collapses {n_em_only}, Lamperti-only "
      f"{n_lam_only} (of n=2000).")
    A("")
    A("## (iii) Feller boundary classification at m = +/-1")
    A("")
    A("sigma(m) = xi sqrt(1-m^2), so sigma^2 = xi^2 (1-m^2) vanishes LINEARLY")
    A("at the boundary: with y = 1 -/+ m (distance to the boundary),")
    A("sigma^2 ~ 2 xi^2 y. The inward drift magnitude at the boundary is")
    A("  c_up = 1 - tanh((J+h)/T)  at m=+1,   c_dn = 1 + tanh((-J+h)/T)  at m=-1,")
    A("both > 0 (drift points inward). Near the boundary the process is the")
    A("CIR-type diffusion dy = c dt + xi sqrt(2y) dW. Its scale density is")
    A("s(y) ~ exp(-int 2c/(2 xi^2 y) dy) = y^(-c/xi^2) and speed density")
    A("m(y) ~ 1/(sigma^2 s) ~ y^(c/xi^2 - 1). Feller test:")
    A("")
    A("- c >= xi^2  (Feller condition 2c >= 2 xi^2, i.e. Bessel-type dimension")
    A("  delta = 2c/xi^2 >= 2): int s diverges at 0 -> boundary UNATTAINABLE;")
    A("  int m converges -> ENTRANCE boundary.")
    A("- 0 < c < xi^2 (delta < 2): int s and int m both converge at 0 ->")
    A("  REGULAR boundary: attainable in finite time, NOT absorbing (positive")
    A("  inward drift, finite speed measure); the clip/reflection choice fixes")
    A("  the boundary behaviour as instantaneous reflection.")
    A("")
    A(f"With T={T_TEMP:g}, xi={XI:g} (xi^2={XI**2:g}), J={HEAD_J:g}:")
    A("")
    A("| h | boundary | inward drift c | delta = 2c/xi^2 | classification |")
    A("|---|---|---|---|---|")
    for h in (0.4, 2.4, 5.0):
        for b in ("+1", "-1"):
            cc, delta, cls = fell[h][b]
            A(f"| {h:g} | {b} | {cc:.3e} | {delta:.3e} | {cls} |")
    A("")
    A("Conclusion: for the headline/passive range (h between 0 and ~2.4 during")
    A("the pre-collapse transient) BOTH boundaries are regular-attainable")
    A("(delta << 2), which is exactly why the EM integrator registers clip hits;")
    A("m=-1 becomes an entrance (unattainable) boundary only for h >~ 3.05")
    A("(where 1 - tanh((J-h)/T) >= xi^2). No boundary is absorbing in any")
    A("studied regime, so clipping approximates the correct reflecting")
    A("behaviour; part (ii) quantifies the residual discretisation effect.")
    A("")
    A(f"Runtime: {time.time()-t0:.0f}s.  Files: clip_stats.csv, "
      "lamperti_paired.csv, lamperti_summary.csv.")

    with open(os.path.join(OUTDIR, "SUMMARY.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"Done in {time.time()-t0:.0f}s -> {OUTDIR}")


if __name__ == "__main__":
    main()
