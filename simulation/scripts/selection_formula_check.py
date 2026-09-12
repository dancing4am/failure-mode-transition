# -*- coding: utf-8 -*-
"""
selection_formula_check.py — selection-formula predictions vs. the measured amplitude surface (SI S2.7).

The S2.7 closed form (confirmed against manuscript/supplementary_information.md
S2.7 and simulation/scripts/phase1a_suzuki_deepen.py):

    P_pred = Phi( -(c/xi) * sqrt(2/lambda) ),
    c      = tanh(h_eff / T),
    lambda = -1 + (J/T) * sech^2(h_eff / T),
    T = 2, xi = 0.5 (committed constants), Phi = standard normal CDF.

Prefactor convention check (S2.7 / phase1a): the transient field is the
W0-scale field h0 = 2*(W0/500) with W0 = INITIAL_WEALTH = 100, i.e. h0 = 0.4
EXACTLY (phase1a: H0 = (A.INITIAL_WEALTH/500)*2).  The "W0-scale field during
the decision window" and "h_eff = 0.4 exactly" therefore coincide; a single
passive prediction column is reported.

E4 — predicted vs measured P(collapse) per stabilizer form (mu = 100):
  (a) passive:   h_eff = 0.4;             measured: results/minimal_model
                 (J=2.5..5, mult=100) + results/extended_sweep (J=6..10, mu=100).
  (b) fixed h:   h_eff = 2.4 / 5.0;       measured: results/active_stabilizer/
                 scaling_control.csv, J in {3,5,7,10,15,20}.  Where lambda <= 0
                 the formula predicts NO instability (P_pred = 0 identically);
                 rows are flagged and the unexplained measured failure reported.
  (c) coupling:  h_eff = 0.4*(1+J)        (h = base*(1+alpha*J), alpha=1, base
                 evaluated at W0 -> 0.4; delivered-field-at-W0 convention);
      quadratic: h_eff = 0.4*(1+(J/5)^2)  (same convention, alpha=1);
                 measured: active_coupling_summary.csv / active_quadratic_summary.csv
                 at mu=100, J in {4, 4.5, 5}.
  (d) stress:    h_eff = realized window-mean field from
                 results/magnitude_matched/ if present; else the stated
                 approximation h_eff = 0.4 + alpha*E[max(0,-m)]*J (alpha=2)
                 with E[max(0,-m)] the mean over the decision window from a
                 quick n=500 rerun of the stress condition (window = first
                 280 steps = 1.4% of horizon, the S2.7 median-collapse window).

E5 — theory-derived h0*(J).  Holding P fixed at target P:
    z_P = Phi^{-1}(P)  =>  c = tanh(h_eff/T) = |z_P| * xi * sqrt(lambda/2)
    =>  h_eff = T * artanh( |z_P| * xi * sqrt(lambda(h_eff, J)/2) )   (implicit)
Solved by brentq on [0, h_lambda0) where h_lambda0 = T*arccosh(sqrt(J/T))
(the lambda=0 boundary), per (xi in {0.35,0.5,0.65}) x (P in {0.10,0.20,0.30}),
J dense in [2.5, 5].

Mapping to the wealth-coupling coefficient h0 (h = h0 * W/500 in h0_scaling.py):
  PRIMARY (window convention, stated assumption): selection happens in the
    initial transient at W = W0 = 100, so h_eff = h0*(W0/500) = 0.2*h0
    =>  h0*_pred = h_eff*/0.2.
  ALTERNATIVE (equilibrium-scale mapping): on the good branch W -> W_eq = 3500
    at mu=100 (income 100 = consumption 30 + W/50), so h_eff = h0*(W_eq/500)
    = 7*h0  =>  h0*_pred = h_eff*/7.  Both are reported; fit quality decides.

Onset locus: lambda = 0 at h_eff  =>  J_onset(h_eff) = T*cosh^2(h_eff/T);
h0* -> 0 as J -> T*cosh^2(0) = T = 2.  Compared against the fitted J_c range
2.16-2.44 (results/h0_scaling/jc_free_fits.csv, unflagged rows).

Outputs (results/selection_formula/):
    e4_predicted_vs_measured.csv
    e5_h0_theory_grid.csv        (dense J grid, both mappings — figure-ready)
    e5_h0_pred_vs_measured.csv   (per measured cell)
    e5_onset_locus.csv
    SUMMARY.md

Run: py -3.13 selection_formula_check.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import active_stabilizer as A  # committed constants + simulator (reused)

RES = HERE.parents[0] / "results"
OUT = RES / "selection_formula"
OUT.mkdir(parents=True, exist_ok=True)

T, XI = A.T_TEMP, A.XI                      # 2.0, 0.5
H0 = (A.INITIAL_WEALTH / 500.0) * 2.0       # 0.4 exactly (S2.7 convention)
W0 = A.INITIAL_WEALTH                       # 100
W_EQ = 3500.0                               # good-branch wealth at mu=100
DECISION_WINDOW_STEPS = 280                 # 1.4% of 20000 (S2.7 median collapse)


def lam(J, h):
    return -1.0 + (J / T) / np.cosh(h / T) ** 2


def p_pred(J, h, xi=XI):
    """S2.7 selection formula; returns (P, lambda). P = 0.0 flagged upstream
    when lambda <= 0 (origin stable -> no instability predicted)."""
    l = lam(J, h)
    if l <= 0:
        return 0.0, l
    c = np.tanh(h / T)
    return float(norm.cdf(-(c / xi) * np.sqrt(2.0 / l))), l


# =========================================================================
# E4 — predicted vs measured table
# =========================================================================
rows = []


def add_row(case, J, h_eff, p_meas, n_seeds, note=""):
    P, l = p_pred(J, h_eff)
    flagged = l <= 0
    ratio = (P / p_meas) if p_meas > 0 else np.nan
    rows.append({
        "case": case, "J": J, "h_eff": h_eff, "lambda": l,
        "lambda_le_0": flagged, "P_pred": P, "P_meas": p_meas,
        "ratio_pred_over_meas": ratio, "n_seeds": n_seeds, "note": note,
    })


# ---- (a) passive baseline, mu = 100 -------------------------------------
mm = pd.read_csv(RES / "minimal_model" / "cell_summary.csv")
mm = mm[(mm["mult"] == 100) & (mm["J"] >= 2.5)].sort_values("J")
for _, r in mm.iterrows():
    add_row("passive", r.J, H0, r.p_collapse, 100,
            "h_eff = 0.4 = 2*(W0/500), W0=100; identical to S2.7 h0")
es = pd.read_csv(RES / "extended_sweep" / "extended_sweep_summary.csv")
es = es[es["mu"] == 100].sort_values("J")
for _, r in es.iterrows():
    add_row("passive", r.J, H0, r.p_collapse, 100, "extended_sweep")

# ---- (b) fixed h = 2.4 and 5.0 ------------------------------------------
sc = pd.read_csv(RES / "active_stabilizer" / "scaling_control.csv")
for cond, h_fix in [("fixed2.4", 2.4), ("fixed5.0", 5.0)]:
    sub = sc[sc["condition"] == cond].sort_values("J")
    for _, r in sub.iterrows():
        add_row(cond, r.J, h_fix, r.p_collapse, 100)

# ---- (c) coupling / quadratic schedules at mu=100, J in {4,4.5,5} -------
HIGH_J = [4.0, 4.5, 5.0]
cp = pd.read_csv(RES / "active_stabilizer" / "active_coupling_summary.csv")
cp = cp[(cp["mu"] == 100) & cp["J"].isin(HIGH_J)].sort_values("J")
for _, r in cp.iterrows():
    h_eff = H0 * (1.0 + r.J)   # delivered field with base evaluated at W0
    add_row("coupling", r.J, h_eff, r.p_collapse, 100,
            "h_eff = 0.4*(1+J), base at W0 (delivered-field convention)")
qd = pd.read_csv(RES / "active_stabilizer" / "active_quadratic_summary.csv")
qd = qd[(qd["mu"] == 100) & qd["J"].isin(HIGH_J)].sort_values("J")
for _, r in qd.iterrows():
    h_eff = H0 * (1.0 + (r.J / 5.0) ** 2)
    add_row("quadratic", r.J, h_eff, r.p_collapse, 100,
            "h_eff = 0.4*(1+(J/5)^2), base at W0")

# ---- (d) stress form ----------------------------------------------------
# PRIMARY: realized window-mean field from results/magnitude_matched (E1 run,
# n=2000, window = t <= 50 i.e. first 1000 steps; that mean ~5.1-5.2 is
# dominated by the passive wealth ramp 0.4 -> 14, not by the stress term,
# which activates in <1.1% of seed-steps — see magnitude_matched/SUMMARY.md).
# SECONDARY (labeled): the W0-base approximation h_eff = 0.4 + 2*E[max(0,-m)]*J
# with E from a quick n=500 rerun over the S2.7 decision window (280 steps).
mmatch = RES / "magnitude_matched"
have_mm = mmatch.exists() and (mmatch / "stress_field_stats.csv").exists()
if have_mm:
    print("[E4d] magnitude_matched present -> h_eff = realized window-mean field")
    sf = pd.read_csv(mmatch / "stress_field_stats.csv")
    cc = pd.read_csv(mmatch / "collapse_comparison.csv")
    for _, r in sf.sort_values("J").iterrows():
        p_meas = float(cc.loc[cc.J == r.J, "p_stress"].iloc[0])
        add_row("stress_realized", float(r.J), float(r.window_mean_h), p_meas,
                2000, "h_eff = realized window-mean field (magnitude_matched, "
                      "window t<=50; ramp-dominated)")
else:
    print("[E4d] results/magnitude_matched ABSENT -> only the stated "
          "approximation h_eff = 0.4 + alpha*E[max(0,-m)]*J is reported")

st = pd.read_csv(RES / "active_stabilizer" / "active_stress_summary.csv")
st = st[(st["mu"] == 100) & st["J"].isin(HIGH_J)].sort_values("J")

ALPHA_STRESS = 2.0
for _, r in st.iterrows():
    # quick n=500 rerun of the stress condition; record window-mean max(0,-m)
    J = float(r.J)
    rng_seed = A.SEED_BASE + 777_000 + int(J * 100)
    n, nw = 500, DECISION_WINDOW_STEPS
    rng = np.random.default_rng(rng_seed)
    m = np.full(n, A.INITIAL_M)
    wealth = np.full(n, A.INITIAL_WEALTH)
    acc = 0.0
    sqrt_dt = np.sqrt(A.DT)
    for t in range(nw):
        m_c = np.clip(m, -1 + A.EPS, 1 - A.EPS)
        h = A.h_active_stress(wealth, m_c, J, ALPHA_STRESS)
        f = -m_c + np.tanh((J * m_c + h) / A.T_TEMP)
        g = A.XI * np.sqrt(1.0 - m_c * m_c)
        m = np.clip(m_c + f * A.DT + g * sqrt_dt * rng.standard_normal(n),
                    -1 + A.EPS, 1 - A.EPS)
        employment = 0.5 + 0.5 * m
        wealth = np.maximum(0.0, wealth + (100.0 * employment
                                           - (30.0 + 10.0 * wealth / 500.0)) * A.DT)
        acc += np.maximum(0.0, -m).mean()
    E_neg = acc / nw
    h_eff = H0 + ALPHA_STRESS * E_neg * J
    add_row("stress_W0approx", J, h_eff, r.p_collapse, 100,
            f"SECONDARY convention: h_eff = 0.4 + 2*E[max(0,-m)]*J, "
            f"E={E_neg:.4f} (n=500 rerun, first {nw} steps; ignores wealth ramp)")

e4 = pd.DataFrame(rows)
e4.to_csv(OUT / "e4_predicted_vs_measured.csv", index=False)

print("\n=== E4 predicted vs measured ===")
print(e4[["case", "J", "h_eff", "lambda", "P_pred", "P_meas",
          "ratio_pred_over_meas"]].to_string(index=False,
          float_format=lambda x: f"{x:.4g}"))

# passive-block agreement stats
pas = e4[e4["case"] == "passive"]
mae = float(np.abs(pas.P_pred - pas.P_meas).mean())
ss_res = float(((pas.P_pred - pas.P_meas) ** 2).sum())
ss_tot = float(((pas.P_meas - pas.P_meas.mean()) ** 2).sum())
r2 = 1 - ss_res / ss_tot
print(f"\npassive block: MAE = {mae:.4f}, R^2 = {r2:.3f} "
      f"(J = {pas.J.min():g}..{pas.J.max():g})")

# =========================================================================
# E5 — theory-derived h0*(J)
# =========================================================================
XI_GRID = (0.35, 0.5, 0.65)
TARGETS = (0.10, 0.20, 0.30)


def h_eff_star(J, xi, target):
    """Solve h = T*artanh(|z_P| * xi * sqrt(lambda(h,J)/2)) for h in
    (0, h_lambda0).  Returns nan when no solution (J <= T: no instability
    at any h; then no positive-lambda window exists)."""
    if J <= T:
        return np.nan
    z = abs(norm.ppf(target))

    def F(h):
        l = lam(J, h)
        if l <= 0:
            return np.tanh(h / T)  # > 0 at the boundary
        return np.tanh(h / T) - z * xi * np.sqrt(l / 2.0)

    h_hi = T * np.arccosh(np.sqrt(J / T))  # lambda = 0 boundary
    f0 = F(0.0)
    if f0 >= 0:
        return 0.0  # target requires no bias at all
    return float(brentq(F, 0.0, h_hi - 1e-12, xtol=1e-12))


# dense theory grid (figure-ready)
Jdense = np.round(np.arange(2.5, 5.0 + 1e-9, 0.01), 4)
grid_rows = []
for xi in XI_GRID:
    for tg in TARGETS:
        for J in Jdense:
            he = h_eff_star(J, xi, tg)
            grid_rows.append({
                "xi": xi, "target": tg, "J": J, "h_eff_star": he,
                "h0_star_pred_window": he / (W0 / 500.0),   # /0.2
                "h0_star_pred_equil": he / (W_EQ / 500.0),  # /7
            })
pd.DataFrame(grid_rows).to_csv(OUT / "e5_h0_theory_grid.csv", index=False)

# per-cell predicted vs measured
sol = pd.read_csv(RES / "h0_scaling" / "h0_solutions.csv")
cell_rows = []
for _, r in sol.iterrows():
    he = h_eff_star(float(r.J), float(r.xi), float(r.target))
    hw = he / 0.2
    hq = he / 7.0
    meas = float(r.h0_star) if np.isfinite(r.h0_star) else np.nan
    in_ci = (np.isfinite(meas) and np.isfinite(r.ci_lo)
             and r.ci_lo <= hw <= r.ci_hi)
    cell_rows.append({
        "J": r.J, "xi": r.xi, "target": r.target,
        "h0_star_meas": meas, "ci_lo": r.ci_lo, "ci_hi": r.ci_hi,
        "h_eff_star_theory": he,
        "h0_star_pred_window": hw, "h0_star_pred_equil": hq,
        "ratio_window_pred_over_meas": hw / meas if meas and np.isfinite(meas) else np.nan,
        "ratio_equil_pred_over_meas": hq / meas if meas and np.isfinite(meas) else np.nan,
        "pred_window_inside_95CI": bool(in_ci),
    })
e5 = pd.DataFrame(cell_rows)
e5.to_csv(OUT / "e5_h0_pred_vs_measured.csv", index=False)

print("\n=== E5 h0*(J): predicted (window mapping h_eff = 0.2*h0) vs measured ===")
print(e5[["J", "xi", "target", "h0_star_meas", "h0_star_pred_window",
          "ratio_window_pred_over_meas", "pred_window_inside_95CI"]]
      .to_string(index=False, float_format=lambda x: f"{x:.4g}"))

ok = e5.dropna(subset=["h0_star_meas"])
mae_w = float(np.abs(ok.h0_star_pred_window - ok.h0_star_meas).mean())
mae_e = float(np.abs(ok.h0_star_pred_equil - ok.h0_star_meas).mean())
med_rw = float(ok.ratio_window_pred_over_meas.median())
med_re = float(ok.ratio_equil_pred_over_meas.median())
ss_res = float(((ok.h0_star_pred_window - ok.h0_star_meas) ** 2).sum())
ss_tot = float(((ok.h0_star_meas - ok.h0_star_meas.mean()) ** 2).sum())
r2_w = 1 - ss_res / ss_tot
n_ci = int(ok.pred_window_inside_95CI.sum())
print(f"\nwindow mapping  (h0* = h_eff*/0.2): MAE = {mae_w:.3f}, R^2 = {r2_w:.3f}, "
      f"median pred/meas = {med_rw:.3f}, inside 95% CI: {n_ci}/{len(ok)}")
print(f"equil  mapping  (h0* = h_eff*/7  ): MAE = {mae_e:.3f}, "
      f"median pred/meas = {med_re:.4f}")
better = "window" if mae_w < mae_e else "equilibrium"
print(f"=> {better} mapping fits better")

# onset locus: J_onset(h_eff) = T*cosh^2(h_eff/T); h0*->0 at J -> 2
onset_rows = [{"h_eff": h, "h0_window": h / 0.2,
               "J_onset": T * np.cosh(h / T) ** 2}
              for h in np.round(np.arange(0.0, 2.0 + 1e-9, 0.05), 4)]
pd.DataFrame(onset_rows).to_csv(OUT / "e5_onset_locus.csv", index=False)

jc = pd.read_csv(RES / "h0_scaling" / "jc_free_fits.csv")
jc_ok = jc[~jc["flagged"]]
print(f"\nonset: theory predicts h0* -> 0 at J = T = 2 exactly "
      f"(J_onset(h_eff) = 2*cosh^2(h_eff/2), = 2 at h_eff = 0).")
print(f"fitted J_c (jc_free_fits, unflagged): "
      f"{jc_ok.J_c.min():.2f}..{jc_ok.J_c.max():.2f} "
      f"(n = {len(jc_ok)}; CI union {jc_ok.ci_Jc_lo.min():.2f}.."
      f"{jc_ok.ci_Jc_hi.max():.2f})")

print(f"\nwrote {OUT}")
