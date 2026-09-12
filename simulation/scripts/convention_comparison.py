# -*- coding: utf-8 -*-
"""
convention_comparison.py — window-convention comparison for the S2.7/S5.8 selection formula.

Computes the SI S2.7 / S5.8 selection-formula prediction under TWO explicit
conventions for h_eff, side by side, for EVERY variant/cell in the existing
predicted-vs-measured table, against the same measured P(collapse).

Formula (identical to selection_formula_check.py, S2.7):
    P_pred = Phi( -(c/xi) * sqrt(2/lambda) ),  c = tanh(h_eff/T),
    lambda = -1 + (J/T)*sech^2(h_eff/T),  T=2, xi=0.5.
Where lambda <= 0 the origin is linearly stable -> no instability predicted
(P_pred = 0 identically); such rows are flagged.

Two conventions for h_eff (Methods paper.md:484-486 defines the operative
field as the SELECTION-TRANSIENT value; SI S5.8 stress row reports the
DECISION-WINDOW value):

(a) SELECTION-TRANSIENT: the field during the O(1) interval just after
    initialization (W ~ W0 = 100, m ~ 0).
      passive          : h_eff = 0.4                     (= 2*(W0/500))
      fixed 2.4 / 5.0  : h_eff = 2.4 / 5.0               (constant)
      coupling (a=1)   : h_eff = 0.4*(1 + J)             (base at W0)
      quadratic (a=1)  : h_eff = 0.4*(1 + (J/5)^2)       (base at W0)
      stress   (a=2)   : h_eff = 0.4 + a*E[max(0,-m)]*J  (base pinned at W0,
                         E[max(0,-m)] from an n>=500 rerun over the O(1)
                         selection transient; base does NOT ramp)

(b) DECISION-WINDOW: realized mean field over t <= 50 (= first 1000 steps).
      fixed 2.4 / 5.0  : h_eff = 2.4 / 5.0               (constant, unchanged)
      passive/coupling/quadratic : realized per-step mean of the delivered
                         field h_fn(W(t),m(t),...) over an n>=500 rerun,
                         steps 0..999 (W ramps).
      stress           : realized window-mean field from
                         results/magnitude_matched/stress_field_stats.csv
                         (n = 2000; window t <= 50), cross-checked by rerun.

Output: results/selection_formula/convention_comparison.csv with columns
    variant, J, mu, h_eff_transient, P_pred_transient, lambda_transient,
    flag_transient_stable, h_eff_window, P_pred_window, lambda_window,
    flag_window_stable, P_measured, n_seeds, note

Run: py -3.13 convention_comparison.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import active_stabilizer as A  # committed constants + field functions

RES = HERE.parents[0] / "results"
OUT = RES / "selection_formula"
OUT.mkdir(parents=True, exist_ok=True)

T, XI = A.T_TEMP, A.XI                  # 2.0, 0.5
H0 = (A.INITIAL_WEALTH / 500.0) * 2.0   # 0.4 exactly (W0 = 100)
WINDOW_STEPS = 1000                     # t <= 50 at dt = 0.05
TRANSIENT_STEPS = 280                   # O(1) selection transient (S2.7 median)
ALPHA_STRESS = 2.0
MU = 100


def lam(J, h):
    return -1.0 + (J / T) / np.cosh(h / T) ** 2


def p_pred(J, h):
    l = lam(J, h)
    if l <= 0:
        return 0.0, l
    c = np.tanh(h / T)
    return float(norm.cdf(-(c / XI) * np.sqrt(2.0 / l))), l


def window_mean_field(h_fn, J, alpha, seed, n=500, nsteps=WINDOW_STEPS):
    """Realized per-step mean of the delivered field over steps 0..nsteps-1
    (wealth ramps from W0). Matches the S2.7 decision window t <= 50."""
    rng = np.random.default_rng(seed)
    m = np.full(n, A.INITIAL_M)
    wealth = np.full(n, A.INITIAL_WEALTH)
    sqrt_dt = np.sqrt(A.DT)
    acc_h = 0.0
    for t in range(nsteps):
        m_c = np.clip(m, -1 + A.EPS, 1 - A.EPS)
        h = h_fn(wealth, m_c, J, alpha)
        acc_h += float(np.mean(h))
        f = -m_c + np.tanh((J * m_c + h) / T)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        m = np.clip(m_c + f * A.DT + g * sqrt_dt * rng.standard_normal(n),
                    -1 + A.EPS, 1 - A.EPS)
        employment = 0.5 + 0.5 * m
        wealth = np.maximum(0.0, wealth + (MU * employment
                            - (30.0 + 10.0 * wealth / 500.0)) * A.DT)
    return acc_h / nsteps


def transient_stress_heff(J, seed, n=500, nsteps=TRANSIENT_STEPS):
    """Selection-transient h_eff for the stress form: base pinned at W0
    (=0.4, no ramp), stress term from realized E[max(0,-m)] over the O(1)
    transient. h_eff = 0.4 + alpha*E[max(0,-m)]*J."""
    rng = np.random.default_rng(seed)
    m = np.full(n, A.INITIAL_M)
    wealth = np.full(n, A.INITIAL_WEALTH)
    sqrt_dt = np.sqrt(A.DT)
    acc = 0.0
    for t in range(nsteps):
        m_c = np.clip(m, -1 + A.EPS, 1 - A.EPS)
        h = A.h_active_stress(wealth, m_c, J, ALPHA_STRESS)
        f = -m_c + np.tanh((J * m_c + h) / T)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        m = np.clip(m_c + f * A.DT + g * sqrt_dt * rng.standard_normal(n),
                    -1 + A.EPS, 1 - A.EPS)
        employment = 0.5 + 0.5 * m
        wealth = np.maximum(0.0, wealth + (MU * employment
                            - (30.0 + 10.0 * wealth / 500.0)) * A.DT)
        acc += float(np.maximum(0.0, -m).mean())
    E_neg = acc / nsteps
    return H0 + ALPHA_STRESS * E_neg * J, E_neg


rows = []


def add_row(variant, J, h_t, h_w, p_meas, n_seeds, note=""):
    Pt, lt = p_pred(J, h_t)
    Pw, lw = p_pred(J, h_w)
    rows.append({
        "variant": variant, "J": J, "mu": MU,
        "h_eff_transient": h_t, "P_pred_transient": Pt,
        "lambda_transient": lt, "flag_transient_stable": lt <= 0,
        "h_eff_window": h_w, "P_pred_window": Pw,
        "lambda_window": lw, "flag_window_stable": lw <= 0,
        "P_measured": p_meas, "n_seeds": n_seeds, "note": note,
    })


# ---- (1) passive baseline ------------------------------------------------
mm = pd.read_csv(RES / "minimal_model" / "cell_summary.csv")
mm = mm[(mm["mult"] == MU) & (mm["J"] >= 2.5)].sort_values("J")
es = pd.read_csv(RES / "extended_sweep" / "extended_sweep_summary.csv")
es = es[es["mu"] == MU].sort_values("J")
passive_cells = [(float(r.J), float(r.p_collapse), "minimal_model")
                 for _, r in mm.iterrows()]
passive_cells += [(float(r.J), float(r.p_collapse), "extended_sweep")
                  for _, r in es.iterrows()]
for J, pmeas, src in passive_cells:
    hw = window_mean_field(A.h_passive, J, 0.0, seed=A.SEED_BASE + 101_000 + int(J * 100))
    add_row("passive", J, H0, hw, pmeas, 100,
            f"transient=0.4 (W0); window=realized mean 2*(W/500) over t<=50 ({src})")

# ---- (2) fixed h = 2.4 and 5.0 ------------------------------------------
sc = pd.read_csv(RES / "active_stabilizer" / "scaling_control.csv")
for cond, h_fix in [("fixed2.4", 2.4), ("fixed5.0", 5.0)]:
    sub = sc[sc["condition"] == cond].sort_values("J")
    for _, r in sub.iterrows():
        add_row(cond, float(r.J), h_fix, h_fix, float(r.p_collapse), 100,
                "fixed field: both conventions equal the constant")

# ---- (3) coupling (alpha=1) & (4) quadratic (alpha=1), J in {4,4.5,5} ----
HIGH_J = [4.0, 4.5, 5.0]
cp = pd.read_csv(RES / "active_stabilizer" / "active_coupling_summary.csv")
cp = cp[(cp["mu"] == MU) & cp["J"].isin(HIGH_J)].sort_values("J")
for _, r in cp.iterrows():
    J = float(r.J)
    h_t = H0 * (1.0 + J)  # base at W0
    h_w = window_mean_field(A.h_active_coupling, J, 1.0,
                            seed=A.SEED_BASE + 202_000 + int(J * 100))
    add_row("coupling", J, h_t, h_w, float(r.p_collapse), 100,
            "transient=0.4*(1+J) at W0; window=realized mean base*(1+J) over t<=50")
qd = pd.read_csv(RES / "active_stabilizer" / "active_quadratic_summary.csv")
qd = qd[(qd["mu"] == MU) & qd["J"].isin(HIGH_J)].sort_values("J")
for _, r in qd.iterrows():
    J = float(r.J)
    h_t = H0 * (1.0 + (J / 5.0) ** 2)  # base at W0
    h_w = window_mean_field(A.h_active_quadratic, J, 1.0,
                            seed=A.SEED_BASE + 303_000 + int(J * 100))
    add_row("quadratic", J, h_t, h_w, float(r.p_collapse), 100,
            "transient=0.4*(1+(J/5)^2) at W0; window=realized mean over t<=50")

# ---- (5) stress form (alpha=2) ------------------------------------------
# window h_eff and measured P from magnitude_matched (n=2000, window t<=50);
# transient h_eff from the W0-base rerun (base pinned at W0, stress from m(t)).
sf = pd.read_csv(RES / "magnitude_matched" / "stress_field_stats.csv").sort_values("J")
cc = pd.read_csv(RES / "magnitude_matched" / "collapse_comparison.csv")
for _, r in sf.iterrows():
    J = float(r.J)
    h_w = float(r.window_mean_h)
    p_meas = float(cc.loc[cc.J == r.J, "p_stress"].iloc[0])
    h_t, E_neg = transient_stress_heff(J, seed=A.SEED_BASE + 404_000 + int(J * 100))
    add_row("stress", J, h_t, h_w, p_meas, 2000,
            f"transient=0.4+2*E[max(0,-m)]*J, E={E_neg:.4f} (base pinned at W0); "
            f"window=realized mean field {h_w:.3f} (magnitude_matched, n=2000)")

df = pd.DataFrame(rows)
df.to_csv(OUT / "convention_comparison.csv", index=False)

pd.set_option("display.width", 200)
show = df[["variant", "J", "h_eff_transient", "P_pred_transient",
           "flag_transient_stable", "h_eff_window", "P_pred_window",
           "flag_window_stable", "P_measured"]]
print(show.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
print(f"\nwrote {OUT / 'convention_comparison.csv'}")

# ---- headline discrepancies ---------------------------------------------
st = df[df.variant == "stress"]
print("\n--- stress-row discrepancy (measured P = 0 for all) ---")
for _, r in st.iterrows():
    print(f"J={r.J:>4}: transient h_eff={r.h_eff_transient:.3f} "
          f"P_pred={r.P_pred_transient:.3f} (|err|={abs(r.P_pred_transient-r.P_measured):.3f}); "
          f"window h_eff={r.h_eff_window:.3f} P_pred={r.P_pred_window:.3f} "
          f"(|err|={abs(r.P_pred_window-r.P_measured):.3f})")

pas = df[df.variant == "passive"]
print("\n--- passive block: MAE(pred vs meas) under each convention ---")
print(f"transient (h_eff=0.4): MAE = {np.abs(pas.P_pred_transient-pas.P_measured).mean():.4f}")
print(f"window   (realized)  : MAE = {np.abs(pas.P_pred_window-pas.P_measured).mean():.4f}")
print("passive window h_eff range: "
      f"{pas.h_eff_window.min():.3f}..{pas.h_eff_window.max():.3f}")
