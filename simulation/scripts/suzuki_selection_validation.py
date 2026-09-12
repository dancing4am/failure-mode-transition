# -*- coding: utf-8 -*-
"""
Suzuki selection-at-the-unstable-origin validation (SI §S2.7).

Parameter-free leading-order prediction of the passive residual:
  P_coll(J) = Phi( -(c/xi) sqrt(2/lambda) ),  c=tanh(h0/T),  lambda=-1+(J/T)sech^2(h0/T),
with c, lambda, xi fixed by the committed constants (no fitting). Compares to the committed
scalar run_cell; also reports collapse timing (early-transient => selection) and the
window-dependent power-law exponent (=> no clean power law). Writes
results/active_stabilizer/suzuki_validation.csv.

Run: py -3.13 suzuki_selection_validation.py
"""
import sys
import csv
from pathlib import Path

import numpy as np
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import active_stabilizer as A  # noqa: E402

MU, NSEEDS = 100, 800
T, XI = A.T_TEMP, A.XI
H0 = (A.INITIAL_WEALTH / 500.0) * 2.0
C = np.tanh(H0 / T)
SECH2 = 1.0 / np.cosh(H0 / T) ** 2
JS = [2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0]
OUT = HERE.parents[0] / "results" / "active_stabilizer" / "suzuki_validation.csv"


def lam(J):
    return -1.0 + (J / T) * SECH2


def p_suzuki(J):
    l = lam(J)
    return norm.cdf(-(C / XI) * np.sqrt(2.0 / l)) if l > 0 else np.nan


print(f"H0={H0} c=f(0)={C:.4f} T={T} xi={XI}  (floor -> 0.5)")
print(f"{'J':>5s} {'lam':>7s} {'P_sim':>8s} {'P_suzuki':>9s} {'med_step/N':>11s}")
Psim, Pthy, rows = [], [], []
for J in JS:
    r = A.run_cell(J, MU, A.h_passive, 0.0, n_seeds=NSEEDS, seed=A.SEED_BASE + int(J * 1000))
    p = float(r.collapsed.mean())
    cs = r.collapse_step[r.collapsed]
    medfrac = float(np.median(cs)) / A.N_STEPS if cs.size else float("nan")
    pt = p_suzuki(J)
    Psim.append(p); Pthy.append(pt)
    rows.append({"J": J, "lambda": round(lam(J), 4), "p_sim": p,
                 "p_suzuki": round(pt, 4) if pt == pt else "nan", "med_step_frac": round(medfrac, 4)})
    print(f"{J:>5.1f} {lam(J):>7.3f} {p:>8.3f} {pt:>9.3f} {medfrac:>11.4f}")

Psim, Pthy = np.array(Psim), np.array(Pthy)
ok = ~np.isnan(Pthy)
mae = float(np.mean(np.abs(Psim[ok] - Pthy[ok])))
r2 = float(1 - np.sum((Psim[ok] - Pthy[ok]) ** 2) / np.sum((Psim[ok] - Psim[ok].mean()) ** 2))
print(f"\nSuzuki vs sim: MAE={mae:.4f}, R^2={r2:.4f}")
Ja = np.array(JS)
exps = {}
for lo, hi in [(2.5, 5.0), (5.0, 12.0), (2.5, 20.0)]:
    m = (Ja >= lo) & (Ja <= hi) & (Psim > 0)
    a1 = float(np.polyfit(np.log(Ja[m]), np.log(Psim[m]), 1)[0])
    exps[f"{lo}-{hi}"] = round(a1, 3)
    print(f"power-law alpha on J in [{lo},{hi}] = {a1:.3f}")

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["J", "lambda", "p_sim", "p_suzuki", "med_step_frac"])
    w.writeheader(); w.writerows(rows)
    f.write(f"# MAE={mae:.4f} R2={r2:.4f} power_law_alpha={exps}\n")
print(f"wrote {OUT}")
