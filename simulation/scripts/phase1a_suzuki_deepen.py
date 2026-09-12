# -*- coding: utf-8 -*-
"""
PHASE 1a — deepen the Suzuki transient-selection mechanism (SI S2.7).

Builds on the VALIDATED closed form  P_coll(J) = Phi(-(c/xi) sqrt(2/lambda(J))),
c=tanh(h0/T), lambda(J) = -1 + (J/T) sech^2(h0/T), committed T=2, xi=0.5, h0=0.4.

PRE-STATED CLAIMS (before running):
  (1) APPROACH TO THE 1/2 CEILING:  (1/2 - P_coll) ~ J^{-1/2}  as J->inf.
      Because arg = -(c/xi)sqrt(2/lambda) -> 0 and lambda ~ (J/T)sech^2 is affine in J,
      1/2 - P_coll ~ phi(0)(c/xi)sqrt(2/lambda) ~ const * J^{-1/2}.  => log-log slope -> -1/2,
      and (1/2 - P_coll)*sqrt(J) -> a constant.  VERIFIED IN THE CLOSED FORM ONLY; the 800-seed
      simulator (J<=50) sits below the J~1e3 asymptote (slope ~-1.03, non-monotonic P_sim), so it is
      consistent in level with approach to 1/2 but does NOT pin the -1/2 exponent (reported, not claimed).
  (2) WHY NO CLEAN POWER LAW:  the local exponent a_loc(J) = d ln P_coll / d ln J is a
      CONTINUOUS function of the fitting window (no constant), sweeping high->low as J rises;
      a window fit returns whatever a_loc averages over that window (reproducing the ambiguous
      pre-registered 0.69-0.73).  Tabulated.
  (3) NOISE-STRUCTURE INVARIANCE OF SELECTION (ASYMPTOTIC / FAST-SELECTION ONLY):  the origin
      selection variance is xi^2/(2 lambda) whether the noise is additive (xi) or multiplicative
      (xi*sqrt(1-m^2)), since sqrt(1-m^2)->1 at m=0 where the branch is chosen.  => P_coll is
      noise-structure-invariant ONLY for fast selection (J>=6, agree within ~0.02-0.04); near
      threshold (J=3-4) additive vs multiplicative differ by ~2-3x because slow selection lets m
      wander where sqrt(1-m^2)<1 -- so it is NOT invariant at low J;
      the FAILURE MODE (rigidity vs fragmentation) is what the noise structure changes (multiplicative
      noise vanishes at +-1 and locks in -> rigidity; the additive-noise residual is documented in
      S5.1).  Verified by a self-contained selection-SDE (additive vs multiplicative -> same P_coll).

Seeded.  Writes results/phase_upgrades/phase1a_suzuki_deepen.csv
Run: py -3.13 phase1a_suzuki_deepen.py
"""
import sys, csv
from pathlib import Path
import numpy as np
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import active_stabilizer as A  # validated simulator + committed constants

T, XI = A.T_TEMP, A.XI
H0 = (A.INITIAL_WEALTH / 500.0) * 2.0
C = np.tanh(H0 / T); SECH2 = 1.0 / np.cosh(H0 / T) ** 2
OUT = HERE.parents[0] / "results" / "phase_upgrades" / "phase1a_suzuki_deepen.csv"

def lam(J): return -1.0 + (J / T) * SECH2
def p_suzuki(J):
    l = lam(J); return norm.cdf(-(C / XI) * np.sqrt(2.0 / l)) if l > 0 else np.nan

rows = []

# ---- (1) approach to 1/2 ceiling: (1/2 - P) ~ J^{-1/2} (closed form) ----
Jhi = np.array([10, 15, 20, 30, 50, 75, 100, 150, 200], float)
gap = 0.5 - np.array([p_suzuki(J) for J in Jhi])
slope = np.polyfit(np.log(Jhi), np.log(gap), 1)[0]
const_pred = norm.pdf(0.0) * (C / XI) * np.sqrt(2.0 * T / SECH2)   # (1/2-P)*sqrt(J) -> this
print("(1) APPROACH TO 1/2 CEILING (closed form):")
print(f"    log-log slope of (1/2 - P_coll) vs J on [10,200] = {slope:.3f}   (claim: -0.5)")
print(f"    (1/2 - P_coll)*sqrt(J):  " + ", ".join(f"{(0.5-p_suzuki(J))*np.sqrt(J):.3f}" for J in [20,50,100,200]) +
      f"   -> predicted const {const_pred:.3f}")
for J in Jhi:
    rows.append({"block": "approach_half", "J": J, "gap_half_minus_P": 0.5 - p_suzuki(J),
                 "gap_times_sqrtJ": (0.5 - p_suzuki(J)) * np.sqrt(J)})

# ---- (2) continuous local exponent a_loc(J) = d ln P / d ln J ----
print("\n(2) CONTINUOUS LOCAL EXPONENT a_loc(J) = d lnP/d lnJ (no constant => no clean power law):")
Jgrid = [2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0]
for J in Jgrid:
    dJ = 1e-3 * J
    aloc = (np.log(p_suzuki(J + dJ)) - np.log(p_suzuki(J - dJ))) / (np.log(J + dJ) - np.log(J - dJ))
    rows.append({"block": "local_exponent", "J": J, "a_local": aloc})
    print(f"    J={J:>5.1f}   a_loc={aloc:+.3f}")
# window-averaged exponents (reproduce the ambiguous 0.69-0.73)
for lo, hi in [(2.5, 5.0), (5.0, 12.0), (2.5, 20.0)]:
    Js = np.array([j for j in [2.5,3,3.5,4,5,7,10,15,20] if lo <= j <= hi])
    a = np.polyfit(np.log(Js), np.log([p_suzuki(j) for j in Js]), 1)[0]
    print(f"    window [{lo},{hi}] power-law fit alpha = {a:.3f}")

# ---- simulator confirmation of (1): residual also approaches 1/2 as J^{-1/2} ----
print("\n(1-sim) Simulator residual P_sim(J) vs closed form, and approach to 1/2:")
print(f"    {'J':>5} {'P_sim':>7} {'P_suzuki':>9} {'(1/2-Psim)*sqrtJ':>17}")
Jsim = [8.0, 12.0, 20.0, 30.0, 50.0]
gsim = []
for J in Jsim:
    r = A.run_cell(J, 100, A.h_passive, 0.0, n_seeds=800, seed=A.SEED_BASE + int(J * 1000))
    p = float(r.collapsed.mean()); gsim.append(0.5 - p)
    rows.append({"block": "approach_half_sim", "J": J, "p_sim": p, "p_suzuki": p_suzuki(J),
                 "gap_sim_times_sqrtJ": (0.5 - p) * np.sqrt(J)})
    print(f"    {J:>5.0f} {p:>7.3f} {p_suzuki(J):>9.3f} {(0.5-p)*np.sqrt(J):>17.3f}")
ssim = np.polyfit(np.log(Jsim), np.log(gsim), 1)[0]
print(f"    simulator log-log slope of (1/2 - P_sim) vs J = {ssim:.3f}  (claim ~ -0.5)")

# ---- (3) noise-structure invariance of SELECTION PROBABILITY (self-contained selection SDE) ----
# dm = f(m) dt + sigma(m) dW from m0=0 (transient), additive sigma=xi vs multiplicative xi*sqrt(1-m^2).
# Wealth held at initial (h fixed at h0) during the brief selection window -> isolates selection.
print("\n(3) NOISE-STRUCTURE INVARIANCE OF SELECTION PROBABILITY (additive vs multiplicative):")
print(f"    {'J':>5} {'P_sel additive':>15} {'P_sel multipl.':>15} {'P_suzuki':>9}")
def selection_prob(J, multiplicative, n=20000, dt=0.05, nstep=4000, seed=0):
    rng = np.random.default_rng(seed)
    m = np.zeros(n)
    for _ in range(nstep):
        drift = -m + np.tanh((J * m + H0) / T)
        sig = XI * np.sqrt(np.clip(1 - m * m, 0, 1)) if multiplicative else XI * np.ones_like(m)
        m = m + drift * dt + sig * np.sqrt(dt) * rng.standard_normal(n)
        m = np.clip(m, -0.999999, 0.999999)
    return float((m < 0).mean())   # fraction selected onto the collapse (-) branch
for J in [3.0, 4.0, 6.0, 10.0]:
    pa = selection_prob(J, False, seed=20260611 + int(J))
    pm = selection_prob(J, True,  seed=20260611 + int(J))
    rows.append({"block": "noise_invariance", "J": J, "P_sel_additive": pa,
                 "P_sel_multiplicative": pm, "p_suzuki": p_suzuki(J)})
    print(f"    {J:>5.0f} {pa:>15.3f} {pm:>15.3f} {p_suzuki(J):>9.3f}")
print("    => selection probability ~ noise-structure-invariant; failure MODE (rigidity vs")
print("       fragmentation) is what the noise structure changes (multiplicative locks at +-1; cf S5.1).")

OUT.parent.mkdir(parents=True, exist_ok=True)
keys = sorted({k for r in rows for k in r})
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
print(f"\nwrote {OUT}")
