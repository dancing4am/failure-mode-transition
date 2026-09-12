# -*- coding: utf-8 -*-
"""
PHASE 1b — rigorous derivation of the multiplicative-noise leg of Theorem 1.

PRE-STATED CLAIM (before running): for the Ito SDE
    dm = f(m) dt + xi*sqrt(1-m^2) dW,   f(m) = -m + tanh((J m + h)/T),
the stationary inter-branch bias inherits the h/J scaling of the deterministic tilt-to-
barrier ratio, NOT via an asserted "prefactor cancellation" (current SI S2.5) but because:
  (A) the Ito effective potential U_eff(m) = \int -f/(1-m^2) dm has a between-well tilt
      DeltaU = U_eff(m-) - U_eff(m+) that equals the deterministic V tilt DeltaV = V(m-) - V(m+)
      to leading order, both ~ 2h/J  (=> J*Delta/h -> 2 as J->inf);
  (B) the wells m_pm(h) are exponentially pinned at +-1 (dm_pm/dh = O(e^{-2J/T})), so the
      multiplicative prefactor / boundary-layer give only exp-suppressed h-corrections, which
      cannot enter the leading h/J term.
Hence the leading log occupation-ratio is (2/xi^2)*DeltaU ~ (4/xi^2)(h/J), prescription-
independent. We verify (A),(B); the full Ito stationary occupation ratio (by quadrature of the
exact P_ss) APPROACHES the exponent-only prediction with an O(1) prefactor correction at moderate
J (lnR/expo drifts 1.13->1.45 over J=3.5-6 as the 1/(1-m^2) factor enters) -- only the LEADING
h/J exponent is established, which is the load-bearing claim (=> Proposition, not Theorem).

Deterministic, parameter-free (committed constants). Writes results/phase_upgrades/phase1b_noise_leg.csv
Run: py -3.13 phase1b_noise_leg.py
"""
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

T, XI = 2.0, 0.5
OUT = Path(__file__).resolve().parents[1] / "results" / "phase_upgrades" / "phase1b_noise_leg.csv"

def f(m, J, h):           # drift  (= -V'(m))
    return -m + np.tanh((J * m + h) / T)

def V(m, J, h):           # deterministic potential, -V' = f
    return 0.5 * m * m - (T / J) * np.log(np.cosh((J * m + h) / T))

def wells(J, h):
    """Return (m_minus, m_plus): the two STABLE roots of f (smallest and largest of the
    three roots; the middle one is the saddle). Robust root-finding by grid sign-change scan."""
    m = np.linspace(-1 + 1e-9, 1 - 1e-9, 20001)
    fv = f(m, J, h)
    roots = []
    for i in range(len(m) - 1):
        if fv[i] == 0.0:
            roots.append(m[i])
        elif fv[i] * fv[i + 1] < 0:
            roots.append(brentq(f, m[i], m[i + 1], args=(J, h)))
    roots = sorted(roots)
    if len(roots) < 2:
        raise RuntimeError(f"expected >=2 roots at J={J},h={h}, got {roots}")
    return roots[0], roots[-1]   # (m_minus, m_plus): the two outer (stable) wells

def dU_tilt(J, h):        # U_eff(m-) - U_eff(m+)  via smooth between-well quadrature
    mm, mp = wells(J, h)
    g = lambda m: -f(m, J, h) / (1.0 - m * m)        # U_eff'(m) = -f/(1-m^2); smooth on (mm,mp)
    val, _ = quad(g, mp, mm, limit=200)
    return val

def dV_tilt(J, h):        # V(m-) - V(m+)
    mm, mp = wells(J, h)
    return V(mm, J, h) - V(mp, J, h)

def occ_ratio_ito(J, h, ngrid=200000, delta=1e-7):
    """Exact Ito stationary occupation ratio R = mass(m>0)/mass(m<0).
    P_ss(m) ∝ (1/(1-m^2)) exp(-(2/xi^2) U_eff(m)), U_eff by cumulative trapezoid from 0."""
    m = np.linspace(-1 + delta, 1 - delta, ngrid)
    g = -f(m, J, h) / (1.0 - m * m)                  # U_eff'
    Ueff = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(m))])
    logp = -np.log(1.0 - m * m) - (2.0 / XI**2) * Ueff
    logp -= logp.max()
    p = np.exp(logp)
    pos = m > 0
    massp = np.trapezoid(p[pos], m[pos]); massm = np.trapezoid(p[~pos], m[~pos])
    return massp / massm

rows = []
print(f"{'J':>5} {'h':>5} {'J*dU/h':>8} {'J*dV/h':>8} {'dU/dV':>7} {'dm+/dh':>10} {'e^-2J/T':>9}")
for J in [4.0, 6.0, 8.0, 12.0, 16.0, 20.0]:   # J>~25: wells pinned within e^{-2J/T}<1e-11 of +-1, below double precision (the claim itself)
    for h in [0.1]:
        dU = dU_tilt(J, h); dVv = dV_tilt(J, h)
        # well shift sensitivity dm+/dh (finite diff)
        _, mp1 = wells(J, h - 1e-4); _, mp2 = wells(J, h + 1e-4)
        dmp_dh = (mp2 - mp1) / 2e-4
        rows.append({"J": J, "h": h, "dU_tilt": dU, "dV_tilt": dVv,
                     "J_dU_over_h": J * dU / h, "J_dV_over_h": J * dVv / h,
                     "dmplus_dh": dmp_dh, "exp_m2JoverT": np.exp(-2 * J / T)})
        print(f"{J:>5.0f} {h:>5.2f} {J*dU/h:>8.3f} {J*dVv/h:>8.3f} {dU/dVv:>7.3f} "
              f"{dmp_dh:>10.2e} {np.exp(-2*J/T):>9.2e}")

print("\n[A] h/J law: J*dU_tilt/h -> 2 and J*dV_tilt/h -> 2 as J->inf (coeff of occ-log-ratio = (2/xi^2)*2 = 4/xi^2 = %.1f)" % (4/XI**2))
print("[A] prescription-independence: dU_tilt/dV_tilt -> 1 (multiplicative leg matches deterministic tilt)")
print("[B] well pinning: dm+/dh tracks e^{-2J/T} (exp-suppressed) -> prefactor h-correction subdominant to h/J")

# full Ito occupation ratio vs exponent-only prediction (moderate J where quadrature is well-conditioned)
print("\nFull Ito stationary occupation ratio vs exponent-only exp((2/xi^2)*dU_tilt):")
print(f"{'J':>5} {'h':>5} {'ln R_ito(quad)':>14} {'(2/xi^2)dU':>12} {'ratio':>7}")
for J in [3.5, 4.0, 5.0, 6.0]:
    h = 0.1
    R = occ_ratio_ito(J, h); lnR = np.log(R)
    expo = (2.0 / XI**2) * dU_tilt(J, h)
    rows.append({"J": J, "h": h, "lnR_ito_quad": lnR, "expo_2_over_xi2_dU": expo})
    print(f"{J:>5.1f} {h:>5.2f} {lnR:>14.3f} {expo:>12.3f} {lnR/expo:>7.3f}")

OUT.parent.mkdir(parents=True, exist_ok=True)
keys = sorted({k for r in rows for k in r})
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
print(f"\nwrote {OUT}")
