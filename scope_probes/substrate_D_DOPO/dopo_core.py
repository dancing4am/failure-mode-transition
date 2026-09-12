"""
dopo_core.py -- Substrate D foundational module.

Model specification below.
Implements the rescaled single-mode degenerate-OPO normal form of
§1.1 and the c-number Langevin integrator; T_eff is computed from
cavity constants per §1.3b (NOT fitted to the transition).

Model:
    dX/dτ = (p-1) X - X^3 + h           [rescaled near-threshold]
    V(X)  = -(p-1) X^2/2 + X^4/4 - h X  [tilted cusp / Landau]
    -∂V/∂X = (p-1) X - X^3 + h          [gradient form]
Overdamped Langevin (Euler-Maruyama):
    X <- X + [-∂V/∂X] dτ + sqrt(2 T_eff dτ) ξ ,  ξ ~ N(0,1)
stationary measure ~ exp(-V(X)/T_eff)  (1-D, detailed balance).

T_eff (§1.3b, derived not assumed):
    g^2 = κ^2 / (2 γ_s γ_p)   (inverse saturation-photon-number;
    fixed by cavity constants ALONE, independent of the
    symmetry-breaking transition -- this independence is what makes
    the §1.6 P2 test non-circular).  T_eff = g^2 · f(p);
    near-threshold locked normalization f(p) = 1.

Z2 (sign) symmetry, h = 0:  X -> -X is a symmetry of the SDE ONLY
when the noise path is also negated, ξ -> -ξ (the noise realization
is not Z2-invariant).  V1 below tests the faithful bitwise statement
[(X0, ξ) -> (-X0, -ξ)] ⇒ exact mirror, plus bitwise drift oddness.
+X* and -X* are the loss-equivalent Z2 wells (signal phase 0 vs π).

Interpreter: Python 3.13 (numpy).  Run: <py3.13> dopo_core.py
"""
from __future__ import annotations

import sys
from dataclasses import dataclass

import numpy as np


# --------------------------------------------------------------------------
# Locked model
# --------------------------------------------------------------------------
@dataclass
class DopoParams:
    p: float = 1.05         # normalized pump (control; pitchfork at p=1)
    h: float = 0.0          # rescaled injected-signal field (§1.3a)
    T_eff: float = 0.05     # effective c-number noise scale (§1.3b)
    dtau: float = 1e-2      # Langevin step (rescaled time)
    seed: int = 0


def drift(X, p: float, h: float):
    """-∂V/∂X = (p-1) X - X^3 + h."""
    return (p - 1.0) * X - X ** 3 + h


def potential(X, p: float, h: float):
    """V(X) = -(p-1) X^2/2 + X^4/4 - h X."""
    return -(p - 1.0) * X ** 2 / 2.0 + X ** 4 / 4.0 - h * X


def t_eff_from_cavity(kappa: float, gamma_s: float, gamma_p: float,
                      f_reg: float = 1.0) -> float:
    """§1.3b: T_eff = g^2 · f(p), g^2 = κ^2/(2 γ_s γ_p).

    Fixed by cavity constants ALONE -- never fitted to the transition
    (the non-circularity guarantee for the §1.6 P2 test).
    Near-threshold locked normalization f_reg = 1.
    """
    g2 = kappa ** 2 / (2.0 * gamma_s * gamma_p)
    return f_reg * g2


def langevin_step(X, p, h, T_eff, dtau, noise):
    """One Euler-Maruyama step (vectorized; X and noise same shape)."""
    return X + drift(X, p, h) * dtau + np.sqrt(2.0 * T_eff * dtau) * noise


def run_ensemble(X0, p, h, T_eff, dtau, n_steps, seed,
                 record_every: int = 0):
    """Vectorized ensemble Langevin (1-D, cheap).

    X0 : scalar or array (n_traj,). Returns (X_final, recorded) where
    recorded is None unless record_every>0 (then array
    (n_rec, n_traj) of the recorded sweep)."""
    rng = np.random.default_rng(seed)
    X = np.array(X0, dtype=float)
    if X.ndim == 0:
        X = X.reshape(1)
    rec = []
    c = np.sqrt(2.0 * T_eff * dtau)
    for t in range(n_steps):
        X = X + drift(X, p, h) * dtau + c * rng.standard_normal(X.shape)
        if record_every and (t % record_every == 0):
            rec.append(X.copy())
    return X, (np.array(rec) if rec else None)


# --------------------------------------------------------------------------
# Analytic stationary distribution (1-D, exact: detailed balance)
# --------------------------------------------------------------------------
def p_ss(Xgrid, p, h, T_eff):
    """Normalized P_ss(X) ∝ exp(-V/T_eff) on a grid (trapezoid norm)."""
    V = potential(Xgrid, p, h)
    w = np.exp(-(V - V.min()) / T_eff)
    _trap = getattr(np, "trapezoid", np.trapz)   # numpy 2.x rename
    Z = _trap(w, Xgrid)
    return w / Z


# --------------------------------------------------------------------------
# V1 -- Z2 exact at h=0  (bitwise)
# --------------------------------------------------------------------------
def _v1_z2_exact():
    rng = np.random.default_rng(1)

    # (a) drift oddness, bitwise, many random X, several p
    worst_odd = 0.0
    for p in (0.3, 0.8, 1.05, 1.5, 3.0):
        X = rng.standard_normal(5000) * 2.0
        d_pos = drift(X, p, 0.0)
        d_neg = drift(-X, p, 0.0)
        if not np.array_equal(d_neg, -d_pos):
            worst_odd = max(worst_odd,
                            float(np.max(np.abs(d_neg + d_pos))))

    # (b) full-trajectory mirror under the TRUE Z2 action (X,ξ)->(-X,-ξ)
    p, T_eff, dtau, n = 1.2, 0.15, 1e-2, 4000
    xi = rng.standard_normal(n)                       # one noise path
    c = np.sqrt(2.0 * T_eff * dtau)
    Xp = 1.0
    Xm = -1.0
    max_mirror_err = 0.0
    for t in range(n):
        Xp = Xp + drift(Xp, p, 0.0) * dtau + c * xi[t]
        Xm = Xm + drift(Xm, p, 0.0) * dtau + c * (-xi[t])   # ξ -> -ξ
        e = abs(Xm - (-Xp))
        if e > max_mirror_err:
            max_mirror_err = e
    return worst_odd, max_mirror_err


# --------------------------------------------------------------------------
# V2 -- Z2 breaking linear in h, below threshold
# --------------------------------------------------------------------------
def _v2_linear_in_h():
    # p = 0 (well below threshold p=1): restoring (1-p)=1, cubic
    # negligible at the small h used (<~0.6% nonlinearity).
    p, T_eff, dtau = 0.0, 1e-2, 1e-2
    hs = np.array([0.02, 0.04, 0.06, 0.08])
    n_traj, burn, meas = 4000, 3000, 9000
    means = []
    for h in hs:
        X = np.zeros(n_traj)
        X, _ = run_ensemble(X, p, h, T_eff, dtau, burn, seed=42)
        # average over an extra measurement window AND the ensemble
        acc = np.zeros(n_traj)
        rng = np.random.default_rng(43)
        cst = np.sqrt(2.0 * T_eff * dtau)
        for _ in range(meas):
            X = X + drift(X, p, h) * dtau + cst * rng.standard_normal(n_traj)
            acc += X
        means.append(acc.mean() / meas)
    means = np.array(means)
    # predicted <X> = h/(1-p) = h  (slope 1, intercept 0)
    A = np.vstack([hs, np.ones_like(hs)]).T
    slope, intercept = np.linalg.lstsq(A, means, rcond=None)[0]
    pred = hs / (1.0 - p)
    lin_ratio = means[-1] / means[0]                # expect ~ 4
    rel_pred = float(np.max(np.abs(means - pred) / pred))
    return slope, intercept, lin_ratio, rel_pred, means


# --------------------------------------------------------------------------
# V3 -- stationary distribution matches exp(-V/T_eff)
# --------------------------------------------------------------------------
def _v3_stationary():
    p, h, T_eff, dtau = 1.2, 0.10, 0.20, 1e-2
    n_traj, n_steps = 40000, 6000      # last X = iid sample of P_ss
    rng = np.random.default_rng(7)
    X = rng.standard_normal(n_traj) * 0.5
    cst = np.sqrt(2.0 * T_eff * dtau)
    for _ in range(n_steps):
        X = X + drift(X, p, h) * dtau + cst * rng.standard_normal(n_traj)
    # empirical vs analytic CDF (KS statistic)
    xs = np.sort(X)
    ecdf = np.arange(1, len(xs) + 1) / len(xs)
    grid = np.linspace(xs.min() - 0.5, xs.max() + 0.5, 4000)
    pdf = p_ss(grid, p, h, T_eff)
    acdf_grid = np.concatenate([[0.0],
                                np.cumsum((pdf[1:] + pdf[:-1]) / 2.0
                                          * np.diff(grid))])
    acdf_grid /= acdf_grid[-1]
    acdf = np.interp(xs, grid, acdf_grid)
    ks = float(np.max(np.abs(ecdf - acdf)))
    return ks


# --------------------------------------------------------------------------
# V4 -- reversibility (adiabatic pump ramp up then down, h=0)
# --------------------------------------------------------------------------
def _v4_reversibility():
    h, T_eff, dtau = 0.0, 2e-2, 1e-2
    n_traj = 2000
    relax, meas = 4000, 2000
    ps_up = np.round(np.arange(0.5, 1.5 + 1e-9, 0.05), 3)
    ps_dn = ps_up[::-1]
    rng = np.random.default_rng(11)
    # small symmetric-ish start; warm-start continuation through ramp
    X = rng.standard_normal(n_traj) * 0.05
    cst = np.sqrt(2.0 * T_eff * dtau)

    def sweep(ps):
        nonlocal X
        out = []
        for p in ps:
            for _ in range(relax):
                X = X + drift(X, p, h) * dtau + cst * rng.standard_normal(n_traj)
            acc = 0.0
            for _ in range(meas):
                X = X + drift(X, p, h) * dtau + cst * rng.standard_normal(n_traj)
                acc += np.abs(X).mean()
            out.append(acc / meas)
        return np.array(out)

    up = sweep(ps_up)
    dn = sweep(ps_dn)[::-1]            # re-order to ascending p
    hyst = float(np.max(np.abs(up - dn)))
    return ps_up, up, dn, hyst


# --------------------------------------------------------------------------
def _tests():
    print("=== dopo_core.py validations (Substrate D, pre-reg 916f74a) ===")

    odd, mirr = _v1_z2_exact()
    ok1 = (odd == 0.0) and (mirr == 0.0)
    print(f"[V1] Z2 exact at h=0      : drift-oddness max|d(-X)+d(X)| "
          f"= {odd:.1e} ; trajectory mirror max|Xm-(-Xp)| = {mirr:.1e} "
          f"-> {'PASS' if ok1 else 'FAIL'} (need both 0, bitwise; "
          f"Z2 acts as (X,ξ)->(-X,-ξ))")

    slope, intc, lin, relp, m = _v2_linear_in_h()
    ok2 = (abs(slope - 1.0) < 0.05 and abs(intc) < 5e-3
           and abs(lin - 4.0) < 0.2 and relp < 0.05)
    print(f"[V2] Z2 breaking ∝ h (p<1): slope={slope:.4f} (pred 1.000) "
          f"intercept={intc:+.4f} ; <X>(.08)/<X>(.02)={lin:.3f} "
          f"(pred 4) ; max rel-err vs h/(1-p)={relp:.3f} "
          f"-> {'PASS' if ok2 else 'FAIL'}")

    ks = _v3_stationary()
    ok3 = ks < 0.03
    print(f"[V3] stationary ~ exp(-V/T): KS statistic = {ks:.4f} "
          f"-> {'PASS' if ok3 else 'FAIL'} (need <0.03)")

    ps, up, dn, hyst = _v4_reversibility()
    ok4 = hyst < 0.05
    print(f"[V4] reversibility (ramp)  : max|<|X|>_up - <|X|>_dn| = "
          f"{hyst:.4f} -> {'PASS' if ok4 else 'FAIL'} (need <0.05; "
          f"Clause 2: expected PASS, unlike Substrate A)")

    allok = ok1 and ok2 and ok3 and ok4
    print(f"=== dopo_core.py: "
          f"{'ALL VALIDATIONS PASSED' if allok else '*** FAILURE ***'} "
          f"===")
    return allok


if __name__ == "__main__":
    try:                                    # Windows console is cp1252
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    sys.exit(0 if _tests() else 1)
