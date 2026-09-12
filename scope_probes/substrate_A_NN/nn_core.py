"""
nn_core.py -- Substrate A cusp/fold probe, foundational module.

Continuous teacher-student perceptron, quadratic (phase-retrieval)
activation, full-batch Langevin training, linear teacher-aligned field.
Model specification below.

Model:
    z(w,x)   = (w . x) / sqrt(N)
    yhat     = g(z) = z^2          (even activation g(x)=x^2)
    teacher  : y = (w* . x)^2 / N  ,  ||w*|| = sqrt(N)  (fixed)
    empirical loss  L0(w) = (1/P) sum_mu 1/2 (yhat_mu - y_mu)^2
    field term      Lh(w) = -(h/N) (w . w*)        [C1: linear, ODD in w]
    optional norm nuisance  Lwd(w) = (lam/2) ||w||^2 / N   [EVEN; NOT h]
    total           L(w)  = L0 + Lh + Lwd

Z2 (sign) symmetry:  yhat depends on w only through (w.x)^2  =>  L0 and
Lwd are EVEN in w; the ONLY Z2-breaking term is the linear field Lh,
with the exact identity
    L(w) - L(-w) = -(2h/N) (w . w*)            (h-linear, odd).

Order parameter:  R = (w . w*) / (||w|| ||w*||) in [-1, 1].
+R* and -R* are the loss-equivalent Z2 wells; R=0 the symmetric saddle.

Langevin (unambiguous temperature T):
    w <- w - grad L * dt + sqrt(2 T dt) * xi ,  xi ~ N(0, I_N)
stationary measure ~ exp(-L(w)/T).

Interpreter: Python 3.13 (numpy present; default 3.14 lacks it).

Self-tests (run: <py3.13> nn_core.py) implement the four required
validations with quantitative evidence:
  1. Z2 at h=0 exact: L(w)=L(-w) to machine precision, many random w.
  2. Z2 breaking at h!=0: L(w)-L(-w) = -(2h/N)(w.w*), linear in h.
  3. Langevin stationary measure: the SAME integrator on an analytic
     1-D OU potential reproduces the EXACT discrete stationary
     variance T/(k(1-k dt/2)); 1-D tilted double well reproduces the
     Boltzmann occupation ratio exp(dV/T).
  4. Quadratic activation + gradient: analytic grad matches central
     finite differences to ~1e-6.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass
class ProbeParams:
    N: int = 1000           # input dimension
    alpha: float = 2.0      # load J = P/N  (P = round(alpha*N))
    T: float = 0.1          # Langevin temperature
    h: float = 0.0          # linear teacher-aligned field (the field)
    lam: float = 0.0        # optional EVEN norm nuisance (NOT h)
    dt: float = 1e-3        # Langevin step
    seed: int = 0

    @property
    def P(self) -> int:
        return int(round(self.alpha * self.N))


def make_teacher(N: int, seed: int = 0) -> np.ndarray:
    """Fixed teacher w*, normalized to ||w*|| = sqrt(N)."""
    rng = np.random.default_rng(seed)
    w = rng.standard_normal(N)
    return w * np.sqrt(N) / np.linalg.norm(w)


def make_data(wstar: np.ndarray, P: int, seed: int = 0):
    """X ~ N(0,I) of shape (P,N); teacher labels y = (X w*)^2 / N."""
    N = wstar.size
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((P, N))
    z = X @ wstar / np.sqrt(N)
    y = z ** 2 / 1.0  # z already O(1); label = g(z) = z^2
    return X, y


def predict(w: np.ndarray, X: np.ndarray) -> np.ndarray:
    """yhat = ((X w)/sqrt(N))^2  (quadratic / phase-retrieval head)."""
    N = w.size
    z = X @ w / np.sqrt(N)
    return z ** 2


def loss(w, X, y, wstar, p: ProbeParams) -> float:
    """Total loss L = L0 (MSE) + Lh (linear field) + Lwd (even nuisance)."""
    N = w.size
    yhat = predict(w, X)
    L0 = 0.5 * np.mean((yhat - y) ** 2)
    Lh = -(p.h / N) * (w @ wstar)
    Lwd = 0.5 * p.lam * (w @ w) / N
    return float(L0 + Lh + Lwd)


def grad(w, X, y, wstar, p: ProbeParams) -> np.ndarray:
    """Analytic gradient of L.

    d/dw [1/2 (yhat-y)^2] , yhat=z^2 , z=(w.x)/sqrt(N) , dz/dw = x/sqrt(N)
      => per-sample = (yhat-y) * 2 z * x/sqrt(N) = (yhat-y)*2*(w.x)/N * x
    """
    N = w.size
    P = X.shape[0]
    z = X @ w / np.sqrt(N)
    yhat = z ** 2
    resid = yhat - y                       # (P,)
    coef = resid * 2.0 * z / np.sqrt(N)     # (P,)
    g0 = (X.T @ coef) / P                   # (N,)  mean over samples
    gh = -(p.h / N) * wstar
    gwd = p.lam * w / N
    return g0 + gh + gwd


def overlap(w: np.ndarray, wstar: np.ndarray) -> float:
    """Order parameter R = (w.w*)/(||w|| ||w*||) in [-1,1]."""
    nw = np.linalg.norm(w)
    ns = np.linalg.norm(wstar)
    if nw == 0 or ns == 0:
        return 0.0
    return float((w @ wstar) / (nw * ns))


def langevin_step(state, grad_fn, T: float, dt: float, rng):
    """One Euler-Maruyama step of overdamped Langevin (generic; the
    SAME routine is used for the perceptron and for the analytic
    validation potentials).  state can be vector or scalar-as-array.
    """
    g = grad_fn(state)
    noise = rng.standard_normal(np.shape(state))
    return state - g * dt + np.sqrt(2.0 * T * dt) * noise


def run_langevin(w0, X, y, wstar, p: ProbeParams, n_steps: int,
                 record_every: int = 10):
    """Langevin trajectory of the perceptron; returns (steps, R_series,
    loss_series, w_final)."""
    rng = np.random.default_rng(p.seed)
    w = np.array(w0, dtype=float)
    Rs, Ls, ts = [], [], []
    gfn = lambda ww: grad(ww, X, y, wstar, p)
    for t in range(n_steps):
        w = langevin_step(w, gfn, p.T, p.dt, rng)
        if t % record_every == 0:
            Rs.append(overlap(w, wstar))
            Ls.append(loss(w, X, y, wstar, p))
            ts.append(t)
    return np.array(ts), np.array(Rs), np.array(Ls), w


# --------------------------------------------------------------------------
# Validations
# --------------------------------------------------------------------------
def _v1_z2_exact_at_h0():
    p = ProbeParams(N=200, alpha=1.5, h=0.0, lam=0.0, seed=1)
    wstar = make_teacher(p.N, seed=2)
    X, y = make_data(wstar, p.P, seed=3)
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(20):
        w = rng.standard_normal(p.N)
        d = abs(loss(w, X, y, wstar, p) - loss(-w, X, y, wstar, p))
        worst = max(worst, d)
    return worst  # expect ~0 (machine precision)


def _v2_z2_breaks_linear_in_h():
    p0 = ProbeParams(N=200, alpha=1.5, lam=0.0, seed=1)
    wstar = make_teacher(p0.N, seed=2)
    X, y = make_data(wstar, p0.P, seed=3)
    rng = np.random.default_rng(11)
    rel_err = []
    ratios = []
    for _ in range(10):
        w = rng.standard_normal(p0.N)
        wdotws = w @ wstar
        for h in (1e-3, 2e-3):
            p = ProbeParams(**{**p0.__dict__, "h": h})
            diff = loss(w, X, y, wstar, p) - loss(-w, X, y, wstar, p)
            analytic = -(2.0 * h / p.N) * wdotws
            rel_err.append(abs(diff - analytic) /
                           (abs(analytic) + 1e-15))
        # linearity: diff(2h)/diff(h) should be exactly 2
        pa = ProbeParams(**{**p0.__dict__, "h": 1e-3})
        pb = ProbeParams(**{**p0.__dict__, "h": 2e-3})
        da = loss(w, X, y, wstar, pa) - loss(-w, X, y, wstar, pa)
        db = loss(w, X, y, wstar, pb) - loss(-w, X, y, wstar, pb)
        ratios.append(db / da)
    return max(rel_err), (np.mean(ratios), np.std(ratios))


def _v3_langevin_stationary():
    # (a) 1-D OU: L(u)=1/2 k u^2 ; EXACT discrete EM stationary variance
    #     Var = T / (k (1 - k dt/2)).
    k, T, dt = 3.0, 0.5, 1e-3
    rng = np.random.default_rng(0)
    u = np.array(0.0)
    gfn = lambda uu: k * uu
    burn, n = 20000, 400000
    samples = np.empty(n)
    for i in range(burn + n):
        u = langevin_step(u, gfn, T, dt, rng)
        if i >= burn:
            samples[i - burn] = u
    var_meas = samples.var()
    var_exact = T / (k * (1.0 - k * dt / 2.0))
    ou_rel = abs(var_meas - var_exact) / var_exact

    # (b) 1-D tilted double well V(u)=1/4 u^4 - 1/2 u^2 - b u ; check the
    #     Boltzmann occupation ratio P(u>0)/P(u<0) = exp(dV/T) using the
    #     two well minima.
    b, Tw, dtw = 0.15, 0.30, 5e-4
    Vp = lambda uu: uu ** 3 - uu - b           # dV/du
    V = lambda uu: 0.25 * uu ** 4 - 0.5 * uu ** 2 - b * uu
    rng2 = np.random.default_rng(1)
    u = np.array(0.0)
    cp = cn = 0
    burn2, n2 = 40000, 1500000
    for i in range(burn2 + n2):
        u = langevin_step(u, Vp, Tw, dtw, rng2)
        if i >= burn2:
            if u > 0:
                cp += 1
            else:
                cn += 1
    # well minima of the symmetric part ~ +/-1 (shifted slightly by b)
    roots = np.sort(np.roots([1, 0, -1, -b]))   # u^3 - u - b = 0
    u_min_neg, _, u_min_pos = roots
    dV = V(u_min_neg) - V(u_min_pos)            # >0 favors +well
    ratio_meas = cp / max(cn, 1)
    ratio_pred = np.exp(dV / Tw)
    dw_rel = abs(np.log(ratio_meas) - np.log(ratio_pred)) / abs(
        np.log(ratio_pred))
    return ou_rel, var_meas, var_exact, dw_rel, ratio_meas, ratio_pred


def _v4_grad_finite_diff():
    p = ProbeParams(N=60, alpha=2.0, h=7e-3, lam=1e-2, seed=4)
    wstar = make_teacher(p.N, seed=5)
    X, y = make_data(wstar, p.P, seed=6)
    rng = np.random.default_rng(9)
    w = rng.standard_normal(p.N)
    g_an = grad(w, X, y, wstar, p)
    eps = 1e-6
    g_fd = np.empty_like(w)
    for i in range(p.N):
        wp = w.copy(); wp[i] += eps
        wm = w.copy(); wm[i] -= eps
        g_fd[i] = (loss(wp, X, y, wstar, p)
                   - loss(wm, X, y, wstar, p)) / (2 * eps)
    rel = np.linalg.norm(g_an - g_fd) / (np.linalg.norm(g_fd) + 1e-15)
    return rel


def _tests():
    print("=== nn_core.py validations ===")

    v1 = _v1_z2_exact_at_h0()
    ok1 = v1 < 1e-12
    print(f"[1] Z2 exact at h=0   : max|L(w)-L(-w)| = {v1:.2e}  "
          f"-> {'PASS' if ok1 else 'FAIL'} (need <1e-12)")

    rel2, (rmean, rstd) = _v2_z2_breaks_linear_in_h()
    ok2 = (rel2 < 1e-9) and (abs(rmean - 2.0) < 1e-6)
    print(f"[2] Z2 breaks ~ h     : max rel-err vs -(2h/N)(w.w*) = "
          f"{rel2:.2e}; diff(2h)/diff(h) = {rmean:.8f}+/-{rstd:.1e} "
          f"-> {'PASS' if ok2 else 'FAIL'} (need rel<1e-9, ratio=2)")

    ou_rel, vm, ve, dw_rel, rm, rp = _v3_langevin_stationary()
    ok3 = (ou_rel < 0.03) and (dw_rel < 0.10)
    print(f"[3] Langevin stationary:")
    print(f"    OU  Var meas={vm:.5f} exact={ve:.5f} rel={ou_rel:.3f} "
          f"(need <0.03)")
    print(f"    DW  occ-ratio meas={rm:.3f} Boltzmann={rp:.3f} "
          f"log-rel={dw_rel:.3f} (need <0.10)")
    print(f"    -> {'PASS' if ok3 else 'FAIL'}")

    rel4 = _v4_grad_finite_diff()
    ok4 = rel4 < 1e-6
    print(f"[4] grad vs finite-diff: rel-err = {rel4:.2e}  "
          f"-> {'PASS' if ok4 else 'FAIL'} (need <1e-6)")

    allok = ok1 and ok2 and ok3 and ok4
    print(f"=== nn_core.py: {'ALL TESTS PASSED' if allok else '*** FAILURES ***'} ===")
    return allok


if __name__ == "__main__":
    import sys
    sys.exit(0 if _tests() else 1)
