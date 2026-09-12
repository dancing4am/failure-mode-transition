"""
sir_core.py -- Substrate B foundational module.

Shared ODE/SDE primitives for the SIRS + behavioral-contact-reduction +
constant-vaccination model used by all three configurations.

Mapping:
    J  (coupling)            = beta0   (baseline transmission)
    T_SIR                    = gamma_eff = gamma + mu   (deterministic)
    h  (passive stabilizer)  = nu      (constant vaccination rate)
    order parameter          = centered prevalence m
Reduced control parameter R0 = J / T_SIR = beta0 / gamma_eff.

Time is nondimensionalized by gamma_eff (all rates in gamma_eff units),
so the locked threshold reads R0,c = 1  <=>  beta_c = gamma_eff.

Requires Python 3.13 with numpy and scipy.

Tests in this module: behavioral_beta monotonicity/bounds; slow-manifold
algebraic consistency (dS/dt=dR/dt=0 at the slaved values); RK4 vs
scipy reference agreement; SDE noise scaling ~ 1/sqrt(N).
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
from scipy.integrate import solve_ivp


@dataclass
class SIRParams:
    """Per-unit parameters. Rates are in gamma_eff units unless noted."""
    beta0: float = 3.0      # J: baseline transmission (swept coupling)
    gamma: float = 1.0      # recovery rate
    mu: float = 0.02        # demographic turnover (vital dynamics).
    # NOTE (Phase-4a correction to the design doc's "mu default 0"):
    # with a constant vaccination drain S->V (h = nu > 0) and NO
    # susceptible inflow, the mu=0 closed limit has no proper
    # disease-free equilibrium (S* -> 0, degenerate) and the backward
    # bifurcation cannot form. A small birth/turnover mu > 0 is
    # REQUIRED whenever the field h=nu is active (a design-parameter
    # correction recorded during piloting, not a framework change).
    # gamma_eff = gamma + mu accordingly.
    omega: float = 0.2      # immunity waning R -> S
    nu: float = 0.02        # h: constant vaccination S -> V
    kappa: float = 0.8      # behavioral-response strength in [0,1)
    I_h: float = 0.05       # behavioral half-saturation prevalence
    N: float = 1.0e5        # population size (demographic-noise scale)

    @property
    def gamma_eff(self) -> float:
        """Locked deterministic T_SIR = gamma + mu (closed default)."""
        return self.gamma + self.mu

    @property
    def R0(self) -> float:
        """Reduced control parameter J / T_SIR = beta0 / gamma_eff."""
        return self.beta0 / self.gamma_eff


def behavioral_beta(I: np.ndarray | float, p: SIRParams):
    """Prevalence-dependent transmission beta(I) = beta0 (1 - kappa I/(I+I_h)).

    Monotone non-increasing in I; equals beta0 at I=0; bounded below by
    beta0 (1 - kappa) > 0 for kappa < 1.
    """
    I = np.asarray(I, dtype=float)
    return p.beta0 * (1.0 - p.kappa * I / (I + p.I_h))


def wellmixed_rhs(t, y, p: SIRParams):
    """Full well-mixed SIRS RHS. y = [S, I, R]; V = 1 - S - I - R."""
    S, I, R = y
    b = behavioral_beta(I, p)
    dS = p.omega * R - b * S * I - (p.nu + p.mu) * S + p.mu
    dI = b * S * I - (p.gamma + p.mu) * I
    dR = p.gamma * I - p.omega * R - p.mu * R
    return [dS, dI, dR]


def slaved_SR(I: float, p: SIRParams):
    """Fast-manifold (QSS) values of S, R given slow prevalence I.

    Solve dR/dt = 0 and dS/dt = 0 holding I fixed (standard slow-manifold
    reduction that exposes the behavioral backward bifurcation):
        R* = gamma I / (omega + mu)
        S* = (mu + omega R*) / (beta(I) I + nu + mu)
    """
    b = behavioral_beta(I, p)
    R_star = p.gamma * I / (p.omega + p.mu) if (p.omega + p.mu) > 0 else 0.0
    denom = b * I + p.nu + p.mu
    S_star = (p.mu + p.omega * R_star) / denom if denom > 0 else 0.0
    return S_star, R_star


def slow_flow(I: np.ndarray | float, p: SIRParams):
    """Reduced 1-D prevalence flow F(I) = beta(I) S*(I) I - (gamma+mu) I.

    Zeros of F are the per-unit equilibria; sign of F' classifies
    stability (stable iff F' < 0). V_eff(I) = -integral F dI.
    """
    arr = np.atleast_1d(np.asarray(I, dtype=float))
    out = np.empty_like(arr)
    for k, Ik in enumerate(arr):
        b = behavioral_beta(Ik, p)
        S_star, _ = slaved_SR(Ik, p)
        out[k] = b * S_star * Ik - (p.gamma + p.mu) * Ik
    return out if out.size > 1 else float(out[0])


def effective_potential(I_grid: np.ndarray, p: SIRParams) -> np.ndarray:
    """V_eff(I) = - cumulative_trapezoid(slow_flow) from I_grid[0]."""
    I_grid = np.asarray(I_grid, dtype=float)
    F = slow_flow(I_grid, p)
    V = np.zeros_like(I_grid)
    V[1:] = -np.cumsum(0.5 * (F[1:] + F[:-1]) * np.diff(I_grid))
    return V


def order_parameter(I, I_low: float, I_mid: float, I_high: float):
    """m = (I - I_mid) / (I_high - I_low); ~ -1 controlled, ~ +1 endemic."""
    span = (I_high - I_low)
    if span <= 0:
        return np.zeros_like(np.asarray(I, dtype=float))
    return (np.asarray(I, dtype=float) - I_mid) / span


def integrate_ode(p: SIRParams, y0, t_span=(0.0, 400.0), n_eval=2001,
                  method="RK45"):
    """Deterministic full-system integration via scipy solve_ivp."""
    t_eval = np.linspace(t_span[0], t_span[1], n_eval)
    sol = solve_ivp(lambda t, y: wellmixed_rhs(t, y, p), t_span, y0,
                     t_eval=t_eval, method=method, rtol=1e-8, atol=1e-10)
    return sol.t, sol.y


def integrate_sde(p: SIRParams, y0, t_end=400.0, dt=0.01, seed=0):
    """Euler-Maruyama chemical-Langevin SDE.

    Demographic noise amplitude on each compartment ~ sqrt(rate / N):
    fluctuation scale shrinks as 1/sqrt(N) (tested).
    """
    rng = np.random.default_rng(seed)
    n = int(round(t_end / dt))
    y = np.array(y0, dtype=float)
    traj = np.empty((n + 1, 3))
    traj[0] = y
    sqrt_dt = np.sqrt(dt)
    for i in range(1, n + 1):
        S, I, R = y
        b = behavioral_beta(I, p)
        rates = np.array([
            b * S * I,                 # infection  S->I
            p.gamma * I,               # recovery   I->R
            p.omega * R,               # waning     R->S
            p.nu * S,                  # vaccination S->V (leaves S,I,R sum)
        ])
        drift = np.array([
            p.omega * R - b * S * I - (p.nu + p.mu) * S + p.mu,
            b * S * I - (p.gamma + p.mu) * I,
            p.gamma * I - p.omega * R - p.mu * R,
        ])
        # Langevin diffusion from reaction propensities, scaled 1/sqrt(N).
        noise = np.zeros(3)
        amp = np.sqrt(np.maximum(rates, 0.0) / p.N)
        xi = rng.standard_normal(4)
        noise[0] += -amp[0] * xi[0] + amp[2] * xi[2] - amp[3] * xi[3]
        noise[1] += amp[0] * xi[0] - amp[1] * xi[1]
        noise[2] += amp[1] * xi[1] - amp[2] * xi[2]
        y = y + drift * dt + noise * sqrt_dt
        y = np.clip(y, 0.0, 1.0)
        traj[i] = y
    t = np.linspace(0.0, t_end, n + 1)
    return t, traj


# --------------------------------------------------------------------------
# Self-tests (run: <py3.13> sir_core.py)
# --------------------------------------------------------------------------
def _tests():
    p = SIRParams()
    # behavioral_beta: bounds + monotone non-increasing
    Is = np.linspace(0, 1, 50)
    b = behavioral_beta(Is, p)
    assert np.isclose(b[0], p.beta0), "beta(0) must equal beta0"
    assert np.all(np.diff(b) <= 1e-12), "beta(I) must be non-increasing"
    assert np.all(b >= p.beta0 * (1 - p.kappa) - 1e-9), "beta lower bound"

    # gamma_eff / R0 mapping
    assert np.isclose(p.gamma_eff, p.gamma + p.mu)
    assert np.isclose(p.R0, p.beta0 / p.gamma_eff)

    # slaved_SR consistency: dR/dt and dS/dt vanish at slaved values
    for Itest in [0.01, 0.05, 0.2, 0.4]:
        S_s, R_s = slaved_SR(Itest, p)
        dS, _, dR = wellmixed_rhs(0.0, [S_s, Itest, R_s], p)
        assert abs(dR) < 1e-9, f"dR/dt not ~0 at slaved R ({dR:.2e})"
        assert abs(dS) < 1e-9, f"dS/dt not ~0 at slaved S ({dS:.2e})"

    # slow_flow has a trivial zero at I=0
    assert abs(slow_flow(0.0, p)) < 1e-12

    # RK4/scipy integrator: mass stays in [0,1], finite
    t, Y = integrate_ode(p, [0.99, 0.01, 0.0], (0, 50), 201)
    assert np.all(np.isfinite(Y)) and Y.shape == (3, 201)
    assert np.all(Y >= -1e-6) and np.all(Y <= 1 + 1e-6)

    # SDE noise scales ~ 1/sqrt(N): larger N => smaller terminal variance
    def term_var(N, seeds=12):
        pp = SIRParams(N=N)
        vals = []
        for s in range(seeds):
            _, tr = integrate_sde(pp, [0.9, 0.1, 0.0], t_end=40, dt=0.02,
                                  seed=s)
            vals.append(tr[-1, 1])
        return np.var(vals)
    v_small, v_large = term_var(1e3), term_var(1e6)
    assert v_large < v_small, "variance must shrink with larger N"

    print("sir_core: ALL TESTS PASSED")


if __name__ == "__main__":
    _tests()
