"""
temperature.py -- Substrate B foundational module.

T_SIR identification:
  PRIMARY (deterministic): T_SIR = gamma_eff = gamma + mu.
  VALIDATION (independent, non-circular): effective temperature measured
  from the stochastic order parameter at the degenerate cusp point, by
  two estimators that do NOT reuse gamma_eff:
    (i)  variance-curvature : T_var = Var(m) * V''(m_min)
    (ii) fluctuation-dissipation: T_fdt = Var(m) / chi,
         chi = d<m>/dh from a small vaccination-field step.

Per P2 (weakened): exact numeric equality of T_eff to gamma_eff is NOT
expected/claimed; we report the ratio descriptively. This module's job
is to produce the two independent estimates and the ratio, not to
enforce J_c = T_SIR.

Interpreter: Python 3.13.
"""
from __future__ import annotations

import numpy as np

from sir_core import (SIRParams, integrate_sde, slow_flow, slaved_SR,
                       order_parameter)
from bistability import cusp_field, wells_and_saddle, is_bistable


def gamma_eff(p: SIRParams) -> float:
    """Locked deterministic T_SIR."""
    return p.gamma_eff


def potential_curvature(p: SIRParams, I_at: float, dI: float = 1e-4):
    """V''(I) = -F'(I) (since V = -integral F). Local stiffness of the
    well, used by the variance-curvature temperature estimator.
    """
    Fp = (float(slow_flow(I_at + dI, p)) - float(slow_flow(I_at - dI, p))
          ) / (2 * dI)
    return -Fp  # V'' = -F'


def _stationary_m(p: SIRParams, I_start, t_end=600.0, dt=0.01, seed=0,
                  burn_frac=0.5):
    """Long SDE run; return centered order-parameter samples after burn-in.

    Order parameter centered on the saddle, scaled by well separation,
    using the deterministic wells of the SAME p.
    """
    I_low, I_mid, I_high, _ = wells_and_saddle(p)
    S0, R0 = slaved_SR(I_start, p)
    t, tr = integrate_sde(p, [S0, I_start, R0], t_end=t_end, dt=dt,
                          seed=seed)
    I = tr[:, 1]
    m = order_parameter(I, I_low, I_mid, I_high)
    k0 = int(burn_frac * len(m))
    return m[k0:], (I_low, I_mid, I_high)


def effective_temperature(p_template: SIRParams, n_seeds: int = 16,
                          t_end: float = 600.0, dt: float = 0.01,
                          dh: float = 5e-3):
    """Measure T_eff at the degenerate cusp point of p_template.

    Returns dict: nu_cusp, T_var, T_fdt, gamma_eff, ratios. If no cusp /
    no bistability is found, returns status='no_cusp' (honest failure;
    not patched).
    """
    nu_star = cusp_field(p_template)
    if nu_star is None:
        return {"status": "no_cusp", "gamma_eff": gamma_eff(p_template)}

    p0 = SIRParams(**{**p_template.__dict__, "nu": float(nu_star)})
    bist, info = is_bistable(p0)
    if not bist:
        return {"status": "no_cusp", "gamma_eff": gamma_eff(p0)}
    I_low, I_mid, I_high, _ = wells_and_saddle(p0)

    # Sample both wells; effective temperature from within-well variance.
    m_lo, m_hi, curv = [], [], []
    for s in range(n_seeds):
        ml, _ = _stationary_m(p0, max(I_low, 1e-4), t_end, dt, seed=s)
        mh, _ = _stationary_m(p0, max(I_high, I_mid + 1e-3), t_end, dt,
                              seed=1000 + s)
        # keep samples that stayed in their starting well (no barrier
        # crossing) so the variance is the *within-well* fluctuation.
        if np.mean(ml) < 0:
            m_lo.append(ml[ml < 0])
        if np.mean(mh) > 0:
            m_hi.append(mh[mh > 0])
    m_lo = np.concatenate(m_lo) if m_lo else np.array([])
    m_hi = np.concatenate(m_hi) if m_hi else np.array([])
    if m_lo.size < 50 or m_hi.size < 50:
        return {"status": "insufficient_well_samples",
                "gamma_eff": gamma_eff(p0), "nu_cusp": nu_star}

    var_well = 0.5 * (np.var(m_lo) + np.var(m_hi))
    kV = abs(0.5 * (potential_curvature(p0, max(I_low, 1e-4))
                    + potential_curvature(p0, max(I_high, I_mid + 1e-3))))
    T_var = var_well * kV

    # Fluctuation-dissipation: small vaccination-field step -> response.
    p_plus = SIRParams(**{**p0.__dict__, "nu": float(nu_star + dh)})
    p_minus = SIRParams(**{**p0.__dict__, "nu": float(nu_star - dh)})
    mp, mm = [], []
    for s in range(n_seeds):
        a, _ = _stationary_m(p_plus, I_high, t_end, dt, seed=7000 + s)
        b, _ = _stationary_m(p_minus, I_high, t_end, dt, seed=8000 + s)
        mp.append(np.mean(a))
        mm.append(np.mean(b))
    chi = (np.mean(mp) - np.mean(mm)) / (2 * dh)
    T_fdt = var_well / chi if chi != 0 and np.isfinite(chi) else np.nan

    g = gamma_eff(p0)
    return {
        "status": "ok",
        "nu_cusp": float(nu_star),
        "var_well": float(var_well),
        "curvature": float(kV),
        "T_var": float(T_var),
        "T_fdt": float(T_fdt),
        "gamma_eff": float(g),
        "ratio_Tvar_over_gamma_eff": float(T_var / g) if g else np.nan,
        "ratio_Tfdt_over_gamma_eff": (float(T_fdt / g)
                                      if g and np.isfinite(T_fdt)
                                      else np.nan),
    }


def _tests():
    p = SIRParams()
    # gamma_eff is the locked deterministic value
    assert np.isclose(gamma_eff(p), p.gamma + p.mu)

    # curvature sign: at a stable well V'' > 0 (F' < 0)
    from bistability import beta_bistable_window
    base = SIRParams(kappa=0.85, I_h=0.05, omega=0.2, nu=0.02)
    win = beta_bistable_window(base, np.linspace(0.5, 8.0, 60))
    assert win.size > 0
    p_b = SIRParams(**{**base.__dict__, "beta0": float(np.median(win))})
    I_low, I_mid, I_high, ok = wells_and_saddle(p_b)
    assert ok
    if I_high > I_mid:
        assert potential_curvature(p_b, I_high) > 0, "well must be convex"
    # saddle is concave (V'' < 0)
    assert potential_curvature(p_b, I_mid) < 0, "saddle must be concave"

    # effective_temperature returns a well-formed result or an honest
    # no_cusp status (small/fast settings for the unit test)
    res = effective_temperature(p_b, n_seeds=4, t_end=120.0, dt=0.02)
    assert res["status"] in (
        "ok", "no_cusp", "insufficient_well_samples")
    assert "gamma_eff" in res
    if res["status"] == "ok":
        assert res["T_var"] > 0 and np.isfinite(res["T_var"])

    print("temperature: ALL TESTS PASSED")
    print(f"  effective_temperature status={res['status']}, "
          f"gamma_eff={res['gamma_eff']:.4f}"
          + (f", T_var={res['T_var']:.4g}, T_fdt={res.get('T_fdt')}"
             if res['status'] == 'ok' else ""))


if __name__ == "__main__":
    _tests()
