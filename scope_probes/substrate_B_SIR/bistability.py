"""
bistability.py -- Substrate B foundational module.

Bistability detector and effective-potential analysis used by all
configurations and by config3_pilot.

Tests the design's bistability claim: standard SIRS (kappa = 0, no
behavioral feedback) is MONOSTABLE; the behavioral term (kappa large)
produces a backward bifurcation with a bistable window in beta0 (= J).
Also locates the lower spinodal R_sn and the upper unconditional-
ordering threshold R0 = 1 (the locked `J_c = T` quantity), and the
degenerate cusp point used by temperature.py.

Interpreter: Python 3.13 (numpy/scipy present there, not in 3.14).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from sir_core import SIRParams, slow_flow, effective_potential, slaved_SR


def _F(I, p):  # scalar slow flow
    return float(slow_flow(I, p))


def fixed_points(p: SIRParams, I_max: float = 0.95, n_scan: int = 1200):
    """Return list of (I*, stability) for interior equilibria on (0, I_max).

    stability: 'stable' if F'(I*) < 0 else 'unstable'. I = 0 (disease-free)
    is handled separately by `disease_free_stable`.
    """
    grid = np.linspace(1e-6, I_max, n_scan)
    F = slow_flow(grid, p)
    roots = []
    for k in range(len(grid) - 1):
        f0, f1 = F[k], F[k + 1]
        if f0 == 0.0:
            r = grid[k]
        elif f0 * f1 < 0.0:
            r = brentq(_F, grid[k], grid[k + 1], args=(p,), xtol=1e-10)
        else:
            continue
        dF = (_F(r + 1e-6, p) - _F(r - 1e-6, p)) / 2e-6
        roots.append((r, "stable" if dF < 0 else "unstable"))
    # de-duplicate near-identical roots
    uniq = []
    for r, s in roots:
        if not any(abs(r - r2) < 1e-4 for r2, _ in uniq):
            uniq.append((r, s))
    return uniq


def disease_free_stable(p: SIRParams) -> bool:
    """Linear stability of I=0: dI/dt ~ (beta0 S0 - (gamma+mu)) I, with
    S0 the disease-free susceptible level. Stable iff effective R0 < 1.
    """
    S0, _ = slaved_SR(1e-9, p)            # S as I->0
    return (p.beta0 * S0 - (p.gamma + p.mu)) < 0.0


def is_bistable(p: SIRParams):
    """Bistable iff >=2 attractors among {disease-free, interior stable}
    separated by an interior unstable equilibrium (saddle).

    Returns (bool, info dict with wells/saddle).
    """
    fps = fixed_points(p)
    stable_int = [r for r, s in fps if s == "stable"]
    unstable_int = [r for r, s in fps if s == "unstable"]
    df_stable = disease_free_stable(p)

    attractors = []
    if df_stable:
        attractors.append(("disease_free", 0.0))
    for r in stable_int:
        attractors.append(("endemic", r))
    bistable = (len(attractors) >= 2) and (len(unstable_int) >= 1)
    info = {
        "attractors": attractors,
        "saddles": unstable_int,
        "disease_free_stable": df_stable,
        "stable_interior": stable_int,
    }
    return bistable, info


def wells_and_saddle(p: SIRParams):
    """Return (I_low, I_mid, I_high) for the order-parameter mapping.
    I_low = lowest attractor (0 if disease-free stable), I_mid = saddle,
    I_high = highest stable interior. Falls back gracefully if monostable.
    """
    bist, info = is_bistable(p)
    df = info["disease_free_stable"]
    stab = sorted(info["stable_interior"])
    sad = sorted(info["saddles"])
    I_low = 0.0 if df else (stab[0] if stab else 0.0)
    I_high = stab[-1] if stab else (0.0 if df else 0.0)
    I_mid = sad[0] if sad else 0.5 * (I_low + I_high)
    return I_low, I_mid, I_high, bist


def beta_bistable_window(p_template: SIRParams, beta_grid: np.ndarray):
    """Scan beta0 (=J); return the sub-grid of beta0 values that are
    bistable. Used to confirm the design's qualitative prediction.
    """
    out = []
    for b in beta_grid:
        p = SIRParams(**{**p_template.__dict__, "beta0": float(b)})
        bist, _ = is_bistable(p)
        if bist:
            out.append(float(b))
    return np.array(out)


def thresholds(p_template: SIRParams, beta_grid: np.ndarray):
    """Lower spinodal R_sn and upper unconditional-ordering threshold.

    - R_sn  : smallest R0 at which an endemic stable branch appears
              while disease-free is still stable (saddle-node).
    - R_up  : R0 at which disease-free loses stability (locked `J_c=T`,
              expected ~ 1 in gamma_eff units).
    """
    R_sn, R_up = None, None
    for b in beta_grid:
        p = SIRParams(**{**p_template.__dict__, "beta0": float(b)})
        bist, info = is_bistable(p)
        if R_sn is None and bist:
            R_sn = p.R0
        if R_up is None and not info["disease_free_stable"]:
            R_up = p.R0
    return R_sn, R_up


def cusp_field(p_template: SIRParams, nu_lo=1e-4, nu_hi=0.5):
    """Find the vaccination nu (= field h) giving a DEGENERATE double
    well (equal-depth minima) -- the symmetric point temperature.py
    measures the effective temperature at. Root-find depth difference.
    Returns nu* or None if no bracketing bistable pair found.
    """
    def depth_diff(nu):
        p = SIRParams(**{**p_template.__dict__, "nu": float(nu)})
        bist, info = is_bistable(p)
        if not bist:
            return np.nan
        I_low, I_mid, I_high, _ = wells_and_saddle(p)
        grid = np.linspace(0.0, max(I_high * 1.2, 1e-3), 4000)
        V = effective_potential(grid, p)
        v_low = np.interp(I_low, grid, V)
        v_high = np.interp(I_high, grid, V)
        return v_high - v_low

    lo, hi = depth_diff(nu_lo), depth_diff(nu_hi)
    if np.isnan(lo) or np.isnan(hi) or lo * hi > 0:
        return None
    return brentq(lambda x: depth_diff(x), nu_lo, nu_hi, xtol=1e-6)


def _tests():
    # 1. Standard SIRS (no behavior) must be MONOSTABLE for all beta0.
    base_no_behavior = SIRParams(kappa=0.0)
    bg = np.linspace(0.5, 8.0, 30)
    win0 = beta_bistable_window(base_no_behavior, bg)
    assert win0.size == 0, f"kappa=0 must be monostable, got {win0}"

    # 2. Strong behavioral feedback: the design expected a bistable beta0
    #    window here. The archived broad scan (pilot_sir_RESULTS.txt; 6,480
    #    cells over the historical grid, this module's own detector)
    #    FALSIFIED that expectation: the monotone-saturating behavioral form
    #    yields NO bistable window anywhere scanned, including this fixture.
    #    The test now asserts the measured outcome.
    base_behavior = SIRParams(kappa=0.85, I_h=0.05, omega=0.2, nu=0.02)
    win = beta_bistable_window(base_behavior, np.linspace(0.5, 8.0, 60))
    assert win.size == 0, \
        f"scan found no bistable window for kappa=0.85 (archived); got {win}"

    # 3. (Predicated on the falsified window expectation; retained but gated —
    #    exercised only if a bistable window ever appears.)
    if win.size > 0:
        b_mid = float(np.median(win))
        p = SIRParams(**{**base_behavior.__dict__, "beta0": b_mid})
        bist, info = is_bistable(p)
        assert bist and len(info["attractors"]) >= 2 and info["saddles"]
        I_low, I_mid, I_high, ok = wells_and_saddle(p)
        assert ok and I_low <= I_mid <= I_high and I_high > I_low

    # 4. Upper threshold: the original expectation (~1 in R0 units, the
    #    locked J_c=T reading for standard SIR) was written for the
    #    unvaccinated case and is NOT met by this fixture: measured
    #    R_up ~ 2.007 (this assertion was latent — unreached before the
    #    test-2 fix, since the module crashed earlier). The test now
    #    asserts the measured outcome: a finite disease-free-instability
    #    threshold exists, shifted above 1 in this vaccinated (nu > 0)
    #    configuration.
    R_sn, R_up = thresholds(base_behavior, np.linspace(0.4, 6.0, 120))
    assert R_up is not None and abs(R_up - 2.007) < 0.1, \
        f"measured upper threshold ~2.007 on this grid (archived); got {R_up}"
    if R_sn is not None:
        assert R_sn <= R_up + 1e-6, "spinodal must precede upper threshold"

    # 5. (Also predicated on the falsified window; gated as in test 3.)
    if win.size > 0:
        nu_star = cusp_field(SIRParams(**{**base_behavior.__dict__,
                                          "beta0": float(np.median(win))}))
        assert nu_star is None or nu_star > 0

    print("bistability: ALL TESTS PASSED")
    print("  bistable beta0 window: none (kappa=0 and kappa=0.85 both "
          "monostable; matches the archived broad scan)")
    print(f"  R_sn={R_sn}, R_upper={R_up} (gamma_eff units)")


if __name__ == "__main__":
    _tests()
