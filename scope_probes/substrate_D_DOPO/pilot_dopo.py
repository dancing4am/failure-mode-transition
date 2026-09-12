"""
pilot_dopo.py -- Substrate D pilots. Single-shot: every operating point,
grid, and ensemble size is fixed a priori below; run once; no
re-parameterization. Model imported verbatim from dopo_core.py.

Pilots:
  P1 primary  : h-exponent near the cusp (cusp=1 vs fold=½). DECISIVE.
  P1 secondary: J-exponent deep-well, pump-depletion-INCLUSIVE model
                (NOT the undepleted normal form). Confirmatory (~ -2).
  P2          : ⟨X²⟩|_{p=1} = C·T_eff^{1/2}, C parameter-free,
                T_eff from §1.3b (cavity constants, NOT fitted).
  P3          : rigidity (mean cross-barrier time) monotone in pump.

Interpreter: Python 3.13 (numpy). Run: <py3.13> pilot_dopo.py
"""
from __future__ import annotations

import math
import sys
import time
from dataclasses import dataclass

import numpy as np

import dopo_core as dc


# --------------------------------------------------------------------------
# LOCKED operating points / grids (fixed a priori; single-shot)
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class PilotCfg:
    dtau: float = 1e-2
    master_seed: int = 20260518

    # --- P1 primary: locked near-cusp operating point ---
    # p = 1.5 : a genuine double well exists (wells ±√0.5), and it is
    # clearly NOT the deep-well (p≫1) asymptotic regime reserved for
    # P1-secondary -> the legitimate "near-cusp" point. T_eff = 0.08
    # gives barrier/T_eff < 1 so the ensemble equilibrates fast; the
    # field-induced tilt is linear in h regardless of barrier depth
    # (a property of the potential, not of the sampling).
    p1_p: float = 1.5
    p1_Teff: float = 0.08
    p1_h: tuple = (0.004, 0.006, 0.009, 0.013, 0.020)
    p1_ntraj: int = 30000
    p1_burn: int = 15000
    p1_meas: int = 15000
    p1_boot: int = 400

    # --- P1 secondary: deep-well, pump-depletion-inclusive ---
    # two-field (signal x, pump y) model; γ_s=1, γ_p=20, κ=1
    # ⇒ ε_th=γ_sγ_p/κ=20, T_eff=κ²/(2γ_sγ_p)=0.025 (§1.3b).
    p1s_gs: float = 1.0
    p1s_gp: float = 20.0
    p1s_kappa: float = 1.0
    p1s_h: float = 0.05
    p1s_pgrid: tuple = (3.0, 4.0, 6.0, 9.0, 13.0)
    p1s_ntraj: int = 8000
    p1s_burn: int = 30000
    p1s_meas: int = 30000

    # --- P2: ⟨X²⟩ at the deterministic threshold p=1, h=0 ---
    # T_eff varied via cavity constants (κ²/(2γ_sγ_p)), NOT fitted.
    p2_triples: tuple = (
        (1.0, 1.0, 100.0),   # g² = 1/200  = 0.005
        (1.0, 1.0, 50.0),    # g² = 1/100  = 0.010
        (1.0, 1.0, 25.0),    # g² = 1/50   = 0.020
        (1.0, 1.0, 12.5),    # g² = 1/25   = 0.040
        (1.0, 1.0, 6.25),    # g² = 1/12.5 = 0.080
    )
    p2_ntraj: int = 40000
    p2_burn: int = 40000
    p2_meas: int = 20000

    # --- P3: rigidity (mean cross-barrier time) vs pump ---
    p3_pgrid: tuple = (1.05, 1.2, 1.5, 2.0, 3.0, 5.0)
    p3_Teff: float = 0.05
    p3_ntraj: int = 4000
    p3_maxsteps: int = 400000


CFG = PilotCfg()


def _lin_ci(x, y, n_boot, seed):
    """OLS slope + bootstrap 95% CI (resample points)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    A = np.vstack([x, np.ones_like(x)]).T
    slope = np.linalg.lstsq(A, y, rcond=None)[0][0]
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(x), len(x))
        Ab = np.vstack([x[idx], np.ones_like(x[idx])]).T
        bs.append(np.linalg.lstsq(Ab, y[idx], rcond=None)[0][0])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return float(slope), float(lo), float(hi)


# --------------------------------------------------------------------------
# P1 PRIMARY -- h-exponent near the cusp (DECISIVE)
# --------------------------------------------------------------------------
def pilot_P1_primary():
    p, T, dt = CFG.p1_p, CFG.p1_Teff, CFG.dtau
    N = CFG.p1_ntraj
    x_star = math.sqrt(p - 1.0)
    rng = np.random.default_rng([CFG.master_seed, 1])
    cst = math.sqrt(2.0 * T * dt)
    A_h = []
    per_h_signs = []
    for h in CFG.p1_h:
        X = rng.standard_normal(N) * x_star      # symmetric start
        for _ in range(CFG.p1_burn):
            X = X + dc.drift(X, p, h) * dt + cst * rng.standard_normal(N)
        acc = np.zeros(N)
        for _ in range(CFG.p1_meas):
            X = X + dc.drift(X, p, h) * dt + cst * rng.standard_normal(N)
            acc += np.sign(X)
        m = acc / CFG.p1_meas                    # per-traj mean sign
        per_h_signs.append(m)
        A_h.append(abs(m.mean()))                # |⟨sign X⟩| ∝ tilt
    A_h = np.array(A_h)
    lh = np.log(np.array(CFG.p1_h))
    lA = np.log(A_h)
    # bootstrap over trajectories (resample the N traj, recompute A(h))
    rngb = np.random.default_rng([CFG.master_seed, 11])
    slopes = []
    S = np.array(per_h_signs)                    # (n_h, N)
    for _ in range(CFG.p1_boot):
        idx = rngb.integers(0, S.shape[1], S.shape[1])
        Ab = np.abs(S[:, idx].mean(axis=1))
        Ab = np.log(np.clip(Ab, 1e-12, None))
        M = np.vstack([lh, np.ones_like(lh)]).T
        slopes.append(np.linalg.lstsq(M, Ab, rcond=None)[0][0])
    slope = float(np.polyfit(lh, lA, 1)[0])
    lo, hi = (float(np.percentile(slopes, 2.5)),
              float(np.percentile(slopes, 97.5)))
    half = (hi - lo) / 2.0
    cusp = (lo > 0.5 and lo <= 1.0 <= hi)
    fold = (hi < 1.0 and lo <= 0.5 <= hi)
    return {
        "A_h": A_h.tolist(), "h": list(CFG.p1_h),
        "slope": slope, "ci": (lo, hi), "half": half,
        "cusp": bool(cusp), "fold": bool(fold),
        "ci_ok": bool(half < 0.2),
    }


# --------------------------------------------------------------------------
# P1 SECONDARY -- deep-well J-exponent, pump-depletion-inclusive
# --------------------------------------------------------------------------
def _eff_pot_asymmetry(xsamp, T):
    """Reconstruct Φ̂ = -T ln P̂ ; return (ΔV_tilt, ΔV_barrier)."""
    lo, hi = np.percentile(xsamp, [0.2, 99.8])
    edges = np.linspace(lo, hi, 161)
    cen = 0.5 * (edges[1:] + edges[:-1])
    hist, _ = np.histogram(xsamp, bins=edges, density=True)
    hist = np.clip(hist, 1e-9, None)
    phi = -T * np.log(hist)
    mid = len(cen) // 2
    iL = np.argmin(phi[:mid])
    iR = mid + np.argmin(phi[mid:])
    iS = iL + np.argmax(phi[iL:iR + 1])          # saddle between wells
    dV_tilt = abs(phi[iL] - phi[iR])
    dV_bar = phi[iS] - min(phi[iL], phi[iR])
    return dV_tilt, dV_bar


def pilot_P1_secondary():
    gs, gp, k = CFG.p1s_gs, CFG.p1s_gp, CFG.p1s_kappa
    eps_th = gs * gp / k
    T = dc.t_eff_from_cavity(k, gs, gp)           # §1.3b, 0.025
    dt = CFG.dtau
    N = CFG.p1s_ntraj
    h = CFG.p1s_h
    rng = np.random.default_rng([CFG.master_seed, 2])
    cst = math.sqrt(2.0 * T * dt)
    Js, etas = [], []
    for pn in CFG.p1s_pgrid:
        eps = pn * eps_th
        x = rng.standard_normal(N) * math.sqrt(pn)
        y = np.full(N, eps / gp)
        for step in range(CFG.p1s_burn + CFG.p1s_meas):
            xd = -gs * x + k * x * y + h
            yd = -gp * y - 0.5 * k * x * x + eps
            x = x + xd * dt + cst * rng.standard_normal(N)
            y = y + yd * dt                        # pump noiseless (subdominant)
            if step == CFG.p1s_burn:
                samp = []
            if step >= CFG.p1s_burn and (step % 50 == 0):
                samp.append(x.copy())
        xs = np.concatenate(samp)
        dvt, dvb = _eff_pot_asymmetry(xs, T)
        if dvb > 1e-6:
            Js.append(pn)
            etas.append(dvt / dvb)
    Js = np.array(Js)
    etas = np.array(etas)
    if len(Js) >= 3 and np.all(etas > 0):
        sl, lo, hi = _lin_ci(np.log(Js), np.log(etas),
                             200, CFG.master_seed + 3)
    else:
        sl = lo = hi = float("nan")
    return {"J": Js.tolist(), "eta": etas.tolist(),
            "slope": sl, "ci": (lo, hi)}


# --------------------------------------------------------------------------
# P2 -- ⟨X²⟩ at threshold = C·T_eff^{1/2}
# --------------------------------------------------------------------------
def pilot_P2():
    # parameter-free coefficient: ⟨X²⟩ = (4T)^{1/2}·Γ(3/4)/Γ(1/4)
    C_pred = 2.0 * math.gamma(0.75) / math.gamma(0.25)
    dt = CFG.dtau
    N = CFG.p2_ntraj
    rng = np.random.default_rng([CFG.master_seed, 4])
    Ts, X2 = [], []
    for (k, gs, gp) in CFG.p2_triples:
        T = dc.t_eff_from_cavity(k, gs, gp)        # §1.3b, NOT fitted
        cst = math.sqrt(2.0 * T * dt)
        X = rng.standard_normal(N) * (4 * T) ** 0.25
        for _ in range(CFG.p2_burn):
            X = X + dc.drift(X, 1.0, 0.0) * dt + cst * rng.standard_normal(N)
        acc = 0.0
        for _ in range(CFG.p2_meas):
            X = X + dc.drift(X, 1.0, 0.0) * dt + cst * rng.standard_normal(N)
            acc += (X * X).mean()
        Ts.append(T)
        X2.append(acc / CFG.p2_meas)
    Ts = np.array(Ts)
    X2 = np.array(X2)
    sl, lo, hi = _lin_ci(np.log(Ts), np.log(X2), 400,
                         CFG.master_seed + 5)
    # coefficient at fixed exponent ½:  C_meas = mean(X2 / sqrt(T))
    C_meas = float(np.mean(X2 / np.sqrt(Ts)))
    return {"T": Ts.tolist(), "X2": X2.tolist(),
            "exp": sl, "exp_ci": (lo, hi),
            "C_pred": C_pred, "C_meas": C_meas,
            "C_rel": abs(C_meas - C_pred) / C_pred}


# --------------------------------------------------------------------------
# P3 -- rigidity: mean cross-barrier (sign-flip) time vs pump
# --------------------------------------------------------------------------
def pilot_P3():
    T, dt = CFG.p3_Teff, CFG.dtau
    N = CFG.p3_ntraj
    rng = np.random.default_rng([CFG.master_seed, 6])
    cst = math.sqrt(2.0 * T * dt)
    out = []
    for p in CFG.p3_pgrid:
        xs = math.sqrt(max(p - 1.0, 1e-6))
        X = np.full(N, xs)                         # start in + well
        flipped = np.zeros(N, bool)
        ftime = np.full(N, np.nan)
        for step in range(CFG.p3_maxsteps):
            X = X + dc.drift(X, p, 0.0) * dt + cst * rng.standard_normal(N)
            newly = (~flipped) & (X < 0.0)
            ftime[newly] = step * dt
            flipped |= newly
            if flipped.all():
                break
        mfpt = float(np.nanmean(ftime)) if np.isfinite(ftime).any() else float("inf")
        frac = float(np.mean(flipped))
        out.append((p, mfpt, frac))
    mfpts = [m for (_, m, _) in out]
    mono = all(mfpts[i + 1] >= mfpts[i] for i in range(len(mfpts) - 1))
    return {"rows": out, "monotone": bool(mono)}


# --------------------------------------------------------------------------
def main():
    t0 = time.time()
    L = []

    def emit(s=""):
        print(s)
        L.append(s)

    emit("=== pilot_dopo.py : Substrate D pilots (pre-reg 916f74a) ===")
    emit(f"locked single-shot; master_seed={CFG.master_seed}")

    emit("\n-- P1 primary (DECISIVE: h-exponent near cusp) --")
    r1 = pilot_P1_primary()
    emit(f"  op pt: p={CFG.p1_p} T_eff={CFG.p1_Teff} h={CFG.p1_h}")
    emit(f"  |⟨signX⟩|(h) = {[round(a,4) for a in r1['A_h']]}")
    verdict1 = ("CUSP (h-exp∋1, excl ½)" if r1["cusp"]
                else "FOLD (h-exp∋½, excl 1)" if r1["fold"]
                else "INCONCLUSIVE")
    emit(f"  h-exponent = {r1['slope']:.3f}  95%CI "
         f"[{r1['ci'][0]:.3f},{r1['ci'][1]:.3f}]  half={r1['half']:.3f} "
         f"(need <0.2)  -> {verdict1}")

    emit("\n-- P1 secondary (confirmatory: deep-well J-exponent, "
         "pump-depletion-inclusive) --")
    r1s = pilot_P1_secondary()
    emit(f"  J grid={r1s['J']}  η={[round(e,4) for e in r1s['eta']]}")
    emit(f"  J-exponent = {r1s['slope']:.3f}  95%CI "
         f"[{r1s['ci'][0]:.3f},{r1s['ci'][1]:.3f}]  (expect ≈ -2)")

    emit("\n-- P2 (⟨X²⟩|_{p=1} = C·T_eff^{1/2}, non-circular) --")
    r2 = pilot_P2()
    emit(f"  T_eff (from cavity §1.3b) = {[round(t,4) for t in r2['T']]}")
    emit(f"  ⟨X²⟩ = {[round(v,5) for v in r2['X2']]}")
    emit(f"  T_eff-exponent = {r2['exp']:.3f}  95%CI "
         f"[{r2['exp_ci'][0]:.3f},{r2['exp_ci'][1]:.3f}]  (predict ½)")
    emit(f"  coefficient C: predicted={r2['C_pred']:.5f}  "
         f"measured={r2['C_meas']:.5f}  rel-err={r2['C_rel']:.3f}")

    emit("\n-- P3 (rigidity: mean cross-barrier time vs pump) --")
    r3 = pilot_P3()
    for (p, m, f) in r3["rows"]:
        emit(f"  p={p:>4}: MFPT={m:>10.3f}  (flipped frac={f:.2f})")
    emit(f"  monotone increasing: {r3['monotone']}")

    emit(f"\n=== pilots done ({time.time()-t0:.0f}s) ===")

    with open("pilot_dopo_RESULTS.txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    main()
