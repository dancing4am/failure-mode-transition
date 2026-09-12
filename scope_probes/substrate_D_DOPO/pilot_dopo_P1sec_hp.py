"""
pilot_dopo_P1sec_hp.py -- Substrate D, P1-secondary high-precision
re-measurement (POWER ONLY). Frozen by
the DOPO high-precision amendment pilot.

Model + estimator imported / copied VERBATIM from pilot_dopo.py
(pre-reg 916f74a). The ONLY changes vs the original are statistical-
power knobs (N, burn, meas, bins, sampling cadence, bootstrap draws),
exposed as CLI args. Operating points (gamma_s=1, gamma_p=20, kappa=1,
T_eff=0.025, h=0.05, p-grid 3/4/6/9/13) and the histogram-Phi-hat
eta = dV_tilt/dV_barrier estimator are UNCHANGED.

Run (Python 3.13, numpy):
    <py> pilot_dopo_P1sec_hp.py [N] [meas] [bins] [boot] [burn] [samp]
Defaults reproduce the original locked settings (N=8000, meas=30000,
bins=160, boot=200, burn=30000, samp=50).
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np

import dopo_core as dc

MASTER_SEED = 20260518

# Locked operating points (identical to PilotCfg, pre-reg 916f74a)
GS, GP, KAPPA = 1.0, 20.0, 1.0
H = 0.05
PGRID = (3.0, 4.0, 6.0, 9.0, 13.0)
DTAU = 1e-2


def eff_pot_asymmetry(xsamp, T, nbins):
    """VERBATIM estimator from pilot_dopo.py::_eff_pot_asymmetry,
    with bin count exposed (160 in the original). Returns
    (dV_tilt, dV_barrier, saddle_count, frac_pos)."""
    lo, hi = np.percentile(xsamp, [0.2, 99.8])
    edges = np.linspace(lo, hi, nbins + 1)
    cen = 0.5 * (edges[1:] + edges[:-1])
    hist, _ = np.histogram(xsamp, bins=edges, density=True)
    counts, _ = np.histogram(xsamp, bins=edges, density=False)
    hist = np.clip(hist, 1e-9, None)
    phi = -T * np.log(hist)
    mid = len(cen) // 2
    iL = np.argmin(phi[:mid])
    iR = mid + np.argmin(phi[mid:])
    iS = iL + np.argmax(phi[iL:iR + 1])
    dV_tilt = abs(phi[iL] - phi[iR])
    dV_bar = phi[iS] - min(phi[iL], phi[iR])
    saddle_count = int(counts[iS])            # diagnostic: clip-floor?
    frac_pos = float(np.mean(xsamp > 0.0))    # diagnostic: occupancy
    return dV_tilt, dV_bar, saddle_count, frac_pos


def lin_ci(x, y, n_boot, seed):
    """VERBATIM from pilot_dopo.py::_lin_ci (OLS slope + point-resample
    bootstrap 95% CI)."""
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


def run(N, meas, nbins, n_boot, burn, samp):
    eps_th = GS * GP / KAPPA
    T = dc.t_eff_from_cavity(KAPPA, GS, GP)       # 0.025
    dt = DTAU
    rng = np.random.default_rng([MASTER_SEED, 2])  # same key path
    cst = math.sqrt(2.0 * T * dt)
    Js, etas, diag = [], [], []
    for pn in PGRID:
        eps = pn * eps_th
        x = rng.standard_normal(N) * math.sqrt(pn)
        y = np.full(N, eps / GP)
        n_snap = meas // samp + 2                  # safe upper bound
        xs = np.empty(N * n_snap, dtype=np.float64)  # pre-alloc, single copy
        k = 0
        for step in range(burn + meas):
            xd = -GS * x + KAPPA * x * y + H
            yd = -GP * y - 0.5 * KAPPA * x * x + eps
            x = x + xd * dt + cst * rng.standard_normal(N)
            y = y + yd * dt                        # pump noiseless
            if step >= burn and (step % samp == 0):
                xs[k * N:(k + 1) * N] = x           # identical to concat
                k += 1
        xs = xs[:k * N]                            # trim to actual snapshots
        dvt, dvb, sc, fp = eff_pot_asymmetry(xs, T, nbins)
        diag.append((pn, dvt, dvb, sc, fp, xs.size))
        if dvb > 1e-6:
            Js.append(pn)
            etas.append(dvt / dvb)
    Js = np.array(Js)
    etas = np.array(etas)
    if len(Js) >= 3 and np.all(etas > 0):
        sl, lo, hi = lin_ci(np.log(Js), np.log(etas),
                            n_boot, MASTER_SEED + 3)
    else:
        sl = lo = hi = float("nan")
    return Js, etas, sl, lo, hi, diag


def main():
    a = sys.argv
    N = int(a[1]) if len(a) > 1 else 8000
    meas = int(a[2]) if len(a) > 2 else 30000
    nbins = int(a[3]) if len(a) > 3 else 160
    boot = int(a[4]) if len(a) > 4 else 200
    burn = int(a[5]) if len(a) > 5 else 30000
    samp = int(a[6]) if len(a) > 6 else 50

    t0 = time.time()
    print("=== P1-secondary high-precision (POWER ONLY; "
          "AMENDMENT1 2026-05-30) ===")
    print(f"N={N} meas={meas} burn={burn} bins={nbins} boot={boot} "
          f"samp_every={samp}  (orig: 8000/30000/30000/160/200/50)")
    Js, etas, sl, lo, hi, diag = run(N, meas, nbins, boot, burn, samp)
    print("\n per-p diagnostics (pump-depletion-inclusive 2-field):")
    print("   p     dV_tilt     dV_barrier   saddle_ct  frac(x>0)  nsamp")
    for (pn, dvt, dvb, sc, fp, ns) in diag:
        print(f"  {pn:>4}  {dvt:>10.4e}  {dvb:>10.4e}  {sc:>8d}  "
              f"{fp:>8.3f}  {ns:>8d}")
    print(f"\n J grid = {Js.tolist()}")
    print(f" eta     = {[round(e, 5) for e in etas.tolist()]}")
    half = (hi - lo) / 2.0 if np.isfinite(lo) else float("nan")
    print(f" J-exponent = {sl:.3f}   95%CI [{lo:.3f}, {hi:.3f}]   "
          f"half-width = {half:.3f}")
    pinned = np.isfinite(half) and half <= 0.1
    near = lambda v: (lo <= v <= hi)
    if pinned:
        if near(-1.0) and not near(-2.0):
            verdict = "PINNED -> corrected -1 (excludes -2)"
        elif near(-2.0) and not near(-1.0):
            verdict = "PINNED -> original -2 (excludes -1) -- SURFACE"
        elif near(-1.0) and near(-2.0):
            verdict = "PINNED width<=0.1 but CI spans both -1 and -2 (check)"
        else:
            verdict = f"PINNED at ~{sl:.2f} (neither -1 nor -2)"
    else:
        verdict = ("INCONCLUSIVE (half-width > 0.1) -- power-only "
                   "insufficient; per AMENDMENT1 do NOT swap estimator "
                   "without authorization")
    print(f" VERDICT: {verdict}")
    print(f"\n=== done ({time.time()-t0:.0f}s) ===")


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    main()
