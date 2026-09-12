"""
h0_scaling.py — Stabilizer-scaling exponent measurement.

Measures h0*(J, xi, target): the stabilizer amplitude h0 at which
P(collapse) equals a target level, on the minimal Curie-Weiss + wealth
model (local generalized copy of minimal_model.run_cell), then fits

    log(h0*) = a + p * log(J)          (primary)
    log(h0*) = a' + p' * log(J - 2)    (secondary; distance to bifurcation J = T = 2)

Design notes
------------
* run_cell_general vectorizes over an EXTRA leading h0 axis: state arrays
  have shape (n_h0, n_seeds). The per-step Wiener increment is drawn ONCE
  with shape (n_seeds,) and broadcast across the h0 axis, so all h0 values
  see IDENTICAL noise paths (matched seeds across h0). Exactly one
  rng.standard_normal(n_seeds) call per step — same RNG stream as
  minimal_model.run_cell, which makes the single-h0 regression check
  bit-identical.
* Seed convention (documented base): seed = 0xB0_5CA1E
      + j_idx * 1_000_003   (j_idx: index into J_GRID)
      + xi_idx * 9973       (xi_idx: index into XI_GRID)
* Per-path collapse outcome matrices (n_h0 x n_seeds bool) are stored to
  .npz incrementally (one file per (J, xi) cell) so nothing is lost and
  runs are resumable.

Usage
-----
    python h0_scaling.py check              # mandatory regression check vs minimal_model
    python h0_scaling.py run --xi-idx K     # run all 6 J cells for XI_GRID[K]
    python h0_scaling.py analyze            # invert, fit, bootstrap, write CSVs + SUMMARY.md
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR.parent / "results" / "h0_scaling"

# ---------------------------------------------------------------------------
# Model constants (mirrors minimal_model.py — do NOT change)
# ---------------------------------------------------------------------------
T_TEMP = 2.0
DT = 0.05
N_STEPS = 20_000
INITIAL_M = 0.0
INITIAL_WEALTH = 100.0
COLLAPSE_WEALTH = 10.0
COLLAPSE_STREAK = 200
RIGIDITY_M = 0.9
FRAGMENTATION_M = 0.3
EPS = 1e-6

# ---------------------------------------------------------------------------
# Amplitude-requirement experiment configuration
# ---------------------------------------------------------------------------
SEED_BASE = 0xB0_5CA1E          # documented base seed for the stabilizer-scaling measurement
J_GRID = (2.5, 3.0, 3.5, 4.0, 4.5, 5.0)
XI_GRID = (0.35, 0.5, 0.65)
TARGETS = (0.10, 0.20, 0.30)
MULT = 100
N_SEEDS_B = 2000
COARSE_H0 = np.linspace(0.0, 6.0, 13)   # coarse bracket pass
N_REFINE = 31                           # refined grid points (>= 25)
REFINE_BAND = (0.05, 0.40)              # refine where P crosses this band
N_BOOT = 1000
BOOT_SEED_BASE = 0xB007


# ---------------------------------------------------------------------------
# Generalized engine
# ---------------------------------------------------------------------------
def run_cell_general(
    J: float,
    mult: float,
    h0_values,
    n_seeds: int = 100,
    n_steps: int = N_STEPS,
    T: float = T_TEMP,
    xi: float = 0.5,
    dt: float = DT,
    seed: int = 0,
    cons_const: float = 30.0,
    cons_coeff: float = 10.0,
    cons_form: str = "affine",
    cons_lin: float = 0.08,
) -> dict:
    """Generalized copy of minimal_model.run_cell.

    Vectorized over a leading axis of h0 values: state shape (n_h0, n_seeds).
    Noise is drawn once per step with shape (n_seeds,) and broadcast across
    the h0 axis (matched seeds across all h0 values; RNG call order identical
    to minimal_model.run_cell: exactly one standard_normal(n_seeds) per step).

    cons_form:
      'affine' : consumption = cons_const + cons_coeff * (wealth / 500)
      'linear' : consumption = cons_lin * wealth
    """
    if cons_form not in ("affine", "linear"):
        raise ValueError(f"unknown cons_form: {cons_form!r}")

    rng = np.random.default_rng(seed)

    h0 = np.atleast_1d(np.asarray(h0_values, dtype=np.float64)).reshape(-1, 1)
    n_h0 = h0.shape[0]

    m = np.full((n_h0, n_seeds), INITIAL_M, dtype=np.float64)
    wealth = np.full((n_h0, n_seeds), INITIAL_WEALTH, dtype=np.float64)

    collapsed = np.zeros((n_h0, n_seeds), dtype=bool)
    collapse_step = np.full((n_h0, n_seeds), -1, dtype=np.int32)
    collapse_type = np.zeros((n_h0, n_seeds), dtype=np.int8)
    wealth_low_streak = np.zeros((n_h0, n_seeds), dtype=np.int32)

    sqrt_dt = np.sqrt(dt)

    for t in range(n_steps):
        # --- Curie-Weiss SDE (identical op order to minimal_model) ---
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = (wealth / 500.0) * h0
        f = -m_c + np.tanh((J * m_c + h) / T)
        g = xi * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)   # ONE draw, broadcast over h0 axis
        m = np.clip(m_c + f * dt + g * dW, -1 + EPS, 1 - EPS)

        # --- Economic feedback ---
        employment = 0.5 + 0.5 * m
        income = mult * employment
        if cons_form == "affine":
            consumption = cons_const + cons_coeff * (wealth / 500.0)
        else:  # linear
            consumption = cons_lin * wealth
        wealth = np.maximum(0.0, wealth + (income - consumption) * dt)

        # --- Collapse detection ---
        is_low = wealth < COLLAPSE_WEALTH
        wealth_low_streak = np.where(is_low, wealth_low_streak + 1, 0)

        new_collapse = (wealth_low_streak >= COLLAPSE_STREAK) & (~collapsed)
        if new_collapse.any():
            collapsed |= new_collapse
            collapse_step[new_collapse] = t
            abs_m = np.abs(m)
            rigid = new_collapse & (abs_m > RIGIDITY_M)
            frag = new_collapse & (abs_m < FRAGMENTATION_M)
            mixed = new_collapse & ~(rigid | frag)
            collapse_type[rigid] = 1
            collapse_type[frag] = 2
            collapse_type[mixed] = 3

    return {
        "h0": h0.ravel(),
        "collapsed": collapsed,
        "collapse_step": collapse_step,
        "collapse_type": collapse_type,
        "final_m": m,
        "final_wealth": wealth,
    }


# ---------------------------------------------------------------------------
# Regression check (mandatory)
# ---------------------------------------------------------------------------
def regression_check(verbose: bool = True) -> bool:
    """Bit-identity check vs minimal_model.run_cell.

    Cell: J=5.0, mult=100, seed = 0xC0DE + 9*1_000_003 + 8*1009, n_seeds=100,
    defaults, single-element h0 axis at h0=2.0.
    """
    sys.path.insert(0, str(SCRIPT_DIR))
    import minimal_model  # noqa: E402

    seed = 0xC0DE + 9 * 1_000_003 + 8 * 1009
    ref = minimal_model.run_cell(J=5.0, mult=100, n_seeds=100, seed=seed)
    gen = run_cell_general(J=5.0, mult=100, h0_values=[2.0], n_seeds=100, seed=seed)

    ok_collapsed = np.array_equal(ref.collapsed, gen["collapsed"][0])
    ok_step = np.array_equal(ref.collapse_step, gen["collapse_step"][0])
    ok_type = np.array_equal(ref.collapse_type, gen["collapse_type"][0])
    ok_m = np.array_equal(ref.final_m, gen["final_m"][0])         # bitwise (bonus)
    ok_w = np.array_equal(ref.final_wealth, gen["final_wealth"][0])  # bitwise (bonus)
    ok = ok_collapsed and ok_step

    if verbose:
        print(f"[regression check] seed={seed} J=5.0 mult=100 n_seeds=100 h0=2.0")
        print(f"  collapsed bit-identical      : {ok_collapsed}")
        print(f"  collapse_step bit-identical  : {ok_step}")
        print(f"  collapse_type bit-identical  : {ok_type}")
        print(f"  final_m bit-identical        : {ok_m}")
        print(f"  final_wealth bit-identical   : {ok_w}")
        print(f"  P(collapse) ref={ref.collapsed.mean():.3f} gen={gen['collapsed'][0].mean():.3f}")
        print(f"  => REGRESSION CHECK {'PASSED' if ok else 'FAILED'}")
    return ok


# ---------------------------------------------------------------------------
# Isotonic regression (PAVA, non-increasing) and inversion
# ---------------------------------------------------------------------------
def pava_nonincreasing(y: np.ndarray) -> np.ndarray:
    """Isotonic regression, non-increasing, unit weights (pool-adjacent-violators)."""
    z = -np.asarray(y, dtype=np.float64)  # solve non-decreasing on -y
    vals: list[float] = []
    wts: list[int] = []
    for v in z:
        vals.append(float(v))
        wts.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            w = wts[-1] + wts[-2]
            vals[-2] = (vals[-1] * wts[-1] + vals[-2] * wts[-2]) / w
            wts[-2] = w
            vals.pop()
            wts.pop()
    out = np.empty_like(z)
    i = 0
    for v, w in zip(vals, wts):
        out[i:i + w] = v
        i += w
    return -out


def invert_h0_star(h0: np.ndarray, P: np.ndarray, target: float) -> float:
    """Invert the (isotonic-regressed, non-increasing) P(h0) curve at `target`.

    Returns np.nan if the target is unreachable on the sampled h0 range
    (no extrapolation): P < target even at h0=0, or P > target at max h0.
    """
    iso = pava_nonincreasing(P)
    if iso.max() < target or iso.min() > target:
        return float("nan")
    # linear interpolation on the isotonic curve; make xp strictly increasing
    xp = iso[::-1]
    fp = h0[::-1]
    xu, inv_idx = np.unique(xp, return_inverse=True)
    fu = np.bincount(inv_idx, weights=fp) / np.bincount(inv_idx)
    return float(np.interp(target, xu, fu))


# ---------------------------------------------------------------------------
# Per-cell measurement (coarse bracket + refined grid), resumable via npz
# ---------------------------------------------------------------------------
def cell_path(out_dir: Path, J: float, xi: float) -> Path:
    return out_dir / f"raw_J{J:.1f}_xi{xi:.2f}.npz"


def measure_cell(
    J: float,
    j_idx: int,
    xi: float,
    xi_idx: int,
    out_dir: Path,
    n_seeds: int = N_SEEDS_B,
    seed_base: int = SEED_BASE,
    log=print,
    **engine_kwargs,
) -> dict:
    """Coarse pass on COARSE_H0, then a refined >=25-point grid spanning the
    region where P crosses REFINE_BAND. Stores combined per-path collapse
    matrix (n_h0 x n_seeds bool) to npz. Resumable: returns cached npz if present.
    """
    path = cell_path(out_dir, J, xi)
    if path.exists():
        with np.load(path) as z:
            return {"h0": z["h0"], "collapsed": z["collapsed"], "seed": int(z["seed"])}

    seed = seed_base + j_idx * 1_000_003 + xi_idx * 9973
    t0 = time.time()
    res_c = run_cell_general(J, MULT, COARSE_H0, n_seeds=n_seeds, xi=xi, seed=seed,
                             **engine_kwargs)
    P_c = res_c["collapsed"].mean(axis=1)
    iso_c = pava_nonincreasing(P_c)
    log(f"  J={J} xi={xi} coarse done in {time.time()-t0:.1f}s  "
        f"P(h0=0)={P_c[0]:.3f} P(h0=6)={P_c[-1]:.3f}", flush=True)

    lo_band, hi_band = REFINE_BAND
    above = COARSE_H0[iso_c >= hi_band]
    below = COARSE_H0[iso_c <= lo_band]
    lo = float(above.max()) if above.size else float(COARSE_H0[0])
    hi = float(below.min()) if below.size else float(COARSE_H0[-1])
    if hi - lo < 0.5:  # degenerate bracket — pad
        lo = max(float(COARSE_H0[0]), lo - 0.5)
        hi = min(float(COARSE_H0[-1]), hi + 0.5)
    if hi <= lo:
        lo, hi = float(COARSE_H0[0]), float(COARSE_H0[-1])

    refined = np.linspace(lo, hi, N_REFINE)
    t1 = time.time()
    res_r = run_cell_general(J, MULT, refined, n_seeds=n_seeds, xi=xi, seed=seed,
                             **engine_kwargs)
    log(f"  J={J} xi={xi} refined [{lo:.2f},{hi:.2f}] x{N_REFINE} done in "
        f"{time.time()-t1:.1f}s", flush=True)

    # Combine, sort by h0, drop exact duplicates (identical seed + h0 => identical rows)
    h0_all = np.concatenate([res_c["h0"], res_r["h0"]])
    C_all = np.concatenate([res_c["collapsed"], res_r["collapsed"]], axis=0)
    h0_u, first = np.unique(h0_all, return_index=True)
    C_u = C_all[first]

    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, h0=h0_u, collapsed=C_u, seed=seed, J=J, xi=xi,
                        mult=MULT, n_seeds=n_seeds,
                        refine_lo=lo, refine_hi=hi)
    return {"h0": h0_u, "collapsed": C_u, "seed": seed}


def cmd_run(xi_idx: int) -> None:
    xi = XI_GRID[xi_idx]
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[run] xi={xi} (xi_idx={xi_idx}), seed base=0x{SEED_BASE:X}", flush=True)
    for j_idx, J in enumerate(J_GRID):
        measure_cell(J, j_idx, xi, xi_idx, RESULTS_DIR)
    print(f"[run] xi={xi} COMPLETE", flush=True)


# ---------------------------------------------------------------------------
# Analysis: monotonicity, inversion, OLS fits, bootstrap
# ---------------------------------------------------------------------------
def monotonicity_violations(h0: np.ndarray, C: np.ndarray) -> list[tuple]:
    """Paired (matched-seed) check that P is non-increasing in h0.

    For adjacent h0 pairs, b = #paths collapsing at the HIGHER h0 only,
    c = #paths collapsing at the lower h0 only. A violation beyond
    Monte-Carlo noise is flagged when b > c and (b - c) > 2*sqrt(b + c)
    (McNemar-style z > 2).
    """
    out = []
    for i in range(len(h0) - 1):
        lo_row, hi_row = C[i], C[i + 1]
        b = int(np.sum(hi_row & ~lo_row))
        c = int(np.sum(lo_row & ~hi_row))
        if b > c and (b - c) > 2.0 * np.sqrt(max(b + c, 1)):
            out.append((float(h0[i]), float(h0[i + 1]), b, c))
    return out


def ols_loglog(J_arr: np.ndarray, h0_star: np.ndarray, x_offset: float = 0.0):
    """OLS fit log(h0*) = a + p*log(J - x_offset). Returns (p, a, R2, n_used)."""
    J_arr = np.asarray(J_arr, dtype=np.float64)
    h0_star = np.asarray(h0_star, dtype=np.float64)
    mask = np.isfinite(h0_star) & (h0_star > 1e-9) & (J_arr - x_offset > 0)
    if mask.sum() < 3:
        return float("nan"), float("nan"), float("nan"), int(mask.sum())
    x = np.log(J_arr[mask] - x_offset)
    y = np.log(h0_star[mask])
    p, a = np.polyfit(x, y, 1)
    yhat = a + p * x
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(p), float(a), r2, int(mask.sum())


def bootstrap_count_matrix(n_seeds: int, n_boot: int, rng: np.random.Generator) -> np.ndarray:
    """B[s, r] = multiplicity of path s in bootstrap replicate r (same resample
    applied across the entire h0 axis — matched seeds make this meaningful)."""
    B = np.zeros((n_seeds, n_boot), dtype=np.float64)
    for r in range(n_boot):
        idx = rng.integers(0, n_seeds, size=n_seeds)
        B[:, r] = np.bincount(idx, minlength=n_seeds)
    return B


def pct_ci(vals: np.ndarray) -> tuple[float, float]:
    v = vals[np.isfinite(vals)]
    if v.size == 0:
        return float("nan"), float("nan")
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def cmd_analyze() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    sol_rows = []       # J, xi, target, h0_star, ci_lo, ci_hi
    fit_rows = []       # xi, target, p, ci_lo, ci_hi, R2, p_alt, ci_alt_lo, ci_alt_hi, R2_alt
    mono_report = []
    na_report = []

    for xi_idx, xi in enumerate(XI_GRID):
        cells = {}
        for j_idx, J in enumerate(J_GRID):
            path = cell_path(RESULTS_DIR, J, xi)
            if not path.exists():
                raise FileNotFoundError(f"missing raw cell: {path} — run `run --xi-idx {xi_idx}` first")
            with np.load(path) as z:
                cells[J] = {"h0": z["h0"].copy(), "C": z["collapsed"].copy()}

        n_seeds = cells[J_GRID[0]]["C"].shape[1]
        rng = np.random.default_rng(BOOT_SEED_BASE + xi_idx)
        B = bootstrap_count_matrix(n_seeds, N_BOOT, rng)  # shared across J within xi

        # point estimates + bootstrap h0* distributions
        h0_star_pt = {tg: np.full(len(J_GRID), np.nan) for tg in TARGETS}
        h0_star_bt = {tg: np.full((len(J_GRID), N_BOOT), np.nan) for tg in TARGETS}

        for j_idx, J in enumerate(J_GRID):
            h0 = cells[J]["h0"]
            C = cells[J]["C"]
            viol = monotonicity_violations(h0, C)
            if viol:
                mono_report.append((J, xi, viol))
            P_pt = C.mean(axis=1)
            P_bt = (C.astype(np.float64) @ B) / n_seeds  # (n_h0, N_BOOT)
            for tg in TARGETS:
                h0s = invert_h0_star(h0, P_pt, tg)
                h0_star_pt[tg][j_idx] = h0s
                if not np.isfinite(h0s):
                    na_report.append((J, xi, tg,
                                      f"P range [{P_pt.min():.3f},{P_pt.max():.3f}] on h0 in "
                                      f"[{h0.min():.2f},{h0.max():.2f}]"))
                for r in range(N_BOOT):
                    h0_star_bt[tg][j_idx, r] = invert_h0_star(h0, P_bt[:, r], tg)

        J_arr = np.asarray(J_GRID)
        for tg in TARGETS:
            for j_idx, J in enumerate(J_GRID):
                lo, hi = pct_ci(h0_star_bt[tg][j_idx])
                sol_rows.append((J, xi, tg, h0_star_pt[tg][j_idx], lo, hi))

            p, a, r2, n_used = ols_loglog(J_arr, h0_star_pt[tg], 0.0)
            p_alt, a_alt, r2_alt, n_alt = ols_loglog(J_arr, h0_star_pt[tg], 2.0)
            p_boot = np.array([ols_loglog(J_arr, h0_star_bt[tg][:, r], 0.0)[0]
                               for r in range(N_BOOT)])
            p_alt_boot = np.array([ols_loglog(J_arr, h0_star_bt[tg][:, r], 2.0)[0]
                                   for r in range(N_BOOT)])
            ci = pct_ci(p_boot)
            ci_alt = pct_ci(p_alt_boot)
            fit_rows.append((xi, tg, p, ci[0], ci[1], r2, n_used,
                             p_alt, ci_alt[0], ci_alt[1], r2_alt))
            print(f"[fit] xi={xi} target={tg}: p={p:.3f} [{ci[0]:.3f},{ci[1]:.3f}] "
                  f"R2={r2:.4f} (n={n_used}) | alt(J-2): p'={p_alt:.3f} "
                  f"[{ci_alt[0]:.3f},{ci_alt[1]:.3f}] R2={r2_alt:.4f}", flush=True)

    # ---- write outputs ----
    with open(RESULTS_DIR / "h0_solutions.csv", "w", encoding="utf-8") as f:
        f.write("J,xi,target,h0_star,ci_lo,ci_hi\n")
        for J, xi, tg, h0s, lo, hi in sol_rows:
            f.write(f"{J},{xi},{tg},{_fmt(h0s)},{_fmt(lo)},{_fmt(hi)}\n")

    with open(RESULTS_DIR / "exponent_fits.csv", "w", encoding="utf-8") as f:
        f.write("xi,target,p,ci_lo,ci_hi,R2,n_J_used,"
                "p_alt_Jminus2,ci_alt_lo,ci_alt_hi,R2_alt\n")
        for row in fit_rows:
            xi, tg, p, lo, hi, r2, n_used, p_alt, alo, ahi, r2a = row
            f.write(f"{xi},{tg},{_fmt(p)},{_fmt(lo)},{_fmt(hi)},{_fmt(r2)},{n_used},"
                    f"{_fmt(p_alt)},{_fmt(alo)},{_fmt(ahi)},{_fmt(r2a)}\n")

    _write_summary(fit_rows, sol_rows, mono_report, na_report)
    print(f"[analyze] wrote outputs to {RESULTS_DIR}", flush=True)


def _fmt(x) -> str:
    return "NA" if (x is None or not np.isfinite(x)) else f"{x:.6g}"


def _write_summary(fit_rows, sol_rows, mono_report, na_report) -> None:
    ps = [r[2] for r in fit_rows if np.isfinite(r[2])]
    lines = []
    lines.append("# Stabilizer-scaling exponent h0*(J): SUMMARY\n")
    lines.append("Model: minimal Curie-Weiss + wealth (local generalized engine, "
                 "bit-identical to `minimal_model.run_cell` at defaults; see "
                 "regression check in h0_scaling.py). mult=100, n_seeds=2000 per "
                 "(J, xi) with matched seeds across the h0 axis "
                 f"(seed base 0x{SEED_BASE:X} + j_idx*1_000_003 + xi_idx*9973). "
                 f"Bootstrap: {N_BOOT} path-level resamples, 95% percentile CIs.\n")
    lines.append("**Scope note:** the log(h0*) = a + p log(J) fit is restricted to the "
                 "six supercritical couplings J in {2.5, ..., 5.0}, i.e. J > T = 2, "
                 "because the Curie-Weiss bifurcation at J = T sits inside the "
                 "manuscript's sweep; a pure power law in J cannot hold across the "
                 "critical point. The secondary descriptor log(h0*) = a' + p' log(J - 2) "
                 "uses the distance to that bifurcation.\n")
    lines.append("## Exponent fits\n")
    lines.append("| xi | target | p (log J) | 95% CI | R^2 | p' (log(J-2)) | 95% CI | R^2 |")
    lines.append("|----|--------|-----------|--------|-----|---------------|--------|-----|")
    for xi, tg, p, lo, hi, r2, n_used, p_alt, alo, ahi, r2a in fit_rows:
        lines.append(f"| {xi} | {tg} | {_fmt(p)} | [{_fmt(lo)}, {_fmt(hi)}] | {_fmt(r2)} "
                     f"| {_fmt(p_alt)} | [{_fmt(alo)}, {_fmt(ahi)}] | {_fmt(r2a)} |")
    if ps:
        pmin, pmax = min(ps), max(ps)
        any_1_in_ci = any(r[3] <= 1.0 <= r[4] for r in fit_rows
                          if np.isfinite(r[3]) and np.isfinite(r[4]))
        all_1_in_ci = all(r[3] <= 1.0 <= r[4] for r in fit_rows
                          if np.isfinite(r[3]) and np.isfinite(r[4]))
        lines.append("\n## Verdict\n")
        lines.append(f"- Measured exponent p ranges from {pmin:.3f} to {pmax:.3f} across the "
                     f"9 (target, xi) combinations.")
        if all_1_in_ci:
            lines.append("- p = 1 lies inside the 95% CI for ALL combinations: the data are "
                         "compatible with h0* scaling proportionally with J (p ~ 1).")
        elif any_1_in_ci:
            lines.append("- p = 1 lies inside the 95% CI for SOME but not all combinations: "
                         "proportional scaling (p ~ 1) is only partially supported; see table.")
        else:
            lines.append("- p = 1 lies OUTSIDE the 95% CI for every combination: proportional "
                         "scaling h0* ~ J is NOT supported; the measured exponent is as tabled.")
        spread = pmax - pmin
        lines.append(f"- Stability: spread of p across xi in {{0.35, 0.5, 0.65}} and targets "
                     f"{{0.10, 0.20, 0.30}} is {spread:.3f}"
                     + (" (stable)." if spread < 0.3 else " (noticeable variation; see table)."))
    lines.append("\n## Monotonicity of P(collapse) in h0\n")
    if mono_report:
        lines.append("Violations beyond Monte-Carlo noise (paired McNemar-style z > 2):")
        for J, xi, viol in mono_report:
            for lo_h, hi_h, b, c in viol:
                lines.append(f"- J={J}, xi={xi}: P rose from h0={lo_h:.3f} to h0={hi_h:.3f} "
                             f"(b={b} vs c={c})")
    else:
        lines.append("No violations beyond Monte-Carlo noise detected (paired matched-seed "
                     "McNemar-style check, z > 2 threshold). P(collapse) is monotone "
                     "non-increasing in h0 up to noise.")
    lines.append("\n## Unreachable targets (recorded NA, no extrapolation)\n")
    if na_report:
        for J, xi, tg, why in na_report:
            lines.append(f"- J={J}, xi={xi}, target={tg}: unreachable ({why})")
    else:
        lines.append("None — every (J, xi, target) was invertible on the sampled h0 range.")
    lines.append("")
    with open(RESULTS_DIR / "SUMMARY.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check")
    p_run = sub.add_parser("run")
    p_run.add_argument("--xi-idx", type=int, required=True, choices=range(len(XI_GRID)))
    sub.add_parser("analyze")
    args = ap.parse_args()

    if args.cmd == "check":
        ok = regression_check()
        sys.exit(0 if ok else 1)
    elif args.cmd == "run":
        if not regression_check(verbose=True):
            print("ABORT: regression check failed", flush=True)
            sys.exit(1)
        cmd_run(args.xi_idx)
    elif args.cmd == "analyze":
        cmd_analyze()


if __name__ == "__main__":
    main()
