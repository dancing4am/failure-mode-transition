"""
jc_free_fit.py — Free-J_c three-parameter fit of the stabilizer scaling.

Context
-------
The amplitude-requirement fit (h0_scaling.py) fitted log(h0*) = a' + p' log(J - 2), hardcoding the deterministic
Curie-Weiss bifurcation J_c = T = 2. But h != 0 makes the pitchfork imperfect
(h(W) spans ~0..14 against T = 2), so there is no sharp bifurcation at exactly
J = T. This script asks the data where J_c actually is:

    h0*(J) = A * (J - J_c)^p'      with J_c FREE.

Method (profile grid search — robust with only 6 points)
--------------------------------------------------------
For each (xi, target) cell with >= 5 valid J points (the n=4 cell
xi=0.65/target=0.10 is fitted too but FLAGGED and excluded from the overall
consistency count):
  * grid J_c in [-1.0, 2.45] step 0.005, enforcing J_c <= min(J) - 0.05;
  * for each J_c: OLS of log(h0*) on log(J - J_c);
  * pick J_c minimizing SSE; report A, p', J_c, R^2.

Bootstrap (N_BOOT = 1000 path-level resamples)
----------------------------------------------
Reuses the stored per-path collapse matrices (raw_J{J}_xi{xi}.npz). The same
resample of the 2000 path indices is applied across the whole h0 axis within a
(J, xi) file (matched seeds make this meaningful), and — mirroring
h0_scaling.cmd_analyze — the same bootstrap count matrix is shared across J
within a xi block (seed BOOT_SEED_BASE + xi_idx), so the h0* bootstrap
distributions are identical to the amplitude-requirement fit's. Per replicate: recompute P(h0),
re-invert h0*(J, target) with the SAME isotonic-PAVA + linear-interp inversion
as h0_scaling.py (imported, not copied), then re-run the profile fit.
  * If a target is unreachable at some J in a replicate, that J is dropped
    from that replicate's fit (frequency recorded).
  * If a replicate has < 4 valid J points, the replicate is dropped
    (count recorded).
95% percentile CIs on J_c and p'.

Outputs (results/h0_scaling/)
-----------------------------
  jc_free_fits.csv  — xi, target, n_points, A, p_prime, J_c, ci_Jc_lo,
                      ci_Jc_hi, ci_pprime_lo, ci_pprime_hi, R2, contains_2,
                      flagged
  B4_SUMMARY.md     — decision verdict up top.

Usage
-----
    python jc_free_fit.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

# Mirror the amplitude-requirement fit's inversion and bootstrap machinery EXACTLY by importing it.
from h0_scaling import (  # noqa: E402
    BOOT_SEED_BASE,
    J_GRID,
    RESULTS_DIR,
    TARGETS,
    XI_GRID,
    bootstrap_count_matrix,
    cell_path,
    invert_h0_star,
    pct_ci,
)

N_BOOT = 1000
JC_LO, JC_HI, JC_STEP = -1.0, 2.45, 0.005
JC_GRID = np.round(np.arange(JC_LO, JC_HI + 1e-9, JC_STEP), 3)
MIN_POINTS = 4          # replicate fits with fewer valid J points are dropped
FLAG_CELL = (0.65, 0.10)  # n=4 cell — fitted but not comparable


# ---------------------------------------------------------------------------
# Profile fit: grid over J_c, OLS of log(h0*) on log(J - J_c) at each J_c
# ---------------------------------------------------------------------------
def profile_fit(J_arr: np.ndarray, h0_star: np.ndarray) -> dict:
    """3-parameter fit h0* = A (J - J_c)^p' via profile grid search over J_c.

    Returns dict with keys: valid, n, Jc, p, A, R2, at_lo, at_hi.
    at_lo/at_hi flag a J_c estimate pinned at a grid boundary.
    """
    J_arr = np.asarray(J_arr, dtype=np.float64)
    h0_star = np.asarray(h0_star, dtype=np.float64)
    mask = np.isfinite(h0_star) & (h0_star > 1e-9)
    n = int(mask.sum())
    if n < MIN_POINTS:
        return {"valid": False, "n": n}

    Jv = J_arr[mask]
    y = np.log(h0_star[mask])
    jc = JC_GRID[JC_GRID <= Jv.min() - 0.05 + 1e-9]  # enforce J_c < min(J) - 0.05
    if jc.size == 0:
        return {"valid": False, "n": n}

    X = np.log(Jv[None, :] - jc[:, None])            # (G, n)
    xm = X.mean(axis=1)
    ym = float(y.mean())
    xc = X - xm[:, None]
    yc = y - ym
    sxx = np.einsum("gi,gi->g", xc, xc)
    sxy = xc @ yc
    slope = sxy / sxx
    ss_tot = float(yc @ yc)
    sse = np.maximum(ss_tot - slope * sxy, 0.0)      # SSE = Syy - b*Sxy
    i = int(np.argmin(sse))
    p = float(slope[i])
    a = ym - p * float(xm[i])
    r2 = 1.0 - sse[i] / ss_tot if ss_tot > 0 else float("nan")
    return {
        "valid": True,
        "n": n,
        "Jc": float(jc[i]),
        "p": p,
        "A": float(np.exp(a)),
        "R2": float(r2),
        "at_lo": i == 0,
        "at_hi": i == jc.size - 1,
    }


# ---------------------------------------------------------------------------
def _fmt(x) -> str:
    if x is None:
        return "NA"
    try:
        return "NA" if not np.isfinite(x) else f"{x:.6g}"
    except TypeError:
        return str(x)


def main() -> None:
    J_arr = np.asarray(J_GRID, dtype=np.float64)
    t_start = time.time()

    rows = []  # one dict per (xi, target)

    for xi_idx, xi in enumerate(XI_GRID):
        # ---- load stored per-path matrices for this xi block ----
        cells = {}
        for J in J_GRID:
            path = cell_path(RESULTS_DIR, J, xi)
            if not path.exists():
                raise FileNotFoundError(f"missing raw cell: {path}")
            with np.load(path) as z:
                cells[J] = {"h0": z["h0"].copy(), "C": z["collapsed"].copy()}
        n_seeds = cells[J_GRID[0]]["C"].shape[1]

        # ---- bootstrap count matrix: identical seeding to h0_scaling ----
        rng = np.random.default_rng(BOOT_SEED_BASE + xi_idx)
        B = bootstrap_count_matrix(n_seeds, N_BOOT, rng)  # shared across J within xi

        # ---- point-estimate P and bootstrap P, then invert h0* ----
        h0_star_pt = {tg: np.full(len(J_GRID), np.nan) for tg in TARGETS}
        h0_star_bt = {tg: np.full((len(J_GRID), N_BOOT), np.nan) for tg in TARGETS}
        for j_idx, J in enumerate(J_GRID):
            h0 = cells[J]["h0"]
            C = cells[J]["C"]
            P_pt = C.mean(axis=1)
            P_bt = (C.astype(np.float64) @ B) / n_seeds   # (n_h0, N_BOOT)
            for tg in TARGETS:
                h0_star_pt[tg][j_idx] = invert_h0_star(h0, P_pt, tg)
                for r in range(N_BOOT):
                    h0_star_bt[tg][j_idx, r] = invert_h0_star(h0, P_bt[:, r], tg)
        print(f"[boot] xi={xi}: inversion of {N_BOOT} replicates done "
              f"({time.time()-t_start:.1f}s elapsed)", flush=True)

        # ---- per (xi, target): point fit + bootstrap fits ----
        for tg in TARGETS:
            pt = profile_fit(J_arr, h0_star_pt[tg])
            flagged = (xi, tg) == FLAG_CELL

            jc_b = np.full(N_BOOT, np.nan)
            pp_b = np.full(N_BOOT, np.nan)
            n_dropped = 0
            n_boundary = 0
            # per-J unreachable frequency across replicates
            unreach = {J: int(np.sum(~np.isfinite(h0_star_bt[tg][j, :])))
                       for j, J in enumerate(J_GRID)}
            for r in range(N_BOOT):
                fb = profile_fit(J_arr, h0_star_bt[tg][:, r])
                if not fb["valid"]:
                    n_dropped += 1
                    continue
                jc_b[r] = fb["Jc"]
                pp_b[r] = fb["p"]
                if fb["at_lo"] or fb["at_hi"]:
                    n_boundary += 1

            ci_jc = pct_ci(jc_b)
            ci_pp = pct_ci(pp_b)
            contains_2 = (np.isfinite(ci_jc[0]) and np.isfinite(ci_jc[1])
                          and ci_jc[0] <= 2.0 <= ci_jc[1])
            rows.append({
                "xi": xi, "target": tg,
                "n_points": pt.get("n", 0),
                "A": pt.get("A", np.nan), "p_prime": pt.get("p", np.nan),
                "J_c": pt.get("Jc", np.nan), "R2": pt.get("R2", np.nan),
                "at_lo": pt.get("at_lo", False), "at_hi": pt.get("at_hi", False),
                "ci_Jc": ci_jc, "ci_pp": ci_pp,
                "contains_2": bool(contains_2), "flagged": flagged,
                "n_dropped": n_dropped, "n_boundary": n_boundary,
                "unreach": unreach,
                "n_boot_used": int(np.sum(np.isfinite(jc_b))),
            })
            print(f"[fit] xi={xi} target={tg}: J_c={_fmt(pt.get('Jc'))} "
                  f"[{_fmt(ci_jc[0])},{_fmt(ci_jc[1])}] p'={_fmt(pt.get('p'))} "
                  f"[{_fmt(ci_pp[0])},{_fmt(ci_pp[1])}] R2={_fmt(pt.get('R2'))} "
                  f"n={pt.get('n', 0)} dropped_reps={n_dropped}"
                  f"{' FLAGGED(n=4)' if flagged else ''}", flush=True)

    # ------------------------------------------------------------------
    # jc_free_fits.csv
    # ------------------------------------------------------------------
    csv_path = RESULTS_DIR / "jc_free_fits.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("xi,target,n_points,A,p_prime,J_c,ci_Jc_lo,ci_Jc_hi,"
                "ci_pprime_lo,ci_pprime_hi,R2,contains_2,flagged\n")
        for r in rows:
            f.write(f"{r['xi']},{r['target']},{r['n_points']},{_fmt(r['A'])},"
                    f"{_fmt(r['p_prime'])},{_fmt(r['J_c'])},"
                    f"{_fmt(r['ci_Jc'][0])},{_fmt(r['ci_Jc'][1])},"
                    f"{_fmt(r['ci_pp'][0])},{_fmt(r['ci_pp'][1])},"
                    f"{_fmt(r['R2'])},{r['contains_2']},{r['flagged']}\n")

    # ------------------------------------------------------------------
    # B4_SUMMARY.md
    # ------------------------------------------------------------------
    comp = [r for r in rows if not r["flagged"]]        # 8 comparable cells
    flg = [r for r in rows if r["flagged"]]
    n_contain = sum(r["contains_2"] for r in comp)
    jc_pts = [r["J_c"] for r in comp if np.isfinite(r["J_c"])]
    consistent = n_contain == len(comp)

    if consistent:
        verdict = ("**VERDICT: J_c consistent with 2 across cells — the (J-T) "
                   "bifurcation-distance language is permitted.**")
    else:
        verdict = ("**VERDICT: J_c departs from 2 — drop (J-T) language; report "
                   "the empirical J-dependence / fitted J_c as a measured "
                   "quantity.**")

    lines = []
    lines.append("# Free-J_c fit h0*(J) = A (J - J_c)^p': SUMMARY\n")
    lines.append(verdict + "\n")
    lines.append(f"- {n_contain} of {len(comp)} comparable cells (n >= 5 J points; "
                 f"the n=4 cell xi=0.65/target=0.10 excluded) have 2 inside the "
                 f"95% bootstrap CI of J_c.")
    if jc_pts:
        lines.append(f"- J_c point estimates (comparable cells) range from "
                     f"{min(jc_pts):.3f} to {max(jc_pts):.3f}.")
    lines.append("")
    lines.append("Method: profile grid search over J_c in "
                 f"[{JC_LO}, {JC_HI}] (step {JC_STEP}; J_c <= min(J) - 0.05); at each "
                 "J_c, OLS of log(h0*) on log(J - J_c); J_c chosen to minimize SSE. "
                 f"Bootstrap: {N_BOOT} path-level resamples reusing the amplitude-requirement fit's stored "
                 "per-path collapse matrices, same resample across the h0 axis within "
                 "each (J, xi) file and shared across J within a xi block (seed "
                 f"0x{BOOT_SEED_BASE:X} + xi_idx, identical to h0_scaling.py); "
                 "inversion of h0* mirrors h0_scaling.invert_h0_star exactly "
                 "(imported). Replicates: unreachable targets drop that J from the "
                 "replicate's fit; replicates with < 4 valid J points are dropped. "
                 "95% percentile CIs.\n")

    lines.append("## Free-J_c fits (comparable cells, n = 6 J points)\n")
    hdr = ("| xi | target | n | A | p' | 95% CI p' | J_c | 95% CI J_c | R^2 | "
           "CI contains 2? |")
    sep = "|----|--------|---|---|----|-----------|-----|------------|-----|----|"
    lines.append(hdr)
    lines.append(sep)
    for r in comp:
        lines.append(
            f"| {r['xi']} | {r['target']} | {r['n_points']} | {_fmt(r['A'])} | "
            f"{_fmt(r['p_prime'])} | [{_fmt(r['ci_pp'][0])}, {_fmt(r['ci_pp'][1])}] | "
            f"{_fmt(r['J_c'])} | [{_fmt(r['ci_Jc'][0])}, {_fmt(r['ci_Jc'][1])}] | "
            f"{_fmt(r['R2'])} | {'YES' if r['contains_2'] else 'NO'} |")
    lines.append("")
    lines.append("## Flagged cell (n=4, not comparable — excluded from the count)\n")
    lines.append(hdr)
    lines.append(sep)
    for r in flg:
        lines.append(
            f"| {r['xi']} | {r['target']} | {r['n_points']} | {_fmt(r['A'])} | "
            f"{_fmt(r['p_prime'])} | [{_fmt(r['ci_pp'][0])}, {_fmt(r['ci_pp'][1])}] | "
            f"{_fmt(r['J_c'])} | [{_fmt(r['ci_Jc'][0])}, {_fmt(r['ci_Jc'][1])}] | "
            f"{_fmt(r['R2'])} | {'YES' if r['contains_2'] else 'NO'} |")
    lines.append("")

    lines.append("## Pattern across cells\n")
    for xi in XI_GRID:
        sub = [r for r in rows if r["xi"] == xi]
        s = ", ".join(f"target {r['target']}: J_c={_fmt(r['J_c'])}"
                      + (" (n=4, flagged)" if r["flagged"] else "")
                      for r in sub)
        lines.append(f"- xi={xi}: {s}")
    for tg in TARGETS:
        sub = [r for r in rows if r["target"] == tg]
        s = ", ".join(f"xi {r['xi']}: J_c={_fmt(r['J_c'])}"
                      + (" (n=4, flagged)" if r["flagged"] else "")
                      for r in sub)
        lines.append(f"- target={tg}: {s}")
    lines.append("")

    lines.append("## Bootstrap pathologies\n")
    lines.append("| xi | target | reps used | reps dropped (<4 valid J) | "
                 "reps with J_c at grid boundary | unreachable-target counts "
                 "per J (of 1000 reps) |")
    lines.append("|----|--------|-----------|----------------------------|"
                 "-------------------------------|------------------------|")
    for r in rows:
        un = "; ".join(f"J={J}: {c}" for J, c in r["unreach"].items() if c > 0)
        lines.append(f"| {r['xi']} | {r['target']} | {r['n_boot_used']} | "
                     f"{r['n_dropped']} | {r['n_boundary']} | {un or 'none'} |")
    lines.append("")
    bnd = [r for r in rows if r.get("at_lo") or r.get("at_hi")]
    if bnd:
        lines.append("Point-estimate J_c pinned at a grid boundary in: "
                     + ", ".join(f"(xi={r['xi']}, target={r['target']}, "
                                 f"{'lower' if r['at_lo'] else 'upper'})"
                                 for r in bnd) + ".")
    else:
        lines.append("No point-estimate J_c pinned at a grid boundary.")
    lines.append("")

    with open(RESULTS_DIR / "B4_SUMMARY.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"[done] wrote {csv_path} and B4_SUMMARY.md "
          f"({time.time()-t_start:.1f}s total)", flush=True)


if __name__ == "__main__":
    main()
