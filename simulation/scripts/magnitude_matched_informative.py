"""
INFORMATIVE magnitude-matched control (stress variant).

Pre-registered in:
  simulation/prereg/PREREG_magnitude_matched_informative_2026-08-21.md

The pre-existing control (magnitude_matched_control.py, alpha=2.0) compares two
all-zero arms and cannot discriminate "magnitude decides" from "form decides".
This script:

  STEP A. SEARCH for an informative operating point: at the headline cell
          (J=5, mu=100, n>=2000), reduce the stress coefficient alpha (grid,
          ascending) with base scale b=2 until P(stress-form collapse) lands
          strictly inside (0.10, 0.50); pick the smallest-index (smallest) alpha
          in band. If no alpha lands in band, lower the base scale b (grid,
          descending) and rescan; pick the smallest-index (b, alpha) in band.

  STEP B. MATCHED CONTROL at the chosen (b*, alpha*): for the headline cell and
          two corroborating cells (J=4, 4.5), compare the stress form against a
          pure constant field equal to the stress form's FULL-TRAJECTORY
          delivered mean (time-average of h over all steps, averaged across
          runs). Matched seeds (identical noise draws) -> paired difference,
          10,000-rep paired bootstrap 95% CI.

Standard constants: dt=0.05, 20000 steps, T=2, xi=0.5, m0=0, W0=100,
collapse = W<10 for 200 consecutive steps, mu=100.
"""

from __future__ import annotations

from pathlib import Path
import time

import numpy as np
import pandas as pd

# ---- Fixed integrator constants (identical to minimal_model.py) ------------
T_TEMP = 2.0
XI = 0.5
DT = 0.05
N_STEPS = 20_000
INITIAL_M = 0.0
INITIAL_WEALTH = 100.0
COLLAPSE_WEALTH = 10.0
COLLAPSE_STREAK = 200
EPS = 1e-6
MU = 100.0

# ---- Pre-registered search grids -------------------------------------------
HEADLINE_J = 5.0
ALPHA_GRID = (0.05, 0.10, 0.15, 0.20, 0.30, 0.50)   # ascending
BASE_GRID = (2.0, 1.5, 1.0, 0.5, 0.25)              # descending; 2.0 = primary
BAND_LO, BAND_HI = 0.10, 0.50                        # strictly inside

# ---- Matched-control design ------------------------------------------------
CELLS_J = (5.0, 4.0, 4.5)          # headline first, then corroborating
N_SEEDS = 2000
SEED_BASE = 0xE1_1B01              # distinct from existing 0xE1_0001 control
N_BOOT = 10_000
BOOT_SEED = 0xB007

# ---- Decision-rule thresholds (pre-registered) -----------------------------
REPRODUCE_ABS = 0.03              # |diff| <= this AND CI contains 0  -> reproduced
REFUTE_ABS = 0.05                # |diff| >  this AND CI excludes 0   -> refuted

OUT = Path(__file__).resolve().parents[1] / "results" / "magnitude_matched"


def _seed_for_J(J: float) -> int:
    """Deterministic per-cell seed keyed to J's index in a canonical J list."""
    canonical = [5.0, 4.0, 4.5, 4.75, 3.5, 3.0]
    j_idx = canonical.index(J) if J in canonical else int(round(J * 10))
    return SEED_BASE + j_idx * 1_000_003


def run_stress(J: float, b: float, alpha: float, seed: int,
               n_seeds: int = N_SEEDS):
    """Stress-form field h = b*(W/500) + alpha*max(0,-m)*J.

    Returns (collapsed[bool array], h_bar_full) where h_bar_full is the
    full-trajectory delivered mean: time-average of h over all steps, averaged
    across all runs.
    """
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M)
    wealth = np.full(n_seeds, INITIAL_WEALTH)
    collapsed = np.zeros(n_seeds, dtype=bool)
    streak = np.zeros(n_seeds, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)

    full_sum = 0.0
    full_n = 0

    for _ in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        base = b * (wealth / 500.0)
        stress = np.maximum(0.0, -m_c)
        h = base + alpha * stress * J

        full_sum += float(h.sum())
        full_n += n_seeds

        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)

        employment = 0.5 + 0.5 * m
        income = MU * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        is_low = wealth < COLLAPSE_WEALTH
        streak = np.where(is_low, streak + 1, 0)
        collapsed |= streak >= COLLAPSE_STREAK

    h_bar_full = full_sum / full_n
    return collapsed, h_bar_full


def run_constant(J: float, h_const: float, seed: int, n_seeds: int = N_SEEDS):
    """Pure constant field h_const: no wealth feedback, no stress term.

    Same seed and same RNG call pattern as run_stress -> identical noise draws
    (paired comparison).
    """
    rng = np.random.default_rng(seed)
    m = np.full(n_seeds, INITIAL_M)
    wealth = np.full(n_seeds, INITIAL_WEALTH)
    collapsed = np.zeros(n_seeds, dtype=bool)
    streak = np.zeros(n_seeds, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)

    for _ in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        f = -m_c + np.tanh((J * m_c + h_const) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n_seeds)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)

        employment = 0.5 + 0.5 * m
        income = MU * employment
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)

        is_low = wealth < COLLAPSE_WEALTH
        streak = np.where(is_low, streak + 1, 0)
        collapsed |= streak >= COLLAPSE_STREAK

    return collapsed


def paired_boot_diff_ci(a: np.ndarray, b: np.ndarray,
                        rng: np.random.Generator, n_boot: int = N_BOOT):
    """Paired (matched-seed) percentile 95% CI on mean(a) - mean(b)."""
    n = len(a)
    d = a.astype(float) - b.astype(float)
    idx = rng.integers(0, n, size=(n_boot, n))
    diffs = d[idx].mean(axis=1)
    return float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def search_operating_point():
    """Pre-registered search at the headline cell (J=5, mu=100).

    Returns (b_star, alpha_star, scan_rows). scan_rows records every (b, alpha)
    tried with its P(stress collapse), in scan order.
    """
    seed = _seed_for_J(HEADLINE_J)
    scan_rows = []

    for b in BASE_GRID:
        for alpha in ALPHA_GRID:
            t0 = time.time()
            coll, h_bar = run_stress(HEADLINE_J, b, alpha, seed)
            p = float(coll.mean())
            in_band = (BAND_LO < p < BAND_HI)
            scan_rows.append({
                "b": b, "alpha": alpha, "J": HEADLINE_J, "mu": MU,
                "p_stress": p, "h_bar_full": h_bar, "in_band": in_band,
                "n": N_SEEDS,
            })
            print(f"  [search] b={b:.2f} alpha={alpha:.2f}  "
                  f"P(stress)={p:.4f}  h_bar={h_bar:.4f}  "
                  f"in_band={in_band}  ({time.time()-t0:.1f}s)")
            if in_band:
                # smallest-index alpha (ascending) at smallest-index b
                # (descending) is the FIRST in-band hit encountered.
                return b, alpha, scan_rows
        # No alpha in band at this b; per prereg, lower b and rescan.
    return None, None, scan_rows


def full_alpha_scan():
    """Full headline-cell (b = 2.0, J = 5, mu = 100) scan over ALPHA_GRID,
    deposited as operating_point_alpha_scan.csv. Unlike search_operating_point(),
    which early-returns at the first in-band alpha, this records P(stress) at
    every alpha so the reported monotone decrease of P with alpha (§S5.6) is
    reproducible from a deposited file rather than asserted."""
    seed = _seed_for_J(HEADLINE_J)
    b = BASE_GRID[0]                       # 2.0, primary base
    rows = []
    for alpha in ALPHA_GRID:
        coll, h_bar = run_stress(HEADLINE_J, b, alpha, seed)
        p = float(coll.mean())
        rows.append({"b": b, "alpha": alpha, "J": HEADLINE_J, "mu": MU,
                     "p_stress": p, "h_bar_full": h_bar,
                     "in_band": bool(BAND_LO < p < BAND_HI), "n": N_SEEDS})
        print(f"  [alpha-scan] alpha={alpha:.2f}  P(stress)={p:.4f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "operating_point_alpha_scan.csv", index=False)
    return df


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Output dir: {OUT}")
    print("\n=== STEP 0: full headline-cell alpha scan (deposited) ===")
    full_alpha_scan()
    print("\n=== STEP A: search for informative operating point (J=5, mu=100) ===")
    b_star, alpha_star, scan_rows = search_operating_point()

    # Search provenance (every (b, alpha) tried) is printed to stdout above and
    # the chosen point is recorded per-cell in informative_control.csv; we do
    # not emit a separate scan file.
    if b_star is None:
        print("\n[SEARCH FAILED] No (b, alpha) on the pre-registered grids gave "
              "P(stress) strictly inside (0.10, 0.50). No verdict on claim 3.")
        # Still write an (empty-verdict) control file for provenance.
        pd.DataFrame(columns=[
            "J", "mu", "b", "alpha", "h_bar_full", "p_stress", "p_const",
            "diff", "diff_lo", "diff_hi", "n", "verdict",
        ]).to_csv(OUT / "informative_control.csv", index=False)
        return

    print(f"\nChosen operating point: b*={b_star}, alpha*={alpha_star}")

    print("\n=== STEP B: matched constant-field control at chosen operating point ===")
    boot_rng = np.random.default_rng(BOOT_SEED)
    rows = []
    for J in CELLS_J:
        seed = _seed_for_J(J)
        t0 = time.time()
        coll_stress, h_bar = run_stress(J, b_star, alpha_star, seed)
        coll_const = run_constant(J, h_bar, seed)   # matched seed -> paired
        p_s = float(coll_stress.mean())
        p_c = float(coll_const.mean())
        diff = p_s - p_c
        lo, hi = paired_boot_diff_ci(coll_stress, coll_const, boot_rng)

        # Pre-registered decision rule
        ci_contains_0 = (lo <= 0.0 <= hi)
        if ci_contains_0 and abs(diff) <= REPRODUCE_ABS:
            verdict = "REPRODUCED"
        elif (not ci_contains_0) and abs(diff) > REFUTE_ABS:
            verdict = "REFUTED"
        else:
            verdict = "INCONCLUSIVE"

        rows.append({
            "J": J, "mu": MU, "b": b_star, "alpha": alpha_star,
            "h_bar_full": h_bar, "p_stress": p_s, "p_const": p_c,
            "diff": diff, "diff_lo": lo, "diff_hi": hi,
            "n": N_SEEDS, "verdict": verdict,
        })
        print(f"J={J:.1f}  h_bar={h_bar:.4f}  P(stress)={p_s:.4f}  "
              f"P(const)={p_c:.4f}  diff={diff:+.4f}  "
              f"CI=[{lo:+.4f}, {hi:+.4f}]  -> {verdict}  "
              f"({time.time()-t0:.1f}s)")

    ctrl_df = pd.DataFrame(rows)
    ctrl_df.to_csv(OUT / "informative_control.csv", index=False)
    print(f"\nWrote {OUT / 'informative_control.csv'}")

    # Overall verdict (headline cell decisive)
    headline = ctrl_df[ctrl_df.J == HEADLINE_J].iloc[0]
    any_refuted = (ctrl_df.verdict == "REFUTED").any()
    if headline.verdict == "REPRODUCED" and not any_refuted:
        overall = "REPRODUCED"
    elif headline.verdict == "REFUTED":
        overall = "REFUTED"
    else:
        overall = "INCONCLUSIVE"
    print(f"\nOVERALL VERDICT (headline J=5 decisive): {overall}")


if __name__ == "__main__":
    main()
