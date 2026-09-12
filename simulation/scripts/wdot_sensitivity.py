"""
wdot_sensitivity.py — Eq (2) sensitivity analysis.

Perturbs the consumption side of the wealth equation (Eq 2 of the manuscript),

    wdot = mult * (0.5 + 0.5*m) - [ cons_const + cons_coeff * (wealth/500) ],
    wealth >= 0 (floor),

around the baseline (cons_const, cons_coeff) = (30, 10), plus one purely
linear alternative consumption = 0.08 * wealth (chosen to match the baseline
total consumption of 40 at W = 500; no constant term). Six configurations:

    baseline   (30, 10)
    const-30%  (21, 10)      const+30%  (39, 10)
    coeff-30%  (30,  7)      coeff+30%  (30, 13)
    linear     consumption = 0.08 * W

For each variant:
  (a) P(collapse | J) at h0 = 2.0, xi = 0.5, mult = 100, J in {2.5..5.0},
      n_seeds = 2000, bootstrap 95% CIs, and the (seed- and bootstrap-paired)
      shift vs baseline per J;
  (b) the amplitude-requirement (h0_scaling) exponent measurement for target 0.30, xi = 0.5 only:
      h0*(J) by coarse+refined h0 grid inversion, OLS fit
      log(h0*) = a + p log(J), bootstrap CI, and delta vs baseline p.

Seeds are IDENTICAL across variants (seed = 0xB0_5CA1E + j_idx*1_000_003
+ 1*9973, i.e. the h0_scaling xi=0.5 seeds), and the same bootstrap resample
matrix is shared, so per-J shifts and delta-p are paired comparisons.

Usage
-----
    python wdot_sensitivity.py run --variant-idx K   # K in 0..5
    python wdot_sensitivity.py analyze
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from h0_scaling import (  # noqa: E402
    BOOT_SEED_BASE,
    J_GRID,
    N_BOOT,
    N_SEEDS_B,
    SEED_BASE,
    bootstrap_count_matrix,
    cell_path,
    measure_cell,
    monotonicity_violations,
    ols_loglog,
    pct_ci,
    _fmt,
)

RESULTS_DIR = SCRIPT_DIR.parent / "results" / "wdot_sensitivity"
XI = 0.5
XI_IDX = 1                 # index of xi=0.5 in h0_scaling.XI_GRID -> same seeds as the amplitude-requirement fit
H0_EVAL = 2.0              # manuscript's hard-coded stabilizer amplitude
TARGET = 0.30              # exponent re-measurement target

VARIANTS = (
    ("baseline",  dict(cons_const=30.0, cons_coeff=10.0, cons_form="affine")),
    ("const-30",  dict(cons_const=21.0, cons_coeff=10.0, cons_form="affine")),
    ("const+30",  dict(cons_const=39.0, cons_coeff=10.0, cons_form="affine")),
    ("coeff-30",  dict(cons_const=30.0, cons_coeff=7.0,  cons_form="affine")),
    ("coeff+30",  dict(cons_const=30.0, cons_coeff=13.0, cons_form="affine")),
    ("linear",    dict(cons_form="linear", cons_lin=0.08)),
)


def variant_dir(name: str) -> Path:
    return RESULTS_DIR / f"raw_{name}"


def cmd_run(variant_idx: int) -> None:
    name, kwargs = VARIANTS[variant_idx]
    vdir = variant_dir(name)
    vdir.mkdir(parents=True, exist_ok=True)
    print(f"[run] variant={name} kwargs={kwargs} xi={XI} (seeds match amplitude-requirement fit xi_idx={XI_IDX})",
          flush=True)
    for j_idx, J in enumerate(J_GRID):
        measure_cell(J, j_idx, XI, XI_IDX, vdir, n_seeds=N_SEEDS_B,
                     seed_base=SEED_BASE, **kwargs)
    print(f"[run] variant={name} COMPLETE", flush=True)


def _from_h0_row(h0: np.ndarray, value: float) -> int:
    idx = np.where(np.isclose(h0, value, atol=1e-12))[0]
    if idx.size == 0:
        raise RuntimeError(f"h0={value} not found in stored grid {h0}")
    return int(idx[0])


def cmd_analyze() -> None:
    from h0_scaling import invert_h0_star

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(BOOT_SEED_BASE + XI_IDX)   # same as amplitude-requirement fit xi=0.5
    B = bootstrap_count_matrix(N_SEEDS_B, N_BOOT, rng)     # shared across variants & J

    J_arr = np.asarray(J_GRID)
    per_variant = {}
    mono_report = []
    na_report = []

    for name, kwargs in VARIANTS:
        vdir = variant_dir(name)
        P_pt = np.full(len(J_GRID), np.nan)
        P_bt = np.full((len(J_GRID), N_BOOT), np.nan)
        h0s_pt = np.full(len(J_GRID), np.nan)
        h0s_bt = np.full((len(J_GRID), N_BOOT), np.nan)

        for j_idx, J in enumerate(J_GRID):
            path = cell_path(vdir, J, XI)
            if not path.exists():
                raise FileNotFoundError(f"missing {path} — run `run --variant-idx` first")
            with np.load(path) as z:
                h0 = z["h0"].copy()
                C = z["collapsed"].copy()
            viol = monotonicity_violations(h0, C)
            if viol:
                mono_report.append((name, J, viol))
            curve_pt = C.mean(axis=1)
            curve_bt = (C.astype(np.float64) @ B) / N_SEEDS_B

            row = _from_h0_row(h0, H0_EVAL)
            P_pt[j_idx] = curve_pt[row]
            P_bt[j_idx] = curve_bt[row]

            h0s = invert_h0_star(h0, curve_pt, TARGET)
            h0s_pt[j_idx] = h0s
            if not np.isfinite(h0s):
                na_report.append((name, J, TARGET,
                                  f"P range [{curve_pt.min():.3f},{curve_pt.max():.3f}]"))
            for r in range(N_BOOT):
                h0s_bt[j_idx, r] = invert_h0_star(h0, curve_bt[:, r], TARGET)

        p, a, r2, n_used = ols_loglog(J_arr, h0s_pt, 0.0)
        p_boot = np.array([ols_loglog(J_arr, h0s_bt[:, r], 0.0)[0] for r in range(N_BOOT)])
        per_variant[name] = dict(P_pt=P_pt, P_bt=P_bt, h0s_pt=h0s_pt, h0s_bt=h0s_bt,
                                 p=p, p_boot=p_boot, r2=r2, n_used=n_used)
        ci = pct_ci(p_boot)
        print(f"[fit] {name}: p={p:.3f} [{ci[0]:.3f},{ci[1]:.3f}] R2={r2:.4f} "
              f"P(h0=2)={np.array2string(P_pt, precision=3)}", flush=True)

    base = per_variant["baseline"]

    # ---- pcollapse_by_variant.csv ----
    with open(RESULTS_DIR / "pcollapse_by_variant.csv", "w", encoding="utf-8") as f:
        f.write("variant,J,P_collapse,ci_lo,ci_hi,delta_vs_baseline,delta_ci_lo,delta_ci_hi\n")
        for name, _ in VARIANTS:
            v = per_variant[name]
            for j_idx, J in enumerate(J_GRID):
                lo, hi = pct_ci(v["P_bt"][j_idx])
                d = v["P_pt"][j_idx] - base["P_pt"][j_idx]
                dlo, dhi = pct_ci(v["P_bt"][j_idx] - base["P_bt"][j_idx])  # paired
                f.write(f"{name},{J},{_fmt(v['P_pt'][j_idx])},{_fmt(lo)},{_fmt(hi)},"
                        f"{_fmt(d)},{_fmt(dlo)},{_fmt(dhi)}\n")

    # ---- exponent_by_variant.csv ----
    with open(RESULTS_DIR / "exponent_by_variant.csv", "w", encoding="utf-8") as f:
        f.write("variant,target,p,ci_lo,ci_hi,R2,n_J_used,"
                "delta_p_vs_baseline,delta_ci_lo,delta_ci_hi\n")
        for name, _ in VARIANTS:
            v = per_variant[name]
            lo, hi = pct_ci(v["p_boot"])
            dp = v["p"] - base["p"]
            dlo, dhi = pct_ci(v["p_boot"] - base["p_boot"])  # paired
            f.write(f"{name},{TARGET},{_fmt(v['p'])},{_fmt(lo)},{_fmt(hi)},"
                    f"{_fmt(v['r2'])},{v['n_used']},{_fmt(dp)},{_fmt(dlo)},{_fmt(dhi)}\n")

    _write_summary(per_variant, base, mono_report, na_report)
    _write_methods_draft(per_variant, base)
    print(f"[analyze] wrote outputs to {RESULTS_DIR}", flush=True)


def _write_summary(per_variant, base, mono_report, na_report) -> None:
    lines = []
    lines.append("# Eq (2) sensitivity: SUMMARY\n")
    lines.append(f"Setup: h0 = {H0_EVAL}, xi = {XI}, mult = 100, J in "
                 f"{{{', '.join(str(j) for j in J_GRID)}}}, n_seeds = {N_SEEDS_B}, "
                 f"seeds and bootstrap resamples identical across variants (paired "
                 f"comparisons). Exponent re-measured at target P = {TARGET}.\n")

    lines.append("## P(collapse | J) at h0 = 2.0 by variant\n")
    header = "| variant | " + " | ".join(f"J={J}" for J in J_GRID) + " |"
    lines.append(header)
    lines.append("|" + "---|" * (len(J_GRID) + 1))
    for name in per_variant:
        v = per_variant[name]
        lines.append(f"| {name} | " + " | ".join(f"{x:.3f}" for x in v["P_pt"]) + " |")

    max_abs_d = 0.0
    max_where = ""
    for name in per_variant:
        if name == "baseline":
            continue
        d = per_variant[name]["P_pt"] - base["P_pt"]
        i = int(np.nanargmax(np.abs(d)))
        if abs(d[i]) > max_abs_d:
            max_abs_d = abs(d[i])
            max_where = f"{name} at J={J_GRID[i]} (delta = {d[i]:+.3f})"
    lines.append(f"\nMaximum |delta P| vs baseline across all variants and J: "
                 f"**{max_abs_d:.3f}** ({max_where}).")

    lines.append("\n## Exponent p (target 0.30, xi = 0.5) by variant\n")
    lines.append("| variant | p | 95% CI | R^2 | delta p vs baseline (paired 95% CI) |")
    lines.append("|---------|---|--------|-----|--------------------------------------|")
    ps = []
    for name in per_variant:
        v = per_variant[name]
        lo, hi = pct_ci(v["p_boot"])
        dlo, dhi = pct_ci(v["p_boot"] - base["p_boot"])
        ps.append(v["p"])
        lines.append(f"| {name} | {_fmt(v['p'])} | [{_fmt(lo)}, {_fmt(hi)}] | {_fmt(v['r2'])} "
                     f"| {_fmt(v['p'] - base['p'])} [{_fmt(dlo)}, {_fmt(dhi)}] |")
    ps = [p for p in ps if np.isfinite(p)]
    lines.append(f"\nRange of p across the six variants: **[{min(ps):.3f}, {max(ps):.3f}]** "
                 f"(spread {max(ps) - min(ps):.3f}).")

    lines.append("\n## Range over which conclusions hold\n")
    lines.append(f"- +/-30% perturbations of the consumption constant and coefficient, and a "
                 f"structurally different purely-linear consumption law, shift P(collapse) at "
                 f"the manuscript's operating point (h0=2, xi=0.5) by at most {max_abs_d:.3f} "
                 f"in absolute probability.")
    lines.append(f"- The stabilizer-scaling exponent stays within [{min(ps):.3f}, {max(ps):.3f}] "
                 f"across all six consumption specifications.")

    lines.append("\n## Monotonicity / NA notes\n")
    if mono_report:
        for name, J, viol in mono_report:
            for lo_h, hi_h, b, c in viol:
                lines.append(f"- {name}, J={J}: P rose from h0={lo_h:.3f} to {hi_h:.3f} "
                             f"(b={b}, c={c})")
    else:
        lines.append("- No monotonicity violations beyond Monte-Carlo noise.")
    if na_report:
        for name, J, tg, why in na_report:
            lines.append(f"- NA: {name}, J={J}, target={tg} unreachable ({why})")
    else:
        lines.append("- No unreachable targets: h0*(J; 0.30) invertible for every variant and J.")
    lines.append("")
    with open(RESULTS_DIR / "SUMMARY.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def _write_methods_draft(per_variant, base) -> None:
    ps = [v["p"] for v in per_variant.values() if np.isfinite(v["p"])]
    max_abs_d = max(
        float(np.nanmax(np.abs(per_variant[name]["P_pt"] - base["P_pt"])))
        for name in per_variant if name != "baseline"
    )
    text = f"""# Draft Methods paragraph — status of Eq (2)

*(Rendered paragraph; the manuscript's §S5.5 text is maintained separately.)*

The wealth dynamics of Eq (2), wdot = mult * (1 + m)/2 - [30 + 10 (W/500)]
with the reflecting floor W >= 0, is a modelling choice: the affine
consumption term is neither derived from microeconomic first principles nor
empirically calibrated. It encodes only two qualitative requirements — a
fixed subsistence cost and a weakly wealth-dependent component — and the
W >= 0 floor prevents unbounded debt. To establish that our conclusions do
not depend on this specific functional form, we repeated the analysis under
+/-30% perturbations of both consumption parameters ((21, 10), (39, 10),
(30, 7), (30, 13)) and under a structurally different, purely linear
consumption law C(W) = 0.08 W (matched to the baseline total consumption of
40 at W = 500, with no subsistence constant). Across all five alternatives,
collapse probabilities at the operating point (h0 = 2, xi = 0.5, mult = 100)
shift by at most {max_abs_d:.2f} in absolute probability over the supercritical
range J in [2.5, 5], and the stabilizer-scaling exponent p (defined by
log h0* = a + p log J at the P = 0.30 iso-probability level) remains within
[{min(ps):.2f}, {max(ps):.2f}] (baseline p = {base['p']:.2f}). Full sensitivity tables are
provided in simulation/results/wdot_sensitivity/.
"""
    with open(RESULTS_DIR / "methods_eq2_draft.md", "w", encoding="utf-8") as f:
        f.write(text)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p_run = sub.add_parser("run")
    p_run.add_argument("--variant-idx", type=int, required=True,
                       choices=range(len(VARIANTS)))
    sub.add_parser("analyze")
    args = ap.parse_args()

    if args.cmd == "run":
        cmd_run(args.variant_idx)
    elif args.cmd == "analyze":
        cmd_analyze()


if __name__ == "__main__":
    main()
