"""Phase-A statistics for the S5.6 ordering test (bootstrap CIs, domain split,
in-domain coupling arm).

Extends `delivered_mean_ordering.py`: same field forms, integrator, seeds, and
five summaries, but (A1) attaches a matched-seed paired bootstrap SE and 95% CI
to every summary and to P(collapse), and evaluates each adjacent-arm ordering
step *statistically* (is the amplitude gap resolved? is the P gap resolved? is a
violation resolved?); (A2) tags each arm in/out of the selection formula's domain
(lambda > 0 at the measured selection-window field) and runs the ordering test
over in-domain arms (PRIMARY) and all arms (SECONDARY); (A3) scans the coupling
coefficient for the strongest arm whose origin is still linearly unstable and
adds it as a sixth arm.

Bootstrap convention matched to `magnitude_matched_informative.paired_boot_diff_ci`:
N_BOOT = 10,000 seed-index resamples, BOOT_SEED = 0xB007, 2.5/97.5 percentile CIs.
Paired: arms at a given J share the same seed, so a single resampled index matrix
applied to two arms is a matched-seed paired bootstrap of their difference.

Outputs (simulation/results/magnitude_matched/):
  delivered_mean_ordering_ci.csv   one row per (J, arm): point/se/ci/n_defined per summary + p_collapse, in_domain
  ordering_verdicts.csv            one row per (J, summary, adjacent pair) with resolution + violation flags
  coupling_domain_scan.csv         (alpha, J, win_mean, lam, in_domain)
  coupling_in_domain_arm.csv       the sixth arm's five summaries + CIs
Post-hoc diagnostic (not part of the two-arm pre-registered plan).
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

from delivered_mean_ordering import (
    T_TEMP, XI, DT, N_STEPS, INITIAL_M, INITIAL_WEALTH,
    COLLAPSE_WEALTH, COLLAPSE_STREAK, EPS, MU, N_SEEDS, SEED_BASE,
    CELLS_J, SELWIN_STEPS, ARMS, OUT, field,
)

N_BOOT = 10_000
BOOT_SEED = 0xB007
SUMMARIES = ["full_mean", "win_mean", "cond_m_neg", "cond_m_neg_collapsers", "amp_over_lam"]


def seed_for_J(J):
    return SEED_BASE + int(round(J * 10)) * 1_000_003


def run_arm_ps(arm, alpha, J, seed, n=N_SEEDS):
    """Identical dynamics to delivered_mean_ordering.run_arm, but returns
    per-seed arrays so every summary can be bootstrapped."""
    rng = np.random.default_rng(seed)
    m = np.full(n, INITIAL_M)
    wealth = np.full(n, INITIAL_WEALTH)
    collapsed = np.zeros(n, dtype=bool)
    streak = np.zeros(n, dtype=np.int32)
    sqrt_dt = np.sqrt(DT)
    full_ps = np.zeros(n)
    win_ps = np.zeros(n)
    neg_sum = np.zeros(n)
    neg_cnt = np.zeros(n, dtype=np.int64)
    for t in range(N_STEPS):
        m_c = np.clip(m, -1 + EPS, 1 - EPS)
        h = field(arm, alpha, wealth, m_c, J)
        full_ps += h
        if t < SELWIN_STEPS:
            win_ps += h
            neg = m_c < 0.0
            neg_sum += np.where(neg, h, 0.0)
            neg_cnt += neg
        f = -m_c + np.tanh((J * m_c + h) / T_TEMP)
        g = XI * np.sqrt(1.0 - m_c * m_c)
        dW = sqrt_dt * rng.standard_normal(n)
        m = np.clip(m_c + f * DT + g * dW, -1 + EPS, 1 - EPS)
        income = MU * (0.5 + 0.5 * m)
        consumption = 30.0 + 10.0 * (wealth / 500.0)
        wealth = np.maximum(0.0, wealth + (income - consumption) * DT)
        streak = np.where(wealth < COLLAPSE_WEALTH, streak + 1, 0)
        collapsed |= streak >= COLLAPSE_STREAK
    return dict(full_mean_ps=full_ps / N_STEPS,
                win_mean_ps=win_ps / SELWIN_STEPS,
                neg_sum=neg_sum, neg_cnt=neg_cnt,
                collapsed=collapsed, J=J)


def lam_of(win_mean, J):
    return -1.0 + (J / T_TEMP) / np.cosh(win_mean / T_TEMP) ** 2


# --- summary evaluated on a resample index matrix (n_boot, n) or the identity ---
def summary_on(name, d, idx=None):
    fm, wm = d["full_mean_ps"], d["win_mean_ps"]
    ns, nc, coll, J = d["neg_sum"], d["neg_cnt"], d["collapsed"], d["J"]
    if idx is None:  # point estimate over all seeds
        if name == "full_mean":
            return float(fm.mean())
        if name == "win_mean":
            return float(wm.mean())
        if name == "cond_m_neg":
            return float(ns.sum() / nc.sum()) if nc.sum() else np.nan
        if name == "cond_m_neg_collapsers":
            s, c = ns[coll].sum(), nc[coll].sum()
            return float(s / c) if c else np.nan
        if name == "amp_over_lam":
            w = float(wm.mean()); L = lam_of(w, J)
            return float(w / L) if L > 0 else np.nan
        raise ValueError(name)
    # bootstrap: vector over resamples
    if name == "full_mean":
        return fm[idx].mean(axis=1)
    if name == "win_mean":
        return wm[idx].mean(axis=1)
    if name == "cond_m_neg":
        num = ns[idx].sum(axis=1); den = nc[idx].sum(axis=1)
        return np.where(den > 0, num / np.where(den == 0, 1, den), np.nan)
    if name == "cond_m_neg_collapsers":
        nsc = ns * coll; ncc = nc * coll
        num = nsc[idx].sum(axis=1); den = ncc[idx].sum(axis=1)
        return np.where(den > 0, num / np.where(den == 0, 1, den), np.nan)
    if name == "amp_over_lam":
        w = wm[idx].mean(axis=1); L = lam_of(w, J)
        return np.where(L > 0, w / np.where(L > 0, L, np.nan), np.nan)
    raise ValueError(name)


def n_defined(name, d):
    nc, coll = d["neg_cnt"], d["collapsed"]
    if name in ("full_mean", "win_mean"):
        return int(d["full_mean_ps"].size)
    if name == "cond_m_neg":
        return int((nc > 0).sum())
    if name == "cond_m_neg_collapsers":
        return int(((nc > 0) & coll).sum())
    if name == "amp_over_lam":
        w = float(d["win_mean_ps"].mean())
        return int(d["full_mean_ps"].size) if lam_of(w, d["J"]) > 0 else 0
    raise ValueError(name)


def ci(vec):
    v = vec[~np.isnan(vec)] if np.ndim(vec) else vec
    if np.size(v) == 0:
        return np.nan, np.nan, np.nan
    return float(np.nanstd(v, ddof=1)), float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def p_collapse_point(d):
    return float(d["collapsed"].mean())


def measure_all(arms):
    """Run each (label, arm_kind, alpha) at every J; return {J: {label: data}}."""
    data = {J: {} for J in CELLS_J}
    for J in CELLS_J:
        seed = seed_for_J(J)
        for label, kind, alpha in arms:
            data[J][label] = run_arm_ps(kind, alpha, J, seed)
    return data


def build_ci_rows(data, arms):
    rows = []
    for J in CELLS_J:
        idx = np.random.default_rng(BOOT_SEED + int(round(J * 10))).integers(0, N_SEEDS, size=(N_BOOT, N_SEEDS))
        for label, kind, alpha in arms:
            d = data[J][label]
            w = float(d["win_mean_ps"].mean())
            L = lam_of(w, J)
            row = dict(J=J, arm=label, alpha=alpha, p_collapse=p_collapse_point(d),
                       lam_at_hwin=L, in_domain=bool(L > 0))
            pcb = d["collapsed"][idx].mean(axis=1)
            se, lo, hi = ci(pcb)
            row.update(p_collapse_se=se, p_collapse_ci_lo=lo, p_collapse_ci_hi=hi)
            for s in SUMMARIES:
                row[s] = summary_on(s, d)
                boot = summary_on(s, d, idx)
                se, lo, hi = ci(boot)
                row[f"{s}_se"], row[f"{s}_ci_lo"], row[f"{s}_ci_hi"] = se, lo, hi
                row[f"{s}_n_defined"] = n_defined(s, d)
            rows.append(row)
    return rows


def ordering_verdicts(data, arms, label_set):
    """For each summary and J, sort the given arms ascending by the summary and
    test each adjacent step statistically."""
    rows = []
    for J in CELLS_J:
        idx = np.random.default_rng(BOOT_SEED + int(round(J * 10))).integers(0, N_SEEDS, size=(N_BOOT, N_SEEDS))
        labels = [l for (l, _, _) in arms if l in label_set]
        for s in SUMMARIES:
            pts = {l: summary_on(s, data[J][l]) for l in labels}
            # only arms on which the summary is defined can be ordered
            usable = [l for l in labels if not np.isnan(pts[l])]
            order = sorted(usable, key=lambda l: pts[l])
            for lo_l, hi_l in zip(order[:-1], order[1:]):
                dlo, dhi = data[J][lo_l], data[J][hi_l]
                d_amp = pts[hi_l] - pts[lo_l]
                amp_boot = summary_on(s, dhi, idx) - summary_on(s, dlo, idx)
                _, a_lo, a_hi = ci(amp_boot)
                amp_res = (a_lo > 0) == (a_hi > 0) and not (a_lo <= 0 <= a_hi)
                p_hi, p_lo = p_collapse_point(dhi), p_collapse_point(dlo)
                d_p = p_hi - p_lo
                p_boot = dhi["collapsed"][idx].mean(axis=1) - dlo["collapsed"][idx].mean(axis=1)
                _, p_clo, p_chi = ci(p_boot)
                p_res = not (p_clo <= 0 <= p_chi)
                is_viol = d_p > 0
                rows.append(dict(J=J, summary=s, arm_lo=lo_l, arm_hi=hi_l,
                                 d_amp=d_amp, d_amp_ci_lo=a_lo, d_amp_ci_hi=a_hi, amp_resolved=amp_res,
                                 p_lo=p_lo, p_hi=p_hi, d_p=d_p, d_p_ci_lo=p_clo, d_p_ci_hi=p_chi,
                                 p_resolved=p_res, is_violation=is_viol,
                                 violation_resolved=bool(is_viol and amp_res and p_res)))
    return rows


def point_ordering(data, labels):
    """Point-estimate ordering PASS/FAIL per summary over the given labels
    (P non-increasing as the summary increases). Returns {summary: (ok, detail)}."""
    out = {}
    for s in SUMMARIES:
        ok_all = True; detail = ""
        for J in CELLS_J:
            pts = {l: summary_on(s, data[J][l]) for l in labels}
            usable = [l for l in labels if not np.isnan(pts[l])]
            order = sorted(usable, key=lambda l: pts[l])
            ps = [p_collapse_point(data[J][l]) for l in order]
            for i in range(len(ps) - 1):
                if ps[i + 1] > ps[i] + 1e-9:
                    ok_all = False
                    detail = (f"J={J}: {order[i]}({pts[order[i]]:.3f}->P={ps[i]:.3f}) "
                              f"< {order[i+1]}({pts[order[i+1]]:.3f}->P={ps[i+1]:.3f})")
                    break
            if not ok_all:
                break
        n_undef = sum(1 for l in labels if np.isnan(summary_on(s, data[CELLS_J[0]][l])))
        out[s] = (ok_all, detail, n_undef)
    return out


def coupling_domain_scan():
    """A3 step 1: scan coupling alpha; record measured win_mean + lambda per J."""
    rows = []
    grid = [round(0.1 * k, 2) for k in range(1, 11)]  # 0.1 .. 1.0
    for a in grid:
        for J in CELLS_J:
            d = run_arm_ps("coupling", a, J, seed_for_J(J))
            w = float(d["win_mean_ps"].mean())
            L = lam_of(w, J)
            rows.append(dict(alpha=a, J=J, win_mean=w, lam=L, in_domain=bool(L > 0)))
    df = pd.DataFrame(rows)
    # largest alpha with lam>0 at ALL three J
    ok = df.groupby("alpha")["in_domain"].all()
    in_domain_alphas = [a for a in grid if ok.get(a, False)]
    astar = max(in_domain_alphas) if in_domain_alphas else None
    return df, astar


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base_arms = [(lbl, lbl, a) for lbl, a in ARMS]  # (label, kind, alpha)
    # kind must be one of passive/stress/coupling/quadratic -> map from label
    kind_of = {"passive": "passive", "stress_a2.0": "stress", "stress_a0.05": "stress",
               "coupling": "coupling", "quadratic": "quadratic"}
    base_arms = [(lbl, kind_of[lbl], a) for lbl, a in ARMS]

    print("=== measuring five base arms with per-seed arrays ===")
    data = measure_all(base_arms)

    # A3: coupling domain scan + in-domain coupling arm
    print("\n=== A3: coupling domain scan (largest alpha with lambda>0 at all J) ===")
    scan_df, astar = coupling_domain_scan()
    scan_df.to_csv(OUT / "coupling_domain_scan.csv", index=False)
    print(scan_df.to_string(index=False))
    print(f"\nchosen in-domain coupling alpha* = {astar}")
    sixth_label = f"coupling_a{astar}"
    for J in CELLS_J:
        data[J][sixth_label] = run_arm_ps("coupling", astar, J, seed_for_J(J))
    sixth_arm = (sixth_label, "coupling", astar)
    all_arms = base_arms + [sixth_arm]

    # A1: CI rows for all six arms
    ci_rows = build_ci_rows(data, all_arms)
    ci_df = pd.DataFrame(ci_rows)
    ci_df.to_csv(OUT / "delivered_mean_ordering_ci.csv", index=False)
    # deposit the sixth arm on its own too
    ci_df[ci_df.arm == sixth_label].to_csv(OUT / "coupling_in_domain_arm.csv", index=False)
    print(f"\nwrote delivered_mean_ordering_ci.csv ({len(ci_df)} rows) and coupling_in_domain_arm.csv")

    # in_domain flags per arm (from J=5 lambda, consistent across J for these arms)
    print("\n=== in_domain per arm (lambda>0 at measured selection-window field) ===")
    for J in CELLS_J:
        flags = {l: bool(lam_of(float(data[J][l]["win_mean_ps"].mean()), J) > 0)
                 for (l, _, _) in all_arms}
        print(f"J={J}: " + ", ".join(f"{k}={v}" for k, v in flags.items()))

    in_domain_labels = [l for (l, _, _) in all_arms
                        if lam_of(float(data[5.0][l]["win_mean_ps"].mean()), 5.0) > 0]
    all_labels = [l for (l, _, _) in all_arms]

    # A2: ordering test PRIMARY (in-domain) and SECONDARY (all)
    print("\n=== A2 PRIMARY ordering test (in-domain arms only:", in_domain_labels, ") ===")
    prim = point_ordering(data, in_domain_labels)
    for s in SUMMARIES:
        ok, det, nun = prim[s]
        print(f"  {s:24s}: {'PASS' if ok else 'FAIL'}  (undefined on {nun} arm(s))" + ("" if ok else f"  {det}"))

    print("\n=== A2 SECONDARY ordering test (all arms:", all_labels, ") ===")
    seco = point_ordering(data, all_labels)
    for s in SUMMARIES:
        ok, det, nun = seco[s]
        print(f"  {s:24s}: {'PASS' if ok else 'FAIL'}  (undefined on {nun} arm(s))" + ("" if ok else f"  {det}"))

    # statistical verdicts: PRIMARY (in-domain adjacency) and SECONDARY (all arms)
    def report(vd_df, tag):
        print(f"\n=== [{tag}] resolved violations (amplitude AND P gaps resolved at 95%) ===")
        rv = vd_df[vd_df.violation_resolved]
        print(rv[["J", "summary", "arm_lo", "arm_hi", "d_amp", "d_amp_ci_lo", "d_amp_ci_hi",
                  "d_p", "d_p_ci_lo", "d_p_ci_hi"]].to_string(index=False) if len(rv) else "  NONE")
        print(f"\n=== [{tag}] unresolved violations (point-estimate violation, not resolved at 95%) ===")
        uv = vd_df[(vd_df.is_violation) & (~vd_df.violation_resolved)]
        print(uv[["J", "summary", "arm_lo", "arm_hi", "d_amp", "d_amp_ci_lo", "d_amp_ci_hi",
                  "amp_resolved", "d_p", "d_p_ci_lo", "d_p_ci_hi", "p_resolved"]].to_string(index=False)
              if len(uv) else "  NONE")
        # per-summary: does a resolved violation exist at any J?
        print(f"\n=== [{tag}] per-summary: refuted (has a resolved violation)? ===")
        for s in SUMMARIES:
            sub = vd_df[vd_df.summary == s]
            res = sub[sub.violation_resolved]
            Js = sorted(set(res.J))
            print(f"  {s:24s}: {'REFUTED (resolved) at J=' + str(Js) if len(res) else 'no resolved violation'}")

    vd_all = pd.DataFrame(ordering_verdicts(data, all_arms, set(all_labels)))
    vd_all.to_csv(OUT / "ordering_verdicts.csv", index=False)
    print(f"\nwrote ordering_verdicts.csv ({len(vd_all)} rows, SECONDARY all-arm adjacency)")
    vd_prim = pd.DataFrame(ordering_verdicts(data, all_arms, set(in_domain_labels)))
    vd_prim.to_csv(OUT / "ordering_verdicts_primary.csv", index=False)
    print(f"wrote ordering_verdicts_primary.csv ({len(vd_prim)} rows, in-domain adjacency)")
    report(vd_prim, "PRIMARY in-domain")
    report(vd_all, "SECONDARY all-arm")


if __name__ == "__main__":
    main()
