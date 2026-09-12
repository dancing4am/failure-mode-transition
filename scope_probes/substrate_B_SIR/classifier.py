"""
classifier.py -- Substrate B foundational module.

Exhaustive Option-C-style failure-mode classifier. Classes:
  controlled | rigidity (endemic lock-in) | fragmentation |
  partial-drift.

FRAGMENTATION requires ALL FOUR locked criteria (designed to exclude
transient heterogeneity, which is demoted to partial-drift):
  C1  >=2-modal final distribution of unit order parameters:
       GMM-BIC favors >=2 components (BIC_2 < BIC_1) AND one component
       mean < 0 (controlled) and one > 0 (endemic)  AND
       Hartigan dip test rejects unimodality (p < 0.05).
  C2  persistence: over the final 50% window, inter-unit Var_i(m_i)(t)
       is not significantly decreasing (OLS slope not signif. < 0).
  C3  stability: the partition is still >=2-modal in an extended segment
       (>=2x horizon + perturbation kick) supplied by the runner.
  C4  not frozen isolation: holds at eps > 0 (caller supplies eps; the
       eps = 0 trivial-independence case is flagged, not counted).

RIGIDITY: strong majority of units terminal m > 0, low field-sensitivity
|dm/dh| (supplied across the h-grid), and no endemic-well escape.

Interpreter: Python 3.13 (numpy/sklearn present there, not in 3.14).
"""
from __future__ import annotations

import numpy as np
from sklearn.mixture import GaussianMixture
import diptest as _diptest  # validated Hartigan dip (Phase-4a fix:
# the hand-rolled GCM/LCM version was incorrect -- it flagged a clean
# normal as multimodal, p=0.003. A pre-registered statistical test must
# use a validated implementation.)


def dip_test(samples: np.ndarray, n_boot: int = 0, seed: int = 0):
    """Hartigan & Hartigan (1985) dip test via the validated `diptest`
    package. Returns (dip, p_value); p is the analytic calibration
    against the unimodal null. `n_boot`/`seed` kept for call-site
    compatibility (unused). n < 4 -> non-informative (dip 0, p 1).
    """
    x = np.asarray(samples, dtype=float).ravel()
    if x.size < 4:
        return 0.0, 1.0
    d, p = _diptest.diptest(x)
    return float(d), float(p)


# --------------------------------------------------------------------------
# Modality of the unit-mean distribution (C1).
# --------------------------------------------------------------------------
def bimodal_split(unit_means: np.ndarray, dip_p_thresh: float = 0.05,
                  seed: int = 0):
    """C1: GMM-BIC(>=2) with one negative- and one positive-mean
    component, AND Hartigan dip rejecting unimodality.
    Returns (is_bimodal: bool, info: dict).
    """
    u = np.asarray(unit_means, dtype=float).reshape(-1, 1)
    n = u.shape[0]
    info = {"n_units": n}
    if n < 4:
        info["reason"] = "too_few_units"
        return False, info
    g1 = GaussianMixture(1, covariance_type="full",
                         random_state=seed).fit(u)
    g2 = GaussianMixture(2, covariance_type="full", n_init=3,
                         random_state=seed).fit(u)
    bic1, bic2 = g1.bic(u), g2.bic(u)
    means = np.sort(g2.means_.ravel())
    gmm_ok = (bic2 < bic1) and (means[0] < 0.0) and (means[-1] > 0.0)
    d, p = dip_test(u.ravel(), seed=seed)
    info.update({"bic1": float(bic1), "bic2": float(bic2),
                 "gmm_means": means.tolist(), "dip": float(d),
                 "dip_p": float(p), "gmm_ok": bool(gmm_ok)})
    return bool(gmm_ok and (p < dip_p_thresh)), info


def _ols_slope_negative(t: np.ndarray, y: np.ndarray) -> bool:
    """True if y is SIGNIFICANTLY decreasing in t (one-sided ~95%)."""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    if t.size < 5 or np.allclose(t, t[0]):
        return False
    A = np.vstack([t, np.ones_like(t)]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    slope = coef[0]
    yhat = A @ coef
    dof = max(t.size - 2, 1)
    s2 = np.sum((y - yhat) ** 2) / dof
    sxx = np.sum((t - t.mean()) ** 2)
    se = np.sqrt(s2 / sxx) if sxx > 0 else np.inf
    if not np.isfinite(se) or se == 0:
        return False
    tstat = slope / se
    return tstat < -1.66  # one-sided, ~95% for moderate dof


def classify(m_units_t: np.ndarray, t: np.ndarray, *,
             eps: float | None = None,
             field_sensitivity: float | None = None,
             extended_unit_means: np.ndarray | None = None,
             final_frac: float = 0.5,
             rigidity_majority: float = 0.8,
             sensitivity_thresh: float = 0.5):
    """Classify one run.

    m_units_t : array [n_units, n_times] of the order parameter
                (n_units = 1 for well-mixed Config 1).
    eps       : homogenizing-coupling value (Config 2) for C4 flag.
    field_sensitivity : |dm*/dh| for this (h,J) cell (rigidity criterion);
                if None, the within-run no-escape proxy is used.
    extended_unit_means : final unit means from the >=2x-horizon +
                perturbation continuation (C3). If None, C3 is marked
                'not_supplied' and fragmentation cannot be asserted.

    Returns (label, diagnostics).
    """
    M = np.atleast_2d(np.asarray(m_units_t, dtype=float))
    n_units, n_t = M.shape
    k0 = int(final_frac * n_t)
    final = M[:, k0:]
    unit_final_mean = final.mean(axis=1)
    diag = {"n_units": n_units,
            "unit_final_mean_summary":
                [float(np.min(unit_final_mean)),
                 float(np.median(unit_final_mean)),
                 float(np.max(unit_final_mean))]}

    # ---- single-unit (Config 1): fragmentation impossible by design ----
    if n_units == 1:
        mbar = float(unit_final_mean[0])
        escaped = (np.sign(final[0]).std() > 0.5)  # crossed the saddle
        if mbar > 0 and not escaped:
            low_sens = (field_sensitivity is None
                        or abs(field_sensitivity) < sensitivity_thresh)
            diag["low_field_sensitivity"] = bool(low_sens)
            return ("rigidity" if low_sens else "partial-drift"), diag
        if mbar < 0 and not escaped:
            return "controlled", diag
        return "partial-drift", diag

    # ---- multi-unit (Config 2/3) ----
    # C1 modality
    c1, c1info = bimodal_split(unit_final_mean)
    diag["C1_bimodal"] = c1info
    # C2 persistence: inter-unit variance over the final window
    var_t = M[:, k0:].var(axis=0)
    tt = np.asarray(t, float)[k0:]
    var_decreasing = _ols_slope_negative(tt, var_t)
    c2 = not var_decreasing
    diag["C2_persistent"] = bool(c2)
    diag["C2_var_decreasing"] = bool(var_decreasing)
    # C3 stability on extended segment
    if extended_unit_means is None:
        c3 = False
        diag["C3_stability"] = "not_supplied"
    else:
        c3, c3info = bimodal_split(np.asarray(extended_unit_means, float))
        diag["C3_stability"] = {"still_bimodal": bool(c3), **c3info}
    # C4 not frozen isolation
    if eps is None:
        c4 = True
        diag["C4_eps"] = "not_supplied(assumed_ok)"
    else:
        c4 = eps > 0.0
        diag["C4_eps"] = float(eps)

    if c1 and c2 and c3 and c4:
        return "fragmentation", diag

    # Not fragmentation: resolve the remaining classes.
    frac_endemic = float(np.mean(unit_final_mean > 0))
    if frac_endemic >= rigidity_majority:
        low_sens = (field_sensitivity is None
                    or abs(field_sensitivity) < sensitivity_thresh)
        diag["low_field_sensitivity"] = bool(low_sens)
        return ("rigidity" if low_sens else "partial-drift"), diag
    if frac_endemic <= (1.0 - rigidity_majority):
        return "controlled", diag
    # heterogeneous but failed C1-C3 (e.g. transient that homogenized,
    # or unstable partition) -> transient heterogeneity bucket
    return "partial-drift", diag


# --------------------------------------------------------------------------
# Self-tests
# --------------------------------------------------------------------------
def _tests():
    rng = np.random.default_rng(0)

    # dip test: unimodal normal -> do NOT reject; clear bimodal -> reject
    uni = rng.normal(0, 1, 200)
    d_u, p_u = dip_test(uni, n_boot=300, seed=1)
    bim = np.concatenate([rng.normal(-3, 0.3, 100),
                          rng.normal(+3, 0.3, 100)])
    d_b, p_b = dip_test(bim, n_boot=300, seed=1)
    assert p_u > 0.05, f"unimodal wrongly rejected (p={p_u:.3f})"
    assert p_b < 0.05, f"bimodal not detected (p={p_b:.3f}, dip={d_b:.4f})"
    assert d_b > d_u, "bimodal dip should exceed unimodal dip"

    # synthetic FRAGMENTATION: K=12 units, half controlled half endemic,
    # persistent variance, stable extended segment, eps>0 -> 'fragmentation'
    nT = 400
    half = np.concatenate([np.full(6, -0.85), np.full(6, +0.85)])
    M = (half[:, None]
         + rng.normal(0, 0.03, (12, nT)))     # flat (persistent) variance
    t = np.linspace(0, 100, nT)
    ext = half + rng.normal(0, 0.03, 12)      # still bimodal after 2x+kick
    lab, dg = classify(M, t, eps=0.05, extended_unit_means=ext)
    assert lab == "fragmentation", f"expected fragmentation, got {lab}"

    # synthetic RIGIDITY: all units locked endemic, low field sensitivity
    Mr = np.full((12, nT), 0.9) + rng.normal(0, 0.01, (12, nT))
    lab_r, _ = classify(Mr, t, eps=0.5, field_sensitivity=0.05,
                        extended_unit_means=np.full(12, 0.9))
    assert lab_r == "rigidity", f"expected rigidity, got {lab_r}"

    # synthetic CONTROLLED: all units disease-free
    Mc = np.full((12, nT), -0.9) + rng.normal(0, 0.01, (12, nT))
    lab_c, _ = classify(Mc, t, eps=0.5,
                        extended_unit_means=np.full(12, -0.9))
    assert lab_c == "controlled", f"expected controlled, got {lab_c}"

    # TRANSIENT heterogeneity must NOT be fragmentation: starts split,
    # variance collapses to one mode over the final window -> partial-drift
    Mt = np.empty((12, nT))
    for i in range(12):
        start = -0.85 if i < 6 else 0.85
        Mt[i] = np.linspace(start, 0.85, nT) + rng.normal(0, 0.02, nT)
    lab_t, dgt = classify(Mt, t, eps=0.05,
                          extended_unit_means=np.full(12, 0.85))
    assert lab_t != "fragmentation", \
        f"transient heterogeneity misclassified as fragmentation ({dgt})"

    # Config-1 single unit: rigidity vs controlled
    s_rig = np.full((1, nT), 0.8)
    assert classify(s_rig, t, field_sensitivity=0.05)[0] == "rigidity"
    s_ctl = np.full((1, nT), -0.8)
    assert classify(s_ctl, t)[0] == "controlled"

    print("classifier: ALL TESTS PASSED")
    print(f"  dip unimodal p={p_u:.3f} (>0.05), bimodal p={p_b:.3f} "
          f"(<0.05); transient->{lab_t}")


if __name__ == "__main__":
    _tests()
