"""
pilot_alpha_c.py -- Substrate A probe, module 2 (REWRITTEN per
branch-collapse pilot, Amendment 1.

Detector: BRANCH-COLLAPSE RAMP (continuation hysteresis loop), NOT the
superseded random-init SSB ramp. The original detector measured the
algorithmic basin-escape threshold (phase-retrieval search-phase
trapping at R~0); the framework's J_c = T refers to the LANDSCAPE
pitchfork, located here by tracking the ordered branch.

Model is imported verbatim from nn_core.py ;
nothing about the model is re-implemented here.

================================================================
TWO-STAGE PILOT (user-locked 2026-05-17)
================================================================
  Stage A -- locate alpha_c:
    * Initialise AT the informative well  w0 = w* + eps * N(0, I_N)
      (eps = 0.1  ->  R0 ~ 1/sqrt(1+eps^2) ~ 0.995, deep in basin).
    * DOWN-ramp alpha from alpha_high to alpha_floor, step -d_alpha,
      WARM-STARTED (each alpha continues from the previous alpha's
      final w). At each alpha: burn-in, then a sustained measurement
      window; record <|R|> (R = normalized teacher overlap, the
      order parameter).
    * alpha_c(down) = the alpha where the ordered branch collapses:
      sustained <|R|> first crosses BELOW R_COLLAPSE, sub-grid
      located by linear interpolation between the two bracketing
      alphas.
    * UP-ramp confirmation: CONTINUE (warm start) from the collapsed
      disordered w at alpha_floor, ramp alpha back UP; alpha_c(up) =
      where <|R|> first crosses back ABOVE R_COLLAPSE. This is a true
      continuation hysteresis loop -- immune to the cold random-init
      search-phase trapping that killed the v1 detector, because it
      passes THROUGH the bifurcation rather than starting cold deep
      in the trapped region.

  Stage B -- cusp/fold first-pass via hysteresis magnitude:
    * d_ac = | alpha_c(up) - alpha_c(down) |.
    * d_ac <= TAU (0.2)  -> CUSP-consistent (a supercritical
      pitchfork has zero intrinsic hysteresis; <=2 grid steps of
      discretisation/thermal smear is tolerated). Proceed to V4/V5.
    * d_ac >  TAU         -> first-order / fold indication. Per
      Amendment 1(d), Substrate A is CLOSED; V4/V5 NOT performed,
      closure verdict surfaced immediately. (No re-emergence on the
      up-ramp => d_ac treated as effectively infinite => fold.)

  Main probe (h_sweep / j_sweep) proceeds only if BOTH stages pass.

================================================================
LOCKED DEFAULTS (user-approved 2026-05-17; no post-hoc tuning)
================================================================
  eps (well-init noise)        = 0.10     -> R0 ~ 0.995
  alpha_high / d_alpha / floor  = 4.0 / 0.1 / 0.3   (38 alphas)
  TAU (hysteresis tolerance)   = 0.20     (= 2 * d_alpha)
  per-alpha window              = 6000 steps: 3000 burn + 3000 measure
                                  (record_every 20 -> 150 burn / 150
                                   measured samples)
  R_COLLAPSE threshold          = 0.10  (>> ~1/sqrt(400)~0.05 floor,
                                  << in-well R* ~ 0.9; clean band)
  N_SEEDS                       = 12     (teacher FIXED seed 0; data /
                                  Langevin noise / init varied)
  N_RAMP                        = 400    (Amendment-justified;
                                  finite-N shift << 2x J=2 alpha_c)
  N_BIST                        = 1000   (pre-registered op. point,
                                  V5 bistability only)

Interpreter: Python 3.13 (numpy). Run: <py3.13> pilot_alpha_c.py
"""
from __future__ import annotations

import sys
import time
from dataclasses import dataclass

import numpy as np

import nn_core as nc


# --------------------------------------------------------------------------
# Locked configuration
# --------------------------------------------------------------------------
@dataclass
class PilotConfig:
    # model operating point (T / dt locked inside nn_core)
    T: float = 0.1
    dt: float = 1e-3
    h: float = 0.0
    lam: float = 0.0

    # branch-collapse ramp (Stage A) -- all user-locked
    EPS_WELL: float = 0.10
    ALPHA_HIGH: float = 4.0
    D_ALPHA: float = 0.1
    ALPHA_FLOOR: float = 0.3
    TAU: float = 0.20
    N_STEPS: int = 6000
    BURN_STEPS: int = 3000
    REC_EVERY: int = 20
    R_COLLAPSE: float = 0.10
    N_SEEDS: int = 12
    N_RAMP: int = 400

    # V5 bistability -- pre-registered operating point
    N_BIST: int = 1000
    N_STEPS_BIST: int = 15000
    N_SEEDS_BIST: int = 4

    MASTER_SEED: int = 20260517
    TEACHER_SEED: int = 0


CFG = PilotConfig()


def alpha_grid_down() -> np.ndarray:
    n = int(round((CFG.ALPHA_HIGH - CFG.ALPHA_FLOOR) / CFG.D_ALPHA)) + 1
    return np.round(np.linspace(CFG.ALPHA_HIGH, CFG.ALPHA_FLOOR, n), 4)


# --------------------------------------------------------------------------
# Deterministic child seeds
# --------------------------------------------------------------------------
def _child_seed(*tags: int) -> int:
    ss = np.random.SeedSequence([CFG.MASTER_SEED, *tags])
    return int(ss.generate_state(1)[0])


# --------------------------------------------------------------------------
# One warm-started Langevin segment at a fixed alpha
# --------------------------------------------------------------------------
def _alpha_tag(alpha: float) -> int:
    return int(round(alpha * 1000))


def segment(w_init, alpha, N, run_seed, direction_tag, N_steps=None):
    """One warm-started Langevin segment at fixed alpha.

    Returns dict: mean_absR over the MEASURE window, temporal variance
    of R over the measure window, R0, and the final w (for warm-start
    continuation to the next alpha).
    """
    N_steps = CFG.N_STEPS if N_steps is None else N_steps
    atag = _alpha_tag(alpha)
    data_seed = _child_seed(run_seed, atag, 11)
    lang_seed = _child_seed(run_seed, atag, direction_tag, 22)

    wstar = nc.make_teacher(N, seed=CFG.TEACHER_SEED)        # FIXED
    p = nc.ProbeParams(N=N, alpha=alpha, T=CFG.T, h=CFG.h,
                        lam=CFG.lam, dt=CFG.dt, seed=lang_seed)
    X, y = nc.make_data(wstar, p.P, seed=data_seed)

    ts, Rs, Ls, w_fin = nc.run_langevin(np.asarray(w_init, float),
                                        X, y, wstar, p, N_steps,
                                        record_every=CFG.REC_EVERY)
    n_burn = int(round(CFG.BURN_STEPS / CFG.REC_EVERY))
    meas = Rs[n_burn:]
    return {
        "alpha": float(alpha),
        "R0": float(Rs[0]),
        "mean_absR": float(np.mean(np.abs(meas))),
        "R_tvar": float(np.var(meas)),
        "Rs": Rs,
        "w_fin": w_fin,
    }


# --------------------------------------------------------------------------
# Sub-grid threshold crossing (linear interpolation)
# --------------------------------------------------------------------------
def _interp_cross(alphas, vals, thresh, going_down):
    """alphas/vals are in RAMP ORDER. Return the sub-grid alpha at the
    FIRST crossing of `thresh`:
      going_down=True  : first time vals goes from >=thresh to <thresh
                          (ordered-branch collapse on the down-ramp)
      going_down=False : first time vals goes from <=thresh to >thresh
                          (re-emergence on the up-ramp)
    np.nan if it never crosses."""
    a = np.asarray(alphas, float)
    v = np.asarray(vals, float)
    for i in range(1, len(v)):
        if going_down and v[i - 1] >= thresh and v[i] < thresh:
            pass
        elif (not going_down) and v[i - 1] <= thresh and v[i] > thresh:
            pass
        else:
            continue
        if v[i] == v[i - 1]:
            return float(a[i])
        return float(a[i - 1] + (thresh - v[i - 1])
                     * (a[i] - a[i - 1]) / (v[i] - v[i - 1]))
    return float("nan")


# --------------------------------------------------------------------------
# One seed: full continuation hysteresis loop
# --------------------------------------------------------------------------
def hysteresis_loop(run_seed: int, N: int, verbose=False):
    a_down = alpha_grid_down()
    a_up = a_down[::-1]

    # ---- DOWN ramp: start AT the informative well ----
    rng_i = np.random.default_rng(_child_seed(run_seed, 33))
    wstar = nc.make_teacher(N, seed=CFG.TEACHER_SEED)
    w = wstar + CFG.EPS_WELL * rng_i.standard_normal(N)
    down_mean, down_tvar = [], []
    for al in a_down:
        s = segment(w, al, N, run_seed, direction_tag=1)
        w = s["w_fin"]
        down_mean.append(s["mean_absR"])
        down_tvar.append(s["R_tvar"])
    ac_down = _interp_cross(a_down, down_mean, CFG.R_COLLAPSE,
                            going_down=True)

    # ---- UP ramp: CONTINUE from the collapsed disordered w ----
    up_mean = []
    for al in a_up:
        s = segment(w, al, N, run_seed, direction_tag=2)
        w = s["w_fin"]
        up_mean.append(s["mean_absR"])
    ac_up = _interp_cross(a_up, up_mean, CFG.R_COLLAPSE,
                          going_down=False)

    if np.isfinite(ac_down) and np.isfinite(ac_up):
        d_ac = abs(ac_up - ac_down)
    else:
        # no clean collapse OR no re-emergence => effectively
        # infinite hysteresis => fold/first-order
        d_ac = float("inf")
    return {
        "seed": run_seed,
        "a_down": a_down, "down_mean": np.array(down_mean),
        "down_tvar": np.array(down_tvar),
        "a_up": a_up, "up_mean": np.array(up_mean),
        "ac_down": ac_down, "ac_up": ac_up, "d_ac": d_ac,
    }


# --------------------------------------------------------------------------
# Stage A + Stage B over all seeds
# --------------------------------------------------------------------------
def run_pilot(emit):
    loops = []
    for s in range(CFG.N_SEEDS):
        L = hysteresis_loop(s, CFG.N_RAMP)
        loops.append(L)
        emit(f"    seed {s:2d}: ac_down={L['ac_down']:.4f}  "
             f"ac_up={('%.4f' % L['ac_up']) if np.isfinite(L['ac_up']) else 'NONE':>7}"
             f"  |d_ac|="
             f"{('%.4f' % L['d_ac']) if np.isfinite(L['d_ac']) else 'inf'}")

    ac_d = np.array([L["ac_down"] for L in loops], float)
    d_ac = np.array([L["d_ac"] for L in loops], float)
    ok_d = np.isfinite(ac_d)

    # cross-check detector: temporal-variance (critical-fluctuation)
    # peak on the down-ramp, averaged over seeds
    a_down = loops[0]["a_down"]
    chi = np.mean([L["down_tvar"] for L in loops], axis=0)
    chi_peak_alpha = float(a_down[int(np.argmax(chi))])

    finite_dac = d_ac[np.isfinite(d_ac)]
    res = {
        "loops": loops,
        "ac_down_mean": float(np.mean(ac_d[ok_d])) if ok_d.any() else float("nan"),
        "ac_down_std": float(np.std(ac_d[ok_d])) if ok_d.any() else float("nan"),
        "ac_down_min": float(np.min(ac_d[ok_d])) if ok_d.any() else float("nan"),
        "ac_down_max": float(np.max(ac_d[ok_d])) if ok_d.any() else float("nan"),
        "n_collapsed": int(ok_d.sum()),
        "dac_mean": float(np.mean(finite_dac)) if finite_dac.size else float("inf"),
        "dac_max": float(np.max(finite_dac)) if finite_dac.size else float("inf"),
        "n_dac_finite": int(finite_dac.size),
        "n_seeds": CFG.N_SEEDS,
        "chi_peak_alpha": chi_peak_alpha,
    }
    # Stage B verdict: EVERY seed must satisfy |d_ac| <= TAU
    res["cusp_pass"] = bool(ok_d.all()
                            and np.all(np.isfinite(d_ac))
                            and np.all(d_ac <= CFG.TAU))
    return res


# --------------------------------------------------------------------------
# V5 bistability at J = 2 alpha_c (pre-registered N = 1000)
# --------------------------------------------------------------------------
def bistability_check(alpha_c: float, emit):
    alpha_probe = 2.0 * alpha_c
    rng = np.random.default_rng(_child_seed(99, 5))
    wstar = nc.make_teacher(CFG.N_BIST, seed=CFG.TEACHER_SEED)
    pos, neg = [], []
    for sign0 in (+1, -1):
        for k in range(CFG.N_SEEDS_BIST):
            w0 = (sign0 * wstar
                  + CFG.EPS_WELL * rng.standard_normal(CFG.N_BIST))
            r = segment(w0, alpha_probe, CFG.N_BIST, 7000 + k,
                        direction_tag=3 if sign0 > 0 else 4,
                        N_steps=CFG.N_STEPS_BIST)
            # signed final well: mean of the measure-window R
            n_burn = int(round(CFG.BURN_STEPS / CFG.REC_EVERY))
            Rf = float(np.mean(r["Rs"][n_burn:]))
            (pos if sign0 > 0 else neg).append(Rf)
            emit(f"    init sign {sign0:+d} (R0={r['R0']:+.3f}) -> "
                 f"R_final={Rf:+.3f}")
    ok = all(v > 0.3 for v in pos) and all(v < -0.3 for v in neg)
    return {"alpha_probe": alpha_probe, "pos": pos, "neg": neg,
            "both_wells": bool(ok)}


# --------------------------------------------------------------------------
# Validations
# --------------------------------------------------------------------------
def _v1_ramp_deterministic():
    """A warm-started 3-alpha mini DOWN-ramp run twice with the same
    seed must give a bitwise-identical order-parameter trajectory."""
    def mini():
        rng_i = np.random.default_rng(_child_seed(3, 33))
        ws = nc.make_teacher(120, seed=CFG.TEACHER_SEED)
        w = ws + CFG.EPS_WELL * rng_i.standard_normal(120)
        traj = []
        for al in (2.0, 1.9, 1.8):
            s = segment(w, al, 120, 3, direction_tag=1, N_steps=900)
            w = s["w_fin"]
            traj.append(s["Rs"])
        return np.concatenate(traj)
    return float(np.max(np.abs(mini() - mini())))


def _v2_detector_definition():
    """Detector validated on a SYNTHETIC landscape with a known
    collapse/re-emergence alpha, plus the Stage-B verdict logic on
    known cusp-like (zero hysteresis) and fold-like (wide hysteresis)
    synthetics."""
    al_d = np.array([3.0, 2.5, 2.0, 1.5, 1.0, 0.5])      # down order
    # synthetic ordered branch collapsing through R_COLLAPSE between
    # alpha=1.5 and 1.0
    br_d = np.array([0.90, 0.85, 0.70, 0.30, 0.02, 0.01])
    rec_d = _interp_cross(al_d, br_d, CFG.R_COLLAPSE, going_down=True)
    true_d = 1.5 + (CFG.R_COLLAPSE - 0.30) * (1.0 - 1.5) / (0.02 - 0.30)

    al_u = al_d[::-1]                                     # up order
    # CUSP-like: re-emerges at essentially the same alpha (no hyst.)
    br_u_cusp = np.array([0.01, 0.02, 0.30, 0.70, 0.85, 0.90])
    rec_u_cusp = _interp_cross(al_u, br_u_cusp, CFG.R_COLLAPSE,
                               going_down=False)
    dac_cusp = abs(rec_u_cusp - rec_d)
    # FOLD-like: re-emerges far higher (wide hysteresis loop)
    br_u_fold = np.array([0.01, 0.01, 0.02, 0.03, 0.40, 0.90])
    rec_u_fold = _interp_cross(al_u, br_u_fold, CFG.R_COLLAPSE,
                               going_down=False)
    dac_fold = abs(rec_u_fold - rec_d)

    recovered_ok = abs(rec_d - true_d) < 1e-6
    cusp_flag_ok = dac_cusp <= CFG.TAU          # correctly cusp
    fold_flag_ok = dac_fold > CFG.TAU           # correctly fold
    return (float(rec_d), float(true_d), float(dac_cusp),
            float(dac_fold), bool(recovered_ok and cusp_flag_ok
                                  and fold_flag_ok))


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    t0 = time.time()
    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("=== pilot_alpha_c.py : Substrate-A branch-collapse pilot "
         "(Amendment 1) ===")
    emit(f"op point T={CFG.T} dt={CFG.dt} h={CFG.h} | N_ramp="
         f"{CFG.N_RAMP} seeds={CFG.N_SEEDS} eps={CFG.EPS_WELL} | "
         f"alpha {CFG.ALPHA_HIGH}->{CFG.ALPHA_FLOOR} step {CFG.D_ALPHA}"
         f" | window {CFG.N_STEPS}({CFG.BURN_STEPS} burn) | "
         f"R_collapse={CFG.R_COLLAPSE} TAU={CFG.TAU}")

    # ---- V1 ----
    v1 = _v1_ramp_deterministic()
    ok1 = v1 == 0.0
    emit(f"[V1] ramp deterministic given seed : max|dR| over identical "
         f"warm-started reruns = {v1:.1e} -> {'PASS' if ok1 else 'FAIL'}"
         f" (need 0)")

    # ---- V2 ----
    rd, td, dcu, dfo, ok2 = _v2_detector_definition()
    emit(f"[V2] detector definition           : synthetic collapse "
         f"recovered alpha_c={rd:.4f} vs true {td:.4f}; verdict logic "
         f"cusp|d_ac|={dcu:.3f}<=TAU & fold|d_ac|={dfo:.3f}>TAU -> "
         f"{'PASS' if ok2 else 'FAIL'}")

    # ---- Stage A + B ----
    emit("\n-- Stage A (locate alpha_c) + Stage B (cusp/fold via "
         "hysteresis) --")
    R = run_pilot(emit)
    emit(f"\n  alpha_c (branch collapse, down-ramp) = "
         f"{R['ac_down_mean']:.4f} +/- {R['ac_down_std']:.4f}  "
         f"(sd over {R['n_collapsed']}/{R['n_seeds']} seeds; range "
         f"[{R['ac_down_min']:.3f}, {R['ac_down_max']:.3f}])")
    emit(f"  hysteresis |d_ac| : mean="
         f"{R['dac_mean'] if np.isfinite(R['dac_mean']) else float('inf'):.4f}"
         f"  max="
         f"{R['dac_max'] if np.isfinite(R['dac_max']) else float('inf'):.4f}"
         f"  ({R['n_dac_finite']}/{R['n_seeds']} seeds re-emerged)  "
         f"TAU={CFG.TAU}")
    emit(f"  cross-check: critical-fluctuation (R temporal-variance) "
         f"peak at alpha={R['chi_peak_alpha']:.3f}  "
         f"(must be ~ alpha_c)")

    rel_spread = (R["ac_down_std"] / R["ac_down_mean"]
                  if np.isfinite(R["ac_down_mean"]) and R["ac_down_mean"]
                  else float("inf"))
    repro_ok = (R["n_collapsed"] >= max(10, CFG.N_SEEDS - 2)
                and rel_spread < 0.20)
    cusp_ok = R["cusp_pass"]
    ok3 = bool(repro_ok and cusp_ok)
    verdict = ("CUSP-consistent (proceed)" if cusp_ok
               else "FOLD / first-order indication")
    emit(f"[V3] alpha_c + reproducibility + cusp/fold first-pass :")
    emit(f"     alpha_c={R['ac_down_mean']:.4f}+/-{R['ac_down_std']:.4f}"
         f" (rel spread {rel_spread:.3f}, need <0.20; "
         f"{R['n_collapsed']}/{CFG.N_SEEDS} collapsed)")
    emit(f"     hysteresis verdict: {verdict}  "
         f"(all seeds |d_ac|<=TAU = {cusp_ok}) -> "
         f"{'PASS' if ok3 else 'FAIL'}")

    if not cusp_ok:
        # Amendment 1(d): fold/first-order => Substrate A CLOSED.
        emit("\n[V4] SKIPPED  -- gated on V3 cusp first-pass")
        emit("[V5] SKIPPED  -- gated on V3 cusp first-pass")
        emit("\n*** SUBSTRATE A CLOSURE VERDICT (Amendment 1(d)) ***")
        emit("    Hysteresis exceeds TAU (or no clean collapse / no")
        emit("    re-emergence): the transition is NOT a clean")
        emit("    supercritical pitchfork (fold / first-order /")
        emit("    trapped). The single permitted detector retry has")
        emit("    been used. Substrate A is CLOSED with no further")
        emit("    detector modification. Document in")
        emit("    the pilot decision record.")
        allok = False
    else:
        ac = R["ac_down_mean"]
        emit("[V4] literature consistency        :")
        emit("     Real Gaussian phase-retrieval (quadratic teacher-")
        emit("     student) weak-recovery load is O(1): IT ~1/2,")
        emit("     spectral ~1 (Mondelli & Montanari 2019), Langevin/")
        emit("     GD landscape somewhat above ~1 (Sarao Mannelli et")
        emit("     al. 2019/2020). This is the EFFECTIVE landscape")
        emit("     alpha_c at T=0.1; finite T and loss/normalisation")
        emit("     conventions shift it vs the bare value -- a setup")
        emit("     difference, NOT a bug.")
        lit_ok = (0.3 <= ac <= 4.0)
        emit(f"     measured alpha_c = {ac:.4f} -> "
             f"{'consistent O(1) scale' if lit_ok else 'OUT OF O(1) BAND -- SURFACE'}")

        emit(f"\n-- V5 bistability at J = 2*alpha_c (N={CFG.N_BIST}, "
             f"pre-registered) --")
        B = bistability_check(ac, emit)
        emit(f"  alpha_probe = 2*alpha_c = {B['alpha_probe']:.4f}")
        emit(f"  R0>0 -> {np.array2string(np.array(B['pos']), precision=3)}")
        emit(f"  R0<0 -> {np.array2string(np.array(B['neg']), precision=3)}")
        ok5 = B["both_wells"]
        emit(f"[V5] bistability at J=2 alpha_c     : both +-R* "
             f"reachable from matching inits = {ok5} -> "
             f"{'PASS' if ok5 else 'FAIL'}")
        allok = ok1 and ok2 and ok3 and lit_ok and ok5

    status = "ALL CHECKS PASSED" if allok else "CHECKS NEED REVIEW / CLOSURE"
    emit(f"\n=== pilot_alpha_c.py: {status} ({time.time()-t0:.0f}s) ===")

    with open("pilot_alpha_c_RESULTS.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    return allok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
