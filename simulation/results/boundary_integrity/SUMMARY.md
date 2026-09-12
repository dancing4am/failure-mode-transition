# Boundary/Clipping Integrity

Base integrator: EM, dt=0.05, 20000 steps, T=2.0, xi=0.5, clip to [-1+1e-6, 1-1e-6]; collapse = W<10 for 200 steps; headline cell J=5, mu=100 (passive), seed 9057473 (0xC0DE + 9*1000003 + 8*1009).
Fixed-field cells use seeds 0xC0DE + j_idx*1000003 + 8*1009 with j_idx = 9,10,11,12 for J = 5,10,15,20 (natural grid extension); the same seed is shared by both h values at a given J. With h fixed the m-SDE decouples from wealth, so mu is irrelevant to clip statistics.

## (i) Clip-hit statistics (n=1000 per cell)

A 'clip hit' = the proposed EM update m + f dt + g dW falls outside [-1+1e-6, 1-1e-6] and is projected back. Fractions of steps are over all n * 20000 trajectory-steps.

| cell | J | h | steps clipped upper | steps clipped lower | traj ever upper | traj ever lower |
|---|---|---|---|---|---|---|
| headline_passive | 5 | passive W-coupled | 0.3526 | 0.0221 | 0.966 | 0.296 |
| fixed_h2.4_J5 | 5 | 2.4 (fixed) | 0.3298 | 0.0000 | 1.000 | 0.006 |
| fixed_h5.0_J5 | 5 | 5 (fixed) | 0.3959 | 0.0000 | 1.000 | 0.000 |
| fixed_h2.4_J10 | 10 | 2.4 (fixed) | 0.3849 | 0.0136 | 0.991 | 0.117 |
| fixed_h5.0_J10 | 10 | 5 (fixed) | 0.4014 | 0.0000 | 1.000 | 0.005 |
| fixed_h2.4_J15 | 15 | 2.4 (fixed) | 0.3598 | 0.0422 | 0.946 | 0.193 |
| fixed_h5.0_J15 | 15 | 5 (fixed) | 0.3989 | 0.0028 | 1.000 | 0.028 |
| fixed_h2.4_J20 | 20 | 2.4 (fixed) | 0.3379 | 0.0642 | 0.907 | 0.226 |
| fixed_h5.0_J20 | 20 | 5 (fixed) | 0.3881 | 0.0140 | 0.990 | 0.076 |

Headline-cell P(collapse) at n=1000: 0.291.

## (ii) Clipped-EM vs boundary-preserving Lamperti (headline cell, n=2000, paired noise)

Transform check (Ito, dm = f dt + xi sqrt(1-m^2) dW):

- u = artanh(m):  u' = 1/(1-m^2), u'' = 2m/(1-m^2)^2  =>
  du = [(f + xi^2 m)/(1-m^2)] dt + [xi/sqrt(1-m^2)] dW.
  The Ito correction is +xi^2 m/(1-m^2), but the diffusion is
  xi*cosh(u), NOT constant xi (a constant-diffusion reading of this
  route rests on a spurious sqrt(1-m^2) factor). Since m=+/-1 is attainable
  (part iii), u = artanh(m) reaches +/-inf in finite time and the scheme
  explodes: empirically, with n=200 paired trajectories, 100.0% exploded past |u|>50 (median explosion step 51). The artanh route is unusable here.
- The correct constant-diffusion (Lamperti) transform is u = arcsin(m)/xi:
  du = [(f + xi^2 m/2)/(xi sqrt(1-m^2))] dt + dW,  m = sin(xi u).
  m stays in [-1,1] identically (no clipping); u is folded back into
  [-pi/(2 xi), pi/(2 xi)] by exact reflection (sin is invariant), which
  realises the instantaneously-reflecting behaviour of the regular
  boundary. The integrable drift singularity at the boundary is tamed by
  a trust-region cap |drift dt| <= 0.5 (engaged on 5.27e-02 of trajectory-steps; folds on 2.00e-01).
  Identical per-step dW draws as the EM run (same generator, same call
  sequence) give a paired comparison.

| metric | clipped EM | Lamperti | paired diff [95% CI] |
|---|---|---|---|
| P(collapse) | 0.3070 | 0.2950 | +0.0120 [+0.0065, +0.0175] |
| rigidity share of collapsed | 0.8811 | 0.7441 |  |
| P(rigidity collapse), uncond. | 0.2705 | 0.2195 | +0.0510 [+0.0390, +0.0630] |
| collapse step median | 292 | 295 | -5.7 [-7.6, -3.8] (both-collapsed) |
| collapse step p99 | 425 | 447 |  |

Discordant pairs: EM-only collapses 28, Lamperti-only 4 (of n=2000).

## (iii) Feller boundary classification at m = +/-1

sigma(m) = xi sqrt(1-m^2), so sigma^2 = xi^2 (1-m^2) vanishes LINEARLY
at the boundary: with y = 1 -/+ m (distance to the boundary),
sigma^2 ~ 2 xi^2 y. The inward drift magnitude at the boundary is
  c_up = 1 - tanh((J+h)/T)  at m=+1,   c_dn = 1 + tanh((-J+h)/T)  at m=-1,
both > 0 (drift points inward). Near the boundary the process is the
CIR-type diffusion dy = c dt + xi sqrt(2y) dW. Its scale density is
s(y) ~ exp(-int 2c/(2 xi^2 y) dy) = y^(-c/xi^2) and speed density
m(y) ~ 1/(sigma^2 s) ~ y^(c/xi^2 - 1). Feller test:

- c >= xi^2  (Feller condition 2c >= 2 xi^2, i.e. Bessel-type dimension
  delta = 2c/xi^2 >= 2): int s diverges at 0 -> boundary UNATTAINABLE;
  int m converges -> ENTRANCE boundary.
- 0 < c < xi^2 (delta < 2): int s and int m both converge at 0 ->
  REGULAR boundary: attainable in finite time, NOT absorbing (positive
  inward drift, finite speed measure); the clip/reflection choice fixes
  the boundary behaviour as instantaneous reflection.

With T=2, xi=0.5 (xi^2=0.25), J=5:

| h | boundary | inward drift c | delta = 2c/xi^2 | classification |
|---|---|---|---|---|
| 0.4 | +1 | 8.993e-03 | 7.194e-02 | regular (attainable, reflecting; not absorbing) |
| 0.4 | -1 | 1.990e-02 | 1.592e-01 | regular (attainable, reflecting; not absorbing) |
| 2.4 | +1 | 1.222e-03 | 9.774e-03 | regular (attainable, reflecting; not absorbing) |
| 2.4 | -1 | 1.383e-01 | 1.106e+00 | regular (attainable, reflecting; not absorbing) |
| 5 | +1 | 9.080e-05 | 7.264e-04 | regular (attainable, reflecting; not absorbing) |
| 5 | -1 | 1.000e+00 | 8.000e+00 | entrance (unattainable) |

Conclusion: for the headline/passive range (h between 0 and ~2.4 during
the pre-collapse transient) BOTH boundaries are regular-attainable
(delta << 2), which is exactly why the EM integrator registers clip hits;
m=-1 becomes an entrance (unattainable) boundary only for h >~ 3.05
(where 1 - tanh((J-h)/T) >= xi^2). No boundary is absorbing in any
studied regime, so clipping approximates the correct reflecting
behaviour; part (ii) quantifies the residual discretisation effect.

Runtime: 11s.  Files: clip_stats.csv, lamperti_paired.csv, lamperti_summary.csv.
