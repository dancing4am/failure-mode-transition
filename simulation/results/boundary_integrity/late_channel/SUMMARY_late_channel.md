# Late collapse channel under a boundary-preserving integrator

**Question.** e3_collapse_times.csv (clipped EM, n=1000) shows that at
constant field h in {2.4, 5.0}, J in {15, 20}, mu=100, 8-28% of wealth
collapses occur after t=100 (p99 up to ~985). Is this late channel a
clipped-EM boundary artifact?

**Method.** Paired-noise comparison, n=2000 per cell: clipped EM (exact
replica of branch_observable.run_cell with fixed h + wealth/collapse
ledger) vs the constant-diffusion Lamperti scheme of
boundary_integrity.lamperti_run (u = arcsin(m)/xi, exact reflecting
fold at |u| = pi/(2 xi), trust cap |drift dt| <= 0.5), with
constant h and the identical wealth/collapse ledger. One dW draw per
(seed, step) feeds BOTH integrators. dt=0.05, 20000 steps, T=2,
xi=0.5, m0=0, W0=100; collapse = W<10 for 200 consecutive steps;
late = collapse time t > 100. Seeds 0xB00001 + i*1000003, cell order
(2.4,15), (2.4,20), (5.0,15), (5.0,20). Script:
simulation/scripts/late_channel_check.py; data:
simulation/results/boundary_integrity/late_channel/.

## Per-cell results (n=2000, paired)

| h | J | arm | P(collapse) | n coll | n late (t>100) | late frac | median t | p90 | p99 | max |
|---|---|---|---|---|---|---|---|---|---|---|
| 2.4 | 15 | EM | 0.1830 | 366 | 75 | 0.205 | 14.5 | 561.7 | 967.0 | 991.1 |
| 2.4 | 15 | Lamperti | 0.2930 | 586 | 302 | 0.515 | 128.7 | 826.0 | 981.2 | 999.5 |
| 2.4 | 20 | EM | 0.2470 | 494 | 127 | 0.257 | 14.5 | 632.6 | 955.9 | 990.4 |
| 2.4 | 20 | Lamperti | 0.3930 | 786 | 409 | 0.520 | 138.6 | 798.8 | 968.0 | 998.6 |
| 5.0 | 15 | EM | 0.0285 | 57 | 11 | 0.193 | 14.1 | 489.9 | 947.9 | 963.7 |
| 5.0 | 15 | Lamperti | 0.0390 | 78 | 35 | 0.449 | 18.7 | 852.4 | 991.5 | 997.8 |
| 5.0 | 20 | EM | 0.0715 | 143 | 45 | 0.315 | 14.5 | 705.9 | 948.3 | 984.5 |
| 5.0 | 20 | Lamperti | 0.1115 | 223 | 117 | 0.525 | 160.5 | 820.7 | 964.9 | 994.6 |

## Paired differences (EM minus Lamperti, 95% CI)

| h | J | dP(collapse) | dP(late collapse) | verdict on late channel |
|---|---|---|---|---|
| 2.4 | 15 | -0.1100 [-0.1259, -0.0941] | -0.1135 [-0.1287, -0.0983] | **survives** |
| 2.4 | 20 | -0.1460 [-0.1630, -0.1290] | -0.1410 [-0.1574, -0.1246] | **survives** |
| 5.0 | 15 | -0.0105 [-0.0164, -0.0046] | -0.0120 [-0.0173, -0.0067] | **survives** |
| 5.0 | 20 | -0.0400 [-0.0506, -0.0294] | -0.0360 [-0.0460, -0.0260] | **survives** |

## Sanity vs e3_collapse_times.csv (EM arm)

Different n (2000 vs 1000) and different seed base (0xB00001 vs
0xC0DE-derived), so only binomial-level agreement is expected.
Two-proportion z-scores (this run vs e3):

| h | J | e3 P | this EM P | z(P) | e3 late frac | this EM late frac | z(late) |
|---|---|---|---|---|---|---|---|
| 2.4 | 15 | 0.178 | 0.183 | +0.34 | 0.281 | 0.205 | -1.98 |
| 2.4 | 20 | 0.269 | 0.247 | -1.30 | 0.197 | 0.257 | +1.87 |
| 5.0 | 15 | 0.035 | 0.029 | -0.97 | 0.229 | 0.193 | -0.41 |
| 5.0 | 20 | 0.084 | 0.071 | -1.22 | 0.155 | 0.315 | +2.67 |

## Integrator diagnostics

| h | J | EM clip-step rate | Lamperti cap rate | Lamperti fold rate |
|---|---|---|---|---|
| 2.4 | 15 | 4.020e-01 | 5.325e-02 | 2.013e-01 |
| 2.4 | 20 | 4.021e-01 | 5.329e-02 | 2.013e-01 |
| 5.0 | 15 | 4.019e-01 | 5.327e-02 | 2.014e-01 |
| 5.0 | 20 | 4.021e-01 | 5.323e-02 | 2.013e-01 |

## Verdict

- (h=2.4, J=15): late channel **survives** under the
  boundary-preserving scheme (EM 75 vs Lamperti 302 late collapses of n=2000).
- (h=2.4, J=20): late channel **survives** under the
  boundary-preserving scheme (EM 127 vs Lamperti 409 late collapses of n=2000).
- (h=5.0, J=15): late channel **survives** under the
  boundary-preserving scheme (EM 11 vs Lamperti 35 late collapses of n=2000).
- (h=5.0, J=20): late channel **survives** under the
  boundary-preserving scheme (EM 45 vs Lamperti 117 late collapses of n=2000).

Runtime: 15s.

---

# Discretization convergence of the late channel (dt 0.05 -> 0.025)

**Question.** Do the late-channel results above survive halving the
step size?  Rerun of the paired EM + Lamperti comparison at dt=0.025,
n_steps=40000 (same physical horizon t_max=1000), n=1000
per cell, paired dW per (seed, step), same four cells, mu=100,
constant h, T=2, xi=0.5, m0=0, W0=100.  The collapse streak is scaled
with 1/dt — 400 consecutive steps at dt=0.025 (= 200 at
dt=0.05) — so collapse remains 'W<10 sustained for 10 physical time
units' at every dt.  Seeds 0xF00001 + i*1000003 (fresh; paths are NOT
shared across dt, so the dt comparison is CROSS-SAMPLE: two-sample
binomial z / unpooled 95% Wald CI on the differences).  Script:
simulation/scripts/late_channel_dt_convergence.py; data:
simulation/results/boundary_integrity/late_channel/dt_convergence.csv.

## Per-cell comparison (scheme x dt)

| h | J | scheme | dt | n | P(collapse) | late frac | median t | p99 t |
|---|---|---|---|---|---|---|---|---|
| 2.4 | 15 | clipped_EM | 0.05 | 2000 | 0.1830 | 0.205 | 14.4 | 967.0 |
| 2.4 | 15 | clipped_EM | 0.025 | 1000 | 0.1800 | 0.228 | 14.4 | 949.9 |
| 2.4 | 15 | lamperti | 0.05 | 2000 | 0.2930 | 0.515 | 128.7 | 981.2 |
| 2.4 | 15 | lamperti | 0.025 | 1000 | 0.2940 | 0.503 | 109.6 | 960.2 |
| 2.4 | 20 | clipped_EM | 0.05 | 2000 | 0.2470 | 0.257 | 14.5 | 955.9 |
| 2.4 | 20 | clipped_EM | 0.025 | 1000 | 0.2440 | 0.180 | 14.3 | 942.0 |
| 2.4 | 20 | lamperti | 0.05 | 2000 | 0.3930 | 0.520 | 138.6 | 968.0 |
| 2.4 | 20 | lamperti | 0.025 | 1000 | 0.4070 | 0.504 | 123.8 | 977.7 |
| 5.0 | 15 | clipped_EM | 0.05 | 2000 | 0.0285 | 0.193 | 14.1 | 947.9 |
| 5.0 | 15 | clipped_EM | 0.025 | 1000 | 0.0300 | 0.200 | 14.5 | 920.7 |
| 5.0 | 15 | lamperti | 0.05 | 2000 | 0.0390 | 0.449 | 18.7 | 991.5 |
| 5.0 | 15 | lamperti | 0.025 | 1000 | 0.0370 | 0.432 | 16.2 | 954.8 |
| 5.0 | 20 | clipped_EM | 0.05 | 2000 | 0.0715 | 0.315 | 14.5 | 948.3 |
| 5.0 | 20 | clipped_EM | 0.025 | 1000 | 0.0910 | 0.187 | 14.4 | 913.3 |
| 5.0 | 20 | lamperti | 0.05 | 2000 | 0.1115 | 0.525 | 160.4 | 964.9 |
| 5.0 | 20 | lamperti | 0.025 | 1000 | 0.1230 | 0.390 | 15.8 | 945.0 |

## dt differences (dt=0.025 minus dt=0.05, two-sample 95% CI)

| h | J | scheme | dP [95% CI] | z(P) | dLateFrac [95% CI] | z(LF) | within CI? |
|---|---|---|---|---|---|---|---|
| 2.4 | 15 | clipped_EM | -0.0030 [-0.0322, +0.0262] | -0.20 | +0.0229 [-0.0511, +0.0968] | +0.61 | P yes, LF yes |
| 2.4 | 15 | lamperti | +0.0010 [-0.0336, +0.0356] | +0.06 | -0.0120 [-0.0820, +0.0581] | -0.33 | P yes, LF yes |
| 2.4 | 20 | clipped_EM | -0.0030 [-0.0356, +0.0296] | -0.18 | -0.0768 [-0.1385, -0.0150] | -2.32 | P yes, LF NO |
| 2.4 | 20 | lamperti | +0.0140 [-0.0232, +0.0512] | +0.74 | -0.0167 [-0.0765, +0.0432] | -0.55 | P yes, LF yes |
| 5.0 | 15 | clipped_EM | +0.0015 [-0.0113, +0.0143] | +0.23 | +0.0070 [-0.1690, +0.1830] | +0.08 | P yes, LF yes |
| 5.0 | 15 | lamperti | -0.0020 [-0.0165, +0.0125] | -0.27 | -0.0163 [-0.2104, +0.1778] | -0.16 | P yes, LF yes |
| 5.0 | 20 | clipped_EM | +0.0195 [-0.0016, +0.0406] | +1.88 | -0.1279 [-0.2384, -0.0174] | -2.16 | P yes, LF NO |
| 5.0 | 20 | lamperti | +0.0115 [-0.0131, +0.0361] | +0.93 | -0.1344 [-0.2427, -0.0261] | -2.40 | P yes, LF NO |

## Verdict

**NOT CONVERGED.** Cells outside the two-sample 95% CI:

- (h=2.4, J=20, clipped_EM): late-frac -0.0768 (dt=0.025 minus dt=0.05)
- (h=5.0, J=20, clipped_EM): late-frac -0.1279 (dt=0.025 minus dt=0.05)
- (h=5.0, J=20, lamperti): late-frac -0.1344 (dt=0.025 minus dt=0.05)

Runtime: 18s.
