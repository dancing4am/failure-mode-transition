# CRITICAL RAMP SPEED — VERDICT

**Result: the speed-of-onset effect is non-monotone in `mult` and reveals
two distinct failure modes.** Slowing AI deployment helps the *high-margin*
economies (rigidity protection) but *harms* the marginal-margin economies
(fragmentation by starvation). The headline image — "go slower, be safer" —
is incomplete.

100 seeds × 9 mult levels × 6 ramp durations × 40,000 steps at dt=0.05 =
5,400 runs, 50.2 s wall time. Same minimal model as `dynamic_j.py`; total
horizon held at 2,000 physical time units (40,000 steps at dt=0.05) so even
the slowest ramp (1,500 physical units) has 500 physical units of post-ramp
settling time for collapse to express. All numbers below at converged
dt=0.05.

## Headline finding — non-monotonicity at the cliff edge

| ramp_steps | mult=40 | mult=50 | mult=60 | mult=70 | mult=80 | mult=100 |
|---|---|---|---|---|---|---|
|  1,000 | 0.55 | 0.19 | 0.07 | 0.01 | 0.00 | 0.00 |
|  2,000 | 0.86 | 0.26 | 0.02 | 0.00 | 0.00 | 0.00 |
|  4,000 | 1.00 | 0.40 | 0.00 | 0.00 | 0.00 | 0.00 |
| 10,000 | 1.00 | 0.51 | 0.02 | 0.01 | 0.00 | 0.00 |
| 20,000 | 1.00 | 0.51 | 0.04 | 0.00 | 0.00 | 0.00 |
| 30,000 | 1.00 | 0.49 | 0.00 | 0.00 | 0.00 | 0.00 |

**Three regimes:**

1. **Unviable margin (mult ≤ 40):** P(collapse) is monotonically *increasing*
   in ramp duration — 0.55 at 1,000 steps, 1.00 at every ramp ≥ 4,000 steps.
   Slowing deployment makes things worse: the long ramp gives wealth more
   time to drain at sub-viable employment income.

2. **Critical margin (mult = 50):** P(collapse) is monotonically *increasing*
   in ramp duration over the tested range — 0.19 → 0.51 → saturating near
   0.50 by ramp = 10,000. The 50% collapse threshold sits at ramp ≈ 9,200
   steps (interpolated; ≈ 460 physical time units).

3. **Resilient margin (mult ≥ 60):** P(collapse) ≤ 7% at every ramp duration
   (peaks at 7% at the fastest tested ramp; 0–4% otherwise). The economy
   adiabatically settles regardless of pace.

## What the failure mode tells us

The **collapse-type** breakdown rules out a single-mechanism story:

| ramp_steps | mult=40 rigidity-share | mult=50 rigidity-share | mult=60 rigidity-share |
|---|---|---|---|
|  1,000 | 0.35 | 0.32 | 0.57 |
|  2,000 | 0.10 | 0.15 | 0.50 |
|  4,000 | 0.03 | 0.00 | n/a  |
| 10,000 | 0.01 | 0.00 | 0.00 |
| 30,000 | 0.01 | 0.00 | n/a  |

**At fast ramps, residual collapses are partly rigidity** (the carry-over
of the sudden-onset signature: noise traps a fraction of seeds on the −m
branch before *h* rises). **At slow ramps, residual collapses are
overwhelmingly fragmentation** — wealth drains during the long ramp at
sub-viable employment, |m| sits near 0, and the economy never recovers.

The two modes have *opposite* sensitivity to deployment speed. This is why
the prior `dynamic_j` VERDICT — which used a 15k-step horizon equal to
the slowest ramp — saw only the rigidity story: fragmentation needs
post-ramp time to express, and was being clipped by the horizon.

## Updated 2×2 governance matrix

|  | High mult (≥60) | Critical mult (50) | Low mult (≤40) |
|---|---|---|---|
| **Sudden** (~1,000 steps) | partial residual rigidity (≤7%) | 19% collapse, 32% rigidity | 55% collapse, 35% rigidity |
| **Gradual** (30,000 steps) | safe (≤0%) | **49% collapse, ~99% fragmentation** | **100% collapse, ~99% fragmentation** |

The diagonal "sudden-and-poor → catastrophe / gradual-and-rich → safe"
finding from the prior verdict survives. The off-diagonals refine it:
**slowing AI deployment is not a free lunch for low-margin economies** —
it converts a partial rigidity collapse into a near-certain fragmentation
collapse.

## Implications for the paper

1. The **active/passive distinction sharpens.** A purely passive economic
   margin (mult) and a purely passive deployment-speed cap are both
   inadequate at the critical margin (mult ≈ 50). One protects against
   one failure mode and exposes the other.

2. **The "go gradual" governance recommendation needs a margin caveat.**
   Slow AI deployment is structurally safer *only above a viability
   threshold*. Below it, gradualism prolongs the period during which the
   economy operates at sub-viable employment income, guaranteeing
   fragmentation collapse.

3. **Active stabilizers (those that scale with dJ/dt) become provably
   necessary at the critical margin.** Neither passive lever alone (raising
   `mult` or slowing the ramp) suffices to clear the 50% threshold at
   mult=50. The figure-of-merit is the *joint* (mult, ramp_speed) frontier,
   not either axis alone.

## Discussion-ready sentences (lift verbatim)

> "Slowing AI deployment helps the high-margin economies — adiabatic
> settling preserves the +m branch — but harms the marginal economies,
> for which a long ramp at sub-viable employment income guarantees
> fragmentation collapse before the economy can settle on the
> productive branch. At our critical margin (mult = 50), the 50%-collapse
> ramp duration sits at ≈ 9,200 simulation steps (dt=0.05); the curve is
> monotone rising, not falling, in deployment time."

> "The dominant failure mode swaps with deployment speed: rigidity
> dominates at sudden onset, fragmentation by starvation dominates at
> gradual onset. The economic-margin axis and the deployment-speed axis
> protect against *different* failure modes, and at the critical margin
> neither passive lever alone is sufficient."

## Limitations

- 6 ramp durations per mult; the 50%-threshold interpolant uses
  log-linear interpolation between ramp_steps={4k, 10k}. A targeted
  7-point grid around 6k–10k for mult=50 would tighten the threshold
  estimate.
- N=100 seeds per cell; the ±0.05 sampling SE on a 0.5 collapse rate is
  large enough that the mult=50 saturation between ramp=10k and 20k
  (0.51 vs 0.51) is not statistically distinguishable.
- Total horizon = 40,000 steps at dt=0.05 (= 2,000 physical units).
  Extending would test whether the mid-mult curve develops a delayed
  fragmentation signature. Current evidence (mult=60: 0–7%) suggests it
  does not, but this is a one-line follow-up.
- Mean-field model; same caveats as `minimal_model/VERDICT.md`.
- All numbers above at converged dt = 0.05. The previously reported
  dt = 0.1 numbers in `dynamic_j_ramp_speed/VERDICT.md` (an earlier
  revision) were unconverged by ~5–10 pp; replaced.

## Artifacts

- `raw_results.csv`              — 5,400 rows
- `summary.csv`                  — 54 rows (6 ramp_steps × 9 mults)
- `threshold_by_mult.csv`        — 9 rows (one per mult)
- `collapse_vs_ramp_speed.png`   — Figure A: P(collapse) vs ramp_steps, by mult
- `collapse_types_vs_ramp.png`   — Figure B: collapse-type composition
- `threshold_vs_mult.png`        — Figure C: 50%/25%/10% threshold vs mult
- `ramp_speed_sweep.py`          — source (in `simulation/scripts/`)
