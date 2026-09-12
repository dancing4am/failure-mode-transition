# Pre-specification — Mechanism of the late collapse channel

> **Editorial amendment (2026-09-11).** The original frozen document is retained in the
> private working repository at git blob `19fb197bd106dcd4af5a5c0655d5c49ea6aa2e32`
> (SHA-256 of the original file content:
> `a3c917bb752f87d6e4691329a6adc6167beaf1b5de89a39f6355a544d87c4f47`). One line (original
> line 5, the Status line) was reworded for public release to remove a version-control
> process remark. No element of the hypotheses, tests, decision rules, or analysis plan
> was changed. The original is available from the author on request.

**Date written:** 2026-08-22
**Author:** Chowon Jung
**Status:** written BEFORE any result in `simulation/results/late_channel_mechanism/` is generated (pre-specification). Committed in its own commit so the ordering is independently verifiable in the repository history.

---

## 1. Background and hypothesis

The manuscript reports two collapse channels in the fixed-field cells (constant field *h*, wealth ledger driven by *m*, collapse = *W* < 10 for 200 consecutive steps):

- an **early** channel — Suzuki transient branch selection near the unstable state, active only when the trajectory starts near the disordered state and must select a branch;
- a **late** channel — collapses occurring long after the transient (reported up to *t* ≈ 985), currently declared open with no mechanism assigned. Prior work excluded the two obvious artifacts: the late channel *strengthens* under the boundary-preserving (Lamperti) integrator rather than being a clipping artifact, and it survives a d*t* convergence check.

**Hypothesis (H1).** The late channel is ordinary noise-activated (Kramers) escape from the metastable high-*m* well over the deterministic drift barrier. This is a distinct process from transient branch selection, which is why a single selection formula cannot describe both.

**Null / alternative (H0).** The late collapses are not activated escape from the well; they too depend on the initial branch-selection transient (e.g. slow stragglers of the same selection process), so they should vanish when branch selection is removed.

All tests use the boundary-preserving Lamperti integrator (`u = arcsin(m)/ξ`, exact reflecting fold at |*u*| = π/(2ξ), trust-region cap |drift·d*t*| ≤ 0.5), constant field *h*, the identical wealth/collapse ledger as `late_channel_check.py`, and base settings d*t* = 0.05, 20,000 steps (*t*_max = 1000), *T* = 2, ξ = 0.5, *μ* = 100, *W*₀ = 100.

Cells: *h* ∈ {2.4, 5.0} × *J* ∈ {7, 10, 15, 20}, *n* = 2,000 trajectories per cell. Per-cell seed = 0xB0_0001 + *i*·1,000,003 (row-major cell index *i*); the two initial conditions within a cell share the same seed so the comparison is paired on the noise stream.

Time windows (physical time): **early** = collapse time *t* ≤ 50; **late** = collapse time *t* > 250. (The transient is complete well before *t* = 250; the gap between 50 and 250 is deliberately excluded from both so the two windows are clean.)

---

## 2. L1 — Initial-condition dependence (decisive test)

For every cell, run two initialisations with the same seeds:

- **(a)** *m*₀ = 0 (disordered start): both channels can operate.
- **(b)** *m*₀ = +0.9 (deep in the protected well): branch selection cannot operate.

**Reported per cell per IC:** P(collapse); collapse-time median, p99, and count with *t* > 250; count of early (*t* ≤ 50) and late (*t* > 250) collapses; and the **late-collapse hazard rate** *r* = (number of collapses in the window (250, 1000]) ÷ (person-time at risk in that window), where person-time at risk = Σ over trajectories uncollapsed at *t* = 250 of (min(collapse time, 1000) − 250). *r* is the maximum-likelihood constant hazard and is the quantity that should be initial-condition-independent for activated escape.

**Decision rule (fixed in advance):**
- **Supported** if, pooled across cells, IC (b) removes at least **80 %** of the early collapses (Σ early(b) ≤ 0.20 · Σ early(a)) **and** the late hazard is preserved within a factor of two: the ratio *r*(b)/*r*(a) ∈ **[0.5, 2.0]**, evaluated per cell where both arms have ≥ 10 late collapses and on the pooled rate.
- **Refuted** if the late collapses disappear together with the early ones: pooled *r*(b) < **0.25** · *r*(a).
- **Inconclusive** otherwise (early not sufficiently removed, too few late events to estimate *r*, or *r*(b)/*r*(a) between 0.25 and 0.5).

---

## 3. L2 — Waiting-time distribution (memorylessness)

Activated escape is a Poisson process: first-passage times into the low basin are exponential, so the collapse times in the late window are exponential up to the fixed wealth-drain lag between *m*-escape and *W* < 10.

For each cell (IC (a)), pool the late collapses (*t* > 250) and form the waiting times *τ* = collapse time − 250. To remove the deterministic drain lag, also form the **m-escape** waiting times: the time at which *m* first crosses below the deterministic saddle *m*_s, restricted to escapes after *t* = 250, minus 250. The escape-time series is the cleaner Poisson observable; the collapse-time series is reported for completeness.

**Reported per cell (where ≥ 20 late events):** coefficient of variation CV = SD/mean of the waiting times; a Q–Q plot against a fitted exponential; the one-sample Kolmogorov–Smirnov statistic and its *p* value against an exponential fitted to the same data (shifted by the observed minimum for the collapse-time series).

**Decision rule (fixed in advance):** consistent with memoryless activated escape if **CV ∈ [0.75, 1.25]** and the KS test does **not** reject the exponential at **p > 0.05**. Cells with < 20 late events are reported but not scored.

---

## 4. L3 — Barrier scaling

For each cell compute the deterministic drift *f*(*m*) = −*m* + tanh((*Jm* + *h*)/*T*) and its fixed points (stable high-*m* well *m*₊, unstable saddle *m*_s). Two barrier measures:

1. **Drift potential (as specified):** *U*(*m*) = ∫[*m* − tanh((*Jm* + *h*)/*T*)] d*m*; barrier Δ*U* = *U*(*m*_s) − *U*(*m*₊).
2. **Multiplicative-noise effective barrier (physically correct for this SDE):** the Itô stationary exponent gives escape rate *r* ∝ exp(−Δψ) with ψ(*m*) = ∫ 2*f*/(ξ²(1 − *m*²)) d*m*; write Δψ = (2/ξ²)·*G* where *G* = [Ĝ(*m*₊) − Ĝ(*m*_s)], Ĝ(*m*) = ∫ *f*/(1 − *m*²) d*m*. This isolates the ξ dependence in the prefactor 2/ξ².

Regress log(late hazard *r*) — measured in L1, IC (a) — on the barrier across the 8 cells.

**Decision rule (fixed in advance):**
- Activated-escape signature confirmed if log *r* is linear in Δ*U* with **R² ≥ 0.80** and negative slope.
- **ξ-consistency:** regress log *r* on *G*; Kramers predicts slope −2/ξ² = **−8** at ξ = 0.5. Treat as consistent if the fitted slope is within **±40 %** of −8, i.e. in [−11.2, −4.8]. (Cells with < 10 late events are excluded from the regression; the multiplicative-noise prefactor and finite-barrier corrections are acknowledged as caveats.)
- **Noise sweep** (if cheap): at one cell, ξ ∈ {0.35, 0.5, 0.65}; the late hazard must increase monotonically with ξ (lower effective barrier), and log *r* must move in the direction predicted by exp(−2Δ*V*/ξ²)-type scaling.

---

## 5. Overall verdict

L1 is decisive; L2 and L3 are corroborating mechanism signatures.

- **Supported (late channel is activated escape):** L1 supported, and at least one of {L2 consistent, L3 R² ≥ 0.80 with correct sign and a physically plausible slope}.
- **Refuted:** L1 refuted (late hazard collapses when branch selection is removed).
- **Inconclusive:** any other combination; the report will state exactly which measurement would settle it.

No manuscript edits follow from this pre-specification. Results go to `simulation/results/late_channel_mechanism/`.
