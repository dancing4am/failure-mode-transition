# α-FIT — collapse-observable scaling exponent: pre-registration

> **Editorial amendment (2026-09-09).** The original frozen document is retained in the
> private working repository at git blob `1edcae3986c285c29e9d90a2d5b716ec1bd85903`
> (SHA-256 of the original file content:
> `db126578610f3c5d539dfa710658d10787ed5b7377627dbbe9e1a288bd572f85`). Thirteen lines at
> original line numbers 49, 219, 238–239, 254–257, and 261–265 were reworded or removed
> for public release because they contained internal workflow/approval phrasing, an
> internal work-stream label, or the name of a non-shipped file; on 2026-09-11 the two
> pointer-map lines at release line numbers 248–249 were further neutralized (an internal
> cluster label and an internal summary-document name). No substantive content
> of §1–§6 — the pipeline lock, fit target and μ-pairing rule, grid, fit/statistics,
> decision rule, or held-cohort resolution rule — was changed. The original is available
> from the author on request.

**Date**: 2026-05-30
**Status**: FROZEN pre-registration. MEASUREMENT NOT performed; NO
code written; held cohort NOT touched. Same discipline as the
Substrate-A probe and the DOPO/Substrate-D re-measurement: freeze the
design (operating point, fit target, grid, statistics, decision rule)
in a timestamped commit BEFORE any measurement, one pre-committed
decision rule, honest-failure path, no post-hoc parameter tuning.

This document timestamps the discriminating re-measurement that will
settle the *observable* collapse-scaling exponent and thereby resolve
the empirical "observed-slope" cohort that was deliberately HELD at
`h/J²` in commit `8a8a345` while the analytical η set was corrected to
`h/J`.

---

## 0. Why this pre-registration exists (settled context, not re-litigated here)

A single label `h/J²` had been carrying **two distinct quantities**:

1. **Analytical η** = ΔV_tilt / ΔV_barrier of the effective potential.
   The barrier does **not** grow as `J`; it approaches a finite ceiling
   (`ΔV_barrier → ½`), so `η = ΔV_tilt/ΔV_barrier ≈ 4h/J → h/J`. This is
   **SETTLED** (Theorem 1; confirmed by the unit-error audit) and was
   **APPLIED** in `8a8a345` (η, Theorem 1, SI §S2.4–§S2.6, finance
   §S6.4, lineage re-anchor, DOPO interim text).

2. **The observable** `P(collapse | J)` and the J-falloff of the
   **passive-stabilizer effectiveness** plotted in Figure 8 Panel B —
   the manuscript's "theoretical `h/J²` envelope". This is a separate,
   model-specific quantity (branch selection under wealth-ruin with an
   endogenous field), NOT η, and SI §S2.6 already states we do **not**
   assert it shares η's leading exponent. These loci were **HELD** at
   `h/J²` pending the present measurement.

**Two facts about the observable are already pre-committed and are NOT
re-opened by this measurement:**

- **`α = 2` is DECISIVELY EXCLUDED a priori.** Two independent routes
  agree: (i) the branch-selection error-function shoulder gives an
  asymptotic exponent of −½, never −2; (ii) Monte-Carlo of the
  constant-h "ended at m<0" toy gives `ΔP·(J−T)²` growing monotonically
  `0.35 → 3.9` across `J = 3…9` (a true `(J−T)^−2` law would hold this
  product flat). Therefore the held cohort's slope-2 / `h/J²` **envelope
  is wrong as a power law.** `α = 2` is **not** a live hypothesis below;
  if a fit CI nonetheless contains 2 that is an UNEXPECTED result to be
  reported as such, not adopted.

- **The manuscript observable ≠ the constant-h toy.** The full pipeline
  (wealth-ruin + endogenous `h(W)` + m↔W feedback) gives
  `P(J=5, μ=100) = 0.30`, matching the manuscript's headline `≈0.33`;
  the constant-h "m<0" toy gives `0.15`. Neither the toy's asymptotic
  −½ nor its pre-asymptotic ≈−1 (cubic-tanh term `c₂ = −(J²/T²)
  sech²(h/T)tanh(h/T)`, an O(1) collapse-promoting **term** that
  steepens the toy) is assumed to be the manuscript exponent. **The
  manuscript observable's exponent must be measured directly**, with the
  manuscript's exact collapse definition and pipeline. That is what this
  pre-registration freezes.

The live hypotheses for the measured exponent are therefore
**α = ½** vs **α ≈ 1** vs **pre-asymptotic drift between them** vs
**no clean power law**. `α = 2` is excluded a priori.

---

## 1. Pipeline (LOCKED — exact manuscript pipeline, NOT the toy)

The measurement runs the **passive-stabilizer condition of the
manuscript's own simulator** (`simulation/scripts/active_stabilizer.py`,
`h_passive`), which is byte-for-byte the minimal-model collapse pipeline
and the Methods protocol of `manuscript/paper.md` (L276–L290). No toy,
no constant-h surrogate.

Per-step update (T = 2, ξ = 0.5, dt = 0.05, ε = 1e−6,
`m ∈ [−1+ε, 1−ε]`):

```
m_c   = clip(m, −1+ε, 1−ε)
h     = 2·(W/500)                      # passive stabilizer (baseline field)
f     = −m_c + tanh((J·m_c + h)/T)     # drift
g     = ξ·sqrt(1 − m_c²)               # multiplicative noise
dW    = sqrt(dt)·𝒩(0,1)
m     = clip(m_c + f·dt + g·dW, −1+ε, 1−ε)

employment  = 0.5 + 0.5·m
income      = μ·employment
consumption = 30 + 10·(W/500)
W           = max(0, W + (income − consumption)·dt)
```

**Fixed initial conditions / constants** (from `active_stabilizer.py`
L28–L46 and paper.md L276–L290):

- `W(0) = 100`, `m(0) = 0`.
- **Horizon = 20,000 steps** at `dt = 0.05` (= 1,000 physical time
  units). This is the manuscript horizon (paper.md L278) and
  `N_STEPS = 20_000`; it is **pinned**, not a tunable.
- **Collapse definition (verbatim manuscript):** a run is *collapsed*
  if `W < 10` (10% of initial) for **≥ 200 consecutive steps**
  (10 physical time units). `COLLAPSE_WEALTH = 10`,
  `COLLAPSE_STREAK = 200`. (Collapse-type typing — rigidity `|m|>0.9`,
  fragmentation `|m|<0.3` — is recorded but is NOT the fit target here.)
- Matched-seed scheme identical to `active_stabilizer.py`
  (`seed_base + j_idx·1_000_003 + m_idx·1009`) so cells are reproducible
  and paired across μ.

**Validation gate (pre-committed, runs BEFORE any α is fit):** the
implementation must reproduce the manuscript headline cell
**`P(collapse | J=5, μ=100) = 0.30 ± 0.03`** (independently confirmed;
manuscript reports ≈0.33). If it does not, **HALT and surface** — do not
tune constants to hit the target.

---

## 2. Fit target (LOCKED)

The held "envelope" lives on the manuscript's **passive-stabilizer
effectiveness** curve (Figure 8 Panel B, paper.md L1255–L1257):

> *passive stabilizer effectiveness `(1 − P(J,μ=100)/P(J,μ=20))` versus
> J … plus the theoretical `h/J²` envelope (purple dotted, normalized).*

The measurement fits the exponent of **this exact quantity**, reporting
**both** a normalized and a raw form:

- **Primary — efficacy ratio** `E(J) = 1 − P(J, μ_high) / P(J, μ_low)`,
  both P under the passive stabilizer. This **is** the Figure 8 Panel B
  curve; the held "`h/J²` envelope" is the `α = 2` line on it.
- **Secondary — raw gap** `ΔP(J) = P(J, μ_low) − P(J, μ_high)` (the
  un-normalized protective reduction). Reported alongside so the result
  does not hinge on the ratio's normalization.

**μ-pairing rule (resolves the saturation defect in the published
figure).** Figure 8 Panel B uses `μ_low = 20`, but at `μ = 20` the
low-margin `P` **saturates at 1.00 across the J grid**, which degenerates
`E(J)` (the ratio becomes `1 − P_high` and its J-dependence collapses
onto `P_high` alone). Therefore:

- **Primary `μ_low = 40`** (and `μ_high = 100`): chosen because both
  `P_low` and `P_high` stay in `(0,1)` across the J grid, so `E(J)`
  carries genuine information at both ends.
- Also report the **manuscript-matching `μ_low = 20`** pair, **explicitly
  flagged as saturating**, for continuity with the published Figure 8.
- Also report `μ_low = 60`.
- **Selection criterion:** the headline α is fit on the μ pairing that
  keeps **both** `P_low, P_high ∈ (0,1)` across the entire J grid
  (expected to be `μ_low = 40`, `μ_high = 100`). If more than one pairing
  qualifies, the `40/100` pairing is the pre-committed headline.

Margins swept: **μ ∈ {40, 60, 100}** (μ = 20 reported only as the
saturating control). Comparing pairings also tests separability /
h-linearity of the effect.

---

## 3. J grid (LOCKED)

Coupling threshold `J_c = T = 2`; the scaling form is fit for `J > J_c`.

- **Full window:** `J ∈ {2.25, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10, 12}`
  (dense near `J_c` to resolve curvature; extends past the manuscript's
  published `J ≤ 5` primary grid and the Figure-8 `J ≤ 10` extension to
  reach the asymptotic window — deliberate, using the identical
  per-run collapse pipeline).
- **High-J subwindow:** `J ≥ 6`, grid `{6, 7, 8, 9, 10, 12}` (6 points).
  This is the asymptotic-window fit that distinguishes a clean limiting
  exponent from pre-asymptotic curvature.

---

## 4. Fit form and statistics (LOCKED)

- **Form:** `E(J) ∝ (J − T)^(−α)` with `T = 2` (likewise `ΔP(J)`). Fit
  by OLS of `log E` on `log(J − T)`:
  - over the **full window** (`J > 2`), and
  - over the **J ≥ 6 subwindow**,
  reporting `α̂`, intercept, and R² for each.
- **Model-free cross-check:** per-J **local slope**
  `d log E / d log(J − T)` between adjacent grid points, reported as a
  sequence so any monotone drift (e.g. ½ → 1 or 1 → ½ with rising J) is
  visible without assuming a single global exponent.
- **Statistics:** seed-resampling **bootstrap 95% CI** on `α̂`. Increase
  seeds per cell until **bootstrap SE(α) ≤ 0.05** (so the ±0.1 target
  resolves ½ and 1, separated by 0.5, cleanly). Floor **≥ 200
  seeds/cell** (manuscript used 100); escalate as needed to meet the SE
  target. Per-cell `P` reported with its own bootstrap CI.

---

## 5. Pre-committed decision rule

Report `α̂ ± 95% CI` for `E(J)` and `ΔP(J)`, over the full window and the
`J ≥ 6` subwindow, plus the local-slope sequence. Classify by the
**headline (40/100) E(J) fit**, with the subwindow and ΔP as
corroboration:

1. **CI contains ½, excludes 1** → **consistent-with-½** (branch-selection
   asymptote dominates). Resolve the held cohort to **`h/√J`**.
2. **CI contains 1, excludes ½** → **consistent-with-1** (cubic-term
   steepening dominates in the accessible window). Resolve the held
   cohort to **`h/J`** (numerically coincident with the analytical η
   exponent, though for a distinct reason — to be stated, not conflated).
3. **`α = 2` is EXCLUDED a priori** (§0). If a CI nonetheless contains 2,
   that is an **UNEXPECTED** result → surface for investigation; do
   **not** auto-resolve the cohort to it.
4. **Full-window and J ≥ 6 disagree** (e.g. full ≈ 1 but subwindow → ½,
   or local slopes drift monotonically) → **pre-asymptotic drift**.
   Resolve to the **asymptotic (J ≥ 6) value** and state the
   pre-asymptotic caveat explicitly in the manuscript.
5. **No clean power law** (poor log-log fit / low R² / non-monotone local
   slopes with no clean limit) → **qualitative reframe**: drop the
   power-law "envelope" language entirely; replace with a described,
   measured monotone decline of effectiveness with J. Held loci are
   **reworded**, not re-exponented.

The held cohort is resolved to the **MEASURED** α (or reworded under
outcome 5) — **not** to either prior derivation guess.

---

## 6. Held-cohort resolution map (NOT executed in this commit)

For the record only. When (and only when) α settles, these loci move
from `h/J²` to the measured value (or are reworded under outcome 5).
**No edit to any of these is made now.**

- Abstract effect-size phrase — `paper.md L15`, `main.tex L148`,
  `paper_only.tex L148`.
- Figure 3 title "Network robustness of the h/J² scaling" and the
  Figure 8 Panel B "theoretical h/J² envelope" — `paper.md L1181`,
  `L1256`; `main.tex L437`, `L482`.
- Slope-attribution prose — SI `L1129–1130`; `main.tex L1288`.
- Robustness / envelope references — `paper.md L345`, `L766`;
  `main.tex L228`, `L318`.
- The observed-residual clause — `main.tex L775`.
- Headline effect-size line in the journal reporting materials.

(SI §S2.4–§S2.6 analytical η loci are already at `h/J` and are NOT in
this map.)

---

## 7. Honest-failure commitment / discipline

- This pre-registration is **frozen on commit**. The operating point
  (§1), fit target and μ rule (§2), grid (§3), fit/statistics (§4), and
  decision rule (§5) are fixed in advance.
- **No post-hoc tuning** of constants, horizon, collapse threshold,
  grid, or μ to manufacture a preferred exponent. The validation gate
  (§1) is a halt-and-surface check, not a tuning target.
- The **held cohort is untouched** until α settles (§6).
- The **measurement run begins only after this pre-registration is frozen.**

---

*Frozen pre-registration commit. NO measurement performed, NO code
written, NO held-cohort edits at the time of freezing. The
measurement run, its implementation, and any held-cohort resolution
all follow the pipeline lock (§1), the fit target + μ-pairing rule
(§2), the grid (§3), and the decision rule (§5).*
