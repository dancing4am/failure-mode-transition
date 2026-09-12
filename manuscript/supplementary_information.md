# Supplementary Information

**Paper:** *Transient Branch Selection and Delivered Stabilizer Amplitude in a Bistable Mean-Field Model: Matched-Seed Measurements of Collapse Probability*

**Author:** Chowon Jung, Department of Computer Science, University of Colorado Boulder, Boulder, Colorado, USA.

This document supplies the technical material referenced in the main
text. §S1 fully specifies the L2, L3, and L4 dynamical systems used in
the ablation cascade. §S2 gives the *a priori* derivation of the
*h*/*J* scaling of Proposition 1 from the Curie–Weiss free-energy
structure, together with the transient-selection derivation of the
residual (§S2.7). §S3 presents Figure S1 (the L0 mean-field validation). §S4
documents the full convergence test that motivates the d*t* = 0.05
choice in Methods, together with the horizon-stability check, the
provenance of the headline-cell values, the boundary/clipping
integrity analysis (including the boundary-preserving integrator
check and the Feller classification), and the collapse-definition
sensitivity grid. §S5 contains sensitivity and
robustness analyses: noise prescription and *h*(*W*) form (§S5.1),
one-at-a-time parameter sensitivity (§S5.2), alternative bistable
model classes (§S5.3), the wealth-coupling coefficient *h*₀ of the
protective field *h*(*W*) = *h*₀·(*W*/500) required to hold
P(collapse) fixed as coupling varies (§S5.4), the sensitivity of
the wealth equation (§S5.5), the delivered-amplitude
(magnitude-matched) controls (§S5.6), the branch-selection observable
and the collapse-time decomposition of fixed-field failures (§S5.7),
and the selection-formula prediction check together with the
theory-derived amplitude surface (§S5.8). §S6 provides implementation details for
the network agent-based model, the exact reduction of the *N*-agent
system to the one-dimensional SDE under the committed configuration,
and the measured validity range of that reduction. §S7 contains
supplementary Figure S2 (sensitivity sweep across noise prescription
and *h*(*W*) form, supporting §S5.1), Figure S3 (scaling-vs-fixed
control, supporting §S5.1), and Figure S4 (alternative-model-class
robustness, supporting §S5.3 — the asymmetry is not
Curie–Weiss-specific). §S8 contains Supplementary Notes 1–8: full
statements and detailed values moved from the main text during
condensation (prior-work comparison; referents, Definition 1, and
noise notation; delivered-amplitude control details;
amplitude-requirement descriptors; late-channel and ramp detail;
boundary and typing provenance; network detail; and the scope-probe
statement). §S9 reports the pre-registered late-channel mechanism
tests (initial-condition dependence, barrier scaling, and the
waiting-time distribution).

---

## §S1 Full model specifications for the ablation cascade

This section expands the Methods description of Levels 2, 3, and 4 of
the ablation cascade with the full dynamical content. Levels 0 and 1
are completely specified in the main text and are not reproduced here.

### S1.1 Level 2 — minimal model + replicator–mutator beliefs

Level 2 augments the L1 dynamics with a population of *K* = 6 belief
types competing inside each simulated system. Let *f*ₖ(*t*) denote the
fraction of agents holding belief type *k* ∈ {0, …, *K* – 1}, with
∑ₖ*f*ₖ = 1. Each belief type sits at a fixed location *b*ₖ ∈ [–1, 1] on
a one-dimensional opinion axis,

$$
b_k = \frac{2k - (K-1)}{K-1}, \qquad k = 0, \ldots, K-1,
$$

so that the belief axis spans the same range as the coordination
order parameter *m*. The replicator-mutator dynamics evolve *f*ₖ
between simulation steps according to

$$
f_k(t + \Delta t) = (1 - \mu_b)\,\frac{f_k(t)\,W_k(m)}{\langle W \rangle} + \frac{\mu_b}{K},
$$

with mutation rate $\mu_{\mathrm{b}}$ = 0.01 and *m*-coupled fitness

$$
W_k(m) = \exp(\gamma\,b_k\,m), \qquad \langle W \rangle = \sum_j f_j(t)\,W_j(m).
$$

The fitness coupling *γ* = 1 implies that beliefs aligned with the
current sign of *m* are favored. The renormalization against ⟨*W*⟩
keeps ∑ₖ*f*ₖ = 1 to numerical precision; an additional unconditional
renormalization by ∑ₖ*f*ₖ at the end of each step absorbs any drift
from finite-precision arithmetic.

The epistemic-coherence factor is the normalized Shannon entropy of
the belief distribution,

$$
\mathrm{EC}(t) = \frac{-\sum_k f_k(t) \ln f_k(t)}{\ln K} \in [0, 1],
$$

with EC = 1 indicating uniform belief diversity and EC = 0 indicating
complete monoculture. The L2 passive stabilizer field is

$$
h_{\mathrm{eff}}(t) = h(W(t)) \cdot \mathrm{EC}(t),
$$

so the operative restoring force in Equation 1 of the main text is
$h_{\mathrm{eff}}$ rather than *h*. All other dynamics — the SDE for *m*, the
wealth equation, and the collapse criteria — are identical to L1.

Mechanism: at high coupling *J*, the *m*-coupled fitness drives the
belief distribution toward whichever single type aligns with the
saturated *m*; EC drops toward its mutation-limited floor ≈ 0.06 (≈ 0
in pure monoculture). The passive stabilizer's effective field is then strongly
attenuated even though the underlying wealth-driven *h* is intact.
At low *J*, *m* never saturates, fitness selection is weak, and EC
remains high enough that $h_{\mathrm{eff}}$ ≈ *h*.

### S1.2 Level 3 — L2 + Minsky–Keen credit cycle

Level 3 adds a four-phase credit machine to the L2 dynamics. Each
simulated system carries an aggregate debt state *D*(*t*) and a phase
indicator *p*(*t*) ∈ {Expansion, Euphoria, MinskyMoment,
Deleveraging}. Debt evolves continuously,

$$
\dot D = \alpha_D \cdot \max(0, m) \cdot W - \beta_D \cdot D,
$$

with debt-growth coefficient *α*ᴅ = 0.5 and amortization rate
*β*ᴅ = 0.05. Optimism (positive *m*) and high wealth fuel debt
accumulation; deleveraging proceeds at a fixed rate.

Phase transitions are gated by the debt-to-wealth ratio
*r* = *D*/max(*W*, 1):

| From | To | Trigger |
|---|---|---|
| Expansion | Euphoria | *r* > 0.5 ∧ *m* > 0.3 |
| Euphoria | MinskyMoment | *r* > 1.0 |
| MinskyMoment | Deleveraging | timer ≥ 100 steps (5 physical time units at d*t* = 0.05) |
| Deleveraging | Expansion | *r* < 0.1 |

When the MinskyMoment trigger fires, the system's wealth is reduced
by 30% (haircut), its debt is reset to zero, and the SDE noise
amplitude is multiplied by 1.5 for the duration of the MinskyMoment
phase. The phase-dependent income factor scales the wealth equation:

| Phase | Income factor |
|---|---|
| Expansion | 1.00 |
| Euphoria | 1.10 |
| MinskyMoment | 0.50 |
| Deleveraging | 0.85 |

The phase machine and the L2 belief layer evolve concurrently with
the *m*-SDE and the wealth equation. The L3 dynamics are completely
specified by these rules together with the L2 dynamics in §S1.1.

### S1.3 Level 4 — full simulator

L4 is the simulator from prior work, archived in a separate
repository, reachable through the Zenodo archive of this repository
(concept DOI 10.5281/zenodo.20061432). L4 is prior work, not part of this deposit; no main-text claim depends on regenerating it. It comprises the L3
belief and credit dynamics together with the following additional
subsystems (TypeScript implementation):

- **Hawkes-process shock layer** [Bacry et al. 2015] — exogenous
  shock events with self-exciting kernel that perturb both *m* and
  wealth at Poisson-cluster times.
- **Lotka–Volterra predator–prey coupling** between an "innovator"
  and an "imitator" sub-population, modulating the effective income
  multiplier.
- **Demographic dynamics** — birth/death rates linked to wealth and
  age structure on the slow layer.
- **Mortality** with Gompertz–Makeham hazard, coupling demographic
  decline to economic stress.
- **Energy and infrastructure layers** with capacity dynamics that
  modulate consumption and noise.
- **Institutional-quality slow variable** with state-dependent
  capacity bounds, which feeds back into *h* through a fixed
  conversion rule (preserving Definition 1 — the conversion is fixed
  even though institutional quality is dynamic).

The 16-subsystem dependency graph is DAG-validated; integration uses
Strang splitting between fast/medium/slow timescales, with the fast
layer integrated by an Euler–Maruyama / Heun RK2 hybrid. The L4
sweep behind the 30–35% band cited in the main text is documented in
the prior repository, reachable through the Zenodo archive of this
repository (concept DOI 10.5281/zenodo.20061432).

L4 was not re-run at d*t* = 0.05 for this submission; the d*t* = 0.1
results from the prior repository (P(collapse) in the 30–35% band,
rigidity dominating) are consistent with the L1 → L2 → L3 trend at
d*t* = 0.05. A d*t* = 0.05 re-run is a one-day compute task that
would refine the L4 number; it is not necessary for the L1–L3
monotone claim presented in the main text.

---

## §S2 Full derivation of Proposition 1

This section derives the *h*/*J* scaling claim of Proposition 1 from the
mean-field free-energy structure of Equation 1. We work in the
deterministic limit first (§S2.1–§S2.4), then add the multiplicative
noise contribution (§S2.5) and state the noise-corrected form
(§S2.6); §S2.7 then derives the observed residual as transient branch
selection at the unstable origin and validates it against the
simulator.

### S2.1 Effective potential

The deterministic part of Equation 1 can be written as a gradient
flow,

$$
\dot m_{\mathrm{det}} = -\frac{\partial V}{\partial m},
$$

with effective potential (free-energy analogue)

$$
V(m;\,J,h,T) = \frac{m^2}{2} - \frac{T}{J}\ln\cosh\!\left(\frac{Jm + h}{T}\right) + C(J, h, T).
$$

The constant *C* fixes the zero of *V* and plays no role in the
dynamics. Taking the derivative reproduces

$$
-\frac{\partial V}{\partial m} = -m + \tanh\!\left(\frac{Jm + h}{T}\right),
$$

confirming the gradient form. The fixed points of the deterministic
dynamics are the critical points of *V*.

### S2.2 Fixed-point structure

Critical points satisfy

$$
m = \tanh\!\left(\frac{Jm + h}{T}\right). \tag{S1}
$$

For *h* = 0, Equation S1 is the standard Curie–Weiss
self-consistency equation. It has a single stable solution *m* = 0
for *J* < *T* and three solutions for *J* > *T*: two stable branches
at *m* = ±*m*\* and an unstable saddle at *m* = 0, with
*m*\*(*J*, *T*) → 1 as *J* / *T* → ∞.

The bifurcation at *J* = *T* is the Curie–Weiss critical point. We
write *J*ᶜ = *T* below.

### S2.3 Asymmetric splitting under *h* > 0

For *h* > 0 and *J* ≫ *T*, all three critical points persist but
shift. Writing *m* = 1 – ε with ε > 0 small for the upper branch,
substituting into Equation S1, and using the asymptotic
tanh(*x*) = 1 – 2 e⁻²ˣ + *O*(e⁻⁴ˣ) yields

$$
\varepsilon_+ \approx 2\,\exp\!\left(-\frac{2(J + h)}{T}\right). \tag{S2}
$$

The same calculation for *m* = –1 + δ gives

$$
\delta_- \approx 2\,\exp\!\left(-\frac{2(J - h)}{T}\right). \tag{S3}
$$

Both branches are exponentially close to ±1 for *J* ≫ *T*, and the
splitting between ε₊ and δ₋ is exponentially small in 2*h*/*T*. The
saddle, by contrast, sits at *m*ₛ = –*h*/(*J* – *T*) ≈ –*h*/*J* for *J* ≫ *T*; the
linear-in-*h* shift of the saddle along the *m* axis is the
basic geometric fact behind the asymmetry.

### S2.4 Barrier height and the *h*/*J* ratio

The barrier between the two branches is

$$
\Delta V_{\text{barrier}} = V(m_s) - \min(V(m_+), V(m_-)).
$$

For *J* ≫ *T* and *h* = 0, evaluating *V* at the symmetric fixed
points and the saddle gives the standard Curie–Weiss expression

$$
\Delta V_{\text{barrier}} = \frac{1}{2} - \frac{T}{J}\ln 2 + \mathcal{O}\!\left(e^{-2J/T}\right) \xrightarrow{J \gg T} \frac{1}{2}. \tag{S4}
$$

The barrier approaches a finite ceiling of ½ for *J* ≫ *T* (in the
same V-units as the tilt); it does not grow with *J*.

The *h* > 0 perturbation tilts the landscape between the two
branches. From §S2.3, the asymmetric energy splitting between the
protective and unprotective branches is

$$
\Delta V_{\text{tilt}} = V(m_-) - V(m_+) \approx \frac{2 h\,m_*}{J} \approx \frac{2 h}{J} \tag{S5}
$$

for *J* ≫ *T*, since *m*\* → 1. Note that the tilt **scales as *h*/*J***
— the tilt itself decreases with coupling, even though *h* is fixed.

The dimensionless quantity controlling the passive stabilizer's
ability to bias the equilibrium probabilities of the two branches is
the ratio of tilt to barrier:

$$
\eta(J, h, T) \equiv \frac{\Delta V_{\text{tilt}}}{\Delta V_{\text{barrier}}} \approx \frac{2h/J}{1/2} = \frac{4h}{J}\bigl[1 + \mathcal{O}(T/J)\bigr] \xrightarrow{J \gg T} \mathcal{O}\!\left(\frac{h}{J}\right). \tag{S6}
$$

This is Proposition 1's central statement: the tilt-to-barrier ratio
— the deterministic preferential pull toward the protective branch
relative to the barrier separating the branches — **scales as
*h*/*J*** at high coupling. Within the committed parameterization,
*h*(*W*) = 2(*W*/500) is bounded above by the equilibrium wealth,
itself bounded by the static parameter *μ* through Equation 2 of the
main text, so *η* → 0 as *J* → ∞ for any fixed *μ*. This is an
asymptotic property of the ratio *η* within the model; what it
implies for collapse probabilities at *finite* coupling is an
empirical question, answered by the measured comparison of the
specified parameterizations (main text; §S5.1): there, fixed fields
*h* = 2.4 and 5.0 hold P(collapse) ≤ 0.001 at *μ* = 100 throughout
the calibrated range *J* ≤ 5 (0.00 at *n* = 100; 0.001 at
*n* = 1,000, §S5.7), so boundedness alone does not determine
finite-*J* outcomes. (In the static sweeps every wealth-fed variant delivers a
field profile built on the base 2(*W*/500), differing in
delivered field-amplitude *profile*; whether any single summary of that
profile orders the outcomes is examined in §S5.6 as an exploratory,
domain-dependent question, not advanced as a result.)

### S2.5 Multiplicative-noise correction

The full SDE in Equation 1 of the main text adds a multiplicative
noise term

$$
\xi\,\sqrt{1 - m^2}\,\zeta(t),
$$

(*ζ*(*t*): the white noise of Equation 1 of the main text, written
*ζ* here — as in §S2.7 — to avoid clash with the tilt-to-barrier
ratio *η* of §S2.4), which vanishes at the saturated branches
*m* = ±1. This has two
effects.

First, the equilibrium measure of the SDE is **not** the Boltzmann
measure of *V* alone. The Itô-form Fokker–Planck equation,

$$
\partial_t P = -\partial_m\!\left[\left(-\partial_m V\right) P\right] + \frac{1}{2}\partial_m^2\!\left[\xi^2 (1 - m^2) P\right],
$$

has steady-state solution proportional to (after multiplying out)

$$
P_{\mathrm{ss}}(m) \propto \frac{1}{\xi^2(1 - m^2)}\exp\!\left(-\frac{2}{\xi^2}\,U_{\mathrm{eff}}(m)\right), \qquad U_{\mathrm{eff}}(m) = \int^m \frac{\partial_m V(m')}{1 - {m'}^2}\,dm'.
$$

Near the saturated branches the prefactor 1/(1 – *m*²) diverges,
reflecting the fact that the noise is suppressed there and the
trajectories spend exponentially long times in the immediate
neighborhood of ±1. The *exponentially weighted* probability mass
on each branch is governed by $U_{\mathrm{eff}}$, not *V*.

Second, the steady-state preference between branches is controlled by
the *difference* of the effective potential, ΔU = $U_{\mathrm{eff}}$(*m*₋) −
$U_{\mathrm{eff}}$(*m*₊). Rather than asserting (as an earlier draft did) that the
multiplicative-noise prefactors "cancel" to leading order, we compute
this difference explicitly (parameter-free; script
`simulation/scripts/phase1b_noise_leg.py`). Two facts hold to leading
order for *J* ≫ *T*. (i) *The multiplicative-noise tilt equals the
deterministic tilt.* Quadrature of ΔU between the interior wells gives
ΔU/ΔV → 1 and *J*·ΔU/*h* → 2 as *J* grows (both ΔU and ΔV ≈ 2*h*/*J*),
so the *h*/*J* scaling is **prescription-independent**, with an
*O*(*T*/*J*) remainder (e.g. ΔU/ΔV = 1.33, 1.08, 1.01 at *J* = 4, 8, 20).
(ii) *The wells are exponentially pinned at ±1.* The field-sensitivity
of the well location d*m*₊/d*h* tracks $e^{-2J/T}$ (verified: 4.4 ×
10⁻² → 3.7 × 10⁻⁹ as *J* = 4 → 20, matching $e^{-2J/T}$), so the
multiplicative prefactor and the boundary-layer concentration are
*h*-insensitive up to $O(e^{-2J/T})$ and **cannot** enter the
leading *h*/*J* tilt — replacing the hand-waved cancellation with this
explicit suppression. What we do **not** establish is a rigorous error
bound on the full Itô stationary *occupation ratio*, whose *O*(1)
Kramers prefactor near ±1 remains numerically delicate (the quadrature
ratio drifts ≈ 13%–45% over *J* = 3.5–6 as the 1/(1 − *m*²) factor
enters); only the leading *h*/*J* exponent of the inter-branch tilt is
secured. For this reason the result below is stated as a **Proposition**
— asymptotic *h*/*J* scaling of the inter-branch tilt, regime *J* ≫ *T*,
leading order — rather than a theorem.

### S2.6 Statement of Proposition 1 and what *h*/*J* governs

Combining §S2.4 and §S2.5: the dimensionless tilt-to-barrier ratio
between the protective branch and the unprotective branch, in the
deterministic mean-field limit,

$$
\eta(J, h, T) = \frac{4 h}{J}\bigl[1 + \mathcal{O}(T/J) + \mathcal{O}(\xi^2)\bigr],
$$

vanishes as *J* → ∞. This is the formal content of Proposition 1.

It is important to specify what physical observable this ratio
governs. *η* is the ratio of the energy *tilt* between branches
($\Delta{}V_{\mathrm{tilt}}$ ≈ 2*h*/*J*) to the *barrier* separating them
($\Delta{}V_{\mathrm{barrier}}$ → ½ for *J* ≫ *T*, in the same V-units as the tilt). Two things follow.

(i) *Landscape topology.* When *η* approaches unity from below, the
tilt is comparable to the barrier and one of the two wells flattens
or disappears. *η* → 0 means the two wells survive and become
energetically symmetric: at fixed *h*, the wrong-branch attractor
persists at sufficiently high coupling, because eliminating the
second well requires a tilt comparable to the barrier (*η* ≳ 1).
This is a statement about the landscape topology of the potential,
not by itself a statement about measured collapse rates: in the
simulations, the fixed fields *h* = 2.4 and 5.0 hold
P(collapse) ≤ 0.001 at *μ* = 100 throughout the calibrated range
*J* ≤ 5 (0.00 at *n* = 100; 0.001 at *n* = 1,000, §S5.7; §S5.1), and the stress form (the one variant whose field responds within a run to *m*) and the
coupling form suppress the residual to ≤ 0.05 at *α* ≥ 0.5 and
eliminate it at the tested *α* = 1.0 and 2.0, while the quadratically coupling-responsive form
retains a 0.14 residual at *α* = 1.0 (§S5.1) — an unbounded
functional form alone is not sufficient. In static sweeps the coupling and quadratic
variants deliver per-*J* constant rescalings of the wealth-fed base
field 2(*W*/500), while the stress form adds a within-run state response
(α·max(0, −*m*)·*J*); §S5.6's exploratory ordering test of whether any
single summary of the delivered field orders the variants is
domain-dependent and yields no refutation under the domain rule the
paper adopts. At the fully-stabilized headline
coefficient *α* = 2.0 the stress variant and a matched constant are
both at P = 0.00 — an all-zero consistency check that cannot
discriminate form from magnitude. At an informative operating point
(*α* = 0.05, §S5.6, control E2) a constant matched to the active form's
full-trajectory delivered mean collapses far *less* than the active
form itself (verdict REFUTED in all three cells), because for constant and rescaled fields, at fixed *J*, the
full-trajectory delivered mean is not a sufficient statistic for the
outcome (E2) — so for the coupling and quadratic variants the
suppression thresholds above are statements about delivered
field-amplitude profiles during the selection window, while the
stress form's suppression is carried by its rectified response
(§S5.8), not by within-run responsiveness of a rescaled field.

(ii) *Stationary branch occupation.* In Kramers theory
[Hänggi et al. 1990] the ratio of
occupation times between branches is exp($\Delta{}V_{\mathrm{tilt}}$/σ²) where σ² is
the effective noise variance. Substituting $\Delta{}V_{\mathrm{tilt}}$ = 2*h*/*J* gives
an exponent ∝ *h*/*J* — the same *h*/*J* scaling carried by the
tilt-to-barrier ratio *η* of §S2.4. Both the deterministic ratio and
the additive-noise Kramers exponent therefore share the same leading
*h*/*J* dependence. For weak noise (*ξ* small) any tilt
above σ² produces strong occupation preference; the effect of high
coupling on the steady-state population is mediated through the
*h*/*J* exponent in this regime.

What our simulation reports — the residual collapse rate at high *J*
— is determined primarily by *initial branch selection* under
multiplicative noise that vanishes at saturation, rather than by
asymptotic stationary occupation. Initial selection is symmetric
when basins are narrow (high *J*), and once an unlucky trajectory
lands on the −*m* branch the noise $\xi \sqrt{1-m^2}$ prevents escape.
The asymmetry shrinks at high coupling because the origin-instability
rate *λ* grows with *J* while the transient bias *c* is pinned at the
selection-window field, leaving the bias less time to act (§S2.7). The scaling of the
*observed* residual collapse rate with *J* is derived in §S2.7 as
transient branch selection at the unstable origin and characterized
empirically in §S5; it does *not* share η's leading exponent.

Both the deterministic tilt-to-barrier ratio and the additive-noise
Kramers occupation exponent carry the same leading *h*/*J* dependence;
the corrected derivation therefore does not support the earlier
suggestion of a prescription-dependent exponent (*h*/*J*² versus
*h*/*J*). The multiplicative-noise prescription shifts the *O*(1)
prefactor and the initial-branch-selection weighting, not this
leading exponent. The qualitative finding — that the inter-branch
tilt of a fixed field shrinks as *h*/*J* while the barrier stays
*O*(1) — is robust to the noise prescription; how this translates
into measured collapse probabilities for specific parameterizations
is quantified in §S5.1 and §S5.4.

### S2.7 The residual as selection at the unstable origin (Suzuki transient selection)

§S2.6 established that the *observed* residual is set by **initial
branch selection** rather than stationary occupation, and deferred its
*J*-dependence to empirics. Here we derive that dependence in closed
form and validate it against the simulator.

**Regime (precondition).** This subsection applies in the common-mode
regime, in which the order parameter carries aggregate *O*(1) noise.
Under purely idiosyncratic noise with large *N* the macroscopic
fluctuation is ∝1/$\sqrt{N}$, no escape occurs, and the residual vanishes
(§S6); the *J*-dependence below is conditional on the common-mode
regime, not independent of it.

**Why selection, not Kramers escape.** At the good well the wealth
feedback drives *W* to its high fixed point (*W* → 3500 at *μ* = 100),
so *h* = 2(*W*/500) → 14 and *η* = 4*h*/*J* is large; steady-state
escape from the good well is negligible. Collapse is decided in the
**initial transient** from (*m* = 0, *W* = 100), where the
transient field is $h_{\mathrm{eff}}$ = *h*(*W*₀) = 0.4.
Numerically the median collapse occurs within 1.4–1.7% of the run horizon
across *J* — collapses are early-transient, confirming selection.

**Linearization and the wrong-branch probability.** For *h* > 0 the
point *m* = 0 is not a fixed point: *f*(0) = tanh($h_{\mathrm{eff}}$/*T*) > 0.
Linearizing the drift,

$$ \dot m = \lambda m + c + \xi\,\zeta(t), \quad c = \tanh(h_{\mathrm{eff}}/T), \quad \lambda = -1 + (J/T)\,\mathrm{sech}^2(h_{\mathrm{eff}}/T), $$

with ζ the common-mode unit white noise (the model's multiplicative
noise $\xi \sqrt{1-m^2}$ reduces to additive ξ at the origin to leading order
in *m*, where the branch is selected). For *J* > *T*·cosh²($h_{\mathrm{eff}}$/*T*)
≈ 2.08 the origin is **linearly unstable** (λ > 0) — Suzuki's
selection at an unstable point [Suzuki 1977; McKane & Tarlie 2001], a
mechanism experimentally confirmed in the transient decay of optically
unstable states [Arecchi & Politi 1980]. Integrating, the eventual
sign of *m* (good + vs collapse − branch) is fixed by the effective
seed *Z* = *c*/λ + 𝒩(0, ξ²/2λ), so the collapse probability is

$$ P_{\mathrm{coll}}(J) = \Phi\!\left(-\frac{c}{\xi}\sqrt{\frac{2}{\lambda(J)}}\right), $$

with Φ the standard normal CDF. Three consequences: (a) **ceiling → ½**
as *J* → ∞ (λ → ∞ ⇒ argument → 0): an infinitely fast instability
leaves the bias no time to act, so selection becomes a coin flip and
the closed-form collapse probability saturates at a ½ ceiling (an
asymptotic statement; see item (b) of the deepened characterization
below for what the simulated range can and cannot pin);
(b) **monotone**
increase in *J*; (c) **no power law** — *P* is a Gaussian CDF of a
function $\propto \lambda(J)^{-1/2}$, λ affine in *J*, so any power-law fit returns
a window-dependent effective exponent.

**Numerical validation (parameter-free).** With *c*, λ(*J*), ξ fixed
by the committed constants (*T* = 2, ξ = 0.5, $h_{\mathrm{eff}}$ = 0.4) — no fitted
parameters — the prediction matches the simulated passive residual at
**R² = 0.93, MAE = 0.024** across *J* ∈ [2.5, 20] (a conservative
window: it includes the low-*J* points that honesty flag (i) below
notes are over-predicted, so the reported R² understates the agreement;
the fit is near-exact for *J* = 4–7). The pre-registered α-fit's target
was the efficacy ratio E(*J*) = 1 − P(*J*, 100)/P(*J*, 40), fitted
against (*J* − 2)^(−*α*): point estimate 0.63 (bootstrap SE 0.10,
95% CI [0.53, 0.93]); its pre-declared SE ≤ 0.05 stopping target was
not reached at the 800-seed cap, and its pre-declared *n* = 5
subsample fit had CI [−0.46, 2.38], spanning both ½ and 1 (protocol
archived at `simulation/prereg/PREREG_alpha_fit_2026-05-30.md`). The
closed form's window-dependent local slope — ≈2.2 in a low-*J* window,
≈0.35 at high *J*, ≈0.69 over the full range — explains why neither
E(*J*) nor P(*J*, 100) admits a single exponent; the full-range ≈0.69
and the measured 0.63 are slopes of *different* functions, not a
numerical reproduction of one by the other (on the measurement grid
the log–log slope of P(*J*, 100) is +1.26). The α-fit thus
measured the local slope of a saturating selection CDF, not a power
law; this **connects** the measured residual to the same model —
accounting for the §S2.6 gap between the theory and its observable —
while making explicit that the residual is governed by the transient
instability rate λ, distinct from the stationary η = *h*/*J* (it does
*not* share η's exponent, §S2.6). This separation between a transient
rate that governs the response to a shock and the stationary barrier
that governs steady-state stability mirrors the stability–resilience
distinction documented for early-warning indicators in bistable systems
[Dai et al. 2015]. Validation script:
`simulation/scripts/suzuki_selection_validation.py`.

**Suzuki↔Kramers crossover.** The two pictures are separated in *time*:
early transient selection (Suzuki) sets the residual, while the late
deep-well regime (Kramers over a large *η*-barrier) contributes
negligibly. This crossover, with the near-symmetric small-*h* limit
where Suzuki scaling carries logarithmic corrections, further smears
any single exponent.

**Honesty flags.** (i) The linearization holds for λ not too small; as
*J* → *T* (λ → 0) the selection time diverges and the leading-order
formula over-predicts (e.g. *J* = 3), so we do not extrapolate below
*J* ≈ 3.5. (ii) The bias uses the transient field $h_{\mathrm{eff}}$; *W* (hence *h*)
grows modestly during the brief selection window, which is why slow
(low-*J*) selection is mispredicted while fast (high-*J*) selection is
near-exact — an approximation we state rather than hide. (iii) This is
Suzuki's 1977 selection theory together with Kramers escape *applied*
to a controller-failure setting, not new mathematics.

**Deepened characterization of the transient-selection mechanism.** Three consequences of the closed form $P_{\mathrm{coll}}$(*J*) = Φ(−(*c*/*ξ*)$\sqrt{2/\lambda}$) sharpen "no clean exponent" from an observation into a mechanistic statement (parameter-free; script `simulation/scripts/phase1a_suzuki_deepen.py`). (a) *Continuous local exponent.* The local slope $a_{\mathrm{loc}}$(*J*) = d ln $P_{\mathrm{coll}}$/d ln *J* is not constant: it sweeps **6.40 → 0.10** across *J* = 2.5 → 20 with no plateau, equaling **0.73** at *J* ≈ 4. A windowed power-law fit therefore returns whatever $a_{\mathrm{loc}}$ averages over its window — which is exactly why the pre-registered fit returned an ambiguous exponent (point estimate **0.63**, 95% CI [0.53, 0.93]; its window centers near *J* ≈ 4): there is no power law to measure, only the local slope of a saturating selection CDF. (b) *Approach to the ½ ceiling* (closed form): (½ − $P_{\mathrm{coll}}$) $\propto J^{-1/2}$ as *J* → ∞ (log-log slope −0.50; $(\tfrac{1}{2} - P_{\mathrm{coll}})\sqrt{J} \to 0.32 = \varphi(0)\,(c/\xi)\,\sqrt{2T/\operatorname{sech}^2(h_{\mathrm{eff}}/T)}$). This is a *closed-form* property: the 800-seed simulator reaches only *J* = 50 — far below the *J* ~ 10³ onset of the asymptote — so the (noisy, non-monotone) simulated residual is consistent with approach to ½ but does **not** pin the −1/2 exponent. (c) *Noise structure sets the failure mode, not the selection probability.* Because the branch is chosen at the origin where $\sqrt{1-m^2}$ → 1, the selection *probability* is noise-structure-invariant **only for fast selection** (*J* ≥ 6: additive vs multiplicative agree to ≈ 0.02–0.04); near threshold (*J* ≈ 3–4) slow selection lets *m* wander where $\sqrt{1-m^2}$ < 1 and the two differ by a factor ≈ 2–3. What the noise structure robustly fixes is the *failure type* (multiplicative noise vanishing at ±1 locks trajectories in → rigidity; additive noise admits fragmentation; §S5.1), not the leading selection probability.

---

## §S3 Supplementary Figure S1 — L0 mean-field validation

**Supplementary Figure S1.** *Mean-field equilibrium curves for the bare Curie–Weiss
SDE (Level 0).* The figure shows the long-time mean of independent
SDE trajectories at each (*J*, *h*) cell, with *h* held fixed at
three values (*h* ∈ {0.0, 0.2, 0.4}) and *J* swept across
{0.5, …, 5.0}. Dashed curves are the analytical mean-field fixed
points obtained by solving Equation S1 numerically; solid colored
curves are the SDE long-time means with shaded bands showing the
interquartile range across seeds. For *h* > 0 (orange and red), the SDE cross-seed mean tracks the upper
mean-field branch at low coupling but falls increasingly below it as
*J* rises (e.g. at *h* = 0.2 the mean is ≈ 0.23 against an upper branch
≈ 0.99 by *J* ≈ 4.75). This is not a breakdown of the mean-field
description but the finite-noise residual-selection effect central to
this paper: a fraction of seeds are noise-selected onto the lower (–*m*)
branch (the lower quartile reaches ≈ –0.95 at high *J* while the upper
quartile continues to track the upper branch), dragging the cross-seed
average down. For *h* = 0 (blue), the SDE mean fluctuates
around zero — an expected finite-sample effect under spontaneous
symmetry breaking, where individual seeds lock onto either the +1
or –1 branch and the cross-seed average cancels — while the
mean-field curves plotted show both symmetric branches.
This validates the analytical backbone of Proposition 1: the mean-field
free-energy structure used to derive the *h*/*J* scaling in §S2 is
the correct effective description of the bare SDE in the absence of
the wealth–field coupling, modulo the finite-sample symmetry-breaking
caveat at *h* = 0 which does not affect Proposition 1 (which applies for
*h* ≠ 0). Source file:
`simulation/results/ablation/level0_equilibrium_curves.png`; rendered
here as Supplementary Figure S1, `manuscript/arxiv/figures/figS1.png`.

The multiplicative noise term $\xi \sqrt{1-m^2}$ vanishes at the saturated
branches and does not distort the *branch structure* (the two wells
remain at ≈ ±1); the high-*J* downward drift of the cross-seed mean is
the §S2.7 residual-selection effect, not a numerical artifact of the
d*t* = 0.1 integration or the clipping rule
*m* ∈ [–1 + ε, 1 – ε] employed in the integrator.

---

## §S4 Convergence-check details

The convergence test reported in Methods used two diagnostic cells —
the residual-rigidity corner (*J* = 5.0, *μ* = 100) and the
fragmentation-cured corner (*J* = 0.5, *μ* = 100) — at three step
sizes d*t* ∈ {0.1, 0.05, 0.025}, with the simulation horizon held
constant at 1,000 physical time units (10,000, 20,000, and 40,000
steps respectively). The full results, archived at
`simulation/results/minimal_model/convergence_check.csv`:

| Cell | *J* | *μ* | d*t* | *N*ₛₜₑₚₛ | P(collapse) | *n*ᶜᵒˡˡ | rigidity-share of collapsed |
|---|---|---|---|---|---|---|---|
| A — residual rigidity | 5.0 | 100 | 0.100 | 10,000 | 0.250 | 25 | 0.880 |
| A — residual rigidity | 5.0 | 100 | 0.050 | 20,000 | 0.350 | 35 | 0.886 |
| A — residual rigidity | 5.0 | 100 | 0.025 | 40,000 | 0.350 | 35 | 0.886 |
| B — fragmentation cured | 0.5 | 100 | 0.100 | 10,000 | 0.000 | 0 | — |
| B — fragmentation cured | 0.5 | 100 | 0.050 | 20,000 | 0.000 | 0 | — |
| B — fragmentation cured | 0.5 | 100 | 0.025 | 40,000 | 0.000 | 0 | — |

Each cell uses 100 independent seeds. P(collapse) at the rigidity
corner shifts upward by ~10 percentage points between d*t* = 0.1 and
d*t* = 0.05, then plateaus at d*t* ≤ 0.05. The Euler–Maruyama scheme
has weak-order $O(\mathrm{d}t)$ bias and strong-order $O(\mathrm{d}t^{1/2})$ bias for
SDEs; the observed ~10-percentage-point shift between d*t* = 0.1 and
d*t* = 0.05 is consistent with the $O(\mathrm{d}t)$ leading-order weak bias; at
that rate the expected further shift between d*t* = 0.05 and
d*t* = 0.025 is ≈5 percentage points, and the absence of an observed
shift indicates the residual $O(\mathrm{d}t)$ bias is at or below the
≈7-percentage-point sampling standard error of the unpaired 100-seed
comparison. We therefore adopt
d*t* = 0.05 as the coarsest step at which the high-*J* residual is
converged within sampling SE.

The rigidity share of the residual is essentially d*t*-stable at
0.88–0.89 across all three step sizes (all shares in this table are
clipped-EM |*m*| = 0.9 gate values; the boundary-integrity block
below shows this statistic is integrator-scheme-sensitive — 0.88
clipped EM versus 0.74 Lamperti at the same gate — while P(collapse)
and timing are robust); the qualitative result —
rigidity dominance of the high-*J* residual — is therefore robust to
integrator refinement, even though the exact P(collapse) headline
required the d*t* = 0.05 rerun.

**Band-vs-cell distinction.** The convergence-check numbers
(0.25 / 0.35 / 0.35) are at the single cell *J* = 5.0, *μ* = 100. The
Results headline of 0.23 for the high-coupling band is the average
of P(collapse) across *J* ∈ {4.0, 4.5, 5.0} at *μ* = 100, which
includes two cells (*J* = 4.0 and *J* = 4.5) where P(collapse) is
0.19 and 0.22 respectively — both below the *J* = 5.0 cell value of
0.28. For comparison of the two numbers: the *single-cell* converged
residual is 0.28; the *band-averaged* converged residual reported as
the paper's headline is 0.23.

The convergence-check script uses an independent seed set from the
main sweep (the convergence script seeds deterministically per cell
and step size, by `1009·J·10 + μ + ⌊1/dt⌋`; because the seed
depends on d*t*, the noise paths are independent — not shared —
across d*t* values, and the comparison across step sizes is unpaired.
The main sweep seeds by
`seed_base + j_idx · 1000003 + m_idx · 1009`). The difference
between the convergence-check single-cell value (0.35) and the
main-sweep single-cell value (0.28) at the same (*J*, *μ*, d*t*) is
within two combined 100-seed standard errors (combined SE ≈ 0.065 at
*p* ≈ 0.3) and reflects normal stochastic variation across independent seed
sets, not a numerical inconsistency.

**Provenance of the headline-cell P(collapse) values.** The
(*J* = 5, *μ* = 100) residual is quoted at several nearby values across
the paper because they are *different measurements*, not one number
reported inconsistently; all single-cell values are mutually
consistent with independent-seed binomial sampling variation
(per-estimate SE ≈ 0.05; all pairwise differences ≤ 2 combined SEs):

| Value | Quantity | Source / configuration |
|---|---|---|
| 0.23 | high-*J* *band* average (not single-cell) | mean over *J* ∈ {4.0, 4.5, 5.0}, *μ* = 100, minimal model L1, d*t* = 0.05 (Figure 1, abstract headline) |
| 0.28 | single cell, scalar model | minimal/active-stabilizer passive, 100 seeds, d*t* = 0.05, main sweep |
| 0.33 | single cell | network ABM 100-seed headline (§S6) and the multiplicative-noise baseline of the §S5.1 sensitivity sweep (independent seed set) |
| 0.35 | single cell | convergence check at d*t* = 0.05/0.025 (independent seed set; §S4 above) |
| 0.41 | single cell | the unperturbed baseline of the §S5.2 one-at-a-time sensitivity analysis (OAT seed set, independent of the main sweep; `oat_sensitivity.csv`) |

The paper's headline is the band-average **0.23**; the converged
single-cell scalar value is **0.28**. The spread among single-cell
values reflects noise prescription and independent seed sets — all
within two combined SEs — not a contradiction. (The network ABM at this
cell is dynamically identical to the scalar model under the shared
per-step shock, §S6; its 0.33 versus the main-sweep 0.28 is
independent-seed sampling variation, binomial SE ≈ 0.046 at
*p* ≈ 0.3.)

**Horizon stability.** To verify that the reported collapse
probabilities are not truncation artifacts of the *t* = 1,000 horizon,
we extended three cells (*J* ∈ {2.5, 3.5, 5.0} at *μ* = 100) to
*t* = 4,000 (80,000 steps at d*t* = 0.05) with *n* = 4,000 independent
trajectories per cell, using the same integrator and collapse
criterion as the main sweep. Before extending, the script reproduces
the archived main-sweep cells per-seed exactly (identical collapse
flags and collapse steps under the standard cell seeds), so the
extension is a continuation of the same dynamics, not a
reimplementation. The cumulative collapse fraction is identical at
horizons *t* = 250, 500, 1,000, 2,000, and 4,000 in every cell: not a
single additional trajectory collapses beyond *t* = 250.

| *J* | *μ* | *n* | P(collapse), *t* = 250 | *t* = 500 | *t* = 1,000 | *t* = 2,000 | *t* = 4,000 | 95% CI (*t* = 4,000) |
|---|---|---|---|---|---|---|---|---|
| 2.5 | 100 | 4,000 | 0.0568 | 0.0568 | 0.0568 | 0.0568 | 0.0568 | [0.0498, 0.0640] |
| 3.5 | 100 | 4,000 | 0.2097 | 0.2097 | 0.2097 | 0.2097 | 0.2097 | [0.1968, 0.2225] |
| 5.0 | 100 | 4,000 | 0.2945 | 0.2945 | 0.2945 | 0.2945 | 0.2945 | [0.2805, 0.3090] |

Collapse times saturate in absolute units rather than tracking the
horizon: the latest collapse across all 12,000 trajectories occurs at
*t* = 34.6 (medians 14.6–16.5, 99th percentiles 20.7–29.8 depending on
*J*), consistent with collapse being an early first-passage event from
the transient following initialization, after which surviving
trajectories are absorbed onto the stabilized branch. Bootstrap 95%
CIs (10,000 path resamples) are horizon-independent to the resolution
shown. One archived main-sweep trajectory collapses at *t* = 943.5
(94% of the *t* = 1,000 horizon); that event belongs to the
(*J* = 5.0, *μ* = 40) cell — a low-*μ* regime outside the *μ* = 100
cells examined here — where late collapses arise from the weaker
income margin, not from horizon truncation of the *μ* = 100
measurement. Within the *μ* = 100 cells of the archived sweep, the
latest collapse occurs at *t* = 33.6. Seeds: cell seed = 0xC0DE +
*j*ᵢ · 1,000,003 + *m*ᵢ · 1,009 + 7,777,777 (horizon-stability offset);
script: `simulation/scripts/horizon_stability.py`; results:
`simulation/results/horizon_stability/`.

**Boundary and clipping integrity.** The multiplicative noise
$\xi \sqrt{1-m^2}$ vanishes at *m* = ±1, and the Euler–Maruyama (EM)
integrator handles the boundary by projecting any proposed update that
leaves [−1 + 10⁻⁶, 1 − 10⁻⁶] back onto that interval ("clipping").
Because the fixed-point analysis of §S2 places the wells exponentially
close to ±1, we quantified how often the clip engages and whether any
reported quantity depends on it (script:
`simulation/scripts/boundary_integrity.py`; results:
`simulation/results/boundary_integrity/`). Clip-hit statistics
(*n* = 1,000 trajectories per cell; a "clip hit" is a step whose
proposed EM update falls outside the interval and is projected back;
fractions are over all trajectory-steps):

| cell | *J* | *h* | steps clipped, upper | steps clipped, lower | trajectories ever clipped, upper | trajectories ever clipped, lower |
|---|---|---|---|---|---|---|
| headline passive | 5 | passive, *W*-coupled | 0.3526 | 0.0221 | 0.966 | 0.296 |
| fixed | 5 | 2.4 | 0.3298 | 0.0000 | 1.000 | 0.006 |
| fixed | 5 | 5.0 | 0.3959 | 0.0000 | 1.000 | 0.000 |
| fixed | 10 | 2.4 | 0.3849 | 0.0136 | 0.991 | 0.117 |
| fixed | 10 | 5.0 | 0.4014 | 0.0000 | 1.000 | 0.005 |
| fixed | 15 | 2.4 | 0.3598 | 0.0422 | 0.946 | 0.193 |
| fixed | 15 | 5.0 | 0.3989 | 0.0028 | 1.000 | 0.028 |
| fixed | 20 | 2.4 | 0.3379 | 0.0642 | 0.907 | 0.226 |
| fixed | 20 | 5.0 | 0.3881 | 0.0140 | 0.990 | 0.076 |

Clipping is pervasive — roughly 33–40% of all steps are projected at
the upper edge — so we checked the reported statistics against an
integrator that never clips. The correct constant-diffusion (Lamperti)
transform for this SDE is *u* = arcsin(*m*)/*ξ*, giving (Itô)
d*u* = [(*f* + *ξ*²*m*/2)/($\xi \sqrt{1-m^2}$)] d*t* + d*W* with
*m* = sin(*ξu*): *m* stays in [−1, 1] identically, and *u* is folded
back into [−π/(2*ξ*), π/(2*ξ*)] by exact reflection (sin is
invariant), realizing the instantaneously reflecting behavior of the
regular boundary (classification below); the integrable drift
singularity at the boundary is tamed by a trust-region cap
|drift · d*t*| ≤ 0.5 (engaged on 5.3% of trajectory-steps; folds on
20.0%). The seemingly natural alternative *u* = artanh(*m*) does
**not** work: its diffusion coefficient is *ξ* cosh(*u*) — not
constant — and because *m* = ±1 is attainable (below), *u* reaches ±∞
in finite time and the scheme explodes (empirically, paired test
trajectories diverged early in the run, so the artanh scheme is
unusable; the archived CSVs retain only the adopted Lamperti scheme's
diagnostics).
Paired comparison at the headline cell (*J* = 5, *μ* = 100, *n* = 2,000,
identical per-step d*W* draws): P(collapse) = 0.3070 (clipped EM)
versus 0.2950 (Lamperti), paired difference +0.0120 [+0.0065, +0.0175];
median collapse step 292 versus 295 (paired difference −5.7 steps =
0.28 physical time units on both-collapsed pairs); p99 collapse step
425 versus 447; discordant pairs 28 (EM-only) versus 4 (Lamperti-only).
P(collapse) and collapse timing are therefore robust to the boundary
treatment. The one sensitive statistic is the rigidity-versus-mixed
split at the |*m*| = 0.9 gate: rigidity share of collapsed 0.881
(clipped EM) versus 0.744 (Lamperti). Both schemes place the
collapse-time |*m*| mass near the ordered state (fragmentation share
≤ 0.2% in both: 0.0000 clipped EM, 0.0017 Lamperti,
`lamperti_summary.csv`), but the clip pins
trajectories at exactly 1 − 10⁻⁶ while the reflected scheme scatters
them over ≈ [0.97, 1], moving ≈ 14% of collapsed trajectories
(rigidity share 0.881 → 0.744) just below the 0.9 gate. The rigidity share is therefore a
scheme-dependent number — 0.88 under clipped EM, 0.74 under Lamperti,
at the same 0.9 gate — and every rigidity-typing statement in this
paper carries that caveat; the qualitative classification (collapse
occurs in an ordered state, essentially never a fragmented one) is
invariant.

The Feller boundary classification explains why the clip engages at
all. σ²(*m*) = *ξ*²(1 − *m*²) vanishes linearly at the boundary; with
*y* the distance to the boundary the near-boundary process is the
CIR-type diffusion d*y* = *c* d*t* + $\xi \sqrt{2y}$ d*W*, with inward
drift *c* = 1 − tanh((*J* + *h*)/*T*) at *m* = +1 and
*c* = 1 + tanh((−*J* + *h*)/*T*) at *m* = −1. For *c* < *ξ*² the
boundary is **regular** — attainable in finite time, not absorbing,
with the clip/reflection choice fixing instantaneous reflection — and
for *c* ≥ *ξ*² it is an **entrance** boundary (unattainable). At
*T* = 2, *ξ* = 0.5, *J* = 5, both boundaries are regular for the
headline/passive field range (e.g. *h* = 0.4: δ = 2*c*/*ξ*² = 0.072 at
*m* = +1, 0.159 at *m* = −1) and for *h* = 2.4; *m* = −1 becomes an
entrance boundary only for *h* ≳ *J* − 1.95 (at *J* = 5 this means
*h* ≳ 3.05, so the fixed *h* = 5.0 cell at *J* = 5 has an unattainable
lower boundary — consistent with its 0.0000 lower-clip fraction in
the table above — whereas at the extended-*J* cells *h* = 5.0,
*J* = 15/20 the lower boundary is again regular, consistent with the
nonzero lower-clip fractions 0.0028 and 0.0140 measured there).
No boundary is absorbing in any studied regime, so
clipping approximates the correct reflecting behavior, with the
paired comparison above quantifying the residual discretization
effect.

Because the extended-*J* fixed-field cells (*J* > 5) have both
substantial lower-boundary clipping and a late-collapse tail (§S5.7),
we re-ran those four cells under both integrators with paired noise
(*n* = 2,000 per cell; script:
`simulation/scripts/late_channel_check.py`; results:
`simulation/results/boundary_integrity/late_channel/`). "Late" means
collapse time *t* > 100:

| *h* | *J* | scheme | P(collapse) | late fraction of collapses | median *t* | p99 *t* |
|---|---|---|---|---|---|---|
| 2.4 | 15 | clipped EM | 0.183 | 0.205 | 14.5 | 967.0 |
| 2.4 | 15 | Lamperti | 0.293 | 0.515 | 128.7 | 981.2 |
| 2.4 | 20 | clipped EM | 0.247 | 0.257 | 14.5 | 955.9 |
| 2.4 | 20 | Lamperti | 0.393 | 0.520 | 138.6 | 968.0 |
| 5.0 | 15 | clipped EM | 0.0285 | 0.193 | 14.1 | 947.9 |
| 5.0 | 15 | Lamperti | 0.0390 | 0.449 | 18.7 | 991.5 |
| 5.0 | 20 | clipped EM | 0.0715 | 0.315 | 14.5 | 948.3 |
| 5.0 | 20 | Lamperti | 0.1115 | 0.525 | 160.5 | 964.9 |

The late channel is not a clipping artifact — it **survives, amplified**,
under the boundary-preserving scheme: the paired differences in
P(collapse) (EM minus Lamperti) are −0.110 [−0.126, −0.094],
−0.146 [−0.163, −0.129], −0.0105 [−0.0164, −0.0046], and
−0.040 [−0.051, −0.029] for the four cells in table order, all
excluding zero, and the paired differences in P(late collapse) are of
the same size. Concretely, the clipped-EM extended-*J* fixed-field
collapse probabilities **understate** the boundary-preserving values by
factors of 1.60, 1.59, 1.37, and 1.56 (≈ 1.5–1.6× at three of the four
cells) — the clip suppresses the noise-driven late-escape route by
pinning trajectories at the boundary. For this reason, every
fixed-field P(collapse) value at *J* > 5 quoted in this SI is labeled
as a clipped-EM value with this check cited; the clipped-EM numbers
are conservative (low-side) readings of the fixed-field failure
probabilities.

**d*t*-convergence of the late channel.** The four cells were re-run
at d*t* = 0.025 (40,000 steps, *n* = 1,000 per scheme×cell) under
both integrators
(`simulation/results/boundary_integrity/late_channel/dt_convergence.csv`).
P(collapse) is d*t*-converged in all 8 scheme×cell checks: every
d*t* = 0.025 value lies inside the paired 95% CI of the d*t* = 0.05
value (|*z*| ≤ 1.9). The *late-fraction share* at *J* = 20 is not:
it drops by ≈ 8–13 percentage points at d*t* = 0.025 in 3 of the 16
convergence tests (8 scheme×cell checks × two statistics) — clipped
EM at (2.4, 20): 0.257 → 0.180 (*z* = −2.3); clipped EM at (5.0, 20):
0.315 → 0.187 (*z* = −2.2); Lamperti at (5.0, 20): 0.525 → 0.390
(*z* = −2.4) — all in the same direction, and the Lamperti median
collapse time at (5.0, 20) drops from 160.5 to 15.8. The existence of
the late channel and its scheme ordering (larger under the
boundary-preserving scheme) are d*t*-robust, but its quantitative
share at the highest coupling is not yet d*t*-converged: every
late-fraction value at *J* = 20 quoted in this paper carries this
caveat.

**Collapse-definition sensitivity.** The collapse criterion
(*W* < 10 for 200 consecutive steps) contains two conventional
constants. We re-scored the *same* *n* = 2,000 headline-cell
trajectories (*J* = 5, *μ* = 100, clipped EM, matched seed) under a
3 × 3 grid of thresholds and durations (script:
`simulation/scripts/definition_sensitivity.py`; results:
`simulation/results/definition_sensitivity/`):

| *W*-threshold \ duration | 100 steps | 200 steps | 400 steps |
|---|---|---|---|
| *W* < 5 | 0.3100 | 0.3070 | 0.2985 |
| *W* < 10 | 0.3100 | **0.3070** (baseline) | 0.2985 |
| *W* < 20 | 0.3100 | 0.3075 | 0.2995 |

P(collapse) spans 0.2985–0.3100 across the nine definitions — a
spread of 1.15 percentage points, smaller than the *n* = 2,000
binomial 95% half-width of ±2.02 points — so the headline probability is
essentially inert to the definition constants. The failure-type gates
(rigidity if |*m*| > *r* at collapse, fragmentation if |*m*| < *f*,
else mixed) were perturbed on the same 614 collapsed trajectories
(baseline definition):

| gates (*r*, *f*) | rigidity | fragmentation | mixed |
|---|---|---|---|
| (0.8, 0.3) | 0.9528 | 0.0000 | 0.0472 |
| (0.9, 0.3) (baseline) | 0.8811 | 0.0000 | 0.1189 |
| (0.95, 0.3) | 0.7997 | 0.0000 | 0.2003 |
| (0.9, 0.2) | 0.8811 | 0.0000 | 0.1189 |
| (0.9, 0.4) | 0.8811 | 0.0000 | 0.1189 |

The rigidity share moves with the rigidity gate (0.80–0.95 across
*r* = 0.95–0.8), consistent with the integrator-scheme sensitivity of
the same statistic documented above; the fragmentation gate is inert
at this cell — the fragmentation share is 0.0000 at every tested gate,
so no collapsed trajectory sits below even the loosest fragmentation
threshold. The typing headline ("rigidity-dominated, fragmentation
absent at the headline cell") is robust; the specific rigidity-share
value is gate- and scheme-dependent, as flagged wherever it is quoted.

---

## §S5 Sensitivity and robustness analyses

The minimal model has six fixed parameters: *T* = 2 (temperature),
*ξ* = 0.5 (noise amplitude), the wealth-conversion coefficient
2/500 entering *h*(*W*), the consumption baseline 30, the
consumption-to-wealth slope 10/500, and the initial wealth
*W*(0) = 100. The sweep parameters are the coupling *J* and the
income multiplier *μ*. Three structural considerations and the
direct measurements in §S5.1–§S5.5 bear on the robustness of the
headline residual.

**First, the ablation cascade is itself a structural sensitivity
test.** L1 → L2 introduces the belief layer with three additional
parameters (*K*, $\mu_{\mathrm{b}}$, *γ*); L2 → L3 introduces the credit machine
with seven additional parameters. The L1–L3 residual point estimates
(23%, 30%, 34%) increase despite the parameter expansion — L1 → L3
+0.11 (unpaired 95% CI [+0.04, +0.18]; the levels use independent
seed bases, and adjacent steps are within sampling error) — so the
asymmetry is not a fine-tuned L1 artifact.

**Second, the Curie–Weiss critical point sits at *J* = *T*.** The
high-coupling regime in which Proposition 1 applies is *J* ≫ *T*; for
*T* = 2 this means *J* ≳ 4. The high-*J* sweep band lies in this
regime by construction. A two-fold change in *T* shifts the
bifurcation point but preserves the qualitative scaling.

**Third, the noise amplitude *ξ* = 0.5 controls the prefactor of
the Kramers-style escape time** but not the *h*/*J* leading-order
scaling of the deterministic landscape ratio. Smaller *ξ* would
sharpen the P(collapse) numbers; larger *ξ* would soften them. The
qualitative finding is insensitive to *ξ* within the
dynamical-stability range.

### S5.1 Direct sensitivity sweep — noise prescription and h(W) form

We ran a focused sensitivity sweep at the high-*J* band
(*J* ∈ {3.5, 4.0, 4.5, 5.0}) × protective-margin band
(*μ* ∈ {60, 80, 100}), 100 seeds per cell, 4,800 total runs. Three
non-baseline variants were tested against the baseline (multiplicative
noise $\xi \sqrt{1-m^2}$; linear *h*(*W*) = 2(*W*/500)):

- **Additive noise:** replace *g* = $\xi \sqrt{1-m^2}$ with *g* = *ξ*.
- **Logarithmic h:** *h*(*W*) = c·log(1 + *W*/100), with c calibrated
  so *h*(*W* = 100) matches the baseline (c = 0.4 / log 2 ≈ 0.577).
- **Constant h:** *h*(*W*) = 0.4 (= baseline value at initial wealth,
  state-independent).

Headline cell (*J* = 5, *μ* = 100), 100 seeds:

| variant | P(collapse) | rigidity share | SE |
|---|---|---|---|
| baseline (multiplicative noise, linear *h*) | 0.33 | 0.94 | 0.047 |
| additive noise, linear *h* | 0.20 | 0.30 | 0.040 |
| multiplicative noise, log *h* | 0.30 | 0.93 | 0.046 |
| multiplicative noise, constant *h* | 0.70 | 0.81 | 0.046 |

High-*J* band average (*J* ∈ {4.0, 4.5, 5.0}, *μ* = 100):

| variant | band-avg P(collapse) |
|---|---|
| baseline | 0.27 |
| additive noise | 0.12 |
| log *h* | 0.27 |
| constant *h* | 0.60 |

(All rigidity shares in this subsection are clipped-EM |*m*| = 0.9
gate values; §S4 shows this statistic is integrator-scheme-sensitive
— 0.88 clipped EM versus 0.74 Lamperti at the headline cell — while
P(collapse) and collapse timing are scheme-robust.)

Three findings. (i) The **log-*h* variant is statistically
indistinguishable from baseline** (0.30 vs. 0.33 at the headline cell;
band averages identical at 0.27; rigidity shares 0.93 vs. 0.94). The
qualitative residual collapse and rigidity-dominance findings are
robust to log-vs-linear *h*(*W*).

(ii) The **constant-*h* variant is substantially worse** (0.70 vs.
0.33 at headline; 0.60 vs. 0.27 band-average). Removing the
wealth-feedback loop produces a *more* pessimistic residual, not a
less pessimistic one. The baseline is therefore not a tuned
optimistic case — the wealth-driven recovery is helping the passive
stabilizer at sub-supercritical *J*, and the bare passive case
(no state-dependent feedback at all) fails substantially worse.

(iii) The **additive-noise variant has a lower P(collapse) but
qualitatively different failure-type composition**: rigidity share
drops from 0.94 to 0.30. This matches the prediction in §S2.6:
under additive noise the system can escape from *m* = ±1 saturation,
so trapped-on-wrong-branch (rigidity) collapses are converted into
disordered-stalling (fragmentation) collapses. The asymmetry
*itself* survives — P(collapse) is still 0.20, well above zero —
but the rigidity-typed signature reported in the main text is
prescription-specific. Main-text framing is consistent with this:
rigidity dominance is reported under the committed multiplicative
prescription, and the robustness of the residual to the noise
prescription is delegated to this section (Limitations).

Sweep script: `simulation/scripts/sensitivity_sweep.py`. Per-cell
CSVs: `simulation/results/sensitivity/`. Visual: Figure S2.

**Scaling-vs-fixed: a larger fixed field only postpones failure (Figure S3).**
The constant-*h* result above shows that *removing* the wealth feedback
worsens the residual; a complementary control probes large fixed
fields at high coupling. The operative variable here is the
**delivered field-amplitude profile over the run** (§S5.6): within any
single static-*J* run the "scaling" field is a constant rescaling of
the wealth-fed base 2(*W*/500), so the scaling-vs-fixed contrast below
is a delivered-magnitude-and-timing contrast between
parameterizations, not a within-run response contrast. We swept *J*
(extending to 20, beyond the calibrated *J* ∈ [0.5, 5], purely to
expose the mechanism) at *μ* = 100, 100 matched seeds, for the passive
field, two *fixed* but large fields (*h* = 2.4 and 5.0, the latter
12.5× the initial-wealth baseline field *h*(*W*₀) = 0.4), and the
coupling-proportional active field *h* = 2(*W*/500)·(1 + *J*):

| *J* | passive (*h*≈0.4) | fixed *h* = 2.4 | fixed *h* = 5.0 | active (*h* = 2(*W*/500)·(1 + *J*)) |
|---|---|---|---|---|
| 3 | 0.17 | 0.00 | 0.00 | 0.00 |
| 5 | 0.37 | 0.00 | 0.00 | 0.00 |
| 7 | 0.33 | 0.02 | 0.00 | 0.00 |
| 10 | 0.33 | 0.08 | 0.00 | 0.00 |
| 15 | 0.44 | 0.23 | 0.05 | 0.01 |
| 20 | 0.50 | 0.26 | 0.12 | 0.02 |

A larger *fixed* field postpones but does not escape failure:
in the 100-seed sweep of the table above the first detected failures
appear once each field's tilt-to-barrier ratio *η* = 4*h*/*J* has
fallen to ≈1.3–1.4 (*h* = 2.4 at *J* ≈ 7, *η* = 1.37; *h* = 5.0 at
*J* ≈ 15, *η* = 1.33); the higher-resolution *n* = 1,000 runs detect
the first nonzero failures earlier, at *η* ≈ 1.9–2.0 (*h* = 2.4 at
*J* = 5, P = 0.001; *h* = 5.0 at *J* = 10, P = 0.002), with failure
becoming appreciable as *η* falls toward ≈0.5–1.0 — the *η*-decay is the
landscape statement that motivated these measurements; whether it sets
the onset or is the mechanism of these failures is open (§S5.8) — while the
coupling-proportional field — which delivers a per-*J* rescaled
wealth-fed magnitude reaching 56–294 at equilibrium (§S5.6) — stays at
P(collapse) ≤ 0.02 throughout; the passive residual itself rises toward
the ½ selection ceiling of §S2.7 as *J* → 20 (P = 0.50 at *n* = 100,
binomial 95% CI ≈ [0.40, 0.60]; the closed form gives 0.425 at
*J* = 20, and item (b) of §S2.7's deepened characterization notes
the simulated range cannot pin the asymptote). The fixed-field
values at *J* > 5 in this table are clipped-EM estimates: the
boundary-preserving Lamperti check (§S4) shows the EM integrator
understates the fixed-field late-escape channel at *J* ∈ {15, 20} by
factors ≈ 1.4–1.6 (e.g. 0.183 EM versus 0.293 Lamperti at *h* = 2.4,
*J* = 15 in the *n* = 2,000 paired rerun), so these entries are
low-side readings. In the delivered-amplitude reading of §S5.6, what
this contrast measures is the delivered field-amplitude profile
required at each *J* — with Proposition 1's *η* = 4*h*/*J* as the
motivating quantity, small for the fixed parameterizations at high *J* and replenished by
the coupling-proportional parameterization — together with
the late-escape channel measured for the small-*η* fixed fields
(§S5.7), identified as noise-activated escape (§S9), its weight at
the highest couplings still open; it is not evidence of a
within-run response to *J*, which no
static sweep can exhibit. (The *J* = 7–20 points lie beyond the
calibrated range and are shown only to expose the mechanism.) Control
script: `simulation/scripts/scaling_control.py`; output:
`simulation/results/active_stabilizer/scaling_control.csv`.

**Response-coefficient (α) sweep.** The table below reports the
band-mean P(collapse) over *J* ∈ {4.0, 4.5, 5.0} at *μ* = 100
(100 seeds per cell) for the three responsive forms; the α ≥ 0.5
suppression threshold in Supplementary Section S2.6 uses the criterion band-mean
P ≤ 0.05.

| α | stress | coupling | quadratic |
|---|---|---|---|
| 0.1 | 0.233 | 0.183 | 0.233 |
| 0.25 | 0.157 | 0.110 | 0.223 |
| 0.5 | 0.013 | 0.040 | 0.197 |
| 1.0 | 0.000 | 0.000 | 0.153 |
| 2.0 | 0.000 | 0.000 | 0.080 |
| 5.0 | 0.000 | 0.000 | 0.000 |

The two linear forms first reach 0.000 at α = 1.0, while the
quadratic form *h* = 2(*W*/500)·(1 + α·(*J*/5)²) retains 0.153 at
α = 1.0 in this independent-seed sweep (0.14 in the matched full-grid
run) and reaches 0.000 only at α = 5. All three forms build on the
wealth-fed base 2(*W*/500) — there is no constant prefactor field —
and within a static-*J* run the coupling and quadratic forms are
constant rescalings of that base, so the α thresholds above are
delivered-magnitude thresholds; §S5.6 gives the form-versus-magnitude
accounting (at an informative operating point a constant matched to the
active form's full-trajectory delivered mean collapses far *less* than
the active form itself — control E2, verdict REFUTED in all three cells — so
the full-trajectory delivered mean is not the decisive quantity; §S5.6's
exploratory ordering test of other single summaries is domain-dependent
and yields no refutation under the domain rule the paper adopts). Data:
`simulation/results/active_stabilizer/alpha_sweep_results.csv`
(α-sweep) and
`simulation/results/active_stabilizer/active_quadratic_summary.csv`
(full-grid quadratic at α = 1.0).

### S5.2 One-at-a-time sensitivity over the six fixed parameters

We ran a one-at-a-time (OAT) sensitivity at the headline cell
(*J* = 5, *μ* = 100), 100 seeds per perturbation, varying each
parameter individually at ±50% from baseline:

| parameter | −50% P(collapse) | +50% P(collapse) | range | rigidity share |
|---|---|---|---|---|
| *T* (temperature) | 0.31 | 0.20 | 0.11 | 0.97 / 0.75 |
| *ξ* (noise amplitude) | 0.15 | 0.26 | 0.11 | 0.93 / 0.85 |
| *h*-conversion (2.0/500) | 0.48 | 0.26 | 0.22 | 0.94 / 0.85 |
| consumption baseline (30) | 0.28 | 0.30 | 0.02 | 0.93 / 0.83 |
| consumption slope (10/500) | 0.31 | 0.32 | 0.01 | 0.84 / 0.94 |
| initial wealth *W*(0) (100) | 0.31 | 0.18 | 0.13 | 0.90 / 0.89 |

Baseline at this cell with the OAT seed set: P(collapse) = 0.41
(SE = 0.05); the 0.41 here vs. 0.33 in §S5.1 reflects independent
seed-set variation — the 0.08 gap is ≈ 1.2 combined SE (per-arm SE
0.05 and 0.047, combined 0.068).

**The qualitative findings are robust across every perturbation.**
Across all 12 perturbed cells of the archived OAT run (hash-seeded;
regeneration is statistically consistent, not bit-identical — individual
cells shift within binomial SE, Reproducibility manifest), P(collapse) stays in the 0.15–0.48
range — no perturbation eliminates the high-coupling residual. The
rigidity-share-of-collapses stays in 0.75–0.97 — rigidity dominance
is preserved across all six parameters at ±50% (all shares in this
table are clipped-EM |*m*| = 0.9 gate values and carry the
integrator-scheme caveat of §S4: 0.88 clipped EM versus 0.74 Lamperti
at the same gate, with P(collapse) and timing robust).

The most consequential single parameter is the *h*-conversion
coefficient (the residual spans 0.22 across ±50%): halving it (to 1.0)
raises the residual to 0.48, and raising it 50% (to 3.0) lowers it to
0.26. This direct dose–response confirms
that *h* is the operative protective field; weakening *h* worsens
the asymmetry, strengthening *h* attenuates it but does not
eliminate it (P(collapse) = 0.26 at +50% *h*-coupling, well above
zero). Consumption-side parameters (baseline, slope) are nearly
inert (range ≤ 0.02).

A formal Sobol decomposition (~10 hours of compute at 100 sample
points × 9,000-run base sweep) would partition variance across
parameter interactions; we did not run it, but the OAT range above
already establishes that no single parameter accounts for the
high-*J* residual being above zero. Sweep script:
`simulation/scripts/oat_sensitivity.py`. Output CSV:
`simulation/results/sensitivity/oat_sensitivity.csv`.

### S5.3 Alternative model classes — is the asymmetry Curie–Weiss-specific?

Within the framework's cusp + reversible scope (Discussion), we additionally
test cross-class robustness across four bistable mean-field SDEs
that share the field-bias structure but differ in their
nonlinearity:

- **Curie–Weiss** (baseline): −*m* + tanh((*Jm* + *h*)/*T*).
- **Voter-like**: −*m* + sign(arg) · min(1, |arg|), where
  arg = (*Jm* + *h*)/*T*. Replaces tanh with a piecewise-linear
  sign function — the same saturation but with a discontinuity at
  zero, a "voter-style" bistable update.
- **Kuramoto-1D**: −sin(*m*π/2) + tanh((*Jm* + *h*)/*T*). Replaces
  the linear restoring force −*m* with a sine-driven analogue
  (Kuramoto-style on a 1D order parameter).
- **Cubic Landau**: −(*m* − *m*³/3) + tanh((*Jm* + *h*)/*T*). Replaces
  the −*m* term with the cubic Landau potential's gradient — a
  generic bistable system without the Curie–Weiss free-energy form.

Wealth feedback, sweep grid (high-*J* band *J* ∈ {3.5, 4.0, 4.5, 5.0}
× *μ* ∈ {60, 80, 100}), seeds (100 per cell), and collapse criteria
are identical to the minimal model. 4,800 runs total.

**Headline cell** (*J* = 5, *μ* = 100):

| variant | P(collapse) | rigidity share | SE |
|---|---|---|---|
| Curie–Weiss | 0.29 | 0.90 | 0.045 |
| voter-like | 0.30 | 0.93 | 0.046 |
| Kuramoto-1D | 0.24 | 0.92 | 0.043 |
| cubic Landau | 0.28 | 1.00 | 0.045 |

**High-*J* band-average** (*J* ∈ {4.0, 4.5, 5.0}, *μ* = 100):
P(collapse) is 0.28 / 0.32 / 0.17 / 0.30 for the four variants;
rigidity share is 0.84 / 0.95 / 0.70 / 1.00. (Rigidity shares in both
tables are clipped-EM |*m*| = 0.9 gate values and carry the
integrator-scheme caveat of §S4.)

The asymmetry is not Curie–Weiss-specific. All four variants show a
nonzero high-*J* residual at maximum *μ*, and rigidity dominates the
residual in three of four (the Kuramoto-1D variant is slightly
softer, plausibly because the sine restoring force allows the system
to escape ±1 saturation more easily under the same noise). The
operative requirement is bistable mean-field structure with field
bias, not the specific Curie–Weiss free-energy form. Sweep script:
`simulation/scripts/alternative_class_sweep.py`. Output CSVs in
`simulation/results/alternative_class/`. Visual: Figure S4.

### S5.4 Wealth-coupling coefficient required to hold P(collapse) fixed

The analyses above perturb the model at the baseline wealth-coupling
coefficient. A complementary measurement asks: what value of the
wealth-coupling coefficient *h*₀ of the protective field
*h*(*W*) = *h*₀·(*W*/500) is required to hold P(collapse) at a fixed
target level as the coupling *J* varies? We measured this on an *h*₀ grid with matched
seeds across the *h*₀ axis, so that comparisons along the grid are
paired: 2,000 paths per (*J*, *ξ*, *h*₀) cell at *μ* = 100, with the
simulation engine regression-checked to be bit-identical to
`minimal_model.run_cell` at the committed defaults. For each
(*J*, *ξ*) the profile P(collapse | *h*₀) — monotone non-increasing
in *h*₀ up to Monte-Carlo noise (paired matched-seed check) — was
inverted at the target levels P ∈ {0.10, 0.20, 0.30} to obtain
$h_0^*(J)$; uncertainty is from a 1,000-replicate path-level
bootstrap (95% percentile CIs). Targets unreachable on the
*h*₀ ∈ [0, 6] grid are recorded NA without extrapolation; this
affects *J* = 4.5 and 5.0 at (*ξ* = 0.65, target 0.10), leaving
*n* = 4 *J*-points for that cell.

Because the Curie–Weiss bifurcation at *J* = *T* = 2 sits inside the
manuscript's sweep, a single power law in *J* cannot hold across it;
the descriptor log $h_0^*$ = *a* + *p* log *J* is therefore fitted
only over the six supercritical couplings *J* ∈ {2.5, 3.0, …, 5.0}.
The nine (*ξ*, target) conditions give:

| *ξ* | target P | *n* (*J*-points) | *p* | 95% CI | *R*² |
|---|---|---|---|---|---|
| 0.35 | 0.10 | 6 | 1.90 | [1.74, 2.05] | 0.958 |
| 0.35 | 0.20 | 6 | 2.39 | [2.12, 2.59] | 0.909 |
| 0.35 | 0.30 | 6 | 2.95 | [2.67, 3.27] | 0.863 |
| 0.50 | 0.10 | 6 | 1.95 | [1.78, 2.13] | 0.919 |
| 0.50 | 0.20 | 6 | 2.64 | [2.36, 2.94] | 0.905 |
| 0.50 | 0.30 | 6 | 3.55 | [3.33, 3.75] | 0.881 |
| 0.65 | 0.10 | 4 | 2.32 | [1.89, 2.65] | 0.982 |
| 0.65 | 0.20 | 6 | 2.72 | [2.42, 3.04] | 0.910 |
| 0.65 | 0.30 | 6 | 3.62 | [3.33, 3.87] | 0.941 |

Two measurement statements follow. First, the wealth-coupling
coefficient *h*₀ required to hold P(collapse) fixed grows superlinearly in *J* over
the measured range: the fitted exponent *p* ranges from 1.90 to 3.62,
and *p* = 1 lies outside the 95% bootstrap CI in all nine
(*ξ*, target) conditions. Second, the exponent is not a single
number: it varies systematically with the target level and the noise
amplitude (spread 1.72 across the nine conditions), so *p* is a local
descriptor of the measured surface, not a substrate constant.

A three-parameter description $h_0^* = A\,(J - J_{\mathrm{c}})^{p'}$, with
$J_{\mathrm{c}}$ fitted by profile grid search (at each candidate $J_{\mathrm{c}}$, OLS of
log $h_0^*$ on log(*J* − $J_{\mathrm{c}}$); $J_{\mathrm{c}}$ chosen to minimize SSE; same
1,000-replicate path-level bootstrap), fits the same data with
*R*² ≥ 0.986 in every comparable cell:

| *ξ* | target P | *n* | *A* | *p*′ | 95% CI *p*′ | $J_{\mathrm{c}}$ | 95% CI $J_{\mathrm{c}}$ | *R*² |
|---|---|---|---|---|---|---|---|---|
| 0.35 | 0.10 | 6 | 1.90 | 0.63 | [0.52, 0.84] | 2.17 | [1.83, 2.31] | 0.9997 |
| 0.35 | 0.20 | 6 | 1.33 | 0.58 | [0.50, 0.73] | 2.37 | [2.21, 2.43] | 0.9988 |
| 0.35 | 0.30 | 6 | 0.81 | 0.60 | [0.53, 0.79] | 2.44 | [2.35, 2.45] | 0.9991 |
| 0.50 | 0.10 | 6 | 3.29 | 0.52 | [0.41, 0.64] | 2.33 | [2.19, 2.42] | 0.9955 |
| 0.50 | 0.20 | 6 | 1.86 | 0.65 | [0.52, 0.80] | 2.37 | [2.23, 2.45] | 0.9984 |
| 0.50 | 0.30 | 6 | 1.01 | 0.80 | [0.66, 0.96] | 2.41 | [2.34, 2.45] | 0.9984 |
| 0.65 | 0.20 | 6 | 2.09 | 0.70 | [0.57, 0.96] | 2.35 | [2.10, 2.44] | 0.9968 |
| 0.65 | 0.30 | 6 | 0.69 | 1.23 | [0.92, 1.56] | 2.16 | [1.92, 2.33] | 0.9870 |

The (*ξ* = 0.65, target 0.10) cell has only *n* = 4 *J*-points and is
flagged not comparable (fitted $J_{\mathrm{c}}$ = 1.79 with 95% CI [−1, 2.24]);
it is excluded from the counts below. Several upper CI ends in the
table truncate at the $J_{\mathrm{c}}$ ≤ 2.45 grid boundary of the profile
search (boundary-pinned bootstrap replicates: 261, 31, and 45 in the
affected cells), so the true upper CI ends lie higher — which
strengthens, not weakens, the departure of $J_{\mathrm{c}}$ from 2. Across the
eight comparable cells the fitted $J_{\mathrm{c}}$ ranges from 2.16 to 2.44, and
the 95% bootstrap CI on $J_{\mathrm{c}}$ excludes 2 in six of the eight; *p*′
ranges from 0.52 to 1.23. $J_{\mathrm{c}}$ is reported as a fitted descriptive quantity
of this three-parameter family: it is not identified with the
deterministic bifurcation point *J* = *T* = 2, and the fits do not
support substituting the bifurcation point for it.

**Alternative descriptor with the onset pinned at the theoretical
value.** Because the theory's onset sits at *J* = 2 exactly (§S5.8),
a natural two-parameter alternative to the log-*J* descriptor pins
the onset there: log $h_0^*$ = *a* + *p*·log(*J* − 2). On the same six
supercritical *J*-points this descriptor fits *better* than the
log-*J* power law in every (*ξ*, target) condition — *R*² =
0.950–0.997 versus 0.863–0.982 — with near-linear exponents
*p* = 0.76–1.45 (`exponent_fits.csv`, columns `p_alt_Jminus2` and
`R2_alt`). We do not adopt it as the headline descriptor because the
free-$J_{\mathrm{c}}$ profile fits above reject the pinned onset: the 95%
bootstrap CI on $J_{\mathrm{c}}$ excludes 2.00 in six of the eight comparable
conditions; the two whose CIs contain 2.00 are (*ξ* = 0.35,
target 0.10: $J_{\mathrm{c}}$ = 2.17, CI [1.83, 2.31]) and (*ξ* = 0.65,
target 0.30: $J_{\mathrm{c}}$ = 2.16, CI [1.92, 2.33]) (`jc_free_fits.csv`).
The better *R*² of the pinned-onset form shows the measured surface
is well described by a near-onset power law; the free-onset fits show
the onset itself sits above 2 — the same unexplained offset discussed
in §S5.8.

Locality caveat: both descriptions are local to the measured range
2.5 ≤ *J* ≤ 5.0 at *μ* = 100, *ξ* ∈ {0.35, 0.5, 0.65}; we make no
claim about the functional form outside this range. Scripts:
`simulation/scripts/h0_scaling.py`,
`simulation/scripts/jc_free_fit.py`. Results:
`simulation/results/h0_scaling/` (`SUMMARY.md`, `exponent_fits.csv`,
`B4_SUMMARY.md`, `jc_free_fits.csv`). §S5.8 overlays a parameter-free
theory-derived surface on these measurements — inverting the §S2.7
selection formula reproduces $h_0^*(J)$ at R² = 0.906 with median
pred/meas 1.12 over *J* ≥ 3.5 — and compares the theory's onset at
*J* = 2.00 with the fitted $J_{\mathrm{c}}$ above (an unexplained 7–18% offset).

### S5.5 Sensitivity of the wealth equation

The wealth dynamics of Equation 2,
Ẇ = *μ*(1 + *m*)/2 − [30 + 10(*W*/500)], with wealth floored
(clamped) at *W* = 0 at each integration step, is a modeling choice: the affine consumption term is
neither derived from microeconomic first principles nor empirically
calibrated. It encodes two qualitative assumptions — a fixed
subsistence cost (the constant 30) and a weakly wealth-dependent
component (the slope 10/500) — and the floor prevents unbounded debt.
To delimit the range over which the conclusions hold, we repeated the
analysis under ±30% perturbations of both consumption parameters
((21, 10), (39, 10), (30, 7), (30, 13)) and under a structurally
different, purely linear consumption law *C*(*W*) = 0.08·*W* (matched
to the baseline total consumption of 40 at *W* = 500, with no
subsistence constant). Setup: *h*₀ = 2.0, *ξ* = 0.5, *μ* = 100,
*J* ∈ {2.5, …, 5.0}, *n* = 2,000 matched-seed paths per cell; seeds
and 1,000-replicate bootstrap resamples are identical across
variants, so all comparisons are paired.

| variant | *J*=2.5 | *J*=3.0 | *J*=3.5 | *J*=4.0 | *J*=4.5 | *J*=5.0 |
|---|---|---|---|---|---|---|
| baseline (30 + 10·*W*/500) | 0.058 | 0.141 | 0.202 | 0.253 | 0.280 | 0.308 |
| constant −30% (21, 10) | 0.018 | 0.093 | 0.165 | 0.224 | 0.258 | 0.292 |
| constant +30% (39, 10) | 0.100 | 0.184 | 0.231 | 0.274 | 0.299 | 0.324 |
| coefficient −30% (30, 7) | 0.057 | 0.139 | 0.200 | 0.251 | 0.278 | 0.308 |
| coefficient +30% (30, 13) | 0.061 | 0.142 | 0.202 | 0.253 | 0.280 | 0.309 |
| purely linear (0.08·*W*) | 0.000 | 0.000 | 0.000 | 0.000 | 0.004 | 0.086 |

Under the four parametric perturbations, P(collapse) at the operating
point shifts by at most 0.047 in absolute probability across the
supercritical range, and the monotone increase of P(collapse) with
*J* is preserved. The subsistence constant dominates this parametric
sensitivity; the wealth-dependent coefficient is nearly inert
(max |ΔP| ≤ 0.0035). Re-measuring the §S5.4 exponent at the P = 0.30
level under each variant gives:

| variant | *p* | 95% CI |
|---|---|---|
| baseline | 3.55 | [3.33, 3.75] |
| constant −30% | 5.19 | [4.98, 5.39] |
| constant +30% | 2.55 | [2.28, 2.77] |
| coefficient −30% | 3.77 | [3.49, 4.07] |
| coefficient +30% | 3.43 | [3.14, 3.60] |
| purely linear | NA (target unreachable) | — |

The exponent stays within [2.55, 5.19] and remains far from 1 in
every parametric variant, so the superlinear growth of the required
wealth-coupling coefficient (§S5.4) is unchanged by the calibration of
Equation 2.

The purely linear consumption law behaves differently in kind:
collapse nearly vanishes (P = 0.000 for all *J* ≤ 4.0, 0.004 at
*J* = 4.5, 0.086 at *J* = 5.0), the P = 0.30 level is unreachable,
and *p* is undefined. The fixed subsistence term — not its precise
value, and not the wealth-dependent slope — is the structural
ingredient of Equation 2 on which the collapse phenomenology depends.
The conclusions therefore hold under ±30% recalibration of Equation 2
and are conditional on the presence of a fixed subsistence cost in
the consumption law. Script: `simulation/scripts/wdot_sensitivity.py`.
Results: `simulation/results/wdot_sensitivity/`. §S5.7 shows that under
the purely linear law the branch-selection pathology itself persists —
suppressed but not eliminated (P(wrong branch at *t* = 50) =
0.0365–0.2355 over *J* = 3.5–5.0, against affine-baseline values of the
same order but larger, e.g. 0.1340 at *J* = 3.5) — while its wealth
signature is muted.

### S5.6 Delivered-amplitude controls (magnitude-matched)

The stabilizer-variant comparisons in §S5.1 and §S2.6 contrast the
passive field, two fixed fields, and three "responsive" functional
forms. This section reads out what each form actually *delivers* as a
field trajectory in the static sweeps, and tests whether anything
beyond the delivered magnitude is at work (script:
`simulation/scripts/magnitude_matched_control.py`; results:
`simulation/results/magnitude_matched/`).

**Specification readout (as implemented).** There is no constant
"*h*ₐ" anywhere in the code: every response-form variant builds on the
wealth-fed base field base = 2·(*W*/500) (the fixed *h* = 2.4 and 5.0
controls are separate state-independent constants) (0.4 at *W*₀ = 100, rising to 14.0 at
the *μ* = 100 equilibrium *W*\* = 3500). The implemented forms are:
passive *h* = base; coupling *h* = base·(1 + *α·J*); quadratic
*h* = base·(1 + *α*·(*J*/5)²); stress
*h* = base + *α*·max(0, −*m*)·*J* (with *α* = 2.0 in the headline
stress runs). Because *J* is constant within every static-sweep run,
the coupling and quadratic fields are **per-*J* constant rescalings of
the wealth-fed passive field** within a run — there is no within-run
*J*-response in any variant; the stress form is the only variant whose field responds within a run
to the order parameter *m* (via max(0, −*m*)). Delivered field per
variant at the initial wealth *W*₀ = 100 and at the equilibrium
*W*\* = 3500 (*α* = 1 for coupling/quadratic):

| *J* | passive (*W*\* / *W*₀) | coupling (*W*\* / *W*₀) | quadratic (*W*\* / *W*₀) | fixed 2.4 | fixed 5.0 |
|---|---|---|---|---|---|
| 3.0 | 14.0 / 0.40 | 56.0 / 1.60 | 19.04 / 0.544 | 2.4 | 5.0 |
| 4.0 | 14.0 / 0.40 | 70.0 / 2.00 | 22.96 / 0.656 | 2.4 | 5.0 |
| 4.5 | 14.0 / 0.40 | 77.0 / 2.20 | 25.34 / 0.724 | 2.4 | 5.0 |
| 5.0 | 14.0 / 0.40 | 84.0 / 2.40 | 28.00 / 0.800 | 2.4 | 5.0 |
| 7.0 | 14.0 / 0.40 | 112.0 / 3.20 | 41.44 / 1.184 | 2.4 | 5.0 |
| 10.0 | 14.0 / 0.40 | 154.0 / 4.40 | 70.00 / 2.000 | 2.4 | 5.0 |
| 15.0 | 14.0 / 0.40 | 224.0 / 6.40 | 140.00 / 4.000 | 2.4 | 5.0 |
| 20.0 | 14.0 / 0.40 | 294.0 / 8.40 | 238.00 / 6.800 | 2.4 | 5.0 |

The coupling/quadratic-versus-fixed contrast at any given *J* is thus,
by construction, a comparison of delivered field-amplitude profiles
(56–294 at equilibrium for the coupling form versus a constant
2.4/5.0), differing additionally only in the wealth-dependence shape
(*h* ∝ *W* versus *h* constant).

**Realized field of the stress form (headline coefficient, *α* = 2.0).**
The stress term activates only when *m* < 0. Measured over *n* = 2,000
trajectories per *J* at *μ* = 100 and *α* = 2.0, the stress term is
active on < 1.1% of seed-steps within the decision window (*t* ≤ 50) — window activation fraction
0.0046–0.0102 across *J* = 4–20 — and on < 0.06% of seed-steps over
the full run: in > 99% of seed-steps the realized stress field equals
the passive field at the same (*W*, *m*), the rectified term being zero
whenever *m* ≥ 0. This is a per-seed-step field equality *within each
trajectory*, not an equality of arm outcomes — the rare active steps,
on the wrong-branch (*m* < 0) excursions, are exactly what makes the
stress arm's P(collapse) differ from the passive arm's. These
activation fractions are over the *t* ≤ 50 decision window; the field
delivered over the O(1) *t* ≤ 1 *selection* window, where the branch is
actually chosen (≈ 0.4–0.85, the passive base at *W*₀ = 100 plus the
rare rectified contribution; §S2.7, and the ordering test below), is a
different quantity and must not be conflated with the decision-window
mean below. The realized *decision-window* (*t* ≤ 50) mean field is ≈
5.13–5.24 (dominated by the passive ramp of base = 2·(*W*/500) from
0.4 toward 14, not by the stress term), and the full-run mean field is
≈ 13.15–13.16.

**Magnitude-matched constant-field control at the headline coefficient
(E1, *α* = 2.0).** We first replaced the stress field with a *pure
constant* field matched to either its realized window-mean (≈ 5.15) or
its realized full-run mean (≈ 13.15), with matched seeds (*n* = 2,000
per cell, paired bootstrap CIs; results `collapse_comparison.csv`):

| *J* | P(stress) | P(const = window mean) | P(const = full mean) | paired diff, stress − window [95% CI] | paired diff, stress − full [95% CI] |
|---|---|---|---|---|---|
| 4.0 | 0.000 | 0.000 | 0.000 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] |
| 4.5 | 0.000 | 0.000 | 0.000 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] |
| 5.0 | 0.000 | 0.000 | 0.000 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] |
| 7.0 | 0.000 | 0.000 | 0.000 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] |
| 10.0 | 0.000 | 0.001 | 0.000 | −0.001 [−0.0025, 0.000] | 0.000 [0.000, 0.000] |
| 15.0 | 0.000 | 0.021 | 0.000 | −0.021 [−0.0275, −0.015] | 0.000 [0.000, 0.000] |
| 20.0 | 0.000 | 0.057 | 0.000 | −0.057 [−0.0675, −0.047] | 0.000 [0.000, 0.000] |

At the headline stress coefficient *α* = 2.0, the stress form does not
collapse at any tested *J* (P = 0.00 throughout), and the constant
matched to its full-trajectory delivered mean is **likewise** at
P = 0.00 at every *J* from 4 to 20 (paired difference 0, degenerate CI
[0, 0], *n* = 2,000). Because both arms are identically all-zero, this
is a **consistency check on two all-zero arms, not a discriminating
test** of form versus magnitude: two arms that never collapse trivially
agree, whichever hypothesis is true. (The window-mean constant ≈ 5.15
does show excess collapse at *J* ≥ 15 — 0.021/0.057 at *J* = 15/20,
paired CIs excluding 0 — but that is a *late-time delivered-magnitude*
effect: freezing the field at the early-window average removes the
wealth ramp for the remaining ≈ 950 time units, so it under-delivers
5.15 versus 13.15 over the full run; it is not evidence about the
max(0, −*m*) response form.)

**Informative magnitude-matched control (E2, *α* = 0.05).** A
discriminating test requires an operating point at which the stress
form collapses with strictly interior probability. Following a
pre-specified analysis plan archived with the code
(`simulation/prereg/PREREG_magnitude_matched_informative_2026-08-21.md`),
we lowered the stress coefficient to base = 2.0, *α* = 0.05, *μ* = 100.
The pre-registered search (`search_operating_point()`) scans *α* in
ascending order and returns the first grid point whose headline-cell P
falls strictly inside the open band (0.10, 0.50); because P decreases
with *α* (0.281, 0.264, 0.242, 0.223, 0.174, 0.049 at
*α* = 0.05–0.50; full headline-cell scan deposited as
`operating_point_alpha_scan.csv`, produced by `full_alpha_scan()`),
this is the *smallest* *α* — the weakest in-band
stabilizer, with the highest in-band P (0.281). Larger in-band *α*
(0.10–0.30) were therefore not re-run in the matched control. (The
pre-registration's own NOTE resolves its "smallest-index alpha" rule to
this smallest-*α* reading.) We re-ran the matched constant control
there (script: `simulation/scripts/magnitude_matched_informative.py`;
results `informative_control.csv`; *n* = 2,000 matched seeds per cell,
95% paired bootstrap CIs):

| *J* | stress full-trajectory delivered mean | P(stress) | P(const = full-traj. mean) | paired diff, stress − const [95% CI] | verdict |
|---|---|---|---|---|---|
| 5.0 | 11.75 | 0.281 | 0.000 | +0.281 [+0.262, +0.301] | REFUTED |
| 4.5 | 12.35 | 0.261 | 0.000 | +0.261 [+0.242, +0.280] | REFUTED |
| 4.0 | 12.67 | 0.234 | 0.000 | +0.234 [+0.216, +0.253] | REFUTED |

The plan fixed its decision rule in advance: **REPRODUCED** if the 95%
CI of the paired difference contains 0 and |diff| ≤ 0.03; **REFUTED**
if the CI excludes 0 and |diff| > 0.05. In all three cells the CI excludes 0
and the difference is ≈ 0.23–0.28 ≫ 0.05: the verdict is **REFUTED in
all three cells**. Matching the full-trajectory delivered mean does *not*
recover the stress outcome — at *J* = 5 the two fields carry the *same*
full-trajectory delivered mean (11.75) yet give P = 0.281 (stress form)
versus P = 0.000 (matched constant). The §S2.7 formula, applied
independently to each arm, predicts both outcomes: the stress arm's
selection-window field ≈ 0.464 gives *λ* = 1.370 and P = 0.291 against
the measured 0.281, while the matched constant *h* = 11.75 gives
*λ* ≤ 0 (linearly stable origin), predicting exactly 0 against the
measured 0.000 — two independent routes to the same contrast.

**Mechanism (why the full-trajectory mean fails).** The full-trajectory
mean is a plain time-average over all 2,000 runs, dominated by the late,
high-wealth portion of the long-surviving trajectories, where the
wealth-fed base has ramped to base = 2·(*W*/500) → 14 as *W* → 3500. A
constant set to that mean therefore delivers ≈ 12 for the *entire* run —
including the early transient, where it strongly over-protects and
drives P to 0. The active field, by contrast, delivers only ≈ 0.4
*during the selection transient* (the base at *W*₀ = 100), which is
precisely when the branch is chosen (§S2.7). Full-trajectory-mean
matching back-fills protective amplitude into the one window where the
active form supplies almost none — which is why control E2 (above)
refutes the full-trajectory mean.

**Ordering test: does any single delivered-amplitude summary predict the
outcomes? (post-hoc; not part of the two-arm pre-registered plan).**
Having ruled out the full-trajectory mean, we asked whether *some* single
delivered-amplitude summary — the amplitude during selection, or a
conditional or rate-relative version of it — orders the static-sweep arms
the same way P(collapse) does, and we evaluated each ordering step
*statistically* rather than by point estimate. For each arm at
*J* ∈ {4.0, 4.5, 5.0}, *μ* = 100, *n* = 2,000 matched seeds
(`delivered_mean_ordering.py`, `delivered_mean_ordering_ci.py`; deposits
`delivered_mean_ordering.csv`, `delivered_mean_ordering_ci.csv`,
`ordering_verdicts.csv`) we attached a matched-seed paired bootstrap 95%
CI (10,000 resamples; the same convention as
`magnitude_matched_informative.paired_boot_diff_ci`) to every quantity,
and measured five candidates against P(collapse): (i) the full-trajectory
delivered mean; (ii) the selection-window mean over *t* ≤ 1; (iii) that
mean conditional on wrong-branch excursions (*m* < 0); (iv) (iii)
restricted to eventually-collapsing runs; and (v) the window mean
relative to the origin-instability rate
$\lambda = -1 + (J/T)\,\mathrm{sech}^2(h_{\mathrm{eff}}/T)$ with
$h_{\mathrm{eff}}$ = (ii). Sorting the arms by a summary, an adjacent pair
is a *violation* when the amplitude rises but P(collapse) also rises; it
**refutes** the summary only when both the amplitude gap and the
P(collapse) gap are resolved at 95%. (This ordering test draws an
independent seed stream from control E2 above — `delivered_mean_ordering.py`
and `magnitude_matched_informative.py` seed from different base constants
— so it re-measures the *α* = 0.05 stress arm at P(collapse) = 0.2865
against E2's 0.281, a difference within binomial sampling error at
*n* = 2,000.)

*Two domain definitions, and which we adopt.* The ordering test's
outcome depends on which arms count as inside the selection formula's
domain, and this SI has used two criteria that do not coincide. §S5.8
excludes an arm on **functional-form** grounds: the formula is an origin
linearization with a *constant* bias *c*, and the rectified stress field
max(0, −*m*)·*J* — identically zero for *m* ≥ 0, growing only on
*m* < 0 — cannot be represented by any constant *c*, so **both stress
arms are outside the domain** irrespective of their λ. A weaker, purely
operational reading instead keeps only arms with λ > 0 at the measured
selection-window field; under it the coupling *α* = 1 arm is excluded
(λ = −0.326 / −0.374 / −0.433 at *J* = 4.0 / 4.5 / 5.0;
`delivered_mean_ordering.csv`, `lam_at_hwin`) while the stress arms are
*re-admitted* (λ = 0.71 / 0.91 / 1.10). **We adopt the §S5.8
functional-form rule** as the paper's domain criterion and apply it
consistently in both sections; we also report the λ > 0 reading in full,
because the ordering-test verdict differs sharply between the two.

*Definedness.* Of the five summaries, only **three — (i), (ii),
(iii) — are defined on every arm.** Summary (iv) is undefined on stress
*α* = 2.0, which has no collapsing runs, so its collapse-conditional mean
does not exist; and it conditions on the very outcome it would predict,
so it is not a predictor in any case. Summary (v) is undefined wherever
*λ* ≤ 0 (`amp_over_lam` blank), so it cannot order an arm whose origin is
linearly stable. A summary undefined on an arm cannot order it in either
direction; we apply this criterion uniformly below.

*Result under the adopted (functional-form) rule.* With both stress arms
outside the domain, the arms that remain — passive, coupling *α* = 1,
quadratic, and the in-domain coupling *α* = 0.7 arm (below) — order
**monotonically** for every summary (i)–(v) at all three *J*: sorted by
any summary, P(collapse) is non-increasing, with no violation resolved or
unresolved (e.g. *J* = 5 full-trajectory mean 11.52 / 24.63 / 58.76 /
78.84 for passive / quadratic / coupling *α* = 0.7 / coupling *α* = 1
against P 0.300 / 0.157 / 0.024 / 0.005). **Under the domain rule the
paper adopts, the ordering test refutes no single delivered-amplitude
summary.** This monotonicity is close to structural rather than
empirical: with both stress arms removed, every surviving arm is the
same wealth-fed base 2(*W*/500) scaled by a per-*J* constant (×1
passive, ×1.6–2.0 quadratic, ×3.8–4.5 coupling *α* = 0.7, ×5–6 coupling
*α* = 1), so each delivered-amplitude summary and P(collapse) are both
monotone in that constant and no summary *could* be refuted. The
ordering test is informative only under the looser λ > 0 reading, where
the state-dependent stress arm re-enters.

*Result under the λ > 0 reading.* Re-admitting the two stress arms (and
excluding coupling *α* = 1) changes the verdict: the full-trajectory mean
(i), the selection-window mean (ii) and the rate-relative summary (v) are
each refuted with a resolved violation at all three *J*, while the
wrong-branch-conditional mean (iii) has **no resolved violation** (its
only such break, at *J* = 5 — stress *α* = 2.0 [cond_m_neg 1.858] versus
the in-domain coupling *α* = 0.7 arm [1.880] — has amplitude-gap CI
[−0.0249, +0.0679] spanning zero; `ordering_verdicts_primary.csv`).
Crucially, **every one of the ten point-estimate violations has the
stress *α* = 2.0 arm on the low side** (`ordering_verdicts_primary.csv`,
all `arm_lo` = `stress_a2.0`). That includes the in-domain coupling
*α* = 0.7 arm we built to probe the effect: scanning the coupling
coefficient (`coupling_domain_scan.csv`), the strongest *α* whose origin
stays linearly unstable is *α* = 0.7 (λ = +0.029 / +0.032 / +0.013; at
*α* = 0.8 the origin turns stable), delivering a selection-window mean of
1.72 / 1.88 / 2.05 — above the stress form's 0.80 / 0.83 / 0.85 — yet
collapsing at P = 0.0155 / 0.0190 / 0.0235, above the stress form's
0.000. But the resolved (ii) violation it produces at *J* = 4.0/4.5 is
the pair stress *α* = 2.0 → coupling *α* = 0.7 — **stress is still the low
side** — so the coupling arm does not establish a break independent of
the stress arm. Under the λ > 0 reading, the "no summary orders" outcome
is produced entirely by the single stress *α* = 2.0 arm that the adopted
functional-form rule excludes.

*Two further facts.* First — a *cross-J* observation, distinct from the
within-*J* ordering test above and not a restatement of the withdrawn
general negative — even the fields with a *J*-independent delivered
amplitude (the passive field, the fixed constants) show P(collapse)
varying across *J*: the passive field delivers 0.44 at *J* = 4.0 / 4.5 /
5.0 while P runs 0.257 → 0.288 → 0.300, and a fixed *h* = 2.4 delivers
2.4 at every *J* while P runs from 0.00 to 0.26; whether the
*η* = 4*h*/*J* decay sets that onset is open (§S5.8), and these same
cells carry the unexplained 5–140× under-prediction. Second, feeding the measured
selection-window amplitudes to the §S2.7 closed form does not reproduce
the stress form: 0.845 predicts P = 0.141, not the measured 0.000 (whereas
0.44 → 0.301 and 0.89 → 0.126 do track the passive and quadratic arms).

**We therefore report the ordering test as an exploratory measurement,
not as a result.** Under the domain rule the paper adopts (§S5.8,
functional-form), no single delivered-amplitude summary is refuted; under
the looser λ > 0 reading three are, but only through the single stress
*α* = 2.0 arm that the adopted rule excludes on functional-form grounds.
The robust finding is the pre-registered magnitude-matched control
(E2, above) — a direct two-arm outcome comparison that invokes no
formula and hence no domain rule: matching a stabilizer's
full-trajectory delivered mean does not reproduce its outcome (REFUTED in
all three cells, paired difference 0.234–0.281). What distinguishes the
parameterizations is left to the per-arm measurements, and the
pre-specified contrast (Figure 2) stands as a measurement.

**The α = 0.05 arm is a wealth-fed field with a weak, measurable state
response.** A post-hoc third arm (not part of the two-arm
pre-registered plan) at the informative operating point (b = 2,
α = 0.05, *μ* = 100, *n* = 2,000; `informative_passive_arm.py`,
`informative_passive_arm.csv`) compares the stress form against the
pure passive field (α = 0) on matched seeds. The rectified term is
active on 22.5–27.8% of decision-window (*t* ≤ 50) seed-steps, adding
a mean field of 0.17–0.23 on those active steps, and it lowers P(collapse) relative to
the passive field by 0.010–0.021 (95% CI excludes 0 at all three cells
*J* = 4.0, 4.5, 5.0) — a small protective effect consistent with the
rectified term suppressing wrong-branch lock-in. The α = 0.05 arm is
therefore not an inert stand-in for the passive field; the
delivered-amplitude contrast reported above (stress versus the matched
constant, same full-trajectory mean, P = 0.281 versus 0.000 at
*J* = 5) does not depend on it being inert.

**Reading.** Within the static sweeps at *μ* = 100, the robust
result is the pre-registered E2 control — a direct outcome comparison,
not a formula-domain question: the
full-trajectory delivered mean is not a sufficient statistic for the
outcome — a constant matched to it collapses far less than the active
form it was matched to (REFUTED in all three cells). Whether any single
delivered-amplitude summary orders the parameterizations is
domain-dependent (ordering test above): under the adopted functional-form
rule none is refuted, and under the λ > 0 reading the refutations are
carried entirely by the excluded stress *α* = 2.0 arm.
Statements elsewhere in this SI that contrast the coupling and quadratic
"responsive" forms with fixed stabilizers in static sweeps are therefore
statements about delivered field-amplitude *profiles* — the coupling and
quadratic forms deliver a uniform per-*J* rescaling of the base, the
stress form a rectified response on wrong-branch (*m* < 0) excursions
(§S5.8) — not about any single summary amplitude, and not about
within-run responsiveness of a rescaled field.

### S5.7 Branch-selection observable and collapse timing

The main text reports P(collapse), a wealth-ledger observable. The
underlying dynamical failure is selection of the wrong (−*m*) branch.
This section reports the branch observable directly, its dissociation
from the wealth observable, and the collapse-time decomposition of the
fixed-field failures (script:
`simulation/scripts/branch_observable.py`; results:
`simulation/results/branch_observable/`; before the runs, the engine
was regression-checked bit-identical to `minimal_model.run_cell` and
to the archived §S5.5 raw-linear results).

**Definition and window.** P(wrong branch) = fraction of runs with
*m* < 0 at *t* = 50 (step 1,000 of 20,000). Justification of the
window: passive-baseline collapse times saturate early (medians
14.7–17.0, maximum 41.9 across *J* = 2.5–5 at *μ* = 100;
collapse-time table below, `e3_collapse_times.csv`), and the deterministic branch-selection
transient is *O*(1) time units for supercritical *J*; *t* = 50 lies
safely after both, yet is early enough (5% of the horizon) to exclude
rare late stochastic escapes from the definition.

**Dissociation table** (*μ* = 100, affine consumption, *n* = 1,000 per
cell; binomial SE ≤ 0.016):

| condition | *J* | P(wrong branch, *t* = 50) | P(*m* < 0 at end) | P(W-collapse) | difference |
|---|---|---|---|---|---|
| passive | 2.5 | 0.015 | 0.000 | 0.061 | −0.046 |
| passive | 3 | 0.048 | 0.000 | 0.113 | −0.065 |
| passive | 3.5 | 0.141 | 0.000 | 0.215 | −0.074 |
| passive | 4 | 0.217 | 0.000 | 0.246 | −0.029 |
| passive | 4.5 | 0.223 | 0.009 | 0.256 | −0.033 |
| passive | 5 | 0.264 | 0.034 | 0.291 | −0.027 |
| fixed *h* = 2.4 | 3 | 0.000 | 0.000 | 0.000 | 0.000 |
| fixed *h* = 2.4 | 5 | 0.000 | 0.000 | 0.001 | −0.001 |
| fixed *h* = 2.4 | 7 | 0.013 | 0.001 | 0.026 | −0.013 |
| fixed *h* = 2.4 | 10 | 0.077 | 0.015 | 0.105 | −0.028 |
| fixed *h* = 2.4 | 15 | 0.122 | 0.101 | 0.178 | −0.056 |
| fixed *h* = 2.4 | 20 | 0.213 | 0.180 | 0.269 | −0.056 |
| fixed *h* = 5.0 | 3–7 | 0.000 | 0.000 | 0.000 | 0.000 |
| fixed *h* = 5.0 | 10 | 0.000 | 0.000 | 0.002 | −0.002 |
| fixed *h* = 5.0 | 15 | 0.026 | 0.005 | 0.035 | −0.009 |
| fixed *h* = 5.0 | 20 | 0.068 | 0.032 | 0.084 | −0.016 |
| stress (*α* = 2) | 4/4.5/5 | 0.000 | 0.000 | 0.000 | 0.000 |
| coupling (*α* = 1) | 4 / 4.5 / 5 | 0.003 | 0.000 | 0.003 / 0.003 / 0.005 | 0.000 to −0.002 |
| quadratic (*α* = 1) | 4 | 0.126 | 0.000 | 0.138 | −0.012 |
| quadratic (*α* = 1) | 4.5 | 0.109 | 0.004 | 0.130 | −0.021 |
| quadratic (*α* = 1) | 5 | 0.144 | 0.021 | 0.153 | −0.009 |

(The fixed-field rows at *J* > 5 are clipped-EM values; the
boundary-preserving check in §S4 shows the EM integrator understates
the fixed-field late channel at *J* ∈ {15, 20} by factors ≈ 1.4–1.6.)
In every affine-consumption cell the dissociation is small and
*negative*: among surviving runs P(*m* < 0 at *t* = 50) ≤ 0.0011 at
every *J*, so wrong-branch-at-*t* = 50 is essentially a subset of the
collapsed set, and the shortfall relative to P(W-collapse) is
post-collapse field destruction (after early wealth collapse the field
*h* ∝ *W* vanishes and a fraction of collapsed runs re-randomizes to
*m* > 0 by *t* = 50). In the affine model the two observables track
the same event. The variant rows also show that the stress and
coupling parameterizations act on branch selection itself
(P(wrong branch) = 0.000–0.003), not merely on wealth bookkeeping —
with the delivered-magnitude reading of §S5.6 applying to what drives
that suppression in static sweeps.

**Purely linear consumption: the observables dissociate.** Under
*C*(*W*) = 0.08·*W* (§S5.5; *h*₀ = 2.0, *n* = 2,000), the wealth
observable goes nearly blind — P(W-collapse) = 0.000 for all
*J* ≤ 4.0, 0.0045 at *J* = 4.5, 0.0865 at *J* = 5.0, bit-identical to
the archived §S5.5 values — while the branch observable remains at the
affine order of magnitude: **P(wrong branch at *t* = 50) = 0.0365,
0.1125, 0.179, 0.2355 at *J* = 3.5, 4, 4.5, 5** (affine baseline at
the same *h*₀ = 2.0, *n* = 2,000: 0.134, 0.2045, 0.2485, 0.2855).
Without the subsistence constant, income ≈ 0 on the wrong branch
decays *W* on a 12.5-time-unit timescale instead of driving it below
threshold; *W* never satisfies the *W* < 10 streak during the
transient, the field *h* = 2(*W*/500) stays large, and nearly every
wrong-branch run is eventually pushed back to *m* > 0 by the end of
the horizon (P(*m* < 0 at end) = 0.000 for *J* ≤ 4.5 and 0.007 at
*J* = 5.0, `e2_dissociation.csv`). The
branch-selection pathology is thus consumption-rule independent; the
main text's collapse-probability observable is a conservative proxy
whose *visibility* depends on the subsistence constant. In the affine
cells the two observables track each other (survivor
P(*m* < 0) ≤ 0.0011), so no affine-cell result is affected.

**Collapse-time decomposition of the fixed-field failures**
(clipped-EM, *n* = 1,000 per cell, physical time units; see the §S4
late-channel block for the boundary-preserving rerun of the *J* ≥ 15
cells):

| condition | *J* | P(collapse) | median *t* | p90 | p99 | max | fraction after *t* = 100 |
|---|---|---|---|---|---|---|---|
| fixed 2.4 | 5 | 0.001 | 13.8 | 13.8 | 13.8 | 13.8 | 0.000 |
| fixed 2.4 | 7 | 0.026 | 14.4 | 39.0 | 778.6 | 801.2 | 0.077 |
| fixed 2.4 | 10 | 0.105 | 14.5 | 503.8 | 816.0 | 975.9 | 0.190 |
| fixed 2.4 | 15 | 0.178 | 14.5 | 718.5 | 985.2 | 998.8 | 0.281 |
| fixed 2.4 | 20 | 0.269 | 14.3 | 562.9 | 918.4 | 986.9 | 0.197 |
| fixed 5.0 | 10 | 0.002 | 13.8 | 13.8 | 13.8 | 13.8 | 0.000 |
| fixed 5.0 | 15 | 0.035 | 14.2 | 481.9 | 837.5 | 858.2 | 0.229 |
| fixed 5.0 | 20 | 0.084 | 14.2 | 298.6 | 825.1 | 865.5 | 0.155 |
| passive | 2.5–5 | 0.061–0.291 | 14.7–17.0 | 17.1–22.8 | 22.2–36.9 | 25.6–41.9 | 0.000 |

The passive residual is purely early-transient: 0.000 of passive
collapses occur after *t* = 100 at any *J*. The fixed-field failures
at *J* = 7–20 are a **mixture** (excluding the fixed *h* = 5.0, *J* = 10
cell, which has only 2 collapses, both early): the bulk shares the
passive early transient (medians 14.2–14.5, versus passive 14.7–17.0),
but 8–28% of them occur after *t* = 100, with p99 collapse times
reaching ≈ 779–985
— a late, noise-driven escape channel with zero counterpart in the
passive baseline. These tail fractions are clipped-EM values at d*t* = 0.05; under
the boundary-preserving Lamperti integrator the late channel is larger
still (late fractions 0.45–0.53 and total P(collapse) higher by
factors ≈ 1.4–1.6 at *h* ∈ {2.4, 5.0}, *J* ∈ {15, 20}; §S4), and the
late-fraction share at *J* = 20 is not yet d*t*-converged — it
shrinks by ≈ 8–13 percentage points at d*t* = 0.025 in three of the
four *J* = 20 scheme×cell tests, while P(collapse) itself is
d*t*-converged in all checks (§S4). A fixed
field therefore thins, but does not change, the dominant
early-transient failure mode, while *adding* a late-escape channel;
the §S2.7 transient-selection account covers the early mode only
(§S5.8).

### S5.8 Selection-formula predictions and the theory-derived amplitude surface

This section stress-tests the §S2.7 closed form
P = Φ(−(*c*/*ξ*)$\sqrt{2/\lambda}$), with *c* = tanh($h_{\mathrm{eff}}$/*T*) and
λ = −1 + (*J*/*T*) sech²($h_{\mathrm{eff}}$/*T*), outside its original passive
validation domain, and inverts it into a parameter-free prediction of
the §S5.4 amplitude surface $h_0^*(J)$ (script:
`simulation/scripts/selection_formula_check.py`; results:
`simulation/results/selection_formula/` —
`e4_predicted_vs_measured.csv`, `e5_h0_theory_grid.csv`,
`e5_h0_pred_vs_measured.csv`, `e5_onset_locus.csv`). Where λ ≤ 0 the
origin is linearly stable and the formula predicts no instability
(P = 0 identically); such rows are flagged, not smoothed.

**Predicted versus measured P(collapse)**, every variant evaluated
under the single **selection-transient convention** Methods defines —
$h_{\mathrm{eff}}$ taken at the selection point *W* ≈ *W*₀ (the early
transient at which the branch is chosen, §S2.7), *not* a full-run or
decision-window average (*μ* = 100; measured values from the archived
static sweeps, *n* = 100 per cell unless noted; fixed-field *J* > 5
measured values are clipped-EM, understated ≈ 1.4–1.6× per the §S4
Lamperti check; source `convention_comparison.csv`):

| case | *J* | $h_{\mathrm{eff}}$ | λ | P(pred) | P(meas) | note |
|---|---|---|---|---|---|---|
| passive | 2.5 | 0.40 | 0.201 | 0.107 | 0.05 | over-predicts (λ small; §S2.7 flag i) |
| passive | 3.0 | 0.40 | 0.442 | 0.200 | 0.14 | |
| passive | 3.5 | 0.40 | 0.682 | 0.249 | 0.21 | |
| passive | 4.0 | 0.40 | 0.922 | 0.281 | 0.19 | |
| passive | 4.5 | 0.40 | 1.162 | 0.302 | 0.22 | |
| passive | 5.0 | 0.40 | 1.403 | 0.319 | 0.28 | |
| passive | 6–10 | 0.40 | 1.88–3.81 | 0.342–0.387 | 0.28–0.43 | extended sweep |
| fixed 2.4 | 3, 5 | 2.4 | −0.54, −0.24 | 0 (λ ≤ 0) | 0.00 | consistent |
| fixed 2.4 | 7 | 2.4 | +0.068 | 6 × 10⁻²⁰ | 0.02 | **measured failure unexplained** (λ marginal: 6.8% above the *J* = 6.56 onset, ≈10× below the §S2.7 calibration scale — hence reported separately from the *J* = 10–20 factors) |
| fixed 2.4 | 10 | 2.4 | 0.525 | 0.0006 | 0.08 | **under-predicts 140×** |
| fixed 2.4 | 15 | 2.4 | 1.288 | 0.019 | 0.23 | **under-predicts 12×** |
| fixed 2.4 | 20 | 2.4 | 2.050 | 0.050 | 0.26 | **under-predicts 5×** |
| fixed 5.0 | 3–10 | 5.0 | −0.96 to −0.87 | 0 (λ ≤ 0) | 0.00 | consistent |
| fixed 5.0 | 15 | 5.0 | −0.80 | 0 (λ ≤ 0) | 0.05 | **measured failure it fails to explain** |
| fixed 5.0 | 20 | 5.0 | −0.73 | 0 (λ ≤ 0) | 0.12 | **measured failure it fails to explain** |
| coupling (*α* = 1) | 4.0 | 2.0 | −0.16 | 0 (λ ≤ 0) | 0.01 | 1 collapse/100 seeds unexplained |
| coupling (*α* = 1) | 4.5, 5.0 | 2.2, 2.4 | −0.19, −0.24 | 0 (λ ≤ 0) | 0.00 | consistent |
| quadratic (*α* = 1) | 4.0 | 0.656 | 0.799 | 0.158 | 0.13 | $h_{\mathrm{eff}}$ = 0.4·(1 + (*J*/5)²) at *W*₀ |
| quadratic (*α* = 1) | 4.5 | 0.724 | 0.979 | 0.161 | 0.11 | |
| quadratic (*α* = 1) | 5.0 | 0.800 | 1.139 | 0.157 | 0.17 | |
| stress (*α* = 2.0) | 4–5 | 0.445–0.453 | 0.90–1.38 | 0.253–0.297 | 0.00 (*n* = 2,000) | **over-predicts**; stress excluded from domain (below) |
| stress (*α* = 2.0) | 7–20 | 0.451–0.474 | 2.33–8.46 | 0.340–0.410 | 0.00 (*n* = 2,000) | **over-predicts, rising to ≈ 0.41**; stress excluded from domain (below) |

Findings, stated as findings: (1) for fixed *h* = 2.4 the formula
under-predicts the *J* = 10, 15, 20 failures by 5×–140×, and
under-predicts the *J* = 7 failure far more severely
(P ≈ 6 × 10⁻²⁰ against a measured 0.02, flagged "measured failure
unexplained" above); (2) for fixed
*h* = 5.0 at *J* = 15/20, λ ≤ 0 — the formula predicts exactly zero —
against measured 0.05/0.12: with a large fixed field the origin is
stable, so these failures are not origin-selection events — they lie
outside the transient-selection mechanism's domain. Their timing is
mixed, not wholly late: median collapse times remain early-transient
(14.1–14.5 clipped EM; 18.7–160.5 Lamperti at d*t* = 0.05, the
latter not d*t*-converged) with late fractions 0.19–0.53 depending on
scheme (§S5.7, §S4), and the late component of these collapses is identified as
noise-activated escape (§S9), its weight at the highest couplings
remaining open; (3) fed the stress form's
selection-transient field, the formula **over-predicts** at every *J* —
0.253/0.281/0.297 at *J* = 4/4.5/5, rising to ≈ 0.41 by *J* = 20 —
against a measured 0.00 (*n* = 2,000) throughout. (The single value fed
here, ≈ 0.45, is a pinned-base window average, not the realized
selection-window delivered mean of ≈ 0.85 measured for the *α* = 2.0
field, §S5.6; feeding any single constant to a state-dependent field is
itself the issue.) This over-prediction is a **domain error, not a
formula failure**. The closed form is the
origin linearization of a *static* tilted double well, in which the
bias enters as a single constant *c* = tanh($h_{\mathrm{eff}}$/*T*); the
stress variant instead carries a *rectified, one-sided* response. Its
field *h* = base + *α*·max(0, −*m*)·*J* supplies a stress term that is
identically zero at the linearization point *m* = 0 and for all
*m* ≥ 0, activating only as the excursion goes negative; a constant
bias *c* — a fixed tilt of the potential independent of the sign of
*m* — cannot represent a bias that is zero on the *m* ≥ 0 side of the
origin and grows with the excursion on the *m* < 0 side. The stress
variant is therefore **excluded from the formula's domain** (its
measured 0.00 at *α* = 2.0 arises because the rectified
max(0, −*m*)·*J* term suppresses wrong-branch lock-in during the
selection transient, which is what gives the wealth-fed base time to
ramp far above the origin bias, §S5.6). (4) *The realized decision-window field is
not a valid rescue convention.* One can instead feed the formula the
realized decision-window mean field over *t* ≤ 50; for the stress
variant this is ≈ 5.1 (λ ≤ 0) so the prediction becomes 0, matching its
measured 0. But the *same* decision-window convention predicts **zero
collapse for the passive baseline** (window-mean field ≈ 3.2–4.6,
λ ≤ 0) against a measured 0.05–0.43, and **zero for the quadratic
variant** (window-mean ≈ 7.3–8.7, λ ≤ 0) against a measured 0.11–0.17.
A convention that zeroes out the very passive baseline the formula was
validated against is therefore not evidence of anything; the
selection-transient convention above — which does validate on the
passive baseline — is the one reported, and the stress over-prediction
is flagged as the domain error it is. The §S2.7 passive-domain
validation is unaffected (here R² = 0.694, MAE = 0.054 over *J* =
2.5–10 at *n* = 100/cell, versus R² = 0.93, MAE = 0.024 at 800 seeds in
§S2.7 — a seed-count difference, not a formula discrepancy).

**Theory-derived amplitude surface.** Inverting the same closed form
at a target probability P gives, per (*ξ*, target), the implicit
relation

$$ h_{\mathrm{eff}}^{*} = 2\,\mathrm{artanh}\!\left(|z_P|\,\xi\,\sqrt{\lambda(h_{\mathrm{eff}}^{*}, J)/2}\right), $$

(with *T* = 2 and *z*ₚ the standard-normal quantile), solved by
bracketed root-finding on a dense *J* grid. Mapping $h_{\mathrm{eff}}$ to the
§S5.4 coefficient *h*₀ requires an assumption about *when* selection
occurs; we state it: selection happens in the early transient at
*W* ≈ *W*₀ = 100, so $h_{\mathrm{eff}}$ = *h*₀·(*W*₀/500) = 0.2 *h*₀. The
alternative equilibrium-scale mapping ($h_{\mathrm{eff}}$ = 7*h*₀ at
*W*\* = 3500) fails by ≈ 30× (median pred/meas 0.034, negative R²);
the window mapping gives MAE 0.450, R² = 0.844, median pred/meas 1.18
over all 52 measured (*J*, *ξ*, target) cells — consistent with the
§S2.7 claim that selection is an early-transient event. Restricted to
the formula's stated domain *J* ≥ 3.5 (34 cells): **MAE = 0.323,
R² = 0.906, median pred/meas = 1.12**, all ratios within 0.86–1.71,
with a systematic over-prediction growing toward the target-0.30 level
(median ratio 1.19). At *J* < 3.5 the formula over-predicts $h_0^*$ by
up to ≈ 6× (worst at *J* = 2.5, target 0.30) — the slow-selection
regime §S2.7 flags. CI coverage is 13 of 52 cells overall (13 of 34
in-domain): the amplitude-requirement CIs are tight (*n* = 2,000), so the theory is a
≈ 10–20% approximation of the surface, not a within-CI predictor.

**Onset comparison.** The theory's onset locus is λ = 0 at
*J* = *T* cosh²($h_{\mathrm{eff}}$/*T*), so the predicted $h_0^*$ → 0 at
*J* → *T* = 2.00 exactly. The §S5.4 free fits give $J_{\mathrm{c}}$ = 2.16–2.44
across the eight unflagged cells (CI union 1.83–2.45): the theoretical
onset sits below every unflagged point estimate by 0.16–0.44 (7–18% of
the fitted value), though inside the union of the 95% CIs. Relatedly,
the theory's log-log slope of $h_0^*(J)$ is 0.89–0.96 over
*J* = 3.5–5, versus measured fitted exponents 1.16–1.93 on the same
restricted window: the measured surface steepens near onset faster
than the theory, which does not reproduce the near-onset behavior.
This offset is stated as an unexplained discrepancy; its sign is
consistent with slow near-threshold selection (finite-noise
trajectories escape the shallow well anyway, so a larger *J* is needed
before a nonzero $h_0^*$ is required), the same regime in which the
forward formula over-predicts.

---

## §S6 Network agent-based model: implementation, exact reduction, and validity range

The network ABM (Results, Figure 5) extends the scalar minimal model
to *N* = 200 agents on a graph. Each agent evolves the Curie–Weiss
SDE (main-text Equation 1) with the global mean-field term *Jm*
replaced by the local mean field *J* · m̄ᵢ, the row-normalized
neighbor average over agent *i*'s neighborhood (computed as the
row-normalized adjacency-matrix product). Wealth dynamics
(Equation 2) are unchanged and depend on the population-averaged
m̄ = (1/N) ∑ *m*ᵢ. Multiplicative noise $\xi \sqrt{1-m_i^2}$ η(t) uses a
single shared η(t) per step (rather than per-agent independent
noise), preserving the aggregate noise scale on m̄; the decomposition
of the two prescriptions, and the measured range over which the
scalar reduction holds, close this section.

Graph statistics for the six tested topologies:

| Topology | Nodes | Edges | Mean degree | Clustering |
|---|---|---|---|---|
| Complete | 200 | 19,900 | 199.0 | 1.00 |
| Erdős–Rényi (*p* = 0.10) | 200 | 1,979 | 19.8 | ~0.10 |
| Watts–Strogatz (*k* = 20, *β* = 0.1) | 200 | 2,000 | 20.0 | ~0.57 |
| Barabási–Albert (*m* = 10) | 200 | 1,900 | 19.0 | ~0.12 |
| Modular (4 × 50, $p_{\mathrm{in}}$ = 0.35, $p_{\mathrm{out}}$ = 0.01) | 200 | 1,867 | 18.7 | ~0.34 |
| Modular-boundary (same graph, *h* on boundary nodes only) | 200 | 1,867 | 18.7 | ~0.34 |

**Exact reduction to the one-dimensional SDE.** Under the committed
configuration — uniform initialization (*m*ᵢ(0) = 0 for all *i*), the
uniform stabilizer field *h*(*W*) applied to every agent,
row-stochastic (row-normalized) adjacency, and a single shared
per-step shock η(t) — the *N*-agent state remains exactly uniform for
all time: if **m** = *c*·**1**, the row-stochastic neighbor average
returns *c*·**1**, the drift is identical across agents, and the
noise increment $\xi \sqrt{1-c^2}$η is identical across agents. The
*N*-agent system therefore reduces *exactly* — not approximately — to
the one-dimensional SDE of the main text, for every row-stochastic
adjacency. The coupling term adj·**m** = **m** is exactly inert for
the five uniform-*h* topologies (complete, Erdős–Rényi,
Watts–Strogatz, Barabási–Albert, modular); the only variant in which
agents can differentiate is modular-boundary, where the stabilizer
field is applied to a subset of agents and the non-uniform field
breaks the symmetry.

The archived data exhibit this reduction to machine precision. In the
100-seed headline re-run (*J* = 5, *μ* = 100), the five uniform-*h*
topologies produce identical per-seed outcome rows — the same 100
trajectories — with P(collapse) = 0.33 (33 collapses = 30 rigidity +
0 fragmentation + 3 mixed; SE 0.047) for each of the five, and the
modular intra- and inter-community trace correlations under the
shared shock equal 1.0 to within 10⁻¹⁴. These are not five
experiments that happen to agree; they are one experiment — the
one-dimensional SDE — run under five topology labels.

| topology | P(collapse), 100 seeds | SE | rigidity share |
|---|---|---|---|
| complete | 0.33 | 0.047 | 0.91 |
| Erdős–Rényi | 0.33 | 0.047 | 0.91 |
| Watts–Strogatz | 0.33 | 0.047 | 0.91 |
| Barabási–Albert | 0.33 | 0.047 | 0.91 |
| modular | 0.33 | 0.047 | 0.91 |
| modular-boundary | 0.37 | 0.048 | 0.92 |

(Rigidity shares here are clipped-EM |*m*| = 0.9 gate values and
carry the integrator-scheme caveat of §S4.)

The five identical rows are identical because of the exact reduction,
not because five independent estimates converged. Pairwise topology
comparisons are therefore degenerate among the uniform-*h* variants
(they compare a sample with itself); the only non-degenerate
comparison is modular-boundary against the uniform-*h* value:
0.37 vs. 0.33, two-proportion *z* = 0.04/$\sqrt{0.047^2+0.048^2}$ ≈ 0.59
(two-sided *P* ≈ 0.55, *n* = 100 seeds per arm) — not separated at
this sample size. Re-run script:
`simulation/scripts/network_abm_100seed.py`. Output:
`simulation/results/network_abm/headline_100seed.csv`.

**Seed-set variation, not a noise-scale or topology effect.** In the
50-seed full sweep
(`simulation/results/network_abm/network_comparison.csv`), the seed
formula includes a per-topology content-hash offset, so dynamically
identical systems are driven by different noise streams; the
headline-cell scatter across the five uniform-*h* topologies
(0.24–0.34 at 50 seeds) is pure independent-seed sampling variation
(binomial SE ≈ 0.06–0.07). The 0.38 endpoint of the full 50-seed
spread belongs to modular-boundary — the symmetry-breaking variant,
which the exact reduction does not cover — not to a uniform-field
topology.
The same accounting applies to the comparison against the scalar main
sweep: the ABM 100-seed headline value (0.33) versus the scalar
main-sweep single-cell value (0.28 at 100 seeds) is a difference
between independent seed sets for identical dynamics — binomial
SE ≈ 0.046–0.047 at *p* ≈ 0.3, a gap of ≈ 0.75 combined SE. An
earlier version of this section attributed that gap to the ABM's
per-step shared η(t) noise scale relative to the scalar model; that
explanation was wrong: under the committed configuration the ABM
dynamics are identical to the scalar dynamics, so no noise-scale or
topology mechanism is available, and the gap is sampling variation.

**Modular-boundary variant.** The modular topology was constructed as
a planted-partition (stochastic block) model with 4 equal-sized
communities. Within-community edge probability $p_{\mathrm{in}}$ = 0.35 yields
~17 within-community neighbors per agent (measured mean 17.2 on the
archived seed-42 graph); between-community edge probability
$p_{\mathrm{out}}$ = 0.01 yields an expected 150 × 0.01 = 1.5 cross-community
neighbors (measured mean 1.47). The boundary-localized stabilizer
variant applies *h* only to the agents with at least one
cross-community edge — on the archived graph, 151 of 200 agents
(75.5%). Because *h* is
then non-uniform, exact uniformity breaks and its trajectories
genuinely differ from the one-dimensional reduction (per-seed trace
correlations 0.79–0.99 rather than 1.0). Its measured headline-cell
P(collapse) of 0.37 is not separated from the uniform-*h* value of
0.33 at these sample sizes (Z ≈ 0.59 above).

**Shared vs. independent noise.** The committed ABM drives all agents
with a single common-mode shock η(t) per step — the configuration in
which the exact reduction holds and the networked system reproduces
the scalar order-parameter dynamics by construction. Replacing the
shared shock with independent per-agent shocks ηᵢ(t) breaks the
uniformity. At the headline cell (*J* = 5, *μ* = 100, 50
deterministic seeds):

| topology | P(collapse), SHARED η(t) | P(collapse), INDEPENDENT ηᵢ(t) |
|---|---|---|
| complete | 0.40 | 0.00 |
| Erdős–Rényi | 0.30 | 0.00 |
| Watts–Strogatz | 0.36 | 0.00 |
| Barabási–Albert | 0.40 | 0.00 |
| modular | 0.32 | 0.00 |
| modular-boundary | 0.32 | 0.00 |

(The SHARED column re-samples the same one-dimensional dynamics under
yet another independent seed base; its 0.30–0.40 spread across
topology labels is again binomial sampling variation.) Under
independent per-agent noise, collapse does not arise on any topology
at *N* = 200, and the modular intra-/inter-community correlations
fall from exactly 1.0 (shared shock: identical trajectories) to
≈0.001. The mechanism is the standard mean-field fluctuation scaling:
with *N* independent shocks the noise on the order parameter averages
to ∝*ξ*/$\sqrt{N}$, and at *N* = 200 the effective amplitude is too small
to drive branch escape within the horizon. An *N*-scan (passive,
independent noise, mean-field coupling) is consistent:
P(collapse) = 0.24, 0.10, 0.02, 0.00 at *N* = 1, 5, 25, 200. This is
the network-level instance of the common-mode scope condition stated
in the main text and in the §S2.7 precondition: the collapse
phenomenology reported in this paper is a property of systems whose
order parameter carries aggregate *O*(1) noise. Under the shared
shock the networked system is the scalar system by construction, so
sparse or modular topology does not alter the common-mode residual;
no emergent synchronization from weak coupling is claimed. Control
script: `simulation/scripts/network_noise_control.py`; output:
`simulation/results/network_abm/noise_control.csv`.

**Validity range of the one-dimensional reduction.** The two noise
prescriptions above are the endpoints of a continuum. Writing the
per-agent shock as a mixture of a common component $\sigma_{\mathrm{com}}$·$\eta_{\mathrm{c}}$(t) and
an independent component $\sigma_{\mathrm{idio}}$·ηᵢ(t), the noise on the order
parameter is $\sigma_{\mathrm{eff}} = \sqrt{\sigma_{\mathrm{idio}}^2/N + \sigma_{\mathrm{com}}^2}$. This expression is an
identity — independent noise averages out in the mean — not a
dynamical result. To document the range of noise compositions over
which the one-dimensional reduction reproduces the ABM, we compared
ABM runs against the scalar model run at matched *ξ* = $\sigma_{\mathrm{eff}}$
(script: `simulation/scripts/reduction_validity.py`; data:
`simulation/results/reduction_validity/reduction_validity.csv`; cell
*J* = 5.0, *μ* = 100; total noise budget $\sigma_{\mathrm{idio}}$² + $\sigma_{\mathrm{com}}$² = 0.25;
bootstrap 95% CIs by path resampling).

Idiosyncratic endpoint ($\sigma_{\mathrm{com}}$ = 0), ABM at *N* agents vs. scalar
model at matched $\xi_{\mathrm{eff}}$ = 0.5/$\sqrt{N}$:

| *N* | ABM P(collapse) [95% CI] (n = 50) | scalar at $\xi_{\mathrm{eff}}$ [95% CI] (n = 1,000) |
|---|---|---|
| 1 | 0.24 [0.12, 0.36] | 0.287 [0.259, 0.316] |
| 5 | 0.10 [0.02, 0.20] | 0.122 [0.102, 0.143] |
| 25 | 0.02 [0.00, 0.06] | 0.005 [0.001, 0.010] |
| 200 | 0.00 [0.00, 0.00] | 0.000 |

Mixed noise at *N* = 200, common-mode share r = $\sigma_{\mathrm{com}}$²/0.25:

| r | ABM P(collapse) [95% CI] (n = 100) | scalar at matched $\xi_{\mathrm{eff}}$ [95% CI] (n = 1,000) |
|---|---|---|
| 0.00 | 0.00 | 0.000 |
| 0.25 | 0.17 [0.10, 0.24] | 0.152 [0.130, 0.175] |
| 0.50 | 0.18 [0.11, 0.26] | 0.212 [0.187, 0.238] |
| 0.75 | 0.32 [0.23, 0.41] | 0.257 [0.230, 0.284] |
| 1.00 | 0.23 [0.15, 0.31] | 0.314 [0.285, 0.344] |

The one-dimensional reduction at $\sigma_{\mathrm{eff}} = \sqrt{\sigma_{\mathrm{idio}}^2/N + \sigma_{\mathrm{com}}^2}$
reproduces the ABM collapse statistics across the full range of noise
compositions, within sampling error at the achieved precision (ABM
n = 50–100, binomial SE ≈ 0.02–0.06; the largest discrepancy is at
the pure common-mode endpoint *r* = 1.00 — ABM 0.23 versus scalar
0.314, a gap of ≈2 SE). This documents the validity range of the reduction
used throughout the paper; because $\sigma_{\mathrm{eff}}$ is an identity, the
agreement characterizes the reduction itself, not a network effect.

The committed sweep seeds use a stable SHA-1 content hash (Methods),
so the network results are bit-reproducible. Scripts:
`simulation/scripts/network_abm.py`,
`simulation/scripts/network_abm_100seed.py`,
`simulation/scripts/network_noise_control.py`,
`simulation/scripts/network_abm_figures.py`. Results:
`simulation/results/network_abm/`,
`simulation/results/reduction_validity/`.

---

## §S7 Supplementary figures

**Supplementary Figure S2.** *Sensitivity sweep — noise prescription and h(W) form
(supporting §S5.1).* Top row: P(collapse) heatmap over the high-*J*
band (*J* ∈ {3.5, 4, 4.5, 5}) × protective-margin band
(*μ* ∈ {60, 80, 100}) for four model variants — baseline
(multiplicative noise $\xi \sqrt{1-m^2}$, linear *h*(*W*) = 2(*W*/500)),
additive noise (constant *ξ*, linear *h*), logarithmic
*h*(*W*) = 0.577·log(1+*W*/100), and constant *h* = 0.4. Bottom row:
rigidity share at *μ* = 100 across the *J* range for each variant.
Headline cell (*J* = 5, *μ* = 100): baseline 0.33, additive 0.20,
log *h* 0.30, constant *h* 0.70. The qualitative residual is robust
to log-vs-linear *h*; constant *h* is substantially worse, ruling
out the baseline as a tuned best-case; under additive noise the
residual persists but rigidity-typed dominance falls (0.94 → 0.30;
shares under the clipped-EM |*m*| = 0.9 gate, integrator-scheme
caveat of §S4), confirming the rigidity-typed failure signature is
multiplicative-noise-specific while the asymmetry itself (and the
*h*/*J* landscape result, §S2.5) is not. 4,800 runs total
(100 seeds per cell × 12 cells × 4 variants). Source files:
`simulation/scripts/sensitivity_sweep.py`,
`simulation/results/sensitivity/figure_s4_sensitivity_sweep.png`
(rendered here as Figure S2, `manuscript/arxiv/figures/figS2.png`).

**Supplementary Figure S3.** *Scaling-vs-fixed control (supporting §S5.1).*
P(collapse) at *μ* = 100 vs. coupling *J* (100 matched seeds) for the
passive field (*h* ≈ 0.4), two *fixed* but large fields (*h* = 2.4 and
5.0, the latter 12.5× the initial-wealth baseline field
*h*(*W*₀) = 0.4), and the coupling-proportional
active field (*h* = 2(*W*/500)·(1 + *J*)). A larger *fixed* field only
*postpones* failure — each rises as its tilt-to-barrier ratio
*η* = 4*h*/*J* falls (first nonzero failures near *η* ≈ 1.9–2.0 at
*n* = 1,000; §S5.1) — while the
coupling-proportional field stays near zero. The shaded region (*J* > 5) is beyond the
calibrated range and shown only to expose the mechanism; the
fixed-field values there are clipped-EM readings, which the
boundary-preserving check of §S4 shows to understate the late-escape
channel at *J* ≥ 15 by factors ≈ 1.4–1.6. The contrast is a
delivered-magnitude/timing contrast between parameterizations
(§S5.6): the coupling-proportional field replenishes *η* by delivering
a per-*J* rescaled wealth-fed magnitude, while the fixed fields' *η*
decays — Proposition 1's landscape statement, which motivated the
measurement; whether the *η*-decay sets the onset of these failures is
open (§S5.8); the late component of the measured fixed-field failures is
identified as noise-activated escape (§S9); its weight at the highest
couplings remains open. Source:
`simulation/scripts/scaling_control.py`; figure
`manuscript/arxiv/figures/figS3.png`.

**Supplementary Figure S4.** *Alternative-model-class robustness (supporting §S5.3).*
Top row: P(collapse) heatmap over the high-*J* band (*J* ∈ {3.5, 4,
4.5, 5}) × protective-margin band (*μ* ∈ {60, 80, 100}) for four
bistable mean-field SDEs sharing the field-bias structure but
differing in nonlinearity — Curie–Weiss tanh, voter-like piecewise
sign, Kuramoto-1D sine, cubic Landau. Bottom row: rigidity share
at *μ* = 100 across the *J* range. At the headline cell
(*J* = 5, *μ* = 100): Curie–Weiss 0.29, voter-like 0.30,
Kuramoto-1D 0.24, cubic Landau 0.28 — a 0.24–0.30 band. The
high-*J* residual and rigidity-dominance pattern survive across
all four model classes (rigidity share at the headline cell is
0.90 / 0.93 / 0.92 / 1.00 respectively; clipped-EM |*m*| = 0.9 gate,
integrator-scheme caveat of §S4), demonstrating that the
asymmetry is generic to bistable mean-field SDEs with field bias
rather than specific to Curie–Weiss. 4,800 runs total (100 seeds ×
12 cells × 4 variants). Source:
`simulation/results/alternative_class/figure_s5_alternative_class.png`
(rendered here as Figure S4, `manuscript/arxiv/figures/figS4.png`).

---

## §S8 Supplementary Notes

Material moved here from the main text during condensation. Each note
carries the full statement whose compressed form appears in the main
text; all numbers are unchanged.

### Note 1 — Detailed comparison with adjacent recent measurements

Adjacent recent measurements differ in what is measured. A recent
preprint studies barrier crossing in a two-dimensional kinetic Ising
lattice under heat-bath Monte Carlo with the coupling held fixed
(*J* = 1), reporting magnetization distributions, first-passage
times, and time-averaged magnetization under bias and stochastic
fields [Oliver-Bonafoux, Toral & Chakrabarti 2026]; it does not
measure a collapse probability as a function of coupling, nor the
field amplitude required to hold one fixed. A Research Square
preprint treats a McKean–Vlasov mean-field double-well for normative
collapse, locating collapse and recovery thresholds and
hysteresis-loop properties under a swept control parameter
[Abreu Filho 2025]; the measurement in the main text is instead a
finite-horizon first-passage probability at fixed parameters under
matched-seed stabilizer contrasts, with neither a sweep protocol nor
a hysteresis loop. In an applied Ising setting, an enforcement field
*h* set against the coupling *J* exhibits the same tension studied
here [Zaklan, Westerhoff & Stauffer 2009], without measuring how the
required field amplitude grows with coupling. The bistable
coordination structure has economics antecedents beyond the
discrete-choice line: multiple search equilibria [Diamond 1982] and
strategic complementarity with multiple equilibria
[Cooper & John 1988]; the mean-field system measured here is an
instance of that coordination-failure structure. None of these works
reports a collapse probability as a function of coupling, or the
field amplitude required to hold one fixed, under matched-seed
stabilizer contrasts.

### Note 2 — Note on referents, Definition 1, and noise notation

**Note on referents.** Equations 1–2 of the main text constitute a
minimal two-variable dynamical system, and the results of the paper
are measurements on that system. The words "wealth", "consumption",
and "subsistence" are mnemonic labels for the slow feedback variable
*W* and its sink terms — chosen for readability, and consistent with
the social-interaction lineage of Equation 1 — not claims about any
empirical economy: no calibration of Equation 2 to data is claimed,
and no real-world unit is attached to *W* (the status and sensitivity
of the wealth equation are quantified in the main-text Methods and
§S5.5). The threshold is chosen as an unambiguously catastrophic
event — a sustained fall to 10% of initial wealth — not a transient
downturn; mapping it to specific real-world event categories is
deferred to future work. The minimal model contains no Hawkes shocks,
Lotka–Volterra terms, demographic dynamics, credit cycle, or
endogenous beliefs; those enter only at ablation Levels 2–4 (§S1).

**Definition 1 (formal statement).** A *passive stabilizer* is a
mechanism entering the system dynamics whose intervention magnitude
is bounded above by a constant independent of the system state and
the coupling level. Specifically: (i) its maximum restoring force
*h*ₘₐₓ is fixed — determined by structural parameters alone, not by
the current coupling *J*, the current stress level, or the time *t*;
(ii) it enters the dynamics as a restoring force favoring a
particular equilibrium. An *active stabilizer* is a mechanism whose
intervention magnitude is itself a function of the system state,
*h* = *g*(*m*, *J*, *t*), with *g* unbounded as stress increases. In
Equation 1, *h*(*W*) = 2·(*W*/500) is passive: the wealth-to-force
conversion rate is fixed at 2/500 and the field is bounded above by
the equilibrium wealth (itself bounded by *μ*). Definition 1 serves
as a labeling device for the parameterizations compared in the main
text; every conclusion is a measured comparison between fully
specified functional forms.

**Noise notation.** The white-noise process *η*(*t*), always written
with its time argument, is distinct from the dimensionless
tilt-to-barrier ratio *η* of Proposition 1; §S2.7 writes the same
white noise *ζ*(*t*).

### Note 3 — Delivered-amplitude controls: full values and caveats

**Binomial bounds for the matched-seed contrast.** The stress form's
P(collapse) = 0.00 across *J* ≥ 4.0 at *μ* = 100 carries a 95% upper
bound of 0.03 per 100-seed cell; where a zero count is reported at
*n* = 2,000 the 95% upper bound is 0.0015, and at *n* = 1,000 it is
0.003 (rule of three; Methods).

**Delivered equilibrium values (constancy identity).** At the
*m* = +1 wealth equilibrium (*W*\* = 3,500, base field 14.0) the
coupling form delivers 56–294 across *J* = 3–20 and the quadratic
form 19.0–238, against the fixed controls' 2.4 and 5.0. There is no
within-run *J*-response in the coupling or quadratic variant, so
their contrast with the fixed fields is definitionally a
delivered-magnitude comparison.

**Stress-form activation and window-mean control.** The stress form
activates its stress term on < 1.1% of decision-window seed-steps at
α = 2.0, *μ* = 100 (the strong term suppresses the very *m* < 0
excursions that would activate it); at the informative α = 0.05 the
weaker term is active on 22.5–27.8% of decision-window steps (§S5.6). A constant field matched to its *decision-window* mean
(≈5.13–5.24) under-delivers relative to the wealth-fed trajectory
over the rest of the run and shows excess collapse at high coupling —
0.001, 0.021, and 0.057 at *J* = 10, 15, 20 (clipped-EM values at
extended *J*, *n* = 2,000). §S4 shows the clipped-EM integrator
understates fixed-field failure at *J* ≥ 15, and this control was run
under clipped EM only — it has no boundary-preserving arm.

**α-sweep summary.** Band-mean residuals are 0.013 (stress) and
0.040 (coupling) at α = 0.5, first reaching 0.00 at the tested
α = 1.0, while the quadratic form retains 0.197 at α = 0.5 and
reaches 0.00 only at α = 5 (criterion disclosed in §S5.1: band-mean
P(collapse) ≤ 0.05 over *J* ∈ {4.0, 4.5, 5.0} at *μ* = 100).

**Extended-*J* per-cell scheme pairs.** For *h* = 2.4, a
boundary-preserving paired re-run at *J* = 15 and 20 gives 0.293 and
0.393 against matched clipped-EM arms of 0.183 and 0.247
(*n* = 2,000); for *h* = 5.0 the boundary-preserving paired values
are 0.039 and 0.1115 against EM arms 0.0285 and 0.0715. These
P(collapse) values are d*t*-converged under both schemes (§S4). The
coupling form's per-*J* rescaled field, delivering 224–294 at
equilibrium by *J* = 15–20, stays at P ≤ 0.02 throughout (clipped-EM,
100 seeds; §S5.1).

### Note 4 — Amplitude requirement: alternative descriptor and overlay detail

An alternative two-parameter descriptor that pins the onset at the
theoretical value, log $h_0^*$ = *a* + *p*·log(*J* − 2), fits the
§S5.4 data with higher *R*² (0.950–0.997) and near-linear exponents
0.76–1.45, but is not adopted because the free-$J_{\mathrm{c}}$ fits reject
$J_{\mathrm{c}}$ = 2 in six of eight comparable conditions (§S5.4).

Overlay detail: the theory's in-domain log-log slope is 0.89–0.96
against a measured 1.16–1.93; of the 34 in-domain ratios, 9 deviate
by more than 20%, all within 0.86–1.71. The onset offset has the
expected sign because near threshold selection is slow and
finite-noise trajectories escape the shallow well anyway, so a larger
*J* is needed before a finite $h_0^*$ is required; no quantitative
account is given (§S5.8).

### Note 5 — Late channel per-cell detail and ramp profile

**Late-channel per-cell counts and scheme ranges.** Paired
late-collapse counts rise from 75, 127, 11, and 45 (clipped-EM) to
302, 409, 35, and 117 under the boundary-preserving scheme at
(*h*, *J*) = (2.4, 15), (2.4, 20), (5.0, 15), (5.0, 20) — factors of
≈2.6–4 (paired schemes at d*t* = 0.05, *n* = 2,000 per cell). Total
P(collapse) in these cells spans the scheme ranges 0.18–0.29,
0.25–0.39, 0.029–0.039, and 0.072–0.112 respectively. Under the
boundary-preserving scheme at d*t* = 0.05 the λ ≤ 0 cell medians are
18.7 and 160.5 (the latter dropping to 15.8 at d*t* = 0.025).

**Ramp-duration profile at the critical margin.** At *μ* = 50 the
collapse probability rises from 0.19 at the shortest (1,000-step) ramp
through 0.26 and 0.40 to a ≈0.50 plateau (0.51, 0.51, 0.49 at 10,000,
20,000, and 30,000 steps), with the interpolated 50% crossing at
≈9,200 steps (`threshold_by_mult.csv`). Under sudden onset the
per-margin P(collapse) values are 0.47, 0.31, 0.11, 0.14, 0.06, 0.02,
0.02, 0.01, and 0.00 at *μ* = 40, 50, 60, 65, 70, 75, 80, 90, and 100
(`simulation/results/dynamic_j/summary.csv`), and the rigidity
share of the residual is 77–100% where *n* ≥ 10 (clipped-EM gate;
§S4).

**Ramp collapse under longer ramps.** Under the dynamic-*J* ramp
(`simulation/scripts/ramp_speed_sweep.py`;
`simulation/results/dynamic_j_ramp_speed/summary.csv` and
`raw_results.csv`), the ramp does not delay collapse near the viability
margin: at multiplier 40 the median collapse time is ≈ 21 time units
(median collapse step ≈ 414–445) regardless of ramp duration, so the
median coupling at the moment of collapse is the schedule read at that
early time — not a ramp effect — falling from 2.363 (multiplier 40) and
2.858 (multiplier 50) at the shortest (1,000-step) ramp to 0.562 and
0.637 at the longest (30,000-step) ramp, below the deterministic
bifurcation *J* = *T* = 2 and barely above the ramp's starting value
*J* = 0.5. (At multiplier 50 the median collapse step grows from ≈ 524
to ≈ 915 while the ramp lengthens 30-fold.) The longest-ramp cell at
multiplier 40 matches the static (*J* = 0.5, *μ* = 40) cell
(P = 1.00). At these longest ramps fragmentation is the largest typed
class of collapsed runs (0.54 and
0.51 of collapses at multipliers 40 and 50; the remainder mixed, with
rigidity ≤ 0.02; clipped-EM |*m*| < 0.3 gate, which carries the scheme
caveat of §S4). What triggers these early low-coupling collapses — as
opposed to wealth exhaustion at high coupling — is left open; the
earlier per-run fragment-before-wealth-floor diagnostic is withdrawn
because, initialized at *m* = 0, its ordering test fires at the first
step and cannot fail.

### Note 6 — Boundary, typing, and collapse-definition provenance

At the headline cell (*J* = 5, *μ* = 100), 0.353 of all
trajectory-steps hit the upper boundary and 0.022 the lower. The
lower boundary becomes an unattainable *entrance* boundary only for
large fixed fields near calibrated *J* (e.g. *h* = 5.0 at *J* = 5);
no boundary is absorbing in any studied regime (§S4). The paired
clipped-EM-vs-Lamperti difference in P(collapse) is +0.012 [+0.0065,
+0.0175]. The clipped scheme pins saturated trajectories at exactly
1 − 10⁻⁶ while the reflecting scheme scatters them over ≈[0.97, 1] —
the origin of the typing-split sensitivity at the |*m*| = 0.9 gate:
0.88 (0.881) in the *n* = 2,000 paired clipped-EM boundary run vs
0.74 (0.744) under the boundary-preserving integrator, with
P(collapse) and timing robust and fragmentation share ~0 under both
(§S4). Additive noise shifts the balance toward fragmentation:
rigidity share 0.94 → 0.30 at the headline cell in the 100-seed
noise-variant sweep (clipped-EM gate; §S5.1); the high-*J*
orthogonal-variations sweep comprises 4,800 runs (§S5.1). The
collapse-definition grid is *W*-threshold ∈ {5, 10, 20} × duration ∈
{100, 200, 400} steps, and the gate-pair perturbation range is
0.95–0.80 (§S4).

### Note 7 — Network extension detail

The six topologies are complete, Erdős–Rényi, Watts–Strogatz,
Barabási–Albert, planted-partition modular, and modular with
boundary-only stabilization, under uniform initial condition, uniform
field, row-stochastic adjacency, and a shared per-step shock. The
apparent spread among the five uniform-field topologies in the
50-seed sweep (0.24–0.34) is independent-seed sampling variation, not
a topology effect; the 0.38 endpoint of the full 50-seed spread
belongs to the symmetry-breaking modular-boundary variant, which the
exact reduction does not cover. ABM sample sizes are *n* = 50–100
(binomial SE ≈ 0.02–0.06); the largest reduction gap, at the pure
common-mode endpoint *r* = 1.00, is ABM 0.23 versus scalar 0.314
(≈2 SE). A platform-dependent seed hash was corrected to a stable
SHA-1 digest for bit-reproducibility (main-text Methods).

### Note 8 — Scope probes: full statement

Three substrates lie outside the cusp/fold + reversibility criterion
for structural reasons, not measurement difficulty: SIR epidemic
dynamics (transcritical at *R*₀ = 1, disease-free equilibrium
privileged; no bistable cells arose in a broad configuration scan)
[van den Driessche & Watmough 2002; Hadeler & van den Driessche
1997]; power-grid cascades (saddle-node, fold not pitchfork)
[Schäfer, Witthaut, Timme & Latora 2018; Manik, Timme & Witthaut
2017; Witthaut et al. 2022]; and phase-retrieval neural-network
training, where the pilot recorded fold / first-order behaviour with
hysteresis (mean hysteresis ≈ 0.75 against the 0.2 threshold) and
most tested seeds failed adiabatic reversibility [Sarao Mannelli
et al. 2019; Mondelli & Montanari 2019]. These do not enter the
positive evidence. The one pilot that failed adiabatic reversibility
(the neural network) is already excluded by the cusp/fold clause, so
the two clauses were not separately tested; the reversibility clause
rests on a heuristic.

To complete the map at the cusp + reversible corner the degenerate
optical parametric oscillator (DOPO) was tested — the canonical
ℤ₂-symmetric pitchfork and physical basis of coherent Ising machines
[Drummond, McNeil & Walls 1980, 1981; Wang et al. 2013; Marandi et
al. 2014]. Pilot measurements, pre-specified (protocol archived in
the code repository and the Zenodo archive), confirmed the near-cusp
field exponent, excluding the fold value, and matched the
parameter-free near-threshold coefficient; the mean cross-barrier
time rises monotonically with the pump above threshold, and the DOPO
passed the same adiabatic-reversibility test that most
neural-network seeds failed. The pre-specified, secondary,
confirmatory-only deep-well coupling-exponent test was inconclusive
at the achieved precision and cannot separate the corrected from the
original exponent. DOPO is a deliberately engineered cusp; its role
is corner-completion of the scope-criterion map — a high-confidence
confirmation of the cusp + reversible corner, with low could-fail
value — not a risky falsification test.

The mean-field idealization is not an approximation for the network
configuration tested: under a common-mode shock with a uniform
stabilizer field, the *N*-agent system on any row-stochastic topology
reduces exactly to the one-dimensional SDE (§S6), so the network
extension tests the noise-composition scope condition, not topology.
The one variant in which agents differentiate (boundary-only
stabilization) is not separated from the uniform-field value at the
tested sample sizes.

## §S9 Late-channel mechanism: noise-activated escape

This section reports the tests that identify the late collapse channel as noise-activated (Kramers) escape from the metastable well, a process distinct from the early transient branch selection. The tests, cells, and pass/fail thresholds were fixed in advance in a pre-registration committed before the results existed (`simulation/prereg/PRESPEC_late_channel_mechanism_2026-08-22.md`); the generating script is `simulation/scripts/late_channel_mechanism.py` and outputs are under `simulation/results/late_channel_mechanism/`.

All runs use the boundary-preserving (constant-diffusion Lamperti) integrator with constant field *h* and the identical wealth/collapse ledger as §S4, at d*t* = 0.05, 20,000 steps ($t_{\max}$ = 1,000), *T* = 2, ξ = 0.5, *μ* = 100, *W*₀ = 100, *n* = 2,000 per cell, over *h* ∈ {2.4, 5.0} × *J* ∈ {7, 10, 15, 20}. Early collapses are those with collapse time *t* ≤ 50; late collapses *t* > 250. The late-escape hazard *r* is the maximum-likelihood constant rate in the window (250, 1,000], i.e. late collapses divided by person-time at risk after *t* = 250 among trajectories uncollapsed at *t* = 250.

The timing separation of the two channels (main text, "The late channel") is: the wealth-fed baseline collapses are purely early (medians 14.7–17.0 across *J* = 2.5–5, all by *t* = 42, zero after *t* = 100), whereas the fixed-field failures are a mixture with passive-like medians 14.2–14.5 alongside the late tail (8–28% after *t* = 100, 99th percentiles ≈779–985; clipped EM, *n* = 1,000).

For reference, under this integrator the late channel strengthens relative to clipped Euler–Maruyama (paired late-collapse counts rise by ≈2.6–4× across the four failing (*h*, *J*) cells and the late fraction of collapses rises from 0.19–0.32 to 0.45–0.53 at d*t* = 0.05); P(collapse) is d*t*-converged in all eight scheme×cell checks, while the late-fraction share at *J* = 20 is not, shrinking ≈8–13 percentage points at d*t* = 0.025 (§S4). The late-channel *weight* at the highest couplings therefore remains open; the mechanism tests below concern the *rate*, which is well defined regardless of the weight.

### §S9.1 Test 1 — initial-condition dependence (decisive)

Each cell is run from two initial conditions on the same seed stream: (a) *m*₀ = 0 (disordered start; both channels can operate), and (b) *m*₀ = +0.9 (deep in the protected well; branch selection cannot operate). If the late channel is activated escape from the well, IC (b) should remove the early collapses while preserving the late-escape rate; if the late collapses are instead slow stragglers of the same selection transient, they should vanish along with the early ones.

| *h* | *J* | IC | P(collapse) | early (*t*≤50) | late (*t*>250) | late rate *r* (per unit *t*) |
|---|---|---|---|---|---|---|
| 2.4 | 7 | *m*₀=0 | 0.024 | 38 | 7 | 4.8×10⁻⁶ |
| 2.4 | 7 | *m*₀=+0.9 | 0.004 | 2 | 5 | 3.3×10⁻⁶ |
| 2.4 | 10 | *m*₀=0 | 0.137 | 141 | 111 | 8.3×10⁻⁵ |
| 2.4 | 10 | *m*₀=+0.9 | 0.078 | 2 | 126 | 8.8×10⁻⁵ |
| 2.4 | 15 | *m*₀=0 | 0.294 | 282 | 245 | 2.13×10⁻⁴ |
| 2.4 | 15 | *m*₀=+0.9 | 0.202 | 11 | 314 | 2.40×10⁻⁴ |
| 2.4 | 20 | *m*₀=0 | 0.413 | 366 | 372 | 3.65×10⁻⁴ |
| 2.4 | 20 | *m*₀=+0.9 | 0.279 | 15 | 437 | 3.52×10⁻⁴ |
| 5.0 | 7 | *m*₀=0 | 0.000 | 0 | 0 | 0 |
| 5.0 | 7 | *m*₀=+0.9 | 0.000 | 0 | 0 | 0 |
| 5.0 | 10 | *m*₀=0 | 0.0005 | 1 | 0 | 0 |
| 5.0 | 10 | *m*₀=+0.9 | 0.000 | 0 | 0 | 0 |
| 5.0 | 15 | *m*₀=0 | 0.040 | 45 | 28 | 1.93×10⁻⁵ |
| 5.0 | 15 | *m*₀=+0.9 | 0.017 | 1 | 28 | 1.89×10⁻⁵ |
| 5.0 | 20 | *m*₀=0 | 0.098 | 82 | 88 | 6.35×10⁻⁵ |
| 5.0 | 20 | *m*₀=+0.9 | 0.059 | 1 | 97 | 6.70×10⁻⁵ |

Pooled across cells, starting deep in the well removes 96.6% of early collapses (955 → 32) while the pooled late hazard is preserved (ratio *r*(b)/*r*(a) = 1.12). The per-cell hazard ratios, for cells with ≥ 10 late collapses in both arms, are 1.06, 1.13, 0.97, 0.98, and 1.05 — all within the pre-declared factor-of-two band and clustered at unity. The late-escape rate is thus independent of how the well was reached, the defining property of activated escape and the pre-declared decisive criterion (supported). The initial-condition-independent constant hazard is also the memoryless (Poisson) signature that Test 3 could not measure directly (below).

### §S9.2 Test 2 — barrier (Arrhenius) scaling and noise dependence

For each cell the deterministic drift *f*(*m*) = −*m* + tanh((*Jm* + *h*)/*T*) has a stable high-*m* well $m_+$ and an unstable saddle $m_s$. Two barriers are computed: the drift-potential barrier $\Delta U = U(m_s) - U(m_+)$ with $U(m) = \int [m - \tanh((Jm+h)/T)]\,\mathrm{d}m$, and the multiplicative-noise geometric barrier $G = \int_{m_s}^{m_+} f/(1-m^2)\,\mathrm{d}m$, for which the Itô stationary theory gives escape rate $r \propto \exp(-(2/\xi^2)\,G)$.

| *h* | *J* | $m_s$ | Δ*U* | *G* | late rate *r* | ln *r* |
|---|---|---|---|---|---|---|
| 2.4 | 10 | −0.302 | 0.638 | 0.833 | 8.3×10⁻⁵ | −9.39 |
| 2.4 | 15 | −0.185 | 0.582 | 0.776 | 2.13×10⁻⁴ | −8.45 |
| 2.4 | 20 | −0.133 | 0.559 | 0.752 | 3.65×10⁻⁴ | −7.92 |
| 5.0 | 15 | −0.388 | 0.805 | 1.009 | 1.93×10⁻⁵ | −10.85 |
| 5.0 | 20 | −0.279 | 0.715 | 0.912 | 6.35×10⁻⁵ | −9.66 |

Across the five cells with ≥ 10 late collapses, ln *r* is linear in Δ*U* (slope −10.98, *R*² = 0.953) and in *G* (slope −10.55, *R*² = 0.953). Kramers with multiplicative noise predicts a *G*-slope of −2/ξ² = −8.0; the fitted −10.55 is within the pre-declared ±40% band ([−11.2, −4.8]). The ≈30% overshoot is consistent with the finite-barrier and prefactor corrections not captured by the leading-exponential estimate. A noise sweep at (*h* = 2.4, *J* = 15) confirms the direction: the late hazard is 0 (0 late collapses), 2.13×10⁻⁴, and 9.3×10⁻⁴ at ξ = 0.35, 0.50, 0.65 — monotonically increasing as the effective barrier 2*G*/ξ² (12.67, 6.21, 3.67) falls.

**Monostable control.** At (*h* = 5.0, *J* = 7) the drift is monostable — the only fixed point is the high-*m* well, with no saddle and hence no barrier to cross — and the observed late-collapse count is exactly zero (P(collapse) = 0.000, both initial conditions). This is consistent with activated escape, though not uniquely diagnostic of it: the bistable (*h* = 5.0, *J* = 10) cell, whose barrier is large (2*G*/*ξ*² ≈ 10.2), also yields zero late collapses in the observation window, as an Arrhenius rate suppressed by that barrier predicts. The monostable zero therefore corroborates, rather than establishes, the escape picture; the initial-condition independence of Test 1 remains the decisive evidence.

### §S9.3 Test 3 — waiting-time distribution (underpowered by truncation)

The pre-registration also specified a memorylessness test: pooled late waiting times should be exponential (coefficient of variation CV ≈ 1, KS not rejecting an exponential at *p* > 0.05). This test failed its threshold — the late collapse-time and *m*-escape-time series have CV ≈ 0.56–0.68 and KS rejects the plain exponential (*p* < 0.05 in 10 of 11 scored series; the remaining series, the (*h* = 5.0, *J* = 15) collapse times at *n* = 28, has *p* = 0.112). The failure is a diagnostic-design limitation, not evidence against escape: the mean escape time 1/*r* is 2,700–52,600 time units while the observation window after *t* = 250 is only 750, so only the first 1–27% of the exponential is visible. Right-truncation of an exponential to a window short compared with its mean drives the CV toward the uniform-front value 1/$\sqrt{3}$ ≈ 0.577, which is what is observed. The test is therefore underpowered rather than contradictory; the constant, initial-condition-independent hazard established in Test 1 is the operative memorylessness signature. Measuring CV ≈ 1 directly would require running the horizon to a multiple of 1/*r* (tens of thousands of time units), which is deferred.

## References for the Supplementary Information

In addition to the references cited in the main text, this
Supplementary Information cites:

- Arecchi, F. T. & Politi, A. (1980). Transient fluctuations in the
  decay of an unstable state. *Physical Review Letters* **45**,
  1219–1222.
- Bacry, E., Mastromatteo, I. & Muzy, J.-F. (2015). Hawkes processes
  in finance. *Market Microstructure and Liquidity* **1**, 1550005.
- Dai, L., Korolev, K. S. & Gore, J. (2015). Relation between
  stability and resilience determines the performance of early
  warning signals under different environmental drivers. *Proceedings
  of the National Academy of Sciences USA* **112**, 10056–10061.
- Hänggi, P., Talkner, P. & Borkovec, M. (1990). Reaction-rate
  theory: fifty years after Kramers. *Reviews of Modern Physics*
  **62**, 251–341.

All four are in `references.bib` as `arecchipoliti1980`, `bacry2015`,
`daikorolevgore2015`, and `hanggi1990`. Bacry et al. is cited in
§S1.3 (the L4 Hawkes shock layer); Hänggi et al. is cited in §S2.6
(Kramers occupation theory); Arecchi & Politi (1980) and Dai et al.
(2015) are cited in §S2.7 (the experimental confirmation of selection
at an unstable point, and the stability–resilience distinction).
McKane & Tarlie (2001), also cited in §S2.7, is a main-text
reference and is not repeated here.
