# Why Fixed Protections Fail Under Rising Coordination: A Structural Failure-Mode Transition in Coupled Systems

## Abstract

Rising coordination is tightening the coupling between
once-independent agents. Governance relies on passive stabilizers
(diversification, constitutional checks, automatic income supports)
that buffer shocks at fixed strength. We identify a structural
prediction unique to coordination dynamics: as coupling rises, the
*type* of failure transitions from fragmentation to rigidity. This
is confirmed in 22 years of S&P 500 data — drawdown-period
correlation rises 0.50 → 0.88 across coupling deciles, rigidity share
0% → 100%, with empirical–model rank correlation ρ = +0.92. The
underlying scaling law: a fixed field *h* opposing coupling *J*
loses effectiveness as *h*/*J*² with an irreducible residual. The
pattern survives across alternative model classes, topologies, and
a second domain (cryptocurrency, 64% supercritical vs ~25%
equities). An active stabilizer eliminates the
residual. A heuristic AI coupling proxy (14 models × 805 prompts)
locates frontier AI in the regime where passive protections decay.
Overcoming the asymmetry requires active stabilizers.

## Introduction

A defining feature of the current decade is rising coordination across
sociotechnical systems. Algorithmic systems are pulling once-independent
decisions into common alignment: high-frequency traders co-execute on
shared signals, recommender feeds homogenize political opinion, large
language models are converging on a small set of evaluative standards,
and global supply chains route through shared infrastructure that links
demand shocks across continents. The empirical signature is a tightening
of cross-domain correlation. Pairwise correlations among S&P 500 sectors
exceeded 0.9 during both the March 2020 selloff and the 2010–12 Eurozone
crisis (see Results). Political identity, media consumption, and
consumer preference are coupling at decadal-record levels [Levitsky &
Ziblatt 2018; Svolik 2019]. Mass-extinction rates exceed background by
tens to hundreds of times [IPBES 2019; Barnosky et al. 2011] despite
high baseline biodiversity, because the dominant ecological stressor
— anthropogenic climate forcing — is correlated across all species
and habitats. In each domain, agents that used to fail independently
are starting to fail together.

The institutional response to this trend has so far extended mechanisms
designed for the previous era. Financial regulators tighten Value-at-Risk
models and circuit-breaker thresholds; democracies bolster constitutional
checks; conservation programs expand species inventories; nascent AI
governance frameworks — the EU AI Act's risk tiers, NIST's AI Risk
Management Framework, and voluntary commitments by the Frontier Model
Forum — add audit requirements and capability thresholds. These
mechanisms share a structural property that distinguishes them from a
different class of intervention. They buffer shocks at a *fixed*
strength baked into the system's architecture, rather than scaling
their response to the magnitude of any particular threat.
Diversification holds the same allocation whether market correlation is
0.3 or 0.9. A circuit breaker triggers at the same 7% threshold whether
one asset is dropping or a thousand are. Constitutional separation of
powers does not become structurally stronger during an authoritarian
push. We term mechanisms of this kind passive stabilizers — fixed
structural protections, automatic and state-independent — to
distinguish them from active stabilizers whose response magnitude
scales with measured stress (central-bank interventions, emergency
fiscal stimulus, military or diplomatic action).

Whether passive stabilization remains adequate under rising coupling
is, we argue, the central open question of AI-era governance. We
identify a structural prediction unique to coordination dynamics —
as coupling rises, the *type* of failure transitions from
fragmentation to rigidity — and confirm it in 22 years of S&P 500
data (Spearman ρ = +0.92 against the model's type-composition).
Underlying the transition is a quantitative scaling law: passive
stabilizer effectiveness scales as *h*/*J*² and becomes structurally
insufficient at high coupling.

### Domain-specific evidence and the missing substrate

The phenomenon of passive-stabilizer failure under high coupling has
been described repeatedly in domain-specific terms: diversification
benefits collapse when correlations rise [Ang & Chen 2002; Longin &
Solnik 2001; Patton 2006]; circuit breakers can accelerate crashes
they are designed to interrupt [Chen et al. 2024; Sifat & Mohamad
2019]; constitutional checks fail under synchronized opinion
[Levitsky & Ziblatt 2018; Svolik 2019]; mass-extinction rates exceed
background despite high biodiversity because the dominant stressor is
correlated across species [IPBES 2019; Barnosky et al. 2011]. These
literatures share a missing common substrate: no quantitative scaling
law connecting coupling to protection effectiveness, and no formal
distinction between mechanisms that scale their response with the
threat and those that do not.

### Closest existing frameworks and what they miss

The closest prior frameworks are Daníelsson's endogenous-risk theory
[2002; 2013] and the parallel liquidity-funding spiral
[Brunnermeier & Pedersen 2009]; Ashby's Law of Requisite Variety
[1956; Muhlert 2026]; Perrow's *Normal Accidents* [1984]; Svolik's
polarization–democracy tradeoff [2019]; and the asymmetric-dependence
literature in finance [Ang & Chen 2002; Longin & Solnik 2001]. Each captures a domain-specific version of the
phenomenon without the cross-domain scaling law, without the *h*/*J*²
quantification, and without the active/passive distinction. A
detailed comparison is in Supplementary Information §S7.

### Contribution

Our contribution is fourfold. *First*, we make the active/passive
stabilizer distinction mathematically precise (Definition 1, Results):
a *passive stabilizer* has intervention magnitude bounded above by a
structural constant independent of *J*, stress level, and time;
an *active stabilizer* has magnitude that scales with the threat.
Diversification, constitutional checks, biodiversity, and the
income-multiplier feedback satisfy the passive definition;
Federal Reserve rate cuts, emergency fiscal stimulus, and military
intervention exemplify the active class. *Second*, we derive a
quantitative scaling law: in any system where a fixed field *h*
opposes a coupling *J*, the field's protective effect scales as
*h*/*J*² and decays asymptotically. The derivation is a standard consequence of Curie–Weiss
free-energy structure, but to our knowledge has not been extracted
as a governance principle. *Third*, we confirm the scaling law in a
9,000-run sweep of a minimal two-equation simulation, and show via
an ablation cascade that adding belief and credit dynamics
monotonically *worsens* the asymmetry — the minimal model is the
most optimistic case, not a toy that overstates. A matched-seed
active-stabilizer experiment confirms the converse: replacing the
bounded passive field with a stress-responsive or coupling-responsive
analogue of finite strength eliminates the high-coupling residual
entirely (Figure 5), demonstrating that the classification has
predictive content. An agent-based extension on six network
topologies (complete, Erdős–Rényi, Watts–Strogatz, Barabási–Albert,
modular, modular-boundary) shows that the scaling is
topology-independent — all six topologies converge to a residual in
the same 0.33–0.37 band at the headline cell, with no pair
statistically separated (Figure 3). A complementary ramp-speed
experiment characterizes a non-monotone speed-of-onset effect that
materially refines the policy implications. *Fourth*, we test the model's mechanism empirically. The headline
finding is the **failure-mode transition** the framework uniquely
predicts: as coupling rises, the *type* of failure changes from
fragmentation to rigidity. In 22 years of S&P 500 sector data,
drawdown-period correlation rises 0.50 → 0.88 across coupling
deciles, rigidity share 0% → 100%, with empirical–model rank
correlation ρ = +0.92 (Figure 7). No generic correlation-risk model
makes this type-of-failure prediction. The frequency test confirms
the same coupling axis: tail-risk and VaR exceedance track the model
at Spearman ρ = 0.81 and 0.87, cross-mechanism ρ = 0.98 (Figure 6).
A heuristic AI proxy on 805 shared prompts [Li et al. 2023] locates
frontier AI above the model's critical threshold *J*ᶜ = 2 — the
substantive claim is regime-placement; the precise *J* depends on a
non-calibrated mapping (Figure 8, SI §S8).

The ingredients are individually known; the contribution lies in
their combination and governance interpretation.

## Results

### The minimal model and a scaling law for passive stabilizers

We study the simplest dynamical system containing the two
ingredients common to every example in the Introduction: a coordinated
state variable subject to tunable coupling, and a fixed structural
feature opposing synchronization. Let *m*(*t*) ∈ (–1, 1) be a
mean-field coordination order parameter and *W*(*t*) the system's
economic wealth. The minimal model is

$$
\dot m = -m + \tanh\!\left(\frac{J\,m + h(W)}{T}\right) + \xi\,\sqrt{1-m^2}\,\eta(t), \tag{1}
$$

$$
\dot W = \mu \cdot \mathrm{employment}(m) - \mathrm{consumption}(W), \tag{2}
$$

with employment(*m*) = (1+*m*)/2, consumption(*W*) = 30 + 10·(*W*/500),
and the passive stabilizer field *h*(*W*) = 2·(*W*/500). The
temperature *T* = 2 sets the noise scale, *ξ* = 0.5 the noise
amplitude, and *η*(*t*) is white noise (Itô convention; see Methods).
The coupling *J* ∈ [0.5, 5] tunes the pull toward synchronized
*m* ≈ ±1, and the income multiplier *μ* ∈ [20, 100] tunes the rate at
which the wealth-driven field *h* grows. Equation 1 is a Curie–Weiss
SDE in the order parameter; Equation 2 is a deterministic wealth
balance — no Hawkes processes, Lotka–Volterra terms, demographic
dynamics, credit cycle, or endogenous beliefs. Two coupled equations,
six parameters, one source of stochasticity.

We identify *h* as the *passive stabilizer*. We define the class
formally on the basis of the *intervention magnitude*, not the
functional form of the rule:

**Definition 1.** *A* passive stabilizer *is a mechanism entering the
system dynamics whose intervention magnitude is bounded above by a
constant independent of the system state and the coupling level.
Specifically:*
- *(i) its maximum restoring force h*ₘₐₓ *is fixed — determined by
  structural parameters alone, not by the current coupling J, the
  current stress level, or the time t;*
- *(ii) it enters the dynamics as a restoring force favoring a
  particular equilibrium.*

*An* active stabilizer *is a mechanism whose intervention magnitude is
itself a function of the system state, h = g(m, J, t), with g
unbounded as stress increases.*

The distinction turns on the intervention's *output* magnitude, not
on whether the rule has fixed coefficients. The Taylor Rule's
coefficients are fixed, but its rate-cut output scales with the
output gap; a circuit breaker halts trading for 15 minutes regardless
of move size, and diversification is bounded by 1/√*N* regardless of
correlation level. The criterion is boundedness of the restoring
force.

In Equation 1, *h*(*W*) = 2·(*W*/500) realizes a passive stabilizer:
the conversion rate from wealth to restoring force is fixed at 2/500,
and the field is bounded above by the equilibrium wealth (itself
bounded by *μ*). Wealth itself responds to the coupling
indirectly — high *J* locks *m*, which shifts employment, which shifts
wealth — but the conversion rule from *W* to *h* does not, and so the
maximum restoring force does not scale with the threat. Diversification,
constitutional checks, biodiversity, and the income-multiplier feedback
in macroeconomic models each satisfy Definition 1 with a different
state variable supplying *h* through a fixed conversion rule. The
minimal model contains no active stabilizer.

The following result is a standard consequence of the Curie–Weiss
free-energy structure; we state it explicitly because, to our
knowledge, it has not been extracted as a governance principle or
connected to the active/passive stabilizer distinction.

**Theorem 1** (*h*/*J*² *scaling of passive stabilizer effectiveness*).
*In the deterministic mean-field limit of Equations 1–2, the
preferential stability of the protective (+m) branch over the
unprotective (–m) branch scales as h/J² for J ≫ T. The passive
stabilizer's ability to steer the system toward the productive
equilibrium therefore vanishes at high coupling.*

Intuitively: at high *J*, the saturated fixed points sit at *m* ≈ ±1;
*h* tilts the landscape between branches by ≤ 2*h* (independent of
*J*) while the protective basin narrows as 1/*J*; the tilt-to-barrier
ratio decays as *h*/*J*². Because *h* is bounded by the equilibrium
wealth (itself bounded by *μ*), passive stabilization vanishes at
high coupling regardless of *μ*. Full derivation in Methods and SI.

Theorem 1 generates two qualitative predictions. First, at high *J*,
no attainable value of *h* fully eliminates collapse — the residual
collapse rate saturates above zero. Second, the residual collapses at
high *J* should be *rigidity-typed* (|*m*| near 1 at the moment of
failure, the system trapped on the unprotective branch) rather than
*fragmentation-typed* (|*m*| near 0, the system disordered and unable
to coordinate productive activity).

We tested both predictions with a 9,000-run sweep over *J* ∈ {0.5, 1.0,
…, 5.0} and *μ* ∈ {20, 30, …, 100}, with 100 independent seeds per
cell and 20,000 Euler–Maruyama steps per run at d*t* = 0.05 (full
numerical protocol in Methods). A run is classified as *collapsed* if
its wealth falls below 10 (10% of the initial value) for 200
consecutive steps (10 physical time units). The collapse criterion is
chosen to be unambiguously catastrophic within the model's dynamics —
representing sustained economic failure below any plausible recovery
margin, not a transient downturn. The mapping from this model-internal
threshold to specific real-world event categories (recession,
financial crisis, civilizational breakdown) is a calibration question
deferred to future work. A collapse is *rigidity-typed* if |*m*| > 0.9
at the collapse instant and *fragmentation-typed* if |*m*| < 0.3. The choice of d*t* = 0.05
follows from a convergence check across d*t* ∈ {0.1, 0.05, 0.025}
(see Methods) showing that P(collapse) saturates at d*t* ≤ 0.05.

Both predictions are confirmed (Figure 1). At the lowest coupling
(*J* = 0.5), the passive stabilizer eliminates collapse entirely as
*μ* rises from 20 to 100: P(collapse) drops from 1.00 to 0.00. In the
highest coupling band (*J* ≥ 4.0), the same stabilizer range only
reduces P(collapse) from 1.00 to 0.23 — a residual that does not
vanish. Of the residual collapses at maximum stabilizer, 90% are
rigidity-typed under multiplicative noise (SI §S5.1: additive noise
shifts the balance toward fragmentation).
The two failure modes are spatially separated in the (*J*, *μ*) plane:
fragmentation occupies the low-*J*, low-*μ* corner, while rigidity
dominates the high-*J* band at every *μ*.

The minimal model exhibits exactly the structural pattern Theorem 1
predicts: passive stabilization protects against fragmentation but
not against rigidity, where synchronized lock-in cannot be overridden
by bounded passive intervention. The 23% high-*J* residual at maximum *μ* is the
finite-noise, finite-*h* signature of the *h*/*J*² → 0 scaling.

### Robustness: an ablation cascade

The minimal model omits dynamical elements commonly invoked in the
socio-economic instability literature — endogenous belief dynamics,
credit cycles, demographic feedback. We tested whether the asymmetry
result is a structural artifact of those omissions with an ablation
cascade: identical (*J*, *μ*) sweep, seeds, and integration protocol,
adding complexity one level at a time.

**Level 0** is the bare Curie–Weiss SDE with fixed exogenous *h*; its
long-time mean matches the mean-field fixed-point analysis underlying
Theorem 1 (Supplementary Figure S1). **Level 1** is the minimal model
(residual 23% at *μ* = 100, 90% rigidity). **Level 2** adds
replicator–mutator dynamics over *K* = 6 belief types with
epistemic-coherence factor *h*ₑff = *h*·EC (Methods): residual
**30%**, rigidity 83%. **Level 3** adds Minsky–Keen credit-cycle
dynamics: residual **34%**, rigidity 81%. **Level 4** is the full simulator (belief + credit + Hawkes shocks +
Lotka–Volterra + demographics + mortality) reported from prior work
at d*t* = 0.1: 30–35% band, rigidity dominates; not re-run at
d*t* = 0.05 (SI §S1.3).

Across L1 → L2 → L3 at consistent d*t* = 0.05 the residual increases
monotonically (Figure 2; 23% → 30% → 34%) with rigidity dominance
preserved (90%, 83%, 81%). L4 (full simulator from prior work, run
at d*t* = 0.1) reports residuals in the 30–35% band, consistent with
the L1–L3 trend; we do not re-establish strict monotonicity through
L4 here, but its residual is within the L1–L3 envelope and adds no
counterevidence. The minimal model is the most optimistic case
across the tested cascade.

### Network robustness: the scaling law beyond mean-field

A natural concern about mean-field is that modular network structure
might restore passive-stabilizer effectiveness by decoupling
communities. We tested this with an agent-based extension: N = 200
agents on six graph topologies — complete (mean-field control),
Erdős–Rényi, Watts–Strogatz, Barabási–Albert, planted-partition
modular (4 communities), and modular with stabilizer applied only at
community boundaries — same Curie–Weiss local dynamics, wealth
feedback, sweep grid, and collapse criteria as the scalar model
(Methods).

The *h*/*J*² scaling survives on every tested topology (Figure 3).
At the headline cell (*J* = 5, *μ* = 100), a 100-seed re-run
(SI §S9) finds all six topologies converge to P(collapse) ∈ [0.33,
0.37], with no pair separated above one standard error (largest gap:
boundary-localized 0.37 vs. complete 0.33, Z ≈ 0.59 in a two-sample
test of proportions). **None** of the six tested topologies — complete,
Erdős–Rényi, Watts–Strogatz, Barabási–Albert, modular, or
boundary-localized — rescues the mean-field residual. Intra-community
and inter-community trajectory correlations both reach 1.00 at high
*J*: the entire network synchronizes despite only ~1% inter-community
edge density. The mean-field prediction holds across every tested
topology.

### Dynamic coupling: speed-of-onset asymmetry

Real-world coupling does not jump between fixed values. The minimal
model accommodates time-varying coupling by replacing the static *J*
with a schedule *J*(*t*). We tested three canonical schedules: a
*gradual* linear ramp from *J* = 0.5 to
*J* = 5.0 over the simulation horizon; an *exponential* S-curve
representing late-saturating adoption; and a *sudden* step function
modeling a discrete capability gain.

Across 100 seeds × 9 *μ* levels × 3 schedules × 30,000 steps (d*t* =
0.05), the gradual and exponential schedules survive at every margin
*μ* ≥ 60 (P(collapse) ≤ 3%), whereas sudden onset triggers rigidity
collapse across the band (P(collapse) = 47% at *μ* = 40, declining to
2% at *μ* = 80; rigidity dominates the residual at every margin
tested — ≥ 77%, 77–100% in cells with *n* ≥ 10 collapses).

The mechanism is *adiabatic settling*: gradual ramps build *h*
before *J* binds, entering Theorem 1's regime on the protective
branch; sudden onset lets noise select the branch.

A ramp-duration sweep over six durations from 1,000 to 30,000 steps
(5,400 runs) reveals three regimes (Figure 4): *resilient* (*μ* ≥ 60,
P(collapse) ≤ 7% at every ramp duration); *critical* (*μ* = 50,
P(collapse) monotonically increasing from 0.19 to 0.49; 50% threshold
near 9,200 steps); *unviable* (*μ* ≤ 40, P(collapse) saturates at
1.00 by ramp ≈ 4,000 steps). The collapse-type composition swaps:
fast ramps produce a partly rigidity-typed residual, while slow ramps
at low margin produce an overwhelmingly fragmentation-typed residual.

Slow deployment protects high-margin economies through adiabatic
settling, but near the viability margin gradualism *prolongs*
operation at high *J* and low wealth, converting partial rigidity
collapse into near-certain fragmentation collapse. Below the viability margin, no deployment speed avoids catastrophic
failure (slow → fragmentation, sudden → rigidity); only
active-stabilizer architecture navigates the transition safely.

### Empirical setting: diversification and coupling in financial markets

Equal-weight diversification is a paradigmatic passive stabilizer
satisfying Definition 1, and *J* is directly observable as the
average pairwise correlation ρ̄. The diversification benefit
satisfies the variance-algebra identity

$$
\mathrm{DB}(\bar\rho, N) = \frac{\bar\sigma_{\mathrm{single}}}{\sigma_{\mathrm{portfolio}}} = \frac{1}{\sqrt{1/N + (1 - 1/N)\,\bar\rho}}. \tag{3}
$$

Across 5,414 rolling 60-day windows of daily 9- to 10-sector-ETF
returns (30 March 2004 – 30 December 2025; Methods), empirical pairs
sit within MAE 0.016 of Equation 3; Spearman ρ(ρ̄, DB) = **–0.993**
(*p* < 10⁻¹⁵). The tight fit confirms that the equal-weight idealization
is quantitatively adequate in this market (Figure S2).
Diversification satisfies Definition 1's boundedness (√*N* cap);
its role is structural, not dynamically directional. The mechanism
tests below test collapse-frequency predictions, not the
diversification-to-*h* correspondence.
The next three subsections test whether the
model's predicted collapse-frequency curve (Figure 1) tracks actual
tail-risk and stabilizer-failure metrics at different coupling
levels — a test of the mechanism, not of the variance-algebra
identity.

### Tail-risk frequency tracks the model's collapse curve

Equation 3 is variance algebra and does not test the model's
mechanism. We binned the 5,414 windows by ρ̄ deciles and
computed the fraction per decile where the equal-weight portfolio's
maximum drawdown exceeded 5% within the 60-day window. Empirical
tail-risk rises from 0.01 at the lowest decile to 0.94 at the
highest, tracking the model's P(collapse) curve with Spearman
ρ = **+0.81** (*p* = 0.004) after normalizing both to [0, 1]
(Figure 6). The result is robust across thresholds (3%, 7%, 10%;
ρ = +0.80, +0.69, +0.75 with *p* = 0.006, 0.026, 0.012
respectively). The model captures rank order
and direction; the quantitative shape is steeper than the data
(log-log slope 2.51 vs. 1.45 at 5%). The empirical curve saturates
toward the P = 1 ceiling at the highest deciles (0.94 at decile 9)
while the model's curve is still rising in the tested *J* range —
the slope mismatch reflects ceiling saturation, not failure of
coupling-dependence. The prediction
holds out of sample: trained on 2004–2014, the model predicts
2015–2025 tail-risk frequency at Spearman ρ = +0.70 (*p* = 0.024,
17% RMSE skill over a constant-frequency baseline); the reverse split
yields ρ = +0.89 (*p* = 4.7×10⁻⁴), with held-out crises including the
COVID-19 selloff and the 2022 inflation shock. This temporal out-of-sample test partially defends against the
endogeneity concern that crises mechanically inflate ρ̄: the
coupling–tail-risk relation is stationary across both halves.
The test does not formally identify causal direction (an instrument
for ρ̄ would do that), and the affine fit has two free parameters
per fold;
we read it as evidence of a stable structural pattern, not
unconfounded causal estimation.

### Cross-mechanism confirmation: VaR exceedance

Value-at-Risk is a second passive stabilizer satisfying Definition 1:
a fixed percentile threshold (here, the historical 1% and 5%
quantiles estimated from the preceding 252 trading days) that does
not scale with coupling. Daníelsson [2002; 2013] argued that VaR
models create the systemic events they are supposed to measure;
the framework predicts that such fixed-threshold mechanisms lose
effectiveness at high coupling.

The VaR 1% exceedance rate — the fraction of days in each window
where the diversified portfolio's loss exceeds the pre-estimated
VaR — tracks the model's P(collapse) curve with Spearman ρ = **+0.87**
(*p* = 0.0009) across coupling deciles. The VaR 5% exceedance rate
yields ρ = +0.65 (*p* = 0.04). The cross-mechanism rank correlation
between the VaR-1% exceedance rate and the tail-risk frequency
P(DD > 5%) is ρ = **+0.98** (*p* = 1.5×10⁻⁶): two independent passive
stabilizers, measured on the same dataset, fail along the same
coupling axis the model predicts.

### Failure-mode transition: rigidity vs. fragmentation

The preceding tests show that passive-stabilizer failure *frequency*
tracks the model. The model makes a stronger, structurally unique
prediction: the *type* of failure changes with coupling. At low *J*
the model's collapses are fragmentation-dominated (|*m*| near 0); at
high *J* they become rigidity-dominated (|*m*| near 1, agents locked
on the wrong branch). No generic correlation-risk model makes this
prediction — it follows specifically from the Curie–Weiss
order-parameter structure.

We tested this by computing the mean pairwise sector correlation
during each crash window's drawdown period. Across coupling deciles,
drawdown-period correlation rises monotonically from **0.50** at the
lowest decile to **0.88** at the highest (Spearman ρ = +1.00 at
decile level, exact permutation *p* = 2/10! ≈ 5.5 × 10⁻⁷;
per-window ρ = +0.80). Classifying
crashes as rigidity-dominated (correlation > 0.8) or
fragmentation-dominated (< 0.3), the rigidity share rises from 0% to
100%, with empirical crossover at *J* ≈ 1.7 — within a factor of 1.8
of the model's predicted crossover at *J* ≈ 3.0 (Figure 7). The rank correlation
between the empirical drawdown-correlation profile and the model's
predicted rigidity-share-of-collapses curve is **+0.92**
(*p* = 1.3 × 10⁻⁴). This is a structural prediction about the nature
of failure, not its probability, and is the strongest defense
against the endogeneity concern: a "crashes inflate ρ̄" confound
predicts frequency only and cannot generate the rigidity-share
crossover or the +0.92 rank correlation.

### Active stabilizers: escaping the *h*/*J*² decay

If the restoring force scales with the threat — an active stabilizer
in the sense of Definition 1 — Theorem 1's decay should vanish. We
replaced *h*(*W*) = 2·(*W*/500) with two active variants:
stress-responsive *h* = *h*ₐ + α·max(0, −*m*)·*J* and
coupling-responsive *h* = *h*ₐ·(1 + α·*J*) (both have unbounded
intervention magnitude). Same 9,000-run sweep with matched seeds at
α = 2.0 (stress) / α = 1.0 (coupling): the passive high-*J* residual
(0.23 band-averaged; single-cell range 0.28–0.41, SI §S5) is
eliminated entirely under both (P(collapse) = 0.00 across *J* ≥ 4.0
at *μ* = 100; Figure 5). An α-sweep places the
minimum-strength threshold at α ≥ 0.5 for both. The active/passive classification has predictive content: the
residual is a structural consequence of bounded *h*, not of the
coordination dynamics.

### AI coupling: a heuristic proxy and its qualitative implication

Embedding the responses of 14 frontier models to 805 shared
AlpacaEval v2 prompts [Li et al. 2023] with a public
sentence-transformer yields a mean pairwise output similarity
ρ = 0.797 (rounded ρ ≈ 0.80 throughout). Through the mean-field
identification ρ̄ = *J*/(1 + *J*) this places the proxy at
*J*ₑff = 3.92 (≈ 4.0; Figure S3; SI §S8). The mapping is heuristic:
cosine similarity of sentence-transformer embeddings is not a direct
measurement of the Curie–Weiss order-parameter correlation. What can
be said is that the proxy sits well above the model's critical
threshold *J*ᶜ = *T* = 2, placing frontier AI in the regime where
passive effectiveness decays. The implication is *qualitative* — the model has no
calibrated time scale, no AI-specific *μ*, and no event-category
mapping (Figure 8 shows the prediction surface). The empirical
question this leaves open is what observable AI failure modes the
regime maps onto.

## Discussion

The failure-mode transition (ρ = +0.92) is the framework's
structurally unique empirical confirmation; the underlying *h*/*J*²
scaling survives auxiliary dynamics, network topology, and four
bistable model classes. The remainder discusses governance design.

### The active/passive distinction as a governance principle

The central conceptual contribution is taxonomic. The systemic-risk
literature has inherited Perrow's [1984] framing — coupling is
dangerous, therefore avoid coupling — but in coordination-intensive
domains coupling is increased, rapidly and largely invisibly, by
deployment decisions regulators do not control. Distinguishing
mechanisms whose intervention magnitude is bounded by structural
parameters (Definition 1) from those whose response scales with
measured stress reframes the question: not whether coupling is too
high, but whether the protective architecture matches the coupling
regime. The distinction itself is not new — fiscal economics contrasts
*automatic* and *discretionary* stabilizers [McKay & Reis 2016].
What is new here is its connection to a specific scaling law
(Theorem 1) and its application beyond fiscal policy.

The framework correctly classifies known financial failures. An
exploratory analysis of Federal Reserve rate-cut episodes (*N* = 12,
2004–2025; FOMC ≥ 10 bp) shows a positive correlation
(*r* = +0.62, *p* = 0.03) between contemporaneous coupling and 90-day
post-cut SPY returns — passive effectiveness *declines* with coupling, while this active
intervention's *rises* with it, motivating (not confirming) the
active/passive distinction.
Financial circuit breakers illustrate the passive side: a fixed 7%
halt regardless of event severity, with the intraday magnet effect
[Chen et al. 2024] consistent with Theorem 1; correlated extreme-
price-movement behavior across HFTs near these thresholds
[Brogaard et al. 2018] is consistent background.

### Cross-domain extensions

In a second financial domain — a 10-token cryptocurrency basket
(2020–2025, 1,864 60-day windows) — 64.4% of windows are
supercritical (ρ̄ > 0.667, *J* > 2) versus ~25% in the S&P 500
panel; the failure-mode transition reproduces, with rigidity share
crossing 50% near the seventh ρ̄-decile (SI §S5.4, Figure S6). In
comparative politics, a panel regression across 176 countries
(V-Dem 2000–2025, country and year fixed effects) finds the
predicted interaction sign (β(*h* × *J*) = −0.209,
*p* = 3.2 × 10⁻⁶): institutional protections become less effective
under rising coupling [Levitsky & Ziblatt 2018; Svolik 2019]. The
result reverses under an alternative proxy, so we treat it as
suggestive. In ecology, mass-extinction rates exceeding background
despite high biodiversity [IPBES 2019; Barnosky et al. 2011] remain
a hypothesis.

### Implications for AI governance

The framework has direct implications for AI governance, which we
develop here as an *application* of the principle rather than a
domain in which the framework has been independently validated. The
prevailing AI-governance portfolio — capability
evaluations, red-teaming, licensing thresholds, Constitutional AI,
open-source mandates — falls almost entirely on the passive side of
Definition 1: each encodes a fixed rule or threshold whose
intervention magnitude does not scale with the coupling level. RLHF
alignment, in particular, applies a fixed training objective against
a fixed preference dataset — a bounded restoring force in the
framework's terms. As frontier models homogenize on shared corpora, benchmarks, and
training, the heuristic proxy places AI above the critical
threshold; the framework qualitatively predicts that passive
effectiveness decays.
*Universal jailbreaks* are consistent with this prediction. Zou et
al. [2023] showed that a single adversarial suffix, optimized on
open-source models, transfers to GPT-3.5 (87%), GPT-4 (47%),
Claude-1 (48%), and PaLM-2 (66%) — bypassing the independent
safety training of each. (The contemporaneous Claude-2, with
additional defensive training, resisted the same attack at 2%; the
heterogeneity is itself informative — passive defenses can
asymmetrically harden against a known attack but still leave the
broader vulnerability surface intact.) The original mechanism is shared training data and gradient
similarity; in the active/passive frame, fixed-strength alignment
fields cannot adjust to a shared vulnerability surface that grows
with cross-model coordination. The two descriptions are compatible
rather than identical, and the present paper offers the second
reading: universal jailbreaks instantiate the structural mismatch
between bounded passive intervention and supercritical coupling.

Any mechanism whose response depends on (*m*, *J*, *t*) is exempt
from the *h*/*J*² decay. Dynamic API rate limiting that scales with
measured cross-model output correlation is one illustrative concept —
the closest AI analogue to Fed responsive intervention. Real active
stabilizers face their own constraints (zero lower bound;
political/debt limits; AI implementation and adversarial costs) —
"active" is necessary, not sufficient. Theorem 1 identifies
the structural class; verifying specific designs is future work. The
dynamic-coupling results refine the "gradualism" recommendation: slow
deployment is protective above the viability margin but
counterproductive below it.

### Limitations and open questions

The mean-field idealization is tested against a six-topology ABM
(*N* = 200, 50 seeds per cell with a 100-seed headline-cell re-run
in SI §S9): the *h*/*J*² scaling survives on every topology, and at
the headline cell all six topologies converge to a 0.33–0.37 band
with no pair separated above one SE — no tested topology, including
boundary-localized stabilization, rescues the residual. Strategic feedback and adaptive
rewiring are untested. A high-*J* sensitivity sweep (SI §S5.1, 4,800
runs) shows the residual is robust to log-vs-linear *h*(*W*) and
substantially *worse* under constant *h* — the baseline is the less
pessimistic case. Under
additive noise at the headline cell (*J* = 5, *μ* = 100), P(collapse)
is 0.20 with rigidity share 0.30: the residual asymmetry persists,
but the rigidity-typed signature is noise-prescription-specific
(SI §S2.6). The asymmetry is confirmed
across voter, Kuramoto, and cubic-Landau variants (SI §S5.3);
compartmental and adaptive-network dynamics remain future work.

Empirical mechanism tests track the model in rank order (Spearman
ρ = 0.81–0.93); the slope mismatch is diagnosed in Results as a
saturation asymmetry.

The active/passive distinction is binary in the analysis but lies on
a spectrum in real institutions. The Federal Reserve correlation
(*N* = 12) is exploratory: the sign and ordering of the J-conditional
effects are consistent with the framework, but the interaction term
is not significant at conventional levels (*p* = 0.30). The V-Dem
interaction reverses under an alternative proxy and the within-decile
scatter is null in both specifications, so this domain is reported
as "consistent with" rather than "confirms". The active variants
tested (two linear forms) are specific; the α ≥ 0.5 threshold is
form-specific. The AI-coupling proxy embeds outputs with a
sentence-transformer, capturing semantic similarity rather than full
behavioral coupling; AlpacaEval is English-only, frontier-only, and
skews toward instruction-following. Critically, the ρ̄ → *J*
mapping via the Curie–Weiss self-consistency identification is
heuristic for embedding similarity; the *J* value indicates
supercritical coupling rather than calibrating a dynamical variable.
The model does not represent any specific AI system.

The empirical rigidity-share crossover (*J* ≈ 1.7) sits below the
model's prediction (*J* ≈ 3.0) by a factor of 1.8 — the same
saturation asymmetry that shifts the log-log slope. Out-of-sample
skill is reported as rank correlation; on RMSE the reverse-split
prediction does not improve on a constant-residual baseline, so the
OOS evidence is for ordering, not level. A separate convexity
prediction for P(tail | *J*) over the tested *J* range is not
supported: the curve is closer to concave before saturating. Network
tests use *N* = 200 at mean degree ≈ 20; sparser regimes and
adaptive rewiring are untested.

Natural next steps: adaptive network rewiring; endogenous *J*–*W*
coupling reflecting crisis-driven correlation spikes; sandboxed
active-stabilizer implementation. The framework makes a concrete forward prediction: in any
system where coupling rises faster than active-stabilizer
architecture, a coordinated multi-agent failure should precede any
single-agent diagnostic. Frontier AI — supercritically coupled with a passive safety
portfolio — is the most exposed test case.

The central lesson the framework suggests is that as coupling rises,
the question is not whether protections are strong enough but whether
they are the right *kind* — and our results indicate that the current
portfolio, in finance demonstrably and in AI plausibly, faces a
regime where passive design is structurally mismatched.

## Methods

### Model specification

The minimal model couples a Curie–Weiss stochastic differential
equation in the order parameter *m*(*t*) ∈ (–1, 1) to a deterministic
wealth equation *W*(*t*):

$$
\dot m = -m + \tanh\!\left(\tfrac{J\,m + h(W)}{T}\right) + \xi\,\sqrt{1-m^2}\,\eta(t), \qquad \dot W = \mu \cdot \mathrm{employment}(m) - \mathrm{consumption}(W),
$$

with employment(*m*) = (1 + *m*)/2, consumption(*W*) = 30 +
10·(*W*/500), and the passive stabilizer field *h*(*W*) = 2·(*W*/500).
We interpret the noise term in the Itô sense throughout. Fixed
parameters: temperature *T* = 2, noise amplitude *ξ* = 0.5,
integration step d*t* = 0.05. Initial conditions *m*(0) = 0,
*W*(0) = 100. A clipping rule *m* ∈ [–1 + ε, 1 – ε] with ε = 10⁻⁶
keeps the multiplicative noise term numerically well-defined at
saturation. The multiplicative form √(1 − *m*²) is the canonical
prescription for bounded order parameters and matches the empirical
fact that financial returns and opinion variables show volatility
compression at extremes. The ε = 10⁻⁶ buffer prevents division-by-zero
artifacts without affecting the dynamics, as confirmed by the Level 0
mean-field validation (Supplementary Figure S1) and the convergence
check (§S4). The sweep
parameters are the coupling *J* and the income multiplier *μ*.

### Theorem 1 derivation

The deterministic fixed-point equation *m* = tanh[(*Jm* + *h*)/*T*]
has saturated solutions at *m* ≈ ±1 for *J* ≫ *T*. Linearization
around each saturated branch and a free-energy analysis of the
symmetry-broken landscape give a branch asymmetry — the preferential
stability of the protective (+m) branch over the unprotective (–m)
branch — that scales as *h*/*J*² to leading order. Because *h* is
bounded above by the equilibrium wealth, itself bounded by *μ*, the
asymmetry vanishes at high coupling regardless of how large *μ* is
made. The full derivation, including the contribution of the
multiplicative noise term and the noise-induced corrections to the
basin geometry, is in Supplementary Information §S2.

### Numerical integration and convergence

We integrate using an Euler–Maruyama scheme at d*t* = 0.05 over
20,000 steps (1,000 physical time units) per run. The choice of d*t*
follows a three-point convergence test at the residual-rigidity
corner (*J* = 5.0, *μ* = 100) and the fragmentation-cured corner
(*J* = 0.5, *μ* = 100), using d*t* ∈ {0.1, 0.05, 0.025} with the
simulation horizon held constant in physical time. P(collapse) at
the rigidity corner takes values 0.25, 0.35, 0.35 across the three
step sizes; the rigidity share of the residual is d*t*-stable at
88–89%; the fragmentation-cured corner
returns 0.00 at every d*t*. We adopt d*t* = 0.05 as the smallest step
at which the high-*J* residual is saturated within the 100-seed
sampling standard error. The convergence script
(`scripts/convergence_check.py`) is in the repository.

### Sweep design and collapse criteria

The static-coupling sweep covers *J* ∈ {0.5, 1.0, 1.5, 2.0, 2.5, 3.0,
3.5, 4.0, 4.5, 5.0} × *μ* ∈ {20, 30, 40, 50, 60, 70, 80, 90, 100}
with 100 independent seeds per cell. Seeds are assigned
deterministically as *s* = *s*₀ + 1000003·*i* + 1009·*j*, where *i*
and *j* index the (*J*, *μ*) grid. A run is *collapsed* if *W* falls
below 10 (1/10 of the initial value *W*(0) = 100) for 200 consecutive
steps (10 physical time units). A collapse is *rigidity-typed* if |*m*| > 0.9
at the collapse instant, *fragmentation-typed* if |*m*| < 0.3, and
*mixed* otherwise.

### Ablation cascade

The cascade preserves the (*J*, *μ*) grid, seeds, and integration
protocol while adding dynamics one level at a time. **L0** is the
bare Curie–Weiss SDE with fixed exogenous *h* (no economic feedback);
its long-time mean tracks the mean-field fixed points and validates
the analytical backbone of Theorem 1. **L1** is the minimal model.
**L2** adds replicator–mutator dynamics over *K* = 6 belief types
(mutation rate 0.01, *m*-coupled fitness strength *γ* = 1) and an
epistemic-coherence factor EC = −∑*f*ₖ log *f*ₖ / log *K* that
modulates the passive field as *h*ₑff = *h* · EC. **L3** adds a
Minsky–Keen four-phase credit cycle (Expansion → Euphoria → Minsky
moment → Deleveraging) gated by a debt-to-wealth ratio (thresholds
0.5, 1.0, 0.1) with phase-dependent income factors (1.00, 1.10,
0.50, 0.85), a 30% wealth haircut at the Minsky moment, and a 1.5×
noise amplitude during the Minsky moment. **L4** is the full simulator from prior work, archived separately
and referenced here rather than re-run; specification in SI §S1.3.

### Dynamic-coupling protocols

The dynamic-J experiment replaces static *J* with *J*(*t*) under
three schedules: *gradual* (linear ramp from 0.5 to 5.0 over the
simulation horizon), *exponential* (S-curve over the same horizon),
and *sudden* (step at *t* = 200, equivalent to 10 physical time
units). The ramp-speed sweep parameterizes the ramp duration over
{1,000, 2,000, 4,000, 10,000, 20,000, 30,000} steps with a fixed
total horizon of 40,000 steps so even the slowest ramp leaves 10,000
steps of post-ramp settling time.

### Diversification-failure analysis

Daily total-return data for nine US sector ETFs (XLB, XLE, XLF, XLI,
XLK, XLP, XLU, XLV, XLY) from 30 March 2004 onward, with XLRE
appended on its 2015-10-08 inception (yielding a 10-sector panel
post-2015), through 30 December 2025. For each 60-day rolling
window we compute the average pairwise correlation ρ̄ across
constituents and the diversification benefit DB = σ̄_single /
σ_portfolio for an equal-weight portfolio of those constituents
(SPY is loaded for separate diagnostic regressions in the
rigidity-vs-fragmentation analysis but is *not* a portfolio
constituent here). The structural curve
DB(ρ̄, *N*) = 1 / √(1/*N* + (1 − 1/*N*) ρ̄) is the analytical variance
ratio under equal-weight, equal-volatility,
equal-pairwise-correlation idealizations; deviations from it are
reported as mean absolute error.

### Tail-risk and VaR exceedance analysis

For the tail-risk test, each 60-day window is assigned a binary
"tail event" indicator for each of four drawdown thresholds (3%, 5%,
7%, 10%): 1 if the maximum peak-to-trough drawdown of the
equal-weight portfolio within the window exceeds the threshold,
0 otherwise. Windows are binned into coupling deciles by ρ̄
(equivalently *J* = ρ̄/(1 − ρ̄)), and the fraction of tail-event
windows per decile is the empirical P(tail | *J*). The model's
P(collapse | *J*) curve, extracted from the extended simulation
sweep (*J* ∈ [0.5, 10]), is normalized to [0, 1] for shape comparison.
Spearman rank correlation quantifies the match. For the VaR test,
the historical 1% and 5% VaR of the equal-weight portfolio is
estimated from the 252 trading days preceding each window. The
exceedance rate is the fraction of days within the window where the
portfolio return falls below the VaR threshold. Windows beginning
before the 253rd trading day are excluded (4,910 of 5,414 windows
retained). Binning and comparison follow the same decile procedure.

### Failure-mode classification

For each window containing a tail event (max drawdown > 5%), the
drawdown period is defined as the trading days from the peak to the
trough of the maximum drawdown within the window. The rigidity score
is the mean pairwise Pearson correlation of daily sector-ETF returns
across drawdown-period days. A crash is classified as
rigidity-dominated if this correlation exceeds 0.8, and
fragmentation-dominated if below 0.3. The coupling-decile analysis
uses the same *J* = ρ̄/(1 − ρ̄) binning as the tail-risk test. The
out-of-sample prediction splits the panel at 2015-01-01 and fits a
two-parameter affine mapping P_emp = *a* · P_model + *b* on each
half's *J*-binned tail-risk frequencies, applied to the held-out
half.

### Active-stabilizer experiment

The active-stabilizer experiment uses the same minimal-model SDE
(Equations 1–2), integration protocol (Euler–Maruyama, d*t* = 0.05,
20,000 steps), (*J*, *μ*) grid, seed formula, and collapse criteria
as the passive baseline. Two active variants replace
*h*(*W*) = 2·(*W*/500): a stress-responsive form
*h* = *h*ₐ + α·max(0, −*m*)·*J*, and a coupling-responsive form
*h* = *h*ₐ·(1 + α·*J*), where *h*ₐ = 2·(*W*/500) is the passive
baseline. The matched-seed design ensures every cell in the active
sweep receives identical noise paths as the passive baseline,
isolating the effect of the stabilizer form from stochastic
variation. The α-sweep covers α ∈ {0.1, 0.25, 0.5, 1.0, 2.0, 5.0}
over the high-*J* band (*J* ∈ {4.0, 4.5, 5.0}) × full *μ* grid with
100 seeds per cell. The extended J-sweep for the AI-coupling overlay
covers *J* ∈ {6, 7, 8, 9, 10} × *μ* ∈ {20, 40, 60, 80, 100} with 100
seeds per cell, appended to the base sweep for a continuous
*J* ∈ [0.5, 10] grid.

### Network agent-based model

The network extension replaces the scalar order parameter *m* with
N = 200 coupled agents, each carrying an individual state
*m*ᵢ ∈ (−1, 1) evolved by the same Curie–Weiss SDE as Equation 1,
with the mean-field term *Jm* replaced by *J* · mean(*m*ⱼ for
*j* ∈ neighbors(*i*)) — the local mean field on the graph.
Multiplicative noise uses a single shared η(t) per step (rationale
in SI §S9). Wealth *W* is shared across agents and depends on the
population-averaged *m*. Six graph topologies were tested: complete
(*N* − 1 = 199 edges per node; mean-field control), Erdős–Rényi
(*p* = 0.10), Watts–Strogatz (*k* = 20, *β* = 0.1), Barabási–Albert
(*m*_attach = 10), planted-partition modular (4 communities of 50,
*p*_in = 0.35, *p*_out = 0.01), and a modular variant with *h*
applied only to agents with cross-community edges (~5% of agents).
The sweep covers *J* ∈ {0.5, 1, 2, 3, 4, 5} × *μ* ∈ {20, 40, 60,
80, 100} with 50 seeds per cell and 20,000 steps at d*t* = 0.05.
Collapse criteria are identical to the scalar model, applied to the
population-averaged *m*.

### AI output-coupling measurement

Pre-generated responses from 14 frontier and near-frontier models to
805 shared instructions (AlpacaEval v2 benchmark [Li et al. 2023])
were embedded with the `all-MiniLM-L6-v2` sentence-transformer
(384 dimensions). For each instruction, the 14 × 14 pairwise
cosine-similarity matrix was computed; the overall coupling ρ is the
mean off-diagonal entry averaged across all 805 instructions. The
mapping *J* = ρ/(1 − ρ) follows the mean-field identification
ρ̄ = *J*/(1 + *J*) used throughout. Per-model, per-category, and
tier-level breakdowns are in SI §S8.

### Statistical inference

Rank-correlation tests use Spearman ρ via `scipy.stats.spearmanr`
(two-sided, mid-rank tie correction). Inference at the decile level
(*N* = 10) is reported with Fisher-*z* 95% confidence intervals
(SE_*z* = 1/√(*N* − 3) = 0.378). For the headline failure-mode
transition test, ρ = +0.92 has 95% CI [+0.71, +0.98]; for tail-risk
vs. model, ρ = +0.81 has CI [+0.37, +0.95]; for VaR-1% vs. model,
ρ = +0.87 has CI [+0.53, +0.97]; for the cross-mechanism
ρ = +0.98, CI [+0.92, +0.995]. Per-window Spearman correlations
(*N* = 5,414 for tail-risk; 2,869 for crash classification) are
reported in §S10 alongside the per-decile values. We declare the
failure-mode transition rank correlation (ρ = +0.92) as the primary
inferential test. The full family of 15 inferential tests is the
primary plus 14 secondary tests: tail-risk vs. model at four
drawdown thresholds (3%, 5%, 7%, 10%); VaR-1% and VaR-5% vs. model;
the cross-mechanism Spearman of VaR-1% against tail-risk frequency;
OOS forward and reverse split Spearman; the rank correlation between
*J*-decile and mean within-crash correlation; the crypto rigidity
Spearman; the panel regression β(*h*×*J*) interaction in the V-Dem
analysis; the Fed-rate-cut Pearson correlation between coupling and
post-cut SPY return; and the Fed-funds-cut response × *J*
interaction regression. The AlpacaEval cosine similarity
ρ̄ = 0.797 is a measured quantity, not an inferential test, and is
reported without an associated p-value. The Bonferroni-adjusted
threshold at family-wise α = 0.05/15 ≈ 0.0033 is survived by the
primary test (*p* = 1.3×10⁻⁴), the cross-mechanism ρ
(*p* = 1.5×10⁻⁶), the VaR-1% ρ (*p* = 9.5×10⁻⁴), the reverse-OOS
ρ = +0.89 (*p* = 4.7×10⁻⁴), and the *J*-decile-vs.-correlation
rank test (exact-permutation *p* = 2/10! ≈ 5.5×10⁻⁷, two-sided
count). The tail-risk vs.-model ρ at the 5% threshold (uncorrected
*p* = 0.004) is borderline; forward-OOS, VaR-5%, the Fed Pearson
*r* = +0.62 (uncorrected *p* = 0.03), and the Fed-regression
interaction (*p* = 0.30) do not survive and are reported as
exploratory, in keeping with the qualifications already stated in
Limitations.

### Data availability

The empirical inputs supporting the findings of this study are of
three types. Daily equity and sector-ETF total-return data
(2004–2025) were obtained from Yahoo Finance via the `yfinance`
Python package and are reproducible by running
`empirical/scripts/fetch_prices.py`. Federal Reserve and macro
series (FEDFUNDS, NFCI) were obtained from FRED via `fredapi`
(API key required from the user; series identifiers documented
in the script). V-Dem indicator series and OWID mirrors are
redistributed in this repository under their respective licenses
(`empirical/data/`). AlpacaEval v2 model-output JSON files for
14 frontier and near-frontier models are redistributed in
`empirical/data/alpacaeval_outputs/` under the AlpacaEval license
(Li et al., 2023). All processed CSVs underlying the figures are
in `empirical/results/` and `simulation/results/`. The complete
dataset and all generated outputs are archived at the repository
listed under "Code availability" below and on Zenodo
(DOI: 10.5281/zenodo.20061432).

### Code availability

All simulation source, sweep CSV outputs, figure-generation
scripts, and the empirical analysis pipeline are released under
the MIT License at the public GitHub repository
(`https://github.com/dancing4am/failure-mode-transition`) and
permanently archived on Zenodo (DOI: 10.5281/zenodo.20061432).
The headline
numerical results are reproducible from a clean checkout by
running, in order, `simulation/scripts/run_sweep.py`,
`simulation/scripts/level2_with_beliefs.py`,
`simulation/scripts/level3_with_credit.py`,
`simulation/scripts/dynamic_j.py`,
`simulation/scripts/ramp_speed_sweep.py`, then the empirical
pipeline beginning with
`empirical/scripts/diversification_failure.py`,
`empirical/scripts/tail_risk_frequency.py`,
`empirical/scripts/var_exceedance.py`,
`empirical/scripts/rigidity_fragmentation.py`,
`empirical/scripts/out_of_sample_test.py`, and
`empirical/scripts/ai_coupling_direct.py`. Random seeds are set
deterministically in every script via the formula
*s* = *s*₀ + 1,000,003·*i* + 1,009·*j* (with *s*₀ = 0xC0DE
for the minimal-model sweep), guaranteeing bit-exact
reproducibility of all reported figures and tables.

### Author contributions

C.J. conceived the study, derived the theoretical results,
implemented the simulation and empirical-analysis code,
performed the analyses, generated the figures, and wrote the
manuscript.

### Acknowledgments

The author thanks the AlpacaEval team for releasing the model-
output JSON files used in the AI output-coupling proxy, and the
V-Dem Institute and Our World in Data for the political indicator
mirrors. No external funding supported this work.

### Competing interests

The author declares no competing interests.

### Corresponding author

Correspondence and material requests should be addressed to
Chowon Jung (dancing4am@gmail.com).

## Figures

**Figure 1.** *Static asymmetry in the minimal model.* P(collapse)
heatmap over the (*J*, *μ*) grid (upper panel), with cell-level
collapse-type composition (lower panels: rigidity share, fragmentation
share). 9,000 runs at d*t* = 0.05; 100 seeds per cell. Fragmentation
occupies the low-*J*, low-*μ* corner and is fully cured by raising
*μ*; rigidity dominates the high-*J* band at every *μ* and saturates
above zero (residual P(collapse) = 0.23 at *μ* = 100, 90%
rigidity-typed).

**Figure 2.** *Ablation cascade L1 → L2 → L3.* Upper panel: P(collapse)
heatmap comparison across the three levels of model complexity, on
the same (*J*, *μ*) grid as Figure 1. Lower panel: rigidity-share
heatmaps for the same three levels. Adding belief dynamics (L2) and
credit-cycle dynamics (L3) monotonically *worsens* the high-*J*
residual collapse rate (23% → 30% → 34% at *μ* = 100), while
preserving rigidity dominance of the residual at every level (90%,
83%, 81%).

**Figure 3.** *Network robustness of the h/J² scaling.* P(collapse)
heatmaps over the (*J*, *μ*) grid for six network topologies with
N = 200 agents: complete graph (mean-field control), Erdős–Rényi,
Watts–Strogatz, Barabási–Albert, planted-partition modular
(4 communities), and modular with passive stabilizer at community
boundaries only. The high-*J* residual survives on every topology. The 50-seed sweep
spans 0.24–0.38 across topologies; a 100-seed re-run at the headline
cell (*J* = 5, *μ* = 100; SI §S9) shows all six topologies converge
to P(collapse) ∈ [0.33, 0.37], with no pair separated above one
standard error. No tested topology — including modular and
boundary-localized stabilization — rescues the mean-field residual.
50 seeds per cell for the full grid; d*t* = 0.05.

**Figure 4.** *Critical ramp speed.* Upper panel: P(collapse) as a
function of ramp duration (log scale), one curve per economic margin
*μ*. Three regimes are visible: resilient (*μ* ≥ 60, ≤ 7% at every
ramp), critical (*μ* = 50, monotonically increasing in ramp duration;
50% threshold ≈ 9,200 steps), unviable (*μ* ≤ 40, saturated at 1.00
beyond ramp ≈ 4,000 steps). Lower panel: collapse-type composition
(rigidity, fragmentation, mixed) versus ramp duration, broken out by
*μ*. The dominant failure mode swaps from rigidity (fast ramps) to
fragmentation (slow ramps) at low margins. 5,400 runs at d*t* = 0.05.

**Figure 5.** *Active vs. passive stabilizer effectiveness.*
P(collapse) heatmaps over the (*J*, *μ*) grid for three conditions
with matched seeds: passive baseline (left), stress-responsive active
stabilizer with α = 2.0 (center), and coupling-responsive active
stabilizer with α = 1.0 (right). The passive baseline reproduces the
high-*J* residual (band-mean 0.23 at *μ* = 100, 90% rigidity-typed);
both active variants eliminate it entirely (P(collapse) = 0.00 across
*J* ≥ 4.0 at *μ* = 100). 9,000 runs per condition at d*t* = 0.05.

**Figure 6.** *Mechanism test: tail-risk frequency and VaR exceedance
vs. coupling.* Empirical P(max drawdown > 5%) (red points) and
historical 1% VaR exceedance rate (purple diamonds) across coupling
deciles, overlaid with the model's P(collapse | *J*, *μ* = 100) curve
(blue line). All three track the same coupling axis. Spearman
ρ = +0.81 (tail-risk vs. model) and +0.87 (VaR-1% vs. model);
cross-mechanism ρ = +0.98 between the two empirical curves.
The model's curve is steeper in log-*J* than the empirical data
(log-log slope 2.51 vs. 1.45 at the 5% threshold), reflecting the
mean-field branch-selection mechanism's overstatement of the rate at
which passive-stabilizer effectiveness decays at intermediate
coupling — the direction and rank order are robust, the absolute
slope is model-specific. 5,414 windows (tail risk), 4,910 windows
(VaR, shorter due to 252-day lookback).

**Figure 7.** *Failure-mode transition: rigidity vs. fragmentation.*
Empirical share of rigidity-typed crashes per decile (red circles;
crash-window mean pairwise correlation > 0.8) and share of
fragmentation-typed crashes (blue squares; correlation < 0.3),
together with the rigidity share among classified crashes (purple
triangles, dashed) and the model's predicted rigidity-share-of-
collapses curve (blue triangles, *μ* = 100). The mean within-crash
pairwise correlation (not plotted, given in the per-decile table)
rises monotonically from 0.50 (lowest coupling) to 0.88 (highest);
the empirical rigidity-share crossover (*J* ≈ 1.7) is within a
factor of 1.8 of the model's predicted crossover (*J* ≈ 3.0).
Spearman ρ(*J*-decile, mean within-crash correlation) = +1.00
(perfect monotone; exact permutation *p* = 2/10! ≈ 5.5 × 10⁻⁷);
Spearman correlation between empirical within-crash correlation and
model rigidity-share = +0.92 (*p* = 1.3 × 10⁻⁴). 2,869 crash windows
(max drawdown > 5% within the 60-day window), ~287 per decile.

**Figure 8.** *Heuristic AI coupling proxy on the model's
prediction surface.* Panel A: P(collapse) versus coupling *J* for
three economic margins (*μ* = 40, 60, 100). The dashed red vertical
line marks the AlpacaEval output-similarity proxy (ρ = 0.797,
*J*ₑff = 3.92; 14 frontier models × 805 prompts; heuristic, not
calibrated). Dotted lines show the two earlier proxy estimates:
adversarial-attack transferability (*J* ≈ 1.2, green) and
benchmark-vector correlation (*J* ≈ 8.7, purple). Black dashed
line: critical threshold *J*ᶜ = *T* = 2. Financial reference lines
at *J* = 1.0 (average) and *J* = 9.0 (crisis). Panel B: passive
stabilizer effectiveness (1 − P(*J*, *μ* = 100)/P(*J*, *μ* = 20))
versus *J* with the same overlays, plus the theoretical *h*/*J*²
envelope (purple dotted, normalized). At the proxy-located coupling,
the model's P(collapse) surface gives 0.19 with passive
effectiveness 81% at *μ* = 100; this is the model's reading of the
proxy, not a calibrated AI prediction (see "AI coupling" subsection
caveats). Sweep extended to *J* = 10 (Methods).
