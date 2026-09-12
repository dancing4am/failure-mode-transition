# Scope probes

Pilot checks of candidate substrates against the scope criterion discussed in the
manuscript: the picture applies to systems whose bistability is a cusp (a ℤ₂-symmetric
double well with a bounded linear field) traversed reversibly, and not to fold /
saddle-node / transcritical bistability with a privileged equilibrium. Each substrate
below was tested for whether it meets that criterion. The code and the raw pilot outputs
are retained here; the manuscript describes these probes qualitatively in the Discussion.

## Substrates

- **Neural network** (`substrate_A_NN/`). Probed for a clean supercritical pitchfork by
  ramping the control parameter up and down and measuring hysteresis. The transition
  showed fold / first-order behaviour — collapse coupling α_c ≈ 2.91, mean hysteresis
  ≈ 0.75 against a 0.2 threshold, with most seeds re-emerging — so it does not meet the
  cusp + reversibility criterion. **Out of scope.**

- **Compartmental epidemic, SIR** (`substrate_B_SIR/`). The bistable transition is
  fold / transcritical, with the disease-free state as a privileged equilibrium.
  **Out of scope** on structural grounds. Code: `sir_core.py`, `bistability.py`,
  `classifier.py`, `temperature.py`.

- **Power grid** (README only, no code). Second-order (swing-equation) Kuramoto
  synchronization; loss of the synchronous state is a fold / saddle-node transition.
  Classified out of scope at the literature level, so no pilot was run. **Out of scope.**

- **Degenerate optical parametric oscillator, DOPO** (`substrate_D_DOPO/`). The canonical
  ℤ₂-symmetric pitchfork, and the confirming corner case. The pilot measured a field
  exponent of 0.995 (95% CI [0.930, 1.070], excluding the fold value ½), the predicted
  variance coefficient to within 0.1% (measured 0.67646 vs predicted 0.67598), and a mean
  cross-barrier time that increases monotonically with pump. **In scope (confirming).**
  Code: `dopo_core.py`, `pilot_dopo.py`, `pilot_dopo_P1sec_hp.py`; outputs in
  `pilot_dopo_RESULTS.txt`.

## Files

Each substrate directory holds its simulation code (`*.py`) and, where a pilot was run,
its raw outcome record (`*_RESULTS.txt`).
