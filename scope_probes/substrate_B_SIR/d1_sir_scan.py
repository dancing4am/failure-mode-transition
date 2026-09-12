"""Archival runner behind pilot_sir_RESULTS.txt: SIR bistability broad scan.

Runs, using ONLY the deposited substrate-B modules in this directory:
  (a) the deposited self-test script verbatim (python bistability.py),
  (b) beta_bistable_window over the two grids coded in _tests(),
  (c) beta_bistable_window over the historical broad-scan grid
      (kappa x I_h x omega x nu x beta0), using the deposited detector.
Decision rule (fixed before running): ANY bistable cell found -> STOP.

NOTE ON PROVENANCE: this is the runner named in pilot_sir_RESULTS.txt.
For archival the module paths were relativized to this directory (the
original audit run invoked the identical modules from a working copy by
absolute path); the scan logic, grids, and settings are unchanged.
Section (a) of the archived RESULTS records the self-test as it behaved
on first run (AssertionError on the module's original kappa=0.85
expectation, since corrected to assert the measured outcome).
"""
import sys, os, itertools, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from sir_core import SIRParams
from bistability import beta_bistable_window, is_bistable

print("=== pilot_sir_RESULTS.txt : Substrate B bistability broad scan ===")
print(f"date: 2026-08-24 | interpreter: {sys.version.split()[0]}")
print("command: python d1_sir_scan.py (wraps deposited bistability.beta_bistable_window)")
print()

print("-- (a) deposited self-test script, verbatim output --")
r = subprocess.run([sys.executable, "bistability.py"],
                   cwd=HERE, capture_output=True, text=True)
print(r.stdout.strip())
print(r.stderr.strip())
print(f"exit code: {r.returncode}")
print()

found = []

print("-- (b) the two grids coded in _tests() --")
bg1 = np.linspace(0.5, 8.0, 30)
w = beta_bistable_window(SIRParams(kappa=0.0), bg1)
print(f"kappa=0.00 (standard SIRS), beta0 in [0.5,8] x30 : window size = {w.size}"
      + (f"  WINDOW={w}" if w.size else "  -> no bistable cells"))
if w.size: found.append(("kappa=0", w))
bg2 = np.linspace(0.5, 8.0, 60)
w = beta_bistable_window(SIRParams(kappa=0.85, I_h=0.05, omega=0.2, nu=0.02), bg2)
print(f"kappa=0.85 I_h=0.05 omega=0.2 nu=0.02, beta0 in [0.5,8] x60 : window size = {w.size}"
      + (f"  WINDOW={w}" if w.size else "  -> no bistable cells"))
if w.size: found.append(("kappa=0.85 fixture", w))
print("NOTE: the deposited module's own test originally asserted the kappa=0.85")
print("fixture should be bistable; that assertion FAILED on first run (see (a)) -")
print("the detector finds no window there. The self-test expectation, not the")
print("scan, was wrong.")
print()

print("-- (c) historical broad-scan grid, deposited detector --")
print("grid: kappa in {0.6,0.85,0.97} x I_h in {0.01,0.05,0.1} x omega in {0.05,0.2,0.5}")
print("      x nu in {0,0.02,0.05} x beta0 in [0.5,10] (80 points)")
bg = np.linspace(0.5, 10.0, 80)
n_cfg = 0
for kappa, I_h, omega, nu in itertools.product((0.6, 0.85, 0.97),
                                               (0.01, 0.05, 0.1),
                                               (0.05, 0.2, 0.5),
                                               (0.0, 0.02, 0.05)):
    n_cfg += 1
    p = SIRParams(kappa=kappa, I_h=I_h, omega=omega, nu=nu)
    w = beta_bistable_window(p, bg)
    if w.size:
        found.append(((kappa, I_h, omega, nu), w))
        print(f"BISTABLE: kappa={kappa} I_h={I_h} omega={omega} nu={nu} window={w}")
print(f"configurations scanned: {n_cfg} x 80 beta0 points = {n_cfg*80} cells")
print(f"bistable cells found: {sum(len(w) for _, w in found)}")
print()
print("VERDICT:", "NO BISTABLE CELLS - the Discussion sentence stands"
      if not found else "BISTABLE CELLS FOUND - STOP THE ROUND")
