"""Full J × mult sweep for the minimal model.

90 cells × 100 seeds × 10,000 steps. Saves long-form raw_results.csv.
"""

from pathlib import Path
import time

import numpy as np
import pandas as pd

from minimal_model import sweep, J_VALUES, MULT_VALUES, N_SEEDS, N_STEPS

OUT = Path(__file__).resolve().parents[1] / "results" / "minimal_model"
OUT.mkdir(parents=True, exist_ok=True)

print(f"Sweep: {len(J_VALUES)} J × {len(MULT_VALUES)} mult × {N_SEEDS} seeds × {N_STEPS} steps")
print(f"Total runs: {len(J_VALUES) * len(MULT_VALUES) * N_SEEDS:,}")

t0 = time.time()
cells = sweep()
elapsed = time.time() - t0
print(f"Sweep done in {elapsed:.1f}s")

rows = []
for cell in cells:
    for s in range(N_SEEDS):
        rows.append({
            "J": cell.J,
            "mult": cell.mult,
            "seed_idx": s,
            "collapsed": int(cell.collapsed[s]),
            "collapse_step": int(cell.collapse_step[s]),
            "collapse_type": int(cell.collapse_type[s]),
            "final_m": float(cell.final_m[s]),
            "final_wealth": float(cell.final_wealth[s]),
        })

df = pd.DataFrame(rows)
df.to_csv(OUT / "raw_results.csv", index=False)
print(f"Wrote {OUT / 'raw_results.csv'}  ({len(df):,} rows)")

# Quick summary table
summary = (
    df.groupby(["J", "mult"])
      .agg(
          p_collapse=("collapsed", "mean"),
          p_rigidity=("collapse_type", lambda s: (s == 1).mean()),
          p_fragmentation=("collapse_type", lambda s: (s == 2).mean()),
          p_mixed=("collapse_type", lambda s: (s == 3).mean()),
          mean_final_m=("final_m", "mean"),
          mean_final_wealth=("final_wealth", "mean"),
      )
      .reset_index()
)
summary.to_csv(OUT / "cell_summary.csv", index=False)
print(f"Wrote {OUT / 'cell_summary.csv'}  ({len(summary)} cells)")

# Compact console digest
print("\nP(collapse) by cell:")
piv_collapse = summary.pivot(index="J", columns="mult", values="p_collapse")
print(piv_collapse.round(2).to_string())

print("\nP(rigidity collapse) by cell:")
piv_rig = summary.pivot(index="J", columns="mult", values="p_rigidity")
print(piv_rig.round(2).to_string())

print("\nP(fragmentation collapse) by cell:")
piv_frag = summary.pivot(index="J", columns="mult", values="p_fragmentation")
print(piv_frag.round(2).to_string())
