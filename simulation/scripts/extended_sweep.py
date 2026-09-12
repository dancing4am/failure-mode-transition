"""Extended J-sweep (J = 6-10) and the combined model prediction curve.

Regenerates
  results/extended_sweep/extended_sweep_summary.csv    (J in {6..10} x mu in {20,40,60,80,100})
  results/extended_sweep/base_sweep_summary.csv         (base 10x9 sweep merged with the extended cells)

The combined curve extends the standard minimal-model sweep beyond the
calibrated range. Its J <= 5 half is the standard minimal-model sweep
(results/minimal_model/cell_summary.csv, produced bit-exactly by run_sweep.py)
rounded to two decimals; its J > 5 half is the extended sweep computed here with
minimal_model.run_cell (same integrator, dt = 0.05, 20,000 steps, 100 seeds per
cell).

RECONSTRUCTION STATUS — statistically consistent, not bit-identical.
The script that originally produced the archived extended-sweep CSVs was lost
before the repository was assembled. This reconstruction uses the repository's
standard deterministic seed formula,

    seed = 0xC0DE + j_idx * 1_000_003 + m_idx * 1009,

with (j_idx, m_idx) enumerating the extended grid. The original script's exact
seed assignment could not be recovered: 28 candidate formulas (index bases and
offsets, alternate seed-base constants, legacy dt = 0.1 / 10,000-step
configuration) were tested against the archived per-cell counts and none
reproduced them exactly. Regenerated values therefore agree with the archived
CSVs to within binomial sampling error (SE ~ 0.05 at 100 seeds per cell) but
are not bit-identical. For that reason the archived CSVs remain authoritative
and this script, by default, writes to results/extended_sweep/regenerated_check/
and prints per-column maximum absolute deviations against the archive. Pass
--overwrite-archived to write to the canonical locations instead. The J <= 5
half of the combined CSV IS exact (it is a deterministic transform of
cell_summary.csv).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from minimal_model import run_cell, N_SEEDS, N_STEPS

ROOT = Path(__file__).resolve().parents[1]
BASE_SUMMARY = ROOT / "results" / "minimal_model" / "cell_summary.csv"
ARCHIVED_EXT = ROOT / "results" / "extended_sweep" / "extended_sweep_summary.csv"
ARCHIVED_COMB = ROOT / "results" / "extended_sweep" / "base_sweep_summary.csv"

EXT_J_VALUES = (6.0, 7.0, 8.0, 9.0, 10.0)
EXT_MU_VALUES = (20, 40, 60, 80, 100)
SEED_BASE = 0xC0DE  # standard repo constant (minimal_model.sweep)


def run_extended() -> pd.DataFrame:
    rows = []
    t0 = time.time()
    for j_idx, J in enumerate(EXT_J_VALUES):
        for m_idx, mu in enumerate(EXT_MU_VALUES):
            seed = SEED_BASE + j_idx * 1_000_003 + m_idx * 1009
            cell = run_cell(J, mu, n_seeds=N_SEEDS, n_steps=N_STEPS, seed=seed)
            rows.append({
                "J": J,
                "mu": mu,
                "p_collapse": round(float(cell.collapsed.mean()), 2),
                "p_rigidity": round(float((cell.collapse_type == 1).mean()), 2),
                "p_fragmentation": round(float((cell.collapse_type == 2).mean()), 2),
                "n_collapsed": int(cell.collapsed.sum()),
            })
            print(f"  J={J:4.1f} mu={mu:3d}  p_collapse={rows[-1]['p_collapse']:.2f}", flush=True)
    print(f"Extended sweep done in {time.time() - t0:.1f}s")
    return pd.DataFrame(rows)


def build_combined(ext: pd.DataFrame) -> pd.DataFrame:
    base = pd.read_csv(BASE_SUMMARY)
    base = base.rename(columns={"mult": "mu"})
    base = base[["J", "mu", "p_collapse", "p_rigidity", "p_fragmentation"]].copy()
    for col in ("p_collapse", "p_rigidity", "p_fragmentation"):
        base[col] = base[col].round(2)
    base["n_collapsed"] = (base["p_collapse"] * 100).round().astype(int)
    combined = pd.concat([base, ext], ignore_index=True)
    combined = combined.sort_values(["mu", "J"], kind="stable").reset_index(drop=True)
    return combined


def compare(regen: pd.DataFrame, archived_path: Path, label: str) -> None:
    if not archived_path.exists():
        print(f"[{label}] no archived file at {archived_path}; skipping comparison")
        return
    arch = pd.read_csv(archived_path)
    a = arch.sort_values(["mu", "J"], kind="stable").reset_index(drop=True)
    r = regen.sort_values(["mu", "J"], kind="stable").reset_index(drop=True)
    if len(a) != len(r):
        print(f"[{label}] ROW-COUNT MISMATCH: archived {len(a)} vs regenerated {len(r)}")
        return
    print(f"[{label}] max |deviation| vs archived, per column:")
    exact = True
    for col in ("p_collapse", "p_rigidity", "p_fragmentation", "n_collapsed"):
        dev = float(np.max(np.abs(a[col].to_numpy(float) - r[col].to_numpy(float))))
        if dev > 0:
            exact = False
        print(f"    {col:16s} {dev:.4g}")
    print(f"[{label}] {'EXACT match' if exact else 'NOT bit-identical (statistical reconstruction)'}")


def main() -> None:
    overwrite = "--overwrite-archived" in sys.argv

    ext = run_extended()
    combined = build_combined(ext)

    if overwrite:
        ext_out = ARCHIVED_EXT
        comb_out = ARCHIVED_COMB
    else:
        check_dir = ROOT / "results" / "extended_sweep" / "regenerated_check"
        check_dir.mkdir(parents=True, exist_ok=True)
        ext_out = check_dir / "extended_sweep_summary.csv"
        comb_out = check_dir / "combined_sweep_summary.csv"

    ext.to_csv(ext_out, index=False)
    print(f"Wrote {ext_out}  ({len(ext)} cells)")
    combined.to_csv(comb_out, index=False)
    print(f"Wrote {comb_out}  ({len(combined)} rows)")

    if not overwrite:
        compare(ext, ARCHIVED_EXT, "extended_sweep_summary")
        compare(combined, ARCHIVED_COMB, "combined_sweep_summary")


if __name__ == "__main__":
    main()
