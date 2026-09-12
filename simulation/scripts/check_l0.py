"""L0 sanity quantification.

At h=0 the SDE is symmetric (±m* both stable), so each seed picks a branch
by noise; ensemble mean → 0 even though every realisation tracks |m*|. Use
|m| comparison for the symmetric case, signed comparison otherwise.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
df = pd.read_csv(ROOT / "results" / "ablation" / "level0_summary.csv")

# For h=0: compare |SDE mean of |m|| ~ |mean-field|; SDE mean is 0 by symmetry.
# Better: re-run h=0 reporting mean(|m|), or just use |mean_field_m| comparison
#  with a separate column of |SDE m|.  Easier: the q25/q75 already tell us the
#  bimodal structure for h=0.

# Signed deviation (good for h>0 where there's a preferred branch):
df["signed_dev"] = (df["mean_m_sde"] - df["mean_field_m"]).abs()
# Branch-agnostic deviation: compare to the closer of ±m*
df["branch_dev"] = np.minimum(
    (df["mean_m_sde"] - df["mean_field_m"]).abs(),
    (df["mean_m_sde"] + df["mean_field_m"]).abs(),
)

print("Signed deviation (|SDE_mean - +branch|) — large at h=0 due to ±branch averaging:")
print(df.groupby("h")[["signed_dev"]].agg(["mean", "max"]))
print()
print("Branch-agnostic deviation (|SDE_mean - closer-of-±branch|):")
print(df.groupby("h")[["branch_dev"]].agg(["mean", "max"]))
print()

# For h=0 specifically, IQR width is the better diagnostic — bimodal SDE
h0 = df[df["h"] == 0.0].copy()
h0["iqr"] = h0["q75_m_sde"] - h0["q25_m_sde"]
print("At h=0, IQR width vs J (wide IQR = bimodal across seeds):")
print(h0[["J", "mean_m_sde", "q25_m_sde", "q75_m_sde", "iqr", "mean_field_m"]].to_string(index=False))
