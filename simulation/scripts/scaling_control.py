# -*- coding: utf-8 -*-
"""
Scaling-vs-fixed control for passive stabilizers (SI §S5, Figure S3).

Confirms the active/passive distinction is SCALING-vs-FIXED, not magnitude: a
larger FIXED field only postpones the high-coupling residual (its eta=4h/J still
erodes), whereas a coupling-scaled ACTIVE field holds eta roughly constant and
escapes. Matched seeds, mu=100; J extended to 20 BEYOND the calibrated J in
[0.5,5] purely to expose the mechanism. Imports the committed scalar run_cell.
Writes results/active_stabilizer/scaling_control.csv and figures/figS3.png.

Run: py -3.13 scaling_control.py
"""
import sys
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import active_stabilizer as A  # noqa: E402

MU, NSEEDS = 100, 100
JS = [3.0, 5.0, 7.0, 10.0, 15.0, 20.0]
CSV_OUT = HERE.parents[0] / "results" / "active_stabilizer" / "scaling_control.csv"
FIG_OUT = HERE.parents[1] / "manuscript" / "arxiv" / "figures" / "figS3.png"


def h_fixed(val):
    def f(W, m, J, alpha):
        return val
    return f


CONDS = [
    ("passive",   "passive (h~0.4)", A.h_passive,         0.0, 0.4),
    ("fixed2.4",  "fixed h=2.4",     h_fixed(2.4),        0.0, 2.4),
    ("fixed5.0",  "fixed h=5.0",     h_fixed(5.0),        0.0, 5.0),
    ("active",    "active (h~J)",    A.h_active_coupling, 1.0, None),
]


def main():
    results = {key: [] for key, *_ in CONDS}
    print(f"=== scaling-vs-fixed control (mu={MU}, {NSEEDS} matched seeds) ===")
    print("  J  " + "".join(f"{lbl:>16s}" for _, lbl, *_ in CONDS))
    rows = []
    for ji, J in enumerate(JS):
        seed = A.SEED_BASE + ji * 1_000_003  # matched across conditions at this J
        line = {}
        for key, lbl, hfn, alpha, _ in CONDS:
            r = A.run_cell(J, MU, hfn, alpha, n_seeds=NSEEDS, seed=seed)
            p = float(r.collapsed.mean())
            results[key].append(p); line[key] = p
            rows.append({"J": J, "condition": key, "p_collapse": p})
        print(f"{J:>5.1f}" + "".join(f"{line[k]:>16.3f}" for k, *_ in CONDS))

    CSV_OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(CSV_OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["J", "condition", "p_collapse"]); w.writeheader(); w.writerows(rows)
    print(f"wrote {CSV_OUT}")

    # ---- Figure S7 ----
    Ja = np.array(JS)
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.axvspan(5.0, 20.0, color="0.92", zorder=0)
    ax.text(12.0, 0.46, "out of calibrated range\n(mechanism illustration)",
            ha="center", va="center", fontsize=8, color="0.4")
    styles = {"passive": ("o-", "#c0392b"), "fixed2.4": ("s--", "#e67e22"),
              "fixed5.0": ("^--", "#8e44ad"), "active": ("D-", "#1b7837")}
    for key, lbl, *_ in CONDS:
        mk, col = styles[key]
        ax.plot(Ja, results[key], mk, color=col, label=lbl, lw=1.6, ms=5)
    ax.axvline(5.0, color="0.6", ls=":", lw=1)
    ax.set_xlabel("coupling  J"); ax.set_ylabel("P(collapse),  μ = 100")
    ax.set_title("Scaling-vs-fixed: a larger fixed field only postpones failure")
    ax.set_ylim(-0.03, 0.55); ax.legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(FIG_OUT, dpi=300)
    print(f"wrote {FIG_OUT}")


if __name__ == "__main__":
    main()
