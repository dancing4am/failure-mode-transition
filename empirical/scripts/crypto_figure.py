"""
Generate crypto comparison figure (SI Figure S6).

Shows: rigidity-share by decile for both S&P 500 and crypto.
Also shows distribution of rho_bar in each domain.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EMP = Path(__file__).resolve().parents[1]
RES = EMP / "results"
OUT = EMP / "figures" / "figure_s6_crypto_vs_sp500.png"

# Crypto data (just computed)
crypto_dec = pd.read_csv(RES / "crypto_decile_summary.csv")
crypto_crash = pd.read_csv(RES / "crypto_crash_summary.csv")
crypto_panel = pd.read_csv(RES / "crypto_panel_full.csv")

# S&P 500 data — load from the rigidity_fragmentation results
sp_rigidity_csv = RES / "rigidity_fragmentation_by_coupling.csv"
sp_panel_csv = RES / "tail_risk_by_coupling.csv"
have_sp = sp_rigidity_csv.exists() and sp_panel_csv.exists()

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

# Panel A: rho_bar distributions (histogram)
axA = axes[0]
axA.hist(crypto_panel["rho_bar"].dropna(), bins=40, alpha=0.6,
         color="#d62728", label=f"Crypto (mean={crypto_panel.rho_bar.mean():.2f})")
if have_sp:
    sp_panel = pd.read_csv(sp_panel_csv)
    if "rho_bar" in sp_panel.columns:
        axA.hist(sp_panel["rho_bar"].dropna(), bins=40, alpha=0.6,
                 color="#1f77b4", label=f"S&P 500 (mean={sp_panel.rho_bar.mean():.2f})")
axA.axvline(0.667, color="black", linestyle="--", linewidth=1.5,
            label="J=2 (critical)")
axA.set_xlabel(r"average pairwise correlation $\bar\rho$")
axA.set_ylabel("number of 60-day windows")
axA.set_title("A. Distribution of coupling")
axA.legend(loc="upper right", fontsize=9)
axA.grid(alpha=0.3)

# Panel B: rigidity-share by decile
axB = axes[1]
axB.plot(crypto_crash["rho_bar_mid"], crypto_crash["rigidity_share"],
         marker="o", color="#d62728", linewidth=2, label="Crypto")
if have_sp:
    sp_rig = pd.read_csv(sp_rigidity_csv)
    if "rho_bar_mid" in sp_rig.columns and "rigidity_share" in sp_rig.columns:
        axB.plot(sp_rig["rho_bar_mid"], sp_rig["rigidity_share"],
                 marker="s", color="#1f77b4", linewidth=2, label="S&P 500")
axB.axhline(0.5, linestyle=":", color="grey", alpha=0.5)
axB.set_xlabel(r"$\bar\rho$ at decile midpoint")
axB.set_ylabel("rigidity share among classified crashes")
axB.set_title("B. Failure-mode transition (both domains)")
axB.legend(loc="upper left", fontsize=9)
axB.grid(alpha=0.3)
axB.set_ylim(-0.02, 1.05)

# Panel C: drawdown-period correlation by decile
axC = axes[2]
axC.plot(crypto_dec["rho_bar_mid"], crypto_dec["drawdown_corr_mean"],
         marker="o", color="#d62728", linewidth=2, label="Crypto")
if have_sp:
    if "drawdown_corr_mean" in sp_rig.columns:
        axC.plot(sp_rig["rho_bar_mid"], sp_rig["drawdown_corr_mean"],
                 marker="s", color="#1f77b4", linewidth=2, label="S&P 500")
axC.set_xlabel(r"$\bar\rho$ at decile midpoint")
axC.set_ylabel("mean drawdown-period correlation")
axC.set_title("C. Crash correlation by coupling")
axC.legend(loc="upper left", fontsize=9)
axC.grid(alpha=0.3)
axC.set_ylim(0, 1.05)

fig.suptitle(
    "Figure S6. Cross-domain replication: cryptocurrency reproduces the "
    "fragmentation->rigidity transition with deeper supercritical regime.",
    fontsize=11, y=1.02
)
fig.tight_layout()
fig.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Wrote {OUT}")
