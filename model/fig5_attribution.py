"""Fig. 5 (paper Fig. S3). Change in annual-mean Coulombic departure D at Lake Karapiro, 2025 to 2100,
under four SSPs, divided exactly among forcing channels (Shapley decomposition), with respiration more
temperature-sensitive than photosynthesis (a; MMRT, 0.65 and 0.32 eV) and with equal activation
energies (b; 0.57 eV). Diamonds: net change. Metabolic parameters are illustrative pending fits to the
ECP metabolism record. Run run_projection.py and attribution_eq.py first."""
import numpy as np
import pandas as pd
from fig_common import *
import run_projection as rp

fig, axs = plt.subplots(1, 2, figsize=(DOUBLE_COL * 0.8, SINGLE_COL * 0.68), sharey=True)
for ax, f, title, letter in zip(axs, ["results/attribution.csv", "results/attribution_eq.csv"],
                                ["Respiration more sensitive than photosynthesis (R > P)",
                                 "Equal sensitivity (R = P)"], "ab"):
    attr = pd.read_csv(f)
    sub = attr[(attr.site == "Karāpiro") & (attr.target == "D")].set_index("ssp").loc[rp.SSPS]
    xs = np.arange(len(rp.SSPS)); pos = np.zeros(len(xs)); neg = np.zeros(len(xs))
    for c in CHANNELS:
        v = sub[c].values
        ax.bar(xs, v, bottom=np.where(v >= 0, pos, neg), width=0.55, color=CHANNEL_COL[c], lw=0,
               label=(CHANNEL_LABEL[c] if c != "met" else "Metabolism") if c != "phys" else None)
        pos += np.where(v >= 0, v, 0); neg += np.where(v < 0, v, 0)
    ax.scatter(xs, sub.total, marker="D", s=9, color=COL_BG_900, zorder=4, label="Net change")
    for xi, t in zip(xs, sub.total):
        ax.text(xi + 0.32, t, f"{t:+.2f}".replace("-", "−"), va="center", fontsize=5.5, color=COL_BG_900)
    ax.axhline(0, color=COL_BG_300, lw=0.4)
    ax.set_xticks(xs, rp.SSPS)
    ax.set_xlim(-0.5, 3.75)
    ax.set_title(title, fontsize=6, color=COL_BG_700, pad=3)
    panel_label(ax, letter, x=0.92, y=0.92)
    style_ax(ax)
axs[0].set_ylabel("Change in D, 2025–2100 (pH units)")
axs[0].legend(loc="upper left", frameon=False)
fig.tight_layout(w_pad=1.2)
save(fig, "fig5_attribution")
