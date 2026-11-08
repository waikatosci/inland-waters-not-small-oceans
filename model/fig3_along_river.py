"""Fig. 3. Annual-mean Coulombic departure, D = PZC(T) - pH, at twelve Waikato mainstem
sites (2025 reference state, from the 2022 site survey). Numbers above bars give the
supersaturation ratio S. Brown bars: Taupo Volcanic Zone reach."""
import numpy as np
import pandas as pd
from fig_common import *
from waikato_sites import load_sites

sites = load_sites(reference="survey")
proj = pd.read_csv("results/projection.csv")
order = list(sites.site)
d25 = (proj[(proj.year == 2025) & (proj["mode"] == "additive") & (proj.ssp == "SSP2-4.5")]
       .set_index("site").loc[order])

fig, ax = plt.subplots(figsize=(DOUBLE_COL * 0.62, SINGLE_COL * 0.62))
x = np.arange(len(order))
cols = [COL_BROWN_400 if s in ("Ohaaki Br.", "Ohakuri") else COL_TEAL_600 for s in order]
ax.bar(x, d25.D, color=cols, width=0.65, lw=0)
for xi, s0, dv in zip(x, sites.S0, d25.D):
    ax.text(xi, dv + 0.015, f"{s0:.1f}", ha="center", va="bottom", fontsize=5, color=COL_BG_700)
ax.set_xticks(x, order, rotation=40, ha="right")
ax.set_ylim(0, 1.0)
ax.set_ylabel("D = PZC − pH (pH units)")
ax.text(0.99, 0.97, "S above bars", transform=ax.transAxes, ha="right", va="top", fontsize=5.5,
        color=COL_BG_600)
style_ax(ax)
save(fig, "fig3_along_river")
