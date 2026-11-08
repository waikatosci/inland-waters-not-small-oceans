"""Within-site coupling of pH and dissolved inorganic phosphorus (DIP) across GLORICH stations,
by alkalinity decile. (a) Spearman rho(pH, DIP): all samples; samples below and above each site's
median discharge; and with oxygen saturation held constant (photosynthesis proxy). (b) The same for
the dissolved fraction DIP/TP, which removes dilution. Points: decile medians. Stations with
>= 30 pairs (>= 60 for the flow split). Run glorich_partial.py, glorich_controls.py,
glorich_ratio.py and glorich_baseflow.py first."""
import numpy as np, pandas as pd
from fig_common import *
fs = pd.read_csv("glorich_flowsplit_by_site.csv")
ct = pd.read_csv("glorich_controls_by_site.csv")
ra = pd.read_csv("glorich_ratio_by_site.csv")
edges = np.quantile(ra.med_Alk, np.linspace(0, 1, 11))
def dec(df, col, alk="med_Alk"):
    b = np.clip(np.searchsorted(edges, df[alk], side="right") - 1, 0, 9)
    g = df.assign(b=b).groupby("b")
    return np.log10(g[alk].median()), g[col].median()
fig, axs = plt.subplots(1, 2, figsize=(DOUBLE_COL * 0.72, SINGLE_COL * 0.62), sharey=True)
series = {
    "a": [(ra, "rho_DIP", None, COL_BG_900, "All samples"),
          (fs[fs.part == "low"], "rho_DIP", None, COL_TEAL_800, "Below median flow"),
          (fs[fs.part == "high"], "rho_DIP", None, COL_TEAL_300, "Above median flow"),
          (ct, "ctl_DO", None, COL_DORANGE, "O$_2$ saturation held constant")],
    "b": [(ra, "rho_frac", None, COL_BG_900, None),
          (fs[fs.part == "low"], "rho_frac", None, COL_TEAL_800, None),
          (fs[fs.part == "high"], "rho_frac", None, COL_TEAL_300, None)]}
for ax, (k, ss) in zip(axs, series.items()):
    ax.axhline(0, color=COL_BG_300, lw=0.4)
    ax.axvline(3, color=COL_GREEN_800, lw=0.6, ls=(0, (4, 2)))
    for df, col, _, c, lab in ss:
        x, y = dec(df.dropna(subset=[col]), col)
        ax.plot(x, y, "o-", ms=2.5, lw=0.8, color=c, mec="none", label=lab)
    ax.set_xlabel("log$_{10}$ alkalinity (µeq L$^{-1}$)")
    panel_label(ax, k, x=0.92)
    style_ax(ax)
axs[0].set_ylabel("Median within-site Spearman ρ")
axs[0].set_title("pH and DIP", fontsize=6.5, color=COL_BG_700, pad=3)
axs[1].set_title("pH and DIP/TP", fontsize=6.5, color=COL_BG_700, pad=3)
axs[1].text(3.02, 0.03, "1 meq L$^{-1}$", transform=axs[1].get_xaxis_transform(), fontsize=5.5, color=COL_GREEN_800, va="bottom")
axs[0].legend(loc="lower left", frameon=False)
fig.tight_layout(w_pad=1.0)
save(fig, "FigS_glorich_controls")
