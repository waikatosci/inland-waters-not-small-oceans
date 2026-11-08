"""Fig. 11 (paper Fig. S6). Carbonate saturation of Waikato River water against CO2 supersaturation S.
Saturation index SI = log10 Omega for calcite and aragonite, for every titrated alkalinity (Fig. 1)
at 16 C and 427.7 ppm, with Ca2+ at the Lake Karapiro median of 5.33 mg/L: median (line) and
5-95% range of alkalinity (band). Dotted: median alkalinity with Ca2+ lowered by 24%, the depletion
attributed to gold clam calcification downstream (Hartland et al. 2026, ES&T). Vertical lines:
S = 1, 2.1 and 4.6, as in Fig. 1b. Ksp from Plummer and Busenberg (1982); Davies activity corrections."""
import numpy as np
import pandas as pd
from fig_common import *
import carbonate_system as cs
from forward_model import si, CA_MGL

m = pd.read_csv("results/titration_meta.csv")
keep = ((m.sheet != "OTAGO FEB2025") & (m.site != "MS3") & (m.gran_r2 >= 0.98)
        & (m.min_pH >= 3.3) & m.alk_meqL.between(0.3, 1.5))
alk = m[keep].alk_meqL.values
T, PA = 16.0, 427.7
S = np.logspace(0, 1, 50)


def curve(a_meq, mineral, ca=CA_MGL):
    return np.array([si(cs.pH_from_Alk_pCO2(a_meq * 50.04, s * PA * 1e-6, T)["CO3"], T, mineral, ca) for s in S])


fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL * 0.78))
for mineral, col in (("calcite", COL_BG_800), ("aragonite", COL_BROWN_400)):
    Y = np.array([curve(a, mineral) for a in alk])
    lo, med, hi = (np.percentile(Y, q, axis=0) for q in (5, 50, 95))
    ax.fill_between(S, lo, hi, color=col, alpha=0.18, lw=0)
    ax.plot(S, med, color=col, lw=0.9, label=mineral.capitalize())
    ax.plot(S, curve(np.median(alk), mineral, CA_MGL * 0.76), color=col, lw=0.8, ls=(0, (1, 1.5)))
    print(mineral, "S=1:", round(med[0], 2), "S=2.1:", round(np.interp(2.1, S, med), 2), "S=4.6:", round(np.interp(4.6, S, med), 2))
for s_ in (1, 2.1, 4.6):
    ax.axvline(s_, color=COL_BG_300, lw=0.5, ls=(0, (3, 2)))
ax.axhline(0, color=COL_BG_300, lw=0.5)
ax.set_xscale("log")
ax.set_xticks([1, 2, 3, 5, 7, 10], ["1", "2", "3", "5", "7", "10"])
ax.set_xlabel("CO$_2$ supersaturation, S")
ax.set_ylabel("Saturation index, log$_{10}$ Ω")
ax.set_ylim(-2.6, 0.2)
from matplotlib.lines import Line2D
h, l = ax.get_legend_handles_labels()
h.append(Line2D([], [], color=COL_BG_700, lw=0.8, ls=(0, (1, 1.5)))); l.append("Ca$^{2+}$ −24% (gold clam)")
ax.legend(h, l, loc="upper right", frameon=False)
style_ax(ax)
fig.tight_layout()
save(fig, "fig11_saturation_S")
