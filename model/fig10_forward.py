"""Fig. 10 (paper Fig. 5). Forward model of Lake Karapiro, 1850 to 2100, with the surface ocean under the
same forcing for comparison.
(a) A labile divalent cation relative to 1850 at the median published fractions of the electrostatic
limit (field surveys, 0.04; diel cycles, 0.15), the limit itself (dashed), and the same limit for the surface ocean (navy).
(b) Dissolved oxyanions relative to 1850: phosphate scaled by its observed fraction (GLORICH), phosphate
and molybdate at the limit (dashed), and the phosphate limit for the surface ocean (navy).
(c) Free CO2 (left axis) and pH (right axis) in the lake and the surface ocean.
(d) Saturation index of calcite and aragonite in the lake and the surface ocean.
Lines to 2025: reconstruction. Colour band: range from SSP2-4.5 to SSP5-8.5, line SSP3-7.0, with
respiration more temperature-sensitive than photosynthesis in the lake. Dotted: lake, SSP3-7.0 with equal
activation energies. Grey band: monitoring period, 1990 to 2025.
Run forward_model.py and ocean_forward.py first."""
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.ticker as mt
from fig_common import *

SC = ("SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
dfa = pd.read_csv("results/forward_karapiro.csv")
EQ = dfa[(dfa.case == "R=P") & ((dfa.ssp == "SSP3-7.0") | (dfa.year <= 2025))].drop_duplicates("year").sort_values("year")
df = dfa[dfa.case == "R>P"]
H = df[(df.year <= 2025) & (df.ssp == "SSP2-4.5")].sort_values("year")
F = {s: df[(df.ssp == s) & (df.year >= 2025)].set_index("year") for s in SC}
oc = pd.read_csv("results/forward_ocean.csv")
OH = oc[(oc.year <= 2025) & (oc.ssp == "SSP2-4.5")].sort_values("year")
OF = {s: oc[(oc.ssp == s) & (oc.year >= 2025)].set_index("year") for s in SC}
OCEAN = COL_TEAL_900


def draw(ax, col, color, ls="-", conv=lambda v: v, lab=None, eq=True, h=None, f=None, alpha=0.18, lw=0.9):
    h = H if h is None else h
    f = F if f is None else f
    ax.plot(h.year, conv(h[col]), color=color, lw=lw, ls=ls)
    lo, mid, hi = (conv(f[s][col]) for s in SC)
    ax.fill_between(lo.index, lo, hi, color=color, alpha=alpha, lw=0)
    ax.plot(mid.index, mid, color=color, lw=lw, ls=ls, label=lab)
    if eq:
        ax.plot(EQ.year, conv(EQ[col]), color=color, lw=0.8, ls=(0, (1, 1.2)))


def ocean(ax, col, ls="-", conv=lambda v: v, lab=None, color=OCEAN, alpha=0.14):
    draw(ax, col, color, ls=ls, conv=conv, lab=lab, eq=False, h=OH, f=OF, alpha=alpha)


rel = lambda v: 1 + v / 100
fig, axs = plt.subplots(2, 2, figsize=(DOUBLE_COL * 0.86, SINGLE_COL * 1.3), sharex=True)
axs = axs.ravel()
for ax in axs:
    ax.axvspan(1990, 2025, color=COL_BG_50, lw=0, zorder=0)

ax = axs[0]
ocean(ax, "Divalent limit", ls=(0, (4, 2)), conv=rel, lab="Surface ocean, limit")
draw(ax, "Divalent limit", COL_BG_500, ls=(0, (4, 2)), conv=rel, lab="Lake, divalent cation limit")
draw(ax, "Divalent diel", COL_BROWN_500, conv=rel, lab="Lake, diel-cycle median (φ = 0.15)")
draw(ax, "Divalent field", COL_TEAL_600, conv=rel, lab="Lake, field-survey median (φ = 0.04)")
ax.set_yscale("log")
ax.set_ylim(0.9, 13)
ax.yaxis.set_major_locator(mt.FixedLocator([1, 1.5, 2, 3, 5, 8]))
ax.yaxis.set_major_formatter(mt.FixedFormatter(["1", "1.5", "2", "3", "5", "8"]))
ax.yaxis.set_minor_locator(mt.NullLocator())
ax.set_ylabel("Labile divalent cation,\nrelative to 1850")
ax.legend(loc="upper left", frameon=False)
panel_label(ax, "a", x=0.97, y=0.97, ha="right")

ax = axs[1]
ocean(ax, "Phosphate limit", ls=(0, (4, 2)), conv=rel, lab="Surface ocean, phosphate limit")
draw(ax, "Phosphate", COL_DORANGE, conv=rel, lab="Lake phosphate, observed fraction")
draw(ax, "Phosphate limit", COL_DORANGE, ls=(0, (4, 2)), conv=rel, lab="Lake phosphate, limit")
draw(ax, "Molybdate limit", COL_GREEN_800, ls=(0, (4, 2)), conv=rel, lab="Lake molybdate, limit")
ax.set_yscale("log")
ax.set_ylim(0.18, 1.45)
ax.yaxis.set_major_locator(mt.FixedLocator([0.2, 0.3, 0.5, 0.75, 1]))
ax.yaxis.set_major_formatter(mt.FixedFormatter(["0.2", "0.3", "0.5", "0.75", "1"]))
ax.yaxis.set_minor_locator(mt.NullLocator())
ax.set_ylabel("Dissolved concentration,\nrelative to 1850")
ax.legend(loc="lower left", frameon=False)
panel_label(ax, "b", x=0.97, y=0.97, ha="right")

ax = axs[2]
LAKE_C = COL_GREEN_800
draw(ax, "co2_umol", LAKE_C)
ocean(ax, "co2_umol")
ax.set_ylim(0, 100)
ax.set_ylabel("Free CO$_2$ (µmol L$^{-1}$), solid")
ax2 = ax.twinx()
draw(ax2, "pH", LAKE_C, ls=(0, (3, 1.5)), eq=False, alpha=0.10)
ocean(ax2, "pH", ls=(0, (3, 1.5)), alpha=0.08)
ax2.set_ylim(7.2, 8.6)
ax2.set_ylabel("pH, dashed")
ax2.spines["top"].set_visible(False)
ax2.tick_params(axis="y", labelsize=6, length=2, width=0.5)
for sp in ("right",):
    ax2.spines[sp].set_linewidth(0.5)
ax.text(1880, 62.5, "Lake CO$_2$", color=LAKE_C, fontsize=5.5, ha="left", va="bottom")
ax.text(1880, 13.5, "Ocean CO$_2$", color=OCEAN, fontsize=5.5, ha="left", va="bottom")
ax2.text(1880, 7.585, "Lake pH", color=LAKE_C, fontsize=5.5, ha="left", va="bottom")
ax2.text(1880, 8.19, "Ocean pH", color=OCEAN, fontsize=5.5, ha="left", va="bottom")
panel_label(ax, "c", x=0.97, y=0.97, ha="right")

ax = axs[3]
ax.axhline(0, color=COL_BG_300, lw=0.6, ls=(0, (2, 2)), zorder=1)
ocean(ax, "si_cal", lab="Surface ocean, calcite", color=COL_TEAL_600)
ocean(ax, "si_arag", lab="Surface ocean, aragonite")
draw(ax, "si_cal", COL_BG_800, lab="Lake, calcite")
draw(ax, "si_arag", COL_BROWN_400, lab="Lake, aragonite")
ax.set_ylim(-2.0, 1.0)
ax.text(1855, 0.04, "Saturation", fontsize=5, color=COL_BG_600, va="bottom")
ax.set_ylabel("Saturation index, log$_{10}$ Ω")
ax.legend(loc="center left", frameon=False, bbox_to_anchor=(0, 0.42))
panel_label(ax, "d", x=0.97, y=0.97, ha="right")

for ax in axs:
    ax.set_xlim(1850, 2100)
    ax.set_xticks([1850, 1900, 1950, 2000, 2050, 2100])
    style_ax(ax)
for ax in axs[2:]:
    ax.set_xlabel("Year")
axs[2].text(2007.5, 0.03, "Monitoring", transform=axs[2].get_xaxis_transform(), ha="center",
            va="bottom", fontsize=5, color=COL_BG_700)
handles = [Line2D([], [], color=COL_BG_500, lw=0.9, label="To 2025: hindcast. After 2025: SSP3-7.0"),
           Patch(color=COL_BG_500, alpha=0.3, lw=0, label="Range SSP2-4.5 to SSP5-8.5"),
           Line2D([], [], color=COL_BG_500, lw=0.8, ls=(0, (1, 1.2)), label="Lake, SSP3-7.0, R = P"),
           Line2D([], [], color=OCEAN, lw=0.9, label="Surface ocean, same forcing")]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(h_pad=0.8, w_pad=1.2, rect=(0, 0.05, 1, 1))
save(fig, "fig10_forward")
