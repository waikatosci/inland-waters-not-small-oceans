"""Fig. 4. (a) All Gran alkalinity titrations of Waikato River water, Oct 2024 to Aug 2026
(fixed sites S1-S4 and mobile-sonde sites), with the median curve and 5-95% envelope.
(b) Equilibrium pH against supersaturation S for each titrated alkalinity at 16 C and
2025 atmospheric CO2 (427.7 ppm), with the ferrihydrite PZC. Top axis: atmospheric CO2
at which water at equilibrium (S = 1) would have the same pH.

Excluded: Otago lakes, the Ruakuri karst control (MS3), and titrations failing QC
(Gran r2 < 0.98, end point below pH 3.3, or alkalinity outside 0.3-1.5 meq/L; this
removes the Jul-Aug 2025 sets analysed in Sep 2025 with a faulty probe)."""
import numpy as np
import pandas as pd
from fig_common import *
import matplotlib.ticker as mt
import carbonate_system as cs
import unified as u

cur = pd.read_csv("results/titration_curves.csv")
m = pd.read_csv("results/titration_meta.csv")
keep = ((m.sheet != "OTAGO FEB2025") & (m.site != "MS3") & (m.gran_r2 >= 0.98)
        & (m.min_pH >= 3.3) & m.alk_meqL.between(0.3, 1.5))
m = m[keep].copy()
m["mobile"] = m.site.str.startswith("MS")
cur = cur[cur.key.isin(m.key)]
mob = m.set_index("key").mobile

fig, axs = plt.subplots(1, 2, figsize=(DOUBLE_COL * 0.72, SINGLE_COL * 0.62))

# a: titration curves
ax = axs[0]
grid = np.linspace(0, 1.2, 241)
Y = []
for k, c in cur.groupby("key"):
    c = c.sort_values("acid_meqL").drop_duplicates("acid_meqL")
    ax.plot(c.acid_meqL, c.pH, lw=0.25, color=COL_BROWN_400 if mob[k] else COL_TEAL_400,
            alpha=0.55 if mob[k] else 0.12, zorder=3 if mob[k] else 2)
    Y.append(np.interp(grid, c.acid_meqL, c.pH, left=np.nan, right=np.nan))
Y = np.array(Y)
ok = np.mean(~np.isnan(Y), axis=0) > 0.5
with np.errstate(all="ignore"):
    lo, med, hi = (np.nanpercentile(Y[:, ok], q, axis=0) for q in (5, 50, 95))
ax.fill_between(grid[ok], lo, hi, color=COL_BG_400, alpha=0.35, lw=0, zorder=4)
ax.plot(grid[ok], med, color=COL_BG_900, lw=0.9, zorder=5)
ax.set_xlim(0, 1.2); ax.set_ylim(3.5, 8.3)
ax.set_xlabel("Acid added (meq L$^{-1}$)")
ax.set_ylabel("pH")
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
ax.legend(handles=[Line2D([], [], color=COL_TEAL_400, lw=0.8, label=f"Sites S1–S4 (n = {(~m.mobile).sum()})"),
                   Line2D([], [], color=COL_BROWN_400, lw=0.8, label=f"Mobile-sonde sites (n = {m.mobile.sum()})"),
                   Line2D([], [], color=COL_BG_900, lw=0.9, label="Median"),
                   Patch(color=COL_BG_400, alpha=0.35, label="5–95% range")],
          loc="lower left", frameon=False)
panel_label(ax, "a", x=0.90)
style_ax(ax)

# b: pH against S
ax = axs[1]
T, PA = 16.0, 427.7
S = np.logspace(0, 1, 60)
P = np.array([[cs.pH_from_Alk_pCO2(a * 50.04, s * PA * 1e-6, T)["pH"] for s in S] for a in m.alk_meqL])
for row, mb in zip(P, m.mobile):
    ax.plot(S, row, lw=0.25, color=COL_BROWN_400 if mb else COL_TEAL_400, alpha=0.55 if mb else 0.12)
plo, pmed, phi = (np.percentile(P, q, axis=0) for q in (5, 50, 95))
ax.fill_between(S, plo, phi, color=COL_BG_400, alpha=0.35, lw=0)
ax.plot(S, pmed, color=COL_BG_900, lw=0.9)
pzc = float(u.pzc(T))
ax.axhline(pzc, color=COL_GREEN_800, lw=0.7, ls=(0, (4, 2)))
ax.text(9.5, pzc + 0.02, f"Ferrihydrite PZC, {T:.0f} °C", ha="right", va="bottom", fontsize=5.5,
        color=COL_GREEN_800)
WAIKATO_S = (2.1, 4.6)   # range of S at 12 main-stem sites, field pH on 25 Oct 2022
ax.axvspan(*WAIKATO_S, color=COL_BG_50, lw=0, zorder=0)
ax.text(np.sqrt(WAIKATO_S[0] * WAIKATO_S[1]), 8.2, "Waikato main stem", ha="center", va="top", fontsize=5.5, color=COL_BG_600,
        bbox=dict(facecolor=COL_BG_50, edgecolor="none", pad=0.6), zorder=4)
s585 = 1135.2 / PA
ax.axvline(s585, color=COL_RED_900, lw=0.6, ls=(0, (1, 1.5)))
label_bubble(ax, "SSP5-8.5, 2100\n(1,135 ppm)", xy=(s585, 7.82), xytext=(1.25, 7.15), ha="left")
for sv in (1, 2.1, 4.6):
    j = np.argmin(abs(S - sv))
    ax.scatter([S[j]], [pmed[j]], s=8, color=COL_BG_900, zorder=6)
    ax.annotate(f"D = {pzc - pmed[j]:.2f}", xy=(S[j], pmed[j]), xytext=(S[j] * (1.06 if sv == 1 else 1.0), pmed[j] + (0.1 if sv == 1 else 0.2)),
                ha="left" if sv == 1 else "center", va="bottom", fontsize=5.5, color=COL_BG_900, zorder=7,
                arrowprops=dict(arrowstyle="-", lw=0.4, color=COL_BG_600, shrinkA=0, shrinkB=2))
ax.set_xscale("log")
ax.set_xticks([1, 2, 3, 5, 7, 10], ["1", "2", "3", "5", "7", "10"])
ax.xaxis.set_minor_formatter(mt.NullFormatter())
ax.set_xlim(1, 10); ax.set_ylim(6.9, 8.45)
ax.set_xlabel("Supersaturation, S")
ax.set_ylabel("pH")
top = ax.secondary_xaxis("top", functions=(lambda s: s * PA, lambda p: p / PA))
top.set_xticks([500, 1000, 2000, 3000])
top.xaxis.set_major_formatter(mt.FixedFormatter(["500", "1,000", "2,000", "3,000"]))
top.xaxis.set_minor_locator(mt.NullLocator())
top.set_xlabel("Equivalent atmospheric CO$_2$ at S = 1 (ppm)")
top.tick_params(width=0.25, length=2.5, direction="in")
panel_label(ax, "b", x=0.02, y=0.10)
style_ax(ax)
fig.tight_layout(w_pad=1.2)
save(fig, "fig4_titration_S")

i1, i3, i7 = (np.argmin(abs(S - v)) for v in (1, 2.1, 4.6))
print(f"n={len(m)}; alk median {m.alk_meqL.median():.3f}; pH S=1 {pmed[i1]:.2f}, 2.1 {pmed[i3]:.2f}, "
      f"7 {pmed[i7]:.2f}; spread at S=1 {phi[i1]-plo[i1]:.3f}")
