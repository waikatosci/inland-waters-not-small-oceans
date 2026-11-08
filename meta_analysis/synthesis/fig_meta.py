"""Composite synthesis figure: pH responses of trace metals and oxyanions in fresh waters, as medians and
interquartile ranges. (a) divalent cations and Fe, (b) oxyanions and arsenite, (c) divalent cations by evidence
class as a fraction of the electrostatic limit. Run meta_ph.py first."""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, "../../model")
import os; os.chdir(os.path.dirname(os.path.abspath(__file__)))
from fig_common import *
from matplotlib.lines import Line2D

P = pd.read_csv("meta_pooled.csv")
import meta_ph as mp

MEAS, LAB = "Published measurements in natural water", "Published laboratory and model results"
# rows within each element band, top to bottom
ROWS = [(MEAS, COL_TEAL_800, 0.27), (LAB, COL_BROWN_400, 0.09), (mp.MHFO, COL_BG_300, -0.09), (mp.MND, COL_BG_600, -0.27)]
LIMIT = {"P": 1.5, "As": 1.5, "As(III)": 0.0, "V": 1, "Mo": 2, "W": 2}
LABEL = {MEAS: "Published measurements in natural water (field, diel, pH experiment)",
         LAB: "Published laboratory sorption, soil solution and free-ion models",
         mp.MHFO: "Speciation model, iron oxide only (dissolved)",
         mp.MND: "Speciation model, iron oxide and organic matter (free ion)"}


def bar(ax, y, r, col, xlim, lw=3.0):
    lo, hi = max(r.q25, xlim[0]), min(r.q75, xlim[1])
    ax.plot([lo, hi], [y, y], color=col, lw=lw, alpha=0.45, solid_capstyle="butt", zorder=2)
    ax.plot(np.clip(r["median"], *xlim), y, "o", ms=3.6, color=col, mfc=col, zorder=3, clip_on=False)


def forest(ax, elements, xlim, limit_line=None, limits=None):
    E = P[P.level == "element"]
    for i, e in enumerate(elements):
        y0 = len(elements) - 1 - i
        if i % 2 == 0:
            ax.axhspan(y0 - 0.5, y0 + 0.5, color=COL_BG_50, lw=0, zorder=0)
        for grp, col, dy in ROWS:
            g = E[(E.element == e) & (E.cls == grp)]
            if len(g):
                bar(ax, y0 + dy, g.iloc[0], col, xlim)
        m = E[(E.element == e) & (E.cls == MEAS)]
        if len(m):
            ax.text(1.02, y0, f"{int(m.iloc[0].k)} ({int(m.iloc[0].studies)})", ha="left", va="center", fontsize=5.2,
                    color=COL_BG_700, transform=ax.get_yaxis_transform())
        if limits and e in limits:
            ax.plot([limits[e]] * 2, [y0 - 0.48, y0 + 0.48], color=COL_RED_900, lw=0.8, ls=(0, (3, 1.5)), zorder=1)
    if limit_line is not None:
        ax.axvline(limit_line, color=COL_RED_900, lw=0.8, ls=(0, (3, 1.5)), zorder=1)
    ax.axvline(0, color=COL_BG_300, lw=0.5, zorder=1)
    ax.set_yticks(range(len(elements)), elements[::-1])
    ax.set_ylim(-0.55, len(elements) - 0.45); ax.set_xlim(*xlim)
    ax.tick_params(axis="y", length=0)
    ax.text(1.02, 1.0, "Values\n(studies)", transform=ax.transAxes, fontsize=5.2, color=COL_BG_700, va="bottom")
    style_ax(ax)


fig = plt.figure(figsize=(DOUBLE_COL, 122 * MM))
gs = fig.add_gridspec(2, 2, width_ratios=[1.12, 1], height_ratios=[1.45, 1], hspace=0.27, wspace=0.30)

ax = fig.add_subplot(gs[0, 0])
forest(ax, ["Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Cd", "Pb"], (-2.3, 0.6), limit_line=-2)
ax.set_xlabel("d log$_{10}$(concentration) / d pH")
ax.set_title("Cations", fontsize=6.5, loc="left", pad=3)
ax.text(-1.94, 5.5, "Electrostatic\nlimit, −2", fontsize=5.2, color=COL_RED_900, va="center")
panel_label(ax, "a", x=0.99, y=1.075, ha="right")

ax = fig.add_subplot(gs[0, 1])
forest(ax, ["P", "As", "As(III)", "V", "Mo", "W", "U"], (-1.0, 2.2), limits=LIMIT)
ax.text(0.62, 4.0, "As(III) uncharged:\nexpected slope 0", fontsize=5.0, color=COL_RED_900, va="center")
ax.set_xlabel("d log$_{10}$(concentration) / d pH")
ax.set_title("Oxyanions", fontsize=6.5, loc="left", pad=3)
panel_label(ax, "b", x=0.99, y=1.075, ha="right")

hand = [Line2D([], [], color=c, lw=3.0, alpha=0.45, marker="o", ms=3.6, mfc=c, mec=c, label=LABEL[g]) for g, c, _ in ROWS]
hand += [Line2D([], [], color=COL_RED_900, lw=0.8, ls=(0, (3, 1.5)), label="Electrostatic limit")]
lax = fig.add_subplot(gs[1, 1]); lax.axis("off")
lax.legend(handles=hand, loc="center left", frameon=False, fontsize=5.4, handlelength=2.0, borderaxespad=0,
           labelspacing=0.9, title="Median (point) and interquartile range (bar)", title_fontsize=5.4,
           alignment="left")

# (c) fraction of the limit by evidence class, divalent cations without Fe
ax = fig.add_subplot(gs[1, 0])
C = P[P.level == "class"].set_index("cls")
cls_order = [mp.MHFO, mp.MND, mp.FRAC, mp.CONC, mp.AMD, mp.DIELT, mp.DIEL, mp.EXP, mp.SOIL, mp.SORB, mp.FREEM]
SHORT = {mp.MHFO: "Speciation model, iron oxide only", mp.MND: "Speciation model, with organic matter",
         mp.FRAC: "Field: labile or free fraction", mp.CONC: "Field: labile or dissolved concentration",
         mp.AMD: "Acid drainage, corrected for dilution", mp.DIELT: "Diel cycle, temperature held constant",
         mp.DIEL: "Diel cycle, pH and temperature together", mp.EXP: "pH experiment, natural water",
         mp.SOIL: "Soil solution at fixed total metal", mp.SORB: "Sorption at fixed total (natural oxides)",
         mp.FREEM: "Published free-ion models (WHAM, NICA-Donnan)"}
GCOL = {mp.MHFO: COL_BG_300, mp.MND: COL_BG_600}
GCOL.update({c: COL_TEAL_800 for c in [mp.FRAC, mp.CONC, mp.AMD, mp.DIELT, mp.DIEL, mp.EXP]})
GCOL.update({c: COL_BROWN_400 for c in [mp.SOIL, mp.SORB, mp.FREEM]})
C = C.loc[[c for c in cls_order if c in C.index]]
for i, (cls, r) in enumerate(C.iterrows()):
    y = len(C) - 1 - i
    col = GCOL[cls]
    ax.plot([-r.q75 / 2, -r.q25 / 2], [y, y], color=col, lw=3.0, alpha=0.45, solid_capstyle="butt")
    ax.plot(-r["median"] / 2, y, "o", ms=3.6, color=col, mfc=col)
    ax.text(1.02, y, f"{int(r.k)} ({int(r.studies)})", va="center", fontsize=5.2, color=COL_BG_700,
            transform=ax.get_yaxis_transform())
for yb in (len(C) - 2.5, len(C) - 8.5):
    ax.axhline(yb, color=COL_BG_100, lw=0.5)
ax.set_yticks(range(len(C)), [SHORT[c] for c in C.index][::-1], fontsize=5.6)
ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.6, len(C) - 0.4)
ax.axvline(0, color=COL_BG_300, lw=0.5); ax.axvline(1, color=COL_RED_900, lw=0.8, ls=(0, (3, 1.5)))
ax.set_xlabel("Fraction of the electrostatic limit, φ = −slope / 2")
ax.text(1.02, 1.0, "Values\n(studies)", transform=ax.transAxes, fontsize=5.2, color=COL_BG_700, va="bottom")
ax.tick_params(axis="y", length=0)
panel_label(ax, "c", x=0.99, y=1.13, ha="right")
style_ax(ax)
save(fig, "Fig_meta_pH_synthesis")
