"""Figure: log10 dissolved fraction (= inorganic labile fraction; no organic matter) against pH, one panel per element.

Colour = water type; solid = HFO 1 mg Fe/L, dashed = HFO 0.1 mg Fe/L (main model, D&M + carbonate adsorption);
dotted (oxyanion panels only) = HFO 1 mg Fe/L without carbonate adsorption, to show that sensitivity.
Grey band = pH 6.5-8.0, the main slope window.
Input : outputs/models_own_curves.csv     Output: outputs/models_own_curves.png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model_setup as M  # noqa: E402

OUT = os.path.join(M.ROOT, "outputs")
COLORS = {"soft": "#2a78d6", "moderate": "#eb6834", "hard": "#1baf7a"}   # validated categorical slots 1-3
WLAB = {"soft": "soft (alk 0.78 meq/L, Ca 5.3 mg/L)", "moderate": "moderate (alk 2, Ca 30)",
        "hard": "hard (alk 4, Ca 60)"}
TITLES = {"Mn": "Mn(II)", "Ni": "Ni", "Co": "Co(II)", "Zn": "Zn", "Cd": "Cd", "Cu": "Cu(II)", "Pb": "Pb",
          "P": "P (phosphate)", "Mo": "Mo (molybdate)", "As": "As(V) (arsenate)", "V": "V(V) (vanadate)",
          "Se": "Se(IV) (selenite)"}
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"


def main():
    usecols = ["model", "water", "Fe_mgL", "pH", "element", "quantity", "value"]
    c = pd.read_csv(os.path.join(OUT, "models_own_curves.csv"), usecols=usecols)
    c = c[(c.quantity == "dissolved_fraction") & (c.Fe_mgL > 0)]
    order = ["Mn", "Ni", "Co", "Zn", "Cd", "Cu", "Pb", "P", "Mo", "As", "V", "Se"]
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.linewidth": 0.6})
    fig, axes = plt.subplots(3, 4, figsize=(13, 9.2), sharex=True)
    for ax, el in zip(axes.flat, order):
        ax.axvspan(6.5, 8.0, color="#f0efec", zorder=0, lw=0)
        ax.grid(True, color=GRID, lw=0.6, zorder=1)
        for water in ["soft", "moderate", "hard"]:
            for fe, ls in [(1.0, "-"), (0.1, "--")]:
                g = c[(c.model == "DLM_HFO_minteqv4+CO3") & (c.water == water) & (c.Fe_mgL == fe) & (c.element == el)]
                g = g.sort_values("pH")
                ax.plot(g.pH, np.log10(g.value), ls=ls, color=COLORS[water], lw=1.5, zorder=3)
            if el in M.OXYANIONS:
                g = c[(c.model == "DLM_HFO_minteqv4_noCO3") & (c.water == water) & (c.Fe_mgL == 1.0) & (c.element == el)]
                g = g.sort_values("pH")
                ax.plot(g.pH, np.log10(g.value), ls=":", color=COLORS[water], lw=1.5, zorder=2)
        ax.set_title(TITLES[el], color=INK, fontsize=10, loc="left", fontweight="bold")
        ax.set_xlim(6.0, 9.0)
        lo, hi = ax.get_ylim()
        if hi - lo < 0.1:      # keep near-flat panels readable but honest
            mid = 0.5 * (lo + hi)
            ax.set_ylim(min(lo, mid - 0.05), max(hi, mid + 0.05))
        if el == "V":
            ax.text(0.03, 0.05, "low confidence: one surface species,\nsource of constant not documented",
                    transform=ax.transAxes, fontsize=7.5, color=INK2, va="bottom")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("log10 dissolved fraction")
    for ax in axes[-1, :]:
        ax.set_xlabel("pH (set by pCO2 at fixed alkalinity)")
    handles = [Line2D([], [], color=COLORS[w], lw=2, label=WLAB[w]) for w in COLORS]
    handles += [Line2D([], [], color=INK2, lw=1.5, ls="-", label="HFO 1 mg Fe/L"),
                Line2D([], [], color=INK2, lw=1.5, ls="--", label="HFO 0.1 mg Fe/L"),
                Line2D([], [], color=INK2, lw=1.5, ls=":", label="HFO 1 mg Fe/L, no carbonate adsorption (oxyanions)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.005))
    fig.suptitle("Dzombak & Morel (1990) hydrous ferric oxide model, no organic matter, 16 °C: "
                 "fraction of each element not sorbed\n(dissolved = free + inorganic complexes; grey band = pH 6.5-8.0 "
                 "slope window; all elements present together)", fontsize=10.5, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    path = os.path.join(OUT, "models_own_curves.png")
    fig.savefig(path, dpi=200, facecolor="white")
    print("wrote", path)


if __name__ == "__main__":
    main()
