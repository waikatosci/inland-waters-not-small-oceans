"""Shared set-up for the talk figures: Nature Geoscience house style (ngeo_style.py),
plus helpers for saving and for the forcing-channel colours used in Figs 5 and 6."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ngeo_style import *  # noqa: F401,F403  (palette, MM, SINGLE_COL, DOUBLE_COL, helpers)
from ngeo_style import apply_style, style_ax, panel_label, label_bubble

apply_style()
# Liberation Sans is metrically identical to Arial; it is only used where Arial is absent.
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = True

OUT = "figures"
os.makedirs(OUT, exist_ok=True)

# Forcing channels (Figs 5-6). Metabolism takes the single accent colour.
CHANNELS = ["atm", "met", "phys", "carb", "pzc"]
CHANNEL_LABEL = {"atm": "Atmospheric CO$_2$", "met": "Metabolism (MMRT)",
                 "phys": "Gas exchange and solubility", "carb": "Carbonate constants",
                 "pzc": "PZC decrease with warming"}
CHANNEL_COL = {"atm": COL_TEAL_800, "met": COL_DORANGE, "phys": COL_BG_100,
               "carb": COL_BG_300, "pzc": COL_BROWN_200}


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=400)
    plt.close(fig)
    print("wrote", name)
