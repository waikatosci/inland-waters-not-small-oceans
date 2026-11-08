"""Figure 4. Bottom-water daily DO at Lake Karapiro (profiler, ~21.5 m) observed in 2024-26 and in 2100 under
three SSPs and three cases (solubility only; plus respiration at 0.62 eV; empirical bottom-water coupling),
observed daily series shifted by the projected change. Boxes: interquartile range; whiskers: 5th-95th
percentiles; dots: daily minimum. Run oxygen.py first."""
import os, sys
sys.path.insert(0, "../model")
sys.path.insert(1, "../model")
os.chdir("../redox/")
import numpy as np, pandas as pd
from fig_common import *   # noqa: F401,F403  (style, palette, save -> ./figures)

daily = pd.read_csv("oxygen_daily_2100.csv", index_col=0)

fig, ax = plt.subplots(figsize=(SINGLE_COL * 1.15, SINGLE_COL * 0.72))

# ---- bottom-water oxygen
SSP = ["SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]
CASES = ["(i) solubility", "(ii) solubility + respiration", "(iii) empirical coupling"]
CLAB = ["Solubility only", "Solubility + respiration (0.62 eV)", "Empirical coupling (upper bound)"]
CCOL = [COL_TEAL_300, COL_TEAL_600, COL_TEAL_900]


def box(ax, x, v, c, w=0.22):
    q5, q25, q50, q75, q95 = np.percentile(v, [5, 25, 50, 75, 95])
    ax.add_patch(plt.Rectangle((x - w / 2, q25), w, q75 - q25, fc=c, ec=COL_BG_900, lw=0.3, zorder=3))
    ax.plot([x - w / 2, x + w / 2], [q50, q50], color="white" if c != COL_BG_100 else COL_BG_900, lw=0.7, zorder=4)
    ax.plot([x, x], [q5, q25], color=COL_BG_900, lw=0.4, zorder=2); ax.plot([x, x], [q75, q95], color=COL_BG_900, lw=0.4, zorder=2)
    ax.scatter([x], [v.min()], s=4, color=COL_BG_900, zorder=4, lw=0)


for y, lab in ((4, "4 mg L$^{-1}$"), (2, "2 mg L$^{-1}$")):
    ax.axhline(y, color=COL_RED_900 if y == 2 else COL_BG_300, lw=0.4, ls=(0, (3, 2)), zorder=1)
    ax.text(-0.45, y + 0.15, lab, fontsize=5, color=COL_STAT, ha="left")
box(ax, 0, daily["observed"].values, COL_BG_100, w=0.3)
for i, s in enumerate(SSP):
    for j, c in enumerate(CASES):
        box(ax, i + 1 + (j - 1) * 0.27, daily[f"{s} {c}"].values, CCOL[j])
ax.set_xticks(range(4)); ax.set_xticklabels(["Observed\n2024-26"] + [s + "\n2100" for s in SSP], fontsize=5.5)
ax.set_xlim(-0.5, 3.55); ax.set_ylim(0, 13.6)
ax.set_ylabel("Bottom-water DO, daily mean (mg L$^{-1}$)")
for j, (c, l) in enumerate(zip(CCOL, CLAB)):
    ax.add_patch(plt.Rectangle((0.56, 12.95 - j * 0.62), 0.13, 0.36, fc=c, ec=COL_BG_900, lw=0.3))
    ax.text(0.75, 13.13 - j * 0.62, l, fontsize=5.2, va="center")
style_ax(ax)
save(fig, "Fig_redox_oxygen")
