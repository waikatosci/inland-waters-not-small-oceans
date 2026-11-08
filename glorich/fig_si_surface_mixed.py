"""Fig. S7: ferrihydrite surface potential, Nernst relation against a one-pK Basic Stern model.
Fig. S8: GLORICH mixed-effects tests of the alkalinity dependence of the within-site pH-DIP slope.
Run surface_potential_bsm.py and glorich_mixed.py first (the latter writes the CSV and JSON read here)."""
import json
import numpy as np, pandas as pd
from fig_common import *
from surface_potential_bsm import psi0, PZC

# ---- Fig. S7
pH = np.linspace(6.0, 9.0, 61)
nern = 59.16 * (PZC - pH)
lo = np.array([min(psi0(p, Ns, C, I) for Ns in (3, 6, 8) for C in (0.9, 1.5) for I in (1.5e-3, 1e-2)) for p in pH]) * 1e3
hi = np.array([max(psi0(p, Ns, C, I) for Ns in (3, 6, 8) for C in (0.9, 1.5) for I in (1.5e-3, 1e-2)) for p in pH]) * 1e3
cen = np.array([psi0(p) for p in pH]) * 1e3
fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL * 0.75))
ax.axvspan(6.5, 7.8, color=COL_BG_50, lw=0, zorder=0)
ax.fill_between(pH, lo, hi, color=COL_TEAL_200, lw=0, label="Basic Stern model, parameter range")
ax.plot(pH, cen, color=COL_TEAL_800, lw=0.9, label="Basic Stern model, central")
ax.plot(pH, nern, color=COL_BG_900, lw=0.8, ls=(0, (4, 2)), label="Nernst relation")
ax.axhline(0, color=COL_BG_300, lw=0.4)
ax.axvline(PZC, color=COL_BG_300, lw=0.4)
ax.text(PZC + 0.04, -80, "PZC", fontsize=5.5, color=COL_BG_700)
ax.text(7.15, 160, "Waikato pH range", ha="center", fontsize=5.5, color=COL_BG_700)
ax.set_xlabel("pH")
ax.set_ylabel("Surface potential, ψ$_0$ (mV)")
ax.set_xlim(6.0, 9.0); ax.set_ylim(-150, 175)
ax.legend(loc="lower left", frameon=False)
style_ax(ax)
fig.tight_layout()
save(fig, "FigS7_surface_potential")

# ---- Fig. S8
r = pd.read_csv("glorich_mixed_results.csv")
D = json.load(open("glorich_mixed_D.json"))
fig, axs = plt.subplots(1, 2, figsize=(DOUBLE_COL * 0.72, SINGLE_COL * 0.62), gridspec_kw=dict(width_ratios=[1, 1.15]))
ax = axs[0]
ab = r[r.term == "pHc:soft"].reset_index(drop=True)
labels = ["No controls", "Discharge and\nseason held constant", "Base flow only"]
y = np.arange(len(ab))[::-1]
for yi, (_, row), lab in zip(y, ab.iterrows(), labels):
    ax.plot([row.beta - 1.645 * row.se, row.beta + 1.645 * row.se], [yi, yi], color=COL_TEAL_800, lw=0.9)
    ax.plot(row.beta, yi, "o", ms=3.2, color=COL_TEAL_800, mec="none")
    ax.text(0.005, yi + 0.18, f"p = {row.p:.0e}".replace("e-", " × 10$^{-").replace("p = 1 × ", "p = ") + "}$",
            fontsize=5.2, color=COL_BG_700, va="bottom")
ax.axvline(0, color=COL_BG_300, lw=0.4)
ax.set_yticks(y, labels)
ax.set_xlim(-0.02, 0.32); ax.set_ylim(-0.6, 2.7)
ax.set_xlabel("Soft minus hard water,\nslope of log$_{10}$ DIP on pH")
panel_label(ax, "a", x=0.92, y=0.08)
style_ax(ax)

ax = axs[1]
la = np.linspace(D["alk_range"][0], D["alk_range"][1], 50)
b0, b1, cv = D["b0"], D["b1"], np.array(D["cov"])
m = b0 + b1 * la
se = np.sqrt(cv[0, 0] + la ** 2 * cv[1, 1] + 2 * la * cv[0, 1])
ax.fill_between(10 ** la / 1000, m - 1.645 * se, m + 1.645 * se, color=COL_TEAL_200, lw=0)
ax.plot(10 ** la / 1000, m, color=COL_TEAL_800, lw=0.9)
ax.axhline(0, color=COL_BG_300, lw=0.4)
ax.axvline(1, color=COL_BG_300, lw=0.6, ls=(0, (4, 2)))
ax.set_xscale("log")
ax.set_xticks([0.1, 0.3, 1, 3], ["0.1", "0.3", "1", "3"])
ax.set_xlabel("Site-median alkalinity (meq L$^{-1}$)")
ax.set_ylabel("Within-site slope of\nlog$_{10}$ DIP on pH")
ax.text(1.05, ax.get_ylim()[0] + 0.02, "1 meq L$^{-1}$", fontsize=5.2, color=COL_BG_700, va="bottom")
panel_label(ax, "b", x=0.92, y=0.88)
style_ax(ax)
fig.tight_layout(w_pad=1.6)
save(fig, "FigS8_glorich_mixed")
x0 = 10 ** (-b0 / b1) / 1000
print(f"slope crosses zero at {x0:.2f} meq/L")
