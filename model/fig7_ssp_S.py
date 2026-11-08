"""Fig. 7. Projected pH and Coulombic departure in 2100 against present-day (2025)
supersaturation S, under four SSPs, for the titrated Waikato alkalinity.
(a) Annual-mean pH. (b) Change in annual-mean D = PZC(T) - pH from 2025 to 2100; right axis gives the
corresponding pH change, exact for the dashed lines. For the solid lines warming also lowers
the PZC, so pH falls less than D rises by 0.01, 0.04, 0.06 and 0.08 (SSP1-2.6 to SSP5-8.5).
Solid lines: atmospheric CO2 and warming propagated together; warming acts on the
metabolic CO2 source (MMRT, GPP/ER = 0.49, the median for rivers; respiration slightly more temperature-sensitive than
photosynthesis: 0.62 and 0.53 eV, medians of whole-ecosystem studies), gas exchange and solubility, the carbonate constants and the
ferrihydrite PZC. Dotted lines: the same with equal activation energies (0.57 eV).
Dashed lines: atmospheric CO2 alone (2025 temperatures).
Shading: 5-95% range of titrated alkalinity. Black: 2025.
The CO2 excess at each S is taken as 70% metabolic and 30% geothermal (the downstream
sites); the geothermal share responds to gas exchange and solubility only. Annual means
over a 16 +/- 5.5 C seasonal cycle. MMRT optima are typical values, not fitted."""
import numpy as np
import pandas as pd
from fig_common import *
import atmosphere as atm
import unified as u
import run_projection as rp

F_MET = 0.7
m = pd.read_csv("results/titration_meta.csv")
keep = ((m.sheet != "OTAGO FEB2025") & (m.site != "MS3") & (m.gran_r2 >= 0.98)
        & (m.min_pH >= 3.3) & m.alk_meqL.between(0.3, 1.5))
ALK = {q: np.percentile(m[keep].alk_meqL, q) * 50.04 for q in (5, 50, 95)}   # mg/L CaCO3

S0 = np.linspace(1, 10, 37)
SSP_COL = {"SSP1-2.6": COL_TEAL_300, "SSP2-4.5": COL_TEAL_800,
           "SSP3-7.0": COL_DORANGE, "SSP5-8.5": COL_RED_900}


def state(s0, alk, ssp, year, warm=True, metab=u.METAB):
    """Annual-mean pH and D for a water at present-day supersaturation s0."""
    a0 = float(atm.co2_atm_ppm(rp.YEAR0, ssp))
    a = float(atm.co2_atm_ppm(year, ssp))
    dT = rp.warming(year, ssp) if warm else 0.0
    dp = (s0 - 1) * a0
    out = [u.site_state(alk, Tm + dT, a, dp_geo0=dp * (1 - F_MET), dp_met0=dp * F_MET,
                        T0_C=rp.T_MEAN, g=rp.G, metab=metab) for Tm in rp.T_CYCLE]
    return np.mean([o["pH"] for o in out]), np.mean([o["D"] for o in out])


def curves(ssp, year, warm=True, metab=u.METAB):
    r = {}
    for q, alk in ALK.items():
        r[q] = np.array([state(s, alk, ssp, year, warm, metab) for s in S0])
    return r


runs = {"2025": curves("SSP2-4.5", 2025), "2025 EQ": curves("SSP2-4.5", 2025, metab=u.METAB_EQ)}
for ssp in rp.SSPS:
    runs[ssp] = curves(ssp, 2100)
    runs[ssp + " CO2"] = curves(ssp, 2100, warm=False)
    runs[ssp + " EQ"] = curves(ssp, 2100, metab=u.METAB_EQ)

rows = []
for k, r in runs.items():
    for i, s in enumerate(S0):
        rows.append(dict(run=k, S0=s, pH_med=r[50][i, 0], D_med=r[50][i, 1],
                         pH_lo=r[5][i, 0], pH_hi=r[95][i, 0]))
pd.DataFrame(rows).to_csv("results/ssp_vs_S.csv", index=False)

WAIKATO_S = (2.1, 4.6)   # range of S at 12 main-stem sites, field pH on 25 Oct 2022
fig, axs = plt.subplots(1, 2, figsize=(DOUBLE_COL * 0.72, SINGLE_COL * 0.62))
base = runs["2025"][50][:, 1]
base_eq = runs["2025 EQ"][50][:, 1]
for j, (ax, lab) in enumerate(zip(axs, ["pH", "Change in D, 2025–2100 (pH units)"])):
    ax.axvspan(*WAIKATO_S, color=COL_BG_50, lw=0, zorder=0)
    r = runs["2025"]
    if j == 0:
        ax.fill_between(S0, r[95][:, 0], r[5][:, 0], color=COL_BG_400, alpha=0.3, lw=0)
        ax.plot(S0, r[50][:, 0], color=COL_BG_900, lw=0.9, label="2025", zorder=5)
    else:
        ax.axhline(0, color=COL_BG_900, lw=0.6)
    for ssp in rp.SSPS:
        r, rc, re = runs[ssp], runs[ssp + " CO2"], runs[ssp + " EQ"]
        c = SSP_COL[ssp]
        if j == 0:
            ax.fill_between(S0, r[95][:, 0], r[5][:, 0], color=c, alpha=0.15, lw=0)
            y, yc, ye = r[50][:, 0], rc[50][:, 0], re[50][:, 0]
        else:
            y, yc, ye = r[50][:, 1] - base, rc[50][:, 1] - base, re[50][:, 1] - base_eq
        ax.plot(S0, y, color=c, lw=0.9, label=f"{ssp}, 2100")
        if j == 1:
            ax.plot(S0, ye, color=c, lw=0.8, ls=(0, (1, 1.2)))
        ax.plot(S0, yc, color=c, lw=0.6, ls=(0, (3, 1.5)))
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 3, 5, 7, 10], ["1", "2", "3", "5", "7", "10"])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlim(1, 10)
    ax.set_xlabel("Present-day supersaturation, S (2025)")
    ax.set_ylabel(lab)
    style_ax(ax)

axs[0].text(np.sqrt(WAIKATO_S[0] * WAIKATO_S[1]), 0.97, "Waikato\nmain stem", transform=axs[0].get_xaxis_transform(),
            ha="center", va="top", fontsize=5.5, color=COL_BG_600)
from matplotlib.lines import Line2D
h, l = axs[0].get_legend_handles_labels()
h += [Line2D([], [], color=COL_BG_600, lw=0.9), Line2D([], [], color=COL_BG_600, lw=0.8, ls=(0, (1, 1.2))),
      Line2D([], [], color=COL_BG_600, lw=0.6, ls=(0, (3, 1.5)))]
l += ["CO$_2$ and warming, R > P", "CO$_2$ and warming, R = P (b)", "CO$_2$ only"]
axs[0].legend(h, l, loc="lower left", frameon=False)
axs[0].set_ylim(6.35, 8.1)
axs[1].set_ylim(-0.02, 0.45)
# Secondary axis: pH change. Exact for the dashed lines (PZC fixed); for the solid lines
# warming also lowers the PZC, so pH falls less than D rises by 0.02 K^-1 x warming.
sec = axs[1].secondary_yaxis("right", functions=(lambda d: -d, lambda p: -p))
sec.set_ylabel("Change in pH, 2025–2100 (PZC fixed)")
sec.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: f"{v:.1f}".replace("-", "−") if abs(v) > 1e-9 else "0.0"))
sec.tick_params(width=0.25, length=2.5, direction="in")
off = ", ".join(f"{-u.DPZC_DT * rp.warming(2100, k):.2f}" for k in rp.SSPS)
fig.text(0.5, 0.01, "Right axis of b is exact for the dashed lines (PZC fixed).\nOn the solid and dotted lines warming also lowers the "
         f"PZC, so pH falls less than D rises, by {off} (SSP1-2.6 to SSP5-8.5).",
         ha="center", va="bottom", fontsize=5.5, color=COL_BG_600)
panel_label(axs[0], "a", x=0.97, y=0.97, ha="right")
panel_label(axs[1], "b", x=0.97, y=0.97, ha="right")
fig.tight_layout(w_pad=2.2, rect=(0, 0.08, 1, 1))
save(fig, "fig7_ssp_S")

# print summary
i32, i7 = (np.argmin(abs(S0 - v)) for v in (2.1, 4.6))
for i, s in ((0, 1), (i32, 2.1), (i7, 4.6)):
    base = runs["2025"][50][i, 1]; beq = runs["2025 EQ"][50][i, 1]
    print(f"S={s}: D2025 {base:.2f}; " + "; ".join(
        f"{k} dD {runs[k][50][i,1]-base:+.3f} (R=P {runs[k+' EQ'][50][i,1]-beq:+.3f}; CO2 only {runs[k+' CO2'][50][i,1]-base:+.3f})" for k in rp.SSPS))
