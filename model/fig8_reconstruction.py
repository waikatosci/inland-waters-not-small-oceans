"""Fig. 8. Reconstruction and projection of annual-mean conditions at Lake Karapiro,
1850 to 2100: (a) pH; (b) ferrihydrite PZC; (c) Coulombic departure D = PZC - pH;
(d) dissolved phosphate relative to 1850 from the electrostatic term alone, in the
strong-sorption limit (dissolved proportional to 1/Kd, Kd proportional to 10^(1.5 D)).
Grey band: the monitoring period, 1990 to 2025, with the change across it given in each panel.

Forcing. Atmospheric CO2: Law Dome and NOAA global means to 2025 (approximate annual
values, interpolated), then SSP2-4.5, SSP3-7.0 and SSP5-8.5 (Meinshausen et al. 2020).
Water temperature: 1.2 C of warming from 1850 to 2025, scaled to the CO2 forcing
(ln CO2; ASSUMED), then NZ projections to 2100 (Ministry for the Environment), with
water warming equal to air warming. Seasonal cycle 16 +/- 5.5 C at 2025.
The site keeps its 2025 alkalinity (51.4 mg/L as CaCO3) and its 2025 CO2 excess, 70%
metabolic and 30% geothermal, both scaled with temperature (MMRT for the metabolic part,
GPP/ER = 0.49, the median for rivers, held constant; gas exchange and solubility for both)."""
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from fig_common import *
import atmosphere as atm
import unified as u
import run_projection as rp
from waikato_sites import load_sites

# Approximate global annual-mean CO2 (ppm): Law Dome ice core to 1958, NOAA thereafter.
HIST = {1850: 285.2, 1870: 287.6, 1890: 293.9, 1910: 299.7, 1930: 307.2, 1950: 311.3,
        1960: 316.9, 1970: 325.7, 1980: 338.8, 1990: 354.4, 2000: 369.7, 2010: 389.9,
        2015: 401.0, 2020: 414.2, 2025: 427.7}
_hist = PchipInterpolator(list(HIST), list(HIST.values()))
DT_HIST = 1.2                       # warming 1850 -> 2025 (C), ASSUMED
SCEN = ["SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]
SSP_COL = {"SSP2-4.5": COL_TEAL_800, "SSP3-7.0": COL_DORANGE, "SSP5-8.5": COL_RED_900}
Z_PO4 = -1.5


def co2(year, ssp):
    if year <= 2025:
        return float(_hist(year))
    # SSP trajectory, offset so it joins the observed 2025 value
    return float(atm.co2_atm_ppm(year, ssp) - atm.co2_atm_ppm(2025, ssp) + HIST[2025])


def dT(year, ssp):
    """Water temperature relative to 2025 (C)."""
    if year <= 2025:
        f = np.log(co2(year, ssp) / HIST[1850]) / np.log(HIST[2025] / HIST[1850])
        return DT_HIST * (f - 1)
    return rp.warming(year, ssp)


site = load_sites().set_index("site").loc["Karāpiro"]
A0 = HIST[2025]


def state(year, ssp):
    a, d = co2(year, ssp), dT(year, ssp)
    out = [u.site_state(site.alk_mgL, Tm + d, a, dp_geo0=site.dp0 * (1 - site.f_met),
                        dp_met0=site.dp0 * site.f_met, T0_C=rp.T_MEAN, g=rp.G)
           for Tm in rp.T_CYCLE]
    pH = np.mean([o["pH"] for o in out])
    D = np.mean([o["D"] for o in out])
    pzc = np.mean([float(u.pzc(Tm + d)) for Tm in rp.T_CYCLE])
    pco2 = np.mean([o["pco2w"] for o in out])
    return dict(year=year, ssp=ssp, co2=a, dT=d, pco2w=pco2, pH=pH, pzc=pzc, D=D)


years = sorted(set(range(1850, 2101, 5)) | {1990, 2025})
rows = [state(y, s) for s in SCEN for y in years if y >= 2025 or s == SCEN[0]]
df = pd.DataFrame(rows)
hist = df[df.year <= 2025]
D1850 = hist.loc[hist.year == 1850, "D"].item()
df["po4_rel"] = 100 * (10 ** (Z_PO4 * (df.D - D1850)) - 1)     # % change in dissolved PO4
df.to_csv("results/reconstruction_karapiro.csv", index=False)
hist = df[df.year <= 2025]

fig, axs = plt.subplots(2, 2, figsize=(DOUBLE_COL * 0.8, SINGLE_COL * 1.05), sharex=True)
axs = axs.ravel()
r90, r25 = (hist[hist.year == y].iloc[0] for y in (1990, 2025))
cols = [("pH", "pH", "a", 0.92, (0.03, 0.05, "left", "bottom"), "{:+.3f}"),
        ("pzc", "Ferrihydrite PZC", "b", 0.92, (0.03, 0.05, "left", "bottom"), "{:+.3f}"),
        ("D", "D = PZC − pH (pH units)", "c", 0.02, (0.97, 0.05, "right", "bottom"), "{:+.3f}"),
        ("po4_rel", "Dissolved phosphate,\nchange from 1850 (%)", "d", 0.92, (0.03, 0.05, "left", "bottom"), "{:+.0f}% (relative to 1990)")]
for ax, (c, lab, letter, lx, (tx, ty, ha, va), fmt) in zip(axs, cols):
    ax.axvspan(1990, 2025, color=COL_BG_100, lw=0, zorder=0)
    ax.plot(hist.year, hist[c], color=COL_BG_900, lw=0.9, label="Hindcast")
    for s_ in SCEN:
        f = df[(df.ssp == s_) & (df.year >= 2025)]
        ax.plot(f.year, f[c], color=SSP_COL[s_], lw=0.9, label=s_)
    ax.set_xlim(1850, 2100)
    ax.set_xticks([1850, 1900, 1950, 2000, 2050, 2100])
    ax.set_ylabel(lab)
    dv = r25[c] - r90[c] if c != "po4_rel" else 100 * (10 ** (Z_PO4 * (r25.D - r90.D)) - 1)
    ax.text(tx, ty, "1990–2025: " + fmt.format(dv).replace("-", "−"),
            transform=ax.transAxes, ha=ha, va=va, fontsize=5.5, color=COL_BG_700)
    panel_label(ax, letter, x=lx)
    style_ax(ax)
for ax in axs[2:]:
    ax.set_xlabel("Year")
axs[0].text(2007.5, 0.97, "Monitoring", transform=axs[0].get_xaxis_transform(),
            ha="center", va="top", fontsize=5, color=COL_BG_700)
axs[0].legend(loc="lower left", frameon=False, bbox_to_anchor=(0.0, 0.12))
fig.tight_layout(h_pad=0.8, w_pad=1.6)
save(fig, "fig8_reconstruction")


def at(y, s="SSP2-4.5"):
    return df[(df.year == y) & ((df.ssp == s) | (df.year <= 2025))].iloc[0]


print("year  co2  dT  pco2w  pH  pzc  D  po4%")
for y in (1850, 1990, 2025):
    r = at(y); print(y, f"{r.co2:.0f} {r.dT:+.2f} {r.pco2w:.0f} {r.pH:.3f} {r.pzc:.3f} {r.D:.3f} {r.po4_rel:+.1f}")
for s in SCEN:
    r = at(2100, s); print(2100, s, f"{r.co2:.0f} {r.dT:+.2f} {r.pco2w:.0f} {r.pH:.3f} {r.pzc:.3f} {r.D:.3f} {r.po4_rel:+.1f}")
a, b, c = at(1850), at(1990), at(2025)
print(f"1850-2025: dpH {c.pH-a.pH:+.3f} dPZC {c.pzc-a.pzc:+.3f} dD {c.D-a.D:+.3f}")
print(f"1990-2025: dpH {c.pH-b.pH:+.3f} dPZC {c.pzc-b.pzc:+.3f} dD {c.D-b.D:+.3f} "
      f"share of 1850-2025 dD {(c.D-b.D)/(c.D-a.D):.0%}")
for s in SCEN:
    e = at(2100, s)
    print(f"{s}: 1850-2100 dD {e.D-a.D:+.3f}; monitoring share {(c.D-b.D)/(e.D-a.D):.0%}; "
          f"PO4 1990-2025 {100*(10**(Z_PO4*(c.D-b.D))-1):+.1f}%")
