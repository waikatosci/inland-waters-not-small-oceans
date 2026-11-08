"""
run_projection.py
-----------------
Waikato mainstem, 2025 -> 2100, four SSPs.

For every site: annual-mean pH and Coulombic departure D = PZC(T) - pH,
computed over a monthly temperature cycle, with
  - additive supersaturation (flux balance; this framework), and
  - multiplicative supersaturation (fixed RS; the earlier convention),
plus an exact Shapley attribution of the 2025 -> 2100 change in D into
five channels: atmosphere, physical (k, K0), metabolic (MMRT), carbonate
constants, and PZC(T).

Writes results/*.csv. Figures are made by make_figures.py.
"""
import os

import numpy as np
import pandas as pd

import atmosphere as atm
import unified as u
from waikato_sites import load_sites

os.makedirs("results", exist_ok=True)

SSPS = ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]

# NZ mean air warming by 2090 relative to 1995-2014 (MfE 2024 projections
# page). Water warming taken as equal to air warming (ASSUMED), with ~0.3 C
# already realised by 2025.
DT_2090 = {"SSP1-2.6": 0.8, "SSP2-4.5": 1.9, "SSP3-7.0": 3.0, "SSP5-8.5": 3.7}
DT_2025 = 0.3


def warming(year, ssp):
    """Water warming (C) relative to 2025."""
    slope = (DT_2090[ssp] - DT_2025) / (2090 - 2025)
    if ssp == "SSP1-2.6":
        year = min(year, 2090)   # SSP1-2.6 stabilises
    return slope * (year - 2025)


T_MEAN, T_AMP = 16.0, 5.5          # ASSUMED annual cycle, ~10.5-21.5 C
MONTHS = np.arange(12)
T_CYCLE = T_MEAN + T_AMP * np.cos(2 * np.pi * (MONTHS - 1) / 12)   # Feb warmest
G = 0.49                           # GPP/ER, median of 222 US rivers (Bernhardt et al. 2022 data; IQR 0.27-0.73)
YEAR0 = 2025


def annual(site, year, ssp, multiplicative=False, switches=None,
           g=G, metab=u.METAB, dpzc=u.DPZC_DT, f_met=None):
    """Annual-mean pH, D and pCO2_w for one site/year/SSP."""
    f = site.f_met if f_met is None else f_met
    atm0 = float(atm.co2_atm_ppm(YEAR0, ssp))
    a = float(atm.co2_atm_ppm(year, ssp))
    dT = warming(year, ssp)
    out = []
    for Tm in T_CYCLE:
        if switches:
            # a frozen channel must see this month's own 2025 temperature
            s = _switched(site, Tm, dT, a, atm0, f, g, metab, dpzc, switches)
        else:
            s = u.site_state(site.alk_mgL, Tm + dT, a,
                             dp_geo0=site.dp0 * (1 - f), dp_met0=site.dp0 * f,
                             T0_C=T_MEAN, g=g, metab=metab,
                             multiplicative=multiplicative, S0=site.S0,
                             dpzc_dT=dpzc, pco2_atm0_ppm=atm0)
        out.append((s["pH"], s["D"], s["pco2w"]))
    out = np.array(out)
    return dict(pH=out[:, 0].mean(), D=out[:, 1].mean(), pco2w=out[:, 2].mean(),
                D_max=out[:, 1].max(), pH_summer=out[0, 0])


def _switched(site, Tm, dT, a, atm0, f, g, metab, dpzc, sw):
    """Channel-frozen state: a frozen channel sees this month's 2025 T."""
    T_now, T_ref = Tm + dT, Tm
    Tp = T_now if sw["phys"] else T_ref
    Tmt = T_now if sw["met"] else T_ref
    Tc = T_now if sw["carb"] else T_ref
    Tz = T_now if sw["pzc"] else T_ref
    A = a if sw["atm"] else atm0
    pr = u.phys_rel(Tp, T_MEAN)
    fr = u.net_source_rel(Tmt, T_MEAN, g, metab)
    pco2w = A + site.dp0 * (1 - f) * pr + site.dp0 * f * fr * pr
    pH = u.cs.pH_from_Alk_pCO2(site.alk_mgL, pco2w * 1e-6, Tc)["pH"]
    return dict(pH=pH, D=u.pzc(Tz, dpzc_dT=dpzc) - pH, pco2w=pco2w)


CHANNELS = ["atm", "phys", "met", "carb", "pzc"]


def attribution(site, ssp, year=2100, target="D"):
    def fn(on):
        return annual(site, year, ssp, switches=on)[target]
    return u.shapley(fn, CHANNELS)


if __name__ == "__main__":
    sites = load_sites()

    # --- main projection table ---------------------------------------------
    rows = []
    for _, st in sites.iterrows():
        for ssp in SSPS:
            for mode in ("additive", "multiplicative"):
                m = mode == "multiplicative"
                for yr in (2025, 2050, 2075, 2100):
                    r = annual(st, yr, ssp, multiplicative=m)
                    rows.append(dict(site=st.site, ssp=ssp, mode=mode, year=yr, **r))
    proj = pd.DataFrame(rows)
    proj.to_csv("results/projection.csv", index=False)

    # --- attribution of the 2025 -> 2100 change -----------------------------
    arows = []
    for _, st in sites.iterrows():
        for ssp in SSPS:
            for target in ("D", "pH"):
                a = attribution(st, ssp, target=target)
                arows.append(dict(site=st.site, ssp=ssp, target=target, **a,
                                  total=sum(a.values())))
    attr = pd.DataFrame(arows)
    attr.to_csv("results/attribution.csv", index=False)

    # --- summary printout ----------------------------------------------------
    p = proj.pivot_table(index=["site", "ssp", "mode"], columns="year", values="D")
    p["dD_2100"] = p[2100] - p[2025]
    pd.set_option("display.width", 200)
    print(p.round(3).to_string())
    print()
    print(attr.round(4).to_string(index=False))
