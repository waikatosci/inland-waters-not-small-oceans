"""Forward model of Lake Karapiro, 1850-2100: carbonate chemistry, surface charge and the
availability of metals, molybdate and phosphate, free CO2 for primary producers, and
carbonate saturation for shell-building animals.

Builds on fig8_reconstruction.py (same forcing, site and metabolic settings) and adds:
  - free CO2 (umol/L) and carbonate ion, annual means over the monthly temperature cycle;
  - saturation indices SI = log10(IAP/Ksp) for calcite and aragonite. Ksp from Plummer and
    Busenberg (1982). Ca2+ = 5.33 mg/L (median of 60 filtered samples, Lake Karapiro top,
    Oct 2024 to Feb 2026), held constant. Activity coefficients for Ca2+ and CO3-- from the
    Davies equation at I = 1.5 mmol/L (from the major-ion medians and alkalinity);
  - the availability of each species relative to 1850, as 10^(f z dD): f is the observed
    fraction of the electrostatic limit (divalent cations at the median of published field labile and
    free fractions, phi = 0.04, and of diel cycles with temperature held constant, phi = 0.15, Data S1;
    phosphate 0.013 from GLORICH) or 1 for the limit itself (divalent cations, phosphate, molybdate).
Writes results/forward_karapiro.csv."""
import numpy as np
import pandas as pd
import carbonate_system as cs
import unified as u
import run_projection as rp
from scipy.interpolate import PchipInterpolator
import atmosphere as atm
from waikato_sites import load_sites

# Forcing exactly as in fig8_reconstruction.py
HIST = {1850: 285.2, 1870: 287.6, 1890: 293.9, 1910: 299.7, 1930: 307.2, 1950: 311.3,
        1960: 316.9, 1970: 325.7, 1980: 338.8, 1990: 354.4, 2000: 369.7, 2010: 389.9,
        2015: 401.0, 2020: 414.2, 2025: 427.7}
_hist = PchipInterpolator(list(HIST), list(HIST.values()))
DT_HIST = 1.2
SCEN = ["SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]


def co2(year, ssp):
    if year <= 2025:
        return float(_hist(year))
    return float(atm.co2_atm_ppm(year, ssp) - atm.co2_atm_ppm(2025, ssp) + HIST[2025])


def dT(year, ssp):
    if year <= 2025:
        f = np.log(co2(year, ssp) / HIST[1850]) / np.log(HIST[2025] / HIST[1850])
        return DT_HIST * (f - 1)
    return rp.warming(year, ssp)


site = load_sites().set_index("site").loc["Karāpiro"]

CA_MGL = 5.33
I_MOLL = 1.5e-3
PHI = {"field": 0.04, "diel": 0.15}   # medians of published responses, divalent cations other than Fe (Fig. 3c)
SPECIES = {  # name: (z, f)
    "Divalent field": (2.0, PHI["field"]), "Divalent diel": (2.0, PHI["diel"]), "Divalent limit": (2.0, 1.0),
    "Phosphate": (-1.5, 0.013), "Phosphate limit": (-1.5, 1.0), "Molybdate limit": (-2.0, 1.0)}


def gamma2(I=I_MOLL, T_C=16.0):
    A = 0.4883 + 8.074e-4 * T_C          # Debye-Huckel A, adequate for 0-30 C
    s = np.sqrt(I)
    return 10 ** (-A * 4 * (s / (1 + s) - 0.3 * I))


def log_ksp(T_C, mineral):
    T = T_C + 273.15
    c = 2839.319 if mineral == "calcite" else 2903.293
    a = -171.9065 if mineral == "calcite" else -171.9773
    return a - 0.077993 * T + c / T + 71.595 * np.log10(T)


def si(co3_molL, T_C, mineral, ca_mgL=CA_MGL):
    g = gamma2(T_C=T_C)
    ca = ca_mgL / 40078.0
    return np.log10(g * ca * g * co3_molL) - log_ksp(T_C, mineral)


def month_state(alk_mgL, T_C, pco2_atm_ppm, dp_geo0, dp_met0, metab=u.METAB):
    st = u.site_state(alk_mgL, T_C, pco2_atm_ppm, dp_geo0=dp_geo0, dp_met0=dp_met0,
                      T0_C=rp.T_MEAN, g=rp.G, metab=metab)
    c = cs.pH_from_Alk_pCO2(alk_mgL, st["pco2w"] * 1e-6, T_C)
    return dict(pH=c["pH"], D=st["D"], pco2w=st["pco2w"], co2_umol=c["CO2"] * 1e6,
                si_cal=si(c["CO3"], T_C, "calcite"), si_arag=si(c["CO3"], T_C, "aragonite"))


CASES = {"R>P": u.METAB, "R=P": u.METAB_EQ}   # MMRT central; equal activation energies (0.57 eV)


def annual(year, ssp, case="R>P"):
    a, d = co2(year, ssp), dT(year, ssp)
    rows = [month_state(site.alk_mgL, Tm + d, a, site.dp0 * (1 - site.f_met), site.dp0 * site.f_met,
                        metab=CASES[case]) for Tm in rp.T_CYCLE]
    out = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
    out.update(year=year, ssp=ssp, co2_atm=a, dT=d, case=case)
    return out


if __name__ == "__main__":
    years = sorted(set(range(1850, 2101, 5)) | {1990, 2025})
    parts = []
    for case in CASES:
        df = pd.DataFrame([annual(y, s, case) for s in SCEN for y in years if y >= 2025 or s == SCEN[0]])
        base = df[df.year == 1850].iloc[0]
        for name, (z, f) in SPECIES.items():
            df[name] = 100 * (10 ** (f * z * (df.D - base.D)) - 1)
        df["co2_rel"] = 100 * (df.co2_umol / base.co2_umol - 1)
        parts.append(df)
    df = pd.concat(parts, ignore_index=True)
    df.to_csv("results/forward_karapiro.csv", index=False)
    show = ["pH", "D", "co2_umol", "si_cal", "si_arag"] + list(SPECIES)
    print(f"gamma2 at 16 C: {gamma2():.3f}")
    for case in CASES:
        print("==", case)
        for y, s in [(1850, "SSP2-4.5"), (1990, "SSP2-4.5"), (2025, "SSP2-4.5")] + [(2100, s) for s in SCEN]:
            r = df[(df.case == case) & (df.year == y) & ((df.ssp == s) | (df.year <= 2025))].iloc[0]
            print(y, s if y == 2100 else "", " ".join(f"{k}={r[k]:.3g}" for k in show + ["co2_rel"]))
