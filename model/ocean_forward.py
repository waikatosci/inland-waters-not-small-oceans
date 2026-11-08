"""Surface ocean under the same forcing as the Lake Karapiro forward model, for comparison (Fig. 5).

Global-mean surface seawater: salinity 35, total alkalinity 2300 umol/kg (held constant), in equilibrium
with the atmosphere (pCO2 = atmospheric CO2, same record and SSP trajectories as forward_model.py).
Sea surface temperature 18 C in 2025, warming at three-quarters of the water warming used for the lake,
which reproduces the CMIP6 multimodel mean change under SSP5-8.5 (+3.47 C, 2080-2099 vs 1870-1899;
Kwiatkowski et al. 2020). Carbonate chemistry from PyCO2SYS (Humphreys et al. 2022) with the constants
of Lueker et al. (2000); pH on the total scale. Free CO2 converted to umol/L with a density of 1.025 kg/L.
D = PZC(T) - pH with the same ferrihydrite PZC(T) as the lake, so the electrostatic limit can be drawn
for the ocean on the same basis (for comparison only: specific adsorption of seawater ions is ignored).
Writes results/forward_ocean.csv."""
import numpy as np
import pandas as pd
import PyCO2SYS as pyco2
import unified as u
import forward_model as fm

TA, SAL, SST2025, RHO, F_SST = 2300.0, 35.0, 18.0, 1.025, 0.75
SCEN = ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]


def state(year, ssp):
    a = fm.co2(year, ssp)
    T = SST2025 + F_SST * fm.dT(year, ssp)
    r = pyco2.sys(par1=TA, par2=a, par1_type=1, par2_type=4, salinity=SAL, temperature=T,
                  opt_k_carbonic=10, opt_pH_scale=1)
    pH = float(r["pH_total"])
    return dict(year=year, ssp=ssp, co2_atm=a, sst=T, pH=pH, co2_umol=float(r["CO2"]) * RHO,
                si_cal=float(np.log10(r["saturation_calcite"])), si_arag=float(np.log10(r["saturation_aragonite"])),
                D=float(u.pzc(T)) - pH)


if __name__ == "__main__":
    years = sorted(set(range(1850, 2101, 5)) | {1990})
    df = pd.DataFrame([state(y, s) for s in SCEN for y in years if y >= 2025 or s == "SSP2-4.5"])
    base = df[df.year == 1850].iloc[0]
    for name, z in (("Divalent limit", 2.0), ("Phosphate limit", -1.5), ("Molybdate limit", -2.0)):
        df[name] = 100 * (10 ** (z * (df.D - base.D)) - 1)
    df["co2_rel"] = 100 * (df.co2_umol / base.co2_umol - 1)
    df.to_csv("results/forward_ocean.csv", index=False)
    # check against CMIP6 (Kwiatkowski et al. 2020): 2080-2099 minus 1870-1899
    h = df[(df.year >= 1870) & (df.year <= 1899)]
    for s in SCEN:
        f = df[(df.ssp == s) & (df.year >= 2080) & (df.year <= 2099)]
        print(s, "dpH %.3f dSST %.2f" % (f.pH.mean() - h.pH.mean(), f.sst.mean() - h.sst.mean()))
    for y, s in [(1850, "SSP2-4.5"), (1990, "SSP2-4.5"), (2025, "SSP2-4.5")] + [(2100, s) for s in SCEN]:
        r = df[(df.year == y) & ((df.ssp == s) | (df.year <= 2025))].iloc[0]
        print(y, s if y == 2100 else "", " ".join(f"{k}={r[k]:.3g}" for k in ["co2_atm", "sst", "pH", "co2_umol", "si_cal", "si_arag", "D", "Divalent limit", "Phosphate limit", "co2_rel"]))
    h1 = df[(df.year >= 1990) & (df.year <= 2025) & (df.ssp == "SSP2-4.5")]
    print("1990-2025 pH slope per yr", np.polyfit(h1.year, h1.pH, 1)[0])
