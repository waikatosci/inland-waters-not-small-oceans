import os
"""Global upscaling of the internal (metabolic) feedback on river acidification, GLORICH rivers.

Per sample: pCO2 of the water from pH, alkalinity and temperature (carbonate_system, Plummer & Busenberg 1982),
S = pCO2_w / atmospheric CO2 in the sample year. Per station (>= 10 valid samples): median S, mean temperature.
Projection 2025 -> 2100 (flux balance, Text S1): the water excess p_w - p_a is held by gas exchange, and its
respiratory part scales with temperature as the net metabolic source (literature MMRT values: E_R 0.62, E_P 0.53 eV,
GPP/ER 0.49); the physical terms (solubility x transfer velocity) nearly cancel and are ignored here.
  atmospheric term:   dp_atm = p_a(2100) - p_a(2025)
  metabolic term:     dp_met = f * (S - 1) * p_a(2025) * (F(T+dT)/F(T) - 1)
f = share of the CO2 excess from respiration that responds to warming: in-water metabolism only (Hotchkiss et al.
2015: about 14% in small streams to 40% in large rivers, 30% on average), or all respiration including soil CO2
delivered by groundwater (f = 1, upper bound). dT: global warming 2025 -> 2100 of 1.4 C (SSP2-4.5) and 3.4 C
(SSP5-8.5), water warming taken equal to air warming (land warms more, so conservative).
Filters: pH 6-9.5, alkalinity >= 200 ueq/L (pCO2 from pH and alkalinity is biased high in soft, acidic,
organic-rich water; Abril et al. 2015), 5-35 C.
Writes station table and summary."""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, "../model")
import carbonate_system as cs, unified as u
from scipy.interpolate import PchipInterpolator

HIST = {1950: 311.3, 1960: 316.9, 1970: 325.7, 1980: 338.8, 1990: 354.4, 2000: 369.7, 2010: 389.9, 2015: 401.0,
        2020: 414.2, 2025: 427.7}
hist = PchipInterpolator(list(HIST), list(HIST.values()))
C = ["STAT_ID", "RESULT_DATETIME", "Discharge_inst", "Temp_water", "pH", "Alkalinity"]
d = pd.read_csv(os.environ.get("GLORICH_CSV", "../external/glorich/hydrochemistry.csv"), usecols=C, dtype=str, encoding="latin-1")
for c in C[2:]:
    d[c] = pd.to_numeric(d[c], errors="coerce")
d["year"] = pd.to_datetime(d.RESULT_DATETIME, errors="coerce", dayfirst=True).dt.year
d = d.dropna(subset=["pH", "Alkalinity", "Temp_water", "year"])
d = d[(d.pH.between(6, 9.5)) & (d.Alkalinity >= 200) & (d.Temp_water.between(5, 35)) & (d.year >= 1950)]
alk_mg = d.Alkalinity.values * 50.04 / 1000
d["pco2"] = np.asarray(cs.pCO2_from_Alk_pH(alk_mg, d.pH.values, d.Temp_water.values)["pCO2_ppm"], float)
d["pa"] = hist(d.year.clip(upper=2025))
d["S"] = d.pco2 / d.pa
n = d.groupby("STAT_ID").size()
st = d[d.STAT_ID.isin(n[n >= 10].index)].groupby("STAT_ID").agg(n=("S", "size"), S=("S", "median"),
        frac_under=("S", lambda s: (s < 1).mean()), T=("Temp_water", "mean"), alk=("Alkalinity", "median"),
        Q=("Discharge_inst", "median"))
PA25 = 427.7
SC = {"SSP2-4.5": (603.0, 1.4), "SSP5-8.5": (1135.0, 3.4)}
for ssp, (pa100, dT) in SC.items():
    r = np.array([float(u.net_source_rel(t + dT, t, 0.49)) for t in st["T"]])
    for lab, f in [("inwater", 0.30), ("allresp", 1.0)]:
        dpm = f * (st.S - 1).clip(lower=0) * PA25 * (r - 1)
        dpa = pa100 - PA25
        st[f"R_{lab}_{ssp}"] = dpm / dpa
        pw0 = st.S * PA25
        st[f"rel_{lab}_{ssp}"] = np.log10((pw0 + dpa + dpm) / pw0) / np.log10(pa100 / PA25)   # 1 = ocean-like
st.to_csv("glorich_stations.csv")
out = dict(stations=len(st), samples=int(st.n.sum()), median_S=st.S.median(), S_q25=st.S.quantile(.25),
           S_q75=st.S.quantile(.75), frac_S_gt1=(st.S > 1).mean(), frac_S_gt3=(st.S > 3).mean(),
           median_frac_samples_under=st.frac_under.median(), frac_stations_some_under=(st.frac_under > 0.1).mean())
for ssp in SC:
    for lab in ["inwater", "allresp"]:
        R = st[f"R_{lab}_{ssp}"]; rel = st[f"rel_{lab}_{ssp}"]
        out[f"{ssp} {lab}: median metab/atm"] = R.median()
        out[f"{ssp} {lab}: frac metab>=atm"] = (R >= 1).mean()
        out[f"{ssp} {lab}: frac metab>=0.4 atm"] = (R >= 0.4).mean()
        out[f"{ssp} {lab}: median pH response vs ocean-like"] = rel.median()
soft = st[st.alk < 1000]
out["soft stations"] = len(soft); out["soft median S"] = soft.S.median()
for ssp in SC:
    out[f"soft {ssp} inwater median R"] = soft[f"R_inwater_{ssp}"].median()
    out[f"soft {ssp} allresp median R"] = soft[f"R_allresp_{ssp}"].median()
pd.Series(out).to_csv("glorich_summary.csv")
print(pd.Series(out).round(3).to_string())
