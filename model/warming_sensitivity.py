"""Sensitivity of the Lake Karapiro projections (2025-2100) to the amount of water warming.
Air warming: each of the six CCAM-downscaled CMIP6 models (MfE/NIWA core public dataset, annual mean change
2080-2099 vs 1995-2014), averaged over the 5 km cells covering Lake Karapiro (data/cmip6_karapiro_T_change.csv).
Water warming = r x air warming, r = 0.6-1.0 (most streams warm 0.6-0.8 C per C of air warming; Morrill et al.
2005). The central case of the paper uses the national ensemble mean (MfE 2024) with r = 1.
Writes results/warming_sensitivity.csv."""
import pandas as pd
import unified as u
import run_projection as rp
from waikato_sites import load_sites

site = load_sites().set_index("site").loc["Karāpiro"].copy(); site["site"] = "Karāpiro"
cm = pd.read_csv("data/cmip6_karapiro_T_change.csv")
SSP = {"ssp245": "SSP2-4.5", "ssp370": "SSP3-7.0", "ssp585": "SSP5-8.5"}
_w = rp.warming


def run(ssp, dt2090, r):
    slope = r * (dt2090 - rp.DT_2025) / (2090 - 2025)
    rp.warming = lambda y, s: slope * (y - 2025)
    a = u.shapley(lambda on: rp.annual(site, 2100, ssp, switches=on)["pH"], rp.CHANNELS)
    rp.warming = _w
    return dict(dpH=sum(a.values()), atm=a["atm"], met=a["met"], met_over_atm=a["met"] / a["atm"],
                dT2100=slope * 75)


rows = []
for key, ssp in SSP.items():
    rows.append(dict(ssp=ssp, model="central (NZ ensemble mean)", r=1.0, air2090=rp.DT_2090[ssp],
                     **run(ssp, rp.DT_2090[ssp], 1.0)))
    for _, m in cm[cm.ssp == key].iterrows():
        for r in (0.6, 1.0):
            rows.append(dict(ssp=ssp, model=m.model, r=r, air2090=m.karapiro, **run(ssp, m.karapiro, r)))
out = pd.DataFrame(rows)
out.to_csv("results/warming_sensitivity.csv", index=False)
gcm = out[~out.model.isin(["MMM"]) & ~out.model.str.startswith("central")]
print(out[out.model.str.startswith("central")].round(3).to_string(index=False))
s = gcm.groupby("ssp").agg(air_min=("air2090", "min"), air_max=("air2090", "max"), dpH_min=("dpH", "min"),
                           dpH_max=("dpH", "max"), share_min=("met_over_atm", "min"),
                           share_max=("met_over_atm", "max"))
print(s.round(3).to_string())
s1 = gcm[gcm.r == 1.0].groupby("ssp").agg(dpH_min=("dpH", "min"), dpH_max=("dpH", "max"),
                                          share_min=("met_over_atm", "min"), share_max=("met_over_atm", "max"))
print("r = 1 only (model spread alone)"); print(s1.round(3).to_string())
