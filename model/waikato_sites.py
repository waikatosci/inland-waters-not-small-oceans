"""
waikato_sites.py
----------------
Reference state for the 12 Waikato mainstem sites: WRC field pH and WRC field
temperature on the October 2022 survey day (WR_Acidification_Data_25-10-22.xlsx,
sheet Student_Data_WR_only), with alkalinity from Waikato Regional Council
monitoring (ALKT, potentiometric titration to pH 4.5, 2015-22 median at each
station; data/wrc_alkalinity_2015_2022.csv). Ohaaki Bridge is not a council
station and uses the median of ECP Gran titrations there (2025-26). The survey's
own bench alkalinities ran 5-90% above the council values and are not used.
Field pCO2 is back-calculated from alkalinity, pH and temperature.

Lake Karapiro, the site used for the forward model, attribution and Figure 5, instead takes its reference state from
long-term council monitoring at Karapiro tailrace (1131-79, monthly, 2015-2023): each grab sample is moved from its
sampling time (median 07:45) to the daily mean with the Hamilton sonde diel curves scaled by 0.28 for a tailrace (as in
the reach balance), and S0 is the median of the sample supersaturations against observed atmospheric CO2. A single
survey day is kept for the along-river snapshot (Figure S2/3), where all sites must share one day.

f_met is the share of each site's CO2 excess assumed to come from local
net heterotrophy (the rest is geothermal or otherwise temperature-
insensitive). ASSUMED from the delta13C-DIC mass balance (geothermal
~71% of the budget in the TVZ reach), not measured per site.
"""
import pandas as pd

import carbonate_system as cs

XLSX = "data/WR_Acidification_Data_25-10-22.xlsx"
PCO2_ATM_2022 = 417.0   # ppm, approximate global mean Sep 2022

WRC_ALK = "data/wrc_alkalinity_2015_2022.csv"
WRC_STATION = {"Taupō gates": "Taupo Control Gates", "Ohaaki Br.": "Ohaaki Bridge",
               "Ohakuri": "Ohakuri Tailrace Br", "Whakamaru": "Whakamaru Tailrace",
               "Waipapa": "Waipapa Tailrace", "Karāpiro": "Karapiro Tailrace",
               "Narrows": "Narrows Boat Ramp", "Horotiu": "Horotiu Br", "Huntly": "Huntly-Tainui Br",
               "Rangiriri": "Rangiriri Br", "Mercer": "Mercer Br", "Tuakau": "Tuakau Br"}

# order down the river, short names
SHORT = {
    "Taupo Control Gate": ("Taupō gates", 0.5),
    "Ohaaki Bridge": ("Ohaaki Br.", 0.3),
    "Ohakuri Tailrac": ("Ohakuri", 0.3),
    "Whak": ("Whakamaru", 0.7),
    "Waip": ("Waipapa", 0.7),
    "Kara": ("Karāpiro", 0.7),
    "Narr": ("Narrows", 0.7),
    "Horo": ("Horotiu", 0.7),
    "Hunt": ("Huntly", 0.7),
    "Rangiriri": ("Rangiriri", 0.7),
    "Merc": ("Mercer", 0.7),
    "Tuakau": ("Tuakau", 0.7),
}


DIEL = "data/diel_climatology_hamilton.csv"
CHEM = "data/ECP_ChemData_JL_norm_to_flow.xlsx"   # council monitoring, held by WRC, not deposited
REF = "data/karapiro_reference.csv"                # derived reference state, deposited
OBS_CO2 = {2014: 397.2, 2015: 400.8, 2016: 404.2, 2017: 406.5, 2018: 408.7, 2019: 411.4, 2020: 414.2,
           2021: 416.4, 2022: 418.5, 2023: 421.1}   # NOAA global annual means (ppm)


def karapiro_reference(alk, years=(2015, 2023), scale=0.28):
    """Median supersaturation of council grab samples at Karapiro tailrace, moved to the daily mean."""
    import numpy as np
    from scipy.interpolate import PchipInterpolator
    air = PchipInterpolator(list(OBS_CO2), list(OBS_CO2.values()))
    clim = pd.read_csv(DIEL)
    ph = clim[clim.iloc[:, 0] == "pH"].set_index("m").drop(columns=clim.columns[0])
    d = pd.read_excel(CHEM, "ECP_ChemData_incl_TE xlsx")
    g = d[(d.DeviceID == "1131_79") & d.Date.dt.year.between(*years)].dropna(subset=["_pH", "_Temp"]).copy()
    b = ((g.Date.dt.hour * 60 + g.Date.dt.minute) // 30).astype(int)
    corr = np.array([-float(ph.loc[m].iloc[i]) for m, i in zip(g.Date.dt.month, b)])
    pH = g._pH.values + scale * corr
    yr = (g.Date.dt.year + g.Date.dt.dayofyear / 365.25).clip(upper=max(OBS_CO2)).values
    S = np.array([float(cs.pCO2_from_Alk_pH(alk, p, t)["pCO2_ppm"]) / float(air(y)) for p, t, y in zip(pH, g._Temp, yr)])
    return dict(n=len(g), S0=float(np.median(S)), pH=float(np.median(pH)), T=float(np.median(g._Temp)),
                y0=int(g.Date.dt.year.min()), y1=int(g.Date.dt.year.max()))


def load_sites(path=XLSX, reference="monitoring"):
    raw = pd.read_excel(path, sheet_name="Student_Data_WR_only", header=0)
    raw = raw.dropna(subset=["Sample Name"])
    walk = pd.read_csv(WRC_ALK).set_index("STATION_NAME")["alk_mgL_CaCO3"]
    rows = []
    for _, r in raw.iterrows():
        name = str(r["Sample Name"])
        key = next(k for k in SHORT if k in name)
        short, fmet = SHORT[key]
        alk = float(walk[WRC_STATION[short]])
        pH, T = float(r["WRC pH"]), float(r["WRC Temperature"])
        pco2 = float(cs.pCO2_from_Alk_pH(alk, pH, T)["pCO2_ppm"])
        rows.append(dict(site=short, alk_mgL=alk, alk_survey=float(r["Experimental Bicarbonate"]), pH_field=pH, T_field=T,
                         pco2_field=pco2, S0=pco2 / PCO2_ATM_2022,
                         dp0=pco2 - PCO2_ATM_2022, f_met=fmet))
    df = pd.DataFrame(rows)
    if reference == "monitoring":
        i = df.index[df.site == "Karāpiro"][0]
        import os
        if os.path.exists(CHEM) and os.path.exists(DIEL):
            k = karapiro_reference(df.at[i, "alk_mgL"])
            pd.DataFrame([k]).to_csv(REF, index=False)
        else:   # council records not available: use the deposited reference state
            k = pd.read_csv(REF).iloc[0].to_dict()
        df.at[i, "pH_field"], df.at[i, "T_field"], df.at[i, "S0"] = k["pH"], k["T"], k["S0"]
        df.at[i, "pco2_field"] = k["S0"] * PCO2_ATM_2022
        df.at[i, "dp0"] = (k["S0"] - 1) * PCO2_ATM_2022
    return df


if __name__ == "__main__":
    print(load_sites().round(2).to_string(index=False))
