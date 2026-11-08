import os
"""Which control removes the soft-water pH-DIP correlation? Partial rho holding each of O2
saturation, temperature and log discharge constant on its own, by alkalinity decile."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata, wilcoxon
COLS = ["STAT_ID","Discharge_inst","Temp_water","pH","DO_mgL","DOSAT","Alkalinity","DIP"]
hc = pd.read_csv(os.environ.get("GLORICH_CSV", "../external/glorich/hydrochemistry.csv"), usecols=COLS, encoding="latin-1", dtype=str)
for c in COLS[1:]:
    hc[c] = pd.to_numeric(hc[c], errors="coerce")
s = hc.dropna(subset=["STAT_ID","pH","Alkalinity","DIP"])
s = s[(s.pH > 3) & (s.pH < 12) & (s.DIP > 0) & (s.Alkalinity > 0)].copy()
T = s.Temp_water
s["DOpct"] = s.DOSAT.where(s.DOSAT.notna(), 100*s.DO_mgL/(14.652-0.41022*T+0.007991*T**2-0.000077774*T**3))
s.loc[(s.DOpct <= 0) | (s.DOpct > 400), "DOpct"] = np.nan
def pr(h, ctrl):
    Z = np.column_stack([np.ones(len(h))] + [rankdata(h[c]) for c in ctrl])
    f = lambda y: y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    return np.corrcoef(f(rankdata(h.pH)), f(rankdata(h.DIP)))[0, 1]
rows = []
for sid, g in s.groupby("STAT_ID"):
    if len(g) < 30: continue
    rec = dict(med_Alk=g.Alkalinity.median())
    for name, ctrl in [("T", ["Temp_water"]), ("Q", ["logQ"]), ("DO", ["DOpct"])]:
        h = g.dropna(subset=[c if c != "logQ" else "Discharge_inst" for c in ctrl])
        if "logQ" in ctrl:
            h = h[h.Discharge_inst > 0].assign(logQ=lambda x: np.log(x.Discharge_inst))
        if len(h) >= 30 and h.pH.nunique() > 2 and h.DIP.nunique() > 2:
            rec["raw_" + name] = spearmanr(h.pH, h.DIP)[0]
            rec["ctl_" + name] = pr(h, ctrl)
            if name == "T": rec["rho_pH_T"] = spearmanr(h.pH, h.Temp_water)[0]; rec["rho_DIP_T"] = spearmanr(h.DIP, h.Temp_water)[0]
    rows.append(rec)
d = pd.DataFrame(rows)
d["dec"] = pd.qcut(d.med_Alk, 10, labels=False)
out = d.groupby("dec").agg(logAlk=("med_Alk", lambda x: np.log10(x.median())),
    nT=("ctl_T","count"), raw_T=("raw_T","median"), ctl_T=("ctl_T","median"),
    nQ=("ctl_Q","count"), raw_Q=("raw_Q","median"), ctl_Q=("ctl_Q","median"),
    nDO=("ctl_DO","count"), raw_DO=("raw_DO","median"), ctl_DO=("ctl_DO","median"),
    pH_T=("rho_pH_T","median"), DIP_T=("rho_DIP_T","median"))
pd.set_option("display.width", 220); print(out.round(3).to_string())
soft = d[d.med_Alk < 1000]
for k in "T", "Q", "DO":
    x = soft.dropna(subset=["ctl_" + k])
    print(f"soft, control {k}: n={len(x)} raw {x['raw_'+k].median():+.3f} -> {x['ctl_'+k].median():+.3f}; "
          f"Wilcoxon raw>0 p={wilcoxon(x['raw_'+k]).pvalue:.2g}, ctl p={wilcoxon(x['ctl_'+k]).pvalue:.2g}")
d.to_csv("glorich_controls_by_site.csv", index=False)
