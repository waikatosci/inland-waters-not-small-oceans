import os
"""Flow-stratified test: within-site rho(pH, DIP/TP) and rho(pH, DIP) using only samples below each
site's median discharge (base flow), where storm dilution and particulate P are minimised."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon
COLS = ["STAT_ID","Discharge_inst","pH","Alkalinity","DIP","TP"]
hc = pd.read_csv(os.environ.get("GLORICH_CSV", "../external/glorich/hydrochemistry.csv"), usecols=COLS, encoding="latin-1", dtype=str)
for c in COLS[1:]:
    hc[c] = pd.to_numeric(hc[c], errors="coerce")
s = hc.dropna(subset=["STAT_ID","pH","Alkalinity","DIP","Discharge_inst"])
s = s[(s.pH > 3) & (s.pH < 12) & (s.DIP > 0) & (s.Alkalinity > 0) & (s.Discharge_inst > 0)].copy()
rows = []
for sid, g in s.groupby("STAT_ID"):
    if len(g) < 60: continue
    for part, h in [("low", g[g.Discharge_inst <= g.Discharge_inst.median()]),
                    ("high", g[g.Discharge_inst > g.Discharge_inst.median()])]:
        rec = dict(STAT_ID=sid, part=part, med_Alk=g.Alkalinity.median(), n=len(h),
                   rho_DIP=spearmanr(h.pH, h.DIP)[0], rho_pH_Q=spearmanr(h.pH, h.Discharge_inst)[0])
        k = h.dropna(subset=["TP"]); k = k[(k.TP > 0) & (k.DIP <= 1.05*k.TP)]
        if len(k) >= 25: rec["rho_frac"] = spearmanr(k.pH, k.DIP/k.TP)[0]
        rows.append(rec)
d = pd.DataFrame(rows)
d["soft"] = np.where(d.med_Alk < 1000, "soft <1 meq/L", np.where(d.med_Alk < 2000, "1-2 meq/L", "hard >2 meq/L"))
pd.set_option("display.width", 200)
out = d.groupby(["soft","part"]).agg(sites=("STAT_ID","nunique"), rho_DIP=("rho_DIP","median"),
                                      n_frac=("rho_frac","count"), rho_frac=("rho_frac","median"),
                                      rho_pH_Q=("rho_pH_Q","median"))
print(out.round(3).to_string())
for part in "low", "high":
    x = d[(d.soft == "soft <1 meq/L") & (d.part == part)]
    print(part, "flow, soft: DIP p=", f"{wilcoxon(x.rho_DIP.dropna()).pvalue:.2g}",
          "frac p=", f"{wilcoxon(x.rho_frac.dropna()).pvalue:.2g}")
d.to_csv("glorich_flowsplit_by_site.csv", index=False)
