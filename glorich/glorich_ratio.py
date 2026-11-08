import os
"""Dilution-free test: within-site rho(pH, DIP/TP). Dilution by high flow lowers DIP and TP together,
so their ratio removes it; the electrostatic effect predicts a higher dissolved fraction at higher pH
in soft water (positive rho). Also partial rho holding log discharge constant."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata, wilcoxon
COLS = ["STAT_ID","Discharge_inst","pH","Alkalinity","DIP","TP","DOSAT","DO_mgL","Temp_water"]
hc = pd.read_csv(os.environ.get("GLORICH_CSV", "../external/glorich/hydrochemistry.csv"), usecols=COLS, encoding="latin-1", dtype=str)
for c in COLS[1:]:
    hc[c] = pd.to_numeric(hc[c], errors="coerce")
s = hc.dropna(subset=["STAT_ID","pH","Alkalinity","DIP","TP"])
s = s[(s.pH > 3) & (s.pH < 12) & (s.DIP > 0) & (s.TP > 0) & (s.Alkalinity > 0) & (s.DIP <= 1.05*s.TP)].copy()
s["frac"] = s.DIP / s.TP
def pr(h, ctrl, a, b):
    Z = np.column_stack([np.ones(len(h))] + [rankdata(h[c]) for c in ctrl])
    f = lambda y: y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    return np.corrcoef(f(rankdata(h[a])), f(rankdata(h[b])))[0, 1]
rows = []
for sid, g in s.groupby("STAT_ID"):
    if len(g) < 30 or g.frac.nunique() < 5: continue
    rec = dict(STAT_ID=sid, med_Alk=g.Alkalinity.median(), n=len(g),
               rho_DIP=spearmanr(g.pH, g.DIP)[0], rho_TP=spearmanr(g.pH, g.TP)[0],
               rho_frac=spearmanr(g.pH, g.frac)[0])
    h = g[g.Discharge_inst > 0].assign(logQ=lambda x: np.log(x.Discharge_inst))
    if len(h) >= 30:
        rec["rho_frac_sameQ"] = spearmanr(h.pH, h.frac)[0]
        rec["rho_frac_Q"] = pr(h, ["logQ"], "pH", "frac")
    rows.append(rec)
d = pd.DataFrame(rows)
d["dec"] = pd.qcut(d.med_Alk, 10, labels=False)
out = d.groupby("dec").agg(logAlk=("med_Alk", lambda x: np.log10(x.median())), n=("rho_frac","size"),
    rho_DIP=("rho_DIP","median"), rho_TP=("rho_TP","median"), rho_frac=("rho_frac","median"),
    nQ=("rho_frac_Q","count"), frac_sameQ=("rho_frac_sameQ","median"), frac_Q=("rho_frac_Q","median"))
pd.set_option("display.width", 200); print(out.round(3).to_string())
soft = d[d.med_Alk < 1000]
print(f"soft n={len(soft)}: rho_DIP {soft.rho_DIP.median():+.3f}, rho_TP {soft.rho_TP.median():+.3f}, "
      f"rho_frac {soft.rho_frac.median():+.3f} (p={wilcoxon(soft.rho_frac.dropna()).pvalue:.2g})")
x = soft.dropna(subset=["rho_frac_Q"])
print(f"soft with Q n={len(x)}: rho_frac {x.rho_frac_sameQ.median():+.3f} -> |Q {x.rho_frac_Q.median():+.3f} (p={wilcoxon(x.rho_frac_Q).pvalue:.2g})")
d.to_csv("glorich_ratio_by_site.csv", index=False)
