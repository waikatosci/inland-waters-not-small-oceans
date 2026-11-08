import os
"""Within-site pH-DIP coupling in GLORICH, raw and with photosynthesis held constant.

For each station with >= 30 paired pH-DIP samples: Spearman rho(pH, DIP), and the partial
rank correlation holding oxygen saturation (photosynthesis proxy), water temperature and
log discharge constant. Also the within-site slope d log10(DIP)/d pH (OLS on residuals),
for comparison with the electrostatic strong-sorption limit of +1.5 (z = -1.5).
Sites are grouped into deciles of median alkalinity, as in Fig. 3b of the Coulombic baseline."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata

COLS = ["STAT_ID","RESULT_DATETIME","Discharge_inst","Temp_water","pH","DO_mgL","DOSAT","Alkalinity","DOC","DIP"]
hc = pd.read_csv(os.environ.get("GLORICH_CSV", "../external/glorich/hydrochemistry.csv"), usecols=COLS, encoding="latin-1", dtype=str)
for c in COLS[2:]:
    hc[c] = pd.to_numeric(hc[c], errors="coerce")
sub = hc.dropna(subset=["STAT_ID","pH","Alkalinity","DIP"])
sub = sub[(sub.pH > 3) & (sub.pH < 12) & (sub.DIP > 0) & (sub.Alkalinity > 0)].copy()
T = sub.Temp_water
o2sat = 14.652 - 0.41022*T + 0.007991*T**2 - 0.000077774*T**3      # mg/L, fresh water, 1 atm
sub["DOpct"] = sub.DOSAT.where(sub.DOSAT.notna(), 100*sub.DO_mgL/o2sat)
sub.loc[(sub.DOpct <= 0) | (sub.DOpct > 400), "DOpct"] = np.nan


def resid(y, Z):
    return y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]


rows = []
for sid, g in sub.groupby("STAT_ID"):
    if len(g) < 30:
        continue
    r_raw = spearmanr(g.pH, g.DIP)[0]
    rec = dict(STAT_ID=sid, n=len(g), med_Alk=g.Alkalinity.median(), med_DOC=g.DOC.median(),
               med_pH=g.pH.median(), rho=r_raw)
    h = g.dropna(subset=["DOpct", "Temp_water", "Discharge_inst"])
    h = h[h.Discharge_inst > 0]
    if len(h) >= 30 and h.pH.std() > 0:
        rec["n_ctrl"] = len(h)
        rec["rho_same"] = spearmanr(h.pH, h.DIP)[0]          # raw rho on the same samples
        for name, ctrl in [("rho_DO", ["DOpct"]), ("rho_DOTQ", ["DOpct", "Temp_water", "logQ"])]:
            hh = h.assign(logQ=np.log(h.Discharge_inst))
            Z = np.column_stack([np.ones(len(hh))] + [rankdata(hh[c]) for c in ctrl])
            ex, ey = resid(rankdata(hh.pH), Z), resid(rankdata(hh.DIP), Z)
            rec[name] = np.corrcoef(ex, ey)[0, 1]
        Z = np.column_stack([np.ones(len(hh)), hh.DOpct, hh.Temp_water, hh.logQ])
        rec["slope_raw"] = np.polyfit(h.pH, np.log10(h.DIP), 1)[0]
        ex, ey = resid(hh.pH.values, Z), resid(np.log10(hh.DIP.values), Z)
        rec["slope_ctrl"] = (ex @ ey) / (ex @ ex)
        rec["rho_pH_DO"] = spearmanr(h.pH, h.DOpct)[0]
    rows.append(rec)
site = pd.DataFrame(rows)
site["logAlk"] = np.log10(site.med_Alk.clip(lower=10))
site["alk_d"] = pd.qcut(site.med_Alk, 10, labels=False)
site.to_csv("glorich_site_rho.csv", index=False)

pd.set_option("display.width", 220)
print(f"stations with >=30 pairs: {len(site)}; with O2, T and Q controls: {site.n_ctrl.notna().sum()}")
agg = site.groupby("alk_d").agg(logAlk=("logAlk", "median"), n=("rho", "size"), rho=("rho", "median"),
                                n_ctrl=("rho_DOTQ", "count"), rho_same=("rho_same", "median"),
                                rho_DO=("rho_DO", "median"), rho_DOTQ=("rho_DOTQ", "median"),
                                slope_raw=("slope_raw", "median"), slope_ctrl=("slope_ctrl", "median"),
                                rho_pH_DO=("rho_pH_DO", "median"))
print(agg.round(3).to_string())
c = site.dropna(subset=["rho_DOTQ"])
soft = c[c.med_Alk < 1000]
print(f"\nsoft (<1 meq/L) sites with controls: {len(soft)}")
print(f"  median rho same samples {soft.rho_same.median():+.3f}; | DO {soft.rho_DO.median():+.3f}; | DO,T,Q {soft.rho_DOTQ.median():+.3f}")
print(f"  share of sites where control makes rho more positive: {(soft.rho_DOTQ > soft.rho_same).mean():.0%}")
print(f"  slope d log DIP/dpH raw {soft.slope_raw.median():+.3f}, controlled {soft.slope_ctrl.median():+.3f} (IQR {soft.slope_ctrl.quantile(.25):+.2f} to {soft.slope_ctrl.quantile(.75):+.2f})")
from scipy.stats import wilcoxon
print("  Wilcoxon controlled rho vs 0:", wilcoxon(soft.rho_DOTQ).pvalue, "; controlled vs raw:", wilcoxon(soft.rho_DOTQ - soft.rho_same).pvalue)
