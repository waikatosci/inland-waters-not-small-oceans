import os
"""Mixed-effects test of the alkalinity dependence of the pH-phosphate slope in GLORICH, with and
without discharge and season. Within-site centred predictors (so the pH slope is a within-site
slope), random intercept and random pH slope per station. Sites with >= 30 paired pH-DIP
records; up to 60 records per site drawn at random. Soft = site-median alkalinity < 1 meq/L.
Models:
  A  logDIP ~ pHc * soft                                  (as in the earlier submission)
  B  logDIP ~ pHc * soft + logQc + sin + cos              (discharge and season held constant)
  C  as B, records below each site's median discharge only (base flow)
  D  as B, with continuous log10 site alkalinity in place of soft
Writes glorich_mixed_results.csv."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import statsmodels.formula.api as smf

COLS = ["STAT_ID", "RESULT_DATETIME", "Discharge_inst", "pH", "Alkalinity", "DIP"]
hc = pd.read_csv(os.environ.get("GLORICH_CSV", "../external/glorich/hydrochemistry.csv"), usecols=COLS, encoding="latin-1", dtype=str)
for c in COLS[2:]:
    hc[c] = pd.to_numeric(hc[c], errors="coerce")
s = hc.dropna(subset=["STAT_ID", "pH", "Alkalinity", "DIP"])
s = s[(s.pH > 3) & (s.pH < 12) & (s.DIP > 0) & (s.Alkalinity > 0)].copy()
s["date"] = pd.to_datetime(s.RESULT_DATETIME, errors="coerce", dayfirst=True)
n = s.groupby("STAT_ID").size()
s = s[s.STAT_ID.isin(n[n >= 30].index)]
alk = s.groupby("STAT_ID").Alkalinity.median()          # ueq/L
s["soft"] = (s.STAT_ID.map(alk) < 1000).astype(int)
s["logAlk"] = np.log10(s.STAT_ID.map(alk))
s["logDIP"] = np.log10(s.DIP)
s["pHc"] = s.pH - s.groupby("STAT_ID").pH.transform("mean")
s["logQ"] = np.log10(s.Discharge_inst.where(s.Discharge_inst > 0))
s["logQc"] = s.logQ - s.groupby("STAT_ID").logQ.transform("mean")
doy = s.date.dt.dayofyear
s["sin"], s["cos"] = np.sin(2 * np.pi * doy / 365.25), np.cos(2 * np.pi * doy / 365.25)
s["base"] = s.logQ < s.groupby("STAT_ID").logQ.transform("median")
rng = np.random.default_rng(1)
samp = s.loc[s.sample(frac=1, random_state=1).groupby("STAT_ID").head(60).index]
print(f"{samp.STAT_ID.nunique()} sites, {len(samp)} records (all); soft sites {samp[samp.soft==1].STAT_ID.nunique()}")


def fit(name, formula, d, term):
    d = d.dropna(subset=[c for c in ["logDIP", "pHc", "logQc", "sin"] if c in formula])
    m = smf.mixedlm(formula, d, groups=d["STAT_ID"], re_formula="~pHc").fit(reml=True, method="lbfgs")
    b, se, p = m.fe_params[term], m.bse[term], m.pvalues[term]
    print(f"{name}: n={len(d)} sites={d.STAT_ID.nunique()} {term} = {b:+.4f} (SE {se:.4f}), p = {p:.1e}; "
          f"pHc (hard/base) = {m.fe_params['pHc']:+.4f}")
    if term == "pHc:logAlk":
        import json
        cv = m.cov_params().loc[["pHc", term], ["pHc", term]].values.tolist()
        json.dump(dict(b0=m.fe_params["pHc"], b1=b, cov=cv,
                       alk_range=[float(d.logAlk.quantile(.02)), float(d.logAlk.quantile(.98))]),
                  open("glorich_mixed_D.json", "w"))
    return dict(model=name, n=len(d), sites=d.STAT_ID.nunique(), term=term, beta=b, se=se, p=p,
                pH_slope_ref=m.fe_params["pHc"])


rows = [
    fit("A no controls", "logDIP ~ pHc * soft", samp, "pHc:soft"),
    fit("B discharge + season", "logDIP ~ pHc * soft + logQc + sin + cos", samp, "pHc:soft"),
    fit("C base flow", "logDIP ~ pHc * soft + logQc + sin + cos", samp[samp.base], "pHc:soft"),
    fit("D continuous alkalinity", "logDIP ~ pHc * logAlk + logQc + sin + cos", samp, "pHc:logAlk"),
]
pd.DataFrame(rows).to_csv("glorich_mixed_results.csv", index=False)
