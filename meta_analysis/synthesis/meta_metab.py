"""Literature synthesis of whole-ecosystem temperature sensitivity (apparent activation energy, eV) of GPP and
ER in fresh waters, and of GPP/ER (g). Freshwater whole systems (lakes, rivers, streams, outdoor mesocosms);
designs that follow temperature within systems (seasonal, diel, warming experiments, geothermal gradients).
Cross-site space-for-time, organism-level, theory and biofilm chamber estimates are kept out of the main
estimate and reported as sensitivity cases. Writes metab_summary.csv and metab_rows.csv."""
import numpy as np, pandas as pd

X = "../extractions/"
Z = 1.959964
d1 = pd.read_csv(X + "metab_ea.csv")
d1 = d1[~((d1.study_id == "Demars2016") & (d1.process == "GPP"))]                       # replaced by full-text row
d1 = d1[~((d1.study_id == "Demars2011KMAE") & (d1.process == "ER") & (d1.Ea_eV == 0.68))]  # replaced by Perkins 2012
d = pd.concat([d1, pd.read_csv(X + "metab2_ea.csv")], ignore_index=True)
d["se2"] = d.se
m = d.se2.isna() & d.ci_low.notna()
d.loc[m, "se2"] = (d.ci_high - d.ci_low)[m] / (2 * Z)
FW = ["lake", "river", "stream", "mesocosm"]
WITHIN = ["within_site_seasonal", "within_site_diel", "experimental_warming", "geothermal_gradient"]
main = d[d.system_type.isin(FW) & d.design.isin(WITHIN) & d.process.isin(["GPP", "ER"])].copy()
main.to_csv("metab_rows.csv", index=False)


def dl(y, se):
    y, v = np.asarray(y, float), np.asarray(se, float) ** 2
    w = 1 / v; mf = (w * y).sum() / w.sum(); Q = (w * (y - mf) ** 2).sum(); k = len(y)
    t2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0
    ws = 1 / (v + t2); mu = (ws * y).sum() / ws.sum()
    return mu, np.sqrt(1 / ws.sum()), np.sqrt(t2), k


rows = []
for proc in ["GPP", "ER"]:
    for name, sub in [("main: freshwater, within-system", main),
                      ("plus biofilm channels", d[d.system_type.isin(FW + ["biofilm"]) & d.design.isin(WITHIN)]),
                      ("cross-site, freshwater", d[d.system_type.isin(FW) & (d.design == "cross_site_space_for_time")])]:
        g = sub[sub.process == proc]
        if g.empty:
            continue
        # one value per study (median within study) for the robust summary
        per = g.groupby("study_id").Ea_eV.median()
        r = dict(process=proc, set=name, rows=len(g), studies=per.size, median_of_studies=per.median(),
                 q25=per.quantile(.25), q75=per.quantile(.75), min=per.min(), max=per.max())
        h = g[g.se2.notna()]
        if len(h) >= 2:
            mu, s, tau, k = dl(h.Ea_eV, h.se2)
            r.update(pooled=mu, lo=mu - Z * s, hi=mu + Z * s, tau=tau, k_pooled=k)
        rows.append(r)

# paired gap ER - GPP within the same study and design
gap = []
for (s, des), g in main.groupby(["study_id", "design"]):
    if {"GPP", "ER"} <= set(g.process):
        gap.append(dict(study=s, design=des, ER=g[g.process == "ER"].Ea_eV.median(),
                        GPP=g[g.process == "GPP"].Ea_eV.median()))
G = pd.DataFrame(gap); G["gap"] = G.ER - G.GPP
rows.append(dict(process="ER minus GPP", set="paired within study", rows=len(G), studies=G.study.nunique(),
                 median_of_studies=G.gap.median(), q25=G.gap.quantile(.25), q75=G.gap.quantile(.75),
                 min=G.gap.min(), max=G.gap.max()))

# GPP/ER
q = pd.read_csv(X + "metab_gpp_er.csv")
print(q[["study_id", "system_type", "statistic", "value", "spread"]].to_string())
S = pd.DataFrame(rows)
S.to_csv("metab_summary.csv", index=False)
pd.set_option("display.width", 250)
print(S.round(3).to_string(index=False))
print(G.round(2).to_string(index=False))


# net-source sensitivity at the literature values
def e_net(ER, gap, g):
    return ER + g / (1 - g) * gap


ER_med = S[(S.process == "ER") & S.set.str.startswith("main")].median_of_studies.iloc[0]
for g in (0.27, 0.49, 0.73, 0.9):
    print(f"g={g}: E_net = {e_net(ER_med, G.gap.median(), g):.2f} eV")
