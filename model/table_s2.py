"""Table S2: one-at-a-time sensitivity of the 2025->2100 change in annual-mean D at Lake Karapiro
(full model, monthly cycle), with phosphate K_d ratio at the electrostatic limit (10^(1.5*dD)).
Central: literature metabolism (unified.METAB), g = 0.49, f = 0.7, dPZC/dT = -0.02. Writes results/table_s2.csv."""
import pandas as pd
import unified as u, run_projection as rp
from waikato_sites import load_sites
site = load_sites().set_index("site").loc["Karāpiro"].copy()
cases = [("Central (GPP/ER = 0.49, E_R = 0.62 eV, E_P = 0.53 eV, f = 0.7, dPZC/dT = −0.02 per °C)", {}),
         ("GPP/ER = 0.27 (river lower quartile)", dict(g=0.27)),
         ("GPP/ER = 0.76 (lake median)", dict(g=0.76)),
         ("GPP/ER = 0.9", dict(g=0.9)),
         ("Optimum of photosynthesis 24 °C", dict(metab=dict(u.METAB, Topt_P=24.0))),
         ("Optimum of photosynthesis 30 °C", dict(metab=dict(u.METAB, Topt_P=30.0))),
         ("E_R = 0.5 eV", dict(metab=dict(u.METAB, Ea_R=0.5))),
         ("E_R = 0.8 eV", dict(metab=dict(u.METAB, Ea_R=0.8))),
         ("E_P = E_R = 0.57 eV", dict(metab=u.METAB_EQ)),
         ("Organism-level theory, E_R = 0.65 eV, E_P = 0.32 eV", dict(metab=u.METAB_THEORY)),
         ("f = 0.3", dict(f_met=0.3)),
         ("f = 1.0", dict(f_met=1.0)),
         ("dPZC/dT = +0.01 per °C", dict(dpzc=0.01)),
         ("dPZC/dT = 0", dict(dpzc=0.0)),
         ("dPZC/dT = −0.007 per °C", dict(dpzc=-0.007)),
         ("dPZC/dT = −0.03 per °C", dict(dpzc=-0.03))]
rows=[]
for name,kw in cases:
    r=dict(case=name)
    for ssp in ("SSP2-4.5","SSP5-8.5"):
        dD=rp.annual(site,2100,ssp,**kw)["D"]-rp.annual(site,2025,ssp,**kw)["D"]
        r[f"dD {ssp}"]=round(dD,3); r[f"Kd ratio {ssp}"]=round(10**(1.5*dD),2)
    rows.append(r)
t=pd.DataFrame(rows); t.to_csv("results/table_s2.csv",index=False); print(t.to_string(index=False))
