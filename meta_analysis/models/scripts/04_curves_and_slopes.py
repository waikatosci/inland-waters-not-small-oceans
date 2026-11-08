"""Turn the raw PHREEQC output into the long-format curves file and the slopes file.

Inputs : outputs/hfo_raw_<model>.csv written by 03_run_hfo_model.py
Outputs: outputs/models_own_curves.csv   (long format: one row per scenario x pH x element x quantity)
         outputs/models_own_slopes.csv   (columns: model, water, sorbent, DOC_mgL, element, quantity,
                                          pH_range, slope, notes)
         outputs/models_own_conditions.csv (one row per scenario x pH: pCO2, alkalinity, SI, surface potential)

Quantities
  dissolved_fraction        dissolved / system total  (dissolved = not sorbed to HFO)
  free_ion                  free-ion concentration, mol/kg water (= mol/L to <0.1 % here). For oxyanions:
                            sum of the protonation forms of the monomeric oxyanion (no metal complexes)
  inorganic_labile          free ion + all dissolved inorganic complexes, mol/kg water. With no organic
                            matter this equals the dissolved concentration.
Slopes are least-squares slopes of log10(quantity) on pH using every grid point (0.05 pH step) inside
the window, end points included. Because the system total is fixed, the slope of log10(dissolved_fraction)
equals the slope of log10(dissolved concentration).
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model_setup as M  # noqa: E402

OUT = os.path.join(M.ROOT, "outputs")
MODELS = {
    "DLM_HFO_minteqv4+CO3": "hfo_raw_DLM_HFO_minteqv4plusCO3.csv",
    "DLM_HFO_minteqv4_noCO3": "hfo_raw_DLM_HFO_minteqv4_noCO3.csv",
}
WINDOWS = [(6.5, 8.0), (7.0, 8.5)]
TOTALS = {el: ug * 1e-6 / M.GFW[el] for el, (name, ug) in M.TRACES.items()}


def sorbent_label(fe):
    return "none" if fe == 0 else f"HFO_{fe:g}mgFe/L"


def long_curves(raw):
    recs = []
    for _, r in raw.iterrows():
        w = r["mass_water_kg"]
        base = dict(model=r["model"], engine=r["engine"], database=r["database"], water=r["water"],
                    sorbent=sorbent_label(r["Fe_mgL"]), Fe_mgL=r["Fe_mgL"], DOC_mgL=r["DOC_mgL"],
                    pH=round(r["pH_calc"], 4), pCO2_uatm=10 ** r["logpCO2_atm"] * 1e6,
                    alkalinity_meqL=r["alk_eq_kgw"] * 1e3, DIC_mmolL=r["dic_molal"] * 1e3,
                    SI_calcite=r["SI_Calcite"])
        for el in M.TRACES:
            diss_mol = r[f"{el}_diss_molal"] * w
            sorb = r[f"{el}_sorbed_mol"]
            total = diss_mol + sorb
            q = {
                "dissolved_fraction": (diss_mol / total, "-"),
                "dissolved": (r[f"{el}_diss_molal"], "mol/kgw"),
                "free_ion": (r[f"{el}_free_molal"], "mol/kgw"),
                "free_fraction": (r[f"{el}_free_molal"] * w / total, "-"),
                "inorganic_labile": (r[f"{el}_diss_molal"], "mol/kgw"),
                "inorganic_labile_fraction": (diss_mol / total, "-"),
                "sorbed_fraction": (sorb / total, "-"),
            }
            if el in M.CATIONS:
                q["free_ion_activity"] = (r[f"{el}_free_act"], "mol/kgw")
            for name, (val, unit) in q.items():
                recs.append(dict(base, element=el, element_class="cation" if el in M.CATIONS else "oxyanion",
                                 quantity=name, value=val, unit=unit))
    return pd.DataFrame(recs)


def conditions(raw):
    keep = ["model", "engine", "water", "Fe_mgL", "pH_target", "pH_calc", "pe_calc", "logpCO2_atm", "alk_eq_kgw",
            "dic_molal", "mu", "init_alk_eq_kgw", "init_Na_molal", "init_Cl_molal", "psi_V", "sigma_eq",
            "Ca_sorbed_mol", "Mg_sorbed_mol", "C_sorbed_mol", "S_sorbed_mol"] + \
           [c for c in raw.columns if c.startswith("SI_") or c.startswith("redoxcheck_")]
    d = raw[[c for c in keep if c in raw.columns]].copy()
    d["pCO2_uatm"] = 10 ** d["logpCO2_atm"] * 1e6
    return d


def slopes(curves, cond):
    out = []
    qmap = {"dissolved_fraction": "dissolved_fraction", "free_ion": "free_ion", "inorganic_labile": "inorganic_labile"}
    sub = curves[curves.quantity.isin(qmap)]
    for key, g in sub.groupby(["model", "water", "sorbent", "Fe_mgL", "DOC_mgL", "element", "quantity"]):
        model, water, sorbent, fe, doc, el, qty = key
        g = g.sort_values("pH")
        for lo, hi in WINDOWS:
            m = (g.pH >= lo - 1e-6) & (g.pH <= hi + 1e-6)
            x = g.pH[m].values
            y = np.log10(g.value[m].values)
            if len(x) < 3 or not np.all(np.isfinite(y)):
                continue
            b, a = np.polyfit(x, y, 1)
            yhat = a + b * x
            ss_res = np.sum((y - yhat) ** 2)
            ss_tot = np.sum((y - y.mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
            local = np.gradient(y, x)
            pc = g.pCO2_uatm[m].values
            sic = g.SI_calcite[m].values
            notes = [f"Fe={fe:g} mg/L", f"n={len(x)}", f"r2={r2:.4f}" if np.isfinite(r2) else "r2=NA (flat)",
                     f"local_slope_min={local.min():.3f}", f"local_slope_max={local.max():.3f}",
                     f"log10_value_at_{lo}={y[0]:.3f}", f"log10_value_at_{hi}={y[-1]:.3f}",
                     f"pCO2_uatm={pc.max():.0f}-{pc.min():.0f}", f"SI_calcite_max={sic.max():.2f}"]
            flags = []
            if pc.max() > 20000:
                flags.append("pCO2 above 20000 uatm needed at low-pH end")
            if sic.max() > 0:
                flags.append("calcite supersaturated at high-pH end (precipitation suppressed)")
            if qty == "inorganic_labile":
                flags.append("no organic matter: inorganic labile = dissolved")
            if el == "V" and sorbent != "none":
                flags.append("V(V)-HFO represented by one species (Hfo_OHVO4-3) in minteq.v4.dat; source not stated; treat with caution")
            if el in ("Mo", "V", "Co") and sorbent != "none":
                flags.append(f"{el}-HFO constants only in minteq.v4.dat (MINTEQA2 v4 DLM set)")
            if el == "Mn" and sorbent != "none":
                flags.append("Mn(II)-HFO constants: D&M (1990) table 10.5 linear free-energy estimates, taken from wateq4f.dat")
            if model.endswith("+CO3") and el in M.OXYANIONS and sorbent != "none":
                flags.append("carbonate competes for HFO sites (Appelo et al. 2002 constants); compare noCO3 variant")
            if sorbent == "none" and qty != "free_ion":
                flags.append("no sorbent: dissolved fraction is 1 by construction")
            out.append(dict(model=model, water=water, sorbent=sorbent, DOC_mgL=doc, element=el, quantity=qty,
                            pH_range=f"{lo}-{hi}", slope=round(b, 4), notes="; ".join(notes + flags)))
    return pd.DataFrame(out)


def main():
    raws = []
    for model, fname in MODELS.items():
        p = os.path.join(OUT, fname)
        raws.append(pd.read_csv(p))
    raw = pd.concat(raws, ignore_index=True)
    # sanity: target pH reached, mass balance closed, redox fixed
    assert (raw.pH_calc - raw.pH_target).abs().max() < 1e-3
    for el in M.TRACES:
        tot = raw[f"{el}_diss_molal"] * raw["mass_water_kg"] + raw[f"{el}_sorbed_mol"]
        err = (tot / TOTALS[el] - 1).abs().max()
        assert err < 2e-3, (el, err)
    red = raw[[c for c in raw.columns if c.startswith("redoxcheck_")]].abs().max().max()
    assert red < 1e-20, red
    curves = long_curves(raw)
    curves.to_csv(os.path.join(OUT, "models_own_curves.csv"), index=False)
    cond = conditions(raw)
    cond.to_csv(os.path.join(OUT, "models_own_conditions.csv"), index=False)
    sl = slopes(curves, cond)
    sl = sl.sort_values(["model", "quantity", "water", "sorbent", "element", "pH_range"])
    sl.to_csv(os.path.join(OUT, "models_own_slopes.csv"), index=False)
    print("curves", curves.shape, "slopes", sl.shape, "max redox-check total", red)


if __name__ == "__main__":
    main()
