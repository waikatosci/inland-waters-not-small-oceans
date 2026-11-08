"""Curves, slopes and a quick figure from outputs/nicad_raw.csv (written by 09_run_nicad_coupled.py).

Outputs
  outputs/models_nicad_curves.csv  long format (one row per scenario x pH x element x quantity)
  outputs/models_nicad_slopes.csv  same columns as models_own_slopes.csv:
                                   model, water, sorbent, DOC_mgL, element, quantity, pH_range, slope, notes
  outputs/models_nicad_curves.png  soft water: log10 free and inorganic-labile fractions vs pH, DOC 0/1/3/10 mg/L,
                                   with and without HFO (1 mg Fe/L)

Quantities (all per kg water; = per L within 0.1 %)
  free_ion            free aquo ion in bulk solution (outside the FA Donnan phase), mol/kgw
  inorganic_labile    free ion + dissolved inorganic complexes in bulk solution, mol/kgw
  dissolved_fraction  (inorganic dissolved + held by FA) / total; FA (colloidal) passes a filter, HFO does not
  plus, in the curves only: free_fraction, inorganic_labile_fraction, FA_bound_fraction, sorbed_fraction, dissolved
Slopes: least squares of log10(quantity) on pH over all grid points in the window (0.05 step, 31 points).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model_setup as M  # noqa: E402

OUT = os.path.join(M.ROOT, "outputs")
WINDOWS = [(6.5, 8.0), (7.0, 8.5)]
CATS = M.CATIONS


def sorbent_label(fe):
    return "none" if fe == 0 else f"HFO_{fe:g}mgFe/L"


def long_curves(raw):
    recs = []
    for _, r in raw.iterrows():
        w = r["mass_water_kg"]
        base = dict(model=r["model"], water=r["water"], sorbent=sorbent_label(r["Fe_mgL"]), Fe_mgL=r["Fe_mgL"],
                    DOC_mgL=r["DOC_mgL"], FA_mgL=r["FA_mgL"], FA_per_DOC=r["FA_per_DOC"],
                    Al_control=r["Al_control"], pH=round(r["pH_calc"], 4),
                    pCO2_uatm=10 ** r["logpCO2_atm"] * 1e6, alkalinity_meqL=r["alk_eq_kgw"] * 1e3,
                    DIC_mmolL=r["dic_molal"] * 1e3, SI_calcite=r["SI_Calcite"],
                    FA_logchi=r.get("FA_logchi", np.nan), FA_q_eqkg=r.get("FA_q_eqkg", np.nan))
        for el in CATS:
            tot = r[f"{el}_total_mol"]
            inorg = r[f"{el}_diss_molal"]
            held = r[f"{el}_FA_held_mol"]
            sorb = r[f"{el}_sorbed_mol"]
            q = {
                "free_ion": (r[f"{el}_free_molal"], "mol/kgw"),
                "inorganic_labile": (inorg, "mol/kgw"),
                "dissolved": (inorg + held / w, "mol/kgw"),
                "dissolved_fraction": ((inorg * w + held) / tot, "-"),
                "free_fraction": (r[f"{el}_free_molal"] * w / tot, "-"),
                "inorganic_labile_fraction": (inorg * w / tot, "-"),
                "FA_bound_fraction": (held / tot, "-"),
                "sorbed_fraction": (sorb / tot, "-"),
            }
            for name, (val, unit) in q.items():
                recs.append(dict(base, element=el, element_class="cation", quantity=name, value=val, unit=unit))
    return pd.DataFrame(recs)


def slopes(curves):
    out = []
    sub = curves[curves.quantity.isin(["dissolved_fraction", "free_ion", "inorganic_labile"])]
    fab = curves[curves.quantity == "FA_bound_fraction"].set_index(
        ["model", "water", "sorbent", "DOC_mgL", "element", "pH"])["value"]
    keys = ["model", "water", "sorbent", "Fe_mgL", "DOC_mgL", "FA_mgL", "FA_per_DOC", "Al_control", "element",
            "quantity"]
    for key, g in sub.groupby(keys):
        model, water, sorbent, fe, doc, fa, fpd, al, el, qty = key
        g = g.sort_values("pH")
        for lo, hi in WINDOWS:
            m = (g.pH >= lo - 1e-6) & (g.pH <= hi + 1e-6)
            x = g.pH[m].values
            y = np.log10(g.value[m].values)
            if len(x) < 3 or not np.all(np.isfinite(y)):
                continue
            b, a = np.polyfit(x, y, 1)
            yhat = a + b * x
            ss_res, ss_tot = np.sum((y - yhat) ** 2), np.sum((y - y.mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
            local = np.gradient(y, x)
            pc = g.pCO2_uatm[m].values
            sic = g.SI_calcite[m].values
            notes = [f"Fe={fe:g} mg/L", f"FA={fa:g} mg/L ({fpd:g} x DOC)" if doc > 0 else "no organic matter",
                     f"n={len(x)}", f"r2={r2:.4f}" if np.isfinite(r2) else "r2=NA (flat)",
                     f"local_slope_min={local.min():.3f}", f"local_slope_max={local.max():.3f}",
                     f"log10_value_at_{lo}={y[0]:.3f}", f"log10_value_at_{hi}={y[-1]:.3f}",
                     f"pCO2_uatm={pc.max():.0f}-{pc.min():.0f}", f"SI_calcite_max={sic.max():.2f}"]
            if doc > 0:
                f_lo = fab[(model, water, sorbent, doc, el, round(x[0], 4))]
                f_hi = fab[(model, water, sorbent, doc, el, round(x[-1], 4))]
                notes.append(f"FA_bound_fraction={f_lo:.3f}@{lo}-{f_hi:.3f}@{hi}")
            flags = []
            if pc.max() > 20000:
                flags.append("pCO2 above 20000 uatm needed at low-pH end")
            if sic.max() > 0:
                flags.append("calcite supersaturated at high-pH end (precipitation suppressed)")
            if doc > 0:
                flags.append("generic NICA-Donnan FA (Milne et al. 2003 Table 3) at 25 C constants")
                if al:
                    flags.append("Al3+ competition on, free Al3+ fixed by Al(OH)3 (log K 8.5, dH -107 kJ/mol; Lofts & Tipping 2011)")
                else:
                    flags.append("no Al or Fe(III) competition for FA sites")
                est = {"Zn": "Zn FA n1,n2 estimated (italic) in Milne 2003 Table 3; generic FA binds Zn weakly",
                       "Mn": "Mn FA parameters from one small dataset with n1, logK2, n2 estimated",
                       "Co": "Co FA n2 estimated", "Pb": "Pb FA n1, n2 fixed at HA values",
                       "Cd": "Cd FA n2 estimated"}
                if el in est:
                    flags.append(est[el])
                if qty == "dissolved_fraction":
                    flags.append("dissolved includes FA-bound metal (FA passes filter)")
            if sorbent == "none" and qty == "dissolved_fraction" and doc == 0:
                flags.append("no sorbent: dissolved fraction is 1 by construction")
            if el == "Mn" and sorbent != "none":
                flags.append("Mn(II)-HFO constants: D&M (1990) table 10.5 estimates (wateq4f.dat)")
            if el == "Co" and sorbent != "none":
                flags.append("Co-HFO constants only in minteq.v4.dat")
            out.append(dict(model=model, water=water, sorbent=sorbent, DOC_mgL=doc, element=el, quantity=qty,
                            pH_range=f"{lo}-{hi}", slope=round(b, 4), notes="; ".join(notes + flags)))
    return pd.DataFrame(out)


def plot(curves):
    INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
    DOC_COL = {0.0: "#86b6ef", 1.0: "#3987e5", 3.0: "#1c5cab", 10.0: "#0d366b"}   # ordinal blue, validated
    c = curves[(curves.water == "soft") & (curves.FA_per_DOC == 1.65) & (~curves.Al_control.astype(bool))]
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.linewidth": 0.6})
    rows = [("free_fraction", "free ion / total"), ("inorganic_labile_fraction", "inorganic labile / total")]
    fig, axes = plt.subplots(2, 7, figsize=(17, 6.4), sharex=True)
    for ri, (qty, lab) in enumerate(rows):
        for ci, el in enumerate(CATS):
            ax = axes[ri, ci]
            ax.axvspan(6.5, 8.0, color="#f0efec", zorder=0, lw=0)
            ax.grid(True, color=GRID, lw=0.6, zorder=1)
            for doc, col in DOC_COL.items():
                for fe, ls in ((0.0, "-"), (1.0, "--")):
                    g = c[(c.DOC_mgL == doc) & (c.Fe_mgL == fe) & (c.element == el) & (c.quantity == qty)]
                    g = g.sort_values("pH")
                    ax.plot(g.pH, np.log10(g.value), color=col, ls=ls, lw=1.6, zorder=3)
            if ri == 0:
                ax.set_title(el, loc="left", fontweight="bold", color=INK, fontsize=10)
            if ci == 0:
                ax.set_ylabel(f"log10 {lab}")
            if ri == 1:
                ax.set_xlabel("pH (set by CO2)")
            ax.set_xlim(6.0, 9.0)
            lo, hi = ax.get_ylim()
            if hi - lo < 0.2:
                mid = 0.5 * (lo + hi)
                ax.set_ylim(mid - 0.1, mid + 0.1)
    handles = [Line2D([], [], color=DOC_COL[d], lw=2, label=f"DOC {d:g} mg/L") for d in DOC_COL]
    handles += [Line2D([], [], color=INK2, lw=1.6, ls="-", label="no oxide"),
                Line2D([], [], color=INK2, lw=1.6, ls="--", label="HFO 1 mg Fe/L")]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Soft water (alk 0.78 meq/L, Ca 5.3 mg/L): generic NICA-Donnan fulvic acid (FA = 1.65 x DOC) "
                 "+ inorganic speciation (+ HFO), 16 C. Grey band: pH 6.5-8.0", x=0.01, ha="left", fontsize=10,
                 color=INK)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(os.path.join(OUT, "models_nicad_curves.png"), dpi=150)


def main():
    raw = pd.read_csv(os.path.join(OUT, "nicad_raw.csv"))
    assert (raw.pH_calc - raw.pH_target).abs().max() < 1e-3
    mb = raw[[c for c in raw.columns if c.endswith("_massbal_relerr")]].abs().max().max()
    assert mb < 1e-5, mb
    red = raw[[c for c in raw.columns if c.startswith("redoxcheck_")]].abs().max().max()
    assert red < 1e-20, red
    curves = long_curves(raw)
    curves.to_csv(os.path.join(OUT, "models_nicad_curves.csv"), index=False)
    sl = slopes(curves).sort_values(["model", "quantity", "water", "sorbent", "DOC_mgL", "element", "pH_range"])
    sl.to_csv(os.path.join(OUT, "models_nicad_slopes.csv"), index=False)
    plot(curves)
    print("curves", curves.shape, "slopes", sl.shape, "max mass-balance error", mb)


if __name__ == "__main__":
    main()
