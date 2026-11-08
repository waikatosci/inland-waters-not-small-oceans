"""Inorganic speciation + Dzombak & Morel (1990) HFO surface complexation across pH, for three waters.

pH is set by changing the CO2 partial pressure at fixed alkalinity: one initial water is defined
(at pH 7.0, alkalinity fixed, charge balanced), and for each target pH PHREEQC adds or removes CO2(g)
until a(H+) = 10^-pH (pseudo-phase Fix_H+ with CO2(g) as the reactant). CO2 exchange does not change
alkalinity; the only other proton sink is the HFO surface itself (reported in the output).

Variants (column 'model'):
  DLM_HFO_minteqv4+CO3   main: minteq.v4.dat + Mn(II) and carbonate surface constants from wateq4f.dat
  DLM_HFO_minteqv4_noCO3 same without carbonate adsorption (closer to Dzombak & Morel's original set)
Sorbent: none, HFO 0.1 mg Fe/L, HFO 1 mg Fe/L. All trace elements are present together (they compete).
Mineral precipitation is not allowed (no EQUILIBRIUM_PHASES other than Fix_H+); saturation indices
are reported so supersaturation can be flagged.

Writes outputs/hfo_raw_<model>.csv (one row per water x sorbent x pH, wide) and
outputs/models_own_conditions.csv.
Usage: python3 scripts/03_run_hfo_model.py [--engine 386|373] [--step 0.05]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engine import Engine, db_path  # noqa: E402
import model_setup as M  # noqa: E402


def mmol(mg_L, key):
    return mg_L / M.GFW[key]


def solution_block(water_key, ref_pH=7.0):
    w = M.WATERS[water_key]
    lines = ["SOLUTION 1", f"   -temp {M.TEMP_C}", "   -units mmol/L", f"   pH {ref_pH}", "   pe 8",
             f"   Alkalinity {w['alk_meq']}    # meq/L (molar units -> equivalents)"]
    for el in ("Ca", "Mg", "Na", "K", "Cl"):
        val = mmol(w[el], el)
        tag = "    charge" if w["balance"] == el else ""
        lines.append(f"   {el} {val:.6e}{tag}")
    lines.append(f"   S(6) {mmol(w['SO4'], 'SO4'):.6e}")
    for el, (name, ug) in M.TRACES.items():
        lines.append(f"   {name} {ug * 1e-3 / M.GFW[el]:.6e}")
    lines.append("END")
    return "\n".join(lines) + "\n"


def surface_block(fe_mg_L):
    a = M.hfo_amounts(fe_mg_L)
    return ("SURFACE 1\n"
            f"   Hfo_wOH {a['weak']:.6e} {M.HFO_AREA_M2_G} {a['g_hfo']:.6e}\n"
            f"   Hfo_sOH {a['strong']:.6e}\n"
            "END\n")


def punch_block(has_surface):
    heads, lines = [], []
    n = [10]

    def add(head, expr):
        heads.append(head)
        lines.append(f"{n[0]} PUNCH {expr}")
        n[0] += 10

    add("pH_calc", '-LA("H+")')
    add("pe_calc", '-LA("e-")')
    add("mu", "MU")
    add("alk_eq_kgw", "ALK")
    add("dic_molal", 'TOT("C(4)")')
    add("logpCO2_atm", 'SI("CO2(g)")')
    add("mass_water_kg", 'TOT("water")')
    add("charge_balance_eq", "CHARGE_BALANCE")
    add("Na_molal", 'TOT("Na")')
    add("Cl_molal", 'TOT("Cl")')
    add("Ca_molal", 'TOT("Ca")')
    for ph in M.SI_PHASES:
        add(f"SI_{ph}", f'SI("{ph}")')
    for chk in ("Mn(3)", "Co(3)", "Cu(1)", "As(3)", "V(4)", "Se(6)", "Se(-2)", "S(-2)"):
        add(f"redoxcheck_{chk}", f'TOT("{chk}")')
    if has_surface:
        add("psi_V", 'EDL("psi", "Hfo")')
        add("sigma_eq", 'EDL("charge", "Hfo")')
        for e in ("Ca", "Mg", "C", "S"):
            add(f"{e}_sorbed_mol", f'SURF("{e}", "Hfo")')
    for el in M.TRACES:
        add(f"{el}_diss_molal", f'TOT("{el}")')
        add(f"{el}_sorbed_mol", f'SURF("{el}", "Hfo")' if has_surface else "0")
        expr = " + ".join(f'MOL("{s}")' for s in M.FREE_SPECIES[el])
        add(f"{el}_free_molal", expr)
        if el in M.CATIONS:
            add(f"{el}_free_act", f'ACT("{M.FREE_SPECIES[el][0]}")')
    blk = "SELECTED_OUTPUT\n   -reset false\n   -state true\nUSER_PUNCH\n   -headings " + " ".join(heads) + "\n   -start\n"
    blk += "\n".join("   " + l for l in lines) + "\n   -end\nEND\n"
    return blk


def build_input(model, water_key, fe_mg_L, ph_grid):
    s = ""
    s += M.REDOX_FIX_MINTEQV4
    s += M.SURF_MN
    if model == "DLM_HFO_minteqv4+CO3":
        s += M.SURF_CO3
    s += M.PHASES_FIX
    has_surface = fe_mg_L > 0
    s += punch_block(has_surface)
    s += solution_block(water_key)
    if has_surface:
        s += surface_block(fe_mg_L)
    for ph in ph_grid:
        s += "USE solution 1\n"
        if has_surface:
            s += "USE surface 1\n"
        s += f"EQUILIBRIUM_PHASES 1\n   Fix_H+ {-ph:.4f} CO2(g) 10.0\nEND\n"
    return s


def system_totals():
    """System totals in mol (per kg water of the initial solution, ~ per L)."""
    return {el: ug * 1e-6 / M.GFW[el] for el, (name, ug) in M.TRACES.items()}


def run(model, engine_which="386", step=0.05):
    ph_grid = np.round(np.arange(6.0, 9.0 + 1e-9, step), 4)
    out = []
    for water in M.WATERS:
        for fe in (0.0, 0.1, 1.0):
            eng = Engine(db_path("minteq.v4.dat"), engine_which)   # fresh IPhreeqc instance per scenario
            txt = build_input(model, water, fe, ph_grid)
            rows = eng.run(txt)
            init = [r for r in rows if str(r["state"]).strip() == "i_soln"]
            react = [r for r in rows if str(r["state"]).strip() == "react"]
            assert len(react) == len(ph_grid), (water, fe, len(react))
            for target, r in zip(ph_grid, react):
                r = dict(r)
                r.update(model=model, engine=eng.version, database="minteq.v4.dat", water=water,
                         sorbent="none" if fe == 0 else "HFO", Fe_mgL=fe, DOC_mgL=0.0, pH_target=target,
                         init_alk_eq_kgw=init[0]["alk_eq_kgw"], init_Na_molal=init[0]["Na_molal"],
                         init_Cl_molal=init[0]["Cl_molal"])
                out.append(r)
    df = pd.DataFrame(out)
    df = df.drop(columns=[c for c in df.columns if c == "state"])
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", default="386")
    ap.add_argument("--step", type=float, default=0.05)
    ap.add_argument("--models", default="DLM_HFO_minteqv4+CO3,DLM_HFO_minteqv4_noCO3")
    a = ap.parse_args()
    outdir = os.path.join(M.ROOT, "outputs")
    os.makedirs(outdir, exist_ok=True)
    for model in a.models.split(","):
        df = run(model, a.engine, a.step)
        tag = model.replace("+", "plus")
        suffix = "" if a.engine == "386" else f"_engine{a.engine}"
        # cross-check runs with the second engine go to validation/, not outputs/
        dest = outdir if a.engine == "386" else os.path.join(M.ROOT, "validation")
        path = os.path.join(dest, f"hfo_raw_{tag}{suffix}.csv")
        df.to_csv(path, index=False)
        print("wrote", path, df.shape)


if __name__ == "__main__":
    main()
