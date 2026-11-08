"""Inorganic speciation + HFO (PHREEQC, as in 03_run_hfo_model.py, main model DLM_HFO_minteqv4+CO3) coupled to
generic NICA-Donnan fulvic acid (nicadonnan.py, Milne et al. 2001/2003) for the three waters of model_setup.py.

Coupling (iterative, self-consistent; all equilibria at 16 C, NICA-Donnan constants used at their 25 C values):
  * PHREEQC handles everything inorganic and the HFO surface; NICA-Donnan handles the fulvic acid (FA) "unit"
    (specifically bound ions + counter-ion excess in the Donnan volume). Both see the same bulk free ions.
  * For every pH point the system totals T_c of Ca, Mg, Na, K and the trace elements are fixed. PHREEQC is given
    bulk totals T_c - H_c, where H_c is the amount held by FA (mol/L). After each PHREEQC run we compute, for every
    component, R_c = (bulk inorganic dissolved + HFO-sorbed) / free and, from NICA-Donnan at the PHREEQC free-ion
    concentrations, K_c = H_c / free. The update is H_c <- T_c K_c / (R_c + K_c). Iterated to |change| < 1e-7
    (relative) for every component and pH.
  * pH is set by CO2 at fixed total alkalinity, as in 03_run_hfo_model.py. With FA present the conserved quantity
    is carbonate + organic alkalinity: FA is taken to enter the water as its salt in equilibrium at pH 7.0 (Na as
    counter-ion), so at pH 7.0 the bulk carbonate alkalinity equals the specified water alkalinity. At other pH
    the bulk carbonate alkalinity is Alk_w - (U(pH) - U(7)), U = cation charge held by FA (specific + Donnan),
    i.e. protons released or taken up by FA are exchanged with the carbonate system. Na_T = Na_w + U(7).
  * FA concentration = f x DOC; main f = 1.65 (Sjostedt et al. 2010 ES&T 44:8587, p. 8589, ADOM/DOC calibrated
    for Swedish surface waters, "All of the active DOM was considered to be fulvic acid"); variant f = 1.27
    (Lofts & Tipping 2011 Environ Chem 8:501, p. 505, F_FA,DOC used with WHAM/Model VII).
  * Ions bound specifically by FA: H, Ca, Mg, Mn, Co, Ni, Cu, Zn, Cd, Pb. Donnan-only: Na, K, Cl, SO4, HCO3, CO3,
    OH. No Al or Fe(III) in these waters (as in the oxide-only runs), so their competition is absent.
  * FA does not sorb to HFO and there are no ternary complexes.

Writes outputs/nicad_raw.csv (one row per scenario x pH, wide).
Usage: python3 scripts/09_run_nicad_coupled.py
"""
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engine import Engine, db_path  # noqa: E402
import model_setup as M  # noqa: E402
import nicadonnan as ND  # noqa: E402

PH_GRID = np.round(np.arange(6.0, 9.0 + 1e-9, 0.05), 4)
I7 = int(np.argmin(np.abs(PH_GRID - 7.0)))
MAJORS = ["Ca", "Mg", "Na", "K"]
COMPS = MAJORS + list(M.TRACES)                 # all components whose bulk totals PHREEQC receives
FA_BINDERS = ["Ca", "Mg"] + M.CATIONS           # bound specifically (plus H)
ZC = {"Ca": 2, "Mg": 2, "Na": 1, "K": 1, "Mn": 2, "Ni": 2, "Co": 2, "Zn": 2, "Cd": 2, "Cu": 2, "Pb": 2}
FREE_NAME = {"Ca": "Ca+2", "Mg": "Mg+2", "Na": "Na+", "K": "K+", "Mn": "Mn+2", "Ni": "Ni+2", "Co": "Co+2",
             "Zn": "Zn+2", "Cd": "Cd+2", "Cu": "Cu+2", "Pb": "Pb+2"}
PHREEQC_NAME = {"Ca": "Ca", "Mg": "Mg", "Na": "Na", "K": "K"}
PHREEQC_NAME.update({el: name for el, (name, ug) in M.TRACES.items()})
MODEL = "DLM_HFO_minteqv4+CO3"
# Variant: free Al3+ fixed by aluminium hydroxide solubility as in Lofts & Tipping 2011 (p. 503: "standard
# solubility product of 8.5 and enthalpy of -107 kJ mol-1 for aluminium(III) hydroxide"). The solid is an
# unlimited Al reservoir; Al bound by FA does not change bulk alkalinity (Al(OH)3 + 3 H+ -> Al3+, which then
# displaces 3 H+ from FA), so Al is excluded from U.
AL_PHASE = """
PHASES
AlOH3_LT2011
     Al(OH)3 + 3 H+ = Al+3 + 3 H2O
     log_k 8.5
     delta_h -107 kJ
END
"""


def initial_totals(eng, water):
    """System totals (mol/kgw) and alkalinity of the water as defined in 03_run_hfo_model.py (charge balanced)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import importlib
    run03 = importlib.import_module("03_run_hfo_model")
    txt = M.REDOX_FIX_MINTEQV4 + M.SURF_MN + M.SURF_CO3 + M.PHASES_FIX
    heads = ["alk", "Cl", "S6", "water"] + COMPS
    lines = ['PUNCH ALK, TOT("Cl"), TOT("S(6)"), TOT("water")']
    lines += [f'PUNCH TOT("{PHREEQC_NAME[c].split("(")[0]}")' for c in COMPS]
    txt += ("SELECTED_OUTPUT\n   -reset false\n   -high_precision true\nUSER_PUNCH\n   -headings " + " ".join(heads) + "\n   -start\n" +
            "\n".join(f"   {10 * (k + 1)} {l}" for k, l in enumerate(lines)) + "\n   -end\nEND\n")
    txt += run03.solution_block(water)
    r = eng.run(txt)[-1]
    return r


def solution_text(k, w, tot, bulk, alk):
    lines = [f"SOLUTION {k}", f"   -temp {M.TEMP_C}", "   -units mol/kgw", "   pH 7.0", "   pe 8",
             f"   Alkalinity {alk:.10e}", f"   Cl {tot['Cl']:.10e}", f"   S(6) {tot['S6']:.10e}"]
    for c in COMPS:
        lines.append(f"   {PHREEQC_NAME[c]} {max(bulk[c], 0.0):.10e}")
    return "\n".join(lines) + "\nEND\n"


def punch_text(has_surface, al=False):
    heads, lines = [], []

    def add(h, e):
        heads.append(h)
        lines.append(e)
    add("pH_calc", '-LA("H+")')
    add("mu", "MU")
    add("alk_eq_kgw", "ALK")
    add("dic_molal", 'TOT("C(4)")')
    add("logpCO2_atm", 'SI("CO2(g)")')
    add("mass_water_kg", 'TOT("water")')
    add("cb_eq", "CHARGE_BALANCE")
    add("SI_Calcite", 'SI("Calcite")')
    for sp, h in [("H+", "m_H"), ("OH-", "m_OH"), ("Cl-", "m_Cl"), ("SO4-2", "m_SO4"), ("HCO3-", "m_HCO3"),
                  ("CO3-2", "m_CO3")]:
        add(h, f'MOL("{sp}")')
    for c in COMPS:
        el = PHREEQC_NAME[c].split("(")[0]
        add(f"{c}_diss_molal", f'TOT("{el}")')
        add(f"{c}_sorbed_mol", f'SURF("{el}", "Hfo")' if has_surface else "0")
        if c in FREE_NAME:
            add(f"{c}_free_molal", f'MOL("{FREE_NAME[c]}")')
    if al:
        add("Al_free_molal", 'MOL("Al+3")')
        add("Al_diss_molal", 'TOT("Al")')
    for chk in ("Mn(3)", "Co(3)", "Cu(1)"):
        add(f"redoxcheck_{chk}", f'TOT("{chk}")')
    if has_surface:
        add("psi_V", 'EDL("psi", "Hfo")')
    blk = ("SELECTED_OUTPUT\n   -reset false\n   -high_precision true\n   -state true\nUSER_PUNCH\n   -headings " + " ".join(heads) +
           "\n   -start\n" + "\n".join(f"   {10 * (k + 1)} PUNCH {l}" for k, l in enumerate(lines)) +
           "\n   -end\nEND\n")
    return blk


def run_phreeqc(eng, w, tot, bulk_list, alk_list, fe, al=False):
    has_surface = fe > 0
    txt = M.REDOX_FIX_MINTEQV4 + M.SURF_MN + M.SURF_CO3 + M.PHASES_FIX + (AL_PHASE if al else "")
    txt += punch_text(has_surface, al)
    for k, (bulk, alk) in enumerate(zip(bulk_list, alk_list)):
        txt += solution_text(k + 1, w, tot, bulk, alk)
    if has_surface:
        a = M.hfo_amounts(fe)
        txt += (f"SURFACE 1\n   Hfo_wOH {a['weak']:.6e} {M.HFO_AREA_M2_G} {a['g_hfo']:.6e}\n"
                f"   Hfo_sOH {a['strong']:.6e}\nEND\n")
    for k, ph in enumerate(PH_GRID):
        txt += f"USE solution {k + 1}\n"
        if has_surface:
            txt += "USE surface 1\n"
        txt += f"EQUILIBRIUM_PHASES 1\n   Fix_H+ {-ph:.4f} CO2(g) 10.0\n"
        txt += "   AlOH3_LT2011 0 1e-3\nEND\n" if al else "END\n"
    rows = eng.run(txt)
    react = [r for r in rows if str(r["state"]).strip() == "react"]
    assert len(react) == len(PH_GRID), len(react)
    return react


def nicad_at(r, mFA, al=False):
    """NICA-Donnan for one PHREEQC result row; returns held (mol/L per component), result dict."""
    bulk = {"H": r["m_H"]}
    for c in FA_BINDERS:
        bulk[c] = r[f"{c}_free_molal"]
    if al:
        bulk["Al"] = r["Al_free_molal"]
    extra = {"Na": (r["Na_free_molal"], 1), "K": (r["K_free_molal"], 1), "Cl": (r["m_Cl"], -1),
             "SO4": (r["m_SO4"], -2), "HCO3": (r["m_HCO3"], -1), "CO3": (r["m_CO3"], -2), "OH": (r["m_OH"], -1)}
    res = ND.solve_donnan(bulk, r["mu"], "FA", extra_ions=extra)
    held = ND.held_per_litre(res, mFA)
    return held, res


def solve_held(c, Tc, Rc, free0, res, mFA):
    """Mass balance for one FA-binding component at fixed Boltzmann factor and fixed other ions:
    free * Rc + mFA * (Q_c(free) + V_D * free * (chi^z - 1)) = Tc. Returns the amount held by FA (mol/L)."""
    chi, VD = res["chi"], res["VD"]
    z = ND.CHARGE[c]
    cD = dict(res["cD"])

    def g(logc):
        f = 10 ** logc
        cD[c] = f * chi ** z
        Q = ND.nica_bound(cD, "FA")[c][0]
        return np.log10(f * Rc + mFA * (Q + VD * f * (chi ** z - 1.0))) - np.log10(Tc)
    hi = np.log10(Tc / Rc) + 1e-9
    lo = hi - 12.0
    while g(lo) > 0:
        lo -= 6.0
    logf = brentq(g, lo, hi, xtol=1e-13, rtol=1e-13, maxiter=300)
    return Tc - 10 ** logf * Rc


def run_scenario(eng, water, fe, doc, fa_factor, tol=1e-7, maxit=200, al=False):
    tot = initial_totals(eng, water)
    T = {c: tot[c] for c in COMPS}
    mFA = fa_factor * doc * 1e-6                       # kg FA per L
    n = len(PH_GRID)
    held = [dict.fromkeys(COMPS, 0.0) for _ in range(n)]
    n_extra = 0.0
    it = 0
    while True:
        it += 1
        U = [sum(ZC.get(c, 0) * h[c] for c in COMPS if c in ZC) for h in held]
        bulk_list, alk_list = [], []
        for k in range(n):
            b = {c: T[c] - held[k][c] for c in COMPS}
            b["Na"] = T["Na"] + n_extra - held[k]["Na"]
            bulk_list.append(b)
            alk_list.append(tot["alk"] + n_extra - U[k])
        react = run_phreeqc(eng, water, tot, bulk_list, alk_list, fe, al)
        if mFA == 0:
            nres = [None] * n
            break
        new_held, nres, dmax = [], [], 0.0
        for k, r in enumerate(react):
            h_nd, res = nicad_at(r, mFA, al)
            nres.append(res)
            wkg = r["mass_water_kg"]
            hk = {}
            for c in COMPS:
                if c not in ZC:          # oxyanions: not bound by FA
                    hk[c] = 0.0
                    continue
                free = r[f"{c}_free_molal"]
                Rc = (r[f"{c}_diss_molal"] * wkg + r[f"{c}_sorbed_mol"]) / free
                Tc = T[c] + (n_extra if c == "Na" else 0.0)
                if c in FA_BINDERS:
                    hk[c] = solve_held(c, Tc, Rc, free, res, mFA)
                else:          # Na, K: Donnan only, linear in the free concentration
                    Kc = h_nd.get(c, 0.0) / free
                    hk[c] = Tc * Kc / (Rc + Kc)
                ref = max(abs(hk[c]), 1e-30)
                dmax = max(dmax, abs(hk[c] - held[k][c]) / max(ref, 1e-6 * Tc))
            new_held.append(hk)
        held = new_held
        U7 = sum(ZC[c] * held[I7][c] for c in ZC)
        dmax = max(dmax, abs(U7 - n_extra) / max(U7, 1e-12))
        n_extra = U7
        if os.environ.get('NICAD_DEBUG'):
            print('   it', it, 'dmax', dmax, flush=True)
        if dmax < tol and it > 2:
            # one more PHREEQC pass with the converged held amounts so outputs are consistent
            U = [sum(ZC[c] * h[c] for c in ZC) for h in held]
            bulk_list = []
            alk_list = []
            for k in range(n):
                b = {c: T[c] - held[k][c] for c in COMPS}
                b["Na"] = T["Na"] + n_extra - held[k]["Na"]
                bulk_list.append(b)
                alk_list.append(tot["alk"] + n_extra - U[k])
            react = run_phreeqc(eng, water, tot, bulk_list, alk_list, fe, al)
            nres = [nicad_at(r, mFA, al)[1] for r in react]
            break
        if it >= maxit:
            raise RuntimeError(f"no convergence {water} Fe={fe} DOC={doc}: dmax={dmax}")
    out = []
    for k, r in enumerate(react):
        rec = dict(r)
        rec.update(water=water, Fe_mgL=fe, DOC_mgL=doc, FA_per_DOC=fa_factor, FA_mgL=fa_factor * doc, Al_control=al,
                   pH_target=PH_GRID[k], iterations=it, Na_extra_eq=n_extra,
                   init_alk_eq_kgw=tot["alk"])
        wkg = r["mass_water_kg"]
        for c in COMPS:
            rec[f"{c}_total_mol"] = T[c] + (n_extra if c == "Na" else 0.0)
            rec[f"{c}_FA_held_mol"] = held[k][c] if mFA > 0 else 0.0
            rec[f"{c}_massbal_relerr"] = ((r[f"{c}_diss_molal"] * wkg + r[f"{c}_sorbed_mol"] +
                                           rec[f"{c}_FA_held_mol"]) / rec[f"{c}_total_mol"] - 1.0)
        if mFA > 0:
            res = nres[k]
            rec.update(FA_logchi=res["logchi"], FA_VD_Lkg=res["VD"], FA_q_eqkg=res["q"],
                       FA_QH_molkg=res["Q"]["H"][0], FA_donnan_resid=res["residual"],
                       FA_U_eqL=sum(ZC[c] * held[k][c] for c in ZC))
            if al:
                rec["Al_FA_specific_molkg"] = res["Q"]["Al"][0]
                rec["Al_FA_held_mol"] = mFA * res["Q"]["Al"][0]
            for c in FA_BINDERS:
                rec[f"{c}_FA_specific_molkg"] = res["Q"][c][0]
                rec[f"{c}_FA_dist2_share"] = res["Q"][c][2] / res["Q"][c][0] if res["Q"][c][0] > 0 else np.nan
        out.append(rec)
    return out


def main():
    eng = Engine(db_path("minteq.v4.dat"))
    rows = []
    t0 = time.time()
    scen = []
    for water in M.WATERS:
        for fe in (0.0, 1.0):
            for doc in (0.0, 1.0, 3.0, 10.0):
                scen.append((water, fe, doc, 1.65))
            for doc in (1.0, 3.0, 10.0):
                scen.append((water, fe, doc, 1.27))
            for doc in (0.0, 1.0, 3.0, 10.0):
                scen.append((water, fe, doc, 1.65, True))
    for sc in scen:
        water, fe, doc, f = sc[:4]
        al = len(sc) > 4
        eng = Engine(db_path("minteq.v4.dat"))     # fresh instance per scenario
        out = run_scenario(eng, water, fe, doc, f, al=al)
        rows.extend(out)
        print(f"{water:9s} Fe={fe:g} DOC={doc:g} f={f} Al={al}: iterations={out[0]['iterations']} "
              f"({time.time() - t0:.0f} s)", flush=True)
    df = pd.DataFrame(rows).drop(columns=["state"])
    df.insert(0, "model", [(f"NICAD_genericFA_{f}xDOC+{MODEL}" if d > 0 else MODEL) + ("+AlOH3" if a else "")
                           for f, d, a in zip(df.FA_per_DOC, df.DOC_mgL, df.Al_control)])
    path = os.path.join(M.ROOT, "outputs", "nicad_raw.csv")
    df.to_csv(path, index=False)
    print("wrote", path, df.shape)


if __name__ == "__main__":
    main()
