"""Validation of the NICA-Donnan implementation (nicadonnan.py) against published numbers.

V1  Milne et al. (2001) ES&T 35:2049, eq 3 (p. 2058): explicit double Langmuir-Freundlich approximation of the
    generic FA and HA proton charge at a fixed ionic strength of 0.1 M. We compute Gamma = Qmax1H + Qmax2H - Q_H
    with the full NICA-Donnan model (0.1 M NaCl, H+ concentration from pH by the Davies equation) and compare.
V2  Rey-Castro et al. (2009) ES&T 43:7184, Table 1 (p. 7186): Visual MINTEQ speciation of a river water (Arve) with
    0.5 mg/L generic fulvic acid at pH 7, I = 6.4e-3 M. Given their bulk free-ion concentrations we solve the Donnan
    balance and compare the Boltzmann factor (their log(cD/c) = 1.17 z) and Donnan concentrations.
V3  Rey-Castro et al. (2009) Table 2 (p. 7188): average conditional affinities <log k'_i,j>_c from their eq 9, using
    their Donnan concentrations (Table 1). Tests the parameter table and the competition terms.
V4  Rey-Castro et al. (2009) Table 1: forward calculation of the free Cu, Zn, Ni, Pb, Cd concentrations from the
    stated totals, with inorganic complexation from PHREEQC (minteq.v4.dat, 25 C, Davies) and our NICA-Donnan.
    Compared with the free concentrations Visual MINTEQ gave them.
V5  Milne et al. (2003) ES&T 37:958, p. 969 text and Fig. 4: generic HA, 10 nM free metal, 1 mM free Ca, binding of
    Zn rises by about 1.85 log units between pH 4 and 8. Ionic strength not stated; we bracket it.

Writes validation/nicad_validation_*.csv and prints a summary.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nicadonnan as ND  # noqa: E402
from engine import Engine, db_path  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VAL = os.path.join(ROOT, "validation")


def davies_gamma(I, z, A=0.509):
    s = np.sqrt(I)
    return 10 ** (-A * z * z * (s / (1 + s) - 0.3 * I))


# ------------------------------------------------------------------------------------------- V1
def v1():
    eq3 = {"FA": (5.74, 5.30e3, 0.36, 2.01, 1.61e10, 0.45), "HA": (2.51, 1.39e4, 0.46, 3.33, 1.00e9, 0.19)}
    rows = []
    I = 0.1
    for humic, (a1, b1, g1, a2, b2, g2) in eq3.items():
        sub = ND.SUBSTRATE[humic]
        for pH in np.round(np.arange(3.0, 10.51, 0.25), 2):
            aH = 10 ** -pH
            cH = aH / davies_gamma(I, 1)
            res = ND.solve_donnan({"H": cH}, I, humic, extra_ions={"Na": (I, 1), "Cl": (I, -1)})
            gam_model = sub["Qmax1"] + sub["Qmax2"] - res["Q"]["H"][0]
            t1 = (b1 * aH) ** g1
            t2 = (b2 * aH) ** g2
            gam_eq3 = a1 * (1 - t1 / (1 + t1)) + a2 * (1 - t2 / (1 + t2))
            rows.append(dict(test="V1_Milne2001_eq3", humic=humic, I=I, pH=pH, Gamma_NICAD=gam_model,
                             Gamma_eq3=gam_eq3, diff=gam_model - gam_eq3, logchi=res["logchi"]))
    df = pd.DataFrame(rows)
    out = []
    for humic, g in df.groupby("humic"):
        g = g.set_index("pH")
        d_model = g.Gamma_NICAD[9.0] - g.Gamma_NICAD[4.0]
        d_eq3 = g.Gamma_eq3[9.0] - g.Gamma_eq3[4.0]
        m = (g.index >= 3.5) & (g.index <= 10.5)
        out.append(dict(humic=humic, max_abs_diff_eqkg_pH3_5_10_5=g["diff"][m].abs().max(),
                        rmse=np.sqrt((g["diff"][m] ** 2).mean()), dGamma_pH4_9_model=d_model,
                        dGamma_pH4_9_eq3=d_eq3, rel_diff_dGamma=(d_model - d_eq3) / d_eq3))
    return df, pd.DataFrame(out)


# ------------------------------------------------------------------------------------------- Rey-Castro data
# Rey-Castro et al. 2009 Table 1 (p. 7186): total, bulk free, Donnan-phase free (mol/L)
RC_T1 = {
    #        total      free bulk    free Donnan
    "Ca": (1.56e-3, 1.46e-3, 3.12e-1),
    "Mg": (3.00e-4, 2.84e-4, 6.08e-2),
    "Na": (2.60e-4, 2.59e-4, 3.80e-3),
    "Cl": (2.14e-4, 2.13e-4, 1.46e-5),
    "NO3": (6.60e-5, 6.58e-5, 4.49e-6),
    "K": (3.90e-5, 3.89e-5, 5.69e-4),
    "FeIII": (1.10e-6, 1.20e-14, 3.75e-11),
    "Al": (5.20e-7, 7.39e-12, 2.32e-8),
    "Cu": (5.00e-8, 6.19e-9, 1.33e-6),
    "Zn": (5.00e-8, 4.51e-8, 9.67e-6),
    "Ni": (1.00e-8, 8.54e-9, 1.83e-6),
    "Pb": (5.00e-9, 3.82e-10, 8.18e-8),
    "Cd": (5.00e-9, 4.28e-9, 9.17e-7),
}
RC_CO3T, RC_CO3_FREE = 1.20e-3, 9.82e-4
RC_SO4T, RC_SO4_FREE = 6.23e-4, 5.21e-4
RC_I = 6.4e-3
RC_FA_KG_L = 0.5e-6
RC_LOGCDH = -5.80
RC_LOGCHI = 1.17
# Rey-Castro Table 2 (p. 7188): <log k'_i,1>_c, <log k'_i,2>_c, <log k'_i>_c (global)
RC_T2 = {"Ca": (-2.2, -5.8, -2.9), "Mg": (-2.2, -6.2, -2.9), "Al": (-4.2, 5.1, -2.6), "H": (2.3, 6.5, 3.3),
         "Cu": (0.1, 2.1, 0.5), "Zn": (-3.9, -4.4, -4.0), "Ni": (-2.2, -2.2, -2.2), "Pb": (-1.3, 3.7, -0.1),
         "Cd": (-1.1, -4.0, -1.6), "Mn": (-1.6, -5.1, -2.2)}


def rc_bulk_H():
    # bulk pH 7.00; H+ concentration via Davies at I = 6.4e-3 (they used SIT; see notes)
    return 10 ** -7.0 / davies_gamma(RC_I, 1)


def v2(include_fe):
    bound = ["H", "Ca", "Mg", "Al", "Cu", "Zn", "Ni", "Pb", "Cd"] + (["FeIII"] if include_fe else [])
    bulk = {i: RC_T1[i][1] for i in bound if i != "H"}
    bulk["H"] = rc_bulk_H()
    # HCO3- is the dominant carbonate species at pH 7; the table's '(CO3)T free' value (9.82e-4) is used as a
    # monovalent anion, sulfate as divalent. Anions contribute < 0.5 % of the Donnan balance here.
    extra = {"Na": (RC_T1["Na"][1], 1), "K": (RC_T1["K"][1], 1), "Cl": (RC_T1["Cl"][1], -1),
             "NO3": (RC_T1["NO3"][1], -1), "HCO3": (RC_CO3_FREE, -1), "SO4": (RC_SO4_FREE, -2)}
    res = ND.solve_donnan(bulk, RC_I, "FA", bound_ions=bound, extra_ions=extra)
    rows = []
    for i in bound:
        rep = RC_LOGCDH if i == "H" else np.log10(RC_T1[i][2])
        rows.append(dict(test="V2_ReyCastro2009_T1_Donnan", Fe_included=include_fe, ion=i,
                         log_cD_ours=np.log10(res["cD"][i]), log_cD_reported=rep,
                         diff=np.log10(res["cD"][i]) - rep))
    rows.append(dict(test="V2_ReyCastro2009_T1_Donnan", Fe_included=include_fe, ion="logchi",
                     log_cD_ours=res["logchi"], log_cD_reported=RC_LOGCHI, diff=res["logchi"] - RC_LOGCHI))
    # implied check with their own cD values: humic charge from NICA at their cD, vs the Donnan excess they imply
    cD_rep = {i: RC_T1[i][2] for i in bound if i != "H"}
    cD_rep["H"] = 10 ** RC_LOGCDH
    Qb = ND.nica_bound(cD_rep, "FA")
    q = ND.humic_charge(Qb, "FA")
    VD = ND.donnan_volume(RC_I, ND.SUBSTRATE["FA"]["b"])
    chi = 10 ** RC_LOGCHI
    don = sum(z * c * (chi ** z - 1) for c, z in [(RC_T1[k][1], 2) for k in ("Ca", "Mg")] +
              [(RC_T1[k][1], 1) for k in ("Na", "K")] + [(rc_bulk_H(), 1)] +
              [(RC_T1[k][1], 3) for k in (["Al"] + (["FeIII"] if include_fe else []))] +
              [(RC_T1[k][1], 2) for k in ("Cu", "Zn", "Ni", "Pb", "Cd")] +
              [(RC_T1["Cl"][1], -1), (RC_T1["NO3"][1], -1), (RC_CO3_FREE, -1), (RC_SO4_FREE, -2)])
    rows.append(dict(test="V2b_charge_at_reported_cD", Fe_included=include_fe, ion="q_NICA_eqkg",
                     log_cD_ours=q, log_cD_reported=-VD * don, diff=q + VD * don))
    return pd.DataFrame(rows), res


def v3():
    """Rey-Castro eq 9 at their Donnan concentrations (Fe(III) omitted: its contribution is stated to be several
    orders of magnitude smaller than the other ions, p. 7189)."""
    par, sub = ND.ION_PARAMS["FA"], ND.SUBSTRATE["FA"]
    cD = {i: RC_T1[i][2] for i in ("Ca", "Mg", "Al", "Cu", "Zn", "Ni", "Pb", "Cd")}
    cD["H"] = 10 ** RC_LOGCDH
    p = (sub["p1"], sub["p2"])
    Qm = (sub["Qmax1"], sub["Qmax2"])
    nH = (par["H"][1], par["H"][3])
    rows = []
    for i in RC_T2:
        per = []
        w = []
        for j in (0, 1):
            logk, n = (par[i][0], par[i][1]) if j == 0 else (par[i][2], par[i][3])
            S = sum(10 ** ((par[m][1] if j == 0 else par[m][3]) * ((par[m][0] if j == 0 else par[m][2])
                                                                    + np.log10(c)))
                    for m, c in cD.items() if m != i)
            per.append(logk - np.log10(1 + S ** p[j]) / (n * p[j]))
            w.append(n / nH[j] * Qm[j])
        glob = (w[0] * per[0] + w[1] * per[1]) / (w[0] + w[1])
        rep = RC_T2[i]
        rows.append(dict(test="V3_ReyCastro2009_T2_eq9", ion=i, carb_ours=per[0], carb_rep=rep[0],
                         phen_ours=per[1], phen_rep=rep[1], global_ours=glob, global_rep=rep[2],
                         max_abs_diff=max(abs(per[0] - rep[0]), abs(per[1] - rep[1]), abs(glob - rep[2]))))
    return pd.DataFrame(rows)


def v4():
    """Forward calculation of free trace metals from totals (Rey-Castro Table 1)."""
    eng = Engine(db_path("minteq.v4.dat"))
    sol = f"""SOLUTION 1
   -temp 25
   -units mol/L
   pH 7.0
   Ca {RC_T1['Ca'][0]}
   Mg {RC_T1['Mg'][0]}
   Na {RC_T1['Na'][0]}
   K {RC_T1['K'][0]}
   Cl {RC_T1['Cl'][0]}
   N(5) {RC_T1['NO3'][0]}
   C(4) {RC_CO3T}
   S(6) {RC_SO4T}
   Cu(2) {RC_T1['Cu'][0]}
   Zn {RC_T1['Zn'][0]}
   Ni {RC_T1['Ni'][0]}
   Pb {RC_T1['Pb'][0]}
   Cd {RC_T1['Cd'][0]}
SELECTED_OUTPUT
   -reset false
USER_PUNCH
   -headings mu Ca_free Mg_free Na_free K_free Cu_T Cu_free Zn_T Zn_free Ni_T Ni_free Pb_T Pb_free Cd_T Cd_free HCO3 SO4free H_molal
   -start
   10 PUNCH MU, MOL("Ca+2"), MOL("Mg+2"), MOL("Na+"), MOL("K+")
   20 PUNCH TOT("Cu"), MOL("Cu+2"), TOT("Zn"), MOL("Zn+2"), TOT("Ni"), MOL("Ni+2"), TOT("Pb"), MOL("Pb+2"), TOT("Cd"), MOL("Cd+2")
   30 PUNCH MOL("HCO3-"), MOL("SO4-2"), MOL("H+")
   -end
END
"""
    r = eng.run(sol)[-1]
    R = {m: r[f"{m}_T"] / r[f"{m}_free"] for m in ("Cu", "Zn", "Ni", "Pb", "Cd")}
    rows = []
    for variant, use_rc_I in (("PHREEQC_I_and_activity", False), ("RC_I_6.4e-3", True)):
        I = RC_I if use_rc_I else r["mu"]
        cH = 10 ** -7.0 / davies_gamma(I, 1)
        # majors: bulk free values as reported (Ca, Mg, Na, K); Al as reported (fixed)
        base = {"H": cH, "Ca": RC_T1["Ca"][1], "Mg": RC_T1["Mg"][1], "Al": RC_T1["Al"][1]}
        extra = {"Na": (RC_T1["Na"][1], 1), "K": (RC_T1["K"][1], 1), "Cl": (RC_T1["Cl"][1], -1),
                 "NO3": (RC_T1["NO3"][1], -1), "HCO3": (RC_CO3_FREE, -1), "SO4": (RC_SO4_FREE, -2)}
        free = {m: RC_T1[m][1] for m in ("Cu", "Zn", "Ni", "Pb", "Cd")}
        for it in range(60):
            bulk = dict(base)
            bulk.update(free)
            res = ND.solve_donnan(bulk, I, "FA", extra_ions=extra)
            new = {}
            for m in free:
                T = RC_T1[m][0]

                def g(logc, m=m):
                    b2 = dict(bulk)
                    b2[m] = 10 ** logc
                    cD = {k: v * res["chi"] ** ND.CHARGE[k] for k, v in b2.items()}
                    Q = ND.nica_bound(cD, "FA")[m][0]
                    held = RC_FA_KG_L * (Q + res["VD"] * 10 ** logc * (res["chi"] ** 2 - 1))
                    return np.log10(10 ** logc * R[m] + held) - np.log10(T)
                new[m] = 10 ** brentq(g, -20, np.log10(T) + 1, xtol=1e-12)
            dmax = max(abs(np.log10(new[m] / free[m])) for m in free)
            free = new
            if dmax < 1e-9:
                break
        for m in free:
            rows.append(dict(test="V4_ReyCastro2009_T1_free_from_total", variant=variant, ion=m,
                             I_used=I, logchi=res["logchi"], R_inorg_PHREEQC=R[m],
                             log_free_ours=np.log10(free[m]), log_free_reported=np.log10(RC_T1[m][1]),
                             diff=np.log10(free[m]) - np.log10(RC_T1[m][1]),
                             frac_FA_ours=1 - free[m] * R[m] / RC_T1[m][0]))
    return pd.DataFrame(rows), r


def v5():
    rows = []
    for I in (0.003, 0.01, 0.03, 0.1):
        vals = {}
        for pH in (4.0, 6.0, 8.0):
            cH = 10 ** -pH / davies_gamma(I, 1)
            bulk = {"H": cH, "Ca": 1e-3, "Zn": 1e-8}
            nacl = I - 3e-3   # background 1:1 salt making up the ionic strength (assumed)
            extra = {"Cl": (2e-3, -1)}
            if nacl > 0:
                extra.update({"Na": (nacl, 1), "NO3": (nacl, -1)})
            res = ND.solve_donnan(bulk, I, "HA", extra_ions=extra)
            vals[pH] = np.log10(res["Q"]["Zn"][0])
        rows.append(dict(test="V5_Milne2003_Fig4_Zn_HA", I=I, logQZn_pH4=vals[4.0], logQZn_pH6=vals[6.0],
                         logQZn_pH8=vals[8.0], rise_pH4_8=vals[8.0] - vals[4.0], reported_rise=1.85,
                         implied_exchange_ratio=(vals[8.0] - vals[4.0]) / 4.0, reported_ratio=0.46))
    return pd.DataFrame(rows)


def main():
    os.makedirs(VAL, exist_ok=True)
    d1, s1 = v1()
    d1.to_csv(os.path.join(VAL, "nicad_validation_V1_milne2001_eq3.csv"), index=False)
    print("V1 Milne 2001 eq 3 (I = 0.1 M)\n", s1.to_string(index=False))
    t2a, resA = v2(False)
    t2b, resB = v2(True)
    t2 = pd.concat([t2a, t2b])
    t2.to_csv(os.path.join(VAL, "nicad_validation_V2_reycastro_donnan.csv"), index=False)
    print("\nV2 Rey-Castro Table 1 Donnan\n", t2.to_string(index=False))
    t3 = v3()
    t3.to_csv(os.path.join(VAL, "nicad_validation_V3_reycastro_eq9.csv"), index=False)
    print("\nV3 Rey-Castro Table 2 (eq 9)\n", t3.round(2).to_string(index=False))
    t4, r = v4()
    t4.to_csv(os.path.join(VAL, "nicad_validation_V4_reycastro_free.csv"), index=False)
    print("\nV4 Rey-Castro Table 1 free metal from total\n", t4.round(3).to_string(index=False))
    print("PHREEQC ionic strength for the Arve water:", r["mu"])
    t5 = v5()
    t5.to_csv(os.path.join(VAL, "nicad_validation_V5_milne2003_fig4.csv"), index=False)
    print("\nV5 Milne 2003 Fig 4 / text (Zn, generic HA)\n", t5.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
