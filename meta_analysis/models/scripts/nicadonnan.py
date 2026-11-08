"""Generic NICA-Donnan model (consistent NICA isotherm, bimodal, with Donnan electrostatics).

Equations
---------
Milne et al. (2003) ES&T 37:958, eq 1 (p. 958): amount of component i bound

    Q_i = sum_j (n_ij / n_Hj) * Qmax_jH * [ (K_ij c_Di)^n_ij / S_j ] * S_j^p_j / (1 + S_j^p_j)
    S_j = sum_k (K_kj c_Dk)^n_kj        (k runs over H+ and every competing metal ion)

where j = 1, 2 are the low-affinity ("carboxylic") and high-affinity ("phenolic") distributions and c_Di is the
concentration of the free ion i in the Donnan phase (mol/L). (Same form: Rey-Castro et al. 2009, eqs 1-3.)

Donnan volume: Milne et al. (2001) eq 2 / Milne et al. (2003) eq 2:  log10 V_D = b (1 - log10 I) - 1   (L/kg)

Donnan phase: c_Di = c_i * chi^z_i (Boltzmann factor chi, the same for every ion; Rey-Castro et al. 2009 Table 1
footnote: log(c_D,i / c_i) = 1.17 z_i in their example). chi is found from electroneutrality of the humic particle
plus its Donnan phase:

    q / V_D + sum_i z_i (c_Di - c_i) = 0,   q = -(Qmax_1H + Qmax_2H) + sum_i z_i Q_i     (eq/kg)

i.e. the residual (negative) charge left after specific binding of protons and metals is neutralised by the excess
of counter-ions in the Donnan volume (Milne et al. 2001, p. 2050; Milne et al. 2003, p. 958-959).

Concentrations, not activities, are used for c_i, as in the data treatment of Milne et al. (2001, p. 2051:
"Proton concentrations have been derived from pH measurements using the Davies equation") and Milne et al.
(2003, p. 959: bound metal "as a function of free metal-ion concentration, proton concentration, and ionic
strength (all in mol L-1)").

Parameters: Milne et al. (2003) Table 3 (p. 968), checked against the PDF page image. Proton parameters in that
table (FA: log K_H1 2.34, n_H1 0.66, log K_H2 8.60, n_H2 0.76) are the Milne et al. (2001) Table 4 generic values
(p. 2057: log K_H1 2.34, m_1 0.38, log K_H2 8.6, m_2 0.53, with m = n_H p; 0.66 x 0.59 = 0.389, 0.76 x 0.70 = 0.532).
"""
import numpy as np
from scipy.optimize import brentq

# ------------------------------------------------------------------ Milne et al. 2003, Table 3, p. 968
# ion: (log K1, n1, log K2, n2). Values printed in italics in the table (estimated, not fitted) are flagged in
# ESTIMATED below.
SUBSTRATE = {
    "FA": dict(b=0.57, Qmax1=5.88, Qmax2=1.86, p1=0.59, p2=0.70),
    "HA": dict(b=0.49, Qmax1=3.15, Qmax2=2.55, p1=0.62, p2=0.41),
}
ION_PARAMS = {
    "FA": {
        "H": (2.34, 0.66, 8.60, 0.76),
        "Al": (-4.11, 0.42, 12.16, 0.31),
        "Ba": (-2.6, 0.90, -4.8, 0.70),
        "Ca": (-2.13, 0.85, -3.0, 0.80),
        "Cd": (-0.97, 0.68, 0.50, 0.50),
        "Co": (-2.64, 0.71, -2.21, 0.65),
        "Cu": (0.26, 0.53, 8.26, 0.36),
        "FeII": (-1.02, 0.30, -1.1, 0.51),
        "FeIII": (6.0, 0.25, 36.0, 0.19),
        "Mg": (-2.1, 0.77, -2.4, 0.59),
        "Mn": (-1.55, 0.72, -1.1, 0.56),
        "Ni": (-2.07, 0.65, 2.03, 0.53),
        "Pb": (-1.16, 0.60, 6.92, 0.69),
        "Sr": (-2.5, 0.85, -4.6, 0.70),
        "Zn": (-3.84, 0.67, -0.73, 0.61),
    },
    "HA": {
        "H": (2.93, 0.81, 8.00, 0.63),
        "Al": (-1.05, 0.40, 8.89, 0.30),
        "Ba": (-1.1, 0.90, -0.7, 0.70),
        "Ca": (-1.37, 0.78, -0.43, 0.75),
        "Cd": (-0.20, 0.73, 2.37, 0.54),
        "Co": (-0.24, 0.79, 1.0, 0.66),
        "Cu": (2.23, 0.56, 6.85, 0.34),
        "FeII": (0.1, 0.30, 2.2, 0.50),
        "FeIII": (3.5, 0.30, 17.5, 0.25),
        "Mg": (-0.6, 0.77, 0.6, 0.59),
        "Mn": (-0.3, 0.72, 1.3, 0.56),
        "Ni": (-0.26, 0.64, 1.0, 0.55),
        "Pb": (1.25, 0.60, 4.84, 0.69),
        "Sr": (-1.36, 0.78, -0.43, 0.75),
        "Zn": (0.11, 0.67, 2.39, 0.27),
    },
}
# Which entries are italic (estimated) in Milne et al. 2003 Table 3 (read from the PDF page image, p. 968)
ESTIMATED = {
    "FA": {"Ba": "all", "Ca": "logK2,n2", "Cd": "n2", "Co": "n2", "FeII": "logK2,n2", "FeIII": "all",
           "Mg": "all", "Mn": "n1,logK2,n2", "Ni": "", "Pb": "n1,n2", "Sr": "all", "Zn": "n1,n2", "Cu": "",
           "Al": "", "H": ""},
    "HA": {"Ba": "all", "Co": "logK2", "FeII": "all", "Mg": "all", "Mn": "all", "Sr": "n1,logK2,n2",
           "Al": "n1,n2", "Ca": "", "Cd": "", "Cu": "", "FeIII": "n1,logK2,n2", "Ni": "logK2", "Pb": "",
           "Zn": "", "H": ""},
}
CHARGE = {"H": 1, "Al": 3, "Ba": 2, "Ca": 2, "Cd": 2, "Co": 2, "Cu": 2, "FeII": 2, "FeIII": 3, "Mg": 2, "Mn": 2,
          "Ni": 2, "Pb": 2, "Sr": 2, "Zn": 2}


def donnan_volume(I, b):
    """V_D in L/kg (Milne et al. 2001 eq 2)."""
    return 10 ** (b * (1.0 - np.log10(I)) - 1.0)


def nica_bound(cD, humic="FA"):
    """Specific (NICA) binding, mol/kg, for every ion in cD (dict ion -> Donnan-phase concentration, mol/L).

    Returns dict ion -> (Q_total, Q_dist1, Q_dist2) and the fractional coverages."""
    sub = SUBSTRATE[humic]
    par = ION_PARAMS[humic]
    nH = (par["H"][1], par["H"][3])
    Qmax = (sub["Qmax1"], sub["Qmax2"])
    p = (sub["p1"], sub["p2"])
    out = {i: [0.0, 0.0, 0.0] for i in cD}
    for j in (0, 1):
        terms = {}
        for i, c in cD.items():
            logK, n = (par[i][0], par[i][1]) if j == 0 else (par[i][2], par[i][3])
            terms[i] = 0.0 if c <= 0 else 10 ** (n * (logK + np.log10(c)))
        S = sum(terms.values())
        if S <= 0:
            continue
        frac_total = S ** p[j] / (1.0 + S ** p[j])
        for i, t in terms.items():
            n = par[i][1] if j == 0 else par[i][3]
            q = (n / nH[j]) * Qmax[j] * (t / S) * frac_total
            out[i][j + 1] = q
            out[i][0] += q
    return out


def humic_charge(Qb, humic="FA"):
    """q in eq/kg: -(Qmax1H + Qmax2H) + sum z_i Q_i."""
    sub = SUBSTRATE[humic]
    return -(sub["Qmax1"] + sub["Qmax2"]) + sum(CHARGE.get(i, 0) * v[0] for i, v in Qb.items()
                                                 if i in CHARGE)


def solve_donnan(bulk, I, humic="FA", bound_ions=None, extra_ions=None, logchi_bounds=(0.0, 8.0)):
    """Solve for the Boltzmann factor.

    bulk       : dict ion -> bulk free concentration (mol/L) for ions that bind specifically (keys of ION_PARAMS)
    extra_ions : dict name -> (concentration mol/L, charge) for ions that take part only in the Donnan balance
                 (Na+, K+, Cl-, SO4-2, HCO3-, ...).
    Returns dict with chi, logchi, VD, cD, Q (NICA bound, mol/kg), q (eq/kg), donnan_excess (mol/L of Donnan
    volume, i.e. c_D - c, per ion), all per kg humic.
    """
    if bound_ions is None:
        bound_ions = [i for i in bulk if i in ION_PARAMS[humic]]
    extra_ions = extra_ions or {}
    VD = donnan_volume(I, SUBSTRATE[humic]["b"])
    allions = {i: (bulk[i], CHARGE[i]) for i in bound_ions}
    allions.update(extra_ions)

    def f(logchi):
        chi = 10 ** logchi
        cD = {i: bulk[i] * chi ** CHARGE[i] for i in bound_ions}
        Qb = nica_bound(cD, humic)
        q = humic_charge(Qb, humic)
        don = sum(z * c * (chi ** z - 1.0) for c, z in allions.values())
        return q / VD + don

    lo, hi = logchi_bounds
    if f(lo) > 0:
        logchi = lo
    else:
        logchi = brentq(f, lo, hi, xtol=1e-12, rtol=1e-12, maxiter=500)
    chi = 10 ** logchi
    cD = {i: bulk[i] * chi ** CHARGE[i] for i in bound_ions}
    Qb = nica_bound(cD, humic)
    q = humic_charge(Qb, humic)
    excess = {i: c * (chi ** z - 1.0) for i, (c, z) in allions.items()}
    return dict(chi=chi, logchi=logchi, VD=VD, cD=cD, Q=Qb, q=q, excess=excess, residual=f(logchi))


def held_per_litre(res, mFA):
    """Moles per litre held by the humic 'unit' (specific + Donnan excess) for each ion, given mFA in kg/L."""
    held = {}
    for i, v in res["Q"].items():
        held[i] = mFA * v[0]
    for i, e in res["excess"].items():
        held[i] = held.get(i, 0.0) + mFA * res["VD"] * e
    return held
