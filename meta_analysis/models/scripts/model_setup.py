"""Shared definitions for the PHREEQC hydrous ferric oxide (HFO) calculations.

Everything that is an assumption is defined here, so the README can point to one place.
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TEMP_C = 16.0

# ---------------------------------------------------------------- waters (mg/L unless stated)
# Alkalinity in meq/L (as HCO3-). Na, K, Cl and SO4 were specified only for the soft water;
# we use the same values for the two harder waters. Charge balance is closed on the ion named in
# 'balance' (PHREEQC "charge" keyword) at the reference pH of the initial solution.
WATERS = {
    "soft": dict(alk_meq=0.78, Ca=5.3, Mg=2.0, Na=10.0, K=2.0, Cl=10.0, SO4=5.0, balance="Na",
                 label="Soft (Waikato-like): alk 0.78 meq/L, Ca 5.3, Mg 2 mg/L"),
    "moderate": dict(alk_meq=2.0, Ca=30.0, Mg=6.0, Na=10.0, K=2.0, Cl=10.0, SO4=5.0, balance="Cl",
                     label="Moderately hard: alk 2 meq/L, Ca 30, Mg 6 mg/L"),
    "hard": dict(alk_meq=4.0, Ca=60.0, Mg=12.0, Na=10.0, K=2.0, Cl=10.0, SO4=5.0, balance="Cl",
                 label="Hard: alk 4 meq/L, Ca 60, Mg 12 mg/L"),
}

# gram formula weights used to convert mg/L to mmol/L (values as in minteq.v4.dat master species)
GFW = {"Ca": 40.078, "Mg": 24.305, "Na": 22.9898, "K": 39.0983, "Cl": 35.453, "SO4": 96.0636,
       "Mn": 54.938, "Ni": 58.69, "Co": 58.9332, "Zn": 65.39, "Cd": 112.41, "Cu": 63.546,
       "Pb": 207.2, "P": 30.9738, "Mo": 95.94, "As": 74.9216, "V": 50.94, "Se": 78.96}

# ---------------------------------------------------------------- trace elements (ug/L, system totals)
# key: element symbol used in outputs; (PHREEQC input name with fixed redox state, ug/L)
TRACES = {
    "Mn": ("Mn(2)", 20.0), "Ni": ("Ni", 1.0), "Co": ("Co(2)", 0.1), "Zn": ("Zn", 2.0),
    "Cd": ("Cd", 0.02), "Cu": ("Cu(2)", 1.0), "Pb": ("Pb", 0.1),
    "P": ("P", 20.0), "Mo": ("Mo", 0.5), "As": ("As(5)", 2.0), "V": ("V(5)", 1.0), "Se": ("Se(4)", 0.2),
}
CATIONS = ["Mn", "Ni", "Co", "Zn", "Cd", "Cu", "Pb"]
OXYANIONS = ["P", "Mo", "As", "V", "Se"]

# "free" species. Cations: the aquo ion. Oxyanions: all protonation forms of the monomeric oxyanion
# (i.e. not complexed with Ca, Mg, Na, K or trace metals). Species names as in minteq.v4.dat.
FREE_SPECIES = {
    "Mn": ["Mn+2"], "Ni": ["Ni+2"], "Co": ["Co+2"], "Zn": ["Zn+2"], "Cd": ["Cd+2"], "Cu": ["Cu+2"], "Pb": ["Pb+2"],
    "P": ["H3PO4", "H2PO4-", "HPO4-2", "PO4-3"],
    "Mo": ["H2MoO4", "HMoO4-", "MoO4-2"],
    "As": ["H3AsO4", "H2AsO4-", "HAsO4-2", "AsO4-3"],
    "V": ["VO2+", "H3VO4", "H2VO4-", "HVO4-2", "VO4-3"],
    "Se": ["H2SeO3", "HSeO3-", "SeO3-2"],
}

# ---------------------------------------------------------------- hydrous ferric oxide (Dzombak & Morel 1990)
# 89 g HFO per mol Fe, 600 m2/g, 0.2 mol weak sites and 0.005 mol strong sites per mol Fe.
# (Same values as the PHREEQC manual's example 8: 0.09 g HFO -> 2e-4 mol weak, 5e-6 mol strong.)
FE_GFW = 55.845
HFO_G_PER_MOL_FE = 89.0
HFO_AREA_M2_G = 600.0
WEAK_SITES_PER_FE = 0.2
STRONG_SITES_PER_FE = 0.005


def hfo_amounts(fe_mg_L):
    mol_fe = fe_mg_L * 1e-3 / FE_GFW
    return dict(mol_fe=mol_fe, g_hfo=mol_fe * HFO_G_PER_MOL_FE,
                weak=mol_fe * WEAK_SITES_PER_FE, strong=mol_fe * STRONG_SITES_PER_FE)


# ---------------------------------------------------------------- database modifications
# 1. Redox held fixed. In a PHREEQC batch reaction every redox element is re-equilibrated with pe.
#    To hold Mn(II), Co(II), Cu(II), As(V), V(V), Se(IV) and S(VI) fixed we redefine the reactions
#    that link the other valence states to the chosen one so that the other states are >=60 orders
#    of magnitude less stable at any pe between -10 and 20. Nothing else in the database changes.
REDOX_FIX_MINTEQV4 = """
SOLUTION_SPECIES
# --- redox held fixed (see model_setup.py); original log_k values in minteq.v4.dat given in comments
e- + Mn+3 = Mn+2
     log_k 100        # original 25.35
     delta_h 0 kJ
Mn+2 + 4 H2O = MnO4- + 8 H+ + 5 e-
     log_k -300       # original -127.794
     delta_h 0 kJ
     -gamma 3 0
Mn+2 + 4 H2O = MnO4-2 + 8 H+ + 4 e-
     log_k -300       # original -118.422
     delta_h 0 kJ
     -gamma 5 0
Co+3 + e- = Co+2
     log_k 100        # original 32.4
     delta_h 0 kJ
Cu+2 + e- = Cu+
     log_k -100       # original 2.69
     delta_h 0 kJ
H3AsO4 + 2 e- + 2 H+ = H3AsO3 + H2O
     log_k -100       # original 18.898
     delta_h 0 kJ
VO2+ + e- + 2 H+ = VO+2 + H2O
     log_k -100       # original 16.903
     delta_h 0 kJ
SeO4-2 + 2 e- + 3 H+ = HSeO3- + H2O
     log_k 100        # original 36.308 (makes Se(VI) negligible relative to Se(IV))
     delta_h 0 kJ
HSeO3- + 6 e- + 6 H+ = HSe- + 3 H2O
     log_k -100       # original 44.86
     delta_h 0 kJ
SO4-2 + 9 H+ + 8 e- = HS- + 4 H2O
     log_k -100       # original 33.66
     delta_h 0 kJ
"""

# 2. Surface constants that minteq.v4.dat lacks, copied verbatim from wateq4f.dat
#    (USGS phreeqc3 GitHub master, database/wateq4f.dat, SURFACE_SPECIES block):
#    Mn(II): "Derived constants table 10.5" of Dzombak & Morel (1990).
#    Carbonate: "Van Geen et al., 1994 reoptimized for HFO" (Appelo et al. 2002, ES&T 36:3096).
SURF_MN = """
SURFACE_SPECIES
# Mn(II): wateq4f.dat, Dzombak & Morel (1990) table 10.5
        Hfo_sOH + Mn+2 = Hfo_sOMn+ + H+
        log_k -0.4
        Hfo_wOH + Mn+2 = Hfo_wOMn+ + H+
        log_k -3.5
"""
SURF_CO3 = """
SURFACE_SPECIES
# Carbonate: wateq4f.dat, Van Geen et al. (1994) re-optimised for HFO (Appelo et al. 2002)
       Hfo_wOH + CO3-2 + H+ = Hfo_wCO3- + H2O
       log_k 12.56
       Hfo_wOH + CO3-2 + 2 H+ = Hfo_wHCO3 + H2O
       log_k 20.62
"""

PHASES_FIX = """
PHASES
Fix_H+
     H+ = H+
     log_k 0.0
END
"""

SI_PHASES = ["Calcite", "Rhodochrosite", "Smithsonite", "Otavite", "Tenorite", "Malachite", "Cerussite",
             "Hydrocerussite", "Ni(OH)2", "CoCO3", "Hydroxylapatite", "CaMoO4", "Zn(OH)2(am)"]
