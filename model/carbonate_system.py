"""
carbonate_system.py

Forward/inverse freshwater carbonate equilibrium engine for the Waikato River
CO2 / alkalinity / pH work.

This extends the single-point calculator in pco2_est.xlsx (Alk, pH, T -> pCO2)
into a general engine that solves between any two of {Alkalinity, DIC, pH, pCO2}
at a given temperature, so it can run forward (e.g. Alk + a forcing pCO2 or a
metabolically-perturbed DIC -> predicted pH) as well as reproduce the sheet's
original inverse direction.

Equilibrium constants
----------------------
K0 (CO2 Henry's constant) and K1 (first dissociation of carbonic acid) are
taken verbatim from pco2_est.xlsx, i.e. Plummer & Busenberg (1982), freshwater,
ionic strength = 0.

K2 (second dissociation, HCO3- <-> H+ + CO3--) is added from the same P&B
(1982) parameterisation. It isn't in the spreadsheet (which assumes alkalinity
is carried by HCO3- alone below pH ~8.3), so it was cross-checked against the
CO3-2 + H+ = HCO3- reaction in the PHREEQC database (Vitens/phreeqpython
phreeqc.dat), which carries the classic 5-term P&B analytical form and gives
log_k = 10.329 at 25 C (pK2 = 10.33), matching the standard freshwater
carbonate-chemistry literature value.

Kw (water dissociation) is included for completeness/robustness at higher pH,
from the same database's H2O = OH- + H+ reaction (Harned & Owen-type fit).
Per the spreadsheet's own assumption #1, the CO3-- and OH- terms are
negligible below pH ~8.3 -- they matter here mainly as headroom for the
photosynthesis-driven daytime pH excursions and any higher-pH reaches.

Units
-----
T          deg C
Alkalinity mg/L as CaCO3 (as in the spreadsheet's yellow input cell)
pCO2       atm (multiply by 1e6 for ppm/uatm)
DIC, species concentrations   mol/L
pH         free-hydrogen-ion scale (-log10[H+]), concentrations not activities
           (spreadsheet assumption #2 -- a few-percent effect at Waikato's low
           ionic strength)

Not included: non-carbonate alkalinity (organic acids/DOM, borate, ammonia).
The Waikato is a "black-ish" lowland river with measurable DOC, so measured
total alkalinity can run a little above carbonate alkalinity -- this is
exactly the kind of offset the equilibration-experiment ground-truthing
should reveal (see model notes at the bottom of this file).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

# Equivalent weight of CaCO3 used for the mg/L <-> eq/L conversion.
# Matches pco2_est.xlsx exactly (cell E10: =B5/50.04/1000).
EQ_WEIGHT_CACO3 = 50.04  # g/eq


# --------------------------------------------------------------------------
# Temperature-dependent equilibrium constants (Plummer & Busenberg, 1982)
# --------------------------------------------------------------------------

def _TK(T_C):
    """deg C -> Kelvin."""
    return np.asarray(T_C, dtype=float) + 273.15


def log_K0(T_C):
    """log10(K0), CO2 Henry's-law solubility constant, mol/L/atm.
    Verbatim from pco2_est.xlsx cell E3."""
    T = _TK(T_C)
    return 108.3865 + 0.01985076 * T - 6919.53 / T - 40.45154 * np.log10(T) + 669365 / T ** 2


def log_K1(T_C):
    """log10(K1), CO2 + H2O <-> H+ + HCO3-.
    Verbatim from pco2_est.xlsx cell E5."""
    T = _TK(T_C)
    return -356.3094 - 0.06091964 * T + 21834.37 / T + 126.8339 * np.log10(T) - 1684915 / T ** 2


def log_K2(T_C):
    """log10(K2), HCO3- <-> H+ + CO3--. Plummer & Busenberg (1982);
    see module docstring for the cross-check."""
    T = _TK(T_C)
    return -107.8871 - 0.03252849 * T + 5151.79 / T + 38.92561 * np.log10(T) - 563713.9 / T ** 2


def log_Kw(T_C):
    """log10(Kw), H2O <-> H+ + OH-."""
    T = _TK(T_C)
    return 293.29227 + 0.1360833 * T - 10576.913 / T - 123.73158 * np.log10(T) - 6.996455e-5 * T ** 2


def constants(T_C):
    """Return {K0, K1, K2, Kw} (not logs) at temperature T_C."""
    return dict(
        K0=10 ** log_K0(T_C),
        K1=10 ** log_K1(T_C),
        K2=10 ** log_K2(T_C),
        Kw=10 ** log_Kw(T_C),
    )


# --------------------------------------------------------------------------
# Unit helpers
# --------------------------------------------------------------------------

def alk_mgL_to_eqL(alk_mgL):
    """Alkalinity mg/L as CaCO3 -> eq/L. Matches the spreadsheet exactly."""
    return np.asarray(alk_mgL, dtype=float) / EQ_WEIGHT_CACO3 / 1000.0


def alk_eqL_to_mgL(alk_eqL):
    return np.asarray(alk_eqL, dtype=float) * EQ_WEIGHT_CACO3 * 1000.0


# --------------------------------------------------------------------------
# Core speciation
# --------------------------------------------------------------------------

def species_from_pH_DIC(pH, DIC, T_C):
    """Full carbonate speciation from pH and DIC (mol/L) at temperature T_C."""
    K = constants(T_C)
    H = 10.0 ** (-np.asarray(pH, dtype=float))
    denom = H ** 2 + K["K1"] * H + K["K1"] * K["K2"]
    a_co2 = H ** 2 / denom
    a_hco3 = K["K1"] * H / denom
    a_co3 = K["K1"] * K["K2"] / denom
    CO2 = DIC * a_co2
    HCO3 = DIC * a_hco3
    CO3 = DIC * a_co3
    OH = K["Kw"] / H
    pCO2 = CO2 / K["K0"]
    Alk_eqL = HCO3 + 2 * CO3 + OH - H
    return dict(H=H, CO2=CO2, HCO3=HCO3, CO3=CO3, OH=OH, pCO2=pCO2, Alk_eqL=Alk_eqL, DIC=DIC)


# --------------------------------------------------------------------------
# Solve directions
# --------------------------------------------------------------------------

def pCO2_from_Alk_pH(alk_mgL, pH, T_C, hco3_only=False):
    """(Alkalinity, pH, T) -> pCO2, DIC.

    This is the spreadsheet's original direction. With hco3_only=True it
    reproduces pco2_est.xlsx to full numerical precision (assumption #1:
    alkalinity carried by HCO3- alone). The default (hco3_only=False) instead
    solves the full charge-balance alkalinity equation, which only starts to
    diverge from the sheet above ~pH 8.3.
    """
    K = constants(T_C)
    H = 10.0 ** (-np.asarray(pH, dtype=float))
    alk_eqL = alk_mgL_to_eqL(alk_mgL)
    if hco3_only:
        HCO3 = alk_eqL
    else:
        HCO3 = (alk_eqL - K["Kw"] / H + H) / (1 + 2 * K["K2"] / H)
    CO2 = HCO3 * H / K["K1"]
    CO3 = K["K2"] * HCO3 / H
    pCO2 = CO2 / K["K0"]
    DIC = CO2 + HCO3 + CO3
    return dict(pCO2_atm=pCO2, pCO2_ppm=pCO2 * 1e6, DIC=DIC, HCO3=HCO3, CO2=CO2, CO3=CO3)


def pH_from_Alk_pCO2(alk_mgL, pCO2_atm, T_C):
    """(Alkalinity, pCO2, T) -> pH, DIC.

    THE forward-model direction: given a conservative alkalinity and a target
    pCO2 (a rising atmospheric background, an experimental headspace level in
    the equilibration rig, or a river pCO2 under a given scenario), solve the
    equilibrium pH. Scalar in pCO2_atm/T_C (loop or vectorise externally for
    scenario sweeps).
    """
    K = constants(T_C)
    alk_eqL = float(alk_mgL_to_eqL(alk_mgL))
    CO2 = K["K0"] * pCO2_atm

    def F(pH):
        H = 10.0 ** (-pH)
        HCO3 = K["K1"] * CO2 / H
        CO3 = K["K2"] * HCO3 / H
        OH = K["Kw"] / H
        return (HCO3 + 2 * CO3 + OH - H) - alk_eqL

    pH = brentq(F, 2.0, 12.0, xtol=1e-10)
    H = 10.0 ** (-pH)
    HCO3 = K["K1"] * CO2 / H
    CO3 = K["K2"] * HCO3 / H
    DIC = CO2 + HCO3 + CO3
    return dict(pH=pH, DIC=DIC, HCO3=HCO3, CO3=CO3, CO2=CO2)


def pH_from_Alk_DIC(alk_mgL, DIC, T_C):
    """(Alkalinity, DIC, T) -> pH.

    The coupling point for the metabolic layer: hold Alk fixed, perturb DIC by
    a net community production/respiration CO2 flux (mol/L over the mixing
    depth/residence time of interest), and re-solve pH. Scalar in DIC/T_C.
    """
    alk_eqL = float(alk_mgL_to_eqL(alk_mgL))

    def F(pH):
        return species_from_pH_DIC(pH, DIC, T_C)["Alk_eqL"] - alk_eqL

    pH = brentq(F, 2.0, 12.0, xtol=1e-10)
    return species_from_pH_DIC(pH, DIC, T_C)


if __name__ == "__main__":
    # Quick self-check against pco2_est.xlsx (Alk=40 mg/L CaCO3, pH=7, T=18C)
    out = pCO2_from_Alk_pH(40, 7, 18, hco3_only=True)
    print("hco3_only=True  (spreadsheet mode):", out["pCO2_ppm"], "ppm  (sheet: 4781.72 ppm)")
    out_full = pCO2_from_Alk_pH(40, 7, 18, hco3_only=False)
    print("hco3_only=False (full alk balance):", out_full["pCO2_ppm"], "ppm")
