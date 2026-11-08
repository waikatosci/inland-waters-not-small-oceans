"""
unified.py
----------
One chain from metabolism to mineral surface charge:

    MMRT metabolism (GPP, ER)  ->  net CO2 source F(T)
    F(T), geothermal load, gas exchange k(T), solubility K0(T)
                               ->  CO2 excess  dp = F / (k K0)          [additive]
    pCO2_w = pCO2_atm(t, SSP) + dp
    (Alk, pCO2_w, T)           ->  pH                                   [carbonate_system.py]
    PZC(T) - pH                ->  Coulombic departure  D               [pH units]
    D                          ->  sorption factor  Kd / Kd(PZC) = 10**(-z_eff * D)

Why the last line has no temperature in it: with a Nernstian surface
potential psi0 = (RT ln10 / F)(PZC - pH), the Boltzmann factor
exp(-z F psi0 / RT) = 10**(-z (PZC - pH)). The RT in the Nernst slope
cancels the RT in the Boltzmann factor exactly, so temperature reaches
the sorption partition only through PZC(T) and pH(T).

Supersaturation is carried as an ADDITIVE excess (a flux balance), not
as a fixed multiple of the atmosphere. A fixed multiple (the RS factor in
WR_Acidification_Data and the S draw in the Coulombic Cenozoic Monte
Carlo) implicitly assumes the internal CO2 source grows in proportion to
atmospheric CO2, which no mechanism supplies. `multiplicative=True`
reproduces the old convention for comparison.

All parameters flagged ILLUSTRATIVE are placeholders chosen to be
literature-typical; they are not fitted to Waikato data.
"""
from __future__ import annotations

import numpy as np

import carbonate_system as cs
import gas_exchange as ge

R = 8.314462618          # J/mol/K
KB_EV = 8.617333262e-5   # eV/K
EV_TO_J_MOL = 96485.33212

# ---------------------------------------------------------------------------
# 1. MMRT rates, parameterised by things you can read off data
# ---------------------------------------------------------------------------

def mmrt_params(Ea_ref_eV, Topt_C, Tref_C=15.0):
    """Convert (apparent activation energy at Tref, optimum temperature)
    into MMRT (dH_ref, dCp) at T0 = Tref.

    MMRT: ln k = ln(kB T/h) - dH(T)/RT + dS(T)/R, with
          dH(T) = dH0 + dCp (T - T0),  dS(T) = dS0 + dCp ln(T/T0).
    Apparent Ea(T) = dH(T) + RT. At Topt, d ln k / dT = 0  ->  dH(Topt) = -R Topt.
    So dH0 = Ea_ref - R T0, and dCp = -(R Topt + dH0) / (Topt - T0).
    """
    T0 = Tref_C + 273.15
    Topt = Topt_C + 273.15
    dH0 = Ea_ref_eV * EV_TO_J_MOL - R * T0
    dCp = -(R * Topt + dH0) / (Topt - T0)
    return dict(T0=T0, dH0=dH0, dCp=dCp)


def mmrt_rel(T_C, Ea_ref_eV, Topt_C, Tref_C=15.0):
    """Rate at T relative to rate at Tref (dimensionless). dS0 cancels."""
    p = mmrt_params(Ea_ref_eV, Topt_C, Tref_C)
    T = np.asarray(T_C, dtype=float) + 273.15
    T0, dH0, dCp = p["T0"], p["dH0"], p["dCp"]
    lnrel = (np.log(T / T0)
             - (dH0 + dCp * (T - T0)) / (R * T) + dH0 / (R * T0)
             + (dCp * np.log(T / T0)) / R)
    return np.exp(lnrel)


# Whole-ecosystem parameters from the literature synthesis (meta/synthesis/meta_metab.py): apparent
# activation energies at low temperature are the medians across studies that follow temperature within
# freshwater systems (seasonal, diel, warming experiments, geothermal gradients): GPP 0.53 eV (9 studies),
# ER 0.62 eV (11 studies). Optima are typical of phytoplankton
# communities (GPP) and ecosystem respiration (ER); freshwater whole-ecosystem optima were not found.
METAB = dict(Ea_P=0.53, Topt_P=27.0, Ea_R=0.62, Topt_R=35.0, Tref=15.0)
# Earlier illustrative set (organism-level metabolic theory: GPP 0.32, ER 0.65 eV), kept for comparison.
METAB_THEORY = dict(METAB, Ea_P=0.32, Ea_R=0.65)
# Alternative: equal apparent activation energies (0.57 eV at 15 C), as estimated for whole
# streams (Demars et al. 2016); optima unchanged. Used for the "R = P" projections.
METAB_EQ = dict(METAB, Ea_P=0.57, Ea_R=0.57)


def net_source_rel(T_C, T0_C, g, metab=METAB):
    """Net internal CO2 source F(T) / F(T0).

    F = ER - GPP (ER includes benthic and Corbicula respiration).
    g = GPP/ER at T0. Written this way,

        F/F0 = (rR(T) - g rP(T)) / (1 - g)

    where rX = X(T)/X(T0). The factor 1/(1-g) is the "net heterotrophy
    leverage": the closer GPP is to ER, the more a small change in either
    moves the net source. Valid for g < 1 (net heterotrophic locally).
    """
    if not g < 1:
        raise ValueError("g = GPP/ER must be < 1 for a net internal CO2 source")
    rP = mmrt_rel(T_C, metab["Ea_P"], metab["Topt_P"], metab["Tref"]) / \
         mmrt_rel(T0_C, metab["Ea_P"], metab["Topt_P"], metab["Tref"])
    rR = mmrt_rel(T_C, metab["Ea_R"], metab["Topt_R"], metab["Tref"]) / \
         mmrt_rel(T0_C, metab["Ea_R"], metab["Topt_R"], metab["Tref"])
    return (rR - g * rP) / (1.0 - g)


def apparent_Ea_net(T_C, g, metab=METAB, dT=0.01):
    """Apparent activation energy (eV) of the net source F at T."""
    T = np.asarray(T_C, dtype=float)
    f1 = net_source_rel(T - dT, T, g, metab)
    f2 = net_source_rel(T + dT, T, g, metab)
    dlnF_dT = (np.log(f2) - np.log(f1)) / (2 * dT)
    return dlnF_dT * KB_EV * (T + 273.15) ** 2


# ---------------------------------------------------------------------------
# 2. Physical conversion of a flux into a partial-pressure excess
# ---------------------------------------------------------------------------

def phys_rel(T_C, T0_C, k600=ge.DEFAULT_K600_M_D):
    """[k(T0) K0(T0)] / [k(T) K0(T)]: how the same CO2 flux maps into a
    different pCO2 excess at T. k600 cancels; kept for clarity."""
    kk = lambda T: ge.k_co2(T, k600) * cs.constants(T)["K0"]
    return kk(T0_C) / kk(T_C)


# ---------------------------------------------------------------------------
# 3. Mineral surface
# ---------------------------------------------------------------------------

PZC25 = 8.1          # ferrihydrite, Dzombak & Morel 1990 (as in the Coulombic paper)
DPZC_DT = -0.02      # pH units per K. ASSUMED: oxide PZCs fall with T; a
                     # ferrihydrite-specific value still needs sourcing.


def pzc(T_C, pzc25=PZC25, dpzc_dT=DPZC_DT):
    return pzc25 + dpzc_dT * (np.asarray(T_C, dtype=float) - 25.0)


def psi0_mV(pH, T_C, pzc25=PZC25, dpzc_dT=DPZC_DT):
    """Nernstian surface potential (mV) at in-situ T."""
    slope = 1000 * R * (np.asarray(T_C) + 273.15) * np.log(10) / 96485.33212
    return slope * (pzc(T_C, pzc25, dpzc_dT) - np.asarray(pH))


def sorption_factor(D, z_eff):
    """Kd / Kd(at PZC) = 10**(-z_eff * D), D = PZC - pH.
    Phosphate z_eff ~ -1.5, molybdate -2 (Coulombic paper Methods)."""
    return 10.0 ** (-z_eff * np.asarray(D))


# ---------------------------------------------------------------------------
# 4. The chain for one site and one state
# ---------------------------------------------------------------------------

def site_state(alk_mgL, T_C, pco2_atm_ppm, dp_geo0, dp_met0, T0_C, g,
               metab=METAB, multiplicative=False, S0=None,
               pzc25=PZC25, dpzc_dT=DPZC_DT, pco2_atm0_ppm=None,
               switches=None):
    """Return pH, pCO2_w and Coulombic departure D for one state.

    dp_geo0, dp_met0 : geothermal and metabolic CO2 excess (ppm) at the
                       reference state (T0_C, reference atmosphere).
    switches         : optional dict to freeze channels for attribution:
                       'atm', 'phys', 'met', 'carb', 'pzc' -> bool (True = let
                       the channel follow T/atm; False = hold at reference).
    """
    sw = dict(atm=True, phys=True, met=True, carb=True, pzc=True)
    if switches:
        sw.update(switches)
    T_phys = T_C if sw["phys"] else T0_C
    T_met = T_C if sw["met"] else T0_C
    T_carb = T_C if sw["carb"] else T0_C
    T_pzc = T_C if sw["pzc"] else T0_C
    atm = pco2_atm_ppm if sw["atm"] else pco2_atm0_ppm

    if multiplicative:
        pco2w = S0 * atm
    else:
        pr = phys_rel(T_phys, T0_C)
        fr = net_source_rel(T_met, T0_C, g, metab)
        pco2w = atm + dp_geo0 * pr + dp_met0 * fr * pr
    pH = cs.pH_from_Alk_pCO2(alk_mgL, pco2w * 1e-6, T_carb)["pH"]
    D = pzc(T_pzc, pzc25, dpzc_dT) - pH
    return dict(pH=pH, pco2w=pco2w, D=D)


def shapley(fn, channels):
    """Exact Shapley attribution of fn(all on) - fn(all off) over channels.
    fn takes a dict {channel: bool}."""
    from itertools import permutations
    n = len(channels)
    contrib = {c: 0.0 for c in channels}
    cache = {}

    def val(on):
        key = tuple(sorted(on))
        if key not in cache:
            cache[key] = fn({c: (c in on) for c in channels})
        return cache[key]

    perms = list(permutations(channels))
    for order in perms:
        on = set()
        prev = val(on)
        for c in order:
            on.add(c)
            cur = val(on)
            contrib[c] += cur - prev
            prev = cur
    return {c: v / len(perms) for c, v in contrib.items()}
