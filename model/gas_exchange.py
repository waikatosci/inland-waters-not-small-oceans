"""
gas_exchange.py
----------------
CO2 outgassing (air-water exchange) for a supersaturated river reach.

Schmidt number for CO2 in fresh water as a function of temperature:

    Sc_CO2(T) = 1742 - 91.24*T + 2.208*T^2 - 0.0219*T^3      (T in degC)

    Freshwater fit tabulated by Raymond et al. (2012), Limnol. Oceanogr. Fluids Environ. 2, 41-53,
    as coded in the GLEON LakeMetabolizer R package (getSchmidt.R). Valid 4-35 degC.

k600 -> k_gas scaling (Jahne et al. 1987):

    k_CO2 = k600 * (Sc_CO2 / 600) ** (-n)

    n = 2/3 for smooth, low-turbulence surfaces (wind speed below about 3.6 m/s at 10 m)
    n = 1/2 for rougher, wind- or wave-affected conditions

k600 is a supplied parameter (m/d). The default is a literature-typical value for a lowland regulated
river (Raymond et al. 2012 report stream and river medians mostly in the range 1-10 m/d). Because
evasion balances the internal supply at steady state, k600 sets the size of the CO2 excess for a given
source but does not enter the 1/S scaling of the atmospheric response.
"""
import numpy as np

DEFAULT_K600_M_D = 2.0  # m/d


def schmidt_co2(T_C):
    T_C = np.asarray(T_C, dtype=float)
    return 1742.0 - 91.24 * T_C + 2.208 * T_C ** 2 - 0.0219 * T_C ** 3


def k_co2(T_C, k600_m_d=DEFAULT_K600_M_D, n=2.0 / 3.0):
    """CO2 gas transfer velocity (m/d) at temperature T, scaled from k600."""
    sc = schmidt_co2(T_C)
    return k600_m_d * (sc / 600.0) ** (-n)


def outgassing_flux(co2aq_M, co2aq_eq_M, T_C, k600_m_d=DEFAULT_K600_M_D, n=2.0 / 3.0):
    """
    Air-water CO2 flux, mol/m2/d, POSITIVE = outgassing (river loses CO2
    to atmosphere), following F = k_CO2 * (CO2_aq - CO2_eq).
    co2aq_M, co2aq_eq_M in mol/L (converted to mol/m3 internally).
    """
    k = k_co2(T_C, k600_m_d, n)  # m/d
    return k * (np.asarray(co2aq_M) - np.asarray(co2aq_eq_M)) * 1000.0  # mol/m3/L->m3 = *1000


def dDIC_dt_gasexchange(co2aq_M, co2aq_eq_M, T_C, depth_m, k600_m_d=DEFAULT_K600_M_D, n=2.0 / 3.0):
    """
    Rate of change of DIC, in mol/L/d (matching the mol/L convention used
    throughout carbonate_system.py), in a well-mixed water column of given
    depth due to air-water exchange alone (negative = losing CO2 i.e.
    outgassing when supersaturated).

    flux is an AREAL rate (mol/m2/d); depth_m*1000 converts a depth in
    metres to a water volume in litres per m2 of river surface, so
    flux / (depth_m * 1000) is mol/L/d.
    """
    flux = outgassing_flux(co2aq_M, co2aq_eq_M, T_C, k600_m_d, n)  # mol/m2/d, +ve = out
    return -flux / (np.asarray(depth_m) * 1000.0)  # mol/L/d
