"""
atmosphere.py
-------------
Rising atmospheric CO2 background used as the "equilibrium target" that
supersaturated river CO2 outgasses towards. Anchored to the CMIP6/IPCC
AR6 Shared Socioeconomic Pathway concentration trajectories from

    Meinshausen, M. et al. (2020) The shared socio-economic pathway (SSP)
    greenhouse gas concentrations and their extensions to 2500.
    Geoscientific Model Development, 13, 3571-3605.
    (Table 5: global-mean CO2, ppm)

Anchor points used (global-mean CO2, ppm), from Table 5 of Meinshausen
et al. (2020):

    scenario     2025    2050    2100
    SSP1-2.6     427.7   469.3   445.6   (peaks in the 2060s, then eases back,
                                          staying above the 2025 level)
    SSP2-4.5     429.0   506.9   602.8
    SSP3-7.0     432.3   540.6   867.2
    SSP5-8.5     431.8   562.8  1135.2

The 2025 values differ slightly between scenarios, although they share the same
recent past; this is an artefact of how each SSP marker run was constructed. The
values are global means, which serve as the gas-exchange equilibrium target
because CO2 (unlike O2) is well mixed, so no Waikato-specific correction is needed
at this precision.

SSP1-2.6 is non-monotonic: concentration peaks in the 2060s and then eases back as
net-negative emissions take effect, ending the century above the 2025 level. The
quadratic fit through the anchors does not capture this downturn precisely, so the
SSP1-2.6 trajectory is indicative only after about 2080; the full Meinshausen time
series would be needed to resolve the late-century tail.
"""

import numpy as np

ANCHORS = {
    # year: ppm -- Meinshausen et al. (2020) Table 5, global-mean CO2
    "SSP1-2.6": {2025: 427.7, 2050: 469.3, 2100: 445.6},
    "SSP2-4.5": {2025: 429.0, 2050: 506.9, 2100: 602.8},
    "SSP3-7.0": {2025: 432.3, 2050: 540.6, 2100: 867.2},
    "SSP5-8.5": {2025: 431.8, 2050: 562.8, 2100: 1135.2},
}


def _quadratic_from_anchors(years, ppms):
    """Fit a quadratic through 3 anchor points for smooth interpolation."""
    coeffs = np.polyfit(years, ppms, 2)
    return np.poly1d(coeffs)


_FITS = {
    scen: _quadratic_from_anchors(list(a.keys()), list(a.values()))
    for scen, a in ANCHORS.items()
}


def co2_atm_ppm(year, scenario="SSP2-4.5"):
    """
    Atmospheric CO2 (ppm) in a given calendar year under a given SSP
    scenario. `year` can be a scalar or array. Quadratic interpolation
    through the 2025/2050/2100 anchors (exact at the anchors); clipped
    to >= 280 ppm (pre-industrial floor) as a sanity bound.
    """
    if scenario not in _FITS:
        raise ValueError(f"scenario must be one of {list(_FITS)}")
    year = np.asarray(year, dtype=float)
    ppm = _FITS[scenario](year)
    return np.clip(ppm, 280.0, None)


def all_scenarios(year):
    return {scen: float(co2_atm_ppm(year, scen)) for scen in ANCHORS}


if __name__ == "__main__":
    for yr in (2025, 2030, 2050, 2075, 2100):
        row = " ".join(f"{s}: {co2_atm_ppm(yr, s):6.1f}" for s in ANCHORS)
        print(f"{yr}  {row}")
