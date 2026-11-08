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

Anchor points used (global-mean CO2, ppm) -- READ VERBATIM FROM TABLE 5 OF
THE MEINSHAUSEN ET AL. (2020) PAPER ITSELF (the first version of this file,
from an earlier session, had these numbers wrong for 6 of the 12 anchors --
right at the two most commonly quoted endpoints, SSP1-2.6 2050 and SSP5-8.5
2100, wrong everywhere else, by up to 23% at SSP3-7.0 2100. Re-verified
against the primary source before use here):

    scenario     2025    2050    2100
    SSP1-2.6     427.7   469.3   445.6   (peaks ~2060s, eases back -- but
                                          stays above today, not below it)
    SSP2-4.5     429.0   506.9   602.8
    SSP3-7.0     432.3   540.6   867.2
    SSP5-8.5     431.8   562.8  1135.2

2025 values differ slightly between scenarios even though they're meant to
share the same near-term past -- a small artefact of how each SSP marker
run was constructed, not a Waikato-relevant effect. These are GLOBAL-MEAN
atmospheric CO2 -- fine as the gas-exchange equilibrium target since CO2
(unlike O2) is well-mixed enough globally that a Waikato-specific
correction is not warranted at this level of precision.

NOTE: SSP1-2.6 is non-monotonic (concentration peaks ~2060s then eases
back slightly as net-negative emissions kick in, ending the century still
well above 2025 levels, not below them) -- the quadratic anchor fit below
will not capture the down-turn shape precisely past 2100; treat this
trajectory as indicative only beyond ~2080, and swap in the full
Meinshausen ppm time series (freely available) if the multi-decade tail
matters for the write-up.
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
