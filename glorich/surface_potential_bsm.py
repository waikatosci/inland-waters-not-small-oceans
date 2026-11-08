"""Surface potential of ferrihydrite against pH: Nernst relation versus a one-pK Basic Stern model
(Hiemstra and Van Riemsdijk). Surface reaction FeOH(-1/2) + H+ = FeOH2(+1/2), log K = PZC = 8.1.
Site density Ns, Stern capacitance C, Gouy-Chapman diffuse layer at ionic strength I (1:1 electrolyte),
charge balance sigma0 + sigma_d = 0 (no ion pairs). Reports the slope d(psi0)/dpH and the departure
from Nernst over the freshwater range, for a range of Ns, C and I."""
import numpy as np
from scipy.optimize import brentq

F, R, T, EPS = 96485.0, 8.314, 298.15, 78.5 * 8.854e-12
PZC = 8.1


def psi0(pH, Ns_nm2=6.0, C=1.15, I=1.5e-3):
    Ns = Ns_nm2 * 1e18 / 6.022e23 * F            # C/m2 per unit charge fraction
    c = I * 1000.0                                 # mol/m3

    def resid(p0):
        theta = 1.0 / (1.0 + 10 ** (pH - PZC) * np.exp(F * p0 / (R * T)))   # fraction protonated
        s0 = Ns * (theta - 0.5)
        pd = p0 - s0 / C
        sd = -np.sqrt(8 * EPS * R * T * c) * np.sinh(F * pd / (2 * R * T))
        return s0 + sd
    return brentq(resid, -0.3, 0.3)


pHs = np.linspace(6.0, 8.1, 22)
nernst = 59.16e-3 * (PZC - pHs)
print("Ns(nm-2)  C(F/m2)  I(mol/L)  slope 6.5-7.8 (mV/pH)  fraction of Nernst  max |dev| 6.5-8.1 (mV)")
for Ns in (3.0, 6.0, 8.0):
    for C in (0.9, 1.15, 1.5):
        for I in (1.5e-3, 1e-2):
            v = np.array([psi0(p, Ns, C, I) for p in pHs])
            m = (pHs >= 6.5) & (pHs <= 7.8)
            sl = -np.polyfit(pHs[m], v[m], 1)[0] * 1e3
            dev = np.max(np.abs(v - nernst)[pHs >= 6.5]) * 1e3
            print(f"{Ns:7.1f}  {C:7.2f}  {I:8.4f}  {sl:21.1f}  {sl/59.16:18.2f}  {dev:10.1f}")
