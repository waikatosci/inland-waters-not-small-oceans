"""Shapley attribution of the 2025 -> 2100 change in D and pH at Lake Karapiro with equal
activation energies for photosynthesis and respiration (unified.METAB_EQ). Same channels and
settings as run_projection.py. Writes results/attribution_eq.csv."""
import pandas as pd
import unified as u
import run_projection as rp
from waikato_sites import load_sites

site = load_sites().set_index("site").loc["Karāpiro"].copy(); site["site"] = "Karāpiro"
rows = []
for ssp in rp.SSPS:
    for target in ("D", "pH"):
        fn = lambda on: rp.annual(site, 2100, ssp, switches=on, metab=u.METAB_EQ)[target]
        a = u.shapley(fn, rp.CHANNELS)
        rows.append(dict(site="Karāpiro", ssp=ssp, target=target, **a, total=sum(a.values())))
out = pd.DataFrame(rows)
out.to_csv("results/attribution_eq.csv", index=False)
print(out.round(4).to_string(index=False))
ref = pd.read_csv("results/attribution.csv")
print(ref[(ref.site == "Karāpiro")].round(4).to_string(index=False))
