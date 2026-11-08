"""Implementation check: reproduce USGS PHREEQC Example 8 (Zn sorption on hydrous ferric oxide).

Reference input and output were fetched on 2026-10-07 from
https://raw.githubusercontent.com/usgs-coupled/phreeqc3/master/examples/ex8, ex8.sel
(the selected-output file shipped with the PHREEQC distribution as the expected result).
Example 8 is the Dzombak and Morel (1990) diffuse-layer model with 0.09 g HFO, 5e-6 mol strong
and 2e-4 mol weak sites, 600 m2/g, 0.1 M NaNO3, total Zn 1e-7 and 1e-4 molal.

This checks that our engines + wrapper reproduce the published reference numbers. It is a
check of the machinery, not an independent test of the Dzombak and Morel constants.
Writes validation/usgs_ex8_comparison.csv
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engine import Engine, db_path, ROOT  # noqa: E402

REF = os.path.join(ROOT, "validation", "usgs_ex8", "ex8.sel")
OUT = os.path.join(ROOT, "validation", "usgs_ex8_comparison.csv")

HEAD = """
SURFACE_SPECIES
     Hfo_sOH  + H+ = Hfo_sOH2+
     log_k  7.18
     Hfo_sOH = Hfo_sO- + H+
     log_k  -8.82
     Hfo_sOH + Zn+2 = Hfo_sOZn+ + H+
     log_k  0.66
     Hfo_wOH  + H+ = Hfo_wOH2+
     log_k  7.18
     Hfo_wOH = Hfo_wO- + H+
     log_k  -8.82
     Hfo_wOH + Zn+2 = Hfo_wOZn+ + H+
     log_k  -2.32
SURFACE 1
     Hfo_sOH        5e-6    600.    0.09
     Hfo_wOH        2e-4
END
SOLUTION 1
     -units  mmol/kgw
     pH      8.0
     Zn      0.0001
     Na      100.    charge
     N(5)    100.
SOLUTION 2
     -units  mmol/kgw
     pH      8.0
     Zn      0.1
     Na      100.    charge
     N(5)    100.
PHASES
     Fix_H+
     H+ = H+
     log_k  0.0
END
SELECTED_OUTPUT
     -reset false
     -pH true
     -solution true
     -state true
     -molalities Zn+2 Hfo_wOZn+ Hfo_sOZn+
END
"""


def build():
    s = HEAD
    for sol in (1, 2):
        for ph in np.arange(5.0, 8.0001, 0.25):
            s += f"USE solution {sol}; USE surface 1\nEQUILIBRIUM_PHASES 1\n   Fix_H+ {-ph:.2f} NaOH 10.0\nEND\n"
    return s


def read_ref():
    df = pd.read_csv(REF, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    df = df[[c for c in df.columns if c and not c.startswith("Unnamed")]]
    return df


def main():
    ref = read_ref()
    ref = ref[ref["state"].str.strip() == "react"].reset_index(drop=True)
    results = []
    engines = [("386", db_path("phreeqc.dat"))]
    # 3.7.3 cannot read the 2024 phreeqc.dat; use the copy bundled with phreeqpython 1.6.2
    import phreeqpython
    engines.append(("373", os.path.join(os.path.dirname(phreeqpython.__file__), "database", "phreeqc.dat")))
    for which, db in engines:
        eng = Engine(db, which)
        rows = eng.run(build())
        rows = [r for r in rows if str(r.get("state", "")).strip() == "react"]
        assert len(rows) == len(ref), (len(rows), len(ref))
        for r, (_, q) in zip(rows, ref.iterrows()):
            for col in ("m_Zn+2", "m_Hfo_wOZn+", "m_Hfo_sOZn+"):
                results.append({
                    "engine": eng.version, "database": os.path.basename(db), "Zn_total_molal": "1e-7" if int(q["soln"]) == 1 else "1e-4",
                    "pH": r["pH"], "quantity": col, "this_study": r[col], "usgs_reference": float(q[col]),
                    "rel_diff": r[col] / float(q[col]) - 1.0,
                })
    df = pd.DataFrame(results)
    df.to_csv(OUT, index=False)
    for eng, g in df.groupby("engine"):
        print(f"{eng}: max |relative difference| vs USGS ex8.sel = {g['rel_diff'].abs().max():.2e} over {len(g)} values")
    # pH of 50 % sorption for each Zn loading, for the README
    for eng, g in df[df.engine == df.engine.iloc[0]].groupby("Zn_total_molal"):
        piv = g.pivot_table(index="pH", columns="quantity", values="this_study")
        tot = piv.sum(axis=1)
        fs = (piv["m_Hfo_wOZn+"] + piv["m_Hfo_sOZn+"]) / tot
        ph50 = np.interp(0.5, fs.values, fs.index.values)
        print(f"  Zn total {eng}: pH50 (50 % sorbed) = {ph50:.2f}")


if __name__ == "__main__":
    main()
