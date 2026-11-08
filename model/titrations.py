"""Parse ECP Gran alkalinity titrations (2024-25 and 2026 workbooks) into long form,
and compute Gran alkalinity for each titration.

Layout: each sheet holds one sampling date. Titrations sit in 5-column blocks
(Sample vol (ml) | ml titrant | pH | sec | degC) under a 'Sample ID:' row, with the
site name two rows above. Titrant concentration is read from each sheet's "Acid conc (M)" cell
(0.0901 M in the 2026 workbook, confirmed Oct 2026); sheets without one use N_ACID.
"""
import re
import numpy as np
import pandas as pd

FILES = ["data/ECP_sampling_Alkalinity.xlsx", "data/ECP_sampling_Alkalinity2026.xlsx"]
SKIP = ("run_records", "BLANK", "BOAT CAL")
N_ACID = 0.1  # eq/L, default where a sheet records no acid concentration
# 2026 workbook: 0.0901 M throughout (one sheet, 24 Mar 26, leaves the cell blank)
ACID_DEFAULT = {"data/ECP_sampling_Alkalinity2026.xlsx": 0.0901}


def sheet_acid(d, default=N_ACID):
    for i in range(min(6, d.shape[0])):
        row = d.iloc[i].tolist()
        for j, v in enumerate(row[:-1]):
            if isinstance(v, str) and "acid conc" in v.lower():
                n = pd.to_numeric(row[j + 1], errors="coerce")
                if np.isfinite(n) and 0.01 < n < 1:
                    return float(n)
    return default


def parse():
    curves, meta = [], []
    for f in FILES:
        x = pd.ExcelFile(f)
        for sh in x.sheet_names:
            if sh in SKIP:
                continue
            d = pd.read_excel(x, sheet_name=sh, header=None)
            nacid = sheet_acid(d, ACID_DEFAULT.get(f, N_ACID))
            col0 = [str(v) for v in d.iloc[:, 0].tolist()]
            id_rows = [i for i, v in enumerate(col0) if v.strip().startswith("Sample ID")]
            for k, r in enumerate(id_rows):
                end = id_rows[k + 1] - 2 if k + 1 < len(id_rows) else d.shape[0]
                for c in range(1, d.shape[1] - 3, 5):
                    sid = d.iloc[r, c]
                    if not isinstance(sid, str) or not sid.strip():
                        continue
                    site = None
                    for back in (1, 2):
                        v = d.iloc[r - back, c] if r - back >= 0 else None
                        if isinstance(v, str) and v.strip():
                            site = v.strip(); break
                    vol = pd.to_numeric(d.iloc[r + 2, c], errors="coerce")
                    ml = pd.to_numeric(d.iloc[r + 2:end, c + 1], errors="coerce")
                    ph = pd.to_numeric(d.iloc[r + 2:end, c + 2], errors="coerce")
                    tc = pd.to_numeric(d.iloc[r + 2:end, c + 4], errors="coerce")
                    ok = ml.notna() & ph.notna()
                    if ok.sum() < 20 or not np.isfinite(vol) or vol <= 0:
                        continue
                    key = f"{sh}|{sid.strip()}"
                    cur = pd.DataFrame(dict(key=key, ml=ml[ok].values, pH=ph[ok].values,
                                            T=tc[ok].values))
                    cur["acid_meqL"] = cur.ml * nacid / (vol / 1000.0)
                    curves.append(cur)
                    meta.append(dict(key=key, sheet=sh, sample=sid.strip(), site_label=site,
                                     vol_ml=vol, acid_M=nacid, T=float(np.nanmedian(cur["T"])),
                                     pH0=float(cur.pH.iloc[0]), min_pH=float(cur.pH.min())))
    return pd.concat(curves, ignore_index=True), pd.DataFrame(meta)


def site_code(sample):
    m = re.match(r"(S\d+[A-Za-z]?|MS\d+[a-z]?)", sample)
    return m.group(1) if m else "other"


def gran_alk(cur, vol_ml, n_acid=N_ACID, lo=3.3, hi=4.3):
    """Gran function F1 = (V0+v)*10^-pH on the acid side; x-intercept -> Ve."""
    s = cur[(cur.pH < hi) & (cur.pH > lo)]
    if len(s) < 5:
        return np.nan, np.nan
    F = (vol_ml + s.ml) * 10 ** (-s.pH)
    b, a = np.polyfit(s.ml, F, 1)
    ve = -a / b
    r2 = np.corrcoef(s.ml, F)[0, 1] ** 2
    return ve * n_acid / (vol_ml / 1000.0), r2   # meq/L


if __name__ == "__main__":
    cur, meta = parse()
    meta["site"] = meta["sample"].map(site_code)
    alk = [gran_alk(cur[cur.key == k], v, n) for k, v, n in zip(meta.key, meta.vol_ml, meta.acid_M)]
    meta["alk_meqL"] = [a for a, _ in alk]
    meta["gran_r2"] = [r for _, r in alk]
    cur.to_csv("results/titration_curves.csv", index=False)
    meta.to_csv("results/titration_meta.csv", index=False)
    print(len(meta), "titrations;", meta.sheet.nunique(), "sampling dates")
    print(meta.site.value_counts().to_string())
    print(meta[["alk_meqL", "gran_r2", "min_pH", "pH0", "T"]].describe().round(3).to_string())
