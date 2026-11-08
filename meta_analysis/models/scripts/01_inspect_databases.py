"""Inventory of hydrous ferric oxide (HFO, Dzombak & Morel 1990 diffuse-layer model) surface constants.

Parses the SURFACE_SPECIES blocks of every database we have and lists, for each sorbate of interest,
the reactions and log K values present. Nothing is added or edited here.

Databases inspected
  databases_downloaded/  : USGS phreeqc3 GitHub master copies (phreeqc.dat, wateq4f.dat, minteq.v4.dat,
                           minteq.dat, llnl.dat), fetched 2026-10-07 (md5 in outputs/environment.txt)
  phreeqpython bundle    : phreeqc.dat, wateq4f_PWN.dat, vitens.dat, pitzer.dat shipped inside phreeqpython 1.6.2

Writes outputs/hfo_constants_inventory.csv (one row per database x reaction)
       outputs/hfo_constants_presence.csv  (sorbate x database presence matrix)
Run: python3 scripts/01_inspect_databases.py
"""
import os
import re

import pandas as pd
import phreeqpython

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "outputs")
DL = os.path.join(ROOT, "databases_downloaded")
BUNDLE = os.path.join(os.path.dirname(phreeqpython.__file__), "database")

DATABASES = [
    ("USGS phreeqc.dat", os.path.join(DL, "phreeqc.dat")),
    ("USGS wateq4f.dat", os.path.join(DL, "wateq4f.dat")),
    ("USGS minteq.v4.dat", os.path.join(DL, "minteq.v4.dat")),
    ("USGS minteq.dat", os.path.join(DL, "minteq.dat")),
    ("USGS llnl.dat", os.path.join(DL, "llnl.dat")),
    ("phreeqpython phreeqc.dat", os.path.join(BUNDLE, "phreeqc.dat")),
    ("phreeqpython wateq4f_PWN.dat", os.path.join(BUNDLE, "wateq4f_PWN.dat")),
    ("phreeqpython vitens.dat", os.path.join(BUNDLE, "vitens.dat")),
]

# sorbate label -> regex matched against the left-hand side of the reaction (aqueous reactant)
SORBATES = {
    "Mn(II)": r"\bMn\+2\b",
    "Ni": r"\bNi\+2\b",
    "Co(II)": r"\bCo\+2\b",
    "Zn": r"\bZn\+2\b",
    "Cd": r"\bCd\+2\b",
    "Cu(II)": r"\bCu\+2\b",
    "Pb": r"\bPb\+2\b",
    "Ca": r"\bCa\+2\b",
    "Mg": r"\bMg\+2\b",
    "phosphate": r"PO4-3|H2PO4|HPO4",
    "molybdate": r"MoO4",
    "arsenate As(V)": r"AsO4|H3AsO4",
    "vanadate V(V)": r"VO4|VO2\+|H2VO4",
    "selenite Se(IV)": r"SeO3",
    "sulfate": r"SO4-2",
    "carbonate": r"CO3-2|HCO3",
    # extra elements reported for completeness
    "arsenite As(III)": r"H3AsO3",
    "selenate Se(VI)": r"SeO4-2",
    "chromate Cr(VI)": r"CrO4",
    "antimonate Sb(V)": r"Sb\(OH\)6",
    "uranyl U(VI)": r"UO2\+2",
    "silicate": r"H4SiO4",
    "borate": r"H3BO3",
    "fluoride": r"\bF-",
    "Fe(II)": r"\bFe\+2\b",
    "Sr": r"\bSr\+2\b",
    "Ba": r"\bBa\+2\b",
}
REQUESTED = ["Mn(II)", "Ni", "Co(II)", "Zn", "Cd", "Cu(II)", "Pb", "Ca", "Mg", "phosphate", "molybdate",
             "arsenate As(V)", "vanadate V(V)", "selenite Se(IV)", "sulfate", "carbonate"]

KEYWORDS = re.compile(r"^\s*(SOLUTION_MASTER_SPECIES|SOLUTION_SPECIES|PHASES|EXCHANGE_MASTER_SPECIES|EXCHANGE_SPECIES|"
                      r"SURFACE_MASTER_SPECIES|SURFACE_SPECIES|RATES|END|MEAN_GAMMAS|LLNL_AQUEOUS_MODEL_PARAMETERS|"
                      r"NAMED_EXPRESSIONS|CALCULATE_VALUES|ISOTOPES|ISOTOPE_RATIOS|ISOTOPE_ALPHAS|PITZER|SIT|"
                      r"GAS_BINARY_PARAMETERS)\b")


def surface_reactions(path):
    """Return list of (reaction, log_k, comment_context) for Hfo_ surface species."""
    with open(path, encoding="latin-1") as f:
        lines = f.read().splitlines()
    out = []
    in_block = False
    last_comment = ""
    cur = None
    for raw in lines:
        line = raw.split("#")[0].rstrip()
        comment = raw.split("#", 1)[1].strip() if "#" in raw else ""
        m = KEYWORDS.match(raw)
        if m:
            in_block = m.group(1) == "SURFACE_SPECIES"
            continue
        if not in_block:
            continue
        if not line.strip():
            if comment:
                last_comment = comment
            continue
        # semicolon-joined forms such as "Hfo_wOH + H4SiO4 = Hfo_wH3SiO4 + H2O ; log_K 4.28"
        parts = [p.strip() for p in line.split(";") if p.strip()]
        for p in parts:
            if "=" in p and "Hfo_" in p and not p.lower().lstrip("-").startswith("log_k"):
                cur = {"reaction": re.sub(r"\s+", " ", p), "log_k": None, "context": last_comment,
                       "inline_comment": comment}
                out.append(cur)
            else:
                mk = re.match(r"-?log_k\s+([-+0-9.eE]+)", p, flags=re.I)
                if mk and cur is not None and cur["log_k"] is None:
                    cur["log_k"] = float(mk.group(1))
        if comment and not any("=" in p for p in parts):
            last_comment = comment
    return out


def main():
    rows = []
    for dbname, path in DATABASES:
        if not os.path.exists(path):
            continue
        for r in surface_reactions(path):
            lhs = r["reaction"].split("=")[0]
            if lhs.replace(" ", "") in ("Hfo_sOH", "Hfo_wOH", "Hfo_sOH+H+", "Hfo_wOH+H+"):
                group = "acid-base"
            else:
                group = ";".join(k for k, rx in SORBATES.items() if re.search(rx, lhs)) or "other"
            site = "strong (Hfo_s)" if "Hfo_s" in lhs else "weak (Hfo_w)"
            rows.append(dict(database=dbname, sorbate=group, site=site, reaction=r["reaction"], log_k=r["log_k"],
                             source_comment_in_database=(r["inline_comment"] or r["context"])))
    inv = pd.DataFrame(rows)
    inv = inv[inv.log_k.notna() | inv.reaction.str.contains("=")]
    inv.to_csv(os.path.join(OUT, "hfo_constants_inventory.csv"), index=False)

    pres = []
    for s in list(SORBATES) + ["acid-base"]:
        rec = {"sorbate": s, "requested": s in REQUESTED or s == "acid-base"}
        for dbname, _ in DATABASES:
            sub = inv[(inv.database == dbname) & (inv.sorbate.str.split(";").apply(lambda L: s in L))]
            if len(sub) == 0:
                rec[dbname] = "-"
            else:
                rec[dbname] = " | ".join(f"{'s' if 'strong' in st else 'w'}: {lk:g}" for st, lk in zip(sub.site, sub.log_k))
        pres.append(rec)
    pres = pd.DataFrame(pres)
    pres.to_csv(os.path.join(OUT, "hfo_constants_presence.csv"), index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 60)
    for dbname, _ in DATABASES:
        print(f"\n== {dbname}")
        for _, r in pres[pres.requested].iterrows():
            print(f"   {r['sorbate']:16s} {r[dbname]}")
    print("\nwrote outputs/hfo_constants_inventory.csv and outputs/hfo_constants_presence.csv")


if __name__ == "__main__":
    main()
