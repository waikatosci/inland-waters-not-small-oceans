"""Published pH responses of labile and dissolved trace elements in fresh waters (round 1 + round 2 extractions).

Harmonises the extracted slopes (d log10 measure / d pH) into evidence classes, adds the
GLORICH river phosphate, our PHREEQC oxide runs and NICA-Donnan organic-matter runs, summarises them by element and evidence class
(median and interquartile range), counts diel directions, and writes
meta_table.csv, meta_pooled.csv, diel_votes.csv."""
import re, json
import numpy as np, pandas as pd

M = "../"
X = M + "extractions/"
Z = 1.959964

FRAC, CONC, AMD = "Field: labile or free fraction", "Field: labile or dissolved concentration", "Acid drainage: dilution-corrected"
DIEL, DIELT, EXP = "Diel cycle (pH and temperature)", "Diel cycle, temperature held constant", "pH experiment"
SOIL, SORB, FREEM = "Soil: dissolved at fixed total", "Sorption at fixed total (lab, natural oxides)", "Published free-ion models"
MHFO, MND = "Oxide model (HFO), dissolved", "Organic matter and oxide model (NICA-Donnan), free ion"

CATIONS = ["Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Cd", "Pb"]
OXY = ["P", "As", "As(III)", "V", "Mo", "W", "U"]
ALIAS = {"AsV": "As", "PO4": "P"}

r1 = pd.concat([pd.read_csv(X + f + "_slopes.csv").assign(file=f, rnd=1) for f in
                ["field_neutral", "field_acid", "lab_kd", "diel", "models_lit"]], ignore_index=True)
r2 = pd.concat([pd.read_csv(X + f + "_slopes.csv").assign(file=f, rnd=2) for f in
                ["field2", "diel2", "lab2", "models2"]], ignore_index=True)
DIEL2_STUDIES = set(r2[r2.file == "diel2"].study_id)
r1 = r1[~((r1.file == "diel") & r1.study_id.isin(DIEL2_STUDIES))]           # superseded by diel2
cat = pd.concat([r1, r2], ignore_index=True)
cat["slope"] = pd.to_numeric(cat.slope, errors="coerce")
cat["slope_se"] = pd.to_numeric(cat.slope_se, errors="coerce")
cat["n"] = pd.to_numeric(cat.n, errors="coerce")
cat = cat[~cat.notes.fillna("").str.contains("SUPERSEDED")]
cat["meas"] = cat.measure.fillna("").str.lower()
cat["el"] = cat.metal.astype(str).replace(ALIAS)
iii = cat.study_id.eq("Stachowicz2008") & cat.notes.fillna("").str.contains(r"As\(III\)|arsenite|As\(OH\)3", regex=True)
cat.loc[iii, "el"] = "As(III)"

rows = []


def add(element, cls, study, slope, se=np.nan, lo=np.nan, hi=np.nan, measure="", design="", pH="", site="", note=""):
    if np.isfinite(se) and not np.isfinite(lo):
        lo, hi = slope - Z * se, slope + Z * se
    rows.append(dict(element=element, cls=cls, study=study, site=site, slope=slope, se=se, lo=lo, hi=hi,
                     measure=measure, design=design, pH_range=pH, note=note))


def prange(r):
    return "" if pd.isna(r.pH_min) else f"{r.pH_min:g}-{r.pH_max:g}"


def classify(r):
    m, s, f = r.meas, r.study_id, r.file
    if f in ("field_neutral", "field_acid", "field2"):
        if "doc held" in m:
            return None                                   # reported in the SI table, not plotted (duplicates the fraction)
        if "so4" in m:
            return AMD
        if m.startswith("model_free_ion"):
            return FREEM                                  # WHAM free ion across natural waters (Han 2013)
        if m.startswith("labile_fraction") or ("free-ion fraction" in m and "measured" in m):
            return FRAC
        if m.startswith(("other", "model")):
            return None
        if r.design == "experimental_pH":
            return EXP
        if m.startswith(("c_dgt", "dissolved", "free_ion")):
            return None if s.startswith(("Xia", "Candeias")) else CONC
        return None
    if f in ("diel", "diel2"):
        if r.design == "experimental_pH":
            return EXP
        if r.design != "diel" or "24-h mean" in m:
            return None
        if "temperature held" in m or "day vs night" in m:
            return DIELT
        if "amplitude" in m or "evening" in m:
            return None                                   # extremes only: used for direction counts
        return DIEL
    if f in ("lab_kd", "lab2"):
        if "soil" in m and "fixed" in m:
            return SOIL
        if "fixed total" in m and "competitive" not in m:
            return SORB
        if s == "Korshin2015" and r.measure == "Kd":
            return "SORB_KD"
        if s == "Tessier1996" and m.startswith("kd"):
            return "SORB_KD"
        if m.startswith("model_free_ion") or "free-ion activity" in m:
            return FREEM
        return None
    return None


num = cat[cat.slope.notna()].copy()
num["cls"] = num.apply(classify, axis=1)
num = num[num.cls.notna() & num.el.isin(CATIONS + OXY)]
# one row per study, site, element, class and measure family (largest n) in the diel classes
num["site"] = num.water_body.fillna("").str[:40]
dd = num[num.cls.isin([DIEL, DIELT])].sort_values("n", ascending=False).drop_duplicates(["study_id", "site", "el", "cls"])
num = pd.concat([num[~num.cls.isin([DIEL, DIELT])], dd])
for _, r in num.iterrows():
    if r.cls == "SORB_KD":
        add(r.el, SORB, r.study_id, -r.slope, r.slope_se, measure="-(d log K / d pH)", design=r.design, pH=prange(r),
            site=r.site, note="binding-constant slope, sign reversed (dissolved at fixed total, strong-sorption limit)")
    else:
        add(r.el, r.cls, r.study_id, r.slope, r.slope_se, measure=r.measure, design=r.design, pH=prange(r), site=r.site,
            note=str(r.quality))

# GLORICH river phosphate
j = json.load(open("../../glorich/glorich_mixed_D.json")); C = np.array(j["cov"])
for alk, lab in [(200, "soft, 200 ueq/L"), (3000, "hard, 3000 ueq/L")]:
    x = np.array([1, np.log10(alk)])
    add("P", CONC, f"GLORICH rivers ({lab})", j["b0"] + j["b1"] * x[1], float(np.sqrt(x @ C @ x)),
        measure="dissolved inorganic P", design="within-site, 1002 rivers", note="mixed model, discharge and season held")

# models: oxide only (dissolved), and NICA-Donnan + oxide (free ion); pH 6.5-8.0
mo = pd.read_csv(M + "models/outputs/models_own_slopes.csv")
mo = mo[(mo.sorbent == "HFO_1mgFe/L") & (mo.quantity == "dissolved_fraction") & (mo.pH_range == "6.5-8.0")]
for _, r in mo.iterrows():
    add(r.element, MHFO, f"PHREEQC HFO, {r.water}, {'no carbonate' if 'noCO3' in r.model else 'with carbonate'}", r.slope,
        measure="dissolved", design="model, pH set by CO2 at fixed alkalinity", pH="6.5-8.0")
mn = pd.read_csv(M + "models/outputs/models_nicad_slopes.csv")
mn = mn[(mn.DOC_mgL > 0) & (mn.quantity == "free_ion") & (mn.pH_range == "6.5-8.0")]
for _, r in mn.iterrows():
    add(r.element, MND, f"NICA-Donnan, {r.water}, DOC {r.DOC_mgL:g}, {r.sorbent}", r.slope, measure="free ion",
        design="model, pH set by CO2 at fixed alkalinity", pH="6.5-8.0")

T = pd.DataFrame(rows)
T.to_csv(M + "synthesis/meta_table.csv", index=False)

# ---------- diel direction counts: one vote per study and element ----------
d = cat[cat.file.isin(["diel", "diel2"]) & (cat.design == "diel")].copy()
d = d[~d.study_id.isin(["Gammons2005a", "Gammons2005b"])]
d = d[~((d.study_id == "Kurz2013") & (d.file == "diel"))]          # round-1 Kurz As, V, Cr rows: not in the paper
d = d[~d.meas.str.contains("temperature held")]


def direction(r):
    if np.isfinite(r.slope):
        return -1 if r.slope < -0.05 else (1 if r.slope > 0.05 else 0)
    t = str(r.notes)
    if re.search(r"NULL", t): return 0
    m = re.search(r"DIRECTION:\s*(probably\s+)?(negative|positive)", t)
    if m: return -1 if m.group(2) == "negative" else 1
    if re.search(r"small and inconsistent|no clear diel", t, re.I): return 0
    return np.nan


d["dir"] = d.apply(direction, axis=1)
d = d.dropna(subset=["dir"])
vote = d.groupby(["study_id", "el"]).dir.mean().apply(lambda v: "neg" if v < -0.33 else ("pos" if v > 0.33 else "null"))
V = vote.reset_index().groupby(["el", "dir"]).study_id.nunique().unstack(fill_value=0)
V.index.name = "metal"
V.to_csv(M + "synthesis/diel_votes.csv")


# ---------- summaries: median and interquartile range ----------
# Slopes are summarised, not pooled: designs differ in how well they isolate pH, so a single pooled
# estimate would hide the spread that is the point of the synthesis.
# Omanovic2015 slopes reflect mixing with mine water between one clean and five affected sites (notes/field2.md):
# listed in Data S1, left out of the summaries.
NOSUM = {"Omanovic2015"}
GROUPS = {"Published measurements in natural water": [FRAC, CONC, AMD, DIELT, DIEL, EXP],
          "Published laboratory and model results": [SOIL, SORB, FREEM],
          MHFO: [MHFO], MND: [MND]}
LIMIT_SLOPE = {"Mn": -2, "Fe": -2, "Co": -2, "Ni": -2, "Cu": -2, "Zn": -2, "Cd": -2, "Pb": -2,
               "P": 1.5, "As": 1.5, "V": 1, "Mo": 2, "W": 2}
S_ = T[~T.study.isin(NOSUM)]
out = []


def summ(g):
    return dict(k=len(g), studies=g.study.nunique(), median=g.slope.median(),
                q25=g.slope.quantile(.25), q75=g.slope.quantile(.75), min=g.slope.min(), max=g.slope.max())


for e in CATIONS + OXY:
    for grp, classes in GROUPS.items():
        g = S_[(S_.element == e) & S_.cls.isin(classes)]
        if len(g):
            out.append(dict(level="element", element=e, cls=grp, **summ(g)))
DIV = [e for e in CATIONS if e != "Fe"]
for cls in [MHFO, MND, CONC, FRAC, AMD, DIELT, DIEL, EXP, SOIL, SORB, FREEM]:
    g = S_[(S_.cls == cls) & S_.element.isin(DIV)]
    if len(g):
        out.append(dict(level="class", element="divalent cations (no Fe)", cls=cls, **summ(g)))
P = pd.DataFrame(out)
lim = P.element.map(LIMIT_SLOPE).fillna(-2)
for c in ["median", "q25", "q75"]:
    P["phi_" + c] = P[c] / lim
P.loc[P.element.isin(["As(III)", "U"]), ["phi_median", "phi_q25", "phi_q75"]] = np.nan
P.to_csv(M + "synthesis/meta_pooled.csv", index=False)
if __name__ == "__main__":
    pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
    print(P.round(3).to_string(index=False))
    print(V)
    print(T.groupby(["cls"]).agg(rows=("slope", "size"), studies=("study", "nunique")))
