"""Draws the random sample of published slopes in Data S1 for checking against the source papers, and writes the
checking sheet. About 10% of the published slopes in each evidence class (rounded up) are drawn with a fixed seed, so
the sample is reproducible and fixed before any checking. Each sampled slope is linked to its extraction record, which
gives where in the source it was taken from and how it was obtained.
Writes check_sample_S1.csv (the sample) and, if a path is given, the checking workbook:
    python check_sample.py [checking_sheet.xlsx]"""
import glob, math, re, sys
import numpy as np, pandas as pd

SEED = 20261119
FRACTION = 0.10
EX = "../meta_analysis/extractions/"

s1 = pd.read_excel("Data_Set_S1_meta_analysis.xlsx", sheet_name="Slopes")
pub = s1[s1.source == "Published"].reset_index(names="s1_row")
pub["s1_row"] += 2                                      # spreadsheet row in the Slopes sheet (header is row 1)
E = pd.concat([pd.read_csv(f).assign(extraction_file=f.split("/")[-1]) for f in sorted(glob.glob(EX + "*_slopes.csv"))],
              ignore_index=True)
E["slope"] = pd.to_numeric(E.slope, errors="coerce")
E["el"] = E.metal.astype(str).replace({"AsV": "As", "PO4": "P"})   # as in meta_analysis/synthesis/meta_ph.py


def link(r):
    el = re.sub(r"\(.*", "", str(r.element))
    c = E[(E.study_id == r.study) & (E.el == el)]
    for sign in (1, -1):                                # Kd slopes enter Data S1 with the sign reversed
        m = c[np.isclose(c.slope * sign, r.slope_dlog10C_dpH, atol=6e-4)]
        if len(m) > 1 and isinstance(r.site, str) and r.site:
            m = m[m.water_body.fillna("").astype(str).str.startswith(r.site[:40])] if \
                m.water_body.fillna("").astype(str).str.startswith(r.site[:40]).any() else m
        if len(m):
            return m.index[0], sign
    raise ValueError(f"no extraction record for {r.study} {r.element} {r.slope_dlog10C_dpH}")


ix = [link(r) for r in pub.itertuples()]
pub["e"] = [i for i, _ in ix]
pub["sign_reversed"] = [s == -1 for _, s in ix]

rng = np.random.default_rng(SEED)
take = []
for cls, g in pub.groupby("evidence_class", sort=True):
    k = max(1, math.ceil(FRACTION * len(g)))
    take += list(rng.choice(g.index.values, size=k, replace=False))
S = pub.loc[sorted(take)].copy()
x = E.loc[S.e]
S["slope_source"] = x.slope_source.values
S["location_in_paper"] = x.location_in_paper.values
S["extraction_notes"] = x.notes.values
S["extraction_file"] = x.extraction_file.values
S["check_id"] = [f"C{i + 1:02d}" for i in range(len(S))]


def raw_rows(r):
    if r.slope_source in ("reported", "qualitative_only"):
        return ""
    f = r.extraction_file.replace("_slopes.csv", "_raw.csv")
    try:
        n = int((pd.read_csv(EX + f).study_id == r.study).sum())
    except FileNotFoundError:
        n = 0
    return f"{f} ({n} rows for this study)" if n else "none: inputs are the source values named in Notes"


S["input_values"] = [raw_rows(r) for r in S.itertuples()]
cols = ["check_id", "s1_row", "study", "authors", "year", "journal", "doi", "element", "evidence_class", "site", "measure",
        "pH_range", "slope_dlog10C_dpH", "sign_reversed", "slope_source", "location_in_paper", "extraction_file",
        "input_values", "extraction_notes"]
S[cols].to_csv("check_sample_S1.csv", index=False)
print(f"{len(S)} of {len(pub)} published slopes sampled (seed {SEED}); by evidence class:")
print(S.evidence_class.value_counts().to_string())
print("by how the slope was obtained:")
print(S.slope_source.value_counts().to_string())

if len(sys.argv) > 1:
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.worksheet.datavalidation import DataValidation
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    F = "Arial"
    thin = Side(style="thin", color="BFBFBF")
    head = Font(name=F, size=10, bold=True)
    body = Font(name=F, size=10)
    wrap = Alignment(wrap_text=True, vertical="top")
    fill_check = PatternFill("solid", fgColor="FFF8E5")   # columns the checker fills in

    ws = wb.active
    ws.title = "How to check"
    text = [
        ("Checking extracted slopes in Data S1 against their sources", True),
        ("", False),
        (f"This sheet lists {len(S)} of the {len(pub)} published slopes in Data S1, drawn at random within each evidence "
         f"class (about 10% of each, rounded up; seed {SEED}, data_sets/check_sample.py). Please check each one against "
         "the source paper and record the result in the shaded columns of the 'Check' sheet.", False),
        ("", False),
        ("For every row", True),
        ("1. Open the paper (DOI column) and go to the place given in 'Where in source'.", False),
        ("2. Confirm it is the right paper, element, site or experiment, and measure (dissolved, labile, free ion, Kd).", False),
        ("3. Check the value, following the rule for 'How obtained' below, and enter Yes, No or Unclear under 'Agrees'.", False),
        ("4. If No, give the value you find in the source and a short reason.", False),
        ("", False),
        ("Rules by 'How obtained'", True),
        ("reported: the slope or regression coefficient in Data S1 equals the value printed in the source, to its printed "
         "precision. Slopes on ln or other units were converted to log10 (divide an ln slope by 2.303); 'Notes' says "
         "where this was done. Where 'Sign reversed' is TRUE the source gives a partition coefficient (Kd), whose slope "
         "enters Data S1 with the opposite sign.", False),
        ("computed_from_table: the slope was fitted in this study to values tabulated in the source. Check about five of "
         "the input values in the table named under 'Input values' (in meta_analysis/extractions/) against the source "
         "table, and that the pH range matches; where 'Input values' says none, check the source values named in "
         "'Notes'. The fit itself is reproduced by the code and need not be redone.", False),
        ("computed_from_reported_points: the slope was fitted to points read from a figure or text. Check that the points "
         "agree with the source to within the reading tolerance stated in 'Notes' (typically 0.05 log units).", False),
        ("qualitative_only: check that the direction of change with pH matches the source.", False),
        ("", False),
        ("The 'Summary' sheet totals the results as you go. Agreement = Yes / (Yes + No); Unclear rows are reported "
         "separately.", False),
    ]
    for i, (t, b) in enumerate(text, 1):
        c = ws.cell(i, 1, t)
        c.font = Font(name=F, size=12 if i == 1 else 10, bold=b)
        c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.column_dimensions["A"].width = 110

    ck = wb.create_sheet("Check")
    hdr = [("Check ID", 9), ("Data S1 row", 9), ("Study", 18), ("Authors", 26), ("Year", 7), ("Journal", 30), ("DOI", 26),
           ("Element", 8), ("Evidence class", 24), ("Site or experiment", 30), ("Measure", 24), ("pH range", 11),
           ("Slope in Data S1", 10), ("Sign reversed", 9), ("How obtained", 18), ("Where in source", 22),
           ("Extraction file", 18), ("Input values", 24), ("Notes", 50),
           ("Found in source?", 10), ("Value in source", 16), ("Agrees", 10), ("If not, correct value", 14),
           ("Comment", 36), ("Checked by", 14), ("Date", 11)]
    first_check_col = 20
    for j, (h, w) in enumerate(hdr, 1):
        c = ck.cell(1, j, h)
        c.font = head
        c.alignment = Alignment(wrap_text=True, vertical="bottom")
        c.border = Border(bottom=Side(style="medium", color="808080"))
        if j >= first_check_col:
            c.fill = fill_check
        ck.column_dimensions[get_column_letter(j)].width = w
    for i, r in enumerate(S.itertuples(), 2):
        vals = [r.check_id, r.s1_row, r.study, r.authors, r.year, r.journal, r.doi, r.element, r.evidence_class,
                r.site if isinstance(r.site, str) else "", r.measure, r.pH_range, r.slope_dlog10C_dpH, r.sign_reversed,
                r.slope_source, r.location_in_paper, r.extraction_file, r.input_values,
                r.extraction_notes if isinstance(r.extraction_notes, str) else ""]
        for j, v in enumerate(vals, 1):
            v = "" if (isinstance(v, float) and np.isnan(v)) else v
            c = ck.cell(i, j, v)
            c.font = body
            c.alignment = wrap
            c.border = Border(bottom=thin)
        doi = str(r.doi)
        if doi.startswith("10."):
            ck.cell(i, 7).hyperlink = "https://doi.org/" + doi
            ck.cell(i, 7).font = Font(name=F, size=10, color="1F4E79", underline="single")
        for j in range(first_check_col, len(hdr) + 1):
            c = ck.cell(i, j)
            c.fill = fill_check
            c.font = body
            c.alignment = wrap
            c.border = Border(bottom=thin)
        ck.cell(i, 13).number_format = "0.000"
    n = len(S) + 1
    yn = DataValidation(type="list", formula1='"Yes,No"', allow_blank=True)
    ag = DataValidation(type="list", formula1='"Yes,No,Unclear"', allow_blank=True)
    ck.add_data_validation(yn)
    ck.add_data_validation(ag)
    yn.add(f"T2:T{n}")
    ag.add(f"V2:V{n}")
    ck.freeze_panes = "D2"
    ck.auto_filter.ref = f"A1:{get_column_letter(len(hdr))}{n}"

    sm = wb.create_sheet("Summary")
    sm.column_dimensions["A"].width = 34
    for col in "BCDEF":
        sm.column_dimensions[col].width = 12
    rows = [["", "Rows", "Checked", "Yes", "No", "Unclear"]]
    groups = ["All"] + sorted(S.slope_source.unique())
    for g in groups:
        crit = "" if g == "All" else f',Check!$O$2:$O${n},"{g}"'
        rows.append([g if g != "All" else "All sampled slopes",
                     f'=COUNTA(Check!$A$2:$A${n})' if g == "All" else f'=COUNTIF(Check!$O$2:$O${n},"{g}")',
                     f'=COUNTIFS(Check!$V$2:$V${n},"<>"{crit})',
                     f'=COUNTIFS(Check!$V$2:$V${n},"Yes"{crit})',
                     f'=COUNTIFS(Check!$V$2:$V${n},"No"{crit})',
                     f'=COUNTIFS(Check!$V$2:$V${n},"Unclear"{crit})'])
    for i, row in enumerate(rows, 1):
        for j, v in enumerate(row, 1):
            c = sm.cell(i, j, v)
            c.font = head if i == 1 or j == 1 else body
    k = len(rows) + 2
    sm.cell(k, 1, "Agreement, Yes / (Yes + No)").font = head
    sm.cell(k, 2, '=IF(D2+E2=0,"",D2/(D2+E2))').number_format = "0%"
    sm.cell(k + 1, 1, "Rows still to check").font = body
    sm.cell(k + 1, 2, "=B2-C2")
    sm.cell(k + 3, 1, "Sampled from").font = body
    sm.cell(k + 3, 2, f"{len(S)} of {len(pub)} published slopes in Data S1, {S.study.nunique()} of "
                      f"{pub.study.nunique()} studies; seed {SEED}").font = body
    for w in (ws, ck, sm):
        w.page_setup.orientation = "landscape"
        w.page_setup.fitToWidth, w.page_setup.fitToHeight = 1, 0
        w.sheet_properties.pageSetUpPr.fitToPage = True
    ck.print_title_rows = "1:1"
    wb.save(sys.argv[1])
    print("wrote", sys.argv[1])
