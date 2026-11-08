# Extraction schema for the pH-lability meta-analysis (read this before extracting)

Purpose: a meta-analysis of how labile (DGT) and dissolved trace-metal concentrations in fresh waters respond to
pH, to separate the pH effect from co-varying season, temperature, DOC and flow. Every number must be traceable
to a specific table, figure caption, text passage or supplementary file of a real, citable paper. Never estimate,
recall from memory or invent values. If a value is not stated, leave the cell empty.

Elements of interest:
- Cations (Chelex DGT; expected to become MORE labile/dissolved as pH falls, i.e. negative slope on pH):
  Mn, Ni, Co, Zn, Cd, Cu, Pb, Fe, Al (others welcome: rare earths, Ba, Sr as controls).
- Oxyanions (DGT with ferrihydrite, Metsorb/TiO2 or zirconium-oxide binding layers; expected to become MORE
  dissolved/labile as pH RISES, i.e. positive slope on pH, because oxide surfaces lose positive charge):
  Mo (molybdate), P (phosphate), V (vanadate), W (tungstate), As (arsenate/arsenite), Sb, Se (selenate/selenite),
  U (uranyl carbonates), Cr(VI) (chromate), Si (silicate, control).
Record both groups. The sign convention is always d log10(measure)/d pH; do not flip signs.

## File 1: slopes  ->  extractions/<your_tag>_slopes.csv
One row per study x water body (or study-wide pooled analysis) x metal x measure. Columns, in this order:
study_id        short key, e.g. Gimpel2003
authors         first author et al.
year
title
journal
doi             (or URL if no DOI)
access          open | SI_only | abstract_only | not_accessed
water_body      name(s)
country
water_type      stream | river | lake | reservoir | wetland | mine_affected | estuary_fresh | lab_natural_water | lab_synthetic | other
design          spatial_gradient | temporal_single_site | diel | experimental_pH | model_comparison | review
metal
measure         C_DGT | labile_fraction (C_DGT / C_filtered) | dissolved (filtered) | Kd | free_ion | model_free_ion | model_labile | model_dissolved | other (say what)
binding_layer   Chelex | ferrihydrite | Metsorb | ZrO2 | other | (empty if not DGT)
pH_min
pH_max
n               number of observations in the relationship
slope           d log10(measure) / d pH  (negative = higher at lower pH). Convert if the paper used ln or another unit and say so in notes.
slope_se        standard error if given or computable; else empty
ci_low, ci_high if given
r               correlation coefficient (Pearson or Spearman; say which in notes)
p_value
slope_source    reported | computed_from_table | computed_by_you_from_reported_points | qualitative_only
DOC_mgL         range or mean
Ca_mgL_or_hardness
alkalinity
temp_C_range
location_in_paper   e.g. "Table 2", "Fig. 4 caption", "SI Table S3", "p. 1945 text"
notes           confounders (season, flow, redox, Fe/Al precipitation), method caveats (resin binding at low pH,
                capacity, diffusion layer), how the slope was obtained, anything a reviewer would question
quality         A = computed by you from tabulated per-sample data | B = regression reported by authors |
                C = only direction or correlation reported

## File 2: raw tabulated data  ->  extractions/<your_tag>_raw.csv
Whenever a paper (or its free supplementary information) tabulates per-sample or per-site values of pH together
with DGT-labile and/or filtered concentrations, copy those rows exactly. Columns:
study_id, site_or_sample, date_or_period, water_type, pH, metal, C_DGT, C_filtered, C_total, unit (ug/L, nM ...),
DOC_mgL, Ca_mgL, temp_C, location_in_paper, notes

## File 3: candidates  ->  extractions/<your_tag>_candidates.csv
Every relevant paper you identified, including those you could not open, so the user can supply PDFs.
study_id, authors, year, title, journal, doi, why_relevant, access, data_expected (e.g. "Table 3 has pH and DGT Ni, Zn at 12 lakes")

## File 4: notes  ->  notes/<your_tag>.md
Short: search strategy, what you could and could not access, key qualitative findings, method caveats.

Tips: ACS, Elsevier and Springer supplementary information is often free even when the article is paywalled.
PubMed Central, author repositories, theses (university repositories) and government reports (USGS, EPA, NIVA,
CEH) often have full text with tables. Record negative and null results too (a slope near zero is data).
