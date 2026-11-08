# Extraction protocol for the pH-lability meta-analysis

Purpose: a meta-analysis of how labile (DGT) and dissolved trace-metal concentrations in fresh waters respond to
pH, separating the pH effect from co-varying season, temperature, DOC and flow. Every number is traceable to a
specific table, figure caption, text passage or supplementary file of a citable paper. Values were taken only
from the paper or its supplementary material; where a value was not stated, the cell is empty.

Elements covered:
- Cations (Chelex DGT; expected to become more labile or dissolved as pH falls, i.e. a negative slope on pH):
  Mn, Ni, Co, Zn, Cd, Cu, Pb, Fe, Al (with rare earths, Ba and Sr as controls where reported).
- Oxyanions (DGT with ferrihydrite, Metsorb/TiO2 or zirconium-oxide binding layers; expected to become more
  dissolved or labile as pH rises, i.e. a positive slope on pH, because oxide surfaces lose positive charge):
  Mo (molybdate), P (phosphate), V (vanadate), W (tungstate), As (arsenate/arsenite), Sb, Se (selenate/selenite),
  U (uranyl carbonates), Cr(VI) (chromate), Si (silicate, control).
Both groups were recorded. The sign convention is always d log10(measure)/d pH.

Extraction was done in topic batches, each with a short tag (for example `diel`, `field_neutral`, `lab_kd`);
a second or third pass on the same topic carries a number (`diel2`, `redox3`). Each batch has the files below.

## File 1: slopes  ->  extractions/<tag>_slopes.csv
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
measure         C_DGT | labile_fraction (C_DGT / C_filtered) | dissolved (filtered) | Kd | free_ion | model_free_ion | model_labile | model_dissolved | other (described)
binding_layer   Chelex | ferrihydrite | Metsorb | ZrO2 | other | (empty if not DGT)
pH_min
pH_max
n               number of observations in the relationship
slope           d log10(measure) / d pH  (negative = higher at lower pH); converted where the paper used ln or another unit, as noted
slope_se        standard error if given or computable; else empty
ci_low, ci_high if given
r               correlation coefficient (Pearson or Spearman, as noted)
p_value
slope_source    reported | computed_from_table | computed_from_reported_points | qualitative_only
DOC_mgL         range or mean
Ca_mgL_or_hardness
alkalinity
temp_C_range
location_in_paper   e.g. "Table 2", "Fig. 4 caption", "SI Table S3", "p. 1945 text"
notes           confounders (season, flow, redox, Fe/Al precipitation), method caveats (resin binding at low pH,
                capacity, diffusion layer), how the slope was obtained
quality         A = computed in this study from tabulated per-sample data | B = regression reported by the authors |
                C = only direction or correlation reported

## File 2: raw tabulated data  ->  extractions/<tag>_raw.csv
Where a paper (or its supplementary information) tabulates per-sample or per-site values of pH together with
DGT-labile and/or filtered concentrations, those rows are copied exactly. Columns:
study_id, site_or_sample, date_or_period, water_type, pH, metal, C_DGT, C_filtered, C_total, unit (ug/L, nM ...),
DOC_mgL, Ca_mgL, temp_C, location_in_paper, notes

## File 3: candidates  ->  extractions/<tag>_candidates.csv
Every relevant paper identified, including those whose full text could not be obtained, as a record of the search.
study_id, authors, year, title, journal, doi, why_relevant, access, data_expected (e.g. "Table 3 has pH and DGT Ni, Zn at 12 lakes")

## File 4: notes  ->  notes/<tag>.md
Search strategy, sources that could and could not be obtained, key qualitative findings and method caveats.
A notes file is kept for the models_lit batch; for the other batches this information is in the notes, access
and why_relevant columns of the extraction tables.

Sources searched included publisher supplementary information (often open where the article is not), PubMed
Central, author and university repositories (including theses) and government reports (USGS, EPA, NIVA, CEH).
Negative and null results were recorded (a slope near zero is data).
