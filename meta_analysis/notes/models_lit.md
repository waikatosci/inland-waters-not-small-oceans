# Published speciation-model outputs vs pH (Part A), tag `models_lit`

Files: `extractions/models_lit_slopes.csv` (6 rows), `extractions/models_lit_candidates.csv` (41 papers).
Built by `models/scripts/07_partA_lit_tables.py`. Our own model runs (Part B) are documented in `models/README.md`.

## How the search was done, and its limits

* The shared web-search budget for this session was used up before Part A started, so no keyword searching was
  possible. Work proceeded from the leads in the brief plus papers already known, retrieved directly.
* Retrieval routes that worked: Crossref metadata look-ups by DOI (for titles, authors, volume and pages), a few
  repository and publisher pages, open PDFs (Xu et al. 2016; a TNO repository PDF), PMC and Frontiers pages, and the
  USGS PHREEQC example files on GitHub.
* Routes that failed: ACS, Elsevier and most publisher hosts (HTTP 403); the NERC (NORA) records for Milne et al.
  2001 and 2003 have no full text; command-line downloads are limited to GitHub and PyPI, so the Visual MINTEQ and
  ORCHESTRA parameter files could not be fetched; the page fetcher refuses URLs that have not appeared in an earlier
  result, so repository pages could not be constructed from titles.
* **Pages that could not be read because of rate limiting (HTTP 429 from api.crossref.org):** look-ups for
  10.1071/EN11016 (Tipping, Lofts and Sonke 2011, Model VII), 10.1023/A:1009627214459 (expected Tipping 1998,
  Model VI) and 10.1021/es102415r (expected Sjöstedt et al. 2010). Those citations are marked "not verified".
* Every candidate whose bibliographic details were recalled rather than checked is marked "not verified" in the
  `why_relevant` or `journal` column. Check them before citing.

## What could be extracted

* **Numbers.** Only one accessible published model output is tabulated per pH: USGS PHREEQC Example 8
  (Parkhurst and Appelo 2013), a Dzombak and Morel-type diffuse-layer model of Zn on hydrous ferric oxide in
  0.1 M NaNO3 (laboratory conditions, not a natural water). Slopes over pH 6.5 to 8.0 computed from the reference
  output: total Zn 1e-7 mol/kg, dissolved -1.64 (SE 0.04), free Zn2+ -1.66; total Zn 1e-4 mol/kg, dissolved -0.35,
  free -0.37. These show how steep a sorption edge is when sorbent is in large excess, and how loading (site
  saturation) flattens it. They are a laboratory reference point, not a freshwater prediction.
* **Qualitative.** Ukotije-Ikwut et al. (2023, Frontiers in Marine Science; Visual MINTEQ with generic
  NICA-Donnan fulvic acid) state that modelled DOM-bound Cu was 96 to 100% and DOM-bound Fe 70 to 100% in blackwater
  rivers and estuaries, variable at pH 3.3 to 6 and "stable at 90-100%" above pH 6 (rows with quality C; no slope).
* **Abstract-level findings.** Lofts and Tipping (2011): WHAM/Model VII agreed best with field free-ion data for Al,
  Co, Ni, Zn and Cd, and voltammetric measurements were systematically higher than predictions. Kalis et al. (2006):
  in situ free-ion concentrations "were difficult to predict using an equilibrium speciation model". Neither abstract
  says anything about pH dependence.

## What the candidates should give once the PDFs are available

* Model free-ion vs pH for fresh waters (most useful for the synthesis figure): Lofts and Tipping 2011 (SI with
  per-sample predictions), Tipping 2005 (Al competition, acid to neutral waters), Tipping et al. 2002 (Al and Fe(III)
  competition), Tipping, Lofts and Lawlor 1998 (Humber rivers), Unsworth et al. 2006 and Sigg et al. 2006 (few sites,
  one pH each), Mueller et al. 2012 and Ahmed et al. 2013 (measured vs predicted free ions).
* Organic-matter model parameters: Milne et al. 2001 and 2003 (generic NICA-Donnan), Rey-Castro et al. 2009,
  Tipping 1998 (Model VI), Tipping, Lofts and Sonke 2011 (Model VII), Gustafsson 2001 (Stockholm Humic Model).
* Oxyanions on ferrihydrite: Gustafsson 2003 (Mo, W), Larsson et al. 2017 (V), Waite et al. 1994 and Mahoney et al.
  2009 (U with carbonate), Appelo et al. 2002, Stachowicz et al. 2008 and Mendez and Hiemstra 2019 (carbonate
  competition). These are laboratory edge studies; figure-only as far as known.
* Biotic Ligand Model pH sensitivity: Di Toro et al. 2001, Santore et al. 2001, De Schamphelaere and Janssen 2002
  (toxicity endpoints; speciation shown in figures).

## Caveats for using model slopes beside field slopes

* Published model curves are mostly for laboratory conditions (sorbent excess, fixed ionic strength, no CO2). Their
  edges are much steeper than what the same constants give in natural waters with about 1 mg/L Fe oxide (compare the
  Example 8 slopes above with `models/outputs/models_own_slopes.csv`).
* In fresh waters pH co-varies with DOC, Ca, alkalinity and season; a "model slope" computed across field samples
  (for example from Lofts and Tipping's SI) mixes these, unlike a single-water pH scan.
* Carbonate competition matters for the oxyanions: if pH is changed by CO2, the surface carbonate load changes with it.
  Model studies that fix pH with acid or base (most laboratory edges) miss this.
* Generic NICA-Donnan parameter tables reproduced in secondary papers disagree (Cu: n1 0.53 and log K2 8.26 in
  Xu et al. 2016 vs 0.63 and 7.26 in Ukotije-Ikwut et al. 2023), so parameters must come from Milne et al. 2003
  itself or from the Visual MINTEQ or ORCHESTRA files.
