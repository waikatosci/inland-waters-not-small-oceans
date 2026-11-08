# Speciation-model estimates of how trace elements respond to pH

Model slopes, d log10(quantity) / d pH, to sit beside the observed slopes and the electrostatic limit in the
synthesis figure. This folder holds Part B (our own calculations). Part A (published model outputs) is in
`../extractions/models_lit_slopes.csv`, `../extractions/models_lit_candidates.csv` and `../notes/models_lit.md`.

Sign convention as in `../SCHEMA.md`: negative slope = more of the element dissolved (or free) at lower pH.

## Headline result: soft (Waikato-like) water, pH 6.5 to 8.0, no organic matter

Least-squares slope of log10(quantity) on pH. "Dissolved" = not sorbed to hydrous ferric oxide (HFO) = free ion
plus inorganic complexes (the inorganic labile pool; with no organic matter these are the same thing).
Main model = Dzombak and Morel diffuse-layer HFO model with carbonate adsorption.

| Element | No oxide: free ion | HFO 0.1 mg Fe/L: dissolved | HFO 1 mg Fe/L: dissolved | HFO 1 mg Fe/L: free ion | HFO 1 mg Fe/L, no carbonate adsorption: dissolved |
|---|---|---|---|---|---|
| Mn(II) | 0.00 | 0.00 | -0.04 | -0.04 | -0.04 |
| Ni | -0.02 | -0.02 | -0.18 | -0.20 | -0.19 |
| Co(II) | -0.01 | 0.00 | -0.05 | -0.06 | -0.05 |
| Zn | -0.04 | -0.06 | -0.41 | -0.45 | -0.43 |
| Cd | -0.01 | -0.01 | -0.17 | -0.19 | -0.18 |
| Cu(II) | -0.72 | -0.43 | -0.51 | -1.23 | -0.54 |
| Pb | -0.54 | -0.33 | -0.32 | -0.86 | -0.40 |
| P (phosphate) | -0.02 | +0.02 | +0.15 | +0.14 | +0.63 |
| Mo (molybdate) | 0.00 | 0.00 | 0.00 | 0.00 | +0.01 |
| As(V) | 0.00 | -0.14 | -0.38 | -0.38 | +0.43 |
| V(V) (low confidence) | 0.00 | -0.07 | -0.55 | -0.55 | -0.14 |
| Se(IV) | 0.00 | 0.00 | -0.01 | -0.01 | +0.09 |

Same quantity (dissolved, HFO 1 mg Fe/L, main model, pH 6.5 to 8.0) across the three waters:

| Element | soft | moderate | hard |
|---|---|---|---|
| Mn(II) | -0.04 | -0.02 | -0.01 |
| Ni | -0.18 | -0.09 | -0.05 |
| Co(II) | -0.05 | -0.02 | -0.01 |
| Zn | -0.41 | -0.23 | -0.15 |
| Cd | -0.17 | -0.08 | -0.05 |
| Cu(II) | -0.51 | -0.40 | -0.31 |
| Pb | -0.32 | -0.15 | -0.11 |
| P (phosphate) | +0.15 | +0.05 | +0.01 |
| Mo (molybdate) | 0.00 | 0.00 | 0.00 |
| As(V) | -0.38 | -0.64 | -0.68 |
| V(V) | -0.55 | -0.58 | -0.49 |
| Se(IV) | -0.01 | -0.02 | -0.02 |

What the numbers say, in plain terms:

* With a realistic amount of iron oxide (1 mg Fe/L) the model gives only shallow cation slopes (about -0.04 for Mn
  to -0.5 for Cu). The oxide is too scarce to hold most of the metal, so the dissolved pool changes little with pH.
  The 0.1 mg Fe/L case is flatter still. Harder water (more Ca and Mg competing for the surface) flattens it further.
* Cu and Pb free-ion slopes are steep even without any oxide (-0.72 and -0.54), because carbonate and hydroxide
  complexes grow as pH rises at fixed alkalinity, and steeper with it (-1.23 and -0.86). Most of this is inorganic
  complexation, not sorption.
* For the oxyanions the sign depends on whether carbonate competes for the oxide surface. Lowering pH by adding CO2
  loads the surface with carbonate (56% of HFO sites at pH 6.5 in the soft water with 1 mg Fe/L, 12% at pH 8). With carbonate
  adsorption on, arsenate and vanadate become less sorbed at low pH (negative slope); with it off, arsenate shows the
  usual positive oxyanion slope (+0.43) and phosphate a steeper one (+0.63 instead of +0.15). This is the single largest
  model-structure uncertainty here, and it bears directly on a CO2-driven mechanism.
* In the main model molybdate is under 1% sorbed and selenite 5 to 10% (pH 6.5 to 8.0, soft water, 1 mg Fe/L) (phosphate, sulfate
  and carbonate take the sites), so their slopes are near zero (+0.01 and +0.09 without carbonate adsorption).
* None of this includes organic matter, which dominates Cu and Pb (and matters for Ni and Zn) in real waters. See
  section 8.

## Reproduce

```
cd meta_analysis/models
python3 scripts/00_check_install.py          # software and database versions -> outputs/environment.txt
python3 scripts/01_inspect_databases.py      # which HFO constants exist where -> outputs/hfo_constants_*.csv
python3 scripts/02_validate_usgs_ex8.py      # USGS Example 8 reproduction -> validation/usgs_ex8_comparison.csv
python3 scripts/03_run_hfo_model.py          # main runs (a few seconds) -> outputs/hfo_raw_*.csv
python3 scripts/03_run_hfo_model.py --engine 373   # same runs on the second engine -> validation/hfo_raw_*_engine373.csv
python3 scripts/04_curves_and_slopes.py      # -> outputs/models_own_curves.csv, models_own_slopes.csv, models_own_conditions.csv
python3 scripts/05_validate_dlm_independent.py   # independent re-implementation check -> validation/*.csv
python3 scripts/06_plot_curves.py            # -> outputs/models_own_curves.png
python3 scripts/07_partA_lit_tables.py       # Part A tables -> ../extractions/models_lit_*.csv
python3 scripts/08_validate_nicad.py         # NICA-Donnan validation -> validation/nicad_validation_V*.csv
python3 scripts/09_run_nicad_coupled.py      # FA + inorganic + HFO coupled runs (about 2 min) -> outputs/nicad_raw.csv
python3 scripts/10_nicad_curves_slopes_plot.py   # -> outputs/models_nicad_curves.csv, models_nicad_slopes.csv, .png
```

## Files

| File | Content |
|---|---|
| `scripts/model_setup.py` | Every assumption in one place: waters, trace totals, HFO properties, redox fix, added constants |
| `scripts/iphreeqc.py`, `scripts/engine.py` | Thin wrappers around two independent builds of the USGS IPhreeqc engine |
| `outputs/models_own_curves.csv` | Long format: model, water, sorbent, Fe_mgL, DOC_mgL, pH, pCO2_uatm, alkalinity, DIC, SI_calcite, element, quantity, value, unit (pH 6.0 to 9.0 in 0.05 steps) |
| `outputs/models_own_slopes.csv` | model, water, sorbent, DOC_mgL, element, quantity, pH_range, slope, notes (notes hold n, r2, local slope range, end values, pCO2 range, flags) |
| `outputs/models_own_conditions.csv` | Per run and pH: pCO2, alkalinity, ionic strength, surface potential, sorbed Ca/Mg/C/S, saturation indices, redox checks |
| `outputs/models_own_curves.png` | Figure: log10 dissolved fraction vs pH, colour = water, solid = HFO 1 mg Fe/L, dashed = 0.1, dotted = 1 mg Fe/L without carbonate adsorption (oxyanions) |
| `outputs/hfo_constants_inventory.csv`, `hfo_constants_presence.csv` | Every HFO surface reaction and log K in each database inspected |
| `outputs/environment.txt` | Versions and md5 sums of the database files |
| `validation/usgs_ex8_comparison.csv` | Our output vs the USGS reference output for Example 8 |
| `validation/dlm_independent_check.csv` | PHREEQC vs our independent implementation, per pH and element |
| `validation/dm_standard_conditions_edges.csv` | Diagnostic single-element sorption edges under Dzombak and Morel's laboratory conditions |
| `scripts/nicadonnan.py`, `scripts/08-10_*.py` | NICA-Donnan model, its validation, coupled runs, slopes and figure (section 10) |
| `outputs/models_nicad_slopes.csv`, `models_nicad_curves.csv`, `models_nicad_curves.png`, `nicad_raw.csv` | Organic-matter runs (section 10) |
| `validation/nicad_validation_V1-V5*.csv` | NICA-Donnan validation against Milne 2001/2003 and Rey-Castro 2009 |
| `databases_downloaded/` | The PHREEQC databases used, as fetched from the USGS phreeqc3 GitHub repository on 2026-10-07 |

## 1. Software (what worked)

* `pip install phreeqpython --break-system-packages` worked: phreeqpython 1.6.2, which bundles IPhreeqc 3.7.3-15968
  and four databases (phreeqc.dat, pitzer.dat, vitens.dat, wateq4f_PWN.dat). It does not bundle minteq.v4.dat or
  the standard wateq4f.dat, its Python layer does not return the selected-output table, and its IPhreeqc cannot read
  the current USGS phreeqc.dat (new keyword GAS_BINARY_PARAMETERS). We therefore call its library directly
  (`iphreeqc.py`) and use it only as a cross-check engine.
* `pip install phreeqc --break-system-packages` (PyPI package phreeqc 1.1.1, IPhreeqc 3.8.6-17100) is the primary
  engine for all results.
* Databases phreeqc.dat, wateq4f.dat, minteq.v4.dat, minteq.dat and llnl.dat were downloaded from
  https://github.com/usgs-coupled/phreeqc3/tree/master/database (md5 sums in `outputs/environment.txt`).
* Python 3.11, numpy 2.4, scipy 1.17, pandas 3.0, matplotlib 3.10.

## 2. Surface constants: what the databases contain

All are Dzombak and Morel (1990) two-site diffuse-layer constants (log K; s = strong site, w = weak site) unless
noted. Full listing with reactions: `outputs/hfo_constants_inventory.csv`.

| Sorbate | phreeqc.dat | wateq4f.dat | minteq.v4.dat |
|---|---|---|---|
| Mn(II) | s -0.4, w -3.5 (D&M Table 10.5, estimated) | same | absent |
| Ni | absent | s 0.37, w -2.5 | s 0.37, w -2.5 |
| Co(II) | absent | absent | s -0.46, w -3.01 |
| Zn | s 0.99, w -1.99 | same | same |
| Cd | s 0.47, w -2.91 | same | s 0.47, w -2.90 |
| Cu(II) | s 2.89, w 0.6 | same | same |
| Pb | s 4.65, w 0.3 | same | same |
| Ca | s 4.97, w -5.85 | same | same |
| Mg | w -4.6 | same | same |
| Phosphate | w 31.29 / 25.39 / 17.72 | same | same, on both site types |
| Molybdate | absent | absent | 9.5 / 2.4, both site types |
| Arsenate | absent | w 29.31 / 23.51 / 10.58 (from AsO4-3) | 8.61 / 2.81 / -10.12 (from H3AsO4; same constants), both site types |
| Vanadate | absent | absent | -16.63, one species Hfo_OHVO4-3 (from VO2+), both site types |
| Selenite | absent | w 12.69 / 5.17 (from SeO3-2) | 4.29 / -3.23 (from HSeO3-; same constants), both site types |
| Sulfate | w 7.78 / 0.79 | same | same, both site types |
| Carbonate | w 12.56 / 20.62 (Van Geen et al. 1994 re-optimised for HFO; Appelo et al. 2002) | same | absent |
| Acid-base | pKa1 7.29, pKa2 8.93 | same | same |

Also present but not used: Sr, Ba, Fe(II), Ag, Be, Hg, Sn, Cr(III), arsenite, borate, selenate, chromate,
antimonate, fluoride, silicate (phreeqc.dat), uranyl (wateq4f.dat). No database has tungstate, so W was not modelled.
The bundled phreeqpython `wateq4f_PWN.dat` lists arsenate constants 31.29 / 25.39 / 10.58, i.e. the phosphate values for
two of the three species; this differs from the USGS wateq4f.dat and was not used.

**What the runs use.** minteq.v4.dat for everything (aqueous species and surface constants), because it is the only
database with Co, Mo and V surface constants, plus two blocks copied verbatim from the USGS wateq4f.dat: the Mn(II)
constants and the carbonate constants (main model only). No constant was taken from anywhere else, and none was
invented or estimated by us. In minteq.v4.dat anions bind to strong as well as weak sites (MINTEQA2 convention); the
strong sites are 2.5% of the total, so this changes little.

**Vanadate is low confidence.** minteq.v4.dat has a single, fully deprotonated surface species (equivalent to
log K 13.57 for Hfo_wOH + VO4-3 = Hfo_wOHVO4-3 with the database's own vanadate hydrolysis constants) and does not
state its source. Under Dzombak and Morel's standard laboratory conditions this makes V(V) sorption peak at pH 6.4
and fall towards lower pH (`validation/dm_standard_conditions_edges.csv`), which is unusual for an oxyanion. A
documented constant set (for example Larsson et al. 2017, CD-MUSIC; listed in the Part A candidates) would be needed.

## 3. Waters, temperature and how pH was set

16 °C. Concentrations in mg/L; alkalinity in meq/L as HCO3-.

| | Alkalinity | Ca | Mg | Na | K | Cl | SO4 | Ionic strength (mmol/L) |
|---|---|---|---|---|---|---|---|---|
| Soft (Waikato-like) | 0.78 | 5.3 | 2 | 15.8 (10 specified) | 2 | 10 | 5 | 1.42 to 1.44 |
| Moderately hard | 2 | 30 | 6 | 10 | 2 | 13.2 (10 assumed) | 5 | 3.39 to 3.47 |
| Hard | 4 | 60 | 12 | 10 | 2 | 12.9 (10 assumed) | 5 | 6.0 to 6.3 |

* Charge balance. The soft water as specified has an anion excess of 0.25 meq/L, closed by raising Na from 10 to
  15.8 mg/L (0.686 mmol/L). The brief gave only alkalinity, Ca and Mg for the harder waters; we used the soft-water
  Na, K, Cl and SO4 for them, which leaves a small cation excess, closed by raising Cl to 13.2 mg/L (moderate) and
  12.9 mg/L (hard).
* pH set by CO2 at fixed alkalinity. Each water is defined once at pH 7.0; for every target pH (6.0 to 9.0 in 0.05
  steps) PHREEQC adds or removes CO2 gas until the target pH is reached (pseudo-phase Fix_H+ with CO2(g) as the
  reactant). CO2 exchange does not change alkalinity; the only other proton exchange is with the oxide surface, which
  shifts alkalinity by at most 0.36% (soft water, 1 mg Fe/L). Each pH point starts from the same initial water and
  oxide, so points are independent equilibria.
* pCO2 needed (µatm, no oxide): soft 43,000 at pH 6.0, 13,600 at 6.5, 4,300 at 7.0, 425 at 8.0, 131 at 8.5, 39 at 9.0;
  moderate 33,900 at 6.5 and 1,060 at 8.0; hard 66,000 at 6.5 and 2,050 at 8.0. The requested 100 to 20,000 µatm
  spans pH 6.33 to 8.61 in the soft water, 6.73 to 8.97 in the moderate water and 7.02 to 9.0 in the hard water. We
  kept one pH grid for all waters so slopes compare like with like; slope notes flag windows that need more than
  20,000 µatm.
* Activity corrections: PHREEQC's treatment as defined species by species in minteq.v4.dat (Debye-Hückel-type or
  Davies). Equilibrium constants are corrected to 16 °C where the database gives an enthalpy;
  the HFO surface constants have none, so they are used at their 25 °C values.

## 4. Hydrous ferric oxide

Dzombak and Morel (1990) standard properties, the same as in USGS PHREEQC Example 8: 89 g HFO per mol Fe, specific
surface area 600 m2/g, 0.2 mol weak sites and 0.005 mol strong sites per mol Fe. Diffuse-layer (Gouy-Chapman)
electrostatics, PHREEQC default (no explicit diffuse-layer composition).

| HFO as Fe | mol Fe/L | g HFO/L | weak sites (mol/L) | strong sites (mol/L) | surface area (m2/L) |
|---|---|---|---|---|---|
| 0.1 mg/L | 1.79e-6 | 1.59e-4 | 3.58e-7 | 8.95e-9 | 0.096 |
| 1 mg/L | 1.79e-5 | 1.59e-3 | 3.58e-6 | 8.95e-8 | 0.96 |

The oxide amount is held fixed at every pH (no Fe dissolution or precipitation).

## 5. Trace elements, redox and precipitation

* Totals (system totals, all present together so they compete for sites): Mn(II) 20, Ni 1, Co(II) 0.1, Zn 2,
  Cd 0.02, Cu(II) 1, Pb 0.1, P 20, Mo 0.5, As(V) 2, V(V) 1, Se(IV) 0.2 µg/L. Note that Mn (0.36 µmol/L) exceeds the
  strong sites at both oxide levels and P (0.65 µmol/L) exceeds the weak sites at 0.1 mg Fe/L, so competition is real.
  Fe and Al were not modelled as trace elements (Fe is the sorbent; Al has no HFO constants here).
* Redox held fixed. In a PHREEQC batch reaction every redox element re-equilibrates with pe. To keep Mn(II), Co(II),
  Cu(II), As(V), V(V), Se(IV) and S(VI) in the stated states, the reactions linking them to their other valence
  states were redefined so the other states are at least 60 orders of magnitude less stable (listed with original
  values in `scripts/model_setup.py`). Check: the other valence states are below 1e-20 mol/kg in every run.
* No mineral was allowed to precipitate. Saturation indices are reported in `models_own_conditions.csv`. Within
  pH 6.5 to 8.0: calcite reaches SI +0.59 (hard) and +0.04 (moderate), and is undersaturated in the soft water;
  hydroxylapatite reaches SI +0.11 (soft), +3.4 (moderate) and +4.4 (hard). Above pH 8 calcite reaches +1.5 (hard,
  pH 9). Everything else (rhodochrosite, smithsonite, otavite, tenorite, malachite, cerussite, hydrocerussite,
  Ni(OH)2, CoCO3, CaMoO4, Zn(OH)2) stays undersaturated except rhodochrosite at +0.04 in the hard water at pH 9.
  Suppressing precipitation is deliberate: calcite and apatite supersaturation is common in surface waters
  (precipitation is slow), and letting calcite precipitate would remove alkalinity and break the fixed-alkalinity
  design.

## 6. Quantities and slopes

* `dissolved_fraction` = (total - sorbed) / total. Because the total is fixed, its slope equals the slope of the
  dissolved concentration.
* `free_ion`: cations, the free aquo ion (Mn+2, Ni+2, ...). Oxyanions, the sum of the protonation forms of the
  monomeric oxyanion (for example H3PO4 + H2PO4- + HPO4-2 + PO4-3), excluding complexes with Ca, Mg, Na, K or metals.
  Concentrations in mol/kg water (equal to mol/L within 0.1% here).
* `inorganic_labile` = free ion + all dissolved inorganic complexes. With no organic matter this equals `dissolved`.
* Slopes: ordinary least squares of log10(quantity) on pH over every grid point in the window, 31 points for
  pH 6.5 to 8.0 and 7.0 to 8.5. Curves are not straight lines, so the notes give r2 and the range of local slopes;
  for Pb (r2 about 0.64) the window slope hides a minimum near pH 7.5.
* Two model variants: `DLM_HFO_minteqv4+CO3` (main, with carbonate adsorption) and `DLM_HFO_minteqv4_noCO3`
  (Dzombak and Morel's original reaction set, no carbonate adsorption). Sorbent levels: none, 0.1 and 1 mg Fe/L.
  `DOC_mgL` is 0 in every row.

## 7. Validation

1. **USGS PHREEQC Example 8 (Zn sorbed on HFO; values reported by another source).** Both engines reproduce the
   reference output shipped with PHREEQC exactly (maximum relative difference 0 over 78 values each; two Zn
   loadings, pH 5 to 8). pH of 50% sorption: 6.42 (Zn 1e-7 mol/kg) and 7.48 (1e-4 mol/kg). This tests our
   engines, wrappers and parsing; Example 8 uses its own constants (pKa 7.18/8.82, Zn 0.66/-2.32).
2. **Independent re-implementation (Cd and phosphate edges, plus the production runs).** `05_validate_dlm_independent.py`
   solves the Dzombak and Morel equations from scratch in Python (surface mass action with the Boltzmann factor,
   site balances, Gouy-Chapman charge-potential relation, water permittivity from Malmberg and Maryott 1956), using
   the aqueous activities and ionic strength that PHREEQC reports and constants parsed from the same database text.
   Agreement with PHREEQC: Cd edge (1e-7 M Cd, 1e-3 M Fe, 0.1 M NaNO3, 25 °C) within 0.03% for sorbed Cd;
   phosphate edge (1e-5 M) within 0.4%; surface potential within 0.03 mV; the production soft water with
   1 mg Fe/L within 0.03% for all 12 elements and the hard water with 0.1 mg Fe/L within 0.14%. pH of 50% sorption:
   Cd 6.62; phosphate 8.69 (falling edge).
3. **Comparison with Dzombak and Morel's own published edges: not done.** The book and the papers with its edges
   could not be accessed from this session (publisher hosts blocked, web search budget exhausted), and no other
   accessible source reported a phosphate or Cd edge with full conditions. To make that check, compare
   `validation/dm_standard_conditions_edges.csv` (pH of 50% sorption under the standard laboratory conditions:
   Pb 4.35, Cu 5.19, Zn 6.31, Cd 6.62, Ni 6.66, Mn 7.09, Co 7.11; Mo 6.99, Se(IV) 8.54, P 8.75 and As(V) 10.13 on the
   falling edge) with the book's figures. The cation order (Pb > Cu > Zn > Cd and Ni > Mn and Co) follows the usual
   affinity sequence for iron oxides.

4. **Second engine.** The full production runs repeated with IPhreeqc 3.7.3 (the phreeqpython build) are identical
   to the IPhreeqc 3.8.6 results at the printed precision (5 significant figures) for every element, water, oxide
   level and pH (`validation/hfo_raw_*_engine373.csv`).

So the numerical implementation is verified three ways; whether the constants themselves reproduce Dzombak and Morel's
published curves is inherited from the USGS databases and not independently checked here.

## 8. Organic matter (NICA-Donnan): not run in round 1 (superseded by section 10)

The brief allowed NICA-Donnan only if the generic fulvic acid parameters (Milne et al. 2001 protons, Milne et al.
2003 metals) could be obtained and verified, and validated against a worked example. That bar was not met, so no
organic-matter calculation was run and `DOC_mgL` is 0 throughout. What was found:

* The primary papers (ES&T) were not accessible (publisher 403; the NERC repository has no full text), the Visual
  MINTEQ and ORCHESTRA parameter files could not be downloaded, and the web search budget for this session ran out.
* Proton parameters for generic fulvic acid agree in two secondary sources (Table A3 of Annals of Botany
  105:1183, doi 10.1093/aob/mcq098; Table 1 of Xu et al. 2016, J Colloid Interface Sci 473:141): Qmax1 5.88 and
  Qmax2 1.86 mol/kg, log KH1 2.34 and log KH2 8.60, nH1 0.66 and nH2 0.76, p1 0.59 and p2 0.70. The Donnan parameter
  b = 0.57 appears only in Xu et al. Ca: log K1 -2.13, n1 0.85, log K2 -3.00, n2 0.80 (Annals of Botany only).
* Cu: Xu et al. give log K1 0.26, n1 0.53, log K2 8.26, n2 0.36 (internally consistent with their nCu/nH ratios);
  a Visual MINTEQ-based table (Ukotije-Ikwut et al. 2023, Front. Mar. Sci., Table 3) prints n1 0.63 and log K2 7.26.
  The two disagree. Cd: log K1 -0.97 and log K2 0.5 in one supplementary table (n not given).
* Zn, Ni, Co, Mn, Pb, Mg and the Al and Fe(III) competitors: no verifiable source found.
* No published worked example with complete inputs was accessible for validation.

To run it later: supply Milne et al. (2003) Tables (or the Visual MINTEQ / ORCHESTRA parameter file) and one worked
example. The brief's DOC conversion (dissolved organic matter = 2 x DOC, 65% active fulvic acid) is the usual WHAM
convention (Bryan et al. 2002); Visual MINTEQ users often calibrate an active-DOM/DOC ratio instead (Sjöstedt et al.
2010). Neither source could be opened here to confirm the exact figures.

## 9. Main limitations

* No organic matter in sections 1-7 (now added in section 10). For Cu and Pb (and partly Ni and Zn) organic complexation controls the free and labile pools in
  real fresh waters, and proton competition on humic sites adds its own pH dependence. The inorganic-only free-ion
  slopes here are therefore not the whole story for those metals.
* One sorbent (HFO). No Mn oxides, Al hydroxides, clays or particulate organic matter. HFO amount independent of pH.
* Equilibrium only; no sorption kinetics, no diffusion-layer effects of the kind that shape DGT measurements.
* Surface constants at 25 °C applied at 16 °C.
* Carbonate adsorption constants (Appelo et al. 2002, as carried in the USGS databases) decide the sign of the As(V)
  slope and the size of the phosphate slope. Report both variants.
* Vanadate constant of unknown provenance (section 2).
* Window slopes summarise curved relationships; use `models_own_curves.csv` for local slopes.

## 10. Organic matter added: generic NICA-Donnan fulvic acid (round 2, tag models2)

Full notes, parameter tables with page references, and WHAM Model VI/VII parameters: `../notes/models2.md`.

* **Model.** `scripts/nicadonnan.py`: consistent bimodal NICA isotherm with Donnan electrostatics (Milne et al.
  2003 eq 1; Donnan volume log V_D = b(1 - log I) - 1), generic fulvic acid parameters from Milne et al. 2003
  Table 3 (p. 968; proton values = Milne et al. 2001 Table 4). Concentrations, not activities, as in Milne's fits.
* **Validation (`08_validate_nicad.py`), all passed.** Milne 2001 eq 3 (generic FA proton charge at 0.1 M):
  within 0.07 eq/kg, charge change pH 4-9 within 0.02%. Rey-Castro et al. 2009 Table 1 (Visual MINTEQ, Arve River
  water + 0.5 mg/L FA): Boltzmann factor 1.167 vs 1.17 per unit charge, Donnan concentrations within 0.005 log;
  free Cu, Zn, Ni, Pb, Cd recomputed from the totals within 0.015 log units. Rey-Castro Table 2 affinities within
  0.05 (Al 0.07). Milne 2003 Fig. 4: Zn binding to generic HA rises 1.82-1.86 log units from pH 4 to 8 (reported
  1.85). With DOC = 0 the coupled code reproduces the round-1 results within 5e-5.
* **Runs (`09_run_nicad_coupled.py`).** Same three waters, pH by CO2 at fixed total alkalinity (FA enters as its
  salt at pH 7; its proton release or uptake is exchanged with the carbonate system), DOC 0, 1, 3, 10 mg/L,
  FA = 1.65 x DOC (Sjostedt et al. 2010; variant 1.27 x DOC, Lofts and Tipping 2011), HFO 0 or 1 mg Fe/L (main
  carbonate-adsorbing DLM). Iterative coupling with PHREEQC; mass balance closed to 1e-8. Variant `+AlOH3` adds Al3+
  competition with Al fixed by Al(OH)3 solubility (changes slopes by <= 0.02). No Fe(III) competition.
* **Outputs.** `outputs/models_nicad_slopes.csv` (same columns as models_own_slopes.csv; quantities free_ion,
  inorganic_labile, dissolved_fraction where dissolved includes FA-bound metal), `outputs/models_nicad_curves.csv`,
  `outputs/models_nicad_curves.png`, `outputs/nicad_raw.csv`.

Soft water, DOC 3 mg/L, pH 6.5-8.0 (main model), d log10 / d pH:

| Element | No oxide: free ion | No oxide: inorganic labile | HFO 1 mg Fe/L: free ion | HFO 1 mg Fe/L: inorganic labile |
|---|---|---|---|---|
| Mn(II) | -0.02 | -0.02 | -0.05 | -0.05 |
| Co(II) | -0.03 | -0.01 | -0.07 | -0.06 |
| Ni | -0.12 | -0.10 | -0.26 | -0.23 |
| Zn | -0.04 | 0.00 | -0.46 | -0.42 |
| Cd | -0.25 | -0.24 | -0.27 | -0.26 |
| Cu(II) | -1.66 | -0.95 | -1.66 | -0.94 |
| Pb | -1.03 | -0.49 | -0.94 | -0.40 |

In plain terms: adding fulvic acid makes the free-ion response to pH much steeper for the metals it binds strongly
(Cu, Pb, Cd, and Ni at higher DOC), because organic binding strengthens as pH rises and the fulvic acid loses
protons. Cu free ion falls by a factor of about 45 per pH unit. Zn, Mn and Co hardly change, because the generic
fulvic acid binds them weakly (these parameters rest on one or two datasets). Free-ion slopes steepen with DOC
(free Cu -1.51, -1.66, -1.78 at DOC 1, 3, 10; free Ni -0.05, -0.12, -0.27, no oxide). The filtered ("dissolved")
concentration of Cu and Pb is almost flat, because fulvic acid keeps them in solution.

Lofts and Tipping (2011) give no per-sample pH in the main paper, so no WHAM free-ion vs pH slope could be extracted
(details in `../notes/models2.md`; the Accessory publication with the inputs is listed in
`../extractions/models2_candidates.csv`).
