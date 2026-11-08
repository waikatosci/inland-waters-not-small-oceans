# Inland waters are not small oceans

Code and data for the analyses behind the manuscript "Inland waters are not small oceans: rising CO₂ and warming act through biology".

Most rivers and lakes hold more CO₂ than the air above them. In this supersaturated water the atmosphere sets only 1/S of the pH response to rising CO₂ (S is the ratio of water to air CO₂), and the remainder follows respiration, which warming accelerates more than photosynthesis. This repository contains the carbonate-system and metabolism model that makes this argument quantitative, the forward model of Lake Karāpiro (Waikato River, New Zealand) from 1850 to 2100, a synthesis of published pH responses of trace metals and oxyanions, the speciation models, the global upscaling across GLORICH river stations, and the oxygen and titration analyses from the Waikato. Every main-text and supplementary display item can be redrawn from the scripts and data here, with the exceptions listed under Data Availability.

The work is the platform paper of the Emerging Climatic Pressures (ECP) research programme on the Waikato River:

> Hartland, A., Allan, M., Arnold, J., Baptista, M. S., Cantalapiedra Mateo, V., Clague, J., Cornwall, C., Duggan, I. C., Durney, P., Esters, L., Fisher, K., Hamilton, D., Hill-Moana, T. P., Holland, A., Kenny, N. J., Lang, J., Lehto, N. J., Marsh, B., Marvin, B., Moorhead, B., Mullarney, J. C., Özkundakci, D., Päffgen, E., Prentice, M. J., Robb, L., Sahlée, E., Stenger, R., Taheri, F., Trifonova, K., Wathen-Smith, G., Wells, N. S. & Wood, S. A. (submitted to Global Change Biology). Inland waters are not small oceans: rising CO₂ and warming act through biology.

---

## Quick Start

```bash
git clone https://github.com/waikatosci/inland-waters-not-small-oceans.git
cd inland-waters-not-small-oceans
pip install -r requirements.txt

cd model
python titrations.py          # parse the Gran titrations -> results/titration_*.csv
python run_projection.py      # projections and Shapley attribution -> results/
python attribution_eq.py      # attribution with equal activation energies
python forward_model.py       # Lake Karāpiro 1850-2100 -> results/forward_karapiro.csv
python ocean_forward.py       # surface ocean under the same forcing -> results/forward_ocean.csv
python fig7_ssp_S.py          # Figure 2
python fig10_forward.py       # Figure 5
```

Each script is run from its own folder and writes figures to a `figures/` subfolder beside it. `run_projection.py` and the along-river figures need the 2022 Waikato survey workbook, which is not yet in the repository (see Data Availability).

---

## Repository Structure

```
inland-waters-not-small-oceans/
│
├── README.md                          This file
├── LICENSE                            MIT License (code)
├── CITATION.cff                       Citation metadata
├── requirements.txt                   pip install path
│
├── model/                             Carbonate system, metabolism, projections and forward model
│   ├── carbonate_system.py            Freshwater carbonate equilibria (Plummer and Busenberg 1982)
│   ├── gas_exchange.py                CO2 transfer velocity and Schmidt-number scaling
│   ├── atmosphere.py                  Atmospheric CO2 1850-2100 (ice core, instrumental, SSPs)
│   ├── unified.py                     Coulombic departure D, metabolic source, Shapley attribution
│   ├── waikato_sites.py               2022 main-stem survey (12 sites)
│   ├── titrations.py                  Gran titrations, Oct 2024 to Aug 2026 (n = 668)
│   ├── run_projection.py              Annual-mean projections and attribution, 2025-2100
│   ├── attribution_eq.py              Attribution with equal activation energies (0.57 eV)
│   ├── forward_model.py               Lake Karāpiro, 1850-2100
│   ├── table_s2.py                    Table S1: one-at-a-time sensitivity
│   ├── fig_concept.py                 Figure 1
│   ├── fig7_ssp_S.py                  Figure 2
│   ├── ocean_forward.py               Surface ocean comparison for Figure 5 (PyCO2SYS)
│   ├── fig10_forward.py               Figure 5
│   ├── fig3_along_river.py            Figure S2
│   ├── fig5_attribution.py            Figure S3
│   ├── fig8_reconstruction.py         Figure S4
│   ├── fig9_drp_observed.py           Figure S5 (needs council monitoring data)
│   ├── fig11_saturation_S.py          Figure S6
│   ├── fig4_titration_S.py            Figure S10
│   ├── warming_sensitivity.py         Sensitivity to water warming across six NZ climate models (Appendix S1)
│   ├── ferrihydrite_concept.png       Figure S1 (illustration)
│   ├── fig_common.py, ngeo_style.py   Figure style
│   ├── data/                          ECP titration workbooks, Karāpiro reference state, NZ CMIP6 warming at Karāpiro
│   └── results/                       Model outputs used in the text and tables
│
├── redox/                             Oxygen at Lake Karāpiro
│   ├── oxygen.py                      Oxygen-temperature coupling and 2100 projections
│   ├── oxygen_daily_layers.csv        Daily profiler means by layer (input to oxygen.py)
│   └── fig_redox_oxygen.py            Figure 4
│
├── meta_analysis/                     Published pH responses, metabolism and speciation
│   ├── SCHEMA.md                      Extraction schema and sign convention
│   ├── extractions/                   Extracted slopes and source records
│   ├── synthesis/
│   │   ├── meta_ph.py                 Evidence classes, medians and interquartile ranges
│   │   ├── meta_metab.py              Table S2: activation energies and GPP/ER
│   │   └── fig_meta.py                Figure 3
│   ├── models/                        PHREEQC diffuse-layer HFO and NICA-Donnan models (Appendix S3)
│   │   ├── README.md                  Model set-up, validation and run order
│   │   ├── scripts/                   00-10, engine and NICA-Donnan implementation
│   │   ├── databases_downloaded/      PHREEQC databases (USGS)
│   │   ├── outputs/                   Slopes and conditions (large curve files omitted)
│   │   └── validation/                USGS Example 8 and independent re-implementation checks
│   └── notes/models_lit.md            Published model outputs used in Figure 3
│
├── glorich/                           pH-phosphate coupling in GLORICH (Appendix S2, Figs S7, S8, S11)
├── upscale/
│   └── glorich_upscale.py             Global upscaling across 2,957 stations (Appendix S3)
│
├── data_sets/
│   ├── build_datasets.py              Builds Data S1 and S2 from meta_analysis/
│   ├── Data_Set_S1_meta_analysis.xlsx
│   └── Data_Set_S2_literature_syntheses.xlsx
│
└── external/glorich/                  Place hydrochemistry.csv from GLORICH here
```

---

## Display Items

| Item | Script | Main inputs |
|---|---|---|
| Figure 1 | `model/fig_concept.py` | none (diagram) |
| Figure 2 | `model/fig7_ssp_S.py` | `model/results/titration_meta.csv` |
| Figure 3 | `meta_analysis/synthesis/fig_meta.py` | `meta_analysis/extractions/`, `meta_analysis/models/outputs/`, `glorich/glorich_mixed_D.json` |
| Figure 4 | `redox/fig_redox_oxygen.py` | `redox/oxygen_daily_2100.csv` |
| Figure 5 | `model/fig10_forward.py` | `model/results/forward_karapiro.csv`, `model/results/forward_ocean.csv` |
| Figure S1 | `model/ferrihydrite_concept.png` | illustration |
| Figure S2 | `model/fig3_along_river.py` | `model/results/projection.csv` |
| Figure S3 | `model/fig5_attribution.py` | `model/results/attribution*.csv` |
| Figure S4 | `model/fig8_reconstruction.py` | model |
| Figure S5 | `model/fig9_drp_observed.py` | council monitoring data (not deposited) |
| Figure S6 | `model/fig11_saturation_S.py` | `model/results/titration_meta.csv` |
| Figs S7, S8 | `glorich/fig_si_surface_mixed.py` | `glorich/glorich_mixed_results.csv` |
| Figure S9 | not deposited | depth-resolved profiler record (held back) |
| Figure S10 | `model/fig4_titration_S.py` | `model/results/titration_*.csv` |
| Figure S11 | `glorich/fig_glorich_controls.py` | `glorich/glorich_*_by_site.csv` |
| Table S1 | `model/table_s2.py` | model |
| Appendix S1, warming sensitivity | `model/warming_sensitivity.py` | `model/data/cmip6_karapiro_T_change.csv` |
| Table S2 | `meta_analysis/synthesis/meta_metab.py` | `meta_analysis/extractions/metab_*.csv` |
| Data S1, S2 | `data_sets/build_datasets.py` | `meta_analysis/` |

The scripts reproduce the values in `results/` and the deposited summary tables exactly. The forward model scales divalent cations at the medians of published field-survey and diel responses (φ = 0.04 and 0.15; Figure 3c), entered in `model/forward_model.py`.

---

## Data Availability

Deposited here: the Gran titrations (`model/data/`), the Lake Karāpiro reference state (`model/data/karapiro_reference.csv`), the downscaled CMIP6 warming over Lake Karāpiro (`model/data/cmip6_karapiro_T_change.csv`), daily bottom-water oxygen at Lake Karāpiro (`redox/`), the published pH-response extractions and Data S1 and S2, and all model outputs.

Not deposited:

- **GLORICH.** The global river chemistry database (Hartmann et al. 2014) is available from its authors. Download `hydrochemistry.csv` and either place it in `external/glorich/` or set `GLORICH_CSV` to its path before running the scripts in `glorich/` and `upscale/`. The station-level results used in the paper are deposited.
- **Council monitoring data and the 2022 survey.** The flow-normalised monitoring records and the 2022 main-stem survey are held by Waikato Regional Council and are available on request. They are needed by `waikato_sites.py`, `run_projection.py` and `fig9_drp_observed.py`; the outputs of those scripts are deposited in `model/results/`. Site alkalinity is the 2015–22 median of council monitoring at each station, deposited as `model/data/wrc_alkalinity_2015_2022.csv`. The CO₂ excess at Lake Karāpiro is the median of 53 council samples at Karāpiro tailrace (2018–2023), each moved to the daily mean with the Hamilton sonde diel curve scaled by 0.28; the result (S = 3.57, pH 7.53, 15.8 °C) is deposited as `model/data/karapiro_reference.csv`, which `waikato_sites.py` uses when the council records are absent. The 2022 survey day is used only for the along-river snapshot (Figure S2). Gran titrations in the 2026 workbook used 0.0901 M acid, read from each sheet by `titrations.py`; earlier titrations used 0.1 M.
- **Climate projections.** Annual mean temperature change, 2080–2099 against 1995–2014, for the six CCAM-downscaled CMIP6 models and four SSPs is from the NIWA core public dataset (Ministry for the Environment, https://climatedata.environment.govt.nz/core-public-dataset.html; CC BY 4.0). The mean over the eight 5 km cells covering Lake Karāpiro is deposited as `model/data/cmip6_karapiro_T_change.csv`.
- **Lake Karāpiro depth profiles.** The depth-resolved profiler record and deployment chemistry at Lake Karāpiro are held back while they are analysed further within the ECP programme, and are available from the corresponding author on request. They underlie Figure S9; the daily bottom-water oxygen means used for Figure 4 are deposited in `redox/`.

---

## Citation

Manuscript:
> Hartland, A. et al. (in preparation). Inland waters are not small oceans: rising CO₂ and warming act through biology. For submission to Global Change Biology.

Global river chemistry:
> Hartmann, J., Lauerwald, R. & Moosdorf, N. (2014). A brief overview of the GLObal RIver CHemistry Database, GLORICH. Procedia Earth and Planetary Science 10, 23-27.

---

## License

Code: MIT License. See [LICENSE](LICENSE). Data: CC BY 4.0.

---

## Acknowledgments

Funded by the Ministry of Business, Innovation and Employment Endeavour Fund through the Emerging Climatic Pressures programme (LVLX2302). We thank Danielle Blackwell for ICP-MS analyses and Waikato Regional Council for monitoring data.

For questions or contributions, open an issue on GitHub or contact the corresponding author: [adam.hartland@lincolnagritech.co.nz](mailto:adam.hartland@lincolnagritech.co.nz)

---

## AI Assistance Statement

The research questions, conceptual framework and interpretations in this repository are the work of the authors. Generative AI (Anthropic's Claude) was used under the authors' direction for literature searches and the extraction of published values for the meta-analysis and literature syntheses (Data S1 and S2), which were checked against the source publications; for writing and refactoring the code for the carbonate, speciation, upscaling and forward models and for the figures; and for checking that values, cross-references and figure numbering agree across the code, data and manuscript. All AI-assisted code and data were reviewed and verified by the authors, including the speciation model against Visual MINTEQ, and the authors take full responsibility for the content and correctness of the code, data and documentation.
