# Spain protected-area prioritization pipeline

Three-script pipeline for national-scale terrestrial conservation prioritization
in Spain under future protected-area (ssp-derived) scenarios, plus a nature-related financial-asset exposure
analysis. 

## Run order

1. `Prioritization_buildlayers_LULC2020_cleaned.R` — builds all input layers
   (cost, protected areas, EUNIS/LULC features, suitability, NCP) and caches
   them under `cache/` and `Input_layers/`. Caches were used to reduce runtime for sequential runs. 
2. `Scenarios_LULC2020_cleaned.R` — solves the SSP1/SSP2/SSP5 prioritizr
   scenarios using the layers from step 1; writes to `solutions/`.
3. `assets_cleaned_v3.R` — overlays financial assets (power, mining,
   industrial zones, agriculture) on the scenario solutions from step 2.


## Large files not included

This repo only keeps small reference/config data. 

You must download these yourself and place each one at the path noted in its
corresponding script 

- **2020 Iberian land cover** (1.2GB) — from **GLC_FCS30D**, Potentially Public ?? -- Confirm with Marta's contact, 
  (https://doi.org/10.5281/zenodo.8239305).
- **2050 SSP5-RCP8.5 land-use projection** (61MB) — a business-as-usual
  land-use projection for 2050 under the SSP5-RCP8.5 scenario, produced by
  the **SEALS model** (Spatial Economic Allocation Landscape Simulator)
  downscaling coarse SSP projections onto the GLC_FCS30D baseline. Public: -- CONFIRM WITH MARTA'S CONTACT
  (https://doi.org/10.5281/zenodo.14506512).
- **Ecosystem type raster** (346MB) — the EEA's **"Ecosystem types of
  Europe"** map (EUNIS level 2, 100m), reprojected for this project. Public:
  (https://www.eea.europa.eu/en/datahub/dataset-redirect/DAT-146-en).
- **Land opportunity cost layer** (65MB) — an ecosystem-services-based cost
  layer. From Spencer et al., 2025 (10.1038/s41597-025-05837-5). 
  confirm the source with whoever originally shared it.
- **WDPA/WDOECM protected areas for Spain** (three files, 70–154MB each) —
  the **World Database on Protected Areas + World Database on OECMs**,
  April 2026 download, filtered to Spain. Public, updated monthly:
  (https://www.protectedplanet.net/en/thematic-areas/wdpa).
- **Industrial zones, mainland/Balearics/Ceuta/Melilla + Canarias,Caso de uso de SIOSE. Suelo Industrial"** 
  (199.95MB and 3.56MB)
  (https://centrodedescargas.cnig.es/CentroDescargas/mapa-suelo-industrial).
- **Spain crop-parcel stack** (5.2GB) — a custom Parquet stack built from
  **EuroCrops**, the EU's harmonized open crop dataset. The underlying
  source is public (for more info: https://github.com/maja601/EuroCrops/wiki/Spain) /
  (Data Source: doi: 10.2905/JRC.FX0BVKR; 10.2905/b9fb9e67-78a9-4327-9d59-39a928d812d3); es_2023_parquet)
  

## What's included here

- The three pipeline scripts
- `GADM/` level-1 and level-2 country-boundary downloads 
- `Input_layers/SuitabilityAverages/scenarioAverages/{Current,ssp585_2050}_mean.tif`
- `input_layers_assets/Canarias.gpkg`
- `data/agri/eurcrops_spain_1km*.tif`
- `solutions/agriculture_province_crop_summary.csv`
- `mask_esp_eck4_2km.tif`, `solution_scenario-A_esp_target-90.tif`
- `Spain_power_mining_28.05.26.xlsx`

