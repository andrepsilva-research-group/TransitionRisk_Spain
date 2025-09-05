# src folder structure

- `libraries.R` - install & load all necessary packages
- `customFunctions2.R` - load created functions necessary throughout the pipeline

**Species Distribution Models**

- `SDM_MIS.R` - complete function to format data, build multiple SDM models, do ensemble and project for current and future conditions. Works for multiple species and is based on biomod2 package's functions
- `SDMfigures.R` - build figures and clean metrics output tables from SDM results for multiple species
      - output 1 - presence maps, using filtered occurrences obtained from GBIF
      - output 2 - Variable importance tables for encsemble models for multiple species
      - output 3 - evalutaion metrics plot for ensemble models for multiple species
      - output 4 - builds polished figure with species distribution maps, per species, from continuous SDM outputs
      - output 5 - builds polished figure of most suitable areas, per species, from binary SDM outputs (using TSS treshold)
      - output 6 - builds final figure of hotspots area maps, by overlapping binary SDM outputs of all species
   
- `mammalTraits_MIS.R` - mines existing databases for trait data of selected mammal species 
