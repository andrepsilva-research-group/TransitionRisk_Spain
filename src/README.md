# src folder structure

## Environment Setup

1.  `libraries.R` - Installs and loads all necessary R packages used throughout the pipeline.

    *Run this script first to make sure all dependencies are in place.*

2.  `customFunctions2.R` - Loads all custom helper functions used across different parts of the pipeline.

## Species Distribution Models

3.  `SDM_MIS.R`- This is the main SDM pipeline script. It:

    -   Formats and prepares species occurrence and environmental data.
    -   Builds multiple SDM models per species (using biomod2 functions).
    -   Creates ensemble models.
    -   Projects SDMs for both current and future environmental conditions.

Requirements: A folder containing training and prediction landscapes (raster files .tif). Download [here](https://drive.google.com/drive/folders/1ik_lyTc5CdcX6HUMIfU7cKGgvi0UqzwH?usp=sharing).

4.  `SDMfigures.R` - This script processes SDM results and produces a range of **summary tables and figures** per species and one final figure combining all.

**Outputs include:**
- **Presence maps** – using filtered GBIF occurrences per species.
- **Variable importance table** – for ensemble models and multiple species
- **Evaluation plots** – performance metrics for ensemble models 
- **Continuous SDM maps** – polished per-species maps of predicted suitability.
- **Binary SDM maps** – binary map of suitable/unsuitable areas, per species, using TSS threshold.
- **Hotspot maps** – overlap of all species’ binary outputs to obtain hotspot areas where threathned species are more likely to occur areas.


5. `mammalTraits_MIS.R` - Fetches and compiles trait data for the selected mammal species from existing trait databases.

Only needed if the metaRange module is being used.
If not required, this script can safely be removed later.

