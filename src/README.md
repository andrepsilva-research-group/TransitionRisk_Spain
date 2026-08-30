# src folder structure
> [!WARNING]
> WE NEED TO ADD FINAL ZENODO LINK HERE
### Setup the environment

1.  `libraries.R`  
Installs and loads all necessary R packages used throughout the pipeline.  
   *Run this script first to make sure all dependencies are in place.*

2.  `customFunctions2.R`  
Loads all custom helper functions used across different parts of the pipeline.

### Get taxa occurrences

3. `taxaOccurrences_mammals.R` & `taxaOccurrences_birds.R`  
Retrieves and processes species occurrence data for a target region (e.g. the Iberian Peninsula) using information from the IUCN Red Lists and the GBIF database.  
Produces a .csv file of cleaned GBIF occurrence records from a list of target species.<br>     

There are two alternative ways to define the list of target species:

**Option 1**: Using IUCN Spatial Data (polygons showing where each species occurs)
- Download spatial data:
   -  Go to: https://www.iucnredlist.org/resources/spatial-data-download
   -  Select your taxonomic group (e.g. Terrestrial Mammals, Birds, etc.)
   -  Log in (free IUCN account required)
   -  Download the ZIP file containing .shp shapefiles
-  Load and filter in R:
   -  Keep only threatened species (CR, EN, VU)
   -  Intersect the ranges with your area of interest (e.g. Spain and Portugal)
   -  Optionally remove or add species manually
-  Feed the list into the custom function to retrieve and clean GBIF occurrence records fro those species.

**Option 2**: Using an External Species List (manually downloaded or curated list).
- Obtain a CSV list:
   - For instance for bird species a list can be obtained from: https://www.birdlife.org/
   - Download a regional or national species list (e.g. “Iberian Peninsula”)
   - Save it under data/ (e.g. data/BirdLife_Iberian_peninsula.csv)
- (Optional) Add taxonomy automatically, using the taxize package:
-  Feed the list into the custom function to retrieve and clean GBIF occurrence records fro those species.
  
### Run species distribution models

4.  `SDM.R`- This is the main pipeline script. It:

    -   Formats and prepares species occurrence and environmental data.
    -   Builds multiple SDM models (Random Forest, Extreme Gradient Boosting, Artificial Neural Network, Maximum Entropy Model - MAXNET) per species (using [biomod2](https://biomodhub.github.io/biomod2/) functions).
    -   Create an ensemble model.
    -   Projects SDMs for both current and future (2030, 2050) environmental conditions.

Requirements: A folder containing training and prediction landscapes (raster files .tif). Download [here](https://drive.google.com/drive/folders/1ik_lyTc5CdcX6HUMIfU7cKGgvi0UqzwH?usp=sharing).

### Analyse outputs and Visualise results

5.  `SDMfigures.R`
Processes SDM results and produces a range of **summary tables and figures** per species and one final figure combining all.

**Outputs include:**
- **Presence maps** – using filtered GBIF occurrences per species.
- **Variable importance table** – for ensemble models and multiple species
- **Evaluation plots** – performance metrics for ensemble models 
- **Continuous SDM maps** – polished per-species maps of predicted suitability.
- **Binary SDM maps** – binary map of suitable/unsuitable areas, per species, using TSS threshold.

6.  `averageSuitability.R`
Average SDMs raster outputs by conservation status and scenarios to achieve a layer with the average suitability for the chosen category (conservation status or future scenario) for available species. Outputs are available in Zenodo as SDM_suitabilityAverages.zip.
