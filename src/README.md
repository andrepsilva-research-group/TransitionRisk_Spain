# 🗂️ src folder structure

### 🧰 Step 1 - Setup the environment

1.  `libraries.R`  
Installs and loads all necessary R packages used throughout the pipeline.  
   *Run this script first to make sure all dependencies are in place.*

2.  `customFunctions2.R`  
Loads all custom helper functions used across different parts of the pipeline.

### 🐾 Step 2 - Get taxa occurrences

3. `taxaOccurrences.R`  
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
  
### 🗺️ Step 3 - Run species distribution models

4.  `SDM_MIS.R`- This is the main pipeline script. It:

    -   Formats and prepares species occurrence and environmental data.
    -   Builds multiple SDM models per species (using biomod2 functions).
    -   Creates ensemble models.
    -   Projects SDMs for both current and future environmental conditions.

Requirements: A folder containing training and prediction landscapes (raster files .tif). Download [here](https://drive.google.com/drive/folders/1ik_lyTc5CdcX6HUMIfU7cKGgvi0UqzwH?usp=sharing).

### 📊 Step 4 - Analyse outputs and Visualise results

5.  `SDMfigures.R`  
Processes SDM results and produces a range of **summary tables and figures** per species and one final figure combining all.

**Outputs include:**
- **Presence maps** – using filtered GBIF occurrences per species.
- **Variable importance table** – for ensemble models and multiple species
- **Evaluation plots** – performance metrics for ensemble models 
- **Continuous SDM maps** – polished per-species maps of predicted suitability.
- **Binary SDM maps** – binary map of suitable/unsuitable areas, per species, using TSS threshold.
- **Hotspot maps** – overlap of all species’ binary outputs to obtain hotspot areas where threathned species are more likely to occur areas.

-------------------

5. `mammalTraits_MIS.R` - Fetches and compiles trait data for the selected mammal species from existing trait databases.

Only needed if the metaRange module is being used.
If not required, this script can safely be removed later.

