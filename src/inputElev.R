## Name: inputElev.R ##
## Authors: Jorinde-M. Rieger ##
## Description: Processes elevation data as input for training and prediction landscapes for SDMRun.R ##
## Date: August 5th 2025 ##

# Input variables -------------------------------------------
# Elevation data underlying WorldClim V2 by Fick and Hijmans (2017)
# (Data availability: https://www.worldclim.org/data/worldclim21.html)
# Define the file paths to the downloaded file
basePathElev <- "data/wc2.1_30s_elev/wc2.1_30s_elev.tif" # Path to the elevation raster downloaded from WorldClim

# Load the elevation raster -------------------------------------------
elevation_raster <- terra::rast(basePathElev)

# Format extent object and crop the baseline raster
extent_crs <- sf::st_transform(extent_sf, crs = crs(elevation_raster))
extent_sp <- terra::vect(extent_crs)
raster_extent <- crop_mask_raster(elevation_raster, extent_sp)

# Aggregate to target resolution
input_resolution <- terra::res(raster_extent)[1]  # Assuming square cells, take the resolution of the first dimension
aggregation_factor <- round(target_resolution / input_resolution)
raster_agg <- terra::aggregate(raster_extent, fact = aggregation_factor, fun = mean)

# Rename the layer to "Elevation"
names(raster_agg) <- "Elevation"

# Create training landscape (baseline) -------------------------------------------
trainingLandscapesElev <- raster_agg
output_file <- file.path(outputPathLandscapes, paste0("trainingLandscapesElev_", baseline_year, "_5km.tif")) # global extent at 5 km resolution
terra::writeRaster(trainingLandscapesElev, output_file, overwrite = TRUE)

# Create prediction landscapes (future scenarios) -------------------------------------------

predictionLandscapesElev <- list()
for (scenario in scenarios) {
  for (year in years) {
    # Use the same elevation raster for all scenarios and years
    predictionLandscapesElev[[paste0(scenario, "_", year)]] <- raster_agg
    
    # Save the elevation raster for this scenario and year
    output_file <- file.path(outputPathLandscapes, paste0("predictionLandscapesElev_", scenario, "_", year, "_5km.tif")) # global extent at 5 km resolution
    terra::writeRaster(predictionLandscapesElev[[paste0(scenario, "_", year)]], output_file, overwrite = TRUE)
  }
}

# Remove all data from memory/global Environment
rm(list = ls())
gc()
