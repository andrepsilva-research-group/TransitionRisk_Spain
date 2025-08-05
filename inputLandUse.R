## Name: inputLandUse.R ##
## Authors: Jorinde-M. Rieger ##
## Description: Applies functions to calculate land-use input data on defined extent for SDMRun.R for the ssp126 and ssp585 scenarios in various years ##
## Date: August 5th 2025 ##

# Input variables -------------------------------------------
# Download high-resolution (300m), global projections of land-use change under SSP-RCP scenarios derived from SEALS (Spatial Economic Allocation Landscape Simulator)
# by Johnson, J. A., & Thakrar, S., (2024) (Data availability: https://doi.org/10.5281/zenodo.14506512)
# Define the file paths to downloaded land-use data
basePathLandUse <- "data/stitched_lulc_esa_scenarios"

# Simplify and define ESA LULC types (39) to the 7 (SEALS) LULC types
# ESA LULC simplification scheme based on  Johnson, J. A., & Thakrar, S., (2024)
# (Code availability: https://github.com/jandrewjohnson/seals)
LULC_Types <- 1:7
LULC_Types_names <- c(
  "Urban",
  "Cropland",
  "Pasture_Grassland",
  "Forest",
  "Nonforest_vegetation",
  "Water",
  "Barren_other")

# To reclassify LandUse classes create a matrix that ranges from and to to define intervals
value_to_landUse <- c(
  190, 1,  # Urban
  10, 2, # Cropland
  11, 2,
  12, 2,
  20, 2,
  30, 2,
  40, 2,
  130, 3,  # Pasture/Grassland
  50, 4, # Forest
  60, 4,
  61, 4,
  62, 4, 
  70, 4, 
  71, 4,
  72, 4, 
  80, 4,
  81, 4, 
  82, 4,
  90, 4, 
  100, 4,
  151, 4,
  160, 4,
  170, 4,
  110, 5, # Non-forest vegetation
  120, 5,
  121, 5,
  122, 5, 
  140, 5,
  150, 5,
  152, 5,
  153, 5,
  180, 5,
  210, 6, # Water
  200, 7, # Barren or Other
  201, 7,
  202, 7,
  220, 7)
landUsematrix <- matrix(value_to_landUse, ncol = 2, byrow = TRUE )

# Create land-use input raster for the trainingLandscape -------------------------------------------
# Process baseline year
baseline_raster <- load_baseline_landUse(baseline_year)

# Format extent object and crop the baseline raster
extent_crs <- sf::st_transform(extent_sf, crs = crs(baseline_raster))
extent_sp <- terra::vect(extent_crs)
baseline_raster_extent <- crop_mask_raster(baseline_raster, extent_sp)

# Apply LULC type mapping and save the mapped baseline raster
mapped_baseline <- terra::classify(baseline_raster_extent, landUsematrix, include.lowest = TRUE)
output_file <- file.path(outputPathLandscapes, paste0("MappedLandUse_base_", baseline_year, "_", gsub(" ", "_", extent), ".tif"))
terra::writeRaster(mapped_baseline, output_file, overwrite = TRUE)
assign(paste0("MappedLandUse_base_", baseline_year, "_", gsub(" ", "_", extent)), mapped_baseline, envir = .GlobalEnv)

# Load the mapped raster stack for the baseline year
mapped_baseline <- load_mapped_baseline_landUse(outputPathLandscapes, baseline_year)

# Apply calculateRasterClass to the baseline raster to create raster classes for the land use types
trainingLandscapesLandUse <- calculateRasterClass(
  OriginalRaster = mapped_baseline,
  extent = extent_sp,
  target_resolution = target_resolution
)

# Replace LULC numbers with names
trainingLandscapesLandUse <- replace_numbers_with_names(trainingLandscapesLandUse, LULC_Types, LULC_Types_names)

# Save the processed baseline raster
output_file <- file.path(outputPathLandscapes, paste0("trainingLandscapesLandUse_", baseline_year, "_5km.tif")) # global extent at 5 km resolution
terra::writeRaster(trainingLandscapesLandUse, output_file, overwrite = TRUE)

# Create land-use input raster for the predictionLandscape -------------------------------------------
# Loop through the scenarios and years crop to the extent (global)
for (scenario_des in scenarios_des) {
  for (year in years) {
    raster <- load_scenario_landUse(scenario_des, year) # Load the raster
    
    # Crop and mask the raster
    raster_extent <- crop_mask_raster(raster, extent_sp)
    
    # Map the scenario name and save the mapped raster
    scenario <- scenario_name_mapping[scenario_des]
    output_file <- file.path(outputPathLandscapes, paste0("LandUse_", scenario, "_", year, "_", gsub(" ", "_", extent), ".tif"))
    terra::writeRaster(raster_extent, output_file, overwrite = TRUE)
    assign(paste0("LandUse_", scenario, "_", year, "_", gsub(" ", "_", extent)), raster_extent)
  }
}

library(foreach)
library(doParallel)
# Register parallel backend
num_cores <- min(parallel::detectCores() - 1, 10)  # use up to 10 cores
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# Loop through the years to create raster stacks and map land-use types
mapped_rasters <- foreach(year = years, .combine = 'c', .packages = c("terra", "sf")) %dopar% {
  # Create raster stacks for years
  raster_stack <- stack_rasters(year, scenarios, extent, outputPathLandscapes)
  
  # Load the raster stack for the year
  raster_stack <- file.path(outputPathLandscapes, paste0("LandUse_scenarioStack_", year, "_", gsub(" ", "_", extent), ".tif"))
  raster_stack <- terra::rast(raster_stack)
  
 # Apply LULC type mapping and save the mapped raster stack
  mapped_scenarios <- terra::classify(raster_stack, landUsematrix, include.lowest = TRUE)
  output_file <- file.path(outputPathLandscapes, paste0("MappedLandUse_scenarios_", year, "_", gsub(" ", "_", extent), ".tif"))
  terra::writeRaster(mapped_scenarios, output_file, overwrite = TRUE)

  # Return the mapped raster stack
  mapped_scenarios
}
stopCluster(cl)

# Load the mapped raster stacks for the target years
LULC_scenarios_list <- list()
for (year in years) {
  LULC_scenarios_list[[as.character(year)]] <- load_mapped_landUse(outputPathLandscapes, year, extent)
}

# Apply calculateRasterClass to the target year rasters
# Loop through the years to process each layer (scenario) (takes approx 8h per landscape on a global scale)
for (year in names(LULC_scenarios_list)) {
  target_raster <- LULC_scenarios_list[[year]] # Get the raster for the year
  
  # Ensure the raster has the same number of layers as the scenarios
  if (nlyr(target_raster) != length(scenarios)) {
    stop(paste("The raster for year", year, "does not have the same number of layers as the scenarios."))
  }
  
  # Dynamically assign scenario names to the raster layers
  scenario_rasters <- setNames(
    lapply(seq_along(scenarios), function(i) target_raster[[i]]),
    scenarios
  )
  
  # Create a list to store the processed rasters for this year
  processed_scenario_rasters <- list()
    for (scenario in names(scenario_rasters)) {
    scenario_raster <- scenario_rasters[[scenario]] # Get the raster for the scenario
    
    # Apply calculateRasterClass to classify the raster
    scenario_raster_classified <- calculateRasterClass(
      OriginalRaster = scenario_raster,  # Single-layer raster
      extent = extent_sp,
      target_resolution = target_resolution
    )
    
    # Save the processed raster
    output_file <- file.path(outputPathLandscapes, paste0("LandUseClass_", scenario, "_", year,"_", gsub(" ", "_", extent), ".tif"))
    terra::writeRaster(scenario_raster_classified, output_file, overwrite = TRUE)
   
    # Store the processed raster in the list
    processed_scenario_rasters[[scenario]] <- scenario_raster_classified
  }
  
  # Update the LULC_scenarios_list with the processed rasters
  LULC_scenarios_list[[year]] <- processed_scenario_rasters
}

# Replace LULC types numbers with names in the raster stacks
predictionLandscapesLandUse <- list()
for (scenario in scenarios) {
  for (year in years) {
    input_file <- file.path(outputPathLandscapes, paste0("LandUseClass_", scenario, "_", year, "_", gsub(" ", "_", extent), ".tif"))
    loaded_raster <- terra::rast(input_file) # Load the raster
    loaded_raster <- replace_numbers_with_names(loaded_raster, LULC_Types, LULC_Types_names)
    
    # Save the processed raster
    output_file <- file.path(outputPathLandscapes, paste0("predictionLandscapesLandUse_", scenario, "_", year, "_5km.tif")) # global extent at 5 km resolution
    terra::writeRaster(loaded_raster, output_file, overwrite = TRUE)
    
    # Store the raster in the list
    predictionLandscapesLandUse[[paste0(scenario, "_", year)]] <- loaded_raster
  }
}

# Remove all data from memory/global Environment
rm(list = ls())
gc()
