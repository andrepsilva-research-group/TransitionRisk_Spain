## Name: inputClimate.R ##
## Authors: Jorinde-M. Rieger ##
## Description: Applies functions to calculate climate input data (temperature and precipitation) for SDMRun.R
## for the ssp126 and ssp585 scenarios in various years##
## Date: August 5th 2025 ##

# Input variables -------------------------------------------
# Climatologies at high resolution for the earth’s land surface areas (CHELSEA) data provided by Karger et al. (2017 & 2021)
# (Data availability: https://chelsa-climate.org/downloads/)
# Define the file paths to the downloaded and processed data
basePathClim <- "data/CHELSA_gfdl-esm4_V.2.1"

yearsOrigin <- c("2011-2040", "2041-2070", "2071-2100") # original in time periods"2011-2040", "2041-2070", "2071-2100"
baseline_yearOrigin <- "1981-2010"

# Climate Models (GCMs)
models <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll")

# Map time periods to adapted years
yearsMapping <- setNames(years, yearsOrigin)

# Define the file paths
basePathClim <- "data/CHELSA_gfdl-esm4_V.2.1"

# Create environmental input Data (climate) as training landscapes (baseline) -------------------------------------------
trainingLandscapesClim <- list()
for (variable in variables) {
  raster <- load_baseline_clim(variable, baseline_yearOrigin)  # Load the raster
  
  # Format extent object and crop the baseline raster
  extent_crs <- sf::st_transform(extent_sf, crs = crs(raster))
  extent_sp <- terra::vect(extent_crs)
  raster_extent <- crop_mask_raster(raster, extent_sp)
  
  # Aggregate to target resolution
  input_resolution <- terra::res(raster_extent)[1]  # Assuming square cells, take the resolution of the first dimension
  aggregation_factor <- round(target_resolution / input_resolution)
  raster_agg <- terra::aggregate(raster_extent, fact = aggregation_factor, fun = mean)

  trainingLandscapesClim[[variable]] <- raster_agg
}

# Convert the list into a SpatRaster stack
trainingLandscapesClim <- terra::rast(trainingLandscapesClim)

# Rename layers to match variable names
names(trainingLandscapesClim) <- variables

# Save the trainingLandscape with the adapted baseline year
output_file <- file.path(outputPathLandscapes, paste0("trainingLandscapesClim_", baseline_year, "_5km.tif")) # global extent at 5 km resolution
terra::writeRaster(trainingLandscapesClim, output_file, overwrite = TRUE)

#plot(trainingLandscapesClim)

# Create environmental input Data (climate) as prediction landscapes (future scenarios) -------------------------------------------
# Register parallel backend
num_cores <- min(parallel::detectCores() - 1, 10)  # Use up to 10 cores
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# Create raster of averaged GCMs in a paralleled loop
foreach(scenario = scenarios, .combine = 'c', .packages = c("terra", "sf")) %:%
  foreach(yearOrigin = yearsOrigin, .combine = 'c') %dopar% {
    raster_list <- list() # Create a list for each scenario-year combination
    
    for (variable in variables) {
      # average climate models and save them
      raster <- average_climate_models(outputPathLandscapes, scenario, yearOrigin, variable)
      raster_list[[variable]] <- raster
    }
  }

# Stop the cluster
stopCluster(cl)

# Create prediction landscapes for climate variables
predictionLandscapesClim <- list()
for (scenario in scenarios) {
  for (yearOrigin in yearsOrigin) {
    raster_list <- list() # Create a list for each scenario-year combination
    for (variable in variables) {
      # Load averaged raster of climate models
      raster <- load_average_scenario_clim(outputPathLandscapes, scenario, yearOrigin, variable)
      
      # crop and mask the raster to the extent
      raster_extent <- crop_mask_raster(raster, extent_sp)
      
      # Aggregate to target resolution
      input_resolution <- terra::res(raster_extent)[1]  # Assuming square cells, take the resolution of the first dimension
      aggregation_factor <- round(target_resolution / input_resolution)
      raster_agg <- aggregate(raster_extent, fact = aggregation_factor, fun = mean)
      
      # Store the processed raster in the list
      raster_list[[variable]] <- raster_agg
    }
    
    # Convert the list of rasters into a SpatRaster stack
    predictionLandscapesClim[[paste0(scenario, "_", yearOrigin)]] <- terra::rast(raster_list)
  }
}

# Rename raster layers within each stack to match variable names
for (i in seq_along(predictionLandscapesClim)) {
  names(predictionLandscapesClim[[i]]) <- variables
}

# Save the climate prediction landscapes  with adapted years
for (scenario in scenarios) {
  for (yearOrigin in yearsOrigin) {
    year <- yearsMapping[yearOrigin] # Get the adapted year
    
    # Save the prediction landscape
    output_file <- file.path(outputPathLandscapes, paste0("predictionLandscapesClim_", scenario, "_", year, "_5km.tif")) # global extent at 5 km resolution
    terra::writeRaster(predictionLandscapesClim[[paste0(scenario, "_", yearOrigin)]], output_file, overwrite = TRUE)
  }
}

# Remove all data from memory/global Environment
rm(list = ls())
gc()
