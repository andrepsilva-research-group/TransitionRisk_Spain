## Name: environemtalScenarios.R ##
## Author: Jorinde-M. Rieger & André P. Silva ##
## Description: prepares calibration and prediction landscapes for SDMs based on selected variables
## and environmental scenarios ##
## Date: August 6th 2025 ##
## Note: this script was not tested and needs to be double checked before running the final SDMs
## the output landscapes of this script are saved at. the variables selected at the end seem to include
## too many climate variables compared to what we initially decided
## https://drive.google.com/drive/folders/1II-rlvB7RQNNez3sSPAy6gKxi1fHQH3k?usp=sharing

# Settings & libraries -----------------------------------------------------------------
# set working directory for maxent.jar file
source("src/libraries.R") # libraries
# When using the pipeline for the first time run taxaOcccurence.R and adapt species and user-login for GBIF Database:
#source("NatPoKe/src/TaxaOccurence.R") # downloads taxa occurences from GBIF Database
source("src/customFunctions2.R") # functions

# Define scenarios and environmental variables -----------------------------------------------------------------
scenarios <- c("ssp126", "ssp585") # define socio-economic pathways
scenarios_des <- c("rcp26_ssp1", "rcp85_ssp5") # scenario names in land-use raster
scenario_names <- c("SSP1-RCP2.6", "SSP5-RCP8.5") # define socio-economic pathways names

# Define a mapping for scenario names
scenario_name_mapping <- c(
  "rcp26_ssp1" = "ssp126",
  "rcp85_ssp5" = "ssp585")

# Define climatologies
variables <- c("bio1", # mean annual air temperature
               "bio10", # mean daily mean air temperatures of the warmest quarter
               "bio11", # mean daily mean air temperatures of the coldest quarter
               "bio12", # mean annual precipitation amount
               "bio16", # mean monthly precipitation amount of the wettest quarter
               "bio17") # mean monthly precipitation amount of the driest quarter

# Define years
years <- c(2030, 2050, 2100)
years <- c(2030, 2050)
baseline_year <- 2015

# Define the target resolution
target_resolution <- 0.04166 # approx. 5km resolution

# Global extent as initial extent for SDMs
extent <- "Global Terrestrial"
extent_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Define output path
outputPathLandscapes <- "output/Landscapes"
if (!dir.exists(outputPathLandscapes)) {
  dir.create(outputPathLandscapes, recursive = TRUE)
}

# When using the pipeline for the first time run on global extent:-----------------------------------------------------------------
# Note: that this is only needed when running the pipeline for the first time, as the landscapes are already prepared. If you want to run it 
# for the first time you need
# Format training and prediction Landscape
source("src/inputClimate.R") # format and reads input climate raster landscapes, adapt: scenarios, years & variables
source("src/inputLandUse.R") # format and reads input land-use raster landscapes, adapt: scenarios, years & variables
source("src/inputElev.R") # format and reads input land-use raster landscapes, adapt: scenarios, years & variables

# Define file paths and load training landscapes
trainingLandscapesClim <- file.path(outputPathLandscapes, paste0("trainingLandscapesClim_", baseline_year, "_5km.tif"))
trainingLandscapesLandUse <- file.path(outputPathLandscapes,paste0("trainingLandscapesLandUse_",  baseline_year, "_5km.tif"))
trainingLandscapesElev <- file.path(outputPathLandscapes, paste0("trainingLandscapesElev_",  baseline_year, "_5km.tif"))
trainingLandscapesClim <- terra::rast(trainingLandscapesClim)
trainingLandscapesLandUse <- terra::rast(trainingLandscapesLandUse)
trainingLandscapesElev <- terra::rast(trainingLandscapesElev)

# Resample the extent of the training landscapes to the land-use training Landscape
trainingLandscapesClim <- terra::resample(trainingLandscapesClim, trainingLandscapesLandUse)
trainingLandscapesElev <- terra::resample(trainingLandscapesElev, trainingLandscapesLandUse)

# Merge the climate, elevation and land-use rasters
trainingLandscapes <- c(trainingLandscapesElev, trainingLandscapesLandUse, trainingLandscapesClim)

# Save the merged training landscape
output_file <- file.path(outputPathLandscapes, paste0("trainingLandscapes_", baseline_year, "_5km.tif"))
writeRaster(trainingLandscapes, output_file, overwrite = TRUE)

# Create a list to store the merged prediction landscapes in parallelization
# Register parallel backend
num_cores <- min(parallel::detectCores() - 1, 10)  # Use up to 10 cores
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# Create an empty list to store merged prediction landscapes
predictionLandscapes <- list()

# Parallelized loop using foreach
foreach::foreach(scenario = scenarios, .combine = 'c', .packages = c("terra", "sf")) %:%
  foreach::foreach(year = years, .combine = 'c') %dopar% {
    # Define file paths for climate and land-use prediction landscapes & load them in the environment, if needed
    predictionLandscapesClim <- file.path(outputPathLandscapes, paste0("predictionLandscapesClim_", scenario, "_", year, "_5km.tif"))
    predictionLandscapesLandUse <- file.path(outputPathLandscapes, paste0("predictionLandscapesLandUse_", scenario, "_", year, "_5km.tif"))
    predictionLandscapesElev <- file.path(outputPathLandscapes, paste0("predictionLandscapesElev_", scenario, "_", year, "_5km.tif"))
    predictionLandscapesClim <- terra::rast(predictionLandscapesClim)
    predictionLandscapesLandUse <- terra::rast(predictionLandscapesLandUse)
    predictionLandscapesElev <- terra::rast(predictionLandscapesElev)
    
    # Ensure CRS, extent, and resolution consistency
    predictionLandscapesClim <- terra::resample(predictionLandscapesClim, predictionLandscapesLandUse)
    predictionLandscapesElev <- terra::resample(predictionLandscapesElev, predictionLandscapesLandUse)
    
    # Merge the climate and land-use rasters
    merged_prediction <- c(predictionLandscapesElev, predictionLandscapesLandUse, predictionLandscapesClim)
    
    # Save the merged prediction landscape
    output_file <- file.path(outputPathLandscapes, paste0("predictionLandscapes_", scenario, "_", year, "_5km.tif"))
    writeRaster(merged_prediction, output_file, overwrite = TRUE)
    
    # Store the merged prediction landscape in the list
    predictionLandscapes[[paste0(scenario, "_", year)]] <- merged_prediction
  }

# Stop the cluster
stopCluster(cl)

# Load formatted training and prediction Landscape -----------------------------------------------------------------
# Load formatted trainingLandscapes at global extent
trainingLandscapes <- file.path(outputPathLandscapes, paste0("trainingLandscapes_", baseline_year, "_5km.tif")) # with Antarctica
trainingLandscapes <- terra::rast(trainingLandscapes)







# Define extent for the predicitionLandscapes
# Put Biomes together in one landscape?
# Tropical Biome
extent <- "Tropical Biome"
extent_name <- "Tropical & Subtropical Moist Broadleaf Forests" # full name of the biome
extent_sf <- load_biome(extent_name)
extent_crs <- sf::st_transform(extent_sf, crs = crs(trainingLandscapes)) # Ensure CRS consistency
extent_sp <- terra::vect(extent_crs) # Convert the sf to a spatial object

# Boreal Biome
extent <- "Boreal Biome"
extent_name <- "Boreal Forests/Taiga" # full name of the biome
extent_sf <- load_biome(extent_name)
extent_crs <- sf::st_transform(extent_sf, crs = crs(trainingLandscapes)) # Ensure CRS consistency
extent_sp <- terra::vect(extent_crs) # Convert the sf to a spatial object

# Load predictionLandscapes and crop to defined extent
predictionLandscapes <- list()
# Loop through scenarios and years
for (scenario in scenarios) {
  for (year in years) {
    raster_path <- file.path(outputPathLandscapes, paste0("predictionLandscapes_", scenario, "_", year, "_5km.tif"))# with Antarctica
    raster <- terra::rast(raster_path)
    raster <- crop_mask_raster(raster, extent_sp)
    predictionLandscapes[[paste0(scenario, "_", year)]] <- raster
  }
}
rm(raster)

# When using the pipeline for the first time: check variable correlation and select suitable variables for the landscapes -----------------------------------------------------------------
# Calculate the correlation matrix for all layers in the trainingLandscapes with Pearson's correlation coefficient
cor_matrix <- terra::layerCor(trainingLandscapes, fun = "cor", use = "complete.obs", maxcell = 0.5*ncell(trainingLandscapes), na.rm = TRUE) # pearson correlation coefficient
write.csv(cor_matrix, file = "output/cor_matrix.csv", row.names = TRUE)
cor_mat <- cor_matrix$correlation

# Plot the correlation coefficients as percentages
png(filename = file.path(outputPathLandscapes, "CorrelationMatrix.png"),  width = 2000, height = 1500, res = 300)
corrplot::corrplot.mixed(
  cor_mat, tl.pos = 'lt', tl.cex = 0.6, number.cex = 0.5, addCoefasPercent = TRUE, tl.col = "black")
dev.off()

# Find highly correlated pairs (absolute correlation > 0.7) and select suitable variables
high_cor_pairs <- which(abs(cor_mat) > 0.7 & upper.tri(cor_mat), arr.ind = TRUE)

# Print pairs for selection process
for(i in seq_len(nrow(high_cor_pairs))) {
  cat(
    rownames(cor_mat)[high_cor_pairs[i, 1]], "and",
    colnames(cor_mat)[high_cor_pairs[i, 2]],
    "correlation:",
    round(cor_mat[high_cor_pairs[i, 1], high_cor_pairs[i, 2]], 2), "\n"
  )
}

# Filter suitable variables for the landscapes -----------------------------------------------------------------
# Remove highly correlated variables from training and prediction Landscapes
vars_to_remove <- c("Elevation", "bio10", "bio11", "bio16", "bio17")

# Subset trainingLandscapes to remove unwanted variables
trainingLandscapes <- trainingLandscapes[[!names(trainingLandscapes) %in% vars_to_remove]]

# Subset each predictionLandscapes raster to remove unwanted variables
for (name in names(predictionLandscapes)) {
  predictionLandscapes[[name]] <- predictionLandscapes[[name]][[!names(predictionLandscapes[[name]]) %in% vars_to_remove]]
}
