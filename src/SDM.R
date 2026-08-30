## Name: SDM.R ##
## Author: Inês Silva & Jorinde-M. Rieger & André P. Silva ##
## Date: UPDATED 21th June 2026 ---- August 23th 2025 ##
## Description: builds and projects ensemble Species Distribution Models (SDMs)
## for multiple target species under current and future climate scenarios, using BIOMOD2



## Script overview:
##   STEP 1 – Load and crop environmental rasters to target biome.
##   STEP 2 – Define and prepare current and future climate scenarios.
##   STEP 3 – For each species, prepare presence/pseudo-absence data.
##   STEP 4 – Train individual models (RF, XGBOOST, ANN, MAXNET/MAXENT).
##   STEP 5 – Build ensemble models and evaluate them.
##   STEP 6 – Project models to current and future conditions.
##   STEP 7 – Save continuous and binary suitability rasters for each scenario.
##   STEP 8 – Record processing time per species and export results.


# !! NAVIGATION WARNINGS !! ----------------------------------------------------

# This is a highly sensible function! Here are some good practices to make sure
# we don't get errors running it:
# (1) do not work inside a One Drive folder;
# (2) avoid saving the outputs using long paths (MAXNET does not deal well with them)
# (3) some minor chnages were done to André's function, namely adding
# the path to the maxent folder as a function argument (might not be necessary in
# the future if we stick with MAXNET but I am still leaving it here)
# (4) function stops if one species has an error SOLUTION? use TryCatch()

# Settings & libraries ---------------------------------------------------------
source("src/libraries.R") # libraries
source("src/customFunctions2.R") # functions

##########
# STEP 1 # Load target species, full environmental training raster & occurrences (done only once!)
##########

# target species
species_table <- read.csv("./data/iberian_threatened_mammals_taxonomy.csv",
                          stringsAsFactors = FALSE)

# environmental variables to keep in landscapes & models
vars_to_keep <- c("Elevation",
                  "Urban",
                  "Cropland",
                  "Pasture_Grassland",
                  "Forest",
                  "Nonforest_vegetation",
                  "Water",
                  "Barren_other",
                  #"bio1", # Annual Mean Temperature
                  #"bio10", # Mean Temp. Warmest Quarter
                  "bio11", # Mean Temp. Coldest Quarter
                  "bio12"#, # Annual Percipitation
                  #"bio16", # Percipitation Wettest Quarter
                  #"bio17" # Percipitation Driest Quarter
)

# training raster
myExpl_full <- rast("./data/sdm/Landscapes/trainingLandscapes_2015_5km.tif")
myExpl_full <- myExpl_full[[vars_to_keep]]
invisible(gc())

# species occurrences
speciesData <- read.csv("./data/GBIF_occurrences_mammals_2026-06-07.csv")
speciesData$species <- gsub(" ", ".", speciesData$species)
invisible(gc())

# crop FUTURE landscapes dynamically 
years <- c(2030, 2050)
future_scenario <- c("ssp126", "ssp585")


# Define extent to crop predicitionLandscapes (Spain and Portugal)
biomes <- c("IberianPeninsula")

target_biome <- "IberianPeninsula"


for (target_biome in biomes) {
  
  message("====================================")
  message("Starting biome: ", target_biome)
  message("====================================")
  
  sdm_output_dir <- "E:/BankOfSpain_june2026/SDMs_mammals"
  dir.create(sdm_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  ############
  # STEP 2.1 # load and crop biome extent
  ############

  extent_sf <- rnaturalearth::ne_countries(scale = "medium", country = c("Spain", "Portugal"), returnclass = "sf")
  extent_crs <- sf::st_transform(extent_sf, crs = crs(myExpl_full)) # Ensure CRS consistency
  extent_sp <- terra::vect(extent_crs) # Convert the sf to a spatial object
  invisible(gc())
  
  ############
  # STEP 2.2 # crop CURRENT landscape to the target biome
  ############
  # crop to biome
  myExpl_current <- crop(
    rast("./data/sdm/Landscapes/trainingLandscapes_2015_5km.tif"),
    extent_sp)
  # select approapriate environmental variables 
  myExpl_current <- myExpl_current[[vars_to_keep]]
  invisible(gc())
  
  ############
  # STEP 2.3 # crop FUTURE landscapes
  ############
  myExpl_future <- list()
  
  # dynamically crop teh landscapes rasters per year and scenario by the biome
  for (scen in future_scenario) {
    for (yr in years) {
      name <- paste0(scen, "_", yr, "_", target_biome)
      message(paste0("Cropping Landscapes for ", name))
      path <- sprintf("./data/sdm/Landscapes/predictionLandscapes_%s_%s_5km.tif", scen, yr)
      # Load & crop
      r <- crop(rast(path), extent_sp)
      # Select same variables as current
      myExpl_future[[name]] <- r[[vars_to_keep]]
      rm(name, path, r)
      invisible(gc())
    }
  }
  
  
  # get results into lists (IMPORTANT STEP! Do not skip!!)
  myExplCurrent <- list(myExpl_current)
  myExplFuture  <- myExpl_future
  
  ############
  # STEP 4.4 # run SDMs for this biome
  ############
  # SDMensembleMultiSpecies() – this function automates the full SDM workflow for
  # each species: filtering and formatting GBIF data, generating pseudo-absences, 
  # fitting multiple algorithms (RF, XGBoost, ANN, MAXNET), building ensemble models, 
  # and projecting them to current and future climate scenarios. 
  # It saves all evaluation metrics, variable importance tables, and raster predictions 
  # (continuous and binary) into the specified output folder.
  
  # store timing info
  species_times <- list()   
  
  # select taregt species for specific biome
  target_species <- species_table %>% 
    # replace spaces with . to match biomod2
    mutate(sci_name = gsub(" ", ".", sci_name)) %>% 
    # keep only names
    dplyr::pull(sci_name) %>% 
    unique()
  
  #sp <- target_species[[2]]  
  
  SDM_NatPoke <- lapply(target_species, function(sp) {
    tryCatch({
      # measure start time
      start_time <- Sys.time()
      
      result <- SDMensembleMultiSpecies(
        targetSpecies = sp,
        speciesData = speciesData,
        myExpl_full = myExpl_full,
        myExplCurrent = myExplCurrent,
        myExplFuture = myExplFuture,
        extent = target_biome,
        output_folder = "E:/BankOfSpain_june2026/SDMs_mammals",
        maxent_source = "C:/Users/maria/Desktop/maxent/maxent/maxent.jar",
        ncoresToUse = 6)
      
      # clean temporary files to avoid filling RAM up
      # message("Cleaning temporary raster files")
      # terra::tmpFiles(remove = TRUE)
      # invisible(gc())
      
      # measure time gone by
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      message("✅ Finished ", sp, " in ", round(elapsed, 2), " minutes")
      # store timing
      species_times[[sp]] <<- elapsed
      return(result)
      
    }, error = function(e) {
      message(paste("⚠️ Skipping", sp, "due to error:", e$message))
      species_times[[sp]] <<- NA   # store NA if failed
      return(NULL)
    })
  })
  
  ############
  # STEP 4.5 # save timming results
  ############
  
  timing_df <- data.frame(Species = names(species_times),
                          Time_difference_mins = as.numeric(species_times),
                          stringsAsFactors = FALSE)
  
  write_xlsx(timing_df, file.path(sdm_output_dir,
                                  paste0(target_biome, "_species_times.xlsx")))
  
  rm(extent_sf, extent_crs, extent_sp, myExpl_current, myExpl_future,
     myExplCurrent, myExplFuture, species_times, SDM_NatPoke)
  gc()
}



################################################################################

##########
# STEP 1 # Load target species, full environmental training raster & occurrences (done only once!)
##########

# target species
species_table <- read.csv("./data/iberian_threatened_birds_taxonomy.csv",
                          stringsAsFactors = FALSE)

# environmental variables to keep in landscapes & models
vars_to_keep <- c("Elevation",
                  "Urban",
                  "Cropland",
                  "Pasture_Grassland",
                  "Forest",
                  "Nonforest_vegetation",
                  "Water",
                  "Barren_other",
                  #"bio1", # Annual Mean Temperature
                  #"bio10", # Mean Temp. Warmest Quarter
                  "bio11", # Mean Temp. Coldest Quarter
                  "bio12"#, # Annual Percipitation
                  #"bio16", # Percipitation Wettest Quarter
                  #"bio17" # Percipitation Driest Quarter
)

# training raster
myExpl_full <- rast("./data/sdm/Landscapes/trainingLandscapes_2015_5km.tif")
myExpl_full <- myExpl_full[[vars_to_keep]]
invisible(gc())

# species occurrences
speciesData <- read.csv("./data/GBIF_occurrences_birds_2026-06-08.csv")
speciesData$species <- gsub(" ", ".", speciesData$species)
invisible(gc())

# crop FUTURE landscapes dynamically 
years <- c(2030, 2050)
future_scenario <- c("ssp126", "ssp585")


# Define extent to crop predicitionLandscapes (Spain and Portugal)
biomes <- c("IberianPeninsula")

target_biome <- "IberianPeninsula"


for (target_biome in biomes) {
  
  message("====================================")
  message("Starting biome: ", target_biome)
  message("====================================")
  
  sdm_output_dir <- "E:/BankOfSpain_june2026/SDMs_birds"
  dir.create(sdm_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  ############
  # STEP 2.1 # load and crop biome extent
  ############
  
  extent_sf <- rnaturalearth::ne_countries(scale = "medium", country = c("Spain", "Portugal"), returnclass = "sf")
  extent_crs <- sf::st_transform(extent_sf, crs = crs(myExpl_full)) # Ensure CRS consistency
  extent_sp <- terra::vect(extent_crs) # Convert the sf to a spatial object
  invisible(gc())
  
  ############
  # STEP 2.2 # crop CURRENT landscape to the target biome
  ############
  # crop to biome
  myExpl_current <- crop(
    rast("./data/sdm/Landscapes/trainingLandscapes_2015_5km.tif"),
    extent_sp)
  # select approapriate environmental variables 
  myExpl_current <- myExpl_current[[vars_to_keep]]
  invisible(gc())
  
  ############
  # STEP 2.3 # crop FUTURE landscapes
  ############
  myExpl_future <- list()
  
  # dynamically crop teh landscapes rasters per year and scenario by the biome
  for (scen in future_scenario) {
    for (yr in years) {
      name <- paste0(scen, "_", yr, "_", target_biome)
      message(paste0("Cropping Landscapes for ", name))
      path <- sprintf("./data/sdm/Landscapes/predictionLandscapes_%s_%s_5km.tif", scen, yr)
      # Load & crop
      r <- crop(rast(path), extent_sp)
      # Select same variables as current
      myExpl_future[[name]] <- r[[vars_to_keep]]
      rm(name, path, r)
      invisible(gc())
    }
  }
  
  
  # get results into lists (IMPORTANT STEP! Do not skip!!)
  myExplCurrent <- list(myExpl_current)
  myExplFuture  <- myExpl_future
  
  ############
  # STEP 4.4 # run SDMs for this biome
  ############
  # SDMensembleMultiSpecies() – this function automates the full SDM workflow for
  # each species: filtering and formatting GBIF data, generating pseudo-absences, 
  # fitting multiple algorithms (RF, XGBoost, ANN, MAXNET), building ensemble models, 
  # and projecting them to current and future climate scenarios. 
  # It saves all evaluation metrics, variable importance tables, and raster predictions 
  # (continuous and binary) into the specified output folder.
  
  # store timing info
  species_times <- list()   
  
  # select taregt species for specific biome
  target_species <- species_table %>% 
    # replace spaces with . to match biomod2
    mutate(sci_name = gsub(" ", ".", sci_name)) %>% 
    # keep only names
    dplyr::pull(sci_name) %>% 
    unique()
  
  #sp <- target_species[[2]]  
  
  SDM_NatPoke <- lapply(target_species, function(sp) {
    tryCatch({
      # measure start time
      start_time <- Sys.time()
      
      result <- SDMensembleMultiSpecies(
        targetSpecies = sp,
        speciesData = speciesData,
        myExpl_full = myExpl_full,
        myExplCurrent = myExplCurrent,
        myExplFuture = myExplFuture,
        extent = target_biome,
        output_folder = "E:/BankOfSpain_june2026/SDMs_birds",
        maxent_source = "C:/Users/maria/Desktop/maxent/maxent/maxent.jar",
        ncoresToUse = 6)
      
      # clean temporary files to avoid filling RAM up
      # message("Cleaning temporary raster files")
      # terra::tmpFiles(remove = TRUE)
      # invisible(gc())
      
      # measure time gone by
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      message("✅ Finished ", sp, " in ", round(elapsed, 2), " minutes")
      # store timing
      species_times[[sp]] <<- elapsed
      return(result)
      
    }, error = function(e) {
      message(paste("⚠️ Skipping", sp, "due to error:", e$message))
      species_times[[sp]] <<- NA   # store NA if failed
      return(NULL)
    })
  })
  
  ############
  # STEP 4.5 # save timming results
  ############
  
  timing_df <- data.frame(Species = names(species_times),
                          Time_difference_mins = as.numeric(species_times),
                          stringsAsFactors = FALSE)
  
  write_xlsx(timing_df, file.path(sdm_output_dir,
                                  paste0(target_biome, "_species_times.xlsx")))
  
  rm(extent_sf, extent_crs, extent_sp, myExpl_current, myExpl_future,
     myExplCurrent, myExplFuture, species_times, SDM_NatPoke)
  gc()
}
