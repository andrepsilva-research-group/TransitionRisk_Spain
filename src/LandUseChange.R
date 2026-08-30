## Name: LandUseChange.R ##
## Authors: Jorinde-M. Rieger ##
## Description: Creates figures: Projected percentage change in land-use types & Spatial distribution of land-use cover change (%)
## for the ssp126 and ssp585 scenarios in various years ##
## Date: August 4th 2025 ##

# Settings & libraries -------------------------------------------
# set working directory to folder
source("src/server/libraries.R") # libraries
source("src/server/customFunctions.R") # functions

# Download high-resolution (300m), global projections of land-use change under SSP-RCP scenarios derived from SEALS (Spatial Economic Allocation Landscape Simulator)
# by Johnson, J. A., & Thakrar, S., (2024) (Data availability: https://doi.org/10.5281/zenodo.14506512)

# Define the file paths to the downloaded and processed data
basePathLandUse <- "data/stitched_lulc_esa_scenarios"
outputPathLandscapes <- "output/Landscapes"
if (!dir.exists(outputPathLandscapes)) {
  dir.create(outputPathLandscapes, recursive = TRUE)
}
outputPathLandUse <- "output/LandUseChange"
if (!dir.exists(outputPathLandUse)) {
  dir.create(outputPathLandUse, recursive = TRUE)
}

# Input variables -------------------------------------------
# Define input variables
scenarios <- c("ssp126", "ssp585") # define socio-economic pathways
scenarios_des <- c("rcp26_ssp1", "rcp85_ssp5") # scenario names in land-use raster
scenario_names <- c("SSP1-RCP2.6", "SSP5-RCP8.5") # define socio-economic pathways names

# Define a mapping for scenario names
scenario_name_mapping <- c(
  "rcp26_ssp1" = "ssp126",
  "rcp85_ssp5" = "ssp585")

# Define years
years <- c(2030, 2050, 2100)
baseline_year <- 2015

# Define the target resolution
#target_resolution <- 0.008333333 # 1km resolution
target_resolution <- 0.04166 # approx. 5km resolution

# Define extent: Spain and Portugal
extent <- "Spain and Portugal"
extent_name <- "Spain and Portugal"
extent_sf <- rnaturalearth::ne_countries(scale = "medium", country = c("Spain", "Portugal"), returnclass = "sf")

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

# To reclassify LULC types create a matrix
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

# Create land-use input raster for the trainingLandscape (baseline) -------------------------------------------
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
output_file <- file.path(outputPathLandscapes, paste0("trainingLandscapesLandUse_", baseline_year, "_", gsub(" ", "_", extent), "_5km.tif"))
terra::writeRaster(trainingLandscapesLandUse, output_file, overwrite = TRUE)
gc()

# Create land-use input raster for the predictionLandscape (future scenarios) -------------------------------------------
# Loop through the scenarios and years crop to the extent
for (scenario_des in scenarios_des) {
  for (year in years) {
    raster <- load_scenario_landUse(scenario_des, year) # Load the raster
    
    # Crop and mask the raster
    raster_extent <- crop_mask_raster(raster, extent_sp) # use extent_sp from baseline raster
    
    # Map the scenario name and save the mapped raster
    scenario <- scenario_name_mapping[scenario_des]
    output_file <- file.path(outputPathLandscapes, paste0("LandUse_", scenario, "_", year, "_", gsub(" ", "_", extent), ".tif"))
    terra::writeRaster(raster_extent, output_file, overwrite = TRUE)
    assign(paste0("LandUse_", scenario, "_", year, "_", gsub(" ", "_", extent)), raster_extent)
  }
}

# Load libraries for parallelization if needed
library(foreach)
library(doParallel)

# Register parallel backend
num_cores <- min(parallel::detectCores() - 1, 10)  # use up to 10 cores
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# Loop through the years to create raster stacks and map LULC types
mapped_rasters <- foreach::foreach(year = years, .combine = 'c', .packages = c("terra", "sf")) %dopar% {
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
    output_file <- file.path(outputPathLandscapes, paste0("LandUseClass_", scenario, "_", year,"_", gsub(" ", "_", extent), "_5km.tif"))
    terra::writeRaster(scenario_raster_classified, output_file, overwrite = TRUE)
    processed_scenario_rasters[[scenario]] <- scenario_raster_classified # store raster in a list
  }
  LULC_scenarios_list[[year]] <- processed_scenario_rasters
}

# Replace LULC types numbers with names in the raster stacks
predictionLandscapesLandUse <- list()
for (scenario in scenarios) {
  for (year in years) {
    input_file <- file.path(outputPathLandscapes, paste0("LandUseClass_", scenario, "_", year, "_", gsub(" ", "_", extent), "_5km.tif"))
    loaded_raster <- terra::rast(input_file) # Load the raster
    loaded_raster <- replace_numbers_with_names(loaded_raster, LULC_Types, LULC_Types_names)
    
    # Save the processed raster
    output_file <- file.path(outputPathLandscapes, paste0("predictionLandscapesLandUse_", scenario, "_", year,"_", gsub(" ", "_", extent), "_5km.tif"))
    terra::writeRaster(loaded_raster, output_file, overwrite = TRUE)
    predictionLandscapesLandUse[[paste0(scenario, "_", year)]] <- loaded_raster # store raster in a list
  }
}

# Example plot
#plot(predictionLandscapesLandUse$ssp126_2100)

# Calculate and create Figure: Projected percentage change in land-use types in extent -------------------------------------------
# Define consistent color palette for the scenarios
scenario_colors <- setNames(
  c("#1f77b4", "#ff7f0e"), scenario_names)

# Load mapped prediction LandscapeLandUse
mapped_scenarios <- list()
for (year in years) {
  mapped_scenarios[[as.character(year)]] <- load_mapped_landUse(outputPathLandscapes, year, extent)
}

# Calculate land use cover percentages and create figure
LandUseCover <- list()
  scenarios_percentages_df_list <- lapply(years, process_and_map_scenarios, use_continent = FALSE)
  scenarios_percentages_df <- do.call(rbind, scenarios_percentages_df_list)
  # safe percentages in the list
  LandUseCover <- scenarios_percentages_df

  # Clean and plot
  scenarios_percentages_df <- scenarios_percentages_df %>%
    dplyr::mutate(Scenario = gsub("^scenario_", "", Scenario)) %>%
    dplyr::mutate(Scenario = gsub("_\\d{4}$", "", Scenario)) %>%
    dplyr::mutate(Scenario = recode(Scenario, !!!setNames(scenario_names, scenarios)))
  
  scenarios_percentages_df_filtered <- scenarios_percentages_df %>%
    filter(landUse != "Water")
  
  LandUseChange_time_plot <- ggplot2::ggplot(scenarios_percentages_df_filtered, aes(x = time, y = value, color = Scenario, group = Scenario)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = scenario_colors, guide = guide_legend(direction = "horizontal")) +
    ggplot2::facet_wrap(~ landUse, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = paste0("Future Land Use Projections for ", extent),
         x = "Year", y = "Total Land Area (%)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = element_text(size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      strip.text = element_text(size = 12),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
    )
  print(LandUseChange_time_plot)
  
  ggplot2::ggsave(
    filename = file.path(outputPathLandUse, paste0("LandUseChange_time_", gsub(" ", "_", extent), ".png")),
    plot = LandUseChange_time_plot,
    width = 10,
    height = 6,
    dpi = 300
  )
gc()

# Load the training and prediction Landscapes -------------------------------------------
# Load training- and predictionLandscapesLandUse for defined extent
trainingLandscapesLandUse <- file.path(outputPathLandscapes, paste0("trainingLandscapesLandUse_", baseline_year,"_", gsub(" ", "_", extent), "_5km.tif"))
trainingLandscapesLandUse <- terra::rast(trainingLandscapesLandUse)

predictionLandscapesLandUse <- list()
for (scenario in scenarios) {
  for (year in years) {
    raster_path <- file.path(outputPathLandscapes, paste0("predictionLandscapesLandUse_", scenario, "_", year,"_", gsub(" ", "_", extent), "_5km.tif"))
    raster <- terra::rast(raster_path)
    predictionLandscapesLandUse[[paste0(scenario, "_", year)]] <- raster
  }
}

# Calculate spatial percentage changes of LULC types -------------------------------------------
percentage_change_rasters_list <- calculate_percentage_changes(
  base_year_raster = trainingLandscapesLandUse,
  target_year_rasters_list = predictionLandscapesLandUse,
  scenarios = scenarios,
  years = years
)
gc()

# Create Figure: Spatial distribution of land-use type cover change (%) in extent -------------------------------------------
# Create a custom color ramp with specified breakpoints
custom_color_ramp <- circlize::colorRamp2(c(-100, 0, 100), c("blue", "yellow", "red"))

# Define years included in the figure
years <- c(2030, 2050)

# Prepare plot grobs for all scenarios and years
for (class in names(trainingLandscapesLandUse)) {
  if (class == "Water") next # filter out water LULC type
  plot_grobs <- list()
  for (i in seq_along(scenarios)) {
    for (j in seq_along(years)) {
      plot_grobs[[paste0("row", i, "_col", j)]] <- gtable::gtable_trim(
        ggplotGrob(
          plot_landUse_spatialChanges(
            raster = percentage_change_rasters_list[[as.character(years[j])]][[scenarios[i]]][[class]],
            biome_geom = extent_sf,
            color_ramp = custom_color_ramp,
            fill_label = "Change in %  ",
            min_value = -100,
            max_value = 100
          ) + theme(legend.position = "none")
        )
      )
    }
  }

  # Column titles (years)
  title_grobs <- c(list(nullGrob()), lapply(years, function(year) textGrob(as.character(year), gp = gpar(fontsize = 15), just = "centre")))

  # Row labels (scenarios)
  row_labels <- lapply(seq_along(scenarios), function(i) textGrob(scenario_names[i], rot = 90, gp = gpar(fontsize = 15)))

  # Build layout matrix
  nrows <- length(scenarios) + 1
  ncols <- length(years) + 1
  layout_matrix <- matrix(NA, nrow = nrows, ncol = ncols)
  layout_matrix[1, ] <- 1:ncols  # Title row

  all_grobs <- title_grobs
  counter <- ncols + 1

  for (i in seq_along(scenarios)) {
    all_grobs <- append(all_grobs, list(row_labels[[i]]))
    for (j in seq_along(years)) {
      all_grobs <- append(all_grobs, list(plot_grobs[[paste0("row", i, "_col", j)]]))
    }
    layout_matrix[i + 1, ] <- counter:(counter + ncols - 1)
    counter <- counter + ncols
  }

  col_widths <- unit.c(unit(1, "cm"), rep(unit(1, "null"), length(years)))
  row_heights <- unit.c(unit(1, "null"), rep(unit(4, "null"), length(scenarios)))

  # Extract the legend from one of the plots
  example_plot <- plot_landUse_spatialChanges(
    raster = percentage_change_rasters_list[[as.character(years[1])]][[scenarios[1]]][[class]],
    biome_geom = extent_sf,
    color_ramp = custom_color_ramp,
    fill_label = "Change in %  ",
    min_value = -100,
    max_value = 100
  )
  shared_legend <- extract_legend(example_plot)

  # Arrange the plots, text and legend
  final_plot <- gridExtra::grid.arrange(
    grobs = all_grobs,
    layout_matrix = layout_matrix,
    widths = col_widths,
    heights = row_heights,
    bottom = shared_legend,
    top = textGrob(
      paste0("Land-Use Change of ", class, " Areas in ", extent_name),
      gp = gpar(fontsize = 20)
    )
  )

  # Save the figure
  ggplot2::ggsave(
    filename = file.path(outputPathLandUse, paste0("LandUseSpatialChanges_", class, "_", gsub(" ", "_", extent), ".png")),
    plot = final_plot,
    width = 10, height = 8, dpi = 300
  )
}
gc()
