## Name: ClimateChange.R ##
## Authors: Jorinde-M. Rieger ##
## Description: Create figures: Projected changes in climatologies (bio1 & 12) for extent & Projected spatial changes in climatologies (bio1 & 12) for extent
## for the ssp126 and ssp585 scenarios for 2030, 2050, (2100) ##
## Date: August 4th 2025 ##

# Settings & libraries -------------------------------------------
# set working directory to folder
source("src/libraries.R") # libraries
source("src/customFunctions.R") # functions

# Climatologies at high resolution for the earth’s land surface areas (CHELSEA) data provided by Karger et al. (2017 & 2021)
# (Data availability: https://chelsa-climate.org/downloads/)
# Define the file paths to the downloaded and processed data
basePathClim <- "data/CHELSA_Data"
outputPathLandscapes <- "output/Landscapes"
if (!dir.exists(outputPathLandscapes)) {
  dir.create(outputPathLandscapes, recursive = TRUE)
}

# Define output path
outputPathClimate <- "output/ClimateChange"
if (!dir.exists(outputPathClimate)) {
  dir.create(outputPathClimate, recursive = TRUE)
  
# Input variables -------------------------------------------
# Define input variables
scenarios <- c("ssp126", "ssp585") # define socio-economic pathways
scenario_names <- c("SSP1-RCP2.6", "SSP5-RCP8.5") # define socio-economic pathways names

# Define climatologies
variables <- c("bio1", # mean annual air temperature
               "bio12") # annual precipitation amount
variable_names <- c("Temperature (bio1)", 
                    "Precipitation (bio12)")
y_labels <- c("Mean annual air temperature (bio1) (°C)", 
              "Annual precipitation amount (bio12)(kg m-2 year-1)")
value_units <- c("°C", 
                "kg m-2 year-1")

# Define years
years <- c(2030, 2050, 2100)
baseline_year <- 2015
yearsOrigin <- c("2011-2040", "2041-2070", "2071-2100") # original in time periods"2011-2040", "2041-2070", "2071-2100"
baseline_yearOrigin <- "1981-2010"
# Map time periods to adapted years
yearsMapping <- setNames(years, yearsOrigin)

# Climate Models (GCMs)
models <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll")

# Spain and Portugal
extent <- "Spain and Portugal"
extent_sf <- rnaturalearth::ne_countries(scale = "medium", country = c("Spain", "Portugal"), returnclass = "sf")

# Define the target resolution
target_resolution <- 0.04166 # approx. 5km resolution

# Define the file paths
basePathClim <- "data/CHELSA_gfdl-esm4_V.2.1"
outputPathLandscapes <- "output/Landscapes"
if (!dir.exists(outputPathLandscapes)) {
  dir.create(outputPathLandscapes, recursive = TRUE)
}

# Define output path
outputPathClimate <- "output/ClimateChange"
if (!dir.exists(outputPathClimate)) {
  dir.create(outputPathClimate, recursive = TRUE)
}

# Create climate rasters as training landscapes (baseline) -------------------------------------------
trainingLandscapesClim <- list()
for (variable in variables) {
  raster <- load_baseline_clim(variable, baseline_yearOrigin) # Load the raster
  
  # Format extent object and crop the baseline raster
  extent_crs <- sf::st_transform(extent_sf, crs = crs(raster))
  extent_sp <- terra::vect(extent_crs)
  raster_extent <- crop_mask_raster(raster, extent_sp)
  
  # Aggregate to target resolution
  input_resolution <- terra::res(raster_extent)[1]  # Assuming square cells, take the resolution of the first dimension
  aggregation_factor <- round(target_resolution / input_resolution)
  raster_agg <- aggregate(raster_extent, fact = aggregation_factor, fun = mean)
  
  trainingLandscapesClim[[variable]] <- raster_agg
}

# Convert the list into a SpatRaster stack
trainingLandscapesClim <- terra::rast(trainingLandscapesClim)

# Rename layers to match variable names
names(trainingLandscapesClim) <- variables

# Save the trainingLandscape with the adapted baseline year
output_file <- file.path(outputPathLandscapes, paste0("trainingLandscapesClim_", baseline_year, "_", gsub(" ", "_", extent), "_5km.tif"))
terra::writeRaster(trainingLandscapesClim, output_file, overwrite = TRUE)

#plot(trainingLandscapesClim)

# Create climate rasters as prediction landscapes (future scenarios) -------------------------------------------
# Register parallel backend
num_cores <- min(parallel::detectCores() - 1, 10)  # Use up to 10 cores
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# Create raster of averaged GCMs in a paralleled loop
foreach::foreach(scenario = scenarios, .combine = 'c', .packages = c("terra", "sf")) %:%
  foreach::foreach(yearOrigin = yearsOrigin, .combine = 'c') %dopar% {
    raster_list <- list()  # Create a list for each scenario-year combination
    
    for (variable in variables) {
      # Average climate models and save them
      raster <- average_climate_models(basePathClim, scenario, yearOrigin, variable)
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
      
      # Format extent object and crop the baseline raster
        extent_crs <- sf::st_transform(extent_sf, crs = crs(raster))
        extent_sp <- terra::vect(extent_crs)
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
    output_file <- file.path(outputPathLandscapes, paste0("predictionLandscapesClim_", scenario, "_", year,"_", gsub(" ", "_", extent), "_5km.tif"))
    terra::writeRaster(predictionLandscapesClim[[paste0(scenario, "_", yearOrigin)]], output_file, overwrite = TRUE)
  }
}

#plot(predictionLandscapesClim[["ssp126_2071-2100"]])

# Load the training and prediction Landscapes for climate variables -------------------------------------------
trainingLandscapesClim <- file.path(outputPathLandscapes, paste0("trainingLandscapesClim_", baseline_year, "_", gsub(" ", "_", extent),"_5km.tif"))
trainingLandscapesClim <- terra::rast(trainingLandscapesClim)

# Load predictionLandscapesClim
predictionLandscapesClim <- list()
# Loop through scenarios and years
for (scenario in scenarios) {
  for (year in years) {
    raster_path <- file.path(outputPathLandscapes, paste0("predictionLandscapesClim_", scenario, "_", year, "_", gsub(" ", "_", extent),"_5km.tif"))
    raster <- terra::rast(raster_path)
    predictionLandscapesClim[[paste0(scenario, "_", year)]] <- raster
  }
}

# Calculate and create Figure: Projected changes in climatologies (bio1 & 12) for extent -------------------------------------------
# Loop through the variables and years to extract mean values
mean_values_list <- list()
for (variable in variables) {
  for (year in years) {
    for (scenario in scenarios) {
      raster_stack <- predictionLandscapesClim[[paste0(scenario, "_", year)]]
      layer <- raster_stack[[variable]]
      mean_val <- as.numeric(terra::global(layer, fun = "mean", na.rm = TRUE))
      mean_values_list[[paste(variable, scenario, year, sep = "_")]] <- data.frame(
        Variable = variable,
        Scenario = scenario,
        Year = year,
        Mean = mean_val,
        Value_Type = variable_names[which(variables == variable)]
      )
    }
  }
}
mean_values_df <- do.call(rbind, mean_values_list)

# Change name of scenarios to readable names
mean_values_df <- mean_values_df %>%
  mutate(Scenario = recode(Scenario,
                           "ssp126" = "SSP1-RCP2.6",
                           "ssp585" = "SSP5-RCP8.5"))


# Define consistent color palette for the scenarios
scenario_colors <- setNames(
  c("#1f77b4", "#ff7f0e"), scenario_names)

# Split data by variable for plotting
variable_dfs <- lapply(variable_names, function(vn) {
  mean_values_df %>% filter(Value_Type == vn)
})
names(variable_dfs) <- variable_names

# Automate plotting for all variables
plots <- lapply(seq_along(variable_names), function(i) {
  plot_timeChanges(variable_dfs[[i]], variable_names[i], y_labels[i])
})
names(plots) <- variable_names

# Extract the shared legend from one of the plots
shared_legend <- extract_legend(plots[["Temperature (bio1)"]])

# Create labels (a) and (b)
label_grobs <- gridExtra::arrangeGrob(
  grobs = list(
    textGrob("(a)", gp = gpar(fontsize = 16), hjust = 0.5),
    textGrob("(b)", gp = gpar(fontsize = 16), hjust = 0.5)
  ),
  ncol = 2
)

# Arrange the 2 plots in a 1x2 grid
plots_grid <- gridExtra::arrangeGrob(
  label_grobs,
  gridExtra::arrangeGrob(
    grobs = lapply(plots, function(p) p + theme(legend.position = "none")),
    ncol = 2,
    nrow = 1
  ),
  ncol = 1,
  heights = unit(c(0.5, 5), "null"),
  top = textGrob(extent_name, gp = gpar(fontsize = 18))
)

# Combine with the shared legend on the right
combined_plot_time <- gridExtra::grid.arrange(
  plots_grid,
  shared_legend,
  ncol = 1,
  heights = unit(c(5, 1), "null")
)

# Save the combined plot
ggplot2::ggsave(
  filename = file.path(outputPathClimate, paste0("ClimateChange_time_", gsub(" ", "_", extent), ".png")),
  plot = combined_plot_time,
  width = 10, height = 6, dpi = 600
)

# Create Figure: Projected spatial changes in climatologies (bio1 & 12) for extent -------------------------------------------

years <- c(2030, 2050) # years to calculate changes for

# Calculate changes
for (variable in variables) {
  baseline_raster <- trainingLandscapesClim[[variable]]
  for (scenario in scenarios) {
    for (year in years) {
      pred_stack <- predictionLandscapesClim[[paste0(scenario, "_", year)]]
      pred_raster <- pred_stack[[variable]]
      change_raster <- pred_raster - baseline_raster
        assign(paste0("Change_", scenario, "_", variable, "_", year, "_", extent), change_raster)
    }
  }
}
# Calculate min/max values for colour scale
min_values <- list()
max_values <- list()
for (variable in variables) {
  min_value <- Inf
  max_value <- -Inf
  for (year in years) {
    for (scenario in scenarios) {
      raster_stack <- get(paste0("Change_", scenario, "_", variable, "_", year, "_", extent))
      min_value <- min(min_value, min(values(raster_stack), na.rm = TRUE))
      max_value <- max(max_value, max(values(raster_stack), na.rm = TRUE))
    }
  }
  min_value <- round(min_value)
  max_value <- ceiling(max_value)

# Create Color ramp
  if (variable == "bio1") {
  color_ramp <- colorRamp2(c(min_value, 0, max_value), c("blue", "white", "red"))
} else if (variable == "bio12") {
  color_ramp <- colorRamp2(c(min_value, 0, max_value), c("saddlebrown", "yellow", "darkgreen"))
}
  fill_label <- paste("Change in", value_units[which(variables == variable)], " ")

  # Prepare plot grobs for all scenarios and years
  plot_grobs <- list()
  for (i in seq_along(scenarios)) {
    for (j in seq_along(years)) {
      plot_grobs[[paste0("row", i, "_col", j)]] <- gtable::gtable_trim(
        ggplot2::ggplotGrob(
          plot_ClimatespatialChanges(
            raster = get(paste0("Change_", scenarios[i], "_", variable, "_", years[j], "_", extent)),
            extent_geom = extent_sf,
            color_ramp = color_ramp,
            fill_label = fill_label,
            min_value = min_value,
            max_value = max_value
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
  example_plot <- plot_ClimatespatialChanges(
    raster = get(paste0("Change_", scenarios[1], "_", variable, "_", years[1], "_", extent)),
    extent_geom = extent_sf,
    color_ramp = color_ramp,
    fill_label = fill_label,
    min_value = min_value,
    max_value = max_value
  )
  shared_legend <- extract_legend(example_plot)

  # Arrange and save
  final_plot <- gridExtra::grid.arrange(
    grobs = all_grobs,
    layout_matrix = layout_matrix,
    widths = col_widths,
    heights = row_heights,
    bottom = shared_legend,
    top = textGrob(
      paste0(variable_names[which(variables == variable)], " ", fill_label, " for the ", extent),
      gp = gpar(fontsize = 20)
    )
  )

  ggplot2::ggsave(
    filename = file.path(outputPathClimate, paste0("ClimateSpatialChange_", variable, "_", gsub(" ", "_", extent), ".png")),
    plot = final_plot,
    width = 12, height = 8, dpi = 600
  )
}
gc()
