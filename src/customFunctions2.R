## Name: CustomFunctions2.R ##
## Authors: Jorinde-M. Rieger ##
## Description: Loads all developed customised functions for LandUseChange.R, ClimateChange.R, inputClimate.R, inputLandUse.R, SDMRun.R ##
## Date: August 5th 2025 ##

#####################################
# General Functions
#####################################
# Functions to load and modify rasters-------------------------------------------
# Function to crop and mask rasters
crop_mask_raster <- function(raster, biome_sp) {
  terra::mask(terra::crop(raster, biome_sp), biome_sp)
}

# Function to extract the legend from a ggplot object
extract_legend <- function(plot) {
  gtable <- ggplotGrob(plot)
  legend <- gtable$grobs[which(sapply(gtable$grobs, function(x) x$name) == "guide-box")][[1]]
  return(legend)
}

# Function to load and select the biome shapefile
load_biome <- function(biome_name) {
  biome_sf <- sf::st_read("data/Ecoregions2017/Ecoregions2017.shp")
  biome_sf[biome_sf$BIOME_NAME == biome_name, ]
}

# Function to load and select continents
load_select_continents <- function(continent_names) {
  continents <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  merged_continents <- list() # Initialize an empty list to store merged continents
  
  # Handle merged "Central & South America"
  if ("Central & South America" %in% continent_names) {
    central_south_america <- continents %>%
      dplyr::filter(subregion %in% c("Central America", "South America", "Caribbean")) %>%
      dplyr::summarise(geometry = st_union(geometry)) %>%
      dplyr::mutate(continent = "Central & South America")
    
    # Ensure column consistency
    missing_columns <- setdiff(names(continents), names(central_south_america))
    for (col in missing_columns) {
      central_south_america[[col]] <- NA  # Add missing columns with NA values
    }
    central_south_america <- central_south_america[names(continents)]  # Reorder columns to match `continents`
    
    # Add to merged continents
    merged_continents[["Central & South America"]] <- central_south_america
    
    # Remove "Central & South America" from the continent names
    continent_names <- setdiff(continent_names, "Central & South America")
  }
  
  # Handle merged "Europe & Asia"
  if ("Europe & Asia" %in% continent_names) {
    europe_asia <- continents %>%
      dplyr::filter(continent %in% c("Europe")) %>%
      dplyr::summarise(geometry = st_union(geometry)) %>%
      dplyr::mutate(continent = "Europe & Asia")
    missing_columns <- setdiff(names(continents), names(europe_asia))
    for (col in missing_columns) {
      europe_asia[[col]] <- NA
    }
    europe_asia <- europe_asia[names(continents)]
    merged_continents[["Europe & Asia"]] <- europe_asia
    continent_names <- setdiff(continent_names, "Europe & Asia")
  }
  
  # Filter the remaining continents
  remaining_continents <- continents %>%
    dplyr::filter(continent %in% continent_names)
  
  # Combine merged and remaining continents
  all_continents <- do.call(rbind, c(merged_continents, list(remaining_continents)))
  
  # Return the final dataset
  return(all_continents)
}

# Function to intersect biome with continents
#intersect_biome_with_continents <- function(biome_sf, continent_geoms) {
intersect_extent_continents <- function(extent_sf, continent_geoms) {
  # Validate and fix geometries
  extent_sf <- st_make_valid(extent_sf)
  continent_geoms <- lapply(continent_geoms, st_make_valid)
  
  # Perform intersection and handle empty geometries
  biome_continents <- setNames(lapply(continent_geoms, function(continent_geom) {
    result <- st_intersection(extent_sf, continent_geom)
    if (is.null(result) || nrow(result) == 0) {
      return(NULL)  # Return NULL if no intersection
    }
    return(result)
  }), names(continent_geoms))
  
  return(biome_continents)
}

# Function to crop the biome boundaries to the continents
crop_biome_to_continent <- function(biome, continent_geom) {
  st_intersection(biome, continent_geom)
}

#####################################
# Climate realated Functions #
#####################################
# Functions to load rasters -------------------------------------------
load_baseline_clim <- function(variable, baseline_year){
  file_path <- file.path(basePathClim, paste0("CHELSA_", variable, "_", baseline_year, "_V.2.1.tif"))
  rast(file_path)}

load_scenario_clim <- function(scenario, variable, year) {
  file_path <- file.path(basePathClim, scenario, paste0("CHELSA_", variable, "_", year, "_gfdl-esm4_", scenario, "_V.2.1.tif"))
  rast(file_path)}

load_average_scenario_clim <- function(outputPath, scenario, year, variable){
  file_path <- file.path(outputPath, paste0("AverageCHELSA",variable, "_", scenario, "_", yearOrigin, ".tif"))
  terra::rast(file_path)}

# Function to calculate the average for a given scenario, year, and variable
# to save computation time this could be saved as netCDF (terra::writeCDF)
average_climate_models <- function(outputPath, scenario, yearOrigin, variable) {
  # Initialize a list to store rasters for all models
  model_rasters <- list()
  
  # Loop through the models
  for (model in models) {
    # Construct the file path for the raster
    raster_file <- file.path(basePathClim, scenario, paste0("CHELSA_", variable, "_", yearOrigin, "_", model, "_", scenario, "_V.2.1.tif"))
    
    # Check if the file exists
    if (!file.exists(raster_file)) {
      warning(paste("File not found:", raster_file))
      next
    }
    
    # Load the raster
    model_rasters[[model]] <- terra::rast(raster_file)
  }
  
  # Combine the rasters into a SpatRaster stack
  model_stack <- terra::rast(model_rasters)
  
  # Calculate the average using terra::app
  averaged_raster <- terra::app(model_stack, fun = mean, na.rm = TRUE)
  
  # Rename the layer
  names(averaged_raster) <- variable
  
  # Save the averaged raster to disk
  output_file <- file.path(outputPath, paste0("AverageCHELSA",variable, "_", scenario, "_", yearOrigin, ".tif"))
  terra::writeRaster(averaged_raster, output_file, overwrite = TRUE)
  
  return(averaged_raster)
}

# Function to extract mean values from a raster stack
extract_mean_values <- function(raster_stack, years, value_type) {
  mean_values <- sapply(1:nlyr(raster_stack), function(i) {
    mean(values(raster_stack[[i]]), na.rm = TRUE)
  })
  data.frame(
    Year = years,
    Scenario = names(raster_stack),
    Mean_Value = mean_values,
    Value_Type = value_type
  )
}

# Functions to vizualize results -------------------------------------------
# Function to create plots
plot_timeChanges <- function(mean_values_df, value_type, y_label) {
  ggplot(mean_values_df, aes(x = Year, y = Mean, color = Scenario, group = Scenario)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = scenario_colors,
                       guide = guide_legend(direction = "horizontal")
    ) +
    labs(
      x = "Year",
      y = y_label
    ) +
    theme_minimal()+
    theme(
      axis.title.x = element_text(size = 14, margin = ggplot2::margin(t = 10)),  # Increase gap for x-axis title
      axis.title.y = element_text(size = 14, margin = ggplot2::margin(r = 10)),   # Increase gap for y-axis title
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "bottom" 
    )
}


plot_ClimatespatialChanges <- function(raster, extent_geom, color_ramp, fill_label, min_value, max_value) {
  # Convert raster to data frame
  raster_df <- as.data.frame(raster, xy = TRUE)
  colnames(raster_df)[3] <- "value"  # Percentage change (%)
  
  # Load country boundaries
  countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  # Exclude sovereign states and keep only relevant territories
  countries <- countries %>%
    dplyr::filter(sovereignt != "France")  # Exclude France as a sovereign state
  
  # Ensure CRS consistency
  extent_geom <- st_transform(extent_geom, crs = st_crs(countries))
  countries <- st_transform(countries, crs = st_crs(extent_geom))
  
  # Validate geometries
  countries <- st_make_valid(countries)
  extent_geom <- st_make_valid(extent_geom)
  
  # Identify countries overlapping with the biome
  overlapping_indices <- st_intersects(countries, extent_geom, sparse = TRUE)
  overlapping_countries <- countries[lengths(overlapping_indices) > 0, ]
  
  # Create the plot
  plot <- ggplot(raster_df) +
    # Add extent boundary
    geom_sf(data = extent_geom, aes(color = "Area"), show.legend = FALSE, fill = "lightgrey", size = 0.2) +
    # Add country boundaries
    geom_sf(data = overlapping_countries, aes(color = "Country Boundaries"), show.legend = FALSE, fill = NA, size = 0.2) +
    # Add raster data
    geom_tile(data = raster_df, aes(x = x, y = y, fill = value)) +
    # Define the color scale for the raster
    scale_fill_gradientn(name = fill_label, 
                         colors = color_ramp(seq(min_value, max_value, length.out = 101)),
                         limits = c(min_value, max_value), 
                         breaks = c(min_value, 0, max_value), 
                         labels = c(min_value, 0, max_value), 
                         na.value = "grey",
                         guide = guide_colorbar(
                           direction = "horizontal",
                           title.position = "left",      # Title to the left of the colorbar
                           label.position = "bottom",    # Labels below the colorbar
                           title.vjust = 0.5,            # Center the title vertically
                           barwidth = unit(6, "cm"),     # Adjust as needed
                           barheight = unit(0.5, "cm")   # Adjust as needed
                         )
    ) +
    # Define the color scale for the extent area and country boundaries
    scale_color_manual(
      name = NULL,
      values = c("Country Boundaries" = "black", "Area" = "lightgrey"),
      breaks = c("Country Boundaries", "Area"),  # Ensure these match the aes(color = ...) values
      labels = c("Country Boundaries", "Area"),
      guide = guide_legend(direction = "horizontal")
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"  # <--- Ensures legends are in a horizontal box
    )+
    # Add labels and theme
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      plot.title = element_blank(),
      legend.title = element_text(size = 22, margin = ggplot2::margin(b = 10)),
      legend.text = element_text(size = 18),
      legend.key.height = unit(1, "cm"),  # Increase the height of the color ramp
      legend.spacing = unit(1, "cm")
    ) +
    coord_sf()  # Use coord_sf() for spatial data
  
  return(plot)
}

#####################################
# Land-Use Change Functions #
#####################################
# Functions to load rasters-------------------------------------------
# Baseline raster
load_baseline_landUse <- function(baseline_year){
  rast("data/stitched_lulc_esa_scenarios/lulc_esa_2015.tif")
}
# Scenario Rasters
load_scenario_landUse <- function(scenario, year) {
  # Construct the file path
  file_path <- file.path(basePathLandUse, scenario, paste0("lulc_esa_gtap1_", scenario, "_", year, "_no_policy.tif"))
  rast(file_path)
}
# Mapped baseline raster
load_mapped_baseline_landUse <- function(outputPath, baseline_year){
  mapped_baseline_path <- file.path(outputPath, paste0("MappedLandUse_base_", baseline_year, "_", gsub(" ", "_", extent), ".tif"))
  rast(mapped_baseline_path)
}
# Mapped raster stacks
load_mapped_landUse <- function(outputPath, year, extent) {
  mapped_file_path <- file.path(outputPath, paste0("MappedLandUse_scenarios_", year, "_", gsub(" ", "_", extent), ".tif"))
  if (file.exists(mapped_file_path)) {
    mapped_raster_stack <- rast(mapped_file_path)
    assign(paste0("MappedLandUse_scenarios_", year, "_", gsub(" ", "_", extent)), mapped_raster_stack, envir = .GlobalEnv)
    return(mapped_raster_stack)
  } else {
    stop(paste("Mapped raster file for year", year, "does not exist."))
  }
}

# Functions to modify rasters-------------------------------------------
# Function to stack rasters
stack_rasters <- function(year, scenarios, extent, outputPathLandscapes) {
  scenarios_list <- list()
  
  for (scenario in scenarios) {
    # Dynamically construct the raster file path
    raster_file <- file.path(outputPathLandscapes, paste0("LandUse_", scenario, "_", year, "_", gsub(" ", "_", extent), ".tif"))
    
    # rasterize file
    scenarios_list[[paste0(scenario, "_", year)]] <- terra::rast(raster_file)
  }
  
  # Create a raster stack from the list of scenarios
  scenarios_stack <- terra::rast(scenarios_list)
  names(scenarios_stack) <- names(scenarios_list)
  
  # Save the raster stack to disk
  stack_output_file <- file.path(outputPathLandscapes, paste0("LandUse_scenarioStack_", year, "_", gsub(" ", "_", extent), ".tif"))
  writeRaster(scenarios_stack, stack_output_file, overwrite = TRUE)
  
  # Return the raster stack
  return(scenarios_stack)
}

# Functions to analyze rasters-------------------------------------------
# Function to calculate the percentages for each land-use type
calculate_landUse_percentages <- function(raster_stack, landUse_types, landUse_names, time) {
  # Create binary maps for each land-use type
  landUse_layers <- lapply(landUse_types, function(cat) {
    terra::app(raster_stack, fun = function(x) {
      return(ifelse(x == cat, 1, 0))
    })
  })
  
  # Calculate the percentages for each land-use type without NAs
  landUse_percentages <- sapply(landUse_layers, function(layer) {
    sum(values(layer), na.rm = TRUE) / sum(!is.na(values(layer))) * 100
  })
  
  # Create a data frame with the results
  percentage_df <- data.frame(
    time = time,
    landUse = landUse_names,
    value = landUse_percentages
  )
  return(percentage_df)
}

# Function to process the mapped scenarios and apply the function to calculate percentages
process_and_map_scenarios <- function(year, use_continent = FALSE, continent = NULL) {
  
  # Load the raster stack for the years
  if (use_continent && !is.null(continent)) {
    mapped_raster_stack <- continent_scenarios[[continent]][[as.character(year)]]
  } else {
    mapped_raster_stack <- mapped_scenarios[[as.character(year)]]
  }
  
  # Calculate land use percentages
  scenarios_percentages_list <- list()
  for (i in 1:nlyr(mapped_raster_stack)) {
    scenario_name <- names(mapped_raster_stack)[i]
    raster_layer <- mapped_raster_stack[[i]]
    percentages_df <- calculate_landUse_percentages(
      raster_stack = raster_layer,
      landUse_types = LULC_Types,
      landUse_names = LULC_Types_names,
      time = year
    )
    percentages_df$Scenario <- scenario_name
    scenarios_percentages_list[[scenario_name]] <- percentages_df
  }
  
  # Combine the results into a single data frame
  scenarios_percentages_df <- do.call(rbind, scenarios_percentages_list)
  return(scenarios_percentages_df)
}

# Function to create binary maps and classes of land-use types
calculateRasterClass <- function(OriginalRaster, extent, target_resolution) {
  # Crop and mask the raster to the biome's boundary
  raster <- mask(crop(OriginalRaster, extent), extent)
  
  # Define the unique land-use classes and remove NAs
  land_use_classes <- terra::freq(raster)[,2]
  land_use_classes <- na.omit(land_use_classes)
  
  # Function to create binary raster for each land-use class
  create_binary_raster <- function(raster, land_use_class) {
    binary_raster <- app(raster, fun = function(x) {
      ifelse(x == land_use_class, 1, 0)
    })
    return(binary_raster)
  }
  
  # Create a list to store binary rasters
  binary_rasters <- list()
  
  # Loop through each land-use class and create binary rasters
  for (class in land_use_classes) {
    binary_rasters[[as.character(class)]] <- create_binary_raster(raster, class)
  }
  
  # Calculate the aggregation factor based on the target resolution
  input_resolution <- res(raster)[1]  # Assuming square cells, take the resolution of the first dimension
  aggregation_factor <- round(target_resolution / input_resolution)
  
  # Aggregate each binary raster by aggregation factor
  aggregated_rasters <- list()
  for (class in names(binary_rasters)) {
    aggregated_raster <- aggregate(binary_rasters[[class]], fact = aggregation_factor, fun = function(x) sum(x > 0, na.rm = TRUE)/length(x))
    masked_raster <- terra::mask(terra::crop(aggregated_raster, extent), extent)
    aggregated_rasters[[class]] <- masked_raster
  }
  
  # Convert the list of rasters to a SpatRaster stack
  aggregated_rasters_stack <- terra::rast(aggregated_rasters)
  return(aggregated_rasters_stack)
}

# Function to replace numbers of LULC Types with names (for baseline raster)
replace_numbers_with_names <- function(raster_list, types, names) {
  # Replace the names of the raster layers with the corresponding land-use names
  names(raster_list) <- names[match(names(raster_list), types)]
  return(raster_list)
}

# Function to replace numbers of LULC Types with names in a nested list (for scenario rasters)
replace_numbers_with_names_nested <- function(nested_list, types, names) {
  for (scenario in names(nested_list)) {
    for (year in names(nested_list[[scenario]])) {
      # Replace the names of the land-use types within each year
      names(nested_list[[scenario]][[year]]) <- names[match(names(nested_list[[scenario]][[year]]), types)]
    }
  }
  return(nested_list)
}

# Function to calculate percentage changes for each land-use class for scenarios
calculate_percentage_changes <- function(base_year_raster, target_year_rasters_list, scenarios, years) {
  percentage_change_rasters_list <- list()
  
  for (year in years) {
    percentage_change_rasters_list[[as.character(year)]] <- list()
    
    for (scenario in scenarios) {
      key <- paste0(scenario, "_", year)
      target_raster <- target_year_rasters_list[[key]]
      percentage_change_classes <- list()
      
      for (class in names(base_year_raster)) {
        base_raster <- base_year_raster[[class]]
        target_raster_class <- target_raster[[class]]
        
        # Define a mask: include only where the class is present in either year
        # Define presence mask
        presence_mask <- (base_raster > 0 | target_raster_class > 0)
        
        # Mask both rasters
        base_masked <- terra::mask(base_raster, presence_mask, maskvalues = FALSE)
        target_masked <- terra::mask(target_raster_class, presence_mask, maskvalues = FALSE)
        
        # Calculate percentage change only in valid areas
        percentage_change <- terra::ifel(
          base_masked == 0,
          target_masked * 100,
          (target_masked - base_masked) / base_masked * 100
        )
        
        percentage_change_classes[[class]] <- percentage_change
      }
      percentage_change_rasters_list[[as.character(year)]][[scenario]] <- percentage_change_classes
    }
  }
  return(percentage_change_rasters_list)
}

# Functions to vizualize results -------------------------------------------
# Function to create individual plots for each scenario, class, and year
plot_landUse_spatialChanges <- function(raster, biome_geom, color_ramp, fill_label, min_value, max_value, coord_limits = NULL) {
  # Convert raster to data frame
  raster_df <- as.data.frame(raster, xy = TRUE)
  colnames(raster_df)[3] <- "value"  # Percentage change (%)
  
  # Load country boundaries
  countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  # Exclude sovereign states and keep only relevant territories
  countries <- countries %>%
    dplyr::filter(sovereignt != "France")  # Exclude France as a sovereign state
  
  # Ensure CRS consistency
  biome_geom <- st_transform(biome_geom, crs = st_crs(countries))
  countries <- st_transform(countries, crs = st_crs(biome_geom))
  
  # Validate geometries
  countries <- st_make_valid(countries)
  biome_geom <- st_make_valid(biome_geom)
  
  # Identify countries overlapping with the biome
  overlapping_indices <- st_intersects(countries, biome_geom, sparse = TRUE)
  overlapping_countries <- countries[lengths(overlapping_indices) > 0, ]
  
  # Create the plot
  plot <- ggplot(raster_df) +
    # Add biome boundary
    geom_sf(data = biome_geom, aes(color = "Area"), fill = "lightgrey", size = 0.2) +
    # Add country boundaries
    geom_sf(data = overlapping_countries, aes(color = "Country Boundaries"), fill = NA, size = 0.2) +
    # Add raster data
    geom_tile(data = raster_df, aes(x = x, y = y, fill = value)) +
    # Define the color scale for the raster
    scale_fill_gradientn(
      name = fill_label,
      colors = color_ramp(seq(-100, 100, length.out = 101)),
      limits = c(min_value, max_value),
      na.value = "grey",
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "left",      # Title to the left of the colorbar
        title.theme = element_text(size = 18,  margin = ggplot2::margin(b = 10)),  # increase size here
        label.position = "bottom",    # Labels below the colorbar
        title.vjust = 0.5,            # Center the title vertically
        barwidth = unit(6, "cm"),     # Adjust as needed
        barheight = unit(1, "cm")   # Adjust as needed
      )
    ) +
    # Define the color scale for the biome and country boundaries
    scale_color_manual(
      name = NULL,
      values = c("Country Boundaries" = "black", "Area" = "lightgrey"),
      breaks = c("Country Boundaries", "Area"),  # Ensure these match the aes(color = ...) values
      labels = c("Country Boundaries", "Area"),
      guide = guide_legend(direction = "horizontal")
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"  # <--- Ensures legends are in a horizontal box
    )+
    # Add labels and theme
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      plot.title = element_blank(),
      legend.title = element_text(size = 22, margin = ggplot2::margin(b = 10)),
      legend.text = element_text(size = 18),
      legend.key.height = unit(1, "cm"),  # Increase the height of the color ramp
      legend.spacing = unit(0.5, "cm")
    ) +
    if (!is.null(coord_limits)) {
      coord_sf(xlim = coord_limits$xlim, ylim = coord_limits$ylim, expand = FALSE)
    } else {
      coord_sf()
    }
  return(plot)
}


#####################################
#  Format species input data Functions for SDMRun.R #
#####################################

# Function to remove species duplicates by cell ID
removeSpeciesDuplicatesbyCellID <- function (dataframe) {
  SpeciesDataOcc <- dataframe %>%
    drop_na(cell) %>%
    group_by(species, cell) %>%
    slice_max(year, with_ties = FALSE) %>%  # Keep most recent record per species-cell
    ungroup() %>%
    distinct(species, cell, year, .keep_all = TRUE) %>%  # Ensure unique species-cell-year
    arrange(species, cell, year)
  # Format numeric columns to have at least 5 decimal digits
  numeric_cols <- sapply(SpeciesDataOcc, is.numeric)
  SpeciesDataOcc[numeric_cols] <- lapply(SpeciesDataOcc[numeric_cols], function(x) round(x, 5))
  return(SpeciesDataOcc)
}