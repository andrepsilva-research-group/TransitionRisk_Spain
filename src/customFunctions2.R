## Name: CustomFunctions2.R ##
## Authors: Jorinde-M. Rieger & Inês Silva ##
## Description: Loads all user-developed functions for LandUseChange.R, ClimateChange.R, inputClimate.R, inputLandUse.R, SDMRun.R ##
## Date: August 5th 2025 ##

#####################
# General Functions #
#####################

# Functions to load and modify rasters-------------------------------------------
# Function to crop and mask rasters
crop_mask_raster <- function(raster, biome_sp) {
  terra::mask(terra::crop(raster, biome_sp), biome_sp)
}

# Function to extract the legend from a ggplot object --------------------------
extract_legend <- function(plot) {
  gtable <- ggplotGrob(plot)
  legend <- gtable$grobs[which(sapply(gtable$grobs, function(x) x$name) == "guide-box")][[1]]
  return(legend)
}

# Function to load and select the biome shapefile ------------------------------
load_biome <- function(biome_name) {
  biome_sf <- sf::st_read("data/Ecoregions2017/Ecoregions2017.shp")
  biome_sf[biome_sf$BIOME_NAME == biome_name, ]
}

# Function to load and select continents ---------------------------------------
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

# Function to intersect biome with continents ----------------------------------
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

# Function to crop the biome boundaries to the continents ----------------------
crop_biome_to_continent <- function(biome, continent_geom) {
  st_intersection(biome, continent_geom)
}

##############################
# Climate realated Functions #
##############################

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

# Function to calculate the average for a given scenario, year, and variable ----
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

# Function to extract mean values from a raster stack --------------------------
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

# Functions to vizualize results -----------------------------------------------
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

# UNDESCRIBED FUNCTION !!!!!!!!!!!! --------------------------------------------
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

#############################
# Land-Use Change Functions #
#############################

# Functions to load rasters ----------------------------------------------------
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

# Functions to modify rasters --------------------------------------------------
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

# Functions to analyze rasters -------------------------------------------------
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

# Function to replace numbers of LULC Types with names (for baseline raster) ---
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

# Functions to vizualize results -----------------------------------------------
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


#######################################
# Format species input data Functions # specific for SDMRun.R 
#######################################

# Function to remove species duplicates by cell ID -----------------------------
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


### Function to retrieve taxa occ from GBIF ------------------------------------

# Downloads GBIF occurrence records for a set of species, applies quality
# filters, and crops occurrences to IUCN range polygons. Retruns: An sf object
# of GBIF occurrences spatially restricted to IUCN ranges.
# Parameters:
#   iucn_shapefile : path to IUCN or BirdLife species range shapefile
#   taxa_filter : optional vector of species names to keep
#   gbif_user : GBIF username
#   gbif_pwd : GBIF password
#   gbif_email : GBIF account email

get_taxa_occurrences <- function(
    iucn_shapefile, # IUCN spatial file
    taxa_filter = NULL, # use in case we want to subset specific species
    gbif_user = "", # GBIF username
    gbif_pwd = "", # GBIF password
    gbif_email = "") # GBIF account email
{
  
  ##########
  # STEP 1 # Read IUCN ranges and species list
  ##########
  
  # progess message
  message("Selecting species names from IUCN ranges file")
  
  IUCN_ranges <- sf::st_read(iucn_shapefile)
  invisible(gc())
  
  IUCN_ranges <- IUCN_ranges %>%
    # fix typo
    mutate(sci_name = if_else(sci_name == "Gyps rueppelli", "Gyps rueppellii", sci_name))
  
  # keep unique species names
  taxa_spp <- unique(IUCN_ranges$sci_name)
  
  # filter taxa if provided
  if (!is.null(taxa_filter)) {
    taxa_spp <- taxa_spp[taxa_spp %in% taxa_filter]} 

  
  ##########
  # STEP 2 # Get GBIF taxon keys
  ##########
  
  gbif_taxon_keys <- name_backbone_checklist(name = taxa_spp) %>%
    filter(matchType == "EXACT") %>%
    pull(usageKey)
  invisible(gc())
  
  ##########
  # STEP 3 # Download occurrences with coordinates only
  ##########
  
  download_key <- occ_download(
    pred_in("taxonKey", gbif_taxon_keys),
    # living species
    pred("OCCURRENCE_STATUS","PRESENT"),
    # type of occurrences
    pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "OBSERVATION", "MACHINE_OBSERVATION", "OCCURRENCE")),
    # with coordinate values
    pred("hasCoordinate", TRUE),
    pred_notnull("decimalLatitude"),
    pred_notnull("decimalLongitude"),
    # from 2015 onwards
    pred_gte("year", 2015),
    # file format to download
    format = "SIMPLE_CSV",
    # GBIF credentials
    user = gbif_user,
    pwd = gbif_pwd,
    email = gbif_email
  )
  
  # wait until finished
  message("Downloading occurrences file from GBIF. Please wait...")
  occ_download_wait(download_key)
  invisible(gc())
  
  # retrieve and import download
  d <- occ_download_get(key = download_key)
  gbif_data <- occ_download_import(d)
  rm(d)
  invisible(gc())
  
  ##########
  # STEP 4 # Filter species with >30 records
  ##########
  
  gbif_data_filtered <- gbif_data %>%
    drop_na(decimalLatitude, decimalLongitude, year) %>% 
    group_by(species) %>% 
    # keep only species with over 30 occ
    dplyr::filter(n() > 30) %>% 
    ungroup()
  
  # are there species with no occurrence records match?
  missing_spp <- setdiff(taxa_spp,  unique(gbif_data_filtered$species))
  
  # Check and report
  if (length(missing_spp) > 0) {
    message("!! WARNING !! The following species do not have enough records:")
    message(paste(missing_spp, collapse = ", "))
  } else {
    message("All taxa_spp species have more than 30 occurrence records! Continuing...")
  }
  
  ##########
  # STEP 5 # Crop GBIF occurrences within IUCN range
  ##########
  
  # progress message
  message("Cropping occurrences by species IUCN range...")
  
  # convert GBIF data to sf points
  gbif_sf <- st_as_sf(
    gbif_data_filtered,
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = st_crs(IUCN_ranges))
  invisible(gc())
  
  # work in a flat earth
  sf_use_s2(FALSE)
  # loop species-by-species to avoid geometry issues
  occurrences_in_range <- bind_rows(lapply(taxa_spp, function(sp) {
    message(sp)
    poly <- IUCN_ranges %>% filter(sci_name == sp)
    pts  <- gbif_sf %>% filter(species == sp)
    if (nrow(poly) == 0 || nrow(pts) == 0) return(NULL)
    return(st_join(pts, poly, join = st_within, left = FALSE))
  }))
  
  # final output
  return(occurrences_in_range)
}

# ####################
# # Taxa Occurrences # 
# ####################
# 
# # This function downloads and clean species occurrence records from GBIF
# # !! NAVIGATION WARNING !! in order to use it, one has to register in the GBIF
# # database website. OR use the api_keys.yml file in the config folder
# 
# get_taxa_occurrences <- function(
#     iucn_shapefile, # (IUCN) spatial file
#     taxa_filter = NULL, # use in case we want to subset specific species
#     gbif_user = "", # GBIF username
#     gbif_pwd = "", # GBIF password
#     gbif_email = "") # GBIF account email
# {
#   
#   ##########
#   # STEP 1 # Read IUCN ranges and species list
#   ##########
#   
#   # progess message
#   message("Selecting species names from IUCN ranges file")
#   
#   IUCN_ranges <- sf::st_read(iucn_shapefile)
#   invisible(gc())
#   
#   # keep unique species names
#   taxa_spp <- unique(IUCN_ranges$sci_name)
#   
#   # filter taxa if provided
#   if (!is.null(taxa_filter)) {
#     taxa_spp <- taxa_spp[taxa_spp %in% taxa_filter]}
#   
#   ##########
#   # STEP 2 # Get GBIF taxon keys
#   ##########
#   
#   gbif_taxon_keys <- name_backbone_checklist(name = taxa_spp) %>%
#     filter(matchType == "EXACT") %>%
#     pull(usageKey)
#   invisible(gc())
#   
#   ##########
#   # STEP 3 # Download occurrences with coordinates only
#   ##########
#   
#   download_key <- occ_download(
#     pred_in("taxonKey", gbif_taxon_keys), 
#     pred("OCCURRENCE_STATUS","PRESENT"),
#     pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "OBSERVATION", "MACHINE_OBSERVATION", "OCCURRENCE")),
#     pred("hasCoordinate", TRUE),
#     pred_notnull("decimalLatitude"),
#     pred_notnull("decimalLongitude"),
#     format = "SIMPLE_CSV",
#     user = gbif_user,
#     pwd = gbif_pwd,
#     email = gbif_email
#   )
#   
#   # wait until finished
#   message("Downloading occurrences file from GBIF. Please wait...")
#   occ_download_wait(download_key)
#   invisible(gc())
#   
#   # retrieve and import download
#   d <- occ_download_get(key = download_key)
#   gbif_data <- occ_download_import(d)
#   rm(d)
#   invisible(gc())
#   
#   ##########
#   # STEP 4 # Filter species with >30 records
#   ##########
#   
#   gbif_data_filtered <- gbif_data %>%
#     drop_na(decimalLatitude, decimalLongitude, year) %>% 
#     dplyr::filter(# keep only species (exclude subspecies for instance)
#       taxonRank == "SPECIES",
#       infraspecificEpithet == "",
#       # only keep occ with year info.
#       year > 0) %>% 
#     group_by(species) %>% 
#     # keep only species with over 30 occ
#     dplyr::filter(n() > 30) %>% 
#     ungroup()
#   
#   ##########
#   # STEP 5 # Crop occurrences within IUCN range
#   ##########
#   
#   # progress message
#   message("Cropping occurrences by species IUCN range...")
#   
#   gbif_sf <- st_as_sf(
#     gbif_data_filtered,
#     coords = c("decimalLongitude", "decimalLatitude"),
#     crs = st_crs(IUCN_ranges))
#   invisible(gc())
#   
#   sf_use_s2(FALSE)
#   occurrences_in_range <- bind_rows(lapply(taxa_spp, function(sp) {
#     message(sp)
#     poly <- IUCN_ranges %>% filter(sci_name == sp)
#     pts  <- gbif_sf %>% filter(species == sp)
#     if (nrow(poly) == 0 || nrow(pts) == 0) return(NULL)
#     st_join(pts, poly, join = st_within, left = FALSE)
#   }))
#   
#   return(occurrences_in_range)
# }

SDMensembleMultiSpecies <- function(targetSpecies, # vector of target species names
                                    speciesData, # target species occurrences file from GBIF
                                    myExpl_full, # training landscape (whole world)
                                    myExplCurrent, # current environment landscape (cropped to biome)
                                    myExplFuture, # future environment landscapes (cropped to biome)
                                    extent, # extent name for files' names (e.g. tropical OR boreal)
                                    output_folder, # folder path to save outputs
                                    maxent_source, # path to maxent.jar file
                                    ncoresToUse # n cores to use in parallelization jobs
) {
  
  # # If changes are required use these args for testing inside the function
  # targetSpecies <- targetSpecies[1]
  # speciesData <- speciesData
  # myExpl_full <- myExpl_full
  # myExplCurrent <- myExplCurrent
  # myExplFuture <- myExplFuture
  # extent <- "GlobalTerrestrial"
  # output_folder <- "./output/28Aug2025"
  # maxent_source <- "C:/Users/maria/Desktop/maxent/maxent/maxent.jar"
  # #"C:/Users/User/OneDrive - Universidade de Lisboa/Ambiente de Trabalho/maxent/maxent/maxent.jar"
  # ncoresToUse <- 6
  
  ##########
  # STEP 1 # Setup & Folder Prep
  ##########
  
  # create output folder
  if(!dir.exists(output_folder)){
    dir.create(output_folder, recursive = TRUE)
  }
  
  if (file.exists(maxent_source)) {
    file.copy(from = maxent_source,
              to = file.path(output_folder, "maxent.jar"),
              overwrite = TRUE)
  } else {
    warning("maxent.jar not found at: ", maxent_source,
            "\nDownload it or place it in this folder before running.")
  }
  
  # set working directory to output folder
  setwd(output_folder) 
  invisible(gc())
  
  ##########
  # STEP 2 # Filter Occurrence Data & Prepare Presence/Pseudo-absence data
  ##########
  
  # print starting message
  message(paste0("Starting for ", targetSpecies))
  
  # Select single species data
  DataSingleSpecies <- speciesData %>%
    dplyr::filter(species == !!targetSpecies)
  invisible(gc())
  
  # Remove NAs and filter out records older than 2015
  # this might reduce the number of presence data to <30 occurences
  DataSingleSpecies <- DataSingleSpecies %>%
    drop_na(decimalLongitude,decimalLatitude, year)
  invisible(gc())
  
  # keep occurrence records after 2015
  DataSingleSpecies <- DataSingleSpecies %>%
    filter(year >= 2015)
  
  # skip sps with less than 30 occ records
  if (nrow(DataSingleSpecies) < 30) {
    message("Skipping ", targetSpecies, " - only ", nrow(DataSingleSpecies), " occurrences >= 2015.")
    next
  }
  
  # assign cell IDs to each occurrence based on myExpl raster
  cellValues <- terra::extract(
    myExpl_full,
    cbind(DataSingleSpecies$decimalLongitude, DataSingleSpecies$decimalLatitude))
  
  cellValues$cell <- terra::cellFromXY(myExpl_full,
                                       cbind(DataSingleSpecies$decimalLongitude, DataSingleSpecies$decimalLatitude))
  
  DataSingleSpecies <- cbind(DataSingleSpecies, cellValues)
  
  # keep one record per cell (to avoid biased occ points)
  DataSingleSpecies_unique <- DataSingleSpecies %>%
    group_by(cell) %>%
    slice_max(year, with_ties = FALSE) %>%  # or slice_head(n = 1) for the first
    ungroup() %>%
    dplyr::filter(complete.cases(.))  # biomod excludes all cells that do not have any data
  
  # subset occurrences
  set.seed(123)
  n_sample <- min(300, nrow(DataSingleSpecies_unique))
  DataSingleSpecies_unique <- DataSingleSpecies_unique %>%
    # use 300 occurrences or all of them if less than 300
    slice_sample(n = n_sample) %>%
    as.data.frame()
  
  # format species occurence data (presence only data)
  myResp <- as.numeric(DataSingleSpecies_unique$species == targetSpecies)
  myRespXY <- DataSingleSpecies_unique[, c("decimalLongitude", "decimalLatitude")]
  
  n.pres <- sum(myResp == 1)
  nb.PA <- c(n.pres, n.pres, n.pres, 10000, 10000, 10000) # number of pseudo-absences per set
  
  # format input data (with initial pseudo-absences set) 
  myBiomodData.PA <- BIOMOD_FormatingData(
    resp.var = myResp,
    expl.var = myExpl_full,
    resp.xy = myRespXY,
    resp.name = targetSpecies,
    PA.nb.rep = 6, # Number of pseudo-absences sets
    PA.nb.absences = nb.PA,  # Adjust as needed. Different for each model
    PA.strategy = 'random', # random PA selection within the given raster
    na.rm = TRUE, # missing values for explanatory variables
    filter.raster = TRUE) # Removes cell duplicates.
  
  # print a message
  message(paste0("Data formatting done for ", targetSpecies))
  
  # save presence points as .csv 
  presence_points <- myBiomodData.PA@coord[myBiomodData.PA@data.species == 1, ]
  presence_df <- as.data.frame(presence_points)
  colnames(presence_df) <- c("Longitude", "Latitude")
  presence_df$type <- "Presence Points"
  presence_df$species <- as.character(targetSpecies)
  presence_df <- na.omit(presence_df)
  
  write.csv(
    presence_df,
    file = file.path(paste0("PresencePoints_", targetSpecies, "_", extent, ".csv")),
    row.names = FALSE)
  
  # Save the presence and pseudo absence points plot
  #png(
  #  filename = file.path(output_folder, paste0("PresencePAPoints_", targetSpecies, "_", extent, ".png")),
  #  width = 2000,
  #  height = 1500,
  #  res = 300
  #)
  #plot(myBiomodData.PA)
  #dev.off()
  
  ##########
  # STEP 3 # Run the models 
  ##########
  
  # selection of models and pseudo-absences set
  models.pa.list <- list(
    RF = c("PA1", "PA2", "PA3"), # Random-forest
    XGBOOST = c("PA1", "PA2", "PA3"), # Extreme Gradient Boosting
    ANN = c("PA1", "PA2", "PA3"), # Artificial Neural Network
    #MAXENT = c("PA4", "PA5", "PA6") # Maximum Entropy Models
    MAXNET = c("PA4", "PA5", "PA6") # replaced MAXENT for MAXNET
  )
  
  # run single models
  myBiomodModelOut <- BIOMOD_Modeling(
    bm.format = myBiomodData.PA,
    modeling.id = paste0("Model_", targetSpecies),
    models = c("RF", "XGBOOST", "ANN", 
               "MAXNET"
               #"MAXENT" # to use MAXENT maxent.jar needs to be inside the working directory
    ), 
    models.pa = models.pa.list,
    CV.strategy = "random",
    CV.nb.rep = 5, # Number of cross-validation runs
    CV.perc = 0.7, # data split, percentage that will be kept for calibration
    OPT.strategy = 'bigboss',
    prevalence = 0.5, # same weight for presences and abs since we have a very inbalanced dataset
    metric.eval = c("TSS", "ROC"), # ADD BOYCE?
    var.import = 3, # could be changed to 1 if we need to save time
    nb.cpu = ncoresToUse, # parallelization
    do.progress = TRUE)
  
  # print progress message
  message(paste0("Single models completed for ", targetSpecies))
  
  # get evaluation scores & variable importance
  eval_scores <- get_evaluations(myBiomodModelOut)
  eval_scores$species <- targetSpecies  # Add species column
  #evaluationScores <- rbind(evaluationScores, eval_scores)  # Combine scores across species
  var_importance <- get_variables_importance(myBiomodModelOut)
  var_importance$species <- targetSpecies  # Add species column
  #variableImportance <- rbind(variableImportance, var_importance)  # Combine importance across species
  
  # save evaluation scores and variable importance to files
  write.csv(eval_scores, file = file.path(paste0("EvalScores_", targetSpecies, "_", extent, ".csv")), row.names = FALSE)
  write.csv(var_importance, file = file.path(paste0("VarImportance_", targetSpecies, "_", extent, "_", ".csv")), row.names = FALSE)
  
  # Save evaluation score boxplots and variables importance
  #png(
  #  filename = file.path(paste0("EvalBoxplot_", targetSpecies, extent, ".png")),
  #  width = 2000,
  #  height = 1500,
  #  res = 300)
  #bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'))
  #dev.off()
  
  # Create a plot for variable importance for all runs
  #varImpData <- bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))$tab
  #filteredData <- varImpData[varImpData$run == "allRun", ]
  #ggplot2::ggplot(filteredData, aes(x = expl.var, y = var.imp, fill = algo)) +
  #  ggplot2::geom_boxplot() +
  #  ggplot2::labs(
  #    title = "Variable Importance for All Runs",
  #    x = "Explanatory Variable",
  #    y = "Variable Importance",
  #    fill = "Model"
  #  ) +
  #  ggplot2::theme_minimal()
  #ggplot2::ggsave(file.path(paste0("VarImpBoxplot_AllRun_", species, extent, ".png")), width = 10, height = 6, dpi = 300)
  
  invisible(gc(rm(myBiomodData.PA)))# clean up to save memory
  invisible(gc(rm(eval_scores, var_importance)))  # clean up to save memory
  
  ##########
  # STEP 4 # Project single models
  ##########
  
  # project single models
  myBiomodProj <- lapply(myExplCurrent, function(env_raster) { # as list to apply to multiple current landscapes 
    BIOMOD_Projection(
      bm.mod = myBiomodModelOut,
      proj.name = 'Current',
      new.env = env_raster,
      models.chosen ='all',
      build.clamping.mask = TRUE,
      nb.cpu = ncoresToUse
    )
  })    
  
  # print progress message
  message(paste0("Single models projections done for ", targetSpecies))
  
  ##########
  # STEP 5 # Do ensemble models
  ##########
  
  # Model ensemble models
  myBiomodEM <- BIOMOD_EnsembleModeling(
    bm.mod = myBiomodModelOut,
    models.chosen ='all',
    em.by ='all',
    em.algo = c('EMmean'),
    metric.select = c('TSS'),
    metric.select.thresh = c(0.6), # threshold will be updated to 0.6
    metric.eval = c('TSS','ROC'),
    nb.cpu = ncoresToUse, #Parallelization
    do.progress = TRUE,
    var.import = 3,
    EMci.alpha = 0.05)
  
  # print progress message
  message(paste0("Ensemble model done for ", targetSpecies))
  
  # get evaluation scores & variable importance for ensemble models
  eval_scoresEM <- get_evaluations(myBiomodEM)
  eval_scoresEM$species <- targetSpecies  # Add species column
  #evaluationScoresEM <- rbind(evaluationScoresEM, eval_scoresEM)  # Combine scores across species
  
  var_importanceEM <- get_variables_importance(myBiomodEM)
  var_importanceEM$species <- targetSpecies  # Add species column
  
  # Save evaluation scores and variable importance to files
  write.csv(eval_scoresEM, file = file.path(paste0("EvalScoresEM_", targetSpecies, "_", extent, ".csv")), row.names = FALSE)
  write.csv(var_importanceEM, file = file.path(paste0("VarImportanceEM_", targetSpecies, "_", extent, ".csv")), row.names = FALSE)
  
  #variableImportanceEM <- rbind(variableImportanceEM, var_importanceEM)  # Combine importance across species
  
  
  # Save evaluation score boxplots and variables importance
  #png(
  #  filename = file.path(paste0("EvalBoxplotEM_", species, extent, ".png")),
  #  width = 2000,
  #  height = 1500,
  #  res = 300)
  #bm_PlotEvalBoxplot(bm.out = myBiomodEM, group.by = c('metric', 'metric'))
  #dev.off()
  
  #png(
  #  filename = file.path(paste0("VarImpBoxplotEM_", species, extent, ".png")),
  #  width = 2000, 
  # height = 1500,
  #  res = 300)
  #bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'merged.by.run'))
  #dev.off()
  
  invisible(gc(rm(eval_scoresEM, var_importanceEM)))  # Clean up to save memory
  
  ##########
  # STEP 6 # Project ensemble models for current conditions
  ##########
  
  # Project ensemble models (from single projections) on current conditions
  myBiomodEMProj <- lapply(myBiomodProj, function(Proj) { # as list to apply to multiple proj 
    BIOMOD_EnsembleForecasting(
      bm.em = myBiomodEM,
      bm.proj = Proj,
      models.chosen ='all',
      metric.binary ='all',
      nb.cpu = ncoresToUse,
      binary.meth = c("TSS"),
      compress = TRUE
    )
  }) 
  
  # print progress message
  message(paste0("Ensemble models' projections for current conditions done for", targetSpecies))
  
  ##########
  # STEP 7 # Project single and ensemble models to future conditions 
  ##########
  
  # Project single models onto future conditions
  myBiomodProjectionFuture <- lapply(myExplFuture, function(future_raster) { # as list to apply to multiple current landscapes 
    BIOMOD_Projection(
      bm.mod = myBiomodModelOut,
      proj.name = "Future",
      new.env = future_raster,
      models.chosen = 'all',
      metric.binary = 'TSS',
      build.clamping.mask = TRUE,
      nb.cpu = ncoresToUse
    )
  }) 
  
  # print progress message
  message(paste0("Single models' projection for future scenarios done for ", targetSpecies))
  
  # Project ensemble-models projections on future variables
  myBiomodEF <- lapply(myBiomodProjectionFuture, function(future_proj) { # as list to apply to multiple future landscapes 
    BIOMOD_EnsembleForecasting(
      bm.em = myBiomodEM,
      bm.proj = future_proj, # not sure if this making the correct correspondence to the layers in  myBiomodProjectionFuture
      models.chosen = 'all',
      #metric.binary = 'all',
      nb.cpu = ncoresToUse,
      binary.meth = c("TSS"),
      compress = "xz"
    )
  })
  
  # print progress message
  message(paste0("Ensemble models' projections for future scenarios done for ", targetSpecies))
  
  ##########
  # STEP 8 # Save ensemble for current and future conditions rasters for each scenario
  ##########
  
  # print progress message
  message(paste0("Saving output rasters for ", targetSpecies))
  
  # Get evaluation results to extract threshold
  evals <- get_evaluations(myBiomodEM)
  th_TSS <- evals$cutoff[evals$metric.eval == "TSS"]
  
  ## Current Conditions Raster ##
  
  EMcurrent <- get_predictions(myBiomodEMProj[[1]], as.data.frame = FALSE)
  
  # save normal suitability (continuous) raster
  EMcurrent_filename <- file.path(
    #output_folder,
    paste0("proj_Current_EM_", gsub(" ", ".", targetSpecies), "_continuous.tif"))
  terra::writeRaster(EMcurrent, EMcurrent_filename, filetype = "GTiff", overwrite = TRUE, gdal = c("COMPRESS=LZW", "PREDICTOR=2", "BIGTIFF=YES"))
  
  # save binary (converted) raster
  bin_rasters <- bm_BinaryTransformation(data = EMcurrent, threshold = th_TSS, do.filtering = FALSE)
  names(bin_rasters) <- paste0("ssp126_2030", names(bin_rasters), "_TSSbin")
  bin_filename <- file.path(
    #output_folder,
    paste0("proj_Current_EM_",gsub(" ", ".", targetSpecies), "_binary.tif"))
  terra::writeRaster(bin_rasters, bin_filename, filetype = "GTiff", overwrite = TRUE, gdal = c("COMPRESS=LZW", "PREDICTOR=2", "BIGTIFF=YES"))
  
  # clean up to save memory
  invisible(gc(rm(EMcurrent, EMcurrent_filename, bin_rasters, bin_filename)))  
  invisible(gc())
  
  ## Future Conditions Rasters ##
  
  # go through each scenario to save it
  lapply(names(myBiomodEF), function(sc){
    
    EFproj <- myBiomodEF[[sc]]
    
    # ---Continuous raster ---
    cont_rasters <- get_predictions(EFproj, as.data.frame = FALSE)
    names(cont_rasters) <- paste0(sc, "_", names(cont_rasters))
    
    # save normal suitability (continuous) raster
    cont_filename <- file.path(
      #output_folder,
      paste0("proj_", sc, "_", gsub(" ", ".", targetSpecies), "_continuous.tif"))
    terra::writeRaster(cont_rasters, cont_filename, filetype = "GTiff", overwrite = TRUE, gdal = c("COMPRESS=LZW", "PREDICTOR=2", "BIGTIFF=YES"))
    
    # --- Binary raster (using TSS threshold & biomod2 function) ---
    bin_rasters <- bm_BinaryTransformation(data = cont_rasters, threshold = th_TSS, do.filtering = FALSE)
    names(bin_rasters) <- paste0("ssp126_2030_", names(bin_rasters), "_TSSbin")
    
    # save a converted (binary) raster
    bin_filename <- file.path(
      #output_folder,
      paste0("proj_", sc, "_", gsub(" ", ".", targetSpecies), "_binary.tif"))
    terra::writeRaster(bin_rasters, bin_filename, filetype = "GTiff", overwrite = TRUE, gdal = c("COMPRESS=LZW", "PREDICTOR=2", "BIGTIFF=YES"))
    #plot(bin_rasters)
    rm(EFproj, cont_rasters, bin_rasters, cont_filename, bin_filename)
    invisible(gc())
  })
  
  ## move up two directories
  setwd("../../")
  
  # Collect metadata for the projection
  #projectionMetadata <- rbind(
  #  projectionMetadata,
  #  data.frame(
  #    species = species,
  #    scenario = scenario,
  #    rasterFile = rasterFilename,
  #    evaluationMetrics = paste(get_evaluations(myBiomodEM), collapse = ";") # Use myBiomodEM here
  #  )
  #)
  #}
  # Return all results as a list
  #return(list(
  #  evaluationScores = evaluationScores,
  #  variableImportance = variableImportance,
  #  evaluationScoresEM = evaluationScoresEM,
  #  variableImportanceEM = variableImportanceEM)
  #)
}


### Fixing species names for prettier plotting ---------------------------------

# This function transforms names WITHOUT spaces into the correct form based on
# the trait dataframe that exists in the data folder of this repo

pretty_species_names <- function(x) {
  library(here)
  # import trait dataframe 
  mammalTraits_2025_03_17 <- read_csv(here("data", "externaldata", "mammalTraits_2025-12-11.csv"),
                                      show_col_types = FALSE)
  
  # pull the species names **WITH SPACES** column 
  with_spaces <- unique(mammalTraits_2025_03_17$sci_name)
  
  # get corresponding names **WITHOUT** spaces
  no_spaces <- gsub(" ", ".", with_spaces)
  
  # match and replace names
  matched <- match(x, no_spaces)
  
  # show warninng if any names are not found
  if (any(is.na(matched))) {
    unmatched <- x[is.na(matched)]
    warning(paste(unmatched, collapse = ", "), "was/were not matched")
  }
  
  return(with_spaces[matched])
}

