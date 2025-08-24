## Name: binaryContinousSDMmaps.R
## Author: Inês Silva
## Description: create presentable plots from continuous and binary SDM rasters
## Date: 24 August 2025

# Settings & libraries ---------------------------------------------------------
source("src/libraries.R") # libraries
source("src/customFunctions2.R") # functions

##########
# STEP 1 # Basemap for Portugal and Spain
##########

extent <- "Spain and Portugal"
ptsp_map <- rnaturalearth::ne_countries(scale = "medium", country = c("Spain", "Portugal"), returnclass = "sf")
invisible(gc())

##########
# STEP 2 # Write function to plot both continuous and binary plots and save
##########

# Load Portugal + Spain as basemap (do this once outside the function)
ptsp_map <- ne_countries(scale = "medium", country = c("Portugal", "Spain"), returnclass = "sf")

species_name = targetSpecies[[1]]
sc = "ssp126"
plot_scenario <- function(sc, species_name, output_folder){
  
  # --- File paths ---
  cont_file <- file.path(output_folder, paste0("EnsembleForecast_", gsub(" ", ".", species_name), "_", sc, "_continuous.tif"))
  bin_file  <- file.path(output_folder, paste0("EnsembleForecast_", gsub(" ", ".", species_name), "_", sc, "_binary.tif"))
  
  if(!file.exists(cont_file) | !file.exists(bin_file)){
    warning(paste("Files not found for", species_name, "scenario", sc))
    return(NULL)
  }
  
  # --- Load rasters ---
  cont_rast <- rast(cont_file)
  bin_rast  <- rast(bin_file)
  
  # --- Convert to data frames ---
  cont_df <- as.data.frame(cont_rast, xy = TRUE, na.rm = TRUE)
  names(cont_df)[3] <- "suitability"
  
  bin_df <- as.data.frame(bin_rast, xy = TRUE, na.rm = TRUE)
  names(bin_df)[3] <- "presence"
  
  # --- Continuous plot ---
  p_cont <- ggplot() +
    geom_raster(data = cont_df, aes(x = x, y = y, fill = suitability)) +
    geom_sf(data = ptsp_map, fill = NA, color = "black", size = 0.5) +
    scale_fill_viridis_c(name = "Suitability", limits = c(0, 1000)) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal()
  
  # --- Binary plot ---
  p_bin <- ggplot() +
    geom_sf(data = ptsp_map, fill = "grey95", color = "black") +
    geom_raster(data = bin_df, aes(x = x, y = y, fill = factor(presence)), alpha = 0.6) +
    scale_fill_manual(
      name = NULL,
      values = c("0" = "transparent", "1" = "green4"),
      labels = c(" ","Presence")
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal()
  
  # --- Combine ---
  combined_plot <- p_cont / p_bin + patchwork::plot_annotation(
    title = bquote(italic(.(species_name)) ~ " Distribution and Presence Maps"),
    subtitle = sc
  )
  return(combined_plot)
  
  # ggsave
  ggsave(plot = combined_plot,
         file = ".tiff",
         bg = 'white', width = 250, height = 230, units = "mm", dpi = 1200, compression ="lzw")
 
}

##########
# STEP 3 # Apply function
##########

species_list <- c("Streptopelia turtur", "Otis tarda", "Chlamydotis undulata",
                     "Neophron percnopterus", "Gyps rueppellii", "Aquila adalberti",
                     "Lanius meridionalis", "Chersophilus duponti", "Fringilla polatzeki",
                     "Pyrrhula murina")

scenarios <- c("ssp126_2030_ptsp", "ssp126_2050_ptsp", "ssp126_2100_ptsp", "ssp585_2030_ptsp", "ssp585_2050_ptsp", "ssp585_2100_ptsp")
output_folder <- ""
  
all_plots <- list()

for(sp in species_list){
  for(sc in scenarios){
    message("Processing ", sp, " - ", sc)
    p <- plot_scenario(sc, sp, output_folder, scale_to = 1)
    if(!is.null(p)){
      all_plots[[paste(sp, sc, sep = "_")]] <- p
    }
  }
}

# Example: show one plot
all_plots[["Streptopelia turtur_ssp126_2030_ptsp"]]
