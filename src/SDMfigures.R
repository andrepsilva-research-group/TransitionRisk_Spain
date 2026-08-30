## Name: SDMfigures.R ##
## Author: Jorinde-M. Rieger & Inês Silva ##
## Date: 30th Aug 2025 ##
## Description: Create all suplementary materials (figures and tables) from
## multiple SDMs ##

source("src/customFunctions2.R")

# paths for outputs
mammalsOutputs <- "E:/BankOfSpain_june2026/SDMs_mammals" #REPLACE WITH YOUR OWN PATH
birdsOutputs <- "E:/BankOfSpain_june2026/SDMs_birds" # REPLACE WITH YOUR OWN PATH
pathSDMoutputs <- c(mammalsOutputs, birdsOutputs)

# folder to save figure and tables on SDM outputs
SDMsFigures_dir <- "E:/BankOfSpain_june2026/SDM_outputfigures" # adapt if needed
if (!dir.exists(SDMsFigures_dir)) {
  dir.create(SDMsFigures_dir, recursive = TRUE)
}
# create subfolder for presence plots
presencePlotBase <- file.path(SDMsFigures_dir, "presencePlots")
if (!dir.exists(presencePlotBase)) {
  dir.create(presencePlotBase, recursive = TRUE)
}


# scenarios names (used for files)
scenarios <- c("Current", "ssp126_2030", "ssp126_2050",
               "ssp585_2030", "ssp585_2050")

# scenarios prettier labels
scenario_labels <- c(
  Current = "Current",
  ssp126_2030 = "SSP1-2.6 (2030)",
  ssp126_2050 = "SSP1-2.6 (2050)",
  ssp585_2030 = "SSP5-8.5 (2030)",
  ssp585_2050 = "SSP5-8.5 (2050)")


iberian_threatened_mammals <- read_csv("data/iberian_threatened_mammals_taxonomy.csv") %>% 
  dplyr::pull(sci_name)

iberian_threatened_birds <- read_csv("data/iberian_threatened_birds_taxonomy.csv") %>% 
  dplyr::pull(sci_name)

targetSpecies <- c(iberian_threatened_birds, iberian_threatened_mammals)

# get basemaps objects
extent <- "IberianPeninsula"
world <- ne_countries(scale = "medium", returnclass = "sf")
iberia <- world[world$name %in% c("Portugal", "Spain"), ]
# vector for masking
iberia_vect <- vect(iberia)   
# iberian peninsula bbox
iberia_bbox <- st_bbox(iberia) 
gc()

############
# OUTPUT 1 # Presence maps
############

# Load GBIF occurrence data (fixed file assignments)
gbif_mammals_all <- read_csv("data/GBIF_occurrences_mammals_2026-06-07.csv",
                             show_col_types = FALSE)

gbif_birds_all <- read_csv("data/GBIF_occurrences_birds_2026-06-08.csv",
                           show_col_types = FALSE)

# Put datasets in a named list
gbif_data <- list(mammals = gbif_mammals_all, birds = gbif_birds_all)

# Standardize column names for each dataset
gbif_data <- lapply(gbif_data, function(df) {
  df <- df[, 2:5]
  
  # careful not to switch axes
  colnames(df)[1:3] <- c("species", "Latitude", "Longitude")
  df$species <- gsub(" ", ".", df$species)
  
  return(df)
})

presencePlots <- list()

for (path in pathSDMoutputs) {
  
  # Detect taxa from path
  taxa <- if (grepl("mammals", path, ignore.case = TRUE)) {
    "mammals"} else {"birds"}
  
  # Select matching GBIF dataset
  gbif_all <- gbif_data[[taxa]]
  
  # Extract region from path (assuming folder structure contains region)
  region <- extent
  
  # Create output folder
  taxa_output_dir <- file.path(presencePlotBase, taxa)
  if (!dir.exists(taxa_output_dir)) {
    dir.create(taxa_output_dir, recursive = TRUE)
  }
  
  for (sps in targetSpecies) {
    
    message(paste0("Building Presence Plot for ", sps, " (", taxa, ")"))
    
    # Presence points used in SDMs
    csv_file <- file.path(path, paste0("PresencePoints_", gsub(" ", ".", sps), "_", region, ".csv"))
    
    if (!file.exists(csv_file)) next
    
    presence_df <- read_csv(csv_file, show_col_types = FALSE)
    colnames(presence_df)[1:2] <- c("Longitude", "Latitude")
    presence_df$Type <- "SDM Presences"
    
    # Total GBIF occurrence records
    gbif_species <- gbif_all %>%
      filter(species == gsub(" ", ".", sps)) %>%
      mutate(Type = "GBIF Occurrences")
    
    # Color-blind palette
    color_map <- c("GBIF" = "#0072B2", "SDM"  = "#E69F00")
    
    # Plot
    p <- ggplot() +
      geom_sf(data = world, fill = "grey95", color = "black") +
      geom_point(data = gbif_species,
                 aes(x = Longitude, y = Latitude, color = "GBIF"), size = 0.4, alpha = 0.6) +
      geom_point(data = presence_df,
                 aes(x = Longitude, y = Latitude, color = "SDM"), size = 0.7, alpha = 0.8) +
      scale_color_manual(values = color_map, breaks = c("GBIF", "SDM"),
        labels = c("GBIF occurrences", "SDM presences"),
        drop = FALSE) +
      theme_minimal(base_size = 12) +
      theme(legend.title = element_blank(), 
            legend.position = "bottom",
            legend.background = element_rect(fill = "white", color = "grey80"),
            panel.grid = element_line(color = "grey90"))
    
    invisible(gc())
    
    # Save plot
    ggsave(plot = p, file = file.path(taxa_output_dir,
                                      paste0("PresencesPlot_", gsub(" ", ".", sps), ".png")),
      bg = "white", width = 300, height = 150, units = "mm", dpi = 1200)
  }
}

############
# OUTPUT 2 # Variable Importance for Ensemble Models Table
############

all_var_importance_em <- data.frame() # create data frame for all variable importance of EM

for (path in pathSMDoutputs) {
  
  # infer group from folder name
  taxon_group <- ifelse(grepl("mammal", path, ignore.case = TRUE),
                        "mammals",
                        ifelse(grepl("bird", path, ignore.case = TRUE),
                               "birds",
                               NA))
  
  for (species in targetSpecies) {
    file_path <- file.path(path,
                           paste0("VarImportanceEM_", gsub(" ", ".", species), "_", extent, ".csv"))
    if (file.exists(file_path)) {
      var_imp <- read.csv(file_path)
      var_imp$species <- species  # Add species column if not present
      var_imp$taxon_group <- taxon_group   # <-- ADD THIS COLUMN
      all_var_importance_em <- rbind(all_var_importance_em, var_imp)
    }
  }
}

# Filter for algo == "EMmean"
var_importance_em_filtered <- all_var_importance_em %>%
  filter(algo == "EMmean")

# Summarize mean and SD for each species and variable
importance_summary <- var_importance_em_filtered %>%
  group_by(taxon_group, species, expl.var) %>%
  dplyr::summarize(
    mean_importance = round(mean(var.imp, na.rm = TRUE), 4), # Calculate mean
    sd_importance = round(sd(var.imp, na.rm = TRUE), 4), # Calculate standard deviation
    .groups = "drop") %>%
  pivot_longer(cols = c(mean_importance, sd_importance), names_to = "metrics", values_to = "importance") %>%
  mutate(metrics = ifelse(metrics == "mean_importance", "Mean", "SD")) %>%
  pivot_wider(names_from = expl.var, values_from = importance) %>%
  dplyr::arrange(species, metrics) %>% 
  mutate(taxon_group = case_when(taxon_group == "birds" ~ "Bird", taxon_group == "mammals" ~ "Mammal"))

# write table to .csv and .xslx (for easy copy paste later)
write.csv(importance_summary,
          file = file.path(SDMsFigures_dir,
                           paste0("VariableImportanceSummaryTable", Sys.Date(), ".csv")), row.names = FALSE)
writexl::write_xlsx(importance_summary,
                    path = file.path(SDMsFigures_dir,
                                     paste0("VariableImportanceSummaryTable", Sys.Date(), ".xlsx")))

############
# OUTPUT 3 # Evaluation plots for Ensemble Models
############

all_EvalScores <- data.frame() 

for (path in pathSDMoutputs) {
  
  taxon_group <- ifelse(grepl("mammal", path, ignore.case = TRUE), "mammals",
                        ifelse(grepl("bird", path, ignore.case = TRUE), "birds", NA))
  
  for (species in targetSpecies) {
    
    file_path <- file.path(path,
                           paste0("EvalScores_", gsub(" ", ".", species), "_", extent, ".csv"))
    
    if (file.exists(file_path)) {
      eval_scores <- read.csv(file_path)
      eval_scores$species <- species
      eval_scores$taxon_group <- taxon_group
      all_EvalScores <- rbind(all_EvalScores, eval_scores)
    }
  }
}

plots_eval <- list()

for (g in unique(all_EvalScores$taxon_group)) {
  
  df_g <- all_EvalScores %>%
    dplyr::filter(taxon_group == g)
  
  p <- ggplot(df_g, aes(x = species, y = calibration, fill = metric.eval)) +
    geom_boxplot() +
    labs(#title = paste("Evaluation Metrics (EM) -", g),
         x = "Species",
         y = "Calibration Score",
         fill = "Metric") +
    scale_fill_viridis_d(name = "Metric") +
    theme_minimal() +
    theme(text=element_text(size=18),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"))
  
  plots_eval[[g]] <- p
}

ggsave(plot = plots_eval$mammals,
       file = file.path(SDMsFigures_dir, "EvaluationMetrics_mammals.png"),
       bg = 'white', width = 300, height = 200, units = "mm", dpi = 1200#, compression ="lzw"
       )

ggsave(plot = plots_eval$birds,
       file = file.path(SDMsFigures_dir, "EvaluationMetrics_birds.png"),
       bg = 'white', width = 300, height = 200, units = "mm", dpi = 1200#, compression ="lzw"
       )

############
# OUTPUT 4 # Continuous landscapes 
############

# Storage for plots
continuous_list <- list()

# first for loop to produce a continuos map per sps and scenario
for (path in pathSMDoutputs) {

  for (species in targetSpecies) {
    
    for (scenario in scenarios) {
      
      # prep rasters' paths
      if (scenario == "Current") {
        tif_path <- file.path(
          path,
          paste0("proj_Current_EM_", gsub(" ", ".", species), "_continuous.tif"))
      } else {
        tif_path <- file.path(
          path,
          paste0("proj_", scenario, "_", extent, "_", gsub(" ", ".", species), "_continuous.tif"))
      }
      
      # only continue if the file exists
      if (!file.exists(tif_path)) {
        #message("File not found: ", tif_path)
        next  # skip to next iteration
      }
      
      # prep the rasters
      rast_obj <- rast(tif_path)
      rast_crop <- crop(rast_obj, iberia_vect)
      rast_mask <- mask(rast_crop, iberia_vect)
      
      cont_df <- as.data.frame(rast_mask, xy = TRUE, na.rm = TRUE)
      names(cont_df)[3] <- "suitability"
      cont_df$suitability <- cont_df$suitability/1000
      
      # plot SDM results
      p <- ggplot() +
        geom_sf(data = world, fill = "grey90", color = "black", linewidth = 0.2) +
        geom_tile(data = cont_df, aes(x = x, y = y, fill = suitability)) +
        tidyterra::scale_fill_terrain_c(name = "Suitability", limits = c(0, 1)) +
        coord_sf(xlim = c(iberia_bbox["xmin"], iberia_bbox["xmax"]),
                 ylim = c(iberia_bbox["ymin"], iberia_bbox["ymax"]),
                 expand = FALSE) +
        labs(x = "Longitude", y = "Latitude", 
             title = scenario_labels[scenario]) +
        theme_minimal()
      
      # Store plot
      continuous_list[[paste(species, scenario, sep = "_")]] <- p
    }
  }
}

# second for loop to organise plots to produce one figure per sps
species_layouts <- list()

for (species in targetSpecies) {
  
  message("Builiding continuous SDM landscapes for ", species)
  
  # Extract plots for this species
  p_current <- continuous_list[[paste(species, "Current", sep = "_")]]
  p_126_2030 <- continuous_list[[paste(species, "ssp126_2030", sep = "_")]]
  p_126_2050 <- continuous_list[[paste(species, "ssp126_2050", sep = "_")]]
  p_585_2030 <- continuous_list[[paste(species, "ssp585_2030", sep = "_")]]
  p_585_2050 <- continuous_list[[paste(species, "ssp585_2050", sep = "_")]]
  
  
  # If any are NULL, skip this species
  if (is.null(p_current)) {
    message("Skipping ", species, " because SDM outputs are missing.")
    next
  }
  
  # Patchwork layout:
  # Current spans two rows, other scenarios stacked in two columns
  species_plot <- 
    (p_current | (p_126_2030 / p_585_2030) | (p_126_2050 / p_585_2050))  +
    plot_layout(widths = c(1, 1, 1),
                guides = "collect") +
    plot_annotation(title = bquote(italic(.(species)))) &
    theme(legend.position = "bottom") 
  
  # Store the per-species patchwork
  species_layouts[[species]] <- species_plot
  
  # ggsave
  ggsave(plot = species_plot,
         file = file.path(SDMsFigures_dir, "ContinuousLandscapes", paste0("SDMlandscapes_", species, ".png")),
         bg = 'white', width = 300, height = 150, units = "mm", dpi = 1200,
         #compression ="lzw"
         )
}

invisible(gc())

############
# OUTPUT 5 # Binary maps from SDM outputs
############

binary_list <- list()

# first for loop to produce a continuos map per sps and scenario

for (path in pathSMDoutputs) {
  
  for (species in targetSpecies) {
    
    for (scenario in scenarios) {
      
      # prep rasters' paths
      if (scenario == "Current") {
        tif_path <- file.path(
          path,
          paste0("proj_Current_EM_", gsub(" ", ".", species), "_binary.tif"))
      } else {
        tif_path <- file.path(
          path,
          paste0("proj_", scenario, "_", extent, "_", gsub(" ", ".", species), "_binary.tif"))
      }
      
      # only continue if the file exists
      if (!file.exists(tif_path)) {
        #message("File not found: ", tif_path)
        next  # skip to next iteration
      }
      
      # prep the rasters
      rast_obj <- rast(tif_path)
      rast_crop <- crop(rast_obj, iberia_vect)
      rast_mask <- mask(rast_crop, iberia_vect)
      
      bin_df <- as.data.frame(rast_mask, xy = TRUE, na.rm = TRUE)
      names(bin_df)[3] <- "suitableArea"
      
      # plot SDM results
      p_bin <- ggplot() +
        geom_sf(data = world, fill = "grey95", color = "black") +
        geom_raster(data = bin_df, aes(x = x, y = y, fill = factor(suitableArea)), alpha = 0.6) +
        coord_sf(xlim = c(iberia_bbox["xmin"], iberia_bbox["xmax"]),
                 ylim = c(iberia_bbox["ymin"], iberia_bbox["ymax"]),
                 expand = FALSE) +
        scale_fill_manual(
          name = NULL,
          values = c("0" = "transparent", "1" = "green4"),
          labels = c(" ","Suitable Area")
        ) +
        labs(x = "Longitude", y = "Latitude", title = scenario_labels[scenario]) +
        theme_minimal() +
        theme(legend.position = "bottom")
    
      
      # Store plot
      binary_list[[paste(species, scenario, sep = "_")]] <- p_bin
    }
  }
}

# second for loop to organise plots to produce one figure per sps
species_layouts <- list()

for (species in targetSpecies) {
  
  message("Builiding binary landscapes for ", species)
  
  # Extract plots for this species
  p_current <- binary_list[[paste(species, "Current", sep = "_")]]
  p_126_2030 <- binary_list[[paste(species, "ssp126_2030", sep = "_")]]
  p_126_2050 <- binary_list[[paste(species, "ssp126_2050", sep = "_")]]
  p_585_2030 <- binary_list[[paste(species, "ssp585_2030", sep = "_")]]
  p_585_2050 <- binary_list[[paste(species, "ssp585_2050", sep = "_")]]
  
  # If any are NULL, skip this species
  if (is.null(p_current)) {
    message("Skipping ", species, " because SDM outputs are missing.")
    next
  }
  
  # Patchwork layout:
  # Current spans two rows, other scenarios stacked in two columns
  species_binary_plot <- 
    (p_current | (p_126_2030 / p_585_2030) | (p_126_2050 / p_585_2050))  +
    plot_layout(widths = c(1, 1, 1),
                guides = "collect") +
    plot_annotation(title = bquote(italic(.(species)))) &
    theme(legend.position = "bottom") 
  
  # Store the per-species patchwork
  species_layouts[[species]] <- species_binary_plot
  
  # ggsave
  ggsave(plot = species_binary_plot,
         file = file.path(SDMsFigures_dir, "BinaryLandscapes", paste0("BinaryLandscapes_", species, ".png")),
         bg = 'white', width = 300, height = 150, units = "mm", dpi = 1200,
         #compression ="lzw"
         )
}

# 
# ############
# # OUTPUT 6 # Hotspot areas maps
# ############
# 
# # all files
# mammal_files <- list.files(mammalsOutputs, pattern = "_binary.tif",
#                            full.names = TRUE, recursive = TRUE)
# bird_files <- list.files(birdsOutputs, pattern = "_binary.tif",
#                          full.names = TRUE, recursive = TRUE)
# binary_files <- c(mammal_files, bird_files)
# 
# 
# # current scenario files
# current_files <- binary_files[grep("proj_Current_EM", binary_files)]
# 
# # future scenario files
# future_ssp126 <- binary_files[grep("ssp126", binary_files)]
# future_ssp585 <- binary_files[grep("ssp585", binary_files)]
# 
# 
# # paths to binary rasters
# raster_paths <- list(
#   Current = current_files,
#   ssp126_2030 = future_ssp126[grepl("2030", future_ssp126)],
#   ssp126_2050 = future_ssp126[grepl("2050", future_ssp126)],
#   ssp585_2030 = future_ssp585[grepl("2030", future_ssp585)],
#   ssp585_2050 = future_ssp585[grepl("2050", future_ssp585)] 
# )
# 
# # empty list to store plots
# hotspots_list <- list()
# 
# # overlap all rasters per scenario and year
# for (scenario in names(raster_paths)) {
#   
#   # stack and sum binary rasters
#   r_sum <- sum(rast(raster_paths[[scenario]]))
#   
#   # crop + mask to Iberia extent (to avoid seeing norther africa in colors)
#   r_crop <- crop(r_sum, iberia_vect)   
#   r_mask <- mask(r_crop, iberia_vect)  
#   
#   # convert to dataframe
#   r_df <- as.data.frame(r_mask, xy = TRUE)
#   
#   # plot
#   p <- ggplot(r_df) +
#     geom_sf(data = world, fill = "grey90", color = "black", linewidth = 0.2) +
#     geom_tile(aes(x = x, y = y, fill = sum)) +
#     coord_sf(xlim = c(iberia_bbox["xmin"], iberia_bbox["xmax"]),
#              ylim = c(iberia_bbox["ymin"], iberia_bbox["ymax"]),
#              expand = FALSE)  +
#     scale_fill_distiller(palette = "Spectral",
#                          direction = -1,
#                          name = "Estimated no of\nthreatened species",
#                          limits = c(0, 16)) +
#     labs(x = "Longitude", y = "Latitude", title = scenario_labels[scenario]) +
#     theme_minimal() +
#     theme(legend.position = "bottom") +
#     guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))
#   
#   # store plots in a list for later
#   hotspots_list[[scenario]] <- p
# }
# 
# # extract plots for each scenario
# hotspots_current <- hotspots_list[["Current"]]
# hotspots_126_2030 <- hotspots_list[["ssp126_2030"]]
# hotspots_126_2050 <- hotspots_list[["ssp126_2050"]]
# hotspots_585_2030 <- hotspots_list[["ssp585_2030"]]
# hotspots_585_2050 <- hotspots_list[["ssp585_2050"]]
# 
# # display plots as needed
# species_plot <- 
#   (hotspots_current | (hotspots_126_2030 / hotspots_585_2030) | (hotspots_126_2050 / hotspots_585_2050))  +
#   plot_layout(widths = c(1, 1, 1),
#               guides = "collect") &
#   theme(legend.position = "bottom") 
# 
# # save final figure
# ggsave(plot = species_plot,
#        file = file.path(outputPathSDMfigures, paste0("HotspotAreasMap_", Sys.Date(), ".tif")),
#        bg = 'white', width = 300, height = 150, units = "mm", dpi = 1200, compression ="lzw")
# 
# ggsave(plot = species_plot,
#        file = file.path(outputPathSDMfigures, paste0("HotspotAreasMap_", Sys.Date(), ".png")),
#        bg = 'white', width = 300, height = 150, units = "mm", dpi = 1200
#        #, compression ="lzw"
#        )
