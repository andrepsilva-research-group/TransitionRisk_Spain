## Name: SDMfigures.R ##
## Author: Jorinde-M. Rieger & Inês Silva ##
## Description: Create all suplementary material's figure from the SDM's runs ##
## Date: 30th Aug 2025 ##


# paths for outputs
mammalsOutputs <- "C:/Users/maria/Desktop/mammalOutputs" #REPLACE WITH YOUR OWN PATH
birdsOutputs <- "C:/Users/maria/Desktop/birdOutputs" # REPLACE WITH YOUR OWN PATH
pathSMDoutputs <- c(mammalsOutputs, birdsOutputs)

#pathSMDoutputs <- "C:/Users/maria/Desktop/birdOutputs"

# new folder to save figures
outputPathSDMfigures <- "C:/Users/maria/Desktop/SDMoutputs31Aug2025" # adapt if needed
if (!dir.exists(outputPathSDMfigures)) {
  dir.create(outputPathSDMfigures, recursive = TRUE)
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

targetSpecies <- c(## BIRDS ##
                   "Streptopelia turtur", "Otis tarda", "Chlamydotis undulata",
                   "Neophron percnopterus", "Gyps rueppellii", "Aquila adalberti",
                   "Lanius meridionalis", "Chersophilus duponti", "Fringilla polatzeki",
                   "Pyrrhula murina", 
                   ## MAMMALS ##
                   "Crocidura canariensis", "Eliomys quercinus", "Galemys pyrenaicus",
                   "Lepus castroviejoi", "Lynx pardinus", "Myotis capaccinii",
                   "Nyctalus lasiopterus", "Nyctalus azoreum", "Pipistrellus maderensis",
                   "Rhinolophus mehelyi", "Oryctolagus cuniculus", "Plecotus teneriffae",
                   "Miniopterus schreibersii")

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

presencePlots <- list()

#species <- "Streptopelia turtur"

for (path in pathSMDoutputs) {

  for (species in targetSpecies) {
    
    message(paste0("Building Presence Plot for ", species))
    # file path for the presence .csv
    csv_file <- file.path(path,
                          paste0("PresencePoints_", gsub(" ", ".", species), "_", extent, ".csv"))
    
    if (file.exists(csv_file)) {
      # read presence points
      presence_df <- readr::read_csv(csv_file, show_col_types = FALSE)
      colnames(presence_df)[1:2] <- c("Longitude", "Latitude")
      presence_df$Type <- "Presence Points"
      
      # Plot occurrences over world map
      p <- ggplot() +
        geom_sf(data = world, fill = "grey95", color = "black") +
        geom_point(
          data = presence_df,
          aes(x = Longitude, y = Latitude),
          color = "red", size = 0.3, alpha = 0.7
        ) +
        ggtitle(bquote("Presence points of " * italic(.(species)) * " since 2015")) +
        theme_minimal()
      
      # save all  species plots inside a list
      presencePlots[[gsub(" ", ".", species)]] <- p
      
      # save plot as image
      ggsave(plot = p,
             file = file.path(outputPathSDMfigures, "presencePlots", paste0("PresencesPlot_", gsub(" ", ".", species), ".png")),
             bg = 'white', width = 300, height = 150, units = "mm", dpi = 1200, 
             #compression ="lzw"
             )
    }
  }
}

############
# OUTPUT 2 # Variable Importance for Ensemble Models Table
############

all_var_importance_em <- data.frame() # create data frame for all variable importance of EM

for (path in pathSMDoutputs) {
  for (species in targetSpecies) {
    file_path <- file.path(path,
                           paste0("VarImportanceEM_", gsub(" ", ".", species), "_", extent, ".csv"))
    if (file.exists(file_path)) {
      var_imp <- read.csv(file_path)
      var_imp$species <- species  # Add species column if not present
      all_var_importance_em <- rbind(all_var_importance_em, var_imp)
    }
  }
}

# Filter for algo == "EMmean"
var_importance_em_filtered <- all_var_importance_em %>%
  filter(algo == "EMmean")

# Summarize mean and SD for each species and variable
importance_summary <- var_importance_em_filtered %>%
  group_by(species, expl.var) %>%
  dplyr::summarize(
    mean_importance = mean(var.imp, na.rm = TRUE), # Calculate mean
    sd_importance = sd(var.imp, na.rm = TRUE), # Calculate standard deviation
    .groups = "drop") %>%
  pivot_longer(cols = c(mean_importance, sd_importance), names_to = "metrics", values_to = "importance") %>%
  mutate(metrics = ifelse(metrics == "mean_importance", "Mean", "SD")) %>%
  pivot_wider(names_from = expl.var, values_from = importance) %>%
  dplyr::arrange(species, metrics)

# write table to .csv and .xslx (for easy copy paste later)
write.csv(importance_summary,
          file = file.path(outputPathSDMfigures,
                           paste0("VariableImportanceSummaryTable", Sys.Date(), ".csv")), row.names = FALSE)
writexl::write_xlsx(importance_summary,
                    path = file.path(outputPathSDMfigures,
                                     paste0("VariableImportanceSummaryTable", Sys.Date(), ".xlsx")))

############
# OUTPUT 3 # Evaluation plots for Ensemble Models
############

all_EvalScoresEM <- data.frame() 

for (path in pathSMDoutputs) {
  
  for (species in targetSpecies) {
    
    file_path <- file.path(path, paste0("EvalScoresEM_", gsub(" ", ".", species), "_", extent, ".csv"))
    
    if (file.exists(file_path)) {
      eval_scores <- read.csv(file_path)
      eval_scores$species <- species  # Add species column if not present
      all_EvalScoresEM <- rbind(all_EvalScoresEM, eval_scores)
    }
  }
}
  
eval_plot_em <- ggplot(all_EvalScoresEM, aes(x = species, y = calibration, fill = metric.eval)) +
  geom_boxplot() +
  labs(title = "Evaluation Metrics for Ensemble Models (EM)",
       x = "Species",
       y = "Calibration Score",
       fill = "Metric") +
  scale_fill_viridis_d(name = "Metric") +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(plot = eval_plot_em,
       file = file.path(outputPathSDMfigures, "EvaluationMetrics_EM_exampleSpecies.tif"),
       bg = 'white', width = 300, height = 200, units = "mm", dpi = 1200, compression ="lzw")


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
          paste0("proj_", scenario, "_", gsub(" ", ".", species), "_continuous.tif"))
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
      
      # plot SDM results
      p <- ggplot() +
        geom_sf(data = world, fill = "grey90", color = "black", linewidth = 0.2) +
        geom_tile(data = cont_df, aes(x = x, y = y, fill = suitability)) +
        tidyterra::scale_fill_terrain_c(name = "Suitability", limits = c(0, 1000)) +
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
         file = file.path(outputPathSDMfigures, "ContinuousLandscapes", paste0("SDMlandscapes_", species, ".png")),
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
          paste0("proj_", scenario, "_", gsub(" ", ".", species), "_binary.tif"))
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
         file = file.path(outputPathSDMfigures, "BinaryLandscapes", paste0("BinaryLandscapes_", species, ".png")),
         bg = 'white', width = 300, height = 150, units = "mm", dpi = 1200,
         #compression ="lzw"
         )
}


############
# OUTPUT 6 # Hotspot areas maps
############

# all files
mammal_files <- list.files(mammalsOutputs, pattern = "_binary.tif",
                           full.names = TRUE, recursive = TRUE)
bird_files <- list.files(birdsOutputs, pattern = "_binary.tif",
                         full.names = TRUE, recursive = TRUE)
binary_files <- c(mammal_files, bird_files)


# current scenario files
current_files <- binary_files[grep("proj_Current_EM", binary_files)]

# future scenario files
future_ssp126 <- binary_files[grep("ssp126", binary_files)]
future_ssp585 <- binary_files[grep("ssp585", binary_files)]


# paths to binary rasters
raster_paths <- list(
  Current = current_files,
  ssp126_2030 = future_ssp126[grepl("2030", future_ssp126)],
  ssp126_2050 = future_ssp126[grepl("2050", future_ssp126)],
  ssp585_2030 = future_ssp585[grepl("2030", future_ssp585)],
  ssp585_2050 = future_ssp585[grepl("2050", future_ssp585)] 
)

# empty list to store plots
hotspots_list <- list()

# overlap all rasters per scenario and year
for (scenario in names(raster_paths)) {
  
  # stack and sum binary rasters
  r_sum <- sum(rast(raster_paths[[scenario]]))
  
  # crop + mask to Iberia extent (to avoid seeing norther africa in colors)
  r_crop <- crop(r_sum, iberia_vect)   
  r_mask <- mask(r_crop, iberia_vect)  
  
  # convert to dataframe
  r_df <- as.data.frame(r_mask, xy = TRUE)
  
  # plot
  p <- ggplot(r_df) +
    geom_sf(data = world, fill = "grey90", color = "black", linewidth = 0.2) +
    geom_tile(aes(x = x, y = y, fill = sum)) +
    coord_sf(xlim = c(iberia_bbox["xmin"], iberia_bbox["xmax"]),
             ylim = c(iberia_bbox["ymin"], iberia_bbox["ymax"]),
             expand = FALSE)  +
    scale_fill_distiller(palette = "Spectral",
                         direction = -1,
                         name = "Estimated no of\nthreatened species",
                         limits = c(0, 16)) +
    labs(x = "Longitude", y = "Latitude", title = scenario_labels[scenario]) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(fill = guide_colorbar(barwidth = 15, barheight = 0.5))
  
  # store plots in a list for later
  hotspots_list[[scenario]] <- p
}

# extract plots for each scenario
hotspots_current <- hotspots_list[["Current"]]
hotspots_126_2030 <- hotspots_list[["ssp126_2030"]]
hotspots_126_2050 <- hotspots_list[["ssp126_2050"]]
hotspots_585_2030 <- hotspots_list[["ssp585_2030"]]
hotspots_585_2050 <- hotspots_list[["ssp585_2050"]]

# display plots as needed
species_plot <- 
  (hotspots_current | (hotspots_126_2030 / hotspots_585_2030) | (hotspots_126_2050 / hotspots_585_2050))  +
  plot_layout(widths = c(1, 1, 1),
              guides = "collect") &
  theme(legend.position = "bottom") 

# save final figure
ggsave(plot = species_plot,
       file = file.path(outputPathSDMfigures, paste0("HotspotAreasMap_", Sys.Date(), ".tif")),
       bg = 'white', width = 300, height = 150, units = "mm", dpi = 1200, compression ="lzw")

ggsave(plot = species_plot,
       file = file.path(outputPathSDMfigures, paste0("HotspotAreasMap_", Sys.Date(), ".png")),
       bg = 'white', width = 300, height = 150, units = "mm", dpi = 1200
       #, compression ="lzw"
       )
