## Name: SDMfigures.R ##
## Author: Jorinde-M. Rieger ##
## Description: Creates figures of SDM results in R ##
## Date: August 7th 2025 ##

# Define the scenarios and species to plot
scenarios <- c("ssp126", "ssp585")
scenario_names <- c("SSP1-RCP2.6", "SSP5-RCP8.5")
years <- c(2030, 2050)

putputPathSDMensemble <- "output/SDMensemble/BoS_Mammals" # adapt to outputPath from SDMs
outputPathSDMfigures <- "output/SDMfigures/BoS_Mammals" # adapt if needed
if (!dir.exists(outputPathSDMfigures)) {
  dir.create(outputPathSDMfigures, recursive = TRUE)
}

# Define Extent
# Spain and Portugal extent
extent <- "Spain and Portugal"
extent_sf <- rnaturalearth::ne_countries(scale = "medium", country = c("Spain", "Portugal"), returnclass = "sf")
extent_crs <- sf::st_transform(extent_sf, crs = crs(trainingLandscapes)) # Ensure CRS consistency
extent_sp <- terra::vect(extent_crs) # Convert the sf to a spatial object

# Target species 
targetSpecies <- c("Nyctalus lasiopterus", "Myotis capaccinii", "Galemys pyrenaicus", "Miniopterus schreibersii", "Oryctolagus cuniculus", "Eliomys quercinus", "Lynx pardinus",
  "Crocidura canariensis","Plecotus teneriffae", "Lepus castroviejoi", "Pipistrellus maderensis","Nyctalus azoreum", "Rhinolophus mehelyi") # all species in this row caused an error (too little presence data?)

# 1. Figures: SDM evaluation metrics ------------------------------------------------
# a. Plot evaluation scores for ensemble models
all_EvalScoresEM <- data.frame() # create data frame for all evaluation scores of EM
for (species in targetSpecies) {
  file_path <- file.path(putputPathSDMensemble, paste0("EvalScoresEM_", species, extent, ".csv"))
  if (file.exists(file_path)) {
    eval_scores <- read.csv(file_path)
    eval_scores$species <- species  # Add species column if not present
    all_EvalScoresEM <- rbind(all_EvalScoresEM, eval_scores)
  }
}

eval_plot_em <- ggplot2::ggplot(all_EvalScoresEM, aes(x = species, y = calibration, fill = metric.eval)) +
  ggplot2::geom_boxplot() +
  ggplot2::labs(
    title = "Evaluation Metrics for Ensemble Models (EM)",
    x = "Species",
    y = "Calibration Score",
    fill = "Metric"
  ) +
  ggplot2::scale_fill_viridis_d(name = "Metric") +  
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

print(eval_plot_em)
ggplot2::ggsave(file.path(outputPathSDMfigures, "EvaluationMetrics_EM_exampleSpecies.png"),
       plot = eval_plot_em, width = 10, height = 6, dpi = 300)

# b. Plot evaluation metrics for RUN 1-6
all_eval_scores <- data.frame() # create data frame for all evaluation scores
for (species in targetSpecies) {
  file_path <- file.path(putputPathSDMensemble, paste0("EvalScores_", species, extent, ".csv"))
  if (file.exists(file_path)) {
    eval_scores <- read.csv(file_path)
    eval_scores$species <- species  # Add species column if not present
    all_eval_scores <- rbind(all_eval_scores, eval_scores)
  }
}

# Filter for ROC and TSS metrics
eval_scores_filtered <- all_eval_scores %>%
  dplyr::filter(metric.eval %in% c("ROC", "TSS"))

eval_plot_models <- ggplot2::ggplot(eval_scores_filtered, aes(x = metric.eval, y = validation, fill = algo)) +
  ggplot2::geom_boxplot(position = position_dodge(1)) +
  ggplot2::facet_wrap(~ species, nrow = 1) +
  ggplot2::labs(
    title = "Validation Scores (ROC and TSS) by Model",
    x = "Metric",
    y = "Validation Score",
    fill = "Algorithm"
  ) +
  ggplot2::scale_fill_viridis_d(name = "Model") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
    strip.text = element_text(size = 12)
  )

print(eval_plot_models)
ggplot2::ggsave(file.path(outputPathSDMfigures, "ValidationScores_Models_exampleSpecies.png"),
       plot = eval_plot_models, width = 14, height = 7, dpi = 300)

# 2. Table: Variable importance ------------------------------------------------  
all_var_importance_em <- data.frame() # create data frame for all variable importance of EM
for (species in targetSpecies) {
  file_path <- file.path(putputPathSDMensemble, paste0("VarImportanceEM_", species, extent, ".csv"))
  if (file.exists(file_path)) {
    var_imp <- read.csv(file_path)
    var_imp$species <- species  # Add species column if not present
    all_var_importance_em <- rbind(all_var_importance_em, var_imp)
  }
}
# Filter for algo == "EMmean"
var_importance_em_filtered <- all_var_importance_em %>%
  filter(algo == "EMmean")

# Summarize mean and SD for each species and variable
importance_summary <- var_importance_em_filtered %>%
  dplyr::group_by(species, expl.var) %>%
  dplyr::summarize(
    mean_importance = mean(var.imp, na.rm = TRUE), # Calculate mean
    sd_importance = sd(var.imp, na.rm = TRUE), # Calculate standard deviation
    .groups = "drop"
  )

# Reshape to wide format with mean and SD rows for each species
importance_table <- importance_summary %>%
  tidyr::pivot_longer(cols = c(mean_importance, sd_importance), names_to = "metrics", values_to = "importance") %>%
  dplyr::mutate(metrics = ifelse(metrics == "mean_importance", "Mean", "SD")) %>%
  tidyr::pivot_wider(names_from = expl.var, values_from = importance) %>%
  dplyr::arrange(species, metrics)
print(importance_table)
write.csv(importance_table, file = file.path(outputPathSDMfigures, "VariableImportanceSummary_exampleSpecies.csv"), row.names = FALSE)
writexl::write_xlsx(importance_table, path = file.path(outputPathSDMfigures, "VariableImportanceSummary_exampleSpecies.xlsx")) # Excel file

# 3. Current Landscapes ------------------------------------------------
currentSuitability_plots <- list()
for (species in targetSpecies) {
  raster_file <- file.path(gsub(" ", ".", species), "proj_Current", paste0("proj_Current_", gsub(" ", ".", species), "_ensemble.tif"))
    emmean_layer <- terra::rast(raster_file)
   # emmean_layer <- terra::raster[[grep("EMmean", names(raster))]]
    raster_df <- as.data.frame(emmean_layer, xy = TRUE, na.rm = TRUE)
    colnames(raster_df)[3] <- "value"
    current_plot <- ggplot2::ggplot(raster_df, aes(x = x, y = y, fill = value)) +
      ggplot2::geom_raster() +
      tidyterra::scale_fill_terrain_c(name = "Suitability") +
      ggplot2::labs(
        title = paste(species, "- Current Ensemble Suitability"),
        x = "Longitude",
        y = "Latitude"
      ) +
      ggplot2::coord_fixed() +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "right"
      )
  currentSuitability_plots[[species]] <- current_plot

  # Optionally save individual PNGs
  ggplot2::ggsave(
     filename = file.path(outputPathSDMfigures, paste0("CurrentSuitability_", species, ".png")),
    plot = current_plot,
    width = 10,
    height = 8,
    dpi = 300)
}

# Combine all current plots into a single grid
if (length(currentSuitability_plots) > 0) {
  combined_current_plot <- gridExtra::grid.arrange(grobs = currentSuitability_plots, ncol = 2)
  ggplot2::ggsave(
    filename = file.path(outputPathSDMfigures, "CurrentSuitability_AllSpecies.png"),
    plot = combined_current_plot,
    width = 16,
    height = 8,
    dpi = 300)
}

# 4. Figure: Presence Points for multiple species ------------------------------------------------
world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

presencePlots <- list()
for (species in targetSpecies) {
  # Build the file path for the presence CSV
  csv_file <- file.path(putputPathSDMensemble, paste0("PresencePoints_", species, "_Spain and Portugal.csv"))
  raster_file <- file.path(gsub(" ", ".", species), "proj_Current", paste0("proj_Current_", gsub(" ", ".", species), "_ensemble.tif"))
  if (file.exists(csv_file) && file.exists(raster_file)) {
    # Read presence points
    presence_df <- readr::read_csv(csv_file, show_col_types = FALSE)
    colnames(presence_df)[1:2] <- c("Longitude", "Latitude")
    presence_df$Type <- "Presence Points"

    # Read current suitability raster
    emmean_layer <- terra::rast(raster_file)
    raster_df <- as.data.frame(emmean_layer, xy = TRUE, na.rm = TRUE)
    colnames(raster_df)[3] <- "value"

    # Plot: raster as background, presence points on top
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = world_sf, fill = NA, color = "black", size = 0.5) + # Add world boundary as background
      ggplot2::geom_raster(data = raster_df, aes(x = x, y = y, fill = value), alpha = 0.8) +
      ggplot2::geom_point(data = presence_df, aes(x = Longitude, y = Latitude), color = "blue", size = 0.7) +
      tidyterra::scale_fill_terrain_c(name = "Suitability") +
      ggplot2::labs(
        title = species,
        x = "Longitude",
        y = "Latitude"
      ) +
       ggplot2::coord_sf(expand = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = "none"
      )
    presencePlots[[species]] <- p
  }
}

# Extract the legend from one of the plots
example_plot <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = world_sf, fill = NA, color = "black", size = 0.5) + # Add world boundary as background
      ggplot2::geom_raster(data = raster_df, aes(x = x, y = y, fill = value), alpha = 0.8) +
      ggplot2::geom_point(data = presence_df, aes(x = Longitude, y = Latitude), color = "blue", size = 0.7) +
      tidyterra::scale_fill_terrain_c(
    name = "Suitability",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.theme = element_text(size = 18),
      label.theme = element_text(size = 18),
      barwidth = unit(6, "cm"),   # Make the colorbar longer
      barheight = unit(0.5, "cm")  # Make the colorbar thicker
    )
  ) +
      ggplot2::labs(
        title = species,
        x = "Longitude",
        y = "Latitude"
      ) +
      ggplot2::coord_sf() +
      ggplot2::theme_minimal() +
      ggplot2::theme(
    legend.direction = "horizontal",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )
shared_legend <- cowplot::get_legend(example_plot)

# Arrange all plots in a grid (max 2 per row)
if (length(presencePlots) > 0) {
  combined_presence_plot <- gridExtra::grid.arrange(
    grobs = presencePlots,
    ncol = 2 # max 2 species per row
  )

  final_plot <- gridExtra::grid.arrange(
    combined_presence_plot,
    shared_legend,
    ncol = 1,  # Legend below the plots
    heights = unit(c(10, 1.5), "null")
  )

  ggplot2::ggsave(
    filename = file.path(outputPathSDMfigures, "PresencePointsWithSuitability_allSpecies.png"),
    plot = final_plot,
    width = 12, height = 6 + 1.5 * ceiling(length(presencePlots)/2), dpi = 300
  )
}

# 5. Figure: Current and Predicted Landscapes ------------------------------------------------
# Current Suitability Landscapes
current_plots <- list()
for (species in targetSpecies) {
  raster_file <- file.path(gsub(" ", ".", species), "proj_Current", paste0("proj_Current_", gsub(" ", ".", species), "_ensemble.tif"))
  if (file.exists(raster_file)) {
    emmean_layer <- terra::rast(raster_file)
    raster_df <- as.data.frame(emmean_layer, xy = TRUE, na.rm = TRUE)
    colnames(raster_df)[3] <- "value"
    currentPlot <- ggplot2::ggplot(raster_df, aes(x = x, y = y, fill = value)) +
      ggplot2::geom_raster() +
      tidyterra::scale_fill_terrain_c(name = "Prediction") +
      ggplot2::labs(title = NULL, x = "Longitude", y = "Latitude") +
      ggplot2::coord_sf(expand = FALSE) + # Ensure correct aspect ratio
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        axis.ticks = element_line(),
        panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_blank(),
        legend.position = "none"
      )
    current_plots[[species]] <- currentPlot
  }
}

# Future Suitability Landscapes
future_plots <- list()
for (year in years) {
  for (scenario in scenarios) {
    for (species in targetSpecies) {
      # Build folder and file names
      folder_name <- paste0("proj_", scenario, "_", year, "_", species)
      file_name <- paste0("proj_", scenario, "_", year, "_", species, "_", gsub(" ", ".", species), ".tif")
      raster_file <- file.path(gsub(" ", ".", species), folder_name, file_name)
      if (file.exists(raster_file)) {
        emmean_layer <- terra::rast(raster_file)
        raster_df <- as.data.frame(emmean_layer, xy = TRUE, na.rm = TRUE)
        colnames(raster_df)[3] <- "value"
        futurePlot <- ggplot2::ggplot(raster_df, aes(x = x, y = y, fill = value)) +
          ggplot2::geom_raster() +
          tidyterra::scale_fill_terrain_c(name = "Suitability",
                               guide = guide_colorbar(
                                 title.position = "top",
                                 title.theme = element_text(size = 22),
                                 label.theme = element_text(size = 20)
                               )) +
          ggplot2::labs(title = NULL, x = "Longitude", y = "Latitude") +
          ggplot2::coord_sf(expand = FALSE) + # Ensure correct aspect ratio
          ggplot2::theme_bw() +
          ggplot2::theme(
            axis.title = element_text(size = 22),
            axis.text = element_text(size = 20),
            axis.ticks = element_line(),
            panel.grid.major = element_line(color = "gray"),
            panel.grid.minor = element_blank(),
            legend.position = "none")
        future_plots[[paste0(scenario, "_", year, "_", species)]] <- futurePlot
      }
    }
  }
}

# Extract the legend from one of the plots
example_plot <- ggplot2::ggplot(raster_df, aes(x = x, y = y, fill = value))+
  ggplot2::geom_raster() +
  tidyterra::scale_fill_terrain_c(name = "Suitability",
                       guide = guide_colorbar(
                         title.position = "top",
                         title.theme = element_text(size = 22),
                         label.theme = element_text(size = 20),
                         barwidth = 20,
                         barheight = 1.5
                       )) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.direction = "horizontal",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20)
  )
shared_legend <- cowplot::get_legend(example_plot)

# Suitability Landscapes for single species and multiple years
for (species in targetSpecies) {
  # Prepare plots for the layout
  current_grob <- gtable::gtable_trim(ggplotGrob(current_plots[[species]]))
  ssp1_2030_grob <- gtable::gtable_trim(ggplotGrob(future_plots[[paste0("ssp126_2030_", species)]]))
  ssp1_2050_grob <- gtable::gtable_trim(ggplotGrob(future_plots[[paste0("ssp126_2050_", species)]]))
  ssp5_2030_grob <- gtable::gtable_trim(ggplotGrob(future_plots[[paste0("ssp585_2030_", species)]]))
  ssp5_2050_grob <- gtable::gtable_trim(ggplotGrob(future_plots[[paste0("ssp585_2050_", species)]]))

  # Titles for columns
  title_grobs <- list(
    grid::nullGrob(),
    grid::textGrob("Current", gp = gpar(fontsize = 20), just = "centre"),
    grid::textGrob("2030", gp = gpar(fontsize = 20), just = "centre"),
    grid::textGrob("2050", gp = gpar(fontsize = 20), just = "centre")
  )

  # Row labels
  row1_label <- gtable_trim::textGrob("SSP1-RCP2.6", rot = 90, gp = gpar(fontsize = 20))
  row2_label <- gtable_trim::textGrob("SSP5-RCP8.5", rot = 90, gp = gpar(fontsize = 20))

  # Layout matrix: 3 columns (Current, 2030, 2050) + 1 for row labels, 3 rows (title, ssp1, ssp5)
  layout_matrix <- matrix(
    c(
      1, 2, 3, 4,      # Title row
      5, 6, 7, 8,      # SSP1 row: label, current, 2030, 2050
      9, 10, 11, 12    # SSP5 row: label, current, 2030, 2050
    ),
    nrow = 3, byrow = TRUE
  )

  # Collect all grobs in order of layout_matrix
  all_grobs <- list(
    title_grobs[[1]], title_grobs[[2]], title_grobs[[3]], title_grobs[[4]],
    row1_label, current_grob, ssp1_2030_grob, ssp1_2050_grob,
    row2_label, current_grob, ssp5_2030_grob, ssp5_2050_grob
  )

  # Set column widths and row heights
  col_widths <- grid::unit.c(unit(1, "cm"), rep(unit(1, "null"), 3))
  row_heights <- grid::unit.c(unit(1, "null"), rep(unit(4, "null"), 2))

   final_plot <- gridExtra::grid.arrange(
    grobs = all_grobs,
    layout_matrix = layout_matrix,
    widths = col_widths,
    heights = row_heights,
    bottom = shared_legend,
    top = textGrob(species, gp = gpar(fontface = "italic", fontsize = 26))
  )

  ggplot2::ggsave(
    filename = file.path(outputPathSDMfigures, paste0("SuitabilityLandscapes_", species, ".png")),
    plot = final_plot,
    width = 24, height = 11, dpi = 300
  )
}
