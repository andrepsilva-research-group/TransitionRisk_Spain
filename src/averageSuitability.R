## Name: averageSuitability.R ##
## Author: Inês Silva ##
## Date: 30th June 2025 ##
## Description: Produces maps of average suitability per pixel to add to an
## overall index for spatial prioritization analysis. ##


source("src/libraries.R")

# paths for outputs
mammalsOutputs <- "E:/BankOfSpain_june2026/SDMs_mammals"
birdsOutputs <- "E:/BankOfSpain_june2026/SDMs_birds"
pathSMDoutputs <- c(mammalsOutputs, birdsOutputs)

#pathSMDoutputs <- "C:/Users/maria/Desktop/birdOutputs"

# new folder to save figures
outputPath_averageSuitbaility <- "./output/averageSuitability_June2026" # adapt if needed
if (!dir.exists(outputPath_averageSuitbaility)) {
  dir.create(outputPath_averageSuitbaility, recursive = TRUE)
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

library(readr)
iberian_threatened_mammals <- read_csv("data/iberian_threatened_mammals_taxonomy.csv")
iberian_threatened_birds <- read_csv("data/iberian_threatened_birds_taxonomy.csv") %>% 
  rename("category" = "RL.Category")

targetSpecies <- c(iberian_threatened_mammals$sci_name, iberian_threatened_birds$sci_name)
# 43 but we only have sdms for less

# get basemaps objects
extent <- "IberianPeninsula"
world <- ne_countries(scale = "medium", returnclass = "sf")
iberia <- world[world$name %in% c("Portugal", "Spain"), ]
# vector for masking
iberia_vect <- vect(iberia)   
# iberian peninsula bbox
iberia_bbox <- st_bbox(iberia) 
gc()


#################################################
# ALL SPECIES TOGETHER AVERAGE SUITABILITY MAPS #
#################################################

# first for loop to produce a continuos map per sps and scenario
# store rasters by scenario
scenario_rasters <- list()
for (path in pathSMDoutputs) {
  
  for (species in targetSpecies) {
    
    for (scenario in scenarios) {
      
      # prep raster path
      if (scenario == "Current") {
        tif_path <- file.path(
          path,
          paste0("proj_Current_EM_", gsub(" ", ".", species), "_continuous.tif")
        )
      } else {
        tif_path <- file.path(
          path,
          paste0("proj_", scenario, "_IberianPeninsula_", gsub(" ", ".", species), "_continuous.tif")
        )
      }
      
      if (!file.exists(tif_path)) next
      
      # read + spatial prep
      r <- rast(tif_path) / 1000
      r <- crop(r, iberia_vect)
      r <- mask(r, iberia_vect)
      
      # initialize list for scenario if needed
      if (is.null(scenario_rasters[[scenario]])) {
        scenario_rasters[[scenario]] <- list()
      }
      
      # store raster
      scenario_rasters[[scenario]][[species]] <- r
    }
  }
}

#plot(scenario_rasters$ssp126_2030$`Eliomys quercinus`)
  
  
# average rasters
mean_rasters <- list()
for (scenario in names(scenario_rasters)) {
  
  # stack rasters
  r_stack <- rast(scenario_rasters[[scenario]])
  
  mean_raster_scenario <- mean(r_stack, na.rm = TRUE)
  
  # pixel-wise mean
  mean_rasters[[scenario]] <- mean_raster_scenario
  
  writeRaster(mean_raster_scenario,
              paste0("output/scenarioAverages/", scenario, "_mean.tif"),
              overwrite = TRUE)
  
}

#mean_rasters$Current
#plot(mean_rasters$Current)

plot_list <- list()

for (scenario in names(mean_rasters)) {
  
  mean_df <- as.data.frame(mean_rasters[[scenario]],
                           xy = TRUE, na.rm = TRUE)
  names(mean_df)[3] <- "meanSuitability"
  mean_df$meanSuitability <- mean_df$meanSuitability
  
  p <- ggplot() +
    geom_sf(data = world, fill = "grey95", color = "black") +
    geom_tile(
      data = mean_df,
      aes(x = x, y = y, fill = meanSuitability)) +
    coord_sf(
      xlim = c(iberia_bbox["xmin"], iberia_bbox["xmax"]),
      ylim = c(iberia_bbox["ymin"], iberia_bbox["ymax"]),
      expand = FALSE) +
    scale_fill_distiller(
      palette = "RdYlGn",
      direction = 1,
      #limits = c(0,1),
      name = "Mean suitability") +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = scenario_labels[scenario]) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  plot_list[[scenario]] <- p
}

# display plots as needed
species_plot <- 
  (plot_list$Current | (plot_list$ssp126_2030 / plot_list$ssp585_2030) | (plot_list$ssp126_2050 / plot_list$ssp585_2050))  +
 plot_annotation(
    title = 'Average Suitability Values',
    subtitle = 'These plots include both bird and mammals species suitability values',
    caption = 'total no of species: 21'  ) +
  plot_layout(widths = c(1, 1, 1),
              guides = "collect") &
  theme(legend.position = "bottom") 

# save plot as image
ggsave(plot = species_plot,
       file = "./output/AverageSuitbaility_MammalBirdSpecies_June2026.png",
       bg = 'white', width = 300, height = 150, units = "mm", dpi = 1200 
       #compression ="lzw"
)

#########################################################
# BASED ON CONSERVATION STATUS AVERAGE SUITABILITY MAPS #
#########################################################

sps_names_all <- bind_rows(iberian_threatened_mammals, iberian_threatened_birds)
# named vector: species -> conservation category
sps_to_cat <- sps_names_all %>%
  dplyr::select(sci_name, category) %>%
  distinct() %>%
  tibble::deframe()

sps_names_all %>% 
  group_by(category) %>% 
  summarise(n = n())

# empty list
status_scenario_rasters <- list()

for (path in pathSMDoutputs) {
  for (species in targetSpecies) {
    # skip species without category
    if (!species %in% names(sps_to_cat)) next
    
    # get all sps with the same conservation status
    status <- sps_to_cat[[species]]
    
    for (scenario in scenarios) {
      # get current raster
      if (scenario == "Current") {
        tif_path <- file.path(
          path,
          paste0("proj_Current_EM_", gsub(" ", ".", species), "_continuous.tif")
        )
      } else {
        # get future scenarios rasters
        tif_path <- file.path(
          path,
          paste0("proj_", scenario, "_IberianPeninsula_", gsub(" ", ".", species), "_continuous.tif")
        )
      }
      
      if (!file.exists(tif_path)) next
      # crop by the iberian peninsula
      r <- rast(tif_path) |>
        crop(iberia_vect) |>
        mask(iberia_vect)
      
      # store in a nested list
      status_scenario_rasters[[status]][[scenario]][[species]] <- r
    }
  }
}

status_scenario_rasters$CR$Current

#plot(status_scenario_rasters$VU$Current$`Lepus castroviejoi`)

mean_status_rasters <- list()

for (status in names(status_scenario_rasters)) {
  
  mean_status_rasters[[status]] <- list()
  
  for (scenario in names(status_scenario_rasters[[status]])) {
    
    r_stack <- rast(status_scenario_rasters[[status]][[scenario]])
    
    # average across species and rescale 0–1
    mean_raster <- mean(r_stack, na.rm = TRUE) / 1000
    
    # store in list
    mean_status_rasters[[status]][[scenario]] <- mean_raster
    
    writeRaster(mean_raster,
                paste0("output/conservationStatusAverages/", status, "_", scenario, "_mean.tif"),
                overwrite = TRUE)
  }
}

# empty list
conservStatus_plots <- list()

# build one plot per conservation status - scenario - year
for (status in names(mean_status_rasters)) {
  
  for (scenario in names(mean_status_rasters[[status]])) {
    
    r <- mean_status_rasters[[status]][[scenario]]
    
    mean_df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
    names(mean_df)[3] <- "meanSuitability"
    
    p <- ggplot() +
      geom_sf(data = world, fill = "grey95", color = "black", linewidth = 0.2) +
      geom_tile(
        data = mean_df,
        aes(x = x, y = y, fill = meanSuitability)) +
      coord_sf(
        xlim = c(iberia_bbox["xmin"], iberia_bbox["xmax"]),
        ylim = c(iberia_bbox["ymin"], iberia_bbox["ymax"]),
        expand = FALSE) +
      scale_fill_distiller(
        palette = "RdYlGn",
        direction = 1,
        limits = c(0, 0.8),
        oob = scales::squish,
        name = "Mean suitability") +
      labs(title = paste(scenario_labels[scenario]),
           x = "Longitude", y = "Latitude") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    conservStatus_plots[[paste(status, scenario, sep = "_")]] <- p
  }
}
conservStatus_plots$VU_Current
conservStatus_plots$VU_ssp126_2030
conservStatus_plots$VU_ssp126_2050
conservStatus_plots$VU_ssp585_2030
conservStatus_plots$VU_ssp585_2050
conservStatus_plots$EN_Current
conservStatus_plots$EN_ssp126_2030
conservStatus_plots$EN_ssp126_2050
conservStatus_plots$EN_ssp585_2030
conservStatus_plots$EN_ssp585_2030
conservStatus_plots$EN_ssp585_2050

############
# SSP1-2.6 #
############

ssp126 <-
  ((conservStatus_plots$VU_Current) |(conservStatus_plots$VU_ssp126_2030) | (conservStatus_plots$VU_ssp126_2050)) /
  (conservStatus_plots$EN_Current | conservStatus_plots$EN_ssp126_2030 | conservStatus_plots$EN_ssp126_2050) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5))

final_ssp126_plot <-
  ssp126 +
  plot_annotation(
    tag_levels = list(c("VUNERABLE\n STATUS (n=13)", "", "", "ENDANGERED\n STATUS (n=7)")),
    tag_suffix = "") &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 12
    )
  )

# save plot as image
ggsave(plot = final_ssp126_plot,
       file = "./output/AverageSuitbaility_MammalBirdSpeciesPerStatus_ssp126_21Jan2026.png",
       bg = 'white', width = 300, height = 150, units = "mm", dpi = 1200 
       #compression ="lzw"
)


############
# SSP5-8.5 #
############

ssp585 <-
  ((conservStatus_plots$VU_Current) |(conservStatus_plots$VU_ssp585_2030) | (conservStatus_plots$VU_ssp585_2050)) /
  (conservStatus_plots$EN_Current | conservStatus_plots$EN_ssp585_2030 | conservStatus_plots$EN_ssp585_2050) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5))

final_ssp585_plot <-
  ssp585 +
  plot_annotation(
    tag_levels = list(c("VUNERABLE\n STATUS (n=13)", "", "", "ENDANGERED\n STATUS (n=7)")),
    tag_suffix = "") &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 12
    )
  )

# save plot as image
ggsave(plot = final_ssp585_plot,
       file = "./output/AverageSuitbaility_MammalBirdSpeciesPerStatus_ssp585_21Jan2026.png",
       bg = 'white', width = 300, height = 150, units = "mm", dpi = 1200 
       #compression ="lzw"
)
