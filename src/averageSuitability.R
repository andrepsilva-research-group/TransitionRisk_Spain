## Name: averageSuitability.R ##
## Authors: Inês Silva
## Description: Produces maps of average suitability per pixel to add to an overall index ##


# paths for outputs
mammalsOutputs <- "./output/mammalOutputs" #REPLACE WITH YOUR OWN PATH
birdsOutputs <- "./output/birdOutputs" # REPLACE WITH YOUR OWN PATH
pathSMDoutputs <- c(mammalsOutputs, birdsOutputs)

#pathSMDoutputs <- "C:/Users/maria/Desktop/birdOutputs"

# new folder to save figures
outputPath_averageSuitbaility <- "./output/averageSuitability_19Jan2026" # adapt if needed
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

# species currently with SDM outputs (19th Jan 2026)
targetSpecies <- c(
  ## BIRDS ##
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


## --- MAMMALS --- ##

# How to obtain an IUCN shapefile ----------------------------------------------
# - Go to the IUCN Red List website: https://www.iucnredlist.org/resources/spatial-data-download
# - Choose the target species group (e.g. "Terrestrial Mammals")
# - Click the "Polygon download" button (you will need to create a free IUCN account)
# - Follow the website instruction for download
# - Once on your computer, extract the downloaded .zip file and move it into your data folder
# ------------------------------------------------------------------------------

sps_group <- "mammals"

# read a IUCN shapefile
sps <- sf::st_read("./data/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")

# filter for threatened species
sps_threatened <- sps %>%
  filter(category %in% c("CR", "EN", "VU"))
invisible(gc(rm(sps)))

# choose a target area
extent_sf <- rnaturalearth::ne_countries(
  scale = "medium",
  country = c("Spain", "Portugal"),
  returnclass = "sf")

# intersect IUCN file with target area
sf_use_s2(FALSE) # work on a flat earth
sps_extent <- st_intersection(sps_threatened, extent_sf)
invisible(gc(rm(extent_sf)))

# inspect species and REMOVE or ADD species manually 
target_sps <- unique(sps_extent$sci_name)
target_sps <- target_sps[target_sps != "Macaca sylvanus"] # REMOVE Macaca sylvanus (non-native species)


# save a .csv file with the taxonomic info
sps_names <- sps_extent %>% 
  dplyr::select(id_no, sci_name, order_, family, genus, legend, category) %>% 
  dplyr::filter(sci_name %in% target_sps) %>% 
  st_drop_geometry() %>% 
  unique()

## --- BIRDS --- ##

# How to obtain a Bird species list from BirdLife ------------------------------
# - Go to BirdLife International website: https://datazone.birdlife.org/search
# - Choose the filters on the left hand side to choose the target area/criteria 
# - Click the "Download" button 
# - Follow the website instruction for download
# - Once on your computer, extract the downloaded .zip file and move it into your data folder
# ------------------------------------------------------------------------------

sps_group <- "birds"

# read csv species list
sps <- read.csv("./data/BIRDLIFE_spsList/species-filter-results.csv")

# filter for threatened species
sps_threatened <- sps %>%
  filter(RL.Category %in% c("CR", "EN", "VU"))
invisible(gc(rm(sps)))

# inspect species and REMOVE or ADD species manually
target_sps <- unique(sps_threatened$Scientific.name)

# get taxonomic info
taxo_info <- taxize::classification(target_sps, db = "gbif")
taxo_df <- do.call(rbind, lapply(names(taxo_info), function(sp) {
  cl <- taxo_info[[sp]]
  if (is.null(cl)) return(NULL)
  data.frame(
    sci_name = sp,
    order_ = cl$name[cl$rank == "order"],
    family = cl$name[cl$rank == "family"],
    genus = cl$name[cl$rank == "genus"],
    stringsAsFactors = FALSE)
}))
invisible(gc(rm(taxo_info)))

sps_names_bird <- left_join(sps_threatened %>% rename(sci_name = Scientific.name),
                       taxo_df, by = "sci_name") %>% 
  dplyr::select(SIS.ID, sci_name, order_, family, genus, Common.name, RL.Category) %>% 
  dplyr::filter(sci_name %in% targetSpecies)


sps_names_std <- sps_names %>%
  dplyr::rename(id_no = id_no,
                sci_name = sci_name,
                order_ = order_,
                family = family,
                genus = genus,
                common_name = legend,
                conservation_status = category)

sps_names_bird_std <- sps_names_bird %>%
  dplyr::rename(id_no = SIS.ID,
                sci_name = sci_name,
                order_ = order_,
                family = family,
                genus = genus,
                common_name = Common.name,
                conservation_status = RL.Category)

sps_names_all <- dplyr::bind_rows(sps_names_std,sps_names_bird_std)


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
          paste0("proj_", scenario, "_", gsub(" ", ".", species), "_continuous.tif")
        )
      }
      
      if (!file.exists(tif_path)) next
      
      # read + spatial prep
      r <- rast(tif_path)
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

# average rasters
mean_rasters <- list()
for (scenario in names(scenario_rasters)) {
  
  # stack rasters
  r_stack <- rast(scenario_rasters[[scenario]])
  
  # pixel-wise mean
  mean_rasters[[scenario]] <- mean(r_stack, na.rm = TRUE)
}

#mean_rasters$Current
#plot(mean_rasters$Current)

plot_list <- list()

for (scenario in names(mean_rasters)) {
  
  mean_df <- as.data.frame(mean_rasters[[scenario]],
                           xy = TRUE, na.rm = TRUE)
  names(mean_df)[3] <- "meanSuitability"
  mean_df$meanSuitability <- mean_df$meanSuitability / 1000
  
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
       file = "./output/AverageSuitbaility_MammalBirdSpecies_21Jan2026.png",
       bg = 'white', width = 300, height = 150, units = "mm", dpi = 1200 
       #compression ="lzw"
)

#########################################################
# BASED ON CONSERVATION STATUS AVERAGE SUITABILITY MAPS #
#########################################################

# named vector: species -> conservation category
sps_to_cat <- sps_names_all %>%
  dplyr::select(sci_name, conservation_status) %>%
  distinct() %>%
  tibble::deframe()

sps_names_all %>% 
  group_by(conservation_status) %>% 
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
          paste0("proj_", scenario, "_", gsub(" ", ".", species), "_continuous.tif")
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

#plot(status_scenario_rasters$VU$Current$`Lepus castroviejoi`)

mean_status_rasters <- list()

for (status in names(status_scenario_rasters)) {
  
  mean_status_rasters[[status]] <- list()
  
  for (scenario in names(status_scenario_rasters[[status]])) {
    
    r_stack <- rast(status_scenario_rasters[[status]][[scenario]])
    
    # average across species and rescale 0–1
    mean_status_rasters[[status]][[scenario]] <-
      mean(r_stack, na.rm = TRUE) / 1000
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
