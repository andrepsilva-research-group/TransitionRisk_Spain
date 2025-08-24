## Name: presenceMaps
## Author: Inês Silva
## Description: build maps with presence points overlapping
## Date: 21 Aug 2025


#setwd("C:/Users/maria/OneDrive - Universidade de Lisboa/bankOfSpain")
source("./src/libraries.R")


# Load data
speciesData <- read.csv("./data/GBIF_birds_30+occurrences_Spain_and_Portugal.csv")
myExpl_full <- rast("./data/landscapes/traininglandscapes_2015_5km.tif")


# Get list of unique species
speciesList <- unique(speciesData$species)

# Load country boundaries (Spain + Portugal)
world <- ne_countries(scale = "medium", returnclass = "sf")

output_folder <- "./output/presenceMaps"
if(!dir.exists(output_folder)){
  dir.create(output_folder, recursive = TRUE)
}

# Loop over each species
for (targetSpecies in speciesList) {
  
  message("Processing: ", targetSpecies)
  
  # Filter single species
  DataSingleSpecies <- speciesData %>%
    dplyr::filter(species == !!targetSpecies) %>%
    drop_na(decimalLongitude, decimalLatitude, year) %>%
    dplyr::filter(year >= 2015)
  
  # Assign cell IDs to each occurrence based on myExpl raster
  cellValues <- terra::extract(
    myExpl_full,
    cbind(DataSingleSpecies$decimalLongitude, DataSingleSpecies$decimalLatitude)
  )
  cellValues$cell <- terra::cellFromXY(myExpl_full,
                                       cbind(DataSingleSpecies$decimalLongitude, DataSingleSpecies$decimalLatitude))
  DataSingleSpecies <- cbind(DataSingleSpecies, cellValues)
  
  set.seed(123) # keep this here as long as we subset only x occ
  
  # for cells with bioclimatic vars values
  DataSingleSpecies_unique <- DataSingleSpecies %>%
    group_by(cell) %>%
    # keep most recent year
    slice_max(year, with_ties = FALSE) %>%  # or slice_head(n = 1) for the first
    ungroup() %>%
    dplyr::filter(complete.cases(.)) 
  #%>% # biomod excludes all cells that do not have any data
    # subset 30 occ
    #slice_sample(n = min(30, nrow(DataSingleSpecies)))
  
  # Plot occurrences on a simple country map
  p <- ggplot() +
    geom_sf(data = world, fill = "grey95", color = "black") +
    geom_point(
      data = DataSingleSpecies_unique,
      aes(x = decimalLongitude, y = decimalLatitude),
      color = "red", size = 0.3, alpha = 0.7
    ) +
    ggtitle(bquote("Presence points of " * italic(.(targetSpecies)) * " since 2015")) +
    theme_minimal()
  
  #print(p)
  
  # save the map
  ggsave(
   filename = file.path(output_folder, paste0("map_", gsub(" ", "_", targetSpecies), ".png")),
   plot = p, width = 6, height = 5
  )
}

