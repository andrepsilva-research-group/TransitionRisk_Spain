## Species selection for Iberian Peninsula ##
## Inês Silva & Jorinde-M. Rieger ##
## 06 June 2026 ##

source("./src/libraries.R")
source("./src/customFunctions2.R")

# choose a target area (Iberian Peninsula)
extent_sf <- rnaturalearth::ne_countries(
  scale = "medium",
  country = c("Spain", "Portugal"),
  returnclass = "sf")


####################
## MAMMAL SPECIES ##
####################

sps_group <- "mammals"

# Read and combine IUCN shapefiles for mammals terrestrial + freshwater
sps_all <- bind_rows(
  st_read("data/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp",
          quiet = TRUE),
  st_read("~/NatPoKe/data/externaldata/MAMMALS_FRESHWATER/MAMMALS_FRESHWATER.shp",
          quiet = TRUE))

# keep threatened species only
sps_all <- sps_all %>%
  filter(category %in% c("CR", "EN", "VU"))

# intersect mammals .shp with iberian peninsula
sf_use_s2(FALSE)
sps_extent <- st_intersection(sps_all, extent_sf)

# get a species list
target_sps <- unique(sps_extent$sci_name)

# inspect species and REMOVE species manually
target_sps <- setdiff(target_sps, "Macaca sylvanus") # REMOVE Macaca sylvanus (non-native species)

# filter final spatial object
sps_extent <- sps_extent %>%
  filter(sci_name %in% target_sps)

# build taxonomy table
sps_names <- sps_extent %>%
  st_drop_geometry() %>%
  select(id_no, sci_name, order_, family, genus, legend) %>%
  distinct(sci_name, .keep_all = TRUE)

# add species occurrences (from GBIF with rgbif pkg

mammalSpecies <- unique(sps_names$sci_name)

message("Querying GBIF for ", length(mammalSpecies), " species...")
pb <- utils::txtProgressBar(min = 0, max = length(mammalSpecies), style = 3)
occ_total <- integer(length(mammalSpecies))
occ_2015  <- integer(length(mammalSpecies))

for (i in seq_along(mammalSpecies)) {
  
  # Total occurrences (with coordinates)
  occ_total[i] <- rgbif::occ_count(
    scientificName = mammalSpecies[i],
    hasCoordinate  = TRUE,
    basisOfRecord = "HUMAN_OBSERVATION;OBSERVATION;OCCURRENCE;MACHINE_OBSERVATION"
  )
  
  # Occurrences from 2015 onwards
  occ_2015[i] <- rgbif::occ_count(
    scientificName = mammalSpecies[i],
    hasCoordinate  = TRUE,
    basisOfRecord = "HUMAN_OBSERVATION;OBSERVATION;OCCURRENCE;MACHINE_OBSERVATION",
    year           = "2015,2024"
  )
  # add progress bar as this step might take a while
  utils::setTxtProgressBar(pb, i)
}
close(pb)

# format gbif output as a dataframe
occ_df <- data.frame(
  sci_name       = mammalSpecies,
  occ_count_all  = occ_total,
  occ_count_2015 = occ_2015,
  row.names      = NULL
)

# left join with taxonomy table
sps_names <- sps_names %>% 
  left_join( occ_df, by = "sci_name")

# write table into a csv file (taxonomy + occ number)
write.csv(sps_names,
          paste0("data/iberian_threatened_", sps_group, "_taxonomy.csv"),
          row.names = FALSE)

##################
## BIRD SPECIES ##
##################

# How to obtain a Bird species list from BirdLife ------------------------------
# - Go to BirdLife International website: https://datazone.birdlife.org/search
# - Choose the filters on the left hand side to choose the target area/criteria 
# - Click the "Download" button 
# - Follow the website instruction for download
# - Once on your computer, extract the downloaded .zip file and move it into your data folder
# ------------------------------------------------------------------------------

sps_group <- "birds"

# read csv species list
bird_sps <- read.csv("data/BIRDLIFE_spsList/species-filter-results.csv")

# filter for threatened species
birds_threatened <- bird_sps %>%
  filter(RL.Category %in% c("CR", "EN", "VU"))
invisible(gc(rm(bird_sps)))

# inspect species and REMOVE or ADD species manually
target_bird_sps <- unique(birds_threatened$Scientific.name)

# get taxonomic info
taxo_info <- taxize::classification(target_bird_sps, db = "gbif")
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

bird_sps_names <- left_join(birds_threatened %>% rename(sci_name = Scientific.name),
                       taxo_df, by = "sci_name") %>% 
  dplyr::select(SIS.ID, Common.name, sci_name, order_, family, genus)

# add species occurrences (from GBIF with rgbif pkg
birdSpecies <- unique(bird_sps_names$sci_name)

message("Querying GBIF for ", length(birdSpecies), " species...")
pb <- utils::txtProgressBar(min = 0, max = length(birdSpecies), style = 3)
occ_total <- integer(length(birdSpecies))
occ_2015  <- integer(length(birdSpecies))

for (i in seq_along(birdSpecies)) {
  
  # Total occurrences (with coordinates)
  occ_total[i] <- rgbif::occ_count(
    scientificName = birdSpecies[i],
    hasCoordinate  = TRUE,
    basisOfRecord = "HUMAN_OBSERVATION;OBSERVATION;OCCURRENCE;MACHINE_OBSERVATION"
  )
  
  # Occurrences from 2015 onwards
  occ_2015[i] <- rgbif::occ_count(
    scientificName = birdSpecies[i],
    hasCoordinate  = TRUE,
    basisOfRecord = "HUMAN_OBSERVATION;OBSERVATION;OCCURRENCE;MACHINE_OBSERVATION",
    year           = "2015,2024"
  )
  # add progress bar as this step might take a while
  utils::setTxtProgressBar(pb, i)
}
close(pb)

# format gbif output as a dataframe
occ_df_birds <- data.frame(
  sci_name       = birdSpecies,
  occ_count_all  = occ_total,
  occ_count_2015 = occ_2015,
  row.names      = NULL
)

# left join with taxonomy
bird_sps_names <- bird_sps_names %>% 
  left_join(occ_df, by = "sci_name")

# write table into a csv file (taxonomy + occ number)
write.csv(bird_sps_names,
          paste0("data/iberian_threatened_", sps_group, "_taxonomy.csv"),
          row.names = FALSE)
