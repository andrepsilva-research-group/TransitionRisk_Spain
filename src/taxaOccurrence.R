## taxaOccurrence.R ##
## Jorinde-M. Rieger & Inês Silva ##
## Description: Downloads taxa occurence for multiple species from GBIF Database ##
## Date: October 13th 2025 ##

source("./src/libraries.R")
source("./src/customFunctions2.R")

##########
# STEP 1 # Get a target species list
##########

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
sps <- sf::st_read("data/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")

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
#target_sps <- c(target_sps, "Galemys pyrenaicus")

# filter spatial object for target species only
sps_extent <- sps_extent %>% 
  dplyr::filter(sci_name %in% target_sps)

# save the result as a shape file
#st_write(sps_extent, "data/iberian_threatened_mammals.shp")

# save a .csv file with the taxonomic info
sps_names <- sps_extent %>% 
  dplyr::select(id_no, sci_name, order_, family, genus, legend) %>% 
  dplyr::filter(sci_name %in% target_sps) %>% 
  st_drop_geometry() %>% 
  unique()
write.csv(sps_names, paste0("data/iberian_threatened_", sps_group, "_taxonomy.csv"), row.names = FALSE)


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
sps <- read.csv("data/BIRDLIFE_spsList/species-filter-results.csv")

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

sps_names <- left_join(sps_threatened %>% rename(sci_name = Scientific.name),
                       taxo_df, by = "sci_name") %>% 
  dplyr::select(SIS.ID, Common.name, sci_name, order_, family, genus)

# save a .csv file with the taxonomic info
write.csv(sps_names, paste0("data/iberian_threatened_", sps_group, "_taxonomy.csv"), row.names = FALSE)
invisible(gc(rm(taxo_df, sps_threatened)))

##########
# STEP 2 # Retrieve occurrences
##########

keys <- yaml::read_yaml("./config/api_keys.yml")

# get_taxa_occurrences() is user defined function that downloads and filters
# species occurrence records from GBIF for a set of taxa defined in an IUCN
# spatial shapefile, returning only records that fall within the species’ range.
# Check out the code behind it in customFunctions2.R script

# GBIF_data <- get_taxa_occurrences(iucn_shapefile = "data/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp",
#                      taxa_filter = target_sps,
#                      gbif_user = keys$gbif_user, 
#                      gbif_pwd = keys$gbif_pwd,
#                      gbif_email = keys$gbif_email)

GBIF_data <- get_taxa_occurrences(iucn_shapefile = "data/species/BOTW_2025.gpkg",
                                  taxa_filter = target_sps[[1]],
                                  gbif_user = keys$gbif_user, 
                                  gbif_pwd = keys$gbif_pwd,
                                  gbif_email = keys$gbif_email)
invisible(gc())

# manually inspect files for weird stuff
species_name_check <- GBIF_data %>%
  # summarise number of records per species and verbatim name
  group_by(species, verbatimScientificName) %>%
  summarise(n_records = n(), .groups = "drop") %>%
  st_drop_geometry() %>% 
  arrange(species, desc(n_records))

# open table for manual check
View(species_name_check)
# STOP script review is done
readline(prompt = "Press [Enter] to continue with filtering once you've reviewed the species name matches...")

# filter the GBIF occurrences
GBIF_data_df <- GBIF_data %>%
  # remove unwanted subspecies, domestic animals, or problematic names
  dplyr::filter(!verbatimScientificName %in% c(
    # weird bat record
    #"Myotis capaccinii, Myotis daubentonii",
    # wrong plover record
    "Charadarius squatarola (Linnaeus, 1758)"
  )) %>%
  # get coordinates again
  mutate(
    decimalLongitude = st_coordinates(.)[,1],
    decimalLatitude = st_coordinates(.)[,2]) %>%
  # drop geometry
  st_set_geometry(NULL)

##########
# STEP 3 # save occurrence records as .csv file
##########

# write species occurrence records to .csv
write.csv(
  GBIF_data_df[, c("species", "decimalLatitude", "decimalLongitude", "year")],
  file = paste0("data/GBIF_occurrences_", sps_group, "_30+occurrences_",".csv"),
  row.names = FALSE)

## QUICK VISUAL CHECK OF OCCURRENCE RECORDS ##
world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot() +
  geom_sf(data = world, fill = "grey95", color = "black") +
  geom_point(
    data = GBIF_data_df,
    aes(x = decimalLongitude, y = decimalLatitude),
    color = "red", alpha = 0.6, size = 1) +
  # facet per species
  facet_wrap(~species) +
  theme_minimal()+
  labs(title = "GBIF Occurrences of Threatened Species in the Iberian Peninsula")

message("✅ GBIF occurrence data cleaned and exported to 'data folder")