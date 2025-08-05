## TaxaOccurrence.R ##
## Jorinde-M. Rieger ##
## Description: Downloads taxa occurence for multiple species from GBIF Database ##
## Date: August 4th 2025 ##

# Define input species -----------------------------------------------------------------
# Use IUCN species names and ranges
# downloaded manually, later find a way to download automatically through R
sps <- sf::st_read("data/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp") 
species_group <- "mammals"

# Bird species BirdLife Iberian peninsula
# downloaded manually from https://www.birdlife.org/ subset of Iberian peninsula
#sps <- read.csv("data/trait_datasets/BirdLife_Iberian_peninsula.csv")
#species_group <- "birds"

# Subset according to IUCN Red list category: vulnerable (VU), endangered (EN), critically endangered (CR)
threatened_status <- c("VU", "EN", "CR") # Check column names to specify the category and specifications
sps_threatened <- sps %>%
  filter(category %in% threatened_status)

# Subset bird species
#sps_threatened <- sps%>%
#  filter(RL.Category %in% threatened_status)
#sps_names <- unique(sps_threatened$Scientific.name)
#length(unique(sps_threatened$Scientific.name))

# Overlay IUCN data with extent of Spain and Portugal, take list of intersection
extent <- "Spain and Portugal"
extent_name <- "Spain and Portugal"
extent_sf <- rnaturalearth::ne_countries(scale = "medium", country = c("Spain", "Portugal"), returnclass = "sf")
sf::sf_use_s2(FALSE) # disable s2
sps_extent <- sf::st_intersection(sps_threatened, extent_sf)
sf::sf_use_s2(TRUE) # anable s2
length(unique(sps_extent$sci_name))
sps_names <- unique(sps_extent$sci_name)

# Inspect species and REMOVE or ADD species manually
sps_names <- sps_names[sps_names != "Macaca sylvanus"] # REMOVE Macaca sylvanus (non-native species)
sps_names <- c(sps_names, "Galemys pyrenaicus") # ADD Galemys pyrenaicus

# Extract occurrences available in GBIF (e.g. mammals)
gbif_taxon_keys <-
  as.data.frame(sps_names) %>%
  pull("sps_names") %>% #use the sps names from the file
  name_backbone_checklist() %>% #match to backbone
  filter(!matchType == "NONE") %>% #get matched names
  pull(usageKey) #get the GBIF taxon keys

# to download datasets from gbif credentials are necessary.
# Register at https://www.gbif.org/user/profile

test <- occ_download( # creates key
  pred_in("taxonKey", gbif_taxon_keys),
  format = "SIMPLE_CSV",
  user = "USERNAME", # ADD USERNAME HERE "jorinde_rgr"
  pwd = "PASSWORD", # ADD PASSWORD HERE "NatPoKe2025"
  email = "EMAIL") # ADD EMAIL ASSOCIATE WITH ACCOUNT HERE "jorinde.rieger@su.se"

# check if download is finished
occ_download_wait('0021647-250717081556266') # ADD download key from "test"

# retrieve the download from GBIF to local computer
d <- occ_download_get(
  key = '0021647-250717081556266', # Exchange with download key "test"
  path = "data/trait_datasets",
  overwrite = TRUE
)

# import download to current session
gbif_data <- occ_download_import(d)

GBIF_sps <- gbif_data %>%
  # remove occurrences without coordinates
  drop_na(c(decimalLatitude, decimalLongitude)) %>% 
  group_by(species) %>%
  # filter individuals with more than 30 occurrences
  dplyr::filter(n() > 30)

# Write species occurences, with the subselection of variables
write.csv(
  GBIF_sps[, c("species", "decimalLatitude", "decimalLongitude", "year")],
  file = paste0("data/trait_datasets/GBIF_",species_group, "_30+occurrences_", gsub(" ", "_", extent),".csv"),
  row.names = FALSE
)

# View the target species and add them to SDMRun.R
targetSpecies <- unique(GBIF_sps$species)

invisible(gc())
