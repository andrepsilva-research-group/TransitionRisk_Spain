## Mammal sps Occurrence Retrieval ##
## Inês Silva & Jorinde-M. Rieger ##
## 07 June 2026 ##

source("./src/libraries.R")
source("./src/customFunctions2.R")


##########
# STEP 1 # import mammals species list
##########

sps_names <- read_csv("data/iberian_threatened_mammals_taxonomy.csv")

##########
# STEP 2 # Retrieve species occurrence records
##########

# load GBIF credentials
keys <- yaml::read_yaml("./config/api_keys.yml")

# run get_taxa_occurrences() function for mammals
GBIF_data1 <- get_taxa_occurrences(
  iucn_shapefile = "data/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp", #"./data/externaldata/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp",
  taxa_filter = sps_names$sci_name,
  gbif_user = keys$gbif_user,
  gbif_pwd = keys$gbif_pwd,
  gbif_email = keys$gbif_email)


# galemis pyrenaicus and mustela lutreola are in a different shape file we need to re-run occ search twice
GBIF_data2 <- get_taxa_occurrences(
  iucn_shapefile = "~/NatPoKe/data/externaldata/MAMMALS_FRESHWATER/MAMMALS_FRESHWATER.shp",
  taxa_filter = sps_names$sci_name,
  gbif_user = keys$gbif_user,
  gbif_pwd = keys$gbif_pwd,
  gbif_email = keys$gbif_email)

GBIF_data <- dplyr::bind_rows(GBIF_data1, GBIF_data2) %>%
  distinct()
# nrows = 51184

##########
# STEP 3 # Remove domestic animals occurrences
##########

# domestic species usually contain these words in names
pattern <- paste(c("domesticus", "familiaris"), collapse = "|")

# filter the GBIF occurrences
GBIF_data_df <- GBIF_data %>%
  # remove unwanted species automatically
  dplyr::filter(
    !grepl(pattern, verbatimScientificName, ignore.case = TRUE)) %>%
  # recheck if species remaining have over 30 occ
  dplyr::filter(n() > 30) %>% 
  # get coordinates again
  mutate(
    decimalLongitude = st_coordinates(.)[,1],
    decimalLatitude = st_coordinates(.)[,2]) %>%
  # drop geometry
  st_set_geometry(NULL)

##########
# STEP 4 # Export screened occurrence dataset to .csv file 
##########

invisible(gc(rm(GBIF_data)))

# occurrences file
write.csv(GBIF_data_df[, c("species", "decimalLatitude", "decimalLongitude", "year")],
          paste0("./data/GBIF_occurrences_", sps_group, "_", Sys.Date(), ".csv"))

message("✅ Screened GBIF  mammals data exported to 'GBIF_occurrences_mammals.csv'.")

message("Following species had range polygons but not enough records:\n",
        paste0(setdiff(unique(sps_names$sci_name),
                       unique(GBIF_data_df$species)), collapse = "\n")
)