## Minning Traits for Mammals Sps (Bank of Spain)
## Inês Silva & Jorinde-M. Rieger
## 31 July 2025


# !! NAVIGATION WARNING !! 
## This script is adapted from the SRIT-Database repo mammalTraits.R script to
## search only for specific species, add galemys pyrenaicus to the terrestrial
## mammals .shp file and use the genus means from homeRange package to fill in
## missing values for species with no values

# vector with specific species
bankSpain_sps <- c("Crocidura canariensis", "Eliomys quercinus", "Galemys pyrenaicus", 
                   "Lepus castroviejoi", "Lynx pardinus", "Myotis capaccinii", "Nyctalus lasiopterus", "Nyctalus azoreum",
                   "Pipistrellus maderensis", "Rhinolophus mehelyi", "Oryctolagus cuniculus",
                   "Plecotus teneriffae", "Miniopterus schreibersii")


if (length(list.files(pattern = "mammalTraits_.*\\.csv$")) != 0) {
  mammalTraits_data <- read_csv(list.files(pattern = "mammalTraits_.*\\.csv$"))
} else {
  
  ### Printing message
  cat(bold(red("No file detected:"),green("creating mammal trait database...")))
  
  ## import IUCN terrestrial mammals .shp
  IUCN_mammals <- sf::st_read("C:/Users/maria/OneDrive - Universidade de Lisboa/ANDRE/SRIT_ANDRE/external_data/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
  
  # import galemys pyrenaicus .shp
  galemys <- sf::st_read("C:/Users/maria/OneDrive - Universidade de Lisboa/ANDRE/galemis_pirenaicos/data_0.shp")
  # galemys is a freshwater mammals and thus is not included in the terrestrial mammals shapefile
  
  ## processing galemys .shp to match IUCN_mammals ##
  
  colnames(galemys) <- tolower(colnames(galemys))
  # get columns in IUCN **but not** in gamelis
  missing_cols <- setdiff(names(IUCN_mammals), names(galemys))
  # add missing columns to galemys (filled with NA)
  galemys[missing_cols] <- NA
  # reorder galemys columns to match IUCN_mammals
  galemys <- galemys[, names(IUCN_mammals)]
  # combine both
  combined_data <- rbind(IUCN_mammals, galemys)
  
  #colnames(IUCN_mammals) == colnames(combined_data)
  #combined_data[combined_data$sci_name == "Galemys pyrenaicus",]
  rm(IUCN_mammals)
  
  ## END OF processing galemys .shp to match IUCN_mammals ##
  invisible(gc())
  
  sf_use_s2(FALSE)
  # filter combined data for target species
  filtered_traits <- combined_data %>%
    distinct(sci_name, .keep_all = TRUE) %>%
    dplyr::filter(sci_name %in% bankSpain_sps) %>%
    dplyr::select(sci_name, family, order_) %>%
    st_drop_geometry()
  invisible(gc())
  
  # import the biome data
  biomes <- read_sf("C:/Users/maria/OneDrive - Universidade de Lisboa/ANDRE/SRIT_ANDRE/external_data/Ecoregions2017/Ecoregions2017.shp") %>% #"ecoregions/Ecoregions2017.shp"
    st_make_valid() %>%
    naniar::replace_with_na(list(BIOME_NAME = c("N/A"))) %>%
    dplyr::select('BIOME_NAME') %>%
    dplyr::distinct() %>%
    drop_na()
  invisible(gc())
  
  # import the continent data (https://figshare.com/articles/dataset/Continent_Polygons/12555170)
  continent <- read_sf("C:/Users/maria/OneDrive - Universidade de Lisboa/ANDRE/SRIT_ANDRE/external_data/continent/Continents.shp")[,2] #"continent/Continents.shp"
  
  # intersect ecoregions and continents shapefiles
  biomes <- st_intersection(st_make_valid(biomes),st_make_valid(continent))
  invisible(gc())
  
  # intersect IUCN and ecoregions and continents shapefiles
  biomes <- st_intersection(st_make_valid(biomes),st_make_valid(combined_data))
  invisible(gc())
  
  # subset of the necessary variables (biome and sps scientific names) 
  sps_biome <- biomes %>%
    dplyr::select('BIOME_NAME', 'sci_name', "CONTINENT") %>% 
    dplyr::filter(sci_name %in% bankSpain_sps) %>%
    st_drop_geometry() %>%
    distinct()
  invisible(gc())
  
  # remove unnecessary objects 
  rm(continent,combined_data)
  invisible(gc())
  
  # import home range data
  HomeRangeData <- GetHomeRangeData()
  # HomeRangeData_pop <- HomeRangeData %>% 
  #   #dplyr::filter(HR_Level %in% "Population (mean)") %>% 
  #   dplyr::filter(Species %in% bankSpain_sps)
  
  # processing data (tacking character values from the isopleth column)
  HomeRangeData <- HomeRangeData[!(HomeRangeData$Isopleth_Size %in% c("Home range","Core area",
                                                                      "70 - 100","75-95","95-99",
                                                                      "60 - 90","65 - 90","50-95",
                                                                      ">95")),]
  
  # processing data (tacking NAs out)
  #HomeRange_sps <- HomeRange_sps[-(which(is.na(HomeRange_sps$Isopleth_Size))),]
  
  # processing data (make isopleth column)
  HomeRangeData$Isopleth_Size <- as.numeric(HomeRangeData$Isopleth_Size)
  
  # processing data (selecting home ranges with isopleth that contain 90% or more)
  HomeRangeData <- HomeRangeData[HomeRangeData$Isopleth_Size >= 90,]
  
  # average homeRanges per genus (because only 6 of our target species have a record for its homeRange size)
  HomeRange_genus <- HomeRangeData %>%
    mutate(Genus = word(Species, 1)) %>% 
    group_by(Genus) %>%
    summarise(
      Max_Home_Range_km2 = max(Home_Range_km2, na.rm = TRUE),
      Min_Home_Range_km2 = min(Home_Range_km2, na.rm = TRUE),
      Mean_HomeRange_km2 = mean(Home_Range_km2, na.rm = TRUE),
      .groups = "drop"
    )
  invisible(gc())
  
  ## IMPORT TRAIT DATASETS ##
  
  # import PHYLACINE1.2 data from: https://megapast2future.github.io/PHYLACINE_1.2/ 
  PHYLACINE1.2_traits <- read_csv("C:/Users/maria/OneDrive - Universidade de Lisboa/ANDRE/SRIT_ANDRE/external_data/PHYLACINE_1.2.1/Data/Traits/Trait_data.csv") 
  PHYLACINE1.2_traits$Binomial.1.2 <- str_replace_all(PHYLACINE1.2_traits$Binomial.1.2, "_", " ") #remove "_" between the generic name and the specific epithet and replace with a "space"
  PHYLACINE1.2_traits$Diet.Meat <- PHYLACINE1.2_traits$Diet.Vertebrate+PHYLACINE1.2_traits$Diet.Invertebrate
  
  # download PANTHERIA data
  #pantheria <- read_tsv("https://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  # naniar::replace_with_na_all(condition = ~.x == -999)
  
  # import SANTINI et al 2022: https://onlinelibrary.wiley.com/doi/10.1111/geb.13476 
  santini2022 <- read_excel("C:/Users/maria/OneDrive - Universidade de Lisboa/ANDRE/SRIT_ANDRE/external_data/geb13476-sup-0002-tables1.xls")  #reports error of RMSE_bvs instead of numeric (it's ok)
  
  # # import ANIMAL DIVERSITY WEB data (deprecated)
  # ADW_202404081457 <- read_excel("C:/Users/maria/OneDrive - Universidade de Lisboa/ANDRE/SRIT_ANDRE/external_data/ADW_202404081457.xls")
  # unique(ADW_202404081457$`Mating System`)
  
  # import COMBINE database
  combine <- read_csv("C:/Users/maria/OneDrive - Universidade de Lisboa/ANDRE/SRIT_ANDRE/external_data/ecy3344-sup-0001-datas1/COMBINE_archives/trait_data_imputed.csv")
  
  
  ## SELECT, MERGE & HARMONIZE TRAITS ##
  
  combined_traits_data <- filtered_traits %>%
    mutate(Genus = word(sci_name, 1)) %>% # add genus column
    left_join(combine, by = c('sci_name' = 'iucn2020_binomial')) %>%
    left_join(PHYLACINE1.2_traits, by = c("sci_name" = "Binomial.1.2")) %>%
    #left_join(pantheria, by = c("sci_name" = "scientificNameStd")) %>% 
    #left_join(ADW_202404081457, by = c('sci_name'= 'Species')) %>% 
    left_join(santini2022, by = c('sci_name'= 'Species')) %>% 
    left_join(HomeRange_genus, by = "Genus") %>%
    left_join(sps_biome, by = c('sci_name' = 'sci_name'))  %>%
    mutate(
      Combined_HomeRange_km2 = coalesce(Mean_HomeRange_km2, home_range_km2)) %>%
    dplyr::select(
      # taxonomic info from IUCN & combine
      sci_name, family.x, order_,
      # Diet info from phylacine 1.2
      Diet.Meat, Diet.Plant, trophic_level,
      # Body size from PHYLACINE 1.2
      Mass.g, max_longevity_d, weaning_mass_g,
      # age & reproduction info from combine
      age_first_reproduction_d, litter_size_n, litters_per_year_n,
      # mating info from ADW
      #'Mating System',
      # abundance info from santini et al 
      PredMd, up75,
      # habitat related info from ecoregions.shp, continents.shp, homeRange pckg and combine
      BIOME_NAME, CONTINENT, Mean_HomeRange_km2, Max_Home_Range_km2, Min_Home_Range_km2, home_range_km2, Combined_HomeRange_km2
    ) %>%
    distinct()
  # %>% drop_na()
  
  # csv file with trait data
  write.csv(combined_traits_data, paste0("mammalTraitsbankOfSpain_",Sys.Date(),".csv"), row.names = FALSE)
  
  # summary info
  complete_data_sps <- data.frame(unique(combined_traits_data$sci_name))
  names(complete_data_sps) <-"Species"
  print(paste("NºSpecies:",nrow(complete_data_sps)))
}