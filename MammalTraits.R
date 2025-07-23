## MammalTraits.R ##
## Afonso Barrocal & Inês Silva, updated by Jorinde-M. Rieger ##
## July 23rd, 2025 ##

##########################################################################
## Things to change: ##
##   1  See what are the limiting factors in the database;
##   2. See mass from HerbiTrait and CarniTrait (do they have new data?);
##   3. Make more comments;
##   4. (...)
##########################################################################

library(easypackages)
easypackages::packages(
  "tidyverse",
  "readr",
  "tidyr",
  "dplyr",
  "pbapply",
  "gridExtra",
  "readxl",
  "ncdf4",
  "sf",
  "terra",
  "DescTools",
  "data.table",
  "ggplot2",
  "tibble",
  "stringr",
  "rnaturalearth",
  "rgbif",
  "DT",
  "HomeRange",
  "naniar",
  "crayon",
  "gdalUtilities",
  "BIEN",
  "dismo",
  "taxize",
  prompt = FALSE)

### trait data is on github
# remotes::install_github("RS-eco/traitdata")
# remotes::install_github("SHoeks/HomeRange", subdir='pkg') # this package also needs to be downloaded from github
# devtools::install_github("RangeShifter/RangeShiftR-package", ref="main", force=TRUE)
# credentials::set_github_pat("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") # last token used for github connection)

# Select mammal species
sps_names <- c("Sus scrofa", "Vulpes vulpes", "Alces alces", "Canis latrans", "Lynx rufus", "Martes americana", "Taxidea taxus", "Ursus americanus", 
"Leontopithecus caissara", "Leopardus pardalis", "Nasua nasua", "Aepyceros melampus", "Colobus angolensis", "Daubentonia madagascariensis",
"Diceros bicornis", "Erythrocebus patas", "Gorilla beringei", "Gorilla gorilla", "Orycteropus afer", "Pan paniscus", "Pan troglodytes",
"Papio anubis", "Papio ursinus", "Cervus nippon", "Cuon alpinus", "Felis chaus", "Macaca fuscata", "Pongo abelii", "Pongo pygmaeus",
"Lynx lynx", "Ursus arctos", "Canis lupus", "Rangifer tarandus", "Puma concolor", "Bison bison", "Ursus arctos", "Panthera onca",
"Crocuta crocuta", "Mandrillus sphinx", "Panthera pardus", "Syncerus caffer", "Acinonyx jubatus", "Panthera leo", "Connochaetes taurinus",
"Loxodonta africana", "Acinonyx jubatus", "Panthera tigris")

if (length(list.files(pattern = "mammalTraits_.*\\.csv$")) != 0) {
  mammalTraits_data <- read_csv(list.files(pattern = "mammalTraits_.*\\.csv$"))
} else {
  
  ### Printing message
  cat(bold(red("No file detected:"),green("creating mammal trait database...")))
  
  ### EXTRACT HOME RANGE DATA ###
  
  # import home range data
  HomeRangeData <- GetHomeRangeData()
  # Add Genus column
HomeRangeData <- HomeRangeData %>%
  mutate(Genus = word(Species, 1))

# Filter to selected species
  mammal_sps <- sps_names
  # filter home range 
  HomeRange_sps <- HomeRangeData %>%
    dplyr::filter(Species %in% mammal_sps) %>%
    dplyr::filter(n() > 30)
  head(HomeRange_sps)

  # create vector with all the species
  target_sps2 <- unique(HomeRange_sps$Species)
  target_sps <- mammal_sps
  tmp <- unique(c(target_sps, target_sps2))
  out <- data.frame(target_sps = target_sps[match(tmp, target_sps)],
                    target_sps2 = target_sps2[match(tmp, target_sps2)])
  sapply(out, function(x)  sum(is.na(x)))
  extra_sps_HomeRange <- out[is.na(out$target_sps), 2]
  total_target_sps <- c(target_sps, extra_sps_HomeRange)

  ### MAX HOMERANGE PER GENUS ###
  
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

HomeRange_species <- HomeRangeData %>%
  group_by(Genus, Species) %>%
  summarise(
    Max_Home_Range_km2 = max(Home_Range_km2, na.rm = TRUE),
    Min_Home_Range_km2 = min(Home_Range_km2, na.rm = TRUE),
    Mean_HomeRange_km2 = mean(Home_Range_km2, na.rm = TRUE),
    .groups = "drop"
  )
# calculate the genus mean maximum and minimum Home Range dimension in km2
  HomeRange_meanGenus <- HomeRange_species %>%
  group_by(Genus) %>%
  summarise(
    Max_Home_Range_km2 = mean(Max_Home_Range_km2, na.rm = TRUE),
    Min_Home_Range_km2 = mean(Min_Home_Range_km2, na.rm = TRUE),
    Mean_HomeRange_km2 = mean(Mean_HomeRange_km2, na.rm = TRUE)
  )
  # Filter the Genus according to target species
  target_genera <- unique(word(mammal_sps, 1))
  HomeRange_meanGenus <- HomeRange_meanGenus %>%
  filter(Genus %in% target_genera)
print(HomeRange_meanGenus)


  ################################################################################
  
  ## Extract trait data and build a R dataframe with extracted data
  IUCN_mammals <- sf::st_read("C:/Users/jori5367/Documents/NatPoKe/data/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
  # import traits from the IUCN shapefile
  sf_use_s2(FALSE)
  sps_traits<- IUCN_mammals %>%
    distinct(sci_name, .keep_all = TRUE) %>%
    dplyr::filter(sci_name %in% total_target_sps) %>%
    dplyr::select(sci_name, family, order_) %>%
    st_drop_geometry()
  invisible(gc())
  
  # import PHYLACINE1.2 data from: https://megapast2future.github.io/PHYLACINE_1.2/ 
  PHYLACINE1.2_traits <- read_csv("trait_datasets/PHYLACINE_1.2.1/Data/Traits/Trait_data.csv") 
  PHYLACINE1.2_traits$Binomial.1.2 <- str_replace_all(PHYLACINE1.2_traits$Binomial.1.2, "_", " ") #remove "_" between the generic name and the specific epithet and replace with a "space"
  PHYLACINE1.2_traits$Diet.Meat <- PHYLACINE1.2_traits$Diet.Vertebrate+PHYLACINE1.2_traits$Diet.Invertebrate
  
  # download PANTHERIA data
  pantheria <- read_tsv("https://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
    naniar::replace_with_na_all(condition = ~.x == -999)
  
  # import SANTINI et al 2022: https://onlinelibrary.wiley.com/doi/10.1111/geb.13476 
  santini2022 <- read_excel("trait_datasets/geb13476-sup-0002-tables1.xls")  #reports error of RMSE_bvs instead of numeric (it's ok)
  
  # import ANIMAL DIVERSITY WEB data (deprecated)
  ADW_202404081457 <- read_excel("trait_datasets/ADW_202404081457.xls")
  unique(ADW_202404081457$`Mating System`)
  
  # import COMBINE database
  combine <- read_csv("trait_datasets/ecy3344-sup-0001-datas1/COMBINE_archives/trait_data_imputed.csv")
  
  # table with info on the n of species per dataset
  datasets_info <- tibble(
    dataset = c("PHYLACINE1.2_traits", "santini2022", "ADW_202404081457", "combine"),
    n_sps = c(length(unique(PHYLACINE1.2_traits$Binomial.1.2)), length(unique(santini2022$Species)),
              length(unique(combine$iucn2020_binomial)), length(unique(combine$iucn2020_binomial))))
  invisible(gc())
  
  # import the biome data
  biomes <- read_sf("C:/Users/jori5367/Documents/NatPoKe/data/Ecoregions2017/Ecoregions2017/Ecoregions2017.shp") %>% #"ecoregions/Ecoregions2017.shp"
    st_make_valid() %>%
    naniar::replace_with_na(list(BIOME_NAME = c("N/A"))) %>%
    select('BIOME_NAME') %>%
    distinct() %>%
    drop_na()
  invisible(gc())
  
  # import the continent data (https://figshare.com/articles/dataset/Continent_Polygons/12555170)
  continent <- read_sf("trait_datasets/continent-poly/Continents.shp")[,2] #"continent/Continents.shp"
  
  # intersect ecoregions and continents shapefiles
  biomes <- st_intersection(st_make_valid(biomes),st_make_valid(continent))
  invisible(gc())
  
  # intersect IUCN and ecoregions and continents shapefiles
  biomes <- st_intersection(st_make_valid(biomes),st_make_valid(IUCN_mammals))
  invisible(gc())
  
  # subset of the necessary variables (biome and sps scientific names) from the mammals_forests dataframe
  sps_biome <- biomes %>%
    dplyr::select('BIOME_NAME', 'sci_name', "CONTINENT") %>% 
    dplyr::filter(sci_name %in% total_target_sps) %>%
    st_drop_geometry() %>%
    distinct()
  invisible(gc())
  
  # remove unnecessary objects to increase internal memory
  rm(biomes,continent,IUCN_mammals)
  invisible(gc())
  
  ################################################################################
  ## FINAL DATASET ##
  
  # combining traits
  combined_traits_data <- sps_traits %>%
    mutate(Genus = word(sci_name, 1)) %>% # add genus column
    left_join(combine, by = c('sci_name' = 'iucn2020_binomial')) %>%
    left_join(PHYLACINE1.2_traits, by = c("sci_name" = "Binomial.1.2")) %>%
    #left_join(pantheria, by = c("sci_name" = "scientificNameStd")) %>% 
    #left_join(ADW_202404081457, by = c('sci_name'= 'Species')) %>% 
    left_join(santini2022, by = c('sci_name'= 'Species')) %>% 
    left_join(HomeRange_meanGenus, by = "Genus") %>%
    left_join(sps_biome, by = c('sci_name' = 'sci_name')) %>%
    dplyr::select('sci_name', 'family.x', 'order_', 'BIOME_NAME', 'CONTINENT',
                  'Diet.Meat', 'Diet.Plant', 'Mass.g', 'max_longevity_d', 'weaning_mass_g','age_first_reproduction_d',
                  'litter_size_n', 'litters_per_year_n', 
                 # 'Mating System',
                  'PredMd','up75',
                  'trophic_level','Mean_HomeRange_km2', 'Max_Home_Range_km2', 'Min_Home_Range_km2') %>%
    distinct() #%>% 
   # drop_na()
View(combined_traits_data)

  # csv file with trait data
  write.csv(combined_traits_data, paste0("mammalTraits+genusRange_",Sys.Date(),".csv"), row.names = FALSE)

  # summary infor
  complete_data_sps <-data.frame(unique(combined_traits_data$sci_name))
  names(complete_data_sps) <-"Species"
  print(paste("NºSpecies:",nrow(complete_data_sps)))
}
