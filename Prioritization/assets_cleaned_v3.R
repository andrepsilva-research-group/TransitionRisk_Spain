#assets_cleaned_v3
setwd("path/to/Spain_TranRisk/Prioritization")
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(readxl)
library(patchwork)
library(ggnewscale)
library(arrow)
library(geodata)
terraOptions(threads = parallel::detectCores() - 1)
# ==============================================================================
# FILE PATHS / SET UP
EXCEL_PATH <- "Spain_power_mining_28.05.26.xlsx"
industrial_main   <- "input_layers_assets/Peninsula_IslasBaleares_Ceuta_Melilla.gpkg"
industrial_canary <- "input_layers_assets/Canarias.gpkg"
AGRI_RAST_PATH   <- "data/agri/eurcrops_spain_1km.tif"
AGRI_FRAC_PATH   <- "data/agri/eurcrops_spain_1km_frac.tif"
PARQUET_PATH     <- "input_layers_assets/es_stack.parquet"
# Coverage Threshold per PU
FRAC_THRESHOLD <- 0.5
# ==============================================================================
# Country Boundaries ----
pt <- ne_countries(scale = "medium", country = "Portugal", returnclass = "sf")
es <- ne_countries(scale = "medium", country = "Spain",    returnclass = "sf")
countries   <- rbind(pt, es) |> st_transform(3035)
countries_v <- vect(countries)
spain       <- countries |> filter(admin == "Spain")
spain_v     <- vect(spain)
# ==============================================================================
# Load Pre-Built Layers from Cache ----
PA_rast     <- rast("cache/PA_rast.tif")
area_raster <- rast("cache/area_raster.tif")
PA_spain   <- mask(PA_rast,     spain_v)
area_spain <- mask(area_raster, spain_v)
# ==============================================================================
# Load Scenario Solutions ----
sol_ssp1 <- rast("solutions/solution_ssp1.tif")
sol_ssp2 <- rast("solutions/solution_ssp2_regional.tif")
sol_ssp5 <- rast("solutions/solution_ssp5_development.tif")
PA_rast <- resample(PA_rast, sol_ssp1, method = "near")
# ==============================================================================
# Set up Asset Analysis ----
REQUIRED_COLS <- c(
  "asset_id", "asset_name", "ISIC_section", "Asset_weight",
  "Industry_production", "Industry_GVA", "Industry_employment",
  "Asset_production_share", "Asset_GVA_share", "Asset_employment_share"
)
coerce_sector <- function(sf_obj, isic_label) {
  result <- sf_obj |>
    mutate(across(any_of(c("Asset_production_share","Asset_GVA_share",
                           "Asset_employment_share","Asset_weight",
                           "Industry_production","Industry_GVA",
                           "Industry_employment")), as.numeric)) |>
    mutate(ISIC_section = isic_label) |>
    select(all_of(REQUIRED_COLS)) |>
    st_transform(4326)
  st_geometry(result) <- "geometry"
  result
}
empty_sector_sf <- function() {
  st_sf(
    asset_id = character(), asset_name = character(),
    ISIC_section = character(), Asset_weight = numeric(),
    Industry_production = numeric(), Industry_GVA = numeric(),
    Industry_employment = numeric(), Asset_production_share = numeric(),
    Asset_GVA_share = numeric(), Asset_employment_share = numeric(),
    geometry = st_sfc(crs = 4326)
  )
}
# ==============================================================================
#SECTOR 1 — POWER ----
power_raw <- read_excel(EXCEL_PATH, sheet = "Power (operating)")
power_sf <- power_raw |>
  mutate(latitude  = as.numeric(Latitude),
         longitude = as.numeric(Longitude)) |>
  filter(!is.na(latitude), !is.na(longitude)) |>
  mutate(asset_id   = `GEM unit/phase ID`,
         asset_name = `Plant / Project name`) |>
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) |>
  coerce_sector("Electricity, gas, steam and air conditioning supply")
# ==============================================================================
#SECTOR 2 — MINING ----
mining_raw <- read_excel(EXCEL_PATH, sheet = "Mining (PRTR23- active)")
mining_sf <- mining_raw |>
  mutate(
    latitude  = as.numeric(Latitude...3),
    longitude = as.numeric(Longitude...4)
  ) |>
  filter(!is.na(latitude), !is.na(longitude)) |>
  mutate(asset_id   = as.character(PRTRcode),
         asset_name = FacilityName) |>
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) |>
  coerce_sector("Mining and quarrying")
# ==============================================================================
#SECTOR 3 — INDUSTRIAL ZONES ----
both_siose_exist <- file.exists(industrial_main) && file.exists(industrial_canary)
if (both_siose_exist) {

  siose_main <- st_read(industrial_main,
                        layer = st_layers(industrial_main)$name[1],
                        quiet = TRUE) |> st_make_valid()
  siose_can  <- st_read(industrial_canary,
                        layer = st_layers(industrial_canary)$name[1],
                        quiet = TRUE) |> st_make_valid()

  cents  <- st_centroid(st_transform(siose_main, 4326))
  coords <- st_coordinates(cents)
  siose_main <- siose_main[!(coords[,"Y"] < 35.95 & coords[,"X"] < -1.6), ]

  if (st_crs(siose_main) != st_crs(siose_can))
    siose_can <- st_transform(siose_can, st_crs(siose_main))

  siose_raw <- bind_rows(siose_main, siose_can)

  # SIOSE AR 2020, as sourced for this project, contains only industrial-use
  # zones -- no land-use filter is applied.
  siose_mfg <- siose_raw

  sf::sf_use_s2(FALSE)
  siose_mfg_ea <- st_transform(siose_mfg, 3035) |>
    st_make_valid() |>
    filter(lengths(st_dimension(st_geometry(siose_mfg))) >= 0) |>
    st_collection_extract("POLYGON")
  siose_mfg_ea$area_m2 <- as.numeric(st_area(siose_mfg_ea))
  siose_mfg_ea$area_ha <- siose_mfg_ea$area_m2 / 1e4
  siose_mfg_ea <- siose_mfg_ea |>
    filter(area_ha > 0) |>
    mutate(.mfg_row_id = as.character(row_number()))
  sf::sf_use_s2(TRUE)

  total_mfg_area_ha <- sum(siose_mfg_ea$area_ha, na.rm = TRUE)

  mfg_sf <- siose_mfg_ea |>
    st_centroid() |> st_transform(4326) |>
    mutate(asset_id = .mfg_row_id, asset_name = paste0("IndustrialZone_",.mfg_row_id),
           Asset_weight = area_ha / total_mfg_area_ha,
           Industry_production = NA_real_, Industry_GVA = NA_real_,
           Industry_employment = NA_real_, Asset_production_share = NA_real_,
           Asset_GVA_share = NA_real_, Asset_employment_share = NA_real_) |>
    coerce_sector("Industrial")

  mfg_poly_4326 <- siose_mfg_ea |>
    st_transform(4326) |>
    select(.mfg_row_id, area_ha,
           any_of(c("cpro","municipio","IGRP_nombre","IGRP_idpob","IGRP_codine")))

} else {
  mfg_sf        <- empty_sector_sf()
  siose_mfg_ea  <- st_sf(.mfg_row_id = character(), area_ha = numeric(),
                         geometry = st_sfc(crs = 3035))
  mfg_poly_4326 <- st_sf(.mfg_row_id = character(), area_ha = numeric(),
                         geometry = st_sfc(crs = 4326))
}
# ==============================================================================
#SECTOR 4 — AGRICULTURE ----
# Builds the cropland rasters 
PROVINCE_CROP_PATH <- "solutions/agriculture_province_crop_summary.csv"

remove_accents <- function(x) {
  x <- gsub("[áàâä]", "a", x); x <- gsub("[ÁÀÂÄ]", "A", x)
  x <- gsub("[éèêë]", "e", x); x <- gsub("[ÉÈÊË]", "E", x)
  x <- gsub("[íìîï]", "i", x); x <- gsub("[ÍÌÎÏ]", "I", x)
  x <- gsub("[óòôö]", "o", x); x <- gsub("[ÓÒÔÖ]", "O", x)
  x <- gsub("[úùûü]", "u", x); x <- gsub("[ÚÙÛÜ]", "U", x)
  x <- gsub("ñ",      "n", x); x <- gsub("Ñ",      "N", x)
  x
}

if (!file.exists(AGRI_RAST_PATH) || !file.exists(AGRI_FRAC_PATH) || !file.exists(PROVINCE_CROP_PATH)) {

  stopifnot(packageVersion("terra") >= "1.7.39")
  dir.create("data/agri", recursive = TRUE, showWarnings = FALSE)
  dir.create("solutions", recursive = TRUE, showWarnings = FALSE)

  rast_stack        <- c(PA_rast, sol_ssp1, sol_ssp2, sol_ssp5)
  names(rast_stack) <- c("pa_baseline","pa_ssp1","pa_ssp2","pa_ssp5")
  template_r        <- sol_ssp1; values(template_r) <- NA
  spain_bbox        <- st_bbox(spain)

  prov_gadm <- geodata::gadm("ESP", level = 2, path = "GADM") |>
    st_as_sf() |>
    st_transform(3035) |>
    mutate(
      prov_id       = row_number(),
      province_name = remove_accents(NAME_2),
      ccaa_name     = remove_accents(NAME_1),
      province_code = CC_2
    )
  prov_rast        <- rasterize(vect(prov_gadm), template_r, field = "prov_id")
  names(prov_rast) <- "prov_id"
  prov_meta <- prov_gadm |>
    st_drop_geometry() |>
    select(prov_id, province_code, province_name, ccaa_name)

  tbl <- open_dataset(PARQUET_PATH) |>
    filter(!is.na(c2023)) |>
    collect() |>
    select(cropfield, c2023, area_ha, geometry)

  es_map <- tryCatch(
    read_csv(
      "https://raw.githubusercontent.com/maja601/EuroCrops/main/csvs/country_mappings/es_2023.csv",
      show_col_types = FALSE
    ),
    error = function(e) NULL
  )

  # --- calculating overlap between agr parcel and PA by CENTROID, running by chunk due to large file size ----
  chunk_size      <- 100000
  n_chunks        <- ceiling(nrow(tbl) / chunk_size)
  result_list     <- vector("list", n_chunks)
  all_parcel_list <- vector("list", n_chunks)

  for (i in seq_len(n_chunks)) {

    idx    <- ((i-1)*chunk_size + 1):min(i*chunk_size, nrow(tbl))
    chunk  <- tbl[idx, ]
    sfc    <- st_as_sfc(chunk$geometry, crs = 3035)
    cents  <- st_centroid(sfc)
    coords <- st_coordinates(cents)

    keep <- coords[,"X"] >= spain_bbox["xmin"] &
      coords[,"X"] <= spain_bbox["xmax"] &
      coords[,"Y"] >= spain_bbox["ymin"] &
      coords[,"Y"] <= spain_bbox["ymax"]
    if (sum(keep) == 0) next

    chunk_keep  <- chunk[keep, ]
    coords_keep <- coords[keep, , drop = FALSE]

    pts       <- vect(coords_keep, crs = "EPSG:3035")
    ext_vals  <- terra::extract(rast_stack, pts)[, -1]
    prov_vals <- terra::extract(prov_rast,  pts)[, 2]

    all_parcel_list[[i]] <- data.frame(
      cropfield = chunk_keep$cropfield,
      crop_code = chunk_keep$c2023,
      area_km2  = chunk_keep$area_ha / 100,
      prov_id   = prov_vals,
      stringsAsFactors = FALSE
    )

    chunk_result <- bind_cols(
      chunk_keep |> select(cropfield, crop_code = c2023, area_km2 = area_ha) |>
        mutate(area_km2 = area_km2 / 100),
      as.data.frame(ext_vals)
    ) |>
      mutate(across(starts_with("pa_"),
                    ~ ifelse(is.na(.), 0L, as.integer(. == 1))))

    in_any_pa <- rowSums(
      chunk_result[, c("pa_baseline","pa_ssp1","pa_ssp2","pa_ssp5")]
    ) > 0
    if (sum(in_any_pa) > 0)
      result_list[[i]] <- chunk_result[in_any_pa, ]
  }

  pa_parcels  <- bind_rows(result_list)
  all_parcels <- bind_rows(all_parcel_list)
  rm(result_list, all_parcel_list); gc()

  # --- Fractional cropland raster for plot -----------------------
  CHUNK_FRAC    <- 500000
  n_chunks_frac <- ceiling(nrow(tbl) / CHUNK_FRAC)
  frac_accum    <- template_r; values(frac_accum) <- 0.0
  cell_area_m2  <- prod(res(template_r))

  for (i in seq_len(n_chunks_frac)) {
    idx <- ((i-1)*CHUNK_FRAC + 1):min(i*CHUNK_FRAC, nrow(tbl))
    sfc <- st_as_sfc(tbl$geometry[idx], crs = 3035)

    cents  <- st_centroid(sfc)
    coords <- st_coordinates(cents)
    keep   <- coords[,"X"] >= spain_bbox["xmin"] &
      coords[,"X"] <= spain_bbox["xmax"] &
      coords[,"Y"] >= spain_bbox["ymin"] &
      coords[,"Y"] <= spain_bbox["ymax"]
    if (sum(keep) == 0) next

    vv <- vect(sfc[keep])
    cf            <- rasterizeGeom(vv, template_r, fun = "area") / cell_area_m2
    cf[is.na(cf)] <- 0
    cf[cf > 1]    <- 1
    frac_accum    <- frac_accum + cf
    frac_accum[frac_accum > 1] <- 1
  }

  frac_accum <- mask(frac_accum, spain_v)
  writeRaster(frac_accum, AGRI_FRAC_PATH, overwrite = TRUE)

  binary_accum <- ifel(frac_accum >= FRAC_THRESHOLD, 1L, NA)
  writeRaster(binary_accum, AGRI_RAST_PATH, overwrite = TRUE)

  rm(tbl); gc()

  # --- Province x crop type summary -- the one agriculture CSV this script writes
  province_crop_summary <- all_parcels |>
    left_join(prov_meta, by = "prov_id") |>
    filter(!is.na(province_code)) |>
    left_join(
      pa_parcels |> select(cropfield, pa_baseline, pa_ssp1, pa_ssp2, pa_ssp5),
      by = "cropfield"
    ) |>
    mutate(across(starts_with("pa_"), ~ replace(., is.na(.), 0L))) |>
    mutate(original_code = as.numeric(crop_code))

  if (!is.null(es_map)) {
    province_crop_summary <- left_join(province_crop_summary,
                                       es_map |> select(original_code, translated_name, HCAT3_name),
                                       by = "original_code")
  } else {
    province_crop_summary <- mutate(province_crop_summary,
                                    translated_name = NA_character_, HCAT3_name = NA_character_)
  }

  province_crop_summary <- province_crop_summary |>
    group_by(province_code, province_name, ccaa_name,
             crop_code, translated_name, HCAT3_name) |>
    summarise(
      n_parcels    = n(),
      total_km2    = round(sum(area_km2,                  na.rm = TRUE), 3),
      km2_baseline = round(sum(area_km2[pa_baseline==1L], na.rm = TRUE), 3),
      km2_ssp1     = round(sum(area_km2[pa_ssp1==1L],     na.rm = TRUE), 3),
      km2_ssp2     = round(sum(area_km2[pa_ssp2==1L],     na.rm = TRUE), 3),
      km2_ssp5     = round(sum(area_km2[pa_ssp5==1L],     na.rm = TRUE), 3),
      .groups = "drop"
    ) |>
    arrange(province_code, desc(total_km2))

  write_csv(province_crop_summary, PROVINCE_CROP_PATH)
}

agri_rast_raw   <- rast(AGRI_RAST_PATH)
agri_rast_raw   <- resample(agri_rast_raw, sol_ssp1, method = "near")
agri_rast_spain <- mask(agri_rast_raw, spain_v)

if (file.exists(AGRI_FRAC_PATH)) {
  agri_frac_raw   <- rast(AGRI_FRAC_PATH)
  agri_frac_raw   <- resample(agri_frac_raw, sol_ssp1, method = "bilinear")
  agri_frac_spain <- mask(agri_frac_raw, spain_v)
} else {
  agri_frac_spain <- agri_rast_spain
}
# ==============================================================================
#Combine Sectors ----
assets_combined <- bind_rows(
  Filter(function(x) nrow(x) > 0, list(power_sf, mining_sf, mfg_sf))
) |>
  mutate(ISIC_section = factor(ISIC_section, levels = c(
    "Agriculture, forestry and fishing",
    "Electricity, gas, steam and air conditioning supply",
    "Industrial", "Mining and quarrying"
  )))
assets_proj <- st_transform(assets_combined, crs = crs(sol_ssp1))
assets_vect <- vect(assets_proj)
# ==============================================================================
#Extract Raster Values at Asset Locations ----
vals_ssp1     <- terra::extract(sol_ssp1, assets_vect)[, 2]
vals_ssp2     <- terra::extract(sol_ssp2, assets_vect)[, 2]
vals_ssp5     <- terra::extract(sol_ssp5, assets_vect)[, 2]
vals_baseline <- terra::extract(PA_rast,  assets_vect)[, 2]
# ==============================================================================
#Assets per Scenario ----
assets_impact <- assets_combined |>
  st_drop_geometry() |>
  mutate(
    baseline_in_pa  = !is.na(vals_baseline) & vals_baseline == 1,
    ssp1_in_pa      = !is.na(vals_ssp1)     & vals_ssp1 == 1,
    ssp2_in_pa      = !is.na(vals_ssp2)     & vals_ssp2 == 1,
    ssp5_in_pa_surv = !is.na(vals_ssp5)     & vals_ssp5 == 1
  )
# ==============================================================================
#Industrial polygon x PA ----
if (nrow(mfg_poly_4326) > 0) {
  mfg_flags <- assets_impact |>
    filter(as.character(ISIC_section) == "Industrial") |>
    select(asset_id, baseline_in_pa, ssp1_in_pa, ssp2_in_pa, ssp5_in_pa_surv) |>
    rename(.mfg_row_id = asset_id)
  mfg_poly_flagged <- mfg_poly_4326 |> left_join(mfg_flags, by = ".mfg_row_id")
} else {
  mfg_poly_flagged <- mfg_poly_4326
}
# ==============================================================================
# Save CSV Power, Industrial, Mining, Scenario Summary Outputs ----

industrial_overlap <- mfg_poly_flagged |>
  st_drop_geometry() |>
  mutate(
    area_km2             = area_ha / 100,
    overlap_baseline_km2 = ifelse(baseline_in_pa,  area_km2, 0),
    overlap_ssp1_km2     = ifelse(ssp1_in_pa,      area_km2, 0),
    overlap_ssp2_km2     = ifelse(ssp2_in_pa,      area_km2, 0),
    overlap_ssp5_km2     = ifelse(ssp5_in_pa_surv, area_km2, 0)
  ) |>
  mutate(across(starts_with("overlap_"), ~replace_na(., 0))) |>
  rename(zone_id = .mfg_row_id) |>
  select(zone_id, any_of(c("cpro","municipio","IGRP_nombre","IGRP_idpob","IGRP_codine")),
         area_ha, area_km2, starts_with("overlap_"))

power_pa_count <- power_raw |>
  mutate(asset_id = `GEM unit/phase ID`) |>
  left_join(
    assets_impact |>
      filter(as.character(ISIC_section) == "Electricity, gas, steam and air conditioning supply") |>
      transmute(asset_id, pa_baseline = as.integer(baseline_in_pa),
                pa_ssp1 = as.integer(ssp1_in_pa), pa_ssp2 = as.integer(ssp2_in_pa),
                pa_ssp5 = as.integer(ssp5_in_pa_surv)),
    by = "asset_id"
  )

mining_pa_count <- mining_raw |>
  mutate(asset_id = as.character(PRTRcode)) |>
  left_join(
    assets_impact |>
      filter(as.character(ISIC_section) == "Mining and quarrying") |>
      transmute(asset_id, pa_baseline = as.integer(baseline_in_pa),
                pa_ssp1 = as.integer(ssp1_in_pa), pa_ssp2 = as.integer(ssp2_in_pa),
                pa_ssp5 = as.integer(ssp5_in_pa_surv)),
    by = "asset_id"
  )

n_total <- nrow(assets_impact)
count_by_scenario <- tibble(
  Scenario     = c("Baseline","SSP1","SSP2","SSP5"),
  Assets_in_PA = c(sum(assets_impact$baseline_in_pa), sum(assets_impact$ssp1_in_pa),
                   sum(assets_impact$ssp2_in_pa),      sum(assets_impact$ssp5_in_pa_surv)),
  Total_assets = n_total,
  Pct_of_total = round(Assets_in_PA / n_total * 100, 1)
)

write.csv(power_pa_count,     "solutions/power_assets_pa_count_update.csv",         row.names = FALSE)
write.csv(mining_pa_count,    "solutions/mining_assets_pa_count_update.csv",         row.names = FALSE)
write.csv(industrial_overlap, "solutions/industrial_per_zone_pa_overlap_update.csv", row.names = FALSE)
write.csv(count_by_scenario,  "solutions/assets_scenario_summary_1km.csv",           row.names = FALSE)
# ==============================================================================
# Plotting setup ----

section_colours <- c(
  "Agriculture, forestry and fishing"                   = "#E69F00",
  "Electricity, gas, steam and air conditioning supply" = "#0072B2",
  "Industrial"                                          = "#5B2C6F",
  "Mining and quarrying"                                = "#117A65"
)
section_label_short <- c(
  "Agriculture, forestry and fishing"                   = "Agriculture",
  "Electricity, gas, steam and air conditioning supply" = "Power",
  "Industrial"                                          = "Industry",
  "Mining and quarrying"                                = "Mining"
)
section_sizes <- c(
  "Agriculture, forestry and fishing"                   = 0.8,
  "Electricity, gas, steam and air conditioning supply" = 0.8,
  "Industrial"                                          = 0.8,
  "Mining and quarrying"                                = 3
)
present_sections <- levels(droplevels(assets_combined$ISIC_section))
if (!is.null(agri_frac_spain))
  present_sections <- c("Agriculture, forestry and fishing",
                        present_sections[present_sections !=
                                           "Agriculture, forestry and fishing"])
section_colours_plot <- section_colours[names(section_colours) %in% present_sections]
section_labels_plot  <- section_label_short[names(section_label_short) %in% present_sections]
section_sizes_plot   <- section_sizes[names(section_sizes) %in% present_sections]
MFG_POLY_FILL  <- "#E8620A"
PA_ALPHA       <- 0.65
AGRI_ALPHA     <- 0.85
AGRI_COLOURS <- c("#FFFDE7","#FFE082","#E69F00","#B35900","#7A3B00")
PA_ONLY_COLOUR <- "#57C84D"
OVERLAP_COLOUR <- "#CC4A1A"
darken_hex <- function(hex, factor=0.4) {
  v <- col2rgb(hex)/255
  rgb(pmax(v[1]*(1-factor),0), pmax(v[2]*(1-factor),0), pmax(v[3]*(1-factor),0))
}
spain_4326     <- st_transform(spain,     4326)
countries_4326 <- st_transform(countries, 4326)
agri_fill_scale <- function(order = 3) {
  scale_fill_gradientn(
    name     = "Cropland cover",
    colours  = AGRI_COLOURS,
    limits   = c(FRAC_THRESHOLD, 1),
    na.value = "transparent",
    labels   = scales::label_percent(accuracy = 1),
    guide    = guide_colorbar(
      order          = order,
      barwidth       = 8,
      barheight      = 0.8,
      direction      = "horizontal",
      title.position = "top"
    )
  )
}
agri_alpha_scale <- function() {
  scale_alpha_continuous(
    range  = c(0.15, AGRI_ALPHA),
    limits = c(FRAC_THRESHOLD, 1),
    guide  = "none"
  )
}
map_theme <- theme_void() +
  theme(legend.position    = "bottom",
        legend.title       = element_text(face="bold", size=24),
        legend.text        = element_text(size=22),
        legend.key.size    = unit(2.4,"lines"),
        legend.spacing.x   = unit(0.6,"cm"),
        plot.title         = element_text(face="bold", size=11),
        plot.background    = element_rect(fill="white", color=NA),
        panel.background   = element_rect(fill="white", color=NA))
common_coord <- coord_sf(xlim=c(-9.5,4.5), ylim=c(35.5,44), expand=FALSE)
make_mfg_poly <- function(in_pa_col) {
  if (nrow(mfg_poly_flagged)==0) return(mfg_poly_flagged)
  mfg_poly_flagged |> mutate(in_pa = .data[[in_pa_col]])
}
mfg_poly_ssp1     <- make_mfg_poly("ssp1_in_pa")
mfg_poly_ssp2     <- make_mfg_poly("ssp2_in_pa")
mfg_poly_ssp5     <- make_mfg_poly("ssp5_in_pa_surv")
mfg_poly_baseline <- make_mfg_poly("baseline_in_pa")
# ==============================================================================
# Overlap layer builder ----

make_overlap_sf <- function(scenario_rast, agri_r = agri_frac_spain) {
  if (is.null(agri_r)) return(NULL)

  scenario_rast <- mask(scenario_rast, spain_v)
  agri_r        <- mask(agri_r,        spain_v)

  agri_bin <- ifel(!is.na(agri_r) & agri_r >= FRAC_THRESHOLD, 1L, 0L)
  pa_bin   <- ifel(!is.na(scenario_rast) & scenario_rast == 1, 1L, 0L)

  agri_masked <- ifel(!is.na(agri_r) & agri_r >= FRAC_THRESHOLD, agri_r, NA)
  agri_agg    <- aggregate(agri_masked, fact = 5, fun = "mean", na.rm = TRUE)
  agri_sf     <- as.polygons(agri_agg, dissolve = FALSE) |>
    st_as_sf() |> st_transform(4326)
  names(agri_sf)[1] <- "crop_fraction"
  agri_sf <- agri_sf[!is.na(agri_sf$crop_fraction), ]

  pa_cat_r <- ifel(agri_bin == 0 & pa_bin == 1, 2L,
                   ifel(agri_bin == 1 & pa_bin == 1, 3L, NA))
  names(pa_cat_r) <- "category"

  pa_sf <- as.polygons(
    aggregate(pa_cat_r, fact = 3, fun = "modal", na.rm = TRUE), dissolve = TRUE
  ) |> st_as_sf() |> st_transform(4326)
  names(pa_sf)[1] <- "category"
  pa_sf$category  <- as.integer(pa_sf$category)
  pa_sf <- pa_sf[!is.na(pa_sf$category), ]

  list(agri_sf = agri_sf, pa_sf = pa_sf)
}
if (!is.null(agri_frac_spain)) {
  overlap_baseline <- make_overlap_sf(PA_rast,                 agri_frac_spain)
  overlap_ssp1     <- make_overlap_sf(sol_ssp1,                agri_frac_spain)
  overlap_ssp2     <- make_overlap_sf(sol_ssp2,                agri_frac_spain)
  overlap_ssp5     <- make_overlap_sf(ifel(sol_ssp5==1,1L,NA), agri_frac_spain)
} else {
  overlap_baseline <- overlap_ssp1 <- overlap_ssp2 <- overlap_ssp5 <- NULL
}
# ==============================================================================
# Per-scenario asset sf layers ----
make_asset_sf <- function(in_pa_flag) {
  assets_proj |>
    mutate(in_pa = in_pa_flag, isic_section = ISIC_section) |>
    st_transform(4326)
}
assets_baseline_sf <- make_asset_sf(!is.na(vals_baseline) & vals_baseline==1)
# ==============================================================================
# Final plotting ----
make_agri_facet_map <- function(overlap_layers, title, subtitle = NULL) {
  if (is.null(overlap_layers))
    return(ggplot() + labs(title=title, subtitle="Agriculture disabled") + map_theme)

  ggplot() +
    geom_sf(data=countries_4326, fill="#f5f5f0", color="grey70", linewidth=0.2) +
    geom_sf(data=spain_4326,     fill="white",   color=NA) +
    geom_sf(data=overlap_layers$agri_sf,
            aes(fill=crop_fraction, alpha=crop_fraction), color=NA) +
    agri_fill_scale() +
    agri_alpha_scale() +
    new_scale_fill() +
    geom_sf(data=overlap_layers$pa_sf |> filter(category==2),
            fill=PA_ONLY_COLOUR, color=NA, alpha=PA_ALPHA) +
    geom_sf(data=overlap_layers$pa_sf |> filter(category==3),
            fill=OVERLAP_COLOUR, color=NA, alpha=0.90) +
    geom_sf(data=overlap_layers$pa_sf |> filter(!is.na(category)),
            aes(fill=factor(category)), color=NA, alpha=0) +
    scale_fill_manual(
      name   = "Protected area",
      values = c("2"=PA_ONLY_COLOUR, "3"=OVERLAP_COLOUR),
      labels = c("2"="PA — no cropland", "3"="Cropland inside PA"),
      guide  = guide_legend(override.aes=list(alpha=1, size=8))
    ) +
    geom_sf(data=spain_4326, fill=NA, color="black", linewidth=0.45) +
    labs(title=title, subtitle=subtitle) +
    common_coord + map_theme +
    theme(plot.title    = element_text(face="bold", size=14),
          plot.subtitle = element_text(size=11))
}

pa_baseline_sf_plot <- PA_rast |> mask(spain_v) |>
  (\(x){ x[x!=1] <- NA; x })() |>
  as.polygons(dissolve=TRUE) |> st_as_sf() |> st_transform(4326)
pa_ssp1_sf_plot <- sol_ssp1 |> mask(spain_v) |>
  (\(x){ x[x!=1] <- NA; x })() |>
  as.polygons(dissolve=TRUE) |> st_as_sf() |> st_transform(4326)
pa_ssp2_sf_plot <- sol_ssp2 |> mask(spain_v) |>
  (\(x){ x[x!=1] <- NA; x })() |>
  as.polygons(dissolve=TRUE) |> st_as_sf() |> st_transform(4326)
pa_ssp5_sf_plot <- ifel(sol_ssp5==1, 1L, NA) |> mask(spain_v) |>
  (\(x){ x[x!=1] <- NA; x })() |>
  as.polygons(dissolve=TRUE) |> st_as_sf() |> st_transform(4326)
scenario_list_22 <- list(
  list(label="Current (Baseline)",                             ovl=overlap_baseline, vals=NULL,
       mfg=mfg_poly_baseline, pa_sf_plot=pa_baseline_sf_plot),
  list(label="Coordinated Nature Transition (SSP1 aligned)",   ovl=overlap_ssp1,     vals=vals_ssp1,
       mfg=mfg_poly_ssp1,     pa_sf_plot=pa_ssp1_sf_plot),
  list(label="Middle of the Road (SSP2 aligned)",               ovl=overlap_ssp2,     vals=vals_ssp2,
       mfg=mfg_poly_ssp2,     pa_sf_plot=pa_ssp2_sf_plot),
  list(label="Production Led (SSP5 aligned)",                   ovl=overlap_ssp5,     vals=vals_ssp5,
       mfg=mfg_poly_ssp5,     pa_sf_plot=pa_ssp5_sf_plot)
)
to_slug <- function(x) tolower(gsub("[^a-zA-Z0-9]+", "_", x))
make_point_sector_panel_22 <- function(sec, vals_vec, mfg_polys, pa_sf) {

  if (!is.null(vals_vec)) {
    mask_sec   <- as.character(assets_impact$ISIC_section) == sec
    in_pa      <- !is.na(vals_vec[mask_sec]) & vals_vec[mask_sec] == 1
    n_sec      <- sum(mask_sec)
    assets_sec <- assets_proj[mask_sec, ] |>
      mutate(in_pa=in_pa, isic_section=sec) |> st_transform(4326)
  } else {
    assets_sec <- assets_baseline_sf |> filter(as.character(isic_section)==sec)
    n_sec      <- nrow(assets_sec)
    in_pa      <- assets_sec$in_pa
  }
  mfg_layer <- if (sec=="Industrial") mfg_polys else NULL

  p <- ggplot() +
    geom_sf(data=countries_4326, fill="#f5f5f0", color="grey70", linewidth=0.2) +
    geom_sf(data=spain_4326,     fill="white",   color=NA) +
    geom_sf(data=pa_sf,          fill=PA_ONLY_COLOUR, color=NA, alpha=PA_ALPHA)

  if (!is.null(mfg_layer) && nrow(mfg_layer) > 0) {
    p <- p +
      geom_sf(data=mfg_layer |> filter(!in_pa), fill=MFG_POLY_FILL,
              color=darken_hex(MFG_POLY_FILL), linewidth=0.15) +
      geom_sf(data=mfg_layer |> filter(in_pa),  fill=MFG_POLY_FILL,
              color=darken_hex(MFG_POLY_FILL), linewidth=0.15, alpha=0.85)
  } else {
    p <- p +
      geom_sf(data=assets_sec |> filter(!in_pa),
              aes(color=isic_section, size=isic_section), alpha=0.8) +
      geom_sf(data=assets_sec |> filter(in_pa),
              aes(color=isic_section, size=isic_section), alpha=0.8)
  }

  p +
    scale_color_manual(
      name   = "ISIC section",
      values = section_colours_plot,
      labels = section_labels_plot,
      drop   = FALSE,
      guide  = guide_legend(override.aes=list(size=8))
    ) +
    scale_size_manual(values=section_sizes_plot, guide="none") +
    geom_sf(data=spain_4326, fill=NA, color="black", linewidth=0.45) +
    labs(title=NULL, subtitle=NULL) +
    common_coord + map_theme
}
for (sec in present_sections) {

  sec_panels <- map(scenario_list_22, function(sc) {
    if (sec == "Agriculture, forestry and fishing") {
      p <- make_agri_facet_map(sc$ovl, title=NULL, subtitle=NULL)
    } else {
      p <- make_point_sector_panel_22(sec, sc$vals, sc$mfg, sc$pa_sf_plot)
    }
    p + labs(title=sc$label) +
      theme(plot.title = element_text(face="bold", size=13))
  })

  sector_plot <- wrap_plots(sec_panels, nrow=1) +
    plot_layout(guides="collect") +
    plot_annotation(
      title = section_label_short[sec],
      theme = theme(
        plot.title      = element_text(face="bold", size=18),
        plot.background = element_rect(fill="white", color=NA)
      )
    ) & theme(
      legend.position = "bottom",
      legend.title    = element_text(face="bold", size=10),
      legend.text     = element_text(size=9),
      legend.key.size = unit(1.2,"lines")
    )

  out_file <- file.path(
    "solutions",
    sprintf("assets_%s_scenarios_compared.png", to_slug(section_label_short[sec]))
  )

  ggsave(out_file, sector_plot, width=16, height=4.5, dpi=300, bg="white")
}
