#assets_animalprod
setwd("path/to/Spain_TranRisk/Prioritization")

library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(readxl)
library(patchwork)
library(ggnewscale)

terraOptions(threads = parallel::detectCores() - 1)

# ==============================================================================
# FILE PATHS
# ==============================================================================

EXCEL_PATH <- "Spain_power_mining_28.05.26.xlsx"

industrial_main   <- "input_layers_assets/Peninsula_IslasBaleares_Ceuta_Melilla.gpkg"
industrial_canary <- "input_layers_assets/Canarias.gpkg"

# NEW: MITECO PRTR intensive livestock (pigs / poultry, IED Annex I 6.6) inventory.
# Place the uploaded workbook here (or update the path) — sheet name expected:
# "Intensive_animal_ag".
ANIMAL_AG_PATH <- "input_layers_assets/Spain_Industrial_animal_agriculture_active only.xlsx"

AGRI_IMPACT_PATH <- "solutions/impact_agriculture.csv"
AGRI_CROP_PATH   <- "solutions/agriculture_crop_type_summary.csv"
AGRI_RAST_PATH   <- "data/agri/eurcrops_spain_1km.tif"
AGRI_FRAC_PATH   <- "data/agri/eurcrops_spain_1km_frac.tif"

SIOSE_INDUSTRIAL_CODES <- NULL
SIOSE_LAND_USE_COL     <- NULL

# Must match 00_prep_eurcrops_agri.R
FRAC_THRESHOLD <- 0.05

# ==============================================================================
# 1. Country Boundaries
# ==============================================================================

pt <- ne_countries(scale = "medium", country = "Portugal", returnclass = "sf")
es <- ne_countries(scale = "medium", country = "Spain",    returnclass = "sf")

countries   <- rbind(pt, es) |> st_transform(3035)
countries_v <- vect(countries)
spain       <- countries |> filter(admin == "Spain")
spain_v     <- vect(spain)

# ==============================================================================
# 2. Load Pre-Built Layers from Cache
# ==============================================================================

message("Loading cache layers...")

PA_rast              <- rast("cache/PA_rast.tif")
area_raster          <- rast("cache/area_raster.tif")
artificial_prop_ssp5 <- rast("cache/artificial_prop_ssp5.tif")

PA_spain   <- mask(PA_rast,     spain_v)
area_spain <- mask(area_raster, spain_v)

# ==============================================================================
# 4. Load Scenario Solutions
# ==============================================================================

message("Loading scenario solutions...")

sol_ssp1 <- rast("solutions/solution_ssp1.tif")
sol_ssp2 <- rast("solutions/solution_ssp2_regional.tif")
sol_ssp5 <- rast("solutions/solution_ssp5_development.tif")

PA_rast <- resample(PA_rast, sol_ssp1, method = "near")

message("Solutions loaded.")

# ==============================================================================
# HELPERS
# ==============================================================================

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
# 5. SECTOR 1 — POWER
# ==============================================================================

message("\n=== SECTOR 1: POWER [Global Energy Monitor] ===")

power_raw <- read_excel(EXCEL_PATH, sheet = "Power (operating)")

power_sf <- power_raw |>
  mutate(latitude  = as.numeric(Latitude),
         longitude = as.numeric(Longitude)) |>
  filter(!is.na(latitude), !is.na(longitude)) |>
  mutate(asset_id   = `GEM unit/phase ID`,
         asset_name = `Plant / Project name`) |>
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) |>
  coerce_sector("Electricity, gas, steam and air conditioning supply")

message(sprintf("  Plants loaded: %d", nrow(power_sf)))

# ==============================================================================
# 6. SECTOR 2 — MINING
# ==============================================================================

message("\n=== SECTOR 2: MINING [MITECO PRTR 2023] ===")

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

message(sprintf("  Facilities loaded: %d", nrow(mining_sf)))

# ==============================================================================
# 7. SECTOR 3 — INDUSTRIAL ZONES
# ==============================================================================

message("\n=== SECTOR 3: INDUSTRIAL ZONES [SIOSE AR 2020] ===")

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
  siose_mfg <- if (!is.null(SIOSE_INDUSTRIAL_CODES)) {
    siose_raw |> filter(.data[[SIOSE_LAND_USE_COL]] %in% SIOSE_INDUSTRIAL_CODES)
  } else siose_raw
  
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
  message(sprintf("  Zones: %d  (%.0f ha total)", nrow(siose_mfg_ea), total_mfg_area_ha))
  
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
  message("  SIOSE files not found — creating empty placeholder.")
  mfg_sf        <- empty_sector_sf()
  mfg_poly_4326 <- st_sf(.mfg_row_id = character(), area_ha = numeric(),
                         geometry = st_sfc(crs = 4326))
}

# ==============================================================================
# 8. SECTOR 4 — AGRICULTURE
# ==============================================================================

message("\n=== SECTOR 4: AGRICULTURE [EuroCrops 2023] ===")


missing_files <- c(AGRI_IMPACT_PATH, AGRI_CROP_PATH, AGRI_RAST_PATH)[
  !file.exists(c(AGRI_IMPACT_PATH, AGRI_CROP_PATH, AGRI_RAST_PATH))
]
if (length(missing_files) > 0)
  stop("Run 03_prep_eurocrops_agriculture.R first. Missing:\n",
       paste(" ", missing_files, collapse = "\n"))

agri_impact <- read.csv(AGRI_IMPACT_PATH)
agri_crops  <- read.csv(AGRI_CROP_PATH)

total_agri_km2    <- agri_impact$total_area_km2
agri_baseline_km2 <- agri_impact$area_baseline_km2
agri_ssp1_km2     <- agri_impact$area_ssp1_km2
agri_ssp2_km2     <- agri_impact$area_ssp2_km2
agri_ssp5_km2     <- agri_impact$area_ssp5_km2

message(sprintf("  Total cropland:  %.0f km²", total_agri_km2))
message(sprintf("  Baseline PA:     %.0f km²  (%.1f%%)",
                agri_baseline_km2, 100*agri_baseline_km2/total_agri_km2))
message(sprintf("  SSP1 PA:         %.0f km²  (%.1f%%)",
                agri_ssp1_km2, 100*agri_ssp1_km2/total_agri_km2))
message(sprintf("  SSP2 PA:         %.0f km²  (%.1f%%)",
                agri_ssp2_km2, 100*agri_ssp2_km2/total_agri_km2))
message(sprintf("  SSP5 PA:         %.0f km²  (%.1f%%)",
                agri_ssp5_km2, 100*agri_ssp5_km2/total_agri_km2))

agri_rast_raw   <- rast(AGRI_RAST_PATH)
agri_rast_raw   <- resample(agri_rast_raw, sol_ssp1, method = "near")
agri_rast_spain <- mask(agri_rast_raw, spain_v)
message("  Binary cropland raster loaded.")

if (file.exists(AGRI_FRAC_PATH)) {
  agri_frac_raw   <- rast(AGRI_FRAC_PATH)
  agri_frac_raw   <- resample(agri_frac_raw, sol_ssp1, method = "bilinear")
  agri_frac_spain <- mask(agri_frac_raw, spain_v)
  message("  Fractional cropland raster loaded — gradient viz enabled.")
} else {
  warning("  AGRI_FRAC_PATH not found — falling back to binary raster for viz.")
  agri_frac_spain <- agri_rast_spain
}


agr_sf <- empty_sector_sf()

# ==============================================================================
# 8b. SECTOR 5 — ANIMAL PRODUCTION [MITECO PRTR — Intensive Livestock]
# ==============================================================================

animal_raw <- read_excel(ANIMAL_AG_PATH, sheet = "Intensive_animal_ag")

animal_filtered <- animal_raw |>
  filter(installationStatus == "functional") |>
  mutate(latitude  = as.numeric(Latitude),
         longitude = as.numeric(Longitude)) |>
  filter(!is.na(latitude), !is.na(longitude)) |>

  filter(!Comunidad_autonoma %in% c("Ceuta","Melilla")) |>
  mutate(asset_id   = as.character(Codigo_PRTR),
         asset_name = installationName)

animal_sf <- animal_filtered |>
  mutate(Asset_weight = NA_real_) |>
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) |>
  coerce_sector("Animal production")

animal_meta <- animal_filtered |>
  transmute(asset_id, asset_name,
            Provincia, Comunidad_autonoma,
            activity_type = IEDMainActivityName,
            latitude, longitude)

message(sprintf("  Installations loaded: %d", nrow(animal_sf)))


# ==============================================================================
# 9. Combine Point/Polygon Sectors
# ==============================================================================

message("\n=== COMBINING POINT/POLYGON SECTORS ===")

assets_combined <- bind_rows(
  Filter(function(x) nrow(x) > 0, list(power_sf, mining_sf, mfg_sf, agr_sf, animal_sf))
) |>
  mutate(ISIC_section = factor(ISIC_section, levels = c(
    "Agriculture, forestry and fishing",
    "Animal production",
    "Electricity, gas, steam and air conditioning supply",
    "Industrial", "Mining and quarrying"
  )))


if (nrow(assets_combined) > 0) {
  coords_4326      <- st_coordinates(st_geometry(assets_combined))
  assets_combined$longitude <- coords_4326[, 1]
  assets_combined$latitude  <- coords_4326[, 2]
}

n_total     <- nrow(assets_combined)
assets_proj <- st_transform(assets_combined, crs = crs(sol_ssp1))
assets_vect <- vect(assets_proj)

message(sprintf("Total point/polygon assets: %d", n_total))
print(table(assets_combined$ISIC_section))

# ==============================================================================
# 10. Extract Raster Values at Asset Locations
# ==============================================================================

message("Extracting raster values...")
vals_ssp1     <- terra::extract(sol_ssp1, assets_vect)[, 2]
vals_ssp2     <- terra::extract(sol_ssp2, assets_vect)[, 2]
vals_ssp5     <- terra::extract(sol_ssp5, assets_vect)[, 2]
vals_baseline <- terra::extract(PA_rast,  assets_vect)[, 2]

# ==============================================================================
# 11. Flag Assets per Scenario
# ==============================================================================

assets_impact <- assets_combined |>
  st_drop_geometry() |>
  mutate(
    baseline_in_pa  = !is.na(vals_baseline) & vals_baseline == 1,
    ssp1_in_pa      = !is.na(vals_ssp1)     & vals_ssp1 == 1,
    ssp2_in_pa      = !is.na(vals_ssp2)     & vals_ssp2 == 1,
    ssp5_in_pa_surv = !is.na(vals_ssp5)     & vals_ssp5 == 1
  )

# ==============================================================================
# 11b. Manufacturing polygon flags
# ==============================================================================

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
# 11c. Animal Production — dedicated per-facility CSV
# ==============================================================================

animal_impact <- assets_impact |>
  filter(as.character(ISIC_section) == "Animal production") |>
  select(asset_id, asset_name, longitude, latitude,
         baseline_in_pa, ssp1_in_pa, ssp2_in_pa, ssp5_in_pa_surv) |>
  left_join(animal_meta |> select(-latitude, -longitude),
            by = c("asset_id", "asset_name"))

write.csv(animal_impact, "solutions/impact_animal_production.csv",
          row.names = FALSE)

message(sprintf(
  "  -> Per-facility animal production PA flags: solutions/impact_animal_production.csv (%d facilities)",
  nrow(animal_impact)))

# ==============================================================================
# 12. Point/Polygon Asset Summary
# ==============================================================================

count_by_scenario <- tibble(
  Scenario     = c("Baseline","SSP1","SSP2","SSP5"),
  Assets_in_PA = c(sum(assets_impact$baseline_in_pa), sum(assets_impact$ssp1_in_pa),
                   sum(assets_impact$ssp2_in_pa),      sum(assets_impact$ssp5_in_pa_surv)),
  Total_assets = n_total,
  Pct_of_total = round(Assets_in_PA / n_total * 100, 1)
)
message("\n=== POINT/POLYGON ASSET IMPACT ==="); print(count_by_scenario)

# ==============================================================================
# 13. Breakdown by ISIC Section
# ==============================================================================

section_breakdown <- assets_impact |>
  group_by(ISIC_section) |>
  summarise(
    n_total      = n(),
    n_baseline   = sum(baseline_in_pa), n_ssp1 = sum(ssp1_in_pa),
    n_ssp2       = sum(ssp2_in_pa),     n_ssp5 = sum(ssp5_in_pa_surv),
    pct_baseline = round(n_baseline/n_total*100,1),
    pct_ssp1     = round(n_ssp1/n_total*100,1),
    pct_ssp2     = round(n_ssp2/n_total*100,1),
    pct_ssp5     = round(n_ssp5/n_total*100,1),
    .groups = "drop"
  ) |> arrange(desc(n_total))

# ==============================================================================
# 14. Save CSV Outputs
# ==============================================================================

total_mfg_ha <- if (nrow(mfg_poly_flagged)>0) sum(mfg_poly_flagged$area_ha,na.rm=TRUE) else NA_real_

area_in_pa <- function(flag_col) {
  if (nrow(mfg_poly_flagged)==0) return(NA_real_)
  sum(mfg_poly_flagged$area_ha[mfg_poly_flagged[[flag_col]]], na.rm=TRUE)
}

# NOTE: Animal production is a point sector like Power/Mining, so it is
# folded into the same impact_point_sectors.csv output rather than getting
# its own file — consistent with how Power and Mining are already handled.
point_summary <- assets_impact |>
  filter(as.character(ISIC_section) %in%
           c("Electricity, gas, steam and air conditioning supply",
             "Mining and quarrying",
             "Animal production")) |>
  group_by(ISIC_section) |>
  summarise(n_total=n(), n_baseline=sum(baseline_in_pa), n_ssp1=sum(ssp1_in_pa),
            n_ssp2=sum(ssp2_in_pa), n_ssp5=sum(ssp5_in_pa_surv),
            pct_baseline=round(n_baseline/n_total*100,1),
            pct_ssp1=round(n_ssp1/n_total*100,1),
            pct_ssp2=round(n_ssp2/n_total*100,1),
            pct_ssp5=round(n_ssp5/n_total*100,1), .groups="drop") |>
  rename(Sector=ISIC_section)

industrial_summary <- tibble(
  Sector             = "Industrial zones (SIOSE AR 2020)",
  n_zones            = nrow(mfg_poly_flagged),
  total_area_km2     = round(total_mfg_ha/100, 2),
  area_baseline_km2  = round(area_in_pa("baseline_in_pa")/100, 2),
  area_ssp1_km2      = round(area_in_pa("ssp1_in_pa")/100, 2),
  area_ssp2_km2      = round(area_in_pa("ssp2_in_pa")/100, 2),
  area_ssp5_km2      = round(area_in_pa("ssp5_in_pa_surv")/100, 2),
  pct_baseline       = round(area_in_pa("baseline_in_pa")/total_mfg_ha*100,1),
  pct_ssp1           = round(area_in_pa("ssp1_in_pa")/total_mfg_ha*100,1),
  pct_ssp2           = round(area_in_pa("ssp2_in_pa")/total_mfg_ha*100,1),
  pct_ssp5           = round(area_in_pa("ssp5_in_pa_surv")/total_mfg_ha*100,1)
)

agriculture_summary <- read.csv(AGRI_IMPACT_PATH)

# --- NEW: fold agricultural parcel counts into section_breakdown -----------

if ("total_n_parcels" %in% names(agriculture_summary)) {
  agri_breakdown_row <- tibble(
    ISIC_section = "Agriculture, forestry and fishing",
    n_total      = agriculture_summary$total_n_parcels,
    n_baseline   = agriculture_summary$n_parcels_baseline,
    n_ssp1       = agriculture_summary$n_parcels_ssp1,
    n_ssp2       = agriculture_summary$n_parcels_ssp2,
    n_ssp5       = agriculture_summary$n_parcels_ssp5,
    pct_baseline = agriculture_summary$pct_n_baseline,
    pct_ssp1     = agriculture_summary$pct_n_ssp1,
    pct_ssp2     = agriculture_summary$pct_n_ssp2,
    pct_ssp5     = agriculture_summary$pct_n_ssp5
  )
  section_breakdown <- bind_rows(section_breakdown, agri_breakdown_row) |>
    arrange(desc(n_total))
} else {
  message("  NOTE: impact_agriculture.csv has no parcel-count columns yet — ",
          "re-run 00_prep_eurcrops_agri.R to regenerate it with n_parcels_* fields.")
}
# -----------------------------------------------------------------------------

section_label_short_csv <- c(
  "Agriculture, forestry and fishing"                   = "Agriculture",
  "Animal production"                                   = "Animal production",
  "Electricity, gas, steam and air conditioning supply" = "Power",
  "Industrial"                                          = "Industry",
  "Mining and quarrying"                                = "Mining"
)
asset_unit_csv <- c(
  "Agriculture, forestry and fishing"                   = "EuroCrops parcels",
  "Animal production"                                   = "PRTR installations",
  "Electricity, gas, steam and air conditioning supply" = "GEM power units",
  "Industrial"                                          = "SIOSE industrial zones",
  "Mining and quarrying"                                = "PRTR facilities"
)

message("\n=== PRIMARY IMPACT SUMMARY ===")
message("Point sectors:");        print(point_summary)
message("\nIndustrial zones:");    print(industrial_summary)
message("\nAgriculture (km²):");  print(agriculture_summary)

write.csv(point_summary,       "solutions/impact_point_sectors.csv",    row.names=FALSE)
write.csv(industrial_summary,  "solutions/impact_industrial_zones.csv", row.names=FALSE)
write.csv(assets_impact,       "solutions/assets_scenario_impact_1km.csv",   row.names=FALSE)
write.csv(count_by_scenario,   "solutions/assets_scenario_summary_1km.csv",  row.names=FALSE)
write.csv(section_breakdown,   "solutions/assets_section_breakdown_1km.csv", row.names=FALSE)

if (!is.null(agri_crops))
  write.csv(agri_crops, "solutions/agriculture_crop_type_summary.csv", row.names=FALSE)

exposure_by_sector <- section_breakdown |>
  mutate(Sector = section_label_short_csv[as.character(ISIC_section)],
         Unit   = asset_unit_csv[as.character(ISIC_section)]) |>
  select(Sector, Unit, n_total, n_baseline, n_ssp1, n_ssp2, n_ssp5) |>
  pivot_longer(c(n_baseline, n_ssp1, n_ssp2, n_ssp5),
               names_to = "Scenario", values_to = "n_exposed") |>
  mutate(Scenario    = recode(Scenario, n_baseline = "Baseline", n_ssp1 = "SSP1",
                              n_ssp2 = "SSP2", n_ssp5 = "SSP5"),
         pct_exposed = round(n_exposed / n_total * 100, 1)) |>
  select(Scenario, Sector, Unit, n_exposed, n_total, pct_exposed) |>
  arrange(Scenario, Sector)

message("\n=== EXPOSED ASSETS BY SCENARIO AND SECTOR ==="); print(exposure_by_sector, n = Inf)
write.csv(exposure_by_sector, "solutions/exposed_assets_by_scenario_sector.csv",
          row.names = FALSE)

message("All CSVs saved.")
message("  -> Exposed assets by scenario x sector: solutions/exposed_assets_by_scenario_sector.csv")
message("  -> Asset counts per scenario: solutions/assets_scenario_summary_1km.csv")
message("  -> Per-sector breakdown:      solutions/assets_section_breakdown_1km.csv")
message("  -> Per-asset PA flags with coordinates (all sectors, baseline/ssp1/ssp2/ssp5): solutions/assets_scenario_impact_1km.csv")
message("  -> Per-facility animal production PA flags (dedicated):        solutions/impact_animal_production.csv")

# ==============================================================================
# 15. Shared Plotting Setup
# ==============================================================================

section_colours <- c(
  "Agriculture, forestry and fishing"                   = "#E69F00",
  "Animal production"                                   = "#D55E00",
  "Electricity, gas, steam and air conditioning supply" = "#0072B2",
  "Industrial"                                          = "#5B2C6F",
  "Mining and quarrying"                                = "#117A65"
)

section_label_short <- c(
  "Agriculture, forestry and fishing"                   = "Agriculture",
  "Animal production"                                   = "Animal production",
  "Electricity, gas, steam and air conditioning supply" = "Power",
  "Industrial"                                          = "Industry",
  "Mining and quarrying"                                = "Mining"
)

# Point sizes — mining doubled relative to power (0.8 → 1.6)
section_sizes <- c(
  "Agriculture, forestry and fishing"                   = 0.8,
  "Animal production"                                   = 0.8,
  "Electricity, gas, steam and air conditioning supply" = 0.8,
  "Industrial"                                          = 0.8,
  "Mining and quarrying"                                = 3
)

present_sections <- levels(droplevels(assets_combined$ISIC_section))
present_sections <- c("Agriculture, forestry and fishing",
                      present_sections[present_sections !=
                                         "Agriculture, forestry and fishing"])

section_colours_plot <- section_colours[names(section_colours) %in% present_sections]
section_labels_plot  <- section_label_short[names(section_label_short) %in% present_sections]
section_sizes_plot   <- section_sizes[names(section_sizes) %in% present_sections]

MFG_POLY_FILL  <- "#E8620A"
MFG_POLY_ALPHA <- 1.0
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

# ==============================================================================
# 15b. Fallback PA sf and rast_to_sf helper
# ==============================================================================

PA_sf_plot <- PA_rast |> mask(spain_v) |>
  (\(x){ x[x!=1] <- NA; x })() |>
  as.polygons(dissolve=TRUE) |> st_as_sf() |> st_transform(4326)

rast_to_sf <- function(r, agg_factor=5) {
  poly <- as.polygons(aggregate(r, fact=agg_factor, fun="modal", na.rm=TRUE),
                      dissolve=TRUE) |>
    st_as_sf() |> st_transform(4326)
  names(poly)[1] <- "selected"; poly
}

sol_ssp1_sf <- rast_to_sf(sol_ssp1)
sol_ssp2_sf <- rast_to_sf(sol_ssp2)
sol_ssp5_sf <- rast_to_sf(sol_ssp5)

map_theme <- theme_void() +
  theme(legend.position    = "bottom",
        legend.title       = element_text(face="bold", size=10),
        legend.text        = element_text(size=9),
        legend.key.size    = unit(1.0,"lines"),
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

n_cols_facet <- min(2, length(present_sections))

# ==============================================================================
# 16. Overlap layer builder
# ==============================================================================

make_overlap_sf <- function(scenario_rast, agri_r = agri_frac_spain) {
  if (is.null(agri_r)) return(NULL)
  
  scenario_rast <- mask(scenario_rast, spain_v)
  agri_r        <- mask(agri_r,        spain_v)
  
  # NA guard: non-agricultural cells (NA in frac raster) → 0, not NA
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

message("Building overlap layers...")
overlap_baseline <- make_overlap_sf(PA_rast,                 agri_frac_spain)
overlap_ssp1     <- make_overlap_sf(sol_ssp1,                agri_frac_spain)
overlap_ssp2     <- make_overlap_sf(sol_ssp2,                agri_frac_spain)
overlap_ssp5     <- make_overlap_sf(ifel(sol_ssp5==1,1L,NA), agri_frac_spain)
message("  Done.")

# ==============================================================================
# 17. Per-scenario asset sf layers
# ==============================================================================

make_asset_sf <- function(in_pa_flag) {
  assets_proj |>
    mutate(in_pa = in_pa_flag, isic_section = ISIC_section) |>
    st_transform(4326)
}

assets_ssp1_sf     <- make_asset_sf(!is.na(vals_ssp1)     & vals_ssp1==1)
assets_ssp2_sf     <- make_asset_sf(!is.na(vals_ssp2)     & vals_ssp2==1)
assets_ssp5_sf     <- make_asset_sf(!is.na(vals_ssp5)     & vals_ssp5==1)
assets_baseline_sf <- make_asset_sf(!is.na(vals_baseline) & vals_baseline==1)

# ==============================================================================
# 18. Map helpers
#
#  Point size is mapped via aes(size=isic_section) + scale_size_manual so
#  mining points render at double the size of power points.
# ==============================================================================

make_scenario_map <- function(overlap_layers, assets_sf_all, title, subtitle = NULL,
                              sol_sf = NULL, mfg_polys = NULL) {
  
  has_layers <- is.list(overlap_layers) && !is.null(overlap_layers)
  has_mfg    <- !is.null(mfg_polys) && nrow(mfg_polys) > 0
  
  outside <- assets_sf_all |>
    filter(!in_pa, !(has_mfg & as.character(isic_section)=="Industrial"))
  inside  <- assets_sf_all |>
    filter( in_pa, !(has_mfg & as.character(isic_section)=="Industrial"))
  
  p <- ggplot() +
    geom_sf(data=countries_4326, fill="#f5f5f0", color="grey70", linewidth=0.2) +
    geom_sf(data=spain_4326,     fill="white",   color=NA)
  
  if (has_layers) {
    p <- p +
      geom_sf(data=overlap_layers$agri_sf,
              aes(fill=crop_fraction, alpha=crop_fraction), color=NA) +
      agri_fill_scale(order=3) +
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
        guide  = guide_legend(order=2, override.aes=list(alpha=1, size=5))
      )
  } else if (!is.null(sol_sf)) {
    p <- p + geom_sf(data=sol_sf |> filter(selected==1), fill=PA_ONLY_COLOUR, color=NA)
  } else {
    p <- p + geom_sf(data=PA_sf_plot, fill=PA_ONLY_COLOUR, color=NA, alpha=PA_ALPHA)
  }
  
  if (has_mfg) {
    p <- p +
      geom_sf(data=mfg_polys, fill=MFG_POLY_FILL,
              color=darken_hex(MFG_POLY_FILL), linewidth=0.15,
              alpha=MFG_POLY_ALPHA, show.legend=FALSE) +
      new_scale_fill() +
      geom_sf(data=mfg_polys[1,], aes(fill="mfg"), color=NA, alpha=0) +
      scale_fill_manual(name="Industry",
                        values=c("mfg"=MFG_POLY_FILL),
                        labels=c("mfg"="Industrial Zones"),
                        guide=guide_legend(order=4, override.aes=list(alpha=1, size=5)))
  }
  
  p +
    geom_sf(data=outside, aes(color=isic_section, size=isic_section), alpha=0.8) +
    geom_sf(data=inside,  aes(color=isic_section, size=isic_section), alpha=0.8) +
    scale_color_manual(
      name   = "ISIC section",
      values = section_colours_plot,
      labels = section_labels_plot,
      drop   = FALSE,
      guide  = guide_legend(order=1, override.aes=list(size=4))
    ) +
    scale_size_manual(values=section_sizes_plot, guide="none") +
    geom_sf(data=spain_4326, fill=NA, color="black", linewidth=0.45) +
    labs(title=title, subtitle=subtitle) +
    common_coord + map_theme +
    theme(legend.key.height=unit(1.0,"lines"),
          legend.key.width =unit(1.2,"lines"),
          legend.spacing.y =unit(0.3,"cm"))
}

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
      guide  = guide_legend(override.aes=list(alpha=1, size=5))
    ) +
    geom_sf(data=spain_4326, fill=NA, color="black", linewidth=0.45) +
    labs(title=title, subtitle=subtitle) +
    common_coord + map_theme +
    theme(plot.title    = element_text(face="bold", size=14),
          plot.subtitle = element_text(size=11))
}

make_sector_facet_panel <- function(sec, overlap_layers, vals_vec = NULL,
                                    mfg_polys = NULL, agri_km2 = NULL) {
  if (sec == "Agriculture, forestry and fishing") {
    make_agri_facet_map(
      overlap_layers,
      title    = section_label_short[sec],
      subtitle = if (!is.null(agri_km2))
        sprintf("Cropland in PA: %.0f km² (%.1f%%)",
                agri_km2, 100*agri_km2/max(total_agri_km2, 1))
      else NULL
    )
  } else {
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
      geom_sf(data=PA_sf_plot,     fill=PA_ONLY_COLOUR, color=NA, alpha=PA_ALPHA)
    
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
        guide  = guide_legend(override.aes=list(size=4))
      ) +
      scale_size_manual(values=section_sizes_plot, guide="none") +
      geom_sf(data=spain_4326, fill=NA, color="black", linewidth=0.45) +
      labs(title    = section_label_short[sec],
           subtitle = sprintf("Assets in PA: %d / %d  (%.1f%%)",
                              sum(in_pa), n_sec,
                              100*sum(in_pa)/max(n_sec,1))) +
      common_coord + map_theme +
      theme(plot.title    = element_text(face="bold", size=14),
            plot.subtitle = element_text(size=14))
  }
}

# ==============================================================================
# 19. Three-Panel Scenario Comparison Map
# ==============================================================================

panel_ssp1 <- make_scenario_map(
  overlap_ssp1, assets_ssp1_sf, mfg_polys=mfg_poly_ssp1,
  title = "SSP1 — 50% protection"
)
panel_ssp2 <- make_scenario_map(
  overlap_ssp2, assets_ssp2_sf, mfg_polys=mfg_poly_ssp2,
  title = "SSP2 — Regional rivalry"
)
panel_ssp5 <- make_scenario_map(
  overlap_ssp5, assets_ssp5_sf, mfg_polys=mfg_poly_ssp5,
  title = "SSP5 — Fossil-fuelled development"
)

combined_map <- (panel_ssp1 | panel_ssp2 | panel_ssp5) +
  plot_layout(guides="collect") +
  plot_annotation(
    title = "Economic assets and agricultural land within protected area boundaries",
    theme = theme(plot.title      = element_text(face="bold", size=13),
                  plot.background = element_rect(fill="white", color=NA))
  ) & theme(legend.position="bottom")

print(combined_map)
ggsave("solutions/assets_all_scenarios_map_1km.png",
       combined_map, width=25, height=9, dpi=300)
message("Combined scenario map saved.")

# ==============================================================================
# 20. Baseline Map
# ==============================================================================

baseline_map <- make_scenario_map(
  overlap_baseline, assets_baseline_sf, mfg_polys=mfg_poly_baseline,
  title = "Current protected areas and assets"
)

print(baseline_map)
ggsave("solutions/assets_baseline_map_1km.png",
       baseline_map, width=9, height=8, dpi=300)

# ==============================================================================
# 21. Facet Maps by Sector — Baseline
# ==============================================================================

baseline_facet_maps <- map(present_sections, function(sec)
  make_sector_facet_panel(sec, overlap_baseline,
                          mfg_polys=mfg_poly_baseline, agri_km2=agri_baseline_km2))

baseline_facet <- wrap_plots(baseline_facet_maps, ncol=n_cols_facet) +
  plot_layout(guides="collect") & theme(legend.position="bottom")

print(baseline_facet)
ggsave("solutions/assets_baseline_map_faced.png",
       baseline_facet,
       width  = 9 * n_cols_facet,
       height = 8 * ceiling(length(present_sections)/n_cols_facet),
       dpi=300, bg="white")

# ==============================================================================
# 22. Per-Scenario Facet Panels
# ==============================================================================

for (sc in list(
  list(vals=vals_ssp1, ovl=overlap_ssp1, mfg=mfg_poly_ssp1,
       agri=agri_ssp1_km2,
       title="Coordinated Nature Transition (SSP1 aligned)",
       file="assets_ssp1_facet_sections.png"),
  list(vals=vals_ssp2, ovl=overlap_ssp2, mfg=mfg_poly_ssp2,
       agri=agri_ssp2_km2,
       title="Falling Short — Middle of the Road (SSP2 aligned)",
       file="assets_ssp2_facet_sections.png"),
  list(vals=vals_ssp5, ovl=overlap_ssp5, mfg=mfg_poly_ssp5,
       agri=agri_ssp5_km2,
       title="Fossil-fuelled Development (SSP5 aligned)",
       file="assets_ssp5_facet_sections.png")
)) {
  panel_maps <- map(present_sections, function(sec)
    make_sector_facet_panel(sec, sc$ovl, vals_vec=sc$vals,
                            mfg_polys=sc$mfg, agri_km2=sc$agri))
  
  panel <- wrap_plots(panel_maps, ncol=n_cols_facet) +
    plot_layout(guides="collect") +
    plot_annotation(
      title = sc$title,
      theme = theme(plot.title      = element_text(face="bold", size=16),
                    plot.background = element_rect(fill="white", color=NA))
    ) & theme(legend.position="bottom")
  
  print(panel)
  ggsave(file.path("solutions", sc$file), panel,
         width  = 9 * n_cols_facet,
         height = 8 * ceiling(length(present_sections)/n_cols_facet),
         dpi=300, bg="white")
  message(sc$title, " facet panel saved.")
}

# ==============================================================================
# 23. Per-Sector Scenario Comparison
# ==============================================================================

message("\n=== SECTION 23: Per-sector scenario comparison ===")

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

scenario_list_23 <- list(
  list(label="Current (Baseline)",       ovl=overlap_baseline, vals=NULL,
       mfg=mfg_poly_baseline, agri_km2=agri_baseline_km2, pa_sf_plot=pa_baseline_sf_plot),
  list(label="SSP1 — Nature Transition", ovl=overlap_ssp1,     vals=vals_ssp1,
       mfg=mfg_poly_ssp1,    agri_km2=agri_ssp1_km2,    pa_sf_plot=pa_ssp1_sf_plot),
  list(label="SSP2 — Middle of the Road",ovl=overlap_ssp2,     vals=vals_ssp2,
       mfg=mfg_poly_ssp2,    agri_km2=agri_ssp2_km2,    pa_sf_plot=pa_ssp2_sf_plot),
  list(label="SSP5 — Fossil Development",ovl=overlap_ssp5,     vals=vals_ssp5,
       mfg=mfg_poly_ssp5,    agri_km2=agri_ssp5_km2,    pa_sf_plot=pa_ssp5_sf_plot)
)

to_slug <- function(x) tolower(gsub("[^a-zA-Z0-9]+", "_", x))

make_point_sector_panel_23 <- function(sec, vals_vec, mfg_polys, pa_sf) {
  
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
      guide  = guide_legend(override.aes=list(size=4))
    ) +
    scale_size_manual(values=section_sizes_plot, guide="none") +
    geom_sf(data=spain_4326, fill=NA, color="black", linewidth=0.45) +
    labs(title=NULL, subtitle=NULL) +
    common_coord + map_theme
}

for (sec in present_sections) {
  
  sec_panels <- map(scenario_list_23, function(sc) {
    if (sec == "Agriculture, forestry and fishing") {
      p <- make_agri_facet_map(sc$ovl, title=NULL, subtitle=NULL)
    } else {
      p <- make_point_sector_panel_23(sec, sc$vals, sc$mfg, sc$pa_sf_plot)
    }
    p + labs(title=sc$label) +
      theme(plot.title = element_text(face="bold", size=14))
  })
  
  sector_plot <- wrap_plots(sec_panels, nrow=1) +
    plot_layout(guides="collect") +
    plot_annotation(
      title = section_label_short[sec],
      theme = theme(
        plot.title      = element_text(face="bold", size=20),
        plot.background = element_rect(fill="white", color=NA)
      )
    ) & theme(legend.position="bottom")
  
  out_file <- file.path(
    "solutions",
    sprintf("assets_%s_scenarios_compared.png", to_slug(section_label_short[sec]))
  )
  
  print(sector_plot)
  ggsave(out_file, sector_plot, width=9*4, height=8, dpi=300, bg="white")
  message(sprintf("  Saved: %s", out_file))
}

message("\n=== ALL OUTPUTS COMPLETE ===")