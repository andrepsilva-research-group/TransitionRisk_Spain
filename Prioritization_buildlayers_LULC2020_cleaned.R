setwd("/Users/fishfeed/Library/CloudStorage/OneDrive-StockholmUniversity/FinBio/R_scripts/NatPoKe-dev/NatPoKe-dev_Ines")

source("src/customFunctions2.R")

library(prioritizr)
library(sf)
library(terra)
library(dplyr)
library(tibble)
library(ggplot2)
library(gurobi)
library(rnaturalearth)

library(foreign)
library(geodata)

set.seed(500)
terraOptions(threads = parallel::detectCores() - 1, memfraction = 0.8)

dir.create("Prioritization/cache",        showWarnings = FALSE, recursive = TRUE)
dir.create("Prioritization/Input_layers", showWarnings = FALSE, recursive = TRUE)

# Country Boundaries ----
pt <- ne_countries(scale = "medium", country = "Portugal", returnclass = "sf")
es <- ne_countries(scale = "medium", country = "Spain",    returnclass = "sf")

countries   <- rbind(pt, es)
countries   <- st_transform(countries, 3035)
countries_v <- vect(countries)

spain   <- countries %>% filter(admin == "Spain")
spain_v <- vect(spain)

# Scenario Definitions ----
scenarios <- list(
  current = list(

    lulc = "Prioritization/Input_layers/SPAIN_SEALS/lulc_glc_fcs30d_2020_pen.tif",
    suit = "~/Library/CloudStorage/OneDrive-StockholmUniversity/FinBio/R_scripts/NatPoKe-dev/NatPoKe-dev_Ines/Prioritization/Input_layers/SuitabilityAverages/scenarioAverages/Current_mean.tif"
  ),
  esa2015 = list(

    lulc = "~/Library/CloudStorage/OneDrive-StockholmUniversity/FinBio/R_scripts/NatPoKe-dev/NatPoKe-dev/data/stitched_lulc_esa_scenarios/lulc_esa_2015.tif",
    suit = "~/Library/CloudStorage/OneDrive-StockholmUniversity/FinBio/R_scripts/NatPoKe-dev/NatPoKe-dev_Ines/Prioritization/Input_layers/SuitabilityAverages/scenarioAverages/Current_mean.tif"
  ),
  ssp585_2050 = list(
    lulc = "Prioritization/Input_layers/SPAIN_SEALS/lulc_glc_fcs30d_seals7_ssp5_rcp85_seals_bau_2050_clipped_1deg_coefs.tif",
    suit = "~/Library/CloudStorage/OneDrive-StockholmUniversity/FinBio/R_scripts/NatPoKe-dev/NatPoKe-dev_Ines/Prioritization/Input_layers/SuitabilityAverages/scenarioAverages/ssp585_2050_mean.tif"
  )
)

# Planning Unit Template ----
pu_template <- rast(
  ext(countries_v),
  resolution = 1000,
  crs        = crs(countries_v)
)

# Cost Layer ----
cost_cache <- "Prioritization/cache/pu_cost_raw.tif"

if (file.exists(cost_cache)) {

  pu_cost <- rast(cost_cache)
} else {

  cost_nc           <- rast("Prioritization/spencer_ecoser_costlayer.nc")
  cost_nc           <- project(cost_nc, crs(pu_template))
  names(cost_nc)[6] <- "combined_cost"
  land_cost         <- cost_nc[["combined_cost"]]

  land_cost_1km   <- resample(land_cost, pu_template, method = "bilinear")
  pu_cost         <- mask(land_cost_1km, spain_v)
  pu_cost         <- trim(pu_cost)
  names(pu_cost)  <- "cost"

  writeRaster(pu_cost, cost_cache, overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))

}

# WDPA Protected Areas ----
wdpa_base <- "Prioritization/WDPA_WDOECM_Apr2026_Public_ESP_shp"

poly_paths <- c(
  file.path(wdpa_base,
            "WDPA_WDOECM_Apr2026_Public_ESP_shp_0",
            "WDPA_WDOECM_Apr2026_Public_ESP_shp-polygons.shp"),
  file.path(wdpa_base,
            "WDPA_WDOECM_Apr2026_Public_ESP_shp_1",
            "WDPA_WDOECM_Apr2026_Public_ESP_shp-polygons.shp"),
  file.path(wdpa_base,
            "WDPA_WDOECM_Apr2026_Public_ESP_shp_2",
            "WDPA_WDOECM_Apr2026_Public_ESP_shp-polygons.shp")
)

valid_status <- c("Designated", "Established", "Inscribed")

wdpa_combined_cache  <- "Prioritization/cache/wdpa_combined.gpkg"
wdpa_dissolved_cache <- "Prioritization/cache/wdpa_spain_dissolved.gpkg"

if (file.exists(wdpa_combined_cache) && file.exists(wdpa_dissolved_cache)) {

  wdpa_combined        <- st_read(wdpa_combined_cache,  quiet = TRUE)
  wdpa_spain_dissolved <- st_read(wdpa_dissolved_cache, quiet = TRUE)

} else {

  poly_list <- lapply(seq_along(poly_paths), function(i) {
    p <- poly_paths[i]
    if (!file.exists(p)) {
       return(NULL)
    }
    raw      <- st_read(p, quiet = TRUE)
    filtered <- raw %>%
      filter(STATUS %in% valid_status) %>%
      filter(REALM != "MARINE") %>%
      st_transform(3035) %>%
      st_make_valid()

    filtered
  })

  poly_list        <- Filter(Negate(is.null), poly_list)
  shared_poly_cols <- Reduce(intersect, lapply(poly_list, names))
  wdpa_combined    <- bind_rows(lapply(poly_list, function(x) x[, shared_poly_cols]))
  wdpa_combined    <- st_make_valid(wdpa_combined)

  st_write(wdpa_combined, wdpa_combined_cache, delete_dsn = TRUE, quiet = TRUE)

  spain_boundary_union <- st_union(st_transform(spain, 3035))
  wdpa_spain_clip      <- suppressWarnings(st_intersection(wdpa_combined, spain_boundary_union))
  wdpa_spain_dissolved <- st_union(wdpa_spain_clip)
  wdpa_spain_dissolved <- st_as_sf(data.frame(geometry = st_geometry(wdpa_spain_dissolved)))

  st_write(wdpa_spain_dissolved, wdpa_dissolved_cache, delete_dsn = TRUE, quiet = TRUE)

}

spain_boundary_union <- st_union(st_transform(spain, 3035))
spain_total_km2      <- as.numeric(st_area(spain_boundary_union)) / 1e6
pa_area_km2          <- as.numeric(st_area(st_union(wdpa_spain_dissolved))) / 1e6

# PA Threshold Rasterization ----
PA_rast_cache <- "Prioritization/cache/PA_rast.tif"
area_cache    <- "Prioritization/cache/area_raster.tif"

if (file.exists(PA_rast_cache) && file.exists(area_cache)) {

  PA_rast     <- resample(rast(PA_rast_cache), pu_cost, method = "near")
  area_raster <- resample(rast(area_cache),    pu_cost, method = "bilinear")
} else {

  wdpa_v_combined <- vect(wdpa_combined)

  pa_cover <- rasterize(
    wdpa_v_combined, pu_cost,
    fun   = "sum",
    cover = TRUE
  )

  PA_rast <- as.logical(pa_cover >= 0.50)

  pu_cost[PA_rast] <- 0
  writeRaster(pu_cost, cost_cache, overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))

  area_raster <- cellSize(pu_cost, unit = "km")

  writeRaster(PA_rast,     PA_rast_cache, overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))
  writeRaster(area_raster, area_cache,    overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))

}

pu_cost[PA_rast] <- 0

PA_area_km2         <- global(area_raster * as.numeric(mask(PA_rast, spain_v)),
                              "sum", na.rm = TRUE)[1, 1]
total_area_km2_rast <- global(mask(area_raster, spain_v), "sum", na.rm = TRUE)[1, 1]

# Renormalize Cost Layer ----
cost_renorm_cache <- "Prioritization/cache/pu_cost_renorm.tif"

if (file.exists(cost_renorm_cache)) {

  pu_cost <- rast(cost_renorm_cache)

} else {

  cost_spain_raw <- mask(pu_cost, spain_v)

  cost_land_only <- cost_spain_raw
  cost_land_only[PA_rast]              <- NA
  cost_land_only[cost_land_only == 0]  <- NA

  cost_min_land <- terra::global(cost_land_only, "min", na.rm = TRUE)[1, 1]
  cost_max_land <- terra::global(cost_land_only, "max", na.rm = TRUE)[1, 1]

  cost_renorm <- (cost_spain_raw - cost_min_land) / (cost_max_land - cost_min_land)
  cost_renorm <- clamp(cost_renorm, lower = 0, upper = 1)

  cost_renorm[PA_rast] <- 0

  cost_renorm        <- mask(cost_renorm, spain_v)
  names(cost_renorm) <- "cost"

  writeRaster(cost_renorm, cost_renorm_cache, overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))

  pu_cost <- cost_renorm
}

PA_rast <- resample(PA_rast, pu_cost, method = "near")
PA_rast <- as.logical(PA_rast)
pu_cost[PA_rast] <- 0

# GADM CCAAs ----
spain_ccaa_all <- geodata::gadm(
  country = "ESP", level = 1,
  path    = "Prioritization/GADM"
) %>%
  st_as_sf() %>%
  st_transform(3035)

# ETM / EUNIS Features ----
etm_save_path <- "Prioritization/Input_layers/EUNIS_features_prop_1km.tif"

if (file.exists(etm_save_path)) {

  EUNIS_features_prop <- rast(etm_save_path)

} else {

  ecotype_rast <- rast("Prioritization/Input_layers/ecotype_rast_projected.tif")
  ecotype_rast <- project(ecotype_rast, crs(pu_cost), method = "near")
  ecotype_rast <- crop(ecotype_rast,  countries_v)
  ecotype_rast <- mask(ecotype_rast,  countries_v)

  value_to_EUNIS1 <- c(
    11, 2,   12, 2,
    13, 3,   14, 3,   15, 3,   16, 3,
    17, 4,   18, 4,   20, 4,   21, 4,   22, 4,
    23, 5,   24, 5,   25, 5,   26, 5,   28, 5,   29, 5,
    31, 6,   32, 6,   33, 6,   34, 6,   35, 6,   36, 6,   37, 6,   38, 6,   40, 6,
    41, 7,   42, 7,   43, 7,   44, 7,   45, 7,
    47, 8,   48, 8,   49, 8,   50, 8,
    52, 9,   53, 9,   54, 9,   55, 9,   56, 9,   57, 9,   58, 9,   59, 9
  )

  EUNIS1matrix    <- matrix(value_to_EUNIS1, ncol = 2, byrow = TRUE)
  EUNIS_1_highres <- terra::classify(ecotype_rast, EUNIS1matrix,
                                     include.lowest = TRUE, others = NA)
  names(EUNIS_1_highres) <- "EUNIS_1"

  all_etm_classes <- sort(unique(values(EUNIS_1_highres)))
  all_etm_classes <- all_etm_classes[!is.na(all_etm_classes)]
  classes_etm     <- all_etm_classes[!all_etm_classes %in% c(2, 3, 9)]

  EUNIS_features_prop <- rast(lapply(classes_etm, function(cl) {
    r             <- as.numeric(EUNIS_1_highres == cl)
    r[is.na(r)]   <- 0
    r_prop        <- resample(r, pu_cost, method = "average")
    r_prop        <- mask(r_prop, countries_v)
    names(r_prop) <- paste0("EUNIS_", cl, "_prop")
    r_prop
  }))

  writeRaster(EUNIS_features_prop, etm_save_path, overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))

}

n_eunis_layers <- nlyr(EUNIS_features_prop)

# LULC Features ----
lulc_to_EUNIS1 <- c(
  10,  9,   11,  9,   12,  9,   20,  9,   30,  9,   190, 9,
  40,  6,
  50,  7,   51,  7,   52,  7,   60,  7,   61,  7,   62,  7,
  70,  7,   71,  7,   72,  7,   80,  7,   81,  7,   82,  7,
  90,  7,   91,  7,   92,  7,   100, 7,
  110, 6,   120, 6,   121, 6,   122, 6,
  130, 5,
  140, 8,
  150, 8,   151, 8,   152, 8,   153, 8,
  170, 4,   180, 4,
  181, 4,   182, 4,   183, 4,   184, 4,   185, 4,   186, 4,   187, 4,
  200, 8,   201, 8,   202, 8,
  210, 3,
  220, 3
)

lulc_EUNIS_matrix    <- matrix(lulc_to_EUNIS1, ncol = 2, byrow = TRUE)
classes_to_keep_lulc <- c(4, 5, 6, 7, 8)

seals7_to_EUNIS1 <- c(
  1, 9,
  2, 9,
  3, 5,
  4, 7,
  5, 6,
  6, 3,
  7, 8
)
seals7_EUNIS_matrix <- matrix(seals7_to_EUNIS1, ncol = 2, byrow = TRUE)
seals7_scenarios    <- c("ssp585_2050")

for (sc in names(scenarios)) {
  save_path <- paste0("Prioritization/Input_layers/lulc_iberia_projected_", sc, "_1km.tif")
  if (file.exists(save_path)) {  next }
  if (!file.exists(scenarios[[sc]]$lulc)) {
    warning("Source LULC file not found for scenario '", sc, "' -- skipping: ", scenarios[[sc]]$lulc)
    next
  }

  lulc_raw <- rast(scenarios[[sc]]$lulc)

  if (crs(lulc_raw) == crs(countries_v)) {
    lulc_crop_native <- crop(lulc_raw, countries_v)
  } else {
    countries_in_raw_crs <- tryCatch(
      project(countries_v, crs(lulc_raw)),
      error = function(e) NULL
    )
    if (is.null(countries_in_raw_crs)) {
      warning("Could not reproject country boundary into the CRS of scenario '",
              sc, "' -- skipping: ", scenarios[[sc]]$lulc)
      next
    }
    lulc_crop_native <- crop(lulc_raw, countries_in_raw_crs)
  }

  lulc_proj <- project(lulc_crop_native, crs(pu_cost), method = "near")
  lulc_proj <- crop(lulc_proj,   countries_v)
  lulc_proj <- mask(lulc_proj,   countries_v)
  writeRaster(lulc_proj, save_path, overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))

}

available_scenarios <- names(scenarios)[
  file.exists(paste0("Prioritization/Input_layers/lulc_iberia_projected_", names(scenarios), "_1km.tif"))
]

lulc_projected <- lapply(available_scenarios, function(sc) {
  rast(paste0("Prioritization/Input_layers/lulc_iberia_projected_", sc, "_1km.tif"))
})
names(lulc_projected) <- available_scenarios

terraOptions(threads = parallel::detectCores() - 1, memfraction = 0.6)
lulc_features_list <- list()

for (sc in available_scenarios) {

  save_path <- paste0("Prioritization/Input_layers/LULC_EUNIS_prop_", sc, "_1km.tif")

  if (file.exists(save_path)) {

    lulc_features_list[[sc]] <- rast(save_path)
    next
  }

  lulc_crop <- lulc_projected[[sc]]

  eunis_matrix_sc   <- if (sc %in% seals7_scenarios) seals7_EUNIS_matrix else lulc_EUNIS_matrix
  lulc_EUNIS        <- terra::classify(lulc_crop, eunis_matrix_sc,
                                       include.lowest = TRUE, others = NA)
  names(lulc_EUNIS) <- "EUNIS_1"

  present <- freq(lulc_EUNIS)$value
  present <- sort(present[!is.na(present)])

  keep_matrix <- cbind(classes_to_keep_lulc, classes_to_keep_lulc)
  lulc_EUNIS_filtered <- terra::classify(lulc_EUNIS, keep_matrix, others = NA)
  all_classes_binary <- segregate(lulc_EUNIS_filtered, other = 0, keep = FALSE)

  agg_factor <- round(res(pu_cost)[1] / res(lulc_EUNIS)[1])
  feat_agg   <- aggregate(all_classes_binary, fact = agg_factor,
                          fun = "mean", na.rm = TRUE)
  feat_prop  <- resample(feat_agg, pu_cost, method = "bilinear")
  feat_prop  <- mask(feat_prop, countries_v)

  detected_classes <- sort(present[present %in% classes_to_keep_lulc])
  names(feat_prop) <- paste0("LULC_EUNIS_", detected_classes, "_prop")

  writeRaster(feat_prop, save_path, overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))

  lulc_features_list[[sc]] <- feat_prop
}

# Suitability Layers ----
suit_layers_norm <- list()

for (sc in available_scenarios) {

  save_path <- paste0("Prioritization/Input_layers/suit_", sc, "_1km.tif")

  if (file.exists(save_path)) {

    r <- rast(save_path)
  } else {

    r           <- rast(scenarios[[sc]]$suit)
    r           <- project(r, crs(pu_cost), method = "bilinear")
    r           <- crop(r,  countries_v)
    r           <- mask(r,  countries_v)
    r           <- resample(r, pu_cost, method = "bilinear")
    r[is.na(r)] <- 0
    names(r)    <- paste0("suit_", sc)
    writeRaster(r, save_path, overwrite = TRUE,
                gdal = c("COMPRESS=LZW", "TILED=YES"))

  }

  max_val            <- terra::global(mask(r, spain_v), "max", na.rm = TRUE)[1, 1]
  r_norm             <- r / max_val
  names(r_norm)      <- paste0("suit_", sc, "_norm")
  suit_layers_norm[[sc]] <- r_norm
}

# NCP Layer ----
ncp_cache <- "Prioritization/cache/ncp_rast.tif"

if (file.exists(ncp_cache)) {

  ncp_rast <- rast(ncp_cache)

} else {

  es_r  <- rast("Prioritization/mask_esp_eck4_2km.tif")
  es_90 <- rast("Prioritization/solution_scenario-A_esp_target-90.tif")

  es_r  <- project(es_r,  "EPSG:3035")
  es_90 <- project(es_90, es_r, method = "near")

  es_90_masked    <- mask(es_90, es_r)
  es_90_3035      <- project(es_90_masked, "EPSG:3035", method = "near")

  ncp_rast        <- resample(es_90_3035, pu_cost, method = "near")
  ncp_rast[is.na(ncp_rast)] <- 0
  ncp_max         <- terra::global(ncp_rast, "max", na.rm = TRUE)[1, 1]
  ncp_rast        <- ncp_rast / ncp_max
  names(ncp_rast) <- "ncp_priority"

  writeRaster(ncp_rast, ncp_cache, overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))

}

# Per-CCAA Baseline PA% ----
normalize_name <- function(x) {
  x %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    tolower() %>%
    trimws() %>%
    gsub("\\s+", " ", .)
}

ssp2_params <- tribble(
  ~NAME_1,                       ~budget_mult,
  "Andalucía",                    1.00,
  "Aragón",                       1.02,
  "Cantabria",                    1.00,
  "Castilla y León",              1.07,
  "Castilla-La Mancha",           1.11,
  "Cataluña",                     1.00,
  "Comunidad de Madrid",          1.00,
  "Comunidad Foral de Navarra",   1.00,
  "Comunidad Valenciana",         1.00,
  "Extremadura",                  1.00,
  "Galicia",                      2.45,
  "Islas Baleares",               1.00,
  "Islas Canarias",               1.00,
  "La Rioja",                     1.00,
  "País Vasco",                   1.27,
  "Principado de Asturias",       1.00,
  "Región de Murcia",             1.00
)

params_norm <- ssp2_params %>% mutate(name_key = normalize_name(NAME_1))

spain_ccaa <- spain_ccaa_all %>%
  mutate(name_key = normalize_name(NAME_1)) %>%
  inner_join(params_norm %>% select(name_key, budget_mult), by = "name_key") %>%
  select(-name_key)

still_missing <- setdiff(
  normalize_name(ssp2_params$NAME_1),
  normalize_name(spain_ccaa$NAME_1)
)

baselines_cache <- "Prioritization/cache/ccaa_baselines_poly.rds"

if (file.exists(baselines_cache)) {

  ccaa_baselines_df <- readRDS(baselines_cache)

} else {

  wdpa_spain_all <- suppressWarnings(
    st_intersection(wdpa_combined, st_union(spain_ccaa_all))
  )

  ccaa_baselines_poly <- lapply(seq_len(nrow(spain_ccaa)), function(i) {
    ccaa_row <- spain_ccaa[i, ]
    name_i   <- ccaa_row$NAME_1
    total_i  <- as.numeric(st_area(ccaa_row)) / 1e6

    pa_clip <- suppressWarnings(
      st_intersection(wdpa_spain_all, st_union(ccaa_row))
    )

    if (nrow(pa_clip) == 0 || all(st_is_empty(pa_clip))) {

      return(tibble(NAME_1 = name_i, total_km2 = total_i,
                    pa_km2_poly = 0, baseline_pct_poly = 0))
    }

    pa_dissolved <- st_union(pa_clip)
    pa_km2       <- as.numeric(st_area(pa_dissolved)) / 1e6

    tibble(
      NAME_1            = name_i,
      total_km2         = total_i,
      pa_km2_poly       = pa_km2,
      baseline_pct_poly = 100 * pa_km2 / total_i
    )
  })

  ccaa_baselines_df <- bind_rows(ccaa_baselines_poly) %>%
    left_join(
      spain_ccaa %>% st_drop_geometry() %>% select(NAME_1, budget_mult),
      by = "NAME_1"
    ) %>%
    mutate(ssp2_pct = pmin(baseline_pct_poly * budget_mult, 95))

  saveRDS(ccaa_baselines_df, baselines_cache)

}

national_poly_pct <- sum(ccaa_baselines_df$pa_km2_poly) /
  sum(ccaa_baselines_df$total_km2) * 100

spain_ccaa <- spain_ccaa %>%
  left_join(
    ccaa_baselines_df %>% select(NAME_1, baseline_pct_poly, ssp2_pct),
    by = "NAME_1"
  )

dir.create("Prioritization/solutions", showWarnings = FALSE, recursive = TRUE)
write.csv(ccaa_baselines_df,
          "Prioritization/solutions/ccaa_baselines_poly.csv",
          row.names = FALSE)

# Artificial Land Cover Proportion Layer (SSP5) ----
artificial_codes <- c(1L, 2L)

artificial_ssp5_cache <- "Prioritization/cache/artificial_prop_ssp5.tif"

if (file.exists(artificial_ssp5_cache)) {

  artificial_prop_ssp5 <- rast(artificial_ssp5_cache)

} else if (file.exists("Prioritization/Input_layers/lulc_iberia_projected_ssp585_2050_1km.tif")) {

  lulc_ssp5_cat <- rast(
    "Prioritization/Input_layers/lulc_iberia_projected_ssp585_2050_1km.tif"
  )

  artificial_ssp5_native    <- ifel(lulc_ssp5_cat %in% artificial_codes, 1L, 0L)
  names(artificial_ssp5_native) <- "artificial"
  artificial_prop_ssp5      <- resample(artificial_ssp5_native, pu_cost, method = "average")
  artificial_prop_ssp5      <- mask(artificial_prop_ssp5, spain_v)
  artificial_prop_ssp5      <- clamp(artificial_prop_ssp5, lower = 0, upper = 1)

  writeRaster(artificial_prop_ssp5, artificial_ssp5_cache, overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))

}
