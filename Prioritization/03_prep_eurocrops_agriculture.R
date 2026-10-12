#prep_eurcrop
setwd("path/to/Spain_TranRisk/Prioritization")

library(terra); library(sf); library(arrow)
library(rnaturalearth); library(dplyr); library(readr); library(geodata)

terraOptions(threads = parallel::detectCores() - 1)

PARQUET_PATH       <- "input_layers_assets/es_stack.parquet"
OUT_DIR            <- "solutions"
RAST_OUT           <- "data/agri/eurcrops_spain_1km.tif"       # binary threshold raster (backward-compatible)
FRAC_RAST_OUT      <- "data/agri/eurcrops_spain_1km_frac.tif"  # NEW: fractional coverage raster
FRAC_THRESHOLD     <- 0.5    # cells with < 50 % cropland coverage are excluded from binary layer
EXCLUDE_CROP_CODES <- NULL    # e.g. c("62.0","63.0") after Run 1

dir.create(OUT_DIR,           recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(RAST_OUT), recursive = TRUE, showWarnings = FALSE)

# terra >= 1.7.39 required for rasterizeGeom(..., fun = "coverage")
stopifnot(packageVersion("terra") >= "1.7.39")

# Replace accented characters with ASCII — works reliably on all platforms
remove_accents <- function(x) {
  x <- gsub("[áàâä]", "a", x); x <- gsub("[ÁÀÂÄ]", "A", x)
  x <- gsub("[éèêë]", "e", x); x <- gsub("[ÉÈÊË]", "E", x)
  x <- gsub("[íìîï]", "i", x); x <- gsub("[ÍÌÎÏ]", "I", x)
  x <- gsub("[óòôö]", "o", x); x <- gsub("[ÓÒÔÖ]", "O", x)
  x <- gsub("[úùûü]", "u", x); x <- gsub("[ÚÙÛÜ]", "U", x)
  x <- gsub("ñ",      "n", x); x <- gsub("Ñ",      "N", x)
  x
}

# ==============================================================================
# 1. Spain boundary
# ==============================================================================

pt <- ne_countries(scale="medium", country="Portugal", returnclass="sf")
es <- ne_countries(scale="medium", country="Spain",    returnclass="sf")
countries  <- rbind(pt, es) |> st_transform(3035)
spain      <- countries |> filter(admin == "Spain")
spain_v    <- vect(spain)
spain_bbox <- st_bbox(spain)

# ==============================================================================
# 2. PA and scenario rasters
# ==============================================================================

PA_rast  <- rast("cache/PA_rast.tif")
sol_ssp1 <- rast("solutions/solution_ssp1.tif")
sol_ssp2 <- rast("solutions/solution_ssp2_regional.tif")
sol_ssp5 <- rast("solutions/solution_ssp5_development.tif")

PA_rast           <- resample(PA_rast, sol_ssp1, method = "near")
rast_stack        <- c(PA_rast, sol_ssp1, sol_ssp2, sol_ssp5)
names(rast_stack) <- c("pa_baseline","pa_ssp1","pa_ssp2","pa_ssp5")
template_r        <- sol_ssp1; values(template_r) <- NA

# ==============================================================================
# 2b. Province boundaries (GADM level 2 — 52 Spanish provinces)
# ==============================================================================

prov_gadm <- geodata::gadm("ESP", level=2, path="GADM") |>
  st_as_sf() |>
  st_transform(3035) |>
  mutate(
    prov_id       = row_number(),
    province_name = remove_accents(NAME_2),
    ccaa_name     = remove_accents(NAME_1),
    province_code = CC_2
  )

prov_rast        <- rasterize(vect(prov_gadm), template_r, field="prov_id")
names(prov_rast) <- "prov_id"

prov_meta <- prov_gadm |>
  st_drop_geometry() |>
  select(prov_id, province_code, province_name, ccaa_name)

# ==============================================================================
# 3. Read parquet
# ==============================================================================

tbl <- open_dataset(PARQUET_PATH) |>
  filter(!is.na(c2023)) |>
  collect() |>
  select(cropfield, c2023, area_ha, geometry)

# ==============================================================================
# 3b. Crop type lookup
# ==============================================================================

es_map <- NULL

tryCatch({
  es_map <- read_csv(
    "https://raw.githubusercontent.com/maja601/EuroCrops/main/csvs/country_mappings/es_2023.csv",
    show_col_types = FALSE
  )
  crop_lookup <- open_dataset(PARQUET_PATH) |>
    filter(!is.na(c2023)) |>
    group_by(c2023) |>
    summarise(total_km2=sum(area_ha, na.rm=TRUE)/100, n_parcels=n()) |>
    collect() |>
    mutate(original_code=as.numeric(c2023)) |>
    left_join(es_map, by="original_code") |>
    arrange(desc(total_km2)) |>
    select(original_code, original_name, translated_name,
           HCAT3_name, HCAT3_code, n_parcels, total_km2)
  write_csv(crop_lookup, file.path(OUT_DIR, "agriculture_crop_lookup.csv"))
}, error = function(e) NULL)

# Apply exclusions
if (!is.null(EXCLUDE_CROP_CODES) && length(EXCLUDE_CROP_CODES) > 0) {
  tbl      <- tbl |> filter(!c2023 %in% EXCLUDE_CROP_CODES)
}

# ==============================================================================
# 4. Chunk loop
#    result_list     — PA parcels only (PA overlap CSVs)
#    all_parcel_list — all Spain parcels with province + cropfield
# ==============================================================================

chunk_size      <- 100000
n_chunks        <- ceiling(nrow(tbl) / chunk_size)
result_list     <- vector("list", n_chunks)
all_parcel_list <- vector("list", n_chunks)

for (i in seq_len(n_chunks)) {
  
  idx        <- ((i-1)*chunk_size + 1):min(i*chunk_size, nrow(tbl))
  chunk      <- tbl[idx, ]
  sfc        <- st_as_sfc(chunk$geometry, crs=3035)
  cents      <- st_centroid(sfc)
  coords     <- st_coordinates(cents)
  
  keep <- coords[,"X"] >= spain_bbox["xmin"] &
    coords[,"X"] <= spain_bbox["xmax"] &
    coords[,"Y"] >= spain_bbox["ymin"] &
    coords[,"Y"] <= spain_bbox["ymax"]
  if (sum(keep) == 0) next
  
  chunk_keep  <- chunk[keep, ]
  coords_keep <- coords[keep, , drop=FALSE]
  
  pts       <- vect(coords_keep, crs="EPSG:3035")
  ext_vals  <- terra::extract(rast_stack, pts)[, -1]
  prov_vals <- terra::extract(prov_rast,  pts)[, 2]
  
  # All parcels — cropfield retained for PA flag join in Section 7
  all_parcel_list[[i]] <- data.frame(
    cropfield = chunk_keep$cropfield,
    crop_code = chunk_keep$c2023,
    area_km2  = chunk_keep$area_ha / 100,
    prov_id   = prov_vals,
    stringsAsFactors = FALSE
  )
  
  # PA parcels only
  chunk_result <- bind_cols(
    chunk_keep |> select(cropfield, crop_code=c2023, area_km2=area_ha) |>
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
# NOTE: tbl is NOT removed here — kept alive for Section 4b rasterizeGeom pass
rm(result_list, all_parcel_list); gc()

# ==============================================================================
# 4b. Fractional cropland raster via rasterizeGeom()
# ==============================================================================

CHUNK_FRAC    <- 500000
n_chunks_frac <- ceiling(nrow(tbl) / CHUNK_FRAC)
frac_accum    <- template_r; values(frac_accum) <- 0.0
cell_area_m2  <- prod(res(template_r))   # computed once — m² per cell

for (i in seq_len(n_chunks_frac)) {
  
  idx <- ((i-1)*CHUNK_FRAC + 1):min(i*CHUNK_FRAC, nrow(tbl))
  sfc <- st_as_sfc(tbl$geometry[idx], crs=3035)
  
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
  cf[cf > 1]    <- 1        # floating-point safety cap
  frac_accum    <- frac_accum + cf
  frac_accum[frac_accum > 1] <- 1   # cumulative cap across chunks
  
}

# --- NEW: save fractional raster, then derive + save binary raster ---------
frac_accum <- mask(frac_accum, spain_v)
writeRaster(frac_accum, FRAC_RAST_OUT, overwrite = TRUE)

binary_accum <- ifel(frac_accum >= FRAC_THRESHOLD, 1L, NA)
writeRaster(binary_accum, RAST_OUT, overwrite = TRUE)
# -----------------------------------------------------------------------------

# ==============================================================================
# 5. Summary statistics
# ==============================================================================

total_agri_km2 <- {
  ds <- open_dataset(PARQUET_PATH) |> filter(!is.na(c2023))
  if (!is.null(EXCLUDE_CROP_CODES) && length(EXCLUDE_CROP_CODES) > 0)
    ds <- ds |> filter(!c2023 %in% EXCLUDE_CROP_CODES)
  ds |> summarise(total=sum(area_ha, na.rm=TRUE)/100) |>
    collect() |> pull(total)
}

km2_in_pa <- function(col)
  sum(pa_parcels$area_km2[pa_parcels[[col]]==1L], na.rm=TRUE)

agri_baseline_km2 <- km2_in_pa("pa_baseline")
agri_ssp1_km2     <- km2_in_pa("pa_ssp1")
agri_ssp2_km2     <- km2_in_pa("pa_ssp2")
agri_ssp5_km2     <- km2_in_pa("pa_ssp5")

n_in_pa <- function(col)
  sum(pa_parcels[[col]] == 1L, na.rm = TRUE)

agri_baseline_n      <- n_in_pa("pa_baseline")
agri_ssp1_n          <- n_in_pa("pa_ssp1")
agri_ssp2_n          <- n_in_pa("pa_ssp2")
agri_ssp5_n          <- n_in_pa("pa_ssp5")
total_agri_n_parcels <- nrow(all_parcels)

# ==============================================================================
# 6. Save CSV outputs
# ==============================================================================

# 6a. Impact summary
agriculture_summary <- tibble(
  Sector             = "Agriculture & Forestry (EuroCrops 2023)",
  excluded_codes     = ifelse(is.null(EXCLUDE_CROP_CODES), "none",
                             paste(EXCLUDE_CROP_CODES, collapse=";")),
  frac_threshold     = FRAC_THRESHOLD,
  total_area_km2     = round(total_agri_km2,    1),
  area_baseline_km2  = round(agri_baseline_km2, 1),
  area_ssp1_km2      = round(agri_ssp1_km2,     1),
  area_ssp2_km2      = round(agri_ssp2_km2,     1),
  area_ssp5_km2      = round(agri_ssp5_km2,     1),
  pct_baseline       = round(agri_baseline_km2/total_agri_km2*100, 1),
  pct_ssp1           = round(agri_ssp1_km2/total_agri_km2*100,     1),
  pct_ssp2           = round(agri_ssp2_km2/total_agri_km2*100,     1),
  pct_ssp5           = round(agri_ssp5_km2/total_agri_km2*100,     1),
  total_n_parcels    = total_agri_n_parcels,
  n_parcels_baseline = agri_baseline_n,
  n_parcels_ssp1     = agri_ssp1_n,
  n_parcels_ssp2     = agri_ssp2_n,
  n_parcels_ssp5     = agri_ssp5_n,
  pct_n_baseline     = round(agri_baseline_n/total_agri_n_parcels*100, 1),
  pct_n_ssp1         = round(agri_ssp1_n/total_agri_n_parcels*100,     1),
  pct_n_ssp2         = round(agri_ssp2_n/total_agri_n_parcels*100,     1),
  pct_n_ssp5         = round(agri_ssp5_n/total_agri_n_parcels*100,     1)
)
write_csv(agriculture_summary, file.path(OUT_DIR, "impact_agriculture.csv"))

# 6b. Per-parcel (PA only)
write_csv(pa_parcels, file.path(OUT_DIR, "agriculture_in_pa_parcels.csv"))

# 6c. National crop type summary
crop_summary <- pa_parcels |>
  group_by(crop_code) |>
  summarise(
    n_parcels    = n(),
    total_km2    = round(sum(area_km2,                  na.rm=TRUE), 2),
    km2_baseline = round(sum(area_km2[pa_baseline==1L], na.rm=TRUE), 2),
    km2_ssp1     = round(sum(area_km2[pa_ssp1==1L],     na.rm=TRUE), 2),
    km2_ssp2     = round(sum(area_km2[pa_ssp2==1L],     na.rm=TRUE), 2),
    km2_ssp5     = round(sum(area_km2[pa_ssp5==1L],     na.rm=TRUE), 2),
    .groups = "drop"
  ) |>
  mutate(original_code = as.numeric(crop_code)) |>
  arrange(desc(total_km2))

if (!is.null(es_map)) {
  crop_summary <- left_join(crop_summary,
                            es_map |> select(original_code, translated_name, HCAT3_name),
                            by = "original_code")
} else {
  crop_summary <- mutate(crop_summary,
                         translated_name=NA_character_, HCAT3_name=NA_character_)
}

crop_summary <- crop_summary |>
  select(crop_code, translated_name, HCAT3_name,
         n_parcels, total_km2, km2_baseline, km2_ssp1, km2_ssp2, km2_ssp5)

write_csv(crop_summary, file.path(OUT_DIR, "agriculture_crop_type_summary.csv"))

# ==============================================================================
# 7. Province x crop type summary
# ==============================================================================

province_crop_summary <- all_parcels |>
  left_join(prov_meta, by="prov_id") |>
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
                                  translated_name=NA_character_, HCAT3_name=NA_character_)
}

province_crop_summary <- province_crop_summary |>
  group_by(province_code, province_name, ccaa_name,
           crop_code, translated_name, HCAT3_name) |>
  summarise(
    n_parcels    = n(),
    total_km2    = round(sum(area_km2,                  na.rm=TRUE), 3),
    km2_baseline = round(sum(area_km2[pa_baseline==1L], na.rm=TRUE), 3),
    km2_ssp1     = round(sum(area_km2[pa_ssp1==1L],     na.rm=TRUE), 3),
    km2_ssp2     = round(sum(area_km2[pa_ssp2==1L],     na.rm=TRUE), 3),
    km2_ssp5     = round(sum(area_km2[pa_ssp5==1L],     na.rm=TRUE), 3),
    .groups = "drop"
  ) |>
  arrange(province_code, desc(total_km2))

write_csv(province_crop_summary,
          file.path(OUT_DIR, "agriculture_province_crop_summary.csv"))

