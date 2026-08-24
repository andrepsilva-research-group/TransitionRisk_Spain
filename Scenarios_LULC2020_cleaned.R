setwd("/Users/fishfeed/Library/CloudStorage/OneDrive-StockholmUniversity/FinBio/R_scripts/NatPoKe-dev/NatPoKe-dev_Ines")

library(prioritizr)
library(sf)
library(terra)
library(dplyr)
library(tibble)
library(ggplot2)
library(gurobi)
library(rnaturalearth)
library(patchwork)
library(foreign)
library(geodata)
#setting up ----
set.seed(500)
terraOptions(threads = parallel::detectCores() - 1, memfraction = 0.8)

dir.create("Prioritization/solutions",               showWarnings = FALSE, recursive = TRUE)
dir.create("Prioritization/solutions/ssp2_regional", showWarnings = FALSE, recursive = TRUE)

normalize_name <- function(x) {
  x %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    tolower() %>%
    trimws() %>%
    gsub("\\s+", " ", .)
}
#establishing feature targets by scenario ----
SSP1_EUNIS_TARGET <- 0.45
SSP1_LULC_TARGET  <- 0.45
SSP1_SUIT_TARGET  <- 0.50
SSP1_NCP_TARGET   <- 0.90

SSP2_EUNIS_TARGET <- 0.20
SSP2_LULC_TARGET  <- 0.20
SSP2_SUIT_TARGET  <- 0.25
SSP2_NCP_TARGET   <- 0.30
#loading PU and input layers ---
pt <- ne_countries(scale = "medium", country = "Portugal", returnclass = "sf")
es <- ne_countries(scale = "medium", country = "Spain",    returnclass = "sf")

countries   <- rbind(pt, es)
countries   <- st_transform(countries, 3035)
countries_v <- vect(countries)

spain   <- countries %>% filter(admin == "Spain")
spain_v <- vect(spain)

pu_cost              <- rast("Prioritization/cache/pu_cost_renorm.tif")
PA_rast              <- rast("Prioritization/cache/PA_rast.tif")
area_raster          <- rast("Prioritization/cache/area_raster.tif")
ncp_rast             <- rast("Prioritization/cache/ncp_rast.tif")
artificial_prop_ssp5 <- rast("Prioritization/cache/artificial_prop_ssp5.tif")

PA_rast     <- resample(PA_rast,     pu_cost, method = "near")
area_raster <- resample(area_raster, pu_cost, method = "bilinear")
PA_rast     <- as.logical(PA_rast)

EUNIS_features_prop <- rast("Prioritization/Input_layers/EUNIS_features_prop_1km.tif")
n_eunis_layers      <- nlyr(EUNIS_features_prop)

lulc_current  <- rast("Prioritization/Input_layers/LULC_EUNIS_prop_current_1km.tif")
n_lulc_layers <- nlyr(lulc_current)

suit_layers_norm <- lapply(c("current", "ssp585_2050"), function(sc) {
  r       <- rast(paste0("Prioritization/Input_layers/suit_", sc, "_1km.tif"))
  max_val <- terra::global(mask(r, spain_v), "max", na.rm = TRUE)[1, 1]
  r_norm  <- r / max_val
  names(r_norm) <- paste0("suit_", sc, "_norm")
  r_norm
})
names(suit_layers_norm) <- c("current", "ssp585_2050")

suit_current_norm <- suit_layers_norm[["current"]]

spain_ccaa_all <- geodata::gadm(
  country = "ESP", level = 1,
  path    = "Prioritization/GADM"
) %>%
  st_as_sf() %>%
  st_transform(3035)

cost_spain <- mask(pu_cost,     spain_v)
area_spain <- mask(area_raster, spain_v)
PA_spain   <- mask(PA_rast,     spain_v)
ncp_spain  <- mask(ncp_rast,    spain_v)

total_area <- terra::global(area_spain, "sum", na.rm = TRUE)[1, 1]
PA_area    <- terra::global(area_spain * as.numeric(PA_spain),
                            "sum", na.rm = TRUE)[1, 1]
#SSP1 ----
eunis_spain <- mask(EUNIS_features_prop, spain_v)
lulc_spain  <- mask(lulc_current,        spain_v)
suit_spain  <- mask(suit_current_norm,   spain_v)

features_ssp1 <- c(eunis_spain, lulc_spain, suit_spain, ncp_spain)

target_vec_ssp1 <- c(
  rep(SSP1_EUNIS_TARGET, n_eunis_layers),
  rep(SSP1_LULC_TARGET,  n_lulc_layers),
  SSP1_SUIT_TARGET,
  SSP1_NCP_TARGET
)

budget_ssp1_upper <- 0.50 * total_area
budget_ssp1_lower <- 0.4995 * total_area
budget_ssp1_lower <- max(budget_ssp1_lower,
                         PA_area + 0.01 * (total_area - PA_area))

p_ssp1 <- problem(cost_spain, features_ssp1) %>%
  add_min_set_objective() %>%
  add_relative_targets(target_vec_ssp1) %>%
  add_binary_decisions() %>%
  add_locked_in_constraints(PA_spain) %>%
  add_linear_constraints(threshold = budget_ssp1_upper, sense = "<=", data = area_spain) %>%
  add_linear_constraints(threshold = budget_ssp1_lower, sense = ">=", data = area_spain) %>%

  add_gurobi_solver(gap = 0.01, verbose = TRUE, threads = 8)

sol_ssp1 <- solve(p_ssp1)

writeRaster(sol_ssp1,
            "Prioritization/solutions/solution_ssp1.tif",
            overwrite = TRUE)

p_eval_ssp1 <- problem(cost_spain, features_ssp1) %>%
  add_min_set_objective() %>%
  add_relative_targets(target_vec_ssp1) %>%
  add_binary_decisions()

tc_ssp1 <- eval_target_coverage_summary(p_eval_ssp1, sol_ssp1)

selected_km2_ssp1 <- terra::global(
  area_spain * as.numeric(sol_ssp1 == 1), "sum", na.rm = TRUE
)[1, 1]
#ssp2 ----
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

params_norm_ssp2 <- ssp2_params %>%
  mutate(name_key = normalize_name(NAME_1))

spain_ccaa <- spain_ccaa_all %>%
  mutate(name_key = normalize_name(NAME_1)) %>%
  inner_join(params_norm_ssp2 %>% select(name_key, budget_mult), by = "name_key") %>%
  select(-name_key)

still_missing <- setdiff(
  normalize_name(ssp2_params$NAME_1),
  normalize_name(spain_ccaa$NAME_1)
)

ccaa_baselines_df <- readRDS("Prioritization/cache/ccaa_baselines_poly.rds")

ccaa_baselines_df <- ccaa_baselines_df %>%
  select(-any_of("budget_mult")) %>%
  mutate(name_key = normalize_name(NAME_1)) %>%
  left_join(
    ssp2_params %>%
      mutate(name_key = normalize_name(NAME_1)) %>%
      select(name_key, budget_mult),
    by = "name_key"
  ) %>%
  select(-name_key) %>%
  mutate(
    ssp2_target_pct = pmin(baseline_pct_poly * budget_mult, 95),
    ssp2_target_pct = pmax(ssp2_target_pct, baseline_pct_poly),
    expansion_pp    = ssp2_target_pct - baseline_pct_poly,
    needs_expansion = expansion_pp > 0.5
  )

spain_ccaa <- spain_ccaa %>%
  left_join(
    ccaa_baselines_df %>% select(NAME_1, baseline_pct_poly, ssp2_target_pct,
                                 expansion_pp, needs_expansion),
    by = "NAME_1"
  )

national_ssp2_pct <- sum(ccaa_baselines_df$ssp2_target_pct * ccaa_baselines_df$total_km2) /
  sum(ccaa_baselines_df$total_km2)
national_base_pct <- sum(ccaa_baselines_df$baseline_pct_poly * ccaa_baselines_df$total_km2) /
  sum(ccaa_baselines_df$total_km2)

regional_solution_paths <- c()
solver_log              <- tibble()

for (i in seq_len(nrow(spain_ccaa))) {

  ccaa_row    <- spain_ccaa[i, ]
  ccaa_name   <- ccaa_row$NAME_1
  ssp2_target <- ccaa_row$ssp2_target_pct
  ssp2_frac   <- ssp2_target / 100
  expand_pp   <- ccaa_row$expansion_pp
  run_solver  <- ccaa_row$needs_expansion

  ccaa_v <- vect(ccaa_row)

  cost_r  <- mask(crop(pu_cost,             ccaa_v), ccaa_v)
  eunis_r <- mask(crop(EUNIS_features_prop, ccaa_v), ccaa_v)
  lulc_r  <- mask(crop(lulc_current,        ccaa_v), ccaa_v)
  suit_r  <- mask(crop(suit_current_norm,   ccaa_v), ccaa_v)
  PA_r    <- mask(crop(PA_rast,             ccaa_v), ccaa_v)
  area_r  <- mask(crop(area_raster,         ccaa_v), ccaa_v)
  ncp_r   <- mask(crop(ncp_rast,            ccaa_v), ccaa_v)

  cost_max_r            <- terra::global(cost_r, "max", na.rm = TRUE)[1, 1]
  cost_r[is.na(cost_r)] <- if (is.na(cost_max_r) || cost_max_r == 0) 1 else cost_max_r
  suit_r[is.na(suit_r)] <- 0
  ncp_r[is.na(ncp_r)]   <- 0
  PA_r[is.na(PA_r)]     <- FALSE
  cost_r[PA_r]          <- 0

  total_area_r <- terra::global(area_r,                    "sum", na.rm = TRUE)[1, 1]
  PA_area_r    <- terra::global(area_r * as.numeric(PA_r), "sum", na.rm = TRUE)[1, 1]
  PA_frac_r    <- PA_area_r / total_area_r

  clean_name <- gsub("[^A-Za-z0-9]", "_",
                     iconv(ccaa_name, from = "UTF-8", to = "ASCII//TRANSLIT"))
  out_path <- paste0(
    "Prioritization/solutions/ssp2_regional/solution_ssp2_", clean_name, ".tif"
  )

  if (!run_solver) {

    sol_ccaa              <- as.int(PA_r)
    sol_ccaa[is.na(sol_ccaa)] <- 0L
    sel_km2_r <- terra::global(area_r * as.numeric(sol_ccaa == 1),
                               "sum", na.rm = TRUE)[1, 1]

    solver_log <- bind_rows(solver_log, tibble(
      NAME_1       = ccaa_name,
      total_km2    = total_area_r,
      PA_km2       = PA_area_r,
      target_pct   = ssp2_target,
      selected_km2 = sel_km2_r,
      selected_pct = 100 * sel_km2_r / total_area_r,
      status       = "PA_ONLY"
    ))
    writeRaster(sol_ccaa, out_path, overwrite = TRUE)
    regional_solution_paths <- c(regional_solution_paths, out_path)
    rm(cost_r, eunis_r, lulc_r, suit_r, PA_r, area_r, ncp_r, sol_ccaa)
    gc(); next
  }

  features_r <- c(eunis_r, lulc_r, suit_r, ncp_r)
  layer_sums <- terra::global(features_r, "sum", na.rm = TRUE)$sum
  keep_idx   <- which(layer_sums > 0)

  if (length(keep_idx) == 0) {

    sol_ccaa              <- as.int(PA_r)
    sol_ccaa[is.na(sol_ccaa)] <- 0L
    sel_km2_r <- terra::global(area_r * as.numeric(sol_ccaa == 1),
                               "sum", na.rm = TRUE)[1, 1]
    solver_log <- bind_rows(solver_log, tibble(
      NAME_1 = ccaa_name, total_km2 = total_area_r, PA_km2 = PA_area_r,
      target_pct = ssp2_target, selected_km2 = sel_km2_r,
      selected_pct = 100 * sel_km2_r / total_area_r,
      status = "FALLBACK: all features zero"
    ))
    writeRaster(sol_ccaa, out_path, overwrite = TRUE)
    regional_solution_paths <- c(regional_solution_paths, out_path)
    rm(cost_r, eunis_r, lulc_r, suit_r, PA_r, area_r, ncp_r, features_r, sol_ccaa)
    gc(); next
  }

  features_r   <- features_r[[keep_idx]]
  n_free_cells <- sum(!values(PA_r) & !is.na(values(cost_r)), na.rm = TRUE)
  budget_headroom <- ssp2_frac - PA_frac_r

  budget_upper_r <- ssp2_frac * total_area_r
  budget_lower_r <- max((ssp2_frac - 0.0005) * total_area_r, PA_area_r)

  kept_names <- names(features_r)
  target_vec <- ifelse(
    kept_names == "ncp_priority",      SSP2_NCP_TARGET,
    ifelse(grepl("suit_", kept_names), SSP2_SUIT_TARGET,
           SSP2_EUNIS_TARGET)
  )

  if (length(target_vec) != nlyr(features_r)) {

    solver_log <- bind_rows(solver_log, tibble(
      NAME_1 = ccaa_name, total_km2 = total_area_r, PA_km2 = PA_area_r,
      target_pct = ssp2_target, selected_km2 = PA_area_r,
      selected_pct = 100 * PA_frac_r,
      status = "SKIPPED: target length mismatch"
    ))
    next
  }

  p_base <- problem(cost_r, features_r) %>%
    add_min_set_objective() %>%
    add_relative_targets(target_vec) %>%
    add_binary_decisions() %>%
    add_linear_constraints(threshold = budget_upper_r, sense = "<=", data = area_r) %>%
    add_linear_constraints(threshold = budget_lower_r, sense = ">=", data = area_r) %>%

    add_gurobi_solver(gap = 0.01, verbose = FALSE)

  p_ccaa <- if (PA_area_r > 0) {
    p_base %>% add_locked_in_constraints(PA_r)
  } else {
    p_base
  }

  solve_result <- tryCatch({
    chk <- presolve_check(p_ccaa)
    if (!isTRUE(chk)) {

      list(solution = solve(p_ccaa, force = TRUE), status = "OK_FORCED")
    } else {
      list(solution = solve(p_ccaa), status = "OK")
    }
  }, error = function(e) {

    fb <- as.int(PA_r); fb[is.na(fb)] <- 0L
    list(solution = fb, status = paste("FALLBACK:", conditionMessage(e)))
  })

  sol_ccaa   <- solve_result$solution
  sol_status <- solve_result$status

  sel_km2_r <- terra::global(
    area_r * as.numeric(sol_ccaa == 1), "sum", na.rm = TRUE
  )[1, 1]

  solver_log <- bind_rows(solver_log, tibble(
    NAME_1       = ccaa_name,
    total_km2    = total_area_r,
    PA_km2       = PA_area_r,
    target_pct   = ssp2_target,
    selected_km2 = sel_km2_r,
    selected_pct = 100 * sel_km2_r / total_area_r,
    status       = sol_status
  ))

  writeRaster(sol_ccaa, out_path, overwrite = TRUE)
  regional_solution_paths <- c(regional_solution_paths, out_path)

  rm(cost_r, eunis_r, lulc_r, suit_r, PA_r, area_r, ncp_r,
     features_r, sol_ccaa, p_base, p_ccaa)
  gc()
}

write.csv(solver_log,
          "Prioritization/solutions/ssp2_solver_log.csv",
          row.names = FALSE)

n_ok     <- sum(solver_log$status %in% c("OK", "OK_FORCED"))
n_paonly <- sum(solver_log$status == "PA_ONLY")
n_failed <- sum(grepl("FALLBACK|SKIPPED", solver_log$status))

national_selected_km2 <- sum(solver_log$selected_km2)
national_total_km2    <- sum(solver_log$total_km2)

if (length(regional_solution_paths) == 0) {
  stop("No regional solution files written — check errors above.")
}

sol_rasters <- lapply(regional_solution_paths, function(p) {
  r <- rast(p); r[is.na(r)] <- 0L; r
})

national_ssp2_mosaic <- mosaic(sprc(sol_rasters), fun = "max")
sol_ssp2             <- resample(national_ssp2_mosaic, cost_spain, method = "near")
sol_ssp2[is.na(cost_spain)]                    <- NA
sol_ssp2[is.na(sol_ssp2) & !is.na(cost_spain)] <- 0L
sol_ssp2             <- as.int(sol_ssp2)
names(sol_ssp2)      <- "solution_ssp2_regional"

writeRaster(sol_ssp2,
            "Prioritization/solutions/solution_ssp2_regional.tif",
            overwrite = TRUE)

selected_km2_ssp2 <- terra::global(
  area_spain * as.numeric(sol_ssp2 == 1), "sum", na.rm = TRUE
)[1, 1]

ssp2_expanded <- solver_log %>%
  left_join(
    ccaa_baselines_df %>% select(NAME_1, baseline_pct_poly),
    by = "NAME_1"
  ) %>%
  mutate(
    expansion_pp         = selected_pct - baseline_pct_poly,
    expansion_pct_growth = ifelse(PA_km2 > 0,
                                  100 * (selected_km2 - PA_km2) / PA_km2,
                                  NA_real_)
  )

write.csv(ssp2_expanded,
          "Prioritization/solutions/ssp2_ccaa_expansion.csv",
          row.names = FALSE)
#SSP5 ----
PA_LOSS_THRESHOLD <- 0.50

converted_mask <- resample(
  artificial_prop_ssp5 > PA_LOSS_THRESHOLD,
  PA_spain,
  method = "near"
)
converted_mask <- mask(converted_mask, spain_v)

pa_lost <- ifel(PA_spain &  converted_mask, TRUE, FALSE)
pa_surv <- ifel(PA_spain & !converted_mask, TRUE, FALSE)

PA_km2_total <- terra::global(
  area_spain * as.numeric(PA_spain), "sum", na.rm = TRUE)[1, 1]
lost_km2 <- terra::global(
  area_spain * as.numeric(pa_lost),  "sum", na.rm = TRUE)[1, 1]
surv_km2 <- terra::global(
  area_spain * as.numeric(pa_surv),  "sum", na.rm = TRUE)[1, 1]

sol_ssp5_dev <- ifel(pa_lost, 2L, ifel(pa_surv, 1L, 0L))

spain_land_mask  <- rasterize(spain_v, sol_ssp5_dev, field = 1L)
sol_ssp5_dev     <- ifel(!is.na(spain_land_mask) & is.na(sol_ssp5_dev), 0L, sol_ssp5_dev)
sol_ssp5_dev[is.na(spain_land_mask)] <- NA
sol_ssp5_dev     <- as.int(sol_ssp5_dev)
names(sol_ssp5_dev) <- "solution_ssp5_development"

writeRaster(
  sol_ssp5_dev,
  "Prioritization/solutions/solution_ssp5_development.tif",
  overwrite = TRUE
)

ccaa_loss_log <- lapply(seq_len(nrow(spain_ccaa)), function(i) {
  ccaa_v_i <- vect(spain_ccaa[i, ])
  name_i   <- spain_ccaa$NAME_1[i]
  area_i   <- mask(crop(area_spain, ccaa_v_i), ccaa_v_i)
  PA_i     <- mask(crop(PA_spain,   ccaa_v_i), ccaa_v_i)
  lost_i   <- mask(crop(pa_lost,    ccaa_v_i), ccaa_v_i)
  surv_i   <- mask(crop(pa_surv,    ccaa_v_i), ccaa_v_i)
  total_i  <- terra::global(area_i,                       "sum", na.rm = TRUE)[1, 1]
  pa_i     <- terra::global(area_i * as.numeric(PA_i),    "sum", na.rm = TRUE)[1, 1]
  lost_i_  <- terra::global(area_i * as.numeric(lost_i),  "sum", na.rm = TRUE)[1, 1]
  surv_i_  <- terra::global(area_i * as.numeric(surv_i),  "sum", na.rm = TRUE)[1, 1]
  tibble(
    NAME_1               = name_i,
    total_km2            = total_i,
    pa_current_km2       = pa_i,
    pa_current_pct       = 100 * pa_i   / total_i,
    pa_lost_km2          = lost_i_,
    pa_lost_pct_of_pa    = if (pa_i > 0) 100 * lost_i_ / pa_i else NA_real_,
    pa_surv_km2          = surv_i_,
    pa_surv_pct_of_spain = 100 * surv_i_ / total_i
  )
})

ccaa_loss_df <- bind_rows(ccaa_loss_log)

write.csv(ccaa_loss_df,
          "Prioritization/solutions/ssp5_ccaa_loss_breakdown.csv",
          row.names = FALSE)
#summary table ----
summary_tbl <- tibble(
  Scenario     = c("SSP1 (ambitious)",
                   "SSP2 (regional rivalry)",
                   "SSP5 (fossil-fuelled development)"),
  Target_pct   = c("50%", "~30% national", "PA erosion by LULC change"),
  Selected_km2 = c(selected_km2_ssp1, selected_km2_ssp2, surv_km2),
  Selected_pct = round(c(selected_km2_ssp1, selected_km2_ssp2, surv_km2) /
                         total_area * 100, 1),
  EUNIS_LULC_target = c(SSP1_EUNIS_TARGET, SSP2_EUNIS_TARGET, NA),
  NCP_target         = c(SSP1_NCP_TARGET,   SSP2_NCP_TARGET,   NA),
  PA_lost_km2        = c(NA, NA, lost_km2),
  PA_lost_pct        = c(NA, NA, round(100 * lost_km2 / PA_km2_total, 1))
)

write.csv(summary_tbl,
          "Prioritization/solutions/scenario_summary.csv",
          row.names = FALSE)

ccaa_ssp1_log <- lapply(seq_len(nrow(spain_ccaa_all)), function(i) {
  ccaa_v_i <- vect(spain_ccaa_all[i, ])
  name_i   <- spain_ccaa_all$NAME_1[i]
  area_i   <- mask(crop(area_spain, ccaa_v_i), ccaa_v_i)
  PA_i     <- mask(crop(PA_spain,   ccaa_v_i), ccaa_v_i)
  sol_i    <- mask(crop(sol_ssp1,   ccaa_v_i), ccaa_v_i)
  total_i  <- terra::global(area_i,                        "sum", na.rm = TRUE)[1, 1]
  pa_i     <- terra::global(area_i * as.numeric(PA_i),     "sum", na.rm = TRUE)[1, 1]
  sel_i    <- terra::global(area_i * as.numeric(sol_i == 1),"sum", na.rm = TRUE)[1, 1]
  tibble(
    NAME_1           = name_i,
    total_km2        = total_i,
    pa_baseline_km2  = pa_i,
    pa_baseline_pct  = 100 * pa_i  / total_i,
    selected_km2     = sel_i,
    selected_pct     = 100 * sel_i / total_i,
    expansion_km2    = sel_i - pa_i,
    expansion_pp     = 100 * (sel_i - pa_i) / total_i,
    expansion_pct_growth = ifelse(pa_i > 0, 100 * (sel_i - pa_i) / pa_i, NA_real_)
  )
})

ccaa_ssp1_df <- bind_rows(ccaa_ssp1_log)
write.csv(ccaa_ssp1_df,
          "Prioritization/solutions/ssp1_ccaa_coverage.csv",
          row.names = FALSE)
#plotting ----
common_theme <- theme_minimal() +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 12, color = "grey40"),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank()
  )

common_coord <- coord_sf(
  xlim = c(-9.5, 4.5), ylim = c(35.5, 44), expand = FALSE
)

countries_4326 <- st_transform(countries, 4326)

rast_to_sf <- function(r) {
  poly           <- as.polygons(r, dissolve = FALSE) %>%
    st_as_sf() %>%
    st_transform(4326)
  names(poly)[1] <- "selected"
  poly
}

PA_sf <- PA_rast %>%
  mask(spain_v) %>%
  (\(x) { x[x != 1] <- NA; x })() %>%
  as.polygons(dissolve = TRUE) %>%
  st_as_sf() %>%
  st_transform(4326)

make_panel <- function(sol_r, title, subtitle = NULL) {
  sol_sf <- rast_to_sf(sol_r)
  ggplot() +
    geom_sf(data = sol_sf,         aes(fill = factor(selected)), color = NA) +
    geom_sf(data = countries_4326, fill = NA, color = "black", linewidth = 0.3) +
    geom_sf(data = PA_sf,          fill = "#228B22", color = "#228B22",
            linewidth = 0.1) +
    scale_fill_manual(
      values = c("0" = "grey90", "1" = "red"),
      labels = c("Not selected", "Selected"),
      name   = NULL
    ) +
    labs(title = title, subtitle = subtitle) +
    common_theme +
    common_coord
}

panel_ssp1 <- make_panel(
  sol_ssp1,
  "SSP1 — ambitious protection")

panel_ssp2 <- make_panel(
  sol_ssp2,
  "SSP2 — regional rivalry",
)

pa_lost_sf <- ifel(mask(pa_lost, spain_v), 1L, NA) %>%
  as.polygons(dissolve = TRUE) %>%
  st_as_sf() %>%
  st_transform(4326)

ssp5_base_sf          <- rast_to_sf(sol_ssp5_dev)
ssp5_base_sf$fill_cat <- "Not selected"

panel_ssp5 <- ggplot() +
  geom_sf(data = ssp5_base_sf,   aes(fill = fill_cat), color = NA) +
  geom_sf(data = countries_4326, fill = NA, color = "black", linewidth = 0.3) +
  geom_sf(data = PA_sf,          fill = "#228B22", color = "#228B22", linewidth = 0.1) +
  geom_sf(data = pa_lost_sf,     aes(fill = "PA lost to development"), color = NA) +
  scale_fill_manual(
    values = c(
      "Not selected"           = "grey90",
      "PA lost to development" = "#CC4A1A"
    ),
    breaks = c("Not selected", "PA lost to development"),
    name   = NULL
  ) +
  labs(
    title    = "SSP5 — fossil-fuelled development",
  ) +
  common_theme +
  common_coord

combined_plot <- (panel_ssp1 | panel_ssp2 | panel_ssp5) +
  plot_annotation(
    title    = "Conservation prioritization scenarios — Spain",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8, color = "grey40")
    )
  )

ggsave(
  "Prioritization/solutions/scenario_comparison_three.png",
  combined_plot,
  width = 18, height = 7, dpi = 300
)
