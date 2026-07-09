# ============================================================
# Build used-available modelling dataset and environmental rasters
# ============================================================
#
# Purpose:
#  - Load spatial and tracking inputs (colony, kelp polygons, dive locations, harvest)
#  - For each year, create used-available datasets via the pipeline helpers:
#      process_dloctare(), generate_available_points(), add_harvest_overlap(), make_year_dmodel()
#  - Combine years, finalise the modelling dataset
#  - Build 10 m environmental rasters (bathymetry, substrate, kelp mask)
#  - Save processed datasets and raster stacks to path_for_processed_data
#
# Assumptions:
#  - The following helper functions are available in the environment:
#      process_dloctare, generate_available_points, add_harvest_overlap,
#      make_year_dmodel, finalise_dmodel, create_10m_raster
#  - The .RData files referenced below (colony, tare, dloc23, dloc24, SK_35km_22_23, SK_35km_23_24)
#    exist at the indicated relative paths and when loaded create objects with the names used here.
#  - path_for_processed_data is defined and points to a writable directory.
#  - Spatial CRS: project uses a projected CRS with meter units (e.g. UTM zone 32N).
# ============================================================

# ----------------------------
# Required packages (simple checks)
# ----------------------------
required_pkgs <- c("sf", "terra", "dplyr", "arrow", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Install required packages before running: ", paste(missing_pkgs, collapse = ", "))
}

# ----------------------------
# 1) Load core spatial/tracking data
#    (these files are expected to create objects with the names below)
# ----------------------------
message("Loading base data (expect pre-saved R objects).")

# Colony location (object: colony)
load("../colony")         # defines 'colony' sf object
study_region <- sf::st_buffer(colony, dist = 33000)  # 33 km study buffer

# Kelp forest polygons (object: tare)
load("../tare")           # defines 'tare' sf object

# Dive locations: 2023 and 2024
load("../dloc23")         # defines 'dloc23' (2023 dives)
load("../dloc24")         # defines 'dloc'  (2024 dives; named 'dloc' after loading)

# Harvest information (two datasets covering different year windows)
load("../SK_35_22_23")    # defines 'SK_35km_22_23'
load("../SK_35_23_24")    # defines 'SK_35km_23_24'

# ----------------------------
# Check that output path exists or will be created
# ----------------------------
if (!exists("path_for_processed_data")) stop("path_for_processed_data must be defined (output directory).")
dir.create(path_for_processed_data, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Only re-run processing if final dmodel file does not already exist
# ----------------------------
dmodel_rdata_path <- file.path(path_for_processed_data, "dmodel.RData")
if (!file.exists(dmodel_rdata_path)) {

  message("Processing year-specific datasets (this may take some time).")

  # ------------------------------------------------------------
  # 2) Process 2023 dive locations: clip to kelp, generate available points, add harvest covariates
  # ------------------------------------------------------------
  result23 <- make_year_dmodel(
    dloc_year   = dloc23,
    tare        = tare,
    colony      = colony,
    harvest     = SK_35km_22_23,
    year        = "2023",
    n_per_point = n_per_point,        # assumed defined earlier in pipeline
    min_dist_m  = 1000,
    ID_buffer_m = ID_buffer_m,        # assumed defined earlier
    group_by_ID = group_by_ID,
    harvest_buffer = harvest_buffer,  # assumed defined earlier
    thin_dist_km = thin_dist_km,
    min_frac_in_bg = min_frac_in_bg,
    thin_occ = thin_occ,
    harvest_stratify = harvest_stratify,
    weight_background_by_area = weight_background_by_area
  )
  dmodel123 <- result23$dmodel

  # ------------------------------------------------------------
  # 3) Process 2024 dive locations
  # ------------------------------------------------------------
  result24 <- make_year_dmodel(
    dloc_year   = dloc,
    tare        = tare,
    colony      = colony,
    harvest     = SK_35km_23_24,
    year        = "2024",
    n_per_point = n_per_point,
    min_dist_m  = 1000,
    ID_buffer_m = ID_buffer_m,        # assumed defined earlier
    group_by_ID = group_by_ID,
    harvest_buffer = harvest_buffer,  # assumed defined earlier
    thin_dist_km = thin_dist_km,
    min_frac_in_bg = min_frac_in_bg,
    thin_occ = thin_occ,
    harvest_stratify = harvest_stratify,
    weight_background_by_area = weight_background_by_area
  )
  dmodel124 <- result24$dmodel

  # ------------------------------------------------------------
  # 4) Quick diagnostics (plots returned by helper functions)
  #    Users can display or save these as needed.
  # ------------------------------------------------------------
  tryCatch({
    print(result23$dloctare$plot)                # dive locations clipped to kelp (2023)
    print(result23$avail_results$plot)           # used vs available (2023)
    print(result24$dloctare$plot)                # (2024)
    print(result24$avail_results$plot)           # (2024)
  }, error = function(e) {
    message("One or more diagnostic plots failed to render: ", conditionMessage(e))
  })

  # ------------------------------------------------------------
  # 5) Combine years and finalise modelling dataframe
  # ------------------------------------------------------------
  dmodel <- finalise_dmodel(
    dmodel23     = dmodel123,
    dmodel24     = dmodel124,
    datetime_col = "UTC_datetime"
  )

  # Save final combined dataset for later runs
  save(dmodel, file = dmodel_rdata_path)
  message("Saved combined dmodel to: ", dmodel_rdata_path)

} else {
  message("Found existing dmodel file: loading ", dmodel_rdata_path)
  load(dmodel_rdata_path)  # loads 'dmodel'
}

# ----------------------------
# Optional checks if 'checks' flag enabled earlier in script
# ----------------------------
if (exists("checks") && isTRUE(checks)) {
  message("Summary of final dmodel:")
  print(summary(dmodel))
  print(str(dmodel))
  print(table(dmodel$Year, dmodel$random))
  print(table(dmodel$Year, dmodel$harvested50m))
  print(table(dmodel$Year, dmodel$harvested50m, dmodel$random))
}

# ============================================================
# Environmental layers: create 10 m rasters and combine into env_stack_10m
# Save to parquet for faster reload
# ============================================================
if(importFile == "parquet"){
env_parquet_path <- file.path(path_for_processed_data, "env_stack_10m.parquet")
} else if(importFile == "raster"){
  env_parquet_path <- file.path(path_for_processed_data, "env_stack_10m.tif")
}

if (!file.exists(env_parquet_path)) {

  study_region <- sf::st_buffer(colony, dist = 33000 + 10000)  # 33 km study buffer

  message("Creating 10 m environmental rasters (bathymetry, substrates, tareskog).")

  # 1) Load GEBCO bathymetry NetCDF and rasterise/reproject to 10 m template
  gebco_nc <- "~/Mounts/R/Prosjekter/Toppskarv/Analyses_Nina_multiple_Sites/Environmental Covariates/gebco_2019_n85.0_s50.0_w0.0_e40.0.nc"
  if (!file.exists(gebco_nc)) stop("GEBCO NetCDF not found at: ", gebco_nc)

  gebco <- terra::rast(gebco_nc)

  gebco_10m <- create_10m_raster(
    input_data = gebco,
    study_region = study_region,
    resolution = 10,
    raster_name = "bathymetry"
  ) |> terra::scale()
  # plot(gebco_10m)  # optional visual check

  # 2) Hard-substrate shapefile -> binary raster
  hard_substr_shp <- "~/Mounts/R/Prosjekter/Toppskarv/Analyses_Nina_multiple_Sites/Environmental Covariates/Substrate_Rock_hard_substr.shp"
  hard_substr <- sf::st_read(hard_substr_shp, quiet = TRUE)

  hard_substr_10m <- create_10m_raster(
    input_data = hard_substr,
    study_region = study_region,
    resolution = 10,
    field = 1,
    background = 0,
    raster_name = "hard_substrate"
  )

  # 3) Sediment substrate shapefile -> binary raster
  sediment_shp <- "~/Mounts/R/Prosjekter/Toppskarv/Analyses_Nina_multiple_Sites/Environmental Covariates/Substrate_Sediment.shp"
  sediment_substr <- sf::st_read(sediment_shp, quiet = TRUE)

  sediment_substr_10m <- create_10m_raster(
    input_data = sediment_substr,
    study_region = study_region,
    resolution = 10,
    field = 1,
    background = 0,
    raster_name = "sediment_substrate"
  )

  # 4) Kelp forest mask shapefile -> binary raster (tareskog)
  tareskog_shp <- "~/Mounts/R/Prosjekter/Toppskarv/Analyses_Nina_multiple_Sites/Environmental Covariates/Tareskog_merge.shp"
  tareskog <- sf::st_read(tareskog_shp, quiet = TRUE)

  tareskog_10m <- create_10m_raster(
    input_data = tareskog,
    study_region = study_region,
    resolution = 10,
    field = 1,
    background = 0,
    raster_name = "tareskog"
  )

  # 5) Stack layers and name consistently
  env_stack_10m <- c(gebco_10m, hard_substr_10m, sediment_substr_10m, tareskog_10m)
  names(env_stack_10m) <- c("bathymetry", "hard_substrate", "sediment_substrate", "tareskog")

  message("Environmental stack built: ", paste(names(env_stack_10m), collapse = ", "))

  terra::writeRaster(
    x = env_stack_10m,
    filename = env_parquet_path,
    overwrite = TRUE,
    datatype = "FLT4S",
    gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES")
  )
  # Optional quick plots for sanity
  # plot(env_stack_10m)

  # 6) Convert raster stack to dataframe (xyz) and save as parquet for quick re-load
  env_stack_df <- terra::as.data.frame(env_stack_10m, xy = TRUE, na.rm = FALSE)
  arrow::write_parquet(env_stack_df, env_parquet_path)
  message("Saved env_stack_10m parquet to: ", env_parquet_path)

} else {
  if(get_env_covs && importFile == "parquet"){
  message("Loading existing env_stack_10m parquet: ", env_parquet_path)
  env_stack_df <- arrow::read_parquet(env_parquet_path)

  env_stack_10m <- terra::rast(env_stack_df, type = "xyz", crs = sf::st_crs(study_region)$wkt)
  names(env_stack_10m) <- c("bathymetry", "hard_substrate", "sediment_substrate", "tareskog")
  }

  if(get_env_covs && importFile == "raster"){
    env_stack_10m <- terra::rast(env_parquet_path, crs = sf::st_crs(study_region)$wkt)
  }
  # Optional: plot(env_stack_10m)
}

# ============================================================
# End of data preparation
# Outputs available in current environment:
#   dmodel, env_stack_10m, env_stack_df
# ============================================================
