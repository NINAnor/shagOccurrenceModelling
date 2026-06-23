# ============================================================
# Create raster predictors from harvest polygons and colony
#
# Outputs:
#   harvested50m  = 1 if within 50 m of harvest polygon, else 0
#   harvested200m = 1 if within 200 m of harvest polygon, else 0
#   DistNest      = distance to colony in meters
# ============================================================

create_harvest_dist_rasters <- function(harvest,
                                        colony,
                                        resolution = 200,
                                        buffer_m = 33000) {
  

  
  # ----------------------------------------------------------
  # Ensure all objects use same CRS
  # ----------------------------------------------------------
  
  harvest <- sf::st_transform(harvest, sf::st_crs(colony))
  #colony  <- sf::st_transform(colony,  sf::st_crs(study_region))
  study_region <- sf::st_buffer(colony, dist = buffer_m)
  
  # ----------------------------------------------------------
  # Convert study region to terra object
  # ----------------------------------------------------------
  
  study_union <- sf::st_union(study_region)
  study_vect  <- terra::vect(study_union)
  
  # ----------------------------------------------------------
  # Create template raster
  # ----------------------------------------------------------
  
  template <- terra::rast(
    extent = terra::ext(study_vect),
    resolution = resolution,
    crs = terra::crs(study_vect)
  )
  
  # ----------------------------------------------------------
  # Create harvest buffers
  # ----------------------------------------------------------
  
  harvest50 <- harvest |>
    sf::st_buffer(dist = 50) |>
    sf::st_union()
  
  harvest200 <- harvest |>
    sf::st_buffer(dist = 200) |>
    sf::st_union()
  
  # ----------------------------------------------------------
  # Rasterize harvest buffers
  # ----------------------------------------------------------
  
  harvest50_vect <- terra::vect(harvest50)
  harvest200_vect <- terra::vect(harvest200)
  
  r_h50 <- terra::rasterize(
    harvest50_vect,
    template,
    field = 1,
    background = 0,
    touches = FALSE
  )
  
  r_h200 <- terra::rasterize(
    harvest200_vect,
    template,
    field = 1,
    background = 0,
    touches = FALSE
  )
  
  names(r_h50) <- "harvested50m"
  names(r_h200) <- "harvested200m"
  
  # ----------------------------------------------------------
  # Create distance-to-colony raster
  # ----------------------------------------------------------
  
  colony_vect <- terra::vect(colony)
  
  r_dist <- terra::distance(
    template,
    colony_vect
  )
  
  names(r_dist) <- "DistNest"
  
  # ----------------------------------------------------------
  # Mask all rasters to study region
  # ----------------------------------------------------------
  
  r_h50 <- terra::mask(r_h50, study_vect)
  r_h200 <- terra::mask(r_h200, study_vect)
  r_dist <- terra::mask(r_dist, study_vect)
  
  # ----------------------------------------------------------
  # Stack rasters
  # ----------------------------------------------------------
  
  raster_stack <- c(r_h50, r_h200, r_dist)
  
  return(list(
    harvested50m = r_h50,
    harvested200m = r_h200,
    DistNest = r_dist,
    raster_stack = raster_stack,
    harvest50_buffer = harvest50,
    harvest200_buffer = harvest200
  ))
}


## -------------------------------------
## Save and load the results
## -------------------------------------

save_raster_results <- function(raster_results,
                                year,
                                path_for_processed_data) {
  
  # ----------------------------------------------------------
  # File paths
  # ----------------------------------------------------------
  
  parquet_file <- file.path(
    path_for_processed_data,
    paste0("raster_stack_", year, ".parquet")
  )
  
  rdata_file <- file.path(
    path_for_processed_data,
    paste0("raster_results_", year, "_metadata.RData")
  )
  
  
  # ----------------------------------------------------------
  # Save raster stack as parquet
  # ----------------------------------------------------------
  
  raster_df <- terra::as.data.frame(
    raster_results$raster_stack,
    xy = TRUE,
    na.rm = FALSE
  )
  
  arrow::write_parquet(
    raster_df,
    parquet_file
  )
  
  
  # ----------------------------------------------------------
  # Save non-raster objects as RData
  # ----------------------------------------------------------
  
  raster_metadata <- raster_results
  
  raster_metadata$harvested50m <- NULL
  raster_metadata$harvested200m <- NULL
  raster_metadata$DistNest <- NULL
  raster_metadata$raster_stack <- NULL
  
  save(
    raster_metadata,
    file = rdata_file
  )
}

load_raster_results <- function(year,
                                path_for_processed_data,
                                crs_reference) {
  
  parquet_file <- file.path(
    path_for_processed_data,
    paste0("raster_stack_", year, ".parquet")
  )
  
  rdata_file <- file.path(
    path_for_processed_data,
    paste0("raster_results_", year, "_metadata.RData")
  )
  
  raster_df <- arrow::read_parquet(parquet_file)
  
  raster_stack <- terra::rast(
    raster_df,
    type = "xyz",
    crs = sf::st_crs(crs_reference)$wkt
  )
  
  names(raster_stack) <- c(
    "harvested50m",
    "harvested200m",
    "DistNest"
  )
  
  load(rdata_file)
  
  raster_results <- raster_metadata
  
  raster_results$harvested50m <- raster_stack[["harvested50m"]]
  raster_results$harvested200m <- raster_stack[["harvested200m"]]
  raster_results$DistNest <- raster_stack[["DistNest"]]
  raster_results$raster_stack <- raster_stack
  
  return(raster_results)
}