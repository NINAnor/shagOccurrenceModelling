#' Create raster predictors describing harvest overlap and distances
#'
#' For a colony-centered study region this creates:
#'  - binary rasters: harvested50m, harvested200m (1 inside buffer, 0 outside)
#'  - continuous distance rasters: DistNest (distance to colony), DistHarvest (distance to harvest polygons)
#'  - raster_stack: a terra SpatRaster stack containing the layers in the above order
#'  - harvest50_buffer / harvest200_buffer: sf buffers (unioned) used to produce binary layers
#'
#' Notes:
#'  - Input CRS should be projected (meter units). The function transforms `harvest`
#'    to colony CRS for consistency.
#'  - Distance rasters are scaled (mean 0, sd 1) using terra::scale — suitable for modelling.
#'
#' @param harvest sf (point/polygon) features of harvest events
#' @param colony  sf geometry (point) for colony location
#' @param resolution raster resolution (meters); default 200
#' @param buffer_m study region buffer around colony (meters); default 33000
#'
#' @return list with elements:
#'   harvested50m, harvested200m, DistNest, DistHarvest (SpatRaster layers),
#'   raster_stack (SpatRaster stack), harvest50_buffer, harvest200_buffer (sf)
#' @export
create_harvest_dist_rasters <- function(harvest,
                                        colony,
                                        resolution = 200,
                                        buffer_m = 33000) {
  # Dependency checks
  needed <- c("sf", "terra")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))

  # Validate inputs
  stopifnot(inherits(harvest, "sf"), inherits(colony, "sf"))
  if (nrow(colony) < 1) stop("`colony` must contain at least one geometry (point).")

  # Ensure harvest and colony share CRS (use colony CRS)
  harvest <- sf::st_transform(harvest, sf::st_crs(colony))

  # Study region: buffer around colony
  study_region <- sf::st_buffer(colony, dist = buffer_m)
  study_union  <- sf::st_union(study_region)
  study_vect   <- terra::vect(study_union)

  # Template raster covering the study region
  template <- terra::rast(
    extent = terra::ext(study_vect),
    resolution = resolution,
    crs = terra::crs(study_vect)
  )

  # Create unioned buffers around harvest polygons for 50m and 200m
  # Use st_union to collapse to single geometry per buffer to speed rasterize
  harvest50  <- harvest |> sf::st_buffer(dist = 50)  |> sf::st_union() |> sf::st_make_valid()
  harvest200 <- harvest |> sf::st_buffer(dist = 200) |> sf::st_union() |> sf::st_make_valid()

  # Convert to terra vect for raster operations
  harvest50_vect  <- terra::vect(harvest50)
  harvest200_vect <- terra::vect(harvest200)

  # Rasterize buffers: 1 inside buffer, 0 outside. Use touches = TRUE to be conservative.
  r_h50 <- terra::rasterize(harvest50_vect, template, field = 1, background = 0, touches = TRUE)
  r_h200 <- terra::rasterize(harvest200_vect, template, field = 1, background = 0, touches = TRUE)
  names(r_h50)  <- "harvested50m"
  names(r_h200) <- "harvested200m"

  # Distance to colony (meters): compute on template, mask to study region, scale
  colony_vect <- terra::vect(colony)
  r_dist_nest <- terra::distance(template, colony_vect)
  r_dist_nest <- terra::mask(r_dist_nest, study_vect)
  r_dist_nest <- terra::scale(r_dist_nest)
  names(r_dist_nest) <- "DistNest"

  # Distance to harvest polygons (meters): compute to unioned harvest geometry
  harvest_union <- sf::st_make_valid(sf::st_union(harvest))
  harvest_vect  <- terra::vect(harvest_union)
  r_dist_harv <- terra::distance(template, harvest_vect)
  r_dist_harv <- terra::mask(r_dist_harv, study_vect)
  r_dist_harv <- terra::scale(r_dist_harv)
  names(r_dist_harv) <- "DistHarvest"

  # Mask binary harvest rasters to study region for consistency
  r_h50  <- terra::mask(r_h50, study_vect)
  r_h200 <- terra::mask(r_h200, study_vect)

  # Combine into a stack (order important for downstream code)
  raster_stack <- c(r_h50, r_h200, r_dist_nest, r_dist_harv)

  list(
    harvested50m = r_h50,
    harvested200m = r_h200,
    DistNest = r_dist_nest,
    DistHarvest = r_dist_harv,
    raster_stack = raster_stack,
    harvest50_buffer = harvest50,
    harvest200_buffer = harvest200
  )
}

#' Save raster results to disk (parquet + RData metadata)
#'
#' Saves the raster stack as a parquet table (x,y,layer-values) and stores
#' non-raster metadata as an RData file. This keeps raster storage portable
#' (but may be large).
#'
#' @param raster_results list returned by create_harvest_dist_rasters()
#' @param year character or numeric year label used in file names
#' @param path_for_processed_data directory to write files (created if missing)
#' @export
save_raster_results <- function(raster_results,
                                year,
                                path_for_processed_data) {
  if (!requireNamespace("arrow", quietly = TRUE)) stop("Install package 'arrow' to write parquet files.")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install package 'terra'.")

  dir.create(path_for_processed_data, recursive = TRUE, showWarnings = FALSE)

  parquet_file <- file.path(path_for_processed_data, paste0("raster_stack_", year, ".parquet"))
  rdata_file <- file.path(path_for_processed_data, paste0("raster_results_", year, "_metadata.RData"))

  # Convert raster stack to data.frame (xyz)
  raster_df <- terra::as.data.frame(raster_results$raster_stack, xy = TRUE, na.rm = FALSE)

  # Write parquet
  arrow::write_parquet(raster_df, parquet_file)

  # Save metadata (non-raster components) as an RData file
  raster_metadata <- raster_results
  # remove large raster objects before saving metadata
  raster_metadata$harvested50m <- NULL
  raster_metadata$harvested200m <- NULL
  raster_metadata$DistNest <- NULL
  raster_metadata$raster_stack <- NULL

  save(raster_metadata, file = rdata_file)
}

#' Load raster results previously saved with save_raster_results()
#'
#' Reconstructs terra SpatRaster layers from the parquet table and attaches
#' metadata loaded from the RData file. The function requires a reference CRS
#' object (e.g. an sf geometry) to set the raster CRS on read.
#'
#' @param year label used when saving
#' @param path_for_processed_data directory where files were written
#' @param crs_reference an sf object (or anything acceptable to sf::st_crs)
#'                      used to set the CRS of the reconstructed raster stack
#' @return list matching the structure returned by create_harvest_dist_rasters()
#' @export
load_raster_results <- function(year,
                                path_for_processed_data,
                                crs_reference) {
  if (!requireNamespace("arrow", quietly = TRUE)) stop("Install package 'arrow' to read parquet files.")
  if (!requireNamespace("terra", quietly = TRUE)) stop("Install package 'terra'.")

  parquet_file <- file.path(path_for_processed_data, paste0("raster_stack_", year, ".parquet"))
  rdata_file <- file.path(path_for_processed_data, paste0("raster_results_", year, "_metadata.RData"))

  if (!file.exists(parquet_file)) stop("Parquet file not found: ", parquet_file)
  if (!file.exists(rdata_file)) stop("RData metadata file not found: ", rdata_file)

  raster_df <- arrow::read_parquet(parquet_file)

  # Recreate a SpatRaster from an xyz table. Ensure CRS is set from crs_reference.
  crs_wkt <- sf::st_crs(crs_reference)$wkt
  raster_stack <- terra::rast(raster_df, type = "xyz", crs = crs_wkt)

  # Ensure layer names match the original contract
  names(raster_stack) <- c("harvested50m", "harvested200m", "DistNest", "DistHarvest")

  # Load metadata and reattach raster layers
  load(rdata_file) # loads raster_metadata
  raster_results <- raster_metadata
  raster_results$harvested50m <- raster_stack[["harvested50m"]]
  raster_results$harvested200m <- raster_stack[["harvested200m"]]
  raster_results$DistNest <- raster_stack[["DistNest"]]
  raster_results$raster_stack <- raster_stack

  raster_results
}
