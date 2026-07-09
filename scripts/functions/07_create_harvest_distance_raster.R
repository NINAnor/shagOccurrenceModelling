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
  #harvest0 <- harvest |> sf::st_union() |> sf::st_make_valid()
  harvest50  <- harvest |> sf::st_buffer(dist = 50)  |> sf::st_union() |> sf::st_make_valid()
  harvest200 <- harvest |> sf::st_buffer(dist = 200) |> sf::st_union() |> sf::st_make_valid()

  # Convert to terra vect for raster operations
  #harvest0_vect  <- terra::vect(harvest0)
  harvest50_vect  <- terra::vect(harvest50)
  harvest200_vect <- terra::vect(harvest200)

  # Rasterize buffers: 1 inside buffer, 0 outside. Use touches = TRUE to be conservative.
 # r_h0 <- terra::rasterize(harvest0_vect, template, field = 1, background = 0, touches = TRUE)
  r_h50 <- terra::rasterize(harvest50_vect, template, field = 1, background = 0, touches = TRUE)
  r_h200 <- terra::rasterize(harvest200_vect, template, field = 1, background = 0, touches = TRUE)
  names(r_h50)  <- "harvested50m"
  names(r_h200) <- "harvested200m"
  #names(r_h0)  <- "harvested0m"

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

#' Save raster results to disk (parquet OR GeoTIFF + RData metadata)
#'
#' Save either a parquet table (x,y,band-values) or a compressed GeoTIFF of the
#' raster stack, and save non-raster metadata as an RData file. The RData file
#' includes a copy of the original list minus the large raster elements and
#' also stores the original layer names so loads reconstruct correctly.
#'
#' @param raster_results list returned by create_harvest_dist_rasters()
#' @param year character or numeric year label used in file names
#' @param path_for_processed_data directory to write files (created if missing)
#' @param importFile character: either "parquet" or "tif" (default "parquet")
#' @param tif_options list with arguments passed to terra::writeRaster when writing TIF:
#'                    e.g. list(datatype = "FLT4S",
#'                                  gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=YES"))
#' @export
save_raster_results <- function(raster_results,
                                year,
                                path_for_processed_data,
                                importFile = c("parquet", "raster"),
                                tif_options = list(datatype = "FLT4S",
                                                   gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))) {
  importFile <- match.arg(importFile)

  if (!requireNamespace("terra", quietly = TRUE)) stop("Install package 'terra'.")
  if (importFile == "parquet" && !requireNamespace("arrow", quietly = TRUE)) {
    stop("Install package 'arrow' to write parquet files.")
  }

  dir.create(path_for_processed_data, recursive = TRUE, showWarnings = FALSE)

  # filenames
  base_name <- paste0("raster_stack_", year)
  parquet_file <- file.path(path_for_processed_data, paste0(base_name, ".parquet"))
  tif_file     <- file.path(path_for_processed_data, paste0(base_name, ".tif"))
  rdata_file   <- file.path(path_for_processed_data, paste0("raster_results_", year, "_metadata.RData"))

  # Store metadata: copy, but record layer names first
  raster_metadata <- raster_results
  raster_metadata$layer_names <- names(raster_results$raster_stack)

  # Remove heavy raster objects from metadata copy before saving
  raster_metadata$harvested50m <- NULL
  raster_metadata$harvested200m <- NULL
  raster_metadata$DistNest <- NULL
  raster_metadata$DistHarvest <- NULL
  raster_metadata$raster_stack <- NULL

  if (importFile == "parquet") {
    # Convert raster stack to data.frame (xyz)
    message("Converting raster stack to dataframe (may be large) ...")
    raster_df <- terra::as.data.frame(raster_results$raster_stack, xy = TRUE, na.rm = FALSE)

    # Write parquet
    message("Writing parquet to: ", parquet_file)
    arrow::write_parquet(raster_df, parquet_file)

    # Save metadata
    message("Saving metadata to: ", rdata_file)
    save(raster_metadata, file = rdata_file)

    invisible(list(parquet = parquet_file, metadata = rdata_file))
  } else {
    # importFile == "tif": write GeoTIFF directly to the target path (no temp file)
    message("Writing GeoTIFF directly to: ", tif_file)

    # call writeRaster with user-provided tif_options (safe defaults exist)
    write_args <- c(list(x = raster_results$raster_stack,
                         filename = tif_file, overwrite = TRUE),
                    tif_options)
    do.call(terra::writeRaster, write_args)

    # Save metadata (record path to tif in metadata)
    raster_metadata$raster_file <- tif_file
    message("Saved GeoTIFF to: ", tif_file)
    message("Saving metadata to: ", rdata_file)
    save(raster_metadata, file = rdata_file)

    invisible(list(tif = tif_file, metadata = rdata_file))
  }
}


#' Load raster results previously saved with save_raster_results()
#'
#' Reconstructs terra SpatRaster layers from a parquet table or a GeoTIFF and
#' reattaches metadata loaded from the RData file. If parquet is used, a reference
#' CRS (e.g. an sf object) is required; GeoTIFF stores CRS internally.
#'
#' @param year label used when saving
#' @param path_for_processed_data directory where files were written
#' @param crs_reference an sf object (or anything acceptable to sf::st_crs);
#'                      required only when importFile == "parquet"
#' @param importFile character: either "parquet" or "tif" (default "parquet")
#' @return list matching the structure returned by create_harvest_dist_rasters()
#' @export
load_raster_results <- function(year,
                                path_for_processed_data,
                                crs_reference = NULL,
                                importFile = c("parquet", "raster")) {
  importFile <- match.arg(importFile)

  if (!requireNamespace("terra", quietly = TRUE)) stop("Install package 'terra'.")
  if (importFile == "parquet" && !requireNamespace("arrow", quietly = TRUE)) {
    stop("Install package 'arrow' to read parquet files.")
  }

  base_name <- paste0("raster_stack_", year)
  parquet_file <- file.path(path_for_processed_data, paste0(base_name, ".parquet"))
  tif_file     <- file.path(path_for_processed_data, paste0(base_name, ".tif"))
  rdata_file   <- file.path(path_for_processed_data, paste0("raster_results_", year, "_metadata.RData"))

  if (!file.exists(rdata_file)) stop("RData metadata file not found: ", rdata_file)

  load(rdata_file) # loads raster_metadata

  if (importFile == "parquet") {
    if (is.null(crs_reference)) stop("crs_reference is required when loading from parquet.")
    if (!file.exists(parquet_file)) stop("Parquet file not found: ", parquet_file)

    message("Reading parquet: ", parquet_file)
    raster_df <- arrow::read_parquet(parquet_file)

    crs_wkt <- sf::st_crs(crs_reference)$wkt
    message("Reconstructing SpatRaster from xyz table (setting CRS from crs_reference).")
    raster_stack <- terra::rast(raster_df, type = "xyz", crs = crs_wkt)

  } else {
    # importFile == "tif"
    if (!file.exists(tif_file)) stop("TIF file not found: ", tif_file)
    message("Reading GeoTIFF: ", tif_file)
    raster_stack <- terra::rast(tif_file)
  }

  # restore layer names if available in metadata; otherwise keep what's in raster_stack
  if (!is.null(raster_metadata$layer_names)) {
    names(raster_stack) <- raster_metadata$layer_names
  }

  # Reattach raster layers into the metadata structure (if the list expects specific elements)
  raster_results <- raster_metadata
  ln <- names(raster_stack)

  # Attach common expected layers if they exist in raster_stack
  if ("harvested50m" %in% ln) raster_results$harvested50m <- raster_stack[["harvested50m"]]
  if ("harvested200m" %in% ln) raster_results$harvested200m <- raster_stack[["harvested200m"]]
  if ("DistNest" %in% ln) raster_results$DistNest <- raster_stack[["DistNest"]]
  if ("DistHarvest" %in% ln) raster_results$DistHarvest <- raster_stack[["DistHarvest"]]

  # Always include the full stack
  raster_results$raster_stack <- raster_stack

  raster_results
}
