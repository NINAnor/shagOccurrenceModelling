# Create or load precomputed harvest-distance raster results for 2023 and 2024
#
# Behavior:
#  - If the parquet file for a given year (raster_stack_<year>.parquet) does not
#    exist in path_for_processed_data, this script builds the rasters using
#    create_harvest_dist_rasters(...) and saves them with save_raster_results(...).
#  - Otherwise it loads previously saved results via load_raster_results(...).
#
# Required objects in the calling environment:
#  - path_for_processed_data : output directory (string)
#  - colony                  : sf geometry used as colony reference CRS
#  - SK_35km_22_23           : harvest sf for 2022/2023
#  - SK_35km_23_24           : harvest sf for 2023/2024
#
# Required helper functions:
#  - create_harvest_dist_rasters()
#  - save_raster_results()
#  - load_raster_results()
#
# The resulting lists raster_results_23 and raster_results_24 are created/loaded in
# the current environment and contain elements such as:
#   $raster_stack, $harvested50m, $harvested200m, $DistNest, $DistHarvest, ...
# ------------------------------------------------------------------------------

# Basic sanity checks for environment objects and helper functions
if (!exists("path_for_processed_data")) stop("path_for_processed_data must be defined (output directory).")
if (!dir.exists(path_for_processed_data)) dir.create(path_for_processed_data, recursive = TRUE)

required_objs <- c("colony", "SK_35km_22_23", "SK_35km_23_24")
missing_objs <- required_objs[!vapply(required_objs, exists, logical(1))]
if (length(missing_objs) > 0) stop("Missing required objects in environment: ", paste(missing_objs, collapse = ", "))

required_fns <- c("create_harvest_dist_rasters", "save_raster_results", "load_raster_results")
missing_fns <- required_fns[!vapply(required_fns, exists, logical(1))]
if (length(missing_fns) > 0) stop("Missing required helper functions: ", paste(missing_fns, collapse = ", "))

# Helper to build/load for one year
build_or_load_year <- function(year_label, harvest_obj) {
  parquet_path <- file.path(path_for_processed_data, paste0("raster_stack_", year_label, ".parquet"))

  if (!file.exists(parquet_path)) {
    message("Raster parquet for ", year_label, " not found. Creating rasters...")

    # Create rasters (uses colony CRS and study buffer inside the function)
    rr <- create_harvest_dist_rasters(
      harvest = harvest_obj,
      colony = colony,
      resolution = 10
    )

    message("Saving raster results for ", year_label, " to ", path_for_processed_data)
    save_raster_results(
      raster_results = rr,
      year = as.character(year_label),
      path_for_processed_data = path_for_processed_data
    )

    # Return the created result
    return(rr)
  } else {
    message("Found existing parquet for ", year_label, " — loading results from disk.")
    rr <- load_raster_results(
      year = as.character(year_label),
      path_for_processed_data = path_for_processed_data,
      crs_reference = colony
    )
    return(rr)
  }
}

# Build or load both years
raster_results_23 <- build_or_load_year("2023", SK_35km_22_23)
raster_results_24 <- build_or_load_year("2024", SK_35km_23_24)

# Summary messages
message("Raster results ready:")
message("  2023 layers: ", paste(names(raster_results_23$raster_stack), collapse = ", "))
message("  2024 layers: ", paste(names(raster_results_24$raster_stack), collapse = ", "))

# Optionally inspect a layer (uncomment to plot)
# plot(raster_results_23$raster_stack)
# plot(raster_results_24$raster_stack)
