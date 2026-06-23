# ============================================================
# Add harvest overlap/timing covariates to used-available data
#
# Purpose:
# For each used/available point, check whether it falls within
# harvest buffers of different radii, then calculate whether
# the point occurred before or after the nearest/earliest
# harvest event.
#
# Inputs:
#   dmodel        = sf object containing used + available points
#   harvest       = sf object containing harvest locations/events
#   radii_m       = buffer radii around harvest locations
#   time_col      = harvest datetime column
#   datetime_col  = point datetime column
#   id_col        = unique point ID column created if missing
#
# Returns:
#   sf object with:
#     - GIS_TIME_<radius>m
#     - DiffTime<radius>m
#     - afterharvest<radius>m
#
# Notes:
#   - Uses st_join(), not st_intersection(), so points with no
#     harvest overlap are retained.
#   - For points overlapping multiple harvest buffers, keeps
#     the earliest GIS_TIME.
#   - CRS units must be meters.
# ============================================================

add_harvest_overlap <- function(dmodel,
                                harvest,
                                radii_m = c(200, 50),
                                time_col = "GIS_TIME",
                                datetime_col = "UTC_datetime",
                                id_col = "seq") {

  # ----------------------------------------------------------
  # Ensure point data has a unique ID
  # Do this is the data has no column called seq
  # ----------------------------------------------------------
  
  if (!id_col %in% names(dmodel)) {
    dmodel[[id_col]] <- seq_len(nrow(dmodel))
  }
  
  # ----------------------------------------------------------
  # Check CRS consistency
  # ----------------------------------------------------------
  
  if (sf::st_crs(dmodel) != sf::st_crs(harvest)) {
    cat("Changing crs of harvest data")
    harvest <- sf::st_transform(harvest, sf::st_crs(dmodel))
  }
  
  # ----------------------------------------------------------
  # Work iteratively through each buffer radius
  # ----------------------------------------------------------
  
  out <- dmodel
  
  for (r in radii_m) {
    
    # Column names to create
    gis_col   <- paste0(time_col, "_", r, "m")
    diff_col  <- paste0("DiffTime", r, "m")
    class_col <- paste0("afterharvest", r, "m")
    
    # --------------------------------------------------------
    # Buffer harvest locations
    # --------------------------------------------------------
    
    harvest_buffer <- harvest %>%
      sf::st_buffer(dist = r) %>%
      dplyr::select(dplyr::all_of(time_col), geometry)
    
    # --------------------------------------------------------
    # Spatial join:
    # keep all points, add harvest time where overlap exists
    # --------------------------------------------------------
    
    joined <- sf::st_join(
      out,
      harvest_buffer,
      join = sf::st_intersects,
      #join = sf::st_nearest_feature,
      left = TRUE
    )
    
    # --------------------------------------------------------
    # If a point overlaps several harvest buffers,
    # keep earliest harvest time for that point.
    # Points with no overlap remain NA.
    # --------------------------------------------------------
    
    harvest_time <- joined %>%
      sf::st_drop_geometry() %>%
      dplyr::group_by(.data[[id_col]]) %>%
      dplyr::summarise(
        !!gis_col := if (all(is.na(.data[[time_col]]))) {
          as.POSIXct(NA)
        } else {
          min(.data[[time_col]], na.rm = TRUE)
        },
        .groups = "drop"
      )
    
    # --------------------------------------------------------
    # Join harvest timing back to main point dataset
    # --------------------------------------------------------
    
    out <- out %>%
      dplyr::left_join(
        harvest_time,
        by = id_col
      )
    
    # --------------------------------------------------------
    # Calculate time difference:
    #
    # harvest time - point time
    #
    # Negative = harvest happened before the point timestamp
    # Positive = harvest happened after the point timestamp
    # --------------------------------------------------------
    
    out[[diff_col]] <- out[[gis_col]] - out[[datetime_col]]
    
    # --------------------------------------------------------
    # Classify point relative to harvest timing
    # --------------------------------------------------------
    
    out[[class_col]] <- dplyr::case_when(
      is.na(out[[diff_col]]) ~ "noharvest",
      out[[diff_col]] < 0    ~ "after",
      out[[diff_col]] > 0    ~ "before",
      out[[diff_col]] == 0   ~ "same_time"
    )
    
    out[[class_col]] <- as.factor(out[[class_col]])
  }
  
  return(out)
}