#' Add harvest overlap and timing covariates to used-available points
#'
#' For each point in `dmodel` this function:
#'  - determines whether the point overlaps buffers around harvest events
#'    at one or more radii (radii_m),
#'  - records the earliest harvest timestamp for overlapping buffers,
#'  - computes the time difference (harvest_time - point_time),
#'  - classifies the point as "before", "after", "same_time" or "noharvest".
#'
#' Behavior notes:
#'  - Uses sf spatial joins (st_join with st_intersects) so points with no
#'    overlapping harvest remain in the output (with NA harvest times).
#'  - If a point overlaps multiple harvest events within the buffer, the
#'    earliest harvest time is used.
#'  - CRS units must be in meters for the buffer distances to be meaningful.
#'
#' Inputs:
#'  - dmodel        : sf object of points (used + available). Must contain a datetime column.
#'  - harvest       : sf object of harvest events (point or polygon) containing a timestamp column.
#'  - radii_m       : numeric vector of buffer radii in meters (default c(200, 50)).
#'  - time_col      : name (string) of timestamp column in `harvest` (default "GIS_TIME").
#'  - datetime_col  : name (string) of timestamp column in `dmodel` (default "UTC_datetime").
#'  - id_col        : name (string) of an identifier column in `dmodel` used to aggregate multiple joins
#'                    (created if missing; default "seq").
#'
#' Returns:
#'  - An sf object (same as dmodel) with additional columns for each radius:
#'     * <time_col>_<r>m   : earliest harvest timestamp within r meters (or NA)
#'     * DiffTime<r>m      : time difference (harvest_time - point_time) (difftime)
#'     * afterharvest<r>m  : factor with levels "noharvest", "after", "before", "same_time"
#'
#' @export
add_harvest_overlap <- function(dmodel,
                                harvest,
                                radii_m = c(200, 50),
                                time_col = "GIS_TIME",
                                datetime_col = "UTC_datetime",
                                id_col = "seq") {

  # ---------------------
  # Namespace checks
  # ---------------------
  needed <- c("sf", "dplyr")
  missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Install missing packages: ", paste(missing_pkgs, collapse = ", "))

  `%>%` <- dplyr::`%>%`

  # ---------------------
  # Ensure unique identifier exists in dmodel
  # ---------------------
  if (!id_col %in% names(dmodel)) {
    dmodel[[id_col]] <- seq_len(nrow(dmodel))
  }

  # ---------------------
  # Ensure datetime columns are present and sensible
  # ---------------------
  if (!datetime_col %in% names(dmodel)) {
    stop("datetime_col not found in dmodel: ", datetime_col)
  }
  if (!time_col %in% names(harvest)) {
    stop("time_col not found in harvest: ", time_col)
  }

  # Convert time columns to POSIXct if they are not already
  if (!inherits(dmodel[[datetime_col]], "POSIXt")) {
    dmodel[[datetime_col]] <- as.POSIXct(dmodel[[datetime_col]])
  }
  if (!inherits(harvest[[time_col]], "POSIXt")) {
    harvest[[time_col]] <- as.POSIXct(harvest[[time_col]])
  }

  # ---------------------
  # Harmonize CRS: transform harvest to dmodel CRS if needed
  # ---------------------
  if (!sf::st_crs(dmodel) == sf::st_crs(harvest)) {
    harvest <- sf::st_transform(harvest, sf::st_crs(dmodel))
  }

  # ---------------------
  # Iteratively build columns for each radius
  # ---------------------
  out <- dmodel

  for (r in radii_m) {
    # Construct names for new columns
    gis_col   <- paste0(time_col, "_", r, "m")   # earliest harvest timestamp within r
    diff_col  <- paste0("DiffTime", r, "m")      # harvest_time - point_time
    class_col <- paste0("afterharvest", r, "m")  # factor classification

    # Buffer harvest events
    harvest_buffer <- harvest %>%
      sf::st_buffer(dist = r) %>%
      dplyr::select(dplyr::all_of(time_col), geometry = geometry)

    # Spatial join: keep all points, attach harvest time where intersects
    joined <- sf::st_join(
      out,
      harvest_buffer,
      join = sf::st_intersects,
      left = TRUE
    )

    # For points that intersect multiple harvest polygons, keep earliest harvest time
    # Group by id_col and summarise. If no harvests overlap, we keep NA.
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

    # Join the earliest harvest time back to the main dataset
    out <- out %>%
      dplyr::left_join(harvest_time, by = id_col)

    # Compute time difference: harvest_time - point_time
    out[[diff_col]] <- out[[gis_col]] - out[[datetime_col]]  # difftime, can be NA

    # Classify relative timing
    out[[class_col]] <- dplyr::case_when(
      is.na(out[[diff_col]]) ~ "noharvest",
      out[[diff_col]] < 0    ~ "after",      # harvest happened before the point
      out[[diff_col]] > 0    ~ "before",     # harvest happened after the point
      out[[diff_col]] == 0   ~ "same_time"
    )

    out[[class_col]] <- factor(out[[class_col]], levels = c("noharvest", "after", "before", "same_time"))
  }

  return(out)
}
