#' Build full modelling dataset for a single year
#'
#' This function prepares a complete used-available modelling dataset for one
#' sampling year by:
#'  1. clipping dive locations to kelp habitat around the colony,
#'  2. resolving duplicate intersections,
#'  3. generating available (background) points (optionally per-ID or grouped),
#'  4. adding harvest overlap / timing covariates.
#'
#' Inputs:
#'  - dloc_year       : sf POINTs of dive locations for the focal year
#'  - tare            : sf POLYGON(s) of kelp habitat
#'  - colony          : sf geometry (point) for the colony location
#'  - harvest         : sf of harvest events (must include a 'year' column if you want per-year filtering)
#'  - year            : character or numeric year to process (e.g. "2024")
#'  - n_per_point     : background points per observed point
#'  - min_dist_m      : minimum known-buffer distance around observations (group mode)
#'  - ID_buffer_m     : buffer radius around each individual's locations (ID mode)
#'  - harvest_buffer  : buffer radius (m) used for harvest stratification when sampling background
#'  - thin_dist_km    : thinning distance (km) used when temporal+spatial thinning is enabled
#'  - min_frac_in_bg  : minimum fraction of background points targeted inside harvest buffer
#'  - group_by_ID     : logical; sample accessible area per ID when TRUE, otherwise sample by group (sex × Logger_model)
#'
#' Returns:
#'  A list with:
#'    - dloctare             : outputs from process_dloctare (including raw intersection)
#'    - dloctare_resolved    : sf of dive locations resolved for duplicates
#'    - avail_results        : output from generate_available_points (available points, regions, etc.)
#'    - dmodel               : final used + available sf with harvest overlap covariates
#'
#' @export
make_year_dmodel <- function(dloc_year,
                             tare,
                             colony,
                             harvest,
                             year,
                             n_per_point = 5,
                             min_dist_m = 500,
                             ID_buffer_m = 10000,
                             harvest_buffer = 50,
                             thin_dist_km = 1,
                             min_frac_in_bg  = 0.5,
                             group_by_ID = TRUE,
                             thin_occ = TRUE,
                             harvest_stratify = TRUE,
                             weight_background_by_area = TRUE
                             ) {

  # -------------------------
  # Dependency checks
  # -------------------------
  required_pkgs <- c("dplyr", "sf", "ggplot2")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Install missing packages: ", paste(missing_pkgs, collapse = ", "))

  `%>%` <- dplyr::`%>%`

  # -------------------------
  # 1) Clip dive locations to kelp habitat around colony and resolve duplicates
  # -------------------------
  dloctare <- process_dloctare(
    dloc = dloc_year,
    tare = tare,
    colony = colony,
    year = as.character(year)
  )

  # Resolved dive locations (ensure naturtype is a factor if present)
  dloctare_resolved <- dloctare$dloctare_resolved
  if ("naturtype" %in% names(dloctare_resolved)) {
    dloctare_resolved <- dloctare_resolved %>%
      dplyr::mutate(naturtype = as.factor(naturtype))
  }

  # -------------------------
  # 2) Prepare harvest subset for the focal year (if harvest has a 'year' column)
  # -------------------------
  harvest_year <- harvest
  if ("year" %in% names(harvest)) {
    harvest_year <- harvest %>% dplyr::filter(.data$year == !!year)
  } else {
    # If no year column is present, use harvest as-is and warn the user
    warning("harvest has no 'year' column: using full harvest dataset for sampling.")
  }

  # -------------------------
  # 3) Generate available (background) points
  # -------------------------
  avail_results <- generate_available_points(
    dloc = dloctare_resolved,
    study_region = dloctare$SK33km,
    tare = tare,
    colony = colony,
    n_per_point = n_per_point,
    min_dist_m = min_dist_m,
    id_buffer_m = ID_buffer_m,
    group_by_id = isTRUE(group_by_ID),
    thin_occ = thin_occ,
    thin_dist_km = thin_dist_km,
    harvest = harvest_year,
    harvest_stratify = harvest_stratify,
    harvest_buffer_m = harvest_buffer,
    min_frac_in_bg = min_frac_in_bg,
    weight_background_by_area = weight_background_by_area
  )

  combined_data <- avail_results$combined_data

  # -------------------------
  # 4) Add harvest overlap / timing covariates to combined dataset
  #    - uses radii 200m and 50m as in original workflow
  # -------------------------
  dmodel_year <- add_harvest_overlap(
    dmodel = combined_data,
    harvest = harvest_year,
    radii_m = c(200, 50),
    time_col = "GIS_TIME",
    datetime_col = "UTC_datetime",
    id_col = "seq"
  )

  # -------------------------
  # 5) Return a structured list of outputs for inspection / downstream use
  # -------------------------
  list(
    dloctare = dloctare,
    dloctare_resolved = dloctare_resolved,
    avail_results = avail_results,
    dmodel = dmodel_year
  )
}
