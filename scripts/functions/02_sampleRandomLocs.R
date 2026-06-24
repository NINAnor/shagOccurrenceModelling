#' Generate available (background) points for habitat selection
#'
#' This function generates available (background) points for use-availability
#' analyses given observed dive locations and spatial constraints (study region,
#' kelp habitat 'tare', colony). It supports per-individual or grouped sampling,
#' optional temporal thinning of occurrences, and optional stratification by a
#' harvest buffer. Background points may be weighted by the sampled-area.
#'
#' Key behavior:
#'  - Observed locations are labelled used = 1; background points used = 0.
#'  - Two sampling modes:
#'      * group_by_id = TRUE  -> sample an accessible area per individual (ID)
#'      * group_by_id = FALSE -> sample per group (sex × Logger_model)
#'  - Optional thinning of occurrences by temporal day + spatial distance
#'    (DBSCAN clustering on coordinates within each individual × day).
#'  - Optional harvest-aware stratification: ensure a minimum fraction of
#'    background points fall inside a buffer around harvest polygons.
#'  - Background points keep metadata replicated from nearby observed points
#'    (timestamp, ID or sampled ID) to facilitate time-matching analyses.
#'
#' Inputs:
#'  - dloc               : sf POINTs of observed dive locations (must contain ID, UTC_datetime, sex, Logger_model)
#'  - study_region       : sf polygon of the study region
#'  - tare               : sf polygons of kelp habitat (available habitat)
#'  - colony             : sf point (or geometry) used to compute DistNest
#'  - n_per_point        : integer; number of background points per observed point
#'  - min_dist_m         : numeric; minimum buffer distance around observations for known-buffer (group mode)
#'  - id_buffer_m        : numeric; buffer distance (meters) around each ID for ID-specific accessible area
#'  - group_by_id        : logical; whether to sample per-ID accessible areas
#'  - thin_occ           : logical; whether to thin observed occurrences
#'  - thin_dist_km       : numeric; thinning distance in km
#'  - harvest            : sf polygons for harvest features (optional)
#'  - harvest_stratify   : logical; whether to stratify background sampling by harvest buffer
#'  - harvest_buffer_m   : numeric; buffer around harvest polygons (meters)
#'  - min_frac_in_bg     : numeric (0-1); ensure at least this fraction of background points fall inside harvest buffer
#'  - weight_background_by_area : logical; assign background weights proportional to sampled-area per stratum
#'
#' Returns:
#'  A list with:
#'    - combined_data: sf of used + available points (with w and DistNest)
#'    - used          : sf of observed (used) points
#'    - available     : sf of background points
#'    - available_regions : list of sampled available regions (sf objects)
#'    - known_buffer  : sf of known-buffer polygons used for sampling
#'    - tare          : unioned tare polygon used (sf)
#'    - plot          : ggplot diagnostic map
#'    - group_by_id   : boolean (same as input)
#'
#' Note:
#'  - Function avoids calling library() inside; it checks for required namespaces.
#'  - CRS of spatial layers is harmonized where necessary.
#'
#' @export
generate_available_points <- function(
    dloc,
    study_region,
    tare,
    colony,

    n_per_point = 5,
    min_dist_m = 500,
    id_buffer_m = 10000,
    group_by_id = FALSE,

    thin_occ = TRUE,
    thin_dist_km = 0.2,

    # harvest-aware background sampling
    harvest = NULL,
    harvest_stratify = TRUE,
    harvest_buffer_m = 200,
    min_frac_in_bg = 0.25,             # ensure at least this fraction inside harvested buffer
    weight_background_by_area = TRUE
) {
  # -----------------------
  # Dependency checks
  # -----------------------
  required_pkgs <- c("dplyr", "sf", "purrr", "ggplot2", "tidyr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Install missing packages: ", paste(missing_pkgs, collapse = ", "))

  `%>%` <- dplyr::`%>%`

  # -----------------------
  # Basic input validations
  # -----------------------
  if (!inherits(dloc, "sf")) stop("dloc must be an sf object")
  if (!inherits(study_region, "sf")) stop("study_region must be an sf object")
  if (!inherits(tare, "sf")) stop("tare must be an sf object")
  if (!inherits(colony, "sf")) stop("colony must be an sf object")

  # Ensure dloc contains essential columns used later
  required_cols <- c("ID", "UTC_datetime", "sex", "Logger_model")
  missing_cols <- setdiff(required_cols, names(dloc))
  if (length(missing_cols) > 0) {
    stop("dloc is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # -----------------------
  # Label used points and report
  # -----------------------
  used <- dloc %>% dplyr::mutate(used = 1)
  cat("Original dataset contains:", nrow(used), "points across",
      length(unique(as.character(used$ID))), "unique individuals.\n")

  # -----------------------
  # Union of tare polygons (available habitat)
  # -----------------------
  tare_union <- sf::st_union(tare)

  # -----------------------
  # Optional temporal + spatial thinning of occurrences
  # - Clusters points by individual × day and keeps first point per cluster
  # -----------------------
  if (isTRUE(thin_occ)) {
    if (!requireNamespace("dbscan", quietly = TRUE)) stop("Package 'dbscan' is required for thinning; install it or set thin_occ = FALSE.")
    # Determine a projected CRS for distance operations; prefer current CRS unless it's lon/lat
    crs_orig <- sf::st_crs(used)
    used_m <- if (sf::st_is_longlat(used)) sf::st_transform(used, 32632) else used
    used_m$.time_key <- as.Date(used_m$UTC_datetime)

    thin_dist_m <- thin_dist_km * 1000
    spec_col <- "ID"
    keys <- interaction(used_m[[spec_col]], used_m$.time_key, drop = TRUE)
    groups <- split(used_m, keys)

    thin_one_group <- function(g) {
      n <- nrow(g)
      if (n <= 1) return(g)
      coords <- sf::st_coordinates(g)
      cl <- dbscan::dbscan(coords, eps = thin_dist_m, minPts = 1)
      ord <- order(g$UTC_datetime)
      keep <- !duplicated(cl$cluster[ord])
      g[ord, ][keep, ]
    }

    used_thin_m <- do.call(rbind, lapply(groups, thin_one_group))
    used_thin_m$.time_key <- NULL

    used <- if (!is.na(crs_orig$epsg) || !is.null(crs_orig$wkt)) {
      sf::st_transform(used_thin_m, crs_orig)
    } else {
      used_thin_m
    }
    cat("After thinning:", nrow(used), "used points remain.\n")
  }

  # -----------------------
  # Build harvest buffer if requested
  # -----------------------
  harv_buf_valid <- NULL
  if (isTRUE(harvest_stratify) && !is.null(harvest) && nrow(harvest) > 0) {
    # Harmonize CRS
    if (sf::st_crs(harvest) != sf::st_crs(dloc)) harvest <- sf::st_transform(harvest, sf::st_crs(dloc))
    harv_union <- sf::st_make_valid(sf::st_union(harvest))
    harv_buf <- sf::st_buffer(harv_union, dist = harvest_buffer_m)
    harv_buf_valid <- sf::st_make_valid(harv_buf)
  }

  # -----------------------
  # Helper: compute area robustly (returns numeric 0 if null/empty)
  # -----------------------
  area_num <- function(g) {
    if (is.null(g)) return(0)
    if (length(g) == 0) return(0)
    if (all(sf::st_is_empty(g))) return(0)
    as.numeric(sf::st_area(g))
  }

  # -----------------------
  # Sampling: two modes
  # -----------------------
  if (isTRUE(group_by_id)) {
    # Per-ID accessible area constrained to tare
    avail_sf <- used %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
      sf::st_buffer(id_buffer_m) %>%
      sf::st_intersection(tare_union)

    available_list <- purrr::map(unique(as.character(used$ID)), function(id_i) {
      dloc_id <- used %>% dplyr::filter(as.character(ID) == id_i)
      n_id <- nrow(dloc_id) * n_per_point

      area_i <- avail_sf %>% dplyr::filter(as.character(ID) == id_i)

      # Harvest stratification if available
      if (!is.null(harv_buf_valid)) {
        inside_i  <- sf::st_intersection(area_i, harv_buf_valid)
        outside_i <- suppressWarnings(sf::st_difference(area_i, harv_buf_valid))

        A_in  <- area_num(inside_i)
        A_out <- area_num(outside_i)
        A_tot <- A_in + A_out

        frac_area_in <- if (A_tot > 0) A_in / A_tot else 0
        target_in <- max(min_frac_in_bg, frac_area_in)
        n_in  <- min(n_id, round(n_id * target_in))
        n_out <- n_id - n_in

        rnd_in  <- if (n_in  > 0 && A_in  > 0) sf::st_sample(inside_i,  size = n_in,  type = "random", exact = TRUE) else sf::st_sfc(crs = sf::st_crs(dloc))
        rnd_out <- if (n_out > 0 && A_out > 0) sf::st_sample(outside_i, size = n_out, type = "random", exact = TRUE) else sf::st_sfc(crs = sf::st_crs(dloc))
        rnd <- c(rnd_in, rnd_out)

        w_in  <- if (n_in  > 0) A_in  / n_in  else 0
        w_out <- if (n_out > 0) A_out / n_out else 0
        w_bg  <- c(rep(w_in, length(rnd_in)), rep(w_out, length(rnd_out)))
      } else {
        rnd <- sf::st_sample(area_i, size = n_id, type = "random", exact = TRUE)
        w_bg <- rep(1, length(rnd))
      }

      # replicate metadata from observed points and time-match randomly
      meta_source <- dloc_id %>% sf::st_drop_geometry()
      meta_rep <- meta_source %>%
        dplyr::slice(rep(dplyr::row_number(), each = n_per_point)) %>%
        dplyr::slice(seq_len(length(rnd)))

      # time-matching by resampling observed timestamps
      if (nrow(meta_rep) > 0) {
        meta_rep$UTC_datetime <- sample(meta_source$UTC_datetime, size = nrow(meta_rep), replace = TRUE)
      }

      available_group <- sf::st_sf(meta_rep, geometry = rnd, crs = sf::st_crs(dloc)) %>%
        dplyr::mutate(used = 0)

      # distance to colony
      if (nrow(available_group) > 0) {
        available_group$DistNest <- as.numeric(sf::st_distance(available_group, colony))
      } else {
        available_group$DistNest <- numeric(0)
      }

      # weights
      if (isTRUE(weight_background_by_area)) {
        available_group$w <- if (length(w_bg) > 0) w_bg else numeric(0)
      } else {
        available_group$w <- 1
      }

      list(
        available = available_group,
        available_region = area_i,
        known_buffer = area_i
      )
    })

    grouping_label <- "ID-specific buffers (constrained to tare; stratified by harvest)"

  } else {
    # Group-level sampling: sex × Logger_model
    group_data <- dloc %>% dplyr::group_by(sex, Logger_model) %>% dplyr::group_split()

    available_list <- purrr::map(group_data, function(dloc_group) {
      group_sex <- unique(dloc_group$sex)
      group_logger <- unique(dloc_group$Logger_model)
      n_group <- nrow(dloc_group) * n_per_point

      known_buffer_group <- sf::st_buffer(dloc_group, dist = min_dist_m)
      known_buffer_union <- sf::st_union(known_buffer_group)

      # available region = study_region minus known buffers, intersected with tare
      available_region_group <- sf::st_difference(sf::st_union(study_region), known_buffer_union) %>%
        sf::st_intersection(tare_union)

      # Harvest stratification
      if (!is.null(harv_buf_valid)) {
        inside_g  <- sf::st_intersection(available_region_group, harv_buf_valid)
        outside_g <- suppressWarnings(sf::st_difference(available_region_group, harv_buf_valid))

        A_in  <- area_num(inside_g)
        A_out <- area_num(outside_g)
        A_tot <- A_in + A_out

        frac_area_in <- if (A_tot > 0) A_in / A_tot else 0
        target_in <- max(min_frac_in_bg, frac_area_in)
        n_in  <- min(n_group, round(n_group * target_in))
        n_out <- n_group - n_in

        rnd_in  <- if (n_in  > 0 && A_in  > 0) sf::st_sample(inside_g,  size = n_in,  type = "random", exact = TRUE) else sf::st_sfc(crs = sf::st_crs(dloc))
        rnd_out <- if (n_out > 0 && A_out > 0) sf::st_sample(outside_g, size = n_out, type = "random", exact = TRUE) else sf::st_sfc(crs = sf::st_crs(dloc))
        rnd <- c(rnd_in, rnd_out)

        w_in  <- if (n_in  > 0) A_in  / n_in  else 0
        w_out <- if (n_out > 0) A_out / n_out else 0
        w_bg  <- c(rep(w_in, length(rnd_in)), rep(w_out, length(rnd_out)))
      } else {
        rnd <- sf::st_sample(available_region_group, size = n_group, type = "random", exact = TRUE)
        w_bg <- rep(1, length(rnd))
      }

      if (length(rnd) < n_group) {
        warning(paste0("Fewer available points sampled for group: sex=", group_sex, ", logger=", group_logger))
      }

      meta_source <- dloc_group %>% sf::st_drop_geometry()
      meta_rep <- meta_source %>%
        dplyr::slice(rep(dplyr::row_number(), each = n_per_point)) %>%
        dplyr::slice(seq_len(length(rnd)))

      # randomize IDs and timestamps to avoid pairing by original order
      if (nrow(meta_rep) > 0) {
        meta_rep$ID <- sample(meta_source$ID, size = nrow(meta_rep), replace = TRUE)
        meta_rep$UTC_datetime <- sample(meta_source$UTC_datetime, size = nrow(meta_rep), replace = TRUE)
      }

      available_group <- sf::st_sf(meta_rep, geometry = rnd, crs = sf::st_crs(dloc)) %>%
        dplyr::mutate(used = 0)

      if (nrow(available_group) > 0) {
        available_group$DistNest <- as.numeric(sf::st_distance(available_group, colony))
      } else {
        available_group$DistNest <- numeric(0)
      }

      if (isTRUE(weight_background_by_area)) {
        available_group$w <- if (length(w_bg) > 0) w_bg else numeric(0)
      } else {
        available_group$w <- 1
      }

      list(
        available = available_group,
        available_region = available_region_group,
        known_buffer = known_buffer_group
      )
    })

    grouping_label <- "sex × Logger_model (constrained to tare; stratified by harvest)"
  }

  # -----------------------
  # Combine results
  # -----------------------
  available <- purrr::map_dfr(available_list, "available")
  available_regions <- purrr::map(available_list, "available_region")
  known_buffers <- purrr::map_dfr(available_list, "known_buffer")

  used$DistNest <- as.numeric(sf::st_distance(used, colony))
  used$w <- 1  # presence weight

  combined_data <- dplyr::bind_rows(used, available)

  # -----------------------
  # Diagnostic plot
  # -----------------------
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = sf::st_union(study_region), fill = "grey95", color = "grey60") +
    ggplot2::geom_sf(data = tare_union, fill = NA, color = "darkgreen", linetype = "dashed") +
    ggplot2::geom_sf(data = known_buffers, fill = "red", color = NA, alpha = 0.2) +
    ggplot2::geom_sf(data = used, color = "red", size = 2, alpha = 0.8) +
    ggplot2::geom_sf(data = available, color = "blue", size = 1, alpha = 0.5) +
    ggplot2::geom_sf(data = colony, color = "black", size = 4) +
    ggplot2::facet_grid(sex ~ Logger_model) +
    ggplot2::labs(
      title = "Used and available locations",
      subtitle = paste0(
        n_per_point, " available per used; grouping = ", grouping_label,
        if (!is.null(harv_buf_valid)) paste0("; harvest_buffer_m = ", harvest_buffer_m) else ""
      )
    ) +
    ggplot2::theme_minimal()

  # -----------------------
  # Return structured output
  # -----------------------
  return(list(
    combined_data = combined_data,     # sf: used + available (columns: w, DistNest, used, etc.)
    used = used,
    available = available,
    available_regions = available_regions,
    known_buffer = known_buffers,
    tare = tare_union,
    plot = p,
    group_by_id = group_by_id
  ))
}
