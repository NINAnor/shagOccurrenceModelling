#' Fit habitat selection GAM and produce predictions & diagnostic plots
#'
#' This wrapper fits a binomial GAM to a spatial use/availability dataset and
#' returns the fitted model, model data, prediction objects, and diagnostic plots.
#'
#' Key behaviour and assumptions:
#'  - Expects dmodel to be an sf object with a POINT geometry column.
#'  - Expects year-specific raster objects (raster_results_23 / raster_results_24)
#'    that contain distance-to-harvest and distance-to-nest layers (e.g. DistHarvest,
#'    DistNest) and harvest buffers (harvest50_buffer, harvest200_buffer).
#'  - Expects env_stack_10m to be a terra SpatRaster or SpatRaster collection of
#'    additional covariates for prediction.
#'  - Uses mgcv::bam() for scalable GAM fitting (discrete = TRUE).
#'  - Produces partial-effect curves for DistHarvest (excluding spatial+ID random
#'    effects) and relative-selection (RS = exp(centered link)).
#'  - Calls maps_by_sex_all_periods(...) (project-specific helper) to produce
#'    spatial prediction rasters/plots by sex × harvest period.
#'
#' @param dmodel sf data.frame with model points (must include Year, geometry, ID, etc.)
#' @param pred_grid data.frame or sf object of prediction locations (must be convertible to sf)
#' @param study_region sf polygon of the study area (used downstream in mapping)
#' @param env_stack_10m terra SpatRaster of 10m environmental covariates
#' @param raster_results_23 list-like object for 2023 rasters (must include $raster_stack, $harvest50_buffer, $harvest200_buffer)
#' @param raster_results_24 list-like object for 2024 rasters (same structure as for 2023)
#' @param pars_to_use named list (must include $Year; value "2023" or "2024")
#' @param response_col string name of column in dmodel to use as binary response (default "used")
#' @param nBasis integer upper-bound for k in the 2D spatial smooth (default 150; reduced heuristically)
#' @param nthreads integer; number of threads passed to mgcv::bam (default 1)
#'
#' @return A list containing:
#'   - method: "gam"
#'   - parameters: the pars_to_use passed in
#'   - model: the fitted mgcv::bam object
#'   - model_data: data.frame used in the fit (with x,y coords and response)
#'   - model_sf: cleaned sf used for modeling (with covariates)
#'   - pred_grid_env: prediction grid with extracted covariates
#'   - dist_harvest_effect_plot: ggplot object of predicted probability vs DistHarvest
#'   - rel_selection_plot: ggplot object of relative selection vs DistHarvest
#'   - res_all: output of maps_by_sex_all_periods(...) (spatial predictions & plots)
#'
#' @keywords habitat selection GAM prediction
#' @export
fit_predict_habitat_model <- function(
    dmodel,
    pred_grid,
    study_region,
    env_stack_10m,
    raster_results_23,
    raster_results_24,
    pars_to_use,
    response_col = "used",
    nBasis = 150,
    nthreads = 1
) {
  # -------------------------
  # 0) Check package dependencies
  # -------------------------
  pkgs <- c("sf", "terra", "dplyr", "mgcv", "ggplot2", "tidyr")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))

  `%>%` <- dplyr::`%>%`  # local pipe binding for readability

  # -------------------------
  # 1) Filter model data for the requested year and deduplicate points
  #    - Keep one point per unique WKT geometry to avoid accidental duplicates
  # -------------------------
  model_sf <- dmodel %>%
    dplyr::filter(Year == pars_to_use$Year) %>%
    dplyr::mutate(
      .rowid = dplyr::row_number(),
      geom_wkt = sf::st_as_text(geometry)
    ) %>%
    dplyr::distinct(geom_wkt, .keep_all = TRUE) %>%
    dplyr::select(-geom_wkt)

  # -------------------------
  # 2) Select year-specific raster stack and harvest buffer polygons
  # -------------------------
  if (pars_to_use$Year == "2023") {
    harvest_stack <- raster_results_23$raster_stack
    harvest_regions <- raster_results_23$harvest50_buffer
    harvest200_regions <- raster_results_23$harvest200_buffer
  } else if (pars_to_use$Year == "2024") {
    harvest_stack <- raster_results_24$raster_stack
    harvest_regions <- raster_results_24$harvest50_buffer
    harvest200_regions <- raster_results_24$raster200_buffer
    # NOTE: original code used harvest200_buffer for 2024 as well; if yours is named differently update accordingly
  } else {
    stop("pars_to_use$Year must be '2023' or '2024'.")
  }

  # Create a simple sf of harvest polygons (used later for mapping)
  harvest_regions_sf <- sf::st_sf(
    harvested50m = 1L,
    geometry = harvest_regions,
    crs = sf::st_crs(harvest_regions)
  )

  # Combine harvest distance rasters with environmental rasters
  # (terra SpatRaster concatenation)
  allRasts <- c(harvest_stack, env_stack_10m)

  # -------------------------
  # 3) Extract raster covariates for model points
  #    - Ensure consistent CRS and that rasters use a projected CRS with meter units
  # -------------------------
  r_crs <- terra::crs(allRasts, proj = TRUE)
  if (is.na(r_crs) || r_crs == "") stop("allRasts has no CRS. Set a projected CRS with meter units.")

  model_sf_r <- sf::st_transform(model_sf, r_crs)

  # terra::extract returns a data.frame with ID column (row number); rename to .rowid for join
  env_values_model <- terra::extract(allRasts, terra::vect(model_sf_r)) %>%
    dplyr::rename(.rowid = ID) %>%
    dplyr::select(-dplyr::any_of(c("harvested50m", "harvested200m")))

  # Remove DistNest column from original sf (will be replaced by raster-extracted values)
  model_sf <- model_sf %>% dplyr::select(-dplyr::any_of("DistNest"))

  # Join extracted raster covariates back to the model sf using .rowid
  model_sf <- model_sf %>% dplyr::left_join(env_values_model, by = ".rowid")

  # -------------------------
  # 4) Extract raster covariates for prediction grid
  # -------------------------
  pred_grid_sf <- sf::st_as_sf(pred_grid)
  pred_grid_r  <- sf::st_transform(pred_grid_sf, r_crs)
  env_values_pred <- terra::extract(allRasts, terra::vect(pred_grid_r)) %>% dplyr::select(-ID)
  pred_grid_env <- pred_grid_r %>% dplyr::bind_cols(env_values_pred)

  # -------------------------
  # 5) Prepare modeling dataframe: response, weights, necessary covariates
  # -------------------------
  if (!"w" %in% names(model_sf)) model_sf$w <- 1  # default weight

  # Verify response column exists and is binary 0/1 (coerce factor if needed)
  if (!response_col %in% names(model_sf)) stop("response_col not found in model data")
  if (is.factor(model_sf[[response_col]])) {
    model_sf[[response_col]] <- as.numeric(as.character(model_sf[[response_col]]))
  }
  if (!all(model_sf[[response_col]] %in% c(0, 1))) stop("response_col must contain only 0/1 values")

  # Ensure required predictors exist after extraction
  need_cols <- c("DistNest", "DistHarvest", "bathymetry")
  miss_cols <- setdiff(need_cols, names(model_sf))
  if (length(miss_cols)) stop("Missing predictor(s) in model data after extraction: ", paste(miss_cols, collapse = ", "))

  # Ensure categorical variables exist (create default levels if absent)
  if (!"sex" %in% names(model_sf)) model_sf$sex <- factor("Unknown")
  if (!"harvesttimeSK24" %in% names(model_sf)) model_sf$harvesttimeSK24 <- factor("all")
  if (!"ID" %in% names(model_sf)) model_sf$ID <- factor("all")

  # Select and drop NA rows
  model_sf_clean <- model_sf %>%
    dplyr::select(ID, sex, harvesttimeSK24, Year,
                  dplyr::all_of(response_col),
                  DistNest, DistHarvest, bathymetry, w, geometry) %>%
    tidyr::drop_na()

  # Normalize weights so that mean(w) == 1 to keep effective sample sizes interpretable
  model_sf_clean$w <- model_sf_clean$w / mean(model_sf_clean$w)

  # Extract coordinates for spatial smooth
  coords_model <- sf::st_coordinates(model_sf_clean)

  model_data <- model_sf_clean %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(
      response = .data[[response_col]],
      x = coords_model[, 1],
      y = coords_model[, 2],
      ID = factor(ID),
      sex = factor(sex),
      harvesttimeSK24 = factor(harvesttimeSK24),
      w = w
    )

  # -------------------------
  # 6) Fit GAM
  #    - Uses a spatial thin-plate spline s(x,y)
  #    - Includes random effect for ID (s(ID, bs="re"))
  #    - Interactions of DistHarvest with sex and harvesttimeSK24
  # -------------------------
  # Heuristic for k: at most nBasis but with a minimum relative to sample size
  nBasis <- min(nBasis, max(40, round(nrow(model_data) / 30)))

  fitted_model <- mgcv::bam(
    response ~
      s(DistNest, k = 5) +
      DistHarvest * sex +
      DistHarvest * harvesttimeSK24 +
      s(bathymetry, k = 5) +
      s(ID, bs = "re") +                 # random effect per ID (cheap & interpretable)
      s(x, y, bs = "tp", k = nBasis),
    data     = model_data,
    family   = binomial(),
    weights  = w,
    method   = "fREML",
    discrete = TRUE,        # speeds up fitting on large datasets
    nthreads = nthreads,    # threading control for bam
    gc.level = 2            # more aggressive memory management during fit
  )

  # -------------------------
  # 7) Partial effect curves for DistHarvest (exclude spatial & ID random effects)
  # -------------------------
  .med <- function(x) if (is.numeric(x)) stats::median(x, na.rm = TRUE) else x[1]

  # Build list of smooth labels to exclude spatial and ID terms from partial predictions
  smooth_labels <- vapply(fitted_model$smooth, function(sm) sm$label, character(1))
  exclude_terms <- smooth_labels[
    startsWith(smooth_labels, "s(x,y") | startsWith(smooth_labels, "s(ID)")
  ]

  # Use the central 2..98% quantiles of DistHarvest for effect plotting
  qDH    <- stats::quantile(model_data$DistHarvest, c(0.02, 0.98), na.rm = TRUE)
  DH_seq <- seq(qDH[1], qDH[2], length.out = 150)

  # Base covariate row: hold other covariates at typical values
  base_row <- data.frame(
    DistNest    = .med(model_data$DistNest),
    bathymetry  = .med(model_data$bathymetry),
    x           = .med(model_data$x),
    y           = .med(model_data$y),
    ID          = model_data$ID[1]
  )

  # Expand base_row across DistHarvest, sexes, and harvesttimeSK24 levels
  lev_period <- levels(model_data$harvesttimeSK24)
  lev_sex    <- levels(model_data$sex)

  effect_grid <- do.call(rbind, lapply(lev_sex, function(sx) {
    do.call(rbind, lapply(lev_period, function(pr) {
      nd <- base_row[rep(1, length(DH_seq)), , drop = FALSE]
      nd$DistHarvest     <- DH_seq
      nd$harvesttimeSK24 <- factor(pr, levels = lev_period)
      nd$sex             <- factor(sx, levels = lev_sex)
      nd
    }))
  }))
  row.names(effect_grid) <- NULL

  # Predict on the link scale excluding spatial & ID random effects
  pred_link <- predict(
    fitted_model,
    newdata = effect_grid,
    type    = "link",
    se.fit  = TRUE,
    exclude = exclude_terms
  )

  effect_grid$eta  <- as.numeric(pred_link$fit)
  effect_grid$se   <- as.numeric(pred_link$se.fit)
  effect_grid$p    <- stats::plogis(effect_grid$eta)
  effect_grid$p_lo <- stats::plogis(effect_grid$eta - 1.96 * effect_grid$se)
  effect_grid$p_hi <- stats::plogis(effect_grid$eta + 1.96 * effect_grid$se)

  # Center the link by sex × harvest period and compute Relative Selection (RS)
  effect_grid <- effect_grid %>%
    dplyr::group_by(sex, harvesttimeSK24) %>%
    dplyr::mutate(
      eta_ctr = eta - mean(eta, na.rm = TRUE),
      RS = exp(eta_ctr)
    ) %>%
    dplyr::ungroup()

  # -------------------------
  # 8) Plotting: partial probability and RS
  # -------------------------
  dist_harvest_effect_plot <- ggplot2::ggplot(
    effect_grid,
    ggplot2::aes(x = DistHarvest, y = p, color = harvesttimeSK24, fill = harvesttimeSK24)
  ) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = p_lo, ymax = p_hi), alpha = 0.15, color = NA) +
    ggplot2::geom_line(size = 1) +
    ggplot2::facet_wrap(~ sex, nrow = 1) +
    ggplot2::labs(
      title = paste("Partial effect of distance to harvest for year:", pars_to_use$Year),
      subtitle = "Predicted probability; excluding s(x,y) and s(ID)",
      x = "Distance to harvest",
      y = "Predicted probability"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::scale_fill_brewer(palette = "Dark2")

  rel_selection_plot <- ggplot2::ggplot(
    effect_grid,
    ggplot2::aes(x = DistHarvest, y = RS, color = harvesttimeSK24)
  ) +
    ggplot2::geom_line(size = 1) +
    ggplot2::facet_wrap(~ sex, nrow = 1) +
    ggplot2::labs(
      title = "Relative selection vs. distance to harvest",
      subtitle = "RS = exp(centered link); excluding s(x,y) and s(ID)",
      x = "Distance to harvest",
      y = "Relative selection (RS)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_brewer(palette = "Dark2")

  # -------------------------
  # 9) Generate spatial predictions & plots by sex × harvest period
  #    - maps_by_sex_all_periods(...) is a project-specific helper that should:
  #       * accept the fitted model, model_data, and pred_grid_env
  #       * return a list containing spatial prediction rasters and plot objects
  #    - Here we pass response_col = "random" in original code; adjust if necessary.
  # -------------------------
  res_all <- maps_by_sex_all_periods(
    fitted_model     = fitted_model,
    model_data       = model_data,
    pred_grid_env    = pred_grid_env,
    study_region     = study_region,
    harvest_polygons = harvest_regions_sf,
    env_stack_10m    = env_stack_10m,
    response_col     = "random",   # keep consistent with upstream usage; change to response_col if desired
    point_labels     = c("0" = "random = 0", "1" = "random = 1"),
    exclude_re       = TRUE
  )

  # -------------------------
  # 10) Return structured results for downstream saving/plotting
  # -------------------------
  list(
    method = "gam",
    parameters = pars_to_use,
    model = fitted_model,
    model_data = model_data,
    model_sf = model_sf_clean,
    pred_grid_env = pred_grid_env,
    dist_harvest_effect_plot = dist_harvest_effect_plot,
    rel_selection_plot = rel_selection_plot,
    res_all = res_all
  )
}
