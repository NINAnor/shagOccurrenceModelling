#' Predict sex-specific maps and compose facetted maps by sex (and period)
#'
#' predict_maps_by_sex:
#'  - Given a fitted GAM-like model and prediction grid with covariates,
#'    produce sex-specific predicted probability rasters, a long prediction
#'    table, and a facetted ggplot (one facet per sex).
#'
#' maps_by_sex_all_periods:
#'  - Runs predict_maps_by_sex for every level of harvesttimeSK24 (periods),
#'    assembles a multi-layer SpatRaster with one layer per sex × period, and
#'    creates a single facetted plot with sex × period facets.
#'
#' Assumptions / requirements:
#'  - pred_grid_env is an sf object with extracted environmental covariates
#'    and point geometry (projected CRS consistent with study_region).
#'  - model_data contains factor columns: sex, harvesttimeSK24, ID and numeric
#'    x/y coordinates used for plotting, and a binary response column.
#'  - fitted_model supports predict(..., type = "response") and accepts an
#'    `exclude` argument (e.g., mgcv::gam/bam objects).
#'
#' Usage:
#'  - Use predict_maps_by_sex() to produce sex-specific maps for a single period.
#'  - Use maps_by_sex_all_periods() to produce maps for all periods and a combined plot.
#'
#' @export
predict_maps_by_sex <- function(fitted_model,
                                model_data,        # data.frame used to fit the model (must include x,y)
                                pred_grid_env,     # sf grid with env covariates (and geometry)
                                study_region,      # sf polygon(s) defining extent
                                harvest_polygons = NULL,   # optional sf polygons for overlay (e.g., harvest)
                                env_stack_10m = NULL,      # optional SpatRaster (e.g., mask layer "tareskog")
                                response_col = NULL,       # e.g., "random", "used", or "response"
                                point_labels = c("0" = "random = 0", "1" = "random = 1"),
                                period_level = NULL,       # e.g., "duringJuneHarvest"; NULL -> first level
                                exclude_re = TRUE) {

  # Basic input assertions
  stopifnot(inherits(pred_grid_env, "sf"), inherits(study_region, "sf"))

  # Determine response column if not provided (prefer common names)
  if (is.null(response_col)) {
    possible <- c("random", "used", "response")
    response_col <- possible[possible %in% names(model_data)][1]
    if (is.na(response_col)) stop("response_col not found. Provide response_col (e.g., 'random'/'used').")
  }
  if (!response_col %in% names(model_data)) stop("response_col '", response_col, "' not in model_data.")

  # Ensure required factor columns exist
  if (!"sex" %in% names(model_data)) stop("model_data must contain factor column 'sex'.")
  if (!"harvesttimeSK24" %in% names(model_data)) stop("model_data must contain factor column 'harvesttimeSK24'.")
  if (!"ID" %in% names(model_data)) stop("model_data must contain 'ID' (factor/character).")

  lev_sex <- levels(model_data$sex)
  if (is.null(lev_sex)) stop("model_data$sex must be a factor with levels.")
  lev_period <- levels(model_data$harvesttimeSK24)
  if (is.null(lev_period)) stop("model_data$harvesttimeSK24 must be a factor with levels.")

  # Select period to map
  if (is.null(period_level)) period_level <- lev_period[1]
  if (!period_level %in% lev_period) stop("period_level '", period_level, "' not found in harvesttimeSK24 levels.")

  # Build base prediction data.frame from pred_grid_env
  coords_pred <- sf::st_coordinates(pred_grid_env)
  base_df <- pred_grid_env %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(
      x = coords_pred[, 1],
      y = coords_pred[, 2],
      ID = model_data$ID[1],  # use a valid ID level; excluded if exclude_re = TRUE in predict
      harvesttimeSK24 = factor(period_level, levels = lev_period)
    )

  # Predict for each sex level and store predictions in a list of data.frames
  pred_list <- lapply(lev_sex, function(sx) {
    nd <- base_df
    nd$sex <- factor(sx, levels = lev_sex)
    nd$pred_prob <- as.numeric(
      predict(fitted_model, newdata = nd, type = "response",
              exclude = if (exclude_re) "s(ID)" else NULL)
    )
    nd
  })
  names(pred_list) <- lev_sex

  # Helper to turn each prediction df into a terra raster layer, masked to tareskog & study region
  make_raster <- function(df, layer_name) {
    r <- terra::rast(df[, c("x", "y", "pred_prob")], type = "xyz", crs = sf::st_crs(study_region)$wkt)
    names(r) <- layer_name

    # If env_stack_10m contains a logical mask layer "tareskog", use it to mask predicted raster
    if (!is.null(env_stack_10m) && "tareskog" %in% names(env_stack_10m)) {
      mask <- env_stack_10m[["tareskog"]] == 1
      mask_res <- terra::resample(mask, r, method = "near")
      r <- terra::mask(r, mask_res, maskvalues = 0)
    }

    # Crop and mask to the study region extent for consistency
    r <- terra::crop(r, terra::vect(study_region))
    r <- terra::mask(r, terra::vect(study_region))
    r
  }

  # Create a SpatRaster layer per sex and combine into a multi-layer SpatRaster
  rast_list <- Map(function(df, sx) make_raster(df, paste0("pred_prob_", sx)), pred_list, names(pred_list))
  pred_raster_by_sex <- do.call(c, rast_list)

  # Build a long data.frame for plotting raster layers with sex labels
  pred_df_plot <- do.call(rbind, lapply(names(rast_list), function(sx) {
    d <- terra::as.data.frame(rast_list[[sx]], xy = TRUE, na.rm = TRUE)
    names(d)[3] <- "pred_prob"
    d$sex <- factor(sx, levels = lev_sex)
    d
  }))

  # Prepare occurrence points (model_data must contain projected numeric x,y)
  if (!all(c("x", "y") %in% names(model_data))) {
    stop("model_data must include numeric columns x and y (projected coordinates) for plotting points.")
  }
  pts_sf <- sf::st_as_sf(model_data, coords = c("x", "y"), crs = sf::st_crs(study_region), remove = FALSE)
  pts_sf <- pts_sf %>% dplyr::select(sex, !!rlang::sym(response_col), x, y, geometry)

  # Convert response to 0/1 and create labeled factor for plotting
  resp_val <- pts_sf[[response_col]]
  if (is.factor(resp_val)) resp_val <- as.numeric(as.character(resp_val))
  pts_sf$resp01 <- as.integer(resp_val)
  pts_sf$point_type <- factor(
    ifelse(pts_sf$resp01 == 1, point_labels["1"], point_labels["0"]),
    levels = unname(point_labels[c("0", "1")])
  )

  # Prepare harvest polygon overlay if provided; ensure CRS match and clip to study region
  harvest_sf <- NULL
  if (!is.null(harvest_polygons)) {
    if (!inherits(harvest_polygons, "sf")) stop("harvest_polygons must be an sf object.")
    harvest_sf <- harvest_polygons
    if (sf::st_crs(harvest_sf) != sf::st_crs(study_region)) harvest_sf <- sf::st_transform(harvest_sf, sf::st_crs(study_region))
    harvest_sf <- suppressWarnings(sf::st_intersection(harvest_sf, sf::st_union(study_region)))
  }

  # Construct facetted plot by sex
  prediction_plot_by_sex <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = pred_df_plot, ggplot2::aes(x = x, y = y, fill = pred_prob)) +
    { if (!is.null(harvest_sf)) ggplot2::geom_sf(data = harvest_sf, fill = NA, color = "red", linewidth = 0.6) else ggplot2::geom_blank() } +
    ggplot2::geom_sf(data = pts_sf, ggplot2::aes(color = point_type), size = 0.35, alpha = 0.1) +
    ggplot2::geom_sf(data = study_region, fill = NA, color = "black", linewidth = 0.4) +
    ggplot2::facet_wrap(~ sex, nrow = 1) +
    ggplot2::coord_sf(crs = sf::st_crs(study_region)) +
    ggplot2::scale_fill_viridis_c(name = "Predicted probability", na.value = NA) +
    ggplot2::scale_color_manual(name = "Points", values = setNames(c("#1f78b4", "#e31a1c"), levels(pts_sf$point_type))) +
    ggplot2::labs(
      title = "Predicted habitat suitability (sex-specific baselines)",
      subtitle = paste("Period:", period_level, "| ID effect excluded in map predictions"),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal()

  # Long table of predictions: bind pred_list into one data.frame with sex column
  pred_df_by_sex <- dplyr::bind_rows(pred_list, .id = "sex") %>% dplyr::mutate(sex = factor(sex, levels = lev_sex))

  list(
    pred_raster_by_sex     = pred_raster_by_sex,    # SpatRaster (one layer per sex)
    pred_df_by_sex         = pred_df_by_sex,        # long table (grid × sex)
    prediction_plot_by_sex = prediction_plot_by_sex # ggplot facetted by sex
  )
}


#' Run sex-specific prediction maps for all harvest periods and compile outputs
#'
#' For each level of model_data$harvesttimeSK24, calls predict_maps_by_sex(),
#' accumulates the returned rasters into a larger SpatRaster with one layer per
#' sex × period, and creates a combined facetted plot (sex × period).
#'
#' @export
maps_by_sex_all_periods <- function(fitted_model,
                                    model_data,
                                    pred_grid_env,
                                    study_region,
                                    harvest_polygons = NULL,
                                    env_stack_10m = NULL,
                                    response_col = "random",
                                    point_labels = c("0" = "random = 0", "1" = "random = 1"),
                                    exclude_re = TRUE) {
  # Dependency checks
  for (pkg in c("sf", "terra", "dplyr", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package: ", pkg)
  }

  # Validate factors and their levels
  if (!("harvesttimeSK24" %in% names(model_data)) || is.null(levels(model_data$harvesttimeSK24))) {
    stop("model_data$harvesttimeSK24 must be a factor with levels.")
  }
  lev_period <- levels(model_data$harvesttimeSK24)
  lev_sex <- levels(model_data$sex)
  if (is.null(lev_sex)) stop("model_data$sex must be a factor with levels.")

  # 1) Run predict_maps_by_sex() for each period and store results
  res_by_period <- setNames(vector("list", length(lev_period)), lev_period)
  for (pr in lev_period) {
    res_by_period[[pr]] <- predict_maps_by_sex(
      fitted_model = fitted_model,
      model_data = model_data,
      pred_grid_env = pred_grid_env,
      study_region = study_region,
      harvest_polygons = harvest_polygons,
      env_stack_10m = env_stack_10m,
      response_col = response_col,
      point_labels = point_labels,
      period_level = pr,
      exclude_re = exclude_re
    )
  }

  # 2) Combine per-period raster stacks into a single SpatRaster
  rast_list <- lapply(names(res_by_period), function(pr) {
    r <- res_by_period[[pr]]$pred_raster_by_sex
    names(r) <- paste0(names(r), "_", pr)  # append period label to layer names
    r
  })
  pred_raster_by_sex_period <- do.call(c, rast_list)

  # 3) Build a single facetted data.frame (long) for sex × period plotting
  pred_df_all <- do.call(rbind, lapply(names(res_by_period), function(pr) {
    r <- res_by_period[[pr]]$pred_raster_by_sex
    sexes <- gsub("^pred_prob_", "", names(r))
    do.call(rbind, lapply(seq_along(sexes), function(i) {
      d <- terra::as.data.frame(r[[i]], xy = TRUE, na.rm = TRUE)
      names(d)[3] <- "pred_prob"
      d$sex <- factor(sexes[i], levels = lev_sex)
      d$harvesttimeSK24 <- factor(pr, levels = lev_period)
      d
    }))
  }))

  # Prepare points for overlay, ensuring x/y columns exist
  if (!all(c("x", "y") %in% names(model_data))) stop("model_data must contain numeric columns x and y for plotting points.")
  pts_sf_all <- sf::st_as_sf(model_data, coords = c("x", "y"), crs = sf::st_crs(study_region), remove = FALSE)

  # Convert response to integer 0/1 and create labeled factor for plotting
  resp_val <- model_data[[response_col]]
  if (is.factor(resp_val)) resp_val <- as.numeric(as.character(resp_val))
  pts_sf_all$resp01 <- as.integer(resp_val)
  pts_sf_all$point_type <- factor(ifelse(pts_sf_all$resp01 == 1, point_labels["1"], point_labels["0"]),
                                  levels = unname(point_labels[c("0", "1")]))

  # Ensure facetting factors are present and properly factored
  pts_sf_all$sex <- factor(pts_sf_all$sex, levels = lev_sex)
  pts_sf_all$harvesttimeSK24 <- factor(pts_sf_all$harvesttimeSK24, levels = lev_period)

  # Harvest polygon overlay handling (CRS harmonization and clipping)
  harvest_sf <- NULL
  if (!is.null(harvest_polygons)) {
    harvest_sf <- harvest_polygons
    if (!inherits(harvest_sf, "sf")) stop("harvest_polygons must be an sf object.")
    if (sf::st_crs(harvest_sf) != sf::st_crs(study_region)) harvest_sf <- sf::st_transform(harvest_sf, sf::st_crs(study_region))
    harvest_sf <- suppressWarnings(sf::st_intersection(harvest_sf, sf::st_union(study_region)))
  }

  # 4) Combined facetted plot (sex × period)
  prediction_plot_sex_period <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = pred_df_all, ggplot2::aes(x = x, y = y, fill = pred_prob)) +
    { if (!is.null(harvest_sf)) ggplot2::geom_sf(data = harvest_sf, fill = NA, color = "red", linewidth = 0.6) else ggplot2::geom_blank() } +
    ggplot2::geom_sf(data = pts_sf_all, ggplot2::aes(color = point_type), size = 0.35, alpha = 0.6) +
    ggplot2::geom_sf(data = study_region, fill = NA, color = "black", linewidth = 0.4) +
    ggplot2::facet_grid(sex ~ harvesttimeSK24) +
    ggplot2::coord_sf(crs = sf::st_crs(study_region)) +
    ggplot2::scale_fill_viridis_c(name = "Predicted probability", na.value = NA) +
    ggplot2::scale_color_manual(name = "Points", values = c("random = 0" = "grey", "random = 1" = "black")[levels(pts_sf_all$point_type)]) +
    ggplot2::labs(
      title = "Predicted habitat suitability (sex × harvest period)",
      subtitle = "Population-level maps (ID random effect excluded in prediction)",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal()

  list(
    by_period = res_by_period,                        # list of per-period results from predict_maps_by_sex
    pred_raster_sex_period = pred_raster_by_sex_period, # SpatRaster with one layer per sex × period
    prediction_plot_sex_period = prediction_plot_sex_period
  )
}
