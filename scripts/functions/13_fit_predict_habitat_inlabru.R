library(INLA)
library(inlabru)
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(tidyr)

fit_predict_habitat_inlabru_mesh <- function(
    dmodel,
    pred_grid,
    env_stack_10m,
    raster_results_23,
    raster_results_24,
    pars_to_use,
    response_col = "used",
    use_spatial = TRUE,
    mesh = NULL,                # optional pre-built mesh (inla.mesh.2d); if NULL function will build one
    mesh_args = list(
      max.edge = c(200, 800),   # outer triangle sizes (meters) — tune for your domain
      cutoff   = 50,            # minimum distance between points
      offset   = c(500, 1000)   # extension of outer boundary
    ),
    spde_priors = list(
      prior.range = c(1000, 0.05),  # P(range < 1000 m) = 0.05 (tune to your system)
      prior.sigma = c(1, 0.05)      # P(sigma > 1) = 0.05
    ),
    n_samples_pred = 500,
    control_compute = list(waic = TRUE, cpo = TRUE),
    verbose = TRUE
) {
  # -------------------------
  # 0) checks
  # -------------------------
  pkgs <- c("sf", "terra", "dplyr", "ggplot2", "tidyr", "INLA", "inlabru")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))

  `%>%` <- dplyr::`%>%`

  # -------------------------
  # 1) Filter model data for year and extract rasters (same extraction logic you provided)
  # -------------------------
  model_sf <- dmodel %>%
    dplyr::filter(Year == pars_to_use$Year) %>%
    dplyr::mutate(
      .rowid = dplyr::row_number(),
      geom_wkt = sf::st_as_text(geometry)
    ) %>%
    dplyr::distinct(geom_wkt, .keep_all = TRUE) %>%
    dplyr::select(-geom_wkt)

  if (pars_to_use$Year == "2023") {
    harvest_stack <- raster_results_23$raster_stack
    harvest_regions <- raster_results_23$harvest50_buffer
    harvest200_regions <- raster_results_23$harvest200_buffer
  } else if (pars_to_use$Year == "2024") {
    harvest_stack <- raster_results_24$raster_stack
    harvest_regions <- raster_results_24$harvest50_buffer
    harvest200_regions <- raster_results_24$raster200_buffer
  } else {
    stop("pars_to_use$Year must be '2023' or '2024'.")
  }

  harvest_regions_sf <- sf::st_sf(
    harvested50m = 1L,
    geometry = harvest_regions,
    crs = sf::st_crs(harvest_regions)
  )

  allRasts <- c(harvest_stack, env_stack_10m)

  r_crs <- terra::crs(allRasts, proj = TRUE)
  if (is.na(r_crs) || r_crs == "") stop("allRasts has no CRS. Set a projected CRS with meter units.")
  model_sf_r <- sf::st_transform(model_sf, r_crs)

  env_values_model <- terra::extract(allRasts, terra::vect(model_sf_r)) %>%
    dplyr::rename(.rowid = ID) %>%
    dplyr::select(-dplyr::any_of(c("harvested50m", "harvested200m")))

  model_sf <- model_sf %>% dplyr::select(-dplyr::any_of("DistNest"))
  model_sf <- model_sf %>% dplyr::left_join(env_values_model, by = ".rowid")

  pred_grid_sf <- sf::st_as_sf(pred_grid)
  pred_grid_r  <- sf::st_transform(pred_grid_sf, r_crs)
  env_values_pred <- terra::extract(allRasts, terra::vect(pred_grid_r)) %>% dplyr::select(-ID)
  pred_grid_env <- pred_grid_r %>% dplyr::bind_cols(env_values_pred)

  # -------------------------
  # 2) Prepare data.frame for modeling (assumes DistNest and bathymetry are already scaled)
  # -------------------------
  if (!"w" %in% names(model_sf)) model_sf$w <- 1

  if (!response_col %in% names(model_sf)) stop("response_col not found in model data")
  if (is.factor(model_sf[[response_col]])) {
    model_sf[[response_col]] <- as.numeric(as.character(model_sf[[response_col]]))
  }
  if (!all(model_sf[[response_col]] %in% c(0, 1))) stop("response_col must contain only 0/1 values")

  need_cols <- c("DistNest", "DistHarvest", "bathymetry")
  miss_cols <- setdiff(need_cols, names(model_sf))
  if (length(miss_cols)) stop("Missing predictor(s) in model data after extraction: ", paste(miss_cols, collapse = ", "))

  if (!"sex" %in% names(model_sf)) model_sf$sex <- factor("Unknown")
  if (!"harvesttimeSK24" %in% names(model_sf)) model_sf$harvesttimeSK24 <- factor("all")
  if (!"ID" %in% names(model_sf)) model_sf$ID <- factor("all")

  model_sf_clean <- model_sf %>%
    dplyr::select(ID, sex, harvesttimeSK24, Year,
                  dplyr::all_of(response_col),
                  DistNest, DistHarvest, bathymetry, w, geometry) %>%
    tidyr::drop_na()

  model_sf_clean$w <- model_sf_clean$w / mean(model_sf_clean$w)

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

  # Check that DistNest & bathymetry are already scaled (user said they are)
  # If not scaled you can scale here, but we respect user's claim and use as-is.
  # Optionally allow pred_grid_env scaled columns: assume they are present and consistent.

  if (!("DistNest" %in% names(pred_grid_env) && "bathymetry" %in% names(pred_grid_env))) {
    warning("Prediction grid missing DistNest/bathymetry columns — full predictions may be incomplete.")
  }

  # -------------------------
  # 3) Build mesh & SPDE (if requested)
  # -------------------------
  spde <- NULL
  if (use_spatial) {
    if (is.null(mesh)) {
      if (verbose) message("Building mesh from model coordinates (tune mesh_args for your domain)...")
      mesh <- inla.mesh.2d(
        loc = cbind(model_data$x, model_data$y),
        max.edge = mesh_args$max.edge,
        cutoff   = mesh_args$cutoff,
        offset   = mesh_args$offset
      )
    } else {
      if (verbose) message("Using provided mesh object.")
    }

    # Build SPDE with PC priors
    spde <- inla.spde2.pcmatern(
      mesh = mesh,
      alpha = 2,
      prior.range = spde_priors$prior.range,
      prior.sigma = spde_priors$prior.sigma
    )
  }

  # -------------------------
  # 4) Build inlabru components
  # -------------------------
  if (use_spatial && !is.null(spde)) {
    components <- ~
      Intercept(1) +
      DistHarvest +
      DistNest +        # continuous, already scaled
      bathymetry +      # continuous, already scaled
      sex +
      harvesttimeSK24 +
      DistHarvest:sex +
      DistHarvest:harvesttimeSK24 +
      f(ID, model = "iid") +
      field(coordinates, model = spde)
  } else {
    components <- ~
      Intercept(1) +
      DistHarvest +
      DistNest +
      bathymetry +
      sex +
      harvesttimeSK24 +
      DistHarvest:sex +
      DistHarvest:harvesttimeSK24 +
      f(ID, model = "iid")
  }

  # Add coordinates column used by field(...) and some predict routines
  model_data$coordinates <- cbind(model_data$x, model_data$y)
  if ("x" %in% names(pred_grid_env) && "y" %in% names(pred_grid_env)) {
    pred_grid_env$coordinates <- cbind(pred_grid_env$x, pred_grid_env$y)
  }

  # -------------------------
  # 5) Define likelihood and fit with bru
  # -------------------------
  like_formula_str <- paste(
    "response ~ Intercept + DistHarvest + DistNest + bathymetry +",
    "sex + harvesttimeSK24 + DistHarvest:sex + DistHarvest:harvesttimeSK24 + f(ID)"
  )

  lnk <- like(
    formula = as.formula(like_formula_str),
    family  = "binomial",
    data    = model_data,
    weights = model_data$w
  )

  if (verbose) message("Fitting with inlabru (this may take some time)...")
  fit_bru <- bru(
    components = components,
    lnk,
    options = list(
      control.compute = control_compute,
      bru_verbose = FALSE
    )
  )

  # -------------------------
  # 6) Predictions
  #    - full prediction including spatial field and ID (for maps): expr_full
  #    - marginal prediction excluding field & ID (for partial effects): expr_marg
  # -------------------------
  pred_out <- NULL
  if ("DistNest" %in% names(pred_grid_env) && "bathymetry" %in% names(pred_grid_env) &&
      "DistHarvest" %in% names(pred_grid_env)) {

    # Full expression (includes field and f(ID) if present). For pred locations ID is usually NA -> f(ID) ignored
    expr_full <- if (use_spatial && !is.null(spde)) {
      ~ link(Intercept + DistHarvest + DistNest + bathymetry +
               sex + harvesttimeSK24 + DistHarvest:sex + DistHarvest:harvesttimeSK24 +
               field(coordinates))
    } else {
      ~ link(Intercept + DistHarvest + DistNest + bathymetry +
               sex + harvesttimeSK24 + DistHarvest:sex + DistHarvest:harvesttimeSK24)
    }

    pred_full <- predict(fit_bru, pred_grid_env, expr_full, n.samples = n_samples_pred)

    pred_grid_env$prob_mean <- pred_full$mean
    pred_grid_env$prob_sd   <- pred_full$sd
    pred_grid_env$prob_q025 <- pred_full$q025
    pred_grid_env$prob_q975 <- pred_full$q975

    pred_out <- pred_grid_env
  }

  # -------------------------
  # 7) Partial effect for DistHarvest excluding spatial field and ID
  # -------------------------
  qDH <- quantile(model_data$DistHarvest, c(0.02, 0.98), na.rm = TRUE)
  DH_seq <- seq(qDH[1], qDH[2], length.out = 150)
  lev_period <- levels(model_data$harvesttimeSK24)
  lev_sex    <- levels(model_data$sex)

  med_DistNest <- median(model_data$DistNest, na.rm = TRUE)
  med_bath     <- median(model_data$bathymetry, na.rm = TRUE)

  effect_grid <- do.call(rbind, lapply(lev_sex, function(sx) {
    do.call(rbind, lapply(lev_period, function(pr) {
      nd <- data.frame(
        DistHarvest = DH_seq,
        DistNest = med_DistNest,
        bathymetry = med_bath,
        sex = factor(sx, levels = lev_sex),
        harvesttimeSK24 = factor(pr, levels = lev_period),
        coordinates = cbind(mean(model_data$x, na.rm = TRUE), mean(model_data$y, na.rm = TRUE))
      )
      nd
    }))
  }))

  # Expression excluding field and ID
  expr_eff <- ~ Intercept + DistHarvest + DistNest + bathymetry +
    sex + harvesttimeSK24 + DistHarvest:sex + DistHarvest:harvesttimeSK24

  pred_eff <- predict(fit_bru, effect_grid, ~ link(expr_eff), n.samples = n_samples_pred)

  effect_grid$eta <- pred_eff$mean
  effect_grid$se  <- pred_eff$sd
  effect_grid$p   <- plogis(effect_grid$eta)
  effect_grid$p_lo <- plogis(effect_grid$eta - 1.96 * effect_grid$se)
  effect_grid$p_hi <- plogis(effect_grid$eta + 1.96 * effect_grid$se)

  effect_grid <- effect_grid %>%
    group_by(sex, harvesttimeSK24) %>%
    mutate(
      eta_ctr = eta - mean(eta, na.rm = TRUE),
      RS = exp(eta_ctr)
    ) %>%
    ungroup()

  # -------------------------
  # 8) Plots
  # -------------------------
  dist_harvest_effect_plot <- ggplot(effect_grid, aes(x = DistHarvest, y = p, color = harvesttimeSK24, fill = harvesttimeSK24)) +
    geom_ribbon(aes(ymin = p_lo, ymax = p_hi), alpha = 0.15, color = NA) +
    geom_line(size = 1) +
    facet_wrap(~ sex, nrow = 1) +
    labs(title = paste0("Partial effect of DistHarvest (year ", pars_to_use$Year, ")"),
         subtitle = "Predicted probability (field & ID excluded)",
         x = "Distance to harvest", y = "Predicted probability") +
    theme_minimal() + scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2")

  rel_selection_plot <- ggplot(effect_grid, aes(x = DistHarvest, y = RS, color = harvesttimeSK24)) +
    geom_line(size = 1) +
    facet_wrap(~ sex, nrow = 1) +
    labs(title = paste0("Relative selection vs DistHarvest (year ", pars_to_use$Year, ")"),
         subtitle = "RS = exp(centered link); field & ID excluded",
         x = "Distance to harvest", y = "Relative selection (RS)") +
    theme_minimal() + scale_color_brewer(palette = "Dark2")

  # -------------------------
  # 9) Return list
  # -------------------------
  out <- list(
    method = "inlabru_spde",
    parameters = pars_to_use,
    model = fit_bru,
    model_data = model_data,
    model_sf = model_sf_clean,
    pred_grid_env = pred_out,
    dist_harvest_effect = effect_grid,
    dist_harvest_effect_plot = dist_harvest_effect_plot,
    rel_selection_plot = rel_selection_plot,
    waic = fit_bru$internal$waic,
    cpo  = fit_bru$internal$cpo,
    mesh = if (use_spatial) mesh else NULL,
    spde = spde,
    harvest_regions = harvest_regions_sf
  )

  return(out)
}
