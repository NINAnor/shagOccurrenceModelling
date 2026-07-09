run_shag_inla_model <- function(
    pars_to_use,
    dmodel,
    study_region,
    path_for_processed_data,
    modelCrs,
    thinning_distance = 1000
) {

  library(sf)
  library(terra)
  library(dplyr)
  library(fmesher)
  library(INLA)
  library(inlabru)
  library(ggplot2)

  message("Loading covariates...")

  if (pars_to_use$Year == "2023") {

    proj_file <- file.path(
      path_for_processed_data,
      "cov_stack_proj_2023.tif"
    )

    if (!file.exists(proj_file)) {

      cov_stack <- terra::rast(c(
        file.path(path_for_processed_data,
                  "raster_stack_2023.tif"),
        file.path(path_for_processed_data,
                  "env_stack_10m.tif")
      ))

      names(cov_stack)

      cov_stack <- cov_stack[[c(
        "DistNest",
        "DistHarvest",
        "bathymetry"
      )]]

      terra::project(
        cov_stack,
        paste0("EPSG:", modelCrs),
        filename = proj_file,
        overwrite = TRUE
      )
    }

    cov_stack <- terra::rast(proj_file)

    load("../SK_35_22_23")

    harvest_region <- SK_35km_22_23 |>
      sf::st_transform(modelCrs) |>
      dplyr::select(geometry) |>
      sf::st_crop(
        study_region |>
          sf::st_transform(modelCrs)
      )

  } else {

    proj_file <- file.path(
      path_for_processed_data,
      "cov_stack_proj_2024.tif"
    )

    if (!file.exists(proj_file)) {

      cov_stack <- terra::rast(c(
        file.path(path_for_processed_data,
                  "raster_stack_2024.tif"),
        file.path(path_for_processed_data,
                  "env_stack_10m.tif")
      ))

      cov_stack <- cov_stack[[c(
        "DistNest",
        "DistHarvest",
        "bathymetry"
      )]]

      terra::project(
        cov_stack,
        paste0("EPSG:", modelCrs),
        filename = proj_file,
        overwrite = TRUE
      )
    }

    cov_stack <- terra::rast(proj_file)

    load("../SK_35_23_24")

    harvest_region <- SK_35km_23_24 |>
      sf::st_transform(modelCrs) |>
      dplyr::select(geometry) |>
      sf::st_crop(
        study_region |>
          sf::st_transform(modelCrs)
      )
  }

  message("Preparing presence data...")

  pres_sf <- dmodel |>
    dplyr::filter(
      Year == pars_to_use$Year,
      sex == pars_to_use$sex,
      harvesttimeSK24 == pars_to_use$harvesttimeSK24,
      random == 1
    )


  pres_sf_unique <- pres_sf

  coords <- sf::st_coordinates(pres_sf_unique)

  pres_sf_unique$X <- coords[, 1]
  pres_sf_unique$Y <- coords[, 2]

  pres_sf_unique <- pres_sf_unique |>
    dplyr::distinct(X, Y, .keep_all = TRUE) |>
    dplyr::select(-X, -Y)


  keep <- rep(TRUE, nrow(pres_sf_unique))

  for (i in seq_len(nrow(pres_sf_unique) - 1)) {

    if (!keep[i]) next

    d <- sf::st_distance(
      pres_sf_unique[i, ],
      pres_sf_unique[(i + 1):nrow(pres_sf_unique), ]
    )

    close <- which(
      as.numeric(d) < thinning_distance
    )

    if (length(close) > 0) {
      keep[(i + 1):nrow(pres_sf_unique)][close] <- FALSE
    }
  }

  pres_sf_thinned <- pres_sf_unique[keep, ]

  message(
    "Retained ",
    nrow(pres_sf_thinned),
    " locations after thinning."
  )

  mesh_builder <- create_fmesher_mesh_from_locations(
    points_sf = pres_sf_thinned,
    study_region_sf = study_region,
    target_crs = modelCrs,
    offset = c(2000, 10000),
    max_edge = c(3000, 10000),
    cutoff = 1000
  )

  mesh <- mesh_builder$mesh
  region_sf <- mesh_builder$study_region
  pres_sf <- mesh_builder$used_points

  spde <- inla.spde2.pcmatern(
    mesh = mesh,
    alpha = 2,
    prior.range = c(3000, 0.5),
    prior.sigma = c(0.5, 0.01)
  )

  cov_env <- list(
    DistNest = cov_stack$DistNest,
    DistHarvest = cov_stack$DistHarvest,
    bathymetry = cov_stack$bathymetry
  )

  expected_intercept_m2 <-
    log(
      nrow(pres_sf) /
        as.numeric(sf::st_area(region_sf))
    )

  pres_sf_model <- pres_sf |>
    dplyr::select(random, geometry)

  cmp_spatial <-
    ~ Intercept(1) +
    DistNest(main = DistNest) +
    DistHarvest(main = DistHarvest) +
    bathymetry(main = bathymetry) +
    spatial(geometry, model = spde)

  cmp_nospatial <-
    ~ Intercept(1) +
    DistNest(main = DistNest) +
    DistHarvest(main = DistHarvest) +
    bathymetry(main = bathymetry)

  fit_spatial <- bru(
    components = cmp_spatial,
    bru_obs(
      formula = geometry ~ .,
      family = "cp",
      data = pres_sf_model,
      domain = list(geometry = mesh),
      samplers = list(geometry = region_sf),
      data_extra = cov_env
    ),
    options = bru_options(
      control.fixed = list(
        mean.intercept = expected_intercept_m2,
        prec.intercept = 1,
        mean = 0,
        prec = 1
      ),
      control.compute = list(
        dic = TRUE,
        waic = TRUE
      )
    )
  )

  fit_nospatial <- bru(
    components = cmp_nospatial,
    bru_obs(
      formula = geometry ~ .,
      family = "cp",
      data = pres_sf_model,
      domain = list(geometry = mesh),
      samplers = list(geometry = region_sf),
      data_extra = cov_env
    ),
    options = bru_options(
      control.compute = list(
        dic = TRUE,
        waic = TRUE
      )
    )
  )

  models <- list(
    spatial = fit_spatial,
    nonspatial = fit_nospatial
  )

  waics <- sapply(
    models,
    function(x) x$waic$waic
  )

  best_model <- models[[which.min(waics)]]

  px <- fm_pixels(
    mesh,
    dim = c(500, 500),
    mask = region_sf
  )

  pred_data <- px

  pred_data$DistNest <-
    terra::extract(
      cov_stack$DistNest,
      terra::vect(pred_data)
    )[, 2]

  pred_data$DistHarvest <-
    terra::extract(
      cov_stack$DistHarvest,
      terra::vect(pred_data)
    )[, 2]

  pred_data$bathymetry <-
    terra::extract(
      cov_stack$bathymetry,
      terra::vect(pred_data)
    )[, 2]

  pred_dh <- predict(
    best_model,
    pred_data,
    ~ exp(Intercept + DistHarvest)
  )

  pred_full <- predict(
    best_model,
    pred_data,
    ~ exp(
      Intercept +
        DistNest +
        DistHarvest +
        bathymetry
    )
  )

  plot_dh <-
    ggplot(pred_dh) +
    geom_sf(aes(colour = mean)) +
    geom_sf(
      data = harvest_region,
      colour = "black"
    ) +
    geom_sf(
      data = pres_sf,
      colour = "#D55E00",
      size = 2
    ) +
    scale_colour_viridis_c(
      option = "cividis",
      name = "Predicted\nrelative\nintensity"
    ) +
    theme_bw()

  plot_full <-
    ggplot(pred_full) +
    geom_sf(aes(colour = mean)) +
    geom_sf(
      data = harvest_region,
      colour = "black"
    ) +
    geom_sf(
      data = pres_sf,
      colour = "#D55E00",
      size = 2
    ) +
    scale_colour_viridis_c(
      option = "cividis",
      name = "Predicted\nrelative\nintensity"
    ) +
    theme_bw()

  return(
    list(
      fit_spatial = fit_spatial,
      fit_nospatial = fit_nospatial,
      best_model = best_model,
      waic = waics,
      mesh = mesh,
      presences = pres_sf,
      pred_dh = pred_dh,
      pred_full = pred_full,
      plot_dh = plot_dh,
      plot_full = plot_full,
      method = "inla"
    )
  )
}



