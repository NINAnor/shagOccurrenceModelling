#' Habitat-selection LGCP (v2): 35 km study region, kelp as covariate,
#' time-specific harvest covariate
#'
#' Differences from `run_shag_inla_model()`:
#'  * Available area = colony buffer (default 35 km) restricted to water, not the
#'    kelp polygons. Kelp enters as a covariate (`tareskog`, fraction of cell).
#'  * Presences = all real dives (kelp and non-kelp), not the kelp-clipped dmodel.
#'  * Integration points come from a raster grid (cell centres, weight = km2),
#'    so 10-100 m scale kelp/harvest structure is not lost on a coarse mesh.
#'  * Harvest covariate is the fraction of each cell within `harvest_buffer_m`
#'    of a harvest polygon with GIS_TIME earlier than a reference time for the
#'    harvest period (i.e. the raster analogue of `harvested200m/50m`).
#'  * No ID / logger iid terms: they are not defined at integration points, so
#'    they would enter the data term but not the integral (see Details).
#'
#' @details
#' Water = GEBCO elevation < 0 OR kelp present in the cell (GEBCO is ~450 m, so
#' a few dives fall on "land" cells; those with no kelp are dropped).
#' Reference time: start of the period (for "before" = first dive of that year).
#'
#' @name shag_inla_v2
NULL

# Period windows used throughout the project (see 04_formatDataForHarvestTime.R)
shag_period_windows <- function(year, tz = "UTC") {
  list(
    start = as.POSIXct(paste0(year, "-06-17 11:00:00"), tz = tz),
    end = as.POSIXct(paste0(year, "-06-27 16:00:00"), tz = tz)
  )
}

#' Build the covariate grid (water mask, kelp, depth, distance to colony)
#' @export
build_inla_grid <- function(
  colony,
  tareskog,
  gebco_path,
  modelCrs = 32633,
  buffer_m = 35000,
  res = 200
) {
  colony_xy <- sf::st_coordinates(sf::st_centroid(sf::st_union(sf::st_geometry(
    colony
  ))))[1, ]
  region <- sf::st_sfc(
    sf::st_buffer(sf::st_point(colony_xy), buffer_m),
    crs = modelCrs
  ) |>
    sf::st_sf()

  bb <- as.numeric(sf::st_bbox(region))
  bb <- c(
    floor(bb[1] / res),
    ceiling(bb[3] / res),
    floor(bb[2] / res),
    ceiling(bb[4] / res)
  ) *
    res
  template <- terra::rast(
    terra::ext(bb),
    resolution = res,
    crs = paste0("EPSG:", modelCrs)
  )

  region_v <- terra::vect(region)
  in_region <- terra::rasterize(region_v, template, field = 1, background = NA)

  # Elevation (m, negative = sea)
  elev <- terra::rast(gebco_path) |>
    terra::crop(terra::ext(terra::project(region_v, "EPSG:4326")) + 0.05) |>
    terra::project(template, method = "bilinear")

  # Kelp fraction per cell
  kelp <- terra::aggregate(
    tareskog,
    fact = res / terra::res(tareskog)[1],
    fun = "mean",
    na.rm = TRUE
  ) |>
    terra::project(template, method = "bilinear")
  kelp <- terra::ifel(is.na(kelp), 0, kelp)
  kelp <- terra::clamp(kelp, 0, 1)

  water <- (elev < 0) | (kelp > 0)
  water <- terra::mask(water, in_region)
  water <- terra::ifel(water, 1, NA)

  # Covariates (z-scored over water cells in the region)
  depth <- log1p(terra::clamp(-elev, lower = 0))
  xy <- terra::xyFromCell(template, seq_len(terra::ncell(template)))
  dist_km <- sqrt((xy[, 1] - colony_xy[1])^2 + (xy[, 2] - colony_xy[2])^2) /
    1000
  dist <- terra::setValues(template, log1p(dist_km))

  zs <- function(r) {
    v <- terra::values(terra::mask(r, water), na.rm = TRUE)
    list(center = mean(v), scale = stats::sd(v))
  }
  sc <- list(depth = zs(depth), dist = zs(dist))
  z <- function(r, s) (r - s$center) / s$scale

  depth_z <- z(depth, sc$depth)
  dist_z <- z(dist, sc$dist)

  grid <- c(kelp, depth_z, depth_z^2, dist_z, water)
  names(grid) <- c("kelp", "depth", "depth2", "dist", "water")

  list(
    region = region,
    grid = grid,
    template = template,
    scaling = sc,
    colony_xy = colony_xy,
    res = res,
    modelCrs = modelCrs
  )
}

#' Fraction of each cell within `buffer_m` of harvest polygons harvested before `t_ref`
#' @export
make_harvest_cover <- function(
  harvest,
  t_ref,
  buffer_m,
  template,
  year = NULL,
  year_only = TRUE
) {
  # Accept the per-year list used by run_shag_inla_v2_all() as well
  if (!inherits(harvest, "sf") && is.list(harvest) && !is.null(year)) {
    harvest <- harvest[[as.character(year)]]
  }
  if (!inherits(harvest, "sf")) {
    stop(
      "`harvest` must be an sf object with a GIS_TIME column (got ",
      paste(class(harvest), collapse = "/"),
      "). Pass e.g. harv_list[[\"2024\"]], not the whole list or NULL."
    )
  }
  h <- harvest
  if (year_only && !is.null(year) && "year" %in% names(h)) {
    h <- h[h$year == as.integer(year), ]
  }
  h <- h[!is.na(h$GIS_TIME) & h$GIS_TIME < t_ref, ]
  if (nrow(h) == 0) {
    return(terra::setValues(template, 0))
  }

  h <- sf::st_transform(h, terra::crs(template)) |>
    sf::st_geometry() |>
    sf::st_buffer(buffer_m) |>
    sf::st_union() |>
    sf::st_make_valid()

  r <- terra::rasterize(
    terra::vect(sf::st_sf(geometry = h)),
    template,
    cover = TRUE,
    background = 0
  )
  terra::clamp(r, 0, 1)
}

#' Prepare raw dives: unique locations per individual, year, period
#' @export
prepare_dives <- function(dives, modelCrs = 32633) {
  dives <- sf::st_transform(dives, modelCrs)
  xy <- sf::st_coordinates(dives)
  dives$X <- xy[, 1]
  dives$Y <- xy[, 2]
  dives |>
    dplyr::mutate(
      Year = as.character(Year),
      .ws = as.POSIXct(paste0(Year, "-06-17 11:00:00"), tz = "UTC"),
      .we = as.POSIXct(paste0(Year, "-06-27 16:00:00"), tz = "UTC"),
      period = dplyr::case_when(
        UTC_datetime < .ws ~ "beforeJuneHarvest",
        UTC_datetime <= .we ~ "duringJuneHarvest",
        TRUE ~ "afterJuneHarvest"
      )
    ) |>
    dplyr::select(-.ws, -.we) |>
    dplyr::distinct(ID, Logger_model, Year, X, Y, .keep_all = TRUE)
}

#' Fit one Year x sex x period LGCP and return the harvest coefficient
#'
#' @param pars_to_use list/data.frame row with Year, sex, harvesttimeSK24
#' @param dives output of `prepare_dives()` (all years)
#' @param harvest sf harvest polygons for the year (needs GIS_TIME, year)
#' @param cov output of `build_inla_grid()`
#' @param harvest_buffer_m 200 (default) or 50 etc.
#' @param spatial include the SPDE field
#' @export
run_shag_inla_model_v2 <- function(
  pars_to_use,
  dives,
  harvest,
  cov,
  harvest_buffer_m = 50,
  harvest_year_only = TRUE,
  kelp_type = c("presence", "fraction"),
  spatial = TRUE,
  mesh_max_edge = c(2500, 10000),
  mesh_offset = c(1000, 5000),
  verbose = FALSE
) {
  library(sf)
  library(terra)
  library(dplyr)
  library(fmesher)
  library(INLA)
  library(inlabru)

  year <- as.character(pars_to_use$Year)
  sex <- as.character(pars_to_use$sex)
  period <- as.character(pars_to_use$harvesttimeSK24)

  # Reference time = start of period
  win <- shag_period_windows(year)
  t_ref <- switch(
    period,
    beforeJuneHarvest = min(dives$UTC_datetime[dives$Year == year]),
    duringJuneHarvest = win$start,
    afterJuneHarvest = win$end
  )

  # Covariates on the grid for this stratum
  harv <- make_harvest_cover(
    harvest,
    t_ref,
    harvest_buffer_m,
    cov$template,
    year = year,
    year_only = harvest_year_only
  )
  grid <- c(cov$grid, harv)
  names(grid)[terra::nlyr(grid)] <- "harvest"

  # Kelp: presence (any kelp in the cell) or fraction of cell covered
  kelp_type <- match.arg(kelp_type)
  if (kelp_type == "presence") {
    grid$kelp <- terra::ifel(grid$kelp > 0, 1, 0)
  }

  # Integration points: water cells inside the region
  cells <- which(!is.na(terra::values(grid$water)[, 1]))
  xy <- terra::xyFromCell(grid, cells)
  vals <- terra::values(grid)[
    cells,
    c("kelp", "depth", "depth2", "dist", "harvest"),
    drop = FALSE
  ]
  cell_km2 <- prod(terra::res(grid)) / 1e6
  ips <- sf::st_as_sf(
    data.frame(x = xy[, 1], y = xy[, 2], vals, weight = cell_km2, .block = 1L),
    coords = c("x", "y"),
    crs = cov$modelCrs
  )

  # Presences (water only)
  pres <- dives |>
    dplyr::filter(Year == year, sex == !!sex, period == !!period)
  n_all <- nrow(pres)
  pv <- terra::extract(grid, cbind(pres$X, pres$Y))
  keep <- !is.na(pv$water)
  pres <- pres[keep, ]
  pv <- pv[keep, c("kelp", "depth", "depth2", "dist", "harvest")]
  pres <- dplyr::bind_cols(
    sf::st_sf(geometry = sf::st_geometry(pres)),
    pv
  )
  message(
    "Dives: ",
    n_all,
    " unique; ",
    nrow(pres),
    " on water cells in region."
  )
  if (nrow(pres) < 50) {
    stop("Too few dives in stratum: ", nrow(pres))
  }

  # Mesh and SPDE
  region_simple <- sf::st_simplify(cov$region, dTolerance = 500)
  mesh <- fmesher::fm_mesh_2d(
    boundary = fmesher::fm_as_segm(region_simple),
    max.edge = mesh_max_edge,
    offset = mesh_offset,
    cutoff = 500,
    crs = sf::st_crs(cov$region)
  )
  priors <- choose_pc_priors_binomial(cov$region)
  spde <- INLA::inla.spde2.pcmatern(
    mesh = mesh,
    alpha = 2,
    prior.range = priors$prior.range,
    prior.sigma = priors$prior.sigma
  )

  area_km2 <- sum(ips$weight)
  cmp <- ~ Intercept(1) +
    kelp(kelp, model = "linear") +
    harvest(harvest, model = "linear") +
    depth(depth, model = "linear") +
    depth2(depth2, model = "linear") +
    dist(dist, model = "linear")
  if (spatial) {
    cmp <- update(cmp, ~ . + field(geometry, model = spde))
  }

  fit <- inlabru::bru(
    components = cmp,
    inlabru::bru_obs(
      formula = geometry ~ .,
      family = "cp",
      data = pres,
      ips = ips
    ),
    options = inlabru::bru_options(
      verbose = verbose,
      control.compute = list(dic = TRUE, waic = TRUE),
      control.fixed = list(
        mean = 0,
        prec = 1,
        mean.intercept = log(nrow(pres) / area_km2),
        prec.intercept = 0.1
      )
    )
  )

  fx <- fit$summary.fixed
  coefs <- data.frame(
    Year = year,
    sex = sex,
    period = period,
    term = rownames(fx),
    mean = fx$mean,
    sd = fx$sd,
    q025 = fx$`0.025quant`,
    q975 = fx$`0.975quant`,
    row.names = NULL
  )
  # Posterior probability that each coefficient is negative
  coefs$p_neg <- vapply(
    coefs$term,
    function(tm) {
      m <- fit$marginals.fixed[[tm]]
      if (is.null(m)) NA_real_ else INLA::inla.pmarginal(0, m)
    },
    numeric(1)
  )
  list(
    fit = fit,
    coefs = coefs,
    mesh = mesh,
    n_dives = nrow(pres),
    t_ref = t_ref,
    harvest_buffer_m = harvest_buffer_m,
    harvest_cover_at_dives = mean(pres$harvest > 0),
    harvest_cover_area = mean(ips$harvest > 0),
    waic = fit$waic$waic
  )
}

#' Fit all Year x sex x period strata and collect harvest coefficients
#' @export
run_shag_inla_v2_all <- function(
  dives,
  harvest_by_year,
  cov,
  harvest_buffer_m = 50,
  out_dir = NULL,
  ...
) {
  strata <- expand.grid(
    Year = c("2023", "2024"),
    sex = c("Female", "Male"),
    harvesttimeSK24 = c(
      "beforeJuneHarvest",
      "duringJuneHarvest",
      "afterJuneHarvest"
    ),
    stringsAsFactors = FALSE
  )
  res <- lapply(seq_len(nrow(strata)), function(i) {
    p <- strata[i, ]
    message("Stratum ", i, "/", nrow(strata), ": ", paste(p, collapse = " / "))
    sdir <- if (!is.null(out_dir)) {
      file.path(out_dir, p$Year, p$sex, p$harvesttimeSK24)
    }
    sfile <- if (!is.null(sdir)) file.path(sdir, "fit_summary.rds")
    # Resume: reuse stratum if already saved
    if (!is.null(sfile) && file.exists(sfile)) {
      message("  already saved, loading.")
      return(readRDS(sfile))
    }
    r <- tryCatch(
      run_shag_inla_model_v2(
        p,
        dives,
        harvest_by_year[[p$Year]],
        cov,
        harvest_buffer_m = harvest_buffer_m,
        ...
      ),
      error = function(e) {
        message("  failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(r) && !is.null(sdir)) {
      dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
      light <- list(
        coefs = r$coefs,
        hyper = r$fit$summary.hyperpar,
        n_dives = r$n_dives,
        t_ref = r$t_ref,
        harvest_buffer_m = r$harvest_buffer_m,
        harvest_cover_at_dives = r$harvest_cover_at_dives,
        harvest_cover_area = r$harvest_cover_area,
        waic = r$waic
      )
      saveRDS(light, sfile)
      utils::write.csv(r$coefs, file.path(sdir, "coefs.csv"), row.names = FALSE)
      saveRDS(r$fit, file.path(sdir, "fit_full.rds"))
      r <- light
    }
    r
  })
  names(res) <- do.call(paste, c(strata, sep = "_"))
  coefs <- dplyr::bind_rows(lapply(res, function(r) if (!is.null(r)) r$coefs))
  list(results = res, coefs = coefs)
}

#' Tabulate and plot the harvest coefficient from `run_shag_inla_v2_all()`
#'
#' Classification: "avoids" if the 95% CrI is entirely below 0, "prefers" if
#' entirely above 0, otherwise "no clear effect". Strata with < min_exposure of
#' dives inside harvested cells are flagged "not estimable".
#' @export
summarise_harvest_effect <- function(
  all_res,
  out_dir,
  buffer_m = 50,
  min_exposure = 0.01
) {
  library(ggplot2)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tab <- dplyr::bind_rows(lapply(names(all_res$results), function(nm) {
    r <- all_res$results[[nm]]
    if (is.null(r)) {
      return(NULL)
    }
    h <- r$coefs[r$coefs$term == "harvest", ]
    h$n_dives <- r$n_dives
    h$pct_dives_in_harvest <- 100 * r$harvest_cover_at_dives
    h$pct_area_in_harvest <- 100 * r$harvest_cover_area
    h
  }))
  tab$period <- factor(
    tab$period,
    levels = c("beforeJuneHarvest", "duringJuneHarvest", "afterJuneHarvest")
  )
  tab$exp_mean <- exp(tab$mean)
  tab$conclusion <- dplyr::case_when(
    tab$pct_dives_in_harvest / 100 <
      min_exposure ~ "not estimable (no exposure)",
    tab$q975 < 0 ~ "avoids harvested areas",
    tab$q025 > 0 ~ "prefers harvested areas",
    TRUE ~ "no clear effect"
  )
  utils::write.csv(
    tab,
    file.path(out_dir, "harvest_effect_table.csv"),
    row.names = FALSE
  )
  saveRDS(tab, file.path(out_dir, "harvest_effect_table.rds"))

  plt <- tab[!grepl("not estimable", tab$conclusion), ]
  skipped <- tab[grepl("not estimable", tab$conclusion), ]
  p <- ggplot(plt, aes(x = mean, y = period, colour = conclusion)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(xmin = q025, xmax = q975)) +
    facet_grid(Year ~ sex) +
    scale_colour_manual(
      values = c(
        "avoids harvested areas" = "#D55E00",
        "prefers harvested areas" = "#0072B2",
        "no clear effect" = "grey40"
      )
    ) +
    labs(
      x = paste0(
        "Harvest coefficient (log scale, ",
        buffer_m,
        " m buffer; 95% credible interval)"
      ),
      y = NULL,
      colour = NULL,
      caption = if (nrow(skipped) > 0) {
        paste(
          strwrap(
            paste0(
              "Not shown, not estimable (<",
              100 * min_exposure,
              "% of dives in harvested areas): ",
              paste(unique(skipped$Year), collapse = ", "),
              " (",
              nrow(skipped),
              " strata)."
            ),
            width = 90
          ),
          collapse = "\n"
        )
      }
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(
    file.path(out_dir, "harvest_effect_forest.png"),
    p,
    width = 9,
    height = 6,
    dpi = 300
  )
  ggsave(
    file.path(out_dir, "harvest_effect_forest.pdf"),
    p,
    width = 9,
    height = 6
  )
  list(table = tab, plot = p)
}
