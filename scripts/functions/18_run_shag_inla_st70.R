#' Spatio-temporal LGCP with a 70-slice daily harvest effect
#'
#' Time: unique dive days of 2023 = 1..n23, of 2024 = n23+1..n23+n24.
#' Harvest (1/0, time-varying): 1 if the 50 m cell lies within `buffer_m` of a
#' cut made from 1 September of the previous year up to the time of the dive
#' (integration points: noon of the day). Implemented through a first-cut
#' raster per year: value = 1 if first_cut <= time.
#' Static: kelp presence, depth, depth^2, distance to colony (RW2 spline).
#' Random: sex (iid intercept). Spatial field on the supplied mesh.
#' Integration points: coarse grid of never-harvested cells (harvest = 0 at all
#' times, one copy per year and sex, weight = area x bird-days) plus a random
#' sample of ever-harvested cells replicated per day and sex.
#' @export
st70_first_cut <- function(harvest, year, template, buffer_m = 50) {
  start <- as.POSIXct(paste0(as.integer(year) - 1, "-09-01"), tz = "UTC")
  h <- harvest[!is.na(harvest$GIS_TIME) & harvest$GIS_TIME >= start, ]
  h <- sf::st_transform(h, terra::crs(template))
  g <- sf::st_buffer(sf::st_geometry(h), buffer_m)
  v <- terra::vect(sf::st_sf(t = as.numeric(h$GIS_TIME), geometry = g))
  terra::rasterize(v, template, field = "t", fun = "min", touches = TRUE)
}

#' @export
st70_prepare <- function(
  dsf,
  dmeta,
  grid,
  template,
  harvest,
  buffer_m = 50,
  static_step = 400,
  n_dyn = 600,
  n_bins = 25,
  seed = 5107,
  tmap = NULL,
  brk = NULL
) {
  set.seed(seed)
  dmeta$date <- as.Date(dmeta$UTC_datetime, tz = "UTC")
  dmeta$Year <- as.character(dmeta$Year)
  yrs <- sort(unique(dmeta$Year))
  if (is.null(tmap)) {
    days <- lapply(yrs, function(y) sort(unique(dmeta$date[dmeta$Year == y])))
    names(days) <- yrs
    tmap <- data.frame(
      Year = rep(yrs, lengths(days)),
      date = do.call(c, days)
    )
    tmap$t <- seq_len(nrow(tmap))
  }
  days <- lapply(yrs, function(y) sort(tmap$date[tmap$Year == y]))
  names(days) <- yrs
  dmeta <- dplyr::left_join(dmeta, tmap, by = c("Year", "date"))
  dmeta$sexi <- ifelse(dmeta$sex == "Female", 1L, 2L)
  Tn <- nrow(tmap)

  fcs <- lapply(yrs, function(y) st70_first_cut(harvest, y, template, buffer_m))
  names(fcs) <- yrs

  ## Dives: harvest status at exact time
  h <- numeric(nrow(dsf))
  xy_d <- sf::st_coordinates(dsf)
  for (y in yrs) {
    i <- which(dmeta$Year == y)
    fcv <- terra::extract(fcs[[y]], xy_d[i, , drop = FALSE])[[1]]
    h[i] <- as.numeric(!is.na(fcv) & fcv <= as.numeric(dmeta$UTC_datetime[i]))
  }
  dsf$h <- h
  dsf$t <- dmeta$t
  dsf$sexi <- dmeta$sexi

  ## Effort: birds with dives, per day and sex
  eff <- dmeta |>
    dplyr::group_by(t, sexi, Year) |>
    dplyr::summarise(E = dplyr::n_distinct(ID), .groups = "drop")

  ## Cells
  xy <- terra::xyFromCell(grid, seq_len(terra::ncell(grid)))
  wat <- !is.na(terra::values(grid$water)[, 1])
  res <- terra::res(grid)[1]
  cell_km2 <- (res / 1000)^2
  gv <- terra::values(grid[[c("kelp", "depth", "depth2", "dist")]])
  gv[, "kelp"] <- as.numeric(gv[, "kelp"] > 0)
  okc <- wat & stats::complete.cases(gv)
  stat_idx <- which(
    okc & ((xy[, 1] %% static_step) < res) & ((xy[, 2] %% static_step) < res)
  )

  if (is.null(brk)) {
    brk <- unique(stats::quantile(
      c(dsf$dist, gv[stat_idx, "dist"]),
      seq(0, 1, length.out = n_bins + 1)
    ))
  }
  nbk <- length(brk) - 1
  bin <- function(x) cut(x, brk, include.lowest = TRUE, labels = FALSE)
  dsf$dbin <- bin(dsf$dist)

  rows <- list()
  for (y in yrs) {
    fcv <- terra::values(fcs[[y]])[, 1]
    end_y <- as.numeric(as.POSIXct(max(days[[y]]) + 1, tz = "UTC"))
    ever <- okc & !is.na(fcv) & fcv < end_y
    sidx <- stat_idx[!ever[stat_idx]]
    ey <- eff[eff$Year == y, ]
    t1 <- min(tmap$t[tmap$Year == y])
    for (s in 1:2) {
      Es <- sum(ey$E[ey$sexi == s])
      rows[[length(rows) + 1]] <- data.frame(
        x = xy[sidx, 1],
        y = xy[sidx, 2],
        gv[sidx, , drop = FALSE],
        h = 0,
        t = t1,
        sexi = s,
        weight = static_step^2 / 1e6 * Es
      )
    }
    eidx <- which(ever)
    n_ever <- length(eidx)
    di <- eidx[sample.int(n_ever, min(n_dyn, n_ever))]
    wcell <- n_ever * cell_km2 / length(di)
    for (r in seq_len(nrow(ey))) {
      noon <- as.numeric(as.POSIXct(tmap$date[tmap$t == ey$t[r]], tz = "UTC")) +
        12 * 3600
      rows[[length(rows) + 1]] <- data.frame(
        x = xy[di, 1],
        y = xy[di, 2],
        gv[di, , drop = FALSE],
        h = as.numeric(fcv[di] <= noon),
        t = ey$t[r],
        sexi = ey$sexi[r],
        weight = wcell * ey$E[r]
      )
    }
    message(
      y,
      ": ",
      n_ever,
      " ever-harvested 50 m cells; ",
      length(sidx),
      " static points per sex."
    )
  }
  ips_df <- do.call(rbind, rows)
  ips_df$dbin <- bin(ips_df$dist)
  ips_df$.block <- 1L
  ips <- sf::st_as_sf(ips_df, coords = c("x", "y"), crs = sf::st_crs(dsf))
  list(
    dsf = dsf,
    ips = ips,
    tmap = tmap,
    brk = brk,
    T = Tn,
    nbk = nbk,
    dmeta = dmeta,
    fcs = fcs
  )
}

#' @export
st70_fit <- function(
  prep,
  mesh,
  field_fixed = c(15000, 1),
  n_samples = 300,
  seed = 4172,
  verbose = FALSE
) {
  library(INLA)
  library(inlabru)
  library(dplyr)
  dsf <- prep$dsf
  ips <- prep$ips
  Tn <- prep$T
  nbk <- prep$nbk
  kappa0 <- sqrt(8) / field_fixed[1]
  tau0 <- 1 / (field_fixed[2] * sqrt(4 * pi) * kappa0)
  spde <- INLA::inla.spde2.matern(
    mesh,
    alpha = 2,
    B.tau = cbind(log(tau0), 1, 0),
    B.kappa = cbind(log(kappa0), 0, 1),
    theta.prior.mean = c(0, 0),
    theta.prior.prec = c(1e6, 1e6)
  )
  pcp <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  cmp <- ~ Intercept(1) +
    kelp(kelp, model = "linear", mean.linear = 0, prec.linear = 0.04) +
    depth(depth, model = "linear", mean.linear = 0, prec.linear = 0.04) +
    depth2(depth2, model = "linear", mean.linear = 0, prec.linear = 0.04) +
    distspl(
      dbin,
      model = "rw2",
      values = seq_len(nbk),
      scale.model = TRUE,
      constr = TRUE,
      hyper = pcp
    ) +
    harvest(h, model = "linear", mean.linear = 0, prec.linear = 0.25) +
    harvest_t(
      main = t,
      weights = h,
      model = "rw1",
      values = seq_len(Tn),
      scale.model = TRUE,
      constr = TRUE,
      hyper = pcp
    ) +
    sex_re(main = sexi, model = "iid", values = 1:2, hyper = pcp) +
    field(geometry, model = spde)
  fit <- inlabru::bru(
    cmp,
    inlabru::bru_obs(geometry ~ ., family = "cp", data = dsf, ips = ips),
    options = inlabru::bru_options(
      verbose = verbose,
      control.fixed = list(
        mean.intercept = log(nrow(dsf) / sum(ips$weight)),
        prec.intercept = 0.01
      )
    )
  )
  smp <- inlabru::generate(
    fit,
    newdata = data.frame(h = 1, t = seq_len(Tn)),
    formula = ~ harvest + harvest_t,
    n.samples = n_samples,
    seed = seed
  )
  q <- t(apply(smp, 1, stats::quantile, probs = c(0.025, 0.5, 0.975)))
  expo <- sf::st_drop_geometry(dsf) |>
    dplyr::group_by(t) |>
    dplyr::summarise(
      dives = dplyr::n(),
      dives_in_harvest = sum(h),
      .groups = "drop"
    )
  curves <- data.frame(
    t = seq_len(Tn),
    mean = rowMeans(smp),
    q025 = q[, 1],
    q975 = q[, 3],
    p_neg = rowMeans(smp < 0)
  ) |>
    dplyr::left_join(prep$tmap, by = "t") |>
    dplyr::left_join(expo, by = "t")
  curves$conclusion <- dplyr::case_when(
    curves$q025 > 0 ~ "prefers harvested areas",
    curves$q975 < 0 ~ "avoids harvested areas",
    TRUE ~ "no clear effect"
  )
  list(fit = fit, curves = curves)
}

#' @export
st70_plot <- function(curves, min_harvest_dives = 20, ylim = c(-6, 5)) {
  library(ggplot2)
  curves$low <- curves$dives_in_harvest < min_harvest_dives
  ggplot(curves, aes(date, mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.2) +
    geom_line(colour = "grey30") +
    geom_point(aes(colour = conclusion, shape = low), size = 2.5) +
    scale_colour_manual(
      values = c(
        "prefers harvested areas" = "#0072B2",
        "avoids harvested areas" = "#D55E00",
        "no clear effect" = "grey60"
      )
    ) +
    scale_shape_manual(
      values = c("FALSE" = 16, "TRUE" = 1),
      name = paste0("< ", min_harvest_dives, " dives\nin harvest")
    ) +
    facet_wrap(~Year, ncol = 1, scales = "free_x") +
    coord_cartesian(ylim = ylim) +
    labs(
      x = NULL,
      y = "Harvest effect (log scale, 95% CrI)\nharvested vs not harvested 50 m cell",
      colour = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}

#' Sex-specific version: separate daily harvest trend per sex (iid groups) plus a
#' sex shift of the harvest level and a sex intercept.
#' @export
st70_fit_sex <- function(
  prep,
  mesh,
  field_fixed = c(15000, 1),
  n_samples = 300,
  seed = 4172,
  verbose = FALSE
) {
  library(INLA)
  library(inlabru)
  library(dplyr)
  dsf <- prep$dsf
  ips <- prep$ips
  Tn <- prep$T
  nbk <- prep$nbk
  kappa0 <- sqrt(8) / field_fixed[1]
  tau0 <- 1 / (field_fixed[2] * sqrt(4 * pi) * kappa0)
  spde <- INLA::inla.spde2.matern(
    mesh,
    alpha = 2,
    B.tau = cbind(log(tau0), 1, 0),
    B.kappa = cbind(log(kappa0), 0, 1),
    theta.prior.mean = c(0, 0),
    theta.prior.prec = c(1e6, 1e6)
  )
  pcp <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
  cmp <- ~ Intercept(1) +
    kelp(kelp, model = "linear", mean.linear = 0, prec.linear = 0.04) +
    depth(depth, model = "linear", mean.linear = 0, prec.linear = 0.04) +
    depth2(depth2, model = "linear", mean.linear = 0, prec.linear = 0.04) +
    distspl(
      dbin,
      model = "rw2",
      values = seq_len(nbk),
      scale.model = TRUE,
      constr = TRUE,
      hyper = pcp
    ) +
    harvest(h, model = "linear", mean.linear = 0, prec.linear = 0.25) +
    sex_int(main = sexi, model = "iid", values = 1:2, hyper = pcp) +
    sex_harv(
      main = sexi,
      weights = h,
      model = "iid",
      values = 1:2,
      hyper = pcp
    ) +
    harvest_t(
      main = t,
      group = sexi,
      weights = h,
      model = "rw1",
      values = seq_len(Tn),
      ngroup = 2,
      control.group = list(model = "iid"),
      scale.model = TRUE,
      constr = TRUE,
      hyper = pcp
    ) +
    field(geometry, model = spde)
  fit <- inlabru::bru(
    cmp,
    inlabru::bru_obs(geometry ~ ., family = "cp", data = dsf, ips = ips),
    options = inlabru::bru_options(
      verbose = verbose,
      control.fixed = list(
        mean.intercept = log(nrow(dsf) / sum(ips$weight)),
        prec.intercept = 0.01
      )
    )
  )
  nd <- data.frame(h = 1, t = rep(seq_len(Tn), 2), sexi = rep(1:2, each = Tn))
  smp <- inlabru::generate(
    fit,
    newdata = nd,
    formula = ~ harvest + sex_harv + harvest_t,
    n.samples = n_samples,
    seed = seed
  )
  q <- t(apply(smp, 1, stats::quantile, probs = c(0.025, 0.5, 0.975)))
  expo <- sf::st_drop_geometry(dsf) |>
    dplyr::group_by(t, sexi) |>
    dplyr::summarise(
      dives = dplyr::n(),
      dives_in_harvest = sum(h),
      .groups = "drop"
    )
  curves <- data.frame(
    t = nd$t,
    sexi = nd$sexi,
    sex = c("Female", "Male")[nd$sexi],
    mean = rowMeans(smp),
    q025 = q[, 1],
    q975 = q[, 3],
    p_neg = rowMeans(smp < 0)
  ) |>
    dplyr::left_join(prep$tmap, by = "t") |>
    dplyr::left_join(expo, by = c("t", "sexi"))
  curves$dives_in_harvest[is.na(curves$dives_in_harvest)] <- 0
  curves$conclusion <- dplyr::case_when(
    curves$q025 > 0 ~ "prefers harvested areas",
    curves$q975 < 0 ~ "avoids harvested areas",
    TRUE ~ "no clear effect"
  )
  list(fit = fit, curves = curves)
}

#' @export
st70_plot_sex <- function(curves, min_harvest_dives = 20, ylim = c(-6, 5)) {
  library(ggplot2)
  curves$low <- curves$dives_in_harvest < min_harvest_dives
  ggplot(curves, aes(date, mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.2) +
    geom_line(colour = "grey30") +
    geom_point(aes(colour = conclusion, shape = low), size = 2.5) +
    scale_colour_manual(
      values = c(
        "prefers harvested areas" = "#0072B2",
        "avoids harvested areas" = "#D55E00",
        "no clear effect" = "grey60"
      )
    ) +
    scale_shape_manual(
      values = c("FALSE" = 16, "TRUE" = 1),
      name = paste0("< ", min_harvest_dives, " dives\nin harvest")
    ) +
    facet_grid(sex ~ Year, scales = "free_x") +
    coord_cartesian(ylim = ylim) +
    labs(
      x = NULL,
      y = "Harvest effect (log scale, 95% CrI)\nharvested vs not harvested 50 m cell",
      colour = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}
