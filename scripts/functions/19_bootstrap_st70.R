#' Bird-level bootstrap of the sex-specific 68-slice model (st70_fit_sex)
#'
#' Birds are resampled with replacement within sex x year; a bird drawn twice
#' contributes its dives twice (unique copy ID, so effort counts both).
#' Time map and distance bins are fixed to the full-data ones. Integration
#' points (same seed) and effort weights are rebuilt from the resample.
#' @export
st70_bootstrap <- function(
  dsf0,
  dmeta,
  grid,
  template,
  harvest,
  mesh,
  full_prep,
  ids,
  save_dir,
  seed_base = 7001,
  n_samples = 100
) {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  dmeta$ID <- as.character(dmeta$ID)
  rows_by_id <- split(seq_len(nrow(dmeta)), dmeta$ID)
  birds <- unique(dmeta[, c("ID", "sex", "Year")])
  grp <- split(birds$ID, paste(birds$sex, birds$Year))
  brk <- full_prep$brk
  for (b in ids) {
    f <- file.path(save_dir, sprintf("boot_%03d.csv", b))
    if (file.exists(f)) {
      next
    }
    set.seed(seed_base + b)
    samp <- unlist(lapply(grp, function(g) {
      sample(g, length(g), replace = TRUE)
    }))
    rows <- unlist(rows_by_id[samp], use.names = FALSE)
    dm <- dmeta[rows, ]
    dm$ID <- rep(paste0(samp, "_b", seq_along(samp)), lengths(rows_by_id[samp]))
    pb <- suppressMessages(st70_prepare(
      dsf0[rows, ],
      dm,
      grid,
      template,
      harvest,
      tmap = full_prep$tmap,
      brk = brk
    ))
    cv <- tryCatch(
      suppressMessages(
        st70_fit_sex(pb, mesh, n_samples = n_samples, seed = b)$curves
      ),
      error = function(e) {
        message("Resample ", b, " failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(cv)) {
      cv$boot <- b
      utils::write.csv(
        cv[, c("boot", "t", "sex", "Year", "date", "mean", "q025", "q975")],
        f,
        row.names = FALSE
      )
    }
  }
  invisible(NULL)
}

#' Combine resamples; classify days (exposure gate + bootstrap agreement)
#' @export
st70_bootstrap_summary <- function(
  full_curves,
  dsf,
  dmeta_full_t,
  save_dir,
  threshold = 0.6,
  min_harvest_dives = 20,
  min_birds = 3
) {
  library(dplyr)
  library(ggplot2)
  bt <- dplyr::bind_rows(lapply(
    list.files(save_dir, pattern = "^boot_.*csv$", full.names = TRUE),
    utils::read.csv
  ))
  bt$Year <- as.character(bt$Year)
  B <- dplyr::n_distinct(bt$boot)
  bs <- bt |>
    dplyr::group_by(t, sex) |>
    dplyr::summarise(
      n_boot = dplyr::n(),
      prop_prefers = mean(q025 > 0),
      prop_avoids = mean(q975 < 0),
      boot_lo = stats::quantile(mean, 0.025),
      boot_hi = stats::quantile(mean, 0.975),
      .groups = "drop"
    )
  expo <- data.frame(
    t = dsf$t,
    sex = ifelse(dsf$sexi == 1, "Female", "Male"),
    h = dsf$h,
    bird = dmeta_full_t$ID
  ) |>
    dplyr::group_by(t, sex) |>
    dplyr::summarise(
      dives = dplyr::n(),
      harvest_dives = sum(h),
      harvest_birds = dplyr::n_distinct(bird[h == 1]),
      .groups = "drop"
    )
  res <- full_curves |>
    dplyr::select(t, sex, Year, date, mean, q025, q975, conclusion) |>
    dplyr::left_join(bs, by = c("t", "sex")) |>
    dplyr::left_join(expo, by = c("t", "sex")) |>
    dplyr::mutate(
      harvest_dives = ifelse(is.na(harvest_dives), 0, harvest_dives),
      harvest_birds = ifelse(is.na(harvest_birds), 0, harvest_birds),
      informative = harvest_dives >= min_harvest_dives &
        harvest_birds >= min_birds,
      status = dplyr::case_when(
        !informative ~ "not informative (low exposure)",
        conclusion == "prefers harvested areas" &
          prop_prefers >= threshold ~ "prefers (robust)",
        conclusion == "avoids harvested areas" &
          prop_avoids >= threshold ~ "avoids (robust)",
        TRUE ~ "no robust effect"
      )
    )
  zones <- res |>
    dplyr::arrange(Year, sex, date) |>
    dplyr::group_by(Year, sex) |>
    dplyr::mutate(
      run = cumsum(status != dplyr::lag(status, default = dplyr::first(status)))
    ) |>
    dplyr::group_by(Year, sex, run, status) |>
    dplyr::summarise(
      from = format(min(date), "%d %b"),
      to = format(max(date), "%d %b"),
      days = dplyr::n(),
      mean_beta = round(mean(mean), 2),
      factor_exp_beta = round(mean(exp(mean)), 2),
      harvest_dives = sum(harvest_dives),
      .groups = "drop"
    ) |>
    dplyr::select(-run) |>
    dplyr::filter(status %in% c("prefers (robust)", "avoids (robust)"))
  plt <- ggplot(res, aes(date, mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_ribbon(aes(ymin = boot_lo, ymax = boot_hi), alpha = 0.2) +
    geom_line(colour = "grey30") +
    geom_point(aes(colour = status), size = 2.5) +
    scale_colour_manual(
      values = c(
        "prefers (robust)" = "#0072B2",
        "avoids (robust)" = "#D55E00",
        "no robust effect" = "grey60",
        "not informative (low exposure)" = "white"
      )
    ) +
    geom_point(
      data = dplyr::filter(res, status == "not informative (low exposure)"),
      shape = 1,
      size = 2.5,
      colour = "grey50"
    ) +
    facet_grid(sex ~ Year, scales = "free_x") +
    coord_cartesian(ylim = c(-6, 5)) +
    labs(
      x = NULL,
      y = "Harvest effect (log scale)\nharvested vs not harvested 50 m cell",
      colour = NULL,
      caption = sprintf(
        "Line: full-data estimate; band: bootstrap 95%% interval (%d resamples). Coloured: full-fit interval excludes 0, >= %d%% of resamples agree,\nand >= %d dives in harvested cells from >= %d birds. Open circles: low exposure.",
        B,
        round(100 * threshold),
        min_harvest_dives,
        min_birds
      )
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  list(results = res, zones = zones, plot = plt, B = B)
}
