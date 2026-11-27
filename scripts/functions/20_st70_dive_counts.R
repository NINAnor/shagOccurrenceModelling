#' Dives per day in harvested and not harvested cells (table + stacked-bar figure)
#'
#' @param prep70 output of st70_prepare() (prep70$dsf has h = harvest status at
#'   the dive time, t = time slice, sexi)
#' @param dmeta dive metadata aligned with prep70$dsf (ID, sex, Year, UTC_datetime)
#' @param n_raw number of raw dive records; n_unique after duplicate removal
#' @param out_dir folder for csv/png files
#' @export
st70_dive_counts <- function(prep70, dmeta, n_raw, n_unique, out_dir) {
  library(dplyr)
  library(ggplot2)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  d <- data.frame(
    ID = as.character(dmeta$ID), sex = as.character(dmeta$sex), Year = as.character(dmeta$Year),
    date = as.Date(dmeta$UTC_datetime, tz = "UTC"), harv = prep70$dsf$h
  )
  daily <- d |>
    group_by(Year, sex, date) |>
    summarise(birds = n_distinct(ID), dives_total = n(), dives_harvested = sum(harv),
              dives_not_harvested = sum(1 - harv), .groups = "drop") |>
    mutate(pct_harvested = round(100 * dives_harvested / dives_total, 1))
  totals <- d |>
    group_by(Year, sex) |>
    summarise(days = n_distinct(date), birds = n_distinct(ID), dives_total = n(),
              dives_harvested = sum(harv), dives_not_harvested = sum(1 - harv), .groups = "drop") |>
    mutate(pct_harvested = round(100 * dives_harvested / dives_total, 1))
  recon <- data.frame(
    Step = c("Raw dive records (2023 + 2024)", "After removing duplicate coordinates within a bird",
             "Dropped: on non-water cells", "Used in the INLA models"),
    Dives = c(n_raw, n_unique, -(n_unique - nrow(d)), nrow(d))
  )
  utils::write.csv(daily, file.path(out_dir, "dives_per_day_harvested.csv"), row.names = FALSE)
  utils::write.csv(totals, file.path(out_dir, "dives_totals_harvested.csv"), row.names = FALSE)
  utils::write.csv(recon, file.path(out_dir, "dives_reconciliation.csv"), row.names = FALSE)

  pl <- daily |>
    select(Year, sex, date, `Harvested cells` = dives_harvested, `Not harvested cells` = dives_not_harvested) |>
    tidyr::pivot_longer(c(`Harvested cells`, `Not harvested cells`), names_to = "type", values_to = "n")
  pl$type <- factor(pl$type, levels = c("Not harvested cells", "Harvested cells"))
  p <- ggplot(pl, aes(date, n, fill = type)) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = c("Harvested cells" = "grey35", "Not harvested cells" = "#E69F00"), name = NULL) +
    facet_grid(sex ~ Year, scales = "free_x") +
    labs(x = NULL, y = "Number of dives per day") +
    theme_bw() + theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "dives_per_day_harvested.png"), p, width = 10, height = 6.5, dpi = 150)
  print(p)
  invisible(list(daily = daily, totals = totals, recon = recon, plot = p))
}
