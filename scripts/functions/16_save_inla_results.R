save_inla_results <- function(results,
                              out_dir,
                              prefix = "inla_model") {

  library(sf)
  library(ggplot2)

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  saved <- character()

  #--------------------------------------------------
  # Save complete object
  #--------------------------------------------------

  rds_file <- file.path(
    out_dir,
    paste0(prefix, "_results.rds")
  )

  saveRDS(results, rds_file)

  saved <- c(saved, rds_file)

  #--------------------------------------------------
  # Save model summaries
  #--------------------------------------------------

  if (!is.null(results$fit_spatial)) {

    txt_file <- file.path(
      out_dir,
      paste0(prefix, "_spatial_summary.txt")
    )

    capture.output(
      summary(results$fit_spatial),
      file = txt_file
    )

    saved <- c(saved, txt_file)

    saveRDS(
      results$fit_spatial,
      file.path(
        out_dir,
        paste0(prefix, "_fit_spatial.rds")
      )
    )
  }

  if (!is.null(results$fit_nospatial)) {

    txt_file <- file.path(
      out_dir,
      paste0(prefix, "_nospatial_summary.txt")
    )

    capture.output(
      summary(results$fit_nospatial),
      file = txt_file
    )

    saved <- c(saved, txt_file)

    saveRDS(
      results$fit_nospatial,
      file.path(
        out_dir,
        paste0(prefix, "_fit_nospatial.rds")
      )
    )
  }

  if (!is.null(results$fit_nospatialindividual)) {

    txt_file <- file.path(
      out_dir,
      paste0(prefix, "_nospatialindividual_summary.txt")
    )

    capture.output(
      summary(results$fit_nospatialindividual),
      file = txt_file
    )

    saved <- c(saved, txt_file)

    saveRDS(
      results$fit_nospatialindividual,
      file.path(
        out_dir,
        paste0(prefix, "_fit_nospatialindividual.rds")
      )
    )
  }

  #--------------------------------------------------
  # Save WAIC table
  #--------------------------------------------------

  if (!is.null(results$waic)) {

    waic_file <- file.path(
      out_dir,
      paste0(prefix, "_waic.csv")
    )

    waic_df <- data.frame(
      model = names(results$waic),
      WAIC = as.numeric(results$waic)
    )

    write.csv(
      waic_df,
      waic_file,
      row.names = FALSE
    )

    saved <- c(saved, waic_file)
  }

  #--------------------------------------------------
  # Save mesh
  #--------------------------------------------------

  if (!is.null(results$mesh)) {

    saveRDS(
      results$mesh,
      file.path(
        out_dir,
        paste0(prefix, "_mesh.rds")
      )
    )
  }

  #--------------------------------------------------
  # Save presence data
  #--------------------------------------------------

  if (!is.null(results$presences) &&
      inherits(results$presences, "sf")) {

    pres_file <- file.path(
      out_dir,
      paste0(prefix, "_presences.gpkg")
    )

    sf::st_write(
      results$presences,
      pres_file,
      delete_dsn = TRUE,
      quiet = TRUE
    )

    saved <- c(saved, pres_file)
  }

  #--------------------------------------------------
  # Save prediction surfaces
  #--------------------------------------------------

  save_sf_prediction <- function(x, stem) {

    if (!inherits(x, "sf")) {
      return(NULL)
    }

    file <- file.path(
      out_dir,
      paste0(prefix, "_", stem, ".gpkg")
    )

    sf::st_write(
      x,
      file,
      delete_dsn = TRUE,
      quiet = TRUE
    )

    file
  }

  if (!is.null(results$pred_dh)) {

    saved <- c(
      saved,
      save_sf_prediction(
        results$pred_dh,
        "pred_distHarvest"
      )
    )
  }

  if (!is.null(results$pred_full)) {

    saved <- c(
      saved,
      save_sf_prediction(
        results$pred_full,
        "pred_full"
      )
    )
  }

  #--------------------------------------------------
  # Save plots
  #--------------------------------------------------

  save_plot <- function(plot_obj, stem) {

    if (!inherits(plot_obj, "ggplot")) {
      return(character())
    }

    png_file <- file.path(
      out_dir,
      paste0(prefix, "_", stem, ".png")
    )

    pdf_file <- file.path(
      out_dir,
      paste0(prefix, "_", stem, ".pdf")
    )

    ggsave(
      png_file,
      plot_obj,
      width = 10,
      height = 8,
      dpi = 300
    )

    ggsave(
      pdf_file,
      plot_obj,
      width = 10,
      height = 8
    )

    c(png_file, pdf_file)
  }

  saved <- c(
    saved,
    save_plot(results$plot_dh, "distHarvest_plot"),
    save_plot(results$plot_full, "full_prediction_plot")
  )

  #--------------------------------------------------
  # Save session info
  #--------------------------------------------------

  session_file <- file.path(
    out_dir,
    paste0(prefix, "_sessionInfo.txt")
  )

  capture.output(
    sessionInfo(),
    file = session_file
  )

  saved <- c(saved, session_file)

  unique(saved)
}
