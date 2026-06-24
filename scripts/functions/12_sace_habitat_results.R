#' Save habitat modelling results in multiple formats (rds, geopkg, GeoTIFF, parquet, plots)
#'
#' This function saves a model results object produced by the modelling pipeline.
#' It supports:
#'  - saving the full results object (RDS/RData),
#'  - extracting and saving fitted model objects,
#'  - storing tabular data as parquet (or CSV fallback if arrow not available),
#'  - writing sf layers as GeoPackage + parquet(WKT),
#'  - writing SpatRaster layers as GeoTIFFs (with GDAL multithreading options) and as xyz parquet,
#'  - saving ggplot objects as PNG + PDF,
#'  - saving basic CV/block tables and a sessionInfo text file.
#'
#' The function also accepts a results list containing multiple named sub-results
#' (for example one entry per harvest variable) and will recurse into each
#' element, creating separate outputs.
#'
#' @param results      A named list or a single results list containing items such as:
#'                     $method, $model, $ensemble_model, $model_data, $model_sf,
#'                     $pred_raster, $pred_raster_by_sex, $pred_raster_sex_period,
#'                     $pred_df, $pred_grid_env, $plot, etc.
#' @param out_dir      Directory where outputs will be written (created if missing)
#' @param prefix       File name prefix for outputs (default "model_fit")
#' @param gdal_threads Number of GDAL threads or "ALL_CPUS" (passed to terra write options)
#'
#' @return character vector of saved file paths (unique)
#' @export
save_habitat_results_pro <- function(results,
                                     out_dir,
                                     prefix = "model_fit",
                                     gdal_threads = "ALL_CPUS") {
  # ---------------------------
  # Dependency checks
  # ---------------------------
  required_pkgs <- c("sf", "terra", "ggplot2")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Install missing packages: ", paste(missing_pkgs, collapse = ", "))
  has_arrow <- requireNamespace("arrow", quietly = TRUE)

  # Ensure output directory exists
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------------------------
  # If results is a named list of sub-results, recurse
  # ---------------------------
  is_single_result <- is.list(results) && all(c("method", "model") %in% names(results))
  if (!is_single_result) {
    all_paths <- list()
    for (nm in names(results)) {
      sub_prefix <- paste0(prefix, "_", nm)
      all_paths[[nm]] <- save_habitat_results_pro(results[[nm]], out_dir, sub_prefix, gdal_threads)
    }
    # flatten and return
    return(unique(unlist(all_paths)))
  }

  # ---------------------------
  # Helper: write tabular object to parquet if possible, otherwise to CSV
  # Returns path written
  # ---------------------------
  safe_write_parquet <- function(df, path_parquet) {
    tryCatch({
      if (has_arrow) {
        # arrow::write_parquet accepts data.frame / tibble
        arrow::write_parquet(df, path_parquet)
        path_parquet
      } else {
        path_csv <- sub("\\.parquet$", ".csv", path_parquet)
        utils::write.csv(df, path_csv, row.names = FALSE)
        path_csv
      }
    }, error = function(e) {
      warning("Failed to write table to ", path_parquet, ": ", conditionMessage(e))
      NA_character_
    })
  }

  # ---------------------------
  # Helper: save an sf object both as GeoPackage and as parquet (WKT geometry column)
  # ---------------------------
  save_sf_parquet_and_gpkg <- function(sf_obj, stem) {
    saved_paths <- character(0)
    if (!inherits(sf_obj, "sf")) return(saved_paths)

    # GeoPackage
    gpkg_path <- file.path(out_dir, paste0(prefix, "_", stem, ".gpkg"))
    tryCatch({
      sf::st_write(sf_obj, gpkg_path, layer = tools::file_path_sans_ext(basename(gpkg_path)),
                   driver = "GPKG", append = FALSE, quiet = TRUE)
      saved_paths <- c(saved_paths, gpkg_path)
    }, error = function(e) {
      warning("Failed to write GeoPackage ", gpkg_path, ": ", conditionMessage(e))
    })

    # Parquet: convert geometry to WKT and drop geometry column
    pq_path <- file.path(out_dir, paste0(prefix, "_", stem, ".parquet"))
    tryCatch({
      dfw <- sf_obj
      dfw$geometry_wkt <- sf::st_as_text(sf::st_geometry(sf_obj))
      dfw <- sf::st_drop_geometry(dfw)
      p <- safe_write_parquet(dfw, pq_path)
      if (!is.na(p)) saved_paths <- c(saved_paths, p)
    }, error = function(e) {
      warning("Failed to write parquet for ", stem, ": ", conditionMessage(e))
    })

    saved_paths
  }

  # ---------------------------
  # Helper: write SpatRaster to GeoTIFF using terra::writeRaster with multithread options
  # Also export an xyz table (parquet or CSV)
  # ---------------------------
  write_raster_and_xyz <- function(spat, stem) {
    saved_paths <- character(0)
    if (!inherits(spat, "SpatRaster")) return(saved_paths)

    tif_path <- file.path(out_dir, paste0(prefix, "_", stem, ".tif"))
    # GDAL write options: compress, predictor, threads; try to use terra::wopt
    wopts <- tryCatch({
      terra::wopt(gdal = c(
        paste0("COMPRESS=LZW"),
        "PREDICTOR=2",
        paste0("NUM_THREADS=", gdal_threads),
        "TILED=YES",
        "BIGTIFF=IF_SAFER"
      ))
    }, error = function(e) {
      NULL
    })

    tryCatch({
      terra::writeRaster(spat, tif_path, overwrite = TRUE, wopt = wopts)
      saved_paths <- c(saved_paths, tif_path)
    }, error = function(e) {
      warning("Failed to write GeoTIFF ", tif_path, ": ", conditionMessage(e))
    })

    # Export xyz table for portability
    xyz_path <- file.path(out_dir, paste0(prefix, "_", stem, "_xyz.parquet"))
    tryCatch({
      df_xyz <- terra::as.data.frame(spat, xy = TRUE, na.rm = FALSE)
      p <- safe_write_parquet(df_xyz, xyz_path)
      if (!is.na(p)) saved_paths <- c(saved_paths, p)
    }, error = function(e) {
      warning("Failed to write raster xyz for ", stem, ": ", conditionMessage(e))
    })

    saved_paths
  }

  # ---------------------------
  # Helper: save ggplot object as PNG + PDF; returns written paths
  # ---------------------------
  save_plot_files <- function(plot_obj, stem, width = 8, height = 6, dpi = 300) {
    if (is.null(plot_obj)) return(character(0))
    pngp <- file.path(out_dir, paste0(prefix, "_", stem, ".png"))
    pdfp <- file.path(out_dir, paste0(prefix, "_", stem, ".pdf"))
    paths <- character(0)
    tryCatch({
      ggplot2::ggsave(pngp, plot = plot_obj, width = width, height = height, dpi = dpi)
      paths <- c(paths, pngp)
    }, error = function(e) {
      warning("Failed to save PNG for ", stem, ": ", conditionMessage(e))
    })
    tryCatch({
      ggplot2::ggsave(pdfp, plot = plot_obj, width = width, height = height)
      paths <- c(paths, pdfp)
    }, error = function(e) {
      warning("Failed to save PDF for ", stem, ": ", conditionMessage(e))
    })
    paths
  }

  # ---------------------------
  # Start saving pieces of the single results object
  # ---------------------------
  saved <- character(0)

  # Save the entire results object (RDS + RData)
  tryCatch({
    rds_path <- file.path(out_dir, paste0(prefix, "_full_results.rds"))
    saveRDS(results, rds_path)
    saved <- c(saved, rds_path)
  }, error = function(e) warning("Failed to save RDS: ", conditionMessage(e)))

  tryCatch({
    rdata_path <- file.path(out_dir, paste0(prefix, "_full_results.RData"))
    save(results, file = rdata_path)
    saved <- c(saved, rdata_path)
  }, error = function(e) warning("Failed to save RData: ", conditionMessage(e)))

  # Save model and model summary (if present)
  if (!is.null(results$model)) {
    tryCatch({
      model_rds <- file.path(out_dir, paste0(prefix, "_model.rds"))
      saveRDS(results$model, model_rds)
      saved <- c(saved, model_rds)

      sum_path <- file.path(out_dir, paste0(prefix, "_model_summary.txt"))
      capture.output({
        cat("Method:", results$method, "\n\n")
        print(summary(results$model))
      }, file = sum_path)
      saved <- c(saved, sum_path)
    }, error = function(e) warning("Failed to save model or summary: ", conditionMessage(e)))
  }

  if (!is.null(results$ensemble_model)) {
    tryCatch({
      ens_rds <- file.path(out_dir, paste0(prefix, "_ensemble_model.rds"))
      saveRDS(results$ensemble_model, ens_rds)
      saved <- c(saved, ens_rds)
    }, error = function(e) warning("Failed to save ensemble model: ", conditionMessage(e)))
  }

  # Save tabular data (model_data, pred_df, block_eval, cv_table)
  if (!is.null(results$model_data)) {
    p <- safe_write_parquet(results$model_data, file.path(out_dir, paste0(prefix, "_model_data.parquet")))
    saved <- c(saved, p)
  }
  if (!is.null(results$pred_df)) {
    p <- safe_write_parquet(results$pred_df, file.path(out_dir, paste0(prefix, "_pred_df.parquet")))
    saved <- c(saved, p)
  }
  if (!is.null(results$block_eval)) {
    p <- safe_write_parquet(results$block_eval, file.path(out_dir, paste0(prefix, "_block_eval.parquet")))
    saved <- c(saved, p)
  }
  if (!is.null(results$cv_table)) {
    tryCatch({
      cv_path <- file.path(out_dir, paste0(prefix, "_cv_table.csv"))
      utils::write.csv(results$cv_table, cv_path, row.names = FALSE)
      saved <- c(saved, cv_path)
    }, error = function(e) warning("Failed to save cv_table: ", conditionMessage(e)))
  }

  # Save sf objects: model_sf, pred_grid_env, harvest_regions (if present)
  if (!is.null(results$model_sf) && inherits(results$model_sf, "sf")) {
    saved <- c(saved, save_sf_parquet_and_gpkg(results$model_sf, "model_sf"))
  }
  if (!is.null(results$pred_grid_env) && inherits(results$pred_grid_env, "sf")) {
    saved <- c(saved, save_sf_parquet_and_gpkg(results$pred_grid_env, "pred_grid_env"))
  }
  if (!is.null(results$harvest_regions) && inherits(results$harvest_regions, "sf")) {
    saved <- c(saved, save_sf_parquet_and_gpkg(results$harvest_regions, "harvest_regions"))
  }

  # Save raster outputs (GeoTIFF + xyz parquet)
  if (!is.null(results$pred_raster) && inherits(results$pred_raster, "SpatRaster")) {
    saved <- c(saved, write_raster_and_xyz(results$pred_raster, "pred_raster"))
  }
  if (!is.null(results$pred_raster_by_sex) && inherits(results$pred_raster_by_sex, "SpatRaster")) {
    saved <- c(saved, write_raster_and_xyz(results$pred_raster_by_sex, "pred_raster_by_sex"))
  }
  if (!is.null(results$pred_raster_sex_period) && inherits(results$pred_raster_sex_period, "SpatRaster")) {
    saved <- c(saved, write_raster_and_xyz(results$pred_raster_sex_period, "pred_raster_sex_period"))
  }

  # Save ggplot objects (plots) if present
  saved <- c(saved,
             save_plot_files(results$plot, "prediction_plot"),
             save_plot_files(results$harvest_conditional_plot, "harvest_conditional_plot"),
             save_plot_files(results$dist_harvest_effect_plot, "dist_harvest_effect_plot"),
             save_plot_files(results$rel_selection_plot, "relative_selection_plot"),
             save_plot_files(results$prediction_plot_by_sex, "prediction_plot_by_sex"),
             save_plot_files(results$prediction_plot_sex_period, "prediction_plot_sex_period"))

  # Save session info for reproducibility
  tryCatch({
    si_path <- file.path(out_dir, paste0(prefix, "_sessionInfo.txt"))
    capture.output(sessionInfo(), file = si_path)
    saved <- c(saved, si_path)
  }, error = function(e) warning("Failed to save sessionInfo: ", conditionMessage(e)))

  # Return unique saved paths, excluding NAs
  unique(saved[!is.na(saved)])
}
