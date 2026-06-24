#' Run exploratory data analysis (EDA) on used-available data by year
#'
#' For each year in `years`, subset the input dataset, coerce/scale covariates,
#' produce simple frequency and cross-tab summaries, and generate diagnostic
#' plots (histogram/boxplot for DistNest, bar charts for harvest status and
#' logger model). Optionally produces a pairwise plot matrix using GGally if
#' available.
#'
#' Inputs:
#'  - dmodel       : sf object containing used/available points and covariates
#'                   (must include Year and datetime/other covariates referenced)
#'  - years        : character vector of years to process (default c("2023","2024"))
#'  - response_col : name of response column in dmodel (default "used")
#'  - cov_cols     : covariates to consider for summaries and plotting
#'  - factor_vars  : subset of cov_cols that should be coerced to factor
#'
#' Returns:
#'  - a named list (one entry per year). Each entry is a list with:
#'      * data_for_model    : sf object subset used for EDA (with response renamed to "response")
#'      * missing_values    : named vector of NA counts per column
#'      * frequency_tables  : named list of frequency tables for categorical vars
#'      * cross_tables      : named list of contingency tables: response by factor
#'      * plots             : list of ggplot (and possibly GGally) objects
#'
#' Notes:
#'  - Function attempts to be robust when some covariates are missing:
#'    it will only operate on the subset of covariates present in dmodel.
#'  - The DistNest variable is scaled (centered & scaled) for plotting.
#'  - Pairwise plots require package GGally; if not installed the pairs plot
#'    is skipped with a warning.
#'
#' @export
run_eda_by_year <- function(dmodel,
                            years = c("2023", "2024"),
                            response_col = "used",
                            cov_cols = c(
                              "sex",
                              "Logger_model",
                              "DistNest",
                              "harvested200m",
                              "harvested50m",
                              "harvesttimeSK24"
                            ),
                            factor_vars = c(
                              "sex",
                              "Logger_model",
                              "harvested200m",
                              "harvested50m",
                              "harvesttimeSK24"
                            )) {

  # -----------------------------
  # Dependency checks
  # -----------------------------
  req_pkgs <- c("dplyr", "sf", "ggplot2", "tidyr")
  missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) stop("Install missing packages: ", paste(missing_pkgs, collapse = ", "))

  has_ggally <- requireNamespace("GGally", quietly = TRUE)

  `%>%` <- dplyr::`%>%`

  eda_results <- list()

  for (yr in years) {
    message("Running EDA for year: ", yr)

    # Keep only covariates that actually exist in the data
    cov_present <- intersect(cov_cols, names(dmodel))
    factor_present <- intersect(factor_vars, cov_present)

    # Subset to year and selected columns (preserve geometry)
    data_for_model <- dmodel %>%
      dplyr::filter(Year == yr) %>%
      dplyr::select(dplyr::all_of(response_col), dplyr::all_of(cov_present), geometry)

    # Standardize response column name to "response"
    if (!response_col %in% names(data_for_model)) {
      stop("response_col '", response_col, "' not found in data for year ", yr)
    }
    names(data_for_model)[names(data_for_model) == response_col] <- "response"

    # Scale DistNest if present (numeric scaling for comparability)
    if ("DistNest" %in% names(data_for_model)) {
      # convert to numeric if sf stored it oddly
      data_for_model <- data_for_model %>%
        dplyr::mutate(DistNest = as.numeric(DistNest)) %>%
        dplyr::mutate(DistNest = as.numeric(scale(DistNest)))
    }

    # Coerce response and factor variables to factor
    data_for_model <- data_for_model %>%
      dplyr::mutate(response = as.factor(response)) %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(factor_present), as.factor))

    # Missing value counts (drop geometry before counting)
    missing_values <- colSums(is.na(sf::st_drop_geometry(data_for_model)))

    # Frequency tables for categorical variables (response + factor_present)
    freq_vars <- c("response", factor_present)
    frequency_tables <- lapply(freq_vars, function(v) {
      if (v %in% names(data_for_model)) table(sf::st_drop_geometry(data_for_model)[[v]])
      else NULL
    })
    names(frequency_tables) <- freq_vars

    # Cross tables: response by each categorical factor
    cross_tables <- lapply(factor_present, function(v) {
      if (v %in% names(data_for_model)) table(sf::st_drop_geometry(data_for_model)$response,
                                              sf::st_drop_geometry(data_for_model)[[v]])
      else NULL
    })
    names(cross_tables) <- paste0("response_by_", factor_present)

    # Plots
    plots <- list()

    # DistNest plots: histogram and boxplot vs response
    if ("DistNest" %in% names(data_for_model)) {
      plots$distnest_hist <- ggplot2::ggplot(data_for_model, ggplot2::aes(x = DistNest)) +
        ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "black") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Distribution of scaled Distance to nest —", yr),
                      x = "Scaled DistNest", y = "Count")

      plots$distnest_box <- ggplot2::ggplot(data_for_model, ggplot2::aes(x = response, y = DistNest, fill = response)) +
        ggplot2::geom_boxplot(alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("DistNest by used/available —", yr),
                      x = "Used / available", y = "Scaled DistNest", fill = "Used")
    }

    # Harvest status plots (proportional bars)
    if ("harvested50m" %in% names(data_for_model)) {
      plots$harvest50 <- ggplot2::ggplot(data_for_model, ggplot2::aes(x = harvested50m, fill = response)) +
        ggplot2::geom_bar(position = "fill") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Harvest status within 50 m —", yr),
                      x = "Harvested 50m", y = "Proportion", fill = "Used")
    }

    if ("harvested200m" %in% names(data_for_model)) {
      plots$harvest200 <- ggplot2::ggplot(data_for_model, ggplot2::aes(x = harvested200m, fill = response)) +
        ggplot2::geom_bar(position = "fill") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Harvest status within 200 m —", yr),
                      x = "Harvested 200m", y = "Proportion", fill = "Used")
    }

    # Logger_model frequency bar
    if ("Logger_model" %in% names(data_for_model)) {
      plots$logger <- ggplot2::ggplot(data_for_model, ggplot2::aes(x = Logger_model, fill = response)) +
        ggplot2::geom_bar(position = "dodge") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Logger model frequencies —", yr),
                      x = "Logger model", y = "Count", fill = "Used")
    }

    # Pairwise plot matrix (if GGally is available)
    # Use only non-geometric, non-NA columns for ggpairs
    data_for_pairs <- sf::st_drop_geometry(data_for_model) %>% tidyr::drop_na()

    if (nrow(data_for_pairs) == 0) {
      warning("No complete rows available for pairwise plots for year ", yr)
      plots$pairs <- NULL
    } else if (!has_ggally) {
      warning("Package 'GGally' not available: skipping pairwise plot for year ", yr)
      plots$pairs <- NULL
    } else {
      # Try to run ggpairs; limit to reasonable number of columns to avoid huge plots
      # Select numeric + factor columns (response included)
      max_cols <- 8
      sel_cols <- names(data_for_pairs)
      # exclude ID-like columns if present
      sel_cols <- setdiff(sel_cols, c("geometry"))
      if (length(sel_cols) > max_cols) sel_cols <- sel_cols[1:max_cols]
      suppressWarnings({
        plots$pairs <- GGally::ggpairs(data_for_pairs[sel_cols], mapping = ggplot2::aes(color = response, alpha = 0.5))
      })
    }

    # Package up results for the year
    eda_results[[as.character(yr)]] <- list(
      data_for_model = data_for_model,
      missing_values = missing_values,
      frequency_tables = frequency_tables,
      cross_tables = cross_tables,
      plots = plots
    )
  }

  eda_results
}
