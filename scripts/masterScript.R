# ============================================================
# MASTER SCRIPT: Fit habitat selection model
# Cleaned, annotated, and formatted for publication
# ============================================================
#
# Purpose:
#   - Prepare data, fit habitat selection models (GAM or ensemble),
#     predict across a study region, and save results and logs.
#
# Usage:
#   - On an interactive/personal machine, set run_env in your .Rprofile
#     (or adjust the flags below).
#   - On a server, pass the parameter-row index (tpr) as the first
#     command-line argument, e.g.:
#       Rscript fit_model.R 3
#
# Requirements / assumptions:
#   - The project defines these variables in .Renviron / .Rprofile:
#       project_dir, temp_dir, data_dir, logs_dir, result_dir, run_env
#   - Project functions are placed under scripts/functions/ and will be sourced.
#   - The pipeline scripts are placed under scripts/pipeline/.
#   - Key functions used here (examples): fit_predict_habitat_model(),
#     save_habitat_results_pro(). These should be defined in sourced files.
#
# Note:
#   - This file is intended for readability and sharing. It keeps the original
#     logic but clarifies variable names, adds documentation, and improves logging.
# ============================================================


# ============================================================
# 0. Global script settings (edit as needed)
# ============================================================
version_number <- 0       # increment for new published runs
checks <- FALSE           # run extra exploratory checks if TRUE
seed <- 1994              # model RNG seed for reproducibility

# default parameter-row index (overridden by command-line argument on server)
param_row <- 1


# ============================================================
# 1. Load environment, command-line args, and basic paths
# ============================================================
# Load project-local environment variables (if present)
if (file.exists(".Renviron")) readRenviron(".Renviron")
if (file.exists(".Rprofile")) source(".Rprofile")

# If run_env is not defined, warn the user and assume "personal"
if (!exists("run_env")) {
  run_env <- "personal"
  warning("run_env not found in .Rprofile/.Renviron: assuming 'personal'.")
}

# On non-personal environments (e.g., a server), read command-line args
if (!run_env %in% "personal") {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    param_row <- as.numeric(args[1])
  }
}

# Print core paths (assumes these are defined in .Rprofile)
message("Working directory: ", ifelse(exists("project_dir"), project_dir, "project_dir not set"))
message("Temporary directory: ", ifelse(exists("temp_dir"), temp_dir, "temp_dir not set"))
message("Data directory: ", ifelse(exists("data_dir"), data_dir, "data_dir not set"))
message("Parameter row (param_row): ", param_row)


# ============================================================
# 2. Output versioning and directories
# ============================================================
version_label <- paste0("version_0.1.", version_number)

path_for_processed_data <- file.path(data_dir, version_label)
dir.create(path_for_processed_data, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 3. Load packages and project functions
#    - All function scripts in scripts/functions/ are sourced here.
# ============================================================
function_files <- list.files(file.path("scripts", "functions"),
                             pattern = "\\.R$", full.names = TRUE)
invisible(lapply(function_files, source))

# Terra/ GDAL threading helper ------------------------------------------------
# Choose thread counts appropriate for your environment.
# Example uses 4 terra threads and 1 OpenMP/GDAL thread each.
set_terra_threads(n_terra = 4, omp_threads = 1, gdal_threads = 1)

# Helper to consistently set terra write options with GDAL multithreading
terra_wopt_mt <- function(comp = "LZW", threads = NULL) {
  if (is.null(threads)) threads <- Sys.getenv("GDAL_NUM_THREADS", unset = "ALL_CPUS")
  terra::wopt(gdal = c(
    paste0("COMPRESS=", comp),
    "PREDICTOR=2",
    paste0("NUM_THREADS=", threads),
    "TILED=YES",
    "BIGTIFF=IF_SAFER"
  ))
}


# ============================================================
# 4. Load and prepare input data (pipeline scripts)
# ============================================================
# These scripts should prepare the following objects (examples):
#  - dmodel: the model dataset (points with covariates and response)
#  - pred_grid / pred_coords: prediction grid for spatial predictions
#  - study_region: polygon for the study area
#  - env_stack_10m, raster_results_23/24: environmental rasters/stacks
# Adjust file paths according to your project layout.

rerunDataFormat <- FALSE
if (rerunDataFormat) unlink(file.path(data_dir, "dmodel.RData"))

# Sampling and buffering parameters used by the preparation scripts:
n_per_point <- 5        # number of random background points per used point
ID_buffer_m <- 10000    # individual ID buffer (meters)
harvest_buffer <- 50    # harvest buffer (meters)

# Source the pipeline preparation scripts (they must exist under scripts/pipeline/)
source(file.path("scripts", "pipeline", "01_dataPreparation", "01_prepData.R"))
source(file.path("scripts", "pipeline", "01_dataPreparation", "02_prepPredictionRaster.R"))
source(file.path("scripts", "pipeline", "01_dataPreparation", "03_predictionCoords.R"))

if (checks) {
  # Optional exploratory analysis (only run when checks == TRUE)
  source(file.path("pipeline", version_label, "02_modelFit", "01_exploration_analysis.R"))
}


# ============================================================
# 5. Parameter grid and strata summaries
# ============================================================
# Create a grid of parameter combinations to run. The script selects one
# row from this grid based on `param_row`. Users can extend the grid with
# more factors as required.

pars <- expand.grid(
  Year = levels(dmodel$Year),
  method = c("gam"),
  stringsAsFactors = FALSE
)

# Which harvest variables we will evaluate (fit_predict_habitat_model
# will loop internally and return results per variable)
harvest_vars_to_use <- c("harvested50m")

# Compute simple odds ratio summaries by strata to help interpret results
or_results <- dmodel %>%
  group_by(sex, Year, harvesttimeSK24) %>%
  summarise(
    used0_h0 = sum(used == 0 & harvested50m == 0),
    used1_h0 = sum(used == 1 & harvested50m == 0),
    used0_h1 = sum(used == 0 & harvested50m == 1),
    used1_h1 = sum(used == 1 & harvested50m == 1),
    used0_h0_200 = sum(used == 0 & harvested200m == 0),
    used1_h0_200 = sum(used == 1 & harvested200m == 0),
    used0_h1_200 = sum(used == 0 & harvested200m == 1),
    used1_h1_200 = sum(used == 1 & harvested200m == 1),
    ids = length(unique(ID)),
    .groups = "drop"
  ) %>%
  mutate(
    OR_50m = (used1_h1 * used0_h0) / (used0_h1 * used1_h0),
    logOR_50m = log(OR_50m),
    OR_200m = (used1_h1_200 * used0_h0_200) / (used0_h1_200 * used1_h0_200),
    logOR_200m = log(OR_200m)
  ) %>%
  sf::st_drop_geometry()

# Validate the requested parameter-row index
if (param_row < 1 || param_row > nrow(pars)) {
  stop("param_row is outside the parameter grid range.")
}

# Select the parameter combination to use in this run
pars_to_use <- pars[param_row, , drop = FALSE]

# Convert to a list for passing into the model-fitting function
pars_to_use_list <- list(
  Year = pars_to_use$Year
)

message("Selected parameters:")
print(pars_to_use)


# ============================================================
# 6. Prepare logging (create log path and logfile)
# ============================================================
path_for_logs <- file.path(logs_dir,
                           version_label,
                           pars_to_use$method,
                           pars_to_use$Year)
dir.create(path_for_logs, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(path_for_logs,
                     paste0("model_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))


# ============================================================
# 7. Fit model with sinked logfile (and robust error handling)
# ============================================================
log_con <- file(logfile, open = "wt")
# Sink stdout and messages to the logfile; split = TRUE keeps console output too
sink(log_con, split = TRUE)
sink(log_con, type = "message")

# Use tryCatch to ensure sinks are reset and the log file closed even on error
tryCatch({

  cat("Started:", as.character(Sys.time()), "\n\n")
  cat("Selected parameters:\n")
  print(pars_to_use)

  cat("\nOdds ratios (50m and 200m harvest buffers):\n")
  rt <- or_results %>%
    dplyr::filter(Year == pars_to_use$Year)
  print(rt)

  # ----------------------------
  # Fit model and make predictions
  # ----------------------------
  # fit_predict_habitat_model should:
  #   - take the prepared data and parameter list
  #   - fit the model(s) (GAM or ensemble) for each harvest variable
  #   - generate predictions across pred_grid or pred_coords
  #   - return a list containing at least:
  #       $method   -- "gam" or "ensemble"
  #       $model    -- the fitted model object (for method == "gam")
  #       $predictions  -- raster or data.frame of predicted values
  #
  # Important arguments used here:
  #   response_col = "random"  -- the outcome column to use
  #   nBasis = 50              -- number of basis functions for GAM splines
  set.seed(seed)
  model_results <- fit_predict_habitat_model(
    dmodel = dmodel,
    pred_grid = pred_grid,
    study_region = study_region,
    env_stack_10m = env_stack_10m,
    raster_results_23 = raster_results_23,
    raster_results_24 = raster_results_24,
    pars_to_use = pars_to_use_list,
    response_col = "random",
    nthreads = 1,
    nBasis = 50
  )

  cat("\n==============================================\n")
  cat("Model summary for year:", pars_to_use$Year, "\n")
  cat("==============================================\n")
  if (!is.null(model_results$method) && model_results$method == "gam") {
    # Print summary for GAM fits (if applicable)
    print(summary(model_results$model))
  } else {
    cat("Model method:", model_results$method, "\n")
  }

  # ----------------------------
  # Save model outputs
  # ----------------------------
  output_dir <- file.path(result_dir, version_label, pars_to_use$method, pars_to_use$Year)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  saved_files <- save_habitat_results_pro(
    results = model_results,
    out_dir = output_dir,
    prefix = "model_fit",
    gdal_threads = Sys.getenv("GDAL_NUM_THREADS", "ALL_CPUS")
  )

  cat("\nSaved files:\n")
  print(saved_files)

  # Record completion time and minimal reproducibility info
  cat("\nFinished:", as.character(Sys.time()), "\n\n")
  cat("R session info:\n")
  print(sessionInfo())

}, error = function(e) {
  # Log the error to the logfile and re-throw it to stop execution
  cat("\nERROR:\n")
  cat(conditionMessage(e), "\n")
  stop(e)
}, finally = {
  # Always reset sinks and close the log connection
  sink(type = "message")
  sink()
  close(log_con)
})

# End of script
# ============================================================
