# ============================================================
# MASTER SCRIPT: Fit habitat selection model
# ============================================================
# This script:
#   1. Sets paths depending on run environment
#   2. Loads packages, functions, and prepared data
#   3. Selects one parameter combination using command-line input
#   4. Creates a logfile
#   5. Fits either GAM or MaxEnt
#   6. Saves model outputs
# ============================================================


# ============================================================
# 1. Run settings
# ============================================================

run_on_personal_computer <- FALSE
run_on_nina_server       <- TRUE

version_number <- 4
checks <- FALSE

# Parameter-row index.
# On the server this is read from command-line arguments.
tpr <- 1


# set model seed
seed <- 1994

# ============================================================
# 2. Set working directory and temporary directory
# ============================================================

if (run_on_personal_computer) {
  
  project_dir <- "R:/Prosjekter/Toppskarv/Analyses_Nina_KelpHarvest/kwaku"
  temp_dir <- tempdir()
  
} else if (run_on_nina_server) {
  
  project_dir <- "~/Mounts/R/Prosjekter/Toppskarv/Analyses_Nina_KelpHarvest/kwaku"
  temp_dir <- "/data/scratch/kwaku"
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) > 0) {
    tpr <- as.numeric(args[1])
  }
  
} else {
  
  stop("No run environment selected.")
}

setwd(project_dir)
#Sys.setenv(TMPDIR = temp_dir)

cat("Working directory:", getwd(), "\n")
cat("Temporary directory:", temp_dir, "\n")
cat("Parameter row:", tpr, "\n")


# ============================================================
# 3. Define output paths
# ============================================================

version_label <- paste0("version_0", version_number)

path_for_processed_data <- file.path(
  "dataset",
  "processed_data",
  version_label
)

dir.create(
  path_for_processed_data,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. Load packages and project functions
# ============================================================

# If current folders do not exist for the current version, duplicate from previous version
# The intention is to edit these new files whist not altering the previous versions
if (run_on_nina_server){
version_label <- paste0("version_0", version_number)

previous_version <- version_number - 1
previous_label <- paste0("version_0", previous_version)

# folders to version
base_dirs <- c("functions", "pipeline")

for (base_dir in base_dirs) {
  
  current_dir <- file.path(base_dir, version_label)
  previous_dir <- file.path(base_dir, previous_label)
  
  # initialise current version from previous version
  if (!dir.exists(current_dir)) {
    
    if (!dir.exists(previous_dir)) {
      stop("Previous version folder does not exist: ", previous_dir)
    }
    
    # create current version folder
    dir.create(current_dir, recursive = TRUE)
    
    # copy all files/folders
    file.copy(
      from = list.files(
        previous_dir,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
      ),
      to = current_dir,
      recursive = TRUE
    )
    
    message("Initialised: ", current_dir)
  }
}

} else {
  # ============================================================
  # 4a. Validate version folders
  # ============================================================
  # Do NOT initialise version folders inside parallel jobs.
  # Parallel jobs should only read existing functions/pipeline code.
  
  required_dirs <- c(
    file.path("functions", version_label),
    file.path("pipeline", version_label)
  )
  
  missing_dirs <- required_dirs[!dir.exists(required_dirs)]
  
  if (length(missing_dirs) > 0) {
    stop(
      "Required version folders are missing:\n",
      paste(missing_dirs, collapse = "\n"),
      "\n\nCreate/copy these folders once before launching parallel jobs."
    )
  }
}

# source package loader
source(
  file.path(
    "functions",
    version_label,
    "00_loadPackages.R"
  )
)


function_files <- list.files(
  file.path("functions", version_label),
  pattern = "\\.R$",
  full.names = TRUE
)

invisible(
  lapply(function_files, source)
)


# ============================================================
# 5. Load and prepare data
# ============================================================
rerunDataFormat <- FALSE
if(rerunDataFormat) unlink(file.path("dataset", "processed_data", version_label, "dmodel.RData"))

n_per_point <- 5 #Number of points to sample for each individual point. 
ID_buffer_m <- 50 # Unit is in kms
source(paste0("pipeline/", version_label,"/01_dataPreparation/01_prepData.R"))
source(paste0("pipeline/", version_label,"/01_dataPreparation/02_prepPredictionRaster.R"))
source(paste0("pipeline/", version_label,"/01_dataPreparation/03_predictionCoords.R"))

if (checks) {
  source(paste0("pipeline/", version_label,"/02_modelFit/01_exploration_analysis.R"))
}

# ============================================================
# 6. Create parameter grid
# ============================================================

pars <- expand.grid(
  sex = levels(dmodel$sex),
  Year = levels(dmodel$Year),
  harvesttimeSK24 = levels(dmodel$harvesttimeSK24),
  method = c("ensemble","gam"),
  stringsAsFactors = FALSE
)

# Harvest variables that fit_predict_habitat_model will loop over
# internally. Both are computed in the same call and returned as
# separate named entries in the results list.
harvest_vars_to_use <- c("harvested50m", "harvested200m")

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
    OR_50m = (used1_h1 * used0_h0) /
      (used0_h1 * used1_h0),
    logOR_50m = log(OR_50m),
    OR_200m = (used1_h1_200 * used0_h0_200) /
      (used0_h1_200 * used1_h0_200),
    logOR_200m = log(OR_200m)
  ) %>%
  sf::st_drop_geometry()



if (tpr < 1 || tpr > nrow(pars)) {
  stop("tpr is outside the parameter grid range.")
}

#for(tpr in 13:24){
pars_to_use <- pars[tpr, ]

pars_to_use_list <- list(
  sex = pars_to_use$sex,
  Year = pars_to_use$Year,
  harvesttimeSK24 = pars_to_use$harvesttimeSK24
)

cat("Selected parameters:\n")
print(pars_to_use)


# ============================================================
# 7. Create logging paths
# ============================================================

path_for_logs <- file.path(
  "logs", 
  version_label,
  pars_to_use$method,
  pars_to_use$Year,
  pars_to_use$harvesttimeSK24,
  pars_to_use$sex
)

dir.create(
  path_for_logs,
  recursive = TRUE,
  showWarnings = FALSE
)

logfile <- file.path(
  path_for_logs,
  paste0(
    "model_run_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".log"
  )
)


# ============================================================
# 8. Run model with logfile
# ============================================================

log_con <- file(logfile, open = "wt")

sink(log_con, split = TRUE)
sink(log_con, type = "message")

tryCatch(
  
  {
    cat("Started:", as.character(Sys.time()), "\n\n")
    cat("Selected parameters:\n")
    print(pars_to_use)
    
    cat("Odds ratios (50m and 200m harvest buffers):\n")
    rt <- or_results %>%
      dplyr::filter(sex == pars_to_use$sex,
                    Year == pars_to_use$Year,
                    harvesttimeSK24 == pars_to_use$harvesttimeSK24)
    print(rt)
    # --------------------------------------------------------
    # Fit model and predict over study region
    #
    # fit_predict_habitat_model now fits one model per entry in
    # harvest_vars_to_use (e.g. "harvested50m" and "harvested200m")
    # and returns a named list of result objects, one per harvest
    # variable: model_results[["harvested50m"]], 
    # model_results[["harvested200m"]].
    # --------------------------------------------------------
    model_results <- fit_predict_habitat_model(
      dmodel = dmodel,
      pred_grid = pred_grid,
      study_region = study_region,
      env_stack_10m = env_stack_10m,
      raster_results_23 = raster_results_23,
      raster_results_24 = raster_results_24,
      pars_to_use = pars_to_use_list,
      method = pars_to_use$method,
      response_col = "random",
      checks = checks,
      nBasis = 50,
      harvest_vars = harvest_vars_to_use
    )
    
    # --------------------------------------------------------
    # Per-harvest-variable summaries
    # --------------------------------------------------------
    
    for (hv in names(model_results)) {
      
      cat("\n==============================================\n")
      cat("Model summary for harvest variable:", hv, "\n")
      cat("==============================================\n")
      
      res_hv <- model_results[[hv]]
      
      model_data_hv <- res_hv$model_data
      print(table(model_data_hv$response, model_data_hv$harvest_var))
      
      if (pars_to_use$method == "gam") {
        print(summary(res_hv$model))
      } else if (pars_to_use$method == "maxent") {
        print(summary(res_hv$model$betas))
      } else if (pars_to_use$method == "ensemble") {
        # No simple coefficient summary for biomod2 ensembles;
        # evaluation metrics are available via biomod2::get_evaluations()
        cat("Ensemble model fitted for", hv, "\n")
      }
    }
    
    # --------------------------------------------------------
    # Save outputs
    #
    # output_dir is shared; save_habitat_results detects the
    # named-list structure and writes harvest-variable-specific
    # prefixes automatically (e.g. model_fit_harvested50m_*,
    # model_fit_harvested200m_*).
    # --------------------------------------------------------
    
    output_dir <- file.path(
      "results",
      version_label,
      pars_to_use$method,
      pars_to_use$Year,
      pars_to_use$harvesttimeSK24,
      pars_to_use$sex
    )
    
    saved_files <- save_habitat_results(
      results = model_results,
      out_dir = output_dir,
      prefix = "model_fit"
    )
    
    cat("\nSaved files:\n")
    print(saved_files)
    
    # --------------------------------------------------------
    # Optional ensemble harvest-effect plots, one per harvest
    # variable
    # --------------------------------------------------------
    
    if (pars_to_use$method == "ensemble") {
      
      for (hv in names(model_results)) {
        
        res_hv <- model_results[[hv]]
        
        if (!is.null(res_hv$harvest_conditional_plot)) {
          
          ggsave(res_hv$harvest_conditional_plot,
                 filename = file.path(
                   output_dir,
                   paste0("harvested_suitability_", hv, ".pdf")
                 ),
                 width = 8,
                 height = 6)
          
          ggsave(res_hv$harvest_conditional_plot,
                 filename = file.path(
                   output_dir,
                   paste0("harvested_suitability_", hv, ".png")
                 ),
                 width = 8,
                 height = 6,
                 dpi = 300)
        }
      }
    }
    
    cat("\nFinished:", as.character(Sys.time()), "\n")
  },
  
  error = function(e) {
    cat("\nERROR:\n")
    cat(conditionMessage(e), "\n")
    stop(e)
  },
  
  finally = {
    sink(type = "message")
    sink()
    close(log_con)
  }
)

# # ============================================================
# # 6. Create parameter grid
# # ============================================================
# 
# pars <- expand.grid(
#   sex = levels(dmodel$sex),
#   Year = levels(dmodel$Year),
#   harvesttimeSK24 = levels(dmodel$harvesttimeSK24),
#   method = c("ensemble","gam"),
#   stringsAsFactors = FALSE
# )
# 
# or_results <- dmodel %>%
#   group_by(sex, Year, harvesttimeSK24) %>%
#   summarise(
#     used0_h0 = sum(used == 0 & harvested50m == 0),
#     used1_h0 = sum(used == 1 & harvested50m == 0),
#     used0_h1 = sum(used == 0 & harvested50m == 1),
#     used1_h1 = sum(used == 1 & harvested50m == 1),
#     ids = length(unique(ID)),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     OR = (used1_h1 * used0_h0) /
#       (used0_h1 * used1_h0),
#     logOR = log(OR)
#   ) %>%
#   sf::st_drop_geometry()
# 
# 
# 
# if (tpr < 1 || tpr > nrow(pars)) {
#   stop("tpr is outside the parameter grid range.")
# }
# 
# #for(tpr in 13:24){
# pars_to_use <- pars[tpr, ]
# 
# pars_to_use_list <- list(
#   sex = pars_to_use$sex,
#   Year = pars_to_use$Year,
#   harvesttimeSK24 = pars_to_use$harvesttimeSK24
# )
# 
# cat("Selected parameters:\n")
# print(pars_to_use)
# 
# 
# # ============================================================
# # 7. Create logging paths
# # ============================================================
# 
# path_for_logs <- file.path(
#   "logs", 
#   version_label,
#   pars_to_use$method,
#   pars_to_use$Year,
#   pars_to_use$harvesttimeSK24,
#   pars_to_use$sex
# )
# 
# dir.create(
#   path_for_logs,
#   recursive = TRUE,
#   showWarnings = FALSE
# )
# 
# logfile <- file.path(
#   path_for_logs,
#   paste0(
#     "model_run_",
#     format(Sys.time(), "%Y%m%d_%H%M%S"),
#     ".log"
#   )
# )
# 
# 
# # ============================================================
# # 8. Run model with logfile
# # ============================================================
# 
# log_con <- file(logfile, open = "wt")
# 
# sink(log_con, split = TRUE)
# sink(log_con, type = "message")
# 
# tryCatch(
#   
#   {
#     cat("Started:", as.character(Sys.time()), "\n\n")
#     cat("Selected parameters:\n")
#     print(pars_to_use)
#     
#     cat("Odd ratio:\n")
#     rt <- or_results %>%
#       dplyr::filter(sex == pars_to_use$sex,
#                     Year == pars_to_use$Year,
#                     harvesttimeSK24 == pars_to_use$harvesttimeSK24)
#     print(rt)
#     # --------------------------------------------------------
#     # Fit model and predict over study region
#     # --------------------------------------------------------
#     
#     model_results <- fit_predict_habitat_model(
#       dmodel = dmodel,
#       pred_grid = pred_grid,
#       study_region = study_region,
#       env_stack_10m = env_stack_10m,
#       raster_results_23 = raster_results_23,
#       raster_results_24 = raster_results_24,
#       pars_to_use = pars_to_use_list,
#       method = pars_to_use$method,
#       response_col = "random",
#       checks = checks,
#       nBasis = 50
#     )
#     
#     cat("\nModel summary:\n")
#     model_data <- model_results$model_data
#     table(model_data$random, model_data$harvested50m)
#     
#     if(pars_to_use$method == "gam"){
#     print(summary(model_results$model))
#     } else if(pars_to_use$method == "maxent"){
#       print(summary(model_results$model$betas))
#     }  
#     
#     # --------------------------------------------------------
#     # Save outputs
#     # --------------------------------------------------------
#     
#     output_dir <- file.path(
#       "results",
#       version_label,
#       pars_to_use$method,
#       pars_to_use$Year,
#       pars_to_use$harvesttimeSK24,
#       pars_to_use$sex
#     )
#     
#     saved_files <- save_habitat_results(
#       results = model_results,
#       out_dir = output_dir,
#       prefix = "model_fit"
#     )
#     
#     cat("\nSaved files:\n")
#     print(saved_files)
#     
#     if(pars_to_use$method == "ensemble"){
#       
#       ggsave(model_results$harvest_conditional_plot,
#              filename = file.path(output_dir, "harvested_suitablility.pdf"),
#              width = 8,
#              height = 6)
#       
#       ggsave(model_results$harvest_conditional_plot,
#              filename = file.path(output_dir, "harvested_suitablility.png"),
#              width = 8,
#              height = 6,
#              dpi = 300)
#     }
#     
#     cat("\nFinished:", as.character(Sys.time()), "\n")
#   },
#   
#   error = function(e) {
#     cat("\nERROR:\n")
#     cat(conditionMessage(e), "\n")
#     stop(e)
#   },
#   
#   finally = {
#     sink(type = "message")
#     sink()
#     close(log_con)
#   }
# )
#}