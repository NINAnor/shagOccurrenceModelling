# .Rprofile (project root)
# Read project-local .Renviron if present
if (file.exists(".Renviron")) readRenviron(".Renviron")

# Runtime flag
run_env <- Sys.getenv("RUN_ENV", unset = NA_character_)

# Determine project root: assume this file is in the project root and WD will be project root
# Use normalizePath to get absolute path (works cross-platform)
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

# Auto-detect RUN_ENV if not set
if (is.na(run_env) || !nzchar(run_env)) {
  host <- tolower(Sys.info()[["nodename"]])
  if (!is.na(host) && grepl("liprstudio", host)) {
    run_env <- "nina"
  } else {
    run_env <- "personal"
  }
}

# Choose temp_dir depending on run_env; allow override via TEMP_DIR_NINA
if (tolower(run_env) == "nina") {
  temp_dir <- Sys.getenv("TEMP_DIR_NINA", unset = "/data/scratch/kwaku")
} else {
  temp_dir <- tempdir()
}

# prefer DATA_DIR environment variable if present, otherwise fall back to proj("dataset")
# path to the results directory
if (tolower(run_env) == "nina") {
  data_dir <- "~/Mounts/R/Prosjekter/Toppskarv/Analyses_Nina_KelpHarvest/kwaku/dataset"
} else {
  data_dir_env <- "C:/ninaProjects/ninaDernhard" #Sys.getenv("DATA_DIR", unset = NA_character_)
  data_dir <- path.expand(data_dir_env)
}



# path to the results directory
if (tolower(run_env) == "nina") {
  result_dir <- "~/Mounts/R/Prosjekter/Toppskarv/Analyses_Nina_KelpHarvest/kwaku/results"
} else {
  result_dir <- "R:/Prosjekter/Toppskarv/Analyses_Nina_KelpHarvest/kwaku/results"
}

# path to log files
if (tolower(run_env) == "nina") {
  logs_dir <- "~/Mounts/R/Prosjekter/Toppskarv/Analyses_Nina_KelpHarvest/kwaku/logs"
} else {
  logs_dir <- "R:/Prosjekter/Toppskarv/Analyses_Nina_KelpHarvest/kwaku/logs"
}

# Convenience helpers
proj <- function(...) {
  if (nargs() == 0L) return(normalizePath(project_dir, winslash = "/", mustWork = FALSE))
  args <- list(...)
  normalizePath(do.call(file.path, c(project_dir, args)), winslash = "/", mustWork = FALSE)
}

goto <- function(...) {
  target <- if (nargs() == 0L) proj() else proj(...)
  setwd(target)
  message("Working directory: ", getwd())
  invisible(target)
}

# Expose variables to the session and print a short startup message
assign("project_dir", project_dir, envir = .GlobalEnv)
assign("temp_dir", temp_dir, envir = .GlobalEnv)
assign("data_dir", data_dir, envir = .GlobalEnv)
assign("result_dir", result_dir, envir = .GlobalEnv)
assign("logs_dir", logs_dir, envir = .GlobalEnv)
message("data_dir: ", data_dir)
message("Project root: ", project_dir, "  |  RUN_ENV: ", run_env, "  |  temp_dir: ", temp_dir)
options(myproject = project_dir)
