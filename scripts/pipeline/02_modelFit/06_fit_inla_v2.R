# Fit LGCP v2 (35 km water region, kelp presence covariate, time-specific
# harvest covariate) for all Year x sex x harvest-period strata and summarise
# the harvest coefficient.
#
# Requires in the environment (from masterScript / 01_prepData.R):
#   colony, dloc23, dloc, SK_35km_22_23, SK_35km_23_24, env_parquet_path,
#   modelCrs, result_dir, version_label
# and the sourced function files (14_build_mesh.R, 17_run_shag_inla_model_v2.R).

inla_v2_buffer_m <- if (exists("inla_v2_buffer_m")) inla_v2_buffer_m else 50
inla_v2_grid_res <- if (exists("inla_v2_grid_res")) inla_v2_grid_res else 50

gebco_nc <- normalizePath(
  file.path(
    here::here(), "..", "..", "Analyses_Nina_multiple_Sites",
    "Environmental Covariates", "gebco_2019_n85.0_s50.0_w0.0_e40.0.nc"
  ),
  winslash = "/"
)

env_stack_10m <- terra::rast(env_parquet_path)

cov35 <- build_inla_grid(
  colony = colony,
  tareskog = env_stack_10m$tareskog,
  gebco_path = gebco_nc,
  modelCrs = modelCrs,
  buffer_m = 35000,
  res = inla_v2_grid_res
)

dives_all <- prepare_dives(
  rbind(
    dloc23[, c("ID", "sex", "Logger_model", "UTC_datetime", "Year")],
    dloc[, c("ID", "sex", "Logger_model", "UTC_datetime", "Year")]
  ),
  modelCrs = modelCrs
)

harv_list <- list("2023" = SK_35km_22_23, "2024" = SK_35km_23_24)

inla_v2_out_dir <- file.path(
  result_dir, version_label, "inla_v2", paste0("buffer", inla_v2_buffer_m, "m")
)
dir.create(inla_v2_out_dir, recursive = TRUE, showWarnings = FALSE)

inla_v2_res <- run_shag_inla_v2_all(
  dives_all, harv_list, cov35,
  harvest_buffer_m = inla_v2_buffer_m,
  out_dir = inla_v2_out_dir,
  kelp_type = "presence"
)

inla_v2_summary <- summarise_harvest_effect(
  inla_v2_res, out_dir = inla_v2_out_dir, buffer_m = inla_v2_buffer_m
)
print(inla_v2_summary$table[, c("Year", "sex", "period", "mean", "q025", "q975",
                                "pct_dives_in_harvest", "conclusion")])
