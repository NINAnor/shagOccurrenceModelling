fit_predict_habitat_model <- function(dmodel,
                                      pred_grid,
                                      study_region,
                                      env_stack_10m,
                                      raster_results_23,
                                      raster_results_24,
                                      pars_to_use,
                                      method = c("gam", "maxent", "ensemble"),
                                      response_col = "used",
                                      checks = FALSE,
                                      block_size_m = 50000,
                                      n_cv_folds = 5,
                                      seed = 123,
                                      nBasis = 10,
                                      harvest_vars = c("harvested50m", "harvested200m")) {
  
  # ==========================================================
  # PURPOSE
  # ==========================================================
  # Fit habitat selection / suitability models using used/random
  # point data, spatial block cross-validation, individual ID effects,
  # and prediction over a spatial grid.
  #
  # Supported methods:
  #   - "gam"      : mgcv GAM with ID random effect
  #   - "maxent"   : maxnet MaxEnt with ID dummy effects
  #   - "ensemble" : biomod2 ensemble with ID dummy effects
  #
  # Important reproducibility choices:
  #   - Continuous predictors are scaled using model data only.
  #   - Spatial CV blocks are created once and reused.
  #   - seed controls block assignment reproducibility.
  #   - ID effects are included during fitting but excluded from predictions.
  #   - Predictions are masked to tareskog == 1.
  #
  # harvest_vars:
  #   - Character vector of harvest covariate names to fit, one model
  #     per element (e.g. c("harvested50m", "harvested200m")).
  #   - For each entry, that variable is used as the "harvest" covariate
  #     in the model formula / predictors. Spatial blocks, scaling
  #     parameters, and ID dummies are computed once and shared across
  #     all entries so results are directly comparable.
  #   - The function returns a named list of result objects, one per
  #     harvest_vars entry, each with the same structure previously
  #     returned for a single model.
  #
  # Based on the modelling workflow you shared. 
  # :contentReference[oaicite:0]{index=0}
  # ==========================================================
  
  method <- match.arg(method)
  
  if (!is.character(harvest_vars) || length(harvest_vars) == 0) {
    stop("harvest_vars must be a non-empty character vector.")
  }
  
  # ==========================================================
  # 0. Package checks
  # ==========================================================
  # Fail early if required packages are missing.
  
  if (!requireNamespace("blockCV", quietly = TRUE)) {
    stop("Package 'blockCV' is required.")
  }
  
  if (method == "maxent" && !requireNamespace("maxnet", quietly = TRUE)) {
    stop("Package 'maxnet' is required for method = 'maxent'.")
  }
  
  if (method == "ensemble" && !requireNamespace("biomod2", quietly = TRUE)) {
    stop("Package 'biomod2' is required for method = 'ensemble'.")
  }
  
  if (!"tareskog" %in% names(env_stack_10m)) {
    stop("env_stack_10m must contain a raster layer named 'tareskog'.")
  }
  
  # ==========================================================
  # 1. Filter model data
  # ==========================================================
  # Restrict analysis to one sex × year × harvest-period combination.
  # pars_to_use must contain:
  #   - sex
  #   - Year
  #   - harvesttimeSK24
  
  cat("Preparing model data...\n")
  
  model_sf <- dmodel %>%
    dplyr::filter(
      sex == pars_to_use$sex,
      Year == pars_to_use$Year,
      harvesttimeSK24 == pars_to_use$harvesttimeSK24
    ) %>%
    dplyr::select(-geom_id)%>%
    dplyr::mutate(geom_wkt = sf::st_as_text(geometry)) %>%
    dplyr::distinct(geom_wkt, .keep_all = TRUE) %>%
    dplyr::select(-geom_wkt)
  
  
  # ==========================================================
  # 2. Select year-specific harvest raster stack and buffer polygons
  # ==========================================================
  # Uses separate harvest layers for 2023 and 2024.
  #
  # harvest50_sf / harvest200_sf are built here from the 50m and 200m
  # harvest buffer polygons (harvest50_buffer, harvest200_buffer),
  # for both years. Each polygon set has its own indicator column
  # (harvested50m / harvested200m), = 1L within that buffer. These
  # are joined to the prediction grid in step 4 to derive
  # harvested50m / harvested200m for pred_grid_env.
  
  if (pars_to_use$Year == "2023") {
    
    harvest_stack <- raster_results_23$raster_stack
    harvest_regions <- raster_results_23$harvest50_buffer
    harvest200_regions <- raster_results_23$harvest200_buffer
    
  } else if (pars_to_use$Year == "2024") {
    
    harvest_stack <- raster_results_24$raster_stack
    harvest_regions <- raster_results_24$harvest50_buffer
    harvest200_regions <- raster_results_24$harvest200_buffer
    
  } else {
    stop("pars_to_use$Year must be either '2023' or '2024'.")
  }
  
  harvest50_sf <- sf::st_sf(
    harvested50m = 1L,
    geometry = harvest_regions,
    crs = sf::st_crs(harvest_regions)
  )
  
  harvest200_sf <- sf::st_sf(
    harvested200m = 1L,
    geometry = harvest200_regions,
    crs = sf::st_crs(harvest200_regions)
  )
  
  # ==========================================================
  # 3. Extract covariates at model points
  # ==========================================================
  # Environmental and harvest rasters are extracted at used/random points.
  # hard_substrate and sediment_substrate are intentionally excluded
  # from later model covariates.
  
  env_values_model <- terra::extract(
    env_stack_10m,
    terra::vect(model_sf)
  ) %>%
    dplyr::select(-ID)
  
  
  
  # harvest_values_model <- terra::extract(
  #   harvest_stack,
  #   terra::vect(model_sf)
  # ) %>%
  #   dplyr::select(
  #     -ID
  #   )
  
  model_sf <- model_sf %>%
    # dplyr::select(-dplyr::any_of(c("DistNest", "harvested50m", "harvested200m"))) %>%
    dplyr::bind_cols(
      env_values_model
    )
  table(model_sf$used, model_sf$harvested50m)
  if ("harvested200m" %in% names(model_sf)) {
    table(model_sf$used, model_sf$harvested200m)
  }
  # ==========================================================
  # 4. Extract covariates at prediction grid
  # ==========================================================
  # Prediction grid receives the same raster covariates as model points.
  
  env_values_pred <- terra::extract(
    env_stack_10m,
    terra::vect(pred_grid)
  ) %>%
    dplyr::select(-ID) %>%
    # Drop any harvested50m/harvested200m layers picked up from
    # env_stack_10m here - these are determined for the prediction
    # grid via harvest_values_pred below, to avoid duplicate-column
    # conflicts in the bind_cols() call.
    dplyr::select(-dplyr::any_of(c("harvested50m", "harvested200m")))
  
  # Ensure prediction grid is sf
  pred_grid_sf <- sf::st_as_sf(pred_grid)
  
  # --------------------------------------------------------
  # Extract harvest buffer indicators (harvested50m, harvested200m)
  # for the prediction grid via spatial join against the buffer
  # polygons (harvest50_sf, harvest200_sf), built in step 2 from
  # raster_results_*$harvest50_buffer and $harvest200_buffer.
  #
  # A grid cell gets harvested50m = 1 if it intersects harvest50_sf,
  # and harvested200m = 1 if it intersects harvest200_sf; 0 otherwise.
  # --------------------------------------------------------
  
  join_harvest_indicator <- function(grid_sf, polygons_sf, col_name) {
    
    joined <- sf::st_join(
      grid_sf,
      polygons_sf[, col_name, drop = FALSE],
      join = sf::st_intersects,
      left = TRUE
    ) %>%
      sf::st_drop_geometry()
    
    dplyr::if_else(
      is.na(joined[[col_name]]),
      0L,
      as.integer(joined[[col_name]])
    )
  }
  
  harvest_values_pred <- tibble::tibble(
    harvested50m = join_harvest_indicator(pred_grid_sf, harvest50_sf, "harvested50m"),
    harvested200m = join_harvest_indicator(pred_grid_sf, harvest200_sf, "harvested200m")
  )
  
  distNest_pred <- terra::extract(
    harvest_stack,
    terra::vect(pred_grid)
  ) %>%
    dplyr::select(-dplyr::any_of(c("ID", "harvested50m", "harvested200m")))
  
  # 
  # harvest_values_model <- sf::st_join(
  #   model_sf %>%
  #     dplyr::select(-c(harvested50m)),
  #   harvest_regions_sf,
  #   join = sf::st_intersects,
  #   left = TRUE
  # ) %>%
  #   dplyr::mutate(
  #     harvested50m = dplyr::if_else(
  #       is.na(harvested50m),
  #       0L,
  #       harvested50m
  #     )
  #   )
  
  pred_grid_env <- pred_grid %>%
    dplyr::select(-dplyr::any_of(c("harvested50m", "harvested200m"))) %>%
    dplyr::bind_cols(
      env_values_pred,
      harvest_values_pred,
      distNest_pred
    )
  
  # ==========================================================
  # 5. Prepare modelling dataframe
  # ==========================================================
  # Keep only the variables used in the model and drop incomplete rows.
  # Response must be binary:
  #   used/random == 1/0, or equivalent response_col.
  #
  # Both harvested50m and harvested200m are retained here (if present)
  # so that either can be selected as the active harvest covariate
  # further down, on a per-harvest_vars-entry basis.
  
  model_sf_clean <- model_sf %>%
    dplyr::select(
      ID,
      Logger_model,
      dplyr::all_of(response_col),
      Year,
      sex,
      harvesttimeSK24,
      DistNest,
      dplyr::any_of(c("harvested50m", "harvested200m")),
      bathymetry,
      tareskog,
      geometry
    ) %>%
    tidyr::drop_na()
  
  coords_model <- sf::st_coordinates(model_sf_clean)
  
  model_data_base <- model_sf_clean %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(
      x = coords_model[, 1],
      y = coords_model[, 2],
      
      # Convert response to numeric 0/1.
      response = as.numeric(as.character(.data[[response_col]])),
      
      # Factor variables retained as categorical predictors.
      Logger_model = factor(Logger_model),
      # tareskog = factor(tareskog),
      ID = factor(ID)
    )
  
  # Convert harvest buffer variables to factors where present, without
  # yet committing to one as "the" active covariate.
  for (hv in intersect(c("harvested50m", "harvested200m"), names(model_data_base))) {
    model_data_base[[hv]] <- factor(model_data_base[[hv]])
  }
  
  stopifnot(all(model_data_base$response %in% c(0, 1)))
  
  # ==========================================================
  # 6. Scale continuous predictors
  # ==========================================================
  # Scaling parameters are estimated from model data only.
  # This avoids leaking prediction-grid information into model fitting.
  # These scaling parameters do not depend on which harvest variable
  # is used, so they are computed once and reused across harvest_vars.
  
  dist_center <- mean(model_data_base$DistNest, na.rm = TRUE)
  dist_scale <- sd(model_data_base$DistNest, na.rm = TRUE)
  
  bathy_center <- mean(model_data_base$bathymetry, na.rm = TRUE)
  bathy_scale <- sd(model_data_base$bathymetry, na.rm = TRUE)
  
  x_center <- mean(model_data_base$x, na.rm = TRUE)
  x_scale <- sd(model_data_base$x, na.rm = TRUE)
  
  y_center <- mean(model_data_base$y, na.rm = TRUE)
  y_scale <- sd(model_data_base$y, na.rm = TRUE)
  
  model_data_base <- model_data_base %>%
    dplyr::mutate(
      DistNest_sc = as.numeric((DistNest - dist_center) / dist_scale),
      bathymetry_sc = as.numeric((bathymetry - bathy_center) / bathy_scale),
      
      # Coordinate scaling for MaxEnt and biomod spatial trend terms.
      x_sc = as.numeric((x - x_center) / x_scale),
      y_sc = as.numeric((y - y_center) / y_scale),
      x2 = x_sc^2,
      y2 = y_sc^2,
      xy = x_sc * y_sc
    )
  
  # ==========================================================
  # 7. Individual ID effect variables
  # ==========================================================
  # GAM uses s(ID, bs = "re").
  # MaxEnt and biomod cannot fit random effects directly, so ID is
  # represented using dummy variables during fitting.
  #
  # For population-level prediction, all ID dummy variables are set to 0.
  # ID dummies do not depend on harvest_vars, so computed once.
  
  id_dummy_model <- stats::model.matrix(~ ID - 1, data = model_data_base)
  colnames(id_dummy_model) <- make.names(colnames(id_dummy_model))
  
  model_data_base <- dplyr::bind_cols(
    model_data_base,
    as.data.frame(id_dummy_model)
  )
  
  id_dummy_cols <- colnames(id_dummy_model)
  
  # ==========================================================
  # 8. Spatial block cross-validation
  # ==========================================================
  # Spatial blocks are created once and reused for ALL harvest_vars
  # entries and all model types, so that the same fold assignment
  # applies regardless of which harvest covariate is used.
  
  cat("Creating spatial blocks...\n")
  
  model_sf_cv <- sf::st_as_sf(
    model_data_base,
    coords = c("x", "y"),
    crs = sf::st_crs(study_region),
    remove = FALSE
  )
  
  set.seed(seed)
  
  spatial_blocks <- blockCV::cv_spatial(
    x = model_sf_cv,
    column = "response",
    r = env_stack_10m,
    #size = block_size_m,
    k = n_cv_folds,
    selection = "random",
    iteration = 100
  )
  
  model_data_base$fold <- spatial_blocks$folds_ids
  
  # biomod2 expects:
  #   TRUE  = training/calibration
  #   FALSE = testing/evaluation
  cv_table <- sapply(seq_len(n_cv_folds), function(k) {
    model_data_base$fold != k
  })
  
  colnames(cv_table) <- paste0("_allData_RUN", seq_len(n_cv_folds))
  cv_table <- as.data.frame(cv_table)
  
  # ==========================================================
  # 9-17. Fit + predict + plot, once per harvest_vars entry
  # ==========================================================
  # The block below was previously a single fit/predict/plot section.
  # It is now wrapped in a loop so it runs once per harvest variable
  # (e.g. "harvested50m", then "harvested200m"), with each iteration
  # producing its own fitted model, predictions, and plots. The
  # harvest variable used in that iteration is referred to generically
  # as `harvest_var_active` inside the loop, and the model formula /
  # predictor sets reference that column directly (internally renamed
  # to `harvest_var`).
  
  results_by_harvest <- list()
  
  for (harvest_var_active in harvest_vars) {
    
    cat("=== Fitting", method, "model using harvest variable:",
        harvest_var_active, "===\n")
    
    if (!harvest_var_active %in% names(model_data_base)) {
      warning(
        paste0(
          "Harvest variable '", harvest_var_active,
          "' not found in model data; skipping this run."
        )
      )
      next
    }
    
    if (!harvest_var_active %in% names(pred_grid_env)) {
      warning(
        paste0(
          "Harvest variable '", harvest_var_active,
          "' not found in pred_grid_env; skipping this run."
        )
      )
      next
    }
    
    # ----------------------------------------------------------
    # Per-iteration model_data: rename the active harvest variable
    # to a fixed internal name 'harvest_var' so the rest of the
    # fitting / prediction code does not need to branch on names.
    # ----------------------------------------------------------
    
    model_data <- model_data_base %>%
      dplyr::rename(harvest_var = dplyr::all_of(harvest_var_active))
    
    block_eval <- NULL
    ensemble_model <- NULL
    biomod_pred_data <- NULL
    fitted_model <- NULL
    
    if (method == "gam") {
      
      # --------------------------------------------------------
      # GAM
      # --------------------------------------------------------
      # Includes:
      #   - harvest effect (harvest_var_active)
      #   - distance to nest
      #   - bathymetry
      #   - tareskog
      #   - spatial smooth
      #   - individual random effect
      
      fitted_model <- mgcv::gam(
        response ~
          harvest_var +
          DistNest_sc +
          bathymetry_sc +
          s(ID, bs = "re")+
          s(x, y, bs = "tp", k = nBasis) ,
        data = model_data,
        family = binomial,
        method = "REML"
      )
      
      # Fold-wise spatial-block evaluation.
      block_eval <- purrr::map_dfr(seq_len(n_cv_folds), function(k) {
        
        train_dat <- model_data[model_data$fold != k, ]
        test_dat <- model_data[model_data$fold == k, ]
        
        m_k <- mgcv::gam(
          response ~
            harvest_var +
            DistNest_sc +
            bathymetry_sc +
            tareskog +
            s(x, y, bs = "tp", k = nBasis) +
            s(ID, bs = "re"),
          data = train_dat,
          family = binomial,
          method = "REML"
        )
        
        # Exclude ID random effect for population-level validation prediction.
        test_dat$pred <- predict(
          m_k,
          newdata = test_dat,
          type = "response",
          exclude = "s(ID)"
        )
        
        test_dat$cv_fold <- k
        test_dat
      })
      
    } else if (method == "maxent") {
      
      # --------------------------------------------------------
      # MaxEnt via maxnet
      # --------------------------------------------------------
      # Includes ID dummy variables during fitting.
      # During prediction, these are set to 0.
      
      maxent_predictors <- model_data %>%
        dplyr::select(
          DistNest_sc,
          bathymetry_sc,
          harvest_var,
          tareskog,
          x_sc,
          y_sc,
          x2,
          y2,
          xy,
          dplyr::all_of(id_dummy_cols)
        )
      
      fitted_model <- maxnet::maxnet(
        p = model_data$response,
        data = maxent_predictors,
        regmult = 0.01,
        f = maxnet::maxnet.formula(
          p = model_data$response,
          data = maxent_predictors,
          classes = "l"
        )
      )
      
      block_eval <- purrr::map_dfr(seq_len(n_cv_folds), function(k) {
        
        train_dat <- model_data[model_data$fold != k, ]
        test_dat <- model_data[model_data$fold == k, ]
        
        train_x <- train_dat %>%
          dplyr::select(
            DistNest_sc,
            bathymetry_sc,
            harvest_var,
            tareskog,
            x_sc,
            y_sc,
            x2,
            y2,
            xy,
            dplyr::all_of(id_dummy_cols)
          )
        
        test_x <- test_dat %>%
          dplyr::select(
            DistNest_sc,
            bathymetry_sc,
            harvest_var,
            tareskog,
            x_sc,
            y_sc,
            x2,
            y2,
            xy,
            dplyr::all_of(id_dummy_cols)
          )
        
        m_k <- maxnet::maxnet(
          p = train_dat$response,
          data = train_x,
          regmult = 0.01,
          f = maxnet::maxnet.formula(
            p = train_dat$response,
            data = train_x,
            classes = "l"
          )
        )
        
        # Remove ID effect from validation prediction.
        for (id_col in id_dummy_cols) {
          test_x[[id_col]] <- 0
        }
        
        test_dat$pred <- as.numeric(
          predict(
            m_k,
            newdata = test_x,
            type = "cloglog"
          )
        )
        
        test_dat$cv_fold <- k
        test_dat
      })
      
    } else if (method == "ensemble") {
      
      # --------------------------------------------------------
      # biomod2 ensemble
      # --------------------------------------------------------
      # Uses:
      #   - MAXNET
      #   - GAM
      #   - RF
      #   - GBM
      #
      # Uses the same spatial-block CV table as other methods.
      # Includes ID dummy variables during fitting.
      
      biomod_df <- model_data %>%
        dplyr::mutate(
          response = as.integer(response),
          harvest_var = as.numeric(as.character(harvest_var)),
          tareskog = as.numeric(as.character(tareskog))
        ) %>%
        dplyr::select(
          response,
          DistNest_sc,
          bathymetry_sc,
          harvest_var,
          tareskog,
          x,
          y,
          x_sc,
          y_sc,
          x2,
          y2,
          xy,
          fold,
          dplyr::all_of(id_dummy_cols)
        ) %>%
        tidyr::drop_na()
      
      biomod_cv_table <- sapply(seq_len(n_cv_folds), function(k) {
        biomod_df$fold != k
      })
      
      colnames(biomod_cv_table) <- paste0("_allData_RUN", seq_len(n_cv_folds))
      
      biomod_cv_table <- as.data.frame(biomod_cv_table)
      biomod_cv_table[] <- lapply(biomod_cv_table, as.logical)
      
      stopifnot(nrow(biomod_cv_table) == nrow(biomod_df))
      stopifnot(!anyNA(biomod_cv_table))
      
      # biomod2 writes model files to working directory.
      # Use a temporary directory so repeated runs do not overwrite project files.
      
      tmp_dir <- tempfile(pattern = "biomod_")
      dir.create(tmp_dir, recursive = TRUE)
      
      old_wd <- getwd()
      setwd(tmp_dir)
      
      on.exit({
        setwd(old_wd)
        if (dir.exists(tmp_dir)) {
          unlink(tmp_dir, recursive = TRUE, force = TRUE)
        }
      }, add = TRUE)
      
      biomod_expl_vars <- biomod_df %>%
        dplyr::select(
          DistNest_sc,
          bathymetry_sc,
          tareskog,
          harvest_var,
          x_sc,
          y_sc,
          x2,
          y2,
          xy,
          dplyr::all_of(id_dummy_cols)
        )
      
      biomod_data <- biomod2::BIOMOD_FormatingData(
        resp.var = biomod_df$response,
        expl.var = biomod_expl_vars,
        resp.xy = biomod_df %>%
          dplyr::select(x, y),
        resp.name = "habitat.model"
      )
      
      biomod_options <- biomod2::bm_ModelingOptions(
        models = c("MAXNET", "GAM", "RF", "GBM"),
        strategy = "user.defined",
        user.val = list(
          MAXNET = list(
            "_allData_allRun" = list(
              feature_classes = "l",
              regmult = 0.01
            )
          ),
          GAM = list(
            "_allData_allRun" = list(
              algo = "GAM_mgcv",
              family = binomial(link = "logit"),
              method = "REML",
              select = TRUE
            )
          ),
          RF = list(
            "_allData_allRun" = list(
              ntree = 150,
              importance = TRUE
            )
          ),
          GBM = list(
            "_allData_allRun" = list(
              n.trees = 250,
              interaction.depth = 3,
              shrinkage = 0.01,
              bag.fraction = 0.5
            )
          )
        )
      )
      
      fitted_model <- biomod2::BIOMOD_Modeling(
        bm.format = biomod_data,
        models = c("MAXNET", "GAM", "RF", "GBM"),
        OPT.user.val = biomod_options,
        
        CV.strategy = "user.defined",
        CV.user.table = biomod_cv_table,
        CV.do.full.models = FALSE,
        
        metric.eval = c("AUCroc", "TSS"),
        var.import = 0,
        nb.cpu = 1,
        modeling.id = paste0("ensemble_model_block_cv_", harvest_var_active)
      )
      
      ensemble_model <- biomod2::BIOMOD_EnsembleModeling(
        bm.mod = fitted_model,
        em.by = "all",
        em.algo = c("EMmean", "EMwmean"),
        metric.select = "AUCroc",
        metric.select.thresh = 0.7,
        metric.eval = c("AUCroc", "TSS")
      )
      
      block_eval <- biomod_df
    }
    
    # ==========================================================
    # Prepare prediction dataframe
    # ==========================================================
    # Prediction data receive same scaling as model data.
    # ID is assigned a valid factor level, but ID effect is excluded later.
    # The active harvest variable is renamed to 'harvest_var' to match
    # the model_data used for fitting.
    
    cat("Preparing prediction data...\n")
    
    coords_pred <- sf::st_coordinates(pred_grid_env)
    
    pred_df <- pred_grid_env %>%
      sf::st_drop_geometry() %>%
      dplyr::rename(harvest_var = dplyr::all_of(harvest_var_active)) %>%
      dplyr::mutate(
        x = coords_pred[, 1],
        y = coords_pred[, 2],
        harvest_var = factor(harvest_var, levels = levels(model_data$harvest_var)),
        tareskog = factor(tareskog, levels = levels(model_data$tareskog)),
        ID = factor(levels(model_data$ID)[1], levels = levels(model_data$ID)),
        DistNest_sc = as.numeric((DistNest - dist_center) / dist_scale),
        bathymetry_sc = as.numeric((bathymetry - bathy_center) / bathy_scale),
        x_sc = as.numeric((x - x_center) / x_scale),
        y_sc = as.numeric((y - y_center) / y_scale),
        x2 = x_sc^2,
        y2 = y_sc^2,
        xy = x_sc * y_sc
      )
    
    # ==========================================================
    # Predict without individual ID effects
    # ==========================================================
    
    cat("Running predictions...\n")
    
    if (method == "gam") {
      
      pred_df$pred_prob <- predict(
        fitted_model,
        newdata = pred_df,
        type = "response",
        exclude = "s(ID)"
      )
      
    } else if (method == "maxent") {
      
      maxent_pred_data <- pred_df %>%
        dplyr::select(
          DistNest_sc,
          bathymetry_sc,
          harvest_var,
          tareskog,
          x_sc,
          y_sc,
          x2,
          y2,
          xy
        ) %>%
        as.data.frame()
      
      for (id_col in id_dummy_cols) {
        maxent_pred_data[[id_col]] <- 0
      }
      
      pred_df$pred_prob <- as.numeric(
        predict(
          fitted_model,
          newdata = maxent_pred_data,
          type = "cloglog"
        )
      )
      
    } else if (method == "ensemble") {
      
      biomod_pred_data <- pred_df %>%
        dplyr::mutate(
          harvest_var = as.numeric(as.character(harvest_var)),
          tareskog = as.numeric(as.character(tareskog))
        ) %>%
        dplyr::select(
          DistNest_sc,
          bathymetry_sc,
          tareskog,
          harvest_var,
          x_sc,
          y_sc,
          x2,
          y2,
          xy
        ) %>%
        as.data.frame()
      
      for (id_col in id_dummy_cols) {
        biomod_pred_data[[id_col]] <- 0
      }
      
      ensemble_projection <- biomod2::BIOMOD_EnsembleForecasting(
        bm.em = ensemble_model,
        new.env = biomod_pred_data,
        proj.name = paste0("current_projection_", harvest_var_active)
      )
      
      pred_vals <- biomod2::get_predictions(ensemble_projection)
      
      pred_vals_one <- pred_vals %>%
        dplyr::filter(algo == "EMwmean") %>%
        dplyr::arrange(points)
      
      pred_df$pred_prob <- pred_vals_one$pred / 1000
    }
    
    # ==========================================================
    # Rasterise predictions and mask to kelp habitat
    # ==========================================================
    # Predictions are restricted to:
    #   - study_region
    #   - tareskog == 1
    
    pred_raster <- terra::rast(
      pred_df[, c("x", "y", "pred_prob")],
      type = "xyz",
      crs = sf::st_crs(study_region)$wkt
    )
    
    names(pred_raster) <- "pred_prob"
    
    tareskog_mask <- env_stack_10m[["tareskog"]] == 1
    
    tareskog_mask <- terra::resample(
      tareskog_mask,
      pred_raster,
      method = "near"
    )
    
    pred_raster <- terra::mask(
      pred_raster,
      tareskog_mask,
      maskvalues = 0
    )
    
    pred_raster <- terra::crop(
      pred_raster,
      terra::vect(study_region)
    )
    
    pred_raster <- terra::mask(
      pred_raster,
      terra::vect(study_region)
    )
    
    # ==========================================================
    # Prepare harvested polygons
    # ==========================================================
    # Use the buffer polygon matching the active harvest variable
    # for this iteration's plot outline (harvest50_buffer for
    # harvested50m, harvest200_buffer for harvested200m; falls back
    # to harvest_regions/50m buffer for any other harvest variable).
    
    harvest_polygons_active <- switch(
      harvest_var_active,
      harvested50m = harvest_regions,
      harvested200m = harvest200_regions,
      harvest_regions
    )
    
    harvest_regions_sf_run <- sf::st_as_sf(harvest_polygons_active)
    
    if (sf::st_crs(harvest_regions_sf_run) != sf::st_crs(study_region)) {
      harvest_regions_sf_run <- sf::st_transform(
        harvest_regions_sf_run,
        sf::st_crs(study_region)
      )
    }
    
    harvest_regions_sf_run <- sf::st_intersection(
      harvest_regions_sf_run,
      sf::st_union(study_region)
    )
    
    # ==========================================================
    # Prediction plot
    # ==========================================================
    # Plot contains:
    #   - masked predicted suitability raster
    #   - harvested polygons
    #   - used/random model points
    #   - study-region boundary
    
    pred_df_plot <- terra::as.data.frame(
      pred_raster,
      xy = TRUE,
      na.rm = TRUE
    )
    
    names(pred_df_plot)[3] <- "pred_prob"
    
    points_plot <- model_sf_clean %>%
      dplyr::mutate(
        point_type = dplyr::case_when(
          .data[[response_col]] == 1 ~ "Used",
          .data[[response_col]] == 0 ~ "Random",
          TRUE ~ NA_character_
        )
      )
    
    prediction_plot <- ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = pred_df_plot,
        ggplot2::aes(x = x, y = y, fill = pred_prob)
      ) +
      ggplot2::geom_sf(
        data = harvest_regions_sf_run,
        fill = NA,
        color = "red",
        linewidth = 0.7
      ) +
      ggplot2::geom_sf(
        data = points_plot,
        ggplot2::aes(color = point_type),
        size = 0.3,
        alpha = 0.1
      ) +
      ggplot2::geom_sf(
        data = study_region,
        fill = NA,
        color = "black",
        linewidth = 0.4
      ) +
      ggplot2::coord_sf(crs = sf::st_crs(study_region)) +
      ggplot2::scale_fill_viridis_c(
        name = "Predicted probability",
        na.value = NA
      ) +
      ggplot2::scale_color_manual(
        name = "Point type",
        values = c(
          "Used" = "black",
          "Random" = "grey60"
        )
      ) +
      ggplot2::labs(
        title = paste("Predicted habitat suitability:", toupper(method)),
        subtitle = paste(
          pars_to_use$Year,
          pars_to_use$sex,
          pars_to_use$harvesttimeSK24,
          "| harvest var:", harvest_var_active,
          "| masked to tareskog == 1"
        ),
        x = "Easting",
        y = "Northing"
      ) +
      ggplot2::theme_minimal()
    
    # ==========================================================
    # Optional diagnostics
    # ==========================================================
    
    model_checks <- NULL
    
    if (checks && method == "gam") {
      model_checks <- list(
        summary = summary(fitted_model),
        gam_check = mgcv::gam.check(fitted_model)
      )
    }
    
    # ==========================================================
    # Optional ensemble harvest-effect plot
    # ==========================================================
    
    harvest_conditional_plot <- NULL
    
    if (method == "ensemble") {
      
      newdat0 <- biomod_pred_data
      newdat1 <- biomod_pred_data
      
      newdat0$harvest_var <- 0
      newdat1$harvest_var <- 1
      
      proj0 <- biomod2::BIOMOD_EnsembleForecasting(
        bm.em = ensemble_model,
        new.env = newdat0,
        proj.name = paste0("harvest0_", harvest_var_active)
      )
      
      proj1 <- biomod2::BIOMOD_EnsembleForecasting(
        bm.em = ensemble_model,
        new.env = newdat1,
        proj.name = paste0("harvest1_", harvest_var_active)
      )
      
      pred0 <- biomod2::get_predictions(proj0) %>%
        dplyr::filter(algo == "EMwmean") %>%
        dplyr::arrange(points)
      
      pred1 <- biomod2::get_predictions(proj1) %>%
        dplyr::filter(algo == "EMwmean") %>%
        dplyr::arrange(points)
      
      effect_df <- data.frame(
        x = pred_df$x,
        y = pred_df$y,
        no_harvest = pred0$pred / 1000,
        harvest = pred1$pred / 1000
      )
      
      effect_df$delta <- effect_df$harvest - effect_df$no_harvest
      
      delta_raster <- terra::rast(
        effect_df[, c("x", "y", "delta")],
        type = "xyz",
        crs = sf::st_crs(study_region)$wkt
      )
      
      delta_raster <- terra::mask(
        delta_raster,
        tareskog_mask,
        maskvalues = 0
      )
      
      delta_raster <- terra::crop(
        delta_raster,
        terra::vect(study_region)
      )
      
      delta_raster <- terra::mask(
        delta_raster,
        terra::vect(study_region)
      )
      
      delta_df <- terra::as.data.frame(
        delta_raster,
        xy = TRUE,
        na.rm = TRUE
      )
      
      harvest_conditional_plot <- ggplot2::ggplot() +
        ggplot2::geom_raster(
          data = delta_df,
          ggplot2::aes(x = x, y = y, fill = delta)
        ) +
        ggplot2::geom_sf(
          data = harvest_regions_sf_run,
          fill = NA,
          color = "red",
          linewidth = 0.7
        ) +
        ggplot2::geom_sf(
          data = points_plot,
          ggplot2::aes(color = point_type),
          size = 0.7,
          alpha = 0.25
        ) +
        ggplot2::geom_sf(
          data = study_region,
          fill = NA,
          color = "black",
          linewidth = 0.4
        ) +
        ggplot2::coord_sf(crs = sf::st_crs(study_region)) +
        ggplot2::scale_fill_viridis_c(
          name = "Harvest effect\nΔ suitability",
          option = "magma"
        ) +
        ggplot2::scale_color_manual(
          name = "Point type",
          values = c(
            "Used" = "black",
            "Random" = "grey60"
          )
        ) +
        ggplot2::theme_minimal()
    }
    
    # ==========================================================
    # Assemble per-harvest-variable result object
    # ==========================================================
    # 'model_data' and 'pred_df' below retain the internal column
    # 'harvest_var'; for clarity, a copy of the active harvest
    # variable name is also stored under 'harvest_var_name'.
    
    results_by_harvest[[harvest_var_active]] <- list(
      method = method,
      harvest_var_name = harvest_var_active,
      parameters = pars_to_use,
      seed = seed,
      
      model = fitted_model,
      ensemble_model = ensemble_model,
      
      model_data = model_data,
      model_sf = model_sf,
      model_sf_clean = model_sf_clean,
      
      pred_grid_env = pred_grid_env,
      pred_df = pred_df,
      pred_raster = pred_raster,
      
      harvest_regions = harvest_regions_sf_run,
      plot = prediction_plot,
      harvest_conditional_plot = harvest_conditional_plot,
      
      checks = model_checks,
      
      spatial_blocks = spatial_blocks,
      cv_table = cv_table,
      block_eval = block_eval,
      block_size_m = block_size_m,
      n_cv_folds = n_cv_folds,
      
      scaling = list(
        DistNest = c(center = dist_center, scale = dist_scale),
        bathymetry = c(center = bathy_center, scale = bathy_scale),
        x = c(center = x_center, scale = x_scale),
        y = c(center = y_center, scale = y_scale)
      ),
      
      id_dummy_cols = id_dummy_cols,
      biomod_pred_data = biomod_pred_data
    )
    
  } # end loop over harvest_vars
  
  cat("Returning results for harvest variables:",
      paste(names(results_by_harvest), collapse = ", "), "\n")
  
  return(results_by_harvest)
}


# ============================================================
# Save habitat modelling outputs
# ============================================================
# save_habitat_results now expects either:
#   - a single result list (as previously returned by
#     fit_predict_habitat_model for one harvest variable), OR
#   - the new named list of result lists returned by
#     fit_predict_habitat_model (one entry per harvest_vars
#     element, e.g. results[["harvested50m"]], results[["harvested200m"]]).
#
# If a named list of result lists is passed (i.e. `results` does not
# itself contain a "method" element but its elements do), this function
# loops over each entry and saves outputs with a harvest-variable-specific
# prefix (e.g. "habitat_model_harvested50m_...",
# "habitat_model_harvested200m_...").

save_habitat_results <- function(results,
                                 out_dir,
                                 prefix = "habitat_model") {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ----------------------------------------------------------
  # Detect whether `results` is a single result object or a
  # named list of result objects (one per harvest variable).
  # ----------------------------------------------------------
  
  is_single_result <- "method" %in% names(results) && "model" %in% names(results)
  
  if (!is_single_result) {
    
    saved_files_all <- list()
    
    for (harvest_var_name in names(results)) {
      
      sub_prefix <- paste0(prefix, "_", harvest_var_name)
      
      saved_files_all[[harvest_var_name]] <- save_habitat_results(
        results = results[[harvest_var_name]],
        out_dir = out_dir,
        prefix = sub_prefix
      )
    }
    
    return(saved_files_all)
  }
  
  # ----------------------------------------------------------
  # 1. Save full result object
  # ----------------------------------------------------------
  
  save(
    results,
    file = file.path(out_dir, paste0(prefix, "_full_results.RData"))
  )
  
  
  # ----------------------------------------------------------
  # 2. Save core model objects
  # ----------------------------------------------------------
  
  method     <- results$method
  parameters <- results$parameters
  model      <- results$model
  checks     <- results$checks
  spatial_blocks <- results$spatial_blocks
  cv_table <- results$cv_table
  block_eval <- results$block_eval
  block_size_m <- results$block_size_m
  n_cv_folds <- results$n_cv_folds
  
  save(
    method,
    parameters,
    model,
    checks,
    spatial_blocks,
    cv_table,
    block_eval,
    block_size_m,
    n_cv_folds,
    file = file.path(out_dir, paste0(prefix, "_model_objects.RData"))
  )
  
  
  # ----------------------------------------------------------
  # 3. Save ordinary dataframes
  # ----------------------------------------------------------
  
  model_data <- results$model_data
  pred_df    <- results$pred_df
  
  save(
    model_data,
    pred_df,
    file = file.path(out_dir, paste0(prefix, "_dataframes.RData"))
  )
  
  
  # ----------------------------------------------------------
  # 4. Helper: save sf objects as parquet using WKT geometry
  # ----------------------------------------------------------
  
  save_sf_parquet <- function(x, filename) {
    
    x_parquet <- x %>%
      dplyr::mutate(
        geometry_wkt = sf::st_as_text(sf::st_geometry(.))
      ) %>%
      sf::st_drop_geometry()
    
    arrow::write_parquet(
      x_parquet,
      file.path(out_dir, filename)
    )
  }
  
  
  # ----------------------------------------------------------
  # 5. Save spatial vector outputs as parquet
  # ----------------------------------------------------------
  
  save_sf_parquet(
    results$model_sf,
    paste0(prefix, "_model_sf.parquet")
  )
  
  save_sf_parquet(
    results$pred_grid_env,
    paste0(prefix, "_pred_grid_env.parquet")
  )
  
  save_sf_parquet(
    results$harvest_regions,
    paste0(prefix, "_harvest_regions.parquet")
  )
  
  
  # ----------------------------------------------------------
  # 6. Save raster prediction as xyz parquet and GeoTIFF
  # ----------------------------------------------------------
  
  pred_raster_df <- terra::as.data.frame(
    results$pred_raster,
    xy = TRUE,
    na.rm = FALSE
  )
  
  arrow::write_parquet(
    pred_raster_df,
    file.path(out_dir, paste0(prefix, "_pred_raster_xyz.parquet"))
  )
  
  terra::writeRaster(
    results$pred_raster,
    file.path(out_dir, paste0(prefix, "_pred_raster.tif")),
    overwrite = TRUE
  )
  
  
  # ----------------------------------------------------------
  # 7. Save prediction plot as PDF and PNG
  # ----------------------------------------------------------
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_prediction_plot.pdf")),
    plot = results$plot,
    width = 8,
    height = 6
  )
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(prefix, "_prediction_plot.png")),
    plot = results$plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  # ----------------------------------------------------------
  # 7. Save harvest conditional prediction plot as PDF and PNG
  # ----------------------------------------------------------
  if(!is.null(results$harvest_conditional_plot)){
    ggplot2::ggsave(
      filename = file.path(out_dir, paste0(prefix, "_harvest_conditional_plot.pdf")),
      plot = results$harvest_conditional_plot,
      width = 8,
      height = 6
    )
    
    ggplot2::ggsave(
      filename = file.path(out_dir, paste0(prefix, "_harvest_conditional_plot.png")),
      plot = results$harvest_conditional_plot,
      width = 8,
      height = 6,
      dpi = 300
    )
  } 
  
  # ----------------------------------------------------------
  # 8. Return saved file paths
  # ----------------------------------------------------------
  
  saved_files <- list.files(
    out_dir,
    pattern = prefix,
    full.names = TRUE
  )
  
  return(saved_files)
}