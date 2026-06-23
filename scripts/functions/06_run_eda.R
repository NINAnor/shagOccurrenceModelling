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
  
  
  eda_results <- list()
  
  for (yr in years) {
    
    message("Running EDA for year: ", yr)
    
    # Keep only columns that exist
    cov_present <- intersect(cov_cols, names(dmodel))
    factor_present <- intersect(factor_vars, cov_present)
    
    data_for_model <- dmodel %>%
      dplyr::filter(Year == yr) %>%
      dplyr::select(
        dplyr::all_of(response_col),
        dplyr::all_of(cov_present),
        geometry
      )
    
    names(data_for_model)[names(data_for_model) == response_col] <- "response"
    
    if ("DistNest" %in% names(data_for_model)) {
      data_for_model <- data_for_model %>%
        dplyr::mutate(
          DistNest = as.numeric(scale(DistNest))
        )
    }
    
    data_for_model <- data_for_model %>%
      dplyr::mutate(
        response = factor(response),
        dplyr::across(
          dplyr::all_of(factor_present),
          as.factor
        )
      )
    
    # Missing values
    missing_values <- colSums(
      is.na(sf::st_drop_geometry(data_for_model))
    )
    
    # Frequency tables for categorical variables
    frequency_tables <- lapply(
      c("response", factor_present),
      function(v) table(data_for_model[[v]])
    )
    
    names(frequency_tables) <- c("response", factor_present)
    
    # Cross tables: response by each factor
    cross_tables <- lapply(
      factor_present,
      function(v) table(data_for_model$response, data_for_model[[v]])
    )
    
    names(cross_tables) <- paste0("response_by_", factor_present)
    
    # Plots
    plots <- list()
    
    if ("DistNest" %in% names(data_for_model)) {
      
      plots$distnest_hist <- ggplot(data_for_model, aes(x = DistNest)) +
        geom_histogram(bins = 30, fill = "steelblue", color = "black") +
        theme_minimal() +
        labs(
          title = paste("Distribution of scaled distance to nest,", yr),
          x = "Scaled DistNest",
          y = "Count"
        )
      
      plots$distnest_box <- ggplot(
        data_for_model,
        aes(x = response, y = DistNest, fill = response)
      ) +
        geom_boxplot(alpha = 0.7) +
        theme_minimal() +
        labs(
          title = paste("Distance to nest by used/available status,", yr),
          x = "Used / available",
          y = "Scaled DistNest",
          fill = "Used"
        )
    }
    
    if ("harvested50m" %in% names(data_for_model)) {
      plots$harvest50 <- ggplot(
        data_for_model,
        aes(x = harvested50m, fill = response)
      ) +
        geom_bar(position = "fill") +
        theme_minimal() +
        labs(
          title = paste("Harvest status within 50 m,", yr),
          x = "Harvest status",
          y = "Proportion",
          fill = "Used"
        )
    }
    
    if ("harvested200m" %in% names(data_for_model)) {
      plots$harvest200 <- ggplot(
        data_for_model,
        aes(x = harvested200m, fill = response)
      ) +
        geom_bar(position = "fill") +
        theme_minimal() +
        labs(
          title = paste("Harvest status within 200 m,", yr),
          x = "Harvest status",
          y = "Proportion",
          fill = "Used"
        )
    }
    
    if ("Logger_model" %in% names(data_for_model)) {
      plots$logger <- ggplot(
        data_for_model,
        aes(x = Logger_model, fill = response)
      ) +
        geom_bar(position = "dodge") +
        theme_minimal() +
        labs(
          title = paste("Logger model frequencies,", yr),
          x = "Logger model",
          y = "Count",
          fill = "Used"
        )
    }
    
    # Pairwise plot
    data_for_pairs <- data_for_model %>%
      sf::st_drop_geometry() %>%
      tidyr::drop_na()
    
    plots$pairs <- GGally::ggpairs(
      data_for_pairs,
      aes(color = response, alpha = 0.5)
    )
    
    eda_results[[yr]] <- list(
      data_for_model = data_for_model,
      missing_values = missing_values,
      frequency_tables = frequency_tables,
      cross_tables = cross_tables,
      plots = plots
    )
  }
  
  return(eda_results)
}