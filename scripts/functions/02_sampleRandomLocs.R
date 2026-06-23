generate_available_points <- function(
    
  # sf object containing observed/tagged locations
  # must contain:
  # - geometry column
  # - ID
  # - sex
  # - Logger_model
  # - UTC_datetime
  dloc,
  
  # sf polygon defining the total study area
  # used when group_by_id = FALSE
  study_region,
  
  # sf polygon defining the area within which
  # available/background points are allowed to fall
  # (e.g. a "tare" mask such as land, habitat, or
  # accessibility layer)
  #
  # applied to the available region AFTER the
  # buffer step, for both group_by_id = TRUE
  # and group_by_id = FALSE
  tare,
  
  # sf point object for colony/nest location
  # used to calculate DistNest
  colony,
  
  # number of available/background points
  # generated per observed location
  # e.g. 5 = 5 background points for every used point
  n_per_point = 5,
  
  # exclusion buffer distance (metres)
  # around observed locations when:
  # group_by_id = FALSE
  #
  # prevents background points being sampled
  # too close to known used locations
  min_dist_m = 500,
  
  # buffer distance (metres)
  # around all locations from each ID when:
  # group_by_id = TRUE
  #
  # defines the accessible area for each individual
  id_buffer_m = 10000,
  
  # controls grouping strategy for background sampling
  #
  # FALSE:
  #   background sampled by:
  #   sex × Logger_model
  #
  # TRUE:
  #   background sampled separately for each ID
  #   using individual-specific buffered areas
  group_by_id = FALSE
  
) {
  
  # --------------------------------------------------
  # Label observed locations as used = 1
  # --------------------------------------------------
  
  used <- dloc %>%
    dplyr::mutate(used = 1)
  
  # --------------------------------------------------
  # Union of tare polygon(s), used to constrain
  # available regions after buffering
  # --------------------------------------------------
  
  tare_union <- sf::st_union(tare)
  
  # ==================================================
  # OPTION 1:
  # Individual-specific background sampling
  # ==================================================
  
  if (isTRUE(group_by_id)) {
    
    # --------------------------------------------------
    # Create accessible area for each ID
    #
    # Steps:
    # 1. Merge all locations from each ID
    # 2. Buffer merged geometry
    # 3. Restrict to the 'tare' region
    #
    # Result:
    # polygon representing individual accessible area,
    # constrained to fall within 'tare'
    # --------------------------------------------------
    
    avail_sf <- used %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(
        geometry = sf::st_union(geometry),
        .groups = "drop"
      ) %>%
      sf::st_buffer(id_buffer_m) %>%
      sf::st_intersection(tare_union)
    
    # --------------------------------------------------
    # Generate available points separately for each ID
    # --------------------------------------------------
    
    available_list <- purrr::map(
      unique(used$ID),
      function(id) {
        
        # observed locations for current ID
        dloc_id <- used %>%
          dplyr::filter(ID == id)
        
        # number of available points to generate
        n_id <- nrow(dloc_id) * n_per_point
        
        # accessible area for current ID
        # (already buffered and constrained to 'tare')
        area_i <- avail_sf %>%
          dplyr::filter(ID == id)
        
        # --------------------------------------------------
        # Random background sampling
        # --------------------------------------------------
        
        rnd <- sf::st_sample(
          area_i,
          size = n_id,
          type = "random",
          exact = TRUE
        )
        
        # --------------------------------------------------
        # Copy metadata from observed points
        # --------------------------------------------------
        
        meta_source <- dloc_id %>%
          sf::st_drop_geometry()
        
        meta <- meta_source %>%
          dplyr::slice(
            rep(dplyr::row_number(), each = n_per_point)
          ) %>%
          dplyr::slice(seq_len(length(rnd)))
        
        # Randomise timestamps
        meta$UTC_datetime <- sample(
          meta_source$UTC_datetime,
          size = nrow(meta),
          replace = TRUE
        )
        
        # --------------------------------------------------
        # Build sf object of available points
        # --------------------------------------------------
        
        available_group <- sf::st_sf(
          meta,
          geometry = rnd,
          crs = sf::st_crs(dloc)
        ) %>%
          dplyr::mutate(used = 0)
        
        # distance from colony/nest
        available_group$DistNest <- as.numeric(
          sf::st_distance(available_group, colony)
        )
        
        list(
          available = available_group,
          available_region = area_i,
          known_buffer = area_i
        )
      }
    )
    
    grouping_label <- "ID-specific buffers (constrained to tare)"
    
  } else {
    
    # ==================================================
    # OPTION 2:
    # Group-level background sampling
    # sex × Logger_model
    # ==================================================
    
    group_data <- dloc %>%
      dplyr::group_by(sex, Logger_model) %>%
      dplyr::group_split()
    
    available_list <- purrr::map(
      group_data,
      function(dloc_group) {
        
        group_sex <- unique(dloc_group$sex)
        group_logger <- unique(dloc_group$Logger_model)
        
        # number of background points
        n_group <- nrow(dloc_group) * n_per_point
        
        # --------------------------------------------------
        # Exclusion buffers around used points
        # --------------------------------------------------
        
        known_buffer_group <- sf::st_buffer(
          dloc_group,
          dist = min_dist_m
        )
        
        known_buffer_group_union <- sf::st_union(
          known_buffer_group
        )
        
        # --------------------------------------------------
        # Available region:
        # study region minus exclusion buffers,
        # further constrained to fall within 'tare'
        # --------------------------------------------------
        
        available_region_group <- sf::st_difference(
          sf::st_union(study_region),
          known_buffer_group_union
        ) %>%
          sf::st_intersection(tare_union)
        
        # --------------------------------------------------
        # Random background sampling
        # --------------------------------------------------
        
        rnd <- sf::st_sample(
          available_region_group,
          size = n_group,
          type = "random",
          exact = TRUE
        )
        
        if (length(rnd) < n_group) {
          warning(
            paste(
              "Fewer available points sampled for:",
              group_sex,
              group_logger
            )
          )
        }
        
        # --------------------------------------------------
        # Copy metadata from observed points
        # --------------------------------------------------
        
        meta_source <- dloc_group %>%
          sf::st_drop_geometry()
        
        meta <- meta_source %>%
          dplyr::slice(
            rep(
              dplyr::row_number(),
              each = n_per_point
            )
          ) %>%
          dplyr::slice(
            seq_len(length(rnd))
          )
        
        # Randomise IDs
        meta$ID <- sample(
          meta_source$ID,
          size = nrow(meta),
          replace = TRUE
        )
        
        # Randomise timestamps
        meta$UTC_datetime <- sample(
          meta_source$UTC_datetime,
          size = nrow(meta),
          replace = TRUE
        )
        
        # --------------------------------------------------
        # Build sf object of available points
        # --------------------------------------------------
        
        available_group <- sf::st_sf(
          meta,
          geometry = rnd,
          crs = sf::st_crs(dloc)
        ) %>%
          dplyr::mutate(used = 0)
        
        # distance from colony/nest
        available_group$DistNest <- as.numeric(
          sf::st_distance(available_group, colony)
        )
        
        list(
          available = available_group,
          available_region = available_region_group,
          known_buffer = known_buffer_group
        )
      }
    )
    
    grouping_label <- "sex × Logger_model (constrained to tare)"
  }
  
  # --------------------------------------------------
  # Combine outputs
  # --------------------------------------------------
  
  available <- purrr::map_dfr(
    available_list,
    "available"
  )
  
  available_regions <- purrr::map(
    available_list,
    "available_region"
  )
  
  known_buffers <- purrr::map_dfr(
    available_list,
    "known_buffer"
  )
  
  used$DistNest <- as.numeric(
    sf::st_distance(used, colony)
  )
  
  combined_data <- dplyr::bind_rows(
    used,
    available
  )
  
  # --------------------------------------------------
  # Diagnostic plot
  # --------------------------------------------------
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = sf::st_union(study_region),
      fill = "grey95",
      color = "grey60"
    ) +
    ggplot2::geom_sf(
      data = tare_union,
      fill = NA,
      color = "darkgreen",
      linetype = "dashed"
    ) +
    ggplot2::geom_sf(
      data = known_buffers,
      fill = "red",
      color = NA,
      alpha = 0.2
    ) +
    ggplot2::geom_sf(
      data = used,
      color = "red",
      size = 2,
      alpha = 0.8
    ) +
    ggplot2::geom_sf(
      data = available,
      color = "blue",
      size = 1,
      alpha = 0.5
    ) +
    ggplot2::geom_sf(
      data = colony,
      color = "black",
      size = 4
    ) +
    ggplot2::facet_grid(
      sex ~ Logger_model
    ) +
    ggplot2::labs(
      title = "Used and available locations",
      subtitle = paste0(
        n_per_point,
        " available points per observed point; grouping = ",
        grouping_label
      )
    ) +
    ggplot2::theme_minimal()
  
  # --------------------------------------------------
  # Return outputs
  # --------------------------------------------------
  
  return(list(
    combined_data = combined_data,
    used = used,
    available = available,
    available_regions = available_regions,
    known_buffers = known_buffers,
    tare = tare_union,
    plot = p,
    group_by_id = group_by_id
  ))
}