# ============================================================
# Process dive locations relative to available kelp habitat
# around a colony location
#
# Workflow:
# 1. Add sampling year to dive locations
# 2. Create a spatial buffer around the colony
# 3. Restrict kelp habitat to within the buffer
# 4. Restrict dive locations to kelp habitat
# 5. Resolve duplicated rows created by st_intersection()
# 6. Check whether duplicates remain
# 7. Produce a diagnostic plot
#
# Inputs:
#   dloc      = sf object containing dive locations (points)
#   tare      = sf object containing kelp polygons
#   colony    = sf object representing colony location
#   year      = character/integer indicating sampling year
#   buffer_m  = buffer radius in meters
#
# Returns:
#   A list containing:
#     - SK33km              : colony buffer polygon
#     - tare_33km           : kelp clipped to buffer
#     - dloctare            : raw intersection output
#     - dloctare_resolved   : duplicate-resolved dive data
#     - duplicate_check     : unresolved duplicate rows
#     - plot                : ggplot map
# ============================================================

process_dloctare <- function(dloc,
                             tare,
                             colony,
                             year = "2024",
                             buffer_m = 33000) {
  
  # ----------------------------------------------------------
  # Helper function:
  # Extract the first non-NA value from a vector
  #
  # Used to merge duplicated rows created during
  # spatial intersection where attributes are split
  # across multiple rows.
  # ----------------------------------------------------------
  
  first_non_na <- function(x) {
    y <- x[!is.na(x)]
    if (length(y) == 0) NA else y[1]
  }
  
  
  # ----------------------------------------------------------
  # Add sampling year to dive location data
  # ----------------------------------------------------------
  
  dloc$Year <- year
  
  
  # ----------------------------------------------------------
  # Create a buffer around the colony
  #
  # Assumes CRS units are meters.
  # Example:
  #   buffer_m = 33000  -> 33 km radius
  # ----------------------------------------------------------
  
  SK33km <- sf::st_buffer(colony, dist = buffer_m)
  
  
  # ----------------------------------------------------------
  # Restrict kelp habitat to within the colony buffer
  #
  # This defines the available habitat area accessible
  # to the birds.
  # ----------------------------------------------------------
  
  tare_33km <- sf::st_intersection(tare, SK33km)
  
  
  # ----------------------------------------------------------
  # Restrict dive locations to available kelp habitat
  #
  # st_intersection() returns only dive locations that
  # overlap kelp polygons.
  #
  # Important:
  # A single dive point may intersect multiple polygons,
  # creating duplicated rows.
  # ----------------------------------------------------------
  
  dloctare <- sf::st_intersection(dloc, tare_33km)
  
  
  # ----------------------------------------------------------
  # Resolve duplicated spatial matches
  #
  # Strategy:
  # 1. Convert geometry to text (geom_id)
  # 2. Group by dive ID + geometry
  # 3. Merge duplicated rows by keeping the first
  #    non-missing value for each column
  #
  # This preserves:
  #   - multiple dives at same location
  #   - one row per actual dive observation
  # ----------------------------------------------------------
  
  dloctare_resolved <- dloctare %>%
    
    # Create a text representation of geometry
    dplyr::mutate(
      geom_id = sf::st_as_text(sf::st_geometry(.))
    ) %>%
    
    # Temporarily remove geometry for tabular operations
    sf::st_drop_geometry() %>%
    
    # Group duplicated spatial matches
    dplyr::group_by(seq, geom_id) %>%
    
    # Merge duplicated rows
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        first_non_na
      ),
      .groups = "drop"
    ) %>%
    
    # Reattach geometry column
    dplyr::left_join(
      
      dloctare %>%
        dplyr::mutate(
          geom_id = sf::st_as_text(sf::st_geometry(.))
        ) %>%
        dplyr::select(seq, geom_id, geometry) %>%
        dplyr::distinct(
          seq,
          geom_id,
          .keep_all = TRUE
        ),
      
      by = c("seq", "geom_id")
    ) %>%
    
    # Convert back to sf object
    sf::st_as_sf()
  
  
  # ----------------------------------------------------------
  # Diagnostic check:
  # Identify any duplicated seq + geometry combinations
  # remaining after cleaning
  # ----------------------------------------------------------
  
  duplicate_check <- dloctare_resolved %>%
    sf::st_drop_geometry() %>%
    dplyr::count(seq, geom_id) %>%
    dplyr::filter(n > 1)
  
  
  # ----------------------------------------------------------
  # Create diagnostic map
  #
  # Grey polygons = available kelp habitat
  # Red outline    = colony buffer
  # Colored points = dive locations by naturtype
  # ----------------------------------------------------------
  
  p <- ggplot() +
    
    # Kelp habitat polygons
    geom_sf(
      data = tare_33km,
      fill = "grey85",
      color = NA
    ) +
    
    # Colony buffer
    geom_sf(
      data = SK33km,
      fill = NA,
      color = "red",
      linewidth = 1
    ) +
    
    # Dive locations
    geom_sf(
      data = dloctare_resolved,
      aes(color = naturtype),
      size = 2,
      alpha = 0.8
    ) +
    
    # Labels
    labs(
      title = paste0(
        "Dive locations in kelp habitat, ",
        year
      ),
      subtitle = paste0(
        "Kelp clipped to ",
        buffer_m / 1000,
        " km around colony"
      ),
      color = "Naturtype"
    ) +
    
    # Theme
    theme_minimal()
  
  
  # ----------------------------------------------------------
  # Return outputs
  # ----------------------------------------------------------
  
  return(list(
    SK33km = SK33km,
    tare_33km = tare_33km,
    dloctare = dloctare,
    dloctare_resolved = dloctare_resolved,
    duplicate_check = duplicate_check,
    plot = p
  ))
}