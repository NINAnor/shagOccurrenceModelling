#' Process dive locations relative to available kelp habitat around a colony
#'
#' This helper restricts kelp habitat to a buffer around a colony, finds dive
#' locations that fall within that kelp, resolves duplicate rows created by
#' spatial intersections, performs a diagnostic duplicate check, and returns a
#' small summary and a ggplot diagnostic map.
#'
#' Workflow:
#'  1. Add sampling year to dive locations
#'  2. Create a spatial buffer around the colony
#'  3. Clip kelp (tare) to the buffer
#'  4. Intersect dive locations with clipped kelp (points within kelp)
#'  5. Resolve duplicated rows created by st_intersection()
#'  6. Check whether duplicates remain
#'  7. Produce a diagnostic plot
#'
#' Inputs:
#'  - dloc     : sf POINTs of dive locations (expected to have a unique id column `seq`)
#'  - tare     : sf POLYGON(S) of kelp habitat
#'  - colony   : sf geometry (point or polygon) for colony location
#'  - year     : character or integer giving the sampling year (default "2024")
#'  - buffer_m : numeric buffer radius in meters (default 33000 -> 33 km)
#'
#' Returns: a list with:
#'  - SK33km            : colony buffer polygon (sf)
#'  - tare_33km         : kelp clipped to buffer (sf)
#'  - dloctare          : raw intersection output (may contain duplicated rows)
#'  - dloctare_resolved : duplicate-resolved dive data (sf)
#'  - duplicate_check   : unresolved duplicate seq+geom combinations (tibble)
#'  - plot              : ggplot diagnostic map
#'
#' Notes and assumptions:
#'  - CRS for inputs must be projected with meter units (st_buffer uses units of CRS)
#'  - The function expects `dloc` to contain a column named `seq` to identify dives.
#'    If `seq` is missing, a temporary unique id called `.seq_tmp` will be created.
#'  - Duplicate resolution strategy: for rows created by st_intersection() where
#'    a single dive point overlapped multiple polygons, attributes are merged by
#'    taking the first non-NA value per column, and geometry is reattached using
#'    a text representation of geometry (geom_id).
#'
#' @export
process_dloctare <- function(dloc,
                             tare,
                             colony,
                             year = "2024",
                             buffer_m = 33000) {
  # ---------------------------
  # Dependency check
  # ---------------------------
  pkgs <- c("sf", "dplyr", "ggplot2", "tidyr")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))

  `%>%` <- dplyr::`%>%`

  # ---------------------------
  # Helper: first_non_na
  # ---------------------------
  # Return the first non-NA value from vector x (or NA if none)
  first_non_na <- function(x) {
    y <- x[!is.na(x)]
    if (length(y) == 0) NA else y[1]
  }

  # ---------------------------
  # Basic input checks & preparation
  # ---------------------------
  if (!inherits(dloc, "sf")) stop("dloc must be an sf object")
  if (!inherits(tare, "sf")) stop("tare must be an sf object")
  if (!inherits(colony, "sf")) stop("colony must be an sf object")

  # Ensure there is an identifier column 'seq'. If not, create a temporary one.
  created_seq_tmp <- FALSE
  if (!"seq" %in% names(dloc)) {
    dloc$.seq_tmp <- seq_len(nrow(dloc))
    dloc <- dloc %>% dplyr::rename(seq = .seq_tmp)
    created_seq_tmp <- TRUE
    warning("Input `dloc` had no 'seq' column: a temporary 'seq' was created.")
  }

  # ---------------------------
  # 1) Add sampling year to dive locations
  # ---------------------------
  dloc$Year <- year

  # ---------------------------
  # 2) Create buffer around colony
  #    - Assumes CRS units are meters
  # ---------------------------
  SK33km <- sf::st_buffer(colony, dist = buffer_m)

  # ---------------------------
  # 3) Clip kelp habitat to buffer (available habitat)
  # ---------------------------
  tare_33km <- sf::st_intersection(tare, SK33km)

  # ---------------------------
  # 4) Intersect dive points with clipped kelp
  #    - st_intersection returns only points that overlap polygons.
  #    - A point may intersect multiple polygons, creating duplicate rows.
  # ---------------------------
  dloctare <- sf::st_intersection(dloc, tare_33km)

  # If no intersections, return early with an informative plot and empty resolved object
  if (nrow(dloctare) == 0) {
    message("No dive locations intersect kelp within buffer.")
    p_empty <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = tare_33km, fill = "grey85", color = NA) +
      ggplot2::geom_sf(data = SK33km, fill = NA, color = "red", linewidth = 1) +
      ggplot2::labs(
        title = paste0("No dive locations in kelp habitat, ", year),
        subtitle = paste0("Kelp clipped to ", buffer_m / 1000, " km around colony")
      ) +
      ggplot2::theme_minimal()
    return(list(
      SK33km = SK33km,
      tare_33km = tare_33km,
      dloctare = dloctare,
      dloctare_resolved = dloctare[0, , drop = FALSE],
      duplicate_check = tibble::tibble(),
      plot = p_empty
    ))
  }

  # ---------------------------
  # 5) Resolve duplicated spatial matches
  #    Strategy:
  #      - Represent geometry as WKT (geom_id) to detect identical spatial rows
  #      - Group by seq + geom_id and summarise across columns keeping first non-NA
  #      - Reattach geometry from a distinct seq+geom_id selection
  # ---------------------------
  # Create geom_id and perform tabular collapse
  dloctare_resolved <- dloctare %>%
    dplyr::mutate(geom_id = sf::st_as_text(sf::st_geometry(.))) %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(seq, geom_id) %>%
    dplyr::summarise(
      dplyr::across(dplyr::everything(), first_non_na),
      .groups = "drop"
    ) %>%
    # Reattach a single geometry per seq + geom_id
    dplyr::left_join(
      dloctare %>%
        dplyr::mutate(geom_id = sf::st_as_text(sf::st_geometry(.))) %>%
        dplyr::select(seq, geom_id, geometry) %>%
        dplyr::distinct(seq, geom_id, .keep_all = TRUE),
      by = c("seq", "geom_id")
    ) %>%
    sf::st_as_sf()

  # ---------------------------
  # 6) Diagnostic duplicate check
  #    - Count any remaining seq + geom_id combinations with >1 row
  # ---------------------------
  duplicate_check <- dloctare_resolved %>%
    sf::st_drop_geometry() %>%
    dplyr::count(seq, geom_id) %>%
    dplyr::filter(n > 1)

  if (nrow(duplicate_check) > 0) {
    warning("Some seq + geometry combinations remain duplicated after resolution. See duplicate_check.")
  }

  # ---------------------------
  # 7) Diagnostic map
  #    - Grey polygons : clipped kelp area
  #    - Red outline   : colony buffer
  #    - Colored points: dive locations (resolved), colored by 'naturtype' if present
  # ---------------------------
  aes_args <- if ("naturtype" %in% names(dloctare_resolved)) {
    aes(color = naturtype)
  } else {
    NULL
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = tare_33km, fill = "grey85", color = NA) +
    ggplot2::geom_sf(data = SK33km, fill = NA, color = "red", linewidth = 1) +
    ggplot2::geom_sf(data = dloctare_resolved, mapping = aes_args, size = 2, alpha = 0.8) +
    ggplot2::labs(
      title = paste0("Dive locations in kelp habitat, ", year),
      subtitle = paste0("Kelp clipped to ", buffer_m / 1000, " km around colony"),
      color = if ("naturtype" %in% names(dloctare_resolved)) "Naturtype" else NULL
    ) +
    ggplot2::theme_minimal()

  # ---------------------------
  # If we created a temporary seq column, revert the name in outputs (keep original data intact)
  # ---------------------------
  if (created_seq_tmp) {
    # rename seq back to .seq_tmp in returned objects to avoid surprising side-effects
    names_to_rename <- c("dloctare", "dloctare_resolved")
    for (nm in names_to_rename) {
      if (exists(nm, inherits = FALSE)) {
        obj <- get(nm)
        if ("seq" %in% names(obj)) names(obj)[names(obj) == "seq"] <- ".seq_tmp"
        assign(nm, obj)
      }
    }
  }

  # ---------------------------
  # Return a list of results
  # ---------------------------
  return(list(
    SK33km = SK33km,
    tare_33km = tare_33km,
    dloctare = dloctare,
    dloctare_resolved = dloctare_resolved,
    duplicate_check = duplicate_check,
    plot = p
  ))
}
