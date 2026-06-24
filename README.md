# Habitat selection modelling pipeline for shags

## Purpose
- Prepare spatial and tracking data, fit habitat-selection models (GAM / ensemble), produce spatial predictions and diagnostic plots, and save results in multiple formats for publication.  
- Designed for a colony-centred study (example: Sklinna). Assumes projected spatial data (meter units, e.g. UTM).

---

## Repository overview
- `scripts/`  
  - `functions/` — helper functions (sourced by the master script)  
  - `pipeline/` — data preparation and pipeline steps  
- master script (entry point) — orchestrates data prep, model fitting, prediction, and saving outputs

---

## Main workflow steps
1. Load project environment and inputs (colony, kelp polygons, dive locations, harvest).  
2. Clip dive locations to kelp habitat around the colony and resolve intersection duplicates.  
3. Generate available (background) points (per-ID or grouped), optional time-matching and spatial thinning.  
4. Annotate used/available points with harvest overlap and timing covariates (multiple buffer radii).  
5. Finalise and combine year-specific modelling datasets.  
6. Build environmental rasters at 10 m (bathymetry, substrate, kelp mask).  
7. Fit habitat-selection models, predict across the study region, make plots.  
8. Save outputs (RDS/RData, GeoPackage, GeoTIFF, parquet, PNG/PDF, `sessionInfo`).

---

## Key helper functions
- `process_dloctare()` — clip dives to kelp, resolve duplicate intersections, produce diagnostic map.  
- `generate_available_points()` — sample background points; optional harvest stratification & weighting.  
- `add_harvest_overlap()` — find earliest overlapping harvest time per point and classify before/after/noharvest.  
- `finalise_dmodel()` — combine years, derive binary harvested flags and harvest-period labels.  
- `make_year_dmodel()` — orchestrates steps for a single year; returns processed objects.  
- `create_10m_raster()` — reproject, crop/mask and rasterise/resample vector or raster input to a common grid.  
- `create_harvest_dist_rasters()` — produce `harvested50m` / `harvested200m` rasters and distance-to-harvest/colony rasters.  
- `fit_predict_habitat_model()` — fit GAM(s), compute partial effects and spatial predictions.  
- `predict_maps_by_sex()`, `maps_by_sex_all_periods()` — produce sex- and period-specific prediction rasters and plots.  
- `save_habitat_results_pro()` — save models, tables, sf layers, rasters and plots in multiple formats.  
- `set_terra_threads()` — configure terra / OpenMP / GDAL / BLAS threading.

---

## Inputs expected (examples)
- Spatial: `colony` (sf), `tare` / `tareskog` (kelp polygons, sf), harvest polygons/events (sf)  
- Tracking: `dloc23`, `dloc24` (sf; columns: `ID`, `UTC_datetime`, `sex`, `Logger_model`, etc.)  
- Environmental: GEBCO NetCDF, substrate shapefiles

---

## Primary outputs
- Combined modelling dataset (`dmodel`) for 2023+2024 (saved as `dmodel.RData`)  
- Environmental raster stack (`env_stack_10m`) saved as parquet + SpatRaster  
- Harvest-distance raster results per year (GeoTIFF/parquet + metadata)  
- Model outputs: `.rds` / `.RData`, `model_summary.txt`, GeoTIFF prediction maps, GeoPackage vectors, parquet tables, PNG/PDF figures, `sessionInfo.txt`

---

## Environment variables and configuration
The scripts expect project variables available via `.Rprofile` or environment:
- `project_dir`
- `temp_dir`
- `data_dir`
- `logs_dir`
- `result_dir`
- `run_env`

Example `.Rprofile` snippet:
```r
project_dir <- "/path/to/project"
temp_dir    <- "/path/to/temp"
data_dir    <- file.path(project_dir, "data")
logs_dir    <- file.path(project_dir, "logs")
result_dir  <- file.path(project_dir, "results")
run_env     <- "personal"  # or "server"
```

## Dependencies

 - R >= 4.0 recommended
 - Main packages (install if missing): ``sf``, ``terra``, ``dplyr``, ``mgcv``, 
 ``ggplot2``, ``tidyr``, ``purrr``, ``arrow``, ``GGally``, ``dbscan``
 - Optional and recommended: ``RhpcBLASctl`` (control BLAS threads)
 - Use ``renv`` or similar to lock package versions.
 
 
## How to run

 - Prepare environment variables and ensure required input files are in expected locations.

 - Install required packages:
 
```r
install.packages(c("sf",
                    "terra",
                    "dplyr",
                    "mgcv",
                    "ggplot2",
                    "tidyr",
                    "purrr",
                    "arrow",
                    "GGally",
                    "dbscan"))
                    
# optional:
install.packages("RhpcBLASctl")
```
 - Run the master script:
    
    * Interactive (R / RStudio): source the master script.
    * Server (batch) with parameter-row index tpr:

```bash
Rscript scripts/fit_model_master.R 3
```
The script reads the command-line argument to select a parameter row for batch runs. Ensure ``run_env`` in ``.Rprofile`` is set appropriately.


## Outputs and files organization

Results are saved under:

``result_dir / version_label / <method> / <Year> / ...``

Typical files:

- ``prefix_full_results.rds`` / ``.RData``
- ``prefix_model.rds``, ``model_summary.txt``
- ``model_data.parquet`` (or CSV fallback)
- ``model_sf.gpkg`` and ``model_sf.parquet`` (WKT)
- ``pred_raster*.tif`` (GeoTIFFs)
- ``pred_raster_xyz.parquet``
- Plots (PNG/PDF)
- ``sessionInfo.txt``
