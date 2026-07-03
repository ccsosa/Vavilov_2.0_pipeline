#' Download, Process, and Harmonize Climate and Soil Predictors for SDM
#'
#' @description
#' Prepares present and future climate layers (WorldClim 2.5 arcmin) plus soil/
#' topography covariates for use in species distribution models covering Latin
#' America. The workflow:
#' \enumerate{
#'   \item Downloads/loads present WorldClim bioclimatic variables, crops and
#'         masks them to a Latin-American extent, and writes one GeoTIFF per layer.
#'   \item Iterates over a set of CMIP6 GCMs x SSP scenarios (2041-2060),
#'         applies the same crop/mask, and writes layers into a
#'         \code{future/<GCM>/<SSP>/} directory tree.
#'   \item Resamples a stack of soil and topography rasters to match the
#'         2.5 arcmin climate grid and writes the harmonized outputs.
#' }
#'
#' @section Inputs:
#' \describe{
#'   \item{\strong{GCMs}}{Character vector of CMIP6 model identifiers:
#'     \code{"ACCESS-CM2"}, \code{"BCC-CSM2-MR"}, \code{"GISS-E2-1-G"},
#'     \code{"INM-CM5-0"}, \code{"MIROC6"}, \code{"MPI-ESM1-2-HR"}.}
#'   \item{\strong{ssp}}{Character vector of emission scenarios:
#'     \code{"ssp245"}, \code{"ssp370"}.}
#'   \item{\strong{extent}}{A \code{SpatVector} (ESRI Shapefile) of Latin-American
#'     country boundaries used for cropping and masking.
#'     Path: \code{<basedir>/Input_data/adm0/adm0_Latam.shp}.}
#'   \item{\strong{soil_vars}}{Character vector of soil/topography variable names
#'     (\code{bdod}, \code{cec}, \code{cfvo}, \code{clay}, \code{nitrogen},
#'     \code{phh2o}, \code{SAGA_TWI}, \code{sand}, \code{silt}, \code{soc}).
#'     Source GeoTIFFs expected under
#'     \code{<basedir>/Input_data/Soil and topography data/<var>.tif}.}
#' }
#'
#' @section Directory structure created:
#' \preformatted{
#' <dir>/2_5min/
#'   present/                      # present bioclim layers (bio_1 to bio_19)
#'   future/
#'     <GCM>/
#'       <SSP>/                    # future bioclim layers per GCM x SSP
#' <basedir>/Input_data/Soil and topography data/2_5min/
#'                                 # resampled soil/topography layers
#' }
#'
#' @section Parameters (inline constants):
#' \describe{
#'   \item{\code{var}}{WorldClim variable type: \code{"bio"} (19 bioclimatic variables).}
#'   \item{\code{res}}{Spatial resolution: \code{2.5} arcmin.}
#'   \item{\code{time}}{CMIP6 time slice: \code{"2041-2060"}.}
#'   \item{\code{download}}{Set to \code{FALSE}; data must already be cached locally at \code{path}.}
#' }
#'
#' @section Output files:
#' Each bioclimatic variable is written as a single-band GeoTIFF named after
#' the original layer (e.g., \code{wc2.1_2.5m_bio_1.tif}). Soil/topography
#' outputs share the same base name as the source file (e.g., \code{clay.tif})
#' and are resampled using bilinear interpolation to match the present climate
#' template raster (\code{clim[[1]]}).
#'
#' @note
#' \itemize{
#'   \item \code{terra} objects (\code{SpatRaster}, \code{SpatVector}) are
#'     \strong{not} safe to export across parallel workers. If this script is
#'     parallelized, reload rasters from file paths inside each worker rather
#'     than exporting \code{clim} or \code{extent}.
#'   \item All \code{writeRaster()} calls use \code{overwrite = TRUE};
#'     existing files will be silently replaced.
#'   \item The commented-out DEM resampling block (1 km to 5 km) is retained
#'     for reference but is not executed.
#' }
#'
#' @seealso
#' \code{\link[geodata]{worldclim_global}},
#' \code{\link[geodata]{cmip6_world}},
#' \code{\link[terra]{crop}},
#' \code{\link[terra]{mask}},
#' \code{\link[terra]{writeRaster}},
#' \code{\link[terra]{resample}}
#'
#' @author Chrystian C. Sosa

# -- Libraries -----------------------------------------------------------------
library(geodata)
library(sf)
library(terra)

# -- Configuration -------------------------------------------------------------

#' CMIP6 General Circulation Models to process
GCMs <- c(
  "ACCESS-CM2",
  "BCC-CSM2-MR",
  "GISS-E2-1-G",
  "INM-CM5-0",
  "MIROC6",
  "MPI-ESM1-2-HR"
)

#' Shared Socioeconomic Pathway emission scenarios
ssp <- c("ssp245", "ssp370")

#' Latin-American administrative boundary used for crop/mask operations.
#' Loaded as a SpatVector; do NOT export this object to parallel workers —
#' pass the file path and reload inside each worker instead.
extent <- terra::vect("D:/PROGRAMAS/Dropbox/VAVILOV_2.0/DATA/Input_data/adm0/adm0_Latam.shp")

#' Root directory for 2.5 arcmin climate data
dir <- "D:/PROGRAMAS/Dropbox/VAVILOV_2.0/DATA/Input_data/climate_data/2_5min"

#' Project base directory (used for soil/topography paths)
basedir <- "D:/PROGRAMAS/Dropbox/VAVILOV_2.0/DATA"

#' Output directory for present-climate layers; created if absent
pre_dir <- paste0(dir, "/", "present")
if (!dir.exists(pre_dir)) dir.create(pre_dir)

#' Output directory for future-climate layers; created if absent
fut_dir <- paste0(dir, "/future")
if (!dir.exists(fut_dir)) dir.create(fut_dir)

# ==============================================================================
# 1. PRESENT CLIMATE
# ==============================================================================

#' Load present WorldClim bioclimatic variables (19 layers, 2.5 arcmin).
#' \code{download = FALSE} assumes tiles are already cached at \code{path}.
clim <- geodata::worldclim_global(var = "bio", res = 2.5, download = FALSE, path = dir)

#' Crop to the Latin-American bounding box, then mask to country boundaries
clim <- terra::crop(clim, extent)
clim <- terra::mask(clim, extent)

#' Preserve layer names from the present stack; reused to name future layers
#' so that present and future GeoTIFFs are consistently labelled across all
#' GCM x SSP combinations.
nlyr_name <- names(clim)

#' Write each bioclimatic variable as an individual single-band GeoTIFF.
#' Layer index \code{i} is used instead of subsetting by name to avoid any
#' mismatch between the SpatRaster internal order and the names vector.
for (i in 1:nlyr(clim)) {
  terra::writeRaster(
    clim[[i]],
    filename  = paste0(pre_dir, "/", names(clim)[[i]], ".tif"),
    overwrite = TRUE
  )
}
rm(i)

# ==============================================================================
# 2. FUTURE CLIMATE  (CMIP6, 2041-2060)
# ==============================================================================

#' Outer loop: iterate over GCMs
for (i in 1:length(GCMs)) {
  
  message("Processing GCM ", i, "/", length(GCMs), ": ", GCMs[[i]])
  
  #' Create a subdirectory for the current GCM if it does not exist
  GCM_dir <- paste0(fut_dir, "/", GCMs[[i]])
  if (!dir.exists(GCM_dir)) dir.create(GCM_dir)
  
  #' Inner loop: iterate over SSP scenarios
  for (j in 1:length(ssp)) {
    
    message("  SSP ", j, "/", length(ssp), ": ", ssp[[j]])
    
    #' Create a subdirectory for the current SSP scenario if it does not exist
    ssp_dir <- paste0(GCM_dir, "/", ssp[[j]])
    if (!dir.exists(ssp_dir)) dir.create(ssp_dir)
    
    #' Load CMIP6 bioclimatic variables for the current GCM x SSP combination.
    #' The SSP string (e.g. "ssp245") is stripped of its prefix and coerced to
    #' integer (245) as required by \code{geodata::cmip6_world()}.
    clim_fut <- geodata::cmip6_world(
      model    = GCMs[[i]],
      ssp      = as.numeric(sub(pattern = "ssp", replacement = "", x = ssp[[j]])),
      time     = "2041-2060",
      var      = "bioc",
      download = FALSE,
      res      = 2.5,
      path     = dir
    )
    
    #' Crop and mask future layers to the same Latin-American extent used for
    #' present climate, ensuring spatial consistency across all stacks.
    clim_fut <- terra::crop(clim_fut, extent)
    clim_fut <- terra::mask(clim_fut, extent)
    
    #' Write future layers using the layer names inherited from the present
    #' stack (\code{nlyr_name}) so filenames are identical across time periods,
    #' GCMs, and SSPs — simplifying downstream file-path construction.
    for (k in 1:length(nlyr_name)) {
      terra::writeRaster(
        clim_fut[[k]],
        filename  = paste0(ssp_dir, "/", nlyr_name[[k]], ".tif"),
        overwrite = TRUE
      )
    }
  }
}

# ==============================================================================
# 3. SOIL AND TOPOGRAPHY PREDICTORS
# ==============================================================================

# NOTE: DEM resampling (1 km -> 5 km) is retained below for reference but
# is currently commented out.
#
# alt <- paste0(basedir, "/Input_data/Other spatial data/merit_DEM_1km_modelling_extent.tif")
# alt <- terra::resample(rast(alt), clim[[1]], method = "near", threads = TRUE)
# writeRaster(alt, paste0(basedir, "/Input_data/Other spatial data/merit_DEM_5km_modelling_extent.tif"))

#' Soil and topography variables to harmonize.
#' Source files expected at:
#'   \code{<basedir>/Input_data/Soil and topography data/<var>.tif}
#' Output files written to:
#'   \code{<basedir>/Input_data/Soil and topography data/2_5min/<var>.tif}
soil_vars <- c("bdod", "cec", "cfvo", "clay", "nitrogen",
               "phh2o", "SAGA_TWI", "sand", "silt", "soc")

#' Resample each soil/topography variable to match the 2.5 arcmin climate grid.
#' \code{clim[[1]]} serves as the spatial template (extent, resolution, CRS).
#' Default resampling method is bilinear, appropriate for continuous variables.
for (i in 1:length(soil_vars)) {
  
  message("Resampling soil variable ", i, "/", length(soil_vars), ": ", soil_vars[[i]])
  
  x <- terra::rast(
    paste0(basedir, "/Input_data/Soil and topography data/", soil_vars[[i]], ".tif")
  )
  
  #' Resample to climate grid; \code{threads = TRUE} enables multi-threaded
  #' processing within terra (does not require a parallel backend).
  x <- terra::resample(x, clim[[1]], threads = TRUE)
  
  terra::writeRaster(
    x,
    filename  = paste0(basedir, "/Input_data/Soil and topography data/2_5min/", soil_vars[[i]], ".tif"),
    overwrite = TRUE
  )
}