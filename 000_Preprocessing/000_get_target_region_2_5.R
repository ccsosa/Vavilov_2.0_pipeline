#' Create a Binary Target-Region Raster at 2.5 arcmin Resolution
#'
#' @description
#' Builds a binary template raster defining the modelling target region for the
#' MultifLandscapes (A1706) pipeline. The workflow:
#' \enumerate{
#'   \item Loads a species concave-hull raster as the spatial reference extent.
#'   \item Crops a WorldClim 2.5 arcmin bioclimatic layer to that extent.
#'   \item Recodes all non-NA cells to 1, producing a binary presence mask.
#'   \item Writes the result as a GeoTIFF to be used as the spatial template
#'         in downstream SDM steps.
#' }
#'
#' @section Inputs:
#' \describe{
#'   \item{\code{adm0}}{Latin-American country boundaries (\code{SpatVector}).
#'     Loaded for spatial context; not directly used in the masking steps shown
#'     but available for optional crop/mask operations.
#'     Path: \code{<dir_2_5>/adm0/adm0_Latam.shp}.}
#'   \item{\code{x}}{Concave-hull raster for a reference species
#'     (\emph{Solanum tuberosum}), used solely to define the target extent.
#'     Path: \code{<dir>/Solanum tuberosum_conc.tif}.}
#'   \item{\code{template}}{WorldClim bioclimatic layer 1 (BIO1, mean annual
#'     temperature) at 2.5 arcmin, used as the resolution/CRS template.
#'     Path: \code{<dir_2_5>/climate_data/2_5min/climate/wc2.1_2.5m/wc2.1_2.5m_bio_1.tif}.}
#' }
#'
#' @section Output:
#' A single-band GeoTIFF (\code{target_region.tif}) where:
#' \itemize{
#'   \item Cells within the target region = \code{1}.
#'   \item Cells outside (originally \code{NA}) remain \code{NA}.
#' }
#' Written to: \code{<dir_2_5>/Raster target region/2_5min/target_region.tif}.
#'
#' @note
#' \itemize{
#'   \item \code{library(raster)} is loaded for compatibility but all spatial
#'     operations use \code{terra}. Consider removing \code{raster} if no
#'     legacy functions are needed downstream.
#'   \item Cell indexing via \code{which(!is.na(template[]))} operates on the
#'     full raster value matrix in memory; for very large extents consider
#'     \code{terra::classify()} as a memory-safer alternative:
#'     \code{classify(template, cbind(NA, NA), others = 1)}.
#'   \item \code{terra} objects (\code{SpatRaster}, \code{SpatVector}) must
#'     \strong{not} be exported to parallel workers — pass file paths and
#'     reload inside each worker.
#' }
#'
#' @seealso
#' \code{\link[terra]{rast}},
#' \code{\link[terra]{crop}},
#' \code{\link[terra]{writeRaster}},
#' \code{\link[terra]{classify}}

# -- Libraries -----------------------------------------------------------------
library(raster) # retained for legacy compatibility; spatial ops use terra
library(terra)

# -- Paths ---------------------------------------------------------------------

#' Root directory for concave-hull species rasters
dir <- "//catalogue/MultifLandscapesA1706/1.Data/Results/concave_hull_rasters_ADM0"

#' Root directory for raw input data (climate, boundaries, templates)
dir_2_5 <- "//catalogue/MultifLandscapesA1706/1.Data/RAW/Input_data"

# -- Load spatial inputs -------------------------------------------------------

#' Latin-American country boundaries.
#' Available for downstream crop/mask steps; not consumed in this script.
#' Do NOT export to parallel workers — pass the path and reload inside workers.
adm0 <- terra::vect(paste0(dir_2_5, "/adm0/", "adm0_Latam.shp"))

#' Concave-hull raster for Solanum tuberosum.
#' Used exclusively to define the spatial extent of the target region.
x <- terra::rast(paste0(dir, "/", "Solanum tuberosum_conc.tif"))

#' WorldClim BIO1 layer at 2.5 arcmin resolution.
#' Provides the CRS, resolution, and cell alignment for the output template.
template <- terra::rast(
  paste0(dir_2_5, "/climate_data/2_5min/climate/wc2.1_2.5m/wc2.1_2.5m_bio_1.tif")
)

# -- Build binary target-region mask -------------------------------------------

#' Crop BIO1 to the extent of the species concave-hull raster.
#' After cropping, non-NA cells define the modelling target region.
template <- terra::crop(template, x)

#' Recode all non-NA cells to 1, producing a binary presence mask.
#' NA cells (outside the hull extent) are left unchanged.
#'
#' @note For large rasters, the memory-safer alternative is:
#'   \code{template <- terra::classify(template, cbind(NA, NA), others = 1)}
template[which(!is.na(template[]))] <- 1

# -- Write output --------------------------------------------------------------

#' Write the binary target-region raster.
#' This file serves as the spatial template (extent, resolution, CRS) for all
#' subsequent SDM raster operations in the MultifLandscapes pipeline.
terra::writeRaster(
  template,
  filename  = paste0(dir_2_5, "/Raster target region/2_5min/target_region.tif"),
  overwrite = TRUE
)