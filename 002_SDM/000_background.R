#' @title Generate Random Background Points for Species Distribution Modeling
#'
#' @description
#' This script generates a large set of random background (pseudo-absence) points
#' across Latin America using a WorldClim 30-arc-second bioclimatic raster as a
#' spatial template. Background points are sampled proportionally to raster cell
#' availability and exported as a CSV file for use in MaxEnt or similar SDM workflows.
#'
#' @details
#' The workflow proceeds as follows:
#' \enumerate{
#'   \item A WorldClim 30s bioclimatic raster (BIO1) is loaded as a spatial template
#'         to define the sampling extent and resolution.
#'   \item Occurrence records are loaded from a pre-cleaned GBIF joined dataset
#'         (not directly used in background sampling here, but loaded for context).
#'   \item \code{dismo::randomPoints()} samples \code{n = 200,000} random geographic
#'         coordinates within non-NA cells of the template raster, with up to 3x
#'         oversampling attempts (\code{tryf = 3}) to reach the target count.
#'   \item Output columns are renamed to \code{lon} and \code{lat}.
#'   \item The resulting background point table is written to
#'         \code{<basedir>/Input_data/background/LATAM_BG.csv}.
#' }
#'
#' A commented-out alternative using \code{dismo::gridSample()} is retained for
#' spatially thinned background sampling if regular grid subsampling is preferred.
#'
#' @section Inputs:
#' \describe{
#'   \item{\code{wc2.1_30s_bio_1.tif}}{WorldClim v2.1 BIO1 raster at 30-arc-second
#'         resolution, used as the spatial template for background sampling.}
#'   \item{\code{GBIF_joined_FINAL_3.csv}}{Cleaned occurrence records joined from
#'         GBIF, loaded via \code{data.table::fread()}.}
#' }
#'
#' @section Outputs:
#' \describe{
#'   \item{\code{LATAM_BG.csv}}{A CSV file with \code{lon} and \code{lat} columns
#'         containing 200,000 random background points. Written to
#'         \code{<basedir>/Input_data/background/}.}
#' }
#'
#' @section Dependencies:
#' \itemize{
#'   \item \pkg{geodata} — spatial data retrieval utilities
#'   \item \pkg{sf} — simple features for vector operations
#'   \item \pkg{terra} — raster and vector spatial data handling
#'   \item \pkg{dismo} — SDM utilities including \code{randomPoints()} and \code{gridSample()}
#'   \item \pkg{data.table} — fast CSV reading via \code{fread()}
#' }
#'
#' @note
#' \itemize{
#'   \item The template raster is loaded with \code{raster::raster()} (i.e., the legacy
#'         \pkg{raster} package) for compatibility with \code{dismo::randomPoints()},
#'         which does not accept \pkg{terra} \code{SpatRaster} objects.
#'   \item \code{tryf = 3} instructs \code{randomPoints()} to attempt sampling up to
#'         \code{3 * n} candidate points internally to compensate for NA-masked cells.
#'   \item The commented-out \code{gridSample()} call can be reactivated for
#'         spatially regularised subsampling of background points on a grid.
#'   \item \code{basedir} must be set correctly before running; all input/output paths
#'         are constructed relative to it.
#' }
#'
#' @seealso
#' \code{\link[dismo]{randomPoints}}, \code{\link[dismo]{gridSample}},
#' \code{\link[raster]{raster}}, \code{\link[data.table]{fread}}
#'
#' @author Chrystian C. Sosa

library(geodata)
library(sf)
library(terra)
library(dismo)
library(data.table)

basedir <- "E:/CSOSA/Dropbox/VAVILOV_2.0/DATA"

#template
template <- raster::raster(paste0(basedir,"/Input_data/climate_data/30s/present/wc2.1_30s_bio_1.tif"))

#occurrences
occ <- fread("E:/CSOSA/Dropbox/VAVILOV_2.0/DATA/DATABASES/occurrences/cleaned/GBIF_joined_FINAL_3.csv")

rp <- dismo::randomPoints(template, n = 200000, tryf = 3)
colnames(rp) <- c("lon","lat")

write.csv(rp,paste0(basedir,"/","Input_data/background","/","LATAM_BG.csv"),na = "",row.names = F)

# background_TG <- data.frame(gridSample(background_TG[,c("lon", "lat")], r = ref))