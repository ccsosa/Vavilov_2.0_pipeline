#' @title Getting soil variables using several cores
#'
#' @description
#' End-to-end preparation script for Species Distribution Modelling (SDM) of
#' D4R cacao-associated species. Loads environmental predictors (soil and
#' topography rasters, terrestrial ecoregions), reads and filters GBIF-derived
#' occurrence records against a curated species list, and extracts soil and
#' ecoregion/biome covariates at each occurrence point in parallel. The
#' resulting per-point covariate table is the input for the downstream
#' modelling stage (target-group background selection, spatial block
#' cross-validation via \pkg{blockCV}, environmental filtering, and Maxent
#' parameter tuning via \pkg{ENMeval}), which is implemented in later parts of
#' the pipeline.
#'
#' @details
#' This script is organized into three conceptual parts, mirrored in the
#' section documentation below:
#' \itemize{
#'   \item \strong{Part 0 -- Preparations}: package loading, user/base-directory
#'     switch, output directory paths, and loading of soil covariate rasters
#'     and terrestrial ecoregions (with an automatic geometry-repair and
#'     GeoPackage-caching step to avoid repeated shapefile validation).
#'   \item \strong{Part 1 -- Preparation of presence data and target group grid}:
#'     reading the curated species list, loading and merging two rounds of
#'     cleaned GBIF occurrence exports, and filtering occurrences down to the
#'     species set of interest, with the merged/filtered table cached to disk
#'     so the (expensive) merge step only runs once.
#'   \item \strong{Part 2 -- Distribution modelling (covariate extraction)}:
#'     parallel extraction of soil covariates and ecoregion/biome membership
#'     at each occurrence point, one species at a time, followed by
#'     row-binding into a single covariate table that is saved as an RDS file
#'     for use by the downstream Maxent modelling script.
#' }
#'
#' \strong{Parallelization note:} \code{parallel::mclapply()} is used with
#' \code{mc.cores = n_cores}, which relies on forking and is therefore
#' Unix/macOS-only (not supported on Windows). \code{terra} SpatRaster and
#' \code{sf} objects are not safely shared across forked/parallel workers, so
#' each worker re-reads the soil rasters and the ecoregions GeoPackage from
#' disk (\code{soil_paths}, \code{eco_path}) rather than relying on objects
#' captured from the parent process. The number of cores (\code{n_cores}) is
#' hardcoded rather than derived from \code{parallel::detectCores()}, in line
#' with the project convention of not saturating shared servers.
#'
#' @section Part 0 - Preparations:
#' Loads all required packages (\pkg{raster}, \pkg{dismo}, \pkg{rgeos},
#' \pkg{rgdal}, \pkg{blockCV}, \pkg{rJava}, \pkg{ENMeval}, \pkg{data.table},
#' \pkg{viridis}, \pkg{sf}, \pkg{terra}, \pkg{usdm}, \pkg{pROC},
#' \pkg{parallel}). Sets \code{basedir} according to the active \code{user}
#' switch (\code{"Chrystian"} or \code{"Tobias"}), and derives
#' \code{outputdir} (SDM results) and \code{spListDir} (buffer-richness
#' species-list summaries) from it.
#'
#' Loads nine 2.5-arcmin soil/topography rasters (bulk density, CEC, coarse
#' fragments, clay, nitrogen, pH, sand, silt, soil organic carbon) into a
#' single multi-layer \code{SpatRaster} (\code{soil_vars}).
#'
#' Loads the Ecoregions2017 layer: if a pre-validated GeoPackage
#' (\code{Ecoregions2017.gpkg}) already exists it is read directly; otherwise
#' the source shapefile is read with \code{sf::sf_use_s2(FALSE)}, repaired
#' with \code{sf::st_make_valid()}, and written out as a GeoPackage so future
#' runs skip the repair step.
#'
#' @section Part 1 - Preparation of presence data and target group grid:
#' Reads the curated species list (\code{Species_20250304_edible_part_curated.xlsx})
#' and derives the unique target species vector \code{sp_list_to}.
#'
#' Loads cleaned GBIF occurrence exports. If the merged/filtered output file
#' (\code{GBIF_joined_FINAL_20260317_spList_filtered.csv}) does not already
#' exist, two raw cleaned exports (\code{GBIF_joined_FINAL_3.csv} and
#' \code{GBIF_data_third_round_cleaned_FINAL_2.csv}) are read, restricted to a
#' common set of columns (species, searched name, country, year, adm1,
#' locality, collection/dataset/institution identifiers, source, coordinates),
#' row-bound into \code{uniq_file}, and filtered to only the species present
#' in \code{sp_list_to} (\code{p1}). Both the unfiltered merge and the
#' filtered subset are written to disk so subsequent runs can load the cached
#' filtered file (\code{p1}) directly instead of repeating the merge.
#'
#' @section Part 2 - Distribution modelling (covariate extraction):
#' For each species in \code{sp_list_to}, runs a parallel worker
#' (\code{parallel::mclapply}, \code{mc.cores = n_cores}) that:
#' \enumerate{
#'   \item Re-reads the soil raster stack and ecoregions GeoPackage from disk
#'     inside the worker (fork-safety for \code{terra}/\code{sf} objects).
#'   \item Subsets occurrence records to the current species.
#'   \item Extracts soil covariate values at each occurrence coordinate with
#'     \code{terra::extract()}.
#'   \item Converts occurrence points to an \code{sf} object and spatially
#'     joins them to the ecoregions layer (\code{sf::st_join()},
#'     \code{st_intersects}) to attach \code{ECO_NAME} and \code{BIOME_NAME}.
#'   \item Returns \code{NULL} for species with zero occurrence records after
#'     filtering (handled downstream by \code{do.call(rbind, ...)}, which
#'     silently drops \code{NULL} list elements).
#' }
#' The per-species results are combined into a single data frame
#' (\code{soil_p_i_2}) and saved as \code{input_data_soil_20260317.RDS} for
#' use in the subsequent Maxent modelling stage.
#'
#' @param user Character switch (\code{"Chrystian"} or \code{"Tobias"}) that
#'   selects the machine-specific \code{basedir} root path. Hardcoded at the
#'   top of the script rather than passed as a function argument.
#' @param n_cores Integer; number of parallel workers passed to
#'   \code{parallel::mclapply()}. Hardcoded (16) rather than derived from
#'   \code{parallel::detectCores()}.
#'
#' @return
#' No R object is returned to the console (this is a script, not a function).
#' Side effects written to disk:
#' \itemize{
#'   \item \code{Ecoregions2017.gpkg} -- validated ecoregions layer (created
#'     on first run if absent).
#'   \item \code{GBIF_joined_FINAL_20260317.csv} -- merged (unfiltered)
#'     occurrence export.
#'   \item \code{GBIF_joined_FINAL_20260317_spList_filtered.csv} -- occurrence
#'     export filtered to the curated species list.
#'   \item \code{input_data_soil_20260317.RDS} -- per-occurrence table of
#'     soil covariates, ecoregion, and biome, for all target species; the
#'     direct input to the downstream Maxent SDM script.
#' }
#'
#' @note
#' This documentation describes the script as written; no logic, control
#' flow, or parameter values have been altered.
NULL

# load packages
library(raster)
library(dismo)
library(rgeos)
library(rgdal)
library(blockCV)
library(rJava)
library(ENMeval)
# library(ecospat)
library(data.table)
library(viridis)
library(sf)
library(terra)
library(usdm)
library(pROC)
library(parallel)
################################################################################
# set user of the script
user <- "Chrystian"
# user <- "Tobias"

# set base directory
if (user == "Chrystian") {
  basedir <- "/catalogue/MultifLandscapesA1706"
}
if (user == "Tobias") {
  basedir <- "C:/Users/tobias/Dropbox/VAVILOV_2.0"
}

################################################################################
# set output directories
outputdir <- paste0(basedir, "/1.Data/Results/SDM results")
spListDir <-
  paste0(basedir, "/1.Data/Results/buffer_richness/summary")

################################################################################
### soil
soil_vars <-
  paste0(
    basedir,
    "/1.Data/RAW/Input_data/Soil and topography data/2_5min/",
    c(
      "bdod.tif",
      "cec.tif",
      "cfvo.tif",
      "clay.tif",
      "nitrogen.tif",
      "phh2o.tif",
      "sand.tif",
      "silt.tif",
      "soc.tif"
    )
  )
#loading all soil variables
soil_vars <- lapply(1:length(soil_vars), function(i) {
  x <- terra::rast(soil_vars[[i]])
  return(x)
})

soil_vars <- terra::rast(soil_vars)
################################################################################
#ecorregions loading
message("Loading ecorregions")

if (file.exists(paste0(
  basedir,
  "/1.Data/RAW/Input_data/Ecoregions/Ecoregions2017.gpkg"
))) {
  message("Loading ecorregions previously fixed!")
  eco <-
    sf::st_read(paste0(
      basedir,
      "/1.Data/RAW/Input_data/Ecoregions/Ecoregions2017.gpkg"
    ))
  
} else {
  sf_use_s2(FALSE)
  eco <-
    sf::st_read(paste0(
      basedir,
      "/1.Data/RAW/Input_data/Ecoregions/Ecoregions2017.shp"
    ))
  eco_fix <- st_make_valid(eco)
  sf::write_sf(
    eco_fix,
    paste0(
      basedir,
      "/1.Data/RAW/Input_data/Ecoregions/Ecoregions2017.gpkg"
    )
  )
  eco <- eco_fix
}
################################################################################
# load presence points

message("loading species list")

sp_list <- readxl::read_xlsx("/catalogue/MultifLandscapesA1706/1.Data/Results/species_lists/species_curated/Species_20250304_edible_part_curated.xlsx")
sp_list_to <- unique(sp_list$species)

message("Loading occurrences")

setwd(dir = paste0(basedir, "/1.Data/RAW/occurrences/cleaned"))

if(!file.exists(paste0(basedir, "/1.Data/RAW/occurrences/cleaned/GBIF_joined_FINAL_20260317_spList_filtered.csv"))){
  p1 <- read.csv("GBIF_joined_FINAL_3.csv")
  length(unique(p1$species_searched))
  p1_plus <- read.csv("GBIF_data_third_round_cleaned_FINAL_2.csv")
  length(unique(p1_plus$species_searched))
  
  p1 <- p1[,c("species","species_searched","country","year","adm1","locality","collection_code",
              "dataset_name","dataset_key","institution_code","source","lon","lat","coords")]
  
  p1_plus <- p1_plus[,c("species","species_searched","country","year","adm1","locality","collection_code",
                        "dataset_name","dataset_key","institution_code","source","lon","lat","coords")]
  uniq_file <- rbind(p1,p1_plus)
  
  p1 <- uniq_file[which(uniq_file$species_searched %in% sp_list_to),]
  write.csv(uniq_file,paste0("./", "GBIF_joined_FINAL_20260317.csv"))
  write.csv(p1,paste0("./", "GBIF_joined_FINAL_20260317_spList_filtered.csv"))
  length(unique(p1$species_searched))
} else {
  p1 <- read.csv(paste0(basedir, "/1.Data/RAW/occurrences/cleaned/GBIF_joined_FINAL_20260317_spList_filtered.csv"),header = T)
}

################################################################################
message("Getting ecorregions and soil information")

n_cores <- 16  # max cores; do not use detectCores() on shared servers

# terra SpatRaster and sf objects are not fork-safe: pass file paths and re-read inside workers
soil_paths <- paste0(
  basedir,
  "/1.Data/RAW/Input_data/Soil and topography data/2_5min/",
  c("bdod.tif","cec.tif","cfvo.tif","clay.tif","nitrogen.tif",
    "phh2o.tif","sand.tif","silt.tif","soc.tif")
)
eco_path <- paste0(basedir, "/1.Data/RAW/Input_data/Ecoregions/Ecoregions2017.gpkg")

soil_p_i <- parallel::mclapply(1:length(sp_list_to), function(i) {
  message(paste0("i: ", i, " / ", length(sp_list_to)))
  
  # re-read inside worker (terra/sf are not fork-safe)
  soil_vars_w <- terra::rast(soil_paths)
  eco_w       <- sf::st_read(eco_path, quiet = TRUE)
  
  p1_sp <- p1[which(p1$species_searched == sp_list_to[[i]]), ]
  
  if (nrow(p1_sp) > 0) {
    p2 <- terra::extract(soil_vars_w, cbind(p1_sp$lon, p1_sp$lat))
    p2 <- cbind(p1_sp[, c("species_searched", "lon", "lat")], p2)
    
    points <- sf::st_as_sf(p1_sp, coords = c("lon", "lat"), crs = sf::st_crs(eco_w))
    points_with_ecoregion <- sf::st_join(x = points, y = eco_w, join = sf::st_intersects)
    
    p2 <- cbind(p2,
                points_with_ecoregion$ECO_NAME,
                points_with_ecoregion$BIOME_NAME)
    colnames(p2)[c(13, 14)] <- c("Ecorregion", "Biome")
  } else {
    p2 <- NULL
  }
  return(p2)
}, mc.cores = n_cores)

################################################################################
#Joining all species in one file

soil_p_i_2 <- do.call(rbind, soil_p_i)

message("saving results (soil + ecorregions)")

saveRDS(
  soil_p_i_2,
  paste0(
    basedir,
    "/1.Data/Process/soil_extreme/",
    "input_data_soil_20260317.RDS"
  )
)
################################################################################