### This script imports GBIF occurrence records for a set of species
# author: Tobias Fremout (tobias.fremout@gmail.com)
#this is a parallelized version to download thousand of records at the same time!

# load packages
library(dismo)
library(terra)
library(parallel)
library(pbapply)
library(readxl)

# set base directory
# basedir <- "D:/PROGRAMAS/Dropbox/VAVILOV_2.0"
basedir <- "E:/CSOSA/Dropbox/VAVILOV_2.0"

# ==============================================================================
# 1. LOAD SPECIES SETS
# ==============================================================================

#' First-round species list (already downloaded in a previous session).
#' Used to identify which species still need GBIF records in this second round.
species_set_first_round <- readxl::read_xlsx(paste0(basedir,"/1.Processing/raw_Species_lists/Species_list_LATAM.xlsx")) #first round

#' Full target species list for the second round of downloads
species_set <- readxl::read_xlsx(paste0(basedir,"/1.Processing/raw_Species_lists/Vavilov_native_species_to_curate.xlsx"))

# get only the species names
species_set <- species_set$species

#' Subset to species not yet downloaded in the first round
sp_missing <- species_set[!species_set %in% species_set_first_round$species]

species_set<- sp_missing

# ==============================================================================
# 2. PARSE GENUS AND SPECIFIC EPITHET
# ==============================================================================

#' Split each binomial into genus and specific epithet columns.
#' Required by \code{dismo::gbif()}, which accepts them as separate arguments.
species_set <- data.frame(species = species_set,
                          genus = NA,
                          specEpithet = NA)
for (i in 1:nrow(species_set)) {
  species_split <- strsplit(species_set$species[i], " ")[[1]]
  species_set$genus[i] <- species_split[1]
  species_set$specEpithet[i] <- species_split[2]
};rm(i)

# ==============================================================================
# 3. DEFINE DOWNLOAD EXTENT FROM WORLDCLIM TEMPLATE
# ==============================================================================

#' Load WorldClim BIO1 at 30 s resolution solely to extract the Latin-American
#' bounding box. The raster object is discarded after extent extraction.
t <- raster::raster(paste0(basedir, "/DATA/Latin America_climate","/","wc2.1_30s_bio_1.tif"))

# get extent of this raster
ex <- raster::extent(t)

#remove raster
rm(t)

# ==============================================================================
# 4. PARALLEL GBIF DOWNLOAD FUNCTION
# ==============================================================================

#' Number of parallel workers for GBIF downloads.
#' \strong{Do not use \code{detectCores() - 1}} — ask the server administrator
#' how many cores are available before changing this value.
numCores <- 6

#' Download GBIF occurrence records for a set of species in parallel.
#'
#' @description
#' Spawns a \code{parallel} socket cluster of \code{numCores} workers and uses
#' \code{parLapplyLB} (load-balanced) to distribute species-level GBIF queries
#' across workers. Each worker calls \code{dismo::gbif()} with up to 300 retries
#' to handle transient API failures, removes duplicate coordinates, and
#' standardises the returned columns. Results are row-bound into a single
#' \code{data.frame}.
#'
#' @param species_set A \code{data.frame} with columns \code{species},
#'   \code{genus}, and \code{specEpithet}. Typically produced by the
#'   binomial-splitting step above.
#' @param ex A \code{raster::Extent} object defining the geographic bounding
#'   box passed to \code{dismo::gbif(ext = )}.
#' @param numCores Integer. Number of parallel socket workers to use.
#'   \strong{Do not use \code{detectCores() - 1}} — confirm available cores
#'   with the server administrator before changing this value.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{\code{species}}{Accepted scientific name from GBIF taxonomy.}
#'   \item{\code{species_searched}}{Name as supplied in \code{species_set$species},
#'     retained for downstream synonym reconciliation.}
#'   \item{\code{lon}, \code{lat}}{Decimal longitude and latitude.}
#'   \item{\code{country}}{Country name from GBIF record.}
#'   \item{\code{year}}{Collection year (\code{NA} if not reported).}
#'   \item{\code{adm1}}{First-level administrative division (\code{NA} if absent).}
#'   \item{\code{locality}}{Locality description (\code{NA} if absent).}
#'   \item{\code{collection_code}}{Collection code (\code{NA} if absent).}
#'   \item{\code{dataset_name}}{Dataset name (\code{NA} if absent).}
#'   \item{\code{dataset_key}}{Dataset key (\code{NA} if absent).}
#'   \item{\code{institution_code}}{Institution code (\code{NA} if absent).}
#'   \item{\code{source}}{Always \code{"GBIF"}.}
#' }
#'
#' @note
#' \itemize{
#'   \item \code{Sys.sleep(3)} is called after each species download to reduce
#'     the risk of hitting GBIF rate limits across workers.
#'   \item If \code{dismo::gbif()} throws an error the species is skipped
#'     (\code{NULL} returned for that element); \code{do.call(rbind, ...)}
#'     silently drops \code{NULL} entries.
#'   \item Optional columns (year, adm1, locality, etc.) default to
#'     \code{rep(NA, nrow(GBIF))} when absent from the API response.
#'   \item Duplicate detection uses \code{duplicated(GBIF$lat, GBIF$lon)};
#'     note this deduplicates on \code{lat} only due to how \code{duplicated()}
#'     handles multiple arguments — verify if strict lat+lon deduplication is
#'     required.
#'   \item \code{terra} objects (\code{SpatRaster}, \code{SpatVector}) must
#'     \strong{not} be exported to parallel workers — pass file paths and
#'     reload inside workers if rasters are needed.
#' }
gbif_parallel_function <- function(species_set,ex,numCores){
  #adding this to allow parallel working for R>=3.6
  rscript_args = c("-e", shQuote("getRversion"))
  #preparing n cores
  cl <- parallel::makeCluster(numCores)
  #exporting objects to each of the n cores used
  parallel::clusterExport(cl, varlist=c("species_set",
                                        "ex"),envir=environment())
  
  # download the GBIF data with parLapplyLB (load-balanced parallel lapply)
  occGBIF <- parallel::parLapplyLB(
    cl=cl,
    X = seq_len(length(species_set$species)),
    fun = function (i){
      
      # keep track of loop progress
      message(paste0("species number ", i, ": ", species_set$species[i], sep=""))
      
      # download GBIF data (t_1 to avoid overlap with R's transpose function)
      # ntries = 300: increased from default 100 to handle intermittent API failures
      t_1 <- try(
        GBIF <- dismo::gbif(
          genus = as.character(species_set$genus[i]),
          species = as.character(species_set$specEpithet[i]),
          ext = ex,
          geo = TRUE,
          removeZeros = TRUE,
          download = TRUE,
          ntries = 300
        ))
      
      if (class(t_1) == "try-error") {
        GBIF <- NULL
      }
      
      # pause between requests to reduce GBIF rate-limit pressure across workers
      Sys.sleep(time = 3)
      
      if (!is.null(GBIF)) {
        
        # delete duplicate coordinates
        j <- which(duplicated(GBIF$lat, GBIF$lon))
        if(length(j) > 0) {
          GBIF <- GBIF[-j,]
        }
        
        # add queried name alongside GBIF-accepted name for downstream
        # species matching and synonym reconciliation
        GBIF$species_searched <- rep(species_set$species[i], nrow(GBIF))
        
        # get the relevant variables
        species <- GBIF$acceptedScientificName
        species_searched <- GBIF$species_searched
        lon <- GBIF$lon
        lat <- GBIF$lat
        country <- GBIF$country
        
        # the following variables may be missing so select only if available
        j <- which(names(GBIF) == "year")
        if (length(j) > 0) {
          year <- GBIF$year
        } else {
          year <- rep(NA, nrow(GBIF))
        }
        j <- which(names(GBIF) == "country")
        if (length(j) > 0) {
          country_name <- GBIF$country
        } else {
          country_name <- rep(NA, nrow(GBIF))
        }
        j <- which(names(GBIF) == "adm1")
        if (length(j) > 0) {
          adm1 <- GBIF$adm1
        } else {
          adm1 <- rep(NA, nrow(GBIF))
        }
        j <- which(names(GBIF) == "locality")
        if (length(j) > 0) {
          locality <- GBIF$locality
        } else {
          locality <- rep(NA, nrow(GBIF))
        }
        j <- which(names(GBIF) == "collectionCode")
        if (length(j) > 0) {
          collection_code <- GBIF$collectionCode
        } else {
          collection_code <- rep(NA, nrow(GBIF))
        }
        j <- which(names(GBIF) == "datasetName")
        if (length(j) > 0) {
          dataset_name <- GBIF$datasetName
        } else {
          dataset_name <- rep(NA, nrow(GBIF))
        }
        j <- which(names(GBIF) == "datasetKey")
        if (length(j) > 0) {
          dataset_key <- GBIF$datasetKey
        } else {
          dataset_key <- rep(NA, nrow(GBIF))
        }
        j <- which(names(GBIF) == "institutionCode")
        if (length(j) > 0) {
          institution_code <- GBIF$institutionCode
        } else {
          institution_code <- rep(NA, nrow(GBIF))
        }
        
        GBIF <- data.frame(species = species,
                           species_searched = species_searched,
                           lon = as.numeric(as.character(lon)),
                           lat = as.numeric(as.character(lat)),
                           country = as.character(country_name),
                           year = as.numeric(as.character(year)),
                           adm1 = as.character(adm1),
                           locality = as.character(locality),
                           collection_code = as.character(collection_code),
                           dataset_name = as.character(dataset_name),
                           dataset_key = as.character(dataset_key),
                           institution_code = as.character(institution_code),
                           source = "GBIF")
      }
      
    })
  parallel::stopCluster(cl)
  # rbind the list; NULL elements (failed downloads) are dropped silently
  occGBIF_rbind <- do.call(rbind, occGBIF)
  return(occGBIF_rbind)
}

rscript_args = c("-e", shQuote("getRversion()"))

# ==============================================================================
# 5. ITERATIVE DOWNLOAD WITH RETRY ROUNDS
# ==============================================================================
# GBIF API failures mean some species return no records in any given attempt.
# The pattern below re-runs gbif_parallel_function() on the remaining missing
# species after each round, accumulating successfully downloaded names in un_A,
# until all species have been attempted or the retry budget is exhausted.
# Rounds 1-11 use numCores = 6; rounds 12-18 use numCores = 12.
# Output CSVs are written after each round so partial results are never lost.
# ==============================================================================
#creating required dirs
original_dir <- paste0(basedir, "/1.Processing/occurrences/", "occurrences", "/original")
if(!dir.exists(original_dir)){
  dir.create(original_dir,recursive = T)
}

cleaned_dir <- paste0(basedir, "/1.Processing/occurrences/", "occurrences", "/cleaned")
if(!dir.exists(cleaned_dir)){
  dir.create(cleaned_dir,recursive = T)
}
# ==============================================================================
setwd(dir = paste0(basedir, "/1.Processing/occurrences/", "occurrences", "/original"))

# -- Round 1 -------------------------------------------------------------------
occGBIF_rbind <- gbif_parallel_function(species_set,ex,numCores)
write.csv(occGBIF_rbind, paste0("GBIF_data_download_", "occurrences_list_second_round", ".csv"))

# -- Round 2: retry species absent from round 1 --------------------------------

#' Check how many species are still missing after round 1
sum(!species_set$species %in% unique(occGBIF_rbind$species_searched))

#' Indices of species not yet downloaded
sp_n <- 1:length(species_set$species)
sp_n <- sp_n[!species_set$species %in% unique(occGBIF_rbind$species_searched)]

missing_sp <- species_set[sp_n,]
occGBIF_rbind2 <- gbif_parallel_function(missing_sp,ex,numCores)
write.csv(occGBIF_rbind2, paste0("GBIF_data_download_", "occurrences_list_second_round_2", ".csv"))

# -- Round 3 -------------------------------------------------------------------

#' Check coverage after round 2
sum(!species_set$species %in% unique(occGBIF_rbind2$species_searched))

sp_n <- 1:length(species_set$species)
sp_n <- sp_n[!species_set$species %in% unique(occGBIF_rbind2$species_searched)]

#' Accumulate all successfully downloaded species names across rounds 1-2
un1 <- unique(occGBIF_rbind$species_searched)
un2 <- unique(occGBIF_rbind2$species_searched)
un_A <- c(un1,un2)

missing_sp2 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind3 <- gbif_parallel_function(missing_sp2,ex,numCores)
write.csv(occGBIF_rbind3, paste0("GBIF_data_download_", "occurrences_list_second_round_3", ".csv"))

# -- Round 4 -------------------------------------------------------------------
un3 <- unique(occGBIF_rbind3$species_searched)
un_A <- c(un_A,un3)
missing_sp3 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind4 <- gbif_parallel_function(missing_sp3,ex,numCores)
write.csv(occGBIF_rbind4, paste0("GBIF_data_download_", "occurrences_list_second_round_4", ".csv"))

# -- Round 5 -------------------------------------------------------------------
un4 <- unique(occGBIF_rbind4$species_searched)
un_A <- c(un_A,un4)
missing_sp4 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind5 <- gbif_parallel_function(missing_sp4,ex,numCores)
write.csv(occGBIF_rbind5, paste0("GBIF_data_download_", "occurrences_list_second_round_5", ".csv"))

# -- Round 6 -------------------------------------------------------------------
un5 <- unique(occGBIF_rbind5$species_searched)
un_A <- c(un_A,un5)
missing_sp5 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind6 <- gbif_parallel_function(missing_sp5,ex,numCores)
write.csv(occGBIF_rbind6, paste0("GBIF_data_download_", "occurrences_list_second_round_6", ".csv"))

# -- Round 7 -------------------------------------------------------------------
un6 <- unique(occGBIF_rbind6$species_searched)
un_A <- c(un_A,un6)
missing_sp6 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind7 <- gbif_parallel_function(missing_sp6,ex,numCores)
write.csv(occGBIF_rbind7, paste0("GBIF_data_download_", "occurrences_list_second_round_7", ".csv"))

# -- Round 8 -------------------------------------------------------------------
un7 <- unique(occGBIF_rbind7$species_searched)
un_A <- c(un_A,un7)
missing_sp7 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind8 <- gbif_parallel_function(missing_sp7,ex,numCores)
write.csv(occGBIF_rbind8, paste0("GBIF_data_download_", "occurrences_list_second_round_8", ".csv"))

# -- Round 9 -------------------------------------------------------------------
un8 <- unique(occGBIF_rbind8$species_searched)
un_A <- c(un_A,un8)
missing_sp8 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind9 <- gbif_parallel_function(missing_sp8,ex,numCores)
write.csv(occGBIF_rbind9, paste0("GBIF_data_download_", "occurrences_list_second_round_9", ".csv"))

# -- Round 10 ------------------------------------------------------------------
un9 <- unique(occGBIF_rbind9$species_searched)
un_A <- c(un_A,un9)
missing_sp9 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind10 <- gbif_parallel_function(missing_sp9,ex,numCores)
write.csv(occGBIF_rbind10, paste0("GBIF_data_download_", "occurrences_list_second_round_10", ".csv"))

# -- Round 11 ------------------------------------------------------------------
un10 <- unique(occGBIF_rbind10$species_searched)
un_A <- c(un_A,un10)
missing_sp10 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind11 <- gbif_parallel_function(missing_sp10,ex,numCores)
write.csv(occGBIF_rbind11, paste0("GBIF_data_download_", "occurrences_list_second_round_11", ".csv"))

# -- Rounds 12-18: increased core count (12) for remaining species -------------

# -- Round 12 ------------------------------------------------------------------
un11 <- unique(occGBIF_rbind11$species_searched)
un_A <- c(un_A,un11)
missing_sp11 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind12 <- gbif_parallel_function(missing_sp11,ex,12)
write.csv(occGBIF_rbind12, paste0("GBIF_data_download_", "occurrences_list_second_round_12", ".csv"))

# -- Round 13 ------------------------------------------------------------------
un12 <- unique(occGBIF_rbind12$species_searched)
un_A <- c(un_A,un12)
missing_sp12 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind13 <- gbif_parallel_function(missing_sp12,ex,12)
write.csv(occGBIF_rbind13, paste0("GBIF_data_download_", "occurrences_list_second_round_13", ".csv"))

# -- Round 14 ------------------------------------------------------------------
un13 <- unique(occGBIF_rbind13$species_searched)
un_A <- c(un_A,un13)
missing_sp13 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind14 <- gbif_parallel_function(missing_sp13,ex,12)
write.csv(occGBIF_rbind14, paste0("GBIF_data_download_", "occurrences_list_second_round_14", ".csv"))

# -- Round 15 ------------------------------------------------------------------
un14 <- unique(occGBIF_rbind14$species_searched)
un_A <- c(un_A,un14)
missing_sp14 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind15 <- gbif_parallel_function(missing_sp14,ex,12)
write.csv(occGBIF_rbind15, paste0("GBIF_data_download_", "occurrences_list_second_round_15", ".csv"))

# -- Round 16 ------------------------------------------------------------------
un15 <- unique(occGBIF_rbind15$species_searched)
un_A <- c(un_A,un15)
missing_sp15 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind16 <- gbif_parallel_function(missing_sp15,ex,12)
write.csv(occGBIF_rbind16, paste0("GBIF_data_download_", "occurrences_list_second_round_16", ".csv"))

# -- Round 17 ------------------------------------------------------------------
un16 <- unique(occGBIF_rbind16$species_searched)
un_A <- c(un_A,un16)
missing_sp16 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind17 <- gbif_parallel_function(missing_sp16,ex,12)
write.csv(occGBIF_rbind17, paste0("GBIF_data_download_", "occurrences_list_second_round_17", ".csv"))

# -- Round 18 ------------------------------------------------------------------
un17 <- unique(occGBIF_rbind17$species_searched)
un_A <- c(un_A,un17)
missing_sp17 <- species_set[!species_set$species %in% un_A,]
occGBIF_rbind18 <- gbif_parallel_function(missing_sp17,ex,12)
write.csv(occGBIF_rbind18, paste0("GBIF_data_download_", "occurrences_list_second_round_18", ".csv"))