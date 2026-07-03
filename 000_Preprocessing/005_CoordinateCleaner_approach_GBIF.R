# Load necessary libraries
require(CoordinateCleaner)  # For automated cleaning of geographic coordinates in species occurrence data
require(readxl)             # For reading Excel files (not used in this code but loaded)
require(data.table)         # For efficient data manipulation and fast file reading

# ==============================================================================
# 1. COORDINATE CLEANING FUNCTION
# ==============================================================================

#' Clean Species Occurrence Coordinates Using CoordinateCleaner
#'
#' @description
#' Applies a battery of automated geographic coordinate tests via
#' \code{CoordinateCleaner::clean_coordinates()} to flag and remove
#' problematic occurrence records. The following tests are applied:
#' \describe{
#'   \item{\code{capitals}}{Flags points within 5000 m of country capitals.}
#'   \item{\code{centroids}}{Flags points within 5000 m of country centroids.}
#'   \item{\code{equal}}{Flags points where lat == lon.}
#'   \item{\code{gbif}}{Flags points outside country borders according to GBIF.}
#'   \item{\code{institutions}}{Flags points within 5000 m of biodiversity
#'     institutions (herbaria, museums).}
#'   \item{\code{outliers}}{Flags geographic outliers within a species'
#'     occurrence point cloud.}
#'   \item{\code{seas}}{Flags points located in the sea.}
#'   \item{\code{zeros}}{Flags points with zero coordinates (0, 0).}
#' }
#'
#' @param x A \code{data.frame} or \code{data.table} of occurrence records
#'   containing at minimum columns \code{lon}, \code{lat}, \code{country},
#'   and \code{species}.
#' @param basedir Character. Project base directory (currently unused inside
#'   the function body but retained as a parameter for future extension).
#' @param cleanDir Character. Path to the cleaned occurrences directory
#'   (currently unused inside the function body but retained for future
#'   extension, e.g. writing flagged records to disk).
#'
#' @return A \code{data.frame} containing only records that passed all
#'   coordinate tests (\code{flags$.summary == TRUE}).
#'   Flagged records (\code{dat_fl}) are computed but not returned; uncomment
#'   the relevant line to save them separately if needed.
#'
#' @note
#' \code{capitals_rad}, \code{inst_rad}, and \code{centroids_rad} are all set
#' to 5000 m. Adjust these radii if a stricter or more lenient spatial buffer
#' around reference features is required.
clean_coordinates <- function(x, basedir, cleanDir) {
  # Apply coordinate cleaning tests to flag problematic records
  # Tests include:
  # - capitals: flags points near country capitals (within 5000m)
  # - centroids: flags points near country centroids (within 5000m)
  # - equal: flags points with equal lat/lon coordinates
  # - gbif: flags points outside country borders according to GBIF
  # - institutions: flags points near biodiversity institutions (within 5000m)
  # - outliers: flags geographic outliers within species occurrence points
  # - seas: flags points located in the sea
  # - zeros: flags points with zero coordinates
  flags <- CoordinateCleaner::clean_coordinates(
    x = x,
    lon = "lon",
    lat = "lat",
    countries = "country",
    capitals_rad = 5000,
    inst_rad = 5000,
    centroids_rad = 5000,
    species = "species",
    tests = c(
      "capitals",
      "centroids",
      "equal",
      "gbif",
      "institutions",
      "outliers",
      "seas",
      "zeros"
    )
  )
  
  # Subset data to exclude flagged problematic records
  dat_cl <- x[flags$.summary, ]   # Good records: passed all tests
  dat_fl <- x[!flags$.summary, ]  # Flagged records: failed one or more tests
  # (not returned; save separately if needed)
  
  # Return cleaned dataset with flagged records removed
  return(dat_cl)
}

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

#' Project base directory
basedir <- "E:/CSOSA/Dropbox/VAVILOV_2.0"

#' Directory containing cleaned occurrence CSVs produced by the previous
#' spatial cleaning step (urban areas, precision, region mismatch filters)
cleanDir <- paste0(basedir, "/", "DATA/DATABASES/occurrences/cleaned")

# ==============================================================================
# 3. APPLY COORDINATE CLEANING PER DOWNLOAD ROUND
# ==============================================================================

# -- Round 1 (commented out; already processed) --------------------------------
# x <- fread(paste0(cleanDir, "/", "GBIF_data_cleaned.csv"))
# dat_cl <- clean_coordinates(x, basedir, cleanDir)
# write.csv(dat_cl, paste0(cleanDir, "/", "GBIF_data_cleaned_FINAL.csv"))

# -- Round 2 (commented out; already processed) --------------------------------
# x <- fread(paste0(cleanDir, "/", "GBIF_data_second_round_cleaned.csv"))
# dat_cl <- clean_coordinates(x, basedir, cleanDir)
# write.csv(dat_cl, paste0(cleanDir, "/", "GBIF_data_second_round_cleaned_FINAL.csv"))

# -- Round 3 -------------------------------------------------------------------

#' Load third-round spatially cleaned occurrence records and apply the
#' CoordinateCleaner battery. Output is written as the FINAL cleaned file
#' for this round before joining with previous rounds.
x <- fread(paste0(cleanDir, "/", "GBIF_data_third_round_cleaned.csv"))
dat_cl <- clean_coordinates(x, basedir, cleanDir)
write.csv(dat_cl, paste0(cleanDir, "/", "GBIF_data_third_round_cleaned_FINAL.csv"))

# Uncomment below to rename columns or run external scripts if needed
# colnames(dat_cl)[1:3]<- c("scientific_name","latitude","longitude")
# write.csv(dat_cl,"D:/BOLDER/bactris_gasipaes/occurences.csv")

# ==============================================================================
# 4. JOIN ROUNDS 1 AND 2 INTO A SINGLE DATASET
# ==============================================================================

#' Load the CoordinateCleaner-cleaned FINAL files for rounds 1 and 2 and
#' row-bind them into a single joined dataset.
x1 <- fread(paste0(cleanDir, "/", "GBIF_data_cleaned_FINAL.csv"))
x2 <- fread(paste0(cleanDir, "/", "GBIF_data_second_round_cleaned_FINAL.csv"))

x3 <- rbind(x1, x2)
write.csv(x3, paste0(cleanDir, "/", "GBIF_joined_FINAL.csv"))

# ==============================================================================
# 5. REMOVE PROBLEMATIC COORDINATES AND HIGH-FREQUENCY DUPLICATE POINTS
# ==============================================================================
# A known problematic coordinate cluster associated with the Field Museum
# (-65.9119, -7.0758) and other high-frequency duplicate points (> 100
# occurrences at the same coordinate pair) are removed in this step.
# The branching logic for GBIF_joined_FINAL_2.csv is retained as a comment
# for reference; the active block processes the third-round FINAL file.

require(data.table)

#' Load third-round FINAL cleaned occurrences for duplicate-coordinate removal.
#' The joined FINAL file path is retained as a comment for reference.
# occs <- fread("E:/CSOSA/Dropbox/VAVILOV_2.0/DATA/DATABASES/occurrences/cleaned/GBIF_joined_FINAL.csv",
occs <- fread("E:/CSOSA/Dropbox/VAVILOV_2.0/DATA/DATABASES/occurrences/cleaned/GBIF_data_third_round_cleaned_FINAL.csv",
              header = T)

#' Drop leading index columns introduced by successive write.csv() calls.
#' The commented line retains an alternative column range for reference.
# occs <- occs[, -c(1, 2,3,4,5)]
occs <- occs[, -c(1, 2, 3, 4)]

#' Rename columns to a consistent schema after dropping the index columns.
#' Column order reflects the structure of the GBIF download + spatial join output.
colnames(occs) <- c(
  "species",
  "species_searched",
  "country",
  "year",
  "adm1",
  "locality",
  "collection_code",
  "dataset_name",
  "dataset_key",
  "institution_code",
  "source",
  "GID_0",
  "NAME_0",
  "geometry_long",
  "geometry_lat",
  "lon",
  "lat"
)

# Commented-out branching block retained for reference:
# handled the Field Museum coordinate cluster and records with > 200 identical
# coordinate pairs across two alternative input files (GBIF_joined_FINAL_2.csv).
# if(file.exists(...)) { ... } else { ... }

#' Build a coordinate string key (lon - lat) for each record to detect
#' high-frequency duplicate coordinate pairs.
occs$coords <- paste(occs$lon, "-", occs$lat)

#' Count occurrences per unique coordinate string
x <- tapply(occs$coords, occs$coords, length)

#' Sort by frequency (descending) to inspect the most duplicated coordinates
x <- x[order(x, decreasing = T)]

#' Identify coordinate pairs with more than 100 occurrences — these are likely
#' georeferencing artefacts (e.g. institution centroids, rounded coordinates)
x2 <- x[x > 100]

#' Remove records matching high-frequency duplicate coordinate pairs
occs3 <- occs[!occs$coords %in% names(x2), ]

write.csv(occs3, paste0(cleanDir, "/", "GBIF_data_third_round_cleaned_FINAL_2.csv"))