### This script cleans occurrence data (GBIF and BIEN data typically) using the following steps:
# remove data in urban areas
# remove records with imprecise coordinates
# remove data older than 1950
# remove data with region mismatch

# author: Tobias Fremout (tobias.fremout@gmail.com)

# load packages
library(terra)
library(sf)

# set base directory
basedir <- "E:/CSOSA/Dropbox/VAVILOV_2.0"

# ==============================================================================
# 1. LOAD SPATIAL REFERENCE LAYERS
# ==============================================================================

#' Latin-American country boundaries used for region-mismatch filtering.
#' CRS is reassigned to WGS84 before spatial intersection.
adm0 <- st_read(paste0(basedir,"/DATA/input_data/","adm0/", "adm0_Latam", ".shp"))

#' WorldClim BIO1 at 30 s resolution used as spatial template to crop the
#' urban areas raster to the modelling extent.
t <- rast(paste0(basedir,"/DATA/input_data/climate_data/present/","wc2.1_30s_bio_1.tif"))

# load urban areas raster
#' GHSL Global Human Settlement Layer — degree of urbanisation (SMOD) built-up
#' estimates at 30 arc-second resolution, 2014. Values > 25 indicate urban
#' cells; records falling in these cells are removed.
urban <- rast(paste0(basedir,"/DATA/input_data/GHSL urban areas/",
                     "ghsl-population-built-up-estimates-degree-urban-smod_built-30ss-2014.tif"))

# crop urban raster to the modelling extent defined by the climate template
urban <- crop(urban, t)

# ==============================================================================
# 2. LOAD OCCURRENCE DATA
# ==============================================================================

setwd(dir = paste0(basedir, "/DATA/DATABASES/occurrences/original"))

#' Load second-round GBIF occurrence records.
#' The first-round file (GBIF_data_download_occurrences_list.csv) is retained
#' as a comment for reference.
# occ <- read.csv(paste0("GBIF_data_download", "_occurrences_list", ".csv"))
occ <- read.csv(paste0("GBIF_data_download_second_round_occurrences_list", ".csv"))

# check structure of loaded data
names(occ)
dim(occ)

# ==============================================================================
# 3. REMOVE RECORDS IN URBAN AREAS
# ==============================================================================

#' Extract GHSL urbanisation values at occurrence coordinates.
#' Records with GHSL value > 25 (urban) are removed.
extr <- terra::extract(urban, cbind(occ$lon, occ$lat))
i <- which(extr > 25)
if(length(i)>0) {
  occ <- occ[-i,]
}
print(paste0(length(i), " points deleted because located in urban areas"))

# ==============================================================================
# 4. COERCE COORDINATE AND YEAR COLUMNS
# ==============================================================================

#' Ensure lon, lat, and year are numeric; coercion may produce NAs if the
#' source CSV stored these as character (e.g. empty strings).
occ$lat <- as.numeric(occ$lat)
occ$lon <- as.numeric(occ$lon)
occ$year <- as.numeric(occ$year)

# ==============================================================================
# 5. REMOVE RECORDS WITH IMPRECISE COORDINATES
# ==============================================================================

#' Count the number of decimal places in a numeric value.
#'
#' @param x A single numeric value.
#' @return Integer. Number of decimal places; 0 for whole numbers.
#' @details Uses \code{format(x, scientific = FALSE)} to avoid scientific
#'   notation before splitting on the decimal point.
decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    strs <- strsplit(as.character(format(x, scientific = F)), "\\.")
    n <- nchar(strs[[1]][2])
  } else {
    n <- 0
  }
  return(n) 
}

# count decimal places of lon and lat coordinates
dec_lon <- sapply(occ$lon, FUN = decimalplaces)
dec_lat <- sapply(occ$lat, FUN = decimalplaces)

#' Remove records with zero decimal places in lon OR lat.
#' Rationale: coordinates rounded to whole degrees have a precision of ~111 km;
#' the probability of such a record genuinely having .00 decimals is 1/100.
i <- which(dec_lon == 0 | dec_lat == 0)
if(length(i)>0) {
  occ <- occ[-i,]
}
print(paste0(length(i), " points deleted because of imprecise coordinates"))

#' Remove records with exactly 1 decimal place in BOTH lon AND lat.
#' Rationale: 1 decimal place gives ~11 km precision; the probability of a
#' record being precise enough under this condition is (1/10)^2 = 1/100.
i <- which(dec_lon == 1 & dec_lat == 1)
if(length(i)>0) {
  occ <- occ[-i,]
}
print(paste0(length(i), " points deleted because of imprecise coordinates"))

# ==============================================================================
# 6. (OPTIONAL) REMOVE RECORDS OLDER THAN 1950
# ==============================================================================
# Records collected before 1950 have a higher probability of imprecise
# georeferencing. This filter is currently disabled; applied in section 8
# on a copy of the cleaned dataset.
#
# i <- which(occ$year < 1950)
# if(length(i)>0) {
#   occ <- occ[-i,]
# }
# print(paste0(length(i), " points deleted because older than 1950"))

# ==============================================================================
# 7. REMOVE RECORDS WITH REGION MISMATCH
# ==============================================================================

#' Project reference CRS (WGS84 geographic)
crs.geo <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#' Convert occurrence data.frame to an \code{sf} point object for spatial
#' intersection with country boundaries.
occSp <- sf::st_as_sf(x = occ, coords = c("lon", "lat"), crs = crs.geo)

#' Assign WGS84 CRS to adm0 boundaries (ensures CRS consistency before
#' intersection; does not reproject).
st_crs(adm0) <- crs.geo

#' Spatial intersection retains only points that fall within a country polygon
#' and appends the polygon attribute \code{NAME_0} (country name from adm0).
overlay <- sf::st_intersection(occSp, adm0)
occ <- overlay

#' Remove records where the GBIF-reported country does not match the country
#' polygon the point falls in (\code{NAME_0}).
j <- which(occ$country != occ$NAME_0)
if(length(j)>0) {
  occ <- occ[-j,]
}
print(paste0(length(j), " points deleted because of region mismatch"))

# ==============================================================================
# 8. RECONSTRUCT DATA FRAME AND WRITE CLEANED OUTPUT
# ==============================================================================

#' Extract coordinates from the \code{sf} geometry column and re-attach them
#' as plain numeric columns after the spatial operations.
coords <- sf::st_coordinates(occ)

occ_df <- as.data.frame(occ)
occ_df$lon <- coords[,1]
occ_df$lat <- coords[,2]
dim(occ_df)

rm(overlay, occSp)
gc()

setwd(dir = paste0(basedir, "/DATA/DATABASES/occurrences/cleaned"))

#' Write the spatially cleaned occurrence dataset (urban, precision, and region
#' filters applied; year filter not yet applied).
# write.csv(occ_df, paste0("GBIF_data", "_cleaned.csv"))
write.csv(occ_df, paste0("GBIF_data_second_round", "_cleaned.csv"))

rm(occ, dec_lat, dec_lon, i, extr, j, coords); gc()

# ==============================================================================
# 9. ADDITIONAL YEAR FILTER  (>= 1950, non-NA)
# ==============================================================================

#' Apply the year filter on a copy of the cleaned dataset so that the
#' all-years version (section 8) is preserved separately.
occ_df2 <- occ_df

#' Remove records collected before 1950.
#' High probability of imprecise georeferencing in pre-1950 herbarium records.
i <- which(occ_df2$year < 1950)
if(length(i)>0) {
  occ_df2 <- occ_df2[-i,]
}

#' Remove records with missing year information.
occ_df2 <- occ_df2[which(!is.na(occ_df2$year)),]

print(paste0(nrow(occ_df2), " points deleted because older than 1950 and NA"))

#' Write the year-filtered occurrence dataset.
write.csv(occ_df2, paste0("GBIF_data_second_round", "_year_cleaned.csv"))

rm(occ, dec_lat, dec_lon, i, extr, j, coords, occ_df2); gc()