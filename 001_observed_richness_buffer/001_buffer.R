#' Species Occurrence Processing and Buffer Rasterization Pipeline
#'
#' @description
#' This script processes geographic occurrence data for multiple species, 
#' projects coordinates, generates individual 10-kilometer buffer zones, 
#' rasterizes the spatial extents using a climate template, and outputs 
#' clean CSV records and GeoTIFF surfaces for subsequent richness analysis.
#'
#' @note 
#' Project: VAVILOV_2.0
#' Dependencies: library(terra), library(sf)

# load packages
library(terra)
library(sf)

# ==============================================================================
# 1. DIRECTORY CONFIGURATION & INITIALIZATION
# ==============================================================================

#' @details Setup base and sub-directory trees for input data and outputs.
basedir <- "E:/CSOSA/Dropbox/VAVILOV_2.0"
#getting occurrences folder
occDir <- paste0(basedir, "/1.Processing/occurrences/cleaned")
#This is the outcome folder
out_dir <- paste0(basedir, "/1.Processing")
#creating outcome folder if it is the  case
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

# richness base directory
out_buffer_dir <- paste0(out_dir, "/buffer_richness")
if (!dir.exists(out_buffer_dir)) {
  dir.create(out_buffer_dir)
}

# species points csv directory
out_buffer_dir_occ <- paste0(out_buffer_dir, "/occs")
if (!dir.exists(out_buffer_dir_occ)) {
  dir.create(out_buffer_dir_occ)
}

# species raster tif directory
out_buffer_dir_raster <- paste0(out_buffer_dir, "/raster")
if (!dir.exists(out_buffer_dir_raster)) {
  dir.create(out_buffer_dir_raster)
}

# summary reports directory
out_buffer_dir_summary <- paste0(out_buffer_dir, "/summary")
if (!dir.exists(out_buffer_dir_summary)) {
  dir.create(out_buffer_dir_summary)
}

# ==============================================================================
# 2. DATA IMPORT AND COLUMN STANDARDIZATION
# ==============================================================================

#' @importFrom utils read.csv
#' @param occDir Character path containing cleaned occurrence CSVs.
#' @return data.frame `occs` containing cleaned spatial coordinates and taxonomy.
occs <-
  read.csv(
    paste0(occDir, "/", "GBIF_data_third_round_cleaned_FINAL_2.csv"),
    header = T,
    row.names = 1
  )

#' @details Structural mapping of raw columns to standardized field names.
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

# ==============================================================================
# 3. TAXONOMIC EXCLUSION & MASTER LIST EXPORT
# ==============================================================================

#' @return data.frame `species_selected` unique list of evaluated target taxa.
species_selected <-
  data.frame(species = unique(occs$species_searched))

#' @importFrom utils write.csv
write.csv(
  species_selected,
  paste0(out_buffer_dir_summary, "/", "species_selected_2.csv"),
  row.names = F
)

# ==============================================================================
# 4. RASTER TEMPLATE & PROGRESS MONITORING INITIALIZATION
# ==============================================================================

#' @param template SpatRaster used as the reference grid for rasterization layout.
template <-
  rast(paste0(
    basedir,
    "/DATA/input_data/climate_data/30s/present/",
    "wc2.1_30s_bio_1.tif"
  ))

#' @param pb txtProgressBar CLI tracker initialized to scale with unique species count.
pb <-
  utils::txtProgressBar(
    min = 0,
    max = nrow(species_selected),
    style = 3)

# ==============================================================================
# 5. SPATIAL BUFFERING AND RASTERIZATION LOOP
# ==============================================================================

#' @description Iterative mapping block running over distinct species units.
#' @section Spatial workflow per iteration:
#' 1. Filter occurrences to target species `i`.
#' 2. Coerce coordinate vectors into a SpatVector object under WGS84 (EPSG:4326).
#' 3. Project to a metric framework (EPSG:3857) to perform spatial distance buffer calculations.
#' 4. Apply a uniform 10,000-meter (10 km) geometric buffer over occurrence vectors.
#' 5. Reproject geometry back into WGS84 geographic system.
#' 6. Discretize vector boundaries into grid pixels using the referenced raster template.
#' 7. Persist metrics as standardized GeoTIFF outputs.
#'
#' @importFrom terra vect project buffer rasterize writeRaster
#' @importFrom utils setTxtProgressBar
#' @return NULL (Writes outputs directly to disk and runs garbage collection (`gc()`) dynamically)
buffer_list <- lapply(1:nrow(species_selected), function(i) {
  utils::setTxtProgressBar(pb, i)
  
  # subsetting for a given species 
  sp_i_occs <-
    occs[which(occs$species_searched == species_selected$species[[i]]), ]
  
  # writing partial csv
  write.csv(
    sp_i_occs,
    paste0(out_buffer_dir_occ, "/", species_selected$species[[i]], ".csv"),
    row.names = F
  )
  
  # vectorizing occurrences into a points shapefile 
  points_i <-
    vect(sp_i_occs, geom = c("lon", "lat"), crs = "EPSG:4326")
  
  # Project to a metric CRS before buffering
  points_i_proj <- project(points_i, "EPSG:3857")  # Web Mercator
  
  # Create 10 km buffer (10,000 meters)
  b_i <- terra::buffer(points_i_proj, width = 10000)
  
  # project back to geographic coordinates
  b_i <- terra::project(b_i, "EPSG:4326")
  
  # rasterize
  r_i <- terra::rasterize(b_i, template)
  
  # save buffer raster as tif
  terra::writeRaster(r_i,
                     paste0(
                       out_buffer_dir_raster,
                       "/",
                       species_selected$species[[i]],
                       ".tif"
                     ),
                     overwrite=T
  )
  
  # memory management: clear temporal spatial arrays
  rm(points_i_proj,r_i,points_i,sp_i_occs);gc()
})

close(pb)