# load packages
library(terra)
library(sf)
library(ggplot2)
library(bit)
require(leaflet)
library(htmlwidgets)
################################################################################
# set base directory
basedir <- "E:/CSOSA/Dropbox/VAVILOV_2.0"
occDir <- paste0(basedir, "/1.Processing/occurrences/cleaned")
out_dir <- paste0(basedir, "1.Processing")
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}
#richness
out_buffer_dir <- paste0(out_dir, "/buffer_richness")
if (!dir.exists(out_buffer_dir)) {
  dir.create(out_buffer_dir)
}

#species_raster_dir
out_buffer_dir_raster <- paste0(out_buffer_dir, "/raster")
if (!dir.exists(out_buffer_dir_raster)) {
  dir.create(out_buffer_dir_raster)
}

#summary_dir
out_buffer_dir_summary <- paste0(out_buffer_dir, "/summary")
if (!dir.exists(out_buffer_dir_summary)) {
  dir.create(out_buffer_dir_summary)
}

# chunks summary_dir
out_buffer_dir_summary_chunks <- paste0(out_buffer_dir_summary, "/chunks")
if (!dir.exists(out_buffer_dir_summary_chunks)) {
  dir.create(out_buffer_dir_summary_chunks)
}
################################################################################

#get species list
species_selected <- read.csv(
  paste0(out_buffer_dir_summary, "/", "species_selected.csv")
)

################################################################################
# species_selected <- species_selected[1:50,]
#load template for rasterize
template <-
  rast(paste0(
    basedir,
    "/DATA/input_data/climate_data/present/",
    "wc2.1_30s_bio_1.tif"
  ))

################################################################################
# #loading first raster to obtain sum step by step
# x <- terra::rast(paste0(out_buffer_dir_raster,"/",species_selected$species[[1]],".tif"))
# 
#   
# pb <-
#   utils::txtProgressBar(min = 2,
#                         max = nrow(species_selected),
#                         style = 3)
# 
# #sum for each new raster object
# for(i in 2:nrow(species_selected)){
#   utils::setTxtProgressBar(pb, i)
#   
#   x_i <- terra::rast(paste0(out_buffer_dir_raster,"/",species_selected$species[[i]],".tif"))
#   x <- sum(c(x, x_i), na.rm = TRUE)
#    rm(x_i)
#   # ;gc(verbose = F)
# }
# close(pb)
################################################################################
#load rasters

# if(!file.exists(paste0(out_buffer_dir_summary, "/richness_map.tif"))){

current<-list.files(out_buffer_dir_raster,pattern =".tif" ,full.names = T)
current_name<-list.files(out_buffer_dir_raster,pattern = ".tif",full.names = F)
current<-lapply(current,rast)

#saving in a stack object
r_stack <- terra::rast(current)

#creating chunks to sum subsets (each 50 layers is created a chunk)
x_chunks <- bit::chunks(1,  nlyr(r_stack), by=50)
message(paste("Total Chunks: ",length(x_chunks)))
#for each chunk

message(paste("Getting raster sum for all chunks"))

for(i in 1:length(x_chunks)){
  # i <- 1
  #obtaining chunk boundaries
  x_i1 <- as.numeric(as.character(x_chunks[[i]])[1])
  x_i2 <- as.numeric(as.character(x_chunks[[i]])[2])
  #sum of each chunks
  x <- sum(r_stack[[x_i1:x_i2]],na.rm = TRUE)
  #writing in a raster
  writeRaster(x,paste0(out_buffer_dir_summary_chunks,"/",i,".tif"),overwrite=T)
};rm(i)

################################################################################
#calling all chunk results
message(paste("Obtaining final summary file"))
x_sum <- list.files(out_buffer_dir_summary_chunks,".tif",full.names = T)
x_sum <- lapply(x_sum,rast)
x_sum <- terra::rast(x_sum)
#getting summary richness
x_sum <- sum(x_sum,na.rm = T)
################################################################################
#save results
message(paste("Saving raster file"))
terra::writeRaster(x_sum,paste0(out_buffer_dir_summary, "/richness_map.tif"),overwrite=T)
# } else {
#   x_sum <- terra::rast(paste0(out_buffer_dir_summary, "/richness_map.tif"))
# }
################################################################################
#calling admin0 for plot
adm0 <- sf::st_read(paste0(basedir,"/DATA/Input_data/adm0/adm0_Latam.shp"))
################################################################################
# # Plot the richness map with administrative boundaries
png(paste0(out_buffer_dir_summary, "/richness_map.png"), 
    width = 10, height = 8, units = "in", res = 300)
plot(x_sum, main = "Species Richness", 
     col = hcl.colors(100, "viridis"),
     axes = TRUE)
plot(st_geometry(adm0), add = TRUE, border = "black", lwd = 1.5)
dev.off()
################################################################################
# occDir <- paste0(basedir, "/DATA/DATABASES/occurrences/cleaned")
# occs <-
#   read.csv(
#     paste0(occDir, "/", "GBIF_data", "_year_cleaned.csv"),
#     header = T,
#     row.names = 1
#   )
# occs <- occs[, -c(1, 2)]
# #fixing column names
# colnames(occs) <- c(
#   "species",
#   "species_searched",
#   "country",
#   "year",
#   "adm1",
#   "locality",
#   "collection_code",
#   "dataset_name",
#   "dataset_key",
#   "institution_code",
#   "source",
#   "GID_0",
#   "NAME_0",
#   "geometry_long",
#   "geometry_lat",
#   "lon",
#   "lat"
# )
# occs <- occs[,c("species_searched","lon","lat")]

pal <- colorNumeric(c("#0C2C84","#41B6C4","#8DEEEE","#FFFFCC","#EEE8CD","#EEB422","red"), values(x_sum),
                    na.color = "transparent")
# 
# %>% addTiles() %>%
# #   addRasterImage(r, colors = pal, opacity = 0.8) %>%
#   addLegend(pal = pal, values = values(x_sum),
#             title = "Surface temp")

x <- leaflet() %>%
  addTiles() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addMiniMap(width = 150, height = 150) %>%
  addRasterImage(x_sum,
                 colors = pal, opacity = 0.3,
                 project = TRUE,
                 maxBytes = 20 * 1024^2) %>%  # Increase to 20MB
  addLegend(pal = pal, values = values(x_sum),
            title = "Species Richness")

saveWidget(x, paste0(out_buffer_dir_summary,"/","buffer_richness.html"), selfcontained = T)
