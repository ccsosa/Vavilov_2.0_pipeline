### this script can be used to obtain change maps, current, and climate change maps for each species!


# author: Chrystian Sosa

library(terra)
library(viridis)
library(sf)
library(foreach)
library(doParallel)
library(bit)
# library(ggplot2)
library(bit)
library(tidyterra)
require(leaflet)
library(htmlwidgets)
library(viridisLite)
library(ggplot2)



changes_excel_conv_obtention <- function(basedir,output_dir,resultdir,species_selected,conversion_table){
  # chunks summary_dir
  output_dir_maps  <-
    paste0(output_dir, "/Summary_maps_folder_CC")
  # if (!dir.exists(output_dir_maps)) {
  #   dir.create(output_dir_maps)
  # }
  
  # #PNG folder (MAP)
  # output_dir_maps_PNG  <-
  #   paste0(output_dir_maps, "/PNG_SP_CC_REALIZED")
  # if (!dir.exists(output_dir_maps_PNG)) {
  #   dir.create(output_dir_maps_PNG)
  # }
  # #TIF folder
  # output_dir_maps_TIF  <-
  #   paste0(output_dir_maps, "/TIF_SP_CC_REALIZED")
  # if (!dir.exists(output_dir_maps_TIF)) {
  #   dir.create(output_dir_maps_TIF)
  # }
  
  #change folder calculation (potential)
  output_dir_maps_AREAS_realized  <-
    paste0(output_dir_maps, "/changes_SP_CC_REALIZED")

  #getting hull method
  hull <- c("ConcaveHull","ConvexHull")
  
  #assigning SSPs
  SSPs <- c("ssp245","ssp370")
  
  
  sp_areas_results <- list()
  for(i in 1:nrow(conversion_table)){
    # i <- 3
    x_j_sp <-list()
    for(j in 1:2){
      # j <- 1
      x_k_sp <-list()
      for(k in 1:2){
        # k <- 1
        #folder to call
        areas_dir_j_k <- paste0(output_dir_maps_AREAS_realized,
                                "/",
                                hull[[j]],
                                "/",
                                SSPs[[k]])
        #getting results for hull and SSP for a given species
        if(file.exists(paste0(areas_dir_j_k,
                              "/",
                              conversion_table$species[[i]],
                              "_areas.csv"))){

          
          x <- read.csv(paste0(areas_dir_j_k,
                               "/",
                               conversion_table$species[[i]],"_areas.csv"),row.names = 1)  
          
           x_k_sp[[k]] <-  
             data.frame(species = conversion_table$species[[i]],
                     family = conversion_table$family[[i]],
                     genus = conversion_table$genus[[i]],
                     species_from_source = conversion_table$species_from_source[[i]],
                     taxonomic_status = conversion_table$taxonomic_status[[i]],
                     ch = hull[[j]],
                     SSP = SSPs[[k]],
                     loss = x$loss, 
                     stable = x$stable,
                     new_areas = x$new_areas)


        } else {
          #In case of no approach
          x_k_sp[[k]] <- data.frame(species = conversion_table$species[[i]],
                                    family = conversion_table$family[[i]],
                                    genus = conversion_table$genus[[i]],
                                    species_from_source = conversion_table$species_from_source[[i]],
                                    taxonomic_status = conversion_table$taxonomic_status[[i]],
                                    ch = hull[[j]],
                                    SSP = SSPs[[k]],
                                    loss = NA, 
                                    stable = NA,
                                    new_areas = NA)
        }
        
      } #k
      #summarizing SSPs
      x_k_sp <- do.call(rbind,x_k_sp)
      x_j_sp[[j]] <- x_k_sp
    } #j
    x_j_sp <- do.call(rbind,x_j_sp)
    sp_areas_results[[i]] <- x_j_sp
  }
  
  sp_areas_results <- do.call(rbind,sp_areas_results)

  writexl::write_xlsx(sp_areas_results,paste0(output_dir_maps_AREAS_realized,"/","sp_areas_results.xlsx"))

  return("DONE!")
}
################################################################################
user <- "Chrystian"
if(user == "Chrystian") {
  basedir <-  "//catalogue/MultifLandscapesA1706"
}

resultdir <- "//catalogue/MultifLandscapesA1706/1.Data/Results/SDM results"

species_selected <- as.data.frame(
  readxl::read_xlsx(
    "//catalogue/MultifLandscapesA1706/1.Data/Results/species_lists/Species_20250304_edible_part_curated.xlsx",
    col_names = T
  )
)

#reading and getting conversion table to use all files togethers
conversion_table <-
  as.data.frame(
    readxl::read_xlsx(
    "E:/CSOSA/Dropbox/VAVILOV_2.0/Results/SUMMARY_FILES/PROJECT_STATUS.xlsx",
    "Hoja1"
  )
)

#template
#load template for rasterize
output_dir <- "//catalogue/MultifLandscapesA1706/1.Data/Results"
changes_map_conv_hull(basedir,output_dir,resultdir,species_selected,conversion_table,template_path,redo_maps)




