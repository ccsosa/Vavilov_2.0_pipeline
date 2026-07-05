### this script can be used to obtain change maps, current, and climate change maps for each species!


# author: Chrystian Sosa

library(terra)
library(viridis)
library(sf)
library(foreach)
library(bit)
# library(ggplot2)
library(bit)
library(tidyterra)
require(leaflet)
library(htmlwidgets)
library(viridisLite)
library(ggplot2)


areas_per_biome <-
  function(resultdir,
           species_selected,
           conversion_table,
           output_dir) {
    # chunks summary_dir
    output_dir_maps  <-
      paste0(output_dir, "/Summary_maps_folder_CC")
    if (!dir.exists(output_dir_maps)) {
      dir.create(output_dir_maps)
    }
    
    #change folder calculation (potential)
    output_dir_maps_AREAS_potential  <-
      paste0(output_dir_maps, "/changes_SP_CC_BIOME")
    if (!dir.exists(output_dir_maps_AREAS_potential)) {
      dir.create(output_dir_maps_AREAS_potential)
    }
    
    ##############################################################################
    ##############################################################################
    ##############################################################################
    #loading biome_shp
    message("loading shapefile and metrics")
    biome_shp <-
      terra::vect(paste0(basedir, "/1.Data/RAW/input_data/Ecoregions/biomes_LATAM_2017", ".shp"))
    #all species
    #adding presence absence results folder
    InDir <-
      paste0(resultdir, "/Distribution maps/Presence-absence")
    InDir_buffer <-
      "//catalogue/MultifLandscapesA1706/1.Data/Results/buffer_richness/raster"
    ##############################################################################
    #getting all species
    species_total <- unique(species_selected$species)
    ##############################################################################
    #getting all unique species
    species_total <- unique(species_selected$species)
    ##############################################################################
    #species with buffer
    sp_buffer_list <- list.files(InDir_buffer, pattern = ".tif")
    #getting specices names for buffer
    sp_buffer_list_name <-
      sub(pattern = ".tif", replacement = "", sp_buffer_list)
    ##############################################################################
    #getting all possible species with results
    sp_modeled_availa <- list.files(InDir, pattern = ".tif")
    sp_modeled_availa <- sub(".tif", "", sp_modeled_availa)
    ##############################################################################
    #load metrics to get modeled species
    metrics <-
      read.csv(paste0(resultdir, "/", "AUC_Maxent_valid_all.csv"))
    #gettting species lists to get results (modeled species with MaxEnt)
    species_list <- metrics$species
    species_list <-
      sub(pattern = "model_evaluation_", replacement = "", species_list)
    species_list <-
      sub(pattern = ".csv", replacement = "", species_list)
    # print(length(species_list))
    # species_list <- species_list[which(species_list!="Cabralea canjerana")]
    # print(length(species_list))
    species_list_df <- data.frame(
      sp = species_total,
      syn = NA,
      sdm = NA,
      buffer = NA,
      sdm_avail_name = NA,
      buffer_avail_name = NA
    )
    ################################################################################
    #species with valid model in the folder
    avail_sp_modeled <-
      species_list[species_list %in% sp_modeled_availa]
    ################################################################################
    #fixing names
    #geting species names accepted
    accp_sp <-
      conversion_table[conversion_table$species %in% avail_sp_modeled, ]
    
    accp_sp$avail_sp <- accp_sp$species
    #getting synonym!
    syn_sp <-
      conversion_table[conversion_table$species_from_source %in% avail_sp_modeled, ]
    #only removing accepting to get the synonyms and accepted names
    syn_sp <- syn_sp[which(syn_sp$taxonomic_status != "Accepted"), ]
    syn_sp$avail_sp <- NA
    for (i in 1:nrow(syn_sp)) {
      # i <- 1
      if (isTRUE(syn_sp$species[[i]] %in% sp_modeled_availa)) {
        syn_sp$avail_sp[[i]] <- syn_sp$species[[i]]
      }  else if (isTRUE(syn_sp$species_from_source[[i]] %in% sp_modeled_availa)) {
        syn_sp$avail_sp[[i]] <- syn_sp$species_from_source[[i]]
      }  else {
        syn_sp$avail_sp[[i]] <- NA
      }
      # syn_sp$species_from_source[[i]] %in% sp_modeled_availa
    }
    #getting available
    avail_sp_modeled_df <- rbind(accp_sp, syn_sp)
    #only getting unique species
    avail_sp_modeled <- unique(avail_sp_modeled_df$avail_sp)
    ##############################################################################
    #species with buffer
    # buffer_sps <- species_total[sp_buffer_list_name %in% species_total]
    #species name that are not modelled!
    buffer_sps <-
      species_total[!species_total %in% avail_sp_modeled]
    #if buffer_sps + avail_sp_modeled =8221 is well done!
    
    #available in the buffer results folder
    buffer_sps <-
      buffer_sps[buffer_sps %in% sp_buffer_list_name] #available to use buffer
    #
    #comparing with conversion table
    buff_sp <-
      conversion_table[conversion_table$species %in% buffer_sps, ]
    #accepted
    buff_sp_acc <- buff_sp[buff_sp$taxonomic_status == "Accepted", ]
    buff_sp_acc$avail_sp <- buff_sp_acc$species
    #rescuing the not accepted (synonyms, invalid, etc...)
    buff_sp_syn <- buff_sp[buff_sp$taxonomic_status != "Accepted", ]
    buff_sp_acc <-
      buff_sp_acc[buff_sp_acc$species %in% sp_buffer_list_name, ]
    buff_sp_syn <-
      buff_sp_syn[buff_sp_syn$species_from_source %in% sp_buffer_list_name, ]
    buff_sp_syn$avail_sp <- NA
    for (i in 1:nrow(buff_sp_syn)) {
      # i <- 1
      if (isTRUE(buff_sp_syn$species[[i]] %in% sp_buffer_list_name)) {
        buff_sp_syn$avail_sp[[i]] <- buff_sp_syn$species[[i]]
      }  else if (isTRUE(buff_sp_syn$species_from_source[[i]] %in% sp_buffer_list_name)) {
        buff_sp_syn$avail_sp[[i]] <- buff_sp_syn$species_from_source[[i]]
      }  else {
        buff_sp_syn$avail_sp[[i]] <- NA
      }
    }
    buffer_sps_df <- rbind(buff_sp_acc, buff_sp_syn)
    buffer_sps_df <-
      buffer_sps_df[!buffer_sps_df$sp_dummy_id %in% avail_sp_modeled_df$sp_dummy_id, ]
    buffer_sps <- unique(buffer_sps_df$avail_sp)
    
    ##############################################################################
    ##############################################################################
    ##############################################################################
    ##############################################################################
    ################################################################################
    #assigning hul
    hull <- c("ConcaveHull")
    #assigning SSPs
    SSPs <- c("ssp245", "ssp370")
    # if(!file.exists( paste0(output_dir_maps_TIF,"/","species_richness_current.tif"))){
    
    #this block is to get SSPs folder names
    for (i in 1:2) {
      SSP_output_dir_maps_AREAS_ConcaveHull <-
        paste0(output_dir_maps_AREAS_potential, "/",
               hull[[1]], "/", SSPs[[i]])
      if (!dir.exists(SSP_output_dir_maps_AREAS_ConcaveHull)) {
        dir.create(SSP_output_dir_maps_AREAS_ConcaveHull, recursive = T)
      }
      # SSP_output_dir_maps_AREAS_ConvHull <-
      #   paste0(output_dir_maps_AREAS_potential, "/",
      #          hull[[2]], "/", SSPs[[i]])
      # if (!dir.exists(SSP_output_dir_maps_AREAS_ConvHull)) {
      #   dir.create(SSP_output_dir_maps_AREAS_ConvHull, recursive = T)
      # }
      
    }
    
    ##############################################################################
    #input folder
    FutDir <-
      "//catalogue/MultifLandscapesA1706/1.Data/Results/Summary_maps_folder_CC/TIF_SP_CC_REALIZED"
    #This block performs the change procedure (FUTURE - CURRENT)
    #hull dirs
    
    
    message("doing summary for future")
    message("Reading all rasters from species list")
    # s <- 1
    #concave hull
    for (d in 1:1) {
      # d <- 1
      #SSPs
      for (s in 1:2) {
        # s <- 1
        #species
        for (i in seq_along(avail_sp_modeled)) {
          # i <- 2
          sp_name <- avail_sp_modeled[[i]]
          message(paste0(hull[[d]], " / ", SSPs[[s]], " / ", sp_name, " / ", i))
          out_path <-
            paste0(FutDir,
                   "/",
                   hull[[d]],
                   "/",
                   SSPs[[s]],
                   "/",
                   sp_name,
                   "_changes",
                   ".tif")
          
          area_j <- list()
          if (file.exists(out_path)) {
            
            #reading raster
            x <- terra::rast(out_path)
            message(paste0("starting ", hull[[d]], " / ", SSPs[[s]], " / ", sp_name, " / ", i))
            for (j in 1:nrow(biome_shp)) {
              
              x_j <- terra::mask(x, biome_shp[j, ])
              x_j <- terra::crop(x_j, biome_shp[j, ])
              x_Area <- terra::expanse(x_j,
                                       unit = "km",
                                       byValue = T,
                                       wide = T)
              #obtaining areas in a CSV file
              area_df <- data.frame(
                sp = sp_name,
                hull = hull[[d]],
                SSP = SSPs[[s]],
                Biome = biome_shp[j, ]$BIOME,
                loss = ifelse("-1" %in% names(x_Area), x_Area$`-1`, 0),
                stable = ifelse("1" %in% names(x_Area), x_Area$`1`, 0),
                new_areas = ifelse("2" %in% names(x_Area), x_Area$`2`, 0)
              )
              area_j[[j]] <- area_df
              rm(x_j, x_Area,area_df)
            }
            
            rm(x);gc() #j
            area_j <- do.call(rbind, area_j)
            ###########################
            #get areas for changes
            
            message(
              paste0(
                "SAVING:",
                output_dir_maps_AREAS_potential,
                "/",
                hull[[d]],
                "/",
                SSPs[[s]],
                "/",
                sp_name,
                "_areas",
                ".csv"
              )
            )
            
            write.csv(
              area_j,
              paste0(
                output_dir_maps_AREAS_potential,
                "/",
                hull[[d]],
                "/",
                SSPs[[s]],
                "/",
                sp_name,
                "_areas",
                ".csv"
              ),
              row.names = F
            )
            rm(area_j)
          }
        }#species
      } #ssp
    } #hull
    return("DONE!")
  }
################################################################################
user <- "Chrystian"
if (user == "Chrystian") {
  basedir <-  "//catalogue/MultifLandscapesA1706"
}

resultdir <-
  "//catalogue/MultifLandscapesA1706/1.Data/Results/SDM results"

species_selected <- as.data.frame(
  readxl::read_xlsx(
    "//catalogue/MultifLandscapesA1706/1.Data/Results/species_lists/Species_20250304_edible_part_curated.xlsx",
    col_names = T
  )
)

#reading and getting conversion table to use all files togethers
conversion_table <-
  readxl::read_xlsx(
    "E:/CSOSA/Dropbox/VAVILOV_2.0/Results/SUMMARY_FILES/PROJECT_STATUS.xlsx",
    "Hoja1"
  )

#config adding a results_dir
output_dir <- "//catalogue/MultifLandscapesA1706/1.Data/Results"



areas_per_biome(resultdir,
           species_selected,
           conversion_table,
           output_dir)
