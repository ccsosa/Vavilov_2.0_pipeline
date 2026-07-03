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



changes_map_conv_hull <- function(basedir,output_dir,resultdir,species_selected,conversion_table,template_path,redo_maps){
  # chunks summary_dir
  output_dir_maps  <-
    paste0(output_dir, "/Summary_maps_folder_CC")
  if (!dir.exists(output_dir_maps)) {
    dir.create(output_dir_maps)
  }
  
  #PNG folder (MAP)
  output_dir_maps_PNG  <-
    paste0(output_dir_maps, "/PNG_SP_CC_REALIZED")
  if (!dir.exists(output_dir_maps_PNG)) {
    dir.create(output_dir_maps_PNG)
  }
  #TIF folder
  output_dir_maps_TIF  <-
    paste0(output_dir_maps, "/TIF_SP_CC_REALIZED")
  if (!dir.exists(output_dir_maps_TIF)) {
    dir.create(output_dir_maps_TIF)
  }
  
  #change folder calculation (potential)
  output_dir_maps_AREAS_realized  <-
    paste0(output_dir_maps, "/changes_SP_CC_REALIZED")
  if (!dir.exists(output_dir_maps_AREAS_realized)) {
    dir.create(output_dir_maps_AREAS_realized)
  }
  #chunks to use transitorelly 
  output_dir_chunks  <-
    paste0(output_dir_maps, "/chunks_p_sp_CC")
  if (!dir.exists(output_dir_chunks)) {
    dir.create(output_dir_chunks)
  }  
  #chunks to use transitorelly for concave
  output_dir_TIF_CH  <-
    paste0(output_dir_maps, "/TIF_CH_sp_CC")
  if (!dir.exists(output_dir_TIF_CH)) {
    dir.create(output_dir_TIF_CH)
  }  
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  #loading adm0
  message("loading shapefile and metrics")
  adm0_path <- paste0(basedir, "/1.Data/RAW/input_data/adm0/adm0_Latam", ".shp")
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
  species_list_df <- data.frame(sp=species_total,
                                syn=NA,
                                sdm=NA,
                                buffer=NA,
                                sdm_avail_name=NA,
                                buffer_avail_name=NA)
  ################################################################################
  #species with valid model in the folder
  avail_sp_modeled <-
    species_list[species_list %in% sp_modeled_availa]
  ################################################################################
  #fixing names
  #geting species names accepted
  accp_sp <-
    conversion_table[conversion_table$species %in% avail_sp_modeled,]
  
  accp_sp$avail_sp <- accp_sp$species
  #getting synonym!
  syn_sp <-
    conversion_table[conversion_table$species_from_source %in% avail_sp_modeled,]
  #only removing accepting to get the synonyms and accepted names
  syn_sp <- syn_sp[which(syn_sp$taxonomic_status != "Accepted"),]
  syn_sp$avail_sp <- NA
  for(i in 1:nrow(syn_sp)){
    # i <- 1
    if(isTRUE(syn_sp$species[[i]] %in% sp_modeled_availa)){
      syn_sp$avail_sp[[i]] <- syn_sp$species[[i]]
    }  else if(isTRUE(syn_sp$species_from_source[[i]] %in% sp_modeled_availa)){
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
    conversion_table[conversion_table$species %in% buffer_sps,]
  #accepted
  buff_sp_acc <- buff_sp[buff_sp$taxonomic_status=="Accepted",]
  buff_sp_acc$avail_sp <- buff_sp_acc$species
  #rescuing the not accepted (synonyms, invalid, etc...)
  buff_sp_syn <- buff_sp[buff_sp$taxonomic_status!="Accepted",]
  buff_sp_acc <- buff_sp_acc[buff_sp_acc$species %in% sp_buffer_list_name,]
  buff_sp_syn <- buff_sp_syn[buff_sp_syn$species_from_source %in% sp_buffer_list_name,]
  buff_sp_syn$avail_sp <- NA
  for(i in 1:nrow(buff_sp_syn)){
    # i <- 1
    if(isTRUE(buff_sp_syn$species[[i]] %in% sp_buffer_list_name)){
      buff_sp_syn$avail_sp[[i]] <- buff_sp_syn$species[[i]]
    }  else if(isTRUE(buff_sp_syn$species_from_source[[i]] %in% sp_buffer_list_name)){
      buff_sp_syn$avail_sp[[i]] <- buff_sp_syn$species_from_source[[i]]
    }  else {
      buff_sp_syn$avail_sp[[i]] <- NA
    }
  }
  buffer_sps_df <- rbind(buff_sp_acc, buff_sp_syn)
  buffer_sps_df <- buffer_sps_df[!buffer_sps_df$sp_dummy_id %in% avail_sp_modeled_df$sp_dummy_id,]
  buffer_sps <- unique(buffer_sps_df$avail_sp)
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ################################################################################
  #assigning SSPs
  SSPs <- c("ssp245","ssp370")
  
  #concave or convex hull folders
  conv_hull_dirs <- c("concave_hull_rasters_ADM0_2_5m","conv_hull_rasters_ADM0_2_5m")
  #concave or convex hull names
  conv_names <- c("ConcaveHull","ConvexHull")
  # if(!file.exists( paste0(output_dir_maps_TIF,"/","species_richness_current.tif"))){
  
  #this block is to get SSPs folder names
  for(i in 1:2){
    #####ADDING CONVEX or CONC
    SSP_dir <- paste0(output_dir_maps_TIF,"/",conv_names[[i]])
    if(!dir.exists(SSP_dir)){
      dir.create(SSP_dir)
    }
    SSP_dir_PNG <- paste0(output_dir_maps_PNG,"/",conv_names[[i]])
    if(!dir.exists(SSP_dir_PNG)){
      dir.create(SSP_dir_PNG)
    }
    SSP_dir_PNG <- paste0(output_dir_maps_PNG,"/",conv_names[[i]])
    if(!dir.exists(SSP_dir_PNG)){
      dir.create(SSP_dir_PNG)
    }
    SSP_dir_Areas_pot <- paste0(output_dir_maps_AREAS_realized,"/",conv_names[[i]])
    if(!dir.exists(SSP_dir_Areas_pot)){
      dir.create(SSP_dir_Areas_pot)
    }
    
    for(j in 1:2){
      SSP_dir <- paste0(output_dir_maps_TIF,"/",conv_names[[i]],"/",SSPs[[j]])
      if(!dir.exists(SSP_dir)){
        dir.create(SSP_dir)
      }
      SSP_dir_PNG <- paste0(output_dir_maps_PNG,"/",conv_names[[i]],"/",SSPs[[j]])
      if(!dir.exists(SSP_dir_PNG)){
        dir.create(SSP_dir_PNG)
      }
      SSP_dir_PNG <- paste0(output_dir_maps_PNG,"/",conv_names[[i]],"/",SSPs[[j]])
      if(!dir.exists(SSP_dir_PNG)){
        dir.create(SSP_dir_PNG)
      }
      
      SSP_dir_Areas_pot <- paste0(output_dir_maps_AREAS_realized,"/",conv_names[[i]],"/",SSPs[[j]])
      if(!dir.exists(SSP_dir_Areas_pot)){
        dir.create(SSP_dir_Areas_pot)
      }
      output_dir_maps_AREAS_realized_s <- paste0(output_dir_maps_AREAS_realized,"/",conv_names[[i]],"/",SSPs[[j]])
      if(!dir.exists(output_dir_maps_AREAS_realized_s)){
        dir.create(output_dir_maps_AREAS_realized_s)
        
      }
    }
  }
  
  ##############################################################################
  #SSPs_dirs
  message("doing summary for future")
  message("Reading all rasters from species list")
  
  inDir <- paste0(basedir,"/","1.Data/Results/Summary_maps_folder_CC/TIF_SP_CC")
  #c <- 1
  for(d in 1:2){
    #defining convex or concave hull folder and names
    # c <- 1
    conv_hull_dir_c <-paste0(basedir,"/1.Data/Results/",conv_hull_dirs[[d]])
    message(conv_names[[d]])
  

    
    for(s in 1:2){
      message(paste0(conv_names[[d]]," / ",SSPs[[s]]))
      
      # s <- 1
      cl <- parallel::makeCluster(4)
      doParallel::registerDoParallel(cl)
      foreach::foreach(i = seq_along(avail_sp_modeled),
                       .packages = c("terra","sf"),
                       .export = c("avail_sp_modeled","conv_hull_dir_c","conv_names",
                                   "SSPs","inDir","output_dir_maps_TIF","adm0_path",
                                   "output_dir_maps_PNG","output_dir_maps_AREAS_realized",
                                   "conversion_table","template_path","adm0_path","redo_maps",
                                   "d","s")) %dopar% {
                                     #loading adm0
                                     adm0 <- terra::vect(adm0_path)
                                     #template
                                     template <- terra::rast(template_path)
                                     template[which(!is.na(template[]))] <- 1
                                     sp_name <- avail_sp_modeled[[i]]
                                     message(paste0(conv_names[[d]]," / ",SSPs[[s]]," / ",sp_name," / ",i))
                                     # out_path <- paste0(FutDir, "/", sp_name,"_",SSPs[[s]],".tif")
                                     
                                     #getting changes file
                                     out_path <- paste0(inDir,"/",SSPs[[s]],"/",sp_name,"_changes",".tif")
                                     final_file_path <- paste0(output_dir_maps_TIF,"/",
                                                               conv_names[[d]],"/",SSPs[[s]],"/",sp_name,"_changes",".tif")
                                     if(!file.exists(final_file_path)){
                                       if(file.exists(out_path)){
                                         change <- terra::rast(out_path)
                                         
                                         #loading C hull
                                         if(d==1){
                                           ch_path <- paste0(conv_hull_dir_c, "/", sp_name, "_conc.tif")
                                         } else if(d==2){
                                           ch_path <- paste0(conv_hull_dir_c, "/", sp_name, "_conv.tif")
                                         }
                                         
                                         if(file.exists(ch_path)){
                                           CH_sp  <- terra::rast(ch_path)
                                           CH_sp <- terra::resample(CH_sp,template,method="near")
                                           
                                           result <- change * CH_sp
                                           terra::writeRaster(result, filename = final_file_path, overwrite = TRUE)
                                           
                                         } else {
                                           sp2 <- conversion_table$species[conversion_table$species_from_source==sp_name]
                                           ch_path2 <- paste0(conv_hull_dir_c, "/", sp2, "_conc.tif")
                                           if(file.exists(ch_path2)){
                                             CH_sp  <- rast(ch_path2)
                                             CH_sp <- terra::resample(CH_sp,template,method="near")
                                             result <- change * CH_sp
                                             terra::writeRaster(result, filename = final_file_path, overwrite = TRUE)
                                             
                                           } else {
                                             print(i)
                                             message(sp_name)
                                             message("NO CONCAVE HULL")
                                             #result <- SDM_sp
                                             result <- NULL
                                           }
                                         }
                                       }
                                       #loading changes
                                       
                                       
                                       
                                       
                                       
                                       if(!is.null(result)){
                                         #get areas for changes
                                         x_Area <- terra::expanse(result,unit="km",byValue=T,wide=T)
                                         #obtaining areas in a CSV file
                                         area_df <- data.frame(sp = sp_name,
                                                               ch = conv_names[[d]],
                                                               SSP = SSPs[[s]],
                                                               loss= ifelse("-1" %in% names(x_Area),x_Area$`-1`, 0),
                                                               stable= ifelse("1" %in% names(x_Area),x_Area$`1`, 0),
                                                               new_areas= ifelse("2" %in% names(x_Area),x_Area$`2`, 0)
                                         )
                                         write.csv(area_df,paste0(output_dir_maps_AREAS_realized,"/",
                                                                  conv_names[[d]],"/",SSPs[[s]],"/",sp_name,
                                                                  "_areas",".csv"))
                                         
                                         ###########################
                                         
                                         
                                         #redo maps?
                                         redo_maps = redo_maps
                                         if(redo_maps==F){
                                           #plot
                                           #definitng colors
                                           change_vals <- c(-1,0,1,2)
                                           change_colrs <- c("#c0392b", #LOSS
                                                             "grey98",#no available
                                                             "#27ae60", #stable
                                                             "#f39c12"# new areas
                                           )
                                           
                                           change_labels <- c("Loss","Not available","Stable","New areas")
                                           setwd(dir = paste0(output_dir_maps_PNG,"/",conv_names[[d]],"/",SSPs[[s]]))
                                           
                                           png(paste0(sp_name,"_changes",".png"), width = 10, height = 10,
                                               units = "in",   # interpret as inches
                                               res = 1000)      # 300 dpi = print quality
                                           par(mfrow = c(1,1))
                                           par(mai = c(0.5, 0.5, 0.5, 0.5))
                                           plot(change,
                                                main="",
                                                legend=F,
                                                col = change_colrs)
                                           plot(sf::st_geometry(sf::st_as_sf(adm0)), 
                                                add = TRUE, 
                                                border = "grey17", 
                                                lwd = 1.25)
                                           legend("bottomleft",
                                                  legend = change_labels,
                                                  fill = change_colrs,
                                                  border = "grey45",
                                                  bty = "o",
                                                  bg="white",
                                                  inset=c(0.01,0.01),
                                                  xpd=NA,
                                                  cex = 0.8)
                                           dev.off()
                                         } else {
                                           "OK"
                                         }
                                       } else {
                                         #obtaining areas in a CSV file
                                         area_df <- data.frame(sp = sp_name,
                                                               ch = conv_names[[d]],
                                                               SSP = SSPs[[s]],
                                                               loss= NA,
                                                               stable= NA,
                                                               new_areas= NA)
                                         
                                         write.csv(area_df,paste0(output_dir_maps_AREAS_realized,"/",
                                                                  conv_names[[d]],"/",SSPs[[s]],"/",sp_name,
                                                                  "_areas",".csv"))
                                       }
                                     }
                                   } #file (i)
      parallel::stopCluster(cl)
    } # SSP (s)
  } #concave or convex hull (d)
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
template_path  <-
  paste0(
    "//catalogue/MultifLandscapesA1706/1.Data/RAW/Input_data/climate_data/2_5min/present/",
    "wc2.1_2.5m_bio_1.tif"
  )

redo_maps <- F
# rast()

# template[which(!is.na(template[]))] <- 1
#config adding a results_dir
output_dir <- "//catalogue/MultifLandscapesA1706/1.Data/Results"
changes_map_conv_hull(basedir,output_dir,resultdir,species_selected,conversion_table,template_path,redo_maps)