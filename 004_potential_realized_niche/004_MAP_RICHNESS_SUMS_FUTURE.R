### this script can be used to summary maps for potential distribution and obtain realized niche summary maps


# author: Chrystian Sosa

library(terra)
library(viridis)
library(sf)
library(foreach)
library(bit)
# library(ggplot2)
# library(bit)
require(leaflet)
library(htmlwidgets)
library(viridisLite)

# set target regionM
summary_richness_function <- function(basedir,resultdir,template,species_selected,num_classes,conversion_table){
  ##############################################################################
  ##############################################################################
  ##############################################################################
  #loading adm0
  message("loading shapefile and metrics")
  adm0_path <- paste0(basedir, "/1.Data/RAW/input_data/adm0/adm0_Latam", ".shp")
  #all species
  #config adding a results_dir
  output_dir <- "//catalogue/MultifLandscapesA1706/1.Data/Results"
  # chunks summary_dir
  output_dir_maps  <-
    paste0(output_dir, "/Summary_maps_folder_CC")
  if (!dir.exists(output_dir_maps)) {
    dir.create(output_dir_maps)
  }
  
  output_dir_maps_PNG  <-
    paste0(output_dir_maps, "/PNG")
  if (!dir.exists(output_dir_maps_PNG)) {
    dir.create(output_dir_maps_PNG)
  }
  output_dir_maps_TIF  <-
    paste0(output_dir_maps, "/TIF")
  if (!dir.exists(output_dir_maps_TIF)) {
    dir.create(output_dir_maps_TIF)
  }
  output_dir_maps_HTML  <-
    paste0(output_dir_maps, "/HTML")
  if (!dir.exists(output_dir_maps_HTML)) {
    dir.create(output_dir_maps_HTML)
  }  
  
  output_dir_chunks  <-
    paste0(output_dir_maps, "/chunks")
  if (!dir.exists(output_dir_chunks)) {
    dir.create(output_dir_chunks)
  }  
  #cutting by convex and concave
  output_dir_TIF_REALIZED  <-
    paste0(output_dir_maps, "/TIF_SP_REALIZED_SDM")
  if (!dir.exists(output_dir_TIF_REALIZED)) {
    dir.create(output_dir_TIF_REALIZED)
  }  
  output_dir_maps_PNG_REALIZED  <-
    paste0(output_dir_maps, "/PNG_SP_REALIZED_SDM")
  if (!dir.exists(output_dir_maps_PNG_REALIZED)) {
    dir.create(output_dir_maps_PNG_REALIZED)
  }
  ##############################################################################
  ##############################################################################
  ##############################################################################
  #loading adm0
  message("loading shapefile and metrics")
  adm0 <- vect(paste0(basedir, "/1.Data/RAW/input_data/adm0/adm0_Latam", ".shp"))
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

  ################################################################################
  
  ##############################################################################
  #adding presence absence results folder
  InDir <- paste0(resultdir,"/Distribution maps/Future/2050/Consensus_maps_folder2/Fut_consensus")
  ##############################################################################
  ##############################################################################
  ################################################################################
  ####First step (summary file) currrent 
  ##############################################################################
  #assigning SSPs
  SSPs <- c("ssp245","ssp370")

################################################################################
################################################################################
################################################################################
################################################################################
#potential niche

    message("doing summary for future")
    message("Reading all rasters from species list")
    # s <- 1
  for(s in 1:2){
  # if(!file.exists(paste0(output_dir_maps_TIF, "/", "species_richness_future_",SSPs[[s]],".tif"))) {
      
    if(!file.exists(paste0(output_dir_maps_TIF, "/", "species_richness_SDM_future_",SSPs[[s]],".tif"))) {
      
      message("Doing summary for current (parallel chunked mode)")
      
      # --- 1. Validate paths sequentially ---------------------------------------
      message("Checking file availability...")
      valid_paths <- vapply(avail_sp_modeled, function(sp) {
        path <- paste0(InDir, "/", sp,"_",SSPs[[s]],".tif")
        if (file.exists(path)) path else NA_character_
      }, character(1))
      
      missing_sps <- avail_sp_modeled[is.na(valid_paths)]
      valid_paths <- valid_paths[!is.na(valid_paths)]
      
      if (length(missing_sps) > 0)
      message(paste0(length(missing_sps), " species raster(s) not found and skipped."))
      message(paste0(length(valid_paths), " valid rasters found."))
      
      # --- 2. Split into chunks of 100 ------------------------------------------
      chunk_size <- 15
      chunk_idx  <- split(seq_along(valid_paths),
                          ceiling(seq_along(valid_paths) / chunk_size))
      n_chunks   <- length(chunk_idx)
      message(paste0("Processing ", n_chunks, " chunks of up to ", chunk_size, " rasters..."))
      
      # --- 3. Setup parallel backend --------------------------------------------
      n_cores <- 8 #min(parallel::detectCores() - 2, n_chunks)
      cl      <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      message(paste0("Using ", n_cores, " cores for ", n_chunks, " chunks"))
      
      # --- Create chunks subfolder ----------------------------------------------
      chunk_dir <- paste0(output_dir_maps_TIF, "/chunks_fut")
      if (!dir.exists(chunk_dir)) dir.create(chunk_dir)
      
      # --- 4. Parallel chunk summation with NA → 0 ------------------------------
      chunk_results <- foreach::foreach(
        chunk          = chunk_idx,
        .packages      = "terra",
        .errorhandling = "pass"
      ) %dopar% {
        
        paths_chunk <- valid_paths[chunk]
        chunk_id    <- chunk[1]  # use first index as unique chunk identifier
        
        stack  <- terra::rast(paths_chunk)
        result <- terra::app(stack, fun = "sum", na.rm = TRUE)
        
        # Save to known folder with predictable name
        out_path <- paste0(chunk_dir, "/chunk_", sprintf("%04d", chunk_id), ".tif")
        terra::writeRaster(result, filename = out_path, overwrite = TRUE)
        out_path
      }
      
      parallel::stopCluster(cl)
      
      # --- 5. Check for chunk errors --------------------------------------------
      failed <- sapply(chunk_results, inherits, "error")
      if (any(failed)) {
        message(paste0(sum(failed), " chunk(s) failed:"))
        print(chunk_results[failed])
        stop("Aborting: fix failed chunks before proceeding.")
      }
      
      # --- 6. Final sum across chunk files --------------------------------------
      message("Combining chunk results into final richness raster...")
      chunk_paths <- unlist(chunk_results)
      chunk_stack <- terra::rast(chunk_paths)
      x           <- terra::app(chunk_stack, fun = "sum", na.rm = TRUE)
      
      # --- 7. Save final raster -------------------------------------------------
      message("Saving species_richness_future.tif ...")
      terra::writeRaster(
        x,
        filename  = paste0(output_dir_maps_TIF, "/", "species_richness_SDM_future_",SSPs[[s]],".tif"),
        overwrite = TRUE
      )
      
      # --- 8. Cleanup chunks folder ---------------------------------------------
      message("Cleaning up chunk files...")
      file.remove(chunk_paths)
      unlink(chunk_dir)
      richness <- x
      message("Done!")
    } else {
      message("loading available raster")
      richness <- terra::rast(paste0(output_dir_maps_TIF, "/",  "species_richness_SDM_future_",SSPs[[s]],".tif"))
      
    }
    ################################################################################
    #plotting current
    
    message("Making a plot")
    
    #Plotting
    par(mar = c(1, 1, 1, 1))
    
    setwd(dir = paste0(output_dir_maps_PNG))
    png(paste0("LATAM_Edible_food_","_species_richness_future_",SSPs[[s]], ".png"), width = 10, height = 10,
        units = "in",   # interpret as inches
        res = 1000)      # 300 dpi = print quality
    par(mfrow = c(1,1))
    par(mai = c(0.5, 0.5, 0.5, 0.5))
    plot(richness,
         main="",
         #main = paste0("Species richness"),
         legend = T,
         col = viridis(length(species_list)))
    plot(sf::st_geometry(sf::st_as_sf(adm0)), add = TRUE, border = "grey17", lwd = 1.25)
    dev.off()
    
  }
  ################################################################################
  ################################################################################
  ################################################################################
  ####Third step (summary file in concave and convex  hull!)
  ################################################################################
    
    # if(!file.exists( paste0(output_dir_maps_TIF,"/","species_richness_current.tif"))){
    #concave or convex hull folders
    conv_hull_dirs <- c("concave_hull_rasters_ADM0_2_5m","conv_hull_rasters_ADM0_2_5m")
    #concave or convex hull names
    conv_names <- c("ConcaveHull","ConvexHull")
    # if(!file.exists( paste0(output_dir_maps_TIF,"/","species_richness_current.tif"))){
    
    #this block is to get SSPs folder names
    for(i in 1:2){
      #####ADDING CONVEX or CONC folder and SSPs
      SSP_dir <- paste0(output_dir_TIF_REALIZED,"/",conv_names[[i]])
      if(!dir.exists(SSP_dir)){
        dir.create(SSP_dir)
      }
      SSP_dir_PNG <- paste0(output_dir_maps_PNG_REALIZED,"/",conv_names[[i]])
      if(!dir.exists(SSP_dir_PNG)){
        dir.create(SSP_dir_PNG)
      }
      
      for(j in 1:2){
        SSP_dir <- paste0(output_dir_TIF_REALIZED,"/",conv_names[[i]],"/",SSPs[[j]])
        if(!dir.exists(SSP_dir)){
          dir.create(SSP_dir)
        }
        SSP_dir_PNG <- paste0(output_dir_maps_PNG_REALIZED,"/",conv_names[[i]],"/",SSPs[[j]])
        if(!dir.exists(SSP_dir_PNG)){
          dir.create(SSP_dir_PNG)
        }
        
      }
    }
}
  #getting SDM * Concave hulls
    #assigning SSPs
    SSPs <- c("ssp245","ssp370")
    
    #concave or convex hull folders
    conv_hull_dirs <- c("concave_hull_rasters_ADM0_2_5m","conv_hull_rasters_ADM0_2_5m")
    #concave or convex hull names
    conv_names <- c("ConcaveHull","ConvexHull")
    
    ############################################################################
    #adding presence absence results folder
    InDir <- paste0(resultdir,"/Distribution maps/Future/2050/Consensus_maps_folder2/Fut_consensus")
    #c <- 1
    for(d in 1:2){
      #defining convex or concave hull folder and names
      # d <- 1
      conv_hull_dir_c <-paste0(basedir,"/1.Data/Results/",conv_hull_dirs[[d]])
      message(conv_names[[d]])
      
      for(s in 1:2){
        message(paste0(conv_names[[d]]," / ",SSPs[[s]]))
        
        # s <- 1
        cl <- parallel::makeCluster(2)
        doParallel::registerDoParallel(cl)
        foreach::foreach(i = seq_along(avail_sp_modeled),
                         .packages = c("terra","sf"),
                         .export = c("avail_sp_modeled","conv_hull_dir_c","conv_names",
                                     "SSPs","InDir","output_dir_maps_TIF","adm0_path",
                                     "output_dir_TIF_REALIZED",
                                     "conversion_table","template_path","adm0_path",
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
                                       out_path <- paste0(InDir,"/","/",sp_name,"_",SSPs[[s]],".tif")
                                       final_file_path <- paste0(output_dir_TIF_REALIZED,"/",
                                                                 conv_names[[d]],"/",SSPs[[s]],"/",sp_name,".tif")
                                       if(!file.exists(final_file_path)){
                                         if(file.exists(out_path)){
                                           SDM_future <- terra::rast(out_path)
                                           
                                           #loading C hull
                                           if(d==1){
                                             ch_path <- paste0(conv_hull_dir_c, "/", sp_name, "_conc.tif")
                                           } else if(d==2){
                                             ch_path <- paste0(conv_hull_dir_c, "/", sp_name, "_conv.tif")
                                           }
                                           
                                           if(file.exists(ch_path)){
                                             CH_sp  <- terra::rast(ch_path)
                                             CH_sp <- terra::resample(CH_sp,template,method="near")
                                             
                                             result <- SDM_future * CH_sp
                                             terra::writeRaster(result, filename = final_file_path, overwrite = TRUE)
                                             
                                           } else {
                                             sp2 <- conversion_table$species[conversion_table$species_from_source==sp_name]
                                             ch_path2 <- paste0(conv_hull_dir_c, "/", sp2, "_conc.tif")
                                             if(file.exists(ch_path2)){
                                               CH_sp  <- rast(ch_path2)
                                               CH_sp <- terra::resample(CH_sp,template,method="near")
                                               result <- SDM_future * CH_sp
                                               terra::writeRaster(result, filename = final_file_path, overwrite = TRUE)
                                               
                                             } else {
                                               print(i)
                                               message(sp_name)
                                               message("NO C HULL")
                                               #result <- SDM_sp
                                               result <- NULL
                                             }
                                           }
                                         }
                                         #loading changes
                                         

                                       }
                                     } #file (i)
        parallel::stopCluster(cl)
      } # SSP (s)
    } #concave or convex hull (d)   
    

  
  
  for(s in 1:2){
    
    out_paths  <- paste0(output_dir_TIF_CH_FUT, "/", avail_sp_modeled,"_",SSPs[[s]], ".tif")
    chunk_size <- 100                           # tune: lower = less RAM, more I/O passes
    
    chunks <- split(out_paths, ceiling(seq_along(out_paths) / chunk_size))
    
    chunk_stack_list <- list()
    for(k in seq_along(chunks)){
      message("Chunk ", k, " / ", length(chunks))
      chunk_test <- data.frame(chunk=chunks[[k]],
                               sp=chunks[[k]],
                               syn=NA,
                               sp_avail=NA,
                               syn_avail=NA)
      
      chunk_test$sp <- sub(paste0(output_dir_TIF_CH,"/"),"",chunk_test$sp)
      chunk_test$sp <- sub(".tif","",chunk_test$sp)
      for(l in 1:nrow(chunk_test)){
        # print(l)
        syn_l <-  conversion_table$species_from_source[which(conversion_table$species==chunk_test$sp[[l]])]
        if(length(syn_l >0)){
          chunk_test$syn[[l]] <- syn_l
        } else {
          chunk_test$syn[[l]] <- NA
        }
        
        
        if(file.exists(chunk_test$chunk[[l]])){
          chunk_test$chunk[[l]] <-chunk_test$chunk[[l]] 
        } else if(!file.exists(chunk_test$chunk[[l]])){
          if(file.exists(paste0(output_dir_TIF_CH,"/",chunk_test$syn[[l]],".tif"))){
            chunk_test$chunk[[l]] <- paste0(output_dir_TIF_CH,"/",chunk_test$syn[[l]],".tif")
          } else {
            chunk_test$chunk[[l]] <- NA
          }
        } else {
          chunk_test$chunk[[l]] <- NA
        }
      }
      chunks[[k]] <- chunk_test$chunk[which(!is.na(chunk_test$chunk))]
      
      
    chunk_stack <- rast(chunks[[k]])          # load chunk_size layers
     chunk_sum   <- app(chunk_stack, fun = "sum", na.rm = TRUE)
    chunk_stack_list[[k]] <- chunk_sum
    #   # sum_raster  <- sum_raster + chunk_sum
    #   # rm(chunk_stack, chunk_sum); gc()
    }
    chunk_stack_list <- do.call(c,chunk_stack_list)
    chunk_stack_list   <- app(chunk_stack_list, fun = "sum", na.rm = TRUE)
    
    # --- Step 3: write final richness raster for SDM * CH---
    terra::writeRaster(
      chunk_stack_list,
      filename  = paste0(output_dir_maps_TIF, "/species_richness_CH_SDM_future","_",SSPs[[s]],".tif"),
      overwrite = TRUE
    )
    }
  
  
  ##############################################################################
  ################################################################################
  # #for buffer (creatying buffer * CH)
  # # --- Step 1: build the masked rasters (write to disk, don't keep in RAM) ---
  # for(i in seq_along(buffer_sps)){
  #   sp_name <- buffer_sps[[i]]
  #   message(sp_name)
  #   out_path <- paste0(output_dir_TIF_CH_FUT, "/", sp_name, ".tif")
  #   
  #   if(!file.exists(out_path)){
  #     SDM_sp <- rast(paste0(InDir_buffer, "/", sp_name, ".tif"))
  #     SDM_sp <- terra::resample(SDM_sp,template2,method="near")
  #     
  #     ch_path <- paste0(conv_hull_dir, "/", sp_name, "_conc.tif")
  #     
  #     if(file.exists(ch_path)){
  #       CH_sp  <- rast(ch_path)
  #       CH_sp <- terra::resample(CH_sp,template2,method="near")
  #       result <- SDM_sp * CH_sp
  #       } else {
  #       message("NO CONCAVE HULL")
  #       result <- SDM_sp
  #     }
  #     
  #     terra::writeRaster(result, filename = out_path, overwrite = TRUE)
  #   }
  # }
  ################################################################################
  # #creating buffers sum chunks
  # 
  # # template <- template2
  # # --- Step 2: chunked sum - never hold more than chunk_size rasters in RAM ---
  # out_paths  <- paste0(output_dir_TIF_CH, "/", buffer_sps, ".tif")
  # # template   <- rast(out_paths[1])
  #  # sum_raster <- setValues(template, 0)        # accumulator, values in RAM but single layer
  # # rm(template); gc()
  # 
  # chunk_size <- 100                           # tune: lower = less RAM, more I/O passes
  # 
  # chunks <- split(out_paths, ceiling(seq_along(out_paths) / chunk_size))
  # chunks_list <- list()
  # for(k in seq_along(chunks)){
  #   message("Chunk ", k, " / ", length(chunks))
  #   chunk_stack <- rast(chunks[[k]])          # load chunk_size layers
  #   chunk_sum   <- sum(chunk_stack, na.rm = TRUE)
  #   chunks_list[[k]] <- chunk_sum
  #   # print(chunk_sum)
  #   # sum_raster  <- sum_raster + chunk_sum
  #   
  #   rm(chunk_stack, chunk_sum); gc()
  # }
  # chunks_list <- do.call(c,chunks_list)
  # chunks_list <- sum(chunks_list,na.rm = T)
  # 
  # ################################################################################
  # #writing buffer summary
  # 
  # # --- Step 3: write final richness raster ---
  # terra::writeRaster(
  #   chunks_list,
  #   filename  = paste0(output_dir_maps_TIF, "/species_richness_CH_BUFFER_future.tif"),
  #   overwrite = TRUE
  # )
  ################################################################################
  ################################################################################
  #getting final summary maps for CH
  
  # CH_SDM <- rast(paste0(output_dir_maps_TIF, "/species_richness_CH_SDM_future.tif"))
  # # CH_BUF <- rast(paste0(output_dir_maps_TIF, "/species_richness_CH_BUFFER_future.tif"))
  # CH_RICH_STACK <- c(CH_SDM,CH_BUF)
  # CH_RICHNESS <- sum(CH_RICH_STACK,na.rm = T)
  # 
  # # --- Step 3: write final richness raster ---
  # terra::writeRaster(
  #   CH_RICHNESS,
  #   filename  = paste0(output_dir_maps_TIF, "/species_richness_CH_future.tif"),
  #   overwrite = TRUE
  # )
  
  ################################################################################
  #plotting current CH*SDm+ CH*BUFFER maps
 for(s in 1:2){ 
   print(s)
   # s <- 1
   x  <- terra::rast(paste0(output_dir_maps_TIF, "/species_richness_CH_SDM_future","_",SSPs[[s]],".tif"))
   
  message("Making a plot")
  
  #Plotting
  par(mar = c(1, 1, 1, 1))
  
  setwd(dir = paste0(output_dir_maps_PNG))
  png(paste0("LATAM_Edible_food_","_species_richness_CH_future_",SSPs[[s]],".png"), width = 10, height = 10,
      units = "in",   # interpret as inches
      res = 1000)      # 300 dpi = print quality
  par(mfrow = c(1,1))
  par(mai = c(0.5, 0.5, 0.5, 0.5))
  plot(x,
       main="",
       #main = paste0("Species richness"),
       legend = T,
       col = viridis(length(species_list)))
  plot(sf::st_geometry(sf::st_as_sf(adm0)), add = TRUE, border = "grey17", lwd = 1.25)
  dev.off()

  ################################################################################
  ################################################################################
  ################################################################################
  ################################################################################
  ################################################################################
 #  message(paste("Step 5.9: HTML file step for richness future"))
 #  richness_total <- terra::rast(paste0(output_dir_maps_TIF, "/", "species_richness_SDM_future_",SSPs[[s]],".tif"))
 # 
 # # Get the range of raster values (excluding NA)
 #  raster_values <- values(richness_total, na.rm = TRUE)
 #  richness_total_rl <- raster::raster(richness_total)  # convert to RasterLayer
 #  
 #  # num_classes <- 8  # Adjust this number as needed
 #  
 #  # Generate breaks using quantiles (0 to 90th percentile, then add max)
 #  breaks <-
 #    quantile(
 #      raster_values,
 #      probs = seq(0, 0.9, length.out = num_classes + 1),
 #      na.rm = TRUE
 #    )
 #  x_b <- quantile(raster_values, probs = 1, na.rm = TRUE)
 #  breaks <- c(breaks, x_b)
 #  breaks <- as.numeric(round(breaks, 0))
 #  breaks <- unique(breaks)
 #  
 #  actual_num_classes <- length(breaks) - 1
 #  colors <- turbo(actual_num_classes)
 #  
 #  # ?????? Automatic label generation ??
 #  # For each bin [breaks[i], breaks[i+1]):
 #  #   . single-value bin  ??? "N"
 #  #   . last bin          ??? ">N"  (open upper end)
 #  #   . all others        ??? "N-M"
 #  
 #  
 #  auto_labels <- sapply(seq_len(actual_num_classes), function(i) {
 #    lo <- breaks[i]
 #    hi <- breaks[i + 1]
 #    if (i == actual_num_classes) {
 #      # last class: open upper bound
 #      paste0(">", lo)
 #    } else if (lo == hi - 1 || lo == hi) {
 #      # single-value class
 #      as.character(lo)
 #    } else {
 #      paste0(lo, "-", hi - 1)              # closed range shown as lo - (hi-1)
 #    }
 #  })
 #  
 #  
 #  
 #  # Create categorical colour palette
 #  pal <- colorBin(
 #    palette   = colors,
 #    domain    = values(richness_total),
 #    bins      = breaks,
 #    na.color  = "transparent"
 #  )
 #  # require(raster)
 #  x_l <- leaflet() %>%
 #    addTiles() %>%
 #    addProviderTiles(providers$CartoDB.Positron) %>%
 #    addMiniMap(width = 150, height = 150) %>%
 #    addRasterImage(
 #      richness_total,
 #      colors   = pal,
 #      opacity  = 0.3,
 #      project  = TRUE,
 #      maxBytes = 20 * 1024 ^ 2
 #    ) %>%
 #    addLegend(
 #      pal      = pal,
 #      values   = values(richness_total),
 #      title    = "Species Richness",
 #      labFormat = function(type, cuts, p)
 #        auto_labels  # fully automatic
 #    )
 # 
 #  htmlwidgets::saveWidget(
 #    x_l,
 #    paste0(
 #      output_dir_maps_HTML,
 #      "/",
 #      "species_richness_total_future_",SSPs[[s]],".html"
 #    ),
 #    selfcontained = T
 #  )
 #  
  
  
  
  message(paste("DONE: "))
  message(paste0("saved in ",output_dir_maps_HTML))
 
  
  ##############################################################################
  ################################################################################
  ################################################################################
  message(paste("Step 5.9: HTML file step for CH richness "))
  richness_total <- rast(paste0(output_dir_maps_TIF, "/species_richness_CH_SDM_future_",SSPs[[s]],".tif"))
  richness_total[richness_total[]==0] <- NA
  # Get the range of raster values (excluding NA)
  raster_values <- values(richness_total, na.rm = TRUE)
  richness_total_rl <- raster::raster(richness_total)  # convert to RasterLayer
  
  # num_classes <- 8  # Adjust this number as needed
  
  # Generate breaks using quantiles (0 to 90th percentile, then add max)
  breaks <-
    quantile(
      raster_values,
      probs = seq(0, 0.9, length.out = num_classes + 1),
      na.rm = TRUE
    )
  x_b <- quantile(raster_values, probs = 1, na.rm = TRUE)
  breaks <- c(breaks, x_b)
  breaks <- as.numeric(round(breaks, 0))
  breaks <- unique(breaks)
  
  actual_num_classes <- length(breaks) - 1
  colors <- turbo(actual_num_classes)
  
  # ?????? Automatic label generation ??
  # For each bin [breaks[i], breaks[i+1]):
  #   . single-value bin  ??? "N"
  #   . last bin          ??? ">N"  (open upper end)
  #   . all others        ??? "N-M"
  
  
  auto_labels <- sapply(seq_len(actual_num_classes), function(i) {
    lo <- breaks[i]
    hi <- breaks[i + 1]
    if (i == actual_num_classes) {
      # last class: open upper bound
      paste0(">", lo)
    } else if (lo == hi - 1 || lo == hi) {
      # single-value class
      as.character(lo)
    } else {
      paste0(lo, "-", hi - 1)              # closed range shown as lo - (hi-1)
    }
  })
  
  
  
  # Create categorical colour palette
  pal <- colorBin(
    palette   = colors,
    domain    = values(richness_total),
    bins      = breaks,
    na.color  = "transparent"
  )
  
  x_l <- leaflet() %>%
    addTiles() %>%
    addProviderTiles(providers$CartoDB.Positron) %>%
    addMiniMap(width = 150, height = 150) %>%
    addRasterImage(
      richness_total,
      colors   = pal,
      opacity  = 0.3,
      project  = TRUE,
      maxBytes = 20 * 1024 ^ 2
    ) %>%
    addLegend(
      pal      = pal,
      values   = values(richness_total),
      title    = "Species Richness",
      labFormat = function(type, cuts, p)
        auto_labels  # fully automatic
    )
  
  htmlwidgets::saveWidget(
    x_l,
    paste0(
      output_dir_maps_HTML,
      "/",
      "species_richness_concave_hull_future_",SSPs[[s]],".html"
    ),
    selfcontained = T
  )
  message(paste("DONE: "))
  message(paste0("saved in ",output_dir_maps_HTML))
 }
  return("DONE!")
}

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
  readxl::read_xlsx(
    "E:/CSOSA/Dropbox/VAVILOV_2.0/Results/SUMMARY_FILES/PROJECT_STATUS.xlsx",
    "Hoja1"
  )

#template
#load template for rasterize
template_path <- paste0(
  "//catalogue/MultifLandscapesA1706/1.Data/RAW/Input_data/climate_data/2_5min/present/",
  "wc2.1_2.5m_bio_1.tif")
  
#   
#   rast(
#   ))
# 
# template[which(!is.na(template[]))] <- 1

num_classes <- 8
summary_richness_function(basedir,resultdir,template,species_selected,num_classes)
