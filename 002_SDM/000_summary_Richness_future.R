### this script can be used to put the AUC values and thresholds of the
# calibrated models in a single table

# author: Chrystian Sosa
require(ggplot2)
library(tidyterra)
library(cowplot)
library(terra)
library(viridis)
library(sf)
library(bit)

user <- "Chrystian"
if(user == "Chrystian") {
  basedir <-  "/catalogue/MultifLandscapesA1706"
}
################################################################################
#get new function

CH_get_Raster_func <- function(A,B,InDir,conv_hull_dir,species_list,current){
  
  files_ch <- lapply(A:B,function(i){
    # message(i)
    # i <- 1
    if(file.exists(paste0(InDir,"/",species_list[[i]],".tif"))& 
       file.exists(paste0(conv_hull_dir,"/","convex hull ",species_list[[i]],".tif"))
    ){
      # get files
      x <- terra::rast(paste0(conv_hull_dir,"/","convex hull ",species_list[[i]],".tif"))
      #resampling
      x <- terra::resample(x,current,"near")
      x1 <- terra::rast(paste0(InDir,"/",species_list[[i]],".tif"))
      x <- x1*x
      # x2 <- x2+x
    } else {
      message(paste0(species_list[[i]]," is not available"))
      # # x1 <- x1
      x <- NULL
    }
    return(x)
  })
  files_ch <- terra::rast(files_ch)
  # #creating summary richness file
  richness <- sum(files_ch ,na.rm = T)
  return(richness)
}
#MaxEnt results
# resultdir <- "/catalogue/MultifLandscapesA1706/1.Data/Results/CWR_results"
resultdir <- "/catalogue/MultifLandscapesA1706/1.Data/Results/SDM results"

# set target region
summary_richness_function <- function(basedir,resultdir,future_proj){
  
  message("loading shapefile and metrics")
  adm0 <- vect(paste0(basedir, "/1.Data/RAW/input_data/adm0/adm0_Latam", ".shp"))
  #load metrics
  metrics <- read.csv(paste0(resultdir,"/","AUC_Maxent_valid.csv"))
  # #gettting speies lists to get results
  species_list <- metrics$species
  species_list <- sub(pattern = "model_evaluation_",replacement = "",species_list)
  species_list <- sub(pattern = ".csv",replacement = "",species_list)
  print(length(species_list))
  species_list <- species_list[which(species_list!="Cabralea canjerana")]
  print(length(species_list))
  # 
  #adding presence absence results folder
  year <- 2050
  InDir <- paste0(resultdir,"/Distribution maps/Presence-absence")
  FutDir <- paste0(resultdir,"/Distribution maps/Future/",year,"/consensus maps")
  ################################################################################
  ####First step (summary file)
  
  message("doing summary for current")
  message("Reading all rasters from species list for future")
  
  
  species_list_filtered <- list.files(FutDir,pattern = ".tif")
  species_list_filtered <- sub(".tif","",species_list_filtered)
  species_list_filtered <- species_list_filtered[species_list_filtered %in% species_list]
  species_list_filtered <- paste0(species_list_filtered,".tif")
  
  
  
  if(!file.exists(paste0(resultdir,"/Distribution maps/Future/",year,"/Future_proj_sum_richness.tif"))){
       # if(!file.exists( paste0(resultdir,"/","species_richness_current.tif"))){
    x <- terra::rast(paste0(FutDir,"/",species_list_filtered))
    x_chunks <- bit::chunks(1,  length(species_list_filtered), by=500) #200000
    sum_i_CH <- list()
    for(i in 1:length(x_chunks)){
      #get limits
      A <- as.numeric(as.character(x_chunks[[i]])[1])
      B <- as.numeric(as.character(x_chunks[[i]])[2])
      
      sum_i_CH[[i]] <- sum(x[[A:B]],na.rm=T)
    }
    
    sum_i_CH <- terra::rast(sum_i_CH)
    x_sum <- sum(sum_i_CH,na.rm = T)
    x_sum[which(x_sum[]==0)] <- NA
    writeRaster(x_sum,paste0(resultdir,"/Distribution maps/Future/",year,"/Future_proj_sum_richness.tif"))
    future <- x_sum
  } else {
    future <- rast(paste0(resultdir,"/Distribution maps/Future/",year,"/Future_proj_sum_richness.tif"))
  }
  
    current <- terra::rast(paste0(resultdir,"/","species_richness_current.tif"))
  
################################################################################
    #PLOTTING
    
    zmin <- min(global(current, "min", na.rm=TRUE)[1,1],
                global(future,  "min", na.rm=TRUE)[1,1])
    zmax <- max(global(current, "max", na.rm=TRUE)[1,1],
                global(future,  "max", na.rm=TRUE)[1,1])
    
    adm0_sf <- sf::st_as_sf(adm0)
    
    border_theme <- theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))
    
    p1 <- ggplot() +
      geom_spatraster(data = current) +
      geom_sf(data = adm0_sf, fill = NA, color = "grey17", linewidth = 0.3) +
      scale_fill_viridis_c(name = "Richness", limits = c(zmin, zmax),
                           na.value = "transparent") +
      labs(title = "Current") +
      theme_void() +
      theme(legend.position = "right") +
      border_theme
    
    p2 <- ggplot() +
      geom_spatraster(data = future) +
      geom_sf(data = adm0_sf, fill = NA, color = "grey17", linewidth = 0.3) +
      scale_fill_viridis_c(name = "Richness", limits = c(zmin, zmax),
                           na.value = "transparent") +
      labs(title = year) +
      theme_void() +
      theme(legend.position = "right") +
      border_theme
    
    # Strip legends from main plots
    p1_noleg <- p1 + theme(legend.position = "none")
    p2_noleg <- p2 + theme(legend.position = "none")
    
    # Extract shared legend from throwaway plot
    legend_plot <- ggplot() +
      geom_spatraster(data = current) +
      scale_fill_viridis_c(name = "Richness", limits = c(zmin, zmax),
                           na.value = "transparent") +
      theme_void() +
      theme(legend.position = "right",
            legend.key.height = unit(3, "cm"))
    
    shared_legend <- cowplot::get_legend(legend_plot)
    
    # Combine
    combined <- cowplot::plot_grid(
      cowplot::plot_grid(p1_noleg, p2_noleg, nrow = 1),
      shared_legend,
      nrow = 1,
      rel_widths = c(1, 0.08)
    )
    
    ggsave(paste0("LATAM_Edible_food_", "_species_richness_total", ".png"),
           combined, width = 20, height = 10, dpi = 1000)
  ################################################################################
  ####Second step (summary file in convex hull!)
  conv_hull_dir <- paste0(resultdir,"/Convex hulls")
  
  # x <- terra::rast(paste0(conv_hull_dir,"/","convex hull ",species_list[[1]],".tif"))
  # x <- terra::resample(x,richness)
  # x1 <- terra::rast(paste0(InDir,"/",species_list[[1]],".tif"))
  # x <- x*x1
  # x2 <- x
  # # for(i in 2:length(species_list)){
  species_list_filtered2 <- sub(".tif","",species_list_filtered)
  
  if(!file.exists(paste0(resultdir,"/Distribution maps/Future/",year,"/Future_proj_CH_sum_richness.tif"))){
    
  x_chunks <- bit::chunks(1,  length(species_list_filtered), by=500) #200000
  sum_i_CH <- list()
  for(i in 1:length(x_chunks)){
    #get limits
    # i <- 1
    print(i)
    A <- as.numeric(as.character(x_chunks[[i]])[1])
    B <- as.numeric(as.character(x_chunks[[i]])[2])
    
    sum_i_CH[[i]] <- CH_get_Raster_func(A=A,
                                        B=B,
                                        InDir = FutDir,
                                        conv_hull_dir = conv_hull_dir,
                                        species_list = species_list_filtered2,
                                        current=current)
  }

  sum_i_CH <- terra::rast(sum_i_CH)
  richness_ch <- terra::sum(richness_ch,na.rm = T)
  richness_ch[which(richness_ch[]==0)] <- NA
  #saving tiff
  richness_fut_ch <- richness_ch
  terra::writeRaster(richness_ch,
                     filename = paste0(resultdir,"/","CH_species_richness_",year,".tif"),
                     overwrite=T)
  
  } else {
  richness_fut_ch <- terra::rast(paste0(resultdir,"/","CH_species_richness_",year,".tif"))  
  }
  
  return("DONE!")
}


summary_richness_function(basedir,resultdir)
