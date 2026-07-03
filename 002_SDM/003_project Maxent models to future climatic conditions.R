# load packages
library(ENMeval)
library(terra)
library(gtools)
library(rJava)
library(parallel)
library(doParallel)

# set user of the script
user <- "Chrystian"
# user <- "Tobias"

#### PART 0: Preparations ####
options(java.parameters = "-Xmx16g") 
# set base directory
if(user == "Chrystian") {
  basedir <-  "//catalogue/MultifLandscapesA1706"
}

project_function <- function(basedir,GCMs,SSPs,year,numCores,n_start,n_stop){
  message("######################")
  message("loading species to run")
  # load species set for which distributions models should be made
  setwd(dir = paste0(basedir, "/1.Data/Results/SDM results"))
  #Using just valid species
  species_set <- read.csv("AUC_Maxent_valid.csv")
  # species_set <- read.csv("AUC_Maxent_warning.csv")
  species_set <- species_set$species
  #
  species_set <-species_set[n_start:n_stop]

  species_set <- 
    sub(pattern = "model_evaluation_",replacement = "",x = species_set)
  species_set <- 
    sub(pattern = ".csv",replacement = "",x = species_set)
  
  # species_set <- species_set[1:50]
  # load raster of target region
  message("loading mask")
  
  setwd(dir = paste0(basedir, "/1.Data/RAW/Input_data/climate_data/2_5min/present"))
  t_r <- terra::rast("wc2.1_2.5m_bio_1.tif")
  t_r[which(!is.na(t_r[]))] <- 1
  # load soil and topography data
  setwd(dir = paste0(basedir, "/1.Data/RAW/Input_data/Soil and topography data/2_5min"))
  soil_top <- terra::rast(list.files())
  soil_top <- terra::crop(soil_top, t_r)
  soil_top <- terra::resample(soil_top,t_r)
  soil_top <- terra::mask(soil_top, t_r)
  nms <- names(soil_top)
  nms <- gsub("-", ".", nms)
  names(soil_top) <- nms
  rm(nms);gc()
  
  #target region wrapped to use in parallel
  t_r_wrap <- terra::wrap(t_r)
  #soil_top layers wrapping to use in parallel
  soil_top_wrap <- terra::wrap(soil_top)
  
  message("creating dirs")
  
  #checking if the folders exists and create them
  #Formula:  Future/year/GCM/Emision_scenario
  future_dir <- paste0(basedir,"/","1.Data/Results/SDM results/","Distribution maps/Future")
  if(!dir.exists(future_dir)){
    dir.create(future_dir)
  }
  #Future/year
  year_dir <- paste0(future_dir,"/",year)
  if(!dir.exists(year_dir)){
    dir.create(year_dir)
  }
  
  #creating dirs
  for(i in 1:length(GCMs)){
    GCM_dir <- paste0(year_dir,"/",GCMs[[i]])
    
    if(!dir.exists(GCM_dir)){
      dir.create(GCM_dir)
    }
    for(j in 1:length(SSPs)){
      ssp_dir <- paste0(GCM_dir,"/",SSPs[[j]])
      if(!dir.exists(ssp_dir)){
        dir.create(ssp_dir)
      } 
    }
  }
  ##############################################################################
  # #closing existing parallel connections if there are some open 
  # if(exists("cl")){
  #   try(stopCluster(cl), silent = TRUE)
  # }
  # gc()
  ##############################################################################
  # #adding this to allow parallel working for R>=3.6
  # rscript_args = c("-e", shQuote("getRversion"))
  # #preparing n cores
  # 
  # cl <- parallel::makeCluster(numCores,
  #                             type = "PSOCK",
  #                             outfile = "")
  # on.exit(stopCluster(cl))  # Asegura que el clúster se detenga, incluso si hay error
  # 
  # # 
  #Loading packages in each socket
  # clusterCall(cl, function() {
  #   library(terra, quietly = TRUE) #to raster operations
  #   library(gtools, quietly = TRUE) #to load filenames 
  #   library(dismo, quietly = TRUE) #to avoid crashes in Maxent model object
  #   return(TRUE)
  # })
  
  # #exporting files to each of the n cores used 
  # parallel::clusterExport(cl, varlist=c(
  #   #Loading environment objects like dirs, GCMS, SSP, region and years    
  #   "basedir","species_set",
  #   "GCMs","SSPs",
  #   "region","year",
  #   "results",
  #   #exporting raster wrapped objects
  #   "t_r_wrap","soil_top_wrap"
  #   
  # ),envir = environment())
  #register parallel in system!
  # doParallel::registerDoParallel(cl)
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  # for-loop
  # parallel::parLapplyLB(
  #   cl=cl,
  #    X = seq_len(length(species_set)),
  #   #X = 1:2, #testing in two species
  #   fun = function (i){
  message("######################")
  message("Starting to run projections for each of the species")
  message("######################")
  
  for(i in 1:length(species_set)){
  
  # for(i in 10){
    # i <- 1
    #Unwrappping target region
    print(paste(species_set[[i]]))
    t_r_unwrapped <- terra::unwrap(t_r_wrap)
    #Unwrapping soil top layers
    soil_top_unwrap <- terra::unwrap(soil_top_wrap)
    # gc()
    # for (i in 1:length(species_set)) {
    # i <- 1
    # load model
    setwd(dir = paste0(basedir,"/","1.Data/results/SDM results", "/Maxent models"))
    
    if(file.exists(paste0("maxent_model_", species_set[i], ".RData"))){
      load(paste0("maxent_model_", species_set[i], ".RData")) # loaded as 'maxent_mod' object
      
      # load thresholds
      setwd(dir = paste0(basedir,"/","1.Data/results/SDM results/","Model thresholds"))
      thresholds <- read.csv(paste0(species_set[i], ".csv"))
      
      # get model threshold
      thr <- thresholds$thr_om
      
      # project model to future climatic conditions for different GCMs and SSPs
      for (g in 1:length(GCMs)) {
        # g <-1
        print(GCMs[g])
        # s <- 1
        for (s in 1:length(SSPs)) {
          
          print(SSPs[s])
          
          # predictor variables
          setwd(dir = paste0(basedir, "/","1.Data/RAW/Input_data/climate_data/2_5min/future/", 
                             GCMs[g], "/", SSPs[s]))       
          env <- terra::rast(gtools::mixedsort(list.files()))
          names(env) <- c(paste0("wc2.1_30s_bio_", 1:19))
          env <- terra::crop(env, t_r_unwrapped)
          env <- terra::resample(env,t_r_unwrapped)
          env <- terra::mask(env, t_r_unwrapped)
          env <- c(env, soil_top_unwrap)
          # rm(soil_top_unwrap,t_r_unwrapped);gc()
          
          
          env2 <- env[[intersect(names(env),names(maxent_mod@presence))]]
          # make model predictions
          pr <- terra::predict(env2, maxent_mod, type = 'cloglog', na.rm=TRUE,cores=numCores) #using one to avoid issues
          # apply threshold
          pr_thr <- pr > thr
          pr_thr <- pr_thr*1
          # plot(pr_thr)
          # save map
          setwd(dir = paste0(basedir,"/","1.Data/Results/SDM results/", 
                             "/Distribution maps/Future/",year,"/", GCMs[g], "/", SSPs[s]))
          terra::writeRaster(pr_thr, paste0(species_set[i],"_",year,".tif"), overwrite = TRUE)
          rm(pr_thr,pr,env);gc()
        }
      }
    } else {
      message(paste0("No results for:", species_set[[i]]))
    }
    
    message("@@@@")
    
  }#)
}

# set target region
# region<- "SA coffee"

# results
# results <- "Maxent 4"

# GCMs and SSPs
GCMs <- c(
  "ACCESS-CM2",
  "GISS-E2-1-G",
  "INM-CM5-0",
  "MIROC6",
  "MPI-ESM1-2-HR"
)
SSPs <- c("ssp245", "ssp370")
year <- "2050"
numCores <- 5

project_function(basedir,GCMs,SSPs,year,numCores,n_start = 0,n_stop = 100)
