## Distribution modelling script for D4R cacao

## summary SDM methodology:
# target-group background points
# spatial block cross-validation (blockCV package)
# environmental filtering
# tuning of Maxent parameter settings (ENMeval)

## content script
# part 0: Preparations
# part 1: Preparation of presence data and target group grid
# part 2: Distribution modelling

#### PART 0: Preparations ####
options(java.parameters = "-Xmx4g")

# load packages
library(raster)
library(dismo)
library(rgeos)
library(rgdal)
library(blockCV)
library(rJava)
# library(ENMeval)
# library(ecospat)
library(data.table)
library(viridis)
library(sf)
library(terra)
library(usdm)
library(pROC)
library(parallel)
library(pbapply)
# set user of the script
user <- "Chrystian"
# user <- "Tobias"

# set base directory
if (user == "Chrystian") {
  basedir <- "//catalogue/MultifLandscapesA1706"
}
if (user == "Tobias") {
  basedir <- "C:/Users/tobia/Dropbox/VAVILOV_2.0"
}

# source help functions
# scriptdir <- paste0(basedir, "/4.Scripts/SDM")
# # source(paste0(scriptdir,"/","env filtering function without rgdal.R"))
# source(paste0(scriptdir, "/", "000_response_curve.R"))

# prepare parallelization
cl <- makeCluster(2)
stopCluster(cl)
gc()

# set output directories
outputdir <- paste0(basedir, "/1.Data/Results/SDM results")
spListDir <-
  paste0(basedir, "/1.Data/Results/buffer_richness/summary")

# create dirs structure
dirs_to <- c(
  "Convex hulls",
  "Distribution maps",
  "Maxent models",
  "Model evaluation",
  "Model thresholds",
  "Response curves",
  "Shapefiles",
  "Training point maps",
  "Tuning results",
  "Variable importance",
  "Distribution maps/Presence-absence",
  "Distribution maps/Suitability",
  "Distribution maps/Pdfs"
)
for (i in 1:length(dirs_to)) {
  x <- paste0(outputdir, "/", dirs_to[[i]])
  if (!dir.exists(x)) {
    dir.create(x)
  }
}

# fix set seed
set.seed(1000)




# load raster of target countries
setwd(dir = paste0(basedir, "/1.Data/RAW/Input_data/Raster target region/30s"))
target_region <- raster("target_region.tif")
target_region[which(!is.na(target_region[]))] <- 1


# load adm0 shapefile
setwd(dir = paste0(basedir, "/1.Data/RAW/Input_data/adm0"))
adm0 <- terra::vect(paste0("adm0_Latam", ".shp"))
adm0_wrap <- terra::wrap(adm0)

if (!file.exists(paste0(basedir, "/envstack_2.5arcmin.RDS"))) {
  # load climate data
  setwd(dir = paste0(
    basedir,
    "/1.Data/RAW/Input_data/climate_data/30s/present"
  ))
  clim <- stack(list.files(pattern = ".tif$"))
  
  # get bio1 and bio12 (for outlier detection)
  bio_1 <- clim[["wc2.1_30s_bio_1"]]
  bio_12 <- clim[["wc2.1_30s_bio_12"]]
  
  # load soil and topography data
  setwd(dir = paste0(
    basedir,
    "/1.Data/RAW/Input_data/Soil and topography data/30s"
  ))
  soil_top <- stack(list.files(pattern = ".tif$", recursive = F))
  
  # put together
  envstack <- stack(clim, soil_top)
  x_names <- names(envstack)
  crs(envstack) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  # remove some climatic variables
  i <- which(
    names(envstack) %in% c(
      "wc2.1_30s_bio_8",
      "wc2.1_30s_bio_9",
      "wc2.1_30s_bio_10",
      "wc2.1_30s_bio_11",
      "wc2.1_30s_bio_13",
      "wc2.1_30s_bio_14",
      "wc2.1_30s_bio_18",
      "wc2.1_30s_bio_19"
    )
  )
  envstack <- envstack[[-i]]
  # prepare ref raster for 2.5 arcmin (ca. 5 km) spatial filtering
  envstack_2.5arcmin <- aggregate(envstack, fact = 5)
  envstack_2.5arcmin <- terra::rast(envstack_2.5arcmin)
  # saveRDS(envstack_2.5arcmin,paste0(basedir,"/envstack_2.5arcmin2.RDS"))
  
} else {
  envstack_2.5arcmin <-
    readRDS(paste0(basedir, "/envstack_2.5arcmin.RDS"))
  # envstack_2.5arcmin <- terra::rast(envstack_2.5arcmin)
  
}

envstack_2.5arcmin_W <- terra::wrap(envstack_2.5arcmin)

# for(i in 1:nlyr(envstack_2.5arcmin)){
#   plot(envstack_2.5arcmin[[i]])
# }
# names(envstack_2.5arcmin) <- c(
#   "wc2.1_30s_bio_1",
#   "wc2.1_30s_bio_12",
#   "wc2.1_30s_bio_15",
#   "wc2.1_30s_bio_16",
#   "wc2.1_30s_bio_17",
#   "wc2.1_30s_bio_2",
#   "wc2.1_30s_bio_3",
#   "wc2.1_30s_bio_4",
#   "wc2.1_30s_bio_5",
#   "wc2.1_30s_bio_6",
#   "wc2.1_30s_bio_7",
#   "bdod_5.15cm_mean_1000",
#   "cec_5.15cm_mean_1000",
#   "cfvo_5.15cm_mean_1000",
#   "clay_5.15cm_mean_1000",
#   "nitrogen_5.15cm_mean_1000",
#   "phh2o_5.15cm_mean_1000",
#   "dtm_twi_merit.dem_m_1km_s0..0cm_2017_v1.0",
#   "sand_5.15cm_mean_1000",
#   "silt_5.15cm_mean_1000",
#   "soc_5.15cm_mean_1000"
# )
################################################################################

#### PART 1: Preparation of presence data and background points ####

# check number of records per species after 5 km spatial filtering
species_set <-
  read.csv(paste0(outputdir, "/", "species_set_SDM_all.csv"))
species_set <- species_set[, 1]


################################################################################
################################################################################
################################################################################
################################################################################
cl <- makeCluster(4)


clusterExport(
  cl,
  varlist = c(
    "envstack_2.5arcmin_W",
    "outputdir",
    "adm0_wrap",
    "species_set"
  ),
  envir = environment()
)
clusterEvalQ(cl, {
  options(java.parameters = "-Xmx4g")
  library(rJava)
  library(dismo)
  library(terra)
  library(raster)
  library(viridis)
  library(sf)
})


#6044 Asclepias subulata
#### PART 2: Suitability modelling ####
for (i in 1:length(species_set)) {
  print(i)
# results <- pblapply(species_set,function(species_name){
# species_set[-c(1:4)]
# results <- pblapply(species_set[-c(1:5900)], function(species_name) {
  # i <- 6647 problematic
  
  species_name <- species_set[i]
  # envstack_2.5arcmin_TO <- terra::unwrap(envstack_2.5arcmin_W)
  # adm0_TO <- terra::unwrap(adm0_wrap)
  message(species_name)
  if (file.exists(paste0(
    outputdir,
    "/Maxent models",
    "/",
    "maxent_model_",
    species_name,
    ".RData"
  ))) {
    if (!file.exists(paste0(
      outputdir,
      "/Distribution maps/Suitability/",
      species_name,
      ".tif"
    ))) {
      setwd(dir = paste0(outputdir, "/Maxent models"))
      load(paste0("maxent_model_", species_name, ".RData"))
      bestmod_eval <- maxent_mod
      # setwd(dir = paste0(outputdir, "/Model evaluation"))
      # bestmod_eval <- read.csv(paste0("model_evaluation_", species_name, ".csv"))
      # make predictions of the best model, only for the target region
      # pr <- predict(envstack_2.5arcmin, bestmod_eval, type = 'cloglog')
      
      env2 <-
        # envstack_2.5arcmin_TO[[intersect(names(envstack_2.5arcmin_TO),
        #                                  names(bestmod_eval@presence))]]
        envstack_2.5arcmin[[intersect(names(envstack_2.5arcmin),
                                         names(bestmod_eval@presence))]]
      
      # env2 <- stack_raster <- stack(env2)
        # pr <- dismo::predict(
        #   object=bestmod_eval,
        #   x= raster(env2))
      # pr <- terra::predict(
      #     env2,
      #   bestmod_eval,
      #   type = 'cloglog',
      #   na.rm = TRUE
      #   #,
      #   # cores = 8
      # ) #using one to avoid issues

      # pr <- predicts::predict(bestmod_eval, env2, args = "outputformat=cloglog")
      pr <- methods::selectMethod("predict", "MaxEnt_model")(
        bestmod_eval, 
        env2, 
        args = "outputformat=cloglog"
      )      
      # plot
      # plot(pr)
      setwd(dir = paste0(outputdir, "/Model thresholds"))
      thr <- read.csv(paste0(species_name, ".csv"))
      # add "10% omission" threshold
      pr_thr <- pr > thr$thr_om
      
      
      # save the maps
      # setwd(dir = paste0(outputdir, "/Distribution maps/Presence-absence"))
      # terra::writeRaster(pr_thr,
      #             species_name,
      #             format = "GTiff",
      #             overwrite = TRUE)
      terra::writeRaster(
        pr_thr,
        paste0(
          outputdir,
          "/Distribution maps/Presence-absence/",
          species_name,
          ".tif"
        ),
        overwrite = TRUE
      )
      # setwd(dir = paste0(outputdir, "/Distribution maps/Suitability"))
      # terra::writeRaster(pr, species_name, format = "GTiff", overwrite = TRUE)
      terra::writeRaster(
        pr,
        paste0(
          outputdir,
          "/Distribution maps/Suitability/",
          species_name,
          ".tif"
        ),
        overwrite = TRUE
      )
      
      presence_species <-
        terra::vect(paste0(outputdir, "/Shapefiles", "/", species_name, ".shp"))
      
      message("Saving model, response curves and variable importance...")
      
      # save the maps as pdfs
      setwd(dir = paste0(outputdir, "/Distribution maps/Pdfs"))
      pdf(paste0(species_name, ".pdf"),
          width = 10,
          height = 5)
      par(mfrow = c(1, 2))
      par(mai = c(0.5, 0.5, 0.5, 0.5))
      plot(
        pr,
        main = paste0("Suitability ", species_name),
        legend = FALSE,
        col = viridis(100)
      )
      plot(
        # st_geometry(st_as_sf(adm0_TO)),
        st_geometry(st_as_sf(adm0)),
        
        add = TRUE,
        border = "grey20",
        lwd = 1.25
      )   # need st_as_sf, st_geometry not working on a SpatVector object
      points(
        presence_species,
        pch = 19,
        cex = 0.8,
        col = "red"
      )
      
      pr_thr[pr_thr == 0] <- NA
      plot(
        pr_thr,
        main = paste0("Presence-absence ", species_name),
        legend = FALSE,
        col = "grey50"
      )
      plot(
        st_geometry(st_as_sf(adm0)),
        
        # st_geometry(st_as_sf(adm0_TO)),
        add = TRUE,
        border = "grey20",
        lwd = 1.25
      )
      # points(presence_species, pch = 19, cex = 0.3, col = "red")
      dev.off()
      par(mfrow = c(1, 1))
      
      message("FIXED!")
      write.csv("FIXED!",paste0(outputdir,"/logs/",species_name,".txt"),row.names = F)
      print(paste0(species_name,"_FIXED!"))
    } else {
      message("DONE!--- SKIPPING")
      write.csv("FIXED!",paste0(outputdir,"/logs/",species_name,".txt"),row.names = F)
      print(paste0(species_name,"_FIXED!-SKIPPED"))
    }
    
  } else {
    message(paste0(species_name, " no modelled!"))
    print(paste0(species_name,"_NO_MODEL!"))
  }
}
# , cl = cl)

stopCluster(cl)
