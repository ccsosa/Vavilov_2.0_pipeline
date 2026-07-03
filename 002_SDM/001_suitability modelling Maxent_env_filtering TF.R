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
library(ENMeval)
# library(ecospat)
library(data.table)
library(viridis)
library(sf)
library(terra)
library(usdm)
library(pROC)
library(parallel)

# set user of the script
user <- "Chrystian"
# user <- "Tobias"

# set base directory
if(user == "Chrystian") {
  basedir <- "//catalogue/MultifLandscapesA1706"
}
if(user == "Tobias") {
   basedir <- "C:/Users/tobia/Dropbox/VAVILOV_2.0"
}

# # source help functions
# scriptdir <- paste0(basedir, "/4.Scripts/SDM")
source(paste0(scriptdir,"/","env filtering function without rgdal.R"))
source(paste0(scriptdir,"/","000_response_curve.R"))

# prepare parallelization
#USE 2 CORES
cl <- makeCluster(2)
stopCluster(cl)
gc()

# set output directories
outputdir <- paste0(basedir, "/1.Data/Results/SDM results")
spListDir <- paste0(basedir, "/1.Data/Results/buffer_richness/summary")

# create dirs structure
dirs_to <- c(
  "Convex hulls","Distribution maps","Maxent models","Model evaluation","Model thresholds",
  "Response curves","Shapefiles","Training point maps","Tuning results","Variable importance",
  "Distribution maps/Presence-absence","Distribution maps/Suitability","Distribution maps/Pdfs")
for(i in 1:length(dirs_to)){
  x <- paste0(outputdir,"/",dirs_to[[i]])
  if(!dir.exists(x)){
    dir.create(x)
  }
}

# fix set seed
set.seed(1000)

# load climate data
setwd(dir = paste0(basedir, "/1.Data/RAW/Input_data/climate_data/30s/present"))
clim <- stack(list.files(pattern = ".tif$"))

# get bio1 and bio12 (for outlier detection)
bio_1 <- clim[["wc2.1_30s_bio_1"]]
bio_12 <- clim[["wc2.1_30s_bio_12"]]

# load soil and topography data
setwd(dir = paste0(basedir, "/1.Data/RAW/Input_data/Soil and topography data/30s"))
soil_top <- stack(list.files(pattern = ".tif$",recursive = F))

# put together
envstack <- stack(clim, soil_top)
x_names <- names(envstack)
crs(envstack) <- "+proj=longlat +datum=WGS84 +no_defs"

# remove some climatic variables
i <- which(names(envstack) %in% c("wc2.1_30s_bio_8",
                                  "wc2.1_30s_bio_9", 
                                  "wc2.1_30s_bio_10", 
                                  "wc2.1_30s_bio_11", 
                                  "wc2.1_30s_bio_13", 
                                  "wc2.1_30s_bio_14", 
                                  "wc2.1_30s_bio_18", 
                                  "wc2.1_30s_bio_19"))
envstack <- envstack[[-i]]

# load raster of target countries
setwd(dir = paste0(basedir, "/1.Data/RAW/Input_data/Raster target region/30s"))
target_region <- raster("target_region.tif")
target_region[which(!is.na(target_region[]))] <- 1

# load altitude raster and crop
setwd(dir = paste0(basedir, "/1.Data/RAW/Input_data/Other spatial data/30s"))
alt <- raster("merit_DEM_1km_modelling_extent.tif")
alt <- crop(alt, envstack)

# load adm0 shapefile
setwd(dir = paste0(basedir, "/1.Data/RAW/Input_data/adm0"))
adm0 <- vect(paste0("adm0_Latam",".shp"))
# plot(adm0)

# load species set 
setwd(dir = paste0(spListDir))
species_set <- read.csv("species_selected.csv")

# prepare ref raster for 30 arcsec spatial filtering
ref_30arcsec <- envstack[[1]]

# prepare ref raster for 2.5 arcmin (ca. 5 km) spatial filtering
ref_2.5arcmin <- aggregate(ref_30arcsec, fact = 5)

# set which one of the two to use
ref <- ref_2.5arcmin

################################################################################
#loading 5 km layers
envstack_2_5m <- c(
paste0(basedir,"/1.Data/RAW/Input_data/climate_data/2_5min/present/",c(
"wc2.1_2.5m_bio_1.tif",
"wc2.1_2.5m_bio_12.tif",
"wc2.1_2.5m_bio_15.tif",
"wc2.1_2.5m_bio_16.tif",
"wc2.1_2.5m_bio_17.tif",
"wc2.1_2.5m_bio_2.tif",
"wc2.1_2.5m_bio_3.tif",
"wc2.1_2.5m_bio_4.tif",
"wc2.1_2.5m_bio_5.tif",
"wc2.1_2.5m_bio_6.tif",
"wc2.1_2.5m_bio_7.tif"
)),

paste0(basedir,"/1.Data/RAW/Input_data/Soil and topography data/2_5min/",c(
  "bdod.tif",
  "cec.tif",
  "cfvo.tif",
  "clay.tif",
  "nitrogen.tif",
  "phh2o.tif",
  "sand.tif",
  "silt.tif",
  "soc.tif"
)),
paste0(basedir,"/1.Data/RAW/Input_data/Other spatial data/2_5min/",c(
  "merit_DEM_5km_modelling_extent.tif"
))
)
  
envstack_2_5m <- envstack_2_5m[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,21,18,19,20)]

envstack_2_5m <- lapply(1:length(envstack_2_5m),function(i){
  x <- raster(envstack_2_5m[[i]])
  return(x)
})

envstack_2_5m <- stack(envstack_2_5m)
names(envstack_2_5m) <- names(envstack)

################################################################################

#### PART 1: Preparation of presence data and background points ####

# load presence points
setwd(dir = paste0(basedir, "/1.Data/RAW/occurrences/cleaned"))
p1 <- read.csv("GBIF_joined_FINAL_3.csv")
p1 <- p1[,c("species_searched", "lon", "lat")]
colnames(p1)[1] <- "species"
presence <- p1

# background points
if(!file.exists(paste0(basedir,"/1.Data/RAW/Input_data/background","/","LATAM_BG_REF.csv"))){
  background_general <- read.csv(paste0(basedir,"/1.Data/RAW/Input_data/background","/","LATAM_BG.csv"))
  background_general <- data.frame(gridSample(background_general[,c("lon", "lat")], r = ref))
  write.csv(background_general,paste0(basedir,"/1.Data/RAW/Input_data/background","/","LATAM_BG_REF.csv")) 
} else {
  background_general <- read.csv(paste0(basedir,"/1.Data/RAW/Input_data/background","/","LATAM_BG_REF.csv"),row.names = 1)
}

# check number of records per species after 5 km spatial filtering
if(!file.exists(paste0(outputdir,"/","species_set_SDM.csv"))){
  results <- data.frame(
    species = species_set$species,
    n = NA
  )
  for (i in 1:length(species_set$species)) {
    print(i)
    j <- which(presence$species == species_set$species[i])
    presence_subset <- presence[j,]
    presence_filtered <- data.frame(gridSample(presence_subset[,c("lon", "lat")], r = ref))
    results$n[i] <- nrow(presence_filtered)
  }
  i <- which(results$n > 29)
  species_set <- results$species[i]
  species_set <- sort(species_set)
  
  # save
  setwd(dir = outputdir)
  write.csv(species_set, "species_set_SDM.csv",row.names = F) 
} else {
  species_set <- read.csv(paste0(outputdir,"/","species_set_SDM.csv"))
  species_set <- species_set[,1]
}

################################################################################
################################################################################
################################################################################
################################################################################

#### PART 2: Suitability modelling ####

# set size of blocks
range <- 200000 

# for-loop
for (i in 1:length(species_set)) {

# i <- 1
  # set species name and select only the occurrence records of this species
  species_name <- species_set[i]
  message(species_name)
  message(i)
  
  j <- which(presence$species == as.character(species_name)) 
  presence_species <- presence[j,]
  presence_species <- data.frame(lon = presence_species$lon,lat = presence_species$lat)
  
  # 5 km gridsample
  presence_species <- data.frame(gridSample(presence_species[,c("lon", "lat")], r = ref))
  names(presence_species) <- c("lon", "lat")
  
  # delete environmental outliers
  # based on bio_1
  e <- raster::extract(bio_1, cbind(presence_species$lon, presence_species$lat))
  q1 <- quantile(e, 0.25, na.rm = TRUE)
  q3 <- quantile(e, 0.75, na.rm = TRUE)
  IQR <- IQR(e, na.rm = TRUE)
  lower_lim <- q1 - 3*IQR
  higher_lim <- q3 + 3*IQR
  l <- which(e < lower_lim)
  h <- which(e > higher_lim)
  if(length(l) > 1 | length(h) > 1) {
    presence_species <- presence_species[-c(l,h),]
  }
  # based on bio_12
  e <- raster::extract(bio_12, cbind(presence_species$lon, presence_species$lat))
  q1 <- quantile(e, 0.25, na.rm = TRUE)
  q3 <- quantile(e, 0.75, na.rm = TRUE)
  IQR <- IQR(e, na.rm = TRUE)
  lower_lim <- q1 - 3*IQR
  higher_lim <- q3 + 3*IQR
  l <- which(e < lower_lim)
  h <- which(e > higher_lim)
  if(length(l) > 1 | length(h) > 1) {
    presence_species <- presence_species[-c(l,h),]
  }
  
  message("fitting pseudoabsences to convex hull...")
  # define a convex hull around the presence points with a buffer corresponding to 10% of the longest axis
  # first convert to SpatialPointsDataFrame
  occ_sp <- presence_species
  coordinates(occ_sp) <- cbind(occ_sp$lon, occ_sp$lat)
  crs.geo <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  proj4string(occ_sp) <- crs.geo
  # calculate distance of longest axis of presence points
  largedist <- max(pointDistance(occ_sp, longlat = TRUE), na.rm = T)
  # define hull around presence points
  occ_coords <- as.matrix(data.frame(occ_sp$lon, occ_sp$lat))
  ch <- chull(occ_coords)
  hull_coords <- occ_coords[c(ch, ch[1]), ]
  hull <- SpatialPolygons(list(Polygons(list(Polygon(hull_coords)), ID = 1)))
  proj4string(hull) <- crs.geo
  # 20% longest axis distance buffer # TF: changed to 20%
  hull_ext <- raster::buffer(hull, width = 0.20 * largedist)
  hullRaster <- rasterize(hull_ext, target_region, field = 1)
  
  # save convex hull as raster
  setwd(paste0(outputdir, "/convex hulls", sep = ""))
  writeRaster(hullRaster, filename = paste0("convex hull ", species_name, sep = ""),
              format = "GTiff", overwrite = TRUE)
  
  # save convex hull as shapefile
  dsn <- paste0(outputdir, "/shapefiles/", sep = "")
  shapefile(x = hull_ext, 
            file = paste0(dsn, "/convex hull 1", species_name, ".shp", sep = ""),
            overwrite = TRUE)
  
  # extract altitude values
  alt_ranges <- extract(alt,occ_sp)
  
  #altitude to use in background
  
  alt_rangs<- c(min(alt_ranges, na.rm = TRUE) - 1000,
                max(alt_ranges, na.rm = TRUE) + 1000)
  
  #getting fixed altitude hullRaster
  altHull_raster  <- crop(alt,hullRaster)
  altHull_raster  <- mask(altHull_raster,hullRaster)
  
  #getting altitude fixed by convex hull
  altHull_raster[which(altHull_raster[]<alt_rangs[1])]<- NA
  altHull_raster[which(altHull_raster[]>alt_rangs[2])] <- NA
  # plot(altHull_raster)
  # obtaining background for convex hull
  nPoints <- 20000
  bgvalues <- extract(altHull_raster, background_general)
  background <- background_general[which(bgvalues==1),]
  
  
  if(nrow(background)>nPoints){
    background <- background[c(1:nPoints),]
  } else {
    message(paste0("Obtaining new points: ",(20000-nrow(background)),
                   " Old: ",nrow(background), " points"))
    bp1 <- dismo::randomPoints(altHull_raster, n = (20000-nrow(background)), 
                               tryf = 3)
    colnames(bp1) <- c("lon","lat")
    background <- rbind(background,bp1)
    # background <- data.frame(gridSample(background[,c("lon", "lat")], r = ref))
  }
  # run PCA for environmental filtering (presence)
  coord <- presence_species
  var_pca <- as.data.frame(raster::extract(envstack, coord))
  var_pca <- cbind(coord, var_pca)
  # var_pca <- na.omit(var_pca)
  var_pca <- var_pca[complete.cases(var_pca),]
  pca <- prcomp(var_pca[,3:ncol(var_pca)], scale.=TRUE)
  scores <- as.data.frame(pca$x)
  
  # environmental filtering presence points
  presence_species <- envSample(coord = var_pca[,1:2], filters=list(scores[,1], scores[,2], scores[,3]), 
                                res=list(((max(scores[,1], na.rm = T)-min(scores[,1], na.rm = T))/30), 
                                         ((max(scores[,2], na.rm = T)-min(scores[,2], na.rm = T))/30),
                                         ((max(scores[,3], na.rm = T)-min(scores[,3], na.rm = T))/30)
                                         ),do.plot = F)
  # environmental filtering background points
  coord<- background
  var_pca_abs <- as.data.frame(raster::extract(envstack, coord))
  var_pca_abs <- cbind(coord, var_pca_abs)
  var_pca_abs <- na.omit(var_pca_abs)
  scores_abs <- predict(pca, var_pca_abs)
  
  background <- envSample(var_pca_abs[,1:2], filters=list(scores_abs[,1], scores_abs[,2], scores_abs[,3]), 
                             res=list(((max(scores_abs[,1], na.rm = T)-min(scores_abs[,1], na.rm = T))/30), 
                                      ((max(scores_abs[,2], na.rm = T)-min(scores_abs[,2], na.rm = T))/30),
                                      ((max(scores_abs[,3], na.rm = T)-min(scores_abs[,3], na.rm = T))/30)
                                      ),do.plot=F)
  
  # plot the presence and background points
  bg <- background
  p <- presence_species
  coordinates(bg) <- cbind(bg$lon, bg$lat)
  coordinates(p) <- cbind(p$lon, p$lat)
  
  # save this map
  message("saving training points as PDF map...")
  
  setwd(dir = paste0(outputdir, "/Training point maps"))
  pdf(paste0(species_set[i],".pdf"), height = 8, width = 8)
  plot(target_region)
  plot(hull,border="black",add = TRUE)
  plot(hull_ext,border="red", add = TRUE)
  points(bg[,c(1,2)], pch=20, cex=0.2, col="gray30")
  points(p[,c(1,2)], pch=20, cex=0.4, col = "red")
  dev.off()
  
  # extract environmental data at presence and background locations
  e_p <- raster::extract(envstack, cbind(presence_species$lon, presence_species$lat))
  presence_species <- cbind(presence_species, e_p)
  e_bg <- raster::extract(envstack, cbind(background$lon, background$lat))
  background <- cbind(background[,1:2], e_bg)
  
  # prepare presence and background data for the spatial block cross-validation
  p <- presence_species
  p$type = 1
  a <- background
  a$type = 0
  dat <- rbind(p,a)
  coordinates(dat)<- cbind(dat$lon, dat$lat)
  crs(dat) <- crs(ref)
  # head(dat)
  # writeOGR(dat[,c(28,1,2)], species_name, driver = "ESRI Shapefile", dsn = "D:/PUL_SDM/data", overwrite = TRUE)
  
  if(nrow(p)>29){
  message(paste0("Using ",nrow(p)," occurrences"))
  message("Spatial blocks...")
    
  # get spatially independent blocks to check how many folds can be used (aim for 8 but use less
  # if necessary: each fold has to have testing presence and absence data)
  # first prepare presence-absence data
  # set the number of folds
  k <- 8
  # set the min. number of testing points minus 1
  min_n <- 3

  # try if 8 folds works
  k <- 8
  t <- try(
    sb <- spatialBlock(speciesData = dat,
                       species = "type",
                       rasterLayer = bio_1,
                       theRange = range,
                       k = k,
                       progress = F,
                       verbose=F,
                       showBlocks=F)
  )
  testpoints <- append(sb$records$test_0, sb$records$test_1)
  if (!all(testpoints > min_n)) {
    class(t) <- "try-error"
  }

  # try with 7 folds
  if (class(t) == "try-error") {
    k <- 7
    t <- try(
      sb <- spatialBlock(speciesData = dat,
                         species = "type",
                         rasterLayer = bio_1,
                         theRange = range,
                         k = k,
                         progress = F,
                         verbose=F,
                         showBlocks=F)
    )
  }
  testpoints <- append(sb$records$test_0, sb$records$test_1)
  if (!all(testpoints > min_n)) {
    class(t) <- "try-error"
  }

  # try with 6 folds
  if (class(t) == "try-error") {
    k <- 6
    t <- try(
      sb <- spatialBlock(speciesData = dat,
                         species = "type",
                         rasterLayer = bio_1,
                         theRange = range,
                         k = k,
                         progress = F,
                         verbose=F,
                         showBlocks=F)
    )
  }
  testpoints <- append(sb$records$test_0, sb$records$test_1)
  if (!all(testpoints > min_n)) {
    class(t) <- "try-error"
  }

  # try with 5 folds
  if (class(t) == "try-error") {
    k <- 5
    t <- try(
      sb <- spatialBlock(speciesData = dat,
                         species = "type",
                         rasterLayer = bio_1,
                         theRange = range,
                         progress = F,
                         k = k,
                         verbose=F,
                         showBlocks=F)
    )
  }
  testpoints <- append(sb$records$test_0, sb$records$test_1)
  if (!all(testpoints > min_n)) {
    class(t) <- "try-error"
  }

  # try with 4 folds
  if (class(t) == "try-error") {
    k <- 4
    t <- try(
      sb <- spatialBlock(speciesData = dat,
                         species = "type",
                         rasterLayer = bio_1,
                         theRange = range,
                         k = k,
                         progress = F,
                         verbose=F,
                         showBlocks=F)
    )
  }
  testpoints <- append(sb$records$test_0, sb$records$test_1)
  if (!all(testpoints > min_n)) {
    class(t) <- "try-error"
  }

  # prepare spatial block partition groups for use in ENMevaluate
  j <- which(dat$type == 1)
  occs.grp <- sb$foldID[j]
  bg.grp <- sb$foldID[-j]
  user_grp <- list(occs.grp, bg.grp)
  names(user_grp) <- c("occs.grp", "bg.grp")
  
  message("ENMEvaluate...")
  
  # tuning arguments
  tune.args <- list(fc = c("LQ","LQH","H"), rm = c(1,3,5))
  
  # to ensure parallel works
  # library(parallel)
  rscript_args = c("-e", shQuote("getRversion()"))
  
  
  # rand <- get.randomkfold(presence_species, background, k = 5)
  
  maxent_tuning <- ENMeval::ENMevaluate(
    occs = presence_species,
    bg = background,
    algorithm = "maxent.jar",
    tune.args = tune.args,
    user.grp = user_grp,
    partitions = "user",
    quiet=TRUE,
    parallel =  T,
    numCores = 9)
  
  # save the tuning results
  setwd(dir = paste0(outputdir, "/Tuning results"))
  write.csv(maxent_tuning@results, paste0(species_name, ".csv"))
  
  message("Best model using CBI...")
  
  # select model with the highest CBI
  bestmod <- which.max(maxent_tuning@results$cbi.val.avg)
  
  # get model evaluation metrics of this model
  bestmod_eval <- maxent_tuning@results[bestmod,]
  # bestmod_eval
  
  # save model evaluation metrics
  setwd(dir = paste0(outputdir, "/Model evaluation"))
  write.csv(bestmod_eval, paste0("model_evaluation_", species_name, ".csv"))
  
  # make predictions of the best model, only for the target region
  pr <- predict(envstack_2_5m, maxent_tuning@models[[bestmod]], type = 'cloglog')
  # plot
  # plot(pr)
  
  message("Evaluating model...")
  
  # make predictions for the presence and background points
  suit_pres <- predict(maxent_tuning@models[[bestmod]], presence_species)
  suit_bg <- predict(maxent_tuning@models[[bestmod]], background)
  
  # evaluate predictive ability of model
  suit_pres <- c(suit_pres)
  suit_bg <- c(suit_bg)
  ev <- dismo::evaluate(p = suit_pres, a = suit_bg)
  
  # get possible thresholds to convert suitability to presence/absence
  thr <- threshold(ev)
  
  # add "10% omission" threshold
  thr$thr_om <- quantile(suit_pres, 0.1, na.rm = TRUE)
  #thr
  
  # apply 10% omission threshold
  pr_thr <- pr
  pr_thr[which(pr_thr[]>thr$thr_om)] <- 1
  pr_thr[which(pr_thr[]<thr$thr_om)] <- 0
  
  # extra evaluation
  if(user == "Chrystian") {
    pred_pres_bin <- as.numeric(suit_pres >= thr$thr_om)
    pred_bg_bin   <- as.numeric(suit_bg >= thr$thr_om)
    
    # observed values (1 for presence, 0 for background)
    observed <- c(rep(1, length(suit_pres)), rep(0, length(suit_bg)))
    predicted <- c(pred_pres_bin, pred_bg_bin)
    
    # confusion matrix
    conf_mat <- table(Predicted = predicted, Observed = observed)
    TP <- conf_mat["1", "1"]
    FN <- conf_mat["0", "1"]
    sensitivity <- TP / (TP + FN)
    TN <- conf_mat["0", "0"]
    FP <- conf_mat["1", "0"]
    specificity <- TN / (TN + FP)
    
    bestmod_eval$sensitivity <- sensitivity
    bestmod_eval$specificity <- specificity
    bestmod_eval$TSS <- (sensitivity + specificity)-1
    bestmod_eval$OPR <- FP/(TP+FP)
    bestmod_eval$UPR <- FN/(TP+FN)
    # bestmod_eval$F1 <- (2*TP)/(FN+(2*TP)+FP)
  }
  
  # save model evaluation metrics
  setwd(dir = paste0(outputdir, "/Model evaluation"))
  write.csv(bestmod_eval, paste0("model_evaluation_", species_name, ".csv"))
  
  # save model threshold
  setwd(dir = paste0(outputdir, "/Model thresholds"))
  write.csv(data.frame(thr), paste0(species_name, ".csv"))
  
  # apply "equal sensitivity and specificity" threshold
  pr_thr <- pr > thr$equal_sens_spec
  pr_thr <- pr > thr$spec_sens
  
  # save the maps
  setwd(dir = paste0(outputdir, "/Distribution maps/Presence-absence")) 
  writeRaster(pr_thr, species_name, format = "GTiff", overwrite = TRUE)
  setwd(dir = paste0(outputdir, "/Distribution maps/Suitability"))
  writeRaster(pr, species_name, format = "GTiff", overwrite = TRUE)
  
  # only select the points inside the target region for plotting
  e <- raster::extract(pr_thr, cbind(presence_species$lon, presence_species$lat))
  j <- which(is.na(e))
  if(length(j)>0) {
    presence_species <- presence_species[-j,]
  }
  
  message("Saving model, response curves and variable importance...")
  
  # save the maps as pdfs
  setwd(dir = paste0(outputdir, "/Distribution maps/Pdfs"))
  pdf(paste0(species_name, ".pdf"), width = 10, height = 5)
  par(mfrow = c(1,2))
  par(mai = c(0.5, 0.5, 0.5, 0.5))
  plot(pr, main = paste0("Suitability ", species_name), 
       legend = FALSE,
       col = viridis(100))
  plot(st_geometry(st_as_sf(adm0)), add = TRUE, border = "grey20", lwd = 1.25)   # need st_as_sf, st_geometry not working on a SpatVector object
  points(presence_species, pch = 19, cex = 0.3, col = "red")
  pr_thr[pr_thr==0] <- NA
  plot(pr_thr, main = paste0("Presence-absence ", species_name), 
       legend = FALSE, col = "grey50")
  plot(st_geometry(st_as_sf(adm0)), add = TRUE, border = "grey20", lwd = 1.25)
  points(presence_species, pch = 19, cex = 0.3, col = "red")
  dev.off()
  par(mfrow = c(1,1))
  
  # save the model
  setwd(dir = paste0(outputdir, "/Maxent models"))
  maxent_mod <- maxent_tuning@models[[bestmod]]
  save(maxent_mod, file = paste0("maxent_model_", species_name, ".RData"))
  
  # variable importance
  varimport <- maxent_tuning@variable.importance[[bestmod]]
  setwd(dir = paste0(outputdir, "/Variable importance"))
  write.csv(varimport, paste0(species_name, ".csv"))

  # x <- respose_curve_function(maxent_mod = maxent_mod,envstack = envstack)
  # ggsave(
  #   paste0(outputdir, "/Response curves", "/",species_name,".pdf"),
  #   x,
  #   width = 12,
  #   height = 10,
  #   units = "in",
  #   dpi = 600
  # )
  
  # # response curves
  # setwd(dir = paste0(outputdir, "/Response curves"))
  # pdf(paste0("response curves ", species_name,".pdf"),
  #     width = 12, height = 10)
  # par(mai = c(0.6,0.6,0.3,0.3))
  # # dismo::response(maxent_mod) # NO LONGER WORKING
  # x
  # # plot(pr, type = "l", las = 1)
  # dev.off()
  
  message("Saving occurences, final step...")
  # save shapefile
  dsn <- paste0(outputdir, "/Shapefiles")
  occ_Sp2 <- sf::st_as_sf(occ_sp)
  #rgdal::writeOGR(occ_sp, species_name, driver = "ESRI Shapefile", dsn = dsn, overwrite = TRUE)
  sf::write_sf(occ_Sp2,paste0(dsn,"/",species_name,".shp"),delete_dsn=TRUE)
  } else {
    message("less than 30 ocurrences, no model!")
  }
    message("DONE!")
}
