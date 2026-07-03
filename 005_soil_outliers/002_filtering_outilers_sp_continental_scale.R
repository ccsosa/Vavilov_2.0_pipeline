## Filtering by quantiles

# capabilities("tcltk")
# load packages# library(ecospat)
library(data.table)
# library(viridis)
library(sf)
library(terra)
library(tidyr)
library(dplyr)
library(dismo)
library(parallel)
library(pbapply)

################################################################################
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
################################################################################
#Defining outcome dir for processing:
temp_dir_eco <- paste0(basedir,
                       "/1.Data/Process/soil_extreme/eco_temp")

if (!dir.exists(temp_dir_eco)) {
  dir.create(temp_dir_eco)
}


################################################################################
#loading template
### soil

soil_vars <-
  paste0(
    basedir,
    "/1.Data/RAW/Input_data/Soil and topography data/30s/",
    c(
      "bdod.tif",
      "cec.tif",
      "cfvo.tif",
      "clay.tif",
      "nitrogen.tif",
      "phh2o.tif",
      "sand.tif",
      "silt.tif",
      "soc.tif"
    )
  )
message("Getting a template for gridSample")

#loading all soil variables
template <- rast(soil_vars[[1]])
template[which(!is.na(template[]))] <- 1

################################################################################
#reading
soil_data <- readRDS(paste0(
  basedir,
  "/1.Data/Process/soil_extreme/",
  "input_data_soil_20260324.RDS"
))

# soil_data_filt <- data.frame(dismo::gridSample(soil_data[, c("lon", "lat")], r = raster::stack(template)))


################################################################################
#get species to obtain ecorregions and soil vars

message("Loading soil and biome variables")


specie_unique <-
  
  # read.csv(paste0(
  readxl::read_xlsx(paste0(
    basedir,
    "/1.Data/Results/species_lists/Species_20250304_edible_part_curated.xlsx"
    # "/1.Data/Results/SDM results/",
    # "species_set_SDM.csv"
  )
  # header = T
  )

specie_unique <- data.frame(x=unique(specie_unique$species))

################################################################################

# 
# if (!file.exists(
#   paste0(
#     basedir,
#     "/1.Data/Process/soil_extreme/",
#     "soil_data_gridSample_input.csv"
#   )
# )) {
#   filt_sp <- list()
#   for (i in 1:nrow(specie_unique)) {
#     print(i)
#     x_i <- subset(soil_data, species_searched == specie_unique$x[[i]])
#     if(nrow(x_i)>0){
#     x_i$coords <- paste0(x_i$lon, "-", x_i$lat)
#     x_i_sp <-
#       data.frame(dismo::gridSample(x_i[, c("lon", "lat")], r = raster::stack(template)))
#     x_i_sp$coords <- paste0(x_i_sp$lon, "-", x_i_sp$lat)
#     x_i <- x_i[x_i$coords %in% x_i_sp$coords,]
#     x_i$coords <- NULL
#     filt_sp[[i]] <- x_i
#     } else {
#       filt_sp[[i]] <- NULL
#     }
#   }
#   
#   filt_sp <- do.call(rbind, filt_sp)
#   write.csv(
#     filt_sp,
#     paste0(
#       basedir,
#       "/1.Data/Process/soil_extreme/",
#       "soil_data_gridSample_input.csv"
#     ),
#     na = "",
#     row.names = F
#   )
# } else {
#   filt_sp <- read.csv(
#     paste0(
#       basedir,
#       "/1.Data/Process/soil_extreme/",
#       "soil_data_gridSample_input.csv"
#     )
#   )
# }
# 
# filt_sp$coords <- paste0(filt_sp$lon, "-", filt_sp$lat)

if (!file.exists(
  paste0(basedir, "/1.Data/Process/soil_extreme/", "soil_data_gridSample_input.csv")
)) {
  library(parallel)
  library(dismo)
  library(raster)
  
  template_stack <- raster::stack(template)
  
  cl <- makeCluster(8)
  clusterExport(cl, varlist = c("soil_data", "specie_unique", "template_stack"),
                envir = environment())
  clusterEvalQ(cl, { library(dismo); library(raster) })
  
  filt_sp <- parLapply(cl, 1:nrow(specie_unique), function(i) {
    x_i <- subset(soil_data, species_searched == specie_unique$x[[i]])
    if (nrow(x_i) > 0) {
      x_i$coords <- paste0(x_i$lon, "-", x_i$lat)
      x_i_sp <- data.frame(dismo::gridSample(x_i[, c("lon", "lat")], r = template_stack))
      x_i_sp$coords <- paste0(x_i_sp$lon, "-", x_i_sp$lat)
      x_i <- x_i[x_i$coords %in% x_i_sp$coords, ]
      x_i$coords <- NULL
      return(x_i)
    } else {
      return(NULL)
    }
  })
  
  stopCluster(cl)
  
  filt_sp <- do.call(rbind, filt_sp)
  write.csv(
    filt_sp,
    paste0(basedir, "/1.Data/Process/soil_extreme/", "soil_data_gridSample_input.csv"),
    na = "", row.names = FALSE
  )
} else {
  filt_sp <- read.csv(
    paste0(basedir, "/1.Data/Process/soil_extreme/", "soil_data_gridSample_input.csv")
  )
}

filt_sp$coords <- paste0(filt_sp$lon, "-", filt_sp$lat)

################################################################################

if (!file.exists(paste0(
  basedir,
  "/1.Data/Process/soil_extreme/",
  "summary_var_per_sp_continental.csv"
))) {
  summary_sp <- list()
  
  for (i in 1:length(specie_unique$x)) {
    # i <- 1
    #subset per species
    x_sp <-
      filt_sp[which(filt_sp$species_searched == specie_unique$x[[i]]), ]
    
    sp_quantile <- list()
    #evaluating per variable
    for (j in 4:12) {
      # j <- 4
      #subset by species and variable
      x_sub <- x_sp[, j]
      #summary file
      df <- data.frame(
        sp = specie_unique$x[[i]],
        var = colnames(filt_sp)[j],
        quantile = c(0.05, 0.1,0.25,0.5,0.75, 0.9, 0.95),
        quant_values = quantile(filt_sp[, j][which(filt_sp$species_searched ==
                                                     specie_unique$x[[i]])], na.rm = T,
                                c(0.05, 0.1,0.25,0.5,0.75, 0.9, 0.95)),
        n_quant = NA,
        n_total = NA,
        mean = NA,
        sd = NA,
        median = NA,
        mad = NA,
        prop = NA
      )
      
      
      #quantile values and getting n
      df$n_quant[[1]] <-
        length(x_sub[na.omit(x_sub) < df$quant_values[[1]]])
      df$n_quant[[2]] <-
        length(x_sub[na.omit(x_sub) < df$quant_values[[2]]])
      df$n_quant[[3]] <-
        length(x_sub[na.omit(x_sub) > df$quant_values[[3]]])
      df$n_quant[[4]] <-
        length(x_sub[na.omit(x_sub) > df$quant_values[[4]]])
      df$n_quant[[5]] <-
        length(x_sub[na.omit(x_sub) > df$quant_values[[5]]])
      df$n_quant[[6]] <-
        length(x_sub[na.omit(x_sub) > df$quant_values[[6]]])
      df$n_quant[[7]] <-
        length(x_sub[na.omit(x_sub) > df$quant_values[[7]]])
      #total points
      df$n_total <-   length(x_sub[na.omit(x_sub)])
      #descriptive statistics
      df$mean  <- mean(x_sub, na.rm = T)
      df$sd  <- sd(x_sub, na.rm = T)
      df$median  <- median(x_sub, na.rm = T)
      df$mad  <- mad(x_sub, na.rm = T)
      sp_quantile[[j]] <- df
      
      
      prop_eco_total
  }
    sp_quantile <- do.call(rbind,sp_quantile)
  #getting final results
  summary_sp[[i]] <- sp_quantile
  }
  summary_sp <- do.call(rbind, summary_sp)
  
  write.csv(
    summary_sp,
    paste0(
      basedir,
      "/1.Data/Process/soil_extreme/",
      "summary_var_per_sp_continental.csv"
    ),
    na = "",
    row.names = F
  )
} else {
  summary_sp <- read.csv(paste0(
    basedir,
    "/1.Data/Process/soil_extreme/",
    "summary_var_per_sp_continental.csv"
  ))
}
################################################################################
# for(i in 1:length(ecosystems_unique)){
# colnames(eco_unique_filtered_data)[4:12] <- 
#   c( "bdod_5-15cm_mean_1000",
#      "cec_5-15cm_mean_1000",
#      "cfvo_5-15cm_mean_1000",   
#      "clay_5-15cm_mean_1000",
#      "nitrogen_5-15cm_mean_1000",
#      "phh2o_5-15cm_mean_1000",   
#      "sand_5-15cm_mean_1000",
#      "silt_5-15cm_mean_1000",
#      "soc_5-15cm_mean_1000")
if (!file.exists(paste0(
  basedir,
  "/1.Data/Process/soil_extreme/",
  "summary_per_var_continental.csv"))) {

colnames(filt_sp)[4:12] <- 
  c( "bdod_5.15cm_mean_1000",
     "cec_5.15cm_mean_1000",
     "cfvo_5.15cm_mean_1000",   
     "clay_5.15cm_mean_1000",
     "nitrogen_5.15cm_mean_1000",
     "phh2o_5.15cm_mean_1000",   
     "sand_5.15cm_mean_1000",
     "silt_5.15cm_mean_1000",
     "soc_5.15cm_mean_1000")


summary_vars <- list()
for(j in 4:12){
  df <- data.frame(
    var = colnames(filt_sp)[j],
    quantile = c(0.05, 0.1,0.25,0.5,0.75, 0.9, 0.95),
    quant_values = quantile(filt_sp[,j], na.rm = T,
                            c(0.05, 0.1,0.25,0.5,0.75, 0.9, 0.95)),
    n_quant = NA,
    n_total = NA,
    mean = NA,
    sd = NA,
    median = NA,
    mad = NA,
    prop = NA
  )
  
  
  #quantile values and getting n
  df$n_quant[[1]] <-
    length(filt_sp[,j][na.omit(filt_sp[,j]) < df$quant_values[[1]]])
  df$n_quant[[2]] <-
    length(filt_sp[,j][na.omit(filt_sp[,j]) < df$quant_values[[2]]])
  df$n_quant[[3]] <-
    length(filt_sp[,j][na.omit(filt_sp[,j]) > df$quant_values[[3]]])
  df$n_quant[[4]] <-
    length(filt_sp[,j][na.omit(filt_sp[,j]) > df$quant_values[[4]]])
  df$n_quant[[5]] <-
    length(filt_sp[,j][na.omit(filt_sp[,j]) > df$quant_values[[5]]])
  df$n_quant[[6]] <-
    length(filt_sp[,j][na.omit(filt_sp[,j]) > df$quant_values[[6]]])
  df$n_quant[[7]] <-
    length(filt_sp[,j][na.omit(filt_sp[,j]) > df$quant_values[[7]]])
  #total points
  df$n_total <-   length(filt_sp[,j][na.omit(filt_sp[,j])])
  #descriptive statistics
  df$mean  <- mean(filt_sp[,j], na.rm = T)
  df$sd  <- sd(filt_sp[,j], na.rm = T)
  df$median  <- median(filt_sp[,j], na.rm = T)
  df$mad  <- mad(filt_sp[,j], na.rm = T)
  df$prop <- df$n_quant/df$n_total *100
  summary_vars[[j]] <- df
}

summary_vars <- do.call(rbind,summary_vars)

write.csv(
  summary_vars,
  paste0(
    basedir,
    "/1.Data/Process/soil_extreme/",
    "summary_per_var_continental.csv"
  ),
  na = "",
  row.names = F
)
} else{
  summary_vars <- read.csv(paste0(
    basedir,
    "/1.Data/Process/soil_extreme/",
    "summary_per_var_continental.csv"
  ))
}






# quantiles <- c(0.05,0.1,0.25,0.5,0.75,0.90,0.95)

if (!file.exists(paste0(
  basedir,
  "/1.Data/Process/soil_extreme/",
  "outliers_list_continental.csv"
))) {


  numCores <- 8
  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(
    cl,
    varlist = c("summary_vars", "summary_sp","filt_sp","specie_unique"),
    envir = environment()
  )


outliers_continental <- pblapply(
      X = seq_len(length(specie_unique$x)),
      FUN = function(i) {

# for(i in seq_along(specie_unique$x)){
  summary_sp_i <- summary_sp[which(summary_sp$sp==specie_unique$x[[i]]),]
  
  result_i <- list()
  
  for(j in 4:12){
    summary_sp_j <- summary_sp_i[which(summary_sp_i$var==colnames(filt_sp)[j]),]
    summary_var_j <- summary_vars[which(summary_vars$var==colnames(filt_sp)[j]),]
    
    result_j <- data.frame(
      species_searched = unique(summary_sp_j$sp),
      var = colnames(filt_sp)[j],
      n_obs_sp = unique(summary_sp_j$n_total),
      sp_median = unique(summary_sp_j$median),
      sp_mad = unique(summary_sp_j$mad),
      sp_mean = unique(summary_sp_j$mean),
      sp_sd = unique(summary_sp_j$sd),
      # n_total = unique(summary_sp_j$n_total),
      sp_Q5 = summary_sp_j$quant_values[which(summary_sp_j$quantile==0.05)],
      sp_Q10 = summary_sp_j$quant_values[which(summary_sp_j$quantile==0.1)],
      sp_Q90 = summary_sp_j$quant_values[which(summary_sp_j$quantile==0.9)],
      sp_Q25 = summary_sp_j$quant_values[which(summary_sp_j$quantile==0.25)],
      sp_Q75 = summary_sp_j$quant_values[which(summary_sp_j$quantile==0.75)],
      sp_Q95 = summary_sp_j$quant_values[which(summary_sp_j$quantile==0.95)],
      # prop = NA,
      Q10_ref = summary_var_j$quant_values[summary_var_j$quantile==0.1],
      Q90_ref = summary_var_j$quant_values[summary_var_j$quantile==0.9],
      Q5_ref = summary_var_j$quant_values[summary_var_j$quantile==0.05],
      Q95_ref = summary_var_j$quant_values[summary_var_j$quantile==0.95],
      Q50_ref = summary_var_j$quant_values[summary_var_j$quantile==0.50],
      Q75_ref = summary_var_j$quant_values[summary_var_j$quantile==0.75],
      Q25_ref = summary_var_j$quant_values[summary_var_j$quantile==0.25],
      Q50sp_Q75eco = NA,
      Q75sp_Q50eco = NA,
      Q5sp_Q95eco = NA,
      Q95sp_Q5eco = NA,
      Q10sp_Q90eco = NA,
      Q90sp_Q10eco = NA,
      Q25sp_Q50eco = NA,
      Q50sp_Q25eco = NA
    )
    # ############################################################################
  
    
    #scenario 0
    result_j$Q5sp_Q95eco[which(result_j$sp_Q5 > result_j$Q95_ref)
    ] <- "Q05(Sp) > Q95(All spp)"
    result_j$Q5sp_Q95eco[which(result_j$sp_Q5 <  result_j$Q5_ref)
    ] <- "Q05(Sp) < Q95 (All spp)"

    result_j$Q95sp_Q5eco[which(result_j$sp_Q95 > result_j$Q5_ref)
    ] <- "Q95(Sp) > Q5 (All spp)"
    result_j$Q95sp_Q5eco[which(result_j$sp_Q95 <  result_j$Q5_ref)
    ] <- "Q95(Sp) < Q5 (All spp)"

    #scenario 1
    result_j$Q10sp_Q90eco[which(result_j$sp_Q10 > result_j$Q90_ref)
    ] <- "Q10(Sp) > Q90 (All spp)"
    result_j$Q10sp_Q90eco[which(result_j$sp_Q10 <  result_j$Q90_ref)
    ] <- "Q10(Sp) < Q90 (All spp)"

    result_j$Q90sp_Q10eco[which(result_j$sp_Q90 > result_j$Q10_ref)
    ] <- "Q90(Sp) > Q10 (All spp)"
    result_j$Q90sp_Q10eco[which(result_j$sp_Q90 <  result_j$Q10_ref)
    ] <- "Q90(Sp) < Q10 (All spp)"

    #scenario 2
    result_j$Q50sp_Q75eco[which(result_j$sp_median > result_j$Q75_ref)
    ] <- "Q50(Sp) > Q75 (All spp)"
    result_j$Q50sp_Q75eco[which(result_j$sp_median < result_j$Q75_ref)
    ] <- "Q50(Sp) < Q75 (All spp)"
    #
    result_j$Q75sp_Q50eco[which(result_j$sp_Q75 > result_j$Q50_ref)
    ] <- "Q75(Sp) > Q50 (All spp)"
    result_j$Q75sp_Q50eco[which(result_j$sp_Q75 < result_j$Q50_ref)
    ] <- "Q75(Sp) < Q50 (All spp)"

    #scenario 3
    result_j$Q25sp_Q50eco[which(result_j$sp_Q25 > result_j$Q50_ref)
    ] <- "Q25(Sp) > Q50 (All spp)"
    result_j$Q25sp_Q50eco[which(result_j$sp_Q25 < result_j$Q50_ref)
    ] <- "Q25(Sp) < Q50 (All spp)"


    result_j$Q50sp_Q25eco[which(result_j$sp_median > result_j$Q25_ref)
    ] <- "Q50(Sp) > Q25(All spp)"
    result_j$Q50sp_Q25eco[which(result_j$sp_median < result_j$Q25_ref)
    ] <- "Q50(Sp) < Q25(All spp)"
    
    result_i[[j]] <- result_j

  }
  result_i <- do.call(rbind,result_i)

})

stopCluster(cl)


outliers_continental <- do.call(rbind,outliers_continental)

  write.csv(
    outliers_continental,
    paste0(
      basedir,
      "/1.Data/Process/soil_extreme/",
      "outliers_list_continental.csv"
    ),
    na = "",
    row.names = F
  )

} else {
  outliers_continental <- read.csv(paste0(
    basedir,
    "/1.Data/Process/soil_extreme/",
    "outliers_list_continental.csv"
  ))
}
