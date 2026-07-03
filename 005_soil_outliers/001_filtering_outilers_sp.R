## Filtering by quantiles

## summary SDM methodology:
# target-group background points
# spatial block cross-validation (blockCV package)
# environmental filtering
# tuning of Maxent parameter settings (ENMeval)

## content script
# part 0: Preparations
# part 1: Preparation of presence data and target group grid
# part 2: Distribution modelling
#install.packages("Rcmdr", dependencies=TRUE) #to install Rcmd
#### PART 0: Preparations ####
# options(Rcmdr.console.output = FALSE)
# library(BiodiversityR, quietly = TRUE)
# options(java.parameters = "-Xmx4g")
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
#get unique ecosystems (biomes)

message("Getting biomes")

ecosystems_unique  <- unique(na.omit(filt_sp$Biome))
ecosystems_unique <-
  ecosystems_unique[which(ecosystems_unique != "N/A")]

#removing
ecosystems_unique  <-
  ecosystems_unique[which(ecosystems_unique != "")]

################################################################################

if (!file.exists(paste0(
  basedir,
  "/1.Data/Process/soil_extreme/",
  "summary_var_per_sp.csv"
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
        mad = NA
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
      
      #creating n per ecosystem
      df[ecosystems_unique] <- NA
      
      
      #evaluating n per ecosystem
      for (k in 11:23) {
        # k <- 11
        eco <- colnames(df)[k]
        df[, eco][[1]] <-   length(na.omit(x_sp[, colnames(filt_sp)[j]][which(x_sp$Biome ==
                                                                                eco)][x_sp[, colnames(filt_sp)[j]][which(x_sp$Biome == eco)] < df$quant_values[[1]]]))
        df[, eco][[2]] <-   length(na.omit(x_sp[, colnames(filt_sp)[j]][which(x_sp$Biome ==
                                                                                eco)][x_sp[, colnames(filt_sp)[j]][which(x_sp$Biome == eco)] < df$quant_values[[2]]]))
        df[, eco][[3]] <-   length(na.omit(x_sp[, colnames(filt_sp)[j]][which(x_sp$Biome ==
                                                                                eco)][x_sp[, colnames(filt_sp)[j]][which(x_sp$Biome == eco)] > df$quant_values[[3]]]))
        df[, eco][[4]] <-   length(na.omit(x_sp[, colnames(filt_sp)[j]][which(x_sp$Biome ==
                                                                                eco)][x_sp[, colnames(filt_sp)[j]][which(x_sp$Biome == eco)] > df$quant_values[[4]]]))
      }
      sp_quantile[[j]] <- df
    }
    sp_quantile <- do.call(rbind, sp_quantile)
    summary_sp[[i]] <- sp_quantile
  }
  
  #getting final results
  summary_sp <- do.call(rbind, summary_sp)
  
  write.csv(
    summary_sp,
    paste0(
      basedir,
      "/1.Data/Process/soil_extreme/",
      "summary_var_per_sp.csv"
    ),
    na = "",
    row.names = F
  )
} else {
  summary_sp <- read.csv(paste0(
    basedir,
    "/1.Data/Process/soil_extreme/",
    "summary_var_per_sp.csv"
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

if (!file.exists(paste0(
  basedir,
  "/1.Data/Process/soil_extreme/",
  "outliers_list.csv"
))) {
  numCores <- 8
  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(
    cl,
    varlist = c("ecosystems_unique", "filt_sp", "summary_sp"),
    envir = environment()
  )
  
  
  eco_outliers <- pblapply(
    X = seq_len(length(ecosystems_unique)),
    FUN = function(i) {
      # i <- 1
      #Subsetting by biome
      eco_i <- filt_sp[which(filt_sp$Biome == ecosystems_unique[[i]]), ]
      
      #values
      sp_eco_list <- list()
      #getting values per ecosystem and variables
      for (j in 4:12) {
        # print(j)
        # j <- 4
        eco_unique_filtered_data <- eco_i
        #quantile 10-90%
        x_q_10_90 <-
          quantile(eco_unique_filtered_data[, j],
                   probs = c(0.1, 0.9),
                   na.rm = TRUE)
         #quantile 5-95%
        x_q_5_95 <-
          quantile(eco_unique_filtered_data[, j],
                   probs = c(0.05, 0.95),
                   na.rm = TRUE)
        
        #quantile 50-75%
        x_q_50_75 <-
          quantile(eco_unique_filtered_data[, j],
                   probs = c(0.5, 0.75),
                   na.rm = TRUE)
        
        x_q_25 <-
          quantile(eco_unique_filtered_data[, j],
                   probs = c(0.25),
                   na.rm = TRUE)
        
        # -	ESCENARIO 1: 
        #   o	Especies que tienen Q10 mayor al Q90 de las especies: total 228 especies para variables de suelo 
        # -	ESCENARIO 2:
        #   o	UP: MEDIANA Q50 (POR CADA ESPECIE) >Q75 (TODAS LAS ESPECIES)
        # -	ESCENARIO 3:
        #   o	DOWN: MEDIANA Q75 (POR CADA ESPECIE) > Q50 (TODAS LAS ESPECIES)
        # 
        
        #copy(just in case)
        eco_unique_filtered_data_sp <- eco_unique_filtered_data
        colnames(eco_unique_filtered_data_sp)[4:12] <- 
          c( "bdod_5.15cm_mean_1000",
             "cec_5.15cm_mean_1000",
             "cfvo_5.15cm_mean_1000",   
             "clay_5.15cm_mean_1000",
             "nitrogen_5.15cm_mean_1000",
             "phh2o_5.15cm_mean_1000",   
             "sand_5.15cm_mean_1000",
             "silt_515cm_mean_1000",
             "soc_5.15cm_mean_1000")
        
        
        #getting a dataframe with species, soil variable, median for biome, and n used
        sp_eco <- eco_unique_filtered_data_sp %>%
          group_by(species_searched) %>%
          summarise(
            biome = ecosystems_unique[[i]],
            var = colnames(eco_i)[j],
            median_sp_value_eco = median(.data[[colnames(eco_unique_filtered_data_sp)[j]]], na.rm = TRUE),
            n_obs_eco = sum(!is.na(.data[[colnames(eco_unique_filtered_data_sp)[j]]]))
            
          )
        
        #species values for the values (overall)
        sp_eco$sp_median <- NA
        sp_eco$sp_mad <- NA
        sp_eco$sp_mean <- NA
        sp_eco$sp_sd <- NA
        sp_eco$n_total <- NA
        sp_eco$sp_Q5 <- NA
        sp_eco$sp_Q10 <- NA
        sp_eco$sp_Q90 <- NA
        sp_eco$sp_Q25 <- NA
        sp_eco$sp_Q75 <- NA
        sp_eco$sp_Q95 <- NA
        #proportion of n eco /total
        sp_eco$prop_eco_total <- NA
        #quantiles per eco (OVERALL)
        sp_eco$Q10_ref_eco <- x_q_10_90[[1]]
        sp_eco$Q90_ref_eco <- x_q_10_90[[2]]
        sp_eco$Q5_ref_eco <- x_q_5_95[[1]]
        sp_eco$Q95_ref_eco <- x_q_5_95[[2]]
        sp_eco$Q50_ref_eco <- x_q_50_75[[1]]
        sp_eco$Q75_ref_eco <- x_q_50_75[[2]]
        sp_eco$Q25_ref_eco <- x_q_25[[1]]
        #Initializing values
        # sp_eco$Q10 <- F
        # sp_eco$Q90 <- F
        # sp_eco$Q5 <- F
        # sp_eco$Q95 <- F
        # sp_eco$Q_10_90_eco <- NA
        sp_eco$Q_50_75_eco <- NA
        #for each species add: median, mad, mean, n (OVERALL)
        for (k in 1:nrow(sp_eco)) {
          # print(k)

          #k<- 1
          #subset per sp and variable
          # x_k_j <-
          #   summary_sp[which(
          #     summary_sp$sp == sp_eco$species_searched[[k]] &
          #       summary_sp$var == colnames(eco_unique_filtered_data_sp)[j]
          #   ), ]
          x_k_j <-
            summary_sp[which(
              summary_sp$sp == sp_eco$species_searched[[k]] &
                summary_sp$var == colnames(eco_i)[j]
            ), ]
          
          
          sp_eco$sp_median[[k]] <- x_k_j$median[[1]]
          sp_eco$sp_mad[[k]] <- x_k_j$mad[[1]]
          sp_eco$sp_mean[[k]] <- x_k_j$mean[[1]]
          sp_eco$sp_sd[[k]] <- x_k_j$sd[[1]]
          sp_eco$n_total[[k]] <- x_k_j$n_total[[1]]
          
          sp_eco$sp_Q5[[k]] <- x_k_j$quant_values[which(x_k_j$quantile==0.05)]
          sp_eco$sp_Q10[[k]] <- x_k_j$quant_values[which(x_k_j$quantile==0.1)]
          sp_eco$sp_Q25[[k]] <- x_k_j$quant_values[which(x_k_j$quantile==0.25)]
          sp_eco$sp_Q75[[k]] <- x_k_j$quant_values[which(x_k_j$quantile==0.75)]
          sp_eco$sp_Q90[[k]] <- x_k_j$quant_values[which(x_k_j$quantile==0.9)]
          sp_eco$sp_Q95[[k]] <-x_k_j$quant_values[which(x_k_j$quantile==0.95)]
          
          #calculating proportion of n/total
          sp_eco$prop_eco_total[[k]] <- (sp_eco$n_obs_eco[[k]] / sp_eco$n_total[[k]]) * 100
        }
        rm(k)
        sp_eco_list[[j]] <- sp_eco
        
      }
      rm(j)
      #getting summary for all variables
      sp_eco_list <- do.call(rbind, sp_eco_list)
      
      
      
      
      sp_eco_list <- sp_eco_list[which(sp_eco_list$n_obs_eco > 0), ]
      
    }
  )
  
  parallel::stopCluster(cl)
  


  eco_outliers <-
    do.call(rbind, eco_outliers)
  #at least five observations per biome
  eco_outliers <- eco_outliers[which(eco_outliers$n_obs_eco>4),]
  #at least 25 records in total
  # eco_outliers <- eco_outliers[which(eco_outliers$n_total>19),]
  eco_outliers <- eco_outliers[which(eco_outliers$n_total>24),]
  
  #backup
  # eco_outliers2 <- eco_outliers
  
  # eco_outliers$Q_5_95_eco <- NA
  # eco_outliers$Q_10_90_eco <- NA
  # eco_outliers$Q_50_75_eco <- NA
  # eco_outliers$Q_25_50_eco <- NA
  eco_outliers$Q50sp_Q75eco <- NA
  eco_outliers$Q75sp_Q50eco <- NA
  eco_outliers$Q5sp_Q95eco <- NA
  eco_outliers$Q95sp_Q5eco <- NA
  eco_outliers$Q10sp_Q90eco <- NA
  eco_outliers$ Q90sp_Q10eco <- NA
  eco_outliers$Q25sp_Q50eco <- NA
  eco_outliers$Q50sp_Q25eco <- NA
  
  #scenario 0
  eco_outliers$Q5sp_Q95eco[which(eco_outliers$sp_Q5 > eco_outliers$Q95_ref_eco)
  ] <- "Q05(Sp) > Q95(All spp in ecosystem)"
  eco_outliers$Q5sp_Q95eco[which(eco_outliers$sp_Q5 <  eco_outliers$Q95_ref_eco)
  ] <- "Q05(Sp) < Q95 (All spp in ecosystem)"
  
  
  eco_outliers$Q95sp_Q5eco[which(eco_outliers$sp_Q95 > eco_outliers$Q5_ref_eco)
  ] <- "Q95(Sp) > Q5(All spp in ecosystem)"
  eco_outliers$Q95sp_Q5eco[which(eco_outliers$sp_Q95 <  eco_outliers$Q5_ref_eco)
  ] <- "Q5(Sp) < Q95 (All spp in ecosystem)"
  
  #scenario 1
  eco_outliers$Q10sp_Q90eco[which(eco_outliers$sp_Q10 > eco_outliers$Q90_ref_eco)
  ] <- "Q10(Sp) > Q90 (All spp in ecosystem)"
  eco_outliers$Q10sp_Q90eco[which(eco_outliers$sp_Q10 <  eco_outliers$Q90_ref_eco)
  ] <- "Q10(Sp) < Q90 (All spp in ecosystem)"
  
  
  eco_outliers$Q90sp_Q10eco[which(eco_outliers$sp_Q90 > eco_outliers$Q10_ref_eco)
  ] <- "Q90(Sp) > Q10 (All spp in ecosystem)"
  eco_outliers$Q90sp_Q10eco[which(eco_outliers$sp_Q90 <  eco_outliers$Q10_ref_eco)
  ] <- "Q10(Sp) < Q90 (All spp in ecosystem)"
  
  #scenario 2
  eco_outliers$Q50sp_Q75eco[which(eco_outliers$sp_median > eco_outliers$Q75_ref_eco)
  ] <- "Q50(Sp) > Q75 (All spp in ecosystem)"
  eco_outliers$Q50sp_Q75eco[which(eco_outliers$sp_median < eco_outliers$Q75_ref_eco)
  ] <- "Q50(Sp) < Q75 (All spp in ecosystem)"
  #
  eco_outliers$Q75sp_Q50eco[which(eco_outliers$sp_Q75 > eco_outliers$Q50_ref_eco)
  ] <- "Q75(Sp) > Q50 (All spp in ecosystem)"
  eco_outliers$Q75sp_Q50eco[which(eco_outliers$sp_Q75 < eco_outliers$Q50_ref_eco)
  ] <- "Q75(Sp) < Q50 (All spp in ecosystem)"
  
  #scenario 3
  eco_outliers$Q25sp_Q50eco[which(eco_outliers$sp_Q25 > eco_outliers$Q50_ref_eco)
  ] <- "Q25(Sp) > Q50 (All spp in ecosystem)"
  eco_outliers$Q25sp_Q50eco[which(eco_outliers$sp_Q25 < eco_outliers$Q50_ref_eco)
  ] <- "Q25(Sp) < Q50 (All spp in ecosystem)"
  
  
  eco_outliers$Q50sp_Q25eco[which(eco_outliers$sp_median > eco_outliers$Q25_ref_eco)
  ] <- "Q50(Sp) > Q25(All spp in ecosystem)"
  eco_outliers$Q50sp_Q25eco[which(eco_outliers$sp_median < eco_outliers$Q25_ref_eco)
  ] <- "Q50(Sp) < Q25 (All spp in ecosystem)"

  message("Saving as a CSV file")
  
  
  write.csv(
    eco_outliers,
    paste0(
      basedir,
      "/1.Data/Process/soil_extreme/",
      "outliers_list2.csv"
    ),
    na = "",
    row.names = F
  )
} else {
  eco_outliers <- read.csv(paste0(
    basedir,
    "/1.Data/Process/soil_extreme/",
    "outliers_list2.csv"
  ))
}

# ################################################################################
# Q95_FINAL <- eco_outliers[which(eco_outliers$Q95==TRUE),c("species_searched","biome","var","median_sp_value_eco","n_obs_eco","sp_median","sp_mad",
#                                                           "sp_mean","sp_sd","n_total","prop_eco_total","Q95_ref_eco","Q95")]
# 
# Q90_FINAL <- eco_outliers[which(eco_outliers$Q90==TRUE),c("species_searched","biome","var","median_sp_value_eco","n_obs_eco","sp_median","sp_mad",
#                                                           "sp_mean","sp_sd","n_total","prop_eco_total","Q90_ref_eco","Q90")]
# 
# Q5_FINAL <- eco_outliers[which(eco_outliers$Q5==TRUE),c("species_searched","biome","var","median_sp_value_eco","n_obs_eco","sp_median","sp_mad",
#                                                           "sp_mean","sp_sd","n_total","prop_eco_total","Q5_ref_eco","Q5")]
# 
# Q10_FINAL <- eco_outliers[which(eco_outliers$Q10==TRUE),c("species_searched","biome","var","median_sp_value_eco","n_obs_eco","sp_median","sp_mad",
#                                                           "sp_mean","sp_sd","n_total","prop_eco_total","Q10_ref_eco","Q10")]
# 
# 
# species_pivot_counts_90 <- Q90_FINAL %>%
#   group_by(species_searched, biome) %>%
#   summarise(n_vars = n_distinct(var),
#             vars_list = paste(unique(var), collapse = "/"), .groups = "drop") 
# 
# 
# species_pivot_counts_95 <- Q95_FINAL %>%
#   group_by(species_searched, biome) %>%
#   summarise(n_vars = n_distinct(var),
#             vars_list = paste(unique(var), collapse = "/"), .groups = "drop") 
# 
# 
# species_pivot_counts_Q10 <- Q10_FINAL %>%
#   group_by(species_searched, biome) %>%
#   summarise(n_vars = n_distinct(var),
#             vars_list = paste(unique(var), collapse = "/"), .groups = "drop") 
# 
# 
# species_pivot_counts_Q5 <- Q5_FINAL %>%
#   group_by(species_searched, biome) %>%
#   summarise(n_vars = n_distinct(var),
#             vars_list = paste(unique(var), collapse = "/"), .groups = "drop") 
# 
# 
# write.csv(
#   species_pivot_counts_90,
#   paste0(
#     basedir,
#     "/1.Data/Process/soil_extreme/",
#     "outliers_list_Q90.csv"
#   ),
#   na = "",
#   row.names = F
# )
# 
# write.csv(
#   species_pivot_counts_95,
#   paste0(
#     basedir,
#     "/1.Data/Process/soil_extreme/",
#     "outliers_list_Q95.csv"
#   ),
#   na = "",
#   row.names = F
# )
# 
# write.csv(
#   species_pivot_counts_Q10,
#   paste0(
#     basedir,
#     "/1.Data/Process/soil_extreme/",
#     "outliers_list_Q10.csv"
#   ),
#   na = "",
#   row.names = F
# )
# 
# write.csv(
#   species_pivot_counts_Q5,
#   paste0(
#     basedir,
#     "/1.Data/Process/soil_extreme/",
#     "outliers_list_Q5.csv"
#   ),
#   na = "",
#   row.names = F
# )
# #%>%
#   # pivot_wider(
  #   names_from = biome,
  #   values_from = n_vars,
  #   values_fill = 0
  # ) %>%
  # rename(species = species_searched)
