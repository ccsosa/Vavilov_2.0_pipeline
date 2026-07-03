# load packages
library(terra)

# set user of the script
user <- 'Chrystian'


if (user == "Chrystian") {
  basedir <- "//catalogue/MultifLandscapesA1706/1.Data/Results"
}


# load species set for which distributions models should be made
setwd(dir = paste0(basedir, "/SDM results"))
species_set <- read.csv("species_set_SDM_all.csv")
species_set <- species_set[, 1]
# j <- which(species_set$incl == "yes")
# species_set <- species_set$species[j]

# GCMs and SSPs
GCMs <- c("ACCESS-CM2",
          "GISS-E2-1-G",
          "INM-CM5-0",
          "MIROC6",
          "MPI-ESM1-2-HR")
SSPs <- c("ssp245", "ssp370")
year <- "2050"

i_list <- list()
for(i in 1:length(GCMs)){
  x_list <- list()
  for(j in 1:2){
    
    x <- list.files(paste0(basedir,"/SDM results","/","Distribution maps/Future/",
                           year,"/",GCMs[i],"/",SSPs[j]),pattern = ".tif")
    x <- sub(pattern = "_2050.tif","",x)
    x_list[[j]] <- data.frame(
               GCM= GCMs[[i]],
               SSP= SSPs[[j]],
               sp=x)
  }
  x_list <- do.call(rbind,x_list)
  i_list[[i]] <- x_list
}

i_list <- do.call(rbind,i_list)


library(dplyr)

# summary_table <- i_list %>%
#   group_by(GCM, SSP) %>%
#   summarise(n_species = n_distinct(sp), .groups = "drop")
# 
# print(summary_table)

library(dplyr)

missing_sp <- expand.grid(GCM = GCMs, SSP = SSPs, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  reframe(
    GCM = GCM,
    SSP = SSP,
    sp_missing = setdiff(species_set, i_list$sp[i_list$GCM == GCM & i_list$SSP == SSP])
  )

# Resumen: cuántas faltan por combinación
missing_sp %>% count(GCM, SSP, name = "n_missing")

# Ver las especies faltantes
print(missing_sp)
sp_miss <- unique(missing_sp$sp_missing)
write.csv(sp_miss,paste0(basedir,"/","missing_sp.csv"),row.names = F)

