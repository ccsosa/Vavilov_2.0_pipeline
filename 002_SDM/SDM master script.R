## master script for species distribution modelling silvopastoral systems

# set user of the script
#user <- "Tobias"
user <- "Chrystian"

if(user == "Chrystian") {
  scriptdir <-  "D:/PROGRAMAS/Dropbox/VAVILOV_2.0/4.Scripts/TO_DOCUMENT/4.Scripts/002_SDM"
}
#
#Step 1. Maxent suitability modelling
setwd(dir = scriptdir)
# file.edit("suitability modelling Maxent.R")
file.edit("001_suitability modelling Maxent_env_filtering TF.R")

#Step 2. Thresholds
setwd(dir = scriptdir)
file.edit("002_AUC and thresholds Maxent models.R")

#Step3.  Project to climate change SSP 245 and SSp375
setwd(dir = scriptdir)
file.edit("003_project Maxent models to future climatic conditions.R")

#check projection status
# setwd(dir = scriptdir)
# file.edit("project Maxent models to future climatic 000_status_project.R")

#Step4. Summarize
setwd(dir = scriptdir)
file.edit("004_consensus maps future climate_C_SSP.R")

# # current richness maps (really slow)
# setwd(dir = paste0(scriptdir,"/summary_maps"))
# file.edit("MAP_RICHNESS_SUMS_CURRENT.R")
# # current richness maps (really slow)
# setwd(dir = paste0(scriptdir,"/summary_maps"))
# file.edit("MAP_RICHNESS_SUMS_FUTURE.R")
