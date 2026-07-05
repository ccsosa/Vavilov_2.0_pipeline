## master script for summary map files

# set user of the script
#user <- "Tobias"
user <- "Chrystian"

if(user == "Chrystian") {
  scriptdir <-  "D:/PROGRAMAS/Dropbox/VAVILOV_2.0/4.Scripts/TO_DOCUMENT/4.Scripts/006_category_maps"
}
basedir <- "D:/PROGRAMAS/Dropbox/VAVILOV_2.0"
setwd(dir = scriptdir)


################################################################################
#Step 1. Getting changes for potential distribution
file.edit("001_summary_Richness_categories_potential.R")
################################################################################
#Step 2. Getting changes for Realized distribution
file.edit("002_summary_Richness_categories_concave.R") 
################################################################################
#Step 3.  Getting changes for potential distribution for future projections
file.edit("003_summary_Richness_categories_final_climate_change.R")
################################################################################
#Step 4. Getting changes for potential distribution for future projections (realized)
file.edit("004_summary_Richness_categories_final_climate_change_realized.R")
################################################################################
