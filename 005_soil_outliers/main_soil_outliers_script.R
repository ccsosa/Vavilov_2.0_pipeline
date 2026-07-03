## master script for soil outliers

# set user of the script
#user <- "Tobias"
user <- "Chrystian"

if(user == "Chrystian") {
  scriptdir <-  "D:/PROGRAMAS/Dropbox/VAVILOV_2.0/4.Scripts/TO_DOCUMENT/4.Scripts/005_soil_outliers"
}
basedir <- "D:/PROGRAMAS/Dropbox/VAVILOV_2.0"
setwd(dir = scriptdir)


################################################################################
#Step 1. Get raster layers value per occurences
file.edit("000_getting_soil_vars_parallel.R")
################################################################################
#Step 2. Find the extreme species and report in an Excel file per ecosystem
file.edit("001_filtering_outilers_sp.R") 
################################################################################
#Step 3.  Find the extreme species and report in an Excel file
file.edit("002_filtering_outilers_sp_continental_scale.R")
################################################################################
#Step 4. Create maps and tif files per ecosystem and layer (potential)
file.edit("003_maps_groups_ecosystems_potential.R")
################################################################################
#Step 5. Create maps and tif files per ecosystem and layer (realized)
file.edit("004_maps_groups_ecosystems_realized.R")