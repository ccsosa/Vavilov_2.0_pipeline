## master script for mapping

# set user of the script
#user <- "Tobias"
user <- "Chrystian"

if(user == "Chrystian") {
  scriptdir <-  "D:/PROGRAMAS/Dropbox/VAVILOV_2.0/4.Scripts/TO_DOCUMENT/4.Scripts/004_realized_niche"
}
basedir <- "D:/PROGRAMAS/Dropbox/VAVILOV_2.0"
setwd(dir = scriptdir)


################################################################################
#Step 1. create TIF for potential distribution model and obtain changes per species!
file.edit("001_GETTING_MAPS_PER_SP_POTENTIAL.R")
################################################################################
#step 2.create TIF for potential distribution model and obtain changes per species!
file.edit("002_GETTING_MAPS_PER_SP_CONV_HULL.R")
################################################################################
#Step 3. create summary maps for current!
file.edit("003_MAP_RICHNESS_SUMS_CURRENT.R")
################################################################################
#step 4.create TIF for potential distribution model and obtain changes per species!
file.edit("004_MAP_RICHNESS_SUMS_FUTURE.R")
################################################################################