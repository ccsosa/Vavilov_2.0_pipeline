## master script for preprocessing

# set user of the script
#user <- "Tobias"
user <- "Chrystian"

if(user == "Chrystian") {
  scriptdir <-  "D:/PROGRAMAS/Dropbox/VAVILOV_2.0/4.Scripts/TO_DOCUMENT/4.Scripts/001_observed_richness_buffer/"
}
################################################################################
#Step 1 Create buffer rasters (Using the GBIF cleaned ata)
basedir <- "D:/PROGRAMAS/Dropbox/VAVILOV_2.0"
setwd(dir = scriptdir)
#INPUTS:
#clean data from GBIF
#template from bio1 to use as template
#this made buffer at 30 seconds (10 km ratio)
file.edit("001_buffer.R")
################################################################################
#step 2. Summary richness for buffer rasters (10 km ratio)
file.edit("002_summary_richness_potential.R")
################################################################################