## master script for preprocessing

# set user of the script
#user <- "Tobias"
user <- "Chrystian"

if(user == "Chrystian") {
  scriptdir <-  "D:/PROGRAMAS/Dropbox/VAVILOV_2.0/4.Scripts/TO_DOCUMENT/4.Scripts/000_Preprocessing"
}
################################################################################
#Step 0. Provide folder to save inputs and results
basedir <- "D:/PROGRAMAS/Dropbox/VAVILOV_2.0"
setwd(dir = scriptdir)
#calling creating folders
source("000_Creating_folders.R")
create_folders_func(basedir = basedir)
# Step 0. Alternative.create folder structure only if you want to edit the source file
# file.edit("000_Creating_folders.R")
################################################################################
#get predictors and template regions rasters
#Step 1. Get bioclimatic layers and cut using soil grids
file.edit("000_get_bioclim.R")

#Step 2.. Obtaining a template target region of Latam (Latam form  as 1 as values for the template)
file.edit("000_get_target_region_2_5.R")
################################################################################
#Step 3. Get species taxonomic accepted names. This code loads a species list,
#use TNRS to get accepted names and uses
#TNRS to get accepted species name. If it is not available then it is used kewr R package. 
#Further, the World Checklist of plants to get the native distribution using the ipni id and saves 
#the file Species_list_20260212_processed.xlsx. This is used for further steps (download occurrences)
file.edit("001_get_sp_status_kew.R")
################################################################################
#Step 4. Obtain occurences from GBIF using a parallelized approach 
#run in several rounds. This is using Rgbif. Thus, the server can fail!
# input data:
#   species_set: Species list in a list
#   ex: #geographical extent of a raster
#   numCores= number of cores (numeric) #never use all cores!

file.edit("002_GBIF data import.R")
################################################################################
#Step 5. Joining GBIF occurrences in one unique file
# remove data in urban areas
# remove records with imprecise coordinates
# remove data older than 1950
# remove data with region mismatch
file.edit("003_join_gbif_data.R")
################################################################################
#Step 6.
#cleaning data
file.edit("004_GBIF data cleaning.R")
################################################################################
#Step 7.
# Apply coordinate cleaning tests to flag problematic records
# Tests include:
# - capitals: flags points near country capitals (within 5000m)
# - centroids: flags points near country centroids (within 5000m)
# - equal: flags points with equal lat/lon coordinates
# - gbif: flags points outside country borders according to GBIF
# - institutions: flags points near biodiversity institutions (within 5000m)
# - outliers: flags geographic outliers within species occurrence points
# - seas: flags points located in the sea
# - zeros: flags points with zero coordinates

#input data:
# x: occurrences from the last step
#basedir: base dir
# cleanDir: Directory to save the clean data
file.edit("005_CoordinateCleaner_approach_GBIF.R")
