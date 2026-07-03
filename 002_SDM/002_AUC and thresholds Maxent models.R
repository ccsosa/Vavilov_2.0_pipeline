### this script can be used to put the AUC values and thresholds of the 
# calibrated models in a single table

# author: Tobias Fremout (tobias.fremout@gmail.com)

# set directory where the results of the modelling are stored # ADAPT FOLDER IF NECESSARY
resultdir <- "//catalogue/MultifLandscapesA1706/1.Data/Results/SDM results"
species_selected <- as.data.frame(
  readxl::read_xlsx(
    "//catalogue/MultifLandscapesA1706/1.Data/Results/species_lists/Species_20250304_edible_part_curated.xlsx",
    col_names = T
  )
)
un_sp <- unique(species_selected$species)
# load evaluation data 
setwd(dir = paste0(resultdir, "/Model evaluation"))
files <- list.files()
files_list <- list()
for (i in 1:length(files)) {
  files_list[[i]] <- read.csv(files[i])
}
AUC_table <- do.call(rbind, files_list)
AUC_table$species <- files
AUC_table$species2 <- sub("model_evaluation_","",AUC_table$species)
AUC_table$species2 <- sub(".csv","",AUC_table$species2)

sp_modeled <- AUC_table$species2[un_sp %in% AUC_table$species2]
others <- AUC_table$species2[!un_sp %in% AUC_table$species2]

# save
setwd(dir = resultdir)
write.csv(AUC_table, "AUC_Maxent_all.csv")

# load thresholds
setwd(dir = paste0(resultdir, "/Model thresholds"))
files <- list.files()
files_list <- list()
for (i in 1:length(files)) {
  files_list[[i]] <- read.csv(files[i])
}
thr_table <- do.call(rbind, files_list)
thr_table$species <- files

# save
setwd(dir = resultdir)
write.csv(thr_table, "thresholds_Maxent.csv")


AUC_table_valid <- AUC_table[which(AUC_table$auc.val.avg>=0.7),]
# save
setwd(dir = resultdir)
write.csv(AUC_table_valid, "AUC_Maxent_valid_all.csv")

#saving warning species
AUC_table_warning <- AUC_table[which(AUC_table$auc.val.avg>=0.6 &
                                       AUC_table$auc.val.avg<0.7),]

# save
setwd(dir = resultdir)
write.csv(AUC_table_warning, "AUC_Maxent_warning.csv")



