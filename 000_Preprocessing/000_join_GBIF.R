require(data.table)


dir <- "E:/CSOSA/Dropbox/VAVILOV_2.0/DATA/DATABASES/occurrences/original"

GBIF1 <- fread(paste0(dir,"/","GBIF_data_download_", "occurrences_list", ".csv"))
GBIF2 <- fread(paste0(dir,"/","GBIF_data_download_second_round_", "occurrences_list", ".csv"))

GBIF_TOTAL <- rbind(GBIF1,GBIF2)
write.csv(GBIF_TOTAL, paste0(dir,"/","GBIF_data_original", ".csv"))
