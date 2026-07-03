list1 <- readxl::read_xlsx("E:/CSOSA/Dropbox/VAVILOV_2.0/DATA/DATABASES/Species_list_LATAM.xlsx")
list1$source <- "First round"
list2 <- readxl::read_xlsx("E:/CSOSA/Dropbox/VAVILOV_2.0/DATA/DATABASES/Vavilov_second_list_LATAM.xlsx")
list2$source <- "FPI_EVERT"
final_list <- rbind(list1,list2)


writexl::write_xlsx(x = final_list,path =
                      "E:/CSOSA/Dropbox/VAVILOV_2.0/DATA/DATABASES/Vavilov_native_species_to_curate.xlsx"
                    # "D:/PROGRAMAS/Dropbox/VAVILOV_2.0/DATA/DATABASES/Vavilov_second_list_LATAM.xlsx"
)


length(unique(final_list$species
)       
)
View(final_list$species[
duplicated(final_list$species
)]
)
