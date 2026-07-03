# Load necessary libraries
require(readxl)
require(data.table)
require(pbapply)
require(openxlsx)
require(TNRS)

# ==============================================================================
# 1. LOAD SPECIES LIST
# ==============================================================================

#' Path to the input species list Excel file for the MultifLandscapes A1706
#' project. Sheet \code{"species_to_kew"} contains the taxa to be resolved.
x <-
  "//catalogue/MultifLandscapesA1706/1.Data/RAW/Input_data/Species_list_202602024.xlsx"
x2 <- readxl::read_xlsx(x, sheet = "species_to_kew")

#' Dummy integer index added to preserve original row order throughout the
#' taxonomic resolution and merging steps.
x2$id_dummy <- 1:nrow(x2)

# ==============================================================================
# 2. LOAD WCVP TAXONOMIC BACKBONE (NON-DARWIN CORE FORMAT)
# ==============================================================================

#' Root directory for the World Checklist of Vascular Plants (WCVP) flat files
#' in non-Darwin Core format, downloaded from Kew / POWO.
db_dir <-
  "//catalogue/MultifLandscapesA1706/1.Data/RAW/Input_data/KEW_POWO/WCVP_no_darwin"

#' WCVP names table. Pipe-delimited, UTF-8 encoded.
#' Key columns used downstream: \code{plant_name_id}, \code{ipni_id},
#' \code{taxon_name}.
db_tax <-
  data.table::fread(paste0(db_dir, "/", "wcvp_names.csv"),
                    sep = "|",
                    encoding = "UTF-8")

# ==============================================================================
# 3. TAXONOMIC NAME RESOLUTION VIA TNRS
# ==============================================================================

#' Unique, whitespace-trimmed species names from the input list.
#' These are submitted to TNRS for resolution against the WCVP source.
species_to_check <- unique(trimws(x2$taxon))

#' Run TNRS resolution if the cached RDS result does not yet exist;
#' otherwise load the cached result to avoid redundant API calls.
#' Settings:
#' \itemize{
#'   \item \code{sources = "wcvp"}: resolve against POWO/WCVP only.
#'   \item \code{classification = "wfo"}: use World Flora Online for family.
#'   \item \code{mode = "resolve"}, \code{matches = "best"}: return the
#'     single best match per name.
#' }
if(!file.exists("//catalogue/MultifLandscapesA1706/1.Data/Process/tax_Status.RDS")){
  tax_status <- TNRS(
    species_to_check,
    sources = c("wcvp"),
    classification = "wfo",
    mode = "resolve",
    matches = "best",
    accuracy = NULL,
    skip_internet_check = FALSE
  )
} else {
  tax_status <- readRDS("//catalogue/MultifLandscapesA1706/1.Data/Process/tax_Status.RDS")
}

# ==============================================================================
# 4. EXTRACT IPNI IDENTIFIER FROM ACCEPTED NAME URL
# ==============================================================================

#' Parse the IPNI numeric ID from the full POWO accepted-name URL.
#' Example URL:
#'   https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:123456-1
#' After stripping the prefix, \code{ipni_id} contains the bare LSID suffix
#' (e.g. \code{"123456-1"}) used to join against \code{db_tax$ipni_id}.
tax_status$ipni_id <- tax_status$Accepted_name_url
tax_status$ipni_id <- trimws(
  sub(pattern = "https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:",
      replacement = "",
      x = tax_status$ipni_id)
)

# ==============================================================================
# 5. APPLY MATCH-SCORE THRESHOLD
# ==============================================================================

#' Binarise \code{Overall_score} into \code{final_status}:
#' \itemize{
#'   \item \code{1} — score >= 0.95 (accepted match).
#'   \item \code{0} — score < 0.95 or \code{NA} (rejected / no match).
#' }
#' Only records with \code{final_status == 1} are carried forward to the
#' distribution lookup step.
tax_status$final_status <- tax_status$Overall_score
tax_status$final_status[which(tax_status$final_status >= 0.95)] <- 1
tax_status$final_status[which(tax_status$final_status <  0.95)] <- 0
tax_status$final_status[which(is.na(tax_status$final_status))]  <- 0

# ==============================================================================
# 6. FALLBACK LOOKUP VIA KEWR FOR UNRESOLVED NAMES
# ==============================================================================
#' Species that failed the 0.95 threshold (\code{final_status == 0}) are
#' candidates for a manual POWO lookup via \pkg{kewr}. The block below
#' illustrates the single-species query pattern; a loop over
#' \code{tax_status_NO_RESULTS} rows would be needed for batch processing.

tax_status_NO_RESULTS <- tax_status[which(tax_status$final_status == 0), ]


################################################################################
#get results from kewr!
################################################################################
tax_no_results <- tax_status[which(tax_status$final_status==0),]


# 
# # for-loop
# pb <- txtProgressBar(min = 0,max =  nrow(results),style=3)
pb <- txtProgressBar(min = 0,max =  nrow(tax_no_results),style=3)
sp_to_test <- list()
for (i in 1:nrow(tax_no_results)) {
  setTxtProgressBar(pb, i)
  # print(i)
  
  # i <- 1
  # search for IPNI id
  s <- search_powo(query = tax_no_results$Name_submitted[i])
  # s
  # only continue if not empty
  if(s[[1]] != 0) {
    
    # get IPNI identifier
    ipni <- s$results[[1]]$fqId
    ipni <- strsplit(ipni, split = "names:")[[1]][2]
    s <- lookup_powo(taxonid = ipni, distribution = F)
    tidied <- tidy(s)
    tidied_s <- tidyr::unnest(tidied, cols=name)
    tidied_s$ID <- tax_no_results$ID[i]
    sp_to_test[[i]] <- tidied_s
  } else {
    sp_to_test[[i]] <- NULL 
  }
}


# Remove NULLs and bind rows (handles mismatched columns automatically)
sp_to_test <- bind_rows(sp_to_test[!sapply(sp_to_test, is.null)])


# Join submitted names efficiently to avoid missmatches
sp_to_test <- sp_to_test %>%
  left_join(
    tax_status %>% select(ID, submitted_name = Name_submitted),
    by = "ID"
  )


# Filter out family and genus level matches
sp_to_test_filtered <- sp_to_test %>%
  filter(!rank %in% c("FAMILY", "GENUS"))



#getting ipni ids to match distributions
sp_to_test_filtered$fqId <- sub(pattern = "urn:lsid:ipni.org:names:",replacement = "",sp_to_test_filtered$fqId)

#fields for final file
tax_status_final <- tax_status[,c("ID","Name_submitted","Overall_score","Taxonomic_status","Accepted_name",
                                  "Accepted_name_id","Accepted_name_rank","Accepted_family","ipni_id",
                                  "final_status")]

#service used
tax_status_final$taxonomic_service <- NA
tax_status_final$taxonomic_service[which(tax_status_final$final_status==1)] <- "TNRS"

#filling out final file fields
if(nrow(sp_to_test_filtered)>0){
  for(i in 1:nrow(sp_to_test_filtered)){
    # i <- 1
    tax_status_final$Overall_score[which(tax_status_final$ID==sp_to_test_filtered$ID[[i]])] <- NA
    tax_status_final$Taxonomic_status[which(tax_status_final$ID==sp_to_test_filtered$ID[[i]])] <- 
      sp_to_test_filtered$taxonomicStatus[[i]]
    tax_status_final$Accepted_name[which(tax_status_final$ID==sp_to_test_filtered$ID[[i]])] <- 
      sp_to_test_filtered$name[[i]]
    tax_status_final$Accepted_name_id[which(tax_status_final$ID==sp_to_test_filtered$ID[[i]])] <- 
      sp_to_test_filtered$fqId[[i]] 
    tax_status_final$Accepted_name_rank[which(tax_status_final$ID==sp_to_test_filtered$ID[[i]])] <- 
      sp_to_test_filtered$rank[[i]] 
    tax_status_final$Accepted_family[which(tax_status_final$ID==sp_to_test_filtered$ID[[i]])] <- 
      sp_to_test_filtered$family[[i]]   
    tax_status_final$final_status[which(tax_status_final$ID==sp_to_test_filtered$ID[[i]])] <- 1
    tax_status_final$taxonomic_service[which(tax_status_final$ID==sp_to_test_filtered$ID[[i]])] <- 
      "kewr"
    tax_status_final$ipni_id[which(tax_status_final$ID==sp_to_test_filtered$ID[[i]])] <- 
      sp_to_test_filtered$fqId[[i]]
  }
  
}

tax_status_final_to_save <- tax_status_final[which(tax_status_final$final_status==1),] 
tax_status_final_to_save$powo_id <- NA

# Fast vectorized lookup instead of loop to obtain POWO id from IPNI
ipni_to_powo <- setNames(db_tax$plant_name_id, db_tax$ipni_id)
tax_status_final_to_save$powo_id <- ipni_to_powo[tax_status_final_to_save$ipni_id]


# ==============================================================================
# 8. LOAD WCVP DISTRIBUTION AND FILTER TO LATIN AMERICA
# ==============================================================================

#' WCVP distribution table. Pipe-delimited, UTF-8 encoded.
#' Key columns: \code{plant_name_id}, \code{continent}, \code{region},
#' \code{area}, \code{introduced}, \code{extinct}, \code{location_doubtful}.
db_dist <-
  data.table::fread(paste0(db_dir, "/", "wcvp_distribution.csv"),
                    sep = "|",
                    encoding = "UTF-8")

#' Subset tax_status to records with an accepted match and a valid POWO ID.
tax_status_final <- tax_status[which(tax_status$final_status == 1), ]

#' Unique POWO plant_name_ids for accepted taxa (NA removed).
unique_taxa <- unique(tax_status_final$powo_id[!is.na(tax_status_final$powo_id)])

#' Subset WCVP names and distribution tables to the accepted taxa only,
#' reducing memory footprint before the continent/region filters.
db_tax_subset  <- db_tax[db_tax$plant_name_id %in% unique_taxa, ]
db_dist_subset <- db_dist[db_dist$plant_name_id %in% db_tax_subset$plant_name_id, ]

#' Restrict to the Americas continents before applying the finer regional filter.
db_dist_subset <- db_dist_subset[db_dist_subset$continent %in%
                                   c("NORTHERN AMERICA", "SOUTHERN AMERICA"), ]

#' Retain only Latin-American WGSRPD Level 2 regions.
db_dist_subset <- db_dist_subset[db_dist_subset$region %in% c(
  "Central America",
  "Caribbean",
  "Brazil",
  "Mexico",
  "Northern South America",
  "Southern South America",
  "Western South America"
), ]

#' Remove doubtful, introduced, and extinct occurrences to retain only
#' confirmed native records for the modelling pipeline.
db_dist_subset <- db_dist_subset[which(db_dist_subset$introduced        == 0), ]
db_dist_subset <- db_dist_subset[which(db_dist_subset$extinct           == 0), ]
db_dist_subset <- db_dist_subset[which(db_dist_subset$location_doubtful == 0), ]

#' Restructure to a tidy output data.frame and resolve taxon names by joining
#' on \code{plant_name_id}.
db_dist_subset <- data.frame(
  plant_name_id = db_dist_subset$plant_name_id,
  taxon_name    = NA,
  continent     = db_dist_subset$continent,
  region        = db_dist_subset$region,
  area          = db_dist_subset$area
)

#' Add accepted taxon names by matching \code{plant_name_id} to the WCVP
#' names subset.
db_dist_subset$taxon_name <- db_tax_subset$taxon_name[match(
  db_dist_subset$plant_name_id,
  db_tax_subset$plant_name_id
)]

# ==============================================================================
# 9. WRITE RESULTS TO A COPY OF THE INPUT WORKBOOK
# ==============================================================================

#' Output path: a copy of the original species list workbook with three new
#' sheets appended — full taxonomic resolution results, the filtered subset,
#' and the Latin-American native distribution records.
copy_file <- "//catalogue/MultifLandscapesA1706/1.Data/Process/Species_list_20260212_processed.xlsx"

# Create a copy of the original input file to preserve the raw data
file.copy(from = x, to = copy_file, overwrite = TRUE)

# Load the copy as an openxlsx workbook object
wb <- loadWorkbook(copy_file)

#' Sheet 1: full TNRS resolution results (all scores, all submitted names).
addWorksheet(wb, "taxonomic_status")
writeData(wb, "taxonomic_status", tax_status)

#' Sheet 2: filtered resolution results (overall_score >= 0.95 only).
addWorksheet(wb, "taxonomic_status_filtered")
writeData(wb, "taxonomic_status_filtered", tax_status_final)

#' Sheet 3: confirmed native Latin-American distribution records for accepted taxa.
addWorksheet(wb, "distribution_latam")
writeData(wb, "distribution_latam", db_dist_subset)

# Save the updated workbook to the copy path
saveWorkbook(wb, copy_file, overwrite = TRUE)

print(paste("Results saved to:", copy_file))