#' Import and Consolidate GBIF Occurrence Records from Second-Round Downloads
#'
#' @description
#' Scans a predefined output directory for CSV files matching the pattern
#' `_list_second_round`, reads and row-binds them into a single data frame,
#' writes the consolidated result to a new CSV file
#' (`GBIF_data_download_second_round_occurrences_list.csv`), and removes the
#' individual input files afterward.
#'
#' @details
#' The script is part of the VAVILOV 2.0 pipeline. It assumes that a prior
#' step has produced per-species (or per-batch) CSV files whose names contain
#' `_list_second_round` inside the directory
#' `<basedir>/DATA/DATABASES/occurrences/Original`.
#'
#' Workflow:
#' \enumerate{
#'   \item Discover all `_list_second_round` CSV files in \code{out_dir}.
#'   \item Read each file with \code{read.csv()} and store in a list.
#'   \item Row-bind all list elements into one data frame with \code{do.call(rbind, ...)}.
#'   \item Write the consolidated data frame to
#'         \code{GBIF_data_download_second_round_occurrences_list.csv}.
#'   \item Delete the individual input files (called twice for safety).
#' }
#'
#' @section Dependencies:
#' \itemize{
#'   \item \pkg{dismo}   — GBIF querying utilities (used in upstream steps).
#'   \item \pkg{terra}   — Spatial raster/vector operations (used in upstream steps).
#'   \item \pkg{parallel} — Parallel backend support (used in upstream steps).
#'   \item \pkg{pbapply} — Progress-bar apply functions (used in upstream steps).
#'   \item \pkg{readxl}  — Excel file reading (used in upstream steps).
#' }
#'
#' @section Input:
#' CSV files matching \code{*_list_second_round*.csv} located in:
#' \preformatted{<basedir>/DATA/DATABASES/occurrences/Original}
#'
#' @section Output:
#' A single consolidated CSV file written to:
#' \preformatted{<basedir>/DATA/DATABASES/occurrences/Original/
#'   GBIF_data_download_second_round_occurrences_list.csv}
#'
#' @section Side Effects:
#' All matched input CSV files are permanently deleted from \code{out_dir}
#' after the consolidated file is written.
#'
#' @note
#' \code{basedir} must be set to the correct local path before running.
#' The active path is \code{"E:/CSOSA/Dropbox/VAVILOV_2.0"}; an alternative
#' path (\code{"D:/PROGRAMAS/Dropbox/VAVILOV_2.0"}) is kept as a comment.
#'
#' @author Tobias Fremout \email{tobias.fremout@@gmail.com}
# load packages
library(dismo)
library(terra)
library(parallel)
library(pbapply)
library(readxl)
# set base directory
# basedir <- "D:/PROGRAMAS/Dropbox/VAVILOV_2.0"
#getting basdir

basedir <- "E:/CSOSA/Dropbox/VAVILOV_2.0"
#providing outcome folder
out_dir = paste0(basedir, "/1.Processing/", "occurrences", "/original")

#detecting files with _list in the name
# x_files <- list.files(pattern = ".csv",path = out_dir,full.names = T)
x_files <- list.files(pattern = "_list",path = out_dir,full.names = T)


#reading csv files
excel_files <- lapply(1:length(x_files),function(i){
  x <- read.csv(x_files[[i]])
  return(x)
})

#putting in one file
gbif_data <- do.call(rbind,excel_files)


#removing files
for(i in 1:length(x_files)){
file.remove(x_files[[i]])
}


#writing up the csv with GBIF data
# write.csv(gbif_data, paste0(out_dir,"/","GBIF_data_download_", "occurrences_list", ".csv"))
write.csv(gbif_data, paste0(out_dir,"/","GBIF_data_download_second_round_", "occurrences_list", ".csv"))

#removing files (just in case)
for(i in 1:length(x_files)){
  file.remove(x_files[[i]])
}