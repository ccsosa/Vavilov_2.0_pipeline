#' Create Vavilov 2.0 Project Folder Structure
#'
#' Sets up the standard directory hierarchy for the Vavilov 2.0 project
#' under a user-defined base directory. Each subdirectory is created only
#' if it does not already exist.
#'
#' @param basedir Character. Absolute path to the root project directory.
#'   Default: \code{"D:/PROGRAMAS/Dropbox/VAVILOV_2.0"}.
#'
#' @return Invisible \code{NULL}. Directories are created as a side effect.
#'
#' @details
#' The following subdirectories are created under \code{basedir}:
#' \itemize{
#'   \item \code{0.ProjectDocuments} — Project planning and administrative documents.
#'   \item \code{1.Data}             — Raw and processed input data.
#'   \item \code{2.Publications}     — Manuscript drafts and related outputs.
#'   \item \code{3.Bibliography}     — Reference management files.
#'   \item \code{4.Scripts}          — Analysis and processing scripts.
#'   \item \code{5.Others}           — Miscellaneous files.
#' }
#'
#' @note
#' Directory creation is conditional: a folder is created only if it does
#' not already exist (\code{!dir.exists()} returns \code{TRUE}).
#'
#' @author Chrystian C. Sosa
#' 
#' 
create_folders_func <- function(basedir){
  if (!dir.exists(paste0(basedir, "/", "0.ProjectDocuments"))) {
    dir.create(paste0(basedir, "/", "0.ProjectDocuments")) }
  if (!dir.exists(paste0(basedir, "/", "0.ProjectDocuments"))) {
    dir.create(paste0(basedir, "/", "1.Processing")) }
  if (!dir.exists(paste0(basedir, "/", "1.Processing"))) {
    dir.create(paste0(basedir, "/", "1.Processing")) }
  if (!dir.exists(paste0(basedir, "/", "2.Publications"))) {
    dir.create(paste0(basedir, "/", "2.Publications")) }
  if (!dir.exists(paste0(basedir, "/", "3.Bibliography"))) {
    dir.create(paste0(basedir, "/", "3.Bibliography")) }
  if (!dir.exists(paste0(basedir, "/", "4.Scripts"))) {
    dir.create(paste0(basedir, "/", "4.Scripts")) }
  if (!dir.exists(paste0(basedir, "/", "5.Others"))) {
    dir.create(paste0(basedir, "/", "5.Others")) }
  if (!dir.exists(paste0(basedir, "/", "1.Data"))) {
    dir.create(paste0(basedir, "/", "1.Data")) }
  if (!dir.exists(paste0(basedir, "/", "1.Data/Results"))) {
    dir.create(paste0(basedir, "/", "1.Data/Results")) }
    
}

