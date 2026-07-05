#' Process Species SDM Consensus and Change Maps for Future Climate Scenarios
#'
#' @description
#' Builds consensus presence-absence rasters across multiple General Circulation
#' Models (GCMs) for each species and Shared Socioeconomic Pathway (SSP),
#' computes gain/loss change maps relative to current distributions, calculates
#' area statistics per change category, and optionally masks future consensus
#' maps with species-specific concave hull rasters. The function is designed to
#' run inside a \code{parallel} cluster via \code{\link[pbapply]{pblapply}}.
#'
#' @param i Integer. Row index into \code{species_set} identifying the species
#'   to process in the current worker.
#'
#' @return Called for its side effects. For each species \eqn{\times} SSP
#'   combination the function writes the following files to disk:
#'   \describe{
#'     \item{\code{output_dir_fut_NOHULL/<species>_<SSP>.tif}}{Raw GCM vote
#'       sum raster (integer 0–\emph{n} GCMs).}
#'     \item{\code{output_dir_fut/<species>_<SSP>.tif}}{Binary consensus
#'       presence-absence raster (majority rule: \eqn{\geq 3} GCMs agreeing
#'       = presence).}
#'     \item{\code{output_changes_dir/<species>_<SSP>.tif}}{Change raster
#'       encoding gain (\code{+2}), loss (\code{-1}), and stable presence
#'       (\code{NA} masked out).}
#'     \item{\code{output_changes_csv_dir/<species>_<SSP>_areas.csv}}{Area
#'       table (km²) per change value with species, year, and SSP metadata
#'       columns appended.}
#'     \item{\code{output_dir_fut_HULL/<species>_<SSP>.tif}}{Concave-hull-
#'       masked future consensus raster (only written when a hull file exists).}
#'     \item{\code{output_dir_logs/<species>.txt}}{Empty sentinel log file
#'       created upon successful completion.}
#'   }
#'   Returns \code{NULL} invisibly.
#'
#' @details
#' \subsection{Consensus logic}{
#'   Binary projections from each GCM are summed into a vote raster
#'   (\code{p_sum}). Cells with \eqn{< 3} votes are reclassified to 0
#'   (absence); cells with \eqn{\geq 3} votes are reclassified to 1
#'   (presence). This majority threshold assumes five GCMs are available;
#'   adjust the hard-coded threshold (\code{3}) if the GCM set changes.
#' }
#'
#' \subsection{Change detection}{
#'   The binary consensus is recoded to 2 (future presence) before
#'   subtracting the current binary SDM, yielding three interpretable values:
#'   \itemize{
#'     \item \code{+2}: gain (present in future, absent currently).
#'     \item \code{-1}: loss (absent in future, present currently).
#'     \item \code{0}: stable (both present or both absent) — set to
#'       \code{NA} and excluded from output.
#'   }
#'   The current SDM is resampled to the future grid with nearest-neighbour
#'   interpolation before subtraction.
#' }
#'
#' \subsection{Concave hull masking}{
#'   When a species-specific concave hull raster exists in \code{chull_dir},
#'   the raw vote-sum raster (before binary reclassification) is multiplied
#'   by the hull mask to restrict projections to the known accessible area.
#' }
#'
#' \subsection{Parallelisation}{
#'   The function relies on the following objects being exported to each
#'   worker via \code{\link[parallel]{clusterExport}} before
#'   \code{\link[pbapply]{pblapply}} is called:
#'   \code{species_set}, \code{basedir}, \code{current_dir}, \code{chull_dir},
#'   \code{GCMs}, \code{SSPs}, \code{year}, \code{output_dir_fut},
#'   \code{output_dir_fut_NOHULL}, \code{output_dir_fut_HULL},
#'   \code{output_changes_dir}, \code{output_changes_csv_dir},
#'   \code{output_dir_logs}.
#' }
#'
#' @section Global variables (must be defined in the calling environment
#'   and exported to the cluster):
#' \describe{
#'   \item{\code{species_set}}{Character vector of species names.}
#'   \item{\code{basedir}}{Root path to the project results directory.}
#'   \item{\code{current_dir}}{Path to the folder containing current
#'     binary SDM \code{.tif} files named \code{<species>.tif}.}
#'   \item{\code{chull_dir}}{Path to the folder containing concave hull
#'     rasters named \code{<species>_conc.tif}.}
#'   \item{\code{GCMs}}{Character vector of GCM identifiers (e.g.,
#'     \code{"ACCESS-CM2"}).}
#'   \item{\code{SSPs}}{Character vector of SSP labels (e.g.,
#'     \code{"ssp245"}, \code{"ssp370"}).}
#'   \item{\code{year}}{Character scalar for the target projection year
#'     (e.g., \code{"2050"}).}
#'   \item{\code{output_dir_fut}}{Output path for binary consensus rasters.}
#'   \item{\code{output_dir_fut_NOHULL}}{Output path for raw vote-sum rasters.}
#'   \item{\code{output_dir_fut_HULL}}{Output path for hull-masked rasters.}
#'   \item{\code{output_changes_dir}}{Output path for change rasters.}
#'   \item{\code{output_changes_csv_dir}}{Output path for area CSV files.}
#'   \item{\code{output_dir_logs}}{Output path for per-species log files.}
#' }
#'
#' @note
#' \itemize{
#'   \item \pkg{terra} is loaded explicitly inside the function body to ensure
#'     availability on each parallel worker.
#'   \item All raster writes use \code{overwrite = TRUE}; existing files are
#'     silently replaced.
#'   \item \code{setwd()} calls inside the function change the working
#'     directory of the worker process; consider replacing them with full
#'     absolute paths in \code{writeRaster} and \code{write.csv} calls for
#'     safer parallel execution.
#'   \item Species with no current SDM file are skipped with a message;
#'     SSP iterations with fewer than one valid GCM projection are also
#'     skipped.
#' }
#'
#' @seealso
#' \code{\link[terra]{rast}}, \code{\link[terra]{resample}},
#' \code{\link[terra]{expanse}}, \code{\link[terra]{writeRaster}},
#' \code{\link[parallel]{makeCluster}},
#' \code{\link[pbapply]{pblapply}}
#'
#' @importFrom terra rast resample expanse writeRaster
#'
#' @examples
#' \dontrun{
#' library(parallel)
#' library(pbapply)
#'
#' # Define all required global objects (basedir, species_set, GCMs, etc.)
#' # then launch the parallel cluster:
#'
#' cl <- makeCluster(8)
#' clusterExport(cl, varlist = c(
#'   "species_set", "basedir", "current_dir", "chull_dir",
#'   "GCMs", "SSPs", "year",
#'   "output_dir_fut", "output_dir_fut_NOHULL", "output_dir_fut_HULL",
#'   "output_changes_dir", "output_changes_csv_dir", "output_dir_logs"
#' ))
#'
#' pblapply(seq_along(species_set), process_species, cl = cl)
#'
#' stopCluster(cl)
#' }
#'
#' @export


# load packages
library(terra)
library(parallel)
library(pbapply)

# set user of the script
user <- 'Chrystian'

if (user == "Chrystian") {
  basedir <- "//catalogue/MultifLandscapesA1706/1.Data/Results"
}

setwd(dir = paste0(basedir, "/SDM results"))
species_set <- read.csv("species_set_SDM_all.csv")
species_set <- species_set[, 1]

GCMs <- c("ACCESS-CM2", "GISS-E2-1-G", "INM-CM5-0", "MIROC6", "MPI-ESM1-2-HR")
SSPs <- c("ssp245", "ssp370")
year <- "2050"

output_dir             <- paste0(basedir, "/SDM results/", "Distribution maps/Future/", year)
output_dir             <- paste0(output_dir, "/", "Consensus_maps_folder2")
output_dir_fut         <- paste0(output_dir, "/", "Fut_consensus")
output_dir_fut_NOHULL  <- paste0(output_dir_fut, "/", "NO_HULL")
output_dir_fut_HULL    <- paste0(output_dir_fut, "/", "CONCAVE_HULL")
output_changes_dir     <- paste0(output_dir, "/", "changes_FC")
output_changes_csv_dir <- paste0(output_dir, "/", "changes_CSV")
output_SUM_dir         <- paste0(output_dir, "/", "sums")

if (!dir.exists(output_dir))             dir.create(output_dir)
if (!dir.exists(output_changes_dir))     dir.create(output_changes_dir)
if (!dir.exists(output_SUM_dir))         dir.create(output_SUM_dir)
if (!dir.exists(output_dir_fut))         dir.create(output_dir_fut)
if (!dir.exists(output_changes_csv_dir)) dir.create(output_changes_csv_dir)
if (!dir.exists(output_dir_fut_NOHULL))  dir.create(output_dir_fut_NOHULL)
if (!dir.exists(output_dir_fut_HULL))    dir.create(output_dir_fut_HULL)

current_dir <- paste0(basedir, "/SDM results/", "Distribution maps", "/", "Presence-absence")
chull_dir   <- "//catalogue/MultifLandscapesA1706/1.Data/Results/concave_hull_rasters_ADM0_2_5m"


process_species <- function(i) {
  library(terra)
  
  print(i)
  print(species_set[i])
  
  # ---- Load current binary SDM ----
  current_sdm_file <- paste0(current_dir, "/", species_set[[i]], ".tif")
  if (file.exists(current_sdm_file)) {
    current_sdm <- terra::rast(current_sdm_file)
  } else {
    current_sdm <- NULL
  }
  
  # ---- Load concave hull raster ----
  c_hull_file <- paste0(chull_dir, "/", species_set[[i]], "_conc.tif")
  if (file.exists(c_hull_file)) {
    c_hull_terra <- terra::rast(c_hull_file)
  } else {
    c_hull_terra <- NULL
  }
  
  if (!is.null(current_sdm)) {
    
    for (s in seq_along(SSPs)) {
      print(SSPs[s])
      p_list <- list()
      
      # ---- Stack GCM projections ----
      for (g in seq_along(GCMs)) {
        GCM_dir <- paste0(basedir, "/SDM results/Distribution maps/Future/",
                          year, "/", GCMs[g], "/", SSPs[s])
        print(GCMs[g])
        
        gcm_file <- paste0(GCM_dir, "/", species_set[i], "_", year, ".tif")
        if (file.exists(gcm_file)) {
          p_list[[g]] <- terra::rast(gcm_file)
        } else {
          p_list[[g]] <- NULL
        }
      }
      
      p_list <- p_list[!sapply(p_list, is.null)]
      
      if (length(p_list) > 0) {
        p     <- terra::rast(p_list)
        message(paste0("summarizing projections for: ", year, "/", GCMs[g], "/", SSPs[s]))
        
        # ---- Consensus: vote sum ----
        p_sum <- sum(p, na.rm = TRUE)
        writeRaster(p_sum,
                    paste0(output_dir_fut_NOHULL, "/", species_set[i], "_", SSPs[s], ".tif"),
                    overwrite = TRUE)
        
        # ---- Binary reclassification (majority threshold = 3) ----
        p_sum[p_sum[] < 3] <- 0
        p_sum[p_sum[] > 0] <- 1
        writeRaster(p_sum,
                    paste0(output_dir_fut, "/", species_set[i], "_", SSPs[s], ".tif"),
                    overwrite = TRUE)
        
        # ---- Change detection ----
        message(paste0("Calculating changes for: ",
                       species_set[i], "/", year, "/", GCMs[g], "/", SSPs[s]))
        
        p_sum2          <- p_sum
        p_sum2[p_sum2 == 1] <- 2            # recode future presence to 2
        
        current_sdm <- terra::resample(current_sdm, p_sum2, "near")
        changes     <- p_sum2 - current_sdm  # +2 = gain, -1 = loss, 0 = stable
        changes[changes == 0] <- NA
        
        # ---- Area statistics ----
        x_expanse        <- terra::expanse(changes, unit = "km", byValue = TRUE, wide = FALSE)
        x_expanse$layer  <- species_set[i]
        x_expanse$year   <- year
        x_expanse$SSP    <- SSPs[s]
        
        writeRaster(changes,
                    paste0(output_changes_dir, "/", species_set[i], "_", SSPs[s], ".tif"),
                    overwrite = TRUE)
        write.csv(x_expanse,
                  paste0(output_changes_csv_dir, "/",
                         species_set[i], "_", SSPs[s], "_areas.csv"),
                  row.names = FALSE)
        
        # ---- Concave hull masking ----
        if (!is.null(c_hull_terra)) {
          p_sum_raw  <- terra::rast(paste0(output_dir_fut_NOHULL, "/",
                                           species_set[i], "_", SSPs[s], ".tif"))
          c_hull_res <- terra::resample(c_hull_terra, p_sum_raw, "near")
          p_sum_CH   <- p_sum_raw * c_hull_res
          writeRaster(p_sum_CH,
                      paste0(output_dir_fut_HULL, "/",
                             species_set[i], "_", SSPs[s], ".tif"),
                      overwrite = TRUE)
        } else {
          message("No CH available!")
        }
        
      } else {
        message("no files to do changes")
      }
    }
    
    # ---- Sentinel log file ----
    file.create(paste0(output_dir_logs, "/", species_set[[i]], ".txt"))
    
  } else {
    message("NO SDM! SKIPPING")
  }
}


# ---------------------------------------------------------------------------
# Parallel execution
# ---------------------------------------------------------------------------

#' @note
#' Adjust the number of cores passed to \code{makeCluster} according to
#' server availability. Do NOT use \code{detectCores() - 1} on shared
#' HPC nodes.

cl <- makeCluster(8)   # <-- confirm core count with server admin before running

clusterExport(cl, varlist = c(
  "species_set", "basedir", "current_dir", "chull_dir",
  "GCMs", "SSPs", "year",
  "output_dir_fut", "output_dir_fut_NOHULL", "output_dir_fut_HULL",
  "output_changes_dir", "output_changes_csv_dir", "output_dir_logs"
))

pbapply::pblapply(seq_along(species_set), process_species, cl = cl)

stopCluster(cl)