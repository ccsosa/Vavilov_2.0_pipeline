require(readxl)
require(data.table)
require(ggplot2)
require(ggpubr)
require(terra)
require(foreach)
require(tidyterra)
require(dplyr)


summmary_raster_var_toile <-
  function(basedir,
           resultdir,
           output_dir,
           species_selected,
           template_file,
           data_ori) {
    out_dir <- output_dir
    out_buffer_dir <- paste0(out_dir, "/soil_extreme")
    if (!dir.exists(out_buffer_dir)) dir.create(out_buffer_dir)
    
    out_buffer_dir_raster <- paste0(out_buffer_dir, "/raster")
    if (!dir.exists(out_buffer_dir_raster)) dir.create(out_buffer_dir_raster)
    
    out_buffer_dir_summary <- paste0(out_buffer_dir, "/summary2")
    if (!dir.exists(out_buffer_dir_summary)) dir.create(out_buffer_dir_summary)
    
    out_buffer_dir_summary_png  <- paste0(out_buffer_dir_summary, "/PNG")
    out_buffer_dir_summary_tif  <- paste0(out_buffer_dir_summary, "/TIF")
    out_buffer_dir_summary_html <- paste0(out_buffer_dir_summary, "/HTML")
    
    message(out_buffer_dir_summary)
    
    if (!dir.exists(out_buffer_dir_summary_png))  dir.create(out_buffer_dir_summary_png)
    if (!dir.exists(out_buffer_dir_summary_tif))  dir.create(out_buffer_dir_summary_tif)
    if (!dir.exists(out_buffer_dir_summary_html)) dir.create(out_buffer_dir_summary_html)
    
    out_buffer_dir_summary_chunks <- paste0(out_buffer_dir_summary, "/chunks")
    if (!dir.exists(out_buffer_dir_summary_chunks)) dir.create(out_buffer_dir_summary_chunks)
    
    ############################################################################
    data      <- data_ori[which(data_ori$prop_eco_total >= 25), ]
    variables <- unique(data$var)
    
    template <- terra::rast(template_file)
    template[which(!is.na(template[]))] <- 1
    
    message("loading shapefile and metrics")
    adm0 <- vect(paste0(basedir, "/1.Data/RAW/input_data/adm0/adm0_Latam.shp"))
    
    ecco <- sf::read_sf(paste0(basedir, "/1.Data/RAW/input_data/Ecoregions/biomes_LATAM_2017.shp"))
    sf::sf_use_s2(FALSE)
    
    species_total <- unique(species_selected$species)
    
    InDir_buffer     <- "/catalogue/MultifLandscapesA1706/1.Data/Results/buffer_richness/raster"
    sp_buffer_list   <- list.files(InDir_buffer, pattern = ".tif")
    sp_buffer_list_name <- sub(pattern = ".tif", replacement = "", sp_buffer_list)
    
    metrics      <- read.csv(paste0(resultdir, "/", "AUC_Maxent_valid_all.csv"))
    species_list <- metrics$species
    species_list <- sub(pattern = "model_evaluation_", replacement = "", species_list)
    species_list <- sub(pattern = ".csv", replacement = "", species_list)
    
    avail_sp_modeled <- species_list[species_list %in% species_total]
    
    buffer_sps <- species_total[!species_total %in% avail_sp_modeled]
    buffer_sps <- buffer_sps[buffer_sps %in% sp_buffer_list_name]
    
    InDir <- paste0(resultdir, "/Distribution maps/Presence-absence")
    
    ############################################################################
    # helper: parallel chunked sum → returns SpatRaster or NULL
    ############################################################################
    chunked_richness_sum <- function(valid_paths, chunk_dir, label = "") {
      
      message(paste0(length(valid_paths), " valid rasters found."))
      
      # --- Early exit ---
      if (length(valid_paths) == 0) {
        message("No valid rasters found, skipping.")
        return(NULL)
      }
      
      chunk_size <- if (length(valid_paths) < 10) length(valid_paths) else 10
      chunk_idx  <- split(seq_along(valid_paths),
                          ceiling(seq_along(valid_paths) / chunk_size))
      n_chunks   <- length(chunk_idx)
      message(paste0("Processing ", n_chunks, " chunks of up to ", chunk_size, " rasters..."))
      
      n_cores <- 8
      cl      <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      message(paste0("Using ", n_cores, " cores for ", n_chunks, " chunks"))
      
      if (!dir.exists(chunk_dir)) dir.create(chunk_dir)
      
      chunk_results <- foreach::foreach(
        chunk          = chunk_idx,
        .packages      = "terra",
        .errorhandling = "pass"
      ) %dopar% {
        paths_chunk <- valid_paths[chunk]
        chunk_id    <- chunk[1]
        stack       <- terra::rast(paths_chunk)
        result      <- terra::app(stack, fun = "sum", na.rm = TRUE)
        out_path    <- paste0(chunk_dir, "/chunk_", sprintf("%04d", chunk_id), ".tif")
        terra::writeRaster(result, filename = out_path, overwrite = TRUE)
        out_path
      }
      
      parallel::stopCluster(cl)
      
      failed <- sapply(chunk_results, inherits, "error")
      if (any(failed)) {
        message(paste0(sum(failed), " chunk(s) failed:"))
        print(chunk_results[failed])
        stop("Aborting: fix failed chunks before proceeding.")
      }
      
      message("Combining chunk results into final richness raster...")
      chunk_paths <- unlist(chunk_results)
      chunk_stack <- terra::rast(chunk_paths)
      x           <- terra::app(chunk_stack, fun = "sum", na.rm = TRUE)
      
      message("Cleaning up chunk files...")
      file.remove(chunk_paths)
      unlink(chunk_dir)
      message("Done!")
      
      return(x)
    }
    
    ############################################################################
    for (i in 1:length(variables)) {
      variable <- variables[[i]]
      
      for (j in 1:1) {
        var_to     <- "Q_10_90_eco"
        data_var_i <- data[which(data$var == variable), ]
        data_var_i_var_toile <-
          data_var_i$species_searched[which(data_var_i[, var_to] == TRUE)]
        
        avail_sp_modeled_i <- data_var_i_var_toile[data_var_i_var_toile %in% avail_sp_modeled]
        buffer_sps_i       <- data_var_i_var_toile[data_var_i_var_toile %in% buffer_sps]
        
        ########################################################################
        # SDM richness
        ########################################################################
        sdm_tif <- paste0(out_buffer_dir_summary_tif, "/", variable, "_", var_to, "_SDM_richness_current.tif")
        
        if (!file.exists(sdm_tif)) {
          message("Checking SDM file availability...")
          valid_paths <- vapply(avail_sp_modeled_i, function(sp) {
            path <- paste0(InDir, "/", sp, ".tif")
            if (file.exists(path)) path else NA_character_
          }, character(1))
          missing_sps <- avail_sp_modeled_i[is.na(valid_paths)]
          valid_paths <- valid_paths[!is.na(valid_paths)]
          if (length(missing_sps) > 0)
            message(paste0(length(missing_sps), " species raster(s) not found and skipped."))
          
          richness <- chunked_richness_sum(
            valid_paths = valid_paths,
            chunk_dir   = paste0(out_buffer_dir_summary_chunks, "/chunks2")
          )
          
          if (!is.null(richness)) {
            message("Saving SDM richness raster...")
            terra::writeRaster(richness, filename = sdm_tif, overwrite = TRUE)
          }
          
        } else {
          richness <- terra::rast(sdm_tif)
        }
        
        ########################################################################
        # Buffer richness
        ########################################################################
        buf_tif <- paste0(out_buffer_dir_summary_tif, "/", variable, "_", var_to, "_buffer_richness_current.tif")
        
        if (!file.exists(buf_tif)) {
          message("Checking buffer file availability...")
          valid_paths <- vapply(buffer_sps_i, function(sp) {
            path <- paste0(InDir_buffer, "/", sp, ".tif")
            if (file.exists(path)) path else NA_character_
          }, character(1))
          missing_sps <- buffer_sps_i[is.na(valid_paths)]
          valid_paths <- valid_paths[!is.na(valid_paths)]
          if (length(missing_sps) > 0)
            message(paste0(length(missing_sps), " species raster(s) not found and skipped."))
          
          richness_buffer <- chunked_richness_sum(
            valid_paths = valid_paths,
            chunk_dir   = paste0(out_buffer_dir_summary_tif, "/chunks2")
          )
          
          if (!is.null(richness_buffer)) {
            message("Saving buffer richness raster...")
            terra::writeRaster(richness_buffer, filename = buf_tif, overwrite = TRUE)
          }
          
        } else {
          richness_buffer <- terra::rast(buf_tif)
        }
        
        ########################################################################
        # Combine SDM + buffer
        ########################################################################
        richness_tif <- paste0(out_buffer_dir_summary_tif, "/", variable, "_", var_to, "_richness_current.tif")
        
        if (!is.null(richness) & !is.null(richness_buffer)) {
          message("SDM AVAIL, BUFFER AVAIL")
          richness_buffer_res <- terra::resample(richness_buffer, template, method = "near")
          richness_total      <- sum(richness, richness_buffer_res, na.rm = TRUE)
        } else if (is.null(richness) & !is.null(richness_buffer)) {
          message("NO SDM, BUFFER AVAIL")
          richness_total <- terra::resample(richness_buffer, template, method = "near")
        } else if (!is.null(richness) & is.null(richness_buffer)) {
          message("SDM AVAIL, NO BUFFER")
          richness_total <- richness
        } else {
          message("NO FILE TO MAP")
          richness_total <- NULL
        }
        
        if (is.null(richness_total)) next
        
        terra::writeRaster(richness_total, filename = richness_tif, overwrite = TRUE)
        richness_total[which(richness_total[] == 0)] <- NA
        
        ########################################################################
        # Plot
        ########################################################################
        message("Making a plot")
        
        biome_centroids <- sf::st_point_on_surface(ecco)
        biome_coords_df <- data.frame(sf::st_coordinates(biome_centroids),
                                      BIOME_NAME = ecco$BIOME)
        biome_coords_df <- biome_coords_df[
          !is.na(biome_coords_df$BIOME_NAME) & biome_coords_df$BIOME_NAME != "N/A", ]
        
        p <- ggplot() +
          scale_fill_gradientn(
            colours  = rev(RColorBrewer::brewer.pal(11, "Spectral")),
            na.value = NA,
            name     = "Species\nRichness"
          ) +
          tidyterra::geom_spatraster(data = richness_total, na.rm = TRUE) +
          geom_sf(data = sf::st_as_sf(adm0), fill = NA, color = "gray", linewidth = 0.5) +
          geom_sf(data = ecco, fill = NA, color = "grey17", linewidth = 0.8) +
          ggrepel::geom_label_repel(
            data           = biome_coords_df,
            aes(x = X, y = Y, label = BIOME_NAME),
            color          = "black",
            na.rm          = TRUE,
            fill           = "gray86",
            size           = 2.1,
            fontface       = "bold",
            box.padding    = 0.5,
            point.padding  = 0.3,
            max.iter       = 50000,
            max.overlaps   = Inf,
            force          = 100,
            force_pull     = 1,
            nudge_x        = ifelse(biome_coords_df$X < -75, -25, 20),
            nudge_y        = ifelse(biome_coords_df$Y > 10, 15, -15),
            arrow          = arrow(length = unit(0.015, "npc"), type = "closed", ends = "last"),
            segment.color  = "grey30",
            segment.size   = 0.4,
            segment.curvature  = 0.2,
            min.segment.length = 0
          ) +
          coord_sf(expand = FALSE) +
          theme_void() +
          theme(
            legend.position = "right",
            legend.title    = element_text(color = "black", size = 12),
            plot.background = element_rect(fill = "white", color = NA)
          )
        
        ggsave(
          paste0(out_buffer_dir_summary_png, "/LATAM_Edible_food_", variable, "_", var_to, "_species_richness_current.png"),
          plot   = p,
          width  = 14,
          height = 14,
          dpi    = 300
        )
      }
    }
    
    return("DONE!")
  }

################################################################################
output_dir <- "/catalogue/MultifLandscapesA1706/1.Data/Results"
user       <- "Chrystian"
if (user == "Chrystian") basedir <- "/catalogue/MultifLandscapesA1706"

resultdir <- "/catalogue/MultifLandscapesA1706/1.Data/Results/SDM results"

species_selected <- as.data.frame(
  readxl::read_xlsx(
    "/catalogue/MultifLandscapesA1706/1.Data/Results/species_lists/species_curated/Species_20250304_edible_part_curated.xlsx",
    col_names = TRUE
  )
)

template_file <- paste0(
  "/catalogue/MultifLandscapesA1706/1.Data/RAW/Input_data/climate_data/2_5min/present/",
  "wc2.1_2.5m_bio_1.tif"
)

data_ori <- read.csv(
  "/catalogue/MultifLandscapesA1706/1.Data/Process/soil_extreme/outliers_list.csv"
)

x <- summmary_raster_var_toile(basedir, resultdir, output_dir, species_selected, template_file, data_ori)