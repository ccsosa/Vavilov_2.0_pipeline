require(readxl)
require(data.table)
require(ggplot2)
require(ggpubr)
require(terra)
require(foreach)
require(tidyterra)
require(dplyr)

################################################################################
# helper: run chunked parallel sum and return SpatRaster or NULL
chunked_sum <- function(valid_paths, chunk_dir, n_cores = 8) {
  if (length(valid_paths) == 0) {
    message("No valid rasters found, skipping.")
    return(NULL)
  }
  
  chunk_size <-
    if (length(valid_paths) < 10)
      length(valid_paths)
  else
    10
  chunk_idx  <- split(seq_along(valid_paths),
                      ceiling(seq_along(valid_paths) / chunk_size))
  n_chunks   <- length(chunk_idx)
  message(paste0(
    "Processing ",
    n_chunks,
    " chunks of up to ",
    chunk_size,
    " rasters..."
  ))
  
  if (!dir.exists(chunk_dir))
    dir.create(chunk_dir, recursive = TRUE)
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  message(paste0("Using ", n_cores, " cores for ", n_chunks, " chunks"))
  
  chunk_results <- foreach::foreach(chunk          = chunk_idx,
                                    .packages      = "terra",
                                    .errorhandling = "pass") %dopar% {
                                      paths_chunk <- valid_paths[chunk]
                                      chunk_id    <- chunk[1]
                                      stack       <-
                                        terra::rast(paths_chunk)
                                      result      <-
                                        terra::app(stack, fun = "sum", na.rm = TRUE)
                                      out_path    <-
                                        paste0(chunk_dir, "/chunk_", sprintf("%04d", chunk_id), ".tif")
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
################################################################################
summmary_raster_var_toile_ecosys <-
  function(basedir,
           resultdir,
           output_dir,
           species_selected,
           template_file,
           data_ori,
           conversion_table) {
    out_dir <- output_dir
    out_buffer_dir <- paste0(out_dir, "/soil_extreme")
    if (!dir.exists(out_buffer_dir))
      dir.create(out_buffer_dir)
    
    out_buffer_dir <- paste0(out_dir, "/soil_extreme/eco_approach")
    if (!dir.exists(out_buffer_dir))
      dir.create(out_buffer_dir)
    
    # out_buffer_dir_raster <- paste0(out_buffer_dir, "/raster_files")
    # if (!dir.exists(out_buffer_dir_raster))
    #   dir.create(out_buffer_dir_raster)
    
    out_buffer_dir_summary <-
      paste0(out_buffer_dir, "/summary_files_potential")
    if (!dir.exists(out_buffer_dir_summary))
      dir.create(out_buffer_dir_summary)
    
    out_buffer_dir_summary_png   <-
      paste0(out_buffer_dir_summary, "/PNG")
    out_buffer_dir_summary_tif   <-
      paste0(out_buffer_dir_summary, "/TIF")
    # out_buffer_dir_summary_html  <-
    #   paste0(out_buffer_dir_summary, "/HTML")
    #
    message(out_buffer_dir_summary)
    
    if (!dir.exists(out_buffer_dir_summary_png))
      dir.create(out_buffer_dir_summary_png)
    if (!dir.exists(out_buffer_dir_summary_tif))
      dir.create(out_buffer_dir_summary_tif)
    # if (!dir.exists(out_buffer_dir_summary_html))
    #   dir.create(out_buffer_dir_summary_html)
    
    out_buffer_dir_summary_chunks <-
      paste0(out_buffer_dir_summary, "/chunks")
    if (!dir.exists(out_buffer_dir_summary_chunks))
      dir.create(out_buffer_dir_summary_chunks)
    
    ############################################################################
    #at least 20 points per ecosystems
    data      <- data_ori[which(data_ori$prop_eco_total >= 20), ]
    #get the variables to map
    variables <- unique(data$var)
    ############################################################################
    #using a template and adding 1 as vaue for the target region
    template <- terra::rast(template_file)
    template[which(!is.na(template[]))] <- 1
    
    ############################################################################
    #loading LATAM shapefile
    message("loading shapefile and metrics")
    adm0 <-
      vect(paste0(basedir, "/1.Data/RAW/input_data/adm0/adm0_Latam.shp"))
    #loading LATAM biomes
    ecco <-
      sf::read_sf(paste0(
        basedir,
        "/1.Data/RAW/input_data/Ecoregions/biomes_LATAM_2017.shp"
      ))
    #avoiding any topology issue
    sf::sf_use_s2(FALSE)
    
    #getting all unique species
    species_total <- unique(species_selected$species)
    ##############################################################################
    #loading buffer species
    InDir_buffer <-
      "//catalogue/MultifLandscapesA1706/1.Data/Results/buffer_richness/raster"
    #species with buffer
    sp_buffer_list <- list.files(InDir_buffer, pattern = ".tif")
    #getting specices names for buffer
    sp_buffer_list_name <-
      sub(pattern = ".tif", replacement = "", sp_buffer_list)
    ##############################################################################
    #adding presence absence results folder
    InDir <-
      paste0(resultdir, "/Distribution maps/Presence-absence")
    
    #getting all possible species with results
    sp_modeled_availa <- list.files(InDir, pattern = ".tif")
    sp_modeled_availa <- sub(".tif", "", sp_modeled_availa)
    ##############################################################################
    #load metrics to get modeled species
    metrics <-
      read.csv(paste0(resultdir, "/", "AUC_Maxent_valid_all.csv"))
    #gettting species lists to get results (modeled species with MaxEnt)
    species_list <- metrics$species
    species_list <-
      sub(pattern = "model_evaluation_", replacement = "", species_list)
    species_list <-
      sub(pattern = ".csv", replacement = "", species_list)
    # print(length(species_list))
    # species_list <- species_list[which(species_list!="Cabralea canjerana")]
    # print(length(species_list))
    ################################################################################
    #species with valid model in the folder
    avail_sp_modeled <-
      species_list[species_list %in% sp_modeled_availa]
    ################################################################################
    #fixing names
    #geting species names accepted
    accp_sp <-
      conversion_table[conversion_table$species %in% avail_sp_modeled, ]
    #getting synonym!
    syn_sp <-
      conversion_table[conversion_table$species_from_source %in% avail_sp_modeled, ]
    #only removing accepting to get the synonyms and accepted names
    syn_sp <- syn_sp[which(syn_sp$taxonomic_status != "Accepted"), ]
    #getting available
    avail_sp_modeled <- c(accp_sp$species, syn_sp$species)
    #only getting unique species
    avail_sp_modeled <- unique(avail_sp_modeled)
    ##############################################################################
    #species with buffer
    # buffer_sps <- species_total[sp_buffer_list_name %in% species_total]
    buffer_sps <-
      species_total[!species_total %in% avail_sp_modeled]
    #if buffer_sps + avail_sp_modeled =8221 is well done!
    buffer_sps <-
      buffer_sps[buffer_sps %in% sp_buffer_list_name] #available to use buffer
    
    buffer_sps <- unique(buffer_sps)
    ################################################################################
    #getting unique biomes to map
    biomes <- unique(data$biome)
    #quantiles
    var_to <- c(
      "Q50sp_Q75eco",
      "Q75sp_Q50eco",
      "Q5sp_Q95eco",
      "Q95sp_Q5eco",
      "Q10sp_Q90eco",
      "Q90sp_Q10eco",
      "Q25sp_Q50eco",
      "Q50sp_Q25eco"
    )
    ############################################################################
    # i <- 1
    # j <- 1
    # k <- 1
    # l <- 1
    ############################################################################
    
    #for each variable
    for (i in 1:length(variables)) {
      print(paste0("i: ", i))
      variable <- variables[[i]]
      print(variable)
      #for each biome
      for (j in 1:length(biomes)) {
        print(paste0("j: ", j))
        biome <- biomes[[j]]
        print(biome)
        ecco_biome_i <- ecco[which(ecco$BIOME ==  biomes[[j]]), ]
        
        # k <- 1
        #subset variable and biome
        for (k in 1:length(var_to)) {
          print(paste0("k: ", k))
          quantile_res <- var_to[[k]]
          print(quantile_res)
          #subset by biome
          data_var_i <-
            data[which(data$var == variable &
                         data$biome == biomes[[j]]), ]
          #getting unique values per variable
          subvar <- unique(data_var_i[, var_to[[k]]])
          
          # data_var_i <- data_var_i$species_searched[which(data_var_i[, var_to] == TRUE)]
          # l <- 1
          #subset by variable unique values
          for (l in 1:length(subvar)) {
            print(paste0("l: ", l))
            #subset
            subvar_l <- subvar[[l]]
            print(subvar_l)
            #subsetting by variable, biome, and subset of variable
            data_subvar_i <-
              data_var_i$species_searched[which(data_var_i[, quantile_res] == subvar[[l]])]
            if (length(data_subvar_i) == 0) {
              message("NO RESULTS")
              next
            }
            
            #getting avaiable models and buffers
            avail_sp_modeled_i <-
              data_subvar_i[data_subvar_i %in% avail_sp_modeled]
            buffer_sps_i       <-
              data_subvar_i[data_subvar_i %in% buffer_sps]
            
            ##########################################################################
            # SDM richness
            #removing () and spaces
            name <- subvar_l
            name <- sub("<", "_LESS_", name, fixed = T)
            name <- sub(">", "_GREATER_", name, fixed = T)
            name <-
              sub("(All spp in ecosystem)", "", name, fixed = T)
            # name <- sub("\\(\\*","",name,fixed = T)
            name <- gsub("[()]", "", name)
            name <- gsub("_+$", "", name)
            name <- gsub(" ", "", name)
            name <- gsub("&", "_and_", name)
            biome_name <- gsub("&", "and", biomes[[j]])
            biome_name <- gsub(" ", "_", biome_name)
            biome_name <- gsub(",", "", biome_name)
            ####################################################################
            
            
            if (!dir.exists(
              paste0(
                out_buffer_dir_summary_tif,
                "/",
                variable,
                "/",
                biome_name,
                "/",
                quantile_res
              )
            )) {
              dir.create(
                paste0(
                  out_buffer_dir_summary_tif,
                  "/",
                  variable,
                  "/",
                  biome_name,
                  "/",
                  quantile_res
                ),
                recursive = T
              )
            }
            dir_1 <-  paste0(
              out_buffer_dir_summary_tif,
              "/",
              variable,
              "/",
              biome_name,
              "/",
              quantile_res
            )
            ####################################################################
            if (!dir.exists(paste0(dir_1,
                                   "/",
                                   name,
                                   "/",
                                   "BUFFER_SDM"))) {
              dir.create(paste0(dir_1,
                                "/",
                                name,
                                "/",
                                "BUFFER_SDM"),
                         recursive = T)
            }
            
            sdm_tif <- paste0(dir_1,
                              "/",
                              name,
                              "/",
                              "BUFFER_SDM",
                              "/",
                              "SDM_richness_current.tif")
            
            if (!file.exists(sdm_tif)) {
              message("Checking SDM file availability...")
              valid_paths <-
                vapply(avail_sp_modeled_i, function(sp) {
                  path <- paste0(InDir, "/", sp, ".tif")
                  if (file.exists(path))
                    path
                  else
                    NA_character_
                }, character(1))
              
              missing_sps <- avail_sp_modeled_i[is.na(valid_paths)]
              valid_paths <- valid_paths[!is.na(valid_paths)]
              if (length(missing_sps) > 0)
                message(paste0(
                  length(missing_sps),
                  " species raster(s) not found and skipped."
                ))
              message(paste0(length(valid_paths), " valid rasters found."))
              
              richness <- chunked_sum(
                valid_paths,
                chunk_dir = paste0(out_buffer_dir_summary_chunks, "/chunks2"),
                n_cores   = 8
              )
              
              if (!is.null(richness)) {
                message("Saving SDM richness raster...")
                terra::writeRaster(richness,
                                   filename = sdm_tif,
                                   overwrite = TRUE)
              }
            } else {
              richness <- terra::rast(sdm_tif)
            }
            
            ##########################################################################
            # Buffer richness
            
            buf_tif <-
              paste0(dir_1,
                     "/",
                     name,
                     "/",
                     "BUFFER_SDM",
                     "/",
                     "buffer_richness_current.tif")
            
            if (!file.exists(buf_tif)) {
              message("Checking buffer file availability...")
              valid_paths <- vapply(buffer_sps_i, function(sp) {
                path <- paste0(InDir_buffer, "/", sp, ".tif")
                if (file.exists(path))
                  path
                else
                  NA_character_
              }, character(1))
              
              missing_sps <- buffer_sps_i[is.na(valid_paths)]
              valid_paths <- valid_paths[!is.na(valid_paths)]
              if (length(missing_sps) > 0)
                message(paste0(
                  length(missing_sps),
                  " species raster(s) not found and skipped."
                ))
              message(paste0(length(valid_paths), " valid rasters found."))
              
              richness_buffer <- chunked_sum(
                valid_paths,
                chunk_dir = paste0(
                  out_buffer_dir_summary_chunks,
                  "/chunks2_buffer"
                ),
                n_cores   = 8
              )
              
              if (!is.null(richness_buffer)) {
                message("Saving buffer richness raster...")
                terra::writeRaster(richness_buffer,
                                   filename = buf_tif,
                                   overwrite = TRUE)
              }
            } else {
              richness_buffer <- terra::rast(buf_tif)
            }
            
            ##########################################################################
            # Combine SDM + buffer
            richness_res <-
              richness  # no resample needed per original logic
            
            total_tif <- paste0(
              out_buffer_dir_summary_tif,
              "/",
              variable,
              "/",
              biome_name,
              "/",
              quantile_res,
              "/",
              name,
              "_",
              "richness_current.tif"
            )
            
            if (!is.null(richness_res) &&
                !is.null(richness_buffer)) {
              message("SDM AVAIL, BUFFER AVAIL")
              richness_res <-
                terra::resample(richness_res, template, method = "near")
              richness_buffer_res <-
                terra::resample(richness_buffer, template, method = "near")
              richness_total      <-
                sum(richness_res, richness_buffer_res, na.rm = TRUE)
              terra::writeRaster(richness_total,
                                 filename = total_tif,
                                 overwrite = TRUE)
              
            } else if (is.null(richness_res) &&
                       !is.null(richness_buffer)) {
              message("NO SDM, BUFFER AVAIL")
              richness_total <-
                terra::resample(richness_buffer, template, method = "near")
              terra::writeRaster(richness_total,
                                 filename = total_tif,
                                 overwrite = TRUE)
              
            } else if (!is.null(richness_res) &&
                       is.null(richness_buffer)) {
              message("SDM AVAIL, NO BUFFER")
              richness_total <- richness_res
              richness_total <-
                terra::resample(richness_res, template, method = "near")
              terra::writeRaster(richness_total,
                                 filename = total_tif,
                                 overwrite = TRUE)
              
            } else {
              message("NO FILE TO MAP")
              richness_total <- NULL
            }
            
            if (is.null(richness_total))
              next
            
            richness_total[which(richness_total[] == 0)] <- NA
            
            ##########################################################################
            
            message("Making a plot")
            
            richness_total_masked <-
              terra::mask(richness_total, terra::vect(ecco_biome_i))
            biome_bbox <- sf::st_bbox(ecco_biome_i)
            biome_centroids <- sf::st_point_on_surface(ecco_biome_i)
            biome_coords_df <-
              data.frame(sf::st_coordinates(biome_centroids),
                         BIOME_NAME = ecco_biome_i$BIOME)
            biome_coords_df <-
              biome_coords_df[!is.na(biome_coords_df$BIOME_NAME) &
                                biome_coords_df$BIOME_NAME != "N/A", ]
            
            
            if (!dir.exists(
              paste0(
                out_buffer_dir_summary_png,
                "/",
                variable,
                "/",
                biome_name,
                "/",
                quantile_res,
                "/",
                name
              )
            )) {
              dir.create(
                paste0(
                  out_buffer_dir_summary_png,
                  "/",
                  variable,
                  "/",
                  biome_name,
                  "/",
                  quantile_res,
                  "/",
                  name
                ),
                recursive = T
              )
            }
            
            
            p <- ggplot() +
              scale_fill_gradientn(
                colours  = rev(RColorBrewer::brewer.pal(11, "Spectral")),
                na.value = NA,
                name     = "Species\nRichness"
              ) +
              tidyterra::geom_spatraster(data = richness_total_masked, na.rm = TRUE) +
              geom_sf(
                data = sf::st_as_sf(adm0),
                fill = NA,
                color = "grey50",
                linewidth = 0.5
              ) +
              geom_sf(
                data = ecco_biome_i,
                fill = NA,
                color = "grey17",
                linewidth = 0.8
              ) +
              # ggrepel::geom_label_repel(
              #   data              = biome_coords_df,
              #   aes(x = X, y = Y, label = BIOME_NAME),
              #   color             = "black",
              #   na.rm             = TRUE,
              #   fill              = "gray86",
              #   size              = 2.1,
              #   fontface          = "bold",
              #   box.padding       = 0.5,
              #   point.padding     = 0.3,
              #   max.iter          = 50000,
            #   max.overlaps      = Inf,
            #   force             = 100,
            #   force_pull        = 1,
            #   nudge_x           = ifelse(biome_coords_df$X < -75, -25, 20),
            #   nudge_y           = ifelse(biome_coords_df$Y > 10,   15, -15),
            #   arrow             = arrow(length = unit(0.015, "npc"), type = "closed", ends = "last"),
            #   segment.color     = "grey30",
            #   segment.size      = 0.4,
            #   segment.curvature = 0.2,
            #   min.segment.length = 0
            # ) +
            coord_sf(expand = FALSE) +
              # labs(title = paste0(variable, " | ", biomes[[j]])) +
              coord_sf(
                xlim   = c(biome_bbox["xmin"], biome_bbox["xmax"]),
                ylim   = c(biome_bbox["ymin"], biome_bbox["ymax"]),
                expand = FALSE
              ) +
              theme_void() +
              theme(
                legend.position = "right",
                legend.title    = element_text(color = "black", size = 12),
                plot.title      = element_text(
                  color = "black",
                  size = 14,
                  face = "bold",
                  hjust = 0.5
                ),
                plot.background = element_rect(fill = "white", color = NA)
              )
            
            
            
            
            ggsave(
              paste0(
                out_buffer_dir_summary_png,
                "/",
                variable,
                "/",
                biome_name,
                "/",
                quantile_res,
                "/",
                name,
                "_",
                "species_richness_current.png"
              ),
              plot = p,
              width = 14,
              height = 14,
              dpi = 1000
            )
          } #end l
        } # end k
      } # end j
    }   # end i
    
    return("DONE!")
  }

################################################################################
#config adding a results_dir
output_dir <- "//catalogue/MultifLandscapesA1706/1.Data/Results"
#user
user <- "Chrystian"
if (user == "Chrystian") {
  basedir <-  "//catalogue/MultifLandscapesA1706"
}

#defining where the SDMs are located
resultdir <-
  "//catalogue/MultifLandscapesA1706/1.Data/Results/SDM results"

#getting species lists
species_selected <- as.data.frame(
  readxl::read_xlsx(
    "E:/CSOSA/Dropbox/VAVILOV_2.0/Results/species_lists/species_curated/TO_CHECK/Species_20260508_edible_part_curated.xlsx",
    col_names = T
  )
)

#defining template raster of latam
template_file <-
  paste0(
    "//catalogue/MultifLandscapesA1706/1.Data/RAW/Input_data/climate_data/2_5min/present/",
    "wc2.1_2.5m_bio_1.tif"
  )

#reading species per quantile and results
data_ori <-
  read.csv(
    "//catalogue/MultifLandscapesA1706/1.Data/Process/soil_extreme/outliers_list2.csv"
  )

#reading and getting conversion table to use all files togethers
conversion_table <-
  readxl::read_xlsx(
    "E:/CSOSA/Dropbox/VAVILOV_2.0/Results/SUMMARY_FILES/PROJECT_STATUS.xlsx",
    "Hoja1"
  )

x <- summmary_raster_var_toile_ecosys(
  basedir,
  resultdir,
  output_dir,
  species_selected,
  template_file,
  data_ori,
  conversion_table
)