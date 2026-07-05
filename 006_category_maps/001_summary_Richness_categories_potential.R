### this script can be used to get summary richness maps using potential niche tif files for categories for current
# author: Chrystian Sosa

library(terra)
library(viridis)
library(sf)
library(foreach)
library(raster)
library(bit)
# library(ggplot2)
library(bit)
require(leaflet)
# library(htmlwidgets)
library(viridisLite)
################################################################################

# set target region
summary_richness_function_potential <-
  function(basedir,
           column_name,
           resultdir,
           template_file,
           species_selected) {
    ##############################################################################
    #config adding a results_dir
    output_dir <- "//catalogue/MultifLandscapesA1706/1.Data/Results"
    out_dir <- output_dir
    #richness
    out_buffer_dir <- paste0(out_dir, "/Summary_maps_categories")
    if (!dir.exists(out_buffer_dir)) {
      dir.create(out_buffer_dir)
    }
    # #species_raster_dir
    # out_buffer_dir_raster <- paste0(out_buffer_dir, "/raster")
    # if (!dir.exists(out_buffer_dir_raster)) {
    #   dir.create(out_buffer_dir_raster)
    # }
    #summary_dir
    out_buffer_dir_summary <- out_buffer_dir
    #paste0(out_buffer_dir, "/summary")
    if (!dir.exists(out_buffer_dir_summary)) {
      dir.create(out_buffer_dir_summary)
    }
    #folder for PNG, TIF, and HTML files
    out_buffer_dir_summary_png <-
      paste0(out_buffer_dir_summary, "/", "PNG")
    out_buffer_dir_summary_tif <-
      paste0(out_buffer_dir_summary, "/", "TIF")
    out_buffer_dir_summary_html <-
      paste0(out_buffer_dir_summary, "/", "HTML")
    
    
    message(out_buffer_dir_summary)
    
    
    if (!dir.exists(out_buffer_dir_summary_png)) {
      dir.create(out_buffer_dir_summary_png)
    }
    if (!dir.exists(out_buffer_dir_summary_tif)) {
      dir.create(out_buffer_dir_summary_tif)
    }
    if (!dir.exists(out_buffer_dir_summary_html)) {
      dir.create(out_buffer_dir_summary_html)
    }
    
    # chunks summary_dir
    out_buffer_dir_summary_chunks <-
      paste0(out_buffer_dir_summary, "/chunks")
    if (!dir.exists(out_buffer_dir_summary_chunks)) {
      dir.create(out_buffer_dir_summary_chunks)
    }
    ################################################################################
    #get species list
    # species_selected <- read.csv(
    #   paste0(out_buffer_dir_summary, "/", "species_selected.csv")
    # )
    #categories to perform analyses
    message(paste0("Step 4... Getting categories for ", column_name))
    
    categories <- unique(species_selected[, c(column_name)])
    
    #removing cats
    categories <-
      categories[!categories %in% c(NA, "TO FILL", "N/A", "?",
                                    "Unclear plant part/food group reported in edible portion",
                                    "Not available plant part")]
    print(categories)
    
    ################################################################################
    template <- terra::rast(template_file)
    template[which(!is.na(template[]))] <- 1
    ##############################################################################
    #loading adm0
    message("loading shapefile and metrics")
    adm0 <-
      vect(paste0(basedir, "/1.Data/RAW/input_data/adm0/adm0_Latam", ".shp"))
    #all species
    species_total <- unique(species_selected$species)
    ##############################################################################
    #loading buffer species
    InDir_buffer <-
      "//catalogue/MultifLandscapesA1706/1.Data/Results/buffer_richness/raster"
    #species with buffer
    sp_buffer_list <- list.files(InDir_buffer, pattern = ".tif")
    sp_buffer_list_name <-
      sub(pattern = ".tif", replacement = "", sp_buffer_list)
    ##############################################################################
    #load metrics to get modeled species
    metrics <-
      read.csv(paste0(resultdir, "/", "AUC_Maxent_valid.csv"))
    #gettting species lists to get results (modeled species with MaxEnt)
    species_list <- metrics$species
    species_list <-
      sub(pattern = "model_evaluation_", replacement = "", species_list)
    species_list <-
      sub(pattern = ".csv", replacement = "", species_list)
    # print(length(species_list))
    # species_list <- species_list[which(species_list!="Cabralea canjerana")]
    # print(length(species_list))
    #species with valid model
    avail_sp_modeled <-
      species_list[species_list %in% species_total]
    ##############################################################################
    #species with buffer
    # buffer_sps <- species_total[sp_buffer_list_name %in% species_total]
    buffer_sps <-
      species_total[!species_total %in% avail_sp_modeled]
    #if buffer_sps + avail_sp_modeled =8221 is well done!
    buffer_sps <-
      buffer_sps[buffer_sps %in% sp_buffer_list_name] #available to use buffer
    ##############################################################################
    #adding presence absence results folder
    InDir <-
      paste0(resultdir, "/Distribution maps/Presence-absence")
    ################################################################################
    ####First step (summary file)
    print(categories)
    print(tapply(species_selected[, c(column_name)],species_selected[, c(column_name)],length))
    #s i <- 15
    # i <- 4
    for (i in 1:length(categories)) {
      #subsetting category for maps
      category <- categories[[i]]
      if (category == "fruits?") {
        category_name <- "fruits_dubius"
      } else if (category == "pulses?") {
        category_name <- "pulses_dubius"
      } else if (category == "oil/fats") {
        category_name <- "oil_fats"
      } else if (category == "gums/mucilages") {
        category_name <- "gums_mucilages"
      } else if (category == "roots/tubers") {
        category_name <- "roots"
      } else if (category == "other food types") {
        category_name <- "other_food_types"
      } else if (category == "food additive") {
        category_name <- "food_additive"
      } else if (category == "sugar and syrups") {
        category_name <- "sugar_and_syrups"
      } else if (category == "nuts/seeds") {
        category_name <- "nuts_seeds"
      } else if(category=="Non-woody epiphyte"){
        category_name <- "NW_epiphyte"
      } else if(category=="unspecified aerial parts"){
        category_name <- "UAP"
      }  else if(category=="seedlings/germinated seeds"){
        category_name <- "SG"
      }  else if(category=="root/tuber vegetables"){
        category_name <- "roots_tubers"
      } else {
        category_name <- category
      }
      message(paste0("starting ",category))
      ################################################################################
      species_subsetted <-
        species_selected[which(species_selected[, column_name] == category),]
      avail_sp_modeled_i <-
        species_subsetted[species_subsetted$species %in% avail_sp_modeled,]
      buffer_sps_i <-
        species_subsetted[species_subsetted$species %in% buffer_sps,]
      ################################################################################
      if (!file.exists(
        paste0(
          out_buffer_dir_summary_tif,
          "/",
          column_name,
          "_",
          category_name,
          "_",
          "richness_current.tif"
        )
      )) {
        if (!file.exists(
          paste0(
            out_buffer_dir_summary_tif,
            "/",
            column_name,
            "_",
            category_name,
            "_",
            "SDM_richness_current.tif"
          )
        )) {
          #message("Doing summary for current (parallel chunked mode)")
          
          # --- 1. Validate paths sequentially ---------------------------------------
          message("Checking file availability...")
          valid_paths <-
            vapply(avail_sp_modeled_i$species, function(sp) {
              path <- paste0(InDir, "/", sp, ".tif")
              if (file.exists(path))
                path
              else
                NA_character_
            }, character(1))
          
          missing_sps <-
            avail_sp_modeled_i$species[is.na(valid_paths)]
          valid_paths <- valid_paths[!is.na(valid_paths)]
          
          if (length(missing_sps) > 0)
            message(paste0(
              length(missing_sps),
              " species raster(s) not found and skipped."
            ))
          message(paste0(length(valid_paths), " valid rasters found."))
          
          if (length(valid_paths) > 0) {
            # --- 2. Split into chunks of 100 ------------------------------------------
            if(length(valid_paths)<10){
              chunk_size <- length(valid_paths) 
            } else {
              chunk_size <- 10 
              }
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
            
            # --- 3. Setup parallel backend --------------------------------------------
            n_cores <- 8 #min(parallel::detectCores() - 2, n_chunks)
            cl      <- parallel::makeCluster(n_cores)
            doParallel::registerDoParallel(cl)
            message(paste0("Using ", n_cores, " cores for ", n_chunks, " chunks"))
            
            # --- Create chunks subfolder ----------------------------------------------
            chunk_dir <-
              paste0(out_buffer_dir_summary_chunks, "/chunks2")
            if (!dir.exists(chunk_dir))
              dir.create(chunk_dir)
            
            # --- 4. Parallel chunk summation with NA → 0 ------------------------------
            chunk_results <- foreach::foreach(
              chunk          = chunk_idx,
              .packages      = "terra",
              .errorhandling = "pass"
            ) %dopar% {
              paths_chunk <- valid_paths[chunk]
              chunk_id    <-
                chunk[1]  # use first index as unique chunk identifier
              
              stack  <- terra::rast(paths_chunk)
              result <- terra::app(stack, fun = "sum", na.rm = TRUE)
              
              # Save to known folder with predictable name
              out_path <-
                paste0(chunk_dir,
                       "/chunk_",
                       sprintf("%04d", chunk_id),
                       ".tif")
              terra::writeRaster(result,
                                 filename = out_path,
                                 overwrite = TRUE)
              out_path
            }
            
            parallel::stopCluster(cl)
            
            # --- 5. Check for chunk errors --------------------------------------------
            failed <- sapply(chunk_results, inherits, "error")
            if (any(failed)) {
              message(paste0(sum(failed), " chunk(s) failed:"))
              print(chunk_results[failed])
              stop("Aborting: fix failed chunks before proceeding.")
            }
            
            # --- 6. Final sum across chunk files --------------------------------------
            message("Combining chunk results into final richness raster...")
            chunk_paths <- unlist(chunk_results)
            chunk_stack <- terra::rast(chunk_paths)
            x           <-
              terra::app(chunk_stack, fun = "sum", na.rm = TRUE)
            
            # --- 7. Save final raster -------------------------------------------------
            message("Saving species_richness_current.tif ...")
            terra::writeRaster(
              x,
              filename  = paste0(
                out_buffer_dir_summary_tif,
                "/",
                column_name,
                "_",
                category_name,
                "_",
                "SDM_richness_current.tif"
              ),
              overwrite = TRUE
            )
            
            # --- 8. Cleanup chunks folder ---------------------------------------------
            message("Cleaning up chunk files...")
            file.remove(chunk_paths)
            unlink(chunk_dir)
            richness <- x
            richness <- terra::resample(richness, template, "near")
            message("Done!")
          } else {
            richness <- NULL
          }
        } else {
          richness <- terra::rast(
            paste0(
              out_buffer_dir_summary_tif,
              "/",
              column_name,
              "_",
              category_name,
              "_",
              "SDM_richness_current.tif"
            )
          )
          
        }
        
        
        ################################################################################
        ################################################################################
        #creating for buffer
        message("Doing for buffer")
        if (!file.exists(
          paste0(
            out_buffer_dir_summary_tif,
            "/",
            column_name,
            "_",
            category_name,
            "_",
            "buffer_richness_current.tif"
          )
        )) {
          message("Doing summary for current (parallel chunked mode)")
          
          # --- 1. Validate paths sequentially ---------------------------------------
          message("Checking file availability...")
          valid_paths <- vapply(buffer_sps_i$species, function(sp) {
            path <- paste0(InDir_buffer, "/", sp, ".tif")
            if (file.exists(path))
              path
            else
              NA_character_
          }, character(1))
          
          missing_sps <- buffer_sps_i$species[is.na(valid_paths)]
          valid_paths <- valid_paths[!is.na(valid_paths)]
          
          if (length(missing_sps) > 0)
            message(paste0(
              length(missing_sps),
              " species raster(s) not found and skipped."
            ))
          message(paste0(length(valid_paths), " valid rasters found."))
          
          
          if (length(valid_paths > 0)) {
            # --- 2. Split into chunks of 100 ------------------------------------------
            
            if(length(valid_paths)<10){
              chunk_size <- length(valid_paths) 
            } else {
              chunk_size <- 10 
            }
            
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
            
            # --- 3. Setup parallel backend --------------------------------------------
            n_cores <- 8 #min(parallel::detectCores() - 2, n_chunks)
            cl      <- parallel::makeCluster(n_cores)
            doParallel::registerDoParallel(cl)
            message(paste0("Using ", n_cores, " cores for ", n_chunks, " chunks"))
            
            # --- Create chunks subfolder ----------------------------------------------
            chunk_dir <-
              paste0(out_buffer_dir_summary_tif, "/chunks2")
            if (!dir.exists(chunk_dir))
              dir.create(chunk_dir)
            
            # --- 4. Parallel chunk summation with NA → 0 ------------------------------
            chunk_results <- foreach::foreach(
              chunk          = chunk_idx,
              .packages      = "terra",
              .errorhandling = "pass"
            ) %dopar% {
              paths_chunk <- valid_paths[chunk]
              chunk_id    <-
                chunk[1]  # use first index as unique chunk identifier
              
              stack  <- terra::rast(paths_chunk)
              result <- terra::app(stack, fun = "sum", na.rm = TRUE)
              
              # Save to known folder with predictable name
              out_path <-
                paste0(chunk_dir,
                       "/chunk_",
                       sprintf("%04d", chunk_id),
                       ".tif")
              terra::writeRaster(result,
                                 filename = out_path,
                                 overwrite = TRUE)
              out_path
            }
            
            parallel::stopCluster(cl)
            
            # --- 5. Check for chunk errors --------------------------------------------
            failed <- sapply(chunk_results, inherits, "error")
            if (any(failed)) {
              message(paste0(sum(failed), " chunk(s) failed:"))
              print(chunk_results[failed])
              stop("Aborting: fix failed chunks before proceeding.")
            }
            
            # --- 6. Final sum across chunk files --------------------------------------
            message("Combining chunk results into final richness raster...")
            chunk_paths <- unlist(chunk_results)
            chunk_stack <- terra::rast(chunk_paths)
            x           <-
              terra::app(chunk_stack, fun = "sum", na.rm = TRUE)
            
            # --- 7. Save final raster -------------------------------------------------
            message("Saving species_richness_buffer_current.tif ...")
            terra::writeRaster(
              x,
              filename  = paste0(
                out_buffer_dir_summary_tif,
                "/",
                column_name,
                "_",
                category_name,
                "_",
                "buffer_richness_current.tif"
              ),
              overwrite = TRUE
            )
            richness_buffer <- x
            # --- 8. Cleanup chunks folder ---------------------------------------------
            message("Cleaning up chunk files...")
            file.remove(chunk_paths)
            unlink(chunk_dir)
            message("Done!")
            
          } else {
            richness_buffer <- NULL
          }
        } else {
          richness_buffer <-
            terra::rast(
              paste0(
                out_buffer_dir_summary_tif,
                "/",
                column_name,
                "_",
                category_name,
                "_",
                "buffer_richness_current.tif"
              )
            )
          
        }
        
        
        if(!is.null(richness)){
          richness_res <-
            terra::resample(richness,template,method="near")
          
        } else {
          richness_res <- NULL
        }
          
          
        
          if (!is.null(richness_res) & !is.null(richness_buffer)) {
          message("SDM AVAIL, BUFFER AVAIL")
          
          richness_buffer_res <-
            terra::resample(richness_buffer, template, method = "near")
          richness_total <- sum(richness_res,richness_buffer_res,na.rm = T)
          
          terra::writeRaster(
            richness_total,
            filename  = paste0(
              out_buffer_dir_summary_tif,
              "/",
              column_name,
              "_",
              category_name,
              "_",
              "richness_current.tif"
            ),
            overwrite = TRUE
          )
          
          
        } else if (is.null(richness_res) &
                   !is.null(richness_buffer)) {
          
          message("SDM AVAIL, BUFFER AVAIL")
          
          richness_buffer_res <-
            terra::resample(richness_buffer, template, method = "near")
          richness_total <- richness_buffer_res
          
          terra::writeRaster(
            richness_total,
            filename  = paste0(
              out_buffer_dir_summary_tif,
              "/",
              column_name,
              "_",
              category_name,
              "_",
              "richness_current.tif"
            ),
            overwrite = TRUE
          )
          
        } else if (!is.null(richness_res) &
                   is.null(richness_buffer)) {
          message("NO SDM, BUFFER AVAIL")
          richness_total <- richness_res
          
          terra::writeRaster(
            richness_total,
            filename  = paste0(
              out_buffer_dir_summary_tif,
              "/",
              column_name,
              "_",
              category_name,
              "_",
              "richness_current.tif"
            ),
            overwrite = TRUE
          )
          
        } else if (is.null(richness_res) &
                   is.null(richness_buffer)) {
          message("NO FILE TO MAP")
          richness_total <- NULL
          
        }
        
      } else {
        richness_total <-
          terra::rast(
            paste0(
              out_buffer_dir_summary_tif,
              "/",
              column_name,
              "_",
              category_name,
              "_",
              "richness_current.tif"
            )
          )
      }
      
      if (!is.null(richness_total)) {
        richness_total[which(richness_total[] == 0)] <- NA
        ################################################################################
        ################################################################################
        ################################################################################
        ################################################################################
        message("Making a plot")
        
        #Plotting
        par(mar = c(1, 1, 1, 1))
        
        setwd(dir = paste0(out_buffer_dir_summary_png))
        png(
          paste0(
            "LATAM_Edible_food_",
            column_name,
            "_",
            category_name,
            "_species_richness_current",
            ".png"
          ),
          width = 10,
          height = 10,
          units = "in",
          # interpret as inches
          res = 1000
        )      # 300 dpi = print quality
        par(mfrow = c(1, 1))
        par(mai = c(0.5, 0.5, 0.5, 0.5))
        plot(
          richness_total,
          main = "",
          #main = paste0("Species richness"),
          legend = T,
          col = viridis(length(species_list))
        )
        plot(
          sf::st_geometry(sf::st_as_sf(adm0)),
          add = TRUE,
          border = "grey17",
          lwd = 1.25
        )
        dev.off()
        ################################################################################
        ################################################################################
        ################################################################################
        ################################################################################
        message(paste("Step 5.9: HTML file step ", category))
        
        # Get the range of raster values (excluding NA)
        raster_values <- values(richness_total, na.rm = TRUE)
        
        # num_classes <- 8  # Adjust this number as needed
        message(paste("Step 5.9.1: HTML breaks ", category))
        
        # Generate breaks using quantiles (0 to 90th percentile, then add max)
        breaks <-
          quantile(
            raster_values,
            probs = seq(0, 0.9, length.out = num_classes + 1),
            na.rm = TRUE
          )
        
        breaks <- unique(breaks)
        
        # Handle single unique value case
        if (length(breaks) < 2) {
          unique_val <- breaks[1]
          breaks <- c(unique_val - 0.5, unique_val + 0.5)
          colors <- "red"
          actual_num_classes <- 1
          auto_labels <- as.character(unique_val)
        } else {
          actual_num_classes <- length(breaks) - 1
          colors <- turbo(actual_num_classes)
          
          auto_labels <-
            sapply(seq_len(actual_num_classes), function(i) {
              lo <- breaks[i]
              hi <- breaks[i + 1]
              if (i == actual_num_classes) {
                paste0(">", lo)
              } else if (lo == hi - 1 || lo == hi) {
                as.character(lo)
              } else {
                paste0(lo, "-", hi - 1)
              }
            })
        }
        
        
        message(paste("Step 5.9.3: HTML colours ", category))
        
        # Create categorical colour palette
        pal <- colorBin(
          palette   = colors,
          domain    = raster_values,
          bins      = breaks,
          na.color  = "transparent"
        )
        
        # richness_total_rl <- terra::project(richness_total,"EPSG:4326")
        # richness_total_rl <- raster::raster(richness_total_rl)
        # 
        x <- leaflet() %>%
          addTiles() %>%
          addProviderTiles(providers$CartoDB.Positron) %>%
          addMiniMap(width = 150, height = 150) %>%
          addRasterImage(
            richness_total,
            colors   = pal,
            opacity  = 0.3,
            project  = TRUE,
            maxBytes = 20 * 1024 ^ 2
          ) %>%
          addLegend(
            pal      = pal,
            values   = raster_values,
            title    = "Species Richness",
            labFormat = function(type, cuts, p)
              auto_labels  # fully automatic
          )
        message(paste("Step 5.9.4: HTML saving ", category))
        
        htmlwidgets::saveWidget(
          x,
          paste0(
            out_buffer_dir_summary_html,
            "/",
            column_name,
            "_",
            category_name,
            "_SDM_buffer_richness.html"
          ),
          selfcontained = F
        )
        message(paste("DONE: ", category))
      } else {
        message(paste0("skipping... ", category))
      }
    }
    return("DONE!")
  }

user <- "Chrystian"
if (user == "Chrystian") {
  basedir <-  "//catalogue/MultifLandscapesA1706"
}

resultdir <-
  "//catalogue/MultifLandscapesA1706/1.Data/Results/SDM results"


species_selected <- as.data.frame(
  readxl::read_xlsx(
    # "/catalogue/MultifLandscapesA1706/1.Data/Results/species_lists/species_curated/Species_20250304_edible_part_curated.xlsx",
    "E:/CSOSA/Dropbox/VAVILOV_2.0/Results/species_lists/species_curated/TO_CHECK/Species_20260508_edible_part_curated.xlsx",
    col_names = T
  )
)

# column_name <- "plant_status"
# colnames(species_selected)[16] <- "growth_form"
# column_name <- "growth_form"
# column_name <- "food_group_to_map"
column_name <- "plant_part"
#template
#load template for rasterize
template_file <-
  paste0(
    "//catalogue/MultifLandscapesA1706/1.Data/RAW/Input_data/climate_data/2_5min/present/",
    "wc2.1_2.5m_bio_1.tif"
  )


num_classes <- 8
summary_richness_function_potential(basedir,
                          column_name,
                          resultdir,
                          template_file,
                          species_selected)

################################################################################
#template
#load template for rasterize
template_file <-
  paste0(
    "//catalogue/MultifLandscapesA1706/1.Data/RAW/Input_data/climate_data/2_5min/present/",
    "wc2.1_2.5m_bio_1.tif"
  )

column_name <- "food_group"

num_classes <- 8
summary_richness_function_potential(basedir,
                                    column_name,
                                    resultdir,
                                    template_file,
                                    species_selected)

column_name <- "cultivated"
summary_richness_function_potential(basedir,
                                   column_name,
                                   resultdir,
                                   template_file,
                                   species_selected)
column_name <- "wild"
summary_richness_function_potential(basedir,
                                   column_name,
                                   resultdir,
                                   template_file,
                                   species_selected)
