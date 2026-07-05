### this script can be used to get summary richness maps using potential niche tif files in several categories for future projections


# author: Chrystian Sosa

library(terra)
library(viridis)
library(sf)
library(foreach)
# library(raster)
library(bit)
# library(ggplot2)
library(bit)
require(leaflet)
# library(htmlwidgets)
library(viridisLite)
################################################################################
#MaxEnt results
# resultdir <- "/catalogue/MultifLandscapesA1706/1.Data/Results/CWR_results"

# set target region
summary_richness_function_FUTURE <-
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
      # paste0(out_buffer_dir, "/summary_future")
    if (!dir.exists(out_buffer_dir_summary)) {
      dir.create(out_buffer_dir_summary)
    }
    #folder for PNG, TIF, and HTML files
    out_buffer_dir_summary_png <-
      paste0(out_buffer_dir_summary, "/", "PNG_future")
    out_buffer_dir_summary_tif <-
      paste0(out_buffer_dir_summary, "/", "TIF_future")
    out_buffer_dir_summary_html <-
      paste0(out_buffer_dir_summary, "/", "HTML_future")
    
    
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
      paste0(out_buffer_dir_summary, "/chunks_F")
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
      categories[!categories %in% c(NA, "TO FILL", "N/A", "?")]
    print(categories)
    
    ################################################################################
    template <- terra::rast(template_file)
    template[which(!is.na(template[]))] <- 1
    ##############################################################################
    #loading adm0
    message("loading shapefile and metrics")
    #all species
    species_total <- unique(species_selected$species)
    ##############################################################################
    #loading buffer species
    # InDir_buffer <-
    #   "//catalogue/MultifLandscapesA1706/1.Data/Results/buffer_richness/raster"
    # #species with buffer
    # sp_buffer_list <- list.files(InDir_buffer, pattern = ".tif")
    # sp_buffer_list_name <-
    #   sub(pattern = ".tif", replacement = "", sp_buffer_list)
    # ##############################################################################
    #load metrics to get modeled species
    metrics1 <-
      read.csv(paste0(resultdir, "/", "AUC_Maxent_valid_ALL.csv"))
    metrics2 <-
      read.csv(paste0(resultdir, "/", "AUC_Maxent_warning.csv"))
    
    metrics <- rbind(metrics1,metrics2)
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
    # buffer_sps <-
    #   species_total[!species_total %in% avail_sp_modeled]
    # #if buffer_sps + avail_sp_modeled =8221 is well done!
    # buffer_sps <-
    #   buffer_sps[buffer_sps %in% sp_buffer_list_name] #available to use buffer
    ##############################################################################
    #adding presence absence results folder
    InDir <-
      paste0(resultdir, "/Distribution maps/Future/2050/Consensus_maps_folder2/Fut_consensus/NO_HULL")
    ################################################################################
    ####First step (summary file)
    
    SSPs <- c("ssp245","ssp370")
    
    
  for(s in 1:2){
    # s <- 2
    for (i in 1:length(categories)) {
      # i <- 1
      adm0 <-
        terra::vect(paste0(basedir, "/1.Data/RAW/input_data/adm0/adm0_Latam", ".shp"))
      
      # i <- 10
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
      } else if (category == "root/tubers") {
        category_name <- "roots"
      } else if (category == "other food types") {
        category_name <- "other_food_types"
      } else if (category == "food additive") {
        category_name <- "food_additive"
      } else if (category == "sugar and syrups") {
        category_name <- "sugar_and_syrups"
      } else if (category == "nuts and seeds") {
        category_name <- "nuts_and_seeds"
      } else if(category=="Non-woody epiphyte"){
        category_name <- "NW_epiphyte"
      } else if(category=="unspecified aerial parts"){
        category_name <- "UAP"
      }  else if(category=="seedlings/germinated seed"){
        category_name <- "SG"
      } else if(category=="cultivated/wild"){
        category_name <- "Cultivated_wild"
      } else if(category=="cultivated/wild/domesticated"){
        category_name <- "Cultivated_wild_domesticated"
      } else if(category=="Gene pool 1  (Cultivated species)"){
        category_name <- "GP1"
      }else if(category=="roots/tubers"){
        category_name <- "root_tubers"
      } else if (category == "nuts/seeds"){
        category_name <- "nuts_and_seeds"
      } else if (category == "seedlings/germinated seeds"){
      category_name <- "SG"
      } else if (category == "root/tuber vegetables"){
        category_name <- "root_tuber"
    } else {
        category_name <- category
      }
      message(paste0("starting ",category))
      ################################################################################
      species_subsetted <-
        species_selected[which(species_selected[, column_name] == category),]
      avail_sp_modeled_i <-
        species_subsetted[species_subsetted$species %in% avail_sp_modeled,]
      avail_sp_modeled_i <- data.frame(species=unique(avail_sp_modeled_i$species))
      # buffer_sps_i <-
      #   species_subsetted[species_subsetted$species %in% buffer_sps,]
      # buffer_sps_i <- data.frame(species=unique(buffer_sps_i$species))
      # ################################################################################

        if (!file.exists(
          paste0(
            out_buffer_dir_summary_tif,
            "/",
            column_name,
            "_",
            category_name,
            "_",
            "SDM_richness_future_",SSPs[[s]],".tif"
          )
        )) {
          #message("Doing summary for current (parallel chunked mode)")
          
          # --- 1. Validate paths sequentially ---------------------------------------
          message("Checking file availability...")
          valid_paths <-
            vapply(avail_sp_modeled_i$species, function(sp) {
              path <- paste0(InDir, "/", sp,"_",SSPs[[s]], ".tif")
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
              paste0(out_buffer_dir_summary_chunks, "/chunks3")
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
              stack  <- terra::app(stack,function(x) ifelse(x>=3,1,0))
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
                "SDM_richness_future_",SSPs[[s]],".tif"
              ),
              overwrite = TRUE
            )
            
            # --- 8. Cleanup chunks folder ---------------------------------------------
            message("Cleaning up chunk files...")
            file.remove(chunk_paths)
            unlink(chunk_dir)
            richness <- x
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
              "SDM_richness_future_",SSPs[[s]],".tif"
            )
          )
          
        }
        
      
      richness_total <- richness
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
            "_species_richness_future_",SSPs[[s]],
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
        
        rm(adm0)
        ################################################################################
        ################################################################################
        ################################################################################
        ################################################################################
        message(paste("Step 5.9: HTML file step ", category))
        
        if(file.exists(paste0(
          out_buffer_dir_summary_tif,
          "/",
          column_name,
          "_",
          category_name,
          "_",
          "SDM_richness_future_",SSPs[[s]],".tif"
        ))){
          richness_total_rl <- 
            terra::rast(
              # raster::raster(            
              paste0(
                out_buffer_dir_summary_tif,
                "/",
                column_name,
                "_",
                category_name,
                "_",
                "SDM_richness_future_",SSPs[[s]],".tif"
              ))
          richness_total_rl[which(richness_total_rl[] == 0)] <- NA
          
          raster_values <- raster::values(richness_total_rl)
          raster_values <- raster_values[!is.na(raster_values)]
          message(paste("Step 5.9.1: HTML breaks ", category))
          
          breaks <- quantile(
            raster_values,
            probs = seq(0, 0.9, length.out = num_classes + 1),
            na.rm = TRUE
          )
          breaks <- unique(breaks)
          
          if (length(breaks) < 2) {
            unique_val <- breaks[1]
            breaks     <- c(unique_val - 0.5, unique_val + 0.5)
            colors     <- "red"
            actual_num_classes <- 1
            auto_labels <- as.character(unique_val)
          } else {
            actual_num_classes <- length(breaks) - 1
            colors <- turbo(actual_num_classes)
            auto_labels <- sapply(seq_len(actual_num_classes), function(i) {
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
          
          pal <- colorBin(
            palette  = colors,
            domain   = range(raster_values),
            bins     = breaks,
            na.color = "transparent"
          )
          
          x <- leaflet() %>%
            addTiles() %>%
            addProviderTiles(providers$CartoDB.Positron) %>%
            addMiniMap(width = 150, height = 150) %>%
            addRasterImage(
              richness_total_rl,
              colors   = pal,
              opacity  = 0.3,
              project  = T,
              maxBytes = 20 * 1024 ^ 2
            ) %>%
            addLegend(
              pal      = pal,
              values   = raster_values,
              title    = "Species Richness",
              labFormat = function(type, cuts, p) auto_labels
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
              
              "_SDM_richness_future_",SSPs[[s]],".html"
            ),
            selfcontained = F
          )
          message(paste("DONE: ", category))
        } else {
          message(paste0("skipping... ", category))
        }
      }
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
# species_selected$all <- 1
# colnames(species_selected)[16] <- "growth_form"
# column_name <- "growth_form"
# column_name <- "food_group_to_map"
# column_name <- "plant_part"
#template
#load template for rasterize
template_file <-
  paste0(
    "//catalogue/MultifLandscapesA1706/1.Data/RAW/Input_data/climate_data/2_5min/present/",
    "wc2.1_2.5m_bio_1.tif"
  )


num_classes <- 8

colnames(species_selected)
column_name <- "wild"
summary_richness_function_FUTURE(basedir,
                          column_name,
                          resultdir,
                          template_file,
                          species_selected)
column_name <- "cultivated"
summary_richness_function_FUTURE(basedir,
                          column_name,
                          resultdir,
                          template_file,
                          species_selected)


column_name <- "CWR"
summary_richness_function_FUTURE(basedir,
                          column_name,
                          resultdir,
                          template_file,
                          species_selected)


column_name <- "domesticated"
summary_richness_function_FUTURE(basedir,
                          column_name,
                          resultdir,
                          template_file,
                          species_selected)


column_name <-  "wild_or_cultivated"
summary_richness_function_FUTURE(basedir,
                          column_name,
                          resultdir,
                          template_file,
                          species_selected)

column_name <- "plant_part"
# unique(species_selected$plant_part)
species_selected$plant_part[which(species_selected$plant_part=="Unclear plant part/food group reported in edible portion")] <- NA
species_selected$plant_part[which(species_selected$plant_part=="Not available plant part")] <- NA
species_selected$plant_part[which(species_selected$plant_part=="Food group reported in edible portion")] <- NA
species_selected$plant_part[which(species_selected$plant_part=="entire plant")] <- "entire_plant"
species_selected$plant_part[which(species_selected$plant_part=="?")] <- NA

summary_richness_function_FUTURE(basedir,
                          column_name,
                          resultdir,
                          template_file,
                          species_selected)


column_name <- "food_group"
unique(species_selected$food_group)
species_selected$food_group[which(species_selected$food_group=="No food group available")] <- NA
species_selected$food_group[which(species_selected$food_group=="unknown")] <- NA
summary_richness_function_FUTURE(basedir,
                          column_name,
                          resultdir,
                          template_file,
                          species_selected)


column_name <- "cultivated"
summary_richness_function_FUTURE(basedir,
                                   column_name,
                                   resultdir,
                                   template_file,
                                   species_selected)
column_name <- "wild"
summary_richness_function_FUTURE(basedir,
                                   column_name,
                                   resultdir,
                                   template_file,
                                   species_selected)
