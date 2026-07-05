basedir <- "//catalogue/MultifLandscapesA1706/1.Data/Results"

#load total species names
species_selected <- as.data.frame(
  readxl::read_xlsx(
    # "/catalogue/MultifLandscapesA1706/1.Data/Results/species_lists/species_curated/Species_20250304_edible_part_curated.xlsx",
    # "E:/CSOSA/Dropbox/VAVILOV_2.0/Results/species_lists/species_curated/TO_CHECK/Species_20260508_edible_part_curated.xlsx",
    "E:/CSOSA/Dropbox/VAVILOV_2.0/Results/species_lists/species_curated/TO_CHECK/Species_20260508_edible_part_curated.xlsx",
    col_names = T
  )
)

#
species <- unique(species_selected$species)
#get valid models (AUC >= 0.7)
valid_models <-
  read.csv(paste0(basedir, "/SDM results/AUC_Maxent_valid_all.csv"))
valid_models$species <- sub(pattern = "model_evaluation_",
                            replacement = "",
                            x = valid_models$species)
valid_models$species <- sub(pattern = ".csv",
                            replacement = "",
                            x = valid_models$species)

################################################################################
#loading conversion table
conversion_table <-
  readxl::read_xlsx(
    "E:/CSOSA/Dropbox/VAVILOV_2.0/Results/SUMMARY_FILES/PROJECT_STATUS.xlsx",
    "Hoja1"
  )


#structure
#main dir
maps_dir <- paste0(basedir, "/", "D4N_data")
if (!dir.exists(maps_dir)) {
  dir.create(maps_dir)
}

#future
maps_Fut <- paste0(maps_dir, "/", "Future")
if (!dir.exists(maps_Fut)) {
  dir.create(maps_Fut)
}#
#SSP2
maps_Fut_SSP2 <- paste0(maps_Fut, "/", "SSP2")
if (!dir.exists(maps_Fut_SSP2)) {
  dir.create(maps_Fut_SSP2)
}
#SSP3
maps_Fut_SSP3 <- paste0(maps_Fut, "/", "SSP3")
if (!dir.exists(maps_Fut_SSP3)) {
  dir.create(maps_Fut_SSP3)
}

#future CH
maps_Fut_CH <- paste0(maps_dir, "/", "Future masked by hull")
if (!dir.exists(maps_Fut_CH)) {
  dir.create(maps_Fut_CH)
}
#SSP2 CH
maps_Fut_SSP2_CH <- paste0(maps_Fut_CH, "/", "SSP2")
if (!dir.exists(maps_Fut_SSP2_CH)) {
  dir.create(maps_Fut_SSP2_CH)
}
#SSP3 CH
maps_Fut_SSP3_CH <- paste0(maps_Fut_CH, "/", "SSP3")
if (!dir.exists(maps_Fut_SSP3_CH)) {
  dir.create(maps_Fut_SSP3_CH)
}

#present
maps_Pres <- paste0(maps_dir, "/", "Presence-absence")
if (!dir.exists(maps_Pres)) {
  dir.create(maps_Pres)
}
#present CH
maps_Pres_CH <-
  paste0(maps_dir, "/", "Presence-absence masked by hull")
if (!dir.exists(maps_Pres_CH)) {
  dir.create(maps_Pres_CH)
}
################################################################################
# folders to get files
#present
SDM_current_dir <-
  paste0(basedir, "/", "SDM results/Distribution maps/Presence-absence")
SDM_current_dir_CH <-
  paste0(basedir, "/", "SDM_buffer_results_all/TIF_CH")
#future
SDM_Fut_dir <-
  paste0(
    basedir,
    "/",
    "SDM results/Distribution maps/Future/2050/Consensus_maps_folder2/Fut_consensus/NO_HULL"
  )
SDM_Fut_dir_CH <-
  paste0(
    basedir,
    "/",
    "SDM results/Distribution maps/Future/2050/Consensus_maps_folder2/Fut_consensus/NO_HULL"
  )
#buffer
buffer_dir <- paste0(basedir, "/", "buffer_richness/raster")
#concave hull files
concdir <-
  "//catalogue/MultifLandscapesA1706/1.Data/Results/concave_hull_rasters_ADM0_2_5m"
# ################################################################################
# ################################################################################
# ################################################################################
for (i in 1:length(species)) {
  # i <- 4
  spp <- species[[i]]
  message(i)
  message(spp)
  ##############################################################################
  spp_syn <- NA_character_
  x_syn <- data.frame()
  ##############################################################################
  #This is to know if the speccies has valid model
  x <- valid_models[valid_models$species %in% spp, ]
  ##############################################################################
  #getting possible synonyms
  if (nrow(x) == 0) {
    message("searching for synonym")
    #obtaining possible synonym
    conv_x <-
      conversion_table[which(conversion_table$species == spp), ]
    #getting synonym species name
    spp_syn <- conv_x$species_from_source
    #observing if the synonym is a valid model!
    x_syn <- valid_models[valid_models$species %in% spp_syn, ]
    #loading concave hull
  } else {
    x_syn <- data.frame()
  }
  ##############################################################################
  #loading concave hull
  if (file.exists(paste0(concdir, "/", spp, "_conc.tif"))) {
    message("available file for accepted")
    conc_hull_sp <-
      terra::rast(paste0(concdir, "/", spp, "_conc.tif"))
  } else if (file.exists(paste0(concdir, "/", spp_syn, "_conc.tif"))) {
    message("available file for synonym")
    conc_hull_sp <-
      terra::rast(paste0(concdir, "/", spp_syn, "_conc.tif"))
  } else {
    conc_hull_sp <- NULL
    message("no concave hull!")
  }
  ##############################################################################
  #copy species for present
  if (nrow(x) > 0 || nrow(x_syn) > 0) {
    #copy accepted name species
    if (file.exists(paste0(SDM_current_dir, "/", spp, ".tif"))) {
      file.copy(
        from = paste0(SDM_current_dir, "/", spp, ".tif"),
        to =
          paste0(maps_Pres, "/", spp, ".tif")
      )
    } else if (file.exists(paste0(SDM_current_dir, "/", spp_syn, ".tif"))) {
      #copy synonym as result for species
      file.copy(
        from = paste0(SDM_current_dir, "/", spp_syn, ".tif"),
        to =
          paste0(maps_Pres, "/", spp, ".tif")
      )
    }
    
    ############################################################################
    #copy realized  for present
    if (!is.null(conc_hull_sp)) {
      if (file.exists(paste0(SDM_current_dir, "/", spp, ".tif"))) {
        SDM_sp <- terra::rast(paste0(SDM_current_dir, "/", spp, ".tif"))
        SDM_sp <- SDM_sp * conc_hull_sp
        
        terra::writeRaster(SDM_sp,
                           paste0(maps_Pres_CH, "/", spp, ".tif"),
                           overwrite = T)
      } else if (file.exists(paste0(SDM_current_dir, "/", spp_syn, ".tif"))) {
        #copy synonym as result for species
        SDM_sp <-
          terra::rast(paste0(SDM_current_dir, "/", spp_syn, ".tif"))
        SDM_sp <- SDM_sp * conc_hull_sp
        terra::writeRaster(SDM_sp,
                           paste0(maps_Pres_CH, "/", spp, ".tif"),
                           overwrite = T)
      }
    }
    ############################################################################
    #copy for future SSP2 potential
    if (file.exists(paste0(SDM_Fut_dir, "/", spp, "_ssp245.tif"))) {
      file.copy(
        from = paste0(SDM_Fut_dir, "/", spp, "_ssp245.tif"),
        to =
          paste0(maps_Fut_SSP2, "/", spp, ".tif")
      )
    } else if (file.exists(paste0(SDM_Fut_dir, "/", spp_syn, "_ssp245.tif"))) {
      file.copy(
        from = paste0(SDM_Fut_dir, "/", spp_syn, "_ssp245.tif"),
        to =
          paste0(maps_Fut_SSP2, "/", spp, ".tif")
      )
    }
    ############################################################################
    #copy for future SSP3 potential
    if (file.exists(paste0(SDM_Fut_dir, "/", spp, "_ssp370.tif"))) {
      file.copy(
        from = paste0(SDM_Fut_dir, "/", spp, "_ssp370.tif"),
        to =
          paste0(maps_Fut_SSP3, "/", spp, ".tif")
      )
    } else if (file.exists(paste0(SDM_Fut_dir, "/", spp_syn, "_ssp370.tif"))) {
      file.copy(
        from = paste0(SDM_Fut_dir, "/", spp_syn, "_ssp370.tif"),
        to =
          paste0(maps_Fut_SSP3, "/", spp, ".tif")
      )
    }
    ############################################################################
    #copy for future SSP2 realized
    if (!is.null(conc_hull_sp)) {
      if (file.exists(paste0(SDM_Fut_dir, "/", spp, "_ssp245.tif"))) {
        SDM_sp <- terra::rast(paste0(SDM_Fut_dir, "/", spp, "_ssp245.tif"))
        SDM_sp <- terra::resample(SDM_sp, conc_hull_sp, "near")
        SDM_sp <- SDM_sp * conc_hull_sp
        terra::writeRaster(SDM_sp,
                           paste0(maps_Fut_SSP2_CH, "/", spp, ".tif"),
                           overwrite = T)
        
      } else if (file.exists(paste0(SDM_Fut_dir, "/", spp_syn, "_ssp245.tif"))) {
        SDM_sp <-
          terra::rast(paste0(SDM_Fut_dir, "/", spp_syn, "_ssp245.tif"))
        SDM_sp <- terra::resample(SDM_sp, conc_hull_sp, "near")
        SDM_sp <- SDM_sp * conc_hull_sp
        terra::writeRaster(SDM_sp,
                           paste0(maps_Fut_SSP2_CH, "/", spp, ".tif"),
                           overwrite = T)
        
      }
    }
    ############################################################################
    #copy for future SSP3 realized
    if (!is.null(conc_hull_sp)) {
      if (file.exists(paste0(SDM_Fut_dir, "/", spp, "_ssp370.tif"))) {
        SDM_sp <- terra::rast(paste0(SDM_Fut_dir, "/", spp, "_ssp370.tif"))
        SDM_sp <- terra::resample(SDM_sp, conc_hull_sp, "near")
        SDM_sp <- SDM_sp * conc_hull_sp
        terra::writeRaster(SDM_sp,
                           paste0(maps_Fut_SSP3_CH, "/", spp, ".tif"),
                           overwrite = T)
        
      } else if (file.exists(paste0(SDM_Fut_dir, "/", spp_syn, "_ssp370.tif"))) {
        SDM_sp <-
          terra::rast(paste0(SDM_Fut_dir, "/", spp_syn, "_ssp370.tif"))
        SDM_sp <- terra::resample(SDM_sp, conc_hull_sp, "near")
        SDM_sp <- SDM_sp * conc_hull_sp
        terra::writeRaster(SDM_sp,
                           paste0(maps_Fut_SSP3_CH, "/", spp, ".tif"),
                           overwrite = T)
      }
    }
  } else {
    ############################################################################
    #copy species for present (potential)
    if (file.exists(paste0(buffer_dir, "/", spp, ".tif"))) {
      file.copy(
        from = paste0(buffer_dir, "/", spp, ".tif"),
        to =
          paste0(maps_Pres, "/", spp, ".tif")
      )
    } else if (file.exists(paste0(buffer_dir, "/", spp_syn, ".tif"))) {
      file.copy(
        from = paste0(buffer_dir, "/", spp_syn, ".tif"),
        to =
          paste0(maps_Pres, "/", spp, ".tif")
      )
    }  else {
      message("no buffer file for potential")
    }
    
    ############################################################################
    #copy species for present (realized)
    if (!is.null(conc_hull_sp)) {
      if (file.exists(paste0(buffer_dir, "/", spp, ".tif"))) {
        SDM_sp <- terra::rast(paste0(buffer_dir, "/", spp, ".tif"))
        SDM_sp <- terra::resample(SDM_sp, conc_hull_sp, "near")
        SDM_sp <- SDM_sp * conc_hull_sp
        terra::writeRaster(SDM_sp,
                           paste0(maps_Pres_CH, "/", spp, ".tif"),
                           overwrite = T)
      } else if (file.exists(paste0(buffer_dir, "/", spp_syn, ".tif"))) {
        SDM_sp <- terra::rast(paste0(buffer_dir, "/", spp_syn, ".tif"))
        SDM_sp <- terra::resample(SDM_sp, conc_hull_sp, "near")
        SDM_sp <- SDM_sp * conc_hull_sp
        terra::writeRaster(SDM_sp,
                           paste0(maps_Pres_CH, "/", spp, ".tif"),
                           overwrite = T)
      }
    } else {
      message("no buffer file for realized")
    }
  }
}
