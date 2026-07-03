library(readxl); require(ggpubr); library(dplyr);library(ggrepel);library(ggsci)
################################################################################
dir <- "E:/CSOSA/Dropbox/VAVILOV_2.0/Results/species_lists/species_curated"
outdir <- "E:/CSOSA/Dropbox/VAVILOV_2.0/Results/SUMMARY_FILES"
data <- readxl::read_xlsx(paste0(dir, "/TO_CHECK/", "Species_20260508_edible_part_curated.xlsx"))
################################################################################
#et unique species
summary_df_part <- data %>%
  distinct(plant_part, family, species) %>%  # unique species per plant part and family
  arrange(plant_part, family, species)
################################################################################
#get unique species
summary_df_part <- summary_df_part[!summary_df_part$plant_part %in% c("TO FILL",
                                                                      "?",
                                                                      NA, 
                                                                      "Not available plant part",
                                                                      "Unclear plant part/food group reported in edible portion"), ]
unique(summary_df_part$plant_part )
#get info only per species and plant part
summary_df_part_n <- data %>%
  distinct(plant_part, species) %>%
  filter(!plant_part %in% c("TO FILL",
                            "?",
                            NA, 
                            "Not available plant part",
                            "Unclear plant part/food group reported in edible portion")) %>%
  group_by(plant_part) %>%
  summarise(n_species = n(), 
            species_list = paste(species, collapse = ", "),
            .groups = "drop") %>%
  arrange(plant_part)

summary_df_part_n <- summary_df_part_n[order(summary_df_part_n$n_species,decreasing = T),]

xbar <- ggbarplot(summary_df_part_n, "plant_part", "n_species",
                  fill = "plant_part", color = "plant_part", palette = "Paired",
                  label = TRUE, 
                  # lab.pos = "in", lab.vjust = 1.7,
                  position = position_dodge(0.5),
                  legend = "none",
                  xlab = "Plant part", ylab = "Species count") +
  rotate_x_text(angle = 45, hjust = 1) +
  labs(fill = "Plant part", color = "Plant part")
# xbar
ggsave(paste0(outdir, "/summary_plant_part.png"), 
       plot = xbar, 
       width = 12, height = 7, 
       dpi = 1200, bg = "white")
################################################################################
################################################################################
################################################################################
#get info only per species and food group
summary_df_food_n <- data %>%
  distinct(food_group, species) %>%
  filter(!food_group %in% c("TO FILL", "?", NA,
                            "No food group available",
                            "food_group"
                            )) %>%
  group_by(food_group) %>%
  summarise(n_species = n(), 
            species_list = paste(species, collapse = ", "),
            .groups = "drop") %>%
  arrange(food_group)

summary_df_food_n <- summary_df_food_n[!summary_df_food_n$food_group %in% c("N/A"),]
summary_df_food_n <- summary_df_food_n[order(summary_df_food_n$n_species,decreasing = T),]

xbar <- ggbarplot(summary_df_food_n, "food_group", "n_species",
                  fill = "food_group", color = "food_group", palette = "Paired",
                  label = TRUE, 
                  # lab.pos = "in", lab.vjust = 1.7,
                  position = position_dodge(0.5),
                  legend = "none",
                  xlab = "Food group", ylab = "Species count") +
  rotate_x_text(angle = 45, hjust = 1) +
  labs(fill = "Food group", color = "food_group")+
  scale_fill_manual(values=pal_d3("category20")(20))+
  scale_color_manual(values=pal_d3("category20")(20))
# xbar
ggsave(paste0(outdir, "/summary_food_group.png"), 
       plot = xbar, 
       width = 12, height = 7, 
       dpi = 1200, bg = "white")
################################################################################
################################################################################
################################################################################

#do both (part and groups)
summary_df_part <- data %>%
  distinct(plant_part, food_group, species) %>%
  filter(!plant_part %in% c("TO FILL", "?", NA, "N/A","Not available plant part",
                            "Unclear plant part/food group reported in edible portion")) %>%
  filter(!food_group %in% c("TO FILL", "?", NA,
                            "No food group available",
                            "food_group"))  %>%
  group_by(plant_part, food_group) %>%
  summarise(n_species = n(), 
            species_list = paste(species, collapse = ", "),
            .groups = "drop") %>%
  arrange(plant_part, food_group )
summary_df_part <- summary_df_part[!summary_df_part$food_group %in% c("N/A"),]



xbar <- ggbarplot(summary_df_part, "food_group", "n_species",
                  fill = "plant_part", color = "plant_part", palette = "Paired",
                  label = FALSE,
                  xlab = "Food group", ylab = "Species count") +
  geom_text_repel(aes(label = ifelse(n_species > 10, n_species, ""),
                      group = plant_part),
                  position = position_stack(vjust = 0.5),
                  size = 3, color = "black",
                  box.padding = 0.2,
                  point.padding = 0.1,
                  segment.color = "grey50",
                  segment.size = 0.3,
                  min.segment.length = 0.2,
                  max.overlaps = Inf) +
  rotate_x_text(angle = 45, hjust = 1)

ggsave(paste0(outdir, "/summary_plant_food.png"), 
       plot = xbar, 
       width = 12, height = 7, 
       dpi = 1200, bg = "white")

################################################################################
################################################################################
################################################################################
library(RColorBrewer)
library(randomcoloR)

# 
# # ggbarplot(summary_df_part, "plant_part","species",
#           # fill = "plant_part", color = "plant_part", palette = "Paired",
#           # label = TRUE,, lab.pos = "in")
summary_df_food_Fam <- data %>%
  distinct(food_group, family, species) %>%
  filter(!food_group %in% c("TO FILL", "?", NA, "N/A",
                            "No food group available",
                            "food_group"
                            )) %>%
  group_by(food_group, family) %>%
  summarise(n_species = n(),
            species_list = paste(species, collapse = ", "),
            .groups = "drop") %>%
  arrange(food_group, family)



top5_families <- summary_df_food_Fam %>%
  group_by(food_group) %>%
  slice_max(order_by = n_species, n = 5) %>%
  arrange(food_group, desc(n_species))


# Generate enough distinct colors for all families
n_families <- length(unique(top5_families$family))
set.seed(42)  # for reproducibility
my_palette <- distinctColorPalette(n_families)

xbar <- ggbarplot(top5_families, "food_group", "n_species",
                  fill = "family",
                  label = TRUE,
                  lab.pos = "out", 
                  # lab.vjust = 1.1,
                  position = position_dodge(0.7),
                  lab.size = 3,           # label size inside bars
                  
                  xlab = "Food group",
                  ylab = "Species count") +
  scale_fill_manual(values = my_palette) +
  rotate_x_text(angle = 45, hjust = 1) +
  labs(fill = "Family")

ggsave(paste0(outdir, "/summary_top5fam_food.png"), 
       plot = xbar, 
       width = 12, height = 7, 
       dpi = 1200, bg = "white")

################################################################################
################################################################################
################################################################################

# get info only per species and food group
summary_df_growth_n <- data %>%
  distinct(growth_form_Engemann, species) %>%
  filter(!growth_form_Engemann %in% c("TO FILL", "?", NA,"growth_form_Engemann")) %>%
  group_by(growth_form_Engemann) %>%
  summarise(n_species = n(),
            species_list = paste(species, collapse = ", "),
            .groups = "drop") %>%
  arrange(growth_form_Engemann)

summary_df_growth_n <- summary_df_growth_n[order(summary_df_growth_n$n_species,decreasing = T),]
colnames(summary_df_growth_n)[1] <- "growth_Form"
xbar <- ggbarplot(summary_df_growth_n,"growth_Form", "n_species",
                  fill = "growth_Form", color = "growth_Form", palette = "Paired",
                  label = TRUE,
                  # lab.pos = "in", lab.vjust = 1.7,
                  position = position_dodge(0.5),
                  legend = "none",
                  xlab = "Plant growth form (Engemann)", ylab = "Species count") +
  rotate_x_text(angle = 45, hjust = 1) +
  labs(fill = "growth_Form", color = "growth_Form")
# xbar
ggsave(paste0(outdir, "/summary_growth_form_Engemann.png"),
       plot = xbar,
       width = 12, height = 7,
       dpi = 1200, bg = "white")

################################################################################
################################################################################
################################################################################

# get info only per species and plant cultivated status
summary_df_growth_n <- data %>%
  distinct(plant_cultivated_status, species) %>%
  filter(!plant_cultivated_status %in% c("TO FILL", "?", NA,"growth_form_Engemann",
                                         "plant_cultivated_status")) %>%
  group_by(plant_cultivated_status) %>%
  summarise(n_species = n(),
            species_list = paste(species, collapse = ", "),
            .groups = "drop") %>%
  arrange(plant_cultivated_status)

summary_df_growth_n <- summary_df_growth_n[order(summary_df_growth_n$n_species,decreasing = T),]
colnames(summary_df_growth_n)[1] <- "plant_cultivated_status"
xbar <- ggbarplot(summary_df_growth_n,"plant_cultivated_status", "n_species",
                  fill = "plant_cultivated_status", color = "plant_cultivated_status", palette = "Paired",
                  label = TRUE,
                  # lab.pos = "in", lab.vjust = 1.7,
                  position = position_dodge(0.5),
                  legend = "none",
                  xlab = "plant cultivated status", ylab = "Species count") +
  rotate_x_text(angle = 45, hjust = 1) +
  labs(fill = "plant_cultivated_status", color = "growth_Form")
# xbar
ggsave(paste0(outdir, "/summary_growth_plant_cultivated_status.png"),
       plot = xbar,
       width = 12, height = 7,
       dpi = 1200, bg = "white")
# # 
# # #get info only per species and food group
# # summary_df_growth_n <- data %>%
# #   distinct(`growth_form (Engemann 2016)`, species) %>%
# #   filter(!`growth_form (Engemann 2016)` %in% c("TO FILL", "?", NA)) %>%
# #   group_by(`growth_form (Engemann 2016)`) %>%
# #   summarise(n_species = n(), 
# #             species_list = paste(species, collapse = ", "),
# #             .groups = "drop") %>%
# #   arrange(`growth_form (Engemann 2016)`)
# # 
# # summary_df_growth_n <- summary_df_growth_n[order(summary_df_growth_n$n_species,decreasing = T),]
# # colnames(summary_df_growth_n)[1] <- "growth_Form"
# # xbar <- ggbarplot(summary_df_growth_n,"growth_Form", "n_species",
# #                   fill = "growth_Form", color = "growth_Form", palette = "Paired",
# #                   label = TRUE, 
# #                   # lab.pos = "in", lab.vjust = 1.7,
# #                   position = position_dodge(0.5),
# #                   legend = "none",
# #                   xlab = "Plant growth form", ylab = "Species count") +
# #   rotate_x_text(angle = 45, hjust = 1) +
# #   labs(fill = "growth_Form", color = "growth_Form")
# # # xbar
# # ggsave(paste0(dir, "/summary_growth_form.png"), 
# #        plot = xbar, 
# #        width = 12, height = 7, 
# #        dpi = 1200, bg = "white")

################################################################################
################################################################################
################################################################################