################################################################################
# SCRIPT: Ecological Quantile Spectrum Analysis (Esquema Conceptual Completo)
# PROYECTO: VAVILOV 2.0
################################################################################

library(ggplot2)
library(dplyr)
library(patchwork)

# Semilla global para reproducibilidad
set.seed(42)

# Definir límites globales del eje X para asegurar alineación perfecta
x_limits <- c(8, 52)

# ==============================================================================
# 1. PANEL SUPERIOR (p1): Ecosystem Q90 > Species Q10
# ==============================================================================

ecosystem_data_A <- data.frame(
  Value = rnorm(1000, mean = 20, sd = 4),
  Group = "Ecosystem (n species)"
)

species_data_A <- data.frame(
  Value = rnorm(300, mean = 38, sd = 4),
  Group = "Species"
)

df1 <- rbind(ecosystem_data_A, species_data_A)

q90_ecosystem <- quantile(ecosystem_data_A$Value, 0.90)
q10_focal     <- quantile(species_data_A$Value, 0.10)

p1 <- ggplot(df1, aes(x = Value, fill = Group)) +
  geom_density(alpha = 0.5, color = "black", size = 0.6) +
  
  # Líneas de los cuantiles
  geom_vline(xintercept = q90_ecosystem, color = "#2c3e50", linetype = "dashed", size = 1) +
  geom_vline(xintercept = q10_focal, color = "#e74c3c", linetype = "dashed", size = 1) +
  
  # ANOTACIONES SEGUROS EN EL TECHO (y = Inf)
  annotate("text", x = q90_ecosystem, y = Inf, label = "Ecosystem Q90 ", 
           color = "#2c3e50", angle = 90, hjust = 1, vjust = 1.5, fontface = "bold") +
  annotate("text", x = q10_focal, y = Inf, label = " Species Q10", 
           color = "#e74c3c", angle = 90, hjust = 1, vjust = -0.5, fontface = "bold") +
  
  scale_x_continuous(breaks = seq(10, 50, by = 5)) +
  scale_fill_manual(values = c("Ecosystem (n species)" = "#95a5a6", "Species" = "#e67e22")) +
  labs(
    title = "Condition: Species Q10 > Ecosystem Q90",
    x = "Trait / Environmental Variable Value",
    y = "Density",
    fill = "Dataset"
  ) +
  # SOLUCIÓN: Fijamos explícitamente límites tanto en Y como en X
  coord_cartesian(ylim = c(0, 0.11), xlim = x_limits) + 
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


# ==============================================================================
# 2. PANEL INFERIOR (p2): Ecosystem Q10 > Species Q90 (CORREGIDO Y SEGURO)
# ==============================================================================

ecosystem_data_B <- data.frame(
  Value = rnorm(1000, mean = 38, sd = 4),
  Group = "Ecosystem (n species)"
)

focal_species_data_B <- data.frame(
  Value = rnorm(300, mean = 20, sd = 4),
  Group = "Species"
)

df2 <- rbind(ecosystem_data_B, focal_species_data_B)

q10_ecosystem <- quantile(ecosystem_data_B$Value, 0.10)
q90_species   <- quantile(focal_species_data_B$Value, 0.90)

p2 <- ggplot(df2, aes(x = Value, fill = Group)) +
  geom_density(alpha = 0.5, color = "black", size = 0.6) +
  
  # Líneas de los cuantiles
  geom_vline(xintercept = q10_ecosystem, color = "#2c3e50", linetype = "dashed", size = 1) +
  geom_vline(xintercept = q90_species, color = "#e74c3c", linetype = "dashed", size = 1) +
  
  # ANOTACIONES SEGUROS EN EL TECHO (y = Inf)
  annotate("text", x = q10_ecosystem, y = Inf, label = " Ecosystem Q10", 
           color = "#2c3e50", angle = 90, hjust = 1, vjust = -0.5, fontface = "bold") +
  annotate("text", x = q90_species, y = Inf, label = "Species Q90 ", 
           color = "#e74c3c", angle = 90, hjust = 1, vjust = 1.5, fontface = "bold") +
  
  scale_x_continuous(breaks = seq(10, 50, by = 5)) +
  scale_fill_manual(values = c("Ecosystem (n species)" = "#95a5a6", "Species" = "#e67e22")) +
  labs(
    title = "Condition: Ecosystem Q10 > Species Q90",
    x = "Trait / Environmental Variable Value",
    y = "Density",
    fill = "Dataset"
  ) +
  # SOLUCIÓN: Mismos límites exactos que el panel superior
  coord_cartesian(ylim = c(0, 0.11), xlim = x_limits) + 
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


# ==============================================================================
# 3. COMPOSICIÓN FINAL CON PATCHWORK Y GUARDADO
# ==============================================================================

combined_panel <- (p1 / p2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "Ecological Quantile Spectrum Analysis",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
  ) & 
  theme(
    legend.position = "top", 
    legend.title = element_text(face = "bold")
  )

# Mostrar en sesión interactiva
combined_panel

# Dirección de destino en Dropbox
setwd("D:/PROGRAMAS/Dropbox/VAVILOV_2.0/Results/SUMMARY_FILES")

# Exportación limpia en Alta Resolución (1000 DPI)
ggsave(
  filename = "quantile_comparison_panel.png", 
  plot = combined_panel, 
  width = 12, 
  height = 10, 
  dpi = 1000
)