# Title : Filtrado y PCA de la matriz VST (formato de figuras de la tesis)
# Author: Keila Adonai Pérez-Fuentes
# Basado en: all_relapse_masters_degree/1_preprocessing/prep_vst.r (solo las figuras)
# Objective: Mismas figuras que el original (histogramas de logCPM antes y
# después del filtrado y PCA por variable clínica), con los datos de
# leukemia_maestria. No recalcula ni guarda matrices: usa la VST del paso 03.
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/):
#   Rscript figuras_formato_tesis/F2_filtrado_pca_vst.R
# Entradas: raw_data/raw_matrix.tsv, preprocess/metadata_NS500560.tsv (paso 01),
#           preprocess/counts_vst.tsv (paso 03)
# Salidas : figuras_formato_tesis/F2_filtrado_pca_vst/
#           filtering_comparison.png, PCA_vst_relapse.png, combined_pca.png

# -------------------------------- PACKAGES ------------------------------------
library(dplyr)
library(tidyverse)
library(edgeR)
library(ggrepel)
library(ggplot2)
library(patchwork)

source("figuras_formato_tesis/fuente_tesis.R")

# ggrepel acomoda las etiquetas al azar: semilla fija para que salgan igual
set.seed(123)

out_dir <- "figuras_formato_tesis/F2_filtrado_pca_vst"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------- DATA ---------------------------------------
counts <- read.table(
  file = "raw_data/raw_matrix.tsv",
  header = T,
  check.names = F)

metadata <- read.table(
  "preprocess/metadata_NS500560.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Matriz VST del paso 03 (símbolos únicos x 35 muestras)
counts_vst_mean <- read.table(
  "preprocess/counts_vst.tsv",
  header = TRUE,
  check.names = FALSE
)

# Editing data -----------------------------------------------------------------
counts <- counts[, colnames(counts) %in% metadata$sample, drop = FALSE]
metadata$sample <- as.character(metadata$sample)
metadata$age_months <- as.numeric(
  gsub(",", ".", metadata$age_months)
)
metadata$age_group_years <- ifelse(
  metadata$age_months >= 120,
  
  "≥10",
  
  "1–9.99"
)
metadata$cl_risk <- ifelse(
  metadata$cl >= 50000,
  
  "≥50000cells/μL",
  
  "<50000cells/μL"
)

# ------------------------ 1. COUNTS MATRIX PREPROCESS -------------------------

# 1.1 Filter low counts genes --------------------------------------------------
counts_filtered <- counts[rowSums(counts >= 10) >= 15, ] #at least the nsample of the smaller condition
dim(counts_filtered)
# [1] 17044    35

# 1.2 Distribution plots -------------------------------------------------------
dge_before <- DGEList(counts = counts)
logCPM_before <- cpm(dge_before, log = T)
mean_logCPM_before <- rowMeans(logCPM_before)
df_before <- data.frame(
  mean_logCPM = mean_logCPM_before
)

distr_before <- ggplot(df_before, aes(x = mean_logCPM)) +
  
  geom_histogram(
    bins = 100,
    fill = "#7B68EE",
    color = "white",
    alpha = 0.9
  ) +
  
  theme_minimal(base_size = 15) +
  
  labs(
    title = "Distribución global de expresión génica",
    subtitle = "Antes del filtrado de genes con conteos bajos",
    x = "Promedio de logCPM",
    y = "Frecuencia"
  ) +
  
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

distr_before

dge_filtered <- DGEList(counts = counts_filtered)
logCPM_filtered <- cpm(dge_filtered, log = TRUE)
mean_logCPM_filtered <- rowMeans(logCPM_filtered)

df_after <- data.frame(
  mean_logCPM = mean_logCPM_filtered
)

distr_after <- ggplot(df_after, aes(x = mean_logCPM)) +
  
  geom_histogram(
    bins = 100,
    fill = "#20B2AA",
    color = "white",
    alpha = 0.9
  ) +
  
  theme_minimal(base_size = 15) +
  
  labs(
    subtitle = "Después del filtrado de genes con conteos bajos",
    x = "Promedio de logCPM",
    y = "Frecuencia"
  ) +
  
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

distr_after

before_after <- distr_before + distr_after
before_after
ggsave(
  filename = file.path(out_dir, "filtering_comparison.png"),
  plot = before_after,
  width = 14,
  height = 6,
  dpi = 300
)

# ------------------------ 5. PCA ANALYSIS -------------------------------------

# Transpose matrix
counts_t <- t(counts_vst_mean)

# PCA
pca <- prcomp(
  counts_t,
  center = TRUE,
  scale. = TRUE
)

# Percentage of explained variance
pca_var <- summary(pca)$importance[2, 1:2] * 100

# PCA dataframe
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  left_join(metadata, by = "sample")

# ------------------------ 6. PCA PLOTS ----------------------------------------

# 6.1 Relapse PCA --------------------------------------------------------------

pca_relapse <- ggplot(
  pca_df,
  aes(PC1, PC2, color = relapse)
) +
  
  geom_point(
    size = 4.5,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    aes(label = sample),
    size = 3.8,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  
  scale_color_manual(
    values = c(
      "si" = "purple",   # turquesa
      "no" = "cadetblue2"    # verde
    ),
    labels = c(
      "si" = "Recaída",
      "no" = "Sin recaída"
    )
  ) +
  
  labs(
    title = "PCA de muestras transformadas con VST",
    
    subtitle = "Agrupación por condición de recaída",
    
    x = paste0(
      "PC1 (",
      round(pca_var[1], 1),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(pca_var[2], 1),
      "%)"
    )
  ) +
  
  theme_minimal(base_size = 15) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 19,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 13,
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 12
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    )
  )

pca_relapse

# 6.2 Sex PCA ------------------------------------------------------------------
pca_sex <- ggplot(
  pca_df,
  aes(PC1, PC2, color = as.factor(sex))
) +
  
  geom_point(
    size = 4.5,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    aes(label = sample),
    size = 3.8,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  
  scale_color_manual(
    values = c(
      "1" = "blue",
      "2" = "pink"  
    ),
    
    labels = c(
      "1" = "Masculino",
      "2" = "Femenino"
    )
  ) +
  
  labs(
    subtitle = "Agrupación por condición de sexo",
    
    x = paste0(
      "PC1 (",
      round(pca_var[1], 1),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(pca_var[2], 1),
      "%)"
    )
  ) +
  
  theme_minimal(base_size = 15) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 19,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 13,
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 12
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    )
  )

pca_sex

# 6.3 Risk PCA ------------------------------------------------------------------
pca_risk <- ggplot(
  pca_df,
  aes(PC1, PC2, color = risk)
) +
  
  geom_point(
    size = 4.5,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    aes(label = sample),
    size = 3.8,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  
  scale_color_manual(
    values = c(
      "RA" = "red",  
      "RH" = "darkseagreen"   
    ),
    
    labels = c(
      "RA" = "Riesgo alto",
      "RH" = "Riesgo habitual"
    )
  ) +
  
  labs(
    subtitle = "Agrupación por condición de riesgo a recaída",
    
    x = paste0(
      "PC1 (",
      round(pca_var[1], 1),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(pca_var[2], 1),
      "%)"
    )
  ) +
  
  theme_minimal(base_size = 15) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 19,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 13,
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 12
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    )
  )

pca_risk

# 6.4 Rearrange PCA ------------------------------------------------------------
pca_rearrange <- ggplot(
  pca_df,
  aes(PC1, PC2, color = rearrange)
) +
  
  geom_point(
    size = 4.5,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    aes(label = sample),
    size = 3.8,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  
  scale_color_manual(
    values = c(
      "BCR::ABL1" = "red",
      "TCF3::PBX1" = "orange",
      "ETV6::RUNX1" = "yellow",   
      "ND" = "darkseagreen"    
    )
  ) +
  
  labs(
    subtitle = "Agrupación por condición de rearreglo génico",
    
    x = paste0(
      "PC1 (",
      round(pca_var[1], 1),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(pca_var[2], 1),
      "%)"
    )
  ) +
  
  theme_minimal(base_size = 15) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 19,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 13,
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 12
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    )
  )

pca_rearrange

# 6.6 Age pca ------------------------------------------------------------------
pca_age <- ggplot(
  pca_df,
  aes(PC1, PC2, color = age_group_years)
) +
  
  geom_point(
    size = 4.5,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    aes(label = sample),
    size = 3.8,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  
  scale_color_manual(
    values = c(
      "≥10" = "red",
      "1–9.99" = "darkseagreen"
    ),
    labels = c(
      "≥10" = "Riesgo alto por edad",
      "1–9.99" = "Riesgo habitual por edad"
    )
  ) +
  
  labs(
    subtitle = "Agrupación por condición de riesgo por edad",
    
    x = paste0(
      "PC1 (",
      round(pca_var[1], 1),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(pca_var[2], 1),
      "%)"
    )
  ) +
  
  theme_minimal(base_size = 15) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 19,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 13,
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 12
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    )
  )

pca_age

# 6.7 WBC pca ------------------------------------------------------------------
pca_wbc <- ggplot(
  pca_df,
  aes(PC1, PC2, color = cl_risk)
) +
  
  geom_point(
    size = 4.5,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    aes(label = sample),
    size = 3.8,
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  
  scale_color_manual(
    
    values = c(
      "≥50000cells/μL" = "red",
      "<50000cells/μL" = "darkseagreen"
    ),
    
    labels = c(
      "≥50000cells/μL" = "Riesgo alto leucocitario",
      
      "<50000cells/μL" = "Riesgo habitual leucocitario"
    )
  )+
  
  labs(
    subtitle = "Agrupación por condición de riesgo por cuenta leucocitaria en sangre periferica",
    
    x = paste0(
      "PC1 (",
      round(pca_var[1], 1),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(pca_var[2], 1),
      "%)"
    )
  ) +
  
  theme_minimal(base_size = 15) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 19,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 13,
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 12
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    )
  )

pca_wbc

# 6.8 Patchwork PCAs -----------------------------------------------------------
combined_pca <- (
  pca_relapse +
    pca_sex+
    pca_rearrange
) / (
  pca_risk +
    pca_wbc +
    pca_age
)

# 6.9 Save PCAs ----------------------------------------------------------------
ggsave(
  filename = file.path(out_dir, "PCA_vst_relapse.png"),
  plot = pca_relapse,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "combined_pca.png"),
  plot = combined_pca,
  width = 22,
  height = 14,
  dpi = 300
)

