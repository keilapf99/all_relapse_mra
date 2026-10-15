# Title : GSEA con Hallmarks de MSigDB (formato de figuras de la tesis)
# Author: Keila Adonai Pérez-Fuentes
# Basado en: all_relapse_masters_degree/4_GSEA_ORA/gsea_hallmarks_plots.r,
#            gsea_hallmarks_plots_centered.r, gsea_hallmarks_top5_publication.r
#            y gsea_hallmarks_top5_publication_es.r (solo las figuras)
# Objective: Mismas figuras que el original, con los resultados del GSEA del
# paso 05a de leukemia_maestria. No recalcula el GSEA.
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/):
#   Rscript figuras_formato_tesis/F5_gsea_msigdb_hallmarks.R
# Entradas: enrichment_analysis/gsea_msigdb_hallmarks/gsea_msigdb_hallmarks.tsv (paso 05a)
# Salidas : figuras_formato_tesis/F5_gsea_msigdb_hallmarks/
#   - Hallmark_relapse_mechanisms.{png,pdf}        : hallmarks agrupados por mecanismo
#   - Hallmark_relapse_mechanisms_centered.png     : misma figura con filas centradas
#   - top5_hallmarks_lolliplot_publication{,_es}.{png,pdf} : top 5 por |NES| (inglés / español)

# -------------------------------- PACKAGES ------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)

source("figuras_formato_tesis/fuente_tesis.R")

out_dir <- "figuras_formato_tesis/F5_gsea_msigdb_hallmarks"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------- DATA ----------------------------------------
# En all_relapse el GSEA se corría con pvalueCutoff = 0.05, así que gsea_res
# solo traía los hallmarks significativos; el paso 05a guarda los 50, por eso
# aquí se filtra por FDR < 0.05.
gsea_res <- read.delim(
  "enrichment_analysis/gsea_msigdb_hallmarks/gsea_msigdb_hallmarks.tsv",
  quote = ""
) %>%
  filter(p.adjust < 0.05)

# ============================================================
# PREPARAR DATOS
# ============================================================

plot_df <- gsea_res %>%
  
  mutate(
    
    Hallmark = gsub("^HALLMARK_", "", Description),
    
    Hallmark_label = case_when(
      
      Hallmark == "EPITHELIAL_MESENCHYMAL_TRANSITION" ~
        "Epithelial Mesenchymal Transition",
      
      Hallmark == "APICAL_JUNCTION" ~
        "Apical Junction",
      
      Hallmark == "IL2_STAT5_SIGNALING" ~
        "IL2 STAT5 Signaling",
      
      Hallmark == "WNT_BETA_CATENIN_SIGNALING" ~
        "WNT Beta Catenin Signaling",
      
      Hallmark == "COMPLEMENT" ~
        "Complement",
      
      Hallmark == "INFLAMMATORY_RESPONSE" ~
        "Inflammatory Response",
      
      Hallmark == "IL6_JAK_STAT3_SIGNALING" ~
        "IL6 JAK STAT3 Signaling",
      
      Hallmark == "XENOBIOTIC_METABOLISM" ~
        "Xenobiotic Metabolism",
      
      Hallmark == "MTORC1_SIGNALING" ~
        "mTORC1 Signaling",
      
      Hallmark == "HEME_METABOLISM" ~
        "Heme Metabolism",
      
      Hallmark == "FATTY_ACID_METABOLISM" ~
        "Fatty Acid Metabolism",
      
      Hallmark == "OXIDATIVE_PHOSPHORYLATION" ~
        "Oxidative Phosphorylation",
      
      Hallmark == "MYC_TARGETS_V1" ~
        "MYC Targets V1",
      
      Hallmark == "KRAS_SIGNALING_UP" ~
        "KRAS Signaling Up",
      
      TRUE ~ Hallmark
      
    ),
    
    Mechanism = case_when(
      
      Hallmark %in% c(
        "EPITHELIAL_MESENCHYMAL_TRANSITION",
        "APICAL_JUNCTION"
      ) ~ "Unlocking Phenotypic\nPlasticity",
      
      Hallmark %in% c(
        "IL2_STAT5_SIGNALING",
        "WNT_BETA_CATENIN_SIGNALING"
      ) ~ "Plasticity / Cellular\nStemness",
      
      Hallmark %in% c(
        "COMPLEMENT",
        "INFLAMMATORY_RESPONSE",
        "IL6_JAK_STAT3_SIGNALING"
      ) ~ "Tumor-Promoting\nInflammation",
      
      Hallmark %in% c(
        "XENOBIOTIC_METABOLISM",
        "MTORC1_SIGNALING"
      ) ~ "Resisting Cell Death",
      
      Hallmark %in% c(
        "HEME_METABOLISM",
        "FATTY_ACID_METABOLISM",
        "OXIDATIVE_PHOSPHORYLATION"
      ) ~ "Deregulating Cellular\nEnergetics",
      
      Hallmark %in% c(
        "MYC_TARGETS_V1",
        "KRAS_SIGNALING_UP"
      ) ~ "Sustaining Proliferative\nSignaling",
      
      TRUE ~ NA_character_
    )
    
  ) %>%
  
  filter(!is.na(Mechanism))

# ============================================================
# ORDEN
# ============================================================

group_order <- c(
  "Sustaining Proliferative\nSignaling",
  "Deregulating Cellular\nEnergetics",
  "Resisting Cell Death",
  "Tumor-Promoting\nInflammation",
  "Plasticity / Cellular\nStemness",
  "Unlocking Phenotypic\nPlasticity"
)

plot_df$Mechanism <- factor(
  plot_df$Mechanism,
  levels = group_order
)

plot_df <- plot_df %>%
  arrange(Mechanism, desc(NES))

plot_df$y <- rev(seq_len(nrow(plot_df)))

# ============================================================
# CENTROS DE LOS BLOQUES
# ============================================================

group_centers <- plot_df %>%
  group_by(Mechanism) %>%
  summarise(
    y = mean(y),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      Mechanism,
      "\n(n = ", n, ")"
    )
  )

# ============================================================
# SEPARADORES
# ============================================================

group_lines <- plot_df %>%
  group_by(Mechanism) %>%
  summarise(
    ymin = min(y),
    .groups = "drop"
  )

# ============================================================
# COLORES
# ============================================================

cols <- c(
  
  "Unlocking Phenotypic\nPlasticity" =
    "#E41A1C",
  
  "Plasticity / Cellular\nStemness" =
    "#7F0000",
  
  "Tumor-Promoting\nInflammation" =
    "#5E3C99",
  
  "Resisting Cell Death" =
    "#E08E00",
  
  "Deregulating Cellular\nEnergetics" =
    "#2166AC",
  
  "Sustaining Proliferative\nSignaling" =
    "#2CA25F"
)

# ============================================================
# PANEL IZQUIERDO
# ============================================================

p_left <- ggplot(
  group_centers,
  aes(x = 1, y = y)
) +
  
  geom_text(
    aes(
      label = label,
      colour = Mechanism
    ),
    fontface = "bold",
    size = 5
  ) +
  
  scale_colour_manual(values = cols) +
  
  xlim(0,2) +
  
  theme_void() +
theme(
  legend.position = "none",
  
  plot.title = element_text(
    family = "Arial",      # tipo de letra
    face = "bold",         # bold, italic, plain
    size = 16,             # tamaño
    hjust = 0.5            # 0=izq, 0.5=centro, 1=derecha
  )
) +
  
  labs(
    title = "Mecanismo biológico"
  )

# ============================================================
# PANEL CENTRAL
# ============================================================

p_mid <- ggplot(
  plot_df,
  aes(y = y)
) +
  
  geom_tile(
    aes(
      x = 0.15,
      fill = Mechanism
    ),
    width = 0.08,
    height = 0.9
  ) +
  
  geom_text(
    aes(
      x = 0.35,
      label = Hallmark_label
    ),
    hjust = 0,
    size = 5
  ) +
  
  scale_fill_manual(values = cols) +
  
  xlim(0,4.5) +
  
  theme_void() +
  
  theme(
    legend.position = "none",
    
    plot.title = element_text(
      family = "Arial",
      face = "bold",
      size = 16,
      hjust = 0.5
    )
  )+
  
  labs(
    title = "Término (Hallmark gene set)"
  )

# ============================================================
# FOREST PLOT
# ============================================================

p_right <- ggplot(
  plot_df,
  aes(
    x = NES,
    y = y,
    colour = Mechanism
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey40"
  ) +
  
  geom_segment(
    aes(
      x = 0,
      xend = NES,
      yend = y
    ),
    linewidth = 1
  ) +
  
  geom_point(
    aes(
      size = -log10(p.adjust)
    )
  ) +
  
  geom_text(
    aes(
      label = sprintf("%.2f", NES)
    ),
    hjust = ifelse(plot_df$NES > 0, -0.8, 1.8),
    size = 4
  ) +
  
  geom_hline(
    data = group_lines,
    aes(yintercept = ymin - 0.5),
    colour = "grey80",
    linetype = 2
  ) +
  
  scale_colour_manual(values = cols) +
  
  scale_size_continuous(
    range = c(2,8),
    name = "p.adjust"
  ) +
  
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(face = "bold")
  ) +
  
  guides(
    colour = "none"
  )

# ============================================================
# FIGURA FINAL
# ============================================================

final_plot <- (p_left | p_mid | p_right) +
  plot_layout(
    widths = c(2.5, 4.5, 6)
  )

final_plot

# ============================================================
# EXPORTAR
# ============================================================

ggsave(
  file.path(out_dir, "Hallmark_relapse_mechanisms.png"),
  final_plot,
  width = 18,
  height = 9,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Hallmark_relapse_mechanisms.pdf"),
  final_plot,
  width = 18,
  height = 9,
  bg = "white"
)

# ============================================================
# VERSIÓN CENTRADA (gsea_hallmarks_plots_centered.r)
# ============================================================

# Usar el mismo rango vertical en los tres paneles para que las filas,
# los nombres de categorías y el forest plot queden centrados entre sí.
common_y_limits <- c(0.5, nrow(plot_df) + 0.5)

# ============================================================
# CENTROS DE LOS BLOQUES
# ============================================================

group_centers <- plot_df %>%
  group_by(Mechanism) %>%
  summarise(
    y = mean(y),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      Mechanism,
      "\n(n = ", n, ")"
    )
  )

# ============================================================
# SEPARADORES
# ============================================================

group_lines <- plot_df %>%
  group_by(Mechanism) %>%
  summarise(
    ymin = min(y),
    .groups = "drop"
  )

# ============================================================
# COLORES
# ============================================================

cols <- c(
  
  "Unlocking Phenotypic\nPlasticity" =
    "#E41A1C",
  
  "Plasticity / Cellular\nStemness" =
    "#7F0000",
  
  "Tumor-Promoting\nInflammation" =
    "#5E3C99",
  
  "Resisting Cell Death" =
    "#E08E00",
  
  "Deregulating Cellular\nEnergetics" =
    "#2166AC",
  
  "Sustaining Proliferative\nSignaling" =
    "#2CA25F"
)

# ============================================================
# PANEL IZQUIERDO
# ============================================================

p_left_centered <- ggplot(
  group_centers,
  aes(x = 1, y = y)
) +
  
  geom_text(
    aes(
      label = label,
      colour = Mechanism
    ),
    fontface = "bold",
    size = 5
  ) +
  
  scale_colour_manual(values = cols) +

  scale_y_continuous(
    limits = common_y_limits,
    expand = expansion(mult = 0)
  ) +
  
  xlim(0,2) +
  
  theme_void() +
theme(
  legend.position = "none",
  
  plot.title = element_text(
    family = "Arial",      # tipo de letra
    face = "bold",         # bold, italic, plain
    size = 16,             # tamaño
    hjust = 0.5            # 0=izq, 0.5=centro, 1=derecha
  )
) +
  
  labs(
    title = "Mecanismo biológico"
  )

# ============================================================
# PANEL CENTRAL
# ============================================================

p_mid_centered <- ggplot(
  plot_df,
  aes(y = y)
) +
  
  geom_tile(
    aes(
      x = 0.15,
      fill = Mechanism
    ),
    width = 0.08,
    height = 0.9
  ) +
  
  geom_text(
    aes(
      x = 0.35,
      label = Hallmark_label
    ),
    hjust = 0,
    size = 5
  ) +
  
  scale_fill_manual(values = cols) +

  scale_y_continuous(
    limits = common_y_limits,
    expand = expansion(mult = 0)
  ) +
  
  xlim(0,4.5) +
  
  theme_void() +
  
  theme(
    legend.position = "none",
    
    plot.title = element_text(
      family = "Arial",
      face = "bold",
      size = 16,
      hjust = 0.5
    )
  )+
  
  labs(
    title = "Término (Hallmark gene set)"
  )

# ============================================================
# FOREST PLOT
# ============================================================

p_right_centered <- ggplot(
  plot_df,
  aes(
    x = NES,
    y = y,
    colour = Mechanism
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey40"
  ) +
  
  geom_segment(
    aes(
      x = 0,
      xend = NES,
      yend = y
    ),
    linewidth = 1
  ) +
  
  geom_point(
    aes(
      size = -log10(p.adjust)
    )
  ) +
  
  geom_text(
    aes(
      label = sprintf("%.2f", NES)
    ),
    hjust = ifelse(plot_df$NES > 0, -0.8, 1.8),
    size = 4
  ) +
  
  geom_hline(
    data = group_lines,
    aes(yintercept = ymin - 0.5),
    colour = "grey80",
    linetype = 2
  ) +
  
  scale_colour_manual(values = cols) +

  scale_y_continuous(
    limits = common_y_limits,
    expand = expansion(mult = 0)
  ) +

  scale_x_continuous(
    expand = expansion(mult = c(0.12, 0.12))
  ) +
  
  scale_size_continuous(
    range = c(2,8),
    name = "p.adjust"
  ) +
  
  labs(
    title = " ",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 16),
    legend.title = element_text(face = "bold")
  ) +
  
  guides(
    colour = "none"
  )

# ============================================================
# FIGURA FINAL
# ============================================================

final_plot_centered <- (p_left_centered | p_mid_centered | p_right_centered) +
  plot_layout(
    widths = c(2.5, 4.5, 6)
  )

final_plot_centered

ggsave(
  file.path(out_dir, "Hallmark_relapse_mechanisms_centered.png"),
  final_plot_centered,
  width = 18,
  height = 9,
  dpi = 600,
  bg = "white"
)

# ============================================================
# TOP 5 HALLMARKS (inglés)
# ============================================================

top5_hallmarks <- gsea_res %>%
  filter(!is.na(NES), !is.na(p.adjust), p.adjust < 0.05) %>%
  mutate(
    abs_NES = abs(NES),
    hallmark = Description %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_title(),
    minus_log10_fdr = -log10(p.adjust),
    nes_label = sprintf("%.2f", NES)
  ) %>%
  slice_max(order_by = abs_NES, n = 5, with_ties = FALSE) %>%
  arrange(NES) %>%
  mutate(hallmark = factor(hallmark, levels = hallmark))

if (nrow(top5_hallmarks) < 5) {
  stop("Fewer than five significant Hallmark gene sets were found.")
}

top5_lolliplot_publication <- ggplot(
  top5_hallmarks,
  aes(y = hallmark)
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.6,
    linetype = "dashed",
    colour = "#727272"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = NES,
      yend = hallmark
    ),
    linewidth = 1.5,
    lineend = "round",
    colour = "#6A1B9A"
  ) +
  geom_point(
    aes(
      x = NES,
      size = minus_log10_fdr
    ),
    shape = 21,
    stroke = 1.1,
    colour = "#4A145F",
    fill = "#8E44AD"
  ) +
  geom_text(
    aes(
      x = NES + if_else(NES >= 0, 0.14, -0.14),
      label = nes_label,
      hjust = if_else(NES >= 0, 0, 1)
    ),
    size = 4.2,
    fontface = "bold",
    colour = "#2B2B2B"
  ) +
  scale_size_continuous(
    range = c(5, 10),
    breaks = scales::pretty_breaks(n = 3),
    name = expression(-log[10]("FDR"))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.08, 0.28))
  ) +
  labs(
    title = "Top Hallmark Gene Sets Enriched in Relapse",
    subtitle = "Five significant gene sets ranked by absolute normalized enrichment score",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    caption = "Point size: adjusted significance (FDR)"
  ) +
  guides(
    size = guide_legend(order = 1)
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_family = "sans", base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      colour = "#25102D",
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 12,
      colour = "#595959",
      margin = margin(b = 18)
    ),
    axis.title.x = element_text(
      face = "bold",
      colour = "#333333",
      margin = margin(t = 10)
    ),
    axis.text.x = element_text(colour = "#4A4A4A"),
    axis.text.y = element_text(
      face = "bold",
      colour = "#333333",
      margin = margin(r = 10)
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#E3E3E3", linewidth = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", colour = "#333333"),
    legend.text = element_text(colour = "#4A4A4A"),
    plot.caption = element_text(
      colour = "#666666",
      hjust = 0,
      margin = margin(t = 14)
    ),
    plot.margin = margin(18, 30, 16, 18)
  )

top5_lolliplot_publication

ggsave(
  file.path(out_dir, "top5_hallmarks_lolliplot_publication.png"),
  top5_lolliplot_publication,
  width = 12,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "top5_hallmarks_lolliplot_publication.pdf"),
  top5_lolliplot_publication,
  width = 12,
  height = 7,
  bg = "white"
)

# ============================================================
# TOP 5 HALLMARKS (español)
# ============================================================

top5_hallmarks_es <- gsea_res %>%
  filter(!is.na(NES), !is.na(p.adjust), p.adjust < 0.05) %>%
  mutate(
    abs_NES = abs(NES),
    hallmark = Description %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_title(),
    minus_log10_fdr = -log10(p.adjust),
    nes_label = sprintf("%.2f", NES)
  ) %>%
  slice_max(order_by = abs_NES, n = 5, with_ties = FALSE) %>%
  arrange(NES) %>%
  mutate(
    hallmark = recode(
      hallmark,
      "Heme Metabolism" = "Metabolismo del hemo",
      "Complement" = "Complemento",
      "Epithelial Mesenchymal Transition" = "Transición epitelio-mesénquima",
      "Protein Secretion" = "Secreción de proteínas",
      "Myc Targets V1" = "Dianas de MYC V1"
    ),
    hallmark = factor(hallmark, levels = hallmark)
  )

if (nrow(top5_hallmarks_es) < 5) {
  stop("Fewer than five significant Hallmark gene sets were found.")
}

top5_lolliplot_publication_es <- ggplot(
  top5_hallmarks_es,
  aes(y = hallmark)
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.6,
    linetype = "dashed",
    colour = "#727272"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = NES,
      yend = hallmark
    ),
    linewidth = 1.5,
    lineend = "round",
    colour = "#6A1B9A"
  ) +
  geom_point(
    aes(
      x = NES,
      size = minus_log10_fdr
    ),
    shape = 21,
    stroke = 1.1,
    colour = "#4A145F",
    fill = "#8E44AD"
  ) +
  geom_text(
    aes(
      x = NES + if_else(NES >= 0, 0.14, -0.14),
      label = nes_label,
      hjust = if_else(NES >= 0, 0, 1)
    ),
    size = 4.2,
    fontface = "bold",
    colour = "#2B2B2B"
  ) +
  scale_size_continuous(
    range = c(5, 10),
    breaks = scales::pretty_breaks(n = 3),
    name = expression(-log[10]("FDR"))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.08, 0.28))
  ) +
  labs(
    title = "Principales conjuntos de genes Hallmark enriquecidos en recaída",
    subtitle = "Cinco conjuntos significativos ordenados por el valor absoluto del NES",
    x = "Puntuación normalizada de enriquecimiento (NES)",
    y = NULL,
    caption = "Tamaño del punto: significancia ajustada (FDR)"
  ) +
  guides(
    size = guide_legend(order = 1)
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_family = "sans", base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      colour = "#25102D",
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 12,
      colour = "#595959",
      margin = margin(b = 18)
    ),
    axis.title.x = element_text(
      face = "bold",
      colour = "#333333",
      margin = margin(t = 10)
    ),
    axis.text.x = element_text(colour = "#4A4A4A"),
    axis.text.y = element_text(
      face = "bold",
      colour = "#333333",
      margin = margin(r = 10)
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#E3E3E3", linewidth = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", colour = "#333333"),
    legend.text = element_text(colour = "#4A4A4A"),
    plot.caption = element_text(
      colour = "#666666",
      hjust = 0,
      margin = margin(t = 14)
    ),
    plot.margin = margin(18, 30, 16, 18)
  )

top5_lolliplot_publication_es

ggsave(
  file.path(out_dir, "top5_hallmarks_lolliplot_publication_es.png"),
  top5_lolliplot_publication_es,
  width = 12,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "top5_hallmarks_lolliplot_publication_es.pdf"),
  top5_lolliplot_publication_es,
  width = 12,
  height = 7,
  bg = "white"
)
