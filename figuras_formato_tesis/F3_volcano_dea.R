# Title : Volcano plot de la expresión diferencial (formato de figuras de la tesis)
# Author: Keila Adonai Pérez-Fuentes
# Basado en: all_relapse_masters_degree/2_DEA/2_DEA_genes_of_interest.r y 2_DEA.r
# Objective: Mismos volcano plots que el original, con la firma de expresión
# diferencial de leukemia_maestria (paso 04). No recalcula el DEA.
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/):
#   Rscript figuras_formato_tesis/F3_volcano_dea.R
# Entradas: DEA/DEA_totales_relapse.tsv (paso 04),
#           figuras_formato_tesis/datos/degs_procesobiologico.tsv (genes etiquetados
#           por proceso biológico; copiado de all_relapse_masters_degree/2_DEA/)
# Salidas : figuras_formato_tesis/F3_volcano_dea/
#   - DEA_relapse_vp.png           : genes etiquetados por proceso biológico (14 x 8 in)
#   - DEA_relapse_vp_top25.png     : 25 genes con mayor |log2FC| etiquetados (11 x 7 in)
#   - DEA_relapse_vp_sinlabels.png : sin etiquetas (11 x 7 in)
#   - top_genes_25_logfc.tsv, genes_of_interest_dea.tsv

# --------------------------- PACKAGES -----------------------------------------
library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(dplyr)

source("figuras_formato_tesis/fuente_tesis.R")

# ggrepel acomoda las etiquetas al azar: semilla fija para que salgan igual
set.seed(123)

out_dir <- "figuras_formato_tesis/F3_volcano_dea"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------- DATA ---------------------------------------
res_order_gsymbols <- read.delim("DEA/DEA_totales_relapse.tsv", quote = "")

# 3.4 Biological process labels -----------------------------------------------
bp_data <- read.table(
  "figuras_formato_tesis/datos/degs_procesobiologico.tsv",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
bp_data <- bp_data[!duplicated(bp_data$DEG), ]

# ----------------------------- STEP 4. Plots ----------------------------------
# 4.1 Label regulation direction -----------------------------------------------
res_order_gsymbols$Type <- ifelse(
  res_order_gsymbols$log2FoldChange > 1 & res_order_gsymbols$padj < 0.05,
  "Sobreexpresados",
  ifelse(
    res_order_gsymbols$log2FoldChange < -1 & res_order_gsymbols$padj < 0.05,
    "Subexpresados",
    "No significativos"
  )
)

# 4.2 Edit design -------------------------------------------------------------
res_order_gsymbols_ns <- res_order_gsymbols %>%
  filter(Type != "Not Significant")

gene_counts <- table(res_order_gsymbols$Type)

color_labels <- c(
  paste0("Sobreexpresados (", gene_counts["Sobreexpresados"], ")"),
  paste0("Subexpresados (", gene_counts["Subexpresados"], ")"),
  paste0("No significativos (", gene_counts["No significativos"], ")")
)

names(color_labels) <- c("Sobreexpresados", "Subexpresados", "No significativos")
top_genes <- res_order_gsymbols %>%
  filter(Type != "No significativos") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 25)

genes_bp_label <- res_order_gsymbols %>%
  filter(gene_symbol %in% bp_data$DEG) %>%
  left_join(bp_data, by = c("gene_symbol" = "DEG"))

bp_colors <- c(
  "Stem_like_Inmadurez_hematopoyetica"    = "#8B4789",   # morado
  "Microambiente_protector_Nicho_medular"  = "#E69500",   # ambar
  "Reprogramacion_metabolica"              = "#2E7D32",   # verde bosque
  "Tolerancia_terapeutica_Resistencia"     = "#6D4C41",   # cafe oscuro
  "Inmunomodulacion_Inmunosupresion"       = "#00838F"    # teal oscuro
)

bp_labels <- c(
  "Stem_like_Inmadurez_hematopoyetica"    = "Stem-like / Inmadurez hematopoyética",
  "Microambiente_protector_Nicho_medular"  = "Microambiente protector / Nicho medular",
  "Reprogramacion_metabolica"             = "Reprogramación metabólica",
  "Tolerancia_terapeutica_Resistencia"    = "Tolerancia terapéutica / Resistencia",
  "Inmunomodulacion_Inmunosupresion"      = "Inmunomodulación / Inmunosupresión"
)

# 4.3 Volcano plot con genes por proceso biológico (DEA_relapse_vp.png) -------------------------------------------------------------
volcano_plot <- ggplot(
  res_order_gsymbols,
  aes(x = log2FoldChange, y = -log10(padj), color = Type)
) +
  geom_point(aes(alpha = Type), size = 2.5) +

  scale_alpha_manual(values = c(
    "No significativos" = 0.25,
    "Sobreexpresados" = 0.9,
    "Subexpresados" = 0.9
  )) +
  guides(alpha = "none") +

  scale_color_manual(
    name = NULL,
    values = c(
      "Sobreexpresados" = "#E64B35",
      "Subexpresados" = "#4DBBD5",
      "No significativos" = "grey"
    ),
    labels = color_labels,
    guide  = guide_legend(order = 2)
  ) +

  new_scale_color() +

  geom_label_repel(
    data = genes_bp_label,
    aes(label = gene_symbol, color = Proceso_biologico_principal),
    size = 4.5,
    fontface = "bold",
    fill = alpha("white", 0.88),
    label.size = 0.5,
    label.padding = unit(0.22, "lines"),
    box.padding = 0.6,
    point.padding = 0.4,
    segment.size = 0.65,
    segment.alpha = 1,
    min.segment.length = 0,
    max.overlaps = 50
  ) +

  scale_color_manual(
    name   = "Proceso biológico",
    values = bp_colors,
    labels = bp_labels,
    guide  = guide_legend(order = 1)
  ) +

  labs(
    title    = "Genes diferencialmente expresados al diagnóstico",
    subtitle = "recaída vs no recaída",
    x        = expression(log[2]~Fold~Change),
    y        = expression(-log[10]~valor~p~ajustado)
  ) +

  theme_minimal(base_size = 15) +
  theme(
    plot.title       = element_text(size = 20, hjust = 0.5, face = "bold.italic"),
    plot.subtitle    = element_text(size = 15, hjust = 0.5, face = "italic"),
    axis.title.x     = element_text(size = 12),
    axis.title.y     = element_text(size = 12),
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(size = 10),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.4, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black")
  ) +

  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey90") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray90")

volcano_plot

# 4.4 Volcano plot con los 25 genes de mayor |log2FC| (DEA_relapse_vp_top25.png) -------------------------------------------------------------
volcano_plot_top25 <- ggplot(
  res_order_gsymbols,
  aes(x = log2FoldChange, y = -log10(padj), color = Type)
) +
  geom_point(aes(alpha = Type), size = 2.5)+
  
  scale_alpha_manual(values = c(
    "No significativos" = 0.25,
    "Sobreexpresados" = 0.9,
    "Subexpresados" = 0.9
  )) + 
  guides(alpha = "none") +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_symbol),
    color = "black",
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = 50
  ) +
  
  scale_color_manual(
    values = c(
      "Sobreexpresados" = "#E64B35",
      "Subexpresados" = "#4DBBD5",
      "No significativos" = "grey"
    ),
    labels = color_labels
  ) +
  
  labs(
    title = "Genes diferencialmente expresados al diagnóstico",
    subtitle = "recaída vs no recaída",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~valor~p~ajustado),
    color = NULL
  ) +
  
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold.italic"),
    plot.subtitle = element_text(size = 15, hjust = 0.5, face = "italic"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey90") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray90")

volcano_plot_top25

# 4.5 Volcano plot sin etiquetas (DEA_relapse_vp_sinlabels.png) ---------------
volcano_plot_sinlabels <- volcano_plot_top25
volcano_plot_sinlabels$layers <- Filter(
  function(l) !inherits(l$geom, "GeomTextRepel"),
  volcano_plot_sinlabels$layers
)
volcano_plot_sinlabels <- volcano_plot_sinlabels +
  theme(legend.position = "right")

# ------------------------------- STEP 5. Save ---------------------------------
write.table(top_genes, file = file.path(out_dir, "top_genes_25_logfc.tsv"),
            sep = "\t",
            quote = F,
            col.names = T)
ggsave(file.path(out_dir, "DEA_relapse_vp.png"), volcano_plot, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "DEA_relapse_vp_top25.png"), volcano_plot_top25, width = 11, height = 7, dpi = 300)
ggsave(file.path(out_dir, "DEA_relapse_vp_sinlabels.png"), volcano_plot_sinlabels, width = 11, height = 7, dpi = 300)

# ------------------------------- STEP 6. Genes of interest --------------------
# Subset res_filter_order_gsymbols (DEGs: |log2FC| >= 1, padj < 0.05) to
# retain only the genes of interest. Genes absent from res_filter_order_gsymbols
# (i.e. not passing the DEG filters) will not appear in the output.

genes_of_interest <- c(
  "CXCL12", "VCAM1",   "THBS1",   "TGFBI",   "NAMPT",
  "CD36",   "PBX1",    "MEIS1",   "TAL1",    "FZD5",
  "WIPI1",  "BNIP3L",  "HAVCR2",  "CD200R1", "LGALS3",
  "CYP1B1", "TTK",     "ICOSLG",  "TNF",     "MLLT3"
)

genes_of_interest_dea <- res_order_gsymbols %>%
  filter(Type != "No significativos") %>%
  filter(gene_symbol %in% genes_of_interest) %>%
  dplyr::select(gene_symbol, log2FoldChange, stat, padj) %>%
  arrange(padj)

write.table(
  genes_of_interest_dea,
  file      = file.path(out_dir, "genes_of_interest_dea.tsv"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

message(
  "Genes of interest found in DEGs: ",
  nrow(genes_of_interest_dea), " / ", length(genes_of_interest), "\n",
  "Saved: genes_of_interest_dea.tsv"
)
