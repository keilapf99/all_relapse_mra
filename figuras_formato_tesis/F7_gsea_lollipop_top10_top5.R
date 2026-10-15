# Title : GSEA con Hallmarks de MSigDB - lollipop top 10 y top 5 (formato de figuras de la tesis)
# Author: Keila Adonai Pérez-Fuentes
# Basado en: los objetos lolli_top10 y lolli_top5 guardados en
#            all_relapse_masters_degree/4_GSEA_ORA/.RData (su script ya no existe;
#            generaron top10_lolliplot_GSEA.png y top5_gsea_nes.png)
# Objective: Reusar esos mismos gráficos (capas, colores, tema, títulos) y
# sustituir sus datos por los del GSEA del paso 05a de leukemia_maestria.
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/):
#   Rscript figuras_formato_tesis/F7_gsea_lollipop_top10_top5.R
# Entradas: enrichment_analysis/gsea_msigdb_hallmarks/gsea_msigdb_hallmarks.tsv (paso 05a),
#           ../../all_relapse_masters_degree/all_relapse_masters_degree/4_GSEA_ORA/.RData
# Salidas : figuras_formato_tesis/F7_gsea_lollipop_top10_top5/
#           top10_lolliplot_GSEA.png, top5_gsea_nes.png

# -------------------------------- PACKAGES ------------------------------------
library(dplyr)
library(ggplot2)

source("figuras_formato_tesis/fuente_tesis.R")

out_dir <- "figuras_formato_tesis/F7_gsea_lollipop_top10_top5"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------- DATA ----------------------------------------
gsea_res <- read.delim(
  "enrichment_analysis/gsea_msigdb_hallmarks/gsea_msigdb_hallmarks.tsv",
  quote = ""
) %>%
  filter(p.adjust < 0.05)

# Gráficos originales (solo se cargan lolli_top10 y lolli_top5)
workspace_all_relapse <- new.env()
load("../../all_relapse_masters_degree/all_relapse_masters_degree/4_GSEA_ORA/.RData",
     envir = workspace_all_relapse)

# Igual que en el original: top N por NES, con `pathway` como factor ordenado
# de mayor a menor NES (el primer nivel queda abajo en el eje y)
top_n_pathways <- function(n) {
  gsea_res %>%
    arrange(desc(NES)) %>%
    slice_head(n = n) %>%
    mutate(
      pathway = gsub("_", " ", gsub("^HALLMARK_", "", Description)),
      pathway = factor(pathway, levels = pathway)
    )
}

# --------------------------------- PLOTS --------------------------------------
lolli_top10 <- workspace_all_relapse$lolli_top10 %+% top_n_pathways(10)
lolli_top5  <- workspace_all_relapse$lolli_top5  %+% top_n_pathways(5)

ggsave(file.path(out_dir, "top10_lolliplot_GSEA.png"), lolli_top10,
       width = 10, height = 8, dpi = 300, bg = "white")
# El original se exportó desde RStudio (2514 x 1238 px, ~144 dpi); mismo tamaño
ggsave(file.path(out_dir, "top5_gsea_nes.png"), lolli_top5,
       width = 17.46, height = 8.6, dpi = 144, bg = "white")
