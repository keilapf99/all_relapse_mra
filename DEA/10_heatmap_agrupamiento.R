# Title : Agrupamiento jerárquico no supervisado de los pacientes
# Author: Keila Adonai Pérez-Fuentes
# Objective: Observar si los pacientes con y sin recaída se separan al
# agruparlos según su expresión génica.
#
# Se hacen dos agrupamientos con la misma metodología:
#   A. Firma de expresión diferencial (237 DEGs, paso 04). DESCRIPTIVO: los genes
#      se eligieron por separar recaída de no recaída en estos mismos pacientes,
#      así que la separación está favorecida por construcción (sesgo de
#      selección o "double dipping"); no es evidencia independiente.
#   B. Control no supervisado: el mismo número de genes, elegidos por ser los
#      de mayor varianza entre los 35 pacientes, sin usar la etiqueta de recaída.
#
# Metodología (ambos): matriz VST (paso 03); cada gen centrado y escalado
# (z-score) entre pacientes; agrupamiento jerárquico de genes y pacientes con
# distancia euclidiana y método de Ward (ward.D2). El dendrograma de pacientes
# se corta en 2 grupos y se compara con recaída (prueba exacta de Fisher e
# índice de Rand ajustado, ARI: 1 = coincidencia perfecta, ~0 = azar).
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/), p. ej.:
#   Rscript DEA/10_heatmap_agrupamiento.R
# Entradas: preprocess/counts_vst.tsv (paso 03), DEA/DEGS_significant_relapse.tsv (paso 04),
#           preprocess/metadata_NS500560.tsv (paso 01)
# Salidas (en DEA/heatmap_agrupamiento/):
#   - heatmap_DEGs.png / .pdf             : A. firma de expresión diferencial
#   - heatmap_genes_variables.png / .pdf  : B. genes de mayor varianza
#   - grupos_pacientes.tsv                : grupo asignado a cada paciente (A y B)
#   - resumen_agrupamiento.tsv            : concordancia de los grupos con recaída

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

out_dir <- "DEA/heatmap_agrupamiento"
dir.create(out_dir, showWarnings = FALSE)

# -------------------------------- DATA ----------------------------------------
vst <- as.matrix(read.table("preprocess/counts_vst.tsv", header = TRUE, check.names = FALSE))
degs <- read.delim("DEA/DEGS_significant_relapse.tsv", quote = "")
meta <- read.delim("preprocess/metadata_NS500560.tsv", colClasses = "character")

meta <- meta %>%
  mutate(Recaída = if_else(relapse == "si", "Con recaída", "Sin recaída"),
         `Riesgo clínico` = if_else(risk == "RA", "Alto", "Habitual"),
         Fusión = if_else(rearrange == "ND", "Sin fusión", rearrange),
         Flowcell = flowcell) %>%
  column_to_rownames("sample")
vst <- vst[, rownames(meta)]
stopifnot(identical(colnames(vst), rownames(meta)))

genes_degs <- intersect(degs$gene_symbol, rownames(vst))
n_genes <- length(genes_degs)
genes_var <- names(sort(apply(vst, 1, var), decreasing = TRUE))[seq_len(n_genes)]

cat(sprintf("%d pacientes | A: %d DEGs | B: %d genes de mayor varianza (%d también son DEGs)\n",
            ncol(vst), n_genes, length(genes_var), length(intersect(genes_var, genes_degs))))

# ------------------------------ FUNCIONES -------------------------------------
# Índice de Rand ajustado (Hubert & Arabie, 1985)
ari <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(x) sum(choose(x, 2))
  esperado <- comb2(rowSums(tab)) * comb2(colSums(tab)) / choose(sum(tab), 2)
  maximo <- (comb2(rowSums(tab)) + comb2(colSums(tab))) / 2
  (comb2(tab) - esperado) / (maximo - esperado)
}

colores_anotacion <- list(
  Recaída = c("Con recaída" = "#B2182B", "Sin recaída" = "#2166AC"),
  `Riesgo clínico` = c("Alto" = "#525252", "Habitual" = "#D9D9D9"),
  Fusión = c("Sin fusión" = "#F0F0F0", "ETV6::RUNX1" = "#E69F00",
             "TCF3::PBX1" = "#56B4E9", "BCR::ABL1" = "#009E73"),
  Flowcell = setNames(c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
                        "#FDB462", "#B3DE69"), sort(unique(meta$Flowcell))),
  Grupo = c("1" = "#7F3C8D", "2" = "#E6A141")
)

agrupar <- function(genes, etiqueta, titulo, archivo) {
  z <- t(scale(t(vst[genes, ])))
  hc_pacientes <- hclust(dist(t(z)), method = "ward.D2")
  hc_genes <- hclust(dist(z), method = "ward.D2")
  grupo <- cutree(hc_pacientes, k = 2)

  tab <- table(Grupo = grupo, Recaída = meta$Recaída)
  ft <- fisher.test(tab)
  resumen <- tibble(
    agrupamiento = etiqueta, n_genes = length(genes),
    grupo1_con_recaida = tab["1", "Con recaída"], grupo1_sin_recaida = tab["1", "Sin recaída"],
    grupo2_con_recaida = tab["2", "Con recaída"], grupo2_sin_recaida = tab["2", "Sin recaída"],
    p_fisher = ft$p.value, ARI = ari(grupo, meta$Recaída)
  )

  anotacion <- meta %>%
    select(Recaída, `Riesgo clínico`, Fusión, Flowcell) %>%
    mutate(Grupo = as.character(grupo[rownames(meta)]), .before = 1)

  subtitulo <- sprintf("Grupo 1: %d con / %d sin recaída | Grupo 2: %d con / %d sin recaída | Fisher p = %s | ARI = %.2f",
                       resumen$grupo1_con_recaida, resumen$grupo1_sin_recaida,
                       resumen$grupo2_con_recaida, resumen$grupo2_sin_recaida,
                       format.pval(resumen$p_fisher, digits = 2), resumen$ARI)

  hm <- pheatmap(z,
                 cluster_rows = hc_genes, cluster_cols = hc_pacientes,
                 cutree_cols = 2,
                 annotation_col = anotacion, annotation_colors = colores_anotacion,
                 color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
                 breaks = seq(-3, 3, length.out = 102),
                 show_rownames = FALSE, fontsize = 9, fontsize_col = 7,
                 border_color = NA, treeheight_row = 30,
                 legend_breaks = c(-3, -1.5, 0, 1.5, 3),
                 legend_labels = c("-3", "-1.5", "0", "1.5", "3"),
                 main = paste0(titulo, "\n", subtitulo, "\nColor: z-score de la expresión VST de cada gen"),
                 silent = TRUE)
  ggsave(file.path(out_dir, paste0(archivo, ".png")), hm$gtable,
         width = 11, height = 9, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(archivo, ".pdf")), hm$gtable,
         width = 11, height = 9, device = cairo_pdf)

  list(resumen = resumen,
       grupos = tibble(sample = names(grupo), !!paste0("grupo_", etiqueta) := grupo))
}

# ------------------------------ ANÁLISIS --------------------------------------
a <- agrupar(genes_degs, "A_DEGs",
             sprintf("A. Firma de expresión diferencial (%d DEGs) - descriptivo", n_genes),
             "heatmap_DEGs")
b <- agrupar(genes_var, "B_mayor_varianza",
             sprintf("B. %d genes de mayor varianza (sin usar la etiqueta de recaída)", n_genes),
             "heatmap_genes_variables")

resumen <- bind_rows(a$resumen, b$resumen)
write_tsv(resumen, file.path(out_dir, "resumen_agrupamiento.tsv"))

meta %>%
  rownames_to_column("sample") %>%
  select(sample, Recaída) %>%
  left_join(a$grupos, by = "sample") %>%
  left_join(b$grupos, by = "sample") %>%
  write_tsv(file.path(out_dir, "grupos_pacientes.tsv"))

print(as.data.frame(resumen), right = FALSE)
