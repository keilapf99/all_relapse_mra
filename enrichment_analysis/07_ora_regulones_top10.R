# Title : Análisis de sobrerrepresentación (ORA) de los regulones de los
#         reguladores maestros
# Author: Keila Adonai Pérez-Fuentes
# Objective: Evaluar si los genes de los regulones de los 10 reguladores
# maestros principales están sobrerrepresentados en la colección Hallmark (H)
# de MSigDB y en los Cancer Hallmarks de Menyhart et al. (2025), conjuntos core
# (genes respaldados por >=2 bases de datos) e integrated (>=1 base de datos).
#
# - Reguladores: los 10 TFs con mayor |NES| entre los significativos
#   (FDR < 0.05) del MRA (paso 06), red de pacientes con recaída.
# - Genes de cada regulón: todos sus genes blanco (activados y reprimidos).
# - Universo: todos los genes de la matriz VST con la que se construyó la red
#   (solo esos genes podían formar parte de un regulón).
# - Prueba hipergeométrica (clusterProfiler::enricher), p ajustada por
#   Benjamini-Hochberg dentro de cada regulón; significativo si FDR < 0.05.
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/), p. ej.:
#   Rscript enrichment_analysis/07_ora_regulones_top10.R
# Entradas: corto_mra/mra_relapse.RData (paso 06), preprocess/matrix_relapse_vst.tsv (paso 03),
#           enrichment_analysis/gsea_msigdb_hallmarks/gene_sets/msigdb_hallmarks_H.tsv (paso 05a),
#           enrichment_analysis/gsea_cancer_hallmarks_menyhart/gene_sets/... (paso 05b)
# Salidas (en enrichment_analysis/ora_regulones_top10/):
#   - top10_reguladores_maestros.tsv : TFs analizados y tamaño de su regulón
#   - ora_regulones_top10.tsv        : resultados completos (ambas colecciones)
#   - ORA_regulones_top10.png / .pdf : conjuntos significativos por regulón

# ------------------------------ PACKAGES --------------------------------------
# libglpk (requerida por igraph -> enrichplot -> clusterProfiler) no está en el
# sistema de este servidor; se carga desde un entorno conda del usuario si existe.
glpk <- "/scratch/home/aperez/.conda_envs/glpk/lib/libglpk.so.40"
if (file.exists(glpk)) dyn.load(glpk, local = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(patchwork)
})

out_dir <- "enrichment_analysis/ora_regulones_top10"
dir.create(out_dir, showWarnings = FALSE)

# -------------------------------- DATA ----------------------------------------
load("corto_mra/mra_relapse.RData")  # alls, regulon_relapse, ...

universo <- rownames(read.table("preprocess/matrix_relapse_vst.tsv", header = TRUE,
                                check.names = FALSE))

leer_menyhart <- function(geneset) {
  readLines(sprintf("enrichment_analysis/gsea_cancer_hallmarks_menyhart/gene_sets/Menyhart_JPA_CancerHallmarks_%s.txt",
                    geneset), warn = FALSE) %>%
    strsplit("\t") %>%
    map_dfr(\(t) { t <- t[nzchar(t)]; tibble(term = t[1], gene = unique(t[-1])) })
}

colecciones <- list(
  "MSigDB Hallmarks" = read_tsv("enrichment_analysis/gsea_msigdb_hallmarks/gene_sets/msigdb_hallmarks_H.tsv",
                                show_col_types = FALSE) %>%
    transmute(term = gs_name, gene = gene_symbol),
  "Cancer Hallmarks (Menyhart, core)" = leer_menyhart("core"),
  "Cancer Hallmarks (Menyhart, integrated)" = leer_menyhart("integrated")
)

# ------------------------ STEP 1. Top 10 reguladores --------------------------
top10 <- alls %>%
  filter(FDR < 0.05) %>%
  slice_max(abs(NES), n = 10) %>%
  mutate(n_genes_regulon = map_int(TF, \(tf) length(regulon_relapse[[tf]]$tfmode)),
         actividad = if_else(NES > 0, "Activo en recaída", "Reprimido en recaída"))

write_tsv(top10, file.path(out_dir, "top10_reguladores_maestros.tsv"))
cat(sprintf("Universo: %d genes\n", length(universo)))
print(top10, row.names = FALSE)

# ------------------------------ STEP 2. ORA -----------------------------------
ora <- function(tf, coleccion) {
  genes <- names(regulon_relapse[[tf]]$tfmode)
  res <- enricher(gene = genes, universe = universo,
                  TERM2GENE = colecciones[[coleccion]],
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                  minGSSize = 10, maxGSSize = 5000)
  if (is.null(res)) return(NULL)
  as_tibble(res@result) %>%
    mutate(TF = tf, coleccion = coleccion, .before = 1)
}

resultados <- expand_grid(TF = top10$TF, coleccion = names(colecciones)) %>%
  pmap_dfr(\(TF, coleccion) ora(TF, coleccion)) %>%
  left_join(select(top10, TF, NES_TF = NES, actividad), by = "TF") %>%
  arrange(coleccion, TF, p.adjust)

write_tsv(resultados, file.path(out_dir, "ora_regulones_top10.tsv"))

resultados %>%
  filter(p.adjust < 0.05) %>%
  count(coleccion, TF, name = "conjuntos_significativos") %>%
  print(n = Inf)

# ------------------------------ STEP 3. PLOT ----------------------------------
siglas <- "\\b(MYC|E2F|G2M|MTORC1|IL2|STAT5|IL6|JAK|STAT3|KRAS|UV|TNFA|NFKB|PI3K|AKT|MTOR|P53|TGF|WNT|DNA|IFN)\\b"
nombre_bonito <- function(x) {
  x |>
    str_remove("^HALLMARK_") |>
    str_replace_all("_", " ") |>
    str_to_sentence() |>
    str_replace_all(regex(siglas, ignore_case = TRUE), toupper)
}

# Regulones ordenados por NES; etiqueta con la dirección de actividad
orden_tf <- top10 %>%
  arrange(desc(NES)) %>%
  mutate(etiqueta = sprintf("%s\n(NES %+.1f)", TF, NES))

# Escalas comunes a ambos paneles, para que compartan una sola leyenda
sig <- filter(resultados, p.adjust < 0.05)
lim_color <- range(-log10(sig$p.adjust))
lim_size  <- range(sig$Count)

grafica_ora <- function(col, titulo, eje_x = TRUE) {
  df <- resultados %>%
    filter(coleccion == col) %>%
    mutate(GeneRatio_num = map_dbl(GeneRatio, \(r) eval(parse(text = r))),
           termino = nombre_bonito(ID))
  sig_terms <- unique(df$termino[df$p.adjust < 0.05])
  df <- df %>%
    filter(termino %in% sig_terms, p.adjust < 0.05) %>%
    mutate(TF = factor(TF, levels = orden_tf$TF, labels = orden_tf$etiqueta),
           termino = fct_reorder(termino, -log10(p.adjust), .fun = sum))
  if (nrow(df) == 0) {
    return(ggplot() + labs(title = titulo, subtitle = "Sin conjuntos con FDR < 0.05") +
             theme_void() + theme(plot.title = element_text(face = "bold")))
  }
  ggplot(df, aes(TF, termino)) +
    geom_point(aes(size = Count, colour = -log10(p.adjust))) +
    scale_colour_gradient(low = "#9ECAE1", high = "#08306B", limits = lim_color,
                          name = expression(-log[10](FDR))) +
    scale_size_continuous(name = "Genes del regulón\nen el conjunto", range = c(2, 8),
                          limits = lim_size) +
    scale_x_discrete(drop = FALSE) +
    labs(title = titulo, x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 12),
          axis.text.x = if (eje_x) element_text(angle = 45, hjust = 1, colour = "black")
                        else element_blank(),
          axis.text.y = element_text(colour = "black"),
          panel.grid.minor = element_blank())
}

n_terms <- resultados %>%
  filter(p.adjust < 0.05) %>%
  distinct(coleccion, ID) %>%
  count(coleccion) %>%
  deframe()

figura <- (grafica_ora("MSigDB Hallmarks", "A. MSigDB Hallmarks", eje_x = FALSE) /
           grafica_ora("Cancer Hallmarks (Menyhart, core)",
                       "B. Cancer Hallmarks (Menyhart et al., 2025; core)", eje_x = FALSE) /
           grafica_ora("Cancer Hallmarks (Menyhart, integrated)",
                       "C. Cancer Hallmarks (Menyhart et al., 2025; integrated)")) +
  plot_layout(heights = pmax(n_terms[names(colecciones)], 2, na.rm = TRUE),
              guides = "collect") +
  plot_annotation(
    title = "ORA de los regulones de los 10 reguladores maestros principales",
    caption = paste0(
      "Prueba hipergeométrica (clusterProfiler::enricher); universo: genes de la matriz VST usada para construir la red.\n",
      "Solo se muestran conjuntos con FDR < 0.05 (Benjamini-Hochberg).\n",
      "NES: actividad del regulador en el MRA (> 0: activo en recaída; < 0: reprimido en recaída)."),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.caption = element_text(size = 8, colour = "#555555", hjust = 0)))

alto <- 3.5 + 0.35 * sum(pmax(n_terms, 2), na.rm = TRUE)
ggsave(file.path(out_dir, "ORA_regulones_top10.png"), figura,
       width = 10, height = alto, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "ORA_regulones_top10.pdf"), figura,
       width = 10, height = alto, device = cairo_pdf, limitsize = FALSE)
