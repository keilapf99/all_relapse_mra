# GSEA de la firma de expresión diferencial completa (recaída vs. sin recaída)
# contra los "Cancer Hallmarks" de Menyhart et al. (NO los Hallmarks de MSigDB).
#
# Basado en gsea_hallmarks_menyhart.R (plantilla genérica). GSEA con
# clusterProfiler; las curvas de enriquecimiento se dibujan con fgsea.
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/), p. ej.:
#   Rscript enrichment_analysis/05b_gsea_cancer_hallmarks_menyhart_relapse.R
# Entradas: DEA/DEA_totales_relapse.tsv (paso 04)
#
# Conjuntos de genes: Menyhart O, Jayasekara Kothalawala W, Gyorffy B. "A gene set
# enrichment analysis for cancer hallmarks." J Pharm Anal. 2025;15(5):101065.
# doi:10.1016/j.jpha.2024.101065 (cancerhallmarks.com). Diez hallmarks de
# Hanahan & Weinberg, en dos versiones:
#   core       genes respaldados por >=2 de 7 bases de datos (consenso; principal)
#   integrated genes respaldados por >=1 base de datos (unión; complementario)
#
# Ranking: estadístico de Wald de DESeq2 (columna `stat`), contraste
# relapse "si" vs. "no". NES > 0 = hallmark enriquecido en pacientes CON recaída.
#
# Salidas (en gsea_cancer_hallmarks_menyhart/):
#   - gene_sets/Menyhart_JPA_CancerHallmarks_<set>.txt : copia local de los conjuntos
#   - gsea_cancer_hallmarks_<set>.tsv                  : resultados completos
#   - GSEA_cancer_hallmarks_NES.png / .pdf             : NES de los 10 hallmarks
#   - GSEA_cancer_hallmarks_core_curvas.png / .pdf     : curvas de enriquecimiento (core)

# libglpk (requerida por igraph -> enrichplot -> clusterProfiler) no está en el
# sistema de este servidor; se carga desde un entorno conda del usuario si existe.
glpk <- "/scratch/home/aperez/.conda_envs/glpk/lib/libglpk.so.40"
if (file.exists(glpk)) dyn.load(glpk, local = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(fgsea)  # plotEnrichment()
  library(patchwork)
})

base_dir <- "."  # raíz del proyecto: ejecutar desde leukemia_maestria/
out_dir  <- file.path(base_dir, "enrichment_analysis/gsea_cancer_hallmarks_menyhart")
gs_dir   <- file.path(out_dir, "gene_sets")
dir.create(gs_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Firma de expresión diferencial ----------------------------------------

degs_totales <- read.delim(file.path(base_dir, "DEA/DEA_totales_relapse.tsv"),
                           quote = "")

# Vector ordenado nombre del gen -> stat. Si un símbolo aparece repetido
# (varios Ensembl ID), se conserva el de mayor |stat|.
ranking <- degs_totales %>%
  filter(!is.na(stat), !is.na(gene_symbol), gene_symbol != "") %>%
  group_by(gene_symbol) %>%
  slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(stat)) %>%
  (\(d) setNames(d$stat, d$gene_symbol))()

cat(sprintf("Firma: %d genes únicos (de %d filas)\n", length(ranking), nrow(degs_totales)))

# ---- Conjuntos de genes de Menyhart -----------------------------------------

# Se descargan una sola vez y se guardan, para que el análisis sea reproducible.
# Devuelve una lista hallmark -> vector de genes (formato de fgsea).
get_cancer_hallmarks <- function(geneset = c("core", "integrated")) {
  geneset <- match.arg(geneset)
  archivo <- file.path(gs_dir, sprintf("Menyhart_JPA_CancerHallmarks_%s.txt", geneset))
  if (!file.exists(archivo)) {
    download.file(sprintf("https://cancerhallmarks.com/download_file/Menyhart_JPA_CancerHallmarks_%s.txt",
                          geneset), archivo, quiet = TRUE)
  }
  tok <- lapply(strsplit(readLines(archivo, warn = FALSE), "\t"), \(t) t[nzchar(t)])
  setNames(lapply(tok, \(t) unique(t[-1])), map_chr(tok, 1))
}

# ---- GSEA -------------------------------------------------------------------

# maxGSSize se sube porque los hallmarks son grandes (hasta ~3,500 genes en
# "integrated"); pvalueCutoff = 1 para guardar los 10 hallmarks.
correr_gsea <- function(geneset) {
  hallmarks <- get_cancer_hallmarks(geneset)
  set.seed(123)
  t2g <- tibble(gs_name = rep(names(hallmarks), lengths(hallmarks)),
                gene_symbol = unlist(hallmarks, use.names = FALSE))
  gs <- GSEA(geneList = ranking, TERM2GENE = t2g, pvalueCutoff = 1, eps = 0,
             seed = TRUE, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 5000,
             verbose = FALSE)
  res <- as_tibble(gs@result) %>%
    mutate(genes_en_set = lengths(hallmarks[ID]), conjunto = geneset) %>%
    arrange(p.adjust)
  write_tsv(res, file.path(out_dir, sprintf("gsea_cancer_hallmarks_%s.tsv", geneset)))
  cat(sprintf("[%s] %d/%d hallmarks significativos (FDR < 0.05)\n",
              geneset, sum(res$p.adjust < 0.05), nrow(res)))
  list(hallmarks = hallmarks, res = res)
}

gsea <- map(c(core = "core", integrated = "integrated"), correr_gsea)

gsea$core$res %>%
  select(ID, setSize, NES, pvalue, p.adjust) %>%
  print(n = Inf)

# ---- Gráfica de NES -----------------------------------------------------------

colores <- c("CON RECAÍDA" = "#B2182B", "SIN RECAÍDA" = "#2166AC")
gris_texto <- "#555555"

a_titulo <- function(x) str_to_sentence(x) |> str_replace("^(.)", toupper)

grafica_nes <- function(res, titulo, leyenda = TRUE) {
  df <- res %>%
    mutate(hallmark = fct_reorder(a_titulo(ID), NES),
           direccion = factor(if_else(NES > 0, "CON RECAÍDA", "SIN RECAÍDA"),
                              levels = names(colores)),
           significativo = p.adjust < 0.05,
           etiqueta = sprintf("FDR %s%s",
                              if_else(p.adjust < 0.001, "< 0.001",
                                      sprintf("= %.3f", p.adjust)),
                              if_else(significativo, " *", "")))
  lim <- max(abs(df$NES)) * 1.55
  ggplot(df, aes(NES, hallmark, fill = direccion, alpha = significativo)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.4) +
    geom_text(aes(label = etiqueta, x = NES + if_else(NES > 0, 0.05, -0.05),
                  hjust = if_else(NES > 0, 0, 1)),
              size = 3.1, colour = "#333333", alpha = 1) +
    scale_fill_manual(values = colores, name = "Enriquecido en", drop = FALSE) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    scale_x_continuous(limits = c(-lim, lim)) +
    labs(title = titulo, x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = if (leyenda) "bottom" else "none")
}

figura_nes <- (grafica_nes(gsea$core$res, "A. Cancer hallmarks: core (consenso)") /
               grafica_nes(gsea$integrated$res, "B. Cancer hallmarks: integrated",
                           leyenda = FALSE)) +
  plot_annotation(
    caption = paste0(
      "GSEA preordenado (clusterProfiler) con el estadístico de Wald de DESeq2 (recaída vs. sin recaída). ",
      "Conjuntos de genes: Menyhart et al., 2025.\n",
      "NES > 0: hallmark enriquecido en pacientes con recaída; NES < 0: en pacientes sin recaída.\n",
      "FDR: p ajustada por Benjamini-Hochberg; * FDR < 0.05 (barras sólidas)."),
    theme = theme(plot.caption = element_text(size = 8.5, colour = gris_texto, hjust = 0))
  )

ggsave(file.path(out_dir, "GSEA_cancer_hallmarks_NES.png"), figura_nes,
       width = 9, height = 10, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "GSEA_cancer_hallmarks_NES.pdf"), figura_nes,
       width = 9, height = 10, device = cairo_pdf)

# ---- Curvas de enriquecimiento (core) ---------------------------------------

# Hallmarks significativos; si no hay, los 4 de menor FDR
ids_curvas <- gsea$core$res %>%
  filter(p.adjust < 0.05) %>%
  pull(ID)
if (length(ids_curvas) == 0) ids_curvas <- head(gsea$core$res$ID, 4)

curva <- function(id) {
  r <- filter(gsea$core$res, ID == id)
  p <- plotEnrichment(gsea$core$hallmarks[[id]], ranking, ticksSize = 0.2)
  p$layers[[1]]$aes_params$colour <- colores[["CON RECAÍDA"]]  # curva de ES
  p +
    labs(title = a_titulo(id),
         subtitle = sprintf("NES = %.2f; FDR = %s", r$NES,
                            if_else(r$p.adjust < 0.001, "< 0.001", sprintf("%.3f", r$p.adjust))),
         x = "Posición del gen en el ranking",
         y = "Enrichment score") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(colour = gris_texto),
          panel.grid.minor = element_blank())
}

curvas <- wrap_plots(map(ids_curvas, curva), ncol = 2) +
  plot_annotation(title = "Cancer hallmarks (core): curvas de enriquecimiento",
                  caption = paste0(
                    "Genes ordenados por el estadístico de Wald de DESeq2: a la izquierda, ",
                    "los más sobreexpresados en recaída; a la derecha, en sin recaída.\n",
                    "Cada línea vertical es un gen del hallmark; la línea punteada marca el ",
                    "enrichment score máximo y mínimo."),
                  theme = theme(plot.title = element_text(face = "bold", size = 13),
                                plot.caption = element_text(size = 8.5, colour = gris_texto,
                                                            hjust = 0)))

ggsave(file.path(out_dir, "GSEA_cancer_hallmarks_core_curvas.png"), curvas,
       width = 10, height = 3.5 * ceiling(length(ids_curvas) / 2) + 0.5,
       dpi = 300, bg = "white")
ggsave(file.path(out_dir, "GSEA_cancer_hallmarks_core_curvas.pdf"), curvas,
       width = 10, height = 3.5 * ceiling(length(ids_curvas) / 2) + 0.5,
       device = cairo_pdf)
