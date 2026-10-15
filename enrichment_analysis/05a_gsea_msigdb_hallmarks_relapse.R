# Title : Gene Set Enrichment Analysis to LLA differential expression genes
# Author: Keila Adonai Pérez-Fuentes
# Date: May 29th 2026
# Objective: Enrich the differential expression genes in ALL relapse in cancer
# hallmarks (colección Hallmark (H) de MSigDB).
#
# Basado en gsea_hallmarks.R. GSEA con clusterProfiler; la firma es el
# estadístico de Wald de DESeq2 (columna `stat`, contraste relapse si vs. no),
# ordenado de mayor a menor. NES > 0 = enriquecido en pacientes CON recaída.
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/), p. ej.:
#   Rscript enrichment_analysis/05a_gsea_msigdb_hallmarks_relapse.R
# Entradas: DEA/DEA_totales_relapse.tsv (paso 04)
# Salidas (en enrichment_analysis/gsea_msigdb_hallmarks/):
#   - gene_sets/msigdb_hallmarks_H.tsv   : copia local de la colección H
#   - gsea_msigdb_hallmarks.tsv          : resultados de los 50 hallmarks
#   - GSEA_msigdb_hallmarks_lollipop.png/.pdf : hallmarks significativos (top 10 por |NES|)

# ------------------------------ PACKAGES --------------------------------------
# libglpk (requerida por igraph -> enrichplot -> clusterProfiler) no está en el
# sistema de este servidor; se carga desde un entorno conda del usuario si existe.
glpk <- "/scratch/home/aperez/.conda_envs/glpk/lib/libglpk.so.40"
if (file.exists(glpk)) dyn.load(glpk, local = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(msigdbr)
})

out_dir <- "enrichment_analysis/gsea_msigdb_hallmarks"
dir.create(out_dir, showWarnings = FALSE)

# -------------------------------- DATA ----------------------------------------
degs_totales <- read.delim("DEA/DEA_totales_relapse.tsv", quote = "")

# Copia local de la colección H: se descarga una sola vez para que el análisis
# no cambie si MSigDB publica una versión nueva
archivo_h <- file.path(out_dir, "gene_sets", "msigdb_hallmarks_H.tsv")
if (!file.exists(archivo_h)) {
  dir.create(dirname(archivo_h), showWarnings = FALSE)
  msigdbr(species = "Homo sapiens", collection = "H") %>%
    dplyr::select(gs_name, gene_symbol, db_version) %>%
    distinct() %>%
    write_tsv(archivo_h)
}
hallmarks <- read_tsv(archivo_h, show_col_types = FALSE)
version_msigdb <- unique(hallmarks$db_version)
hallmarks <- hallmarks %>% dplyr::select(gs_name, gene_symbol)

# ---------------------- STEP 1. Vector needed to GSEA -------------------------
# Nombre del gen -> stat, ordenado de mayor a menor. Si un símbolo aparece
# repetido (varios Ensembl ID), se conserva el de mayor |stat|.
gene_list <- degs_totales %>%
  filter(!is.na(stat)) %>%
  group_by(gene_symbol) %>%
  slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(stat)) %>%
  (\(d) setNames(d$stat, d$gene_symbol))()

cat(sprintf("Firma: %d genes únicos | %d hallmarks de MSigDB (versión %s)\n",
            length(gene_list), n_distinct(hallmarks$gs_name), version_msigdb))

# 2. Run GSEA analysis ---------------------------------------------------------
# pvalueCutoff = 1 para guardar los 50 hallmarks; la significancia se define
# después con FDR (p.adjust, Benjamini-Hochberg) < 0.05.
set.seed(123)
gsea_hallmark <- GSEA(geneList = gene_list,
                      TERM2GENE = hallmarks,
                      pvalueCutoff = 1,
                      pAdjustMethod = "BH",
                      minGSSize = 10,
                      maxGSSize = 500,
                      eps = 0,
                      seed = TRUE,
                      verbose = FALSE)

gsea_res <- as_tibble(gsea_hallmark@result) %>%
  arrange(p.adjust)

write_tsv(gsea_res, file.path(out_dir, "gsea_msigdb_hallmarks.tsv"))

cat(sprintf("%d/%d hallmarks significativos (FDR < 0.05): %d en recaída, %d sin recaída\n",
            sum(gsea_res$p.adjust < 0.05), nrow(gsea_res),
            sum(gsea_res$p.adjust < 0.05 & gsea_res$NES > 0),
            sum(gsea_res$p.adjust < 0.05 & gsea_res$NES < 0)))
gsea_res %>%
  filter(p.adjust < 0.05) %>%
  select(ID, setSize, NES, p.adjust) %>%
  print(n = Inf)

# 3. Data plots ----------------------------------------------------------------
colores <- c("CON RECAÍDA" = "#B2182B", "SIN RECAÍDA" = "#2166AC")

# Siglas que deben quedar en mayúsculas en los nombres de los hallmarks
siglas <- "\\b(MYC|E2F|G2M|MTORC1|IL2|STAT5|IL6|JAK|STAT3|KRAS|UV|TNFA|NFKB|PI3K|AKT|MTOR|P53|TGF|WNT|DNA|IFN)\\b"

# Top 10 hallmarks significativos por |NES|
top10 <- gsea_res %>%
  filter(p.adjust < 0.05) %>%
  slice_max(abs(NES), n = 10) %>%
  mutate(
    pathway = Description |>
      str_remove("HALLMARK_") |>
      str_replace_all("_", " ") |>
      str_to_sentence() |>
      str_replace_all(regex(siglas, ignore_case = TRUE), toupper),
    pathway = fct_reorder(pathway, NES),
    direccion = factor(if_else(NES > 0, "CON RECAÍDA", "SIN RECAÍDA"),
                       levels = names(colores))
  )

lolli_top10 <- ggplot(top10, aes(x = NES, y = pathway)) +
  geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = NES, yend = pathway),
               linewidth = 1.6, colour = "grey80") +
  geom_point(aes(size = -log10(p.adjust), colour = direccion)) +
  geom_text(aes(label = sprintf("%.2f", NES),
                hjust = if_else(NES > 0, -0.6, 1.6)),
            size = 3.6, fontface = "bold") +
  scale_colour_manual(values = colores, name = "Enriquecido en") +
  scale_size_continuous(name = expression(-log[10](FDR)), range = c(4, 10)) +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  labs(title = "Hallmarks de MSigDB enriquecidos (FDR < 0.05)",
       subtitle = sprintf("Top %d por |NES| de %d significativos",
                          nrow(top10), sum(gsea_res$p.adjust < 0.05)),
       x = "Normalized Enrichment Score (NES)", y = NULL,
       caption = paste0("GSEA preordenado (clusterProfiler) con el estadístico de Wald de DESeq2 ",
                        "(recaída vs. sin recaída).\nNES > 0: enriquecido en pacientes con recaída; ",
                        "NES < 0: en pacientes sin recaída.")) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = "#555555"),
        plot.caption = element_text(colour = "#555555", hjust = 0, size = 8.5),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "GSEA_msigdb_hallmarks_lollipop.png"), lolli_top10,
       width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "GSEA_msigdb_hallmarks_lollipop.pdf"), lolli_top10,
       width = 10, height = 7, device = cairo_pdf)
