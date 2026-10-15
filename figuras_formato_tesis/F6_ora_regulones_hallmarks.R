# Title : ORA de los regulones del top 10 de reguladores maestros contra los
#         Hallmarks de MSigDB (formato de figuras de la tesis)
# Author: Keila Adonai Pérez-Fuentes
# Basado en: all_relapse_masters_degree/4_GSEA_ORA/ora_top10_master_regulators_hallmarks.r
#            y ora_top10_master_regulators_publication_plots.r
# Objective: Mismas figuras que el original, con los resultados del ORA del
# paso 07 de leukemia_maestria (colección MSigDB Hallmarks). No recalcula el ORA:
# convierte la tabla de clusterProfiler::enricher al formato de columnas que usaba
# all_relapse y llama al mismo código de las figuras.
#
# Diferencia de método con all_relapse (por eso hay menos asociaciones):
#   - all_relapse: universo = todos los símbolos de la matriz cruda (39,301 genes).
#   - leukemia_maestria (paso 07): universo = genes de la matriz VST con la que
#     se construyó la red (14,208), que son los únicos que podían estar en un regulón.
#   En ambos la p se ajusta por BH dentro de cada regulón y se usa FDR < 0.05.
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/):
#   Rscript figuras_formato_tesis/F6_ora_regulones_hallmarks.R
# Entradas: enrichment_analysis/ora_regulones_top10/ora_regulones_top10.tsv (paso 07),
#           enrichment_analysis/ora_regulones_top10/top10_reguladores_maestros.tsv (paso 07),
#           preprocess/matrix_relapse_vst.tsv (paso 03, universo)
# Salidas : figuras_formato_tesis/F6_ora_regulones_hallmarks/
#   - ora_top10_master_regulators.csv
#   - ora_regulon_universe_genes.txt
#   - ora_hallmarks_top10_master_regulators_all_results.csv
#   - ora_hallmarks_top10_master_regulators_significant.csv
#   - ora_hallmarks_pathways_by_regulon{,_3col}.{png,pdf}
#   - ora_hallmarks_pathways_by_regulon_canva_slide.png

# -------------------------------- PACKAGES ------------------------------------
library(dplyr)
library(stringr)

source("figuras_formato_tesis/fuente_tesis.R")

out_dir <- "figuras_formato_tesis/F6_ora_regulones_hallmarks"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------- DATA ----------------------------------------
ora_l <- read.delim("enrichment_analysis/ora_regulones_top10/ora_regulones_top10.tsv",
                    quote = "") %>%
  filter(coleccion == "MSigDB Hallmarks")

top_master_regulators <- read.delim(
  "enrichment_analysis/ora_regulones_top10/top10_reguladores_maestros.tsv"
) %>%
  mutate(
    abs_NES   = abs(NES),
    direction = if_else(NES >= 0, "NES positivo", "NES negativo")
  ) %>%
  arrange(FDR, desc(abs_NES))

universe_genes <- rownames(read.table("preprocess/matrix_relapse_vst.tsv",
                                      header = TRUE, check.names = FALSE))

# ── Hallmark label formatter (igual que en all_relapse) ───────────────────────
format_hallmark <- function(x) {
  x %>%
    str_remove("^HALLMARK_") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    # Protein / pathway acronyms ─ order matters (longest first)
    str_replace_all("\\bMtorc1\\b",  "mTORC1") %>%
    str_replace_all("\\bMtor\\b",    "mTOR")   %>%
    str_replace_all("\\bStat3\\b",   "STAT3")  %>%
    str_replace_all("\\bStat5\\b",   "STAT5")  %>%
    str_replace_all("\\bIl2\\b",     "IL2")    %>%
    str_replace_all("\\bIl6\\b",     "IL6")    %>%
    str_replace_all("\\bJak\\b",     "JAK")    %>%
    str_replace_all("\\bKras\\b",    "KRAS")   %>%
    str_replace_all("\\bMyc\\b",     "MYC")    %>%
    str_replace_all("\\bPi3k\\b",    "PI3K")   %>%
    str_replace_all("\\bAkt\\b",     "AKT")    %>%
    str_replace_all("\\bP53\\b",     "p53")    %>%
    str_replace_all("\\bDna\\b",     "DNA")    %>%
    str_replace_all("\\bE2f\\b",     "E2F")    %>%
    str_replace_all("\\bG2m\\b",     "G2M")    %>%
    str_replace_all("\\bUv\\b",      "UV")
}

# ── Tabla del paso 07 -> columnas de all_relapse ──────────────────────────────
# Las proporciones se calculan como en all_relapse: gene_ratio = overlap / genes
# del regulón y background_ratio = genes del hallmark / universo. (En la tabla de
# enricher, GeneRatio y BgRatio se dividen solo entre los genes con alguna
# anotación Hallmark, por eso no se usan directamente.)
ora_results <- ora_l %>%
  left_join(top_master_regulators %>% select(TF, n_genes_regulon), by = "TF") %>%
  transmute(
    TF,
    hallmark         = ID,
    hallmark_label   = format_hallmark(ID),
    regulon_size     = n_genes_regulon,
    hallmark_size    = as.integer(sub("/.*", "", BgRatio)),
    overlap          = Count,
    gene_ratio       = overlap / regulon_size,
    background_ratio = hallmark_size / length(universe_genes),
    fold_enrichment  = gene_ratio / background_ratio,
    odds_ratio       = oddsRatio,
    pvalue,
    overlap_genes    = str_replace_all(geneID, "/", ";"),
    p.adjust,
    significant      = p.adjust < 0.05
  ) %>%
  mutate(
    # Global BH adjustment across all TF × hallmark pairs (reference only).
    p.adjust.global = p.adjust(pvalue, method = "BH"),
    minus_log10_fdr = -log10(pmax(p.adjust, .Machine$double.xmin))
  ) %>%
  left_join(
    top_master_regulators %>% transmute(TF, MRA_NES = NES, MRA_FDR = FDR),
    by = "TF"
  )

ora_significant <- ora_results %>%
  filter(significant) %>%
  arrange(TF, p.adjust, desc(fold_enrichment))

# ── Write outputs ─────────────────────────────────────────────────────────────
write.csv(
  top_master_regulators %>% select(TF, NES, pvalue, FDR, abs_NES, direction),
  file.path(out_dir, "ora_top10_master_regulators.csv"),
  row.names = FALSE
)

writeLines(
  sort(universe_genes),
  file.path(out_dir, "ora_regulon_universe_genes.txt")
)

write.csv(
  ora_results,
  file.path(out_dir, "ora_hallmarks_top10_master_regulators_all_results.csv"),
  row.names = FALSE
)

write.csv(
  ora_significant,
  file.path(out_dir, "ora_hallmarks_top10_master_regulators_significant.csv"),
  row.names = FALSE
)

# ── Figuras (mismo código que all_relapse) ────────────────────────────────────
source("figuras_formato_tesis/ora_publication_plots_tesis.R")

cat("Universe genes      :", length(universe_genes), "\n")
cat("Top master regulators:", paste(top_master_regulators$TF, collapse = ", "), "\n")
cat("Significant pairs   :", nrow(ora_significant), "\n")
