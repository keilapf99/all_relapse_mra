# Title : Reguladores maestros de los genes leading edge de las vías enriquecidas
# Author: Keila Adonai Pérez-Fuentes
# Objective: Para cada vía enriquecida en el GSEA de la firma de expresión
# diferencial (MSigDB Hallmarks, paso 05a; Cancer Hallmarks de Menyhart core,
# paso 05b), responder:
#   1. ¿Algún regulón contiene todos los genes leading edge (core enrichment)
#      de la vía? -> cobertura = genes del leading edge en el regulón / total
#   2. ¿Qué reguladores maestros controlan esos genes? -> prueba hipergeométrica
#      de sobrerrepresentación del leading edge en cada regulón.
#
# - Vías: FDR < 0.05 en el GSEA.
# - Reguladores: los reguladores maestros con FDR < 0.05 en el MRA (paso 06),
#   con su regulón en la red de pacientes con recaída (todos sus genes blanco).
# - Universo: genes de la matriz VST con la que se construyó la red; el leading
#   edge se restringe a esos genes.
# - Benjamini-Hochberg global sobre todas las pruebas (vías x reguladores).
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/), p. ej.:
#   Rscript enrichment_analysis/09_leading_edge_reguladores.R
# Entradas: corto_mra/mra_relapse.RData (paso 06), preprocess/matrix_relapse_vst.tsv (paso 03),
#           enrichment_analysis/gsea_msigdb_hallmarks/gsea_msigdb_hallmarks.tsv (paso 05a),
#           enrichment_analysis/gsea_cancer_hallmarks_menyhart/gsea_cancer_hallmarks_core.tsv (paso 05b)
# Salidas (en enrichment_analysis/leading_edge_reguladores/):
#   - leading_edge_vs_regulones.tsv       : todas las combinaciones vía x regulador
#   - cobertura_por_via.tsv               : regulón que más cubre el leading edge de cada vía
#   - reguladores_significativos_por_via.tsv : pares vía - regulador con FDR < 0.05

suppressPackageStartupMessages({
  library(tidyverse)
})

out_dir <- "enrichment_analysis/leading_edge_reguladores"
dir.create(out_dir, showWarnings = FALSE)

# -------------------------------- DATA ----------------------------------------
load("corto_mra/mra_relapse.RData")  # alls, regulon_relapse, ...

universo <- rownames(read.table("preprocess/matrix_relapse_vst.tsv", header = TRUE,
                                check.names = FALSE))

reguladores <- alls %>% filter(FDR < 0.05)
regulones <- map(set_names(reguladores$TF), \(tf) names(regulon_relapse[[tf]]$tfmode))
# Todos los regulones de la red (no solo los maestros), para la cobertura máxima
regulones_red <- map(regulon_relapse, \(r) names(r$tfmode))

vias <- bind_rows(
  read_tsv("enrichment_analysis/gsea_msigdb_hallmarks/gsea_msigdb_hallmarks.tsv",
           show_col_types = FALSE) %>% mutate(coleccion = "MSigDB Hallmarks"),
  read_tsv("enrichment_analysis/gsea_cancer_hallmarks_menyhart/gsea_cancer_hallmarks_core.tsv",
           show_col_types = FALSE) %>% mutate(coleccion = "Cancer Hallmarks (Menyhart, core)")
) %>%
  filter(p.adjust < 0.05) %>%
  transmute(coleccion, via = ID, NES_via = NES, FDR_via = p.adjust,
            leading_edge = map(strsplit(core_enrichment, "/"), \(g) intersect(g, universo)),
            n_leading_edge = lengths(leading_edge))

cat(sprintf("%d vías enriquecidas | %d reguladores maestros | universo: %d genes\n",
            nrow(vias), nrow(reguladores), length(universo)))

# ------------------- STEP 1. Leading edge vs. cada regulón --------------------
N <- length(universo)

resultados <- vias %>%
  select(coleccion, via, NES_via, FDR_via, leading_edge, n_leading_edge) %>%
  expand_grid(TF = reguladores$TF) %>%
  mutate(
    n_genes_regulon = lengths(regulones[TF]),
    genes_compartidos = map2(leading_edge, TF, \(le, tf) sort(intersect(le, regulones[[tf]]))),
    k = lengths(genes_compartidos),
    cobertura = k / n_leading_edge,
    # P(X >= k), X ~ hipergeométrica: leading edge (n) tomado de N genes, de
    # los cuales n_genes_regulon pertenecen al regulón
    pvalue = phyper(k - 1, n_genes_regulon, N - n_genes_regulon, n_leading_edge,
                    lower.tail = FALSE)
  ) %>%
  mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
  left_join(select(reguladores, TF, NES_TF = NES, FDR_TF = FDR), by = "TF") %>%
  mutate(actividad_TF = if_else(NES_TF > 0, "Activo en recaída", "Reprimido en recaída"),
         genes_compartidos = map_chr(genes_compartidos, paste, collapse = ",")) %>%
  select(coleccion, via, NES_via, FDR_via, n_leading_edge, TF, NES_TF, actividad_TF,
         n_genes_regulon, k, cobertura, pvalue, FDR, genes_compartidos) %>%
  arrange(FDR)

write_tsv(resultados, file.path(out_dir, "leading_edge_vs_regulones.tsv"))

# ------------------- STEP 2. Cobertura máxima por vía -------------------------
cobertura <- vias %>%
  mutate(
    cob_union_maestros = map_dbl(leading_edge, \(le) mean(le %in% unlist(regulones))),
    cob_red = map(leading_edge, \(le) map_dbl(regulones_red, \(r) mean(le %in% r))),
    TF_max_red = map_chr(cob_red, \(x) names(which.max(x))),
    cobertura_max_red = map_dbl(cob_red, max)
  ) %>%
  select(coleccion, via, NES_via, n_leading_edge, cob_union_maestros,
         TF_max_red, cobertura_max_red) %>%
  left_join(resultados %>%
              group_by(via) %>%
              slice_max(cobertura, n = 1, with_ties = FALSE) %>%
              ungroup() %>%
              select(via, TF_max_maestro = TF, cobertura_max_maestro = cobertura),
            by = "via") %>%
  relocate(TF_max_maestro, cobertura_max_maestro, .after = n_leading_edge) %>%
  arrange(coleccion, desc(NES_via))

write_tsv(cobertura, file.path(out_dir, "cobertura_por_via.tsv"))

cat("\n1. Cobertura del leading edge por el mejor regulón de cada vía:\n")
cobertura %>%
  transmute(via = str_remove(via, "HALLMARK_"), n_LE = n_leading_edge,
            maestro = TF_max_maestro, cob = round(cobertura_max_maestro, 2),
            union_157 = round(cob_union_maestros, 2)) %>%
  as.data.frame() %>%
  print(right = FALSE)
cat(sprintf("Cobertura máxima de un solo regulón: %.0f%% | unión de los reguladores maestros (mediana): %.0f%%\n",
            100 * max(cobertura$cobertura_max_red), 100 * median(cobertura$cob_union_maestros)))

# ------------------- STEP 3. Reguladores significativos -----------------------
significativos <- resultados %>% filter(FDR < 0.05)
write_tsv(significativos, file.path(out_dir, "reguladores_significativos_por_via.tsv"))

cat(sprintf("\n2. Pares vía - regulador maestro con FDR < 0.05: %d\n", nrow(significativos)))
significativos %>%
  transmute(via = str_remove(via, "HALLMARK_"), TF, NES_TF = round(NES_TF, 1),
            genes = sprintf("%d/%d", k, n_leading_edge), FDR = signif(FDR, 2),
            genes_compartidos) %>%
  as.data.frame() %>%
  print(right = FALSE)
