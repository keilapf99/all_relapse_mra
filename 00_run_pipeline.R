# Title : Ejecuta toda la tubería de análisis en orden
# Author: Keila Adonai Pérez-Fuentes
#
# Uso (desde la carpeta raíz del proyecto, leukemia_maestria/):
#   Rscript 00_run_pipeline.R          # corre todos los pasos (01 a 10)
#   Rscript 00_run_pipeline.R 04 05a   # corre solo los pasos indicados
#
# Cada paso se ejecuta en un proceso de R independiente, de modo que un paso
# solo puede usar lo que los anteriores guardaron en disco (no objetos que
# hayan quedado en memoria). Si un paso falla, la tubería se detiene.
# La salida de cada paso queda en logs/<paso>.log

pasos <- c(
  "01" = "preprocess/01_first_explo_pca.R",
  "02" = "statistics_metadata/02_comparacion_variables_fisher.R",
  "03" = "preprocess/03_preprocessing_vst_pca.R",
  "04" = "DEA/04_DEA_relapse.R",
  "05a" = "enrichment_analysis/05a_gsea_msigdb_hallmarks_relapse.R",
  "05b" = "enrichment_analysis/05b_gsea_cancer_hallmarks_menyhart_relapse.R",
  "06" = "corto_mra/06_mra_relapse.R",
  "07" = "enrichment_analysis/07_ora_regulones_top10.R",
  "08" = "enrichment_analysis/08_red_regulones_hallmarks.R",
  "09" = "enrichment_analysis/09_leading_edge_reguladores.R",
  "10" = "DEA/10_heatmap_agrupamiento.R"
)

if (!file.exists("raw_data/raw_matrix.tsv")) {
  stop("Ejecuta este script desde la carpeta raíz del proyecto (leukemia_maestria/)")
}

seleccion <- commandArgs(trailingOnly = TRUE)
if (length(seleccion) > 0) pasos <- pasos[seleccion]
dir.create("logs", showWarnings = FALSE)

for (paso in names(pasos)) {
  log <- file.path("logs", sprintf("%s.log", tools::file_path_sans_ext(basename(pasos[[paso]]))))
  cat(sprintf("[%s] Paso %s: %s ... ", format(Sys.time(), "%H:%M:%S"), paso, pasos[[paso]]))
  estado <- system2("Rscript", pasos[[paso]], stdout = log, stderr = log)
  if (estado != 0) {
    cat("ERROR\n")
    stop(sprintf("El paso %s falló; revisa %s", paso, log))
  }
  cat("OK\n")
}

cat("Tubería completa.\n")
