# Title : Genera todas las figuras con el formato de la tesis
# Author: Keila Adonai Pérez-Fuentes
#
# Uso (desde la carpeta raíz del proyecto, leukemia_maestria/), después de
# correr la tubería (00_run_pipeline.R):
#   Rscript figuras_formato_tesis/00_run_figuras.R          # todas
#   Rscript figuras_formato_tesis/00_run_figuras.R F3 F6    # solo algunas
#
# Cada script usa los resultados ya guardados por la tubería; no recalcula
# ningún análisis. La salida de cada uno queda en logs/<script>.log

figuras <- c(
  "F1" = "figuras_formato_tesis/F1_estadistica_poblacion.R",
  "F2" = "figuras_formato_tesis/F2_filtrado_pca_vst.R",
  "F3" = "figuras_formato_tesis/F3_volcano_dea.R",
  "F4" = "figuras_formato_tesis/F4_mra_corto.R",
  "F5" = "figuras_formato_tesis/F5_gsea_msigdb_hallmarks.R",
  "F6" = "figuras_formato_tesis/F6_ora_regulones_hallmarks.R",
  "F7" = "figuras_formato_tesis/F7_gsea_lollipop_top10_top5.R"
)

if (!file.exists("raw_data/raw_matrix.tsv")) {
  stop("Ejecuta este script desde la carpeta raíz del proyecto (leukemia_maestria/)")
}

seleccion <- commandArgs(trailingOnly = TRUE)
if (length(seleccion) > 0) figuras <- figuras[seleccion]
dir.create("logs", showWarnings = FALSE)

for (fig in names(figuras)) {
  log <- file.path("logs", sprintf("%s.log", tools::file_path_sans_ext(basename(figuras[[fig]]))))
  cat(sprintf("[%s] Figura %s: %s ... ", format(Sys.time(), "%H:%M:%S"), fig, figuras[[fig]]))
  estado <- system2("Rscript", figuras[[fig]], stdout = log, stderr = log)
  if (estado != 0) {
    cat("ERROR\n")
    stop(sprintf("La figura %s falló; revisa %s", fig, log))
  }
  cat("OK\n")
}

cat("Figuras completas.\n")
