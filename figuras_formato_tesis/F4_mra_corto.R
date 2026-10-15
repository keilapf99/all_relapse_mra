# Title : Master Regulator Analysis - CORTO (formato de figuras de la tesis)
# Author: Keila Adonai Pérez-Fuentes
# Basado en: all_relapse_masters_degree/3_MRA_corto/mra_corto_adapted.r (STEP 3)
#            y mra_tfs_interes.r
# Objective: Mismos mraplot que el original, con el objeto mrs del paso 06 de
# leukemia_maestria. No recalcula el MRA.
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/):
#   Rscript figuras_formato_tesis/F4_mra_corto.R
# Entradas: corto_mra/mra_relapse.RData (paso 06)
# Salidas : figuras_formato_tesis/F4_mra_corto/
#   - cortoMRStop1.pdf    : reguladores maestros 1 a 5
#   - cortoMRStop2.pdf    : reguladores maestros 6 a 10
#   - mra_tfs_interes.pdf : TFs de interés biológico

# ----------------------------- PACKAGES ----------------------------------------
library(corto)

out_dir <- "figuras_formato_tesis/F4_mra_corto"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------- DATA -----------------------------------------
load("corto_mra/mra_relapse.RData")  # mrs, alls, tops, regulon_relapse, ...

# STEP 3. PLOTS ----------------------------------------------------------------
pdf(file = file.path(out_dir, "cortoMRStop1.pdf"), width = 14, height = 12)
mraplot(mrs, 
        mrs = 5,
        title = "Reguladores maestros en la recaída de la LLA pediátrica")
dev.off()
ordered_tfs <- names(sort(mrs$pvalue))
top2 <- ordered_tfs[6:10]
pdf(file = file.path(out_dir, "cortoMRStop2.pdf"), width = 14, height = 12)
mraplot(
  mrs,
  mrs = top2,
  title = "Reguladores maestros en la recaída de la LLA pediátrica"
)
dev.off()

# ----------------------- TFs DE INTERES ---------------------------------------
tfs_interes <- c("TCF3", "PBX1", "ETV6", "RUNX1", "EBF1",
                 "PAX5", "EA2", "IKZF1", "IKZF2", "IKZF3", "HLF")

# Verificar cuales estan presentes en el objeto mrs
tfs_encontrados <- tfs_interes[tfs_interes %in% names(mrs$nes)]
tfs_faltantes   <- tfs_interes[!tfs_interes %in% names(mrs$nes)]

message("TFs encontrados en mrs: ", paste(tfs_encontrados, collapse = ", "))
if (length(tfs_faltantes) > 0)
  message("TFs NO encontrados en mrs: ", paste(tfs_faltantes, collapse = ", "))

pdf(file = file.path(out_dir, "mra_tfs_interes.pdf"),
    width = 14, height = 12)
mraplot(mrs,
        mrs   = tfs_encontrados,
        title = "Reguladores maestros de interes en la recaida de la LLA pediatrica")
dev.off()

message("Plot guardado: mra_tfs_interes.pdf")
