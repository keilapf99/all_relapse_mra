# Title : Master Regulator Analysis - CORTO
# Author: Keila Adonai Pérez-Fuentes
# Date: January 20th 2026
# Objective: 
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/), p. ej.:
#   Rscript corto_mra/06_mra_relapse.R
# Entradas: preprocess/matrix_relapse_vst.tsv, preprocess/matrix_control_vst.tsv (paso 03),
#           raw_data/tfs_motif_yestf.txt, DEA/DEA_totales_relapse.tsv (paso 04)
# Salidas : corto_mra/mra_relapse.RData, corto_mra/cortoMRS.pdf,
#           corto_mra/mra_relapse_todos_TFs.tsv, corto_mra/mra_relapse_TFs_FDR005.tsv

# Packages ---------------------------------------------------------------------
library(corto)

# corto (bootstraps) y el MRA (permutaciones) son aleatorios: semilla fija
set.seed(123)

# Data -------------------------------------------------------------------------
counts_relapse <- read.table(file = "preprocess/matrix_relapse_vst.tsv",
                     header = T,
                     check.names = F)
counts_control <- read.table(file = "preprocess/matrix_control_vst.tsv",
                             header = T,
                             check.names = F)
tfs<- readLines("raw_data/tfs_motif_yestf.txt")
degs_total <- read.delim(file = "DEA/DEA_totales_relapse.tsv", quote = "")

# STEP 1. REGULON --------------------------------------------------------------
centroids_relapse <- tfs[tfs %in% rownames(counts_relapse)]
regulon_relapse <- corto(counts_relapse, 
                 centroids=centroids_relapse, 
                 nbootstraps=100, 
                 p=1e-5, 
                 nthreads=6)

centroids_control <- tfs[tfs %in% rownames(counts_control)]
regulon_control <- corto(counts_control, 
                         centroids=centroids_control, 
                         nbootstraps=100, 
                         p=1e-5, 
                         nthreads=6)

# STEP 2. MOLECULAR SIGNATURE --------------------------------------------------
# This ms was previously calculated by deseq dea
molecular_signature <- degs_total$stat
names(molecular_signature) <- degs_total$gene_symbol
head(molecular_signature)

# Delete NAs
molecular_signature <- molecular_signature[!is.na(molecular_signature) & is.finite(molecular_signature)]

# Filtrar por genes presentes en tu regulón
genes_regulon <- unique(unlist(sapply(regulon_relapse, function(x) names(x$tfmode))))
molecular_signature <- molecular_signature[names(molecular_signature) %in% genes_regulon]

# STEP 2. MASTER REGULATOR ANALYSIS --------------------------------------------
source("corto_mra/mra_mod.R")
mrs<- mra_mod(molecular_signature, 
           regulon = regulon_relapse, 
           minsize = 20, 
           nperm = 1000, 
           nthreads = 6, 
           verbose = TRUE, 
           atacseq = NULL)

# STEP 3. PLOTS ----------------------------------------------------------------
pdf(file = "corto_mra/cortoMRS.pdf", width = 15, height = 10)
mraplot(mrs, mrs = 5,10)
dev.off()

alls <- data.frame(TF= names(mrs$nes),
                   NES = mrs$nes,
                   pvalue = mrs$pvalue)
alls$FDR <-  p.adjust(alls$pvalue, method = "fdr", n = length(alls$pvalue))

tops <- alls[alls$FDR < 0.05, ]

# STEP 4. SAVE FILES -----------------------------------------------------------
save(alls, mrs, tops, regulon_relapse, regulon_control, file = "corto_mra/mra_relapse.RData")
write.table(alls[order(alls$FDR), ], file = "corto_mra/mra_relapse_todos_TFs.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tops[order(tops$FDR), ], file = "corto_mra/mra_relapse_TFs_FDR005.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
