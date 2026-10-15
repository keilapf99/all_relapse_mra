# Title : Preprocessing and normalization of RNA-seq counts for network analysis
# Author: Keila Adonai Pérez-Fuentes
# Date: May 23rd, 2026
# Objective:
# To preprocess and normalize RNA-seq count data from leukemia samples in order
# to generate high-quality expression matrices suitable for downstream network
# analysis. 
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/), p. ej.:
#   Rscript preprocess/03_preprocessing_vst_pca.R
# Entradas: raw_data/raw_matrix.tsv, preprocess/metadata_NS500560.tsv (paso 01)
# Salidas : preprocess/counts_vst.tsv, matrix_relapse_vst.tsv, matrix_control_vst.tsv,
#           metadata_relapse.tsv, metadata_control.tsv, PCA_vst_*.pdf

# -------------------------------- PACKAGES ------------------------------------
library(DESeq2)  # vst()
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(edgeR)
library(ggrepel)

# --------------------------------- DATA ---------------------------------------

counts <- read.table(
  file = "raw_data/raw_matrix.tsv",
  header = T,
  check.names = F)

metadata <- read.table(
  "preprocess/metadata_NS500560.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- counts[, colnames(counts) %in% metadata$sample, drop = FALSE]

# ------------------------ 1. COUNTS MATRIX PREPROCESS -------------------------

# 1.1 Filter low counts genes --------------------------------------------------
dge <- DGEList(counts = counts)
cpm_counts <- cpm(dge)
logCPM <- cpm(dge, log = TRUE)
pdf("preprocess/densidades_logCPM.pdf", width = 7, height = 6)
plotDensities(logCPM, legend = FALSE)
dev.off()


# Al menos 10 lecturas en al menos 15 muestras (mismo filtro que el DEA, paso 04)
counts <- counts[rowSums(counts >= 10) >= 15, ]
dim(counts)

# 1.2 Verification of input data ------------------------------------------------
# 1.2.1 Numeric matrix
counts <- as.matrix(counts)
class(counts)
unique(sapply(counts, class))

# 1.2.2 Verify identical samples order
rownames(metadata) <- metadata$sample
metadata <- metadata[colnames(counts), , drop = FALSE]
stopifnot(identical(colnames(counts), rownames(metadata)))
identical(colnames(counts), rownames(metadata))

############################# 2. Normalization #################################
counts_vst<- vst(counts)

# 2.1 Chunk genes/Duplicated genes ---------------------------------------------

# 2.1.1 Convert to gene symbol

gsymbols <- mapIds(org.Hs.eg.db, keys = rownames(counts_vst), 
                   keytype = "ENSEMBL", 
                   column = "SYMBOL")

counts_vst <-bind_cols(gsymbols, counts_vst)
dim(counts_vst)
# [1] 14953    37
counts_vst <- counts_vst[!is.na(counts_vst$...1),]
dim(counts_vst)
# [1] 12878    37

# 2.2 Media by sample ----------------------------------------------------------

# Define the function to obtain the media of gene expression by sample ----
# Function to calculate the media by column

colMean <- function(v) {
  apply(v, 2, mean)
}

counts <- do.call(rbind, lapply(split(counts_vst, counts_vst[,1]), 
                                function(chunk){colMean(chunk[,-1])}))

dim(counts)
# [1] 12891    35


################################# 3. PCA Plots #################################
counts_t <- t(counts)
pca <- prcomp(counts_t, center = TRUE, scale. = TRUE)
# Varianza explicada
pca_var <- summary(pca)$importance[2, 1:2] * 100

# 3.1 PCA + metadata -----------------------------------------------------------
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  left_join(metadata, by = "sample")

# 3.2 Format for facet --------------------------------------------------------
pca_long <- pca_df %>%
  dplyr::select(sample, PC1, PC2, instrument, relapse, flowcell) %>%
  pivot_longer(
    cols = c(instrument, relapse, flowcell),
    names_to = "variable",
    values_to = "group"
  )

# 3.3 Visualization ------------------------------------------------------------

pdf("preprocess/PCA_vst_instrument.pdf", width = 7, height = 6)

ggplot(pca_df, aes(PC1, PC2, color = instrument)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_brewer(palette = "Set3") +
  labs(
    title = "PCA (VST) – Instrument",
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)")
  ) +
  theme_bw()

dev.off()


pdf("preprocess/PCA_vst_flowcell.pdf", width = 7, height = 6)

ggplot(pca_df, aes(PC1, PC2, color = flowcell)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "PCA (VST) – Flowcell",
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)")
  ) +
  theme_bw()

dev.off()


pdf("preprocess/PCA_vst_relapse.pdf", width = 7, height = 6)

ggplot(pca_df, aes(PC1, PC2, color = relapse)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(
    values = c(
      "si" = "#D73027",   # rojo fuerte
      "no" = "#4575B4"    # azul fuerte
    )
  ) +
  labs(
    title = "PCA (VST) – Relapse",
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)")
  ) +
  theme_bw()

dev.off()
######################## 4. Classify by phenotype ##############################
metadata_norelapse <- metadata %>% 
  filter(relapse == "no")
metadata_relapse <- metadata %>% 
  filter(relapse == "si")
matrix_norelapse <- counts[, colnames(counts)%in%metadata_norelapse$sample]
matrix_relapse <- counts[, colnames(counts)%in%metadata_relapse$sample]

############################### 5. Save files ##################################

write.table(counts, file = "preprocess/counts_vst.tsv", 
            row.names = T, 
            col.names = T,
            sep = "\t",
            quote = F)
write.table(matrix_norelapse, file = "preprocess/matrix_control_vst.tsv", 
            row.names = T, 
            col.names = T,
            sep = "\t",
            quote = F)
write.table(matrix_relapse, file = "preprocess/matrix_relapse_vst.tsv", 
            row.names = T, 
            col.names = T,
            sep = "\t",
            quote = F)
write.table(metadata_relapse, file = "preprocess/metadata_relapse.tsv", 
            row.names = T, 
            col.names = T,
            sep = "\t",
            quote = F)
write.table(metadata_norelapse, file = "preprocess/metadata_control.tsv", 
            row.names = T, 
            col.names = T,
            sep = "\t",
            quote = F)
