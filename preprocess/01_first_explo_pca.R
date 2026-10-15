# Title : PCA for batches effect
# Author: Keila Adonai Pérez-Fuentes
# Date: January 14th, 2026
# Objective: To generate PCA plots for the identification and assessment of batch effects.
# This script normalizes RNA-seq data and removes unwanted sources of variation,
# including unknown effects, prior to downstream analyses.

# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/), p. ej.:
#   Rscript preprocess/01_first_explo_pca.R
# Entradas: raw_data/raw_matrix.tsv, raw_data/raw_metadata.tsv, raw_data/batches_metadata.tsv
# Salidas : preprocess/metadata_NS500560.tsv, preprocess/PCA_TMM_*.pdf

# Packages ---------------------------------------------------------------------
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(edgeR)
library(ggrepel)

# Data -------------------------------------------------------------------------

counts <- read.table(
 file = "raw_data/raw_matrix.tsv",
 header = T,
 check.names = F)

metadata  <- read.table(
  file = "raw_data/raw_metadata.tsv",
  header = T,
  check.names = F)

metadata_batches <- read.table(
  file = "raw_data/batches_metadata.tsv",
  header = T,
  check.names = F)

####################### 1. Counts matrix preprocess #############################

# 1.1 Chunk genes/Duplicated genes ----------------------------------------------

# ------------------------------------------------- 1.1.1 Convert to gene symbol

gsymbols <- mapIds(org.Hs.eg.db, keys = rownames(counts), 
                   keytype = "ENSEMBL", 
                   column = "SYMBOL")

counts <-bind_cols(gsymbols, counts)
dim(counts)
# [1] 78932    52
counts <- counts[!is.na(counts$...1),]
dim(counts)
#[1] 36922    52


# ---------------------------------------------------------- 1.2 Media by sample

# Define the function to obtain the media of gene expression by sample ----
# Function to calculate the media by column

colMean <- function(v) {
  apply(v, 2, mean)
}

counts <- do.call(rbind, lapply(split(counts, counts[,1]), 
                                function(chunk){colMean(chunk[,-1])}))

dim(counts)
# [1] 38367    50


# 1.2.1 Filter low counts genes ------------------------------------------------
counts <- counts[rowSums(counts >= 10) >= 38, ] #at least the 70% of your samples
dim(counts)
#[1]  12143    50

# 1.3. Verification of input data ------------------------------------------------
# ---------------------------------------------------------- 1.3.1 Numeric matrix
class(counts)
unique(sapply(counts, class))

############################ 2. Metadata format ################################

# -------------------- 2.1 Delete samples_id format to match with input_metadata
metadata_batches <- metadata_batches %>%
  mutate(
    sample = sample %>%
      str_remove("_S[0-9]+$") %>%        # quitar _Sx
      str_remove("^[0-9]{6}_") %>%       # quitar prefijos fecha
      toupper()
  )

metadata <- metadata %>%
  rownames_to_column("sample")

metadata_full <- metadata %>%
  left_join(metadata_batches, by = "sample")


# ------------------------------------------- 2.2 Verify identical samples order
setdiff(metadata_batches$sample, metadata$sample)
# [1] "145"    "15"     "169"    "197"    "283"    "327-16" "461-16" "540" 
metadata_batches %>%
  filter(sample %in% setdiff(metadata_batches$sample, metadata$sample)) %>%
  dplyr::select(sample, instrument)

# sample instrument
# 1    145    VH00706
# 2     15    VH00706
# 3    169    VH00706
# 4    197    VH00706
# 5    283   NS500560
# 6 327-16   NB552515
# 7 461-16   NB552515
# 8    540   NS500560

common_samples <- intersect(colnames(counts), metadata_full$sample)

stopifnot(length(common_samples) > 2)

counts <- counts[, common_samples, drop = FALSE]

metadata_full <- metadata_full %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, colnames(counts)))

identical(colnames(counts), metadata_full$sample)
# [1] TRUE

############################ 3. TMM normalization ##############################

counts <- as.matrix(counts)
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")
counts_tmm <- cpm(
  dge,
  log = F,
  prior.count = 1,
  normalized.lib.sizes = TRUE
)

# write.table(counts_tmm, 
#             file = "/home/aperez/maestria_leukemia/raw_data/matrix_forPCA_tmm.tsv",
#             sep = "\t",
#             col.names = T)

################################# 4. PCA Plots #################################

counts_t <- t(counts_tmm)
pca <- prcomp(counts_t, center = TRUE, scale. = TRUE)
# Varianza explicada
pca_var <- summary(pca)$importance[2, 1:2] * 100

# ----------------------------------------------------------- 4.1 PCA + metadata
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  left_join(metadata_full, by = "sample")

## -------------------------------------------------------- 4.2 Format for facet
pca_long <- pca_df %>%
  dplyr::select(sample, PC1, PC2, instrument, relapse, flowcell) %>%
  pivot_longer(
    cols = c(instrument, relapse, flowcell),
    names_to = "variable",
    values_to = "group"
  )

# ------------------------------------------------------------ 4.3 Visualization

pdf("preprocess/PCA_TMM_instrument.pdf", width = 7, height = 6)

ggplot(pca_df, aes(PC1, PC2, color = instrument)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_brewer(palette = "Set3") +
  labs(
    title = "PCA (TMM) – Instrument",
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)")
  ) +
  theme_bw()

dev.off()


pdf("preprocess/PCA_TMM_flowcell.pdf", width = 7, height = 6)

ggplot(pca_df, aes(PC1, PC2, color = flowcell)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "PCA (TMM) – Flowcell",
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)")
  ) +
  theme_bw()

dev.off()


pdf("preprocess/PCA_TMM_relapse.pdf", width = 7, height = 6)

ggplot(pca_df, aes(PC1, PC2, color = relapse)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(
    values = c(
      "si" = "#D73027",   # rojo fuerte
      "no" = "#4575B4"    # azul fuerte
    )
  ) +
  labs(
    title = "PCA (TMM) – Relapse",
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)")
  ) +
  theme_bw()

dev.off()

########################## Relapse and instruments #############################

table(metadata_full$instrument, metadata_full$relapse)
# no si
# NB552515  9  0
# NS500560 15 21
# VH00706   3  0
# VH01014   3  0

pdf("preprocess/PCA_TMM_relapse_instrument_shape.pdf", width = 7, height = 6)

ggplot(pca_df, aes(PC1, PC2, color = relapse, shape = instrument)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(
    values = c("si" = "#D73027", "no" = "#4575B4")
  ) +
  labs(
    title = "PCA (TMM) – Relapse e instrumento",
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)")
  ) +
  theme_bw()

dev.off()

metadata_full_NS500560 <- metadata_full %>%
  filter(instrument == "NS500560")

############################## SAVE DATA #######################################

write.table(
  metadata_full_NS500560,
  file = "preprocess/metadata_NS500560.tsv",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)




