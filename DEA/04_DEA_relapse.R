# Title : Differential Expression Analysis in childhood ALL relapse
# Author: Keila Adonai Pérez-Fuentes
# Date: January 20th 2025
# Objective: This script performs a differential expression analysis using DESeq2
# between two different phenotypes, ALL no relapse vs ALL relapse and also explore the 
# probably sex effect, female and male, to identify any significant differences
# that could help explain the mechanisms of relapse

# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/), p. ej.:
#   Rscript DEA/04_DEA_relapse.R
# Entradas: raw_data/raw_matrix.tsv, preprocess/metadata_NS500560.tsv (paso 01)
# Salidas : DEA/DEA_totales_relapse.tsv (firma completa), DEA/DEGS_significant_relapse.tsv,
#           DEA/DEA_relapse_vp.png, DEA/DEA_relapse_vp2.png

# Packages ---------------------------------------------------------------------
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("org.Hs.eg.db")
library(ggrepel)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Data ------------------------------------------------------------------------ 
counts_lla <- read.table("raw_data/raw_matrix.tsv",
                         sep = "\t",
                         header = T,
                         check.names = F)
metadata <- read.table("preprocess/metadata_NS500560.tsv", 
                       header = T, 
                       sep = "\t")

# Comparing relapse vs no relapse ~ sex

# STEP 1. Input files ------------------------------------------------------------------
#Prepare the data into the input files needed

# The counts data need to match with the metadata names
# 1.1 Numeric matrix
counts <- as.matrix(counts_lla)
class(counts)
unique(sapply(counts, class))

# 1.2 Verify identical samples order
rownames(metadata) <- metadata$sample
counts <- counts[, rownames(metadata), drop = FALSE]
stopifnot(identical(colnames(counts), rownames(metadata)))
identical(colnames(counts), rownames(metadata))

#Convert the metadata into factors to log2fc be calculating correctly 
# Level 1 = Reference
# Level 2 = Comparison 
metadata$relapse <- factor(metadata$relapse, levels = c("no", "si")) 

# STEP 2. Create dds object ----------------------------------------------------

dds_relapse <- DESeqDataSetFromMatrix(countData = counts,
                                      colData = metadata, 
                                      design = ~relapse)

# Preprocessing ----------------------------------------------------------------

# Filtering out lowly expressed genes
# Al menos 10 lecturas en al menos 15 muestras (mismo filtro que la VST, paso 03)
keep <- rowSums(counts(dds_relapse) >= 10) >= 15
dds_relapse <- dds_relapse[keep, ]

# Normalization and differential expression ------------------------------------
dds_relapse <- DESeq(dds_relapse)

# STEP 3. Obtain the results ---------------------------------------------------

# Note: it is important the order of the comparison 
# Example = res <- results(dds, contrast=c("condition", "comparison", "reference"))

res <- results(dds_relapse, contrast=c("relapse", "si", "no"))
res_padj_order <- res[order(res$padj),]
res_order <- data.frame(gene = rownames(res_padj_order), res_padj_order) # SAVE
head(res_order)

# Filter by |log2FC| > 1 and FDR < 0.05
res_subset <- subset(res, abs(res$log2FoldChange) >= 1)
res_subset <- subset(res_subset, res_subset$padj < 0.05)
res_subset <- res_subset[order(res_subset$log2FoldChange),]
res_filter_order <- data.frame(gene = rownames(res_subset), res_subset) # SAVE

# Change to gene_symbol
gsymbols <- mapIds(org.Hs.eg.db, keys = rownames(res_filter_order), 
                   keytype = "ENSEMBL", 
                   column = "SYMBOL")

gsymbols_results <-bind_cols(gsymbols, res_filter_order)
dim(gsymbols_results)
gsymbols_results <- gsymbols_results[!is.na(gsymbols_results$...1),]
dim(gsymbols_results)
# [1] 177   8

colnames(gsymbols_results)[1] <- "gene_symbol" # save

# STEP 4. Plots ----------------------------------------------------------------

# Convertir res a data.frame
df_genes <- as.data.frame(res)
df_genes$ensembl_id <- rownames(df_genes)

# Eliminar genes sin padj
df_genes <- df_genes[!is.na(df_genes$padj), ]

df_genes$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = df_genes$ensembl_id,
  keytype = "ENSEMBL",
  column = "SYMBOL"
)
df_genes <- df_genes[!is.na(df_genes$gene_symbol), ] # save

df_genes$Type <- ifelse(
  df_genes$log2FoldChange > 1 & df_genes$padj < 0.05, "Upregulated",
  ifelse(
    df_genes$log2FoldChange < -1 & df_genes$padj < 0.05, "Downregulated",
    "Not Significant"
  )
)

df_genes_filtered <- df_genes %>%
  filter(Type != "Not Significant")

gene_counts <- table(df_genes$Type)

color_labels <- c(
  paste0("Upregulated (", gene_counts["Upregulated"], ")"),
  paste0("Downregulated (", gene_counts["Downregulated"], ")"),
  paste0("Not Significant (", gene_counts["Not Significant"], ")")
)

names(color_labels) <- c("Upregulated", "Downregulated", "Not Significant")

volcano_plot <- ggplot(
  df_genes,
  aes(x = log2FoldChange, y = -log10(padj), color = Type)
) +
  geom_point(size = 3, alpha = 0.8) +
  
  geom_text_repel(
    data = subset(df_genes, Type != "Not Significant"),
    aes(label = gene_symbol),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = 50
  ) +
  
  scale_color_manual(
    values = c(
      "Upregulated" = "brown1",
      "Downregulated" = "blue",
      "Not Significant" = "gray"
    ),
    labels = color_labels
  ) +
  
  labs(
    title = "DEGs at Diagnosis: Relapse vs Non-relapse",
    x = "Log2 Fold Change",
    y = "-log10 adjusted p-value"
  ) +
  
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")



ggsave("DEA/DEA_relapse_vp.png", volcano_plot, width = 8, height = 6, dpi = 300)



# Genes de interés a etiquetar en el segundo volcano plot
genes_volcano <- c("BNIP3L", "CD200R1", "CD36", "CXCL12", "HAVCR2", "ICOSLG",
                   "LGALS3", "MEIS1", "MLLT3", "NAMPT", "PBX1", "TAL1", "TGFBI",
                   "THBS1", "TNF", "TTK", "VCAM1", "WIPI1", "CYP1B1")

volcano_plot <- ggplot(
  df_genes,
  aes(x = log2FoldChange, y = -log10(padj), color = Type)
) +
  geom_point(size = 3, alpha = 0.8) +
  
  geom_text_repel(
    data = subset(df_genes, gene_symbol %in% genes_volcano),
    aes(label = gene_symbol),
    color = "black",
    size = 4,
    box.padding = 0.4,      # 🔥 MÁS espacio alrededor del texto
    point.padding = 0.3,    # 🔥 empuja lejos del punto
    segment.color = "black",
    segment.size = 0.3,
    force = 15,              # 🔥 separa más las etiquetas
    nudge_y = 0.3,          # 🔥 empujón vertical inicial
    arrow = arrow(
      length = unit(0.01, "npc"),
      type = "closed",
      ends = "last"
    ),
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = c(
      "Upregulated" = "brown1",
      "Downregulated" = "blue",
      "Not Significant" = "gray"
    ),
    labels = color_labels
  ) +
  
  labs(
    title = "DEGs al diagnóstico: recaída vs no recaída",
    x = "Log2 Fold Change",
    y = "-log10 adjusted p-value"
  ) +
  
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")

ggsave("DEA/DEA_relapse_vp2.png", volcano_plot, width = 8, height = 6, dpi = 300)
# STEP 5. Save -----------------------------------------------------------------
write.table(df_genes, file = "DEA/DEA_totales_relapse.tsv",
            sep = "\t", 
            quote = F,
            col.names = T)
write.table(gsymbols_results, file = "DEA/DEGS_significant_relapse.tsv",
            sep = "\t", 
            quote = F,
            col.names = T)
