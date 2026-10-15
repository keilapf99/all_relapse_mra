# Title : Gene Set Enrichment Analysis to LLA differential expression genes
# Author: Keila Adonai Pérez-Fuentes
# Date: May 29th 2026
# Objective: Enrich the differential expression genes in ALL relapse in cancer 
# hallmarks 

# ------------------------------ PACKAGES --------------------------------------
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot) 
library(ggplot2)
# install.packages("msigdbr")
library(msigdbr)
library(stringr)
library(viridis)
library(ggnewscale)
# install.packages("tidytext")
library(tidytext)
library(genekitr)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(tidytext)  # reorder_within()


# -------------------------------- DATA ----------------------------------------
degs_totales <- read.table(
  "~/Desktop/aperez/all_relapse_masters_degree/2_DEA/DEA_totales_relapse.tsv",
  header = TRUE,
  sep = "\t",
  quote = "")

hallmark_genes <- msigdbr(species = "Homo sapiens", collection = "H")
head(hallmark_genes)
library(msigdbr)

hallmarks <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

unique(hallmarks[, c("gs_name", "gs_description")])

# ---------------------- STEP 1. Vector needed to GSEA -------------------------
gene_list <- degs_totales$stat
names(gene_list) <- degs_totales$gene_symbol

# Quitar NA
gene_list <- na.omit(gene_list)

# Convertir SYMBOL → ENTREZ
gene_ids <- bitr(names(gene_list),
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

# Reconstruir vector con ENTREZ
gene_list_entrez <- gene_list[gene_ids$SYMBOL]
names(gene_list_entrez) <- gene_ids$ENTREZID

# Ordenar (MUY importante)
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)


# 2. Run GSEA analysis ---------------------------------------------------------

# Preparar TERM2GENE correctamente
hallmark_df <- hallmark_genes[, c("gs_name", "ncbi_gene")]
hallmark_df$ncbi_gene <- as.character(hallmark_df$ncbi_gene)

# Ejecutar GSEA
gsea_hallmark <- GSEA(geneList = gene_list_entrez, 
                      TERM2GENE = hallmark_df,
                      pvalueCutoff = 0.05,
                      verbose = TRUE)


# 3. Data plots ----------------------------------------------------------------

# Extraer resultados como data.frame
gsea_res <- gsea_hallmark@result

# Filtrar vías significativas y top por NES
top_gsea <- gsea_res %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(abs(NES)))


# Top 5 pathways por NES
top5 <- top_gsea %>%
  slice_max(NES, n = 10) %>%
  mutate(
    pathway = Description |>
      gsub("HALLMARK_", "", x = _) |>
      gsub("_", " ", x = _) |>
      tools::toTitleCase(),
    
    # ordenar de mayor a menor NES
    pathway = fct_reorder(pathway, NES, .desc = TRUE)
  )

# Lollipop plot top 5
lolli_top5 <- ggplot(
  top5,
  aes(
    x = NES,
    y = pathway
  )
) +
  
  geom_segment(
    aes(
      x = 0,
      xend = NES,
      y = pathway,
      yend = pathway
    ),
    linewidth = 1.8,
    color = "grey80"
  ) +
  
  geom_point(
    aes(size = -log10(p.adjust)),
    color = "#8E44AD"
  ) +
  
  geom_text(
    aes(label = sprintf("%.2f", NES)),
    nudge_x = 0.15,   # más separación respecto al punto
    hjust = 0,
    size = 4,
    fontface = "bold"
  ) +
  
  scale_size_continuous(
    name = expression(-log[10](FDR)),
    range = c(5, 12)
  ) +
  
  coord_cartesian(
    xlim = c(0, max(top5$NES) + 0.5)
  ) +
  
  labs(
    title = "Top 10 Hallmark Gene Sets según NES",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    
    axis.text.y = element_text(
      face = "bold",
      color = "black",
      size = 11
    ),
    
    axis.text.x = element_text(
      color = "black"
    ),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.position = "right"
  )

ggsave(
  filename = "top10_lolliplot_GSEA.png",
  plot = lolli_top10,
  width = 10,   # recomiendo agrandar un poco
  height = 8,
  dpi = 300
)


