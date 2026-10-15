# Title : Red de regulones y hallmarks
# Author: Keila Adonai Pérez-Fuentes
# Objective: Observar globalmente en qué hallmarks participan los regulones de
# todos los reguladores maestros significativos y cómo se conectan entre sí.
#
# 1. ORA (clusterProfiler::enricher) del regulón de cada regulador maestro con
#    FDR < 0.05 en el MRA (paso 06), contra la colección Hallmark (H) de MSigDB y
#    los Cancer Hallmarks de Menyhart et al. (2025), conjunto core.
#    - Genes del regulón: todos sus genes blanco (activados y reprimidos).
#    - Universo: genes de la matriz VST con la que se construyó la red.
# 2. Corrección de Benjamini-Hochberg GLOBAL: sobre todas las pruebas
#    (regulones x hallmarks) de cada colección, para que el FDR < 0.05 aplique a
#    la red completa y no a cada regulón por separado.
# 3. Red bipartita regulador -> hallmark: una arista por cada enriquecimiento
#    significativo (FDR global < 0.05).
# 4. Proyección hallmark - hallmark: dos hallmarks se conectan si comparten
#    reguladores; peso = índice de Jaccard entre sus conjuntos de reguladores.
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/), p. ej.:
#   Rscript enrichment_analysis/08_red_regulones_hallmarks.R
# Entradas: corto_mra/mra_relapse.RData (paso 06), preprocess/matrix_relapse_vst.tsv (paso 03),
#           conjuntos de genes guardados por los pasos 05a y 05b
# Salidas (en enrichment_analysis/red_regulones_hallmarks/):
#   - ora_regulones_todos.tsv           : todas las pruebas, con FDR global
#   - red_aristas.tsv / red_nodos.tsv   : red bipartita (importables en Cytoscape)
#   - red_hallmarks_aristas.tsv         : proyección hallmark - hallmark
#   - red_regulones_hallmarks.png/.pdf  : figura de la red bipartita
#   - participacion_hallmarks.png/.pdf  : número de regulones por hallmark
#   - red_hallmarks.png/.pdf            : figura de la proyección hallmark - hallmark

# ------------------------------ PACKAGES --------------------------------------
# libglpk (requerida por igraph -> enrichplot -> clusterProfiler) no está en el
# sistema de este servidor; se carga desde un entorno conda del usuario si existe.
glpk <- "/scratch/home/aperez/.conda_envs/glpk/lib/libglpk.so.40"
if (file.exists(glpk)) dyn.load(glpk, local = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(tidygraph)
  library(ggraph)
})

out_dir <- "enrichment_analysis/red_regulones_hallmarks"
dir.create(out_dir, showWarnings = FALSE)

# -------------------------------- DATA ----------------------------------------
load("corto_mra/mra_relapse.RData")  # alls, regulon_relapse, ...

universo <- rownames(read.table("preprocess/matrix_relapse_vst.tsv", header = TRUE,
                                check.names = FALSE))

colecciones <- list(
  "MSigDB Hallmarks" = read_tsv("enrichment_analysis/gsea_msigdb_hallmarks/gene_sets/msigdb_hallmarks_H.tsv",
                                show_col_types = FALSE) %>%
    transmute(term = gs_name, gene = gene_symbol),
  "Cancer Hallmarks (Menyhart, core)" = readLines(
    "enrichment_analysis/gsea_cancer_hallmarks_menyhart/gene_sets/Menyhart_JPA_CancerHallmarks_core.txt",
    warn = FALSE) %>%
    strsplit("\t") %>%
    map_dfr(\(t) { t <- t[nzchar(t)]; tibble(term = t[1], gene = unique(t[-1])) })
)

reguladores <- alls %>%
  filter(FDR < 0.05) %>%
  mutate(n_genes_regulon = map_int(TF, \(tf) length(regulon_relapse[[tf]]$tfmode)),
         actividad = if_else(NES > 0, "Activo en recaída", "Reprimido en recaída"))

cat(sprintf("Reguladores maestros (FDR < 0.05): %d | universo: %d genes\n",
            nrow(reguladores), length(universo)))

# ------------------------------ STEP 1. ORA -----------------------------------
# pvalueCutoff = 1 para conservar todas las pruebas; la corrección se hace
# después de forma global. enricher devuelve NULL si el regulón no comparte
# genes con la colección.
ora <- function(tf, coleccion) {
  res <- suppressWarnings(
    enricher(gene = names(regulon_relapse[[tf]]$tfmode), universe = universo,
             TERM2GENE = colecciones[[coleccion]], pAdjustMethod = "BH",
             pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 5000))
  if (is.null(res)) return(NULL)
  as_tibble(res@result) %>%
    mutate(TF = tf, coleccion = coleccion, .before = 1)
}

razon <- function(r) map_dbl(strsplit(r, "/"), \(x) as.numeric(x[1]) / as.numeric(x[2]))

resultados <- expand_grid(TF = reguladores$TF, coleccion = names(colecciones)) %>%
  pmap_dfr(\(TF, coleccion) ora(TF, coleccion)) %>%
  rename(p.adjust_por_regulon = p.adjust) %>%
  group_by(coleccion) %>%
  mutate(FDR_global = p.adjust(pvalue, method = "BH")) %>%
  ungroup() %>%
  mutate(odds_ratio = razon(GeneRatio) / razon(BgRatio)) %>%
  select(-qvalue) %>%
  left_join(select(reguladores, TF, NES_TF = NES, actividad), by = "TF") %>%
  arrange(coleccion, FDR_global)

write_tsv(resultados, file.path(out_dir, "ora_regulones_todos.tsv"))

# ------------------------------ STEP 2. RED -----------------------------------
siglas <- "\\b(MYC|E2F|G2M|MTORC1|IL2|STAT5|IL6|JAK|STAT3|KRAS|UV|TNFA|NFKB|PI3K|AKT|MTOR|P53|TGF|WNT|DNA|IFN)\\b"
nombre_bonito <- function(x) {
  x |>
    str_remove("^HALLMARK_") |>
    str_replace_all("_", " ") |>
    str_to_sentence() |>
    str_replace_all(regex(siglas, ignore_case = TRUE), toupper) |>
    str_replace("NFKB", "NF-κB") |>
    str_replace("TNFA", "TNFα")
}

aristas <- resultados %>%
  filter(FDR_global < 0.05) %>%
  transmute(from = TF,
            to = paste0(nombre_bonito(ID),
                        if_else(str_detect(coleccion, "Menyhart"), " (M)", "")),
            hallmark_id = ID, coleccion, genes_compartidos = Count,
            odds_ratio, FDR_global, peso = -log10(FDR_global))

nodos <- bind_rows(
  reguladores %>%
    filter(TF %in% aristas$from) %>%
    transmute(name = TF, tipo = "Regulador maestro", grupo = actividad, NES_MRA = NES,
              n_genes = n_genes_regulon),
  aristas %>%
    distinct(name = to, coleccion) %>%
    mutate(tipo = "Hallmark", grupo = coleccion,
           n_genes = map_int(name, \(h) NA_integer_))
) %>%
  left_join(count(aristas, from, name = "grado_tf"), by = c("name" = "from")) %>%
  left_join(count(aristas, to, name = "grado_hallmark"), by = c("name" = "to")) %>%
  mutate(grado = coalesce(grado_tf, grado_hallmark)) %>%
  select(-grado_tf, -grado_hallmark)

write_tsv(aristas, file.path(out_dir, "red_aristas.tsv"))
write_tsv(nodos, file.path(out_dir, "red_nodos.tsv"))

cat(sprintf("\nRed: %d aristas | %d reguladores con al menos un hallmark (de %d) | %d hallmarks\n",
            nrow(aristas), n_distinct(aristas$from), nrow(reguladores), n_distinct(aristas$to)))
aristas %>% count(coleccion, to, sort = TRUE, name = "reguladores") %>% print(n = Inf)

# ------------------- STEP 3. Proyección hallmark - hallmark --------------------
tfs_por_hallmark <- split(aristas$from, aristas$to)
pares <- combn(names(tfs_por_hallmark), 2, simplify = FALSE)
red_hallmarks <- map_dfr(pares, \(p) {
  a <- tfs_por_hallmark[[p[1]]]; b <- tfs_por_hallmark[[p[2]]]
  comunes <- intersect(a, b)
  tibble(from = p[1], to = p[2], reguladores_compartidos = length(comunes),
         jaccard = length(comunes) / length(union(a, b)),
         reguladores = paste(sort(comunes), collapse = ", "))
}) %>%
  filter(reguladores_compartidos > 0) %>%
  arrange(desc(jaccard))

write_tsv(red_hallmarks, file.path(out_dir, "red_hallmarks_aristas.tsv"))
cat("\nHallmarks conectados por reguladores compartidos:\n")
print(red_hallmarks, n = Inf)

# ------------------------------ STEP 4. FIGURAS -------------------------------
colores_nodo <- c("Activo en recaída" = "#B2182B", "Reprimido en recaída" = "#2166AC",
                  "MSigDB Hallmarks" = "#F2B701",
                  "Cancer Hallmarks (Menyhart, core)" = "#3C9D5D")
gris_texto <- "#555555"

# 4.1 Red bipartita
grafo <- tbl_graph(nodes = nodos, edges = select(aristas, from, to, peso), directed = FALSE)

set.seed(123)
fig_red <- ggraph(grafo, layout = "fr", niter = 2000) +
  geom_edge_link(aes(width = peso), colour = "grey70", alpha = 0.7) +
  geom_node_point(aes(colour = grupo, shape = tipo, size = grado)) +
  geom_node_text(aes(label = name, fontface = if_else(tipo == "Hallmark", "bold", "plain")),
                 repel = TRUE, size = 3, max.overlaps = Inf, seed = 123,
                 bg.colour = "white", bg.r = 0.12) +
  scale_edge_width_continuous(range = c(0.3, 2.2), name = expression(-log[10](FDR))) +
  scale_colour_manual(values = colores_nodo, name = NULL) +
  scale_shape_manual(values = c("Regulador maestro" = 16, "Hallmark" = 15), name = NULL) +
  scale_size_continuous(range = c(2.5, 9), name = "Conexiones") +
  labs(title = "Red de reguladores maestros y hallmarks",
       subtitle = sprintf("%d de %d reguladores maestros con enriquecimiento significativo; %d aristas (FDR global < 0.05)",
                          n_distinct(aristas$from), nrow(reguladores), nrow(aristas)),
       caption = paste0(
         "ORA del regulón de cada regulador maestro (FDR < 0.05 en el MRA). (M): Cancer Hallmarks de Menyhart et al. (2025, core); ",
         "sin marca: MSigDB Hallmarks.\nBenjamini-Hochberg sobre todas las pruebas de cada colección. ",
         "Universo: genes de la matriz VST usada para construir la red de regulación.")) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = gris_texto),
        plot.caption = element_text(colour = gris_texto, size = 8, hjust = 0),
        legend.position = "right",
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 10, 10, 10))

ggsave(file.path(out_dir, "red_regulones_hallmarks.png"), fig_red,
       width = 13, height = 10, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "red_regulones_hallmarks.pdf"), fig_red,
       width = 13, height = 10, device = cairo_pdf)

# 4.2 Participación: número de reguladores por hallmark, según su actividad
participacion <- aristas %>%
  left_join(select(reguladores, TF, actividad), by = c("from" = "TF")) %>%
  count(to, coleccion, actividad) %>%
  mutate(to = fct_reorder(to, n, .fun = sum))

fig_part <- ggplot(participacion, aes(n, to, fill = actividad)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = colores_nodo[c("Activo en recaída", "Reprimido en recaída")],
                    name = "Regulador maestro") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), breaks = scales::breaks_pretty()) +
  facet_grid(coleccion ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(coleccion = \(x) str_replace(x, " \\(", "\n("))) +
  labs(title = "Hallmarks en los que participan los regulones",
       x = "Número de reguladores maestros cuyo regulón está enriquecido", y = NULL,
       caption = "FDR global < 0.05 (Benjamini-Hochberg sobre todas las pruebas de cada colección).") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.caption = element_text(colour = gris_texto, size = 8.5, hjust = 0),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
        legend.position = "bottom")

ggsave(file.path(out_dir, "participacion_hallmarks.png"), fig_part,
       width = 9, height = 1.8 + 0.35 * n_distinct(aristas$to), dpi = 300, bg = "white")
ggsave(file.path(out_dir, "participacion_hallmarks.pdf"), fig_part,
       width = 9, height = 1.8 + 0.35 * n_distinct(aristas$to), device = cairo_pdf)

# 4.3 Proyección hallmark - hallmark
nodos_h <- nodos %>% filter(tipo == "Hallmark")
grafo_h <- tbl_graph(nodes = nodos_h,
                     edges = select(red_hallmarks, from, to, jaccard, reguladores_compartidos),
                     directed = FALSE)

set.seed(123)
fig_h <- ggraph(grafo_h, layout = "fr", niter = 2000) +
  geom_edge_link(aes(width = jaccard, label = reguladores_compartidos),
                 colour = "grey60", alpha = 0.8, angle_calc = "along",
                 label_dodge = unit(2.5, "mm"), label_size = 3) +
  geom_node_point(aes(colour = grupo, size = grado)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.4, fontface = "bold",
                 max.overlaps = Inf, seed = 123, bg.colour = "white", bg.r = 0.12) +
  scale_edge_width_continuous(range = c(0.5, 3), name = "Jaccard (reguladores\ncompartidos)") +
  scale_colour_manual(values = colores_nodo, name = NULL) +
  scale_size_continuous(range = c(3, 10), name = "Reguladores") +
  labs(title = "Hallmarks conectados por reguladores maestros compartidos",
       caption = paste0("Dos hallmarks se conectan si sus regulones enriquecidos comparten reguladores maestros; ",
                        "el número sobre la arista indica cuántos.\n",
                        "Hallmarks sin reguladores compartidos aparecen aislados. (M): Menyhart et al. (2025, core).")) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.caption = element_text(colour = gris_texto, size = 8, hjust = 0),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 10, 10, 10))

ggsave(file.path(out_dir, "red_hallmarks.png"), fig_h,
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "red_hallmarks.pdf"), fig_h,
       width = 10, height = 8, device = cairo_pdf)
