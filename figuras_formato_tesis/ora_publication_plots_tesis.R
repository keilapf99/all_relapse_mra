# ============================================================================
# ora_top10_master_regulators_publication_plots.r
#
# Horizontal bar chart: one panel per master regulator (stacked vertically).
# aes(x = gene_ratio, y = hallmark, fill = TF, alpha = FDR)
# Fill color differs by panel; alpha encodes FDR significance.
#
# Output: ora_hallmarks_pathways_by_regulon.{png,pdf}
#
# Copia de all_relapse_masters_degree/4_GSEA_ORA/ora_top10_master_regulators_publication_plots.r
# Cambios: lee y escribe en `out_dir` (lo define F6_ora_regulones_hallmarks.R,
# que hace source() de este archivo) y cada regulador tiene un color fijo, el
# mismo que en la figura de la tesis (Dark2 de RColorBrewer), aunque cambie
# el número de reguladores con asociaciones significativas.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(forcats)
})

HAS_GGSCI   <- requireNamespace("ggsci",        quietly = TRUE)
HAS_RBREWER <- requireNamespace("RColorBrewer", quietly = TRUE)

# ---- Remove previous outputs ------------------------------------------------
old_files <- file.path(out_dir, c(
  "ora_hallmarks_pathways_by_regulon.png",
  "ora_hallmarks_pathways_by_regulon.pdf"
))
invisible(file.remove(old_files[file.exists(old_files)]))

# ---- Label correction -------------------------------------------------------
fix_label <- function(x) {
  x %>%
    str_replace_all("\\bMtorc1\\b",  "mTORC1") %>%
    str_replace_all("\\bMtor\\b",    "mTOR")   %>%
    str_replace_all("\\bStat3\\b",   "STAT3")  %>%
    str_replace_all("\\bIl6\\b",     "IL6")    %>%
    str_replace_all("\\bJak\\b",     "JAK")    %>%
    str_replace_all("\\bKras\\b",    "KRAS")   %>%
    str_replace_all("\\bMyc\\b",     "MYC")    %>%
    str_replace_all("\\bPi3k\\b",    "PI3K")   %>%
    str_replace_all("\\bAkt\\b",     "AKT")    %>%
    str_replace_all("\\bP53\\b",     "p53")    %>%
    str_replace_all("\\bDna\\b",     "DNA")    %>%
    str_replace_all("\\bE2f\\b",     "E2F")    %>%
    str_replace_all("\\bG2m\\b",     "G2M")    %>%
    str_replace_all("\\bUv\\b",      "UV")
}

# ---- Load data --------------------------------------------------------------
ora_sig <- read.csv(
  file.path(out_dir, "ora_hallmarks_top10_master_regulators_significant.csv"),
  check.names = FALSE
) %>%
  mutate(hallmark_label = fix_label(hallmark_label))

universe_size <- length(readLines(file.path(out_dir, "ora_regulon_universe_genes.txt")))

master_regs <- read.csv(
  file.path(out_dir, "ora_top10_master_regulators.csv"),
  check.names = FALSE
) %>%
  transmute(TF, MRA_NES = NES)

if (nrow(ora_sig) == 0) stop("No significant ORA associations (FDR < 0.05).")

# TFs with >= 1 significant association, ordered by |MRA_NES| descending.
# This defines both the panel order (top = highest |NES|) and the fill-palette
# assignment (first TF gets palette color 1, etc.).
tf_order <- master_regs %>%
  filter(TF %in% unique(ora_sig$TF)) %>%
  arrange(desc(abs(MRA_NES))) %>%
  pull(TF)

# ---- Prepare plot data ------------------------------------------------------
# Keep only the 5 most significant hallmarks per TF (lowest p.adjust).
# Within each panel hallmarks are ordered by gene_ratio ascending so the
# highest-ratio pathway lands at the top after horizontal rendering.
fig_data <- ora_sig %>%
  mutate(TF = factor(TF, levels = tf_order)) %>%
  group_by(TF) %>%
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
  mutate(hallmark_label = fct_reorder(hallmark_label, gene_ratio)) %>%
  ungroup()

# ---- Fill colour scale (one distinct colour per TF panel) -------------------
# Colores de la figura de la tesis (all_relapse se generó sin ggsci, con Dark2)
tf_colors_tesis <- c(
  LTF    = "#1B9E77", E4F1   = "#D95F02", HMG20B = "#7570B3",
  KLF16  = "#E7298A", ZNF282 = "#66A61E", ZNF581 = "#E6AB02",
  ZNF324 = "#A6761D", CEBPD  = "#666666", ATF6B  = "#7D3C98"
)

if (all(tf_order %in% names(tf_colors_tesis))) {
  fill_scale <- scale_fill_manual(
    values = tf_colors_tesis[tf_order],
    guide  = guide_legend(
      title          = "Master Regulator",
      title.position = "top",
      nrow           = 3,
      keywidth       = unit(0.85, "lines"),
      keyheight      = unit(0.70, "lines")
    )
  )
} else if (HAS_GGSCI) {
  fill_scale <- ggsci::scale_fill_d3(
    palette = "category10",
    guide   = guide_legend(
      title          = "Master Regulator",
      title.position = "top",
      nrow           = 3,
      keywidth       = unit(0.85, "lines"),
      keyheight      = unit(0.70, "lines")
    )
  )
} else if (HAS_RBREWER) {
  pal_vec <- c(RColorBrewer::brewer.pal(8, "Dark2"), "#7D3C98")
  fill_scale <- scale_fill_manual(
    values = setNames(pal_vec[seq_along(tf_order)], tf_order),
    guide  = guide_legend(
      title          = "Master Regulator",
      title.position = "top",
      nrow           = 3,
      keywidth       = unit(0.85, "lines"),
      keyheight      = unit(0.70, "lines")
    )
  )
} else {
  # Built-in fallback: D3 category10 hex values
  d3_10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
             "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")
  fill_scale <- scale_fill_manual(
    values = setNames(d3_10[seq_along(tf_order)], tf_order),
    guide  = guide_legend(
      title          = "Master Regulator",
      title.position = "top",
      nrow           = 3,
      keywidth       = unit(0.85, "lines"),
      keyheight      = unit(0.70, "lines")
    )
  )
}

# ---- Alpha scale: continuous FDR (trans = "reverse" so low FDR = high alpha) -
# Requires scales >= 1.3.0.  Falls back to swapped range if unavailable.
alpha_scale <- tryCatch(
  scale_alpha_continuous(
    range  = c(0.3, 1),
    trans  = "reverse",
    name   = "FDR",
    guide  = guide_legend(
      title.position = "top",
      title.hjust    = 0,
      keywidth       = unit(1.2, "lines"),
      keyheight      = unit(0.75, "lines"),
      override.aes   = list(fill = "grey35", colour = NA)
    )
  ),
  error = function(e) {
    scale_alpha_continuous(
      range  = c(1, 0.3),
      name   = "FDR",
      guide  = guide_legend(
        title.position = "top",
        title.hjust    = 0,
        keywidth       = unit(1.2, "lines"),
        keyheight      = unit(0.75, "lines"),
        override.aes   = list(fill = "grey35", colour = NA)
      )
    )
  }
)

# ---- Build plot -------------------------------------------------------------
# aes(x = gene_ratio, y = hallmark, fill = TF, alpha = p.adjust) exactly as
# requested. ggplot2 3.3+ auto-detects horizontal orientation when x is
# continuous and y is discrete, so no coord_flip() is needed.

p_ora <- ggplot(
  fig_data,
  aes(
    x     = gene_ratio,
    y     = hallmark_label,
    fill  = TF,
    alpha = p.adjust
  )
) +
  geom_col(
    colour = NA,      # no bar borders
    width  = 0.70     # moderately thin bars with natural spacing
  ) +
  fill_scale +
  alpha_scale +
  scale_x_continuous(
    name   = "Gene Ratio",
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.10))
  ) +
  # One panel per TF, stacked vertically; panel height proportional to
  # number of hallmarks per TF
  facet_grid(
    TF ~ .,
    scales = "free_y",
    space  = "free_y"
  ) +
  labs(
    title    = "Programas transcripcionales enriquecidos en los regulones de los reguladores maestros de recaida",
    subtitle = paste0(
      "Top 5 vias Hallmark mas significativas por regulador (FDR ajustado BH < 0.05)  |  ",
      "Universo: ", format(universe_size, big.mark = ","), " genes"
    ),
    y        = NULL,
    caption  = "Bar length = gene ratio. Opacity = FDR: solid (low) = most significant."
  ) +
  theme_bw(base_size = 11) +
  theme(
    # --- Title block ---
    plot.title         = element_text(
      face   = "bold", size = 13, colour = "#1C2B3A",
      margin = margin(b = 4)
    ),
    plot.subtitle      = element_text(
      size   = 9, colour = "#546E7A",
      margin = margin(b = 10)
    ),
    plot.caption       = element_text(
      size   = 7.5, colour = "#90A4AE",
      hjust  = 0, margin = margin(t = 8)
    ),
    plot.margin        = margin(16, 10, 12, 22),
    plot.background    = element_rect(fill = "white", colour = NA),
    # --- Facet strip: TF name on the right, rotated 90° ---
    strip.background   = element_rect(fill = "#D8D8D8", colour = "#555555", linewidth = 0.6),
    strip.text.y       = element_text(
      angle  = -90,
      face   = "bold",
      size   = 9,
      colour = "#1C2B3A"
    ),
    strip.placement    = "outside",
    # --- Axes ---
    axis.title.x       = element_text(
      face   = "bold", size = 11, colour = "#1C2B3A",
      margin = margin(t = 6)
    ),
    axis.text.x        = element_text(size = 9,  colour = "#546E7A"),
    axis.text.y        = element_text(size = 10, colour = "#0D1117", face = "plain"),
    axis.ticks.y       = element_blank(),
    axis.ticks.x       = element_line(colour = "grey80", linewidth = 0.3),
    # --- Grid: horizontal lines only ---
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.30),
    # --- Panel border ---
    panel.border       = element_rect(colour = "#555555", fill = NA, linewidth = 0.65),
    panel.background   = element_rect(fill = "white", colour = NA),
    panel.spacing.y    = unit(0.04, "lines"),
    # --- Legends ---
    legend.position    = "right",
    legend.box         = "vertical",
    legend.title       = element_text(face = "bold", size = 9.5, colour = "#1C2B3A"),
    legend.text        = element_text(size = 9,  colour = "#1C2B3A"),
    legend.key.size    = unit(0.80, "lines"),
    legend.background  = element_rect(fill = "white", colour = NA),
    legend.margin      = margin(0, 4, 0, 4)
  )

# ---- Dynamic figure dimensions ----------------------------------------------
# Height: each of the 101 TF-hallmark bars gets ~0.22 in + facet strip overhead
n_bars     <- nrow(fig_data)
n_panels   <- length(tf_order)
fig_height <- n_bars * 0.22 + n_panels * 0.18 + 3.8
fig_width  <- 12

ggsave(
  file.path(out_dir, "ora_hallmarks_pathways_by_regulon.png"),
  p_ora,
  width  = fig_width,
  height = fig_height,
  dpi    = 300,
  bg     = "white"
)
ggsave(
  file.path(out_dir, "ora_hallmarks_pathways_by_regulon.pdf"),
  p_ora,
  width  = fig_width,
  height = fig_height,
  bg     = "white"
)

message(sprintf(
  "\nFigure saved: ora_hallmarks_pathways_by_regulon.{png,pdf}  [%g x %.1f in]",
  fig_width, fig_height
))

# ============================================================================
# LANDSCAPE VERSION — facet_wrap 3 columns, 14 x 10 in
# For thesis figures and presentation slides.
# Reuses fig_data, fill_scale and alpha_scale built above.
# ============================================================================

p_ora_3col <- ggplot(
  fig_data,
  aes(
    x     = gene_ratio,
    y     = hallmark_label,
    fill  = TF,
    alpha = p.adjust
  )
) +
  geom_col(
    colour = NA,
    width  = 0.62     # slightly thinner to absorb the fixed-height panels
  ) +
  fill_scale +
  alpha_scale +
  # Override fill legend for right-side placement: single column
  guides(
    fill = guide_legend(
      title          = "Master Regulator",
      title.position = "top",
      ncol           = 1,
      keywidth       = unit(0.85, "lines"),
      keyheight      = unit(0.65, "lines")
    )
  ) +
  scale_x_continuous(
    name   = "Gene Ratio",
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.12))
  ) +
  facet_wrap(
    ~ TF,
    scales = "free_y",
    ncol   = 3
  ) +
  labs(
    title    = "Programas transcripcionales enriquecidos en los regulones de los reguladores maestros de recaida",
    subtitle = paste0(
      "Top 5 vias Hallmark mas significativas por regulador (FDR ajustado BH < 0.05)  |  ",
      "Universo: ", format(universe_size, big.mark = ","), " genes"
    ),
    y        = NULL,
    caption  = "Bar length = gene ratio. Opacity = FDR: solid (low) = most significant."
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    # --- Title block ---
    plot.title         = element_text(
      face   = "bold", size = 13.5, colour = "#1C2B3A",
      margin = margin(b = 4)
    ),
    plot.subtitle      = element_text(
      size   = 8.5, colour = "#546E7A",
      margin = margin(b = 10)
    ),
    plot.caption       = element_text(
      size   = 7.5, colour = "#90A4AE",
      hjust  = 0, margin = margin(t = 8)
    ),
    plot.margin        = margin(14, 12, 12, 14),
    plot.background    = element_rect(fill = "white", colour = NA),
    # --- Strip at TOP of each panel (facet_wrap default) ---
    strip.background   = element_rect(fill = "#D8D8D8", colour = "#555555", linewidth = 0.6),
    strip.text         = element_text(
      face   = "bold",
      size   = 10.5,
      colour = "#1C2B3A",
      margin = margin(t = 3, b = 3)
    ),
    # --- Axes ---
    axis.title.x       = element_text(
      face   = "bold", size = 10.5, colour = "#1C2B3A",
      margin = margin(t = 6)
    ),
    axis.text.x        = element_text(size = 8.5, colour = "#546E7A"),
    axis.text.y        = element_text(size = 9,   colour = "#0D1117"),
    axis.ticks.y       = element_blank(),
    axis.ticks.x       = element_line(colour = "grey80", linewidth = 0.3),
    # --- Grid: horizontal lines only ---
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.30),
    # --- Panel borders ---
    panel.border       = element_rect(colour = "#555555", fill = NA, linewidth = 0.65),
    panel.background   = element_rect(fill = "white", colour = NA),
    panel.spacing.x    = unit(0.80, "lines"),
    panel.spacing.y    = unit(0.80, "lines"),
    # --- Legend ---
    legend.position    = "right",
    legend.box         = "vertical",
    legend.title       = element_text(face = "bold", size = 9.5, colour = "#1C2B3A"),
    legend.text        = element_text(size = 9,   colour = "#1C2B3A"),
    legend.key.size    = unit(0.78, "lines"),
    legend.background  = element_rect(fill = "white", colour = NA),
    legend.margin      = margin(0, 4, 0, 4),
    legend.spacing.y   = unit(0.6, "lines")
  )

ggsave(
  file.path(out_dir, "ora_hallmarks_pathways_by_regulon_3col.png"),
  p_ora_3col,
  width  = 14,
  height = 10,
  dpi    = 300,
  bg     = "white"
)
ggsave(
  file.path(out_dir, "ora_hallmarks_pathways_by_regulon_3col.pdf"),
  p_ora_3col,
  width  = 14,
  height = 10,
  bg     = "white"
)

message("\nLandscape version saved: ora_hallmarks_pathways_by_regulon_3col.{png,pdf}  [14 x 10 in]")

# ============================================================================
# CANVA SLIDE VERSION — 16:9 widescreen (13.33 x 7.5 in), facet_wrap 3x3,
# top 3 hallmarks per TF, larger fonts for legibility on a slide.
# ============================================================================

# Keep only the 3 most significant hallmarks per TF for this compact layout.
fig_data_canva <- ora_sig %>%
  mutate(TF = factor(TF, levels = tf_order)) %>%
  group_by(TF) %>%
  slice_min(order_by = p.adjust, n = 3, with_ties = FALSE) %>%
  mutate(hallmark_label = fct_reorder(hallmark_label, gene_ratio)) %>%
  ungroup()

p_ora_canva <- ggplot(
  fig_data_canva,
  aes(
    x     = gene_ratio,
    y     = hallmark_label,
    fill  = TF,
    alpha = p.adjust
  )
) +
  geom_col(
    colour = NA,
    width  = 0.65
  ) +
  fill_scale +
  alpha_scale +
  guides(
    fill = guide_legend(
      title          = "Master Regulator",
      title.position = "top",
      ncol           = 1,
      keywidth       = unit(1.0, "lines"),
      keyheight      = unit(0.85, "lines")
    )
  ) +
  scale_x_continuous(
    name   = "Gene Ratio",
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.12))
  ) +
  facet_wrap(
    ~ TF,
    scales = "free_y",
    ncol   = 3
  ) +
  labs(
    title    = str_wrap(
      "Programas transcripcionales enriquecidos en los regulones de los reguladores maestros de recaida",
      width = 70
    ),
    subtitle = paste0(
      "Top 3 vias Hallmark mas significativas por regulador (FDR ajustado BH < 0.05)  |  ",
      "Universo: ", format(universe_size, big.mark = ","), " genes"
    ),
    y        = NULL,
    caption  = "Bar length = gene ratio. Opacity = FDR: solid (low) = most significant."
  ) +
  theme_bw(base_size = 14) +
  theme(
    # --- Title block ---
    plot.title         = element_text(
      face   = "bold", size = 18, colour = "#1C2B3A",
      margin = margin(b = 4)
    ),
    plot.subtitle      = element_text(
      size   = 12, colour = "#546E7A",
      margin = margin(b = 10)
    ),
    plot.caption       = element_text(
      size   = 10, colour = "#90A4AE",
      hjust  = 0, margin = margin(t = 8)
    ),
    plot.margin        = margin(16, 14, 12, 14),
    plot.background    = element_rect(fill = "white", colour = NA),
    # --- Strip at TOP of each panel (facet_wrap default) ---
    strip.background   = element_rect(fill = "#D8D8D8", colour = "#555555", linewidth = 0.6),
    strip.text         = element_text(
      face   = "bold",
      size   = 13,
      colour = "#1C2B3A",
      margin = margin(t = 4, b = 4)
    ),
    # --- Axes ---
    axis.title.x       = element_text(
      face   = "bold", size = 13, colour = "#1C2B3A",
      margin = margin(t = 6)
    ),
    axis.text.x        = element_text(size = 10.5, colour = "#546E7A"),
    axis.text.y        = element_text(size = 11,   colour = "#0D1117"),
    axis.ticks.y       = element_blank(),
    axis.ticks.x       = element_line(colour = "grey80", linewidth = 0.3),
    # --- Grid: horizontal lines only ---
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.30),
    # --- Panel borders ---
    panel.border       = element_rect(colour = "#555555", fill = NA, linewidth = 0.65),
    panel.background   = element_rect(fill = "white", colour = NA),
    panel.spacing.x    = unit(0.90, "lines"),
    panel.spacing.y    = unit(0.90, "lines"),
    # --- Legend ---
    legend.position    = "right",
    legend.box         = "vertical",
    legend.title       = element_text(face = "bold", size = 12, colour = "#1C2B3A"),
    legend.text        = element_text(size = 11,  colour = "#1C2B3A"),
    legend.key.size    = unit(1.0, "lines"),
    legend.background  = element_rect(fill = "white", colour = NA),
    legend.margin      = margin(0, 4, 0, 4),
    legend.spacing.y   = unit(0.6, "lines")
  )

ggsave(
  file.path(out_dir, "ora_hallmarks_pathways_by_regulon_canva_slide.png"),
  p_ora_canva,
  width  = 13.33,
  height = 7.5,
  dpi    = 300,
  bg     = "white"
)

message("\nCanva slide version saved: ora_hallmarks_pathways_by_regulon_canva_slide.png  [13.33 x 7.5 in, 16:9]")
