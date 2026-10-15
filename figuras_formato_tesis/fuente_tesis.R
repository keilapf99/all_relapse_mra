# Fuente de las figuras con formato de la tesis
#
# Las figuras originales (all_relapse_masters_degree) se hicieron en una Mac,
# donde "sans" es Helvetica/Arial. En este servidor "sans" es DejaVu Sans, que
# es más ancha y recorta etiquetas. Aquí se registra FreeSans (clon libre de
# Helvetica) como fuente por defecto ("" y "sans") y bajo "Arial" y "Helvetica" para los PNG
# (ggsave usa ragg). Si FreeSans no existe, no se cambia nada.

free_dir <- "/usr/share/fonts/truetype/freefont"
if (requireNamespace("systemfonts", quietly = TRUE) &&
    file.exists(file.path(free_dir, "FreeSans.ttf"))) {
  for (familia in c("", "sans", "Arial", "Helvetica")) {
    systemfonts::register_font(
      name       = familia,
      plain      = file.path(free_dir, "FreeSans.ttf"),
      bold       = file.path(free_dir, "FreeSansBold.ttf"),
      italic     = file.path(free_dir, "FreeSansOblique.ttf"),
      bolditalic = file.path(free_dir, "FreeSansBoldOblique.ttf")
    )
  }
}
rm(free_dir)
