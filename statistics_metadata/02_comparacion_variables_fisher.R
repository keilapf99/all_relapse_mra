# Comparación de variables clínicas entre pacientes con y sin recaída (LLA)
# con prueba exacta de Fisher para cada variable categórica.
#
# Salidas (en esta misma carpeta):
#   - comparacion_variables_LLA_fisher_A4.png / .pdf : figura con p-valores
#   - pruebas_fisher_variables_LLA.tsv               : tabla de pruebas

library(tidyverse)
library(patchwork)

base_dir <- "."  # raíz del proyecto: ejecutar desde leukemia_maestria/
out_dir  <- file.path(base_dir, "statistics_metadata")

# ---- Datos ------------------------------------------------------------------

meta <- read.delim(file.path(base_dir, "preprocess/metadata_NS500560.tsv"),
                   colClasses = "character")

num <- function(x) as.numeric(sub(",", ".", x))

# El orden de los niveles es el orden de la leyenda (de arriba hacia abajo).
# El OR se reporta para la 2a categoría frente a la 1a (CON vs. SIN recaída).
meta <- meta %>%
  mutate(
    grupo = factor(if_else(relapse == "si", "CON RECAÍDA", "SIN RECAÍDA"),
                   levels = c("CON RECAÍDA", "SIN RECAÍDA")),
    edad_anios = num(age_months) / 12,
    cl_uL = num(cl),  # células/µL
    sexo = factor(if_else(sex == "1", "Masculino", "Femenino"),
                  levels = c("Femenino", "Masculino")),
    edad = factor(if_else(edad_anios >= 10, "≥10", "1-9.99"),
                  levels = c("1-9.99", "≥10")),
    leucocitos = factor(if_else(cl_uL >= 50000, "≥50,000/µL", "<50,000/µL"),
                        levels = c("<50,000/µL", "≥50,000/µL")),
    # Agrupación según St. Jude: riesgo alto = BCR::ABL1+;
    # riesgo habitual = ETV6::RUNX1+, TCF3::PBX1+ y sin fusiones detectadas.
    fusion = factor(if_else(rearrange == "BCR::ABL1",
                            "Riesgo alto", "Riesgo habitual"),
                    levels = c("Riesgo habitual", "Riesgo alto")),
    riesgo = factor(if_else(risk == "RA", "Riesgo alto", "Riesgo habitual"),
                    levels = c("Riesgo habitual", "Riesgo alto"))
  )

variables <- c(
  sexo       = "SEXO",
  edad       = "EDAD (AÑOS)",
  leucocitos = "CL EN SANGRE PERIFÉRICA AL DX",
  fusion     = "FUSIONES GÉNICAS",
  riesgo     = "CLASIFICACIÓN DE RIESGO"
)

# ---- Pruebas exactas de Fisher ---------------------------------------------

formatea_p <- function(p) if_else(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))

prueba_fisher <- function(var) {
  tab <- table(meta[[var]], meta$grupo)
  tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  ft  <- fisher.test(tab)
  es_2x2 <- nrow(tab) == 2
  # Tabla invertida para que el OR sea de la 2a categoría frente a la 1a
  ft_or <- if (es_2x2) fisher.test(tab[2:1, ]) else NULL
  tibble(
    variable   = var,
    tabla      = paste(apply(tab, 1, paste, collapse = "/"),
                       rownames(tab), sep = " ", collapse = "; "),
    dimension  = paste(dim(tab), collapse = "x"),
    p_valor    = ft$p.value,
    comparacion_or = if (es_2x2) paste(rownames(tab)[2], "vs.", rownames(tab)[1]) else NA_character_,
    odds_ratio = if (es_2x2) unname(ft_or$estimate) else NA_real_,
    ic95_inf   = if (es_2x2) ft_or$conf.int[1] else NA_real_,
    ic95_sup   = if (es_2x2) ft_or$conf.int[2] else NA_real_
  )
}

resultados <- map_dfr(names(variables), prueba_fisher) %>%
  mutate(significativo = p_valor < 0.05)

# Complemento: las variables continuas sin categorizar (Mann-Whitney)
continuas <- map_dfr(c("edad_anios", "cl_uL"), \(v) {
  wt <- wilcox.test(meta[[v]] ~ meta$grupo, exact = FALSE)
  tibble(variable = v, prueba = "Mann-Whitney (continua)", p_valor = wt$p.value)
})

write_tsv(bind_rows(mutate(resultados, prueba = "Fisher exacta"), continuas),
          file.path(out_dir, "pruebas_fisher_variables_LLA.tsv"))

print(resultados)
print(continuas)

# ---- Gráficos ---------------------------------------------------------------

paleta <- c(
  "azul_claro" = "#9ECAE1",
  "azul_medio" = "#4292C6",
  "gris"       = "#6B7280",
  "azul_oscuro" = "#08519C"
)

# Texto blanco sobre colores oscuros y negro sobre claros
color_texto <- function(hex) {
  rgb <- col2rgb(hex) / 255
  if_else(colSums(rgb * c(0.299, 0.587, 0.114)) < 0.5, "white", "black")
}
gris_texto <- "#555555"

tema <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5, colour = gris_texto,
                                     margin = margin(b = 4)),
        axis.title.y = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "#333333"),
        legend.key.size = unit(0.45, "cm"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

# Colores por categoría, en el orden de la leyenda; "SIN FUSIONES" va en gris
colores_para <- function(niveles) {
  if ("SIN FUSIONES" %in% niveles) {
    cols <- paleta[c("azul_claro", "azul_medio", "azul_oscuro", "gris")]
  } else {
    cols <- paleta[seq_along(niveles)]
  }
  setNames(unname(cols), niveles)
}

barras <- function(df, var, titulo) {
  colores <- colores_para(levels(df[[var]]))
  conteos <- df %>%
    count(grupo, valor = .data[[var]]) %>%
    group_by(grupo) %>%
    mutate(etiqueta = sprintf("%d (%s%%)", n, round(100 * n / sum(n), 1))) %>%
    ungroup() %>%
    mutate(texto = color_texto(colores[as.character(valor)]))

  ggplot(conteos, aes(grupo, n, fill = valor)) +
    geom_col(width = 0.9) +
    geom_text(aes(label = etiqueta, colour = texto),
              position = position_stack(vjust = 0.5), size = 3.2) +
    scale_colour_identity() +
    scale_fill_manual(values = colores) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(title = titulo, y = "Número de pacientes") +
    tema
}

p_muestra <- barras(mutate(meta, muestra = grupo), "muestra", "MUESTRA POBLACIONAL")
paneles <- imap(variables, \(titulo, var) barras(meta, var, titulo))

# ---- Tabla de resultados de Fisher bajo la figura ----------------------------

nombres_tabla <- c(
  sexo       = "Sexo",
  edad       = "Edad (años)",
  leucocitos = "CL en sangre periférica al Dx",
  fusion     = "Fusiones génicas",
  riesgo     = "Clasificación de riesgo"
)

tabla_fisher <- resultados %>%
  transmute(
    Variable = nombres_tabla[variable],
    `Categoría evaluada` = coalesce(str_remove(comparacion_or, " vs\\. .*"),
                                    "Distribución global (4×2)"),
    `Referencia` = coalesce(str_remove(comparacion_or, ".* vs\\. "), "—"),
    `OR (IC 95%)` = if_else(is.na(odds_ratio), "—",
                            sprintf("%.2f (%.2f–%.2f)", odds_ratio, ic95_inf, ic95_sup)),
    `Valor de p` = sub("p = |p ", "", formatea_p(p_valor))
  )

tema_tabla <- gridExtra::ttheme_minimal(
  base_size = 11,
  core = list(fg_params = list(hjust = 0, x = 0.04),
              bg_params = list(fill = c("#F3F6FA", "white"), col = NA)),
  colhead = list(fg_params = list(hjust = 0, x = 0.04, fontface = "bold"),
                 bg_params = list(fill = paleta[["azul_claro"]], col = NA))
)

grob_tabla <- gridExtra::tableGrob(tabla_fisher, rows = NULL, theme = tema_tabla)
# Ancho de la tabla: 95% de la página, repartido entre las columnas
grob_tabla$widths <- unit(c(0.26, 0.22, 0.20, 0.17, 0.10), "npc")

figura <- wrap_plots(c(list(p_muestra), paneles), ncol = 2) /
  wrap_elements(full = grob_tabla) +
  plot_layout(heights = c(3, 0.6)) +
  plot_annotation(
    caption = paste0(
      "Prueba exacta de Fisher bilateral; significativo si p < 0.05.\n",
      "Fusiones génicas (St. Jude): riesgo alto = BCR::ABL1+; ",
      "riesgo habitual = ETV6::RUNX1+, TCF3::PBX1+ o sin fusiones detectadas.\n",
      "OR: razón de momios de recaída de la categoría evaluada frente a la de referencia ",
      "(la primera de cada leyenda); ",
      "IC 95%: intervalo de confianza al 95%."),
    theme = theme(plot.caption = element_text(size = 8.5, colour = gris_texto, hjust = 0))
  )

ggsave(file.path(out_dir, "comparacion_variables_LLA_fisher_A4.png"), figura,
       width = 8.27 * 1.25, height = 11.69 * 1.25, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "comparacion_variables_LLA_fisher_A4.pdf"), figura,
       width = 8.27 * 1.25, height = 11.69 * 1.25, device = cairo_pdf)
