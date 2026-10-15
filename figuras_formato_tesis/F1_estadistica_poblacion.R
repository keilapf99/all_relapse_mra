# Title : Clinical statistics ALL relapse cohort (formato de figuras de la tesis)
# Author: Keila Adonai Pérez-Fuentes
# Basado en: all_relapse_masters_degree/0_1_statistics_poblation/statistics.r
# Objective: Mismo gráfico que el original, con la metadata de leukemia_maestria.
#
# Ejecutar desde la carpeta raíz del proyecto (leukemia_maestria/):
#   Rscript figuras_formato_tesis/F1_estadistica_poblacion.R
# Entradas: preprocess/metadata_control.tsv, preprocess/metadata_relapse.tsv (paso 03)
# Salidas : figuras_formato_tesis/F1_estadistica_poblacion/comparacion_variables_LLA_{A4,SLIDES}.png

# =============================================================================
# 0. PACKAGES
# =============================================================================

library(tidyverse)
library(patchwork)
library(dplyr)

source("figuras_formato_tesis/fuente_tesis.R")

# =============================================================================
# 1. DATA
# =============================================================================
out_dir <- "figuras_formato_tesis/F1_estadistica_poblacion"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

metadata_control <- read.table(
  file = "preprocess/metadata_control.tsv",
  header = T,
  sep = "\t")

metadata_relapse <- read.table(
  file = "preprocess/metadata_relapse.tsv",
  header = T,
  sep = "\t")

# Join metadatas
metadata <- bind_rows(
  metadata_control,
  metadata_relapse
)

# Grupos de riesgo por edad y cuenta leucocitaria (en all_relapse se agregaban
# en prep_vst.r; la metadata de leukemia_maestria no los trae)
metadata$age_months <- as.numeric(
  gsub(",", ".", metadata$age_months)
)
metadata$age_group_years <- ifelse(
  metadata$age_months >= 120,
  "≥10",
  "1–9.99"
)
metadata$cl_risk <- ifelse(
  metadata$cl >= 50000,
  "≥50000cells/μL",
  "<50000cells/μL"
)

# =============================================================================
# 3. Cleaning and clinical variables
# =============================================================================

metadata <- metadata %>%
  
  mutate(
    
    # -------------------------------------------------------------------------
    # SEXO
    # -------------------------------------------------------------------------
    sex = recode(
      as.character(sex),
      "1" = "Masculino",
      "2" = "Femenino"
    ),
    
    # -------------------------------------------------------------------------
    # FUSIONES
    # -------------------------------------------------------------------------
    rearrange = case_when(
      rearrange == "ND" ~ "SIN FUSIONES",
      TRUE ~ rearrange
    ),
    
    # -------------------------------------------------------------------------
    # RIESGO
    # -------------------------------------------------------------------------
    risk = recode(
      risk,
      "RA" = "Riesgo alto",
      "RH" = "Riesgo habitual"
    ),
    
    # -------------------------------------------------------------------------
    # GRUPOS DE RECAÍDA
    # -------------------------------------------------------------------------
    Grupo = case_when(
      relapse == "si" ~ "CON RECAÍDA",
      relapse == "no" ~ "SIN RECAÍDA"
    )
  )

# =============================================================================
# 4. Factor orders
# =============================================================================
# Esto controla:
# - orden de barras
# - orden de colores
# - orden en leyendas

metadata <- metadata %>%
  
  mutate(
    
    # -------------------------------------------------------------------------
    # EDAD
    # -------------------------------------------------------------------------
    age_group_years = factor(
      age_group_years,
      levels = c(
        "<1",
        "1–9.99",
        "≥10"
      )
    ),
    
    # -------------------------------------------------------------------------
    # CUENTA LEUCOCITARIA
    # -------------------------------------------------------------------------
    
    cl_risk = case_when(
      cl_risk == "<50000cells/μL" ~ "<50,000/μL",
      cl_risk == "≥50000cells/μL" ~ "≥50,000/μL",
      TRUE ~ cl_risk
    ),
    
    cl_risk = factor(
      cl_risk,
      levels = c(
        "<50,000/μL",
        "≥50,000/μL"
      )
    ),
    
    # -------------------------------------------------------------------------
    # FUSIONES
    # SIN FUSIONES queda hasta abajo
    # -------------------------------------------------------------------------
    rearrange = factor(
      rearrange,
      levels = c(
        "ETV6::RUNX1",
        "TCF3::PBX1",
        "BCR::ABL1",
        "SIN FUSIONES"
      )
    )
  )

# =============================================================================
# 5. Table functions
# =============================================================================
# Automatiza los conteos

make_table <- function(data, var, varname) {
  
  data %>%
    
    group_by(
      Grupo,
      Nivel = .data[[var]]
    ) %>%
    
    summarise(
      n = n(),
      .groups = "drop"
    ) %>%
    
    mutate(
      Variable = varname
    )
}

# =============================================================================
# 6. Plots dataframes
# =============================================================================

df <- bind_rows(
  
  # ---------------------------------------------------------------------------
  # COHORTE TOTAL
  # ---------------------------------------------------------------------------
  
  metadata %>%
    
    dplyr::count(Grupo) %>%
    
    dplyr::mutate(
      Variable = "Muestra poblacional",
      Nivel = Grupo
    ) %>%
    
    dplyr::select(
      Variable,
      Nivel,
      Grupo,
      n
    ),
  
  # ---------------------------------------------------------------------------
  # SEXO
  # ---------------------------------------------------------------------------
  
  make_table(
    metadata,
    "sex",
    "Sexo"
  ),
  
  # ---------------------------------------------------------------------------
  # EDAD
  # ---------------------------------------------------------------------------
  
  make_table(
    metadata,
    "age_group_years",
    "Edad"
  ),
  
  # ---------------------------------------------------------------------------
  # CUENTA LEUCOCITARIA
  # ---------------------------------------------------------------------------
  
  make_table(
    metadata,
    "cl_risk",
    "CL en sangre periférica al diagnóstico"
  ),
  
  # ---------------------------------------------------------------------------
  # FUSIONES
  # ---------------------------------------------------------------------------
  
  make_table(
    metadata,
    "rearrange",
    "Fusiones génicas"
  ),
  
  # ---------------------------------------------------------------------------
  # RIESGO
  # ---------------------------------------------------------------------------
  
  make_table(
    metadata,
    "risk",
    "Clasificación de riesgo"
  )
)

# =============================================================================
# 7. Plots functions
# =============================================================================

plot_var <- function(data, varname, paleta) {
  
  # ---------------------------------------------------------------------------
  # FILTRAR VARIABLE
  # ---------------------------------------------------------------------------
  
  sub <- data %>%
    
    filter(
      Variable == varname
    ) %>%
    
    group_by(Grupo) %>%
    
    mutate(
      total = sum(n),
      pct = round(100 * n / total, 1)
    )
  
  # ---------------------------------------------------------------------------
  # ORDEN ESPECÍFICO PARA ALGUNAS VARIABLES
  # ---------------------------------------------------------------------------
  
  sub <- sub %>%
    
    mutate(
      Nivel = case_when(
        
        # -------------------------------------------------
        # CUENTA LEUCOCITARIA
        # -------------------------------------------------
        
        varname == "Cuenta leucocitaria en sangre periférica" ~ factor(
          Nivel,
          levels = c(
            "<50000cells/μL",
            "≥50000cells/μL"
          )
        ),
        
        # -------------------------------------------------
        # FUSIONES
        # -------------------------------------------------
        
        varname == "Fusiones génicas" ~ factor(
          Nivel,
          levels = c(
            "ETV6::RUNX1",
            "TCF3::PBX1",
            "BCR::ABL1",
            "SIN FUSIONES"
          )
        ),
        
        TRUE ~ factor(Nivel)
      )
    )
  
  # ---------------------------------------------------------------------------
  # TÍTULOS
  # ---------------------------------------------------------------------------
  
  title_text <- if (
    varname == "CL en sangre periférica al diagnóstico"
  ) {
    
    "CL EN SANGRE PERIFÉRICA AL DX"
    
  } else {
    
    toupper(varname)
  }
  
  # ---------------------------------------------------------------------------
  # PLOT
  # ---------------------------------------------------------------------------
  
  ggplot(
    sub,
    aes(
      x = Grupo,
      y = n,
      fill = Nivel
    )
  ) +
    
    # -------------------------------------------------------------------------
  # BARRAS
  # -------------------------------------------------------------------------
  
  geom_bar(
    stat = "identity"
  ) +
    
    # -------------------------------------------------------------------------
  # LABELS
  # -------------------------------------------------------------------------
  
  geom_text(
    aes(
      label = paste0(
        n,
        " (",
        pct,
        "%)"
      )
    ),
    
    position = position_stack(
      vjust = 0.5
    ),
    
    size = 3,
    color = "black"
  ) +
    
    # -------------------------------------------------------------------------
  # COLORES
  # -------------------------------------------------------------------------
  
  # Paleta: azul en la versión A4 (la de la tesis), morado en la de diapositivas
  scale_fill_manual(
    values = paleta
  ) +
    
    # -------------------------------------------------------------------------
  # LABELS GENERALES
  # -------------------------------------------------------------------------
  
  labs(
    title = title_text,
    x = "",
    y = "Número de pacientes"
  ) +
    
    # -------------------------------------------------------------------------
  # THEME
  # -------------------------------------------------------------------------
  
  theme_minimal() +
    
    theme(
      
      plot.title = element_text(
        face = "bold",
        size = 12
      ),
      
      axis.text = element_text(
        color = "black"
      ),
      
      legend.title = element_blank()
    )
}

# =============================================================================
# 8. Plots
# =============================================================================

vars <- unique(df$Variable)

paleta_azul <- c(
  "#9ECAE1", # azul claro
  "#4292C6", # azul medio
  "#08519C", # azul oscuro
  "#6B7280"  # gris acero
)

paleta_morada <- c(
  "#D4B9DA", # lila claro
  "#C994C7", # morado medio
  "#9E4EA3", # púrpura fuerte
  "#6A1B6A"  # morado oscuro elegante
)

# =============================================================================
# 9. Join plots
# =============================================================================
# A4 vertical (figura de la tesis): azul, 2 columnas
big_plot_A4 <- wrap_plots(
  map(vars, ~plot_var(df, .x, paleta_azul)),
  ncol = 2
)

# Diapositivas: morado, 3 columnas
big_plot_slides <- wrap_plots(
  map(vars, ~plot_var(df, .x, paleta_morada)),
  ncol = 3
)

# =============================================================================
# 10. Save plots
# =============================================================================
# FORMATO A4 VERTICAL

ggsave(
  
  file.path(out_dir, "comparacion_variables_LLA_A4.png"),
  
  big_plot_A4,
  
  width = 8.27,
  height = 11.69,
  
  dpi = 300
)

ggsave(
  
  file.path(out_dir, "comparacion_variables_LLA_SLIDES.png"),
  
  big_plot_slides,
  
  width = 13.33,
  height = 7.5,
  
  dpi = 300
)
