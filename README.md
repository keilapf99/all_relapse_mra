# Recaída en LLA pediátrica: análisis transcriptómico

Tesis de maestría. Autora: Keila Adonai Pérez-Fuentes.

Comparación de pacientes con LLA **con recaída** vs. **sin recaída** a partir de
datos de RNA-seq (35 pacientes secuenciados en el instrumento NS500560).

## Cómo correr la tubería

Desde la carpeta raíz del proyecto (`leukemia_maestria/`):

```bash
Rscript 00_run_pipeline.R            # todos los pasos, del 01 al 10 (~5 min)
Rscript 00_run_pipeline.R 04 05a 05b # solo algunos pasos
```

Todas las rutas son relativas a la raíz del proyecto, así que la carpeta puede
moverse de servidor o de computadora sin cambiar los scripts. Para correr un
script suelto en RStudio, primero fija el directorio de trabajo en la raíz
(`setwd(".../leukemia_maestria")`). La salida de cada paso queda en `logs/`.

## Orden de ejecución

Los scripts están numerados en el orden en que deben correrse. Un script
depende de los archivos que generan los anteriores.

| # | Script | Qué hace | Lee | Genera |
|---|---|---|---|---|
| 01 | `preprocess/01_first_explo_pca.R` | PCA exploratorio para evaluar efectos de lote (instrumento, flowcell). Se conservan solo las muestras del instrumento NS500560. | `raw_data/raw_matrix.tsv`, `raw_data/raw_metadata.tsv`, `raw_data/batches_metadata.tsv` | `preprocess/metadata_NS500560.tsv`, PCA (TMM) |
| 02 | `statistics_metadata/02_comparacion_variables_fisher.R` | Asociación de las variables clínicas de riesgo (NCI / St. Jude) con la recaída: prueba exacta de Fisher (y Mann-Whitney para edad y CL continuas). | `preprocess/metadata_NS500560.tsv` | `comparacion_variables_LLA_fisher_A4.png/.pdf`, `pruebas_fisher_variables_LLA.tsv` |
| 03 | `preprocess/03_preprocessing_vst_pca.R` | Filtrado de genes poco expresados (≥10 lecturas en ≥15 muestras), normalización VST, promedio de símbolos duplicados y separación en matrices por grupo (para el análisis de redes). | `raw_matrix.tsv`, `metadata_NS500560.tsv` | `preprocess/counts_vst.tsv`, `matrix_relapse_vst.tsv`, `matrix_control_vst.tsv`, PCA (VST) |
| 04 | `DEA/04_DEA_relapse.R` | Expresión diferencial con DESeq2 (`~ relapse`, contraste si vs. no), mismo filtro que el paso 03 (≥10 lecturas en ≥15 muestras). | `raw_data/raw_matrix.tsv`, `metadata_NS500560.tsv` | `DEA/DEA_totales_relapse.tsv` (firma completa), `DEA/DEGS_significant_relapse.tsv` (\|log2FC\| ≥ 1, FDR < 0.05), volcano plot |
| 05a | `enrichment_analysis/05a_gsea_msigdb_hallmarks_relapse.R` | GSEA (clusterProfiler) de la firma completa, ordenada por el `stat` de DESeq2, contra la colección Hallmark (H) de MSigDB. | `DEA/DEA_totales_relapse.tsv` | `enrichment_analysis/gsea_msigdb_hallmarks/` |
| 05b | `enrichment_analysis/05b_gsea_cancer_hallmarks_menyhart_relapse.R` | GSEA (clusterProfiler) de la misma firma contra los Cancer Hallmarks de Menyhart et al. (2025), conjuntos *core* (principal) e *integrated*. | `DEA/DEA_totales_relapse.tsv` | `enrichment_analysis/gsea_cancer_hallmarks_menyhart/` |
| 06 | `corto_mra/06_mra_relapse.R` | Análisis de reguladores maestros con corto (redes por grupo + firma de expresión diferencial). Usa `mra_mod.R` vía `source()`. | `matrix_relapse_vst.tsv`, `matrix_control_vst.tsv`, `raw_data/tfs_motif_yestf.txt`, `DEA/DEA_totales_relapse.tsv` | `corto_mra/mra_relapse.RData`, `cortoMRS.pdf`, `mra_relapse_todos_TFs.tsv`, `mra_relapse_TFs_FDR005.tsv` |
| 07 | `enrichment_analysis/07_ora_regulones_top10.R` | ORA (clusterProfiler::enricher) de los regulones de los 10 reguladores maestros con mayor \|NES\| (FDR < 0.05), contra MSigDB H y Menyhart *core* e *integrated*. Universo: genes de la matriz VST. | `corto_mra/mra_relapse.RData`, `preprocess/matrix_relapse_vst.tsv`, conjuntos de genes de 05a y 05b | `enrichment_analysis/ora_regulones_top10/` |
| 08 | `enrichment_analysis/08_red_regulones_hallmarks.R` | ORA de los regulones de los 157 reguladores maestros significativos contra MSigDB H y Menyhart *core*, con BH global por colección; red bipartita regulador–hallmark y su proyección hallmark–hallmark (Jaccard de reguladores compartidos). | Mismas que el 07 | `enrichment_analysis/red_regulones_hallmarks/` (tablas de nodos y aristas para Cytoscape, figuras) |
| 09 | `enrichment_analysis/09_leading_edge_reguladores.R` | Para cada vía enriquecida en 05a y 05b (*core*): cobertura de sus genes leading edge por cada regulón y sobrerrepresentación del leading edge en los regulones de los 157 reguladores maestros (hipergeométrica, BH global). | `corto_mra/mra_relapse.RData`, `preprocess/matrix_relapse_vst.tsv`, resultados de 05a y 05b | `enrichment_analysis/leading_edge_reguladores/` (tablas) |
| 10 | `DEA/10_heatmap_agrupamiento.R` | Agrupamiento jerárquico no supervisado de los pacientes (z-score de la VST, distancia euclidiana, Ward.D2): A. con los 237 DEGs (descriptivo; sesgado a favor de separar los grupos porque los genes se eligieron con esa etiqueta) y B. con los 237 genes de mayor varianza (control sin usar la etiqueta). Concordancia con recaída: Fisher y ARI. | `preprocess/counts_vst.tsv`, `DEA/DEGS_significant_relapse.tsv`, `metadata_NS500560.tsv` | `DEA/heatmap_agrupamiento/` |

```
raw_data/
   │
   ▼
01 PCA exploratorio ──► metadata_NS500560.tsv
   │
   ├──► 02 Variables clínicas (Fisher)
   ├──► 03 VST ─────────────────────────────────┐
   └──► 04 DEA (DESeq2) ──┬──► 05a GSEA MSigDB  │
                          ├──► 05b GSEA Menyhart│
                          └──► 06 MRA ◄─────────┘
                                 │
                                 ▼
                          10 Heatmap / agrupamiento de pacientes (usa 03 y 04)
                          07 ORA de los regulones del top 10
                          08 ORA de los 157 regulones → red regulador–hallmark
                          09 Leading edge de las vías (05a, 05b) vs. regulones
                             (07, 08 y 09 usan los resultados de 05a y 05b)
```

Los pasos 02, 03 y 04 son independientes entre sí; solo requieren el 01.
El 06 requiere el 03 y el 04; el 07, 08 y 09 requieren el 03, 05a, 05b y 06.

## Scripts auxiliares (no forman parte del flujo numerado)

- `corto_mra/mra_mod.R`: versión modificada de `corto::mra` que acepta un vector
  como firma; la carga `06_mra_relapse.R`. No cambiar su nombre.
- `enrichment_analysis/gsea_hallmarks.R`: versión original del GSEA con MSigDB;
  reemplazada por el paso 05a. Sus figuras (`enrichment_analysis/GSEA_Hallmarks_*.png`,
  `barplot_GSEA_A4.png`, `DEA_relapse_vp.png`) son de la firma anterior.
- `enrichment_analysis/gsea_hallmarks_menyhart.R`: plantilla genérica en la que
  se basa el paso 05b.

## Números de referencia (corrida del 13 de octubre de 2026)

| Paso | Resultado |
|---|---|
| 01 | 35 pacientes del NS500560 (20 con recaída, 15 sin recaída). La muestra 325 está en la metadata pero no en `raw_matrix.tsv`. |
| 03 | 17,044 genes pasan el filtro → 14,221 con símbolo génico → 14,208 símbolos únicos (matrices VST: 14,208 × 35). |
| 04 | Firma completa: 14,221 genes; DEGs (\|log2FC\| ≥ 1, FDR < 0.05): 237 (204 ↑, 33 ↓). |
| 05a | MSigDB H (versión 2026.1.Hs): 26/50 hallmarks con FDR < 0.05 (25 enriquecidos en recaída, 1 sin recaída: WNT/β-catenin). |
| 05b | Core: 3/10 hallmarks con FDR < 0.05 (Tissue invasion and metastasis, Resisting cell death, Genome instability), todos enriquecidos en recaída. Integrated: 9/10. |
| 06 | 243 TFs evaluados; 157 con FDR < 0.05. |
| 07 | 8 de los 10 regulones con al menos un conjunto sobrerrepresentado (FDR < 0.05) en alguna colección; ZNF444 y KLF16 sin ninguno. |
| 08 | 27 de 157 reguladores con al menos un hallmark (FDR global < 0.05); 37 aristas, 13 hallmarks. Módulos: inflamatorio (CEBPB/D/E, MAFB, LTF, NFKB2, MXD1) y metabólico (oxidative phosphorylation / reprogramming energy metabolism; HMG20B, ZNF581, ZNF282, BPTF, NR2F6…). |
| 09 | Ningún regulón contiene el leading edge completo de una vía (máximo 26%: CEBPD en IL6/JAK/STAT3); la unión de los 157 regulones cubre en mediana el 28%. 13 pares vía–regulador con FDR < 0.05, casi todos del módulo inflamatorio (CEBPD, CEBPB, CEBPE, LTF, MAFB, CREB5); además ZNF680 → mitotic spindle y VDR → coagulation. |
| 10 | Ni con los DEGs se separan recaída y no recaída (ARI = 0.04): aparece un subgrupo de 9 pacientes (8 con recaída, todos de riesgo alto) con genes de megacariocito, eritroides y mieloides elevados y CD19/PAX5 disminuidos, compatible con menor proporción de blastos en la muestra. Con los genes de mayor varianza: ARI = −0.02 (Fisher p = 0.21). |

## Figuras con el formato de la tesis (`figuras_formato_tesis/`)

Las figuras de la tesis se habían generado en `all_relapse_masters_degree/`.
Estos scripts reproducen ese mismo formato (mismo código de ggplot, colores,
títulos y tamaños) pero con los resultados de esta tubería. No recalculan
ningún análisis: solo leen lo que guardan los pasos 01–07.

```bash
Rscript figuras_formato_tesis/00_run_figuras.R        # todas (~1 min)
Rscript figuras_formato_tesis/00_run_figuras.R F3 F6  # solo algunas
```

| Script | Basado en (all_relapse_masters_degree) | Usa | Genera (en `figuras_formato_tesis/<script>/`) |
|---|---|---|---|
| `F1_estadistica_poblacion.R` | `0_1_statistics_poblation/statistics.r` | metadata del paso 03 | `comparacion_variables_LLA_A4.png` (azul, 2 columnas, como en la tesis), `_SLIDES.png` (morado, 3 columnas) |
| `F2_filtrado_pca_vst.R` | `1_preprocessing/prep_vst.r` (figuras) | pasos 01 y 03 | `filtering_comparison.png`, `PCA_vst_relapse.png`, `combined_pca.png` |
| `F3_volcano_dea.R` | `2_DEA/2_DEA_genes_of_interest.r` y `2_DEA.r` | paso 04 | `DEA_relapse_vp.png` (genes por proceso biológico), `_top25.png`, `_sinlabels.png` |
| `F4_mra_corto.R` | `3_MRA_corto/mra_corto_adapted.r` y `mra_tfs_interes.r` | paso 06 | `cortoMRStop1.pdf`, `cortoMRStop2.pdf`, `mra_tfs_interes.pdf` |
| `F5_gsea_msigdb_hallmarks.R` | `4_GSEA_ORA/gsea_hallmarks_plots*.r`, `gsea_hallmarks_top5_publication*.r` | paso 05a | `Hallmark_relapse_mechanisms{,_centered}.png`, `top5_hallmarks_lolliplot_publication{,_es}.png` |
| `F6_ora_regulones_hallmarks.R` | `4_GSEA_ORA/ora_top10_master_regulators_*.r` | pasos 03 y 07 (solo MSigDB H) | `ora_hallmarks_pathways_by_regulon{,_3col,_canva_slide}.png` y tablas `.csv` |
| `F7_gsea_lollipop_top10_top5.R` | objetos `lolli_top10` / `lolli_top5` de `4_GSEA_ORA/.RData` | paso 05a | `top10_lolliplot_GSEA.png`, `top5_gsea_nes.png` |

Auxiliares: `fuente_tesis.R` (usa FreeSans, clon de Helvetica, para que el texto
se vea como en las figuras originales hechas en Mac), `ora_publication_plots_tesis.R`
(copia del código de figuras del ORA de all_relapse) y
`datos/degs_procesobiologico.tsv` (genes etiquetados en el volcano).

Diferencias de contenido con las figuras de la tesis (vienen de los análisis,
no del formato):

- **ORA (F6):** con el universo correcto (14,208 genes de la matriz VST, en vez
  de los 39,301 símbolos de la matriz cruda) quedan 12 asociaciones
  regulador–hallmark significativas en 7 reguladores, en lugar de 101 en 9.
  Cada regulador conserva el color que tenía en la tesis.
- **GSEA (F5, F7):** los NES cambian en el segundo decimal. El tamaño de los
  puntos (−log10 FDR) llega a ~20 porque el paso 05a calcula p exactas
  (`eps = 0`); en la tesis el mínimo estaba truncado en 1e−10.
- **MRA (F4):** mismos 10 reguladores principales, con el orden ligeramente
  distinto (KLF16 pasa al 3.º lugar, ZNF581 al 5.º, CEBPD al 8.º).
- **Estadística, PCA y volcano (F1–F3):** mismos números que en la tesis.
- No se reprodujeron `barplot_GSEA(_A4).png`, `Hallmarks_relapse_mechanisms.png`
  ni `recaida_gsea.png` de all_relapse: eran borradores cuyo código ya no existe.

## Notas

- Semilla fija (`set.seed(123)`) en los pasos con componentes aleatorios
  (GSEA y MRA); dos corridas consecutivas dan resultados idénticos.
- GSEA con el número de permutaciones por defecto de clusterProfiler (1000) y
  semilla fija. Genome instability (core, FDR = 0.045) es un resultado frontera:
  con otras semillas su FDR varió entre 0.011 y 0.115 (significativo en ~30% de
  20 semillas) y con 100,000 permutaciones queda en ≈ 0.066.
- `clusterProfiler`/`enrichplot` dependen de `igraph`, que requiere la librería
  del sistema `libglpk`, ausente en este servidor. Se instaló en un entorno conda
  del usuario (`/scratch/home/aperez/.conda_envs/glpk`) y los scripts 05a, 05b y
  07 la cargan con `dyn.load()` si existe. En otra computadora con `libglpk`
  instalada esa línea simplemente no hace nada.
- Las listas de genes de MSigDB y Menyhart se guardan localmente en
  `gene_sets/` dentro de cada carpeta de resultados, para que el análisis no
  cambie si las bases publican versiones nuevas.
- `corto_mra/mra_relapse_030426.RData` y `mra_relapse_110126.RData` son
  resultados de versiones anteriores (filtro distinto); los vigentes son los de
  `mra_relapse.RData`.
- `corto_mra/results.rds` no lo genera ningún script de esta tubería (contiene
  metadatos de casos TARGET).
- Respaldo del proyecto antes de esta reestructuración:
  `../leukemia_maestria_respaldo_2026-10-12/`.
