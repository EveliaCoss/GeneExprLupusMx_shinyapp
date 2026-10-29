# Textos de la interfaz. Para agregar inglés: copiar TXT$es como TXT$en y traducir.
# tr("clave") busca en el idioma activo (APP_LANG) y, si falta, usa español.

LANG <- Sys.getenv("APP_LANG", "es")

URL_PKG      <- "https://github.com/EveliaCoss/GeneExprLupusMx"
URL_APP_REPO <- "https://github.com/EveliaCoss/GeneExprLupusMx_shinyapp"
URL_SITE     <- "https://eveliacoss.github.io/GeneExprLupusMx"
URL_ARTICULO <- "https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2026.1771959/full"
URL_DOI      <- "https://doi.org/10.3389/fimmu.2026.1771959"
URL_ZENODO   <- "https://doi.org/10.5281/zenodo.17742206"
URL_CODIGO   <- "https://github.com/NeuroGenomicsMX/SLE_Transcriptomic_analysis_tolerogenic"

# Figuras del sitio del paquete (se actualizan junto con el paquete)
img_url <- function(f) sprintf("%s/articles/img/%s", URL_SITE, f)

# Lecciones del paquete (sufijo _es / _en según idioma)
leccion_url <- function(id) sprintf("%s/articles/%s_%s.html", URL_SITE, id, LANG)

TXT <- list(es = list(
  titulo = "GeneExprLupusMx · Explorador",
  github_app = "Código de esta app en GitHub",
  tab_inicio = "1. Inicio",
  tab_participantes = "2. Participantes",
  tab_exploracion = "3. Exploración",
  tab_deg = "4. Expresión diferencial",
  tab_go = "5. Funciones biológicas",
  tab_descargas = "6. Descarga de los datos",
  tab_aprende = "7. Aprende a hacer estos análisis",

  # Pie de página (visible en todas las secciones)
  pie = paste0(
    "Datos: <a href='", URL_ARTICULO, "' target='_blank'>Hernández-Ledesma, Coss-Navarrete <i>et al.</i> (2026) ",
    "<i>Front. Immunol.</i> 17:1771959</a> y <a href='", URL_ZENODO, "' target='_blank'>Zenodo 10.5281/zenodo.17742206</a>, ",
    "distribuidos en el paquete de R <a href='", URL_PKG, "' target='_blank'><code>GeneExprLupusMx</code></a>. ",
    "¿Quieres hacer estas gráficas tú misma(o)? Sigue las <a href='", URL_SITE, "' target='_blank'>lecciones del paquete</a>."
  ),

  # ---- Inicio ----
  inicio_bienvenida = "Bienvenida(o)",
  inicio_intro = paste0(
    "<p>Esta aplicación te permite explorar de forma interactiva los resultados del artículo ",
    "<a href='", URL_ARTICULO, "' target='_blank'><i>Gene expression profiling of dendritic cell tolerance dysfunction ",
    "in women with systemic lupus erythematosus</i></a> (<i>Front. Immunol.</i> 2026), ",
    "usando los datos y resultados incluidos en el paquete de R <a href='", URL_PKG, "' target='_blank'><code>GeneExprLupusMx</code></a>.</p>",
    "<p>Recorre las secciones en orden: cada una se apoya en la anterior, igual que las lecciones del paquete.</p>"
  ),
  inicio_estudio_t = "Paso 1 · ¿Qué se estudió?",
  inicio_estudio = paste0(
    "<p>El <b>lupus eritematoso sistémico</b> (LES, o <i>SLE</i> en inglés) es una enfermedad autoinmune que afecta sobre todo a mujeres. ",
    "El estudio pregunta si las células que deberían <i>calmar</i> la respuesta inmune funcionan distinto en las pacientes.</p>",
    "<ul><li><b>33 mujeres mexicanas:</b> 10 controles sanas (<code>Ctrl</code>) y 23 pacientes con lupus (<code>SLE</code>), reclutadas entre junio de 2022 y febrero de 2023.</li>",
    "<li><b>3 tipos celulares por participante:</b> monocitos de sangre (<code>monocyte</code>), células dendríticas derivadas de monocitos ",
    "(<code>moDC</code>, cultivadas con GM-CSF + IL-4; activadoras) y células dendríticas tolerogénicas (<code>tolDC</code>, cultivadas con GM-CSF + IL-10; reguladoras).</li>",
    "<li><b>99 muestras</b> de RNA-seq (33 × 3).</li></ul>"
  ),
  inicio_datos_t = "Paso 2 · ¿De dónde vienen los datos?",
  inicio_datos = paste0(
    "<p>Los datos provienen del estudio publicado en <b>2026</b> en <i>Frontiers in Immunology</i> ",
    "(Hernández-Ledesma, Coss-Navarrete <i>et al.</i>, 17:1771959) y se obtuvieron de tres fuentes:</p>",
    "<ul><li><a href='", URL_ARTICULO, "' target='_blank'><b>Artículo original (2026)</b></a>: diseño del estudio, métodos, figuras y material suplementario.</li>",
    "<li><a href='", URL_ZENODO, "' target='_blank'><b>Zenodo, registro 17742206 (2025)</b></a>: análisis depositado por las autoras; ",
    "cuantificación con nf-core/rnaseq (<b>Salmon</b>) y análisis con <b>tximport</b> y <b>DESeq2</b>.</li>",
    "<li><a href='", URL_CODIGO, "' target='_blank'><i class='fab fa-github' role='img' aria-label='GitHub'></i> <b>SLE_Transcriptomic_analysis_tolerogenic</b></a>: código del análisis original y sus figuras.</li></ul>",
    "<p>El paquete <code>GeneExprLupusMx</code> reorganiza estas fuentes en objetos de R listos para usar ",
    "(cuentas en <code>filtered_txi</code>, expresión diferencial con DESeq2 en <code>de_results</code>, enriquecimiento con g:Profiler ",
    "<code>e113_eg59_p19</code> en <code>go_results</code>…). La app muestra esos <b>objetos ya calculados</b>; no recalcula los análisis. ",
    "Los scripts que los generan están en <a href='", URL_PKG, "/tree/main/data-raw' target='_blank'><code>data-raw/</code></a>.</p>",
    "<p>El paquete tiene un <b>propósito educativo y de investigación reproducible</b>, y promueve el <b>software libre y el código abierto</b>: ",
    "es software libre con licencia MIT, su código es público en <a href='", URL_PKG, "' target='_blank'>GitHub</a>, ",
    "está construido con paquetes de R de libre uso y sus datos provienen de fuentes de acceso abierto (CC BY).</p>",
    "<p><b>¿Quieres los datos?</b> Descarga cualquier tabla en CSV, con sus licencias y cómo citarla, en la sección ",
    "<a href='#' onclick=\"Shiny.setInputValue('ir_a', 'descargas', {priority: 'event'}); window.scrollTo(0, 0); return false;\">",
    "<b>6. Descarga de los datos</b></a>.</p>",
    "<p class='small text-muted'>Los datos crudos (FASTQ) están en el European Genome-phenome Archive (EGA) con acceso controlado; ",
    "los identificadores están en la Tabla suplementaria 2 del artículo.</p>"
  ),
  inicio_ruta_t = "Paso 3 · ¿Cómo usar esta app?",
  inicio_ruta = paste0(
    "<ol><li><b>Participantes:</b> conoce a las voluntarias (edad, tratamiento, actividad de la enfermedad).</li>",
    "<li><b>Exploración:</b> mira cómo se agrupan las 99 muestras y la expresión de cualquier gen.</li>",
    "<li><b>Expresión diferencial:</b> elige una comparación y descubre qué genes cambian.</li>",
    "<li><b>Funciones biológicas:</b> ¿en qué procesos participan esos genes?</li>",
    "<li><b>Descarga de los datos:</b> lleva cualquier tabla a Excel, Python o R.</li>",
    "<li><b>Aprende a hacer estos análisis:</b> reproduce cada gráfica con las lecciones del paquete.</li></ol>"
  ),
  inicio_datos_img = "Datos incluidos en el paquete: 10 controles y 23 pacientes, cada una con monocitos, moDC y tolDC (figura del paquete GeneExprLupusMx).",
  inicio_dif_img = paste0(
    "Diferenciación <i>in vitro</i>: a partir de los monocitos de cada voluntaria se generaron moDC y tolDC agregando citocinas al cultivo ",
    "(figura del paquete GeneExprLupusMx). Si necesitas más información sobre los métodos, consulta el ",
    "<a href='", URL_ARTICULO, "' target='_blank'>artículo original</a>."
  ),
  stat_deg_detalle = "monocito %d↑ %d↓ · moDC %d↑ %d↓ · tolDC %d↑ %d↓",
  sec_t = "Explora la app",
  sec_intro = "Haz clic en un panel para ir a la sección.",
  sec_participantes_t = "2 · Participantes",
  sec_participantes_txt = "Compara edad, diagnóstico y tratamiento entre controles y pacientes.",
  sec_exploracion_t = "3 · Exploración",
  sec_exploracion_txt = "Mira cómo se agrupan las 99 muestras (PCA) y busca la expresión de un gen.",
  sec_deg_t = "4 · Expresión diferencial",
  sec_deg_txt = "Elige una comparación y descubre qué genes están sobreexpresados o reprimidos en el volcano plot.",
  sec_go_t = "5 · Funciones biológicas",
  sec_go_txt = "Descubre en qué procesos participan los genes que cambian (g:Profiler).",
  sec_descargas_t = "6 · Descarga de los datos",
  sec_descargas_txt = "Lleva cualquier tabla del paquete a un archivo CSV.",
  sec_aprende_t = "7 · Aprende a hacer estos análisis",
  sec_aprende_txt = "Reproduce cada gráfica con las lecciones y cita el paquete y los datos.",
  conoce_t = "Conoce el paquete GeneExprLupusMx",
  conoce_txt = paste0(
    "<p><code>GeneExprLupusMx</code> es un paquete de datos de R, bilingüe (español e inglés), con los perfiles de RNA-seq, ",
    "los metadatos clínicos y los resultados intermedios de cada análisis de este estudio. Sus <b>lecciones</b> reproducen paso a paso ",
    "las figuras del artículo, de la estadística descriptiva a la expresión diferencial, el enriquecimiento funcional y los modelos con interacción.</p>",
    "<p>Está pensado para la enseñanza y la investigación reproducible en bioinformática, inmunología y ciencias biomédicas, ",
    "y para hacer visibles datos biomédicos generados en América Latina.</p>"
  ),
  conoce_sitio = "Sitio web y lecciones",
  conoce_github = "Código en GitHub",
  conoce_descarga = "Guía de descarga de datos",
  stat_mujeres = "mujeres mexicanas (10 Ctrl · 23 SLE)",
  stat_muestras = "transcriptomas de RNA-seq",
  stat_celulas = "tipos celulares: monocito · moDC · tolDC",
  stat_deg = "genes diferencialmente expresados en lupus",
  stat_isg = "genes estimulados por interferón (ISG) entre los DEG",
  stat_lecciones = "lecciones que reproducen las figuras del artículo",

  # ---- Recuadro "aprende a hacerla" ----
  aprende_caja = "¿Cómo se hace esta gráfica?",
  aprende_caja_txt = "Aprende a construirla paso a paso, con código en R, en la lección:",

  # ---- Participantes ----
  part_t = "Las participantes del estudio",
  part_intro = "Información clínica de las 33 voluntarias (objeto <code>volunteer_info</code>). Elige una variable para ver su distribución en cada grupo.",
  part_variable = "Variable",
  part_vars = c("Edad (años)" = "Age_years", "Edad al diagnóstico (años)" = "Age_at_diagnosis_years",
                "Años con el diagnóstico" = "Years_with_diagnosis", "Dosis de glucocorticoides (mg)" = "GCC_dose_mg"),
  part_resumen = "Resumen por grupo",
  part_cols = c("n", "Media", "Mediana", "Desv. estándar", "Mín.", "Máx."),
  part_n = "Participantes",
  part_tabla = "Tabla de participantes",
  part_nota = "Las variables de diagnóstico y tratamiento solo aplican a las pacientes (SLE).",

  # ---- Exploración ----
  pca_t = "¿Cómo se agrupan las 99 muestras? (PCA)",
  pca_intro = paste0(
    "Análisis de componentes principales (PCA) de la expresión transformada con VST (<code>vst(dds, blind = FALSE)</code>), ",
    "con los 13,596 genes de <code>dds</code>, igual que en la lección 2. Cada punto es una muestra; las muestras parecidas quedan cerca. ",
    "A la izquierda, cuánta variación guarda cada componente; a la derecha, elige qué componentes graficar."
  ),
  scree_t = "Varianza explicada por cada PC (scree plot)",
  scree_x = "Componente principal",
  scree_y = "Porcentaje de varianza explicada",
  scree_sel = "PC seleccionados",
  scree_nota = paste("Los círculos rojos marcan los PCs que estás graficando. Al cambiar los tipos celulares o los genes se recalcula el PCA",
                     "(como en la sección 11.4 de la lección); los genes más variables se eligen dentro de las muestras seleccionadas, como hace DESeq2::plotPCA()."),
  pca_card_t = "PCA de las muestras",
  pca_ngenes = "Genes usados",
  pca_ngenes_op = c("Todos (13,596)" = "todos", "500 más variables" = "500",
                    "250 más variables" = "250", "100 más variables" = "100"),
  heat_t = "Expresión de varios genes (heatmap)",
  heat_intro = paste0(
    "Elige un conjunto de genes para ver su expresión (z-score de la VST, <code>zscore_mat</code>) en las 99 muestras, ",
    "ordenadas por tipo celular y grupo. Los valores se limitan a [−2, 2] para que los extremos no dominen el color, como en la lección 2. ",
    "Prueba con los ISG: la <b>firma de interferón</b> aparece como columnas rojas en los monocitos de algunas pacientes."
  ),
  heat_set = "Conjunto de genes",
  heat_isg10 = "10 ISG de la lección 2",
  heat_isg_deg = "ISG diferencialmente expresados",
  heat_deg = "DEG SLE vs Ctrl",
  heat_isg = "Todos los ISG (isg_genes)",
  heat_sle = "Genes asociados al lupus (sle_genes)",
  heat_toldc = "Firma de las tolDC (tolDC_genes)",
  heat_dc = "Diferenciación de células dendríticas (dc_genes)",
  heat_modc = "Diferenciación a moDC (moDC_genes)",
  heat_custom = "Elegir genes…",
  heat_genes = "Genes",
  heat_cluster = "Agrupar genes con patrón parecido (clustering)",
  heat_n = "%d genes con expresión en los datos.",
  heat_min = "Elige al menos 2 genes.",
  heat_dl = "Descargar z-scores (CSV)",
  heat_x = "Muestras (por tipo celular y grupo)",

  pca_x = "Eje X",
  pca_y = "Eje Y",
  pca_color = "Colorear por",
  pca_celulas = "Tipos celulares",
  genes_t = "Expresión de genes",
  gen_t = "Expresión de un gen",
  agrupar_por = "Comparar por",
  agrupar_op = c("Grupo (Ctrl / SLE)" = "grupo", "Edad (< 30 / ≥ 30 años)" = "edad",
                 "Actividad de la enfermedad (SLEDAI)" = "actividad", "Uso de glucocorticoides" = "gc"),
  edad_menor = "< 30 años", edad_mayor = "≥ 30 años",
  gc_no = "SLE sin GC", gc_si = "SLE con GC",
  dot_t = "Biomarcadores y conjuntos de genes (dot plot)",
  dot_intro = paste0(
    "Para cada tipo celular y grupo: el <b>color</b> es la expresión media (z-score) y el <b>tamaño</b>, el porcentaje de muestras por encima del promedio del gen. ",
    "Con los biomarcadores de la lección 3 comprobamos que cada muestra es del tipo celular que dice ser (Figura 2B)."
  ),
  dot_marcadores = "Biomarcadores de la lección 3",
  dot_por = "Filas: tipo celular ×",
  dot_tamano = "Tamaño del punto: % de muestras por encima del promedio del gen (z-score > 0). En actividad y glucocorticoides, las controles aparecen como «Ctrl».",
  dot_color = "Expresión media (z-score)", dot_color_corto = "Expresión<br>media<br>(z-score)",
  dot_pct = "% de muestras por encima del promedio",
  dot_min = "Elige al menos 1 gen.",
  heat_por = "Separar muestras por",
  heat_pistas = "Barras de anotación",
  pistas_op = c("Tipo celular" = "celula", "Grupo" = "grupo", "Edad" = "edad", "Actividad (SLEDAI)" = "actividad", "Glucocorticoides" = "gc"),
  heat_aviso = "Los conjuntos con muchos genes (p. ej., todos los ISG) tardan más en dibujarse; espera unos segundos.",
  gen_intro = "Escribe un gen para ver su expresión (z-score de la VST) en cada tipo celular. Prueba con IFI27 o IFI44L (firma de interferón) o con CD14 y CD1C (marcadores celulares).",
  gen_sel = "Gen",

  # ---- DEG ----
  deg_t = "Genes diferencialmente expresados (DEG)",
  deg_intro = "Elige una comparación y ajusta los umbrales. El volcano plot y la tabla se actualizan al instante.",
  deg_ajustes = "Ajustes (figuras y tabla)",
  interp_t = "Interpretación del volcano plot",
  png_boton = "Descargar figura (PNG)",
  bio_t = "Biotipos de los DEG",
  bio_interp_t = "Interpretación de la distribución",
  bio_up = "DEG sobreexpresados", bio_down = "DEG reprimidos", bio_fondo = "Todos los genes analizados",
  bio_na = "Sin anotación",
  bio_x = "% de genes",
  bio_interp_1 = "De los <b>%d DEG</b> de <b>%s</b> (con los umbrales actuales), así se reparten los biotipos:",
  bio_interp_mas = "El biotipo con mayor diferencia es <b>%s</b>: representa el %.1f%% de los DEG, frente al %.1f%% de todos los genes analizados, es decir, está <b>sobrerrepresentado</b> entre los genes que cambian.",
  bio_interp_menos = "El biotipo con mayor diferencia es <b>%s</b>: representa el %.1f%% de los DEG, frente al %.1f%% de todos los genes analizados, es decir, está <b>subrepresentado</b> entre los genes que cambian.",
  bio_interp_nota = "Biotipos de <code>biotype_summary</code> (anotación de Ensembl). «Todos los genes analizados» son los genes con padj en esta comparación. Los biotipos PCG y lncRNA se usan en las lecciones 4 y 6.",
  prop_t = "¿Cuántos genes cambian? (proporción UP/DOWN)",
  prop_nota = "Como la Figura 4B (lección 5): proporción de DEG sobreexpresados y reprimidos en cada comparación de este tipo, con el número de genes dentro de la barra. En negritas, la comparación del volcano plot.",
  lfch_t = "log2FC de los DEG en todas las comparaciones (heatmap; escala de −2 a 2)",
  lfch_n = "Genes a mostrar (los de menor padj)",
  lfch_bio = "Mostrar el biotipo de cada gen",
  lfch_dir = "Dirección en la comparación elegida",
  lfch_interp_t = "Interpretación del heatmap",
  lfch_interp_1 = "El heatmap muestra los <b>%d DEG</b> con menor padj de <b>%s</b> (<span class='txt-up'><b>%d sobreexpresados</b></span> y <span class='txt-down'><b>%d reprimidos</b></span>) y su log2FC en todas las comparaciones de este tipo: %s. Rojo = el gen está sobreexpresado respecto a la referencia de esa columna; azul = reprimido.",
  lfch_interp_2 = "<b>%d de %d genes (%d%%)</b> cambian en la misma dirección en todas las columnas (cambio compartido); %d cambian en dirección opuesta en alguna columna.",
  lfch_interp_sle = "Aquí las columnas son los tres tipos celulares: un gen rojo en todas las columnas está sobreexpresado en las pacientes en monocitos, moDC y tolDC (por ejemplo, los genes de la firma de interferón).",
  lfch_interp_dif = "Aquí las columnas son la diferenciación en controles y en pacientes: si los colores se parecen entre Ctrl y SLE, el programa de diferenciación es compartido (lección 5).",
  lfch_interp_int = "Aquí las columnas son SLEvsCtrl_moDC y SLEvsCtrl_tolDC: un color igual en ambas indica que el efecto del lupus sobre la diferenciación es parecido en moDC y tolDC (lección 7).",
  diag_ctrl = "Controles sanas (Ctrl) · referencia",
  diag_sle = "Pacientes con lupus (SLE)",
  diag_referencia = "referencia",
  diag_grupo_Ctrl = "Ctrl · referencia",
  diag_grupo_SLE = "SLE · referencia",
  diag_citocinas_moDC = "GM-CSF + IL-4",          # como en la lección 5
  diag_citocinas_tolDC = "GM-CSF + IL-10",
  diag_desc_moDC = "célula dendrítica derivada de monocito",
  diag_desc_tolDC = "célula dendrítica tolerogénica",
  diag_pie_sle = "Mismo tipo celular, distinto grupo: se compara a las pacientes contra las controles.",
  diag_pie_dif = "Mismo grupo, distinto tipo celular: se compara la célula dendrítica contra los monocitos del mismo grupo (lección 5, Figura 4A). LFC &gt; 0 = sobreexpresado en la célula dendrítica.",
  diag_pie_int = "Se compara el <i>cambio</i> de la diferenciación entre grupos: LFC = (cambio en SLE) − (cambio en Ctrl).",
  cel_monocyte = "monocitos", cel_moDC = "moDC", cel_tolDC = "tolDC",
  grupo_Ctrl = "controles sanas (Ctrl)", grupo_SLE = "pacientes con lupus (SLE)",
  interp_sle = paste0(
    "En los <b>%s</b> de las pacientes con lupus (SLE), <b class='txt-up'>%d genes están sobreexpresados</b> y ",
    "<b class='txt-down'>%d genes están reprimidos</b> respecto a los %s de las controles sanas (Ctrl, referencia). ",
    "Los puntos de la derecha (LFC &gt; 0) están sobreexpresados en SLE; los de la izquierda, reprimidos en SLE."
  ),
  interp_dif = paste0(
    "En las <b>%s</b>, al diferenciarse de monocitos a <b>%s</b>, <b class='txt-up'>%d genes están sobreexpresados</b> y ",
    "<b class='txt-down'>%d genes están reprimidos</b> respecto a los monocitos (referencia). ",
    "Los puntos de la derecha (LFC &gt; 0) están sobreexpresados en las %s; los de la izquierda, reprimidos respecto a los monocitos."
  ),
  interp_int = paste0(
    "<b class='txt-up'>%d genes</b> tienen mayor sobreexpresión (o menor represión) al diferenciarse de monocitos a <b>%s</b> en las pacientes con lupus que en las controles, y ",
    "<b class='txt-down'>%d genes</b> tienen mayor represión (o menor sobreexpresión) en las pacientes al diferenciarse a %s. ",
    "Aquí el LFC compara cambios: (diferenciación en SLE) − (diferenciación en Ctrl), por lo que la referencia es el cambio en las controles."
  ),
  interp_umbral = "Umbrales actuales: |LFC| ≥ %s y padj &lt; %s.",
  interp_mejor = "El gen más significativo es <b>%s</b>, que %s (LFC = %.2f; padj = %s).",
  interp_sube = "está sobreexpresado", interp_baja = "está reprimido",
  deg_familia = "Tipo de comparación",
  deg_cmp = "Comparación",
  deg_lfc = "|log2 fold change| mínimo",
  deg_padj = "padj máximo",
  deg_buscar = "Resaltar genes",
  deg_buscar_ph = "Escribe uno o más genes…",
  deg_buscar_ayuda = paste0(
    "<b>Resaltar genes:</b> escribe uno o más genes para marcarlos con un círculo morado. Puedes buscar <b>cualquier gen</b> de la comparación ",
    "(21,718), sea o no DEG: si no pasa los umbrales aparecerá entre los puntos grises."
  ),
  deg_genes_fuera = "No aparecen en el volcano porque DESeq2 no les asignó padj (filtrado independiente): %s.",
  deg_resaltados = "Genes resaltados",
  deg_etiquetas = "Mostrar nombres de los DEG más significativos",
  deg_top = "Número de genes a etiquetar",
  deg_top_info = "Se etiquetan los %d DEG con menor padj (de %d que pasan los umbrales). Si dos nombres se enciman, arrastra uno con el mouse.",
  deg_top_cero = "Ningún gen pasa los umbrales actuales.",
  deg_todos = "Mostrar todos los genes en la tabla (no solo DEG)",
  deg_resumen_t = "Número de DEG en todas las comparaciones (con los umbrales actuales)",
  deg_resumen_nota = "En negritas, se resalta la comparación mostrada en el Volcano plot.",
  deg_volcano_t = "Volcano plot",
  deg_tabla_t = "Tabla de resultados",
  deg_up = "Sobreexpresados", deg_down = "Reprimidos", deg_ns = "Sin cambio",
  deg_dl_todos = "Descargar comparación completa (CSV)",
  deg_dl_deg = "Descargar solo DEG (CSV)",
  deg_nota_interaccion = paste0(
    "<b>Ojo:</b> en la interacción el log2 fold change es una <i>diferencia de diferencias</i>: ",
    "(célula dendrítica vs monocito en SLE) − (célula dendrítica vs monocito en Ctrl). ",
    "En el artículo estos genes se seleccionaron solo con padj &lt; 0.05 (sin umbral de LFC)."
  ),
  deg_nota_umbral = "En el artículo, los DEG se definieron con |LFC| ≥ 1 y padj < 0.05.",
  fam_sle = "SLE vs Ctrl, en cada tipo celular",
  fam_dif = "Diferenciación: célula dendrítica vs monocito, en cada grupo",
  fam_int = "Efecto de interacción del SLE en la diferenciación (SLEvsCtrl)",
  fam_sle_desc = "Objeto <code>de_results</code> · ¿Qué genes cambian en las pacientes respecto a las controles? (Figura 2 C–E del artículo)",
  fam_dif_desc = "Objeto <code>de_celltype_results</code> · ¿Qué genes cambian al diferenciar los monocitos a moDC o tolDC dentro de cada grupo? (Figuras 4 y 5)",
  fam_int_desc = paste0("Objeto <code>de_interaction_results</code> · <b>Efecto de interacción del SLE</b> (<i>SLE interaction effect on the in vitro moDCs and tolDCs differentiation</i>, ",
                        "): ¿qué genes cambian de forma distinta al diferenciarse en las pacientes que en las controles? ",
                        "Comparaciones <span style='color:#4daf4a'><b>SLEvsCtrl_moDC</b></span> y <span style='color:#ff7f00'><b>SLEvsCtrl_tolDC</b></span>; ",
                        "en la lección 7 se presentan como <i>between-group comparison of contrasts</i> (Figura 7)."),
  # Modelo de DESeq2 de cada comparación (documentado en ?de_results, ?de_celltype_results y ?de_interaction_results)
  modelo_sle = paste0(
    "<b>Modelo (DESeq2):</b> un modelo por tipo celular, solo con las muestras de %s: <code>~ Group</code> (referencia: Ctrl). ",
    "Resultado: <code>results(name = \"Group_SLE_vs_Ctrl\")</code>; prueba de Wald y corrección de Benjamini–Hochberg."
  ),
  modelo_dif_Ctrl = paste0(
    "<b>Modelo (DESeq2):</b> solo las muestras de las controles: <code>~ Cell_type</code> (referencia: monocyte). ",
    "Resultado: <code>results(name = \"Cell_type_%s_vs_monocyte\")</code>; prueba de Wald y corrección de Benjamini–Hochberg."
  ),
  modelo_dif_SLE = paste0(
    "<b>Modelo (DESeq2):</b> solo las muestras de las pacientes: <code>~ SLEDAI_category + Cell_type + Dose_category</code> (referencia: monocyte), ",
    "que descuenta el efecto de la actividad de la enfermedad (SLEDAI) y de la dosis de glucocorticoides. ",
    "Resultado: <code>results(name = \"Cell_type_%s_vs_monocyte\")</code>; prueba de Wald y corrección de Benjamini–Hochberg."
  ),
  modelo_int = paste0(
    "<b>Modelo (DESeq2):</b> las 99 muestras juntas: <code>~ Group + Cell_type + Group:Cell_type</code> (referencias: Ctrl y monocyte). ",
    "El término de interacción <code>Group:Cell_type</code> mide si el cambio de monocito a %s es distinto en SLE. ",
    "Resultado: <code>results(name = \"GroupSLE.Cell_type%s\")</code>; prueba de Wald y corrección de Benjamini–Hochberg."
  ),
  cmp_monocyte = "Monocitos: SLE vs Ctrl",
  cmp_moDC = "moDC: SLE vs Ctrl",
  cmp_tolDC = "tolDC: SLE vs Ctrl",
  cmp_Ctrl_moDC = "Controles: moDC vs monocitos",
  cmp_Ctrl_tolDC = "Controles: tolDC vs monocitos",
  cmp_SLE_moDC = "SLE: moDC vs monocitos",
  cmp_SLE_tolDC = "SLE: tolDC vs monocitos",
  cmp_int_moDC = "SLEvsCtrl_moDC: diferenciación a moDC en SLE vs en Ctrl",
  cmp_int_tolDC = "SLEvsCtrl_tolDC: diferenciación a tolDC en SLE vs en Ctrl",
  col_gen = "Gen", col_dir = "Dirección", col_biotipo = "Biotipo",
  col_isg = "ISG", col_sle = "Asociado a lupus", col_toldc = "Firma tolDC",
  si = "Sí",

  # ---- GO ----
  go_t = "¿En qué procesos participan los DEG? (enriquecimiento funcional)",
  go_intro = "Términos significativos de g:Profiler (FDR < 0.05) para cada lista de DEG. Cada punto es un término; más a la derecha, más significativo.",
  go_conjunto = "Conjunto de resultados",
  go_lista = "Lista de genes",
  go_fuente = "Base de datos",
  go_top = "Número de términos (en total o por base de datos)",
  go_ajustes = "Ajustes (figuras y tabla)",
  go_manhattan_t = "Panorama: todos los términos significativos por base de datos (Manhattan plot)",
  go_manhattan_nota = paste0(
    "Como en g:Profiler: cada punto es un término significativo (FDR &lt; 0.05) de la lista elegida, agrupado por base de datos; ",
    "más arriba, más significativo (tope en 16). El tamaño indica cuántos genes tiene el término. ",
    "Muestra <b>todas</b> las bases de datos; las casillas de la izquierda filtran el dot plot y la tabla."
  ),
  go_manhattan_y = "−log10(padj)",
  go_modo = "Términos a mostrar en el dot plot",
  go_modo_op = c("Los más significativos de cada base de datos" = "fuente", "Los más significativos en total" = "total"),
  go_fuentes_nota = "Para este conjunto, el paquete consultó solo GO:BP, KEGG y Reactome (ver ?go_celltype_results y ?go_interaction_results).",
  fuentes_nombres = c("GO:MF" = "GO: función molecular", "GO:BP" = "GO: proceso biológico", "GO:CC" = "GO: componente celular",
                      KEGG = "KEGG", REAC = "Reactome", WP = "WikiPathways", TF = "Factores de transcripción",
                      MIRNA = "miRNA", HPA = "Human Protein Atlas", CORUM = "Complejos proteicos", HP = "Fenotipo humano"),
  go_fig_t = "Términos más significativos (dot plot)",
  go_heat_t = "Términos compartidos entre listas (heatmap)",
  go_heat_nota = paste0(
    "Filas: los términos más significativos de <b>cada lista</b> (el número se elige en el panel); columnas: todas las listas del conjunto. ",
    "El color del círculo es el número de genes del término; en blanco, el término no es significativo en esa lista. ",
    "Los términos se agrupan según en qué listas aparecen."
  ),
  go_heat_fuente = "Base de datos del heatmap",
  go_heat_expr = "Expresión",
  go_heat_comparacion = "Comparación",
  go_heat_color = "Número de genes",
  heat_pista_celula2 = "Tipo celular",
  heat_pista_grupo2 = "Grupo",
  go_tabla_t = "Tabla de términos",
  go_dl = "Descargar términos (CSV)",
  go_set_sle = "SLE vs Ctrl por tipo celular (go_results)",
  go_set_dif = "Diferenciación (go_celltype_results)",
  go_set_int = "Efecto de interacción del SLE (go_interaction_results)",
  go_listas_ayuda = "Las listas se nombran <code>&lt;grupo&gt;_&lt;UP|DOWN&gt;_&lt;tipo celular&gt;</code>: UP = genes sobreexpresados, DOWN = genes reprimidos.",
  go_x = "−log10(padj)", go_tam = "Genes en el término",
  go_vacio = "No hay términos significativos para esta combinación.",

  # ---- Descargas ----
  dl_t = "Descarga los datos",
  dl_intro = paste0(
    "Todas las tablas provienen del paquete <code>GeneExprLupusMx</code>, generadas a partir de los datos del artículo ",
    "<a href='", URL_ARTICULO, "' target='_blank'>Hernández-Ledesma, Coss-Navarrete <i>et al.</i> (2026)</a>. ",
    "Si las usas, <b>cita el paquete y los datos</b> (artículo y Zenodo; ver sección 7)."
  ),
  dl_objeto = "Tabla",
  dl_boton = "Descargar CSV",
  dl_dim = "Dimensiones",
  dl_vista = "Vista previa (primeras filas y columnas)",
  dl_r = "¿Prefieres R? Todo está disponible al cargar el paquete:",
  dl_dds = "El objeto <code>dds</code> (DESeqDataSet) solo puede usarse en R con DESeq2.",
  leccion_descarga = "Guía de descarga de datos",
  dl_origen_t = "¿De dónde vienen los datos y cómo citarlos?",
  dl_cita_datos_t = "Cómo citar los datos (artículo y Zenodo)",
  dl_cita_paquete_t = "Si además usas el paquete o sus lecciones",
  dl_citar_aviso = paste0(
    "<b>Importante:</b> si usas estos datos o resultados, <b>cítalos con base en el artículo original</b> ",
    "(Hernández-Ledesma, Coss-Navarrete <i>et al.</i>, 2026, <i>Front. Immunol.</i> 17:1771959) y el análisis depositado en Zenodo. ",
    "Si además usas el paquete o sus lecciones, cita también <code>GeneExprLupusMx</code>."
  ),
  dl_origen_tabla = paste0(
    "<table class='table table-sm'><thead><tr><th>Fuente</th><th>¿Qué contiene?</th><th>Licencia</th></tr></thead><tbody>",
    "<tr><td><a href='", URL_ARTICULO, "' target='_blank'>Artículo original</a> (<i>Front. Immunol.</i> 2026)</td>",
    "<td>Diseño del estudio, métodos, figuras y material suplementario (incluidos los identificadores de EGA).</td>",
    "<td><a href='https://creativecommons.org/licenses/by/4.0/' target='_blank'>CC BY</a> (acceso abierto)</td></tr>",
    "<tr><td><a href='", URL_ZENODO, "' target='_blank'>Zenodo, registro 17742206</a></td>",
    "<td>Análisis de expresión génica depositado por las autoras: cuantificación con nf-core/rnaseq (Salmon) y análisis con tximport y DESeq2.</td>",
    "<td><a href='https://creativecommons.org/licenses/by/4.0/' target='_blank'>CC BY 4.0</a></td></tr>",
    "<tr><td><a href='", URL_CODIGO, "' target='_blank'><i class='fab fa-github'></i> SLE_Transcriptomic_analysis_tolerogenic</a></td>",
    "<td>Código del análisis original y sus figuras.</td>",
    "<td>Sin licencia declarada en el repositorio</td></tr>",
    "<tr><td><a href='", URL_PKG, "' target='_blank'><i class='fab fa-github'></i> GeneExprLupusMx</a></td>",
    "<td>Paquete de R: los datos reorganizados en objetos de R, las lecciones y las funciones.</td>",
    "<td>MIT</td></tr>",
    "<tr><td>European Genome-phenome Archive (EGA)</td>",
    "<td>Datos crudos de secuenciación (FASTQ); los identificadores están en la Tabla suplementaria 2 del artículo.</td>",
    "<td>Acceso controlado (requiere solicitud y aprobación)</td></tr>",
    "</tbody></table>"
  ),
  dl_licencia_nota = paste0(
    "CC BY permite usar, distribuir y adaptar el material siempre que se dé crédito a las autoras y se cite la publicación original. ",
    "El repositorio de código no declara una licencia; para reutilizar ese código, consulta a las autoras. ",
    "Licencias verificadas en el artículo, en Zenodo y en el paquete (octubre de 2026)."
  ),

  # ---- Aprende ----
  aprende_t = "Aprende a hacer estos análisis",
  aprende_intro = "Cada gráfica de esta app se construye paso a paso, con código en R, en las lecciones del paquete. Están pensadas para estudiantes y docentes; no necesitas experiencia previa en RNA-seq.",
  aprende_instalar = "Instala el paquete",
  aprende_paquete_txt = "Visita el sitio del paquete para explorar todas las lecciones, la referencia de cada objeto y las funciones. Si usas las lecciones, las funciones o los objetos de R en tu trabajo o en tus clases, <b>cita el paquete</b>:",
  aprende_datos_cita = "Para citar los datos (artículo original y Zenodo) y conocer sus licencias, ve a la sección ",
  cita_t = "¿Cómo citar?",
  cita_txt = "Hay dos citas: una para el <b>paquete</b> y otra para los <b>datos</b>. Si usas este recurso, por favor incluye ambas.",
  cita_paquete_t = "1) El paquete (lecciones, funciones y objetos de R)",
  cita_paquete = paste0("Coss E (2026) <i>GeneExprLupusMx: Datos transcriptómicos de lupus en mujeres mexicanas con lecciones de aprendizaje</i>. ",
                        "Paquete de R. <a href='", URL_PKG, "' target='_blank'>", URL_PKG, "</a>."),
  cita_pendiente = "<b>Pendiente:</b> el DOI de Zenodo del paquete se agregará próximamente; mientras tanto, cita el repositorio de GitHub y la versión que usaste (<code>packageVersion(\"GeneExprLupusMx\")</code>). En R, <code>citation(\"GeneExprLupusMx\")</code> muestra la cita.",
  cita_datos_t = "2) Los datos (artículo original y análisis en Zenodo)",
  cita_articulo = paste0("Hernández-Ledesma AL, Coss-Navarrete EL, <i>et al</i> (2026) <a href='", URL_ARTICULO, "' target='_blank'>Gene expression profiling of dendritic cell ",
                         "tolerance dysfunction in women with systemic lupus erythematosus</a>. <i>Front. Immunol.</i> 17:1771959. ",
                         "<a href='", URL_DOI, "' target='_blank'>doi: 10.3389/fimmu.2026.1771959</a>."),
  cita_zenodo = paste0("Coss-Navarrete EL, Hernández-Ledesma AL, <i>et al</i> (2025) <a href='", URL_ZENODO, "' target='_blank'>Analysis of gene expression profiles during ",
                       "tolerance induction in dendritic cells from women with Systemic Lupus Erythematosus</a> [Data set]. <i>Zenodo</i>. ",
                       "<a href='", URL_ZENODO, "' target='_blank'>doi: 10.5281/zenodo.17742206</a>."),
  cita_codigo = paste0("Código del análisis original: <a href='", URL_CODIGO, "' target='_blank'>NeuroGenomicsMX/SLE_Transcriptomic_analysis_tolerogenic</a>."),
  cita_app = "Esta aplicación es un complemento interactivo del paquete; los resultados son los del artículo original.",
  lecciones_titulo = c("1 · Conociendo a las participantes",
                       "2 · De las cuentas de RNA-seq a la exploración de datos",
                       "3 · Genes diferencialmente expresados entre pacientes y controles",
                       "4 · Funciones biológicas de los genes diferencialmente expresados",
                       "5 · Diferenciación de monocitos a células dendríticas",
                       "6 · La firma transcripcional de las células dendríticas tolerogénicas",
                       "7 · ¿Se diferencian igual las células de pacientes y controles?"),
  lecciones_figura = c("Figura S3", "PCA y heatmaps", "Figura 2", "Figuras 3 y S4", "Figuras 4 y 5", "Figura 6", "Figura 7"),
  funciones_t = "Cómo usar las funciones del paquete",
  funciones_txt = "Heatmaps con create_heatmap() y get_heatmap_mat(), y reducción de términos GO con processGO() y Reduce_GOTerm().",
  lecciones_claves = list(
    c("estadística descriptiva", "Tabla 1", "frecuencias", "histograma y densidad", "Q–Q y Shapiro–Wilk", "prueba t, Wilcoxon y Fisher"),
    c("tximport", "DESeqDataSet", "filtrado de genes", "normalización", "VST", "heatmap de z-scores", "PCA"),
    c("biomarcadores celulares", "DESeq2", "prueba de Wald", "|LFC| ≥ 1 y padj < 0.05", "volcano plots", "ComplexHeatmap"),
    c("g:Profiler", "GO y Reactome", "genes asociados al lupus", "ISG", "prueba hipergeométrica", "dot plots"),
    c("diferenciación in vitro", "diseño con covariables", "proporción de DEG", "matriz de LFC", "rutas de Reactome", "ácidos grasos"),
    c("firma tolDC (Robertson 2021)", "heatmaps de expresión y LFC", "biotipos PCG y lncRNA", "nivel de referencia"),
    c("modelo con interacción", "Group:Cell_type", "diferencia de diferencias", "genes compartidos y exclusivos", "procesos de Reactome")
  ),
  funciones_claves = c("create_heatmap()", "get_heatmap_mat()", "processGO()", "Reduce_GOTerm()", "rrvgo")
))

tr <- function(key) {
  x <- TXT[[LANG]][[key]]
  if (is.null(x)) TXT$es[[key]] else x
}

# Lecciones: id del archivo en el sitio pkgdown; títulos y figuras en TXT
LECCIONES <- data.frame(
  id = c("clinical", "intro_RNAseq", "l3_DEGs", "l4_GOterms", "l5_DEGs_monocytes_diff", "l6_DEGs_tolDcs", "l7_DEGs_between_groups"),
  titulo = tr("lecciones_titulo"),
  figura = tr("lecciones_figura"),
  img = sprintf("lesson%d_overview_%s.png", 1:7, LANG)
)
LECCIONES$claves <- tr("lecciones_claves")
