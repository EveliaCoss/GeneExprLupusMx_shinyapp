# Explorador interactivo de GeneExprLupusMx
# Ejecutar: shiny::runApp()  (idioma: Sys.setenv(APP_LANG = "es"))
library(shiny)
library(bslib)
library(DT)
library(plotly)

source("textos.R", encoding = "UTF-8")

# ---- Datos ----------------------------------------------------------------
# Orden: 1) carpeta en GENEEXPRLUPUSMX_DATA (útil si el paquete instalado está
# desactualizado), 2) paquete instalado, 3) copia local del paquete en desarrollo.
DATA_DIR <- Sys.getenv("GENEEXPRLUPUSMX_DATA", "")
DATA_LOCAL <- "~/Documents/GeneExprLupusMx/data"

.cache <- new.env()
obj <- function(name) {
  if (is.null(.cache[[name]])) {
    e <- new.env()
    if (nzchar(DATA_DIR)) {
      load(file.path(DATA_DIR, paste0(name, ".rda")), envir = e)
    } else if (requireNamespace("GeneExprLupusMx", quietly = TRUE)) {
      data(list = name, package = "GeneExprLupusMx", envir = e)
    } else {
      load(file.path(DATA_LOCAL, paste0(name, ".rda")), envir = e)
    }
    x <- e[[name]]
    # Orden de los tipos celulares como en las lecciones
    if (name == "metadata_cleaned") x$Cell_type <- factor(x$Cell_type, c("monocyte", "moDC", "tolDC"))
    .cache[[name]] <- x
  }
  .cache[[name]]
}

# PCA de la lección 2 (VST de dds), precalculados con tools/pca.R por combinación de tipos celulares
PCA <- readRDS("datos/pca.rds")

# Comparaciones de expresión diferencial disponibles en el paquete
CMP <- data.frame(
  id = c("monocyte", "moDC", "tolDC", "Ctrl_moDC", "Ctrl_tolDC", "SLE_moDC", "SLE_tolDC", "int_moDC", "int_tolDC"),
  familia = rep(c("sle", "dif", "int"), c(3, 4, 2)),
  objeto = rep(c("de_results", "de_celltype_results", "de_interaction_results"), c(3, 4, 2)),
  elemento = c("monocyte", "moDC", "tolDC", "Ctrl_moDC", "Ctrl_tolDC", "SLE_moDC", "SLE_tolDC", "moDC", "tolDC"),
  leccion = rep(c("l3_DEGs", "l5_DEGs_monocytes_diff", "l7_DEGs_between_groups"), c(3, 4, 2))
)
CMP$etiqueta <- vapply(CMP$id, function(i) tr(paste0("cmp_", i)), "")
FAMILIAS <- setNames(c("sle", "dif", "int"), c(tr("fam_sle"), tr("fam_dif"), tr("fam_int")))

GO_SETS <- data.frame(
  id = c("sle", "dif", "int"),
  objeto = c("go_results", "go_celltype_results", "go_interaction_results"),
  leccion = c("l4_GOterms", "l5_DEGs_monocytes_diff", "l7_DEGs_between_groups"),
  etiqueta = c(tr("go_set_sle"), tr("go_set_dif"), tr("go_set_int"))
)

DESCARGABLES <- c(
  "metadata_cleaned", "volunteer_info", "biotype_summary",
  paste0("de_results$", c("monocyte", "moDC", "tolDC")),
  paste0("de_celltype_results$", c("Ctrl_moDC", "Ctrl_tolDC", "SLE_moDC", "SLE_tolDC")),
  paste0("de_interaction_results$", c("moDC", "tolDC")),
  "go_results", "go_celltype_results", "go_interaction_results",
  "sle_genes", "isg_genes", "dc_genes", "moDC_genes", "tolDC_genes",
  "zscore_mat", paste0("filtered_txi$", c("counts", "abundance", "length"))
)
tabla_descarga <- function(key) {
  p <- strsplit(key, "$", fixed = TRUE)[[1]]
  x <- obj(p[1])
  if (length(p) == 2) x[[p[2]]] else x
}
escribir_csv <- function(x, file) write.csv(x, file, row.names = is.matrix(x))

anotaciones <- function() {
  if (is.null(.cache$anot)) {
    b <- obj("biotype_summary"); t <- obj("tolDC_genes")
    .cache$anot <- list(
      biotipo = setNames(b$Group_biotype, b$gene),
      isg = unique(obj("isg_genes")$genes),
      sle = unique(obj("sle_genes")$Gene_name),
      toldc = setNames(t$tolDC_signature, t$gene)
    )
  }
  .cache$anot
}

direccion <- function(d, lfc, padj) {
  sig <- !is.na(d$padj) & d$padj < padj
  ifelse(sig & d$log2FoldChange >= lfc, "UP", ifelse(sig & d$log2FoldChange <= -lfc, "DOWN", "NS"))
}

# Colores de las lecciones del paquete (group_color, cell_color y volcano de la lección 3)
COL_GRUPO   <- c(Ctrl = "#96D4CC", SLE = "#B493B4")
COL_CELULA  <- c(monocyte = "black", moDC = "palevioletred", tolDC = "gold")
COLORES_DIR <- c(UP = "#CD2626", DOWN = "#1874CD", NS = "#7F7F7F")   # firebrick3, dodgerblue3 y gray50 en hexadecimal (válidos también en CSS)
COL_PALETA  <- c("#3D3F8F", "#007F99", "#B493B4", "#E0A100", "#96D4CC", "#7F62B3", "#E05252", "#14164D", "#F2D77A", "#5DA9E9")

# Agrupaciones de muestras: las controles quedan como "Ctrl" en actividad y glucocorticoides
# (SLEDAI_category = "Zero" incluye a las 10 controles y ninguna control recibe glucocorticoides)
PALETAS <- list(
  grupo = COL_GRUPO,
  edad = setNames(c("#C9ECE7", "#5FB3A8", "#E3D3E3", "#8E5F8E"),
                  paste(rep(c("Ctrl", "SLE"), each = 2), c(tr("edad_menor"), tr("edad_mayor")))),
  actividad = c(Ctrl = "#96D4CC", "SLE · SLEDAI 0" = "#E3D3E3", "SLE · SLEDAI 1–3" = "#B493B4", "SLE · SLEDAI ≥ 4" = "#6E3F6E"),
  gc = setNames(c("#96D4CC", "#E3D3E3", "#7F62B3"), c("Ctrl", tr("gc_no"), tr("gc_si")))
)
agrupar <- function(m, var) {
  sle <- m$Group == "SLE"
  x <- switch(var,
    grupo = as.character(m$Group),
    edad = paste(m$Group, ifelse(as.numeric(as.character(m$Age)) < 30, tr("edad_menor"), tr("edad_mayor"))),
    actividad = ifelse(sle, c(Zero = "SLE · SLEDAI 0", No_activity = "SLE · SLEDAI 1–3", Activity = "SLE · SLEDAI ≥ 4")[as.character(m$SLEDAI_category)], "Ctrl"),
    gc = ifelse(sle, ifelse(m$GC == "Si", tr("gc_si"), tr("gc_no")), "Ctrl"))
  factor(x, names(PALETAS[[var]]))
}
# Biomarcadores de la lección 3 (Figura 2B)
MARCADORES <- list(
  monocyte = c("CD36", "CSF3R", "DPYD", "FCN1", "VCAN"),
  moDC     = c("CD1C", "IRF4", "NDRG2", "CCRL2", "SPARC"),
  tolDC    = c("IGF2BP3", "CD163", "FCGR2A", "DTX1", "GK", "EPB41L3", "TNFSF9", "CXCL1", "EFEMP2")
)
# Tipo celular y grupo de cada muestra en el orden de las lecciones (monocyte, moDC, tolDC)
condicion <- function(m, var) {
  g <- agrupar(m, var)
  lv <- as.vector(t(outer(levels(m$Cell_type), levels(g), paste, sep = " · ")))
  factor(paste(m$Cell_type, g, sep = " · "), lv)
}
leyenda <- function(colores, titulo) tags$span(tags$b(paste0(titulo, ":")),
  lapply(names(colores), function(k) tagList(tags$span(class = "chip", style = paste0("background:", colores[[k]])), k)))

# Colores de las comparaciones de contrastes entre grupos (comparison_color de la lección 7)
COL_COMPARACION <- c(SLEvsCtrl_moDC = "#4daf4a", SLEvsCtrl_tolDC = "#ff7f00")
color_texto <- function(x, color) sprintf("<span style='color:%s'><b>%s</b></span>", color, x)

# Biotipos (biotype_color de las lecciones 4 y 6 para PCG y lncRNA; el resto, de la misma paleta Dark2)
COL_BIOTIPO <- c(PCG = "#7570B3", lncRNA = "#E7298A", Pseudogene = "#66A61E", "small ncRNAs" = "#E6AB02", TEC = "#A6761D", "NA" = "#BDBDBD")

# Colores de g:Profiler (gprofiler2::gostplot) para cada base de datos, en su orden
COL_FUENTES <- c("GO:MF" = "#dc3912", "GO:BP" = "#ff9900", "GO:CC" = "#109618", KEGG = "#dd4477", REAC = "#3366cc",
                 WP = "#0099c6", TF = "#5574a6", MIRNA = "#22aa99", HPA = "#6633cc", CORUM = "#66aa00", HP = "#990099")

TEMA <- bs_theme(
  bootswatch = "flatly", primary = "#3D3F8F", "navbar-bg" = "#3D3F8F",
  "body-color" = "#14164D", "link-color" = "#007F99", "link-hover-color" = "#005F73", "border-color" = "#DFE2FB",
  heading_font = font_link("News Cycle", "https://fonts.googleapis.com/css2?family=News+Cycle:wght@400;700&display=swap")
)

caja_leccion <- function(id) {
  div(class = "alert alert-info mt-2 mb-0",
      strong(icon("graduation-cap"), " ", tr("aprende_caja")), " ", tr("aprende_caja_txt"), " ",
      a(href = leccion_url(id), target = "_blank", LECCIONES$titulo[LECCIONES$id == id]),
      sprintf("(%s).", LECCIONES$figura[LECCIONES$id == id]))
}

# Indicadores del inicio, calculados de los objetos del paquete
.deg_sle <- unique(unlist(lapply(obj("de_results"), function(d) d$gene[direccion(d, 1, 0.05) != "NS"])))
ESTADISTICAS <- c(
  stat_mujeres = nrow(obj("volunteer_info")), stat_muestras = nrow(obj("metadata_cleaned")), stat_celulas = 3,
  stat_deg = length(.deg_sle), stat_isg = sum(.deg_sle %in% obj("isg_genes")$genes), stat_lecciones = nrow(LECCIONES)
)
# Sobreexpresados (↑) y reprimidos (↓) en cada tipo celular, debajo del número de DEG
.n_dir <- sapply(obj("de_results"), function(d) table(factor(direccion(d, 1, 0.05), c("UP", "DOWN"))))
DEG_DETALLE <- do.call(sprintf, c(list(tr("stat_deg_detalle")), as.list(.n_dir[, c("monocyte", "moDC", "tolDC")])))

caja_valor <- function(k) {
  extra <- if (k == "stat_deg") p(class = "small mb-0", DEG_DETALLE)
  value_box(title = tr(k), value = ESTADISTICAS[[k]], extra, class = if (k == "stat_lecciones") "destacado")
}

# Panel del inicio que lleva a una sección de la app
panel_seccion <- function(id) {
  tags$a(class = "panel-seccion", href = "#",
         onclick = sprintf("Shiny.setInputValue('ir_a', '%s', {priority: 'event'}); window.scrollTo(0, 0); return false;", id),
         card(class = "h-100 mb-0",
              img(src = sprintf("sec_%s.png", id), alt = tr(paste0("sec_", id, "_t")), style = "width:100%"),
              card_body(strong(tr(paste0("sec_", id, "_t"))), p(class = "small mb-0", tr(paste0("sec_", id, "_txt"))))))
}

instalar_ui <- function() tagList(
  pre(style = "white-space: pre-wrap; overflow-wrap: anywhere;", "# install.packages(\"remotes\")\nremotes::install_github(\"EveliaCoss/GeneExprLupusMx\")\nlibrary(GeneExprLupusMx)\nvignette(package = \"GeneExprLupusMx\")"),
  a(href = URL_SITE, target = "_blank", URL_SITE)
)

figura <- function(src, pie) div(class = "figura",
  img(src = src, alt = gsub("<[^>]+>", "", pie)), p(class = "text-muted small mb-0", HTML(pie)))

# Descarga la figura plotly tal como se ve, en PNG a 3× de resolución (Plotly.downloadImage en el navegador)
boton_png <- function(id, nombre, entrada = "deg_cmp") tags$button(
  type = "button", class = "btn btn-outline-primary btn-sm align-self-start mt-1",
  onclick = sprintf("Plotly.downloadImage(document.getElementById('%s'), {format: 'png', scale: 3, filename: 'GeneExprLupusMx_%s_' + (Shiny.shinyapp.$inputValues.%s || '')})", id, nombre, entrada),
  icon("image"), " ", tr("png_boton"))

claves <- function(x) div(class = "claves", lapply(x, function(k) span(class = "badge clave", k)))

# ---- UI -------------------------------------------------------------------
ui <- page_navbar(
  id = "nav",
  title = tr("titulo"),
  theme = TEMA,
  navbar_options = navbar_options(bg = "#3D3F8F", theme = "dark"),
  header = tags$head(
    tags$link(rel = "stylesheet", href = "estilo.css"),
    # Mermaid, como en las lecciones del paquete (diagramas de la interpretación del volcano)
    tags$script(type = "module", HTML('import mermaid from "https://cdn.jsdelivr.net/npm/mermaid@11/dist/mermaid.esm.min.mjs";
      mermaid.initialize({ startOnLoad: false, themeVariables: { fontSize: "15px" } }); window.mermaid = mermaid;'))),
  fillable = FALSE,
  lang = LANG,
  footer = div(class = "pie small text-center p-3 mt-4", HTML(tr("pie"))),

  nav_panel(tr("tab_inicio"), value = "inicio",
    h3(tr("inicio_bienvenida")), HTML(tr("inicio_intro")),
    do.call(layout_column_wrap, c(list(width = "180px", class = "mb-3"), unname(lapply(names(ESTADISTICAS), caja_valor)))),
    card(card_header(tr("inicio_estudio_t")),
      layout_columns(col_widths = c(4, 4, 4),
        div(HTML(tr("inicio_estudio"))),
        figura("datos_paquete.png", tr("inicio_datos_img")),
        figura("diferenciacion.png", tr("inicio_dif_img"))
      )),
    layout_columns(col_widths = c(7, 5),
      card(card_header(tr("inicio_datos_t")), HTML(tr("inicio_datos"))),
      card(card_header(tr("inicio_ruta_t")), HTML(tr("inicio_ruta")))
    ),
    h4(tr("sec_t")), p(class = "text-muted", tr("sec_intro")),
    do.call(layout_column_wrap, c(list(width = 1/3, class = "mb-3"),
                                  lapply(c("participantes", "exploracion", "deg", "go", "descargas", "aprende"), panel_seccion))),
    card(card_header(tr("conoce_t")),
      layout_columns(col_widths = c(6, 6),
        div(HTML(tr("conoce_txt")),
            div(class = "d-flex gap-3 flex-wrap",
                a(href = URL_SITE, target = "_blank", icon("book-open"), " ", tr("conoce_sitio")),
                a(href = URL_PKG, target = "_blank", icon("github"), " ", tr("conoce_github")),
                a(href = leccion_url("download"), target = "_blank", icon("download"), " ", tr("conoce_descarga")))),
        div(h5(tr("aprende_instalar")), instalar_ui())
      ))
  ),

  nav_panel(tr("tab_participantes"), value = "participantes",
    h3(tr("part_t")), p(HTML(tr("part_intro"))),
    layout_columns(col_widths = c(3, 9),
      card(selectInput("part_var", tr("part_variable"), tr("part_vars")), p(class = "small text-muted", tr("part_nota"))),
      card(plotlyOutput("part_plot", height = 350), card_footer(caja_leccion("clinical")))
    ),
    card(card_header(tr("part_resumen")), tableOutput("part_resumen")),
    card(card_header(tr("part_tabla")), DTOutput("part_tabla"))
  ),

  nav_panel(tr("tab_exploracion"), value = "exploracion",
    h3(tr("pca_t")), p(HTML(tr("pca_intro"))),
    layout_columns(col_widths = c(5, 7),
      card(card_header(tr("scree_t")),
        plotlyOutput("scree_plot", height = 430),
        p(class = "small text-muted mb-0", tr("scree_nota"))),
      card(card_header(tr("pca_card_t")),
        layout_columns(col_widths = c(4, 4, 4),
          selectInput("pca_x", tr("pca_x"), c("PC1" = 1)),
          selectInput("pca_y", tr("pca_y"), c("PC2" = 2)),
          selectInput("pca_ngenes", tr("pca_ngenes"), tr("pca_ngenes_op"))
        ),
        layout_columns(col_widths = c(4, 8),
          selectInput("pca_color", tr("pca_color"), c("Cell_type", "Group", "SLEDAI_category", "Group_age", "GC")),
          checkboxGroupInput("pca_celulas", tr("pca_celulas"), c("monocyte", "moDC", "tolDC"), c("monocyte", "moDC", "tolDC"), inline = TRUE)
        ),
        plotlyOutput("pca_plot", height = 380))
    ),
    caja_leccion("intro_RNAseq"),
    h3(tr("genes_t")),
    layout_columns(col_widths = c(5, 7),
      card(card_header(tr("gen_t")), p(class = "small", tr("gen_intro")),
        layout_columns(col_widths = c(6, 6),
          selectizeInput("gen", tr("gen_sel"), choices = NULL),
          selectInput("gen_por", tr("agrupar_por"), tr("agrupar_op"))),
        plotlyOutput("gen_plot", height = 420)),
      card(card_header(tr("dot_t")), p(class = "small", HTML(tr("dot_intro"))),
        layout_columns(col_widths = c(6, 6),
          selectInput("dot_set", tr("heat_set"), NULL),
          selectInput("dot_por", tr("dot_por"), tr("agrupar_op"))),
        conditionalPanel("input.dot_set == 'custom'", selectizeInput("dot_genes", tr("heat_genes"), NULL, multiple = TRUE)),
        uiOutput("dot_ui"),
        p(class = "small text-muted mb-0", tr("dot_tamano")),
        card_footer(caja_leccion("l3_DEGs")))
    ),
    h3(tr("heat_t")), p(HTML(tr("heat_intro"))),
    layout_columns(col_widths = c(3, 9),
      card(
        selectInput("heat_set", tr("heat_set"), NULL),
        conditionalPanel("input.heat_set == 'custom'", selectizeInput("heat_genes", tr("heat_genes"), NULL, multiple = TRUE)),
        selectInput("heat_por", tr("heat_por"), tr("agrupar_op")),
        checkboxGroupInput("heat_pistas", tr("heat_pistas"), tr("pistas_op"), c("celula", "grupo")),
        checkboxInput("heat_cluster", tr("heat_cluster"), TRUE),
        uiOutput("heat_n"),
        div(class = "alert alert-warning small py-2", icon("hourglass-half"), " ", tr("heat_aviso")),
        downloadButton("dl_heat", tr("heat_dl"), class = "btn-primary btn-sm")
      ),
      card(
        uiOutput("heat_leyenda", class = "small mb-1"),
        uiOutput("heat_ui"),
        card_footer(caja_leccion("intro_RNAseq")))
    )
  ),

  nav_panel(tr("tab_deg"), value = "deg",
    h3(tr("deg_t")), p(tr("deg_intro")),
    # El panel de ajustes ocupa todo el alto (afecta a las figuras y a la tabla) y acompaña al desplazarse
    layout_columns(col_widths = c(3, 9),
      card(class = "panel-ajustes", div(class = "ajustes-fijos",
        h5(icon("sliders"), " ", tr("deg_ajustes")),
        uiOutput("deg_leccion"),
        selectInput("deg_familia", tr("deg_familia"), FAMILIAS),
        selectInput("deg_cmp", tr("deg_cmp"), NULL),
        sliderInput("deg_lfc", tr("deg_lfc"), 0, 5, 1, step = 0.25),
        selectInput("deg_padj", tr("deg_padj"), c(0.05, 0.01, 0.001, 0.1), 0.05),
        p(class = "small text-muted", tr("deg_nota_umbral")),
        hr(class = "my-2"),
        checkboxInput("deg_etiquetas", tr("deg_etiquetas"), TRUE),
        conditionalPanel("input.deg_etiquetas",
          sliderInput("deg_top", tr("deg_top"), 5, 30, 10, step = 5),
          uiOutput("deg_top_info")),
        hr(class = "my-2"),
        p(class = "small mb-1", HTML(tr("deg_buscar_ayuda"))),
        selectizeInput("deg_genes", tr("deg_buscar"), NULL, multiple = TRUE, options = list(placeholder = tr("deg_buscar_ph"))),
        uiOutput("deg_genes_info")
      )),
      div(
        uiOutput("deg_desc"),
        layout_columns(col_widths = c(7, 5),
          card(card_header(tr("deg_volcano_t")), plotlyOutput("deg_volcano", height = 500), boton_png("deg_volcano", "volcano")),
          card(card_header(tr("deg_resumen_t")), p(class = "small text-muted mb-0", tr("deg_resumen_nota")),
               plotlyOutput("deg_resumen", height = 480), boton_png("deg_resumen", "numero_DEG"))),
        card(class = "interpretacion", card_header(icon("lightbulb"), " ", tr("interp_t")),
             layout_columns(col_widths = c(6, 6), uiOutput("deg_diagrama"), uiOutput("deg_interp"))),
        # Dos columnas: ¿cuántos genes cambian? (Figura 4B) + interpretación | heatmap de log2FC
        layout_columns(col_widths = c(5, 7),
          div(
            card(card_header(tr("prop_t")), p(class = "small text-muted mb-0", HTML(tr("prop_nota"))),
                 uiOutput("deg_prop_leyenda", class = "small"),
                 plotlyOutput("deg_prop", height = 300), boton_png("deg_prop", "proporcion_DEG")),
            card(class = "interpretacion", card_header(icon("lightbulb"), " ", tr("lfch_interp_t")), uiOutput("deg_lfch_interp"))
          ),
          card(card_header(tr("lfch_t")), full_screen = TRUE,
            layout_columns(col_widths = c(6, 6),
              selectInput("lfch_n", tr("lfch_n"), c(25, 50, 100, 200, 500), 50),
              checkboxInput("lfch_bio", tr("lfch_bio"), FALSE)),
            uiOutput("deg_lfch_leyenda", class = "small mb-1"),
            uiOutput("deg_lfch_ui"),
            boton_png("deg_lfch", "heatmap_log2FC"))
        ),
        card(card_header(tr("bio_t")),
          layout_columns(col_widths = c(7, 5),
            div(uiOutput("deg_bio_leyenda", class = "small"), plotlyOutput("deg_bio", height = 300), boton_png("deg_bio", "biotipos")),
            div(class = "interpretacion ps-3", h6(icon("lightbulb"), " ", tr("bio_interp_t")), uiOutput("deg_bio_interp")))),
        card(card_header(tr("deg_tabla_t")), full_screen = TRUE, fill = FALSE,
          div(class = "d-flex gap-2 flex-wrap align-items-center",
              downloadButton("dl_cmp", tr("deg_dl_todos"), class = "btn-primary btn-sm"),
              downloadButton("dl_deg", tr("deg_dl_deg"), class = "btn-primary btn-sm"),
              checkboxInput("deg_todos", tr("deg_todos"), width = "auto")),
          card_body(fillable = FALSE, class = "pt-0", DTOutput("deg_tabla")))
      )
    )
  ),

  nav_panel(tr("tab_go"), value = "go",
    h3(tr("go_t")), p(tr("go_intro")),
    # Mismo esquema que Expresión diferencial: panel de ajustes a todo lo alto y fijo al desplazarse; figura y tabla a la derecha
    layout_columns(col_widths = c(3, 9),
      card(class = "panel-ajustes", div(class = "ajustes-fijos",
        h5(icon("sliders"), " ", tr("go_ajustes")),
        uiOutput("go_leccion"),
        selectInput("go_set", tr("go_conjunto"), setNames(GO_SETS$id, GO_SETS$etiqueta)),
        selectInput("go_lista", tr("go_lista"), NULL),
        p(class = "small text-muted", HTML(tr("go_listas_ayuda"))),
        checkboxGroupInput("go_fuente", tr("go_fuente"), NULL),
        uiOutput("go_fuentes_nota"),
        radioButtons("go_modo", tr("go_modo"), tr("go_modo_op")),
        sliderInput("go_top", tr("go_top"), 5, 50, 10, step = 5)
      )),
      div(
        card(card_header(tr("go_manhattan_t")), p(class = "small text-muted mb-0", HTML(tr("go_manhattan_nota"))),
             plotlyOutput("go_manhattan", height = 420), boton_png("go_manhattan", "manhattan", "go_lista")),
        layout_columns(col_widths = c(5, 7),
          card(card_header(tr("go_fig_t")), uiOutput("go_dot_leyenda", class = "small"), uiOutput("go_plot_ui"), boton_png("go_plot", "dotplot", "go_lista")),
          card(card_header(tr("go_heat_t")), full_screen = TRUE,
            p(class = "small text-muted mb-1", HTML(tr("go_heat_nota"))),
            selectInput("go_heat_fuente", tr("go_heat_fuente"), NULL, width = "100%"),
            uiOutput("go_heat_leyenda", class = "small mb-1"),
            uiOutput("go_heat_ui"), boton_png("go_heat", "terminos_compartidos", "go_set"))),
        card(card_header(tr("go_tabla_t")), full_screen = TRUE, fill = FALSE,
          div(downloadButton("dl_go", tr("go_dl"), class = "btn-primary btn-sm")),
          card_body(fillable = FALSE, class = "pt-0", DTOutput("go_tabla")))
      )
    )
  ),

  nav_panel(tr("tab_descargas"), value = "descargas",
    h3(tr("dl_t")), p(HTML(tr("dl_intro"))),
    layout_columns(col_widths = c(4, 8),
      card(
        selectInput("dl_obj", tr("dl_objeto"), DESCARGABLES),
        downloadButton("dl_csv", tr("dl_boton"), class = "btn-primary"),
        uiOutput("dl_info"),
        hr(), p(tr("dl_r")),
        pre("# install.packages(\"remotes\")\nremotes::install_github(\"EveliaCoss/GeneExprLupusMx\")\nlibrary(GeneExprLupusMx)\ndata(package = \"GeneExprLupusMx\")"),
        p(class = "small", HTML(tr("dl_dds"))),
        a(href = leccion_url("download"), target = "_blank", icon("book"), " ", tr("leccion_descarga"))
      ),
      card(card_header(tr("dl_vista")), tableOutput("dl_vista"))
    ),
    card(card_header(icon("scale-balanced"), " ", tr("dl_origen_t")),
      div(class = "callout-pending", HTML(tr("dl_citar_aviso"))),
      HTML(tr("dl_origen_tabla")),
      p(class = "small text-muted", HTML(tr("dl_licencia_nota"))),
      h5(tr("dl_cita_datos_t")), tags$blockquote(HTML(tr("cita_articulo"))), tags$blockquote(HTML(tr("cita_zenodo"))),
      h5(tr("dl_cita_paquete_t")), tags$blockquote(HTML(tr("cita_paquete"))))
  ),

  nav_panel(tr("tab_aprende"), value = "aprende",
    h3(tr("aprende_t")), p(tr("aprende_intro")),
    # Cuadrícula de 3 × 3: las 7 lecciones, las funciones y la instalación
    do.call(layout_column_wrap, c(list(width = 1/3, class = "mb-3"), lapply(seq_len(nrow(LECCIONES)), function(i) card(class = "leccion-card",
      card_header(class = "d-flex justify-content-between align-items-start gap-2",
        a(href = leccion_url(LECCIONES$id[i]), target = "_blank", LECCIONES$titulo[i]),
        span(class = "small text-muted text-nowrap", LECCIONES$figura[i])),
      a(href = leccion_url(LECCIONES$id[i]), target = "_blank",
        img(src = img_url(LECCIONES$img[i]), style = "width:100%", alt = LECCIONES$titulo[i], loading = "lazy")),
      claves(LECCIONES$claves[[i]])
    )), list(
      card(class = "leccion-card",
        card_header(a(href = leccion_url("functions"), target = "_blank", tr("funciones_t"))),
        p(tr("funciones_txt")), claves(tr("funciones_claves")),
        pre(class = "mb-0", "library(GeneExprLupusMx)\n?create_heatmap\n?Reduce_GOTerm")),
      card(card_header(tr("aprende_instalar")), instalar_ui())
    ))),
    # Visitar y citar el paquete (las citas de los datos están en "Descarga de los datos")
    card(card_header(tr("cita_paquete_t")),
      layout_columns(col_widths = c(6, 6),
        div(
          p(HTML(tr("aprende_paquete_txt"))),
          div(class = "d-flex gap-3 flex-wrap mb-2",
              a(href = URL_SITE, target = "_blank", icon("book-open"), " ", tr("conoce_sitio")),
              a(href = URL_PKG, target = "_blank", icon("github"), " ", tr("conoce_github"))),
          p(class = "small mb-0", HTML(tr("aprende_datos_cita")),
            a(href = "#", onclick = "Shiny.setInputValue('ir_a', 'descargas', {priority: 'event'}); window.scrollTo(0, 0); return false;",
              tr("tab_descargas")), ".")),
        div(
          tags$blockquote(class = "mt-0", HTML(tr("cita_paquete"))),
          div(class = "callout-pending mb-0", HTML(tr("cita_pendiente"))))))
  ),

  # Enlace al repositorio de la app, a la derecha del menú
  nav_spacer(),
  nav_item(tags$a(href = URL_APP_REPO, target = "_blank", class = "nav-link github-link", title = tr("github_app"),
                  `aria-label` = tr("github_app"), icon("github", class = "fa-lg")))
)

# ---- Server ---------------------------------------------------------------
server <- function(input, output, session) {

  # Paneles del inicio
  observeEvent(input$ir_a, nav_select("nav", input$ir_a))

  # Participantes
  output$part_plot <- renderPlotly({
    v <- obj("volunteer_info"); v <- v[!is.na(v[[input$part_var]]), ]
    plot_ly(v, x = ~Group, y = v[[input$part_var]], color = ~Group, colors = COL_GRUPO,
            type = "box", boxpoints = "all", jitter = 0.4, pointpos = 0, text = ~Volunteer_ID) |>
      layout(yaxis = list(title = names(tr("part_vars"))[tr("part_vars") == input$part_var]), xaxis = list(title = ""), showlegend = FALSE)
  })
  output$part_resumen <- renderTable({
    v <- obj("volunteer_info")
    r <- do.call(rbind, lapply(split(v[[input$part_var]], v$Group), function(g) {
      g <- g[!is.na(g)]
      if (length(g)) c(length(g), mean(g), median(g), sd(g), min(g), max(g))
    }))
    colnames(r) <- tr("part_cols")
    r
  }, rownames = TRUE, digits = 1)
  output$part_tabla <- renderDT(datatable(obj("volunteer_info"), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)))

  # Exploración: PCA y gen
  pca <- reactive({
    tipos <- paste(intersect(c("monocyte", "moDC", "tolDC"), input$pca_celulas), collapse = "+")
    req(nzchar(tipos), input$pca_ngenes)
    key <- paste(tipos, input$pca_ngenes, sep = "|")
    r <- PCA[[key]]
    list(x = cbind(as.data.frame(r$x), obj("metadata_cleaned")[rownames(r$x), ]), var = r$variance)
  })
  # Opciones de PC según las muestras elegidas; se conserva la selección si sigue siendo válida
  observeEvent(pca(), {
    v <- pca()$var
    ch <- setNames(seq_along(v), sprintf("PC%d (%s%%)", seq_along(v), v))
    sel <- function(x, def) if (isTruthy(x) && as.integer(x) <= length(v)) x else def
    updateSelectInput(session, "pca_x", choices = ch, selected = sel(input$pca_x, 1))
    updateSelectInput(session, "pca_y", choices = ch, selected = sel(input$pca_y, 2))
  })
  output$scree_plot <- renderPlotly({
    v <- pca()$var; d <- data.frame(pc = seq_along(v), var = v)
    elegidos <- d[d$pc %in% as.integer(c(input$pca_x, input$pca_y)), ]
    plot_ly(d, x = ~pc, y = ~var, type = "scatter", mode = "lines+markers", name = "PC",
            line = list(color = "steelblue"), marker = list(color = "steelblue", size = 7),
            text = ~sprintf("PC%d: %s%%", pc, var), hoverinfo = "text") |>
      add_markers(data = elegidos, x = ~pc, y = ~var, inherit = FALSE, name = tr("scree_sel"),
                  marker = list(color = "darkred", size = 13, symbol = "circle-open", line = list(width = 3)), hoverinfo = "skip") |>
      # etiquetas solo en PC1-PC4, en columna a la derecha del codo para que no se encimen
      add_annotations(data = head(d, 4), x = ~pc, y = ~var, text = ~sprintf("<b>PC%d: %s%%</b>", pc, var),
                      font = list(color = "darkred", size = 13), showarrow = TRUE, arrowcolor = "darkred", arrowhead = 0,
                      axref = "x", ayref = "y", ax = 3 + nrow(d) * 0.12, ay = max(v) * c(0.95, 0.78, 0.61, 0.44)[seq_len(min(4, nrow(d)))],
                      xanchor = "left", inherit = FALSE) |>
      layout(xaxis = list(title = tr("scree_x")), yaxis = list(title = tr("scree_y")), showlegend = FALSE)
  })
  output$pca_plot <- renderPlotly({
    d <- pca(); req(input$pca_x, input$pca_y)
    px <- paste0("PC", input$pca_x); py <- paste0("PC", input$pca_y)
    req(px %in% names(d$x), py %in% names(d$x))
    col <- switch(input$pca_color, Cell_type = COL_CELULA, Group = COL_GRUPO, COL_PALETA)
    plot_ly(d$x, x = d$x[[px]], y = d$x[[py]], color = d$x[[input$pca_color]], colors = col, type = "scatter", mode = "markers",
            marker = list(size = 10, opacity = 0.8), text = ~paste(sample_ID, Group, Cell_type, sep = " · "), hoverinfo = "text") |>
      layout(xaxis = list(title = sprintf("%s (%s%%)", px, d$var[as.integer(input$pca_x)])),
             yaxis = list(title = sprintf("%s (%s%%)", py, d$var[as.integer(input$pca_y)])))
  })
  updateSelectizeInput(session, "gen", choices = rownames(obj("zscore_mat")), selected = "IFI27", server = TRUE)
  output$gen_plot <- renderPlotly({
    req(input$gen %in% rownames(obj("zscore_mat")))
    m <- obj("metadata_cleaned")
    d <- data.frame(z = obj("zscore_mat")[input$gen, ], g = agrupar(m, input$gen_por), Cell_type = m$Cell_type, id = m$sample_ID)
    plot_ly(d, x = ~Cell_type, y = ~z, color = ~g, colors = PALETAS[[input$gen_por]],
            type = "box", boxpoints = "all", jitter = 0.4, pointpos = 0, text = ~id) |>
      layout(boxmode = "group", title = input$gen, yaxis = list(title = "z-score (VST)"), xaxis = list(title = ""))
  })

  # Heatmap de varios genes (z-scores de la VST, como en la lección 2)
  conjuntos <- local({
    r <- rownames(obj("zscore_mat"))
    isg <- unique(obj("isg_genes")$genes)
    l <- list(
      isg10 = c("IFI27", "IFI44L", "IFI44", "USP18", "IFI6", "ISG15", "MX1", "RSAD2", "IFIT1", "OAS1"),  # isg_extended de la lección
      isg_deg = intersect(.deg_sle, isg), deg = .deg_sle, isg = isg,
      sle = unique(obj("sle_genes")$Gene_name), toldc = unique(obj("tolDC_genes")$gene),
      dc = unique(obj("dc_genes")$gene), modc = unique(obj("moDC_genes")$gene))
    lapply(l, function(g) g[g %in% r])
  })
  updateSelectInput(session, "heat_set", choices = c(
    setNames(names(conjuntos), sprintf("%s (%d)", vapply(paste0("heat_", names(conjuntos)), tr, ""), lengths(conjuntos))),
    setNames("custom", tr("heat_custom"))))
  updateSelectizeInput(session, "heat_genes", choices = rownames(obj("zscore_mat")), selected = conjuntos$isg10, server = TRUE)

  # Dot plot de biomarcadores (lección 3, Figura 2B)
  updateSelectInput(session, "dot_set", choices = c(
    setNames("marcadores", sprintf("%s (%d)", tr("dot_marcadores"), length(unlist(MARCADORES)))),
    setNames(names(conjuntos), sprintf("%s (%d)", vapply(paste0("heat_", names(conjuntos)), tr, ""), lengths(conjuntos))),
    setNames("custom", tr("heat_custom"))))
  updateSelectizeInput(session, "dot_genes", choices = rownames(obj("zscore_mat")), selected = c("CD14", "CD1C", "CD163", "IFI27"), server = TRUE)
  dot_df <- reactive({
    g <- switch(req(input$dot_set), marcadores = unlist(MARCADORES, use.names = FALSE), custom = input$dot_genes, conjuntos[[input$dot_set]])
    validate(need(length(g) >= 1, tr("dot_min")))
    m <- obj("metadata_cleaned"); cond <- condicion(m, req(input$dot_por))
    z <- obj("zscore_mat")[g, , drop = FALSE]
    do.call(rbind, lapply(levels(droplevels(cond)), function(k) {
      zk <- z[, cond == k, drop = FALSE]
      data.frame(gene = g, cond = sprintf("%s (n=%d)", k, ncol(zk)), media = rowMeans(zk), pct = 100 * rowMeans(zk > 0))
    }))
  })
  output$dot_ui <- renderUI(plotlyOutput("dot_plot", height = sprintf("%dpx", 170 + 34 * length(unique(dot_df()$cond)))))
  output$dot_plot <- renderPlotly({
    d <- dot_df(); g <- unique(d$gene); filas <- unique(d$cond)
    p <- plot_ly(d, x = ~gene, y = ~cond, type = "scatter", mode = "markers",
                 text = ~sprintf("%s · %s<br>%s: %.2f<br>%s: %.0f%%", gene, cond, tr("dot_color"), media, tr("dot_pct"), pct), hoverinfo = "text",
                 marker = list(color = ~pmin(pmax(media, -2), 2), cmin = -2, cmax = 2, sizemode = "diameter", size = ~4 + 0.16 * pct,
                               line = list(width = 0), opacity = 0.85, showscale = TRUE,
                               colorscale = list(c(0, "darkblue"), c(0.25, "dodgerblue"), c(0.5, "white"), c(0.75, "red"), c(1, "darkred")),
                               colorbar = list(title = tr("dot_color_corto"), len = 0.6))) |>
      layout(xaxis = list(title = "", categoryorder = "array", categoryarray = g, tickangle = -45),
             yaxis = list(title = "", categoryorder = "array", categoryarray = rev(filas)))
    if (identical(input$dot_set, "marcadores")) {   # franjas de color con el tipo celular de cada grupo de biomarcadores
      fin <- cumsum(lengths(MARCADORES)); ini <- fin - lengths(MARCADORES) + 1
      p <- layout(p,
        shapes = lapply(seq_along(MARCADORES), function(i) list(type = "rect", xref = "x", yref = "paper", x0 = ini[i] - 1.45, x1 = fin[i] - 0.55,
                                                                y0 = 1.02, y1 = 1.12, fillcolor = COL_CELULA[[names(MARCADORES)[i]]], line = list(width = 0))),
        annotations = lapply(seq_along(MARCADORES), function(i) list(x = (ini[i] + fin[i]) / 2 - 1, y = 1.07, xref = "x", yref = "paper", showarrow = FALSE, yanchor = "middle",
                                                                     text = sprintf("<b>%s</b>", names(MARCADORES)[i]),
                                                                     font = list(color = if (names(MARCADORES)[i] == "tolDC") "black" else "white"))),
        margin = list(t = 70))
    }
    p
  })

  heat_mat <- reactive({
    g <- if (identical(input$heat_set, "custom")) input$heat_genes else conjuntos[[req(input$heat_set)]]
    validate(need(length(g) >= 2, tr("heat_min")))
    z <- obj("zscore_mat")[g, , drop = FALSE]
    if (isTRUE(input$heat_cluster)) z <- z[hclust(dist(z))$order, , drop = FALSE]
    z
  })
  output$heat_n <- renderUI(p(class = "small text-muted", sprintf(tr("heat_n"), nrow(heat_mat()))))
  alto_genes <- function(n) max(260, min(1400, 14 * n))               # px para las filas de genes
  output$heat_ui <- renderUI(plotlyOutput("heat_plot", height = sprintf("%dpx", 100 + 18 * length(input$heat_pistas) + alto_genes(nrow(heat_mat())))))
  # Barras de anotación elegidas: valores y colores de cada una
  pistas <- reactive({
    m <- obj("metadata_cleaned")
    lapply(setNames(nm = input$heat_pistas), function(k)
      if (k == "celula") list(v = m$Cell_type, col = COL_CELULA) else list(v = agrupar(m, k), col = PALETAS[[k]]))
  })
  output$heat_leyenda <- renderUI({
    nombres <- setNames(names(tr("pistas_op")), tr("pistas_op"))
    div(class = "d-flex flex-wrap column-gap-4 row-gap-1", lapply(names(pistas()), function(k) leyenda(pistas()[[k]]$col, nombres[[k]])))
  })
  output$heat_plot <- renderPlotly({
    z <- heat_mat(); m <- obj("metadata_cleaned")
    # muestras ordenadas por tipo celular y la agrupación elegida, con una columna vacía entre bloques (como facet_grid de la lección)
    bloques <- split(seq_len(ncol(z)), droplevels(condicion(m, req(input$heat_por))))
    idx <- unlist(lapply(bloques, function(b) c(b, NA)))
    idx <- idx[-length(idx)]
    x <- ifelse(is.na(idx), strrep(" ", cumsum(is.na(idx))), as.character(m$sample_ID[idx]))
    zz <- pmin(pmax(z[, idx], -2), 2)                                  # limitar a [-2, 2] como en la lección
    grp <- agrupar(m, input$heat_por)
    txt <- matrix(sprintf("%s · %s<br>%s · %s<br>z = %.2f", rep(rownames(z), length(idx)), rep(m$sample_ID[idx], each = nrow(z)),
                          rep(m$Cell_type[idx], each = nrow(z)), rep(grp[idx], each = nrow(z)), z[, idx]), nrow(z))
    nombres <- setNames(names(tr("pistas_op")), tr("pistas_op"))
    pista <- function(k) {
      v <- pistas()[[k]]$v[idx]; colores <- pistas()[[k]]$col
      plot_ly(x = x, y = nombres[[k]], z = matrix(as.integer(v), 1), type = "heatmap",
        colorscale = do.call(rbind, lapply(seq_along(colores), function(i) rbind(c((i - 1) / length(colores), colores[[i]]), c(i / length(colores), colores[[i]])))),
        zmin = 0.5, zmax = length(colores) + 0.5, showscale = FALSE, text = matrix(as.character(v), 1), hoverinfo = "text")
    }
    principal <- plot_ly(x = x, y = rownames(z), z = zz, type = "heatmap", zmin = -2, zmax = 2,
                         colorscale = list(c(0, "steelblue"), c(0.5, "white"), c(1, "firebrick")),
                         text = txt, hoverinfo = "text", colorbar = list(title = "z-score", len = 0.5))
    k <- length(pistas())
    p <- if (k == 0) principal else
      do.call(subplot, c(lapply(names(pistas()), pista), list(principal), list(nrows = k + 1, shareX = TRUE, margin = 0,
                         heights = c(rep(18, k), alto_genes(nrow(z))) / (18 * k + alto_genes(nrow(z))))))
    ejes <- setNames(list(list(showticklabels = nrow(z) <= 120, tickfont = list(size = 10), autorange = "reversed")),
                     if (k == 0) "yaxis" else paste0("yaxis", k + 1))
    do.call(layout, c(list(p, xaxis = list(showticklabels = FALSE, title = tr("heat_x"))), ejes))
  })
  output$dl_heat <- downloadHandler(
    filename = function() sprintf("GeneExprLupusMx_zscore_%s.csv", input$heat_set),
    content = function(file) escribir_csv(heat_mat(), file)
  )

  # Expresión diferencial
  observeEvent(input$deg_familia, {
    c2 <- CMP[CMP$familia == input$deg_familia, ]
    updateSelectInput(session, "deg_cmp", choices = setNames(c2$id, c2$etiqueta))
    updateSliderInput(session, "deg_lfc", value = if (input$deg_familia == "int") 0 else 1)
  })
  cmp <- reactive({ req(input$deg_cmp %in% CMP$id); CMP[CMP$id == input$deg_cmp, ] })
  deg_raw <- reactive(obj(cmp()$objeto)[[cmp()$elemento]])
  updateSelectizeInput(session, "deg_genes", choices = obj("de_results")$monocyte$gene, server = TRUE)

  deg <- reactive({
    d <- deg_raw(); a <- anotaciones()
    d$direccion <- direccion(d, input$deg_lfc, as.numeric(input$deg_padj))
    d$biotipo <- unname(a$biotipo[d$gene])
    d$ISG <- ifelse(d$gene %in% a$isg, tr("si"), "")
    d$SLE <- ifelse(d$gene %in% a$sle, tr("si"), "")
    d$tolDC <- unname(a$toldc[d$gene])
    d
  })

  output$deg_desc <- renderUI({
    r <- cmp()
    modelo <- switch(r$familia,
      sle = sprintf(tr("modelo_sle"), tr(paste0("cel_", r$elemento))),
      dif = { p <- strsplit(r$elemento, "_")[[1]]; sprintf(tr(paste0("modelo_dif_", p[1])), p[2]) },
      int = sprintf(tr("modelo_int"), r$elemento, r$elemento))
    div(class = "alert alert-secondary",
      HTML(tr(paste0("fam_", r$familia, "_desc"))),
      div(class = "mt-2", HTML(modelo)),
      if (r$familia == "int") div(class = "mt-2", HTML(tr("deg_nota_interaccion"))))
  })
  output$deg_leccion <- renderUI(div(class = "small mb-3", caja_leccion(cmp()$leccion)))

  # Comparaciones del tipo elegido, con sus anotaciones (grupo / tipo celular / comparación) y conteos UP/DOWN
  familia_info <- reactive({
    c2 <- CMP[CMP$familia == cmp()$familia, ]
    if (c2$familia[1] == "dif") c2 <- c2[order(sub("_.*", "", c2$elemento), match(sub(".*_", "", c2$elemento), c("moDC", "tolDC"))), ]
    lfc <- input$deg_lfc; padj <- as.numeric(input$deg_padj)
    c2$celula <- sub(".*_", "", c2$elemento)
    c2$corto <- switch(c2$familia[1], sle = c2$elemento, dif = c2$elemento, int = paste0("SLEvsCtrl_", c2$elemento))
    c2$grupo <- switch(c2$familia[1], sle = "SLE vs Ctrl", dif = sub("_.*", "", c2$elemento), int = c2$corto)
    dirs <- lapply(seq_len(nrow(c2)), function(i) direccion(obj(c2$objeto[i])[[c2$elemento[i]]], lfc, padj))
    c2$up <- vapply(dirs, function(d) sum(d == "UP"), 1L); c2$down <- vapply(dirs, function(d) sum(d == "DOWN"), 1L)
    c2
  })
  colores_grupo <- function(f) switch(f$familia[1], sle = c("SLE vs Ctrl" = "#B493B4"), dif = COL_GRUPO, int = COL_COMPARACION)
  nombre_grupo <- function(f) if (f$familia[1] == "int") tr("go_heat_comparacion") else tr("heat_pista_grupo2")
  # pista de color categórica (una fila u horizontal)
  pista_cat <- function(x, y, v, colores, texto, vertical = FALSE) {
    codigo <- match(v, names(colores))
    plot_ly(x = x, y = y, z = if (vertical) matrix(codigo, ncol = 1) else matrix(codigo, 1), type = "heatmap",
            colorscale = do.call(rbind, lapply(seq_along(colores), function(i) rbind(c((i - 1) / length(colores), colores[[i]]), c(i / length(colores), colores[[i]])))),
            zmin = 0.5, zmax = length(colores) + 0.5, showscale = FALSE, xgap = 1, ygap = 1,
            text = if (vertical) matrix(texto, ncol = 1) else matrix(texto, 1), hoverinfo = "text")
  }

  # ¿Cuántos genes cambian? Proporción UP/DOWN con el número dentro (lección 5, Figura 4B)
  output$deg_prop <- renderPlotly({
    f <- familia_info(); tot <- pmax(f$up + f$down, 1)
    et <- ifelse(f$id == input$deg_cmp, sprintf("<b>%s</b>", f$corto), f$corto)
    lv <- rev(et)
    barras <- plot_ly(y = factor(et, lv), x = f$up / tot, type = "bar", orientation = "h", name = "UP",
                      marker = list(color = COLORES_DIR[["UP"]]), text = f$up, textposition = "inside", insidetextanchor = "middle",
                      textfont = list(color = "white"), hovertext = sprintf("%s · UP: %d", f$corto, f$up), hoverinfo = "text") |>
      add_trace(x = f$down / tot, name = "DOWN", marker = list(color = COLORES_DIR[["DOWN"]]), text = f$down,
                hovertext = sprintf("%s · DOWN: %d", f$corto, f$down))
    subplot(
      pista_cat(nombre_grupo(f), factor(et, lv), f$grupo, colores_grupo(f), f$grupo, vertical = TRUE),
      pista_cat(tr("heat_pista_celula2"), factor(et, lv), f$celula, COL_CELULA, f$celula, vertical = TRUE),
      barras, widths = c(0.08, 0.08, 0.84), shareY = TRUE, margin = 0.005) |>
      layout(barmode = "stack", showlegend = TRUE, legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.12),
             xaxis = list(showticklabels = FALSE), xaxis2 = list(showticklabels = FALSE),
             xaxis3 = list(tickformat = ".0%", range = c(0, 1), title = ""), yaxis = list(title = ""))
  })

  # Heatmap de log2FC: DEG de la comparación elegida (filas) en todas las comparaciones del tipo (columnas)
  lfch <- reactive({
    f <- familia_info(); d <- deg()
    sig <- d[d$direccion != "NS", ]
    sig <- head(sig[order(sig$padj), ], as.integer(input$lfch_n))
    validate(need(nrow(sig) > 0, tr("deg_top_cero")))
    m <- sapply(seq_len(nrow(f)), function(i) { x <- obj(f$objeto[i])[[f$elemento[i]]]; x$log2FoldChange[match(sig$gene, x$gene)] })
    m <- matrix(m, nrow(sig), dimnames = list(sig$gene, f$corto))
    dir <- factor(sig$direccion, c("UP", "DOWN"))
    o <- order(dir, -abs(sig$log2FoldChange)); m <- m[o, , drop = FALSE]
    list(m = m, dir = dir[o], f = f, biotipo = ifelse(is.na(sig$biotipo[o]), "NA", sig$biotipo[o]))
  })
  # Distribución de biotipos: DEG sobreexpresados, DEG reprimidos y todos los genes analizados (con padj) de la comparación
  bio <- reactive({
    d <- deg(); d <- d[!is.na(d$padj), ]
    b <- ifelse(is.na(d$biotipo), "NA", d$biotipo)
    grupos <- list(UP = b[d$direccion == "UP"], DOWN = b[d$direccion == "DOWN"], fondo = b)
    do.call(rbind, lapply(names(grupos), function(k) {
      n <- table(factor(grupos[[k]], names(COL_BIOTIPO)))
      data.frame(grupo = k, biotipo = names(n), n = as.vector(n), pct = if (sum(n)) 100 * as.vector(n) / sum(n) else 0, total = sum(n))
    }))
  })
  output$deg_bio_leyenda <- renderUI({
    b <- bio(); presentes <- unique(b$biotipo[b$n > 0])
    leyenda(setNames(COL_BIOTIPO[presentes], ifelse(presentes == "NA", tr("bio_na"), presentes)), tr("col_biotipo"))
  })
  output$deg_bio <- renderPlotly({
    b <- bio(); b <- b[b$biotipo %in% unique(b$biotipo[b$n > 0]), ]
    et <- c(UP = tr("bio_up"), DOWN = tr("bio_down"), fondo = tr("bio_fondo"))
    b$fila <- factor(sprintf("%s (n=%s)", et[b$grupo], format(b$total, big.mark = ",", trim = TRUE)), rev(unique(sprintf("%s (n=%s)", et[b$grupo], format(b$total, big.mark = ",", trim = TRUE)))))
    b$nombre <- ifelse(b$biotipo == "NA", tr("bio_na"), b$biotipo)
    plot_ly(b, y = ~fila, x = ~pct, color = ~factor(biotipo, names(COL_BIOTIPO)), colors = COL_BIOTIPO, type = "bar", orientation = "h",
            text = ~ifelse(pct >= 6, sprintf("%.0f%%", pct), ""), textposition = "inside", insidetextanchor = "middle", textfont = list(color = "white"),
            hovertext = ~sprintf("%s · %s<br>%s genes (%.1f%%)", fila, nombre, format(n, big.mark = ","), pct), hoverinfo = "text") |>
      layout(barmode = "stack", showlegend = FALSE, xaxis = list(title = tr("bio_x"), range = c(0, 100), ticksuffix = "%"), yaxis = list(title = ""))
  })
  output$deg_bio_interp <- renderUI({
    b <- bio(); r <- cmp()
    deg_b <- aggregate(n ~ biotipo, b[b$grupo != "fondo", ], sum); tot <- sum(deg_b$n)
    if (tot == 0) return(p(tr("deg_top_cero")))
    deg_b$pct <- 100 * deg_b$n / tot
    fondo <- b[b$grupo == "fondo", ]; deg_b$fondo <- fondo$pct[match(deg_b$biotipo, fondo$biotipo)]
    deg_b <- deg_b[deg_b$n > 0, ]; deg_b <- deg_b[order(-deg_b$n), ]
    nombre <- function(x) ifelse(x == "NA", tr("bio_na"), x)
    filas <- sprintf("<li><b>%s</b>: %d DEG (%.1f%%; %.1f%% entre todos los genes analizados)</li>", nombre(deg_b$biotipo), deg_b$n, deg_b$pct, deg_b$fondo)
    dif <- deg_b$pct - deg_b$fondo; k <- which.max(abs(dif))
    tagList(
      p(HTML(sprintf(tr("bio_interp_1"), tot, r$etiqueta))),
      tags$ul(class = "small", HTML(paste(filas, collapse = ""))),
      p(class = "small", HTML(sprintf(if (dif[k] > 0) tr("bio_interp_mas") else tr("bio_interp_menos"),
                                      nombre(deg_b$biotipo[k]), deg_b$pct[k], deg_b$fondo[k]))),
      p(class = "small mb-0", HTML(tr("bio_interp_nota")))
    )
  })

  output$deg_prop_leyenda <- renderUI({
    f <- familia_info()
    div(class = "d-flex flex-wrap column-gap-3", leyenda(colores_grupo(f), nombre_grupo(f)), leyenda(COL_CELULA[unique(f$celula)], tr("heat_pista_celula2")))
  })
  output$deg_lfch_leyenda <- renderUI({
    f <- familia_info()
    div(class = "d-flex flex-wrap column-gap-3 row-gap-1",
        leyenda(colores_grupo(f), nombre_grupo(f)), leyenda(COL_CELULA[unique(f$celula)], tr("heat_pista_celula2")),
        leyenda(COLORES_DIR[c("UP", "DOWN")], tr("lfch_dir")),
        if (isTRUE(input$lfch_bio)) {
          b <- intersect(names(COL_BIOTIPO), lfch()$biotipo)
          leyenda(setNames(COL_BIOTIPO[b], ifelse(b == "NA", tr("bio_na"), b)), tr("col_biotipo"))
        })
  })
  output$deg_lfch_ui <- renderUI(plotlyOutput("deg_lfch", height = sprintf("%dpx", 230 + min(1300, max(260, 12 * nrow(lfch()$m))))))
  output$deg_lfch <- renderPlotly({
    h <- lfch(); m <- h$m; f <- h$f; L <- 2   # escala fija de -2 a 2
    # huecos entre grupos (column_split) y entre UP/DOWN (row_split)
    xs <- seq_len(ncol(m)) + cumsum(c(0, f$grupo[-1] != f$grupo[-nrow(f)]))
    ys <- seq_len(nrow(m)) + cumsum(c(0, as.character(h$dir[-1]) != as.character(h$dir[-length(h$dir)])))
    X <- seq_len(max(xs)); Y <- seq_len(max(ys))
    z <- matrix(NA_real_, length(Y), length(X)); z[ys, xs] <- pmin(pmax(m, -L), L)
    txt <- matrix("", length(Y), length(X))
    txt[ys, xs] <- sprintf("%s · %s<br>log2FC = %.2f", rep(rownames(m), ncol(m)), rep(colnames(m), each = nrow(m)), m)
    # colores de la lección 5: darkblue, dodgerblue, white, red, darkred en -2, -1, 0, 1, 2
    p1 <- (L - 1) / (2 * L)
    escala <- list(c(0, "darkblue"), c(p1, "dodgerblue"), c(0.5, "white"), c(1 - p1, "red"), c(1, "darkred"))
    principal <- plot_ly(x = X, y = Y, z = z, type = "heatmap", zmin = -L, zmax = L, colorscale = escala, xgap = 2, ygap = if (nrow(m) <= 100) 1 else 0,
                         text = txt, hoverinfo = "text", colorbar = list(title = "log2FC", len = 0.35, thickness = 12, y = 1, yanchor = "top"))
    xbio <- max(X) + 1.5   # columna de biotipo junto al heatmap, con un pequeño hueco
    if (isTRUE(input$lfch_bio)) {   # columna de biotipo a la derecha, como anotación de filas
      cb <- COL_BIOTIPO[intersect(names(COL_BIOTIPO), h$biotipo)]
      zb <- rep(NA_real_, length(Y)); zb[ys] <- match(h$biotipo, names(cb))
      tb <- rep("", length(Y)); tb[ys] <- sprintf("%s · %s", rownames(m), ifelse(h$biotipo == "NA", tr("bio_na"), h$biotipo))
      principal <- principal |> add_heatmap(x = xbio, y = Y, z = matrix(zb, ncol = 1), inherit = FALSE, showscale = FALSE, ygap = if (nrow(m) <= 100) 1 else 0,
        colorscale = do.call(rbind, lapply(seq_along(cb), function(i) rbind(c((i - 1) / length(cb), cb[[i]]), c(i / length(cb), cb[[i]])))),
        zmin = 0.5, zmax = length(cb) + 0.5, text = matrix(tb, ncol = 1), hoverinfo = "text")
    }
    xpos <- function(v) { out <- rep(NA, length(X)); out[xs] <- v; out }
    conteo <- plot_ly(x = xs, y = f$up, type = "bar", name = "UP", marker = list(color = COLORES_DIR[["UP"]]),
                      hovertext = sprintf("%s · UP: %d", f$corto, f$up), hoverinfo = "text") |>
      add_trace(y = f$down, name = "DOWN", marker = list(color = COLORES_DIR[["DOWN"]]), hovertext = sprintf("%s · DOWN: %d", f$corto, f$down))
    alto <- min(1300, max(260, 12 * nrow(m)))
    pg <- pista_cat(X, nombre_grupo(f), xpos(f$grupo), colores_grupo(f), xpos(f$corto))
    pc <- pista_cat(X, tr("heat_pista_celula2"), xpos(f$celula), COL_CELULA, xpos(f$corto))
    ndir <- table(h$dir)
    subplot(conteo, pg, pc, principal, nrows = 4, shareX = TRUE, margin = 0.004,
            heights = c(60, 16, 16, alto) / (92 + alto)) |>
      layout(barmode = "stack", showlegend = FALSE, plot_bgcolor = "white",
             xaxis = list(range = c(0.5, if (isTRUE(input$lfch_bio)) xbio + 0.5 else max(X) + 0.5),
                          tickvals = c(xs, if (isTRUE(input$lfch_bio)) xbio), ticktext = c(f$corto, if (isTRUE(input$lfch_bio)) tr("col_biotipo")),
                          tickangle = -30, side = "bottom", zeroline = FALSE),
             yaxis = list(title = "DEG", tickfont = list(size = 9)),
             yaxis4 = list(tickvals = ys, ticktext = rownames(m), showticklabels = nrow(m) <= 100, tickfont = list(size = 9),
                           range = c(max(Y) + 0.5, 0.5), zeroline = FALSE, showgrid = FALSE),
             annotations = lapply(levels(h$dir)[ndir[levels(h$dir)] > 0], function(k) {
               yy <- ys[h$dir == k]
               list(text = sprintf("<b>%s (%d)</b>", k, ndir[[k]]), x = 1.01, xref = "paper", xanchor = "left",
                    y = mean(range(yy)), yref = "y4", showarrow = FALSE, font = list(color = COLORES_DIR[[k]], size = 11), textangle = 0)
             }))
  })
  output$deg_lfch_interp <- renderUI({
    h <- lfch(); m <- h$m; f <- h$f; r <- cmp()
    signo <- sign(m); misma <- sum(apply(signo, 1, function(x) length(unique(x[!is.na(x)])) == 1))
    tagList(
      p(HTML(sprintf(tr("lfch_interp_1"), nrow(m), r$etiqueta, sum(h$dir == "UP"), sum(h$dir == "DOWN"), paste(f$corto, collapse = ", ")))),
      p(HTML(sprintf(tr("lfch_interp_2"), misma, nrow(m), round(100 * misma / nrow(m)), nrow(m) - misma))),
      p(class = "small mb-0", HTML(tr(paste0("lfch_interp_", r$familia))))
    )
  })

  # Diagrama de la comparación, con la configuración de la lección 5 (Figura 4A): flowchart LR y colores por tipo celular
  output$deg_diagrama <- renderUI({
    r <- cmp()
    fondo <- c(monocyte = "#000000", moDC = "#db7093", tolDC = "#ffd700")
    letra <- c(monocyte = "#ffffff", moDC = "#ffffff", tolDC = "#000000")
    estilo <- function(id, cel) sprintf("  style %s fill:%s,color:%s,stroke:%s", id, fondo[[cel]], letra[[cel]], fondo[[cel]])
    texto <- switch(r$familia,
      sle = {
        cel <- r$elemento
        c("flowchart LR",
          sprintf('  S(["%s<br/>SLE"]) -- "vs" --- R(["%s<br/>Ctrl · %s"])', cel, cel, tr("diag_referencia")),
          sprintf("  style R fill:#96D4CC,color:#14164D,stroke:%s,stroke-width:4px", fondo[[cel]]),
          sprintf("  style S fill:#B493B4,color:#14164D,stroke:%s,stroke-width:4px", fondo[[cel]]))
      },
      dif = {
        p <- strsplit(r$elemento, "_")[[1]]; cel <- p[2]
        # RL: el monocito (referencia) queda a la derecha y la flecha de diferenciación apunta a la izquierda
        c("flowchart RL",
          sprintf('  M(["Monocito<br/>%s"]) -- "%s" --> D(["%s<br/>%s"])', tr(paste0("diag_grupo_", p[1])), tr(paste0("diag_citocinas_", cel)),
                  cel, tr(paste0("diag_desc_", cel))),
          estilo("M", "monocyte"), estilo("D", cel))
      },
      int = {
        cel <- r$elemento
        # pacientes arriba y controles (referencia) abajo; en cada una, el monocito (referencia) a la derecha
        c("flowchart TB",
          sprintf('  subgraph S["%s"]', tr("diag_sle")), "    direction RL",
          sprintf('    M2(["Monocito"]) -- "%s" --> X2(["%s"])', tr(paste0("diag_citocinas_", cel)), cel),
          "  end",
          sprintf('  subgraph C["%s"]', tr("diag_ctrl")), "    direction RL",
          sprintf('    M1(["Monocito"]) -- "%s" --> X1(["%s"])', tr(paste0("diag_citocinas_", cel)), cel),
          "  end",
          "  S ~~~ C",  # S antes que C: pacientes arriba
          estilo("M1", "monocyte"), estilo("X1", cel), estilo("M2", "monocyte"), estilo("X2", cel),
          "  style C fill:#E8F6F4,stroke:#96D4CC", "  style S fill:#F3ECF3,stroke:#B493B4",
          sprintf("  linkStyle 0,1 stroke:%s,stroke-width:3px", COL_COMPARACION[[paste0("SLEvsCtrl_", cel)]]))  # no la liga invisible S ~~~ C
      })
    tagList(
      tags$pre(class = "mermaid", paste(texto, collapse = "\n")),
      p(class = "small text-muted text-center mb-0", HTML(tr(paste0("diag_pie_", r$familia)))),
      # dibujar en cuanto Mermaid termine de cargar
      tags$script(HTML("(function dibujar() { if (window.mermaid) mermaid.run({ querySelector: '#deg_diagrama pre.mermaid' }); else setTimeout(dibujar, 200); })();"))
    )
  })

  # Interpretación del volcano: cuántos genes están sobreexpresados y reprimidos, y respecto a qué referencia
  output$deg_interp <- renderUI({
    d <- deg(); r <- cmp()
    up <- sum(d$direccion == "UP"); down <- sum(d$direccion == "DOWN")
    celula <- function(x) tr(paste0("cel_", x))
    texto <- switch(r$familia,
      sle = sprintf(tr("interp_sle"), celula(r$elemento), up, down, celula(r$elemento)),
      dif = {
        p <- strsplit(r$elemento, "_")[[1]]
        sprintf(tr("interp_dif"), tr(paste0("grupo_", p[1])), p[2], up, down, p[2])
      },
      int = sprintf(tr("interp_int"), up, r$elemento, down, r$elemento))
    sig <- d[d$direccion != "NS", ]
    mejor <- if (nrow(sig)) sig[which.min(sig$padj), ]
    tagList(
      p(HTML(texto)),
      p(class = "small mb-0", HTML(sprintf(tr("interp_umbral"), input$deg_lfc, input$deg_padj))),
      if (!is.null(mejor)) p(class = "small mb-0", HTML(sprintf(tr("interp_mejor"), mejor$gene,
        if (mejor$log2FoldChange > 0) tr("interp_sube") else tr("interp_baja"), mejor$log2FoldChange, signif(mejor$padj, 3)))),
      if (up + down == 0) p(class = "small mb-0", tr("deg_top_cero"))
    )
  })

  output$deg_volcano <- renderPlotly({
    d <- deg(); d <- d[!is.na(d$padj), ]
    d$y <- -log10(pmax(d$padj, .Machine$double.xmin))
    n <- table(factor(d$direccion, c("UP", "DOWN", "NS")))
    d$direccion <- factor(d$direccion, c("UP", "DOWN", "NS"),
                          sprintf("%s (%d)", c(tr("deg_up"), tr("deg_down"), tr("deg_ns")), n))
    lfc <- input$deg_lfc; py <- -log10(as.numeric(input$deg_padj))
    resaltar <- d[d$gene %in% input$deg_genes, ]
    # Nombres de los DEG con menor padj (sin repetir los resaltados)
    top <- if (isTRUE(input$deg_etiquetas)) {
      sig <- d[d$direccion != levels(d$direccion)[3] & !d$gene %in% resaltar$gene, ]
      head(sig[order(sig$padj), ], input$deg_top)
    }
    p <- plot_ly(d, x = ~log2FoldChange, y = ~y, color = ~direccion, colors = setNames(COLORES_DIR, levels(d$direccion)),
                 type = "scattergl", mode = "markers", marker = list(size = 5, opacity = 0.6),
                 text = ~sprintf("%s<br>LFC: %.2f<br>padj: %.2g", gene, log2FoldChange, padj), hoverinfo = "text") |>
      layout(title = list(text = if (cmp()$familia == "int") color_texto(cmp()$etiqueta, COL_COMPARACION[[paste0("SLEvsCtrl_", cmp()$elemento)]]) else cmp()$etiqueta,
                          font = list(size = 14)),
             xaxis = list(title = "log2 fold change"), yaxis = list(title = "−log10(padj)"),
             legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.22, yanchor = "top"),
             shapes = list(
               list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = py, y1 = py, line = list(dash = "dot", color = "grey")),
               list(type = "line", x0 = -lfc, x1 = -lfc, y0 = 0, y1 = 1, yref = "paper", line = list(dash = "dot", color = "grey")),
               list(type = "line", x0 = lfc, x1 = lfc, y0 = 0, y1 = 1, yref = "paper", line = list(dash = "dot", color = "grey"))))
    if (NROW(top)) p <- p |> add_annotations(data = top, x = ~log2FoldChange, y = ~y, text = ~gene, inherit = FALSE,
                                             font = list(size = 10, color = "#14164D"), showarrow = TRUE, arrowhead = 0, arrowcolor = "grey",
                                             # en abanico: cada etiqueta con distinta distancia para que no se encimen
                                             ax = sign(top$log2FoldChange + 1e-9) * (30 + 18 * ((seq_len(nrow(top)) - 1) %/% 6)),
                                             ay = -12 - 15 * ((seq_len(nrow(top)) - 1) %% 6))
    if (nrow(resaltar)) p <- p |>
      add_markers(data = resaltar, x = ~log2FoldChange, y = ~y, inherit = FALSE, name = tr("deg_resaltados"), hoverinfo = "skip",
                  marker = list(size = 12, symbol = "circle-open", color = "#7F62B3", line = list(width = 2.5))) |>
      add_annotations(data = resaltar, x = ~log2FoldChange, y = ~y, text = ~sprintf("<b>%s</b>", gene), inherit = FALSE,
                      font = list(size = 12, color = "#7F62B3"), showarrow = TRUE, arrowhead = 0, arrowcolor = "#7F62B3", ax = 30, ay = -30)
    config(p, edits = list(annotationTail = TRUE))   # las etiquetas se pueden arrastrar si se enciman
  })

  # Etiquetas: el máximo del control depende de cuántos DEG pasan los umbrales
  n_sig <- reactive(sum(deg()$direccion != "NS"))
  observeEvent(n_sig(), {
    tope <- max(5, min(30, 5 * ceiling(n_sig() / 5)))
    updateSliderInput(session, "deg_top", max = tope, value = min(isolate(input$deg_top) %||% 10, tope))
  })
  output$deg_top_info <- renderUI(p(class = "small text-muted",
    if (n_sig() == 0) tr("deg_top_cero") else sprintf(tr("deg_top_info"), min(input$deg_top, n_sig()), n_sig())))
  output$deg_genes_info <- renderUI({
    fuera <- setdiff(input$deg_genes, deg()$gene[!is.na(deg()$padj)])
    if (length(fuera)) p(class = "small text-muted", sprintf(tr("deg_genes_fuera"), paste(fuera, collapse = ", ")))
  })

  deg_mostrar <- reactive({
    d <- deg()
    if (!input$deg_todos) d <- d[d$direccion != "NS", ]
    d <- d[order(d$padj), c("gene", "direccion", "log2FoldChange", "lfcSE", "baseMean", "pvalue", "padj", "biotipo", "ISG", "SLE", "tolDC")]
    names(d)[c(1, 2, 8:11)] <- c(tr("col_gen"), tr("col_dir"), tr("col_biotipo"), tr("col_isg"), tr("col_sle"), tr("col_toldc"))
    d
  })
  output$deg_tabla <- renderDT(
    datatable(deg_mostrar(), rownames = FALSE, filter = "top", fillContainer = FALSE,
              options = list(pageLength = 10, lengthMenu = c(10, 25, 50, 100), scrollX = TRUE)) |>
      formatSignif(c("log2FoldChange", "lfcSE", "baseMean", "pvalue", "padj"), 3)
  )
  output$dl_cmp <- downloadHandler(
    filename = function() sprintf("GeneExprLupusMx_%s_%s.csv", cmp()$objeto, cmp()$elemento),
    content = function(file) escribir_csv(deg_raw(), file)
  )
  output$dl_deg <- downloadHandler(
    filename = function() sprintf("GeneExprLupusMx_%s_%s_DEG_lfc%s_padj%s.csv", cmp()$objeto, cmp()$elemento, input$deg_lfc, input$deg_padj),
    content = function(file) { d <- deg(); escribir_csv(d[d$direccion != "NS", ], file) }
  )

  output$deg_resumen <- renderPlotly({
    lfc <- input$deg_lfc; padj <- as.numeric(input$deg_padj)
    n <- do.call(rbind, lapply(seq_len(nrow(CMP)), function(i) {
      dir <- direccion(obj(CMP$objeto[i])[[CMP$elemento[i]]], lfc, padj)
      data.frame(cmp = CMP$etiqueta[i], dir = c("UP", "DOWN"), n = c(sum(dir == "UP"), sum(dir == "DOWN")))
    }))
    # etiquetas en dos líneas para que quepan junto al volcano; la comparación elegida va en negritas
    et <- sub(": ", "<br>", CMP$etiqueta)
    int <- CMP$familia == "int"
    k <- paste0("SLEvsCtrl_", CMP$elemento[int])
    et[int] <- mapply(color_texto, k, COL_COMPARACION[k])   # solo el nombre corto, con su color
    et[CMP$id == input$deg_cmp] <- sprintf("<b>%s</b>", et[CMP$id == input$deg_cmp])
    n$cmp <- factor(et[match(n$cmp, CMP$etiqueta)], rev(et))
    plot_ly(n, y = ~cmp, x = ~n, color = ~dir, colors = COLORES_DIR[c("UP", "DOWN")], type = "bar", orientation = "h",
            text = ~n, textposition = "auto") |>
      layout(barmode = "group", xaxis = list(title = "DEG"), yaxis = list(title = "", tickfont = list(size = 11)),
             legend = list(title = list(text = tr("col_dir")), orientation = "h", y = 1.08))
  })

  # Funciones biológicas
  go_set <- reactive(GO_SETS[GO_SETS$id == input$go_set, ])
  observeEvent(input$go_set, {
    g <- obj(go_set()$objeto)
    updateSelectInput(session, "go_lista", choices = sort(unique(g$query)))
    fuentes <- intersect(names(COL_FUENTES), unique(g$source))
    nombres <- tr("fuentes_nombres")
    etiquetas <- ifelse(nombres[fuentes] == fuentes, fuentes, sprintf("%s (%s)", nombres[fuentes], fuentes))
    updateCheckboxGroupInput(session, "go_fuente", choices = setNames(fuentes, etiquetas),
                             selected = intersect(fuentes, c("GO:BP", "REAC", "KEGG")))
  })
  observeEvent(input$go_set, {
    fuentes <- intersect(names(COL_FUENTES), unique(obj(go_set()$objeto)$source))
    nombres <- tr("fuentes_nombres")
    updateSelectInput(session, "go_heat_fuente", choices = setNames(fuentes, ifelse(nombres[fuentes] == fuentes, fuentes, sprintf("%s (%s)", nombres[fuentes], fuentes))),
                      selected = if ("REAC" %in% fuentes) "REAC" else fuentes[1])
  })
  output$go_fuentes_nota <- renderUI(if (input$go_set != "sle") p(class = "small text-muted", tr("go_fuentes_nota")))
  go_filtrado <- reactive({
    g <- obj(go_set()$objeto)
    g <- g[g$query == input$go_lista & g$source %in% input$go_fuente, ]
    g[order(g$p_value), ]
  })
  output$go_leccion <- renderUI(div(class = "small mb-3", caja_leccion(go_set()$leccion)))
  # Nombre de la lista; en las comparaciones de contrastes, con su color
  titulo_lista <- function(q) {
    k <- sub("_(UP|DOWN)", "", q)
    if (k %in% names(COL_COMPARACION)) color_texto(q, COL_COMPARACION[[k]]) else q
  }

  # Manhattan plot al estilo de g:Profiler: todos los términos significativos de la lista, por base de datos
  output$go_manhattan <- renderPlotly({
    g <- obj(go_set()$objeto); g <- g[g$query == req(input$go_lista), ]
    validate(need(nrow(g) > 0, tr("go_vacio")))
    fuentes <- intersect(names(COL_FUENTES), unique(g$source))
    g$source <- factor(g$source, fuentes)
    g <- g[order(g$source, g$term_id), ]
    # cada base de datos ocupa un bloque de ancho 1; sus términos se reparten dentro del bloque
    g$x <- as.integer(g$source) - 1 + ave(seq_len(nrow(g)), g$source, FUN = function(i) (seq_along(i) - 0.5) / length(i))
    g$y <- pmin(-log10(g$p_value), 16)                                  # tope en 16, como gostplot
    plot_ly(g, x = ~x, y = ~y, color = ~source, colors = COL_FUENTES[fuentes], type = "scatter", mode = "markers",
            marker = list(size = ~5 + 2.2 * log10(term_size), opacity = 0.75, line = list(width = 0)),
            text = ~sprintf("%s<br>%s<br>padj: %.2g<br>%s: %d", term_id, term_name, p_value, tr("go_tam"), intersection_size),
            hoverinfo = "text") |>
      layout(title = list(text = titulo_lista(input$go_lista), font = list(size = 14)), showlegend = FALSE,
             xaxis = list(title = "", tickvals = seq_along(fuentes) - 0.5, ticktext = fuentes, range = c(0, length(fuentes)),
                          zeroline = FALSE, showgrid = FALSE),
             yaxis = list(title = tr("go_manhattan_y"), range = c(0, max(g$y) * 1.1 + 0.5)),
             shapes = c(
               lapply(seq_along(fuentes)[seq_along(fuentes) %% 2 == 0], function(i)
                 list(type = "rect", x0 = i - 1, x1 = i, y0 = 0, y1 = 1, yref = "paper", fillcolor = "#ECEEFE", line = list(width = 0), layer = "below")),
               list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = -log10(0.05), y1 = -log10(0.05), line = list(dash = "dot", color = "grey")))))
  })

  output$go_dot_leyenda <- renderUI({
    fuentes <- intersect(names(COL_FUENTES), unique(go_dot()$source))
    leyenda(COL_FUENTES[fuentes], tr("go_fuente"))
  })
  output$go_plot_ui <- renderUI({
    n <- nrow(go_dot())
    plotlyOutput("go_plot", height = sprintf("%dpx", 160 + 24 * max(n, 5)))
  })
  # Dot plot: los N más significativos en total o de cada base de datos
  go_dot <- reactive({
    g <- go_filtrado()
    if (identical(input$go_modo, "fuente")) g <- do.call(rbind, lapply(split(g, factor(g$source, names(COL_FUENTES))), head, input$go_top))
    else g <- head(g, input$go_top)
    g
  })
  output$go_plot <- renderPlotly({
    g <- go_dot()
    validate(need(nrow(g) > 0, tr("go_vacio")))
    fuentes <- intersect(names(COL_FUENTES), unique(g$source))
    g$term_name <- factor(g$term_name, rev(unique(g$term_name)))
    g$etiqueta <- factor(ifelse(nchar(as.character(g$term_name)) > 32, paste0(substr(g$term_name, 1, 31), "…"), as.character(g$term_name)),
                         rev(unique(ifelse(nchar(as.character(g$term_name)) > 32, paste0(substr(g$term_name, 1, 31), "…"), as.character(g$term_name)))))
    plot_ly(g, x = ~-log10(p_value), y = ~etiqueta, color = ~factor(source, fuentes), colors = COL_FUENTES[fuentes],
            type = "scatter", mode = "markers",
            marker = list(sizemode = "diameter", size = ~6 + 16 * sqrt(intersection_size / max(intersection_size)), line = list(width = 0), opacity = 0.85),
            text = ~sprintf("%s<br>%s<br>padj: %.2g<br>%s: %d", term_id, term_name, p_value, tr("go_tam"), intersection_size),
            hoverinfo = "text") |>
      # la leyenda de bases de datos va arriba, en HTML (go_dot_leyenda), para que quepa en una fila
      layout(title = list(text = titulo_lista(input$go_lista), font = list(size = 14)), showlegend = FALSE,
             xaxis = list(title = tr("go_x")), yaxis = list(title = "", automargin = TRUE))
  })
  # Heatmap de términos compartidos entre listas (como el heatmap de ComplexHeatmap de la lección):
  # filas = los N términos más significativos de cada lista; columnas = todas las listas del conjunto
  go_heat <- reactive({
    g <- obj(go_set()$objeto); g <- g[g$source == req(input$go_heat_fuente), ]
    validate(need(nrow(g) > 0, tr("go_vacio")))
    # columnas: grupo -> tipo celular -> UP/DOWN, como samples_order
    q <- unique(obj(go_set()$objeto)$query); p <- do.call(rbind, strsplit(q, "_"))
    ord <- order(match(p[, 1], c("Ctrl", "SLE", "SLEvsCtrl")), match(p[, 3], c("monocyte", "moDC", "tolDC")), match(p[, 2], c("UP", "DOWN")))
    cols <- data.frame(query = q[ord], grupo = p[ord, 1], dir = p[ord, 2], celula = p[ord, 3])
    top <- unique(unlist(lapply(split(g, g$query), function(x) head(x$term_name[order(x$p_value)], input$go_top))))
    m <- matrix(NA_integer_, length(top), nrow(cols), dimnames = list(top, cols$query))
    gt <- g[g$term_name %in% top, ]
    m[cbind(gt$term_name, gt$query)] <- gt$intersection_size
    # patrón de presencia (row_split): en qué listas aparece cada término
    patron <- apply(!is.na(m), 1, function(x) paste(ifelse(x, cols$dir, "NA"), collapse = "_"))
    m <- m[order(factor(patron, unique(patron[order(-rowSums(!is.na(m)), patron)])), apply(m, 1, max, na.rm = TRUE) * -1), , drop = FALSE]
    list(m = m, cols = cols, patron = patron[rownames(m)])
  })
  output$go_heat_leyenda <- renderUI({
    cols <- go_heat()$cols
    div(class = "d-flex flex-wrap column-gap-3 row-gap-1",
        leyenda(COLORES_DIR[c("UP", "DOWN")], tr("go_heat_expr")),
        leyenda(COL_CELULA[intersect(names(COL_CELULA), cols$celula)], tr("heat_pista_celula2")),
        if (all(cols$grupo == "SLEvsCtrl")) leyenda(COL_COMPARACION, tr("go_heat_comparacion"))
        else leyenda(c(Ctrl = "#96D4CC", SLE = "#B493B4")[unique(cols$grupo)], tr("heat_pista_grupo2")))
  })
  output$go_heat_ui <- renderUI(plotlyOutput("go_heat", height = sprintf("%dpx", 170 + 20 * max(nrow(go_heat()$m) + length(unique(go_heat()$patron)), 6))))
  output$go_heat <- renderPlotly({
    h <- go_heat(); m <- h$m; cols <- h$cols
    # posiciones con hueco entre grupos (column_split) y entre patrones (row_split)
    xs <- seq_len(ncol(m)) + cumsum(c(0, cols$grupo[-1] != cols$grupo[-nrow(cols)]))
    ys <- seq_len(nrow(m)) + cumsum(c(0, h$patron[-1] != h$patron[-length(h$patron)]))
    d <- data.frame(x = rep(xs, each = nrow(m)), y = rep(ys, ncol(m)), n = as.vector(m),
                    term = rep(rownames(m), ncol(m)), lista = rep(cols$query, each = nrow(m)))
    cmax <- max(100, d$n, na.rm = TRUE); cmin <- min(d$n, na.rm = TRUE)
    pos <- function(b) (b - cmin) / (cmax - cmin)
    escala <- list(c(0, "#FDE725"), c(pos(max(10, cmin + 1)), "#5DC863"), c(pos(50), "#3B528B"), c(pos(100), "#440154"), c(1, "#440154"))
    # cuadrícula completa: las columnas/filas de hueco quedan vacías (NA) y se ven como bandas grises
    X <- seq_len(max(xs)); Y <- seq_len(max(ys))
    fondo <- matrix(NA_real_, length(Y), length(X)); fondo[ys, xs] <- 0
    principal <- plot_ly() |>
      add_heatmap(x = X, y = Y, z = fondo, colorscale = list(c(0, "white"), c(1, "white")),
                  showscale = FALSE, xgap = 1, ygap = 1, hoverinfo = "skip") |>
      add_markers(data = d[!is.na(d$n), ], x = ~x, y = ~y, text = ~sprintf("%s<br>%s<br>%s: %d", term, lista, tr("go_tam"), n), hoverinfo = "text",
                  marker = list(size = 11, color = ~n, cmin = cmin, cmax = cmax, colorscale = escala, line = list(width = 0),
                                showscale = TRUE, colorbar = list(title = tr("go_heat_color"), len = 0.35, thickness = 12, y = 1, yanchor = "top")))
    pista <- function(v, colores, nombre) {
      z <- rep(NA_real_, length(X)); z[xs] <- match(v, names(colores))
      txt <- rep("", length(X)); txt[xs] <- paste(cols$query, v, sep = " · ")
      plot_ly(x = X, y = nombre, z = matrix(z, 1), type = "heatmap",
      colorscale = do.call(rbind, lapply(seq_along(colores), function(i) rbind(c((i - 1) / length(colores), colores[[i]]), c(i / length(colores), colores[[i]])))),
      zmin = 0.5, zmax = length(colores) + 0.5, showscale = FALSE, xgap = 1, text = matrix(txt, 1), hoverinfo = "text")
    }
    alto <- 20 * (max(ys) + 1)
    subplot(principal,
            pista(cols$dir, COLORES_DIR[c("UP", "DOWN")], tr("go_heat_expr")),
            pista(cols$celula, COL_CELULA, tr("heat_pista_celula2")),
            # en las comparaciones de contrastes, la última pista es la comparación (como en la lección 7)
            if (all(cols$grupo == "SLEvsCtrl")) pista(paste(cols$grupo, cols$celula, sep = "_"), COL_COMPARACION, tr("go_heat_comparacion"))
            else pista(cols$grupo, c(Ctrl = "#96D4CC", SLE = "#B493B4"), tr("heat_pista_grupo2")),
            nrows = 4, shareX = TRUE, margin = 0.003, heights = c(alto, 16, 16, 16) / (alto + 48)) |>
      layout(showlegend = FALSE, plot_bgcolor = "#D9D9D9",
             xaxis = list(showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE, range = c(0.5, max(xs) + 0.5)),
             yaxis = list(tickvals = ys, range = c(max(ys) + 0.5, 0.5), ticktext = ifelse(nchar(rownames(m)) > 36, paste0(substr(rownames(m), 1, 35), "…"), rownames(m)),
                          showgrid = FALSE, zeroline = FALSE, tickfont = list(size = 10)))
  })

  output$go_tabla <- renderDT(
    datatable(go_filtrado()[, c("source", "term_id", "term_name", "p_value", "term_size", "intersection_size", "intersection")],
              rownames = FALSE, filter = "top", fillContainer = FALSE,
              options = list(pageLength = 10, lengthMenu = c(10, 25, 50, 100), scrollX = TRUE)) |>
      formatSignif("p_value", 3)
  )
  output$dl_go <- downloadHandler(
    filename = function() sprintf("GeneExprLupusMx_%s_%s.csv", go_set()$objeto, input$go_lista),
    content = function(file) escribir_csv(go_filtrado(), file)
  )

  # Descargas
  dl_x <- reactive(tabla_descarga(input$dl_obj))
  output$dl_info <- renderUI({
    obj_name <- sub("\\$.*", "", input$dl_obj)
    tagList(p(class = "mt-2 small", strong(tr("dl_dim")), ": ", paste(format(dim(dl_x()), big.mark = ","), collapse = " × ")),
            a(href = sprintf("%s/reference/%s.html", URL_SITE, obj_name), target = "_blank", icon("circle-info"), sprintf(" ?%s", obj_name)))
  })
  output$dl_vista <- renderTable({
    x <- dl_x(); x <- x[seq_len(min(8, nrow(x))), seq_len(min(8, ncol(x))), drop = FALSE]
    if (is.matrix(x)) x <- cbind(gene = rownames(x), as.data.frame(x))
    x[] <- lapply(x, function(c) if (inherits(c, "Date")) format(c) else c)
    x
  }, digits = 3)
  output$dl_csv <- downloadHandler(
    filename = function() sprintf("GeneExprLupusMx_%s.csv", sub("$", "_", input$dl_obj, fixed = TRUE)),
    content = function(file) escribir_csv(dl_x(), file)
  )
}

shinyApp(ui, server)
