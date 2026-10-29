# Verificación rápida de la lógica del servidor: Rscript tests/check.R (desde la raíz de la app)
app <- source("app.R")$value
shiny::testServer(app, {
  # Conteos de DEG documentados en el paquete (|LFC| >= 1, padj < 0.05; interacción solo padj)
  esperado <- list(monocyte = c(37, 0), moDC = c(4, 40), tolDC = c(4, 12),
                   Ctrl_moDC = c(2688, 2947), SLE_tolDC = c(2804, 3596),
                   int_moDC = c(1, 28), int_tolDC = c(3, 38))
  for (i in names(esperado)) {
    f <- CMP$familia[CMP$id == i]
    session$setInputs(deg_familia = f, deg_cmp = i, deg_lfc = if (f == "int") 0 else 1, deg_padj = "0.05", deg_todos = FALSE)
    n <- c(sum(deg()$direccion == "UP"), sum(deg()$direccion == "DOWN"))
    stopifnot(identical(n, as.integer(esperado[[i]])), nrow(deg_mostrar()) == sum(n))
    session$setInputs(deg_etiquetas = TRUE, deg_top = 30, deg_genes = c("IFI27", "CD14")); invisible(output$deg_volcano)
  }
  stopifnot(nrow(read.csv(output$dl_deg)) == 41, nrow(read.csv(output$dl_cmp)) == 21718)
  # Proporción UP/DOWN y heatmap de log2FC (diferenciación: 4 columnas, 50 genes)
  session$setInputs(deg_familia = "dif", deg_cmp = "SLE_tolDC", deg_lfc = 1, deg_padj = "0.05", lfch_n = "50", lfch_bio = TRUE)
  stopifnot(nrow(familia_info()) == 4, ncol(lfch()$m) == 4, nrow(lfch()$m) == 50)
  invisible(output$deg_prop); invisible(output$deg_lfch); invisible(output$deg_lfch_interp)
  # Biotipos: los porcentajes de cada grupo suman 100 y los DEG coinciden con el conteo de la comparación
  b <- bio(); stopifnot(all(abs(tapply(b$pct, b$grupo, sum) - 100) < 1e-6), sum(b$n[b$grupo != "fondo"]) == n_sig())
  invisible(output$deg_bio); invisible(output$deg_bio_interp)
  # PCA como en la lección 2 (PC1 63.44 %, PC2 9.51 %)
  session$setInputs(pca_celulas = c("monocyte", "moDC", "tolDC"), pca_ngenes = "todos", pca_x = "1", pca_y = "3", pca_color = "Cell_type")
  stopifnot(nrow(pca()$x) == 99, identical(head(pca()$var, 2), c(63.44, 9.51)))
  invisible(output$pca_plot); invisible(output$scree_plot)
  session$setInputs(pca_ngenes = "500"); stopifnot(pca()$var[1] != 63.44)
  session$setInputs(pca_celulas = "tolDC", pca_ngenes = "100"); stopifnot(nrow(pca()$x) == 33)
  # Heatmap: los 10 ISG de la lección y genes elegidos a mano
  session$setInputs(heat_set = "isg10", heat_cluster = TRUE, heat_por = "grupo", heat_pistas = c("celula", "grupo"))
  stopifnot(nrow(heat_mat()) == 10, ncol(heat_mat()) == 99); invisible(output$heat_plot)
  session$setInputs(heat_set = "custom", heat_genes = c("CD14", "CD1C", "IFI27"), heat_cluster = FALSE)
  stopifnot(identical(rownames(heat_mat()), c("CD14", "CD1C", "IFI27")))
  stopifnot(nrow(read.csv(output$dl_heat)) == 3)
  # Agrupaciones: sin NA y con las controles como "Ctrl"
  m <- obj("metadata_cleaned")
  for (v in c("grupo", "edad", "actividad", "gc")) stopifnot(!anyNA(agrupar(m, v)))
  stopifnot(all(agrupar(m, "gc")[m$Group == "Ctrl"] == "Ctrl"), sum(agrupar(m, "actividad") == "Ctrl") == 30)
  # Heatmap con barras de edad y actividad, separado por glucocorticoides
  session$setInputs(heat_set = "isg10", heat_por = "gc", heat_pistas = c("celula", "edad", "actividad")); invisible(output$heat_plot)
  session$setInputs(heat_pistas = character(0)); invisible(output$heat_plot)
  # Dot plot de biomarcadores: 19 genes × 6 filas (tipo celular × grupo), como la Figura 2B
  session$setInputs(dot_set = "marcadores", dot_por = "grupo")
  stopifnot(nrow(dot_df()) == 19 * 6); invisible(output$dot_plot)
  session$setInputs(dot_por = "actividad"); stopifnot(length(unique(dot_df()$cond)) == 12)
  session$setInputs(gen = "IFI27", gen_por = "edad"); invisible(output$gen_plot)
  session$setInputs(go_set = "sle", go_lista = "SLE_UP_monocyte", go_fuente = "GO:BP", go_top = 10)
  stopifnot(nrow(go_filtrado()) > 0, all(go_filtrado()$source == "GO:BP"))
  # Reactome visible: con "por base de datos" el dot plot incluye términos de cada fuente elegida
  session$setInputs(go_set = "sle", go_lista = "SLE_UP_monocyte", go_fuente = c("GO:BP", "REAC"), go_modo = "fuente", go_top = 5)
  stopifnot(setequal(unique(go_dot()$source), c("GO:BP", "REAC"))); invisible(output$go_plot); invisible(output$go_manhattan)
  # Heatmap de términos compartidos: 8 listas en diferenciación, ordenadas Ctrl -> SLE
  session$setInputs(go_set = "dif", go_lista = "Ctrl_UP_moDC", go_heat_fuente = "REAC", go_top = 5)
  stopifnot(ncol(go_heat()$m) == 8, colnames(go_heat()$m)[1] == "Ctrl_UP_moDC"); invisible(output$go_heat)
  session$setInputs(dl_obj = "volunteer_info"); stopifnot(nrow(read.csv(output$dl_csv)) == 33)
})
cat("OK\n")
