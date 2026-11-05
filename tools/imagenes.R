# Genera las imágenes de www/ a partir de las figuras del paquete GeneExprLupusMx.
# Uso (desde la raíz de la app): Rscript tools/imagenes.R [ruta_al_paquete]
library(magick)
pkg <- commandArgs(TRUE)[1]
if (is.na(pkg)) pkg <- "~/Documents/GeneExprLupusMx"
leer <- function(f) image_read(file.path(pkg, f))
recortar <- function(img, x0, y0, x1, y1) image_crop(img, sprintf("%dx%d+%d+%d", x1 - x0, y1 - y0, x0, y0))
miniatura <- function(img, f) {
  img <- image_scale(image_trim(img), "560x300")
  image_write(image_extent(image_background(img, "white"), "560x300", color = "white"), file.path("www", f))
}

promo <- leer("img/overview_promo.png")      # 4800 x 2884
l1 <- leer("vignettes/articles/img/lesson1_overview_es.png")
l4 <- leer("vignettes/articles/img/lesson4_overview_es.png")

# Células de la app (img/): monocito, moDC y tolDC. Reemplazan a las del póster del paquete.
celulas <- lapply(c("monocito", "modc", "toldc"), function(x) image_trim(image_read(sprintf("img/celula_%s.png", x))))

# Tapa la célula vieja dentro de la caja c(x0, y0, x1, y1) y pega la nueva centrada. El fondo se rehace
# con la columna a la izquierda y la fila arriba de la caja (blanco, o las líneas de la cuadrícula del PCA).
cambiar_celula <- function(img, caja, cel) {
  x0 <- caja[1]; y0 <- caja[2]; w <- caja[3] - x0; h <- caja[4] - y0
  estirar <- function(x) image_resize(x, sprintf("%dx%d!", w, h), filter = "Point")
  fondo <- image_composite(estirar(recortar(img, x0 - 1, y0, x0, y0 + h)), estirar(recortar(img, x0, y0 - 1, x0 + w, y0)), "darken")
  fondo <- image_composite(fondo, image_scale(cel, sprintf("%dx%d", w, h)), "over", gravity = "center")
  image_composite(img, fondo, offset = sprintf("+%d+%d", x0, y0))
}
cajas <- rbind(                         # en overview_promo.png: x0, y0, x1, y1, célula
  c(105, 2060, 330, 2285, 1), c(770, 1965, 965, 2170, 2), c(770, 2172, 965, 2370, 3),    # Diferenciación
  c(1512, 862, 1603, 953, 1), c(1740, 795, 1875, 925, 2), c(1707, 1172, 1828, 1292, 3),  # PCA
  c(2615, 809, 2705, 898, 1), c(2945, 802, 3053, 915, 2), c(3282, 797, 3388, 916, 3))    # volcanos
for (i in seq_len(nrow(cajas))) promo <- cambiar_celula(promo, cajas[i, 1:4], celulas[[cajas[i, 5]]])

# Paso 1: sección "Datos incluidos en el paquete" con las células debajo de su columna.
# Se parte el panel bajo la tabla (y = 620) y se inserta una franja (la fila 620 estirada,
# que conserva los bordes laterales) con monocito, moDC y tolDC centrados en su columna.
datos <- recortar(promo, 55, 655, 1185, 1312)                      # 1130 x 657
centros <- c(480, 762, 982)                                        # centro de Monocito, moDC y tolDC
franja <- image_scale(recortar(datos, 0, 620, 1130, 621), "1130x170!")
for (i in 1:3) {
  alto <- c(110, 140, 140)[i]                                      # el monocito es más pequeño
  ci <- image_scale(celulas[[i]], paste0("x", alto))
  franja <- image_composite(franja, ci, offset = sprintf("+%d+%d", centros[i] - image_info(ci)$width %/% 2, 15 + (140 - alto) %/% 2))
}
datos <- image_append(c(recortar(datos, 0, 0, 1130, 620), franja, recortar(datos, 0, 620, 1130, 657)), stack = TRUE)
image_write(image_scale(datos, "800x"), "www/datos_paquete.png")

# Paso 1: sección "Diferenciación in vitro"
image_write(image_scale(recortar(promo, 55, 1840, 1185, 2560), "800x"), "www/diferenciacion.png")

# Miniaturas de las secciones de la app
miniatura(recortar(l1, 1788, 672, 2364, 948), "sec_participantes.png")
miniatura(recortar(promo, 1248, 763, 2352, 1440), "sec_exploracion.png")
miniatura(recortar(promo, 2436, 763, 3564, 1428), "sec_deg.png")
miniatura(recortar(l4, 564, 238, 1134, 552), "sec_go.png")
miniatura(recortar(promo, 72, 1464, 1164, 1788), "sec_descargas.png")
miniatura(recortar(promo, 3636, 1728, 4728, 2196), "sec_aprende.png")
