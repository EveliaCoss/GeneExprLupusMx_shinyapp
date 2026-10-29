# Precalcula los PCA de la lección 2 ("De las cuentas de RNA-seq a la exploración de datos"):
# vst(dds, blind = FALSE) y prcomp() sin escalar, como calc_pca() de la lección.
# Uno por cada combinación de tipos celulares y por cada número de genes: todos (como en la lección)
# o los 500, 250 y 100 más variables dentro de esas muestras (como DESeq2::plotPCA()).
# Uso (desde la raíz de la app): Rscript tools/pca.R [ruta_al_paquete]
suppressMessages(library(DESeq2))
pkg <- commandArgs(TRUE)[1]
if (is.na(pkg)) pkg <- "~/Documents/GeneExprLupusMx"
load(file.path(pkg, "data", "dds.rda"))

vsd <- assay(vst(dds, blind = FALSE))
tipos <- c("monocyte", "moDC", "tolDC")
combos <- unlist(lapply(1:3, function(k) combn(tipos, k, simplify = FALSE)), recursive = FALSE)
n_genes <- c(todos = Inf, "500" = 500, "250" = 250, "100" = 100)
pca <- list()
for (ct in combos) for (ng in names(n_genes)) {
  m <- vsd[, dds$Cell_type %in% ct]
  rv <- apply(m, 1, var)
  m <- m[order(rv, decreasing = TRUE)[seq_len(min(n_genes[[ng]], nrow(m)))], ]
  p <- prcomp(t(m))
  pca[[paste(paste(ct, collapse = "+"), ng, sep = "|")]] <- list(
    x = p$x, variance = round(100 * p$sdev^2 / sum(p$sdev^2), 2))  # % de varianza por PC
}
saveRDS(pca, "datos/pca.rds")
