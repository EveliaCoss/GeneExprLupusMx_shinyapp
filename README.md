# GeneExprLupusMx · Explorador (Shiny)

Aplicación Shiny para explorar de forma interactiva los resultados del artículo
[Hernández-Ledesma, Coss-Navarrete *et al.* (2026) *Front. Immunol.* 17:1771959](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2026.1771959/full),
usando los datos del paquete de R [`GeneExprLupusMx`](https://github.com/EveliaCoss/GeneExprLupusMx).
El aspecto (colores, tipografía, figuras) sigue al [sitio del paquete](https://eveliacoss.github.io/GeneExprLupusMx/);
las infografías se cargan desde ese sitio, así que se necesita conexión a internet para verlas.

Secciones (en orden progresivo): Inicio · Participantes · Exploración (PCA y expresión por gen) ·
Expresión diferencial (9 comparaciones, volcano, tabla, descargas) · Funciones biológicas (g:Profiler) ·
Descargas (cualquier tabla del paquete en CSV) · Aprende y cita (enlaces a las 7 lecciones).

## Ejecutar

```r
install.packages(c("shiny", "bslib", "DT", "plotly"))
shiny::runApp()
```

### ¿De dónde lee los datos?

1. Si existe la variable de entorno `GENEEXPRLUPUSMX_DATA`, carga los `.rda` de esa carpeta.
2. Si no, usa el paquete instalado (`remotes::install_github("EveliaCoss/GeneExprLupusMx")`).
3. Si no está instalado, usa `~/Documents/GeneExprLupusMx/data` (copia local en desarrollo).

Los objetos pesados (`filtered_txi`, `zscore_mat`) solo se cargan cuando se necesitan.

### PCA

`datos/pca.rds` (≈800 KB) guarda los PCA de la lección 2 (`vst(dds, blind = FALSE)` + `prcomp()`, sin escalar)
para cada combinación de tipos celulares y con todos los genes o los 500, 250 y 100 más variables,
así la app no necesita DESeq2. Para regenerarlo:
`Rscript tools/pca.R ~/Documents/GeneExprLupusMx` (requiere DESeq2).

### Imágenes

Las imágenes de `www/` ("Datos incluidos en el paquete" con las células, "Diferenciación in vitro" y miniaturas de las secciones) son recortes de
las figuras del paquete. Para regenerarlas tras actualizar el paquete: `Rscript tools/imagenes.R ~/Documents/GeneExprLupusMx`.

## Idiomas

Todos los textos están en `textos.R` (`TXT$es`). Para agregar inglés, copia `TXT$es` como `TXT$en`,
tradúcelo y ejecuta con `Sys.setenv(APP_LANG = "en")`. Los enlaces a las lecciones ya cambian
automáticamente a la versión `_en`.

## Cómo citar

Cita el **paquete** y los **datos**:

- Coss E (2026) *GeneExprLupusMx: Datos transcriptómicos de lupus en mujeres mexicanas con lecciones de aprendizaje*. Paquete de R. <https://github.com/EveliaCoss/GeneExprLupusMx> (DOI de Zenodo pendiente).
- Hernández-Ledesma AL, Coss-Navarrete EL, *et al* (2026) Gene expression profiling of dendritic cell tolerance dysfunction in women with systemic lupus erythematosus. *Front. Immunol.* 17:1771959. doi: [10.3389/fimmu.2026.1771959](https://doi.org/10.3389/fimmu.2026.1771959).
- Coss-Navarrete EL, Hernández-Ledesma AL, *et al* (2025) Analysis of gene expression profiles during tolerance induction in dendritic cells from women with Systemic Lupus Erythematosus [Data set]. *Zenodo*. doi: [10.5281/zenodo.17742206](https://doi.org/10.5281/zenodo.17742206).

## Verificación

```sh
Rscript tests/check.R
```
