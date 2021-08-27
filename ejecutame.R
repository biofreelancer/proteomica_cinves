# BioFreelancer 2021
# Script para procesar datos resultado de analisis proteomicos de la Unidad de Genómica, Proteómica y Metabolómica, CINVESTAV-IPN
# Si quieres un script como este para automatizar tus analisis, solicitalo a;
# cursos.biofreelancer@gmail.com

# Lista de paquetes necesarios
# library("purrr")
# library("dplyr")
# library("ggplot2")
# library("scales")
# library("ggsci")
# library("stringr")
# library("matrixStats")
# library("openxlsx")
# library("ggrepel")
# library("cowplot")
# library("tidyr")
# library("scales")
# library("factoextra")
# library("pheatmap")

library("dplyr")

# DECLARACION DE PARAMETROS ====
# GENERALES
# Archivos de entrada
sample_sheet <- "data/sample_sheet.csv"               # Ejemplo "data/sample_sheet.csv"  
archivo_proteinas <- "data/dataset_proteomica.csv"    # "data/dataset_proteomica.csv"

# dale un set de ids para hacer boxplot enfocado a ellos
set_interes <- "data/set1.csv"           # Archivo que contiene el set de Accession a destacar en el grafico de caja y violin
subtitulo_interes <- "Proteinas UP"      # El nombre que aparecera como subtitulo en el boxplot-violin de interes

carpeta_peptidos <- "data/archivos_peptidos/"         # "data/archivos_peptidos/"

# PCA
pc_enx <- 1   # que componente quieres graficar en el eje X
pc_eny <- 2   # que componente quieres graficar en el eje Y

# FILTRADO GENERAL PARA ARCHIVO_PROTEINAS
# limite pvalue < X       # para definir lo UP y DOWN, del volcano, heatmap, etc
pvalue_limite <- 0.05 

# limite o threshold del fold change para destacar proteinas (tanto up, como down)
# limite ratio > log2(X) para lo up
# limite ratio < log2(X) para lo down
ratio_limite <- 2 %>%               # un fold change de dos significa un cambio de el doble de expresion
  log2()

# limte de CV > X
cv_limite <- 0

# limite peptide count >= X
peptide_count_limite <- 2
# limite unique peptide >= X
unique_peptide_limite <- 1
# limite de counts en cond1 y cond2
# limite >= X
count_limite <- 2

# para el VOLCANO etiquetamos los nombres de los top N ====
n_best <- 3    # numero de proteinas top up y down a resaltar con nombre

# FILTRADO PARA CARPETA_PEPTIDOS
# sera el +- error ppm a resaltar en el histograma
corte_error <- 5

# GENERALES FLUJO DE DATOS
# crear la carpeta de salida por si no existe
ruta_salida <- "results/"               # asi se va a llamar el directorio donde todo aparecera
# crear el directorio de salida, si ya existe no se queja ( showWarnings = FALSE )
dir.create( path = ruta_salida,
            showWarnings = FALSE )

# DEFINICION DE FUNCIONES QUE SE USARAN MUCHO
guardar_imagen <- function( nombre_imagen, el_plot, ancho, alto ) {
  
  ggsave( filename = paste(ruta_salida, nombre_imagen, sep = "/"),
          plot = el_plot,
          width = ancho,
          height = alto,
          dpi = 300 )
  
}

# EJECUCION DE SCRIPTS ====

# llamamos al primer script que hace TODO el manejo de datos
source( file = "scripts/dataclean.R" )

# graficamos el rango dinamico
source( file = "scripts/rangodinam.R" )

# graficamos el volcan
source( file = "scripts/volcan.R" )

# graficamos el violin
source( file = "scripts/violin.R" )

# graficamos el PCa y heatmap
source( file = "scripts/heat_pca.R" )

# Operamos sobre los archivos individuales de peptios
# para hacer histogramas, dispersion y pastel
source( file = "scripts/graficos_peptido.R" )

# FIN DEL CURSO!! GRACIAS POR SU CONFIANZA
# BioFreelancer 2021
# www.biofreelancer.com
