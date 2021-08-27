# BioFreelancer 2021
# Creado para la UGPM LaNSE CINVESTAV

# Cargamos librerias base ====
library("dplyr")
library("stringr")
library("matrixStats")
library("openxlsx")

# Limpieza general de filas. Eliminar infinitos y asi ====
data_toda <- read.csv( file = archivo_proteinas )

# separamos los onlyids
only_ids <- data_toda %>% 
  filter( Unique.peptides == 0 )

# separamos proteinas exclusivas
exclusivas <- data_toda %>%
  filter( Max.fold.change == Inf )

# encuentra las celdas que se llaman REVERSE
indice_reverse <- str_detect( string = data_toda$Accession,
                              pattern = "REVERSE")

# nomas por morbo sacamos las reverse
reverse <- data_toda %>% 
  filter( indice_reverse )

# Sacamos las chidas
cuantificadas <- data_toda %>% 
  filter( Unique.peptides != 0,
          Max.fold.change != Inf, 
          !indice_reverse )

# Calcular valores por Condicion (mean, cv, etc ) ====
# Leemos la hoja de muestras
muestras <- read.csv( file = sample_sheet )

# encontrar las condiciones posibles
condiciones <- muestras %>%
  pull( condition ) %>% 
  unique()

# Revisamos las condiciones
condiciones

# Sacamos los nombres de las condiciones...
condicion1 <- condiciones[1]
condicion2 <- condiciones[2]

# Revisamos objetos
condicion1
condicion2

# Calculamos para condicion 1 ====

# Sacamos los identificadores de las muestras
ids_condicion1 <- muestras %>% 
  filter( condition == condicion1 ) %>% 
  pull( muestra )

# Revisamos
ids_condicion1

# creamos una expresion regular para buscar todas las muestras de la condicion 1
expresion_1 <- paste( ids_condicion1,
                      collapse = "|" )

# Revisamos
expresion_1

# Vamos a buscar el numero de columnas que corresponden a las muestras de la condicion 1
posicion_muestras_1 <- str_detect( string = colnames( cuantificadas ),
                                   pattern = expresion_1 ) %>%
  which()

# Revisamos
posicion_muestras_1

# vamos a extraer un dataframe que contenga solo las columnas de las inyecciones
data_condicion_1 <- cuantificadas %>% 
  select( posicion_muestras_1 )

# Calculamos valores de interes para las inyecciones de la condicion 1
condicion1_final <- data_condicion_1 %>% 
  mutate( average.cond1 = rowMeans( . ) ) %>% 
  mutate( log10.cond1 = log10( average.cond1 ) ) %>%
  mutate( sd.cond1 = as.matrix( data_condicion_1 ) %>% rowSds( ) ) %>% 
  mutate( CV.cond1 = sd.cond1 / average.cond1 ) %>% 
  mutate( Count.cond1 = rowSums( data_condicion_1 > 0 ) )

# Calculamos para condicion 2 ====
# Sacamos los identificadores de las muestras para la condicion 2
ids_condicion2 <- muestras %>%
  filter( condition == condicion2 ) %>%
  pull( muestra )

# Revisamos
ids_condicion2

# creamos una expresion regular para buscar todas las muestras de la condicion 2
expresion_2 <- paste( ids_condicion2,
                      collapse = "|" )

# Revisamos
expresion_2

# Vamos a buscar el numero de columnas que corresponden a las muestras de la condicion 1
posicion_muestras_2 <- str_detect( string = colnames( cuantificadas ),
                                   pattern = expresion_2 ) %>%
  which()

# Revisamos
posicion_muestras_2

# hacemos lo mismo pero para la condicion 2
data_condicion_2 <- cuantificadas %>% 
  select( posicion_muestras_2 )

# Calculamos valores de interes para las inyecciones de la condicion 1
condicion2_final <- data_condicion_2 %>% 
  mutate( average.cond2 = rowMeans( . ) ) %>% 
  mutate( log10.cond2 = log10( average.cond2) ) %>%
  mutate( sd.cond2 = as.matrix( data_condicion_2 ) %>% rowSds( ) ) %>% 
  mutate( CV.cond2 = sd.cond2 / average.cond2 ) %>% 
  mutate( Count.cond2 = rowSums( data_condicion_2 > 0 ) )

# juntamos toda la info calculada ====
# preparamos las columnas base de la data original
base_data <- cuantificadas %>% 
  select( Accession, Peptide.count, Unique.peptides,
          Anova..p.,
          Max.fold.change,
          Description )

# pegaremos lado a lado los dataframes
# Esto funciona bien porque NO reordenamos ningun DF previamente
juntado_data <- bind_cols( base_data,
                           condicion1_final,
                           condicion2_final )

# Nos tardamos pero... chulada :) 

# filtramos por el n de muestras con datos en las inyecciones ====
# cual es el maximo de muestras en la condicion 1
n_muestras_1 <- muestras %>% 
  filter( condition == condicion1 ) %>% 
  nrow( )

# cual es el maximo de muestras en la condicion 2
n_muestras_2 <- muestras %>% 
  filter( condition == condicion2 ) %>% 
  nrow( )

# eliminar counts < al numero de muestras
filtrado_data <- juntado_data %>% 
  filter( Count.cond1 == n_muestras_1,
          Count.cond2 == n_muestras_2 )

# guardamos cuales se perdieron por counts
sacado_counts <- filtrado_data <- juntado_data %>% 
  filter( Count.cond1 == 0 |
          Count.cond2 == 0 )

# Ultimos calculos ====
final_data <- juntado_data %>% 
  mutate( nlog10.Anova = -log10( Anova..p. ) ) %>% 
  mutate( ratio = average.cond2 / average.cond1 ) %>%
  mutate( log2.Ratio = log2( ratio )  )

# arrange y ordenar ids por abundancia en control
final_data_a <- final_data %>% 
  arrange( desc( log10.cond1 ) ) %>% 
  mutate( id.cond1 = 1:nrow( . ) )

# arrange y ordenar ids por abundancia en casos
final_data_b <- final_data_a %>% 
  arrange( desc( log10.cond2 ) ) %>% 
  mutate( id.cond2 = 1:nrow( . ) )

# exportar la data a un archivo excel...
# guardamos en una lista para el excel final
# mas info aqui: https://stackoverflow.com/questions/27713310/easy-way-to-export-multiple-data-frame-to-multiple-excel-worksheets
# la sintaxis es: list( "nombre de hoja" = dataframe )
para_excel <- list( "todos" = data_toda,
                    "only ids" = only_ids,
                    "exclusivas" = exclusivas,
                    "reverse" = reverse,
                    "cuantificadas" = cuantificadas,
                    "sacadas counts == 0" = sacado_counts,
                    "final para volcano" = final_data_b )

# guardamos
write.xlsx( x = para_excel,
            file = paste( ruta_salida, "resultado_proteinas.xlsx", sep = "/" ),
            overwrite = TRUE )

## FIN DEL SCRIPT