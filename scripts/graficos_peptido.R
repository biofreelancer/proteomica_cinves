# Cargamos librerias
library("purrr")
library("dplyr")
library("ggplot2")
library("scales")
library("ggsci")

# Primero encontrar las rutas de todos los archivos que se necesitan
rutas_archivos.v <- list.files( path = carpeta_peptidos,                       # buscar dentro del directorio de trabajo
                                pattern = ".csv" )                               # queremos la ruta completa. P. ej. "criptos_cierre_mensual/BTC.csv"

# crear una funcion que lee un df y pone el nombre del archivo
leer_peptidos <- function( el_archivo ){
  
  read.csv( file = paste0(carpeta_peptidos, el_archivo) ) %>% 
    mutate( injection = el_archivo )
  
}

# Cargamos todos los archivos
todas_las_cargas.df <- rutas_archivos.v %>%                                   # a partir de todos los archivos encontrados (de sus rutas, en realidad)
  map_df( ~ leer_peptidos( el_archivo = . ) ) %>% 
  mutate( charge = case_when( precursor.z == 1 ~ "1",
                              precursor.z == 2 ~ "2",
                              precursor.z > 2 ~ "3+" ) )

# y si hacemos un histograma?
histo1 <- ggplot( data = todas_las_cargas.df,
                  mapping = aes( x = peptidePrecursor.deltaMhpPPM ) ) +
  geom_histogram( binwidth = 1,
                  fill = "cornflowerblue",
                  alpha = 0.5 )

# vis
histo1

# Vamos a resaltar en el histograma algunos bins
filtrado_error <- todas_las_cargas.df %>% 
  filter( peptidePrecursor.deltaMhpPPM >= -corte_error,
          peptidePrecursor.deltaMhpPPM <= corte_error   )

# resaltamos algunos bins
histo2 <- histo1 +
  geom_histogram( data = filtrado_error,
                  binwidth = 1,
                  fill = "cornflowerblue",
                  color = "black" )

# Vis
histo2

# areglemos el plot
minx <- min( todas_las_cargas.df$peptidePrecursor.deltaMhpPPM,
             na.rm = TRUE )
maxx <- max( todas_las_cargas.df$peptidePrecursor.deltaMhpPPM,
             na.rm = TRUE )
# redondeamos los numeros de los ejes
marcas_x <- seq( from = minx,
                 to = maxx,
                 by = 5 ) %>%
  round()

# mejoramos eje, agregamos titulos, y aligeramos tema
histo3 <- histo2 +
  scale_x_continuous( breaks = marcas_x ) +
  labs( x = "ppm Error",
        y = "Raw counts" ) +
  theme_light( base_size = 15 )

# vis
histo3

# guardamos plot
guardar_imagen( nombre_imagen = "histograma_completo.png", el_plot = histo3, ancho = 7, alto = 5 )

# a ver con facets...
histo4 <- histo3 + facet_wrap( ~ injection ) +
  theme( strip.text.x = element_text(size = 10) )

# Vis
histo4

# guardamos plot
guardar_imagen( nombre_imagen = "histo_muestras.png",
                el_plot = histo4, ancho = 10, alto = 5 )


# Para un pie chart de peptide.matchType
# precalculamos la tabla con los valores para el pastel
pastel_d1 <- todas_las_cargas.df %>% 
  group_by( peptide.matchType ) %>% 
  summarize( counts = n( ) ) %>% 
  arrange( -counts )

# calculamos porcentajes
pastel_d2 <- pastel_d1 %>% 
  mutate( proporcion = counts / sum(counts),
          porcentaje = percent( proporcion, accuracy = 0.1 ),
          etiqueta = paste( porcentaje, peptide.matchType ) ) %>% 
  mutate( etiqueta = factor( etiqueta,
                             levels = etiqueta ) )

# graficamos
pastel1 <- ggplot( data = pastel_d2,
                   mapping = aes( x = 1,
                                  y = proporcion,
                                  fill = etiqueta ) ) +
  geom_col( width = 1, color = "black" )

# vis
pastel1

# cambiamos colores y polarizamos
pastel2 <-  pastel1 +
  scale_fill_npg( ) +
  coord_polar( theta = "y" )

# Vis
pastel2

# agregamos titulos y limpiamos tem
pastel3 <- pastel2 +
  labs( fill = "Match Type" ) +
  theme_void( base_size = 15 )

# Vis
pastel3

# guardamos plot
guardar_imagen( nombre_imagen = "pastel.png",
                el_plot = pastel3, ancho = 5, alto = 5 )

# hagamos un scater de movimiento (se ve chido)
# Nota: quiza se veria bien como un histograma 2d
puntos_mz <- ggplot( data = todas_las_cargas.df,
                     mapping = aes( x = precursor.mz,
                                    y = peptidePrecursor.deltaMhpPPM,
                                    color = peptide.matchType ) ) +
  geom_point( alpha = 0.3 ) +
  labs( y = "Error PPM",
        x = "m/z") +
  scale_color_npg( ) +
  theme_classic( base_size = 15 )

# vis
puntos_mz

# guardamos plot
guardar_imagen( nombre_imagen = "puntosmz.png",
                el_plot = puntos_mz, ancho = 10, alto = 5 )

# probemos el plot de movilidad
movilidad <- ggplot( data = todas_las_cargas.df,
                     mapping = aes( x = precursor.Mobility,
                                    y = precursor.mz,
                                    color = charge ) ) +
  geom_point( size = 0.5, alpha = 0.5 ) +
  theme_classic( base_size = 15 )

# Vis
movilidad

# guardamos plot
guardar_imagen( nombre_imagen = "movilidad.png",
                el_plot = movilidad, ancho = 7, alto = 5 )

# Fin del Script