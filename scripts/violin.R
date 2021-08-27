# cargamos libs
library( "dplyr" )
library( "ggplot2" )
library( "tidyr" )
library( "stringr" )
library( "ggsci" )
library( "scales" )

# Vamos a construir un boxplot basico
# Necesitamos datos de accession y las mediciones
data_box <- base_data %>% 
  select( Accession ) %>% 
  bind_cols( data_condicion_1, data_condicion_2 )

# pasamos a formato long
largo1 <- pivot_longer( data = data_box,
                        cols = -Accession,
                        names_to = "Xmuestra",
                        values_to = "medida" )

# vamos a quitarle la espantosa X a los nombres de muestras
largo2 <- largo1 %>% 
  mutate( muestra = str_sub( Xmuestra, start = 2)  )

# unimos con la hoja de muestras
largo3 <- left_join( x = largo2,
                     y = muestras,
                     by = "muestra" )

# Todo chido para hacer el boxplot
# metemos todo en una funcion que hace las cajas
graficadora <- function( la_data, el_nombre_del_set ){
  
  ggplot( data = la_data,
          mapping = aes( x = condition ,
                         y = medida,
                         fill = condition) ) +
    geom_violin( alpha = 0.7 ) +
    geom_boxplot( width = 0.05, fill = "white", alpha = 1 ) +
    scale_y_continuous( labels = comma ) +
    scale_fill_d3() +
    labs( title = "Peptidoma Diferencial",
          subtitle = el_nombre_del_set,
          x = "",
          y = "intensidad" ) +
    theme_bw() +
    theme( legend.position = "none",
           text = element_text( size = 15) )
  
}

# probamos
todos_box <- graficadora( la_data = largo3, el_nombre_del_set = "All proteins" )

# Guardamos el plot
guardar_imagen( nombre_imagen = "boxplot_todos_los_datos.png",
                el_plot = todos_box, ancho = 5, alto = 5 )

# creamos un dataset con Accession de interes
# cargamos accession de interes
de_interes <- read.csv( file = set_interes )

# veamos en que posicion estan los accession de interes
todos_accession <- largo3 %>% pull( Accession )

interes_accession <- de_interes %>% pull( Accession )

# Cremoas un vector con las posiciones de interes
posiciones_interes <- todos_accession %in% interes_accession

# usamos las posiciones de interes para filtrar filas
largo_interes <- largo3 %>% 
  filter( posiciones_interes )

# probamos graficar los de interes
boxinteres <- graficadora( la_data = largo_interes,
                      el_nombre_del_set = subtitulo_interes )

# Vis
boxinteres

# Guardamos el plot
guardar_imagen( nombre_imagen = "boxplot_de_interes.png",
                el_plot = boxinteres, ancho = 5, alto = 5 )

# FIN DEL SCRIPT